// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include "PlanarProfile.hpp"
#include "MeshDiagnostics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <utility>

#include <earcut/earcut.hpp>

namespace Slic3r {
namespace {

// Picking is interactive and this API has no cancellation callback. Bound the
// quadratic contour validation rather than blocking indefinitely on hostile data.
constexpr size_t max_profile_vertices = 4096;
struct Loop { size_t begin; size_t end; };
using Edge = std::pair<int, int>;

PlanarProfileResult failure(FaceExtrusionError error, const char *message)
{
    PlanarProfileResult result;
    result.error = error;
    result.message = message;
    return result;
}

long double cross2(const Vec2d &a, const Vec2d &b)
{
    return static_cast<long double>(a.x()) * b.y() - static_cast<long double>(a.y()) * b.x();
}

long double area2(const std::vector<Vec2d> &points, const Loop &loop)
{
    // Shift separately for each hole: a small hole far from the profile origin
    // must not disappear through cancellation of large shoelace terms.
    const Vec2d origin = points[loop.begin];
    long double area = 0.0L;
    for (size_t i = loop.begin; i < loop.end; ++i) {
        const size_t next = i + 1 == loop.end ? loop.begin : i + 1;
        area += cross2(points[i] - origin, points[next] - origin);
    }
    return area;
}

bool on_segment(const Vec2d &p, const Vec2d &a, const Vec2d &b, double epsilon)
{
    const Vec2d ab = b - a;
    const double threshold = epsilon * ab.stableNorm();
    return std::abs(cross2(ab, p - a)) <= threshold &&
           (p - a).dot(ab) >= -threshold && (p - b).dot(ab) <= threshold;
}

bool edges_conflict(const Edge &a, const Edge &b, const std::vector<Vec2d> &points, double epsilon)
{
    int common = -1;
    if (a.first == b.first || a.first == b.second) common = a.first;
    else if (a.second == b.first || a.second == b.second) common = a.second;
    if (common >= 0) {
        const Vec2d u = points[a.first == common ? a.second : a.first] - points[common];
        const Vec2d v = points[b.first == common ? b.second : b.first] - points[common];
        // Adjacent straight-through segments are legal. Backtracking overlaps
        // beyond the common endpoint and must not be "repaired" by earcut.
        return std::abs(cross2(u, v)) <= epsilon * std::max(u.stableNorm(), v.stableNorm()) && u.dot(v) > 0.0;
    }
    const Vec2d &p = points[a.first], &q = points[a.second];
    const Vec2d &r = points[b.first], &s = points[b.second];
    if (on_segment(p, r, s, epsilon) || on_segment(q, r, s, epsilon) ||
        on_segment(r, p, q, epsilon) || on_segment(s, p, q, epsilon))
        return true;
    const auto c1 = cross2(q - p, r - p), c2 = cross2(q - p, s - p);
    const auto c3 = cross2(s - r, p - r), c4 = cross2(s - r, q - r);
    return ((c1 > 0 && c2 < 0) || (c1 < 0 && c2 > 0)) &&
           ((c3 > 0 && c4 < 0) || (c3 < 0 && c4 > 0));
}

bool inside_loop(const Vec2d &point, const Loop &loop, const std::vector<Vec2d> &points)
{
    // Boundary contact was rejected separately; a strict winding test suffices.
    int winding = 0;
    for (size_t i = loop.begin; i < loop.end; ++i) {
        const Vec2d &a = points[i], &b = points[i + 1 == loop.end ? loop.begin : i + 1];
        const auto orientation = cross2(b - a, point - a);
        if (a.y() <= point.y() && b.y() > point.y() && orientation > 0) ++winding;
        if (a.y() > point.y() && b.y() <= point.y() && orientation < 0) --winding;
    }
    return winding != 0;
}

FaceExtrusionError validate_contours(const std::vector<Vec2d> &points, const std::vector<Loop> &loops,
                                    const std::vector<Edge> &edges, double epsilon)
{
    for (size_t i = 0; i < points.size(); ++i)
        for (size_t j = i + 1; j < points.size(); ++j)
            if ((points[i] - points[j]).stableNorm() <= epsilon)
                return FaceExtrusionError::InvalidVertex;
    for (size_t i = 0; i < edges.size(); ++i)
        for (size_t j = i + 1; j < edges.size(); ++j)
            if (edges_conflict(edges[i], edges[j], points, epsilon))
                return FaceExtrusionError::SelfIntersection;
    for (const Loop &loop : loops) {
        double perimeter = 0.0;
        for (size_t i = loop.begin; i < loop.end; ++i)
            perimeter += (points[i] - points[i + 1 == loop.end ? loop.begin : i + 1]).stableNorm();
        if (std::abs(area2(points, loop)) <= epsilon * perimeter)
            return FaceExtrusionError::DegenerateTriangle;
    }
    for (size_t i = 1; i < loops.size(); ++i) {
        if (!inside_loop(points[loops[i].begin], loops.front(), points))
            return FaceExtrusionError::SelfIntersection;
        for (size_t j = 1; j < i; ++j)
            if (inside_loop(points[loops[i].begin], loops[j], points) ||
                inside_loop(points[loops[j].begin], loops[i], points))
                return FaceExtrusionError::SelfIntersection; // Nested holes are not islands.
    }
    return FaceExtrusionError::None;
}

// Earcut may filter exact collinear boundary vertices. Retain their original
// positions and topology by subdividing only the corresponding boundary edge;
// never invent a point, move a coordinate or replace the requested outline.
bool restore_boundary_vertices(std::vector<Vec3i> &triangles, const std::vector<Vec2d> &points,
                               const std::vector<int> &previous, const std::vector<int> &next, double epsilon)
{
    std::vector<bool> used(points.size(), false);
    for (const Vec3i &triangle : triangles)
        for (int i : {triangle[0], triangle[1], triangle[2]}) used[size_t(i)] = true;
    for (size_t p = 0; p < points.size(); ++p) {
        if (used[p]) continue;
        int a = previous[p], b = next[p];
        while (!used[size_t(a)] && a != int(p)) a = previous[size_t(a)];
        while (!used[size_t(b)] && b != int(p)) b = next[size_t(b)];
        if (a == int(p) || b == int(p) || a == b || !on_segment(points[p], points[a], points[b], epsilon))
            return false;
        bool inserted = false;
        for (size_t f = 0; f < triangles.size() && !inserted; ++f) {
            const Vec3i old = triangles[f];
            for (int side = 0; side < 3; ++side) {
                const int from = old[side], to = old[(side + 1) % 3], third = old[(side + 2) % 3];
                if ((from == a && to == b) || (from == b && to == a)) {
                    triangles[f] = Vec3i(from, int(p), third);
                    triangles.emplace_back(int(p), to, third);
                    used[p] = inserted = true;
                    break;
                }
            }
        }
        if (!inserted) return false;
    }
    return true;
}

// Roundoff in a rotated projection can make mathematically collinear points
// look like a microscopic ear. Choose a different INTERNAL diagonal instead of
// emitting a zero-area triangle or moving/dropping the original 3D points. Both
// replacements must be nondegenerate and consistently oriented; all original
// boundary directions and area are independently checked afterwards.
bool choose_nonflat_diagonals(std::vector<Vec3i> &triangles, const std::vector<Vec2d> &points, double epsilon)
{
    auto signed_area = [&](const Vec3i &t) { return cross2(points[t[1]] - points[t[0]], points[t[2]] - points[t[0]]); };
    auto usable = [&](const Vec3i &t) {
        const double span = std::max({(points[t[1]] - points[t[0]]).stableNorm(),
                                      (points[t[2]] - points[t[0]]).stableNorm(),
                                      (points[t[2]] - points[t[1]]).stableNorm()});
        return std::abs(signed_area(t)) > epsilon * span;
    };
    for (size_t f = 0; f < triangles.size(); ++f) {
        if (usable(triangles[f])) continue;
        const Vec3i flat = triangles[f];
        std::array<int, 3> sides{{0, 1, 2}};
        std::stable_sort(sides.begin(), sides.end(), [&](int a, int b) {
            return (points[flat[a]] - points[flat[(a + 1) % 3]]).squaredNorm() >
                   (points[flat[b]] - points[flat[(b + 1) % 3]]).squaredNorm();
        });
        bool replaced = false;
        for (int side : sides) {
            const int a = flat[side], b = flat[(side + 1) % 3], middle = flat[(side + 2) % 3];
            for (size_t g = 0; g < triangles.size(); ++g) {
                if (g == f) continue;
                const Vec3i &neighbor = triangles[g];
                bool has_a = false, has_b = false;
                int other = -1;
                for (int i : {neighbor[0], neighbor[1], neighbor[2]}) {
                    has_a |= i == a;
                    has_b |= i == b;
                    if (i != a && i != b) other = i;
                }
                if (!has_a || !has_b || other < 0 || other == middle) continue;
                const Vec3i first(a, middle, other), second(middle, b, other);
                if (usable(first) && usable(second) && ((signed_area(first) > 0) == (signed_area(second) > 0))) {
                    triangles[f] = first;
                    triangles[g] = second;
                    replaced = true;
                    break;
                }
            }
            if (replaced) break;
        }
        if (!replaced) return false;
    }
    return true;
}

} // namespace

PlanarProfileResult triangulate_planar_profile(const std::vector<Vec3d> &outer,
                                               const std::vector<std::vector<Vec3d>> &holes,
                                               const FaceExtrusionOptions &options)
{
    if (!std::isfinite(options.planarity_tolerance_mm) || options.planarity_tolerance_mm < 0.0 ||
        !std::isfinite(options.length_tolerance_mm) || options.length_tolerance_mm <= 0.0 ||
        !std::isfinite(options.angular_tolerance_rad) || options.angular_tolerance_rad <= 0.0 ||
        options.angular_tolerance_rad >= 0.5 * std::acos(-1.0))
        return failure(FaceExtrusionError::InvalidTolerance, "Invalid planar-profile tolerances.");
    if (outer.size() < 3)
        return failure(FaceExtrusionError::EmptySelection, "Pick at least three ordered profile points.");
    size_t count = outer.size();
    if (count > max_profile_vertices)
        return failure(FaceExtrusionError::TooLarge, "The profile exceeds the interactive limit of 4096 vertices.");
    for (const auto &hole : holes) {
        if (hole.size() < 3)
            return failure(FaceExtrusionError::EmptySelection, "Every hole requires at least three ordered points.");
        if (hole.size() > max_profile_vertices - count)
            return failure(FaceExtrusionError::TooLarge, "The profile exceeds the interactive limit of 4096 vertices.");
        count += hole.size();
    }

    std::vector<Vec3d> original = outer;
    original.reserve(count);
    std::vector<Loop> loops{{0, outer.size()}};
    for (const auto &hole : holes) {
        const size_t begin = original.size();
        original.insert(original.end(), hole.begin(), hole.end());
        loops.push_back({begin, original.size()});
    }
    for (const Vec3d &p : original) {
        if (!p.allFinite())
            return failure(FaceExtrusionError::InvalidVertex, "Profile points must be finite.");
        if (p.cwiseAbs().maxCoeff() > std::numeric_limits<float>::max())
            return failure(FaceExtrusionError::PrecisionLoss, "Profile coordinates exceed mesh storage limits.");
    }
    const double epsilon = options.length_tolerance_mm;
    const Vec3d origin = outer.front();
    Vec3d baseline = Vec3d::Zero();
    double longest = 0.0;
    for (const Vec3d &p : outer) {
        const Vec3d delta = p - origin;
        const double length = delta.stableNorm();
        if (length > longest) { longest = length; baseline = delta; }
    }
    Vec3d normal = Vec3d::Zero();
    double largest_cross = 0.0;
    for (const Vec3d &p : outer) {
        const Vec3d candidate = baseline.cross(p - origin);
        const double magnitude = candidate.stableNorm();
        if (magnitude > largest_cross) { largest_cross = magnitude; normal = candidate; }
    }
    if (longest <= epsilon || largest_cross <= epsilon * longest)
        return failure(FaceExtrusionError::DegenerateTriangle, "The outer profile points are collinear or too close.");
    normal /= largest_cross;
    Vec3d u = normal.unitOrthogonal(), v = normal.cross(u);
    std::vector<Vec2d> projected;
    projected.reserve(count);
    double extent = 0.0;
    for (const Vec3d &p : original) {
        const Vec3d delta = p - origin;
        if (std::abs(delta.dot(normal)) > options.planarity_tolerance_mm)
            return failure(FaceExtrusionError::NonPlanar, "The profile and its holes must lie in one plane within the millimeter tolerance.");
        projected.emplace_back(delta.dot(u), delta.dot(v));
        extent = std::max(extent, projected.back().stableNorm());
    }
    const double topology_epsilon = std::max(epsilon, 128.0 * std::numeric_limits<double>::epsilon() * extent);
    std::vector<Edge> edges;
    std::vector<int> previous(count), next(count);
    edges.reserve(count);
    for (const Loop &loop : loops)
        for (size_t i = loop.begin; i < loop.end; ++i) {
            next[i] = int(i + 1 == loop.end ? loop.begin : i + 1);
            previous[i] = int(i == loop.begin ? loop.end - 1 : i - 1);
            edges.emplace_back(int(i), next[i]);
        }
    const auto topology = validate_contours(projected, loops, edges, topology_epsilon);
    if (topology != FaceExtrusionError::None)
        return failure(topology, "Invalid profile: repeated/degenerate points, crossing or touching contours, or a hole outside/inside another hole.");
    if (area2(projected, loops.front()) < 0.0L) {
        normal = -normal;
        v = -v;
        for (Vec2d &p : projected) p.y() = -p.y();
    }
    long double expected_area2 = area2(projected, loops.front());
    for (size_t i = 1; i < loops.size(); ++i) expected_area2 -= std::abs(area2(projected, loops[i]));
    if (!std::isfinite(expected_area2) || expected_area2 <= 0.0L)
        return failure(FaceExtrusionError::DegenerateTriangle, "The profile has no representable positive area.");

    indexed_triangle_set cap;
    cap.vertices.reserve(count);
    std::vector<Vec2d> stored_projection;
    stored_projection.reserve(count);
    const double position_budget = std::max(epsilon, options.planarity_tolerance_mm);
    for (const Vec3d &p : original) {
        const Vec3f encoded = p.cast<float>();
        const Vec3d actual = encoded.cast<double>();
        if ((actual - p).stableNorm() > position_budget ||
            std::abs((actual - origin).dot(normal)) > options.planarity_tolerance_mm)
            return failure(FaceExtrusionError::PrecisionLoss, "The original 3D profile positions cannot be represented within the millimeter tolerance.");
        cap.vertices.push_back(encoded);
        stored_projection.emplace_back((actual - origin).dot(u), (actual - origin).dot(v));
    }
    if (validate_contours(stored_projection, loops, edges, topology_epsilon) != FaceExtrusionError::None)
        return failure(FaceExtrusionError::PrecisionLoss, "Float storage changes the profile topology or closes a gap/hole.");
    for (const Loop &loop : loops)
        if ((area2(projected, loop) > 0.0L) != (area2(stored_projection, loop) > 0.0L))
            return failure(FaceExtrusionError::PrecisionLoss, "Float storage reverses a profile contour.");

    std::vector<std::vector<std::array<double, 2>>> polygon;
    polygon.reserve(loops.size());
    for (const Loop &loop : loops) {
        polygon.emplace_back();
        polygon.back().reserve(loop.end - loop.begin);
        for (size_t i = loop.begin; i < loop.end; ++i)
            polygon.back().push_back({{projected[i].x(), projected[i].y()}});
    }
    const std::vector<uint32_t> indices = mapbox::earcut<uint32_t>(polygon);
    if (indices.empty() || indices.size() % 3 != 0)
        return failure(FaceExtrusionError::NonManifold, "Profile triangulation failed.");
    cap.indices.reserve(count + 2 * holes.size() - 2);
    for (size_t i = 0; i < indices.size(); i += 3) {
        if (indices[i] >= count || indices[i + 1] >= count || indices[i + 2] >= count)
            return failure(FaceExtrusionError::InvalidTriangleIndex, "Profile triangulation returned an invalid index.");
        cap.indices.emplace_back(int(indices[i]), int(indices[i + 1]), int(indices[i + 2]));
    }
    if (!restore_boundary_vertices(cap.indices, projected, previous, next, topology_epsilon) ||
        cap.indices.size() != count + 2 * holes.size() - 2)
        return failure(FaceExtrusionError::NonManifold, "Triangulation did not retain every requested boundary vertex and hole.");

    if (!choose_nonflat_diagonals(cap.indices, projected, topology_epsilon))
        return failure(FaceExtrusionError::DegenerateTriangle, "The requested profile cannot be triangulated without unresolvable thin triangles.");

    struct Use { int count = 0; int direction = 0; int face = -1; int side = -1; };
    std::map<Edge, Use> usage;
    std::vector<Vec3i> neighbors(cap.indices.size(), Vec3i(-1, -1, -1));
    long double triangle_area2 = 0.0L;
    const double cosine = std::cos(options.angular_tolerance_rad);
    for (size_t f = 0; f < cap.indices.size(); ++f) {
        Vec3i &face = cap.indices[f];
        auto area = cross2(projected[face[1]] - projected[face[0]], projected[face[2]] - projected[face[0]]);
        if (area < 0.0L) { std::swap(face[1], face[2]); area = -area; }
        const double span = std::max({(projected[face[1]] - projected[face[0]]).stableNorm(),
                                      (projected[face[2]] - projected[face[0]]).stableNorm(),
                                      (projected[face[2]] - projected[face[1]]).stableNorm()});
        if (area <= topology_epsilon * span)
            return failure(FaceExtrusionError::DegenerateTriangle, "Triangulation contains a degenerate or unresolvable triangle.");
        const Vec3d original_normal = (original[face[1]] - original[face[0]]).cross(original[face[2]] - original[face[0]]);
        if (original_normal.normalized().dot(normal) < cosine)
            return failure(FaceExtrusionError::NonPlanar, "The profile triangulation does not define one plane within the angular tolerance.");
        const Vec3d a = cap.vertices[face[0]].cast<double>(), b = cap.vertices[face[1]].cast<double>(), c = cap.vertices[face[2]].cast<double>();
        const Vec3d stored_normal = (b - a).cross(c - a);
        const Vec3f float_normal = (cap.vertices[face[1]] - cap.vertices[face[0]]).cross(cap.vertices[face[2]] - cap.vertices[face[0]]);
        const double stored_span = std::max({(b - a).stableNorm(), (c - a).stableNorm(), (c - b).stableNorm()});
        if (stored_normal.stableNorm() <= epsilon * stored_span || stored_normal.normalized().dot(normal) < cosine ||
            !float_normal.allFinite() || !std::isfinite(float_normal.norm()) || float_normal.norm() == 0.f)
            return failure(FaceExtrusionError::PrecisionLoss, "Float storage collapses, reverses or makes a cap triangle unrepresentable.");
        triangle_area2 += area;
        for (int side = 0; side < 3; ++side) {
            const int from = face[side], to = face[(side + 1) % 3];
            Use &edge = usage[std::minmax(from, to)];
            ++edge.count;
            edge.direction += from < to ? 1 : -1;
            if (edge.count == 1) { edge.face = int(f); edge.side = side; }
            else if (edge.count == 2) {
                neighbors[f][side] = edge.face;
                neighbors[edge.face][edge.side] = int(f);
            } else
                return failure(FaceExtrusionError::NonManifold, "The triangulated profile contains a non-manifold edge.");
        }
    }
    const long double area_budget = std::max(static_cast<long double>(epsilon) * epsilon,
        128.0L * std::numeric_limits<double>::epsilon() * count * expected_area2);
    if (std::abs(triangle_area2 - expected_area2) > area_budget)
        return failure(FaceExtrusionError::NonManifold, "Triangulation does not cover exactly the requested profile area.");
    std::set<Edge> wanted_boundary, actual_boundary;
    std::map<Edge, int> wanted_direction;
    for (size_t ring = 0; ring < loops.size(); ++ring) {
        const bool reverse = ring == 0 ? area2(projected, loops[ring]) < 0.0L : area2(projected, loops[ring]) > 0.0L;
        for (size_t i = loops[ring].begin; i < loops[ring].end; ++i) {
            const Edge key = std::minmax(int(i), next[i]);
            wanted_boundary.insert(key);
            wanted_direction[key] = (int(i) < next[i] ? 1 : -1) * (reverse ? -1 : 1);
        }
    }
    for (const auto &entry : usage) {
        if (entry.second.count == 1) {
            actual_boundary.insert(entry.first);
            const auto expected = wanted_direction.find(entry.first);
            if (expected == wanted_direction.end() || expected->second != entry.second.direction)
                return failure(FaceExtrusionError::InconsistentWinding, "Triangulation reversed an outline or filled a hole boundary.");
        } else if (entry.second.direction != 0)
            return failure(FaceExtrusionError::InconsistentWinding, "The triangulated cap has inconsistent winding.");
    }
    if (actual_boundary != wanted_boundary)
        return failure(FaceExtrusionError::NonManifold, "Triangulation changed an outline or hole boundary.");
    std::vector<bool> seen(cap.indices.size(), false);
    std::vector<int> pending{0};
    seen[0] = true;
    for (size_t i = 0; i < pending.size(); ++i)
        for (int side = 0; side < 3; ++side) {
            const int neighbor = neighbors[pending[i]][side];
            if (neighbor >= 0 && !seen[neighbor]) { seen[neighbor] = true; pending.push_back(neighbor); }
        }
    const auto diagnostics = its_mesh_diagnostics(cap); // Open cap: no closed-mesh ray test.
    if (pending.size() != cap.indices.size() || diagnostics.non_manifold_vertices != 0 ||
        diagnostics.non_manifold_edges != 0 || diagnostics.same_direction_edges != 0 || diagnostics.open_edges != count)
        return failure(FaceExtrusionError::NonManifold, "Triangulation is not one manifold profile with the requested boundaries.");

    cap.properties.assign(cap.indices.size(), FaceProperty{eNormal, 0.0});
    PlanarProfileResult result;
    result.cap = std::move(cap);
    result.origin = origin;
    result.normal = normal;
    result.area_mm2 = double(expected_area2 * 0.5L);
    result.boundary_edges = count;
    return result;
}

FaceExtrusionResult extrude_planar_profile(const std::vector<Vec3d> &outer, const Vec3d &displacement,
                                          const std::vector<std::vector<Vec3d>> &holes,
                                          const FaceExtrusionOptions &options)
{
    const PlanarProfileResult profile = triangulate_planar_profile(outer, holes, options);
    if (!profile) {
        FaceExtrusionResult result;
        result.error = profile.error;
        result.message = profile.message;
        return result;
    }
    std::vector<int> faces(profile.cap.indices.size());
    std::iota(faces.begin(), faces.end(), 0);
    return extrude_planar_face(profile.cap, faces, displacement, options);
}

} // namespace Slic3r
