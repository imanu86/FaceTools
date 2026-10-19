// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include "FaceExtrusion.hpp"
#include "MeshDiagnostics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <utility>

namespace Slic3r {
namespace {

FaceExtrusionResult failure(FaceExtrusionError error, const char *message)
{
    FaceExtrusionResult result;
    result.error = error;
    result.message = message;
    return result;
}

Vec3d triangle_normal(const indexed_triangle_set &mesh, const Vec3i &face)
{
    const Vec3d a = mesh.vertices[face[0]].cast<double>();
    const Vec3d b = mesh.vertices[face[1]].cast<double>();
    const Vec3d c = mesh.vertices[face[2]].cast<double>();
    return (b - a).cross(c - a);
}

bool degenerate(const indexed_triangle_set &mesh, const Vec3i &face, double epsilon)
{
    const Vec3d a = mesh.vertices[face[0]].cast<double>();
    const Vec3d b = mesh.vertices[face[1]].cast<double>();
    const Vec3d c = mesh.vertices[face[2]].cast<double>();
    const double longest = std::max({(b - a).norm(), (c - a).norm(), (c - b).norm()});
    return longest <= epsilon || (b - a).cross(c - a).norm() <= epsilon * longest;
}

struct HalfEdge {
    int from;
    int to;
    int face;
    int side;

    std::pair<int, int> key() const { return std::minmax(from, to); }
};

double cross2(const Vec2d &a, const Vec2d &b) { return a.x() * b.y() - a.y() * b.x(); }

bool on_segment(const Vec2d &p, const Vec2d &a, const Vec2d &b, double epsilon)
{
    const Vec2d ab = b - a;
    const double length = ab.norm();
    return std::abs(cross2(ab, p - a)) <= epsilon * length &&
           (p - a).dot(ab) >= -epsilon * length &&
           (p - b).dot(ab) <= epsilon * length;
}

// Shared endpoints and a shared full edge are legal; a T-junction, crossing,
// projected duplicate vertex or collinear overlap beyond an endpoint is not.
bool edges_conflict(int a, int b, int c, int d, const std::vector<Vec2d> &points, double epsilon)
{
    if (std::minmax(a, b) == std::minmax(c, d))
        return false;
    int common = -1;
    if (a == c || a == d) common = a;
    else if (b == c || b == d) common = b;
    if (common >= 0) {
        const Vec2d u = points[a == common ? b : a] - points[common];
        const Vec2d v = points[c == common ? d : c] - points[common];
        return std::abs(cross2(u, v)) <= epsilon * std::max(u.norm(), v.norm()) && u.dot(v) > 0.0;
    }

    const Vec2d &pa = points[a], &pb = points[b], &pc = points[c], &pd = points[d];
    if (on_segment(pa, pc, pd, epsilon) || on_segment(pb, pc, pd, epsilon) ||
        on_segment(pc, pa, pb, epsilon) || on_segment(pd, pa, pb, epsilon))
        return true;
    const double c1 = cross2(pb - pa, pc - pa), c2 = cross2(pb - pa, pd - pa);
    const double c3 = cross2(pd - pc, pa - pc), c4 = cross2(pd - pc, pb - pc);
    return ((c1 > 0.0 && c2 < 0.0) || (c1 < 0.0 && c2 > 0.0)) &&
           ((c3 > 0.0 && c4 < 0.0) || (c3 < 0.0 && c4 > 0.0));
}

bool inside_triangle(const Vec2d &p, const Vec3i &triangle, const std::vector<Vec2d> &points, double epsilon)
{
    bool positive = true, negative = true;
    for (int side = 0; side < 3; ++side) {
        const Vec2d &a = points[triangle[side]], &b = points[triangle[(side + 1) % 3]];
        const double value = cross2(b - a, p - a);
        const double threshold = epsilon * (b - a).norm();
        positive &= value > threshold;
        negative &= value < -threshold;
    }
    return positive || negative;
}

bool triangles_conflict(const Vec3i &a, const Vec3i &b, const std::vector<Vec2d> &points, double epsilon)
{
    int shared = 0;
    for (int v : {a[0], a[1], a[2]})
        shared += v == b[0] || v == b[1] || v == b[2];
    // Valid, consistently wound triangles on opposite sides of their shared
    // edge cannot overlap. Winding and duplicate faces were checked earlier.
    if (shared == 2)
        return false;
    for (int i = 0; i < 3; ++i) {
        if (inside_triangle(points[a[i]], b, points, epsilon) ||
            inside_triangle(points[b[i]], a, points, epsilon))
            return true;
        for (int j = 0; j < 3; ++j)
            if (edges_conflict(a[i], a[(i + 1) % 3], b[j], b[(j + 1) % 3], points, epsilon))
                return true;
    }
    return false;
}

// Sweep bounding boxes first: no polygon quantization, union or healing that
// could close a small hole. Worst case is quadratic for overlapping boxes.
FaceExtrusionError validate_projection(const indexed_triangle_set &patch, const Vec3d &normal,
                                      const Vec3d &origin, double epsilon)
{
    const Vec3d u = normal.unitOrthogonal();
    const Vec3d v = normal.cross(u);
    std::vector<Vec2d> points;
    points.reserve(patch.vertices.size());
    for (const Vec3f &p : patch.vertices) {
        const Vec3d relative = p.cast<double>() - origin;
        points.emplace_back(relative.dot(u), relative.dot(v));
    }
    struct Box {
        Vec2d min;
        Vec2d max;
        size_t face;
    };
    std::vector<Box> boxes;
    boxes.reserve(patch.indices.size());
    for (size_t i = 0; i < patch.indices.size(); ++i) {
        const Vec3i &face = patch.indices[i];
        const Vec2d &a = points[face[0]], &b = points[face[1]], &c = points[face[2]];
        boxes.push_back({a.cwiseMin(b).cwiseMin(c), a.cwiseMax(b).cwiseMax(c), i});
    }
    std::sort(boxes.begin(), boxes.end(), [](const Box &a, const Box &b) {
        return a.min.x() < b.min.x() || (a.min.x() == b.min.x() && a.face < b.face);
    });
    // Bound synchronous GUI work even for an adversarial needle tessellation.
    constexpr size_t max_pair_checks = 10000000;
    size_t checks = 0;
    for (size_t i = 0; i < boxes.size(); ++i) {
        for (size_t j = i + 1; j < boxes.size() && boxes[j].min.x() <= boxes[i].max.x() + epsilon; ++j) {
            if (++checks > max_pair_checks)
                return FaceExtrusionError::TooLarge;
            if (boxes[j].min.y() > boxes[i].max.y() + epsilon || boxes[j].max.y() < boxes[i].min.y() - epsilon)
                continue;
            if (triangles_conflict(patch.indices[boxes[i].face], patch.indices[boxes[j].face], points, epsilon))
                return FaceExtrusionError::SelfIntersection;
        }
    }
    return FaceExtrusionError::None;
}

} // namespace

FaceExtrusionResult extrude_planar_face(const indexed_triangle_set &source,
                                       const std::vector<int> &triangle_ids,
                                       const Vec3d &displacement,
                                       const FaceExtrusionOptions &options)
{
    if (!std::isfinite(options.planarity_tolerance_mm) || options.planarity_tolerance_mm < 0.0 ||
        !std::isfinite(options.length_tolerance_mm) || options.length_tolerance_mm <= 0.0 ||
        !std::isfinite(options.angular_tolerance_rad) || options.angular_tolerance_rad <= 0.0 ||
        options.angular_tolerance_rad >= 0.5 * std::acos(-1.0))
        return failure(FaceExtrusionError::InvalidTolerance, "Invalid extrusion tolerances.");
    const double epsilon = options.length_tolerance_mm;
    const double distance = displacement.stableNorm();
    if (!displacement.allFinite() || !std::isfinite(distance) || distance <= epsilon)
        return failure(FaceExtrusionError::InvalidDisplacement, "Extrusion distance must be finite and non-zero.");
    if (triangle_ids.empty())
        return failure(FaceExtrusionError::EmptySelection, "Select a planar face first.");
    // Each input triangle contributes at most 6 vertices and 8 output faces.
    if (triangle_ids.size() > size_t(std::numeric_limits<int>::max()) / 8)
        return failure(FaceExtrusionError::TooLarge, "The selected patch exceeds mesh index limits.");

    std::vector<int> ids = triangle_ids;
    std::sort(ids.begin(), ids.end());
    if (std::adjacent_find(ids.begin(), ids.end()) != ids.end())
        return failure(FaceExtrusionError::DuplicateTriangle, "The selection contains a duplicate triangle.");

    indexed_triangle_set patch;
    patch.indices.reserve(ids.size());
    std::map<std::array<float, 3>, int> welded_vertices;
    std::set<std::array<int, 3>> unique_faces;
    // Welding is EXACT and confined to this copy. A fuzzy weld can erase holes.
    for (int id : ids) {
        if (id < 0 || size_t(id) >= source.indices.size())
            return failure(FaceExtrusionError::InvalidTriangleIndex, "A selected triangle no longer exists.");
        const Vec3i &src = source.indices[id];
        Vec3i face;
        for (int i = 0; i < 3; ++i) {
            if (src[i] < 0 || size_t(src[i]) >= source.vertices.size() || !source.vertices[src[i]].allFinite())
                return failure(FaceExtrusionError::InvalidVertex, "A selected triangle has an invalid vertex.");
            const Vec3f &p = source.vertices[src[i]];
            const std::array<float, 3> key{{p.x(), p.y(), p.z()}};
            const auto inserted = welded_vertices.emplace(key, int(patch.vertices.size()));
            if (inserted.second)
                patch.vertices.push_back(p);
            face[i] = inserted.first->second;
        }
        if (face[0] == face[1] || face[1] == face[2] || face[2] == face[0] || degenerate(patch, face, epsilon))
            return failure(FaceExtrusionError::DegenerateTriangle, "The selected patch contains a degenerate triangle.");
        std::array<int, 3> key{{face[0], face[1], face[2]}};
        std::sort(key.begin(), key.end());
        if (!unique_faces.insert(key).second)
            return failure(FaceExtrusionError::DuplicateTriangle, "The selected patch contains coincident triangles.");
        patch.indices.push_back(face);
    }

    size_t seed = 0;
    double greatest_area = 0.0;
    for (size_t i = 0; i < patch.indices.size(); ++i) {
        const double area = triangle_normal(patch, patch.indices[i]).norm();
        if (area > greatest_area) { greatest_area = area; seed = i; }
    }
    const Vec3d normal = triangle_normal(patch, patch.indices[seed]).normalized();
    const Vec3d origin = patch.vertices[patch.indices[seed][0]].cast<double>();
    const double cos_angle = std::cos(options.angular_tolerance_rad);
    double max_planarity_error = 0.0;
    for (const Vec3f &p : patch.vertices)
        max_planarity_error = std::max(max_planarity_error, std::abs((p.cast<double>() - origin).dot(normal)));
    if (max_planarity_error > options.planarity_tolerance_mm)
        return failure(FaceExtrusionError::NonPlanar, "The selected surface is not planar within the millimeter tolerance.");
    for (const Vec3i &face : patch.indices) {
        const double alignment = triangle_normal(patch, face).normalized().dot(normal);
        if (alignment < 0.0)
            return failure(FaceExtrusionError::InconsistentWinding, "The selected surface contains reversed triangles.");
        if (alignment < cos_angle)
            return failure(FaceExtrusionError::NonPlanar, "The selected triangle normals do not define one plane.");
    }
    const double signed_distance = normal.dot(displacement);
    if (std::abs(signed_distance / distance) < cos_angle)
        return failure(FaceExtrusionError::InvalidDisplacement, "Extrusion must be perpendicular to the selected plane.");

    // Orient the patch towards the translated cap, including negative depths.
    if (signed_distance < 0.0)
        for (Vec3i &face : patch.indices)
            std::swap(face[1], face[2]);
    std::vector<HalfEdge> edges;
    edges.reserve(3 * patch.indices.size());
    for (size_t i = 0; i < patch.indices.size(); ++i)
        for (int side = 0; side < 3; ++side)
            edges.push_back({patch.indices[i][side], patch.indices[i][(side + 1) % 3], int(i), side});
    std::sort(edges.begin(), edges.end(), [](const HalfEdge &a, const HalfEdge &b) {
        return a.key() < b.key() || (a.key() == b.key() && a.face < b.face);
    });
    std::vector<HalfEdge> boundary;
    std::vector<Vec3i> neighbors(patch.indices.size(), Vec3i(-1, -1, -1));
    std::vector<int> boundary_in(patch.vertices.size(), 0), boundary_out(patch.vertices.size(), 0);
    for (size_t i = 0; i < edges.size();) {
        size_t end = i + 1;
        while (end < edges.size() && edges[end].key() == edges[i].key()) ++end;
        if (end - i > 2)
            return failure(FaceExtrusionError::NonManifold, "More than two selected triangles share an edge.");
        if (end - i == 1) {
            boundary.push_back(edges[i]);
            ++boundary_out[edges[i].from];
            ++boundary_in[edges[i].to];
        } else {
            const HalfEdge &a = edges[i], &b = edges[i + 1];
            if (a.from != b.to || a.to != b.from)
                return failure(FaceExtrusionError::InconsistentWinding, "Two triangles share an edge in the same direction.");
            neighbors[a.face][a.side] = b.face;
            neighbors[b.face][b.side] = a.face;
        }
        i = end;
    }
    if (boundary.empty())
        return failure(FaceExtrusionError::NonManifold, "The selected surface has no planar boundary.");
    for (size_t i = 0; i < patch.vertices.size(); ++i)
        if (boundary_in[i] != boundary_out[i] || boundary_in[i] > 1)
            return failure(FaceExtrusionError::NonManifold, "The selected boundary branches or touches itself.");
    // Unlike TriangleMesh's quick stats, this also rejects bow-tie vertices.
    // This is an OPEN patch, so the diagnostic does not perform ray tests.
    const MeshDiagnosticStats diagnostics = its_mesh_diagnostics(patch);
    if (diagnostics.non_manifold_vertices != 0 || diagnostics.non_manifold_edges != 0)
        return failure(FaceExtrusionError::NonManifold, "The selected surface contains a non-manifold vertex or edge.");

    std::vector<bool> seen(patch.indices.size(), false);
    std::vector<int> queue{0};
    seen[0] = true;
    for (size_t i = 0; i < queue.size(); ++i)
        for (int side = 0; side < 3; ++side) {
            const int next = neighbors[queue[i]][side];
            if (next >= 0 && !seen[next]) { seen[next] = true; queue.push_back(next); }
        }
    if (queue.size() != patch.indices.size())
        return failure(FaceExtrusionError::Disconnected, "Select one edge-connected planar region.");
    const FaceExtrusionError projected = validate_projection(patch, normal, origin, epsilon);
    if (projected != FaceExtrusionError::None)
        return failure(projected, projected == FaceExtrusionError::TooLarge ?
                       "The selected patch is too complex to validate interactively." :
                       "The selected triangles overlap, cross or contain a T-junction.");

    indexed_triangle_set output;
    const int count = int(patch.vertices.size());
    output.vertices = patch.vertices;
    output.vertices.reserve(2 * patch.vertices.size());
    const double rounding_tolerance = std::max(epsilon, distance * 1e-4);
    for (const Vec3f &p : patch.vertices) {
        const Vec3d translated = p.cast<double>() + displacement;
        // Guard the narrowing conversion itself, not just the resulting float.
        if (!translated.allFinite() || translated.cwiseAbs().maxCoeff() > std::numeric_limits<float>::max())
            return failure(FaceExtrusionError::PrecisionLoss, "Extrusion coordinates exceed mesh storage limits.");
        const Vec3f top = translated.cast<float>();
        const Vec3d actual = top.cast<double>() - p.cast<double>();
        if ((actual - displacement).norm() > rounding_tolerance || actual.dot(displacement / distance) <= epsilon)
            return failure(FaceExtrusionError::PrecisionLoss, "The requested thickness cannot be represented at these coordinates.");
        output.vertices.push_back(top);
    }
    // Float storage can collapse a narrow opening at a large translation even
    // when the thickness itself is representable. Revalidate the translated cap.
    indexed_triangle_set top_patch;
    top_patch.vertices.assign(output.vertices.begin() + count, output.vertices.end());
    top_patch.indices = patch.indices;
    for (const Vec3i &face : top_patch.indices)
        if (degenerate(top_patch, face, epsilon) ||
            triangle_normal(top_patch, face).dot(triangle_normal(patch, face)) <= 0.0)
            return failure(FaceExtrusionError::PrecisionLoss, "Float conversion collapses or reverses the translated cap.");
    if (validate_projection(top_patch, normal, origin + displacement, epsilon) != FaceExtrusionError::None)
        return failure(FaceExtrusionError::PrecisionLoss, "Float conversion changes the translated cap topology.");
    output.indices.reserve(2 * patch.indices.size() + 2 * boundary.size());
    double expected_volume = 0.0;
    for (const Vec3i &face : patch.indices) {
        output.indices.emplace_back(face[0], face[2], face[1]);
        output.indices.emplace_back(face[0] + count, face[1] + count, face[2] + count);
        expected_volume += 0.5 * triangle_normal(patch, face).dot(displacement);
    }
    for (const HalfEdge &edge : boundary) {
        output.indices.emplace_back(edge.from, edge.to, edge.to + count);
        output.indices.emplace_back(edge.from, edge.to + count, edge.from + count);
    }
    double volume = 0.0;
    for (const Vec3i &face : output.indices) {
        if (degenerate(output, face, epsilon))
            return failure(FaceExtrusionError::PrecisionLoss, "Extrusion produced a degenerate triangle after float conversion.");
        const Vec3d a = output.vertices[face[0]].cast<double>() - origin;
        const Vec3d b = output.vertices[face[1]].cast<double>() - origin;
        const Vec3d c = output.vertices[face[2]].cast<double>() - origin;
        volume += a.dot(b.cross(c)) / 6.0;
        // TriangleMesh's volume/normal APIs also perform these operations in
        // float. Reject extreme coordinates that overflow those downstream.
        const Vec3f float_normal = (output.vertices[face[1]] - output.vertices[face[0]]).cross(
                                    output.vertices[face[2]] - output.vertices[face[0]]);
        if (!float_normal.allFinite() || !std::isfinite(float_normal.norm()) || float_normal.norm() == 0.f)
            return failure(FaceExtrusionError::PrecisionLoss, "Extrusion exceeds mesh arithmetic limits.");
    }
    if (!std::isfinite(volume) || volume <= 0.0 || volume > std::numeric_limits<float>::max() ||
        std::abs(volume - expected_volume) > std::max(expected_volume * 1e-3, epsilon * epsilon * epsilon))
        return failure(FaceExtrusionError::PrecisionLoss, "Extrusion volume failed the numerical consistency check.");
    const MeshDiagnosticStats closed = its_edge_diagnostics(output);
    if (closed.open_edges != 0 || closed.non_manifold_edges != 0 || closed.same_direction_edges != 0)
        return failure(FaceExtrusionError::NonManifold, "Extrusion failed the closed-surface check.");
    output.properties.assign(output.indices.size(), FaceProperty{eNormal, 0.0});
    FaceExtrusionResult result;
    result.mesh = std::move(output);
    result.boundary_edges = boundary.size();
    result.max_planarity_error_mm = max_planarity_error;
    return result;
}

} // namespace Slic3r
