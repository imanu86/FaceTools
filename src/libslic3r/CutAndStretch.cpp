// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include "CutAndStretch.hpp"

#include "FaceExtrusion.hpp"
#include "MeshBoolean.hpp"
#include "MeshDiagnostics.hpp"
#include "TriangleMeshSlicer.hpp"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <Eigen/LU>
#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Slic3r {

struct CutAndStretchSourceFrame {
    Matrix3d source_to_cut = Matrix3d::Identity();
    Matrix3d cut_to_source = Matrix3d::Identity();
    Vec3d source_anchor = Vec3d::Zero();
    Vec3d cut_anchor = Vec3d::Zero();
    Transform3d cut_to_world = Transform3d::Identity();
    double determinant = 1.0;
    std::array<Vec3d, 3> source_plane;
    Vec3d source_plane_metric = Vec3d::Zero();
    bool has_source_plane = false;
    Transform3d original_placement = Transform3d::Identity();
    Vec3d requested_origin = Vec3d::Zero();
    Vec3d requested_normal = Vec3d::UnitZ();

    CGAL::Exact_predicates_exact_constructions_kernel::Plane_3 exact_plane() const
    {
        using K = CGAL::Exact_predicates_exact_constructions_kernel;
        const auto point = [](const Vec3d &v) { return K::Point_3(v.x(), v.y(), v.z()); };
        if (has_source_plane) return K::Plane_3(point(source_plane[0]), point(source_plane[1]), point(source_plane[2]));
        // The generic path certifies ONLY the supplied mathematical world plane,
        // not a nearby source face or an interval interpreted as user intent.
        std::array<K::FT, 4> c = {{0, 0, 0, 0}};
        for (int row = 0; row < 3; ++row) {
            const K::FT n(requested_normal[row]);
            for (int col = 0; col < 3; ++col) c[col] += n * K::FT(original_placement.linear()(row, col));
            c[3] += n * (K::FT(original_placement.translation()[row]) - K::FT(requested_origin[row]));
        }
        return K::Plane_3(c[0], c[1], c[2], c[3]);
    }

    bool on_source_plane(const Vec3d &point) const
    {
        if (!has_source_plane) return false;
        using K = CGAL::Exact_predicates_exact_constructions_kernel;
        const auto p = [](const Vec3d &v) { return K::Point_3(v.x(), v.y(), v.z()); };
        return CGAL::coplanar(p(source_plane[0]), p(source_plane[1]), p(source_plane[2]), p(point));
    }
    Vec3d to_cut(const Vec3d &point) const
    {
        Vec3d out = source_to_cut * (point - source_anchor) + cut_anchor;
        if (has_source_plane) {
            out.z() = source_plane_metric.dot(point - source_plane[0]);
            // Explicit original-facet provenance only: this computes a true zero
            // coordinate, never projects an off-plane source vertex onto a face.
            if (on_source_plane(point)) out.z() = 0.0;
        } else {
            out.z() = cut_to_world.linear().col(2).dot(original_placement.linear() * point +
                       (original_placement.translation() - requested_origin));
        }
        return out;
    }
    Vec3d to_source(const Vec3d &point) const { return source_anchor + cut_to_source * (point - cut_anchor); }
};

struct CutAndStretchGuidedState {
    // Read-only, source-frame operands retained only to reject collisions with
    // regions which the user did NOT authorize to move. Never Model pointers.
    std::vector<indexed_triangle_set> excluded_upper;
    CutAndStretchSelectionSummary selection;
};

namespace {

// This is the existing cut_mesh classifier/welder's band, not a geometry offset.
constexpr double cut_band = double(1e-3f);
constexpr double max_rounding_mm = 1e-4;
constexpr size_t max_faces = 2000000;
constexpr size_t max_cap_faces = 20000;
// Bound new component/Boolean work without changing legacy resource limits.
constexpr size_t max_guided_components = 64;
struct Canceled {};

struct Work {
    const std::function<bool()> &cancel;
    const std::function<void(int)> &progress;
    void check() const { if (cancel && cancel()) throw Canceled{}; }
    void report(int percent) const { check(); if (progress) progress(percent); check(); }
};

void require(bool condition, const char *message)
{
    if (!condition) throw std::runtime_error(message);
}

struct SolidInfo {
    double volume = 0.0;
    double area = 0.0;
    Vec3d min = Vec3d::Constant(std::numeric_limits<double>::infinity());
    Vec3d max = Vec3d::Constant(-std::numeric_limits<double>::infinity());
};

SolidInfo inspect(const indexed_triangle_set &mesh, const Work &work, bool closed = true,
                  bool single_component = true)
{
    work.check();
    require(!mesh.indices.empty() && !mesh.vertices.empty(), "The cut produced an empty surface or solid.");
    require(mesh.indices.size() <= max_faces && mesh.vertices.size() <= size_t(std::numeric_limits<int>::max()),
            "This mesh is too large for a safely validated cut and extension.");
    for (size_t i = 0; i < mesh.vertices.size(); ++i) {
        if ((i & 1023) == 0) work.check();
        require(mesh.vertices[i].allFinite(), "The mesh contains a non-finite vertex.");
    }
    SolidInfo info;
    const Vec3d origin = mesh.vertices.front().cast<double>();
    double compensation = 0.0;
    for (size_t i = 0; i < mesh.indices.size(); ++i) {
        if ((i & 1023) == 0) work.check();
        const Vec3i &f = mesh.indices[i];
        for (int index : {f[0], f[1], f[2]}) {
            require(index >= 0 && size_t(index) < mesh.vertices.size(), "The mesh contains an invalid triangle index.");
            const Vec3d p = mesh.vertices[index].cast<double>();
            info.min = info.min.cwiseMin(p);
            info.max = info.max.cwiseMax(p);
        }
        const Vec3d a = mesh.vertices[f[0]].cast<double>() - origin;
        const Vec3d b = mesh.vertices[f[1]].cast<double>() - origin;
        const Vec3d c = mesh.vertices[f[2]].cast<double>() - origin;
        const Vec3d cross = (b - a).cross(c - a);
        const Vec3f stored_cross = (mesh.vertices[f[1]] - mesh.vertices[f[0]]).cross(mesh.vertices[f[2]] - mesh.vertices[f[0]]);
        require(cross.allFinite() && cross.norm() > 0.0 && stored_cross.allFinite() &&
                std::isfinite(stored_cross.squaredNorm()) && stored_cross.squaredNorm() > 0.f,
                "The mesh contains a degenerate or unrepresentable triangle.");
        info.area += 0.5 * cross.norm();
        const double term = a.dot(b.cross(c)) / 6.0 - compensation;
        const double sum = info.volume + term;
        compensation = (sum - info.volume) - term;
        info.volume = sum;
    }
    require(std::isfinite(info.area) && info.area > 0.0, "The mesh has invalid surface area.");
    if (!closed) return info;

    std::vector<Vec3i> neighbors;
    const auto edge_stats = its_edge_diagnostics(mesh, &neighbors);
    require(edge_stats.open_edges == 0 && edge_stats.non_manifold_edges == 0 && edge_stats.same_direction_edges == 0,
            "Cut and extend requires closed, consistently oriented manifold meshes.");
    // Source, legacy operands and each guided component require one boundary.
    // The guided upper COLLECTION alone may be disconnected, before splitting;
    // all its other manifold/orientation/intersection/solid checks still run.
    std::vector<unsigned char> seen(mesh.indices.size(), 0);
    std::vector<int> queue{0};
    seen[0] = 1;
    for (size_t i = 0; i < queue.size(); ++i) {
        if ((i & 1023) == 0) work.check();
        for (int next : {neighbors[queue[i]][0], neighbors[queue[i]][1], neighbors[queue[i]][2]}) {
            require(next >= 0 && size_t(next) < seen.size(), "The cut surface has an incomplete seam.");
            if (!seen[next]) { seen[next] = 1; queue.push_back(next); }
        }
    }
    require(!single_component || queue.size() == mesh.indices.size(),
            "Disconnected solids or separate cavity shells are not supported by this cut and extension.");
    require(std::isfinite(info.volume) && info.volume > 0.0 && std::isfinite(its_volume(mesh)) && its_volume(mesh) > 0.f,
            "The mesh does not enclose a representable positive volume.");
    const auto diagnostics = its_mesh_diagnostics(mesh);
    work.check();
    require(!diagnostics.non_manifold_vertices && !diagnostics.has_reversed_faces,
            "The mesh has a non-manifold vertex or reversed faces.");
    auto solid = MeshBoolean::cgal::triangle_mesh_to_cgal(mesh);
    work.check();
    require(!MeshBoolean::cgal::does_self_intersect(*solid), "The mesh contains self-intersections.");
    work.check();
    require(MeshBoolean::cgal::does_bound_a_volume(*solid), "The mesh does not bound a consistently oriented solid.");
    work.check();
    return info;
}

double rounding_budget(const SolidInfo &info)
{
    const double magnitude = std::max({1.0, info.min.cwiseAbs().maxCoeff(), info.max.cwiseAbs().maxCoeff()});
    const double budget = std::max(1e-6, 2.0 * std::numeric_limits<float>::epsilon() * magnitude);
    require(budget <= max_rounding_mm, "These coordinates exceed the precision budget for a seam-preserving cut. Move the model nearer the origin or use a larger-scale feature.");
    return budget;
}

SolidInfo transformed_bounds(const SolidInfo &info, const Transform3d &transform)
{
    SolidInfo out;
    for (int corner = 0; corner < 8; ++corner) {
        Vec3d p;
        for (int axis = 0; axis < 3; ++axis) p[axis] = (corner & (1 << axis)) ? info.max[axis] : info.min[axis];
        p = transform * p;
        require(p.allFinite(), "The transformed bounds exceed coordinate storage limits.");
        out.min = out.min.cwiseMin(p);
        out.max = out.max.cwiseMax(p);
    }
    return out;
}

void volume_matches(double actual, double expected, double area, double budget, const char *message)
{
    const double allowance = area * budget + 64.0 * std::numeric_limits<double>::epsilon() * std::abs(expected);
    require(std::isfinite(actual) && std::isfinite(expected) && std::abs(actual - expected) <= allowance, message);
}

indexed_triangle_set transformed(const indexed_triangle_set &source, const Transform3d &transform,
                                 double budget, const Work &work)
{
    indexed_triangle_set out = source;
    require(transform.matrix().allFinite(), "The plane transformation is not finite.");
    for (size_t i = 0; i < out.vertices.size(); ++i) {
        if ((i & 1023) == 0) work.check();
        const Vec3d exact = transform * source.vertices[i].cast<double>();
        require(exact.allFinite() && exact.cwiseAbs().maxCoeff() <= std::numeric_limits<float>::max(),
                "The transformed mesh exceeds coordinate storage limits.");
        out.vertices[i] = exact.cast<float>();
        require((out.vertices[i].cast<double>() - exact).norm() <= budget,
                "The transformation loses too much precision at the cut seam.");
    }
    // prepare/apply use right-handed rigid frames only. No shell flipping or
    // mesh repair is performed. The caller supplies a correctly wound world mesh.
    out.properties.assign(out.indices.size(), FaceProperty{eNormal, 0.0});
    return out;
}

double physical_area(const indexed_triangle_set &mesh, const Matrix3d &linear, const Work &work)
{
    double area = 0.0;
    for (size_t i = 0; i < mesh.indices.size(); ++i) {
        if ((i & 1023) == 0) work.check();
        const Vec3i &f = mesh.indices[i];
        const Vec3d a = mesh.vertices[f[0]].cast<double>();
        const Vec3d ab = linear * (mesh.vertices[f[1]].cast<double>() - a);
        const Vec3d ac = linear * (mesh.vertices[f[2]].cast<double>() - a);
        area += 0.5 * ab.cross(ac).norm();
    }
    require(std::isfinite(area) && area > 0.0, "The physical surface area is not representable.");
    return area;
}

void set_source_facet(const indexed_triangle_set &source, size_t facet, const Matrix3d &inverse,
                      Vec3d &normal, CutAndStretchSourceFrame &frame)
{
    require(facet < source.indices.size(), "The selected plane facet no longer exists in the source mesh.");
    const Vec3i &f = source.indices[facet]; // inspect() already validated all source indices
    for (int i = 0; i < 3; ++i) frame.source_plane[i] = source.vertices[f[i]].cast<double>();
    frame.has_source_plane = true;
    require(!frame.exact_plane().is_degenerate(), "The selected plane facet is degenerate.");
    Vec3d local_normal = (frame.source_plane[1] - frame.source_plane[0]).cross(frame.source_plane[2] - frame.source_plane[0]).normalized();
    const Vec3d covector = inverse.transpose() * local_normal;
    require(covector.allFinite() && std::isfinite(covector.norm()) && covector.norm() > 0.0,
            "The selected facet normal is not representable in physical millimeters.");
    Vec3d candidate = covector.normalized();
    if (candidate.dot(normal) < 0.0) { candidate = -candidate; local_normal = -local_normal; }
    require(candidate.cross(normal).norm() <= 128.0 * std::numeric_limits<double>::epsilon(),
            "The moving direction must be parallel to the selected source-facet normal.");
    normal = candidate;
    frame.source_plane_metric = local_normal / covector.norm();
}

SolidInfo source_frame_bounds(const SolidInfo &local, const CutAndStretchSourceFrame &frame)
{
    SolidInfo out;
    for (int corner = 0; corner < 8; ++corner) {
        Vec3d p;
        for (int axis = 0; axis < 3; ++axis) p[axis] = (corner & (1 << axis)) ? local.max[axis] : local.min[axis];
        p = frame.to_cut(p);
        require(p.allFinite(), "The physical cut-frame bounds exceed numeric limits.");
        out.min = out.min.cwiseMin(p);
        out.max = out.max.cwiseMax(p);
    }
    return out;
}

indexed_triangle_set source_frame_transform(const indexed_triangle_set &mesh, const CutAndStretchSourceFrame &frame,
                                            bool to_source, double budget, const Work &work)
{
    indexed_triangle_set out = mesh;
    using K = CGAL::Exact_predicates_exact_constructions_kernel;
    const auto plane = frame.exact_plane();
    for (size_t i = 0; i < mesh.vertices.size(); ++i) {
        if ((i & 1023) == 0) work.check();
        const Vec3d original = mesh.vertices[i].cast<double>();
        Vec3d exact = to_source ? frame.to_source(original) : frame.to_cut(original);
        if (!to_source && std::abs(exact.z()) < cut_band) {
            require(plane.has_on(K::Point_3(original.x(), original.y(), original.z())),
                    "A noncoplanar source feature lies inside the native cut uncertainty band; the plane must not be snapped to a face.");
            exact.z() = 0.0; // certified zero in the supplied plane, never a tolerance weld
        }
        require(exact.allFinite() && exact.cwiseAbs().maxCoeff() <= std::numeric_limits<float>::max(),
                "The source-local conversion exceeds coordinate storage limits.");
        out.vertices[i] = exact.cast<float>();
        const Vec3d narrowed = out.vertices[i].cast<double>();
        // Local-coordinate error alone is meaningless under instance/volume
        // scale. Validate the actual float backcast in PHYSICAL millimeters.
        const Vec3d error = to_source ? Vec3d(frame.source_to_cut * (narrowed - exact)) : Vec3d(narrowed - exact);
        require(error.allFinite() && error.norm() <= budget,
                "Source-local storage loses too much physical precision at the cut seam.");
        const Vec3d restored = to_source ? frame.to_cut(narrowed) : frame.to_cut(frame.to_source(narrowed));
        require(restored.allFinite() && (restored - (to_source ? original : narrowed)).norm() <= budget,
                "The source-local physical roundtrip cannot preserve the cut geometry.");
    }
    if (frame.determinant < 0.0)
        for (Vec3i &face : out.indices) std::swap(face[1], face[2]);
    out.properties.assign(out.indices.size(), FaceProperty{eNormal, 0.0});
    return out;
}

// Fail before cut_mesh can silently collapse distinct section vertices in its
// fixed 1e-3 XY weld / Z classification band. Exact coincident candidates are
// the same geometric point; nearby distinct candidates are NOT merged here.
void guard_cut_band(const indexed_triangle_set &mesh, double budget, const Work &work)
{
    // Include both endpoints' narrowing uncertainty in the reject band; this
    // does not move geometry. In particular, native epson is a FLOAT constant.
    const double separation = cut_band + 2.0 * budget;
    std::vector<Vec2d> points;
    std::set<std::pair<int, int>> edges;
    for (size_t i = 0; i < mesh.vertices.size(); ++i) {
        if ((i & 1023) == 0) work.check();
        const Vec3d p = mesh.vertices[i].cast<double>();
        const double z = std::abs(p.z());
        require(z <= budget || z > cut_band + budget,
                "A feature lies in the native cutter's 0.001 mm uncertainty band. Choose another plane.");
        if (z <= budget) points.emplace_back(p.x(), p.y());
    }
    for (size_t i = 0; i < mesh.indices.size(); ++i) {
        if ((i & 1023) == 0) work.check();
        const Vec3i &f = mesh.indices[i];
        for (int side = 0; side < 3; ++side) {
            int a = f[side], b = f[(side + 1) % 3];
            const auto key = std::minmax(a, b);
            if (!edges.insert(key).second) continue;
            const Vec3d p = mesh.vertices[a].cast<double>(), q = mesh.vertices[b].cast<double>();
            if ((p.z() < -cut_band && q.z() > cut_band) || (q.z() < -cut_band && p.z() > cut_band)) {
                const Vec3d hit = p + (q - p) * (-p.z() / (q.z() - p.z()));
                points.emplace_back(hit.x(), hit.y());
            }
        }
    }
    std::sort(points.begin(), points.end(), [](const Vec2d &a, const Vec2d &b) {
        return a.x() < b.x() || (a.x() == b.x() && a.y() < b.y());
    });
    points.erase(std::unique(points.begin(), points.end(), [](const Vec2d &a, const Vec2d &b) { return a == b; }), points.end());
    size_t checks = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        work.check();
        for (size_t j = i + 1; j < points.size() && points[j].x() - points[i].x() <= separation; ++j) {
            require(++checks <= 10000000, "The cut section is too complex to validate safely.");
            require(std::abs(points[j].y() - points[i].y()) > separation,
                    "Distinct section vertices are closer than the native cutter's 0.001 mm welding band.");
        }
    }
}

using ExactKernel = CGAL::Exact_predicates_exact_constructions_kernel;
using ExactPolygon = CGAL::Polygon_2<ExactKernel>;
using ExactSet = CGAL::Polygon_set_2<ExactKernel>;
struct CapVertexInfo { int index = -1; };
struct CapFaceInfo { int nesting = -1; };
using CapVertexBase = CGAL::Triangulation_vertex_base_with_info_2<CapVertexInfo, ExactKernel>;
using CapFaceBase = CGAL::Triangulation_face_base_with_info_2<CapFaceInfo, ExactKernel,
    CGAL::Constrained_triangulation_face_base_2<ExactKernel>>;
using CapTDS = CGAL::Triangulation_data_structure_2<CapVertexBase, CapFaceBase>;
using CapCDT = CGAL::Constrained_Delaunay_triangulation_2<ExactKernel, CapTDS, CGAL::Exact_predicates_tag>;

// The native cap triangulator emits zero-area ears at collinear boundary
// vertices. Dropping such ears leaves T-junctions. Instead, cap only the native
// cut's open boundary with a constrained triangulation that retains EVERY seam
// vertex, including collinear ones. No vertex is moved, merged or synthesized.
void cap_cut_boundary(indexed_triangle_set &mesh, std::vector<int> &source_faces,
                      bool normals_down, double budget, const Work &work)
{
    inspect(mesh, work, false);
    require(source_faces.size() == mesh.indices.size(), "The cutter did not return complete face provenance.");
    std::vector<Vec3i> neighbors;
    const auto edges = its_edge_diagnostics(mesh, &neighbors);
    require(!edges.non_manifold_edges && !edges.same_direction_edges,
            "The native cut produced inconsistent shell edges.");
    std::vector<std::pair<int, int>> boundary;
    std::map<int, int> incoming, outgoing;
    std::set<int> ids;
    for (size_t i = 0; i < mesh.indices.size(); ++i) {
        if ((i & 1023) == 0) work.check();
        for (int side = 0; side < 3; ++side) {
            if (neighbors[i][side] >= 0) continue;
            int a = mesh.indices[i][side], b = mesh.indices[i][(side + 1) % 3];
            require(std::abs(mesh.vertices[a].z()) <= budget && std::abs(mesh.vertices[b].z()) <= budget,
                    "The native cut left an opening away from the cut plane.");
            require(++outgoing[a] == 1 && ++incoming[b] == 1,
                    "The cut section branches or touches itself.");
            boundary.emplace_back(a, b);
            ids.insert(a); ids.insert(b);
        }
    }
    require(!boundary.empty() && boundary.size() <= max_cap_faces,
            "The cut has no bounded section, or its section is too complex.");
    for (int id : ids)
        require(incoming[id] == 1 && outgoing[id] == 1, "The native cut boundary does not form closed loops.");

    CapCDT cdt;
    std::map<int, CapCDT::Vertex_handle> handles;
    for (int id : ids) {
        work.check();
        const Vec3f &p = mesh.vertices[id];
        auto handle = cdt.insert(ExactKernel::Point_2(double(p.x()), double(p.y())));
        require(handle->info().index < 0,
                "Distinct cut vertices have the same plane projection; this section cannot be capped without welding.");
        handle->info().index = id;
        handles.emplace(id, handle);
    }
    for (const auto &edge : boundary) {
        work.check();
        cdt.insert_constraint(handles.at(edge.first), handles.at(edge.second));
    }
    require(cdt.number_of_vertices() == ids.size(),
            "The cut boundary crosses itself and would require new intersection vertices.");
    for (auto v = cdt.finite_vertices_begin(); v != cdt.finite_vertices_end(); ++v)
        require(v->info().index >= 0, "The cap triangulator unexpectedly created a seam vertex.");

    std::vector<CapCDT::Edge> frontier;
    auto flood = [&](CapCDT::Face_handle start, int nesting) {
        std::vector<CapCDT::Face_handle> queue{start};
        for (size_t i = 0; i < queue.size(); ++i) {
            if ((i & 1023) == 0) work.check();
            auto face = queue[i];
            if (face->info().nesting >= 0) continue;
            face->info().nesting = nesting;
            for (int side = 0; side < 3; ++side) {
                const CapCDT::Edge edge(face, side);
                auto next = face->neighbor(side);
                if (next->info().nesting >= 0) continue;
                if (cdt.is_constrained(edge)) frontier.push_back(edge);
                else queue.push_back(next);
            }
        }
    };
    flood(cdt.infinite_face(), 0);
    for (size_t i = 0; i < frontier.size(); ++i) {
        work.check();
        const auto edge = frontier[i];
        auto neighbor = edge.first->neighbor(edge.second);
        if (neighbor->info().nesting < 0) flood(neighbor, edge.first->info().nesting + 1);
    }
    const size_t first_cap = mesh.indices.size();
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
        if (face->info().nesting % 2 != 1) continue; // preserve holes, no hull fill
        Vec3i triangle(face->vertex(0)->info().index, face->vertex(1)->info().index, face->vertex(2)->info().index);
        if (normals_down) std::swap(triangle[1], triangle[2]);
        mesh.indices.push_back(triangle);
        source_faces.push_back(-1);
        require(mesh.indices.size() - first_cap <= max_cap_faces, "The cut cap is too complex to validate safely.");
    }
    require(mesh.indices.size() > first_cap, "The cut contact section has zero area.");
    mesh.properties.assign(mesh.indices.size(), FaceProperty{eNormal, 0.0});
    work.check();
}

indexed_triangle_set cap_from_ids(const indexed_triangle_set &mesh, const std::vector<int> &source_faces,
                                  double budget, const Work &work)
{
    require(source_faces.size() == mesh.indices.size(), "The cutter did not return complete cap provenance.");
    indexed_triangle_set cap;
    std::map<int, int> remap;
    for (size_t i = 0; i < mesh.indices.size(); ++i) {
        if ((i & 1023) == 0) work.check();
        if (source_faces[i] != -1) continue;
        Vec3i f;
        for (int side = 0; side < 3; ++side) {
            const int old = mesh.indices[i][side];
            auto inserted = remap.emplace(old, int(cap.vertices.size()));
            if (inserted.second) {
                require(std::abs(mesh.vertices[old].z()) <= budget, "A cut cap is not on the requested plane within the precision budget.");
                cap.vertices.push_back(mesh.vertices[old]); // never flatten the cap
            }
            f[side] = inserted.first->second;
        }
        cap.indices.push_back(f);
        require(cap.indices.size() <= max_cap_faces, "The cut cap is too complex for exact contact validation.");
    }
    inspect(cap, work, false);
    return cap;
}

ExactSet cap_footprint(const indexed_triangle_set &cap, const Work &work)
{
    ExactSet out;
    // Double constructors represent the stored float coordinates exactly. No
    // Clipper quantization, expanded boundary or area-tolerance subset test.
    for (size_t i = 0; i < cap.indices.size(); ++i) {
        if ((i & 63) == 0) work.check();
        ExactPolygon triangle;
        for (int index : {cap.indices[i][0], cap.indices[i][1], cap.indices[i][2]}) {
            const Vec3f &p = cap.vertices[index];
            triangle.push_back(ExactKernel::Point_2(double(p.x()), double(p.y())));
        }
        require(triangle.orientation() != CGAL::COLLINEAR, "A cut cap contains a degenerate projected triangle.");
        if (triangle.orientation() == CGAL::CLOCKWISE) triangle.reverse_orientation();
        out.join(triangle);
    }
    work.check();
    return out;
}

indexed_triangle_set base_cap(const indexed_triangle_set &unit_bridge, const Work &work)
{
    // Re-extrude the untouched base instead of scaling its tiny representational
    // Z residual along with the height: scaling would move the stationary seam.
    require(unit_bridge.vertices.size() % 2 == 0, "The prepared bridge is invalid.");
    const int n = int(unit_bridge.vertices.size() / 2);
    indexed_triangle_set cap;
    cap.vertices.assign(unit_bridge.vertices.begin(), unit_bridge.vertices.begin() + n);
    for (const Vec3i &f : unit_bridge.indices)
        if (f[0] < n && f[1] < n && f[2] < n) cap.indices.emplace_back(f[0], f[2], f[1]);
    inspect(cap, work, false);
    return cap;
}

FaceExtrusionResult extrude(const indexed_triangle_set &cap, double height)
{
    std::vector<int> ids(cap.indices.size());
    std::iota(ids.begin(), ids.end(), 0);
    return extrude_planar_face(cap, ids, Vec3d(0.0, 0.0, height));
}

indexed_triangle_set union_solids(const indexed_triangle_set &a, const indexed_triangle_set &b,
                                  const Work &work, int begin, int end)
{
    work.report(begin);
    indexed_triangle_set output;
    if (MeshBoolean::openmeshcraft::is_enabled()) {
        std::vector<TriangleMesh> results;
        bool failed = false;
        MeshBoolean::openmeshcraft::make_boolean(TriangleMesh(a), TriangleMesh(b), results, "UNION",
            [&work]() { return work.cancel && work.cancel(); },
            [&work, begin, end](float p) { work.report(begin + int((end - begin) * std::clamp(p, 0.f, 100.f) / 100.f)); },
            [&failed]() { failed = true; });
        work.check();
        require(!failed && results.size() == 1 && !results.front().empty(),
                "The Boolean backend could not join the exact cut interfaces safely.");
        output = std::move(results.front().its);
    } else {
        // Unlike legacy MCUT's wrapper, this path never flips independent
        // shells to hide a failed orientation. Validate the actual CSG result.
        output = a;
        MeshBoolean::cgal::plus(output, b);
        work.check();
    }
    inspect(output, work);
    output.properties.assign(output.indices.size(), FaceProperty{eNormal, 0.0});
    work.report(end);
    return output;
}

struct GuidedComponent {
    indexed_triangle_set mesh;
    std::vector<int> source_faces;
};

std::vector<GuidedComponent> split_guided_upper(const indexed_triangle_set &mesh,
    const std::vector<int> &source_faces, const Work &work)
{
    require(source_faces.size() == mesh.indices.size(), "The upper cut lost source-facet provenance.");
    std::vector<Vec3i> neighbors;
    const auto edges = its_edge_diagnostics(mesh, &neighbors);
    require(!edges.open_edges && !edges.non_manifold_edges && !edges.same_direction_edges,
            "The upper cut components do not have closed, consistent seams.");
    std::vector<unsigned char> seen(mesh.indices.size(), 0);
    std::vector<int> remap(mesh.vertices.size(), -1);
    std::vector<GuidedComponent> components;
    for (size_t root = 0; root < mesh.indices.size(); ++root) {
        if ((root & 1023) == 0) work.check();
        if (seen[root]) continue;
        require(components.size() < max_guided_components,
                "Guided cut supports at most 64 upper components; choose another plane.");
        GuidedComponent component;
        std::vector<int> queue{int(root)}, touched;
        seen[root] = 1;
        for (size_t head = 0; head < queue.size(); ++head) {
            if ((head & 1023) == 0) work.check();
            const int face_id = queue[head];
            Vec3i face;
            for (int side = 0; side < 3; ++side) {
                const int old = mesh.indices[face_id][side];
                if (remap[old] < 0) {
                    remap[old] = int(component.mesh.vertices.size());
                    component.mesh.vertices.push_back(mesh.vertices[old]);
                    touched.push_back(old);
                }
                face[side] = remap[old];
                const int next = neighbors[face_id][side];
                require(next >= 0 && size_t(next) < seen.size(), "An upper component has an incomplete seam.");
                if (!seen[next]) { seen[next] = 1; queue.push_back(next); }
            }
            component.mesh.indices.push_back(face);
            component.source_faces.push_back(source_faces[face_id]);
        }
        for (int old : touched) remap[old] = -1;
        component.mesh.properties.assign(component.mesh.indices.size(), FaceProperty{eNormal, 0.0});
        inspect(component.mesh, work); // single closed boundary, including vertex links and CGAL
        components.push_back(std::move(component));
    }
    require(!components.empty(), "The cut has no upper component to associate with the guides.");
    return components;
}

size_t resolve_guided_component(const indexed_triangle_set &source,
    const CutAndStretchSourceFrame &frame, const std::vector<GuidedComponent> &components,
    const std::vector<int> &guide_facets, const Work &work)
{
    require(!guide_facets.empty() && guide_facets.size() <= source.indices.size(),
            "Supply the nonempty original facet mask of every selected guide.");
    std::vector<int> guides = guide_facets;
    std::sort(guides.begin(), guides.end());
    require(guides.front() >= 0 && size_t(guides.back()) < source.indices.size(),
            "A guide facet no longer exists in the original source.");
    require(std::adjacent_find(guides.begin(), guides.end()) == guides.end(),
            "The guided facet mask contains duplicate original IDs; supply their unique union.");
    work.check();
    std::vector<int> source_component(source.indices.size(), -1);
    for (size_t component = 0; component < components.size(); ++component) {
        const auto &ids = components[component].source_faces;
        for (size_t i = 0; i < ids.size(); ++i) {
            if ((i & 1023) == 0) work.check();
            const int id = ids[i];
            if (id == -1) continue; // generated cap, never an original guide
            require(id >= 0 && size_t(id) < source.indices.size(), "An upper fragment has invalid source provenance.");
            require(source_component[id] < 0 || source_component[id] == int(component),
                    "Fragments of one original facet belong to different upper components.");
            source_component[id] = int(component);
        }
    }
    std::vector<Vec3i> source_neighbors;
    its_edge_diagnostics(source, &source_neighbors); // source was fully inspected before the native cut
    const auto plane = frame.exact_plane();
    std::vector<signed char> coplanar(source.indices.size(), -1);
    const auto is_coplanar = [&](int face) {
        signed char &cached = coplanar[face];
        if (cached < 0) {
            cached = 1;
            for (int vertex : source.indices[face]) {
                const Vec3d p = source.vertices[vertex].cast<double>();
                if (!plane.has_on(ExactKernel::Point_3(p.x(), p.y(), p.z()))) { cached = 0; break; }
            }
        }
        return cached != 0;
    };
    // -2 unknown; -4 queued; -1 no owner; -3 multiple owners; >=0 unique owner.
    std::vector<int> coplanar_owner(source.indices.size(), -2);
    int selected = -1;
    for (size_t i = 0; i < guides.size(); ++i) {
        if ((i & 1023) == 0) work.check();
        const int guide = guides[i];
        int owner = source_component[guide];
        if (owner < 0) {
            require(is_coplanar(guide),
                    "A selected guide facet is wholly fixed or has no moving fragment; restrict the guide mask or choose another plane.");
            // Adjacency alone could associate a stationary shelf with a column
            // rising beside it. Material is OPPOSITE the outward facet normal;
            // certify that this original coplanar guide really belongs to the
            // moving material side, including mirrored/nonuniform placement.
            const auto point = [&](int vertex) {
                const Vec3d p = source.vertices[vertex].cast<double>();
                return ExactKernel::Point_3(p.x(), p.y(), p.z());
            };
            const Vec3i &facet = source.indices[guide];
            const ExactKernel::Plane_3 support(point(facet[0]), point(facet[1]), point(facet[2]));
            const Vec3d positive_direction = frame.cut_to_source.col(2);
            const ExactKernel::FT outward_dot = support.a() * ExactKernel::FT(positive_direction.x()) +
                support.b() * ExactKernel::FT(positive_direction.y()) + support.c() * ExactKernel::FT(positive_direction.z());
            require(outward_dot < 0,
                    "A coplanar guide is on the fixed material side; restrict the mask or choose a moving-side guide.");
            if (coplanar_owner[guide] == -2) {
                std::set<int> owners;
                std::vector<int> patch{guide};
                coplanar_owner[guide] = -4;
                for (size_t head = 0; head < patch.size(); ++head) {
                    if ((head & 1023) == 0) work.check();
                    for (int adjacent : source_neighbors[patch[head]]) {
                        require(adjacent >= 0 && size_t(adjacent) < source.indices.size(),
                                "The original coplanar guide patch has an invalid adjacency.");
                        if (is_coplanar(adjacent)) {
                            if (coplanar_owner[adjacent] == -2) {
                                coplanar_owner[adjacent] = -4;
                                patch.push_back(adjacent);
                            }
                        } else if (source_component[adjacent] >= 0) {
                            owners.insert(source_component[adjacent]);
                        }
                    }
                }
                const int resolved = owners.size() == 1 ? *owners.begin() : owners.empty() ? -1 : -3;
                for (int face : patch) coplanar_owner[face] = resolved;
            }
            owner = coplanar_owner[guide];
            require(owner >= 0,
                    "A coplanar guide has no unique upper component; add a guide away from the plane or restrict the mask.");
        }
        require(selected < 0 || selected == owner,
                "The selected guide facets identify different upper components; select one moving component.");
        selected = owner;
    }
    require(selected >= 0 && size_t(selected) < components.size(), "The guide mask does not identify a moving component.");
    return size_t(selected);
}

indexed_triangle_set select_guided_contact(const indexed_triangle_set &lower_cap,
    const ExactSet &selected_footprint, const Work &work)
{
    indexed_triangle_set contact;
    contact.vertices = lower_cap.vertices; // FaceExtrusion copies only referenced vertices, exactly
    for (size_t i = 0; i < lower_cap.indices.size(); ++i) {
        if ((i & 63) == 0) work.check();
        const Vec3i &face = lower_cap.indices[i];
        ExactPolygon triangle;
        for (int index : face) {
            const Vec3f &p = lower_cap.vertices[index];
            triangle.push_back(ExactKernel::Point_2(double(p.x()), double(p.y())));
        }
        require(triangle.orientation() != CGAL::COLLINEAR, "The stationary contact cap has a degenerate triangle.");
        if (triangle.orientation() == CGAL::CLOCKWISE) triangle.reverse_orientation();
        ExactSet triangle_set;
        triangle_set.join(triangle);
        ExactSet overlap = triangle_set;
        overlap.intersection(selected_footprint);
        if (overlap.is_empty()) continue; // boundary-only contact has no material area
        ExactSet outside = triangle_set;
        outside.difference(selected_footprint);
        require(outside.is_empty(),
                "A stationary cap triangle partially overlaps the guided component; this interface cannot be selected without retriangulation.");
        contact.indices.push_back(face); // entire ORIGINAL cap triangle; no clipping, snapping or epsilon
    }
    require(!contact.indices.empty(), "The guided component has no positive-area contact with the stationary half.");
    contact.properties.assign(contact.indices.size(), FaceProperty{eNormal, 0.0});
    inspect(contact, work, false);
    return contact;
}

bool same_selection(const CutAndStretchSelectionSummary &a, const CutAndStretchSelectionSummary &b)
{
    return a.upper_component_count == b.upper_component_count &&
           a.selected_upper_component == b.selected_upper_component &&
           a.guide_facet_count == b.guide_facet_count &&
           a.selected_source_facet_count == b.selected_source_facet_count;
}

void require_strict_disjoint(const indexed_triangle_set &a, const SolidInfo &a_info,
    const indexed_triangle_set &b, const SolidInfo &b_info, const Work &work)
{
    work.check();
    for (int axis = 0; axis < 3; ++axis)
        if (a_info.max[axis] < b_info.min[axis] || b_info.max[axis] < a_info.min[axis]) return;
    require(b.indices.size() <= max_faces && a.indices.size() <= max_faces - b.indices.size() &&
            b.vertices.size() <= size_t(std::numeric_limits<int>::max()) &&
            a.vertices.size() <= size_t(std::numeric_limits<int>::max()) - b.vertices.size(),
            "Collision validation of the guided component exceeds mesh limits.");
    // Aggregation is ONLY an intersection diagnostic, NEVER a delivered mesh or
    // a replacement for solid union. Separate indices let CGAL detect any touch
    // between the two operands, including coplanar/edge/vertex contact.
    indexed_triangle_set diagnostic = a;
    its_merge(diagnostic, b);
    auto combined = MeshBoolean::cgal::triangle_mesh_to_cgal(diagnostic);
    work.check();
    require(!MeshBoolean::cgal::does_self_intersect(*combined),
            "The guided moving component or its bridge touches an excluded upper component; nothing was changed.");
    work.check();
    // Surface intersection alone does not detect one solid contained in another.
    auto left = MeshBoolean::cgal::triangle_mesh_to_cgal(a);
    auto right = MeshBoolean::cgal::triangle_mesh_to_cgal(b);
    work.check();
    MeshBoolean::cgal::intersect(*left, *right); // checked backend: failure throws, not 'empty'
    work.check();
    require(MeshBoolean::cgal::empty(*left),
            "The guided moving component or its bridge overlaps an excluded upper component; nothing was changed.");
}

} // namespace

static CutAndStretchPlan prepare_cut_and_stretch_impl(const indexed_triangle_set &source,
    const Vec3d &plane_origin, const Vec3d &moving_normal,
    const std::function<bool()> &cancel, const std::function<void(int)> &progress,
    const Transform3d *source_to_world, const size_t *plane_facet,
    const std::vector<int> *guide_facets = nullptr)
{
    CutAndStretchPlan plan;
    const Work work{cancel, progress};
    const char *stage = "Source validation";
    try {
        work.report(0);
        require(plane_origin.allFinite() && moving_normal.allFinite() && std::isfinite(moving_normal.norm()) && moving_normal.norm() > 1e-12,
                "Choose a finite cut plane and a nonzero moving-side normal.");
        SolidInfo source_info = inspect(source, work);
        Vec3d n = moving_normal.normalized();
        Eigen::Quaterniond q;
        plan.cut_to_world.linear() = q.setFromTwoVectors(Vec3d::UnitZ(), n).toRotationMatrix();
        plan.cut_to_world.translation() = plane_origin;
        stage = "Cut-frame conversion";
        indexed_triangle_set cut_source;
        double input_budget;
        if (!source_to_world) {
            // Preserve the legacy WORLD-float contract and its original guard.
            input_budget = rounding_budget(source_info);
            cut_source = transformed(source, plan.cut_to_world.inverse(), input_budget, work);
        } else {
            const Transform3d &placement = *source_to_world;
            require(placement.matrix().allFinite() && placement.matrix().row(3) == Eigen::RowVector4d(0, 0, 0, 1),
                    "The source placement must be a finite affine transformation.");
            const Matrix3d linear = placement.linear();
            const double determinant = linear.determinant();
            require(std::isfinite(determinant) && determinant != 0.0 && linear.fullPivLu().isInvertible(),
                    "The source placement is singular or outside the supported numeric range.");
            const Matrix3d inverse = linear.inverse();
            require(inverse.allFinite() && std::isfinite(linear.norm() * inverse.norm()) && linear.norm() * inverse.norm() <= 1e12,
                    "The source placement is too ill-conditioned for a precise physical cut.");
            auto frame = std::make_shared<CutAndStretchSourceFrame>();
            frame->source_anchor = 0.5 * source_info.min + 0.5 * source_info.max;
            frame->original_placement = placement;
            frame->requested_origin = plane_origin;
            frame->requested_normal = moving_normal;
            if (plane_facet) {
                set_source_facet(source, *plane_facet, inverse, n, *frame);
                frame->requested_origin = placement * frame->source_plane[0];
            }
            // Follow a stable source tangent instead of introducing an unrelated
            // world-XY rotation before storing the physical cut mesh as float.
            Vec3d tangent = Vec3d::Zero();
            for (int axis = 0; axis < 3; ++axis) {
                const Vec3d projected = linear.col(axis) - n.dot(linear.col(axis)) * n;
                if (projected.squaredNorm() > tangent.squaredNorm()) tangent = projected;
            }
            require(tangent.allFinite() && tangent.norm() > 0.0, "The source has no stable physical cut-plane tangent.");
            plan.cut_to_world.linear().col(0) = tangent.normalized();
            plan.cut_to_world.linear().col(1) = n.cross(plan.cut_to_world.linear().col(0)).normalized();
            plan.cut_to_world.linear().col(2) = n;
            frame->source_to_cut = plan.cut_to_world.linear().transpose() * linear;
            // Algebraically n_world^T A = n_source^T / |A^-T n_source|.
            // Use that covector directly instead of reintroducing cancellation
            // into exact source-facet distances through matrix products.
            if (frame->has_source_plane) frame->source_to_cut.row(2) = frame->source_plane_metric.transpose();
            frame->cut_to_source = frame->source_to_cut.inverse();
            frame->determinant = determinant;
            const Vec3d anchor_relative_to_plane = linear * frame->source_anchor + (placement.translation() - frame->requested_origin);
            require(anchor_relative_to_plane.allFinite(), "The plane and source placement exceed numeric limits.");
            frame->cut_anchor.z() = frame->has_source_plane ?
                frame->source_to_cut.row(2).dot(frame->source_anchor - frame->source_plane[0]) : n.dot(anchor_relative_to_plane);
            // Shift ONLY the coordinate origin tangentially onto the same plane.
            // All evaluations below are anchored; no source vertex ever becomes
            // a world float, and no model/instance placement is changed.
            plan.cut_to_world.translation() = placement.translation() + linear * frame->source_anchor - n * frame->cut_anchor.z();
            require(plan.cut_to_world.matrix().allFinite() && frame->cut_to_source.allFinite(),
                    "The near-source cut frame exceeds numeric limits.");
            frame->cut_to_world = plan.cut_to_world;
            input_budget = rounding_budget(source_frame_bounds(source_info, *frame));
            source_info.area = physical_area(source, frame->source_to_cut, work);
            source_info.volume *= std::abs(determinant);
            require(std::isfinite(source_info.volume) && source_info.volume > 0.0,
                    "The physical source volume is not representable.");
            cut_source = source_frame_transform(source, *frame, false, input_budget, work);
            plan.source_frame = std::move(frame);
        }
        const SolidInfo cut_info = inspect(cut_source, work);
        const double budget = std::max(input_budget, rounding_budget(cut_info));
        volume_matches(cut_info.volume, source_info.volume, cut_info.area + source_info.area, budget,
                       "The plane transformation changes the source volume beyond its precision budget.");
        require(cut_info.min.z() < -cut_band && cut_info.max.z() > cut_band,
                "The cut plane must cross the solid, not merely touch an exterior face.");
        guard_cut_band(cut_source, budget, work);
        work.report(20);
        std::vector<int> upper_source, lower_source;
        // cut_mesh owns a raw seam map until it returns. Never unwind through
        // a user callback there: request cooperative cancellation and defer any
        // callback exception until the native routine has released its map.
        std::exception_ptr callback_error;
        bool native_canceled = false;
        auto native_cancel = [&]() {
            try { if (cancel && cancel()) native_canceled = true; }
            catch (...) { callback_error = std::current_exception(); native_canceled = true; }
            return native_canceled;
        };
        cut_mesh(cut_source, 0.0f, &plan.moving, &plan.stationary, false, &upper_source, &lower_source,
                 [&](int p) {
                     if (native_cancel()) return;
                     try { if (progress) progress(20 + std::clamp(p, 0, 100) / 4); }
                     catch (...) { callback_error = std::current_exception(); native_canceled = true; }
                     native_cancel();
                 }, native_cancel);
        if (callback_error) std::rethrow_exception(callback_error);
        if (native_canceled) throw Canceled{};
        work.check();
        stage = "Stationary cut cap";
        cap_cut_boundary(plan.stationary, lower_source, false, budget, work);
        const SolidInfo lower_info = inspect(plan.stationary, work);
        stage = "Moving cut cap";
        cap_cut_boundary(plan.moving, upper_source, true, budget, work);
        // Legacy still requires one whole moving half. Guided mode validates the
        // closed upper collection first, then EACH positive component separately.
        const SolidInfo upper_info = inspect(plan.moving, work, true, guide_facets == nullptr);
        require(lower_info.max.z() <= budget && upper_info.min.z() >= -budget,
                "The cut halves do not respect the selected stationary and moving sides.");
        volume_matches(lower_info.volume + upper_info.volume, source_info.volume,
                       lower_info.area + upper_info.area + source_info.area, budget,
                       "The native cut failed solid-volume conservation; the source is unchanged.");
        stage = "Contact section";
        auto lower_cap = cap_from_ids(plan.stationary, lower_source, budget, work);
        auto upper_cap = cap_from_ids(plan.moving, upper_source, budget, work);
        work.report(55);
        if (guide_facets == nullptr) {
            ExactSet uncovered = cap_footprint(lower_cap, work);
            const ExactSet upper_set = cap_footprint(upper_cap, work);
            uncovered.difference(upper_set);
            work.check();
            require(uncovered.is_empty(),
                    "The stationary cap contains a coplanar shelf outside the moving cap. This contact shape is not supported; choose another cut plane.");
            // Legacy: proven B subset A => actual contact A intersect B is B.
        } else {
            require(source_to_world != nullptr && bool(plan.source_frame),
                    "Guided preparation requires the unchanged local source and its placement.");
            stage = "Guided component selection";
            auto components = split_guided_upper(plan.moving, upper_source, work);
            const size_t chosen = resolve_guided_component(source, *plan.source_frame, components, *guide_facets, work);
            stage = "Guided contact section";
            const auto chosen_cap = cap_from_ids(components[chosen].mesh, components[chosen].source_faces, budget, work);
            lower_cap = select_guided_contact(lower_cap, cap_footprint(chosen_cap, work), work);
            auto guided = std::make_shared<CutAndStretchGuidedState>();
            guided->selection.upper_component_count = components.size();
            guided->selection.selected_upper_component = chosen;
            guided->selection.guide_facet_count = guide_facets->size();
            std::set<int> selected_source_ids;
            for (size_t i = 0; i < components[chosen].source_faces.size(); ++i) {
                if ((i & 1023) == 0) work.check();
                const int id = components[chosen].source_faces[i];
                if (id >= 0) selected_source_ids.insert(id);
            }
            guided->selection.selected_source_facet_count = selected_source_ids.size();
            plan.moving = std::move(components[chosen].mesh);
            guided->excluded_upper.reserve(components.size() - 1);
            stage = "Fixed remainder union";
            size_t completed = 0;
            for (size_t i = 0; i < components.size(); ++i) {
                if (i == chosen) continue;
                work.check();
                const int begin = 55 + int(25 * completed / (components.size() - 1));
                const int end = 55 + int(25 * (completed + 1) / (components.size() - 1));
                // A REAL Boolean cancels the internal caps. Never concatenate
                // closed lower + upper shells, and never move the excluded body.
                plan.stationary = union_solids(plan.stationary, components[i].mesh, work, begin, end);
                guided->excluded_upper.push_back(std::move(components[i].mesh));
                ++completed;
            }
            const SolidInfo stationary_info = inspect(plan.stationary, work);
            const SolidInfo chosen_info = inspect(plan.moving, work);
            volume_matches(stationary_info.volume + chosen_info.volume, source_info.volume,
                           stationary_info.area + chosen_info.area + source_info.area, budget,
                           "Reconstructing the fixed remainder changed the source volume.");
            plan.scope = CutAndStretchScope::GuidedComponent;
            plan.selection = guided->selection;
            plan.guided_state = std::move(guided);
            stage = "Guided contact bridge";
        }
        auto bridge = extrude(lower_cap, 1.0);
        require(bool(bridge), bridge.message.c_str());
        plan.bridge_unit = std::move(bridge.mesh);
        const SolidInfo bridge_info = inspect(plan.bridge_unit, work);
        for (const Vec3i &f : lower_cap.indices) {
            const Vec3d a = lower_cap.vertices[f[0]].cast<double>();
            const Vec3d b = lower_cap.vertices[f[1]].cast<double>();
            const Vec3d c = lower_cap.vertices[f[2]].cast<double>();
            const double projected_area = 0.5 * (b - a).cross(c - a).z();
            require(projected_area > 0.0, "The contact cap has an inconsistent projected orientation.");
            plan.contact_area_mm2 += projected_area;
        }
        require(std::isfinite(plan.contact_area_mm2) && plan.contact_area_mm2 > 0.0, "The cut has no positive-area contact section.");
        volume_matches(bridge_info.volume, plan.contact_area_mm2, bridge_info.area, budget,
                       "The unit bridge does not preserve the contact section area.");
        plan.source_volume_mm3 = source_info.volume;
        plan.stationary.properties.assign(plan.stationary.indices.size(), FaceProperty{eNormal, 0.0});
        plan.moving.properties.assign(plan.moving.indices.size(), FaceProperty{eNormal, 0.0});
        work.report(100);
        plan.valid = true;
        return plan;
    } catch (const Canceled &) {
        plan = CutAndStretchPlan{};
        plan.message = "Cut and extend preparation was canceled. The source is unchanged.";
    } catch (const std::exception &e) {
        plan = CutAndStretchPlan{};
        plan.message = std::string(stage) + ": " + e.what();
    } catch (...) {
        plan = CutAndStretchPlan{};
        plan.message = "Unable to prepare a safe cut and extension. The source is unchanged.";
    }
    return plan;
}

CutAndStretchPlan prepare_cut_and_stretch(const indexed_triangle_set &source_world,
    const Vec3d &plane_origin, const Vec3d &moving_normal,
    const std::function<bool()> &cancel, const std::function<void(int)> &progress)
{
    return prepare_cut_and_stretch_impl(source_world, plane_origin, moving_normal, cancel, progress, nullptr, nullptr);
}

CutAndStretchPlan prepare_cut_and_stretch(const indexed_triangle_set &source_local, const Transform3d &source_to_world,
    const Vec3d &plane_origin_world, const Vec3d &moving_normal_world,
    const std::function<bool()> &cancel, const std::function<void(int)> &progress)
{
    return prepare_cut_and_stretch_impl(source_local, plane_origin_world, moving_normal_world, cancel, progress, &source_to_world, nullptr);
}

CutAndStretchPlan prepare_cut_and_stretch_from_face(const indexed_triangle_set &source_local, const Transform3d &source_to_world,
    size_t plane_facet, const Vec3d &moving_normal_world,
    const std::function<bool()> &cancel, const std::function<void(int)> &progress)
{
    return prepare_cut_and_stretch_impl(source_local, Vec3d::Zero(), moving_normal_world, cancel, progress, &source_to_world, &plane_facet);
}

CutAndStretchPlan prepare_cut_and_stretch_guided(
    const indexed_triangle_set &source_local, const Transform3d &source_to_world,
    const Vec3d &plane_origin_world, const Vec3d &moving_normal_world,
    const std::vector<int> &guide_facets,
    const std::function<bool()> &cancel, const std::function<void(int)> &progress)
{
    return prepare_cut_and_stretch_impl(source_local, plane_origin_world, moving_normal_world,
                                       cancel, progress, &source_to_world, nullptr, &guide_facets);
}

CutAndStretchPlan prepare_cut_and_stretch_guided_from_face(
    const indexed_triangle_set &source_local, const Transform3d &source_to_world,
    size_t plane_facet, const Vec3d &moving_normal_world,
    const std::vector<int> &guide_facets,
    const std::function<bool()> &cancel, const std::function<void(int)> &progress)
{
    return prepare_cut_and_stretch_impl(source_local, Vec3d::Zero(), moving_normal_world,
                                       cancel, progress, &source_to_world, &plane_facet, &guide_facets);
}

static CutAndStretchResult apply_cut_and_stretch_impl(const CutAndStretchPlan &plan, double distance_mm,
    const std::function<bool()> &cancel, const std::function<void(int)> &progress, bool to_source)
{
    CutAndStretchResult result;
    const Work work{cancel, progress};
    try {
        work.report(0);
        require(bool(plan), "Prepare a valid cut plane before extending the solid.");
        require(plan.scope == CutAndStretchScope::WholeHalfSpace || plan.scope == CutAndStretchScope::GuidedComponent,
                "The prepared cut scope is invalid.");
        const bool guided = plan.scope == CutAndStretchScope::GuidedComponent;
        require(guided ? bool(plan.guided_state) && bool(plan.source_frame) : !plan.guided_state,
                "The prepared cut scope does not match its immutable provenance.");
        if (guided) {
            require(same_selection(plan.selection, plan.guided_state->selection) &&
                    plan.selection.upper_component_count > 0 &&
                    plan.selection.upper_component_count <= max_guided_components &&
                    plan.selection.selected_upper_component < plan.selection.upper_component_count &&
                    plan.guided_state->excluded_upper.size() + 1 == plan.selection.upper_component_count,
                    "The guided component selection changed after preparation.");
        }
        require(std::isfinite(distance_mm) && distance_mm > 0.0,
                "Only positive extension is supported. Zero makes no edit; shortening requires a separate remove-band operation.");
        require(plan.cut_to_world.matrix().allFinite() &&
                (plan.cut_to_world.linear().transpose() * plan.cut_to_world.linear() - Matrix3d::Identity()).norm() <= 1e-10 &&
                std::abs(plan.cut_to_world.linear().determinant() - 1.0) <= 1e-10,
                "The prepared cut frame is not a valid rigid transformation.");
        require(std::isfinite(plan.contact_area_mm2) && plan.contact_area_mm2 > 0.0 &&
                std::isfinite(plan.source_volume_mm3) && plan.source_volume_mm3 > 0.0,
                "The prepared contact area or source volume is invalid.");
        if (to_source) {
            require(bool(plan.source_frame), "Prepare the cut from the original local mesh and its placement before requesting a local result.");
            require((plan.cut_to_world.matrix().array() == plan.source_frame->cut_to_world.matrix().array()).all(),
                    "The prepared world preview frame changed after source-local preparation.");
        }
        const SolidInfo lower_info = inspect(plan.stationary, work);
        const SolidInfo upper_info = inspect(plan.moving, work);
        const SolidInfo unit_info = inspect(plan.bridge_unit, work);
        double budget = std::max({rounding_budget(lower_info), rounding_budget(upper_info), rounding_budget(unit_info)});
        if (!to_source)
            budget = std::max({budget, rounding_budget(transformed_bounds(lower_info, plan.cut_to_world)),
                              rounding_budget(transformed_bounds(upper_info, plan.cut_to_world))});
        require(distance_mm > 4.0 * budget, "The extension is too small to represent without disturbing the cut seam.");
        require((guided || lower_info.max.z() <= budget) && upper_info.min.z() >= -budget &&
                std::abs(unit_info.min.z()) <= budget && std::abs(unit_info.max.z() - 1.0) <= budget,
                "The prepared operands no longer match the cut plane.");
        volume_matches(lower_info.volume + upper_info.volume, plan.source_volume_mm3,
                       lower_info.area + upper_info.area, budget, "The prepared source volume changed.");
        volume_matches(unit_info.volume, plan.contact_area_mm2, unit_info.area, budget, "The prepared contact section changed.");
        const auto cap = base_cap(plan.bridge_unit, work);
        auto bridge = extrude(cap, distance_mm);
        require(bool(bridge), bridge.message.c_str());
        const SolidInfo bridge_info = inspect(bridge.mesh, work);
        Transform3d translation = Transform3d::Identity();
        translation.translation().z() = distance_mm;
        budget = std::max(budget, rounding_budget(transformed_bounds(upper_info, translation)));
        auto moving = transformed(plan.moving, translation, budget, work);
        const SolidInfo moved_info = inspect(moving, work);
        budget = std::max({budget, rounding_budget(moved_info), rounding_budget(bridge_info)});
        volume_matches(moved_info.volume, upper_info.volume, moved_info.area + upper_info.area, budget,
                       "Translation changed the moving part's volume beyond the precision budget.");
        if (guided) {
            for (const auto &excluded : plan.guided_state->excluded_upper) {
                const SolidInfo excluded_info = inspect(excluded, work);
                // The global lower half cannot intersect material strictly above
                // its cap. Excluded UPPER components can: reject all new touches,
                // overlaps and containment before any candidate union is built.
                require_strict_disjoint(moving, moved_info, excluded, excluded_info, work);
                require_strict_disjoint(bridge.mesh, bridge_info, excluded, excluded_info, work);
            }
        }
        work.report(20);
        auto candidate = union_solids(plan.stationary, bridge.mesh, work, 20, 50);
        candidate = union_solids(candidate, moving, work, 50, 80);
        const SolidInfo candidate_info = inspect(candidate, work);
        const double expected = plan.source_volume_mm3 + plan.contact_area_mm2 * distance_mm;
        volume_matches(candidate_info.volume, expected, candidate_info.area + lower_info.area + moved_info.area + bridge_info.area,
                       budget, "The joined solid failed the expected extension-volume check.");
        indexed_triangle_set output;
        if (to_source) {
            budget = std::max(budget, rounding_budget(candidate_info));
            const auto &frame = *plan.source_frame;
            output = source_frame_transform(candidate, frame, true, budget, work);
            const SolidInfo final_info = inspect(output, work); // positive SOURCE-LOCAL winding, including mirrors
            const double area_mm2 = physical_area(output, frame.source_to_cut, work);
            volume_matches(final_info.volume * std::abs(frame.determinant), expected, area_mm2 + candidate_info.area, budget,
                           "Source-local conversion failed the physical extension-volume check.");
        } else {
            budget = std::max(budget, rounding_budget(transformed_bounds(candidate_info, plan.cut_to_world)));
            output = transformed(candidate, plan.cut_to_world, budget, work);
            const SolidInfo final_info = inspect(output, work);
            budget = std::max(budget, rounding_budget(final_info));
            volume_matches(final_info.volume, expected, final_info.area + candidate_info.area, budget,
                           "World-coordinate conversion failed the extension-volume check.");
        }
        work.report(100);
        result.mesh = std::move(output);
        result.valid = true;
    } catch (const Canceled &) {
        result = CutAndStretchResult{};
        result.message = "Cut and extend was canceled. The source is unchanged.";
    } catch (const std::exception &e) {
        result = CutAndStretchResult{};
        result.message = e.what();
    } catch (...) {
        result = CutAndStretchResult{};
        result.message = "Unable to join the cut and extended solid safely. The source is unchanged.";
    }
    return result;
}

CutAndStretchResult apply_cut_and_stretch(const CutAndStretchPlan &plan, double distance_mm,
    const std::function<bool()> &cancel, const std::function<void(int)> &progress)
{
    return apply_cut_and_stretch_impl(plan, distance_mm, cancel, progress, false);
}

CutAndStretchResult apply_cut_and_stretch_to_source(const CutAndStretchPlan &plan, double distance_mm,
    const std::function<bool()> &cancel, const std::function<void(int)> &progress)
{
    return apply_cut_and_stretch_impl(plan, distance_mm, cancel, progress, true);
}

} // namespace Slic3r
