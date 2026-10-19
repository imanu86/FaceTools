// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include "MeshSolidValidation.hpp"
#include "MeshBoolean.hpp"
#include "MeshDiagnostics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Slic3r {

void validate_solid_mesh(const indexed_triangle_set &mesh, const std::function<void()> &check_cancel)
{
    const auto checkpoint = [&] { if (check_cancel) check_cancel(); };
    checkpoint();
    if (mesh.vertices.empty() || mesh.indices.empty())
        throw std::runtime_error("Mesh does not contain a non-empty solid.");
    if (mesh.vertices.size() > size_t(std::numeric_limits<int>::max()) ||
        mesh.indices.size() > size_t(std::numeric_limits<int>::max()))
        throw std::runtime_error("Mesh exceeds the indexed mesh size limit.");

    Vec3d minimum = Vec3d::Constant(std::numeric_limits<double>::infinity());
    Vec3d maximum = Vec3d::Constant(-std::numeric_limits<double>::infinity());
    for (size_t i = 0; i < mesh.vertices.size(); ++i) {
        if ((i & 1023) == 0) checkpoint();
        if (!mesh.vertices[i].allFinite())
            throw std::runtime_error("Mesh has a non-finite vertex.");
        const Vec3d vertex = mesh.vertices[i].cast<double>();
        minimum = minimum.cwiseMin(vertex);
        maximum = maximum.cwiseMax(vertex);
    }

    // Index and arithmetic checks precede every consumer that assumes valid
    // triangles. Canonical face keys also reject opposite-winding duplicates;
    // they never change the source face order, winding or properties.
    std::vector<unsigned char> referenced(mesh.vertices.size(), 0);
    std::vector<std::array<int, 3>> face_keys;
    face_keys.reserve(mesh.indices.size());
    for (size_t i = 0; i < mesh.indices.size(); ++i) {
        if ((i & 1023) == 0) checkpoint();
        const Vec3i &face = mesh.indices[i];
        for (int index : {face[0], face[1], face[2]}) {
            if (index < 0 || size_t(index) >= mesh.vertices.size())
                throw std::runtime_error("Mesh has an invalid triangle index.");
            referenced[size_t(index)] = 1;
        }
        if (face[0] == face[1] || face[1] == face[2] || face[2] == face[0])
            throw std::runtime_error("Mesh has a degenerate triangle.");
        const Vec3d a = mesh.vertices[face[0]].cast<double>();
        const Vec3d b = mesh.vertices[face[1]].cast<double>();
        const Vec3d c = mesh.vertices[face[2]].cast<double>();
        const double area = (b - a).cross(c - a).norm() * 0.5;
        const Vec3f float_cross = (mesh.vertices[face[1]] - mesh.vertices[face[0]]).cross(
                                  mesh.vertices[face[2]] - mesh.vertices[face[0]]);
        const float float_squared_norm = float_cross.squaredNorm();
        if (!std::isfinite(area) || area <= 0.0 ||
            !std::isfinite(float_squared_norm) || float_squared_norm <= 0.f)
            throw std::runtime_error("Mesh has a degenerate or unrepresentable triangle.");
        std::array<int, 3> key{face[0], face[1], face[2]};
        std::sort(key.begin(), key.end());
        face_keys.push_back(key);
    }
    for (size_t i = 0; i < referenced.size(); ++i) {
        if ((i & 1023) == 0) checkpoint();
        if (!referenced[i])
            throw std::runtime_error("Mesh contains an unused vertex.");
    }
    checkpoint();
    std::sort(face_keys.begin(), face_keys.end());
    checkpoint();
    if (std::adjacent_find(face_keys.begin(), face_keys.end()) != face_keys.end())
        throw std::runtime_error("Mesh contains duplicate triangles, including possible opposite winding.");
    if (!minimum.allFinite() || !maximum.allFinite())
        throw std::runtime_error("Mesh has invalid bounds.");

    // Preserve the GUI's translation-invariant compensated volume check. Only
    // the total solid volume must be positive: inward cavity shells stay inward.
    const Vec3d origin = 0.5 * (minimum + maximum);
    double volume = 0.0, compensation = 0.0;
    for (size_t i = 0; i < mesh.indices.size(); ++i) {
        if ((i & 1023) == 0) checkpoint();
        const Vec3i &face = mesh.indices[i];
        const Vec3d a = mesh.vertices[face[0]].cast<double>() - origin;
        const Vec3d b = mesh.vertices[face[1]].cast<double>() - origin;
        const Vec3d c = mesh.vertices[face[2]].cast<double>() - origin;
        const double term = a.dot(b.cross(c)) / 6.0 - compensation;
        const double sum = volume + term;
        compensation = (sum - volume) - term;
        volume = sum;
    }
    checkpoint();
    const float stored_volume = its_volume(mesh);
    checkpoint();
    if (!std::isfinite(volume) || volume <= 0.0 || volume > std::numeric_limits<float>::max() ||
        !std::isfinite(stored_volume) || stored_volume <= 0.f)
        throw std::runtime_error("Mesh has no representable positive solid volume.");

    const MeshDiagnosticStats diagnostics = its_mesh_diagnostics(mesh);
    checkpoint();
    if (diagnostics.open_edges || diagnostics.non_manifold_edges || diagnostics.non_manifold_vertices ||
        diagnostics.same_direction_edges || diagnostics.has_reversed_faces)
        throw std::runtime_error("Mesh is not closed, manifold and consistently oriented.");
    auto solid = MeshBoolean::cgal::triangle_mesh_to_cgal(mesh);
    checkpoint();
    const bool self_intersects = MeshBoolean::cgal::does_self_intersect(*solid);
    checkpoint();
    if (self_intersects)
        throw std::runtime_error("Mesh has self-intersections.");
    const bool bounds_volume = MeshBoolean::cgal::does_bound_a_volume(*solid);
    checkpoint();
    if (!bounds_volume)
        throw std::runtime_error("Mesh shells do not bound a consistently oriented solid.");
}

} // namespace Slic3r
