// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include <catch2/catch.hpp>

#include "libslic3r/CutPlaneSelection.hpp"
#include "libslic3r/CutAndStretch.hpp"
#include "libslic3r/PlanarProfile.hpp"
#include "libslic3r/MeshDiagnostics.hpp"
#include "libslic3r/TriangleMeshSlicer.hpp"

#include <algorithm>
#include <limits>

using namespace Slic3r;

namespace {

double signed_volume(const indexed_triangle_set &mesh)
{
    REQUIRE_FALSE(mesh.vertices.empty());
    const Vec3d origin = mesh.vertices.front().cast<double>();
    double volume = 0.0;
    for (const Vec3i &face : mesh.indices) {
        const Vec3d a = mesh.vertices[face[0]].cast<double>() - origin;
        const Vec3d b = mesh.vertices[face[1]].cast<double>() - origin;
        const Vec3d c = mesh.vertices[face[2]].cast<double>() - origin;
        volume += a.dot(b.cross(c)) / 6.0;
    }
    return volume;
}

void unchanged_mesh(const indexed_triangle_set &before, const indexed_triangle_set &after)
{
    REQUIRE(after.vertices.size() == before.vertices.size());
    REQUIRE(after.indices.size() == before.indices.size());
    REQUIRE(after.properties.size() == before.properties.size());
    for (size_t i = 0; i < before.vertices.size(); ++i) CHECK(after.vertices[i] == before.vertices[i]);
    for (size_t i = 0; i < before.indices.size(); ++i) CHECK(after.indices[i] == before.indices[i]);
    for (size_t i = 0; i < before.properties.size(); ++i) {
        CHECK(after.properties[i].type == before.properties[i].type);
        CHECK(after.properties[i].area == before.properties[i].area);
    }
}

void require_closed(const indexed_triangle_set &mesh)
{
    REQUIRE_FALSE(mesh.empty());
    const auto diagnostics = its_mesh_diagnostics(mesh);
    CHECK(diagnostics.open_edges == 0);
    CHECK(diagnostics.non_manifold_edges == 0);
    CHECK(diagnostics.non_manifold_vertices == 0);
    CHECK(diagnostics.same_direction_edges == 0);
    CHECK_FALSE(diagnostics.has_reversed_faces);
}

indexed_triangle_set translated(indexed_triangle_set mesh, const Vec3d &offset)
{
    for (Vec3f &vertex : mesh.vertices) vertex = (vertex.cast<double>() + offset).cast<float>();
    return mesh;
}

} // namespace

TEST_CASE("Face-first L workflow chooses an independent plane and keeps the full arm rigid", "[cut_plane_selection]")
{
    // Same dimensions as the delivered L: 40x10x24, arm thickness 4, volume 3600.
    const auto generated = extrude_planar_profile(
        {Vec3d(0,0,0), Vec3d(10,0,0), Vec3d(10,0,20),
         Vec3d(40,0,20), Vec3d(40,0,24), Vec3d(0,0,24)}, Vec3d(0,10,0));
    INFO(generated.message);
    REQUIRE(bool(generated));
    for (const Vec3d &offset : {Vec3d(0,0,0), Vec3d(128,128,0)}) {
        const auto source = translated(generated.mesh, offset);
        const auto before = source;
        // The FIRST guide is lateral, not Z-normal. It must not orient the cut.
        // The underside guide exactly on the plane identifies material above it.
        const std::vector<CutPlaneFaceGuide> guides = {
            {offset + Vec3d(20,10,22), Vec3d::UnitY()},
            {offset + Vec3d(30,5,24), Vec3d::UnitZ()},
            {offset + Vec3d(30,5,20), -Vec3d::UnitZ()}
        };
        const Vec3d plane_origin = offset + Vec3d(20,5,20);
        for (double sign : {1.0, -1.0}) {
            const auto selection = select_cut_moving_side(guides, plane_origin, sign * Vec3d::UnitZ());
            INFO(selection.message);
            REQUIRE(bool(selection));
            CHECK(selection.side == int(sign));
            CHECK(selection.moving_normal.isApprox(Vec3d::UnitZ(), 1e-12));
            const auto plan = prepare_cut_and_stretch(source, plane_origin, selection.moving_normal);
            INFO(plan.message);
            REQUIRE(bool(plan));
            CHECK(plan.contact_area_mm2 == Approx(100.0));
            CHECK(plan.source_volume_mm3 == Approx(3600.0));
            const auto moving_before = plan.moving;
            const auto stationary_before = plan.stationary;
            const auto result = apply_cut_and_stretch(plan, 10.0);
            INFO(result.message);
            REQUIRE(bool(result));
            require_closed(result.mesh);
            CHECK(signed_volume(result.mesh) == Approx(4600.0).margin(1e-4));
            double low = std::numeric_limits<double>::infinity(), high = -low;
            double bottom = low;
            for (const Vec3f &vertex : result.mesh.vertices) {
                const Vec3d point = vertex.cast<double>() - offset;
                bottom = std::min(bottom, point.z());
                if (point.x() > 15.0) {
                    low = std::min(low, point.z());
                    high = std::max(high, point.z());
                }
            }
            CHECK(bottom == Approx(0.0).margin(1e-5));
            CHECK(low == Approx(30.0).margin(1e-5)); // no old arm left at Z=20
            CHECK(high == Approx(34.0).margin(1e-5));
            CHECK(high - low == Approx(4.0).margin(1e-5));
            unchanged_mesh(moving_before, plan.moving);
            unchanged_mesh(stationary_before, plan.stationary);
            unchanged_mesh(before, source);
        }
    }
}

TEST_CASE("Face-first tube workflow preserves a through hole across the stationary bridge and moving regions", "[cut_plane_selection]")
{
    const std::vector<Vec3d> outer = {{-10,-10,0}, {10,-10,0}, {10,10,0}, {-10,10,0}};
    const std::vector<Vec3d> hole = {{-5,-5,0}, {-5,5,0}, {5,5,0}, {5,-5,0}};
    const auto generated = extrude_planar_profile(outer, Vec3d(0,0,30), {hole});
    INFO(generated.message);
    REQUIRE(bool(generated));
    const Vec3d offset(128,128,0);
    const auto source = translated(generated.mesh, offset);
    const auto before = source;
    const std::vector<CutPlaneFaceGuide> guides = {
        {offset + Vec3d(10,0,25), Vec3d::UnitX()},
        {offset + Vec3d(7.5,0,30), Vec3d::UnitZ()},
        {offset + Vec3d(5,0,25), -Vec3d::UnitX()}
    };
    const Vec3d plane_origin = offset + Vec3d(0,0,15);
    const auto selection = select_cut_moving_side(guides, plane_origin, Vec3d::UnitZ());
    INFO(selection.message);
    REQUIRE(bool(selection));
    CHECK(selection.moving_normal.isApprox(Vec3d::UnitZ(), 1e-12));
    const auto plan = prepare_cut_and_stretch(source, plane_origin, selection.moving_normal);
    INFO(plan.message);
    REQUIRE(bool(plan));
    CHECK(plan.contact_area_mm2 == Approx(300.0));
    const auto result = apply_cut_and_stretch(plan, 10.0);
    INFO(result.message);
    REQUIRE(bool(result));
    require_closed(result.mesh);
    CHECK(signed_volume(result.mesh) == Approx(12000.0).margin(1e-4));
    double top = -std::numeric_limits<double>::infinity();
    for (const Vec3f &vertex : result.mesh.vertices) top = std::max(top, double(vertex.z()));
    CHECK(top == Approx(40.0).margin(1e-5));
    const auto sections = slice_mesh_ex(result.mesh, std::vector<float>{10.f,20.f,35.f});
    REQUIRE(sections.size() == 3);
    for (const auto &section : sections) {
        REQUIRE(section.size() == 1);
        CHECK(section.front().holes.size() == 1);
    }
    unchanged_mesh(before, source);
}

TEST_CASE("Face-first guides cannot silently move a second disconnected branch beyond the cut", "[cut_plane_selection]")
{
    // A connected U below the plane, but two separate upper legs. A guide on
    // only one leg must not silently cause the other upper leg to move as well.
    const auto generated = extrude_planar_profile(
        {Vec3d(0,0,0), Vec3d(40,0,0), Vec3d(40,0,30), Vec3d(30,0,30),
         Vec3d(30,0,5), Vec3d(10,0,5), Vec3d(10,0,30), Vec3d(0,0,30)}, Vec3d(0,10,0));
    INFO(generated.message);
    REQUIRE(bool(generated));
    require_closed(generated.mesh);
    const auto before = generated.mesh;
    const Vec3d origin(20,5,15);
    const auto selection = select_cut_moving_side(
        {{Vec3d(5,5,30), Vec3d::UnitZ()}}, origin, Vec3d::UnitZ());
    INFO(selection.message);
    REQUIRE(bool(selection));
    const auto plan = prepare_cut_and_stretch(generated.mesh, origin, selection.moving_normal);
    INFO(plan.message);
    CHECK_FALSE(bool(plan));
    CHECK(plan.message.find("Disconnected") != std::string::npos);
    CHECK(plan.stationary.empty());
    CHECK(plan.moving.empty());
    CHECK(plan.bridge_unit.empty());
    unchanged_mesh(before, generated.mesh);
}

TEST_CASE("A moving guide face also defines the exact cut plane on later plates", "[cut_plane_selection]")
{
    const auto generated = extrude_planar_profile(
        {Vec3d(0,0,0), Vec3d(10,0,0), Vec3d(10,0,20),
         Vec3d(40,0,20), Vec3d(40,0,24), Vec3d(0,0,24)}, Vec3d(0,10,0));
    INFO(generated.message);
    REQUIRE(bool(generated));
    const auto source_before = generated.mesh;
    size_t underside_facet = generated.mesh.indices.size();
    for (size_t i = 0; i < generated.mesh.indices.size(); ++i) {
        const Vec3i &f = generated.mesh.indices[i];
        const Vec3d a = generated.mesh.vertices[f[0]].cast<double>();
        const Vec3d b = generated.mesh.vertices[f[1]].cast<double>();
        const Vec3d c = generated.mesh.vertices[f[2]].cast<double>();
        if (a.z() == 20.0 && b.z() == 20.0 && c.z() == 20.0 && (b-a).cross(c-a).z() < 0.0) {
            underside_facet = i;
            break;
        }
    }
    REQUIRE(underside_facet < generated.mesh.indices.size());
    const std::vector<Vec3d> offsets = {
        Vec3d(0,0,0), Vec3d(450,100,0), Vec3d(1024,-1024,0), Vec3d(15000,24000,0)
    };
    for (size_t placement = 0; placement < offsets.size(); ++placement) {
        Transform3d world = Transform3d::Identity();
        const double scale_z = placement == 1 ? 1.5 : placement == 2 ? 1.25 : 1.0;
        const Vec3d scales(placement == 2 ? -1.0 : 1.25, 0.8, scale_z);
        world.linear() = Eigen::AngleAxisd(0.17 * double(placement), Vec3d::UnitZ()).toRotationMatrix() * scales.asDiagonal();
        world.translation() = offsets[placement];
        const Transform3d world_before = world;
        const auto guide = [&](const Vec3d &point, const Vec3d &normal) {
            CutPlaneFaceGuide face;
            face.point_world = world * point;
            face.outward_normal_world = (world.linear().inverse().transpose() * normal).normalized();
            return face;
        };
        // The underside is simultaneously a moving guide AND the exact plane.
        // It can be any item in the selection, not just the first guide.
        const auto underside = guide(Vec3d(30,5,20), -Vec3d::UnitZ());
        for (size_t chosen_index = 0; chosen_index < 3; ++chosen_index) {
            INFO("placement=" << placement << ", guide index=" << chosen_index);
            std::vector<CutPlaneFaceGuide> guides = {
                guide(Vec3d(30,5,24), Vec3d::UnitZ()),
                guide(Vec3d(20,10,22), Vec3d::UnitY())
            };
            guides.insert(guides.begin() + chosen_index, underside);
            const auto &chosen = guides[chosen_index];
            const auto selection = select_cut_moving_side(guides, chosen.point_world, chosen.outward_normal_world);
            INFO(selection.message);
            REQUIRE(bool(selection));
            CHECK(selection.side == -1);
            CHECK(selection.moving_normal.isApprox(Vec3d::UnitZ(), 1e-12));
            const auto plan = prepare_cut_and_stretch_from_face(generated.mesh, world,
                underside_facet, selection.moving_normal);
            INFO(plan.message);
            REQUIRE(bool(plan));
            // Preview stays on the original plate, with a frame on the chosen
            // infinite plane. A tangential change of frame origin is immaterial.
            CHECK(std::abs(selection.moving_normal.dot(plan.cut_to_world.translation() - chosen.point_world)) < 1e-9);
            CHECK(plan.cut_to_world.linear().col(2).isApprox(selection.moving_normal, 1e-12));
            const double distance = 10.0;
            const auto result = apply_cut_and_stretch_to_source(plan, distance);
            INFO(result.message);
            REQUIRE(bool(result));
            require_closed(result.mesh);
            const double local_shift = distance / scale_z;
            CHECK(signed_volume(result.mesh) == Approx(3600.0 + 100.0 * local_shift).margin(0.01));
            Vec3d minimum = result.mesh.vertices.front().cast<double>();
            Vec3d maximum = minimum;
            for (const Vec3f &vertex : result.mesh.vertices) {
                minimum = minimum.cwiseMin(vertex.cast<double>());
                maximum = maximum.cwiseMax(vertex.cast<double>());
            }
            CHECK((minimum - Vec3d::Zero()).cwiseAbs().maxCoeff() < 1e-4);
            CHECK((maximum - Vec3d(40,10,24 + local_shift)).cwiseAbs().maxCoeff() < 1e-4);
            for (const Vec3d &corner : {Vec3d(40,0,20), Vec3d(40,10,20), Vec3d(40,0,24), Vec3d(40,10,24)}) {
                const Vec3d expected_local = corner + Vec3d(0,0,local_shift);
                CHECK(std::any_of(result.mesh.vertices.begin(), result.mesh.vertices.end(), [&](const Vec3f &vertex) {
                    return (vertex.cast<double>() - expected_local).norm() < 1e-4;
                }));
                const Vec3d expected_world = world * corner + distance * selection.moving_normal;
                CHECK(std::any_of(plan.moving.vertices.begin(), plan.moving.vertices.end(), [&](const Vec3f &vertex) {
                    const Vec3d rendered = plan.cut_to_world * (vertex.cast<double>() + Vec3d(0,0,distance));
                    return (rendered - expected_world).norm() < 1e-4;
                }));
            }
            unchanged_mesh(source_before, generated.mesh);
            CHECK((world.matrix().array() == world_before.matrix().array()).all());
        }
    }
}
