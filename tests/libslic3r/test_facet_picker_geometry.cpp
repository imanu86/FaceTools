// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include <catch2/catch.hpp>

#include "slic3r/GUI/Gizmos/FacetPicker.hpp"
#include "slic3r/GUI/MeshUtils.hpp"
#include "libslic3r/Model.hpp"

#include <Eigen/LU>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>

using namespace Slic3r;
using Slic3r::GUI::FacetPicker;

namespace {

indexed_triangle_set original_facets()
{
    indexed_triangle_set its;
    its.vertices = {{0.125f, 1.25f, 2.5f}, {17.375f, 3.625f, 7.875f},
                    {-2.25f, 19.5f, 4.25f}, {0.125f, 1.25f, 3.5f}};
    its.indices = {{0, 1, 2}, {0, 2, 3}};
    return its;
}

struct Fixture {
    Model model;
    ModelObject *object = model.add_object();
    ModelVolume *volume = object->add_volume(TriangleMesh(original_facets()), ModelVolumeType::MODEL_PART, false);
};

void unchanged(const indexed_triangle_set &actual, const indexed_triangle_set &before)
{
    REQUIRE(actual.vertices.size() == before.vertices.size());
    REQUIRE(actual.indices.size() == before.indices.size());
    CHECK(std::memcmp(actual.vertices.data(), before.vertices.data(), actual.vertices.size() * sizeof(Vec3f)) == 0);
    CHECK(std::memcmp(actual.indices.data(), before.indices.data(), actual.indices.size() * sizeof(Vec3i)) == 0);
}

void invalid(const FacetPicker::Hit &hit)
{
    CHECK_FALSE(hit.valid());
    CHECK(hit.facet == -1);
    CHECK(hit.volume == nullptr);
    CHECK(hit.mesh == nullptr);
    CHECK(hit.world_normal == Vec3d::Zero());
}

} // namespace

TEST_CASE("Facet picker derives the original facet plane in double across physical poses", "[facet_picker_geometry]")
{
    Fixture fixture;
    const auto mesh = fixture.volume->get_mesh_shared_ptr();
    const indexed_triangle_set before = mesh->its;
    const Vec3d a = before.vertices[0].cast<double>();
    const Vec3d edge1 = before.vertices[1].cast<double>() - a;
    const Vec3d edge2 = before.vertices[2].cast<double>() - a;
    Matrix3d shear;
    shear << 1.0, 0.25, -0.125, 0.0, 1.5, 0.5, 0.0, 0.0, 0.75;
    const Matrix3d rotation = Eigen::AngleAxisd(0.37, Vec3d(1, 2, 3).normalized()).toRotationMatrix();
    const std::array<Matrix3d, 4> linear{{Matrix3d::Identity(), shear,
        rotation * Vec3d(2.0, 0.75, 3.0).asDiagonal(),
        rotation * Vec3d(-2.0, 0.75, 3.0).asDiagonal()}};
    for (const Matrix3d &basis : linear) {
        for (double plate : {0.0, 450.0, 1024.0, 15000.0}) {
            Transform3d world = Transform3d::Identity();
            world.linear() = basis;
            world.translation() = Vec3d(plate + 0.125, -plate + 0.25, 20.5);
            const auto hit = FacetPicker::from_source_facet(*fixture.volume, 0, world);
            REQUIRE(hit.valid());
            CHECK(hit.volume == fixture.volume);
            CHECK(hit.mesh == mesh);
            CHECK(hit.facet == 0);
            CHECK(hit.world_plane_point == world * a);
            CHECK(hit.world_pos == hit.world_plane_point);
            CHECK(hit.world_normal.allFinite());
            CHECK(hit.world_normal.norm() == Approx(1.0).margin(2e-15));
            const Vec3d expected = (basis.inverse().transpose() * edge1.cross(edge2)).normalized();
            CHECK((hit.world_normal - expected).norm() < 2e-15);
            CHECK(std::abs(hit.world_normal.dot(basis * edge1)) < 1e-12);
            CHECK(std::abs(hit.world_normal.dot(basis * edge2)) < 1e-12);
            // Inverse-transpose retains the original outward side under a mirror.
            CHECK(hit.world_normal.dot((basis * edge1).cross(basis * edge2)) * basis.determinant() > 0.0);
            for (int corner : {0, 1, 2}) {
                const Vec3d point = world * before.vertices[corner].cast<double>();
                CHECK(std::abs(hit.world_normal.dot(point - hit.world_plane_point)) < 5e-12);
            }
        }
    }
    CHECK(fixture.volume->get_mesh_shared_ptr() == mesh);
    unchanged(fixture.volume->mesh().its, before);
}

TEST_CASE("Facet picker does not inherit the production raycaster float normal or rounded point", "[facet_picker_geometry]")
{
    Fixture fixture;
    Slic3r::GUI::MeshRaycaster raycaster(fixture.volume->mesh());
    CHECK(&raycaster.get_aabb_mesh().vertices() == &fixture.volume->mesh().its.vertices);
    CHECK(&raycaster.get_aabb_mesh().indices() == &fixture.volume->mesh().its.indices);
    Transform3d world = Transform3d::Identity();
    world.linear() << 1.0, 0.25, -0.125, 0.0, 1.5, 0.5, 0.0, 0.0, 0.75;
    world.translation() = Vec3d(15000.125, -1024.25, 20.5);
    auto hit = FacetPicker::from_source_facet(*fixture.volume, 0, world);
    REQUIRE(hit.valid());
    const Vec3d float_derived = FacetPicker::world_normal(world, raycaster.get_triangle_normal(0).cast<double>());
    // This is the actual production raycaster normal, not a synthetic normal shim.
    REQUIRE((float_derived - hit.world_normal).norm() > 1e-10);
    const Vec3d original_support = hit.world_plane_point;
    const auto original_hit = hit;
    hit.world_pos += Vec3d(0.001, -0.002, 0.003); // Display/ray metadata must not define the plane.
    CHECK(hit.world_plane_point == original_support);
    CHECK(hit.same_geometry(original_hit));
}

TEST_CASE("Facet picker rejects invalid source facets and affine arithmetic without mutation", "[facet_picker_geometry]")
{
    Fixture fixture;
    const Transform3d identity = Transform3d::Identity();
    invalid(FacetPicker::from_source_facet(*fixture.volume, 2, identity));
    invalid(FacetPicker::from_source_facet(*fixture.volume, std::numeric_limits<size_t>::max(), identity));
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (const Vec3d &normal : {Vec3d(0, 0, 0), Vec3d(nan, 1, 0), Vec3d(0, inf, 1)})
        CHECK(FacetPicker::world_normal(identity, normal) == Vec3d::Zero());
    for (int kind = 0; kind < 7; ++kind) {
        Transform3d world = identity;
        if (kind == 0) world.linear().col(2).setZero();
        if (kind == 1) world.linear().col(2) = world.linear().col(1);
        if (kind == 2) world.matrix()(3, 0) = 0.125;
        if (kind == 3) world.matrix()(3, 3) = 2.0;
        if (kind == 4) world.translation().x() = nan;
        if (kind == 5) world.linear()(0, 0) = inf;
        if (kind == 6) world.translation().z() = 2.0 * double(std::numeric_limits<float>::max());
        const indexed_triangle_set before = fixture.volume->mesh().its;
        invalid(FacetPicker::from_source_facet(*fixture.volume, 0, world));
        unchanged(fixture.volume->mesh().its, before);
    }
    Transform3d finite_scale = identity;
    finite_scale.linear() = Vec3d(1e-5, 2e3, 4.0).asDiagonal();
    CHECK(FacetPicker::from_source_facet(*fixture.volume, 0, finite_scale).valid());
    for (int kind = 0; kind < 6; ++kind) {
        TriangleMesh mesh;
        mesh.its = original_facets(); // Malformed fixtures must not enter TriangleMesh's stats constructor.
        if (kind == 0) mesh.its.indices[0][1] = -1;
        if (kind == 1) mesh.its.indices[0][2] = int(mesh.its.vertices.size());
        if (kind == 2) mesh.its.vertices[1] = mesh.its.vertices[0];
        if (kind == 3) mesh.its.vertices[2] = 2.0f * mesh.its.vertices[1] - mesh.its.vertices[0];
        if (kind == 4) mesh.its.vertices[0].x() = float(nan);
        if (kind == 5) mesh.its.vertices[2].z() = float(inf);
        fixture.volume->set_mesh(std::move(mesh));
        const indexed_triangle_set before = fixture.volume->mesh().its;
        invalid(FacetPicker::from_source_facet(*fixture.volume, 0, identity));
        unchanged(fixture.volume->mesh().its, before);
    }
    fixture.volume->set_mesh(TriangleMesh(original_facets()));
    fixture.volume->set_type(ModelVolumeType::PARAMETER_MODIFIER);
    invalid(FacetPicker::from_source_facet(*fixture.volume, 0, identity));
    fixture.volume->set_type(ModelVolumeType::MODEL_PART);
    std::shared_ptr<const TriangleMesh> empty;
    fixture.volume->set_mesh(empty);
    invalid(FacetPicker::from_source_facet(*fixture.volume, 0, identity));
}

TEST_CASE("Facet highlight cache includes immutable mesh volume facet and exact pose", "[facet_picker_geometry]")
{
    Fixture fixture;
    const Transform3d identity = Transform3d::Identity();
    const auto first = FacetPicker::from_source_facet(*fixture.volume, 0, identity);
    REQUIRE(first.valid());
    CHECK(first.same_geometry(FacetPicker::from_source_facet(*fixture.volume, 0, identity)));
    CHECK_FALSE(first.same_geometry(FacetPicker::Hit{}));
    CHECK_FALSE(first.same_geometry(FacetPicker::from_source_facet(*fixture.volume, 1, identity)));
    Model copy(fixture.model); // Real Model clone deliberately preserves IDs and shares immutable meshes.
    const auto *other = copy.objects.front()->volumes.front();
    REQUIRE(other->id() == fixture.volume->id());
    REQUIRE(other->mesh_ptr() == fixture.volume->mesh_ptr());
    CHECK_FALSE(first.same_geometry(FacetPicker::from_source_facet(*other, 0, identity)));
    for (int kind = 0; kind < 3; ++kind) {
        Transform3d pose = identity;
        if (kind == 0) pose.translation() = Vec3d(450, 1024, 0);
        if (kind == 1) pose.linear() = Eigen::AngleAxisd(0.25, Vec3d::UnitY()).toRotationMatrix();
        if (kind == 2) pose.linear()(0, 0) = std::nextafter(1.0, 2.0);
        CHECK_FALSE(first.same_geometry(FacetPicker::from_source_facet(*fixture.volume, 0, pose)));
    }
    const auto old_id = fixture.volume->id();
    const indexed_triangle_set old_mesh = first.mesh->its;
    fixture.volume->set_mesh(TriangleMesh(original_facets())); // Equal contents, new immutable identity.
    REQUIRE(fixture.volume->id() == old_id);
    CHECK_FALSE(first.same_geometry(FacetPicker::from_source_facet(*fixture.volume, 0, identity)));
    unchanged(first.mesh->its, old_mesh); // Cached snapshots retain their source lifetime.
}
