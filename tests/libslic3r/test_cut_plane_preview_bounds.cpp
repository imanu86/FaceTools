// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include <catch2/catch.hpp>

#include "slic3r/GUI/GLModel.hpp"
#include "slic3r/GUI/Gizmos/FaceToolsPreviewBounds.hpp"

#include <array>
#include <cmath>
#include <limits>

using namespace Slic3r;
using Slic3r::GUI::GLModel;
using Slic3r::GUI::translated_face_tools_preview_bounds;

namespace {

// Keep the production caller's defined, double-finite and float-finite guard.
// Planarity is legitimate here; a positive-volume predicate is not.
bool accepted_preview_bounds(const BoundingBoxf3 &bounds)
{
    return bounds.defined && bounds.min.allFinite() && bounds.max.allFinite() &&
        bounds.min.cast<float>().allFinite() && bounds.max.cast<float>().allFinite();
}

std::array<Vec3f, 4> initialize_plane(GLModel &model, const Vec3d &normal, double radius)
{
    const Vec3d u = normal.unitOrthogonal(), v = normal.cross(u);
    GLModel::Geometry geometry;
    geometry.format = {GLModel::PrimitiveType::Triangles, GLModel::Geometry::EVertexLayout::P3N3};
    for (const Vec2d &corner : {Vec2d(-1,-1), Vec2d(1,-1), Vec2d(1,1), Vec2d(-1,1)}) {
        const Vec3f point = (radius * (corner.x() * u + corner.y() * v)).cast<float>();
        REQUIRE(point.allFinite());
        geometry.add_vertex(point, Vec3f(normal.cast<float>()));
    }
    geometry.add_triangle(0,1,2);
    geometry.add_triangle(0,2,3);
    REQUIRE(geometry.vertices_count() == 4);
    std::array<Vec3f, 4> corners;
    for (size_t i = 0; i < corners.size(); ++i)
        corners[i] = geometry.extract_position_3(i);
    // This overload populates bounds on the CPU. No wx app, GL context or render.
    model.init_from(std::move(geometry));
    REQUIRE(model.is_initialized());
    REQUIRE(model.get_indices_count() == 6);
    REQUIRE(model.mesh == nullptr);
    REQUIRE(model.get_bounding_box().defined);
    return corners;
}

} // namespace

TEST_CASE("Cut preview bounds preserve real GLModel planes before GPU upload", "[cut_plane_preview_bounds]")
{
    const std::array<Vec3d, 9> normals{{
        Vec3d::UnitZ(), Vec3d::UnitY(), Vec3d::UnitX(),
        -Vec3d::UnitZ(), -Vec3d::UnitY(), -Vec3d::UnitX(),
        Vec3d(1, 2, 3).normalized(), Vec3d(1, 1, 0).normalized(),
        Vec3d(1e-8, 0, 1).normalized()
    }};
    for (const Vec3d &normal : normals) {
        for (double radius : {1.0, 37.25}) {
            INFO("normal=" << normal.transpose() << " radius=" << radius);
            GLModel model;
            const auto corners = initialize_plane(model, normal, radius);
            const BoundingBoxf3 local = model.get_bounding_box();
            for (double position : {0.0, 450.0, 1024.0, 15000.0, 450.0, 0.0}) {
                INFO("position=" << position);
                const Vec3d translation(position, -0.5 * position, 0.25 * position);
                const BoundingBoxf3 bounds = translated_face_tools_preview_bounds(model.get_bounding_box(), translation);
                REQUIRE(accepted_preview_bounds(bounds));
                BoundingBoxf3 expected;
                for (const Vec3f &corner : corners) {
                    const Vec3d world_corner = corner.cast<double>() + translation;
                    expected.merge(world_corner);
                    CHECK(bounds.contains(world_corner));
                    CHECK(std::abs(normal.dot(corner.cast<double>())) <=
                          8.0 * std::numeric_limits<float>::epsilon() * radius);
                }
                // Exact bounds of the actual float vertices, not padded/idealized geometry.
                CHECK(bounds.min == expected.min);
                CHECK(bounds.max == expected.max);
                for (int axis = 0; axis < 3; ++axis)
                    if (local.min[axis] == local.max[axis])
                        CHECK(bounds.min[axis] == bounds.max[axis]);
                // Absolute positioning reuses the model; no accumulated translation.
                CHECK(model.get_bounding_box().defined == local.defined);
                CHECK(model.get_bounding_box().min == local.min);
                CHECK(model.get_bounding_box().max == local.max);
            }
        }
    }
}

TEST_CASE("Cut preview translation never defines empty or invalid bounds", "[cut_plane_preview_bounds]")
{
    GLModel model;
    REQUIRE_FALSE(model.is_initialized());
    for (const Vec3d &translation : {Vec3d(0, 0, 0), Vec3d(450, 1024, 15000)}) {
        const BoundingBoxf3 bounds = translated_face_tools_preview_bounds(model.get_bounding_box(), translation);
        CHECK_FALSE(bounds.defined);
        CHECK_FALSE(accepted_preview_bounds(bounds));
    }
    initialize_plane(model, Vec3d::UnitZ(), 10.0);
    model.reset(); // No GPU buffers were created.
    CHECK_FALSE(accepted_preview_bounds(translated_face_tools_preview_bounds(model.get_bounding_box(), Vec3d::Zero())));

    BoundingBoxf3 undefined;
    undefined.min = Vec3d(-1, -2, -3);
    undefined.max = Vec3d(1, 2, 3);
    const BoundingBoxf3 bounds = translated_face_tools_preview_bounds(undefined, Vec3d(10, 20, 30));
    CHECK_FALSE(bounds.defined); // Finite min/max must not resurrect an undefined input.
    CHECK_FALSE(accepted_preview_bounds(bounds));
    CHECK_FALSE(undefined.defined);
    CHECK(undefined.min == Vec3d(-1, -2, -3));
    CHECK(undefined.max == Vec3d(1, 2, 3));
}

TEST_CASE("Cut preview bounds retain NaN infinity and float overflow rejection", "[cut_plane_preview_bounds]")
{
    GLModel model;
    initialize_plane(model, Vec3d::UnitZ(), 10.0);
    const BoundingBoxf3 local = model.get_bounding_box();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const double float_overflow = 2.0 * double(std::numeric_limits<float>::max());
    for (double bad : {nan, inf, -inf, float_overflow, -float_overflow}) {
        for (int axis = 0; axis < 3; ++axis) {
            INFO("axis=" << axis << " value=" << bad);
            Vec3d translation = Vec3d::Zero();
            translation[axis] = bad;
            CHECK_FALSE(accepted_preview_bounds(translated_face_tools_preview_bounds(local, translation)));
            for (bool bad_minimum : {false, true}) {
                BoundingBoxf3 bad_bounds = local;
                (bad_minimum ? bad_bounds.min : bad_bounds.max)[axis] = bad;
                CHECK_FALSE(accepted_preview_bounds(translated_face_tools_preview_bounds(bad_bounds, Vec3d::Zero())));
            }
        }
    }
    CHECK(model.get_bounding_box().min == local.min);
    CHECK(model.get_bounding_box().max == local.max);
}

TEST_CASE("Cut preview bounds reject overflow produced by otherwise finite addition", "[cut_plane_preview_bounds]")
{
    const double float_max = double(std::numeric_limits<float>::max());
    GLModel model;
    initialize_plane(model, Vec3d::UnitZ(), 0.75 * float_max);
    const BoundingBoxf3 local = model.get_bounding_box();
    const Vec3d translation(0.75 * float_max, 0.0, 0.0);
    REQUIRE(accepted_preview_bounds(local));
    REQUIRE(translation.cast<float>().allFinite());
    const BoundingBoxf3 bounds = translated_face_tools_preview_bounds(local, translation);
    REQUIRE(bounds.min.allFinite());
    REQUIRE(bounds.max.allFinite());
    CHECK_FALSE(bounds.max.cast<float>().allFinite());
    CHECK_FALSE(accepted_preview_bounds(bounds));

    BoundingBoxf3 large;
    const double double_max = std::numeric_limits<double>::max();
    large.merge(Vec3d(0.5 * double_max, 0, 0));
    large.merge(Vec3d(double_max, 1, 0));
    REQUIRE(large.min.allFinite());
    REQUIRE(large.max.allFinite());
    const BoundingBoxf3 overflowed = translated_face_tools_preview_bounds(large, Vec3d(double_max, 0, 0));
    CHECK_FALSE(overflowed.max.allFinite());
    CHECK_FALSE(accepted_preview_bounds(overflowed));
}
