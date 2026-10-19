// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include <catch2/catch.hpp>

#include "libslic3r/PlanarProfile.hpp"
#include "libslic3r/MeshDiagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>

using namespace Slic3r;

namespace {
using Contour = std::vector<Vec3d>;
using Holes = std::vector<Contour>;

Contour rectangle(double x0, double y0, double x1, double y1, double z = 0.0)
{
    return {{x0, y0, z}, {x1, y0, z}, {x1, y1, z}, {x0, y1, z}};
}

void unchanged(const Contour &actual, const Contour &original)
{
    REQUIRE(actual.size() == original.size());
    for (size_t i = 0; i < actual.size(); ++i) CHECK(actual[i] == original[i]);
}

void same_mesh(const indexed_triangle_set &actual, const indexed_triangle_set &expected)
{
    REQUIRE(actual.vertices.size() == expected.vertices.size());
    REQUIRE(actual.indices.size() == expected.indices.size());
    REQUIRE(actual.properties.size() == expected.properties.size());
    for (size_t i = 0; i < actual.vertices.size(); ++i) CHECK(actual.vertices[i] == expected.vertices[i]);
    for (size_t i = 0; i < actual.indices.size(); ++i) CHECK(actual.indices[i] == expected.indices[i]);
    for (size_t i = 0; i < actual.properties.size(); ++i) {
        CHECK(actual.properties[i].type == expected.properties[i].type);
        CHECK(actual.properties[i].area == expected.properties[i].area);
    }
}

void require_error(const PlanarProfileResult &result, FaceExtrusionError error)
{
    CHECK_FALSE(bool(result));
    CHECK(result.error == error);
    CHECK_FALSE(result.message.empty());
    CHECK(result.cap.vertices.empty());
    CHECK(result.cap.indices.empty());
    CHECK(result.cap.properties.empty());
}

void require_error(const FaceExtrusionResult &result, FaceExtrusionError error)
{
    CHECK_FALSE(bool(result));
    CHECK(result.error == error);
    CHECK_FALSE(result.message.empty());
    CHECK(result.mesh.vertices.empty());
    CHECK(result.mesh.indices.empty());
    CHECK(result.mesh.properties.empty());
}

double volume(const indexed_triangle_set &mesh)
{
    const Vec3d origin = mesh.vertices.front().cast<double>();
    double result = 0.0;
    for (const Vec3i &f : mesh.indices) {
        const Vec3d a = mesh.vertices[f[0]].cast<double>() - origin;
        const Vec3d b = mesh.vertices[f[1]].cast<double>() - origin;
        const Vec3d c = mesh.vertices[f[2]].cast<double>() - origin;
        result += a.dot(b.cross(c)) / 6.0;
    }
    return result;
}

void require_cap(const PlanarProfileResult &result, const Contour &outer, const Holes &holes, double area)
{
    INFO(result.message);
    REQUIRE(bool(result));
    Contour all = outer;
    for (const auto &hole : holes) all.insert(all.end(), hole.begin(), hole.end());
    CHECK(result.origin == outer.front());
    CHECK(result.normal.norm() == Approx(1.0));
    CHECK(result.area_mm2 == Approx(area).epsilon(1e-9).margin(1e-9));
    CHECK(result.boundary_edges == all.size());
    REQUIRE(result.cap.vertices.size() == all.size());
    REQUIRE(result.cap.indices.size() == all.size() + 2 * holes.size() - 2);
    for (size_t i = 0; i < all.size(); ++i)
        CHECK(result.cap.vertices[i] == all[i].cast<float>()); // No flattening/recentering/hull.
    std::set<int> referenced;
    for (const Vec3i &face : result.cap.indices) {
        for (int index : {face[0], face[1], face[2]}) {
            REQUIRE(index >= 0);
            REQUIRE(size_t(index) < all.size());
            referenced.insert(index);
        }
        const Vec3d n = (result.cap.vertices[face[1]].cast<double>() - result.cap.vertices[face[0]].cast<double>()).cross(
                        result.cap.vertices[face[2]].cast<double>() - result.cap.vertices[face[0]].cast<double>());
        CHECK(n.dot(result.normal) > 0.0);
    }
    CHECK(referenced.size() == all.size());
    const auto stats = its_mesh_diagnostics(result.cap);
    CHECK(stats.open_edges == all.size());
    CHECK(stats.non_manifold_edges == 0);
    CHECK(stats.non_manifold_vertices == 0);
    CHECK(stats.same_direction_edges == 0);
}

void require_solid(const FaceExtrusionResult &result, size_t vertices, size_t holes, double expected_volume)
{
    INFO(result.message);
    REQUIRE(bool(result));
    REQUIRE_FALSE(result.mesh.empty());
    REQUIRE(result.mesh.vertices.size() == 2 * vertices);
    CHECK(result.boundary_edges == vertices);
    const auto stats = its_mesh_diagnostics(result.mesh);
    CHECK(stats.open_edges == 0);
    CHECK(stats.non_manifold_edges == 0);
    CHECK(stats.non_manifold_vertices == 0);
    CHECK(stats.same_direction_edges == 0);
    CHECK_FALSE(stats.has_reversed_faces);
    std::set<std::pair<int, int>> edges;
    for (const Vec3i &face : result.mesh.indices)
        for (int side = 0; side < 3; ++side)
            edges.insert(std::minmax(face[side], face[(side + 1) % 3]));
    CHECK(int(result.mesh.vertices.size()) - int(edges.size()) + int(result.mesh.indices.size()) == 2 - 2 * int(holes));
    CHECK(volume(result.mesh) == Approx(expected_volume).epsilon(1e-4).margin(1e-6));
    CHECK(volume(result.mesh) > 0.0);
    REQUIRE(result.mesh.properties.size() == result.mesh.indices.size());
    for (const FaceProperty &property : result.mesh.properties) {
        CHECK(property.type == eNormal);
        CHECK(property.area == 0.0);
    }
}

bool covered_xy(const indexed_triangle_set &cap, const Vec2d &point)
{
    auto cross = [](const Vec2d &a, const Vec2d &b) { return a.x() * b.y() - a.y() * b.x(); };
    for (const Vec3i &face : cap.indices) {
        bool positive = true, negative = true;
        for (int side = 0; side < 3; ++side) {
            const Vec2d a = cap.vertices[face[side]].head<2>().cast<double>();
            const Vec2d b = cap.vertices[face[(side + 1) % 3]].head<2>().cast<double>();
            const double signed_area = cross(b - a, point - a);
            positive &= signed_area >= -1e-10;
            negative &= signed_area <= 1e-10;
        }
        if (positive || negative) return true;
    }
    return false;
}
} // namespace

TEST_CASE("Planar profile quad preserves positions and delegates signed extrusion", "[planar_profile]")
{
    const Contour triangle{{0, 0, 0}, {6, 0, 0}, {0, 4, 0}};
    require_cap(triangulate_planar_profile(triangle), triangle, {}, 12.0);
    require_solid(extrude_planar_profile(triangle, Vec3d(0, 0, 2)), 3, 0, 24.0);
    auto outer = rectangle(1, 2, 11, 10, 7);
    const auto original = outer;
    for (bool reversed : {false, true}) {
        if (reversed) std::reverse(outer.begin(), outer.end());
        const auto cap = triangulate_planar_profile(outer);
        require_cap(cap, outer, {}, 80.0);
        CHECK(cap.normal.dot(Vec3d::UnitZ()) == Approx(reversed ? -1.0 : 1.0));
        for (double sign : {-1.0, 1.0}) {
            const Vec3d displacement = sign * 2.5 * cap.normal;
            const auto result = extrude_planar_profile(outer, displacement);
            require_solid(result, 4, 0, 200.0);
            std::vector<int> faces(cap.cap.indices.size());
            std::iota(faces.begin(), faces.end(), 0);
            const auto directly = extrude_planar_face(cap.cap, faces, displacement);
            REQUIRE(bool(directly));
            same_mesh(result.mesh, directly.mesh);
            for (size_t i = 0; i < 4; ++i)
                CHECK((result.mesh.vertices[i + 4].cast<double>() - result.mesh.vertices[i].cast<double>() - displacement).norm() < 1e-5);
        }
    }
    std::reverse(outer.begin(), outer.end());
    unchanged(outer, original);
}

TEST_CASE("Planar profile triangulates concave outlines without a convex hull", "[planar_profile]")
{
    const Contour outer{{0, 0, 0}, {4, 0, 0}, {4, 1, 0}, {1, 1, 0}, {1, 4, 0}, {0, 4, 0}};
    const auto before = outer;
    const auto result = triangulate_planar_profile(outer);
    require_cap(result, outer, {}, 7.0);
    CHECK(covered_xy(result.cap, Vec2d(3.4, 0.4)));
    CHECK(covered_xy(result.cap, Vec2d(0.3, 3.2)));
    CHECK_FALSE(covered_xy(result.cap, Vec2d(2.3, 2.1)));
    for (double depth : {-10.0, 10.0})
        require_solid(extrude_planar_profile(outer, Vec3d(0, 0, depth)), 6, 0, 70.0);
    unchanged(outer, before);
}

TEST_CASE("Planar profile preserves holes of either winding and multiple holes", "[planar_profile]")
{
    const auto outer = rectangle(0, 0, 20, 20);
    Holes holes{rectangle(2, 2, 7, 8), rectangle(11, 12, 17, 18)};
    const auto original = holes;
    for (unsigned mask = 0; mask < 4; ++mask) {
        holes = original;
        for (size_t i = 0; i < holes.size(); ++i)
            if (mask & (1u << i)) std::reverse(holes[i].begin(), holes[i].end());
        const auto before = holes;
        const auto cap = triangulate_planar_profile(outer, holes);
        require_cap(cap, outer, holes, 334.0);
        CHECK_FALSE(covered_xy(cap.cap, Vec2d(4, 5)));
        CHECK_FALSE(covered_xy(cap.cap, Vec2d(14, 15)));
        CHECK(covered_xy(cap.cap, Vec2d(9, 9)));
        for (double depth : {-2.0, 2.0})
            require_solid(extrude_planar_profile(outer, Vec3d(0, 0, depth), holes), 12, 2, 668.0);
        for (size_t i = 0; i < holes.size(); ++i) unchanged(holes[i], before[i]);
    }
    const Contour concave_hole{{3, 3, 0}, {7, 3, 0}, {7, 4, 0}, {4, 4, 0}, {4, 7, 0}, {3, 7, 0}};
    const auto cap = triangulate_planar_profile(outer, {concave_hole});
    require_cap(cap, outer, {concave_hole}, 393.0);
    CHECK_FALSE(covered_xy(cap.cap, Vec2d(6.3, 3.4)));
    CHECK_FALSE(covered_xy(cap.cap, Vec2d(3.4, 6.2)));
    CHECK(covered_xy(cap.cap, Vec2d(5.3, 5.1)));
    require_solid(extrude_planar_profile(outer, Vec3d(0, 0, 3), {concave_hole}), 10, 1, 1179.0);
    const Holes small_hole{rectangle(9.999, 9.999, 10.001, 10.001)};
    const auto narrow = triangulate_planar_profile(outer, small_hole);
    require_cap(narrow, outer, small_hole, 400.0 - 0.002 * 0.002);
    CHECK_FALSE(covered_xy(narrow.cap, Vec2d(10, 10)));
    require_solid(extrude_planar_profile(outer, Vec3d(0, 0, 2), small_hole), 8, 1, 2.0 * (400.0 - 0.002 * 0.002));
}

TEST_CASE("Planar profile retains collinear boundary points instead of welding or dropping them", "[planar_profile]")
{
    Contour outer{{0, 0, 0}, {5, 0, 0}, {10, 0, 0}, {10, 5, 0}, {10, 10, 0}, {5, 10, 0}, {0, 10, 0}, {0, 5, 0}};
    Holes holes{{{2, 2, 0}, {5, 2, 0}, {8, 2, 0}, {8, 5, 0}, {8, 8, 0}, {5, 8, 0}, {2, 8, 0}, {2, 5, 0}}};
    for (bool transformed : {false, true}) {
        if (transformed) {
            Transform3d transform = Transform3d::Identity();
            transform.translate(Vec3d(30, 40, 50));
            transform.rotate(Eigen::AngleAxisd(0.4, Vec3d(1, 2, 3).normalized()));
            for (Vec3d &p : outer) p = transform * p;
            for (Vec3d &p : holes.front()) p = transform * p;
        }
        const auto cap = triangulate_planar_profile(outer, holes);
        require_cap(cap, outer, holes, 64.0);
        require_solid(extrude_planar_profile(outer, cap.normal * 2.0, holes), 16, 1, 128.0);
    }
}

TEST_CASE("Planar profile supports transformed world coordinates and reflected ordering", "[planar_profile]")
{
    const auto local = rectangle(0, 0, 10, 8);
    const auto hole = rectangle(2, 2, 5, 4);
    for (unsigned mask = 0; mask < 4; ++mask) {
        Transform3d volume_transform = Transform3d::Identity();
        volume_transform.translate(Vec3d(8, -5, 6));
        volume_transform.rotate(Eigen::AngleAxisd(0.4, Vec3d(1, 2, 3).normalized()));
        volume_transform.scale(Vec3d(mask & 1 ? -1.2 : 1.2, 0.8, 1.6));
        Transform3d instance = Transform3d::Identity();
        instance.translate(Vec3d(60, 70, 50));
        instance.rotate(Eigen::AngleAxisd(0.2, Vec3d::UnitY()));
        instance.scale(Vec3d(0.75, mask & 2 ? -1.3 : 1.3, 1.8));
        const Transform3d world = instance * volume_transform;
        auto outer = local;
        Holes holes{hole};
        for (Vec3d &p : outer) p = world * p;
        for (Vec3d &p : holes.front()) p = world * p;
        const auto before = outer;
        const auto holes_before = holes;
        const Vec3d x = world.linear() * Vec3d::UnitX(), y = world.linear() * Vec3d::UnitY();
        const Vec3d expected_normal = x.cross(y).normalized();
        const double expected_area = 74.0 * x.cross(y).norm();
        const auto cap = triangulate_planar_profile(outer, holes);
        require_cap(cap, outer, holes, expected_area);
        CHECK(cap.normal.dot(expected_normal) == Approx(1.0));
        for (double sign : {-1.0, 1.0}) {
            const Vec3d displacement = sign * 3.0 * expected_normal;
            const auto result = extrude_planar_profile(outer, displacement, holes);
            require_solid(result, 8, 1, expected_area * 3.0);
            for (size_t i = 0; i < 8; ++i) {
                const Vec3d actual = result.mesh.vertices[i + 8].cast<double>() - result.mesh.vertices[i].cast<double>();
                CHECK((actual - displacement).norm() < 1e-4);
                CHECK(actual.norm() == Approx(3.0).margin(1e-4));
            }
        }
        unchanged(outer, before);
        unchanged(holes.front(), holes_before.front());
    }
}

TEST_CASE("Planar profile honors planarity tolerance without flattening accepted points", "[planar_profile]")
{
    auto outer = rectangle(0, 0, 10, 10);
    outer[2].z() = double(1e-5f);
    const auto original = outer;
    const auto accepted = triangulate_planar_profile(outer);
    require_cap(accepted, outer, {}, 100.0);
    CHECK(accepted.cap.vertices[2].z() == 1e-5f);
    REQUIRE(bool(extrude_planar_profile(outer, accepted.normal * 2.0)));
    FaceExtrusionOptions tight;
    tight.planarity_tolerance_mm = 1e-7;
    require_error(triangulate_planar_profile(outer, {}, tight), FaceExtrusionError::NonPlanar);
    unchanged(outer, original);
    outer[2].z() = 0.1;
    require_error(triangulate_planar_profile(outer), FaceExtrusionError::NonPlanar);
    const auto flat = rectangle(0, 0, 10, 10);
    require_error(triangulate_planar_profile(flat, {rectangle(2, 2, 4, 4, 1)}), FaceExtrusionError::NonPlanar);
}

TEST_CASE("Planar profile rejects empty duplicate degenerate and nonfinite point input", "[planar_profile]")
{
    require_error(triangulate_planar_profile({}), FaceExtrusionError::EmptySelection);
    require_error(triangulate_planar_profile({{0, 0, 0}, {1, 0, 0}}), FaceExtrusionError::EmptySelection);
    require_error(triangulate_planar_profile({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}}), FaceExtrusionError::DegenerateTriangle);
    const auto valid = rectangle(0, 0, 10, 10);
    require_error(triangulate_planar_profile(valid, {Contour{}}), FaceExtrusionError::EmptySelection);
    auto repeated_close = valid;
    repeated_close.push_back(repeated_close.front());
    require_error(triangulate_planar_profile(repeated_close), FaceExtrusionError::InvalidVertex);
    auto duplicate = valid;
    duplicate.insert(duplicate.begin() + 1, duplicate.front());
    require_error(triangulate_planar_profile(duplicate), FaceExtrusionError::InvalidVertex);
    duplicate = valid;
    duplicate.insert(duplicate.begin() + 1, Vec3d(5e-7, 0, 0));
    require_error(triangulate_planar_profile(duplicate), FaceExtrusionError::InvalidVertex);
    for (double invalid : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        auto bad = valid;
        bad[2].x() = invalid;
        require_error(triangulate_planar_profile(bad), FaceExtrusionError::InvalidVertex);
        require_error(extrude_planar_profile(bad, Vec3d::UnitZ()), FaceExtrusionError::InvalidVertex);
    }
    unchanged(valid, rectangle(0, 0, 10, 10));
}

TEST_CASE("Planar profile rejects crossing touching backtracking and invalid hole containment", "[planar_profile]")
{
    const auto valid = rectangle(0, 0, 10, 10);
    const Contour bowtie{{0, 0, 0}, {10, 10, 0}, {0, 10, 0}, {10, 0, 0}};
    require_error(triangulate_planar_profile(bowtie), FaceExtrusionError::SelfIntersection);
    const Contour backtrack{{0, 0, 0}, {5, 0, 0}, {2, 0, 0}, {5, 5, 0}, {0, 5, 0}};
    require_error(triangulate_planar_profile(backtrack), FaceExtrusionError::SelfIntersection);
    const Contour touch{{0, 0, 0}, {10, 0, 0}, {10, 10, 0}, {5, 0, 0}, {0, 10, 0}};
    require_error(triangulate_planar_profile(touch), FaceExtrusionError::SelfIntersection);
    for (const auto &holes : {Holes{rectangle(12, 2, 14, 4)},
                              Holes{rectangle(9, 2, 11, 4)},
                              Holes{rectangle(0, 2, 2, 4)},
                              Holes{rectangle(2, 2, 6, 6), rectangle(4, 4, 8, 8)},
                              Holes{rectangle(1, 1, 9, 9), rectangle(3, 3, 5, 5)},
                              Holes{rectangle(-1, -1, 11, 11)},
                              Holes{bowtie}}) {
        const auto before = holes;
        const auto result = triangulate_planar_profile(valid, holes);
        CHECK_FALSE(bool(result));
        CHECK(result.cap.empty());
        CHECK_FALSE(result.message.empty());
        for (size_t i = 0; i < holes.size(); ++i) unchanged(holes[i], before[i]);
    }
}

TEST_CASE("Planar profile rejects invalid displacement and invalid tolerances", "[planar_profile]")
{
    const auto outer = rectangle(0, 0, 10, 10);
    for (const Vec3d &d : {Vec3d(0, 0, 0), Vec3d(1, 0, 0), Vec3d(1, 0, 1), Vec3d(0, 0, 5e-7),
                           Vec3d(0, 0, std::numeric_limits<double>::infinity()),
                           Vec3d(0, 0, std::numeric_limits<double>::quiet_NaN())})
        require_error(extrude_planar_profile(outer, d), FaceExtrusionError::InvalidDisplacement);
    for (size_t field = 0; field < 3; ++field) {
        FaceExtrusionOptions options;
        double *value = field == 0 ? &options.planarity_tolerance_mm :
                        field == 1 ? &options.length_tolerance_mm : &options.angular_tolerance_rad;
        for (double invalid : {-1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
            *value = invalid;
            require_error(triangulate_planar_profile(outer, {}, options), FaceExtrusionError::InvalidTolerance);
            require_error(extrude_planar_profile(outer, Vec3d::UnitZ(), {}, options), FaceExtrusionError::InvalidTolerance);
        }
    }
    FaceExtrusionOptions options;
    options.angular_tolerance_rad = 0.5 * std::acos(-1.0);
    require_error(triangulate_planar_profile(outer, {}, options), FaceExtrusionError::InvalidTolerance);
    options.angular_tolerance_rad = 0.0;
    require_error(triangulate_planar_profile(outer, {}, options), FaceExtrusionError::InvalidTolerance);
    options = {};
    options.length_tolerance_mm = 0.0;
    require_error(triangulate_planar_profile(outer, {}, options), FaceExtrusionError::InvalidTolerance);
    options = {};
    options.planarity_tolerance_mm = 0.0;
    require_cap(triangulate_planar_profile(outer, {}, options), outer, {}, 100.0);
}

TEST_CASE("Planar profile rejects lossy coordinates collapsed holes and unrepresentable thickness", "[planar_profile]")
{
    const auto far = rectangle(1e8 + 1, 1e8 + 1, 1e8 + 11, 1e8 + 11);
    const auto before = far;
    require_error(triangulate_planar_profile(far), FaceExtrusionError::PrecisionLoss);
    unchanged(far, before);
    const auto representable = rectangle(1e7, 1e7, 1e7 + 16, 1e7 + 16, 1e7);
    require_cap(triangulate_planar_profile(representable), representable, {}, 256.0);
    require_error(extrude_planar_profile(representable, Vec3d(0, 0, 0.001)), FaceExtrusionError::PrecisionLoss);
    auto enormous = rectangle(0, 0, 10, 10);
    enormous[0].x() = 2.0 * std::numeric_limits<float>::max();
    require_error(triangulate_planar_profile(enormous), FaceExtrusionError::PrecisionLoss);
    const double width = double(1e20f);
    require_error(triangulate_planar_profile(rectangle(0, 0, width, width)), FaceExtrusionError::PrecisionLoss);
    FaceExtrusionOptions tiny;
    tiny.length_tolerance_mm = 1e-9;
    const Holes hole{rectangle(10, 10, 10 + 2e-7, 10 + 2e-7)};
    require_error(triangulate_planar_profile(rectangle(0, 0, 20, 20), hole, tiny), FaceExtrusionError::PrecisionLoss);
}

TEST_CASE("Planar profile results are deterministic and all calls leave contours immutable", "[planar_profile]")
{
    const auto outer = rectangle(-10, -10, 20, 20, 3);
    const Holes holes{rectangle(-5, -5, -1, 1, 3), rectangle(7, 8, 12, 13, 3)};
    const auto outer_before = outer;
    const auto holes_before = holes;
    const auto reference = triangulate_planar_profile(outer, holes);
    REQUIRE(bool(reference));
    const auto solid = extrude_planar_profile(outer, Vec3d(0, 0, -2), holes);
    REQUIRE(bool(solid));
    for (int i = 0; i < 8; ++i) {
        const auto repeated = triangulate_planar_profile(outer, holes);
        REQUIRE(bool(repeated));
        same_mesh(repeated.cap, reference.cap);
        CHECK(repeated.origin == reference.origin);
        CHECK(repeated.normal == reference.normal);
        CHECK(repeated.area_mm2 == reference.area_mm2);
        const auto extruded = extrude_planar_profile(outer, Vec3d(0, 0, -2), holes);
        REQUIRE(bool(extruded));
        same_mesh(extruded.mesh, solid.mesh);
        unchanged(outer, outer_before);
        for (size_t h = 0; h < holes.size(); ++h) unchanged(holes[h], holes_before[h]);
    }
}

TEST_CASE("Planar profile bounds synchronous work before copying or triangulating huge input", "[planar_profile]")
{
    const Contour huge(4097, Vec3d::Zero());
    require_error(triangulate_planar_profile(huge), FaceExtrusionError::TooLarge);
    require_error(triangulate_planar_profile(rectangle(0, 0, 10, 10), {huge}), FaceExtrusionError::TooLarge);
}
