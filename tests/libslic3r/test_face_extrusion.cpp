// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include <catch2/catch.hpp>

#include "libslic3r/FaceExtrusion.hpp"
#include "libslic3r/MeshDiagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>

using namespace Slic3r;

namespace {

indexed_triangle_set square()
{
    indexed_triangle_set mesh;
    mesh.vertices = {{0.f, 0.f, 0.f}, {10.f, 0.f, 0.f}, {10.f, 10.f, 0.f}, {0.f, 10.f, 0.f}};
    mesh.indices = {{0, 1, 2}, {0, 2, 3}};
    return mesh;
}

indexed_triangle_set annulus(float inner_half_size = 5.f)
{
    indexed_triangle_set mesh;
    const float low = 10.f - inner_half_size, high = 10.f + inner_half_size;
    mesh.vertices = {{0.f, 0.f, 0.f}, {20.f, 0.f, 0.f}, {20.f, 20.f, 0.f}, {0.f, 20.f, 0.f},
                     {low, low, 0.f}, {high, low, 0.f}, {high, high, 0.f}, {low, high, 0.f}};
    mesh.indices = {{0, 1, 5}, {0, 5, 4}, {1, 2, 6}, {1, 6, 5},
                    {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7}};
    return mesh;
}

std::vector<int> all_faces(const indexed_triangle_set &mesh)
{
    std::vector<int> ids(mesh.indices.size());
    std::iota(ids.begin(), ids.end(), 0);
    return ids;
}

double signed_volume(const indexed_triangle_set &mesh)
{
    const Vec3d origin = mesh.vertices.front().cast<double>();
    double volume = 0.0;
    for (const Vec3i &f : mesh.indices) {
        const Vec3d a = mesh.vertices[f[0]].cast<double>() - origin;
        const Vec3d b = mesh.vertices[f[1]].cast<double>() - origin;
        const Vec3d c = mesh.vertices[f[2]].cast<double>() - origin;
        volume += a.dot(b.cross(c)) / 6.0;
    }
    return volume;
}

void require_closed(const indexed_triangle_set &mesh, int euler_characteristic)
{
    REQUIRE_FALSE(mesh.empty());
    const auto stats = its_mesh_diagnostics(mesh);
    CHECK(stats.open_edges == 0);
    CHECK(stats.non_manifold_edges == 0);
    CHECK(stats.non_manifold_vertices == 0);
    CHECK(stats.same_direction_edges == 0);
    CHECK_FALSE(stats.has_reversed_faces);
    std::set<std::pair<int, int>> edges;
    for (const Vec3i &f : mesh.indices)
        for (int side = 0; side < 3; ++side)
            edges.insert(std::minmax(f[side], f[(side + 1) % 3]));
    CHECK(int(mesh.vertices.size()) - int(edges.size()) + int(mesh.indices.size()) == euler_characteristic);
    CHECK(signed_volume(mesh) > 0.0);
    CHECK(mesh.properties.size() == mesh.indices.size());
    for (const FaceProperty &property : mesh.properties) {
        CHECK(property.type == eNormal);
        CHECK(property.area == 0.0);
    }
}

void require_error(const FaceExtrusionResult &result, FaceExtrusionError error)
{
    CHECK_FALSE(bool(result));
    CHECK(result.error == error);
    CHECK(result.mesh.empty());
    CHECK_FALSE(result.message.empty());
}

void require_unchanged(const indexed_triangle_set &a, const indexed_triangle_set &b)
{
    REQUIRE(a.vertices.size() == b.vertices.size());
    REQUIRE(a.indices.size() == b.indices.size());
    REQUIRE(a.properties.size() == b.properties.size());
    for (size_t i = 0; i < a.vertices.size(); ++i) CHECK(a.vertices[i] == b.vertices[i]);
    for (size_t i = 0; i < a.indices.size(); ++i) CHECK(a.indices[i] == b.indices[i]);
    for (size_t i = 0; i < a.properties.size(); ++i) {
        CHECK(a.properties[i].type == b.properties[i].type);
        CHECK(a.properties[i].area == b.properties[i].area);
    }
}

} // namespace

TEST_CASE("Face extrusion closes only boundary edges", "[face_extrusion]")
{
    const auto input = square();
    for (double depth : {3.0, -3.0}) {
        const auto result = extrude_planar_face(input, {0, 1}, Vec3d(0.0, 0.0, depth));
        INFO(result.message);
        REQUIRE(bool(result));
        CHECK(result.boundary_edges == 4);
        CHECK(result.mesh.vertices.size() == 8);
        CHECK(result.mesh.indices.size() == 12); // no walls along the cap diagonal
        CHECK(signed_volume(result.mesh) == Approx(300.0));
        require_closed(result.mesh, 2);
    }
}

TEST_CASE("Face extrusion supports one triangle and reversed source orientation", "[face_extrusion]")
{
    auto input = square();
    std::swap(input.indices[0][1], input.indices[0][2]);
    const auto result = extrude_planar_face(input, {0}, Vec3d(0.0, 0.0, 2.0));
    REQUIRE(bool(result));
    CHECK(result.mesh.vertices.size() == 6);
    CHECK(result.mesh.indices.size() == 8);
    CHECK(result.boundary_edges == 3);
    CHECK(signed_volume(result.mesh) == Approx(100.0));
    require_closed(result.mesh, 2);
}

TEST_CASE("Face extrusion preserves an annular through hole", "[face_extrusion]")
{
    for (float inner_half_size : {5.f, 0.001f}) {
        const auto input = annulus(inner_half_size);
        const auto result = extrude_planar_face(input, all_faces(input), Vec3d(0.0, 0.0, 2.0));
        INFO(result.message);
        REQUIRE(bool(result));
        CHECK(result.boundary_edges == 8);
        CHECK(result.mesh.vertices.size() == 16);
        CHECK(result.mesh.indices.size() == 32);
        const double inner_side = double(input.vertices[5].x()) - input.vertices[4].x();
        CHECK(signed_volume(result.mesh) == Approx((400.0 - inner_side * inner_side) * 2.0));
        require_closed(result.mesh, 0); // genus 1, not a filled box
        // Each original triangle appears in both caps; no cap covers the hole.
        for (size_t i = 0; i < 2 * input.indices.size(); ++i) {
            const Vec3i &f = result.mesh.indices[i];
            const Vec3d center = (result.mesh.vertices[f[0]].cast<double>() +
                                 result.mesh.vertices[f[1]].cast<double>() +
                                 result.mesh.vertices[f[2]].cast<double>()) / 3.0;
            CHECK_FALSE((std::abs(center.x() - 10.0) < inner_half_size &&
                         std::abs(center.y() - 10.0) < inner_half_size));
        }
    }
}

TEST_CASE("Face extrusion preserves multiple separate holes", "[face_extrusion]")
{
    indexed_triangle_set input;
    for (int y = 0; y <= 3; ++y)
        for (int x = 0; x <= 6; ++x)
            input.vertices.emplace_back(float(x), float(y), 0.f);
    for (int y = 0; y < 3; ++y)
        for (int x = 0; x < 6; ++x) {
            if (y == 1 && (x == 1 || x == 4)) continue;
            const int a = 7 * y + x, b = a + 1, c = a + 8, d = a + 7;
            input.indices.emplace_back(a, b, c);
            input.indices.emplace_back(a, c, d);
        }
    const auto result = extrude_planar_face(input, all_faces(input), Vec3d::UnitZ());
    INFO(result.message);
    REQUIRE(bool(result));
    CHECK(result.boundary_edges == 26);
    CHECK(result.mesh.indices.size() == 116);
    CHECK(signed_volume(result.mesh) == Approx(16.0));
    require_closed(result.mesh, -2); // genus 2
}

TEST_CASE("Face extrusion preserves a concave outline", "[face_extrusion]")
{
    indexed_triangle_set input;
    input.vertices = {{0.f, 0.f, 0.f}, {4.f, 0.f, 0.f}, {4.f, 1.f, 0.f},
                      {1.f, 1.f, 0.f}, {1.f, 4.f, 0.f}, {0.f, 4.f, 0.f}};
    input.indices = {{0, 1, 3}, {1, 2, 3}, {0, 3, 5}, {3, 4, 5}};
    const auto result = extrude_planar_face(input, all_faces(input), Vec3d(0.0, 0.0, 2.0));
    INFO(result.message);
    REQUIRE(bool(result));
    CHECK(result.boundary_edges == 6);
    CHECK(result.mesh.indices.size() == 20);
    CHECK(signed_volume(result.mesh) == Approx(14.0));
    require_closed(result.mesh, 2);
}

TEST_CASE("Face extrusion welds exact triangle-soup vertices without changing input", "[face_extrusion]")
{
    const auto indexed = square();
    indexed_triangle_set soup;
    for (const Vec3i &face : indexed.indices) {
        const int base = int(soup.vertices.size());
        for (int index : {face[0], face[1], face[2]}) soup.vertices.push_back(indexed.vertices[index]);
        soup.indices.emplace_back(base, base + 1, base + 2);
    }
    soup.properties.assign(2, FaceProperty{eSmallHole, 42.0});
    const auto original = soup;
    const auto result = extrude_planar_face(soup, {1, 0}, Vec3d(0.0, 0.0, 1.0));
    REQUIRE(bool(result));
    CHECK(result.mesh.vertices.size() == 8);
    CHECK(result.mesh.indices.size() == 12);
    require_unchanged(soup, original);
    require_closed(result.mesh, 2);
    const auto failed = extrude_planar_face(soup, {0, 1}, Vec3d::Zero());
    CHECK_FALSE(bool(failed));
    require_unchanged(soup, original);
}

TEST_CASE("Face extrusion world thickness survives nonuniform scale and reflection", "[face_extrusion]")
{
    auto local = annulus();
    Transform3d tilt = Transform3d::Identity();
    tilt.rotate(Eigen::AngleAxisd(0.6, Vec3d(1.0, 1.0, 0.0).normalized()));
    its_transform(local, tilt, true);
    const Vec3i seed = local.indices.front();
    const Vec3d local_normal = (local.vertices[seed[1]].cast<double>() - local.vertices[seed[0]].cast<double>()).cross(
                              local.vertices[seed[2]].cast<double>() - local.vertices[seed[0]].cast<double>()).normalized();
    for (double mirror : {1.0, -1.0}) {
        Transform3d transform = Transform3d::Identity();
        transform.translate(Vec3d(100.0, -50.0, 30.0));
        transform.rotate(Eigen::AngleAxisd(0.4, Vec3d::UnitZ()));
        transform.scale(Vec3d(mirror * 2.0, 3.0, 0.5));
        auto world = local;
        its_transform(world, transform, true);
        const Vec3d normal = (transform.linear().inverse().transpose() * local_normal).normalized();
        for (double sign : {1.0, -1.0}) {
            const Vec3d displacement = sign * 2.5 * normal;
            const auto result = extrude_planar_face(world, all_faces(world), displacement);
            INFO(result.message);
            REQUIRE(bool(result));
            require_closed(result.mesh, 0);
            const size_t count = result.mesh.vertices.size() / 2;
            for (size_t i = 0; i < count; ++i) {
                const Vec3d actual = result.mesh.vertices[i + count].cast<double>() - result.mesh.vertices[i].cast<double>();
                CHECK((actual - displacement).norm() < 1e-4);
                CHECK(actual.norm() == Approx(2.5).margin(1e-4));
            }
            auto stored = result.mesh;
            const Transform3d inverse = transform.inverse();
            its_transform(stored, inverse, true);
            CHECK(signed_volume(stored) > 0.0);
        }
    }
}

TEST_CASE("Face extrusion rejects malformed selections and geometry", "[face_extrusion]")
{
    auto input = square();
    const Vec3d displacement(0.0, 0.0, 1.0);
    SECTION("empty") { require_error(extrude_planar_face(input, {}, displacement), FaceExtrusionError::EmptySelection); }
    SECTION("out of range") { require_error(extrude_planar_face(input, {2}, displacement), FaceExtrusionError::InvalidTriangleIndex); }
    SECTION("negative id") { require_error(extrude_planar_face(input, {-1}, displacement), FaceExtrusionError::InvalidTriangleIndex); }
    SECTION("duplicate selection") { require_error(extrude_planar_face(input, {0, 0}, displacement), FaceExtrusionError::DuplicateTriangle); }
    SECTION("bad vertex") {
        input.indices[0][0] = -1;
        require_error(extrude_planar_face(input, {0}, displacement), FaceExtrusionError::InvalidVertex);
    }
    SECTION("repeated vertex") {
        input.indices[0][1] = input.indices[0][0];
        require_error(extrude_planar_face(input, {0}, displacement), FaceExtrusionError::DegenerateTriangle);
    }
    SECTION("collinear") {
        input.vertices[2] = Vec3f(20.f, 0.f, 0.f);
        require_error(extrude_planar_face(input, {0}, displacement), FaceExtrusionError::DegenerateTriangle);
    }
    SECTION("duplicate face reversed") {
        input.indices[1] = Vec3i(2, 1, 0);
        require_error(extrude_planar_face(input, {0, 1}, displacement), FaceExtrusionError::DuplicateTriangle);
    }
    SECTION("one reversed face") {
        std::swap(input.indices[1][0], input.indices[1][1]);
        require_error(extrude_planar_face(input, {0, 1}, displacement), FaceExtrusionError::InconsistentWinding);
    }
    SECTION("outside planarity tolerance") {
        input.vertices[3].z() = 0.01f;
        require_error(extrude_planar_face(input, {0, 1}, displacement), FaceExtrusionError::NonPlanar);
    }
    SECTION("within planarity tolerance, no silent flattening") {
        input.vertices[3].z() = 1e-5f;
        const auto result = extrude_planar_face(input, {0, 1}, displacement);
        REQUIRE(bool(result));
        CHECK(result.max_planarity_error_mm <= FaceExtrusionOptions{}.planarity_tolerance_mm);
        REQUIRE(result.mesh.vertices.size() == 8);
        for (size_t i = 0; i < 4; ++i) CHECK(result.mesh.vertices[i] == input.vertices[i]);
    }
    SECTION("NaN") {
        input.vertices[0].x() = std::numeric_limits<float>::quiet_NaN();
        require_error(extrude_planar_face(input, {0}, displacement), FaceExtrusionError::InvalidVertex);
    }
    SECTION("infinity") {
        input.vertices[0].y() = std::numeric_limits<float>::infinity();
        require_error(extrude_planar_face(input, {0}, displacement), FaceExtrusionError::InvalidVertex);
    }
    SECTION("unselected invalid face is irrelevant") {
        input.indices[1] = Vec3i(-1, -1, -1);
        CHECK(bool(extrude_planar_face(input, {0}, displacement)));
    }
}

TEST_CASE("Face extrusion rejects disconnected and non-manifold patches", "[face_extrusion]")
{
    const Vec3d displacement(0.0, 0.0, 1.0);
    SECTION("disjoint regions") {
        auto input = square();
        input.vertices.insert(input.vertices.end(), {Vec3f(20.f, 0.f, 0.f), Vec3f(30.f, 0.f, 0.f), Vec3f(20.f, 10.f, 0.f)});
        input.indices.emplace_back(4, 5, 6);
        require_error(extrude_planar_face(input, all_faces(input), displacement), FaceExtrusionError::Disconnected);
    }
    SECTION("three triangles on an edge") {
        indexed_triangle_set input;
        input.vertices = {{0.f, 0.f, 0.f}, {4.f, 0.f, 0.f}, {0.f, 4.f, 0.f}, {0.f, -4.f, 0.f}, {1.f, 3.f, 0.f}};
        input.indices = {{0, 1, 2}, {1, 0, 3}, {0, 1, 4}};
        require_error(extrude_planar_face(input, all_faces(input), displacement), FaceExtrusionError::NonManifold);
    }
    SECTION("bow-tie vertex") {
        indexed_triangle_set input;
        input.vertices = {{0.f, 0.f, 0.f}, {2.f, 0.f, 0.f}, {0.f, 2.f, 0.f}, {-2.f, 0.f, 0.f}, {0.f, -2.f, 0.f}};
        input.indices = {{0, 1, 2}, {0, 3, 4}};
        require_error(extrude_planar_face(input, all_faces(input), displacement), FaceExtrusionError::NonManifold);
    }
    SECTION("overlapping connected fan") {
        indexed_triangle_set input;
        input.vertices = {{0.f, 0.f, 0.f}, {2.f, 0.f, 0.f}, {0.f, 2.f, 0.f}, {-2.f, 0.f, 0.f},
                          {0.f, -2.f, 0.f}, {3.f, 0.f, 0.f}, {0.f, 3.f, 0.f}};
        input.indices = {{0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 5}, {0, 5, 6}};
        require_error(extrude_planar_face(input, all_faces(input), displacement), FaceExtrusionError::SelfIntersection);
    }
    SECTION("T junction is not automatically repaired") {
        auto input = square();
        input.vertices.emplace_back(5.f, 5.f, 0.f);
        input.indices = {{0, 1, 2}, {0, 4, 3}, {4, 2, 3}};
        CHECK_FALSE(bool(extrude_planar_face(input, all_faces(input), displacement)));
    }
    SECTION("nearby unwelded seam is not snapped") {
        auto input = square();
        input.vertices.emplace_back(0.f, 0.0001f, 0.f);
        input.vertices.emplace_back(9.9999f, 10.f, 0.f);
        input.indices[1] = Vec3i(4, 5, 3);
        require_error(extrude_planar_face(input, all_faces(input), displacement), FaceExtrusionError::Disconnected);
    }
}

TEST_CASE("Face extrusion rejects invalid displacement and unrepresentable output", "[face_extrusion]")
{
    auto input = square();
    for (const Vec3d &d : {Vec3d(0.0, 0.0, 0.0), Vec3d(0.0, 0.0, 1e-9), Vec3d(1.0, 0.0, 0.0),
                          Vec3d(0.0, 0.0, std::numeric_limits<double>::infinity()),
                          Vec3d(0.0, 0.0, std::numeric_limits<double>::quiet_NaN())})
        require_error(extrude_planar_face(input, {0, 1}, d), FaceExtrusionError::InvalidDisplacement);
    FaceExtrusionOptions bad;
    bad.planarity_tolerance_mm = -1.0;
    require_error(extrude_planar_face(input, {0, 1}, Vec3d::UnitZ(), bad), FaceExtrusionError::InvalidTolerance);
    bad = {};
    bad.angular_tolerance_rad = std::numeric_limits<double>::quiet_NaN();
    require_error(extrude_planar_face(input, {0, 1}, Vec3d::UnitZ(), bad), FaceExtrusionError::InvalidTolerance);
    SECTION("overflow") {
        require_error(extrude_planar_face(input, {0, 1}, Vec3d(0.0, 0.0, 1e300)), FaceExtrusionError::PrecisionLoss);
    }
    SECTION("float loses thickness far from origin") {
        for (Vec3f &p : input.vertices) p.z() = 1e8f;
        require_error(extrude_planar_face(input, {0, 1}, Vec3d::UnitZ()), FaceExtrusionError::PrecisionLoss);
    }
    SECTION("downstream float cross product overflow") {
        for (Vec3f &p : input.vertices) p *= 1e20f;
        require_error(extrude_planar_face(input, {0, 1}, Vec3d::UnitZ()), FaceExtrusionError::PrecisionLoss);
    }
}
