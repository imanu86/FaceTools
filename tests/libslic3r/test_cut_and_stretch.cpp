// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include <catch2/catch.hpp>

#include "libslic3r/CutAndStretch.hpp"
#include "libslic3r/FaceExtrusion.hpp"
#include "libslic3r/MeshBoolean.hpp"
#include "libslic3r/MeshDiagnostics.hpp"
#include "libslic3r/TriangleMeshSlicer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

using namespace Slic3r;

namespace {

std::vector<int> faces(const indexed_triangle_set &mesh)
{
    std::vector<int> out(mesh.indices.size());
    std::iota(out.begin(), out.end(), 0);
    return out;
}

indexed_triangle_set l_arm()
{
    // Shaft x=0..10,z=0..20; arm x=0..30,z=20..24; depth y=10.
    indexed_triangle_set cap;
    cap.vertices = {{0, 0, 0}, {10, 0, 0}, {10, 0, 20}, {30, 0, 20}, {30, 0, 24}, {0, 0, 24}};
    cap.indices = {{0, 1, 2}, {0, 2, 5}, {2, 3, 4}, {2, 4, 5}};
    auto source = extrude_planar_face(cap, faces(cap), Vec3d(0, 10, 0));
    INFO(source.message);
    REQUIRE(bool(source));
    return source.mesh;
}

indexed_triangle_set hollow_shaft()
{
    indexed_triangle_set cap;
    cap.vertices = {{0, 0, 0}, {20, 0, 0}, {20, 20, 0}, {0, 20, 0},
                    {5, 5, 0}, {15, 5, 0}, {15, 15, 0}, {5, 15, 0}};
    cap.indices = {{0, 1, 5}, {0, 5, 4}, {1, 2, 6}, {1, 6, 5},
                   {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7}};
    auto source = extrude_planar_face(cap, faces(cap), Vec3d(0, 0, 30));
    INFO(source.message);
    REQUIRE(bool(source));
    return source.mesh;
}

double volume(const indexed_triangle_set &mesh)
{
    const Vec3d o = mesh.vertices.front().cast<double>();
    double out = 0;
    for (const Vec3i &f : mesh.indices) {
        const Vec3d a = mesh.vertices[f[0]].cast<double>() - o;
        const Vec3d b = mesh.vertices[f[1]].cast<double>() - o;
        const Vec3d c = mesh.vertices[f[2]].cast<double>() - o;
        out += a.dot(b.cross(c)) / 6.0;
    }
    return out;
}

void require_solid(const indexed_triangle_set &mesh, int euler = 2)
{
    REQUIRE_FALSE(mesh.empty());
    const auto d = its_mesh_diagnostics(mesh);
    CHECK(d.open_edges == 0);
    CHECK(d.non_manifold_edges == 0);
    CHECK(d.non_manifold_vertices == 0);
    CHECK(d.same_direction_edges == 0);
    CHECK_FALSE(d.has_reversed_faces);
    const auto solid = MeshBoolean::cgal::triangle_mesh_to_cgal(mesh);
    CHECK_FALSE(MeshBoolean::cgal::does_self_intersect(*solid));
    CHECK(MeshBoolean::cgal::does_bound_a_volume(*solid));
    CHECK(volume(mesh) > 0);
    std::set<int> vertices;
    std::set<std::pair<int, int>> edges;
    for (const Vec3i &f : mesh.indices)
        for (int side = 0; side < 3; ++side) {
            vertices.insert(f[side]);
            edges.insert(std::minmax(f[side], f[(side + 1) % 3]));
        }
    CHECK(int(vertices.size()) - int(edges.size()) + int(mesh.indices.size()) == euler);
    std::vector<Vec3i> adjacent;
    its_edge_diagnostics(mesh, &adjacent);
    std::set<int> seen{0};
    std::vector<int> queue{0};
    for (size_t i = 0; i < queue.size(); ++i)
        for (int next : {adjacent[queue[i]][0], adjacent[queue[i]][1], adjacent[queue[i]][2]}) {
            REQUIRE(next >= 0);
            if (seen.insert(next).second) queue.push_back(next);
        }
    CHECK(queue.size() == mesh.indices.size());
    CHECK(mesh.properties.size() == mesh.indices.size());
}

void unchanged(const indexed_triangle_set &before, const indexed_triangle_set &after)
{
    REQUIRE(before.vertices.size() == after.vertices.size());
    REQUIRE(before.indices.size() == after.indices.size());
    REQUIRE(before.properties.size() == after.properties.size());
    for (size_t i = 0; i < before.vertices.size(); ++i) CHECK(before.vertices[i] == after.vertices[i]);
    for (size_t i = 0; i < before.indices.size(); ++i) CHECK(before.indices[i] == after.indices[i]);
    for (size_t i = 0; i < before.properties.size(); ++i) {
        CHECK(before.properties[i].type == after.properties[i].type);
        CHECK(before.properties[i].area == after.properties[i].area);
    }
}

void failed_plan(const CutAndStretchPlan &plan)
{
    CHECK_FALSE(bool(plan));
    CHECK_FALSE(plan.message.empty());
    CHECK(plan.stationary.empty());
    CHECK(plan.moving.empty());
    CHECK(plan.bridge_unit.empty());
}

void failed_result(const CutAndStretchResult &result)
{
    CHECK_FALSE(bool(result));
    CHECK_FALSE(result.message.empty());
    CHECK(result.mesh.empty());
}

std::pair<double, double> arm_extents(const indexed_triangle_set &mesh, const Transform3d &local_to_world,
                                    const Vec3d &normal, const Vec3d &plane_origin)
{
    const Transform3d inverse = local_to_world.inverse();
    double low = std::numeric_limits<double>::infinity(), high = -low;
    for (const Vec3f &vertex : mesh.vertices) {
        const Vec3d p = vertex.cast<double>();
        if ((inverse * p).x() > 15.0) {
            const double projection = (p - plane_origin).dot(normal);
            low = std::min(low, projection);
            high = std::max(high, projection);
        }
    }
    REQUIRE(std::isfinite(low));
    REQUIRE(std::isfinite(high));
    return {low, high};
}

} // namespace

TEST_CASE("Cut and extend moves an L arm rigidly and extends only its shaft", "[cut_and_stretch]")
{
    const auto source = l_arm();
    const auto before = source;
    const auto plan = prepare_cut_and_stretch(source, Vec3d(0, 0, 20), Vec3d::UnitZ());
    INFO(plan.message);
    REQUIRE(bool(plan));
    CHECK(plan.contact_area_mm2 == Approx(100.0));
    CHECK(plan.source_volume_mm3 == Approx(3200.0));
    const auto stationary_before = plan.stationary;
    const auto moving_before = plan.moving;
    const auto bridge_before = plan.bridge_unit;
    for (double extension : {1.0, 10.0, 17.5, 100.125}) {
        const auto result = apply_cut_and_stretch(plan, extension);
        INFO(result.message);
        REQUIRE(bool(result));
        require_solid(result.mesh);
        CHECK(volume(result.mesh) == Approx(3200.0 + 100.0 * extension).epsilon(1e-5));
        const auto arm = arm_extents(result.mesh, Transform3d::Identity(), Vec3d::UnitZ(), Vec3d(0, 0, 20));
        CHECK(arm.first == Approx(extension).margin(1e-5));
        CHECK(arm.second == Approx(extension + 4.0).margin(1e-5));
        CHECK(arm.second - arm.first == Approx(4.0).margin(1e-5));
    }
    unchanged(before, source);
    unchanged(stationary_before, plan.stationary);
    unchanged(moving_before, plan.moving);
    unchanged(bridge_before, plan.bridge_unit);
}

TEST_CASE("Cut and extend preserves a hollow shaft and its contact hole", "[cut_and_stretch]")
{
    const auto source = hollow_shaft();
    const auto before = source;
    const auto plan = prepare_cut_and_stretch(source, Vec3d(0, 0, 20), Vec3d::UnitZ());
    INFO(plan.message);
    REQUIRE(bool(plan));
    CHECK(plan.contact_area_mm2 == Approx(300.0));
    const auto result = apply_cut_and_stretch(plan, 10.0);
    INFO(result.message);
    REQUIRE(bool(result));
    require_solid(result.mesh, 0);
    CHECK(volume(result.mesh) == Approx(12000.0));
    const auto sections = slice_mesh_ex(result.mesh, std::vector<float>{10.f, 25.f, 35.f});
    REQUIRE(sections.size() == 3);
    for (const auto &section : sections) {
        REQUIRE(section.size() == 1);
        CHECK(section.front().holes.size() == 1);
    }
    unchanged(before, source);
}

TEST_CASE("Cut and extend respects arbitrary planes, nonuniform scale and mirrors", "[cut_and_stretch]")
{
    for (const Vec3d &scale : {Vec3d(1, 1, 1), Vec3d(1.5, 0.75, 1.25), Vec3d(-1.5, 0.75, 1.25)}) {
        Transform3d transform = Transform3d::Identity();
        const Matrix3d rotation = Eigen::AngleAxisd(0.57, Vec3d(1, 2, 3).normalized()).toRotationMatrix();
        transform.linear() = rotation * scale.asDiagonal();
        transform.translation() = Vec3d(3, -4, 2);
        TriangleMesh source(l_arm());
        source.transform(transform, true);
        const auto before = source.its;
        const Vec3d origin = transform * Vec3d(0, 0, 20);
        const Vec3d normal = (transform.linear().inverse().transpose() * Vec3d::UnitZ()).normalized();
        const auto plan = prepare_cut_and_stretch(source.its, origin, normal);
        INFO(scale.transpose());
        INFO(plan.message);
        REQUIRE(bool(plan));
        const double area = 100.0 * std::abs(scale.x() * scale.y());
        CHECK(plan.contact_area_mm2 == Approx(area).epsilon(1e-5));
        const auto result = apply_cut_and_stretch(plan, 10.0);
        INFO(result.message);
        REQUIRE(bool(result));
        require_solid(result.mesh);
        CHECK(volume(result.mesh) == Approx(volume(source.its) + 10.0 * area).epsilon(1e-5));
        const auto old_arm = arm_extents(source.its, transform, normal, origin);
        const auto new_arm = arm_extents(result.mesh, transform, normal, origin);
        CHECK(new_arm.first - old_arm.first == Approx(10.0).margin(1e-4));
        CHECK(new_arm.second - old_arm.second == Approx(10.0).margin(1e-4));
        CHECK(new_arm.second - new_arm.first == Approx(old_arm.second - old_arm.first).margin(1e-4));
        unchanged(before, source.its);
    }
}

TEST_CASE("Cut and extend rejects a coplanar stationary shelf instead of stretching it", "[cut_and_stretch]")
{
    const auto source = l_arm();
    const auto before = source;
    // Reversing the moving side makes the entire arm the stationary cap, but
    // only the shaft covers it from the other side: B is not a subset of A.
    const auto plan = prepare_cut_and_stretch(source, Vec3d(0, 0, 20), -Vec3d::UnitZ());
    INFO(plan.message);
    failed_plan(plan);
    CHECK(plan.message.find("coplanar shelf") != std::string::npos);
    unchanged(before, source);
}

TEST_CASE("Cut and extend rejects shortening, no-op, malformed plans and tangent cuts", "[cut_and_stretch]")
{
    const auto source = l_arm();
    const auto plan = prepare_cut_and_stretch(source, Vec3d(0, 0, 20), Vec3d::UnitZ());
    REQUIRE(bool(plan));
    for (double distance : {0.0, -10.0, 1e-9, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
        failed_result(apply_cut_and_stretch(plan, distance));
    failed_result(apply_cut_and_stretch(CutAndStretchPlan{}, 10.0));
    auto corrupt = plan;
    corrupt.source_volume_mm3 *= 2;
    failed_result(apply_cut_and_stretch(corrupt, 10.0));
    corrupt = plan;
    corrupt.cut_to_world.linear()(0, 0) = 2;
    failed_result(apply_cut_and_stretch(corrupt, 10.0));
    failed_plan(prepare_cut_and_stretch(source, Vec3d(0, 0, 24), Vec3d::UnitZ()));
    failed_plan(prepare_cut_and_stretch(source, Vec3d::Zero(), Vec3d::Zero()));
}

TEST_CASE("Cut and extend refuses uncertain native-cut features and invalid topology", "[cut_and_stretch]")
{
    const auto source = l_arm();
    // A deliberately non-coplanar near-plane feature must not be swallowed by
    // the native 0.001 mm classifier, even if the final cap could be closed.
    failed_plan(prepare_cut_and_stretch(source, Vec3d(0, 0, 20.0005), Vec3d::UnitZ()));
    auto open = source;
    open.indices.pop_back();
    failed_plan(prepare_cut_and_stretch(open, Vec3d(0, 0, 20), Vec3d::UnitZ()));
    auto broken = source;
    broken.indices[0][0] = -1;
    failed_plan(prepare_cut_and_stretch(broken, Vec3d(0, 0, 20), Vec3d::UnitZ()));
    broken = source;
    broken.vertices[0].x() = std::numeric_limits<float>::quiet_NaN();
    failed_plan(prepare_cut_and_stretch(broken, Vec3d(0, 0, 20), Vec3d::UnitZ()));
    // Borderline native float welding band plus scaled-coordinate rounding:
    // the mid-edge cut point is about 0.00100005 mm from its neighbor.
    indexed_triangle_set narrow_cap;
    narrow_cap.vertices = {{0, 0, 0}, {0.0020001f, 0, 0}, {0.0020001f, 10, 0}, {0, 10, 0}};
    narrow_cap.indices = {{0, 1, 2}, {0, 2, 3}};
    const auto narrow = extrude_planar_face(narrow_cap, {0, 1}, Vec3d(0, 0, 10));
    REQUIRE(bool(narrow));
    const auto near_weld = prepare_cut_and_stretch(narrow.mesh, Vec3d(0, 0, 5), Vec3d::UnitZ());
    failed_plan(near_weld);
    CHECK(near_weld.message.find("welding band") != std::string::npos);
    TriangleMesh distant(source);
    distant.translate(Vec3f(10000, 10000, 0));
    failed_plan(prepare_cut_and_stretch(distant.its, Vec3d(10000, 10000, 20), Vec3d::UnitZ()));
}

TEST_CASE("Cut and extend supports ordinary build-plate world positions", "[cut_and_stretch]")
{
    TriangleMesh source(l_arm());
    source.translate(Vec3f(128, 128, 0));
    const auto plan = prepare_cut_and_stretch(source.its, Vec3d(128, 128, 20), Vec3d::UnitZ());
    INFO(plan.message);
    REQUIRE(bool(plan));
    const auto result = apply_cut_and_stretch(plan, 10.0);
    INFO(result.message);
    REQUIRE(bool(result));
    CHECK(volume(result.mesh) == Approx(4200.0).epsilon(1e-5));
    require_solid(result.mesh);
}

TEST_CASE("Cut and extend positive distance can move the negative world side", "[cut_and_stretch]")
{
    const auto source = hollow_shaft();
    const auto plan = prepare_cut_and_stretch(source, Vec3d(0, 0, 20), -Vec3d::UnitZ());
    INFO(plan.message);
    REQUIRE(bool(plan));
    const auto result = apply_cut_and_stretch(plan, 10.0);
    INFO(result.message);
    REQUIRE(bool(result));
    require_solid(result.mesh, 0);
    CHECK(volume(result.mesh) == Approx(12000.0));
    float low = 1000, high = -1000;
    for (const auto &v : result.mesh.vertices) { low = std::min(low, v.z()); high = std::max(high, v.z()); }
    CHECK(low == Approx(-10.0));
    CHECK(high == Approx(30.0));
}

TEST_CASE("Cut and extend rejects disconnected and self-intersecting source solids", "[cut_and_stretch]")
{
    auto source = l_arm();
    auto second = source;
    for (auto &v : second.vertices) v.x() += 50.f;
    const int offset = int(source.vertices.size());
    source.vertices.insert(source.vertices.end(), second.vertices.begin(), second.vertices.end());
    for (const Vec3i &f : second.indices) source.indices.push_back((f.array() + offset).matrix());
    source.properties.resize(source.indices.size());
    const auto disconnected = prepare_cut_and_stretch(source, Vec3d(0, 0, 20), Vec3d::UnitZ());
    failed_plan(disconnected);
    CHECK(disconnected.message.find("Disconnected") != std::string::npos);

    auto crossing = l_arm();
    crossing.vertices[3] = Vec3f(5, 5, -5);
    const auto solid = MeshBoolean::cgal::triangle_mesh_to_cgal(crossing);
    REQUIRE(MeshBoolean::cgal::does_self_intersect(*solid));
    const auto before = crossing;
    failed_plan(prepare_cut_and_stretch(crossing, Vec3d(0, 0, 20), Vec3d::UnitZ()));
    unchanged(before, crossing);
}

TEST_CASE("Cut and extend preserves tiny sources at normal world coordinates", "[cut_and_stretch]")
{
    TriangleMesh source(l_arm());
    Transform3d transform = Transform3d::Identity();
    transform.linear() *= 0.1;
    transform.translation() = Vec3d(128.123, 128.567, 0.0);
    source.transform(transform, true);
    const auto plan = prepare_cut_and_stretch(source.its, transform * Vec3d(0, 0, 20), Vec3d::UnitZ());
    INFO(plan.message);
    REQUIRE(bool(plan));
    const auto result = apply_cut_and_stretch(plan, 10.125);
    INFO(result.message);
    REQUIRE(bool(result));
    require_solid(result.mesh);
    CHECK(volume(result.mesh) == Approx(plan.source_volume_mm3 + plan.contact_area_mm2 * 10.125).epsilon(1e-4));
}

TEST_CASE("Cut and extend cancellation never publishes partial geometry", "[cut_and_stretch]")
{
    const auto source = l_arm();
    const auto before = source;
    failed_plan(prepare_cut_and_stretch(source, Vec3d(0, 0, 20), Vec3d::UnitZ(), [] { return true; }));
    int progress = 0;
    failed_plan(prepare_cut_and_stretch(source, Vec3d(0, 0, 20), Vec3d::UnitZ(),
        [&] { return progress >= 30; }, [&](int p) { progress = p; }));
    const auto plan = prepare_cut_and_stretch(source, Vec3d(0, 0, 20), Vec3d::UnitZ());
    REQUIRE(bool(plan));
    for (int stop : {0, 20, 50, 80, 100}) {
        progress = 0;
        failed_result(apply_cut_and_stretch(plan, 10.0, [&] { return progress >= stop; }, [&](int p) { progress = p; }));
    }
    failed_plan(prepare_cut_and_stretch(source, Vec3d(0, 0, 20), Vec3d::UnitZ(), {},
        [](int p) { if (p >= 40) throw std::runtime_error("test callback failure inside native cut"); }));
    failed_result(apply_cut_and_stretch(plan, 10.0, {},
        [](int p) { if (p == 100) throw std::runtime_error("test terminal callback failure"); }));
    unchanged(before, source);
}
