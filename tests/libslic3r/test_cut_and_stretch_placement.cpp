// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include <catch2/catch.hpp>
#include "libslic3r/CutAndStretch.hpp"
#include "libslic3r/FaceExtrusion.hpp"
#include "libslic3r/MeshBoolean.hpp"
#include "libslic3r/MeshDiagnostics.hpp"
#include "libslic3r/TriangleMeshSlicer.hpp"
#include <Eigen/LU>
#include <array>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <tuple>

using namespace Slic3r;
namespace {
indexed_triangle_set make_l()
{
    indexed_triangle_set cap;
    cap.vertices = {{0,0,0}, {10,0,0}, {10,0,20}, {40,0,20}, {40,0,24}, {0,0,24}};
    cap.indices = {{0,1,2}, {0,2,5}, {2,3,4}, {2,4,5}};
    const auto out = extrude_planar_face(cap, {0,1,2,3}, Vec3d(0,10,0));
    REQUIRE(out);
    return out.mesh;
}
indexed_triangle_set make_tube()
{
    indexed_triangle_set cap;
    cap.vertices = {{0,0,0}, {20,0,0}, {20,20,0}, {0,20,0}, {5,5,0}, {15,5,0}, {15,15,0}, {5,15,0}};
    cap.indices = {{0,1,5}, {0,5,4}, {1,2,6}, {1,6,5}, {2,3,7}, {2,7,6}, {3,0,4}, {3,4,7}};
    const auto out = extrude_planar_face(cap, {0,1,2,3,4,5,6,7}, Vec3d(0,0,30));
    REQUIRE(out);
    return out.mesh;
}
size_t underside_facet(const indexed_triangle_set &mesh)
{
    size_t result = mesh.indices.size();
    float highest = -std::numeric_limits<float>::infinity();
    for (size_t i = 0; i < mesh.indices.size(); ++i) {
        const auto &f = mesh.indices[i];
        const auto &a = mesh.vertices[f[0]], &b = mesh.vertices[f[1]], &c = mesh.vertices[f[2]];
        if (a.z() == b.z() && a.z() == c.z() && (b-a).cross(c-a).z() < 0 && a.z() > highest) {
            highest = a.z(); result = i;
        }
    }
    REQUIRE(result < mesh.indices.size());
    return result;
}
double volume(const indexed_triangle_set &mesh)
{
    const Vec3d origin = mesh.vertices.front().cast<double>();
    double result = 0;
    for (const auto &f : mesh.indices) {
        const Vec3d a = mesh.vertices[f[0]].cast<double>() - origin;
        const Vec3d b = mesh.vertices[f[1]].cast<double>() - origin;
        const Vec3d c = mesh.vertices[f[2]].cast<double>() - origin;
        result += a.dot(b.cross(c)) / 6.0;
    }
    return result;
}
void solid(const indexed_triangle_set &mesh, int euler)
{
    REQUIRE_FALSE(mesh.empty());
    const auto d = its_mesh_diagnostics(mesh);
    CHECK(d.open_edges == 0); CHECK(d.non_manifold_edges == 0);
    CHECK(d.non_manifold_vertices == 0); CHECK(d.same_direction_edges == 0);
    CHECK_FALSE(d.has_reversed_faces);
    const auto cgal = MeshBoolean::cgal::triangle_mesh_to_cgal(mesh);
    CHECK_FALSE(MeshBoolean::cgal::does_self_intersect(*cgal));
    CHECK(MeshBoolean::cgal::does_bound_a_volume(*cgal));
    CHECK(volume(mesh) > 0.0); // ORIGINAL local winding, even for mirrored placement
    std::set<int> vertices;
    std::set<std::pair<int,int>> edges;
    for (const auto &f : mesh.indices)
        for (int i = 0; i < 3; ++i) {
            vertices.insert(f[i]);
            edges.insert(std::minmax(f[i], f[(i+1)%3]));
        }
    CHECK(int(vertices.size()) - int(edges.size()) + int(mesh.indices.size()) == euler);
}
void unchanged(const indexed_triangle_set &a, const indexed_triangle_set &b)
{
    CHECK(a.vertices == b.vertices);
    CHECK(a.indices == b.indices);
    REQUIRE(a.properties.size() == b.properties.size());
    for (size_t i = 0; i < a.properties.size(); ++i) {
        CHECK(a.properties[i].type == b.properties[i].type);
        CHECK(a.properties[i].area == b.properties[i].area);
    }
}
std::set<std::tuple<float,float,float>> point_set(const indexed_triangle_set &mesh)
{
    std::set<std::tuple<float,float,float>> result;
    for (const auto &p : mesh.vertices) result.emplace(p.x(), p.y(), p.z());
    return result;
}
void failed(const CutAndStretchPlan &p)
{
    CHECK_FALSE(p); CHECK_FALSE(p.message.empty());
    CHECK(p.stationary.empty()); CHECK(p.moving.empty()); CHECK(p.bridge_unit.empty());
    CHECK_FALSE(p.source_frame);
}
void failed(const CutAndStretchResult &r)
{
    CHECK_FALSE(r); CHECK_FALSE(r.message.empty()); CHECK(r.mesh.empty());
}
const std::array<Vec3d,4> placements = {{Vec3d(0,0,0), Vec3d(450,100,0), Vec3d(1024,-1024,0), Vec3d(10000.125,20000.375,23.625)}};
}

TEST_CASE("Placement-aware L cut from its selected underside is invariant on later print plates", "[cut_and_stretch_placement]")
{
    const auto source = make_l(), before = source;
    std::set<std::tuple<float,float,float>> reference;
    for (const auto &offset : placements) {
        Transform3d pose = Transform3d::Identity();
        pose.translation() = offset;
        const Transform3d saved_pose = pose;
        const Vec3d selected_underside = pose * Vec3d(30,5,20);
        const auto plan = prepare_cut_and_stretch_from_face(source, pose, underside_facet(source), Vec3d::UnitZ());
        INFO(offset.transpose()); INFO(plan.message);
        REQUIRE(plan);
        REQUIRE(plan.source_frame);
        CHECK(plan.contact_area_mm2 == Approx(100));
        CHECK(plan.source_volume_mm3 == Approx(3600));
        CHECK(std::abs((plan.cut_to_world.translation() - selected_underside).z()) < 1e-10);
        const auto moving = plan.moving, stationary = plan.stationary, bridge = plan.bridge_unit;
        const auto result = apply_cut_and_stretch_to_source(plan, 10.0);
        INFO(result.message);
        REQUIRE(result);
        solid(result.mesh, 2);
        CHECK(volume(result.mesh) == Approx(4600).margin(1e-5));
        double bottom = 1e10, arm_low = 1e10, arm_high = -1e10;
        for (const auto &p : result.mesh.vertices) {
            bottom = std::min(bottom, double(p.z()));
            if (p.x() > 15) { arm_low = std::min(arm_low, double(p.z())); arm_high = std::max(arm_high, double(p.z())); }
        }
        CHECK(bottom == Approx(0).margin(1e-6));
        CHECK(arm_low == Approx(30).margin(1e-6));
        CHECK(arm_high == Approx(34).margin(1e-6));
        CHECK(arm_high - arm_low == Approx(4).margin(1e-6));
        if (reference.empty()) reference = point_set(result.mesh);
        else CHECK(point_set(result.mesh) == reference); // no world-float placement roundtrip
        CHECK((pose.matrix().array() == saved_pose.matrix().array()).all());
        unchanged(before, source); unchanged(moving, plan.moving); unchanged(stationary, plan.stationary); unchanged(bridge, plan.bridge_unit);
    }
}

TEST_CASE("Placement-aware tube keeps its through hole in the original local frame", "[cut_and_stretch_placement]")
{
    const auto source = make_tube(), before = source;
    for (const auto &offset : placements) for (bool mirrored_pose : {false, true}) {
        Transform3d pose = Transform3d::Identity(); pose.translation() = offset;
        if (mirrored_pose)
            pose.linear() = Eigen::AngleAxisd(0.33, Vec3d(1,2,3).normalized()).toRotationMatrix() * Vec3d(-1.25,0.75,1.5).asDiagonal();
        const auto pose_before = pose;
        // XY of a plane point is immaterial: even an absolute-Z plane at distant
        // world XY=0 must use a computation frame near this source.
        const Vec3d plane = mirrored_pose ? Vec3d(pose * Vec3d(10,10,15)) : Vec3d(0,0,offset.z()+15.0);
        const Vec3d normal = (pose.linear().inverse().transpose() * Vec3d::UnitZ()).normalized();
        const auto plan = prepare_cut_and_stretch(source, pose, plane, normal);
        INFO(plan.message); REQUIRE(plan);
        CHECK(plan.contact_area_mm2 == Approx(mirrored_pose ? 281.25 : 300).epsilon(1e-5));
        const auto out = apply_cut_and_stretch_to_source(plan, 10.125);
        INFO(out.message); REQUIRE(out);
        solid(out.mesh, 0);
        CHECK(volume(out.mesh) == Approx(mirrored_pose ? 11025.0 : 12037.5).margin(0.01));
        const auto sections = slice_mesh_ex(out.mesh, std::vector<float>{5.f,20.f,35.f});
        REQUIRE(sections.size() == 3);
        for (const auto &s : sections) { REQUIRE(s.size() == 1); CHECK(s.front().holes.size() == 1); }
        CHECK((pose.matrix().array() == pose_before.matrix().array()).all());
        unchanged(before, source);
    }
}

TEST_CASE("Placement-aware physical millimeters preserve scaled mirrored and rotated local source poses", "[cut_and_stretch_placement]")
{
    const auto source = make_l(), before = source;
    for (const Vec3d &scale : {Vec3d(1,1,1), Vec3d(1.5,0.75,1.25), Vec3d(-1.5,0.75,1.25), Vec3d(1.5,0.75,-1.25)}) {
        Transform3d volume_pose = Transform3d::Identity();
        volume_pose.linear() = Eigen::AngleAxisd(0.57, Vec3d(1,2,3).normalized()).toRotationMatrix() * scale.asDiagonal();
        volume_pose.translation() = Vec3d(0.375,-0.625,0.125);
        for (const Vec3d &offset : placements) {
            Transform3d instance_pose = Transform3d::Identity(); instance_pose.translation() = offset;
            const auto instance_before = instance_pose, volume_before = volume_pose;
            const Transform3d pose = instance_pose * volume_pose;
            const Vec3d normal = (pose.linear().inverse().transpose() * Vec3d::UnitZ()).normalized();
            const Vec3d origin = pose * Vec3d(30,5,20);
            const auto plan = prepare_cut_and_stretch_from_face(source, pose, underside_facet(source), normal);
            INFO(scale.transpose()); INFO(offset.transpose()); INFO(plan.message); REQUIRE(plan);
            const double determinant = std::abs(pose.linear().determinant());
            const double area = 100.0 * std::abs(scale.x()*scale.y());
            CHECK(plan.contact_area_mm2 == Approx(area).epsilon(2e-5));
            CHECK(plan.source_volume_mm3 == Approx(3600*determinant).epsilon(1e-10));
            CHECK((plan.cut_to_world.linear().col(2) - normal).norm() < 1e-12);
            CHECK(std::abs(normal.dot(plan.cut_to_world.translation() - origin)) < 1e-8);
            const auto out = apply_cut_and_stretch_to_source(plan, 10.125);
            INFO(out.message); REQUIRE(out);
            solid(out.mesh, 2);
            CHECK(volume(out.mesh)*determinant == Approx(3600*determinant + area*10.125).epsilon(2e-5));
            double arm_low = 1e10, arm_high = -1e10, bottom = 1e10;
            for (const auto &v : out.mesh.vertices) {
                const Vec3d p = v.cast<double>();
                bottom = std::min(bottom, p.z()*std::abs(scale.z()));
                if (p.x() > 15) {
                    const double height = normal.dot(pose.linear() * (p - Vec3d(0,0,20)));
                    arm_low = std::min(arm_low, height); arm_high = std::max(arm_high, height);
                }
            }
            CHECK(bottom == Approx(0).margin(1e-4));
            CHECK(arm_low == Approx(10.125).margin(1e-4));
            CHECK(arm_high - arm_low == Approx(4*std::abs(scale.z())).margin(1e-4));
            CHECK((instance_pose.matrix().array() == instance_before.matrix().array()).all());
            CHECK((volume_pose.matrix().array() == volume_before.matrix().array()).all());
            unchanged(before, source);
        }
    }
}

TEST_CASE("Explicit source-facet provenance does not authorize snapping generic near-world planes", "[cut_and_stretch_placement]")
{
    const auto source = make_l(), before = source;
    Transform3d pose = Transform3d::Identity();
    pose.linear() = Eigen::AngleAxisd(0.57, Vec3d(1,2,3).normalized()).toRotationMatrix() * Vec3d(1.5,0.75,1.25).asDiagonal();
    pose.translation() = Vec3d(1024.375,-1024.625,0.125);
    const Vec3d normal = (pose.linear().inverse().transpose() * Vec3d::UnitZ()).normalized();
    const size_t facet = underside_facet(source);
    const auto exact = prepare_cut_and_stretch_from_face(source, pose, facet, normal);
    INFO(exact.message); REQUIRE(exact);
    const Vec3d approximate_world_point = pose * Vec3d(30,5,20);
    failed(prepare_cut_and_stretch(source, pose, approximate_world_point, normal));
    for (double offset : {1e-12, 1e-9, 1e-7, 0.0005})
        failed(prepare_cut_and_stretch(source, pose, approximate_world_point + offset*normal, normal));
    failed(prepare_cut_and_stretch_from_face(source, pose, source.indices.size(), normal));
    failed(prepare_cut_and_stretch_from_face(source, pose, facet, Vec3d::UnitX()));
    auto broken = source; broken.indices[facet][1] = broken.indices[facet][0];
    failed(prepare_cut_and_stretch_from_face(broken, pose, facet, normal));
    unchanged(before, source);
}

TEST_CASE("Placement-aware local backcast uses physical error and never silently loses a fractional extension", "[cut_and_stretch_placement]")
{
    auto source = make_l();
    for (auto &p : source.vertices) {
        p *= 1000.0f;
        p.z() += 100000000.0f; // local Z float ULP=8, physical ULP=0.008 mm
    }
    const auto before = source;
    Transform3d pose = Transform3d::Identity();
    pose.linear() *= 0.001;
    pose.translation() = Vec3d(10000,20000,-100000);
    const auto plan = prepare_cut_and_stretch_from_face(source, pose, underside_facet(source), Vec3d::UnitZ());
    INFO(plan.message); REQUIRE(plan);
    // Exact storage is permitted; a fractional displacement is not rounded away.
    const auto exact = apply_cut_and_stretch_to_source(plan, 10.0);
    INFO(exact.message); REQUIRE(exact);
    solid(exact.mesh, 2);
    CHECK(volume(exact.mesh)*1e-9 == Approx(4600).margin(1e-4));
    const auto rejected = apply_cut_and_stretch_to_source(plan, 10.125);
    failed(rejected);
    CHECK(rejected.message.find("precision") != std::string::npos);
    unchanged(before, source);
}

TEST_CASE("Placement-aware cut retains intrinsic precision topology and provenance guards", "[cut_and_stretch_placement]")
{
    const auto source = make_l();
    Transform3d pose = Transform3d::Identity(); pose.translation() = Vec3d(10000,20000,0);
    const Vec3d plane = pose * Vec3d(30,5,20);
    const auto plan = prepare_cut_and_stretch(source, pose, plane, Vec3d::UnitZ());
    REQUIRE(plan);
    for (double distance : {0.0, -1.0, 1e-10, 1e8, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()})
        failed(apply_cut_and_stretch_to_source(plan, distance));
    auto changed = plan; changed.cut_to_world.translation().x() += 1;
    failed(apply_cut_and_stretch_to_source(changed, 10));
    changed = plan; changed.source_frame.reset();
    failed(apply_cut_and_stretch_to_source(changed, 10));
    failed(apply_cut_and_stretch_to_source(CutAndStretchPlan{}, 10));
    const auto old_plan = prepare_cut_and_stretch(source, Vec3d(30,5,20), Vec3d::UnitZ());
    REQUIRE(old_plan);
    failed(apply_cut_and_stretch_to_source(old_plan, 10));
    // The old world API still rejects distant world floats: it has NOT simply
    // had its global precision guard removed to make the new path pass.
    TriangleMesh world(source); world.transform(pose, true);
    failed(prepare_cut_and_stretch(world.its, plane, Vec3d::UnitZ()));
    failed(apply_cut_and_stretch(plan, 10));
    auto large = source; for (auto &p : large.vertices) p *= 100.0f;
    failed(prepare_cut_and_stretch(large, pose, pose * Vec3d(3000,500,2000), Vec3d::UnitZ()));
    failed(prepare_cut_and_stretch(source, pose, plane, -Vec3d::UnitZ())); // unsupported stationary shelf
    failed(prepare_cut_and_stretch(source, pose, plane + Vec3d(0,0,0.0005), Vec3d::UnitZ()));
    failed(prepare_cut_and_stretch(source, pose, pose * Vec3d(30,5,24), Vec3d::UnitZ())); // tangent, no CAD fallback
    auto invalid_pose = pose; invalid_pose.linear().col(0).setZero();
    failed(prepare_cut_and_stretch(source, invalid_pose, plane, Vec3d::UnitZ()));
    invalid_pose = pose; invalid_pose.translation().x() = std::numeric_limits<double>::infinity();
    failed(prepare_cut_and_stretch(source, invalid_pose, plane, Vec3d::UnitZ()));
    invalid_pose = pose; invalid_pose.matrix()(3,0) = 1;
    failed(prepare_cut_and_stretch(source, invalid_pose, plane, Vec3d::UnitZ()));
    failed(prepare_cut_and_stretch(source, pose, plane, Vec3d::Zero()));
    failed(prepare_cut_and_stretch(source, pose, Vec3d(0,0,std::numeric_limits<double>::quiet_NaN()), Vec3d::UnitZ()));
    auto broken = source; broken.indices.pop_back();
    failed(prepare_cut_and_stretch(broken, pose, plane, Vec3d::UnitZ()));
}

TEST_CASE("Placement-aware cancellation clears operands and local results without mutating source or plan", "[cut_and_stretch_placement]")
{
    const auto source = make_l(), before = source;
    Transform3d pose = Transform3d::Identity(); pose.translation() = Vec3d(1024,-1024,0);
    const auto pose_before = pose;
    const Vec3d plane = pose * Vec3d(30,5,20);
    failed(prepare_cut_and_stretch(source, pose, plane, Vec3d::UnitZ(), [] { return true; }));
    int progress = 0;
    failed(prepare_cut_and_stretch(source, pose, plane, Vec3d::UnitZ(), [&] { return progress >= 30; }, [&](int p) { progress = p; }));
    const auto plan = prepare_cut_and_stretch(source, pose, plane, Vec3d::UnitZ());
    REQUIRE(plan);
    const auto lower = plan.stationary, upper = plan.moving, bridge = plan.bridge_unit;
    for (int stop : {0,20,80,100}) {
        progress = 0;
        failed(apply_cut_and_stretch_to_source(plan, 10, [&] { return progress >= stop; }, [&](int p) { progress = p; }));
    }
    failed(prepare_cut_and_stretch(source, pose, plane, Vec3d::UnitZ(), {},
        [](int p) { if (p >= 30) throw std::runtime_error("synthetic preparation callback error"); }));
    failed(apply_cut_and_stretch_to_source(plan, 10, {},
        [](int p) { if (p == 100) throw std::runtime_error("synthetic local-result callback error"); }));
    unchanged(before, source); unchanged(lower, plan.stationary); unchanged(upper, plan.moving); unchanged(bridge, plan.bridge_unit);
    CHECK((pose.matrix().array() == pose_before.matrix().array()).all());
}

TEST_CASE("An intrinsic oblique source facet requires its actual double normal rather than a promoted float measurement normal", "[cut_and_stretch_placement]")
{
    auto source = make_l();
    const size_t facet = underside_facet(source); // preserve the actual facet identity before baking the shear
    for (Vec3f &p : source.vertices) p.z() += 0.25f * p.x() + 0.75f * p.y();
    const auto before = source; // all sheared vertex coordinates remain binary-exact
    Transform3d pose = Transform3d::Identity();
    pose.translation() = Vec3d(1024.125,-1024.375,0);
    const auto pose_before = pose;
    const Vec3i &f = source.indices[facet];
    const Vec3f &a = source.vertices[f[0]], &b = source.vertices[f[1]], &c = source.vertices[f[2]];
    // Measurement caches normalize a Vec3f cross product, then promote it to
    // double. Promotion cannot recover its lost angular precision.
    const Vec3f measured = -(b-a).cross(c-a).normalized();
    const Vec3d measured_world = (pose.linear().inverse().transpose() * measured.cast<double>()).normalized();
    const Vec3d precise_source = -(b.cast<double>()-a.cast<double>()).cross(c.cast<double>()-a.cast<double>()).normalized();
    const Vec3d precise_world = (pose.linear().inverse().transpose() * precise_source).normalized();
    CHECK(precise_source.isApprox(Vec3d(-0.25,-0.75,1).normalized(), 1e-12));
    CHECK(measured_world.cross(precise_world).norm() > 1e-10);
    const auto rejected = prepare_cut_and_stretch_from_face(source, pose, facet, measured_world);
    failed(rejected);
    CHECK(rejected.message.find("parallel") != std::string::npos);
    const auto plan = prepare_cut_and_stretch_from_face(source, pose, facet, precise_world);
    INFO(plan.message); REQUIRE(plan);
    const double physical_contact = Vec3d(10,0,2.5).cross(Vec3d(0,10,7.5)).norm();
    CHECK(plan.contact_area_mm2 == Approx(physical_contact).margin(1e-4));
    CHECK(plan.source_volume_mm3 == Approx(3600).margin(1e-5));
    const double distance = 10.0;
    const auto out = apply_cut_and_stretch_to_source(plan, distance);
    INFO(out.message); REQUIRE(out);
    solid(out.mesh, 2);
    CHECK(volume(out.mesh) == Approx(3600 + physical_contact*distance/std::abs(pose.linear().determinant())).margin(0.01));
    unchanged(before, source);
    CHECK((pose.matrix().array() == pose_before.matrix().array()).all());
}
