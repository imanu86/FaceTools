// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

// DRAFT ONLY. Not registered in CMake and not executed by the authoring agent.
// Compile this SAME file against the frozen neutral wrapper for native RED,
// then against the candidate. No user data, archive, app or filesystem fixture.
#include <catch2/catch.hpp>

#include "libslic3r/CutAndStretch.hpp"
#include "libslic3r/FaceExtrusion.hpp"
#include "libslic3r/MeshBoolean.hpp"
#include "libslic3r/MeshDiagnostics.hpp"
#include "libslic3r/PlanarProfile.hpp"
#include "slic3r/GUI/Gizmos/CutAndStretchScopeConfirmation.hpp" // parent-owned production helper, installed before native RED

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <vector>

using namespace Slic3r;

namespace {

indexed_triangle_set profile_solid(const std::vector<Vec2d> &xz)
{
    std::vector<Vec3d> outline;
    for (const Vec2d &p : xz) outline.emplace_back(p.x(), 0.0, p.y());
    auto generated = extrude_planar_profile(outline, Vec3d(0, 10, 0));
    INFO(generated.message);
    REQUIRE(bool(generated));
    return std::move(generated.mesh);
}

indexed_triangle_set u_source()
{
    // Outer 40x30, two 10-wide legs, base height 5, depth 10: volume 7000.
    return profile_solid({{0,0},{40,0},{40,30},{30,30},{30,5},{10,5},{10,30},{0,30}});
}

indexed_triangle_set l_source()
{
    // Stem 10x10x20 plus 40x10x4 arm: volume 3600.
    return profile_solid({{0,0},{10,0},{10,20},{40,20},{40,24},{0,24}});
}

indexed_triangle_set t_source()
{
    return profile_solid({{-5,0},{5,0},{5,20},{20,20},{20,24},{-20,24},{-20,20},{-5,20}});
}

indexed_triangle_set roof_source()
{
    // Left leg ends at 20. Right leg supports a roof at 30..35 over the left
    // leg. Above Z15 there are TWO disjoint components, one enclosing L-shaped
    // AABB: +5 is safe, +10 touches roof, +12 overlaps, +30 hits it with bridge.
    return profile_solid({{0,0},{40,0},{40,35},{0,35},{0,30},{30,30},
                          {30,5},{10,5},{10,20},{0,20}});
}

double volume(const indexed_triangle_set &mesh)
{
    REQUIRE_FALSE(mesh.vertices.empty());
    const Vec3d origin = mesh.vertices.front().cast<double>();
    double out = 0.0;
    for (const Vec3i &f : mesh.indices) {
        const Vec3d a = mesh.vertices[f[0]].cast<double>() - origin;
        const Vec3d b = mesh.vertices[f[1]].cast<double>() - origin;
        const Vec3d c = mesh.vertices[f[2]].cast<double>() - origin;
        out += a.dot(b.cross(c)) / 6.0;
    }
    return out;
}

void unchanged(const indexed_triangle_set &a, const indexed_triangle_set &b)
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

void failed(const CutAndStretchPlan &plan)
{
    CHECK_FALSE(bool(plan));
    CHECK_FALSE(plan.message.empty());
    CHECK(plan.stationary.empty());
    CHECK(plan.moving.empty());
    CHECK(plan.bridge_unit.empty());
    CHECK_FALSE(bool(plan.guided_state));
}

void failed(const CutAndStretchResult &result)
{
    CHECK_FALSE(bool(result));
    CHECK_FALSE(result.message.empty());
    CHECK(result.mesh.empty());
}

void quality(const indexed_triangle_set &mesh)
{
    REQUIRE_FALSE(mesh.empty());
    indexed_triangle_set welded;
    std::map<std::array<float,3>, int> vertex_ids;
    std::set<std::array<int,3>> face_keys;
    for (const Vec3i &f : mesh.indices) {
        Vec3i face;
        for (int i = 0; i < 3; ++i) {
            REQUIRE(f[i] >= 0);
            REQUIRE(size_t(f[i]) < mesh.vertices.size());
            const Vec3f &p = mesh.vertices[f[i]];
            REQUIRE(p.allFinite());
            const std::array<float,3> key{{p.x(), p.y(), p.z()}};
            const auto found = vertex_ids.emplace(key, int(welded.vertices.size()));
            if (found.second) welded.vertices.push_back(p);
            face[i] = found.first->second;
        }
        REQUIRE(face[0] != face[1]); REQUIRE(face[1] != face[2]); REQUIRE(face[2] != face[0]);
        const Vec3d a = welded.vertices[face[0]].cast<double>();
        const Vec3d b = welded.vertices[face[1]].cast<double>();
        const Vec3d c = welded.vertices[face[2]].cast<double>();
        REQUIRE((b-a).cross(c-a).squaredNorm() > 0.0);
        std::array<int,3> key{{face[0],face[1],face[2]}};
        std::sort(key.begin(), key.end());
        REQUIRE(face_keys.insert(key).second); // duplicate faces AFTER exact diagnostic weld
        welded.indices.push_back(face);
    }
    const auto d = its_mesh_diagnostics(welded);
    CHECK(d.open_edges == 0); CHECK(d.non_manifold_edges == 0);
    CHECK(d.non_manifold_vertices == 0); CHECK(d.same_direction_edges == 0);
    CHECK_FALSE(d.has_reversed_faces);
    const auto solid = MeshBoolean::cgal::triangle_mesh_to_cgal(welded);
    CHECK_FALSE(MeshBoolean::cgal::does_self_intersect(*solid));
    CHECK(MeshBoolean::cgal::does_bound_a_volume(*solid));
    CHECK(volume(welded) > 0.0);
}

template<class Predicate>
std::vector<int> facets(const indexed_triangle_set &mesh, Predicate predicate)
{
    std::vector<int> out;
    for (size_t i = 0; i < mesh.indices.size(); ++i) {
        const Vec3i &f = mesh.indices[i];
        if (predicate(mesh.vertices[f[0]].cast<double>(), mesh.vertices[f[1]].cast<double>(),
                      mesh.vertices[f[2]].cast<double>())) out.push_back(int(i));
    }
    REQUIRE_FALSE(out.empty());
    return out;
}

std::vector<int> top_guides(const indexed_triangle_set &source, bool left, double height = 30.0)
{
    return facets(source, [=](const Vec3d &a, const Vec3d &b, const Vec3d &c) {
        return a.z() == height && b.z() == height && c.z() == height &&
               (left ? std::max({a.x(),b.x(),c.x()}) <= 10.0 : std::min({a.x(),b.x(),c.x()}) >= 30.0);
    });
}

std::vector<int> plane_facets(const indexed_triangle_set &source, double z)
{
    return facets(source, [=](const Vec3d &a, const Vec3d &b, const Vec3d &c) {
        return a.z() == z && b.z() == z && c.z() == z;
    });
}

double top(const indexed_triangle_set &mesh, bool left)
{
    double z = -std::numeric_limits<double>::infinity();
    for (const Vec3f &p : mesh.vertices)
        if (left ? p.x() <= 10.f : p.x() >= 30.f) z = std::max(z, double(p.z()));
    REQUIRE(std::isfinite(z));
    return z;
}

void same_solid_geometry(indexed_triangle_set a, indexed_triangle_set b)
{
    auto a_minus_b = a;
    MeshBoolean::cgal::minus(a_minus_b, b);
    CHECK(a_minus_b.indices.empty());
    auto b_minus_a = b;
    MeshBoolean::cgal::minus(b_minus_a, a);
    CHECK(b_minus_a.indices.empty());
}

void unchanged_other_leg(const indexed_triangle_set &before, const indexed_triangle_set &after, bool moved_left)
{
    // Independent spatial CSG oracle: includes the other leg's ORIGINAL Z15
    // crossing, so mere unchanged bounds/volume or reused cap triangulation is
    // insufficient. Integer fixture coordinates permit exact set comparison.
    TriangleMesh window = make_cube(20, 12, 35);
    window.translate(moved_left ? 25.0 : -5.0, -1.0, 14.0);
    auto a = before, b = after;
    MeshBoolean::cgal::intersect(a, window.its);
    MeshBoolean::cgal::intersect(b, window.its);
    REQUIRE_FALSE(a.empty()); REQUIRE_FALSE(b.empty());
    CHECK(volume(a) == Approx(1600.0));
    CHECK(volume(b) == Approx(1600.0));
    same_solid_geometry(std::move(a), std::move(b));
}

void confirmation_contract(const std::function<CutAndStretchPlan()> &prepare,
                           CutAndStretchScope scope)
{
    const auto real_plan = [&]() -> std::shared_ptr<const CutAndStretchPlan> {
        auto prepared = prepare(); // actual kernel geometry, never a fabricated valid flag
        INFO(prepared.message);
        REQUIRE(bool(prepared));
        REQUIRE(prepared.scope == scope);
        REQUIRE_FALSE(prepared.moving.empty());
        REQUIRE_FALSE(prepared.stationary.empty());
        REQUIRE_FALSE(prepared.bridge_unit.empty());
        return std::make_shared<const CutAndStretchPlan>(std::move(prepared));
    };
    Slic3r::GUI::CutAndStretchScopeConfirmation consent;
    const auto wrong_scope = scope == CutAndStretchScope::GuidedComponent ?
        CutAndStretchScope::WholeHalfSpace : CutAndStretchScope::GuidedComponent;
    auto current = real_plan();
    CHECK_FALSE(consent.matches(current, scope)); // plane preparation is NOT region consent
    CHECK_FALSE(consent.confirm(current, wrong_scope));
    CHECK_FALSE(consent.matches(current, scope));
    CHECK(consent.confirm(current, scope));
    CHECK(consent.matches(current, scope));
    CHECK_FALSE(consent.matches(current, wrong_scope));
    CHECK_FALSE(consent.confirm(current, wrong_scope)); // failed confirmation clears old consent
    CHECK_FALSE(consent.matches(current, scope));
    REQUIRE(consent.confirm(current, scope));
    auto next = real_plan(); // same geometry, a genuinely rebuilt immutable region
    REQUIRE(current != next);
    CHECK_FALSE(consent.matches(next, scope));
    CHECK(consent.matches(current, scope));
    consent.clear();
    CHECK_FALSE(consent.matches(current, scope));
    REQUIRE(consent.confirm(current, scope));
    const std::weak_ptr<const CutAndStretchPlan> lifetime = current;
    current.reset();
    CHECK(lifetime.expired()); // helper must not keep old prepared geometry alive
    CHECK_FALSE(consent.matches(lifetime.lock(), scope));
    CHECK_FALSE(consent.matches(next, scope));
    CHECK(consent.confirm(next, scope));
    CHECK(consent.matches(next, scope));
    consent.clear();
    CHECK_FALSE(consent.matches(next, scope));
}

} // namespace

TEST_CASE("Guided U extends only the guide-associated upper leg", "[cut_guided_scope][cut_guided_red]")
{
    for (bool left : {true, false}) {
        const auto source = u_source(), before = source;
        const auto guides = top_guides(source, left);
        const auto plan = prepare_cut_and_stretch_guided(source, Transform3d::Identity(),
            Vec3d(20,5,15), Vec3d::UnitZ(), guides);
        INFO(plan.message);
        // The frozen neutral adapter fails HERE with the ACTUAL legacy
        // Disconnected result. Do not weaken this into an expected rejection.
        REQUIRE(bool(plan));
        CHECK(plan.scope == CutAndStretchScope::GuidedComponent);
        CHECK(plan.selection.upper_component_count == 2);
        CHECK(plan.selection.selected_upper_component < 2);
        CHECK(plan.selection.guide_facet_count == guides.size());
        CHECK(plan.selection.selected_source_facet_count > 0);
        CHECK(plan.source_volume_mm3 == Approx(7000));
        CHECK(plan.contact_area_mm2 == Approx(100));
        CHECK(volume(plan.stationary) == Approx(5500));
        CHECK(volume(plan.moving) == Approx(1500));
        const auto stationary_before = plan.stationary, moving_before = plan.moving, bridge_before = plan.bridge_unit;
        const auto result = apply_cut_and_stretch_to_source(plan, 10.0);
        INFO(result.message); REQUIRE(bool(result));
        quality(result.mesh);
        CHECK(volume(result.mesh) == Approx(8000).epsilon(1e-6));
        CHECK(top(result.mesh, left) == Approx(40));
        CHECK(top(result.mesh, !left) == Approx(30));
        float bottom = std::numeric_limits<float>::infinity();
        for (const Vec3f &p : result.mesh.vertices) bottom = std::min(bottom, p.z());
        CHECK(bottom == 0.f);
        unchanged_other_leg(source, result.mesh, left);
        unchanged(before, source);
        unchanged(stationary_before, plan.stationary); unchanged(moving_before, plan.moving);
        unchanged(bridge_before, plan.bridge_unit);
    }
}

TEST_CASE("Guided resolution consumes every original facet not only the first seed", "[cut_guided_scope]")
{
    const auto source = u_source(), before = source;
    auto guides = top_guides(source, true);
    const auto right = top_guides(source, false);
    guides.insert(guides.end(), right.begin(), right.end());
    auto plan = prepare_cut_and_stretch_guided(source, Transform3d::Identity(), Vec3d(20,5,15), Vec3d::UnitZ(), guides);
    failed(plan);
    CHECK(plan.message.find("different upper components") != std::string::npos);
    // Same first ID/seed, but ONE extra source facet must change acceptance.
    guides = top_guides(source, true);
    guides.push_back(right.back());
    plan = prepare_cut_and_stretch_guided(source, Transform3d::Identity(), Vec3d(20,5,15), Vec3d::UnitZ(), guides);
    failed(plan);
    CHECK(plan.message.find("different upper components") != std::string::npos);
    unchanged(before, source);
}

TEST_CASE("Guided straddling facets resolve their upper fragments uniquely", "[cut_guided_scope]")
{
    const auto source = u_source();
    const auto guides = facets(source, [](const Vec3d &a, const Vec3d &b, const Vec3d &c) {
        return a.x() == 0.0 && b.x() == 0.0 && c.x() == 0.0 &&
               std::min({a.z(),b.z(),c.z()}) < 15 && std::max({a.z(),b.z(),c.z()}) > 15;
    });
    const auto plan = prepare_cut_and_stretch_guided(source, Transform3d::Identity(), Vec3d(20,5,15), Vec3d::UnitZ(), guides);
    INFO(plan.message); REQUIRE(bool(plan));
    CHECK(plan.scope == CutAndStretchScope::GuidedComponent);
    const auto result = apply_cut_and_stretch_to_source(plan, 10.0);
    INFO(result.message); REQUIRE(bool(result));
    CHECK(top(result.mesh, true) == Approx(40)); CHECK(top(result.mesh, false) == Approx(30));
}

TEST_CASE("Guided invalid fixed and duplicate masks fail closed", "[cut_guided_scope]")
{
    const auto source = u_source(), before = source;
    const auto left = top_guides(source, true);
    const auto bottom = plane_facets(source, 0);
    for (const auto &ids : std::vector<std::vector<int>>{
            {}, {-1}, {int(source.indices.size())}, {left.front(), left.front()},
            {left.front(), bottom.front()}, bottom}) {
        const auto plan = prepare_cut_and_stretch_guided(source, Transform3d::Identity(),
            Vec3d(20,5,15), Vec3d::UnitZ(), ids);
        failed(plan);
    }
    unchanged(before, source);
}

TEST_CASE("Guided original coplanar L facets use adjacency rather than nearest geometry", "[cut_guided_scope]")
{
    const auto source = l_source(), before = source;
    const auto guides = plane_facets(source, 20.0); // original underside is omitted by native uncapped cut
    const auto plan = prepare_cut_and_stretch_guided_from_face(source, Transform3d::Identity(),
        size_t(guides.front()), Vec3d::UnitZ(), guides);
    INFO(plan.message); REQUIRE(bool(plan));
    CHECK(plan.scope == CutAndStretchScope::GuidedComponent);
    CHECK(plan.selection.upper_component_count == 1);
    CHECK(plan.selection.guide_facet_count == guides.size());
    CHECK(plan.contact_area_mm2 == Approx(100));
    const auto out = apply_cut_and_stretch_to_source(plan, 10.0);
    INFO(out.message); REQUIRE(bool(out)); quality(out.mesh);
    CHECK(volume(out.mesh) == Approx(4600).epsilon(1e-6));
    CHECK(top(out.mesh, false) == Approx(34));
    const auto fixed_guide = prepare_cut_and_stretch_guided_from_face(source, Transform3d::Identity(),
        size_t(guides.front()), -Vec3d::UnitZ(), guides);
    failed(fixed_guide);
    CHECK(fixed_guide.message.find("fixed material side") != std::string::npos);
    unchanged(before, source);
}

TEST_CASE("Guided partial lower-cap overlap rejects without clipping contact triangles", "[cut_guided_scope]")
{
    // Centered pedestal under a wide plate. At Z10, moving DOWN selects the
    // narrow pedestal, whose 4x4 footprint is strictly inside the 10x10 plate
    // cap and only partly covers its large triangles. No cap clipping is allowed.
    TriangleMesh pedestal = make_cube(4,4,10);
    pedestal.translate(3,3,0);
    TriangleMesh plate = make_cube(10,10,2);
    plate.translate(0,0,10);
    auto source = pedestal.its;
    MeshBoolean::cgal::plus(source, plate.its);
    quality(source);
    const auto before = source;
    const auto guides = plane_facets(source, 0);
    const auto plan = prepare_cut_and_stretch_guided(source, Transform3d::Identity(),
        Vec3d(5,5,10), -Vec3d::UnitZ(), guides);
    INFO(plan.message); failed(plan);
    CHECK(plan.message.find("partially overlaps") != std::string::npos);
    unchanged(before, source);
}

TEST_CASE("Guided T identifies the whole connected bar not an invented CAD branch", "[cut_guided_scope]")
{
    const auto source = t_source();
    const auto top_faces = plane_facets(source, 24);
    // One original facet is only a guide; the geometry result MUST expose both
    // arms. Fresh scope confirmation is a parent-owned GUI gate, not this API.
    const auto plan = prepare_cut_and_stretch_guided(source, Transform3d::Identity(),
        Vec3d(0,5,20), Vec3d::UnitZ(), {top_faces.front()});
    INFO(plan.message); REQUIRE(bool(plan));
    CHECK(plan.scope == CutAndStretchScope::GuidedComponent);
    CHECK(plan.selection.upper_component_count == 1);
    double min_x = 1000, max_x = -1000;
    for (const Vec3f &p : plan.moving.vertices) {
        const Vec3d world = plan.cut_to_world * p.cast<double>();
        min_x = std::min(min_x, world.x()); max_x = std::max(max_x, world.x());
    }
    CHECK(min_x == Approx(-20)); CHECK(max_x == Approx(20));
    const auto out = apply_cut_and_stretch_to_source(plan, 10);
    INFO(out.message); REQUIRE(bool(out)); quality(out.mesh);
    for (double x : {-20.0,20.0}) {
        double maximum = -1000;
        for (const Vec3f &p : out.mesh.vertices)
            if (std::abs(double(p.x()) - x) < 1e-6) maximum = std::max(maximum, double(p.z()));
        CHECK(maximum == Approx(34));
    }
}

TEST_CASE("Actual guided T region requires fresh identity-bound scope confirmation", "[cut_guided_scope][cut_scope_confirmation]")
{
    const auto source = t_source(), before = source;
    const auto top_faces = plane_facets(source, 24);
    confirmation_contract([&] {
        return prepare_cut_and_stretch_guided(source, Transform3d::Identity(),
            Vec3d(0,5,20), Vec3d::UnitZ(), {top_faces.front()});
    }, CutAndStretchScope::GuidedComponent);
    unchanged(before, source);
}

TEST_CASE("Actual legacy L region requires explicit whole-half scope confirmation", "[cut_guided_scope][cut_scope_confirmation]")
{
    const auto source = l_source(), before = source;
    confirmation_contract([&] {
        return prepare_cut_and_stretch(source, Transform3d::Identity(),
            Vec3d(0,5,20), Vec3d::UnitZ());
    }, CutAndStretchScope::WholeHalfSpace);
    unchanged(before, source);
}

TEST_CASE("Guided collision guard refuses moving and bridge contact with excluded regions", "[cut_guided_scope]")
{
    const auto source = roof_source(), before = source;
    quality(source);
    const auto plan = prepare_cut_and_stretch_guided(source, Transform3d::Identity(),
        Vec3d(20,5,15), Vec3d::UnitZ(), top_guides(source, true, 20));
    INFO(plan.message); REQUIRE(bool(plan));
    CHECK(plan.selection.upper_component_count == 2);
    const auto stationary = plan.stationary, moving = plan.moving, bridge = plan.bridge_unit;
    // Overlapping AABBs alone must NOT reject this geometrically disjoint move.
    const auto safe = apply_cut_and_stretch_to_source(plan, 5);
    INFO(safe.message); REQUIRE(bool(safe)); quality(safe.mesh);
    CHECK(volume(safe.mesh) == Approx(volume(source) + 500).epsilon(1e-6));
    for (double distance : {10.0,12.0,30.0}) {
        const auto out = apply_cut_and_stretch_to_source(plan, distance);
        INFO(distance); INFO(out.message); failed(out);
        CHECK(out.message.find("excluded upper component") != std::string::npos);
    }
    unchanged(before, source); unchanged(stationary, plan.stationary);
    unchanged(moving, plan.moving); unchanged(bridge, plan.bridge_unit);
}

TEST_CASE("Guided placement preserves world millimeters source pose and the other leg", "[cut_guided_scope]")
{
    const auto source = u_source(), before = source;
    for (double mirror : {1.0,-1.0}) {
        Transform3d pose = Transform3d::Identity();
        pose.linear() = Eigen::AngleAxisd(0.37, Vec3d(1,2,3).normalized()).toRotationMatrix() *
                        Vec3d(mirror*1.25,0.8,1.5).asDiagonal();
        pose.translation() = Vec3d(15000,-24000,12.5);
        const auto pose_before = pose;
        const Vec3d n = (pose.linear().inverse().transpose() * Vec3d::UnitZ()).normalized();
        const auto plan = prepare_cut_and_stretch_guided(source, pose, pose*Vec3d(20,5,15), n,
                                                        top_guides(source, true));
        INFO(plan.message); REQUIRE(bool(plan));
        CHECK(plan.scope == CutAndStretchScope::GuidedComponent);
        CHECK(plan.contact_area_mm2 == Approx(100.0).epsilon(1e-5));
        const auto out = apply_cut_and_stretch_to_source(plan, 10.0);
        INFO(out.message); REQUIRE(bool(out)); quality(out.mesh);
        CHECK(top(out.mesh, true) == Approx(30.0+10.0/1.5).margin(1e-4));
        CHECK(top(out.mesh, false) == Approx(30.0).margin(1e-4));
        CHECK(volume(out.mesh)*std::abs(pose.linear().determinant()) ==
              Approx(volume(source)*std::abs(pose.linear().determinant()) + 1000).epsilon(2e-5));
        CHECK((pose.matrix().array() == pose_before.matrix().array()).all());
        unchanged(before, source);
    }
}

TEST_CASE("Guided cancellation and corrupt scope never publish partial results", "[cut_guided_scope]")
{
    const auto source = u_source(), before = source;
    const auto ids = top_guides(source, true);
    failed(prepare_cut_and_stretch_guided(source, Transform3d::Identity(), Vec3d(20,5,15), Vec3d::UnitZ(), ids,
                                         [] { return true; }));
    int progress = 0;
    failed(prepare_cut_and_stretch_guided(source, Transform3d::Identity(), Vec3d(20,5,15), Vec3d::UnitZ(), ids,
        [&] { return progress >= 60; }, [&](int p) { progress = p; }));
    failed(prepare_cut_and_stretch_guided(source, Transform3d::Identity(), Vec3d(20,5,15), Vec3d::UnitZ(), ids, {},
        [](int p) { if (p >= 55) throw std::runtime_error("synthetic guided callback failure"); }));
    const auto plan = prepare_cut_and_stretch_guided(source, Transform3d::Identity(), Vec3d(20,5,15), Vec3d::UnitZ(), ids);
    INFO(plan.message); REQUIRE(bool(plan));
    for (int stop : {0,20,50,80,100}) {
        progress = 0;
        failed(apply_cut_and_stretch_to_source(plan, 10, [&] { return progress >= stop; }, [&](int p) { progress = p; }));
    }
    auto corrupt = plan;
    corrupt.scope = CutAndStretchScope::WholeHalfSpace;
    failed(apply_cut_and_stretch_to_source(corrupt, 10));
    corrupt = plan; corrupt.guided_state.reset();
    failed(apply_cut_and_stretch_to_source(corrupt, 10));
    corrupt = plan; ++corrupt.selection.upper_component_count;
    failed(apply_cut_and_stretch_to_source(corrupt, 10));
    unchanged(before, source);
}

TEST_CASE("Explicit legacy whole-half contracts are unchanged by guided APIs", "[cut_guided_scope][cut_and_stretch]")
{
    const auto u = u_source();
    const auto legacy_u = prepare_cut_and_stretch(u, Transform3d::Identity(), Vec3d(20,5,15), Vec3d::UnitZ());
    failed(legacy_u);
    CHECK(legacy_u.message.find("Disconnected") != std::string::npos);
    const auto l = l_source();
    const auto legacy_l = prepare_cut_and_stretch(l, Transform3d::Identity(), Vec3d(0,0,20), Vec3d::UnitZ());
    INFO(legacy_l.message); REQUIRE(bool(legacy_l));
    CHECK(legacy_l.scope == CutAndStretchScope::WholeHalfSpace);
    CHECK_FALSE(bool(legacy_l.guided_state));
    const auto out = apply_cut_and_stretch_to_source(legacy_l, 10);
    INFO(out.message); REQUIRE(bool(out)); quality(out.mesh);
    CHECK(volume(out.mesh) == Approx(4600).epsilon(1e-6));
    const auto reverse = prepare_cut_and_stretch(l, Transform3d::Identity(), Vec3d(0,0,20), -Vec3d::UnitZ());
    failed(reverse);
    CHECK(reverse.message.find("coplanar shelf") != std::string::npos);
    auto multishell = u;
    auto second = u;
    for (Vec3f &p : second.vertices) p.x() += 60.f;
    its_merge(multishell, second);
    const auto guided_multi = prepare_cut_and_stretch_guided(multishell, Transform3d::Identity(),
        Vec3d(20,5,15), Vec3d::UnitZ(), top_guides(u,true));
    failed(guided_multi);
    CHECK(guided_multi.message.find("Disconnected") != std::string::npos);
}
