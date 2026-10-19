// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include <catch2/catch.hpp>

#include "libslic3r/CutGuideSelection.hpp"
#include "libslic3r/CutAndStretch.hpp"
#include "libslic3r/CutPlaneSelection.hpp"
#include "libslic3r/MeshDiagnostics.hpp"
#include "libslic3r/PlanarProfile.hpp"
#include "libslic3r/Measure.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/TriangleSelector.hpp"

#include <Eigen/LU>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

using namespace Slic3r;

namespace {

using MeshPtr = std::shared_ptr<const TriangleMesh>;
using Patches = std::vector<CutGuidePatch>;
using PaintData = std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>;

MeshPtr mesh_from(indexed_triangle_set its)
{
    return std::make_shared<TriangleMesh>(std::move(its));
}

MeshPtr hinge(double angle_degrees)
{
    const double angle = angle_degrees * PI / 180.;
    indexed_triangle_set its;
    its.vertices = {Vec3f(0, 0, 0), Vec3f(1, 0, 0), Vec3f(0, 1, 0),
                    Vec3f(0, float(-std::cos(angle)), float(std::sin(angle)))};
    its.indices = {Vec3i(0, 1, 2), Vec3i(1, 0, 3)};
    return mesh_from(std::move(its));
}

MeshPtr strip(const std::vector<double> &segment_angles)
{
    indexed_triangle_set its;
    Vec3d position = Vec3d::Zero();
    its.vertices.emplace_back(position.cast<float>());
    its.vertices.emplace_back((position + Vec3d::UnitY()).cast<float>());
    for (size_t i = 0; i < segment_angles.size(); ++i) {
        const double angle = segment_angles[i] * PI / 180.;
        position += Vec3d(std::cos(angle), 0., std::sin(angle));
        its.vertices.emplace_back(position.cast<float>());
        its.vertices.emplace_back((position + Vec3d::UnitY()).cast<float>());
        const int a = int(2 * i);
        its.indices.emplace_back(a, a + 2, a + 1);
        its.indices.emplace_back(a + 2, a + 3, a + 1);
    }
    return mesh_from(std::move(its));
}

MeshPtr ring()
{
    indexed_triangle_set its;
    its.vertices = {{-2,-2,0}, {2,-2,0}, {2,2,0}, {-2,2,0},
                    {-1,-1,0}, {1,-1,0}, {1,1,0}, {-1,1,0}};
    its.indices = {{0,1,5}, {0,5,4}, {1,2,6}, {1,6,5},
                   {2,3,7}, {2,7,6}, {3,0,4}, {3,4,7}};
    return mesh_from(std::move(its));
}

MeshPtr islands(size_t count, bool coincident = false)
{
    indexed_triangle_set its;
    for (size_t i = 0; i < count; ++i) {
        const int a = int(its.vertices.size());
        const float x = coincident ? 0.f : float(3 * i);
        its.vertices.emplace_back(x, 0.f, 0.f);
        its.vertices.emplace_back(x + 1.f, 0.f, 0.f);
        its.vertices.emplace_back(x, 1.f, 0.f);
        its.indices.emplace_back(a, a + 1, a + 2);
    }
    return mesh_from(std::move(its));
}

std::vector<int> measured_planes(const TriangleMesh &mesh)
{
    // Use the real measurement partition, including its connected-component rule.
    // No feature extraction is needed for open meshes or for the adapter's mask.
    Measure::Measuring measuring(mesh.its);
    std::vector<int> planes(mesh.its.indices.size(), -1);
    for (int plane = 0; plane < measuring.get_num_of_planes(); ++plane) {
        for (int facet : measuring.get_plane_triangle_indices(plane)) {
            REQUIRE(facet >= 0);
            REQUIRE(size_t(facet) < planes.size());
            REQUIRE((planes[size_t(facet)] == -1 || planes[size_t(facet)] == plane));
            planes[size_t(facet)] = plane;
        }
    }
    REQUIRE(std::all_of(planes.begin(), planes.end(), [](int plane) { return plane >= 0; }));
    return planes;
}

std::vector<int> all_facets(const TriangleMesh &mesh)
{
    std::vector<int> facets(mesh.its.indices.size());
    std::iota(facets.begin(), facets.end(), 0);
    return facets;
}

Vec3d vertex_zero(const TriangleMesh &mesh, size_t facet)
{
    return mesh.its.vertices[size_t(mesh.its.indices[facet][0])].cast<double>();
}

Vec3d hit_sample(const TriangleMesh &mesh, size_t facet)
{
    const Vec3i &f = mesh.its.indices[facet];
    // Deliberately retain a non-vertex, double-precision click witness.
    return 0.23 * mesh.its.vertices[size_t(f[0])].cast<double>() +
           0.31 * mesh.its.vertices[size_t(f[1])].cast<double>() +
           0.46 * mesh.its.vertices[size_t(f[2])].cast<double>();
}

Vec3d precise_normal(const TriangleMesh &mesh, size_t facet)
{
    const Vec3i &f = mesh.its.indices[facet];
    const Vec3d a = mesh.its.vertices[size_t(f[0])].cast<double>();
    const Vec3d b = mesh.its.vertices[size_t(f[1])].cast<double>();
    const Vec3d c = mesh.its.vertices[size_t(f[2])].cast<double>();
    return (b - a).cross(c - a).normalized();
}

std::vector<int> native_mask(TriangleSelector &selector)
{
    std::vector<int> mask;
    for (const auto &triangle : selector.get_triangles()) {
        if (triangle.valid() && !triangle.is_split() && triangle.is_selected_by_seed_fill())
            mask.push_back(triangle.source_triangle);
    }
    std::sort(mask.begin(), mask.end());
    mask.erase(std::unique(mask.begin(), mask.end()), mask.end());
    return mask;
}

std::vector<int> native_region(const TriangleMesh &mesh, size_t seed, const Vec3d &hit,
                               float angle, const Transform3d &pose = Transform3d::Identity())
{
    TriangleSelector selector(mesh);
    selector.seed_fill_select_triangles(hit.cast<float>(), int(seed), pose, {}, angle, 0.f, true);
    return native_mask(selector);
}

std::vector<int> patch_mask(const Patches &patches)
{
    std::vector<int> mask;
    for (const auto &patch : patches)
        mask.insert(mask.end(), patch.triangles.begin(), patch.triangles.end());
    std::sort(mask.begin(), mask.end());
    REQUIRE(std::adjacent_find(mask.begin(), mask.end()) == mask.end());
    return mask;
}

void check_patches(const TriangleMesh &mesh, const std::vector<int> &planes,
                   const Patches &patches, const std::vector<int> &expected_mask)
{
    CHECK(patch_mask(patches) == expected_mask);
    for (const auto &patch : patches) {
        REQUIRE_FALSE(patch.triangles.empty());
        CHECK(std::is_sorted(patch.triangles.begin(), patch.triangles.end()));
        REQUIRE(patch.seed_facet < mesh.its.indices.size());
        CHECK(std::binary_search(patch.triangles.begin(), patch.triangles.end(), int(patch.seed_facet)));
        CHECK(patch.hit_local.allFinite());
        CHECK(patch.normal_local.allFinite());
        CHECK((patch.normal_local - precise_normal(mesh, patch.seed_facet)).norm() < 1e-12);
        for (int facet : patch.triangles) {
            REQUIRE(facet >= 0);
            REQUIRE(size_t(facet) < planes.size());
            CHECK(patch.plane_id == planes[size_t(facet)]);
        }
    }
}

void check_identical_patches(const Patches &a, const Patches &b)
{
    REQUIRE(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        CHECK(a[i].plane_id == b[i].plane_id);
        CHECK(a[i].triangles == b[i].triangles);
        CHECK(a[i].seed_facet == b[i].seed_facet);
        CHECK(a[i].hit_local == b[i].hit_local);
        CHECK(a[i].normal_local == b[i].normal_local);
    }
}

void check_unchanged_mesh(const indexed_triangle_set &before, const indexed_triangle_set &after)
{
    CHECK(before.vertices == after.vertices);
    CHECK(before.indices == after.indices);
    REQUIRE(before.properties.size() == after.properties.size());
    for (size_t i = 0; i < before.properties.size(); ++i) {
        CHECK(before.properties[i].type == after.properties[i].type);
        CHECK(before.properties[i].area == after.properties[i].area);
    }
}

} // namespace

TEST_CASE("Cut guide angular regions enumerate actual native selection flags without neighboring paint-state expansion", "[cut_guide_selection]")
{
    const MeshPtr mesh = hinge(30.);
    const auto before = mesh->its;
    const auto planes = measured_planes(*mesh);
    CutGuideSelection selection(mesh, planes);
    const Vec3d hit = hit_sample(*mesh, 0);
    const auto region = selection.angular_region(0, hit, 10.f);
    REQUIRE(region);
    CHECK(*region == std::vector<int>{0});
    CHECK(*region == native_region(*mesh, 0, hit, 10.f));

    TriangleSelector native(*mesh);
    const auto paint_before = native.serialize();
    native.seed_fill_select_triangles(hit.cast<float>(), 0, Transform3d::Identity(), {}, 10.f, 0.f, true);
    CHECK(native_mask(native) == *region);
    CHECK(native.serialize() == paint_before);
    // This public mesh API is intentionally NOT the mask oracle: same-state neighbors leak into it.
    int state = -1;
    CHECK(native.get_seed_fill_mesh(state).indices.size() > region->size());
    const Patches patches = selection.toggle({}, *region, 0, hit);
    REQUIRE(patches.size() == 1);
    CHECK(patches.front().hit_local == hit);
    check_patches(*mesh, planes, patches, {0});
    check_unchanged_mesh(before, mesh->its);
}

TEST_CASE("Cut guide angular thresholds retain native epsilon clamp and float boundary behavior", "[cut_guide_selection]")
{
    for (double bend : {0.5, 1.0}) {
        const MeshPtr mesh = hinge(bend);
        CutGuideSelection selection(mesh, measured_planes(*mesh));
        const auto region = selection.angular_region(0, hit_sample(*mesh, 0), 0.f);
        REQUIRE(region);
        CHECK(*region == native_region(*mesh, 0, hit_sample(*mesh, 0), 0.f));
        CHECK(region->size() == (bend == 0.5 ? size_t(2) : size_t(1)));
    }

    const MeshPtr mesh = hinge(30.);
    CutGuideSelection selection(mesh, measured_planes(*mesh));
    const auto normals = its_face_normals(mesh->its);
    const double dot = std::clamp(double(normals[0].dot(normals[1])), 0., 1.);
    const float boundary = float(std::acos(std::min(1., dot + EPSILON)) * 180. / PI);
    const std::array<float, 5> angles{boundary - 0.001f, std::nextafter(boundary, 0.f), boundary,
                                    std::nextafter(boundary, 90.f), boundary + 0.001f};
    for (float angle : angles) {
        INFO("Angle near native float threshold: " << angle);
        const auto region = selection.angular_region(0, hit_sample(*mesh, 0), angle);
        REQUIRE(region);
        CHECK(*region == native_region(*mesh, 0, hit_sample(*mesh, 0), angle));
    }
    CHECK(selection.angular_region(0, hit_sample(*mesh, 0), angles.front())->size() == 1);
    CHECK(selection.angular_region(0, hit_sample(*mesh, 0), angles.back())->size() == 2);

    const MeshPtr opposing = hinge(120.);
    CutGuideSelection opposite_selection(opposing, measured_planes(*opposing));
    CHECK(its_face_normals(opposing->its)[0].dot(its_face_normals(opposing->its)[1]) < 0.f);
    CHECK(opposite_selection.angular_region(0, hit_sample(*opposing, 0), 89.99f)->size() == 1);
    for (float angle : {89.999f, 90.f}) {
        const auto region = opposite_selection.angular_region(0, hit_sample(*opposing, 0), angle);
        CHECK(*region == native_region(*opposing, 0, hit_sample(*opposing, 0), angle));
        CHECK(*region == std::vector<int>{0, 1});
    }
}

TEST_CASE("Cut guide fill follows local adjacent angles through a smooth chain under mirrored nonuniform poses", "[cut_guide_selection]")
{
    const MeshPtr mesh = strip({0., 20., 40., 60.});
    const auto before = mesh->its;
    CutGuideSelection selection(mesh, measured_planes(*mesh));
    const Vec3d hit = hit_sample(*mesh, 0);
    const auto region = selection.angular_region(0, hit, 21.f);
    REQUIRE(region);
    CHECK(*region == all_facets(*mesh));
    CHECK(*selection.angular_region(0, hit, 19.f) == std::vector<int>{0, 1});
    CHECK(precise_normal(*mesh, 0).dot(precise_normal(*mesh, 6)) < std::cos(21. * PI / 180.));

    for (double mirror : {1., -1.}) {
        Transform3d pose = Transform3d::Identity();
        pose.linear() = Eigen::AngleAxisd(0.43, Vec3d(1., 2., 3.).normalized()).toRotationMatrix() *
                        Vec3d(mirror, 1., 4.).asDiagonal();
        const Vec3d first_world = (pose.linear().inverse().transpose() * precise_normal(*mesh, 0)).normalized();
        const Vec3d next_world = (pose.linear().inverse().transpose() * precise_normal(*mesh, 2)).normalized();
        CHECK(first_world.dot(next_world) < std::cos(21. * PI / 180.));
        // Native seed-fill's pose is used by overhang filtering, not its adjacent-angle metric.
        CHECK(native_region(*mesh, 0, hit, 21.f, pose) == *region);
        CHECK(*selection.angular_region(0, hit, 21.f) == *region);
    }
    check_unchanged_mesh(before, mesh->its);
}

TEST_CASE("Cut guide connectivity uses native indexed edges rather than coincident positions or same-directed edges", "[cut_guide_selection]")
{
    const MeshPtr disconnected = islands(2, true);
    const auto planes = measured_planes(*disconnected);
    REQUIRE(planes[0] != planes[1]);
    CutGuideSelection selection(disconnected, planes);
    for (size_t seed : {size_t(0), size_t(1)}) {
        const auto region = selection.angular_region(seed, hit_sample(*disconnected, seed), 90.f);
        CHECK(*region == std::vector<int>{int(seed)});
        CHECK(*region == native_region(*disconnected, seed, hit_sample(*disconnected, seed), 90.f));
    }

    indexed_triangle_set its;
    its.vertices = {{0,0,0}, {1,0,0}, {0,1,0}, {0,-1,0}};
    its.indices = {{0,1,2}, {0,1,3}}; // The shared edge has the same direction: not a native neighbor.
    const MeshPtr same_direction = mesh_from(std::move(its));
    CutGuideSelection separated(same_direction, measured_planes(*same_direction));
    CHECK(*separated.angular_region(0, hit_sample(*same_direction, 0), 90.f) == std::vector<int>{0});
}

TEST_CASE("Cut guide region cache shares immutable masks only for the same mesh angle and connected region", "[cut_guide_selection]")
{
    const MeshPtr mesh = hinge(30.);
    const auto planes = measured_planes(*mesh);
    CutGuideSelection selection(mesh, planes);
    const auto first = selection.angular_region(0, hit_sample(*mesh, 0), 40.f);
    REQUIRE(first);
    const auto first_snapshot = *first;
    const auto inside = selection.angular_region(1, hit_sample(*mesh, 1), 40.f);
    CHECK(inside.get() == first.get());
    const Vec3d rounded_hit = hit_sample(*mesh, 1).cast<float>().cast<double>();
    CHECK(selection.angular_region(1, rounded_hit, 40.f).get() == first.get());

    const auto changed_angle = selection.angular_region(1, rounded_hit, 45.f);
    CHECK(changed_angle.get() != first.get());
    CHECK(*changed_angle == *first);
    const auto narrowed = selection.angular_region(1, rounded_hit, 10.f);
    CHECK(*narrowed == std::vector<int>{1});
    CHECK(*first == first_snapshot);

    const MeshPtr same_geometry = std::make_shared<TriangleMesh>(*mesh);
    CutGuideSelection other_mesh(same_geometry, planes);
    const auto other_region = other_mesh.angular_region(0, hit_sample(*same_geometry, 0), 40.f);
    CHECK(other_region.get() != first.get());
    CHECK(*other_region == first_snapshot);

    const MeshPtr separate_mesh = islands(2);
    CutGuideSelection separate(separate_mesh, measured_planes(*separate_mesh));
    const auto left = separate.angular_region(0, hit_sample(*separate_mesh, 0), 30.f);
    const auto right = separate.angular_region(1, hit_sample(*separate_mesh, 1), 30.f);
    CHECK(left.get() != right.get());
    CHECK(*left == std::vector<int>{0});
    CHECK(*right == std::vector<int>{1});
}

TEST_CASE("Cut guide patches preserve exact partial measured-plane masks and split their disconnected components", "[cut_guide_selection]")
{
    const MeshPtr mesh = strip({0., 0., 0.});
    const auto planes = measured_planes(*mesh);
    REQUIRE(std::all_of(planes.begin(), planes.end(), [&planes](int id) { return id == planes.front(); }));
    CutGuideSelection selection(mesh, planes);
    const Vec3d clicked_hit = hit_sample(*mesh, 4);
    const Patches partial = selection.toggle({}, {0, 1, 4, 5}, 4, clicked_hit);
    REQUIRE(partial.size() == 2);
    check_patches(*mesh, planes, partial, {0, 1, 4, 5});
    const auto clicked = std::find_if(partial.begin(), partial.end(), [](const CutGuidePatch &patch) { return patch.seed_facet == 4; });
    const auto fallback = std::find_if(partial.begin(), partial.end(), [](const CutGuidePatch &patch) { return patch.seed_facet == 0; });
    REQUIRE(clicked != partial.end());
    REQUIRE(fallback != partial.end());
    CHECK(clicked->triangles == std::vector<int>{4, 5});
    CHECK(clicked->hit_local == clicked_hit);
    CHECK(fallback->triangles == std::vector<int>{0, 1});
    CHECK(fallback->hit_local == vertex_zero(*mesh, 0));
    check_identical_patches(partial, selection.toggle({}, {5, 1, 4, 0}, 4, clicked_hit));

    const MeshPtr cube = std::make_shared<TriangleMesh>(its_make_cube(2., 3., 4.));
    const auto cube_planes = measured_planes(*cube);
    CutGuideSelection cube_selection(cube, cube_planes);
    const auto region = cube_selection.angular_region(0, hit_sample(*cube, 0), 90.f);
    const Patches faces = cube_selection.toggle({}, *region, 0, hit_sample(*cube, 0));
    REQUIRE(faces.size() == 6);
    check_patches(*cube, cube_planes, faces, all_facets(*cube));
    CHECK(faces.front().seed_facet == 0);
}

TEST_CASE("Cut guide atomic toggles union partial overlaps and subtract only an already fully selected batch", "[cut_guide_selection]")
{
    const MeshPtr mesh = strip({0., 0., 0.});
    const auto planes = measured_planes(*mesh);
    CutGuideSelection selection(mesh, planes);
    const Vec3d original_hit = hit_sample(*mesh, 1);
    const Patches current = selection.toggle({}, {0, 1, 2, 3}, 1, original_hit);
    const Patches current_before = current;
    const Patches united = selection.toggle(current, {2, 3, 4, 5}, 4, hit_sample(*mesh, 4));
    REQUIRE(united.size() == 1);
    check_patches(*mesh, planes, united, {0, 1, 2, 3, 4, 5});
    CHECK(united.front().seed_facet == 1);
    CHECK(united.front().hit_local == original_hit);
    check_identical_patches(current, current_before);

    const Patches subtracted = selection.toggle(united, {2, 3, 4, 5}, 4, hit_sample(*mesh, 4));
    REQUIRE(subtracted.size() == 1);
    check_patches(*mesh, planes, subtracted, {0, 1});
    CHECK(subtracted.front().seed_facet == 1);
    CHECK(subtracted.front().hit_local == original_hit);
    CHECK(selection.toggle(subtracted, {0, 1}, 0, hit_sample(*mesh, 0)).empty());

    const Patches split = selection.toggle(united, {2, 3}, 2, hit_sample(*mesh, 2));
    REQUIRE(split.size() == 2);
    check_patches(*mesh, planes, split, {0, 1, 4, 5});
    CHECK(split.front().seed_facet == 1);
    CHECK(split.front().hit_local == original_hit);
    CHECK(split.back().seed_facet == 4);
    CHECK(split.back().hit_local == vertex_zero(*mesh, 4));
}

TEST_CASE("Cut guide removed witnesses are replaced by a real source vertex rather than a patch or hole centroid", "[cut_guide_selection]")
{
    const MeshPtr mesh = ring();
    const auto planes = measured_planes(*mesh);
    REQUIRE(std::all_of(planes.begin(), planes.end(), [&planes](int id) { return id == planes.front(); }));
    CutGuideSelection selection(mesh, planes);
    const Patches current = selection.toggle({}, all_facets(*mesh), 0, hit_sample(*mesh, 0));
    REQUIRE(current.size() == 1);
    const Patches remaining = selection.toggle(current, {0}, 0, hit_sample(*mesh, 0));
    REQUIRE(remaining.size() == 1);
    check_patches(*mesh, planes, remaining, {1, 2, 3, 4, 5, 6, 7});
    CHECK(remaining.front().seed_facet == 1);
    CHECK(remaining.front().hit_local == vertex_zero(*mesh, 1));
    CHECK((std::abs(remaining.front().hit_local.x()) >= 1. || std::abs(remaining.front().hit_local.y()) >= 1.));
    const Patches repeated = selection.toggle(current, {0}, 0, hit_sample(*mesh, 0));
    check_identical_patches(remaining, repeated);
}

TEST_CASE("Cut guide normals promote source vertices before subtraction and cross product and tolerate legacy witness rounding", "[cut_guide_selection]")
{
    indexed_triangle_set its;
    its.vertices = {{-2.1253f, 1.7539f, 3.9781f}, {6.0347f, -4.4519f, 2.3331f},
                    {0.3337f, 8.8273f, -1.1669f}};
    its.indices = {{0, 1, 2}};
    const MeshPtr oblique = mesh_from(std::move(its));
    const Vec3f &a = oblique->its.vertices[0], &b = oblique->its.vertices[1], &c = oblique->its.vertices[2];
    const Vec3d float_cross_promoted = (b - a).cross(c - a).cast<double>().normalized();
    const Vec3d expected_normal = precise_normal(*oblique, 0);
    REQUIRE((float_cross_promoted - expected_normal).norm() > 1e-10);
    CutGuideSelection precise(oblique, measured_planes(*oblique));
    const Vec3d original_hit = hit_sample(*oblique, 0).cast<float>().cast<double>();
    const Patches selected = precise.toggle({}, {0}, 0, original_hit);
    REQUIRE(selected.size() == 1);
    CHECK(selected.front().hit_local == original_hit);
    CHECK((selected.front().normal_local - expected_normal).norm() < 1e-12);

    const MeshPtr planar = strip({0.});
    const auto planes = measured_planes(*planar);
    CutGuideSelection selection(planar, planes);
    Patches legacy = selection.toggle({}, {0}, 0, hit_sample(*planar, 0));
    legacy.front().normal_local.z() = std::nextafter(legacy.front().normal_local.z(), 2.);
    const Vec3d legacy_hit = legacy.front().hit_local;
    const Patches joined = selection.toggle(legacy, {1}, 1, hit_sample(*planar, 1));
    REQUIRE(joined.size() == 1);
    check_patches(*planar, planes, joined, {0, 1});
    CHECK(joined.front().seed_facet == 0);
    CHECK(joined.front().hit_local == legacy_hit);
    CHECK(joined.front().normal_local == precise_normal(*planar, 0));
}

TEST_CASE("Cut guide source preflight rejects bounded invalid geometry before native selector allocation", "[cut_guide_selection]")
{
    const MeshPtr valid = hinge(30.);
    const auto planes = measured_planes(*valid);
    REQUIRE_NOTHROW(CutGuideSelection::validate_source(*valid));
    CHECK_THROWS(CutGuideSelection(MeshPtr{}, {}));
    CHECK_THROWS(CutGuideSelection(valid, {}));
    auto bad_planes = planes;
    bad_planes[0] = -1;
    CHECK_THROWS(CutGuideSelection(valid, bad_planes));
    bad_planes[0] = int(valid->its.indices.size());
    CHECK_THROWS(CutGuideSelection(valid, bad_planes));

    CutGuideSelectionLimits small_limit;
    small_limit.max_source_facets = 1;
    CHECK_THROWS(CutGuideSelection::validate_source(*valid, small_limit));
    CHECK_THROWS(CutGuideSelection(valid, planes, small_limit));
    TriangleMesh unused_vertices;
    unused_vertices.its.vertices = {{0,0,0}, {1,0,0}, {0,1,0}, {7,8,9}};
    unused_vertices.its.indices = {{0,1,2}};
    CHECK_THROWS(CutGuideSelection::validate_source(unused_vertices, small_limit));
    small_limit.max_source_facets = 2;
    REQUIRE_NOTHROW(CutGuideSelection::validate_source(unused_vertices, small_limit));
    small_limit.max_source_facets = 1;
    CutGuideSelectionLimits zero_limit;
    zero_limit.max_source_facets = 0;
    CHECK_THROWS(CutGuideSelection::validate_source(*valid, zero_limit));
    zero_limit = {};
    zero_limit.max_guides = 0;
    CHECK_THROWS(CutGuideSelection::validate_source(*valid, zero_limit));

    TriangleMesh empty;
    CHECK_THROWS(CutGuideSelection::validate_source(empty));
    TriangleMesh invalid_index;
    invalid_index.its = valid->its;
    invalid_index.its.indices[0][0] = std::numeric_limits<int>::max();
    CHECK_THROWS(CutGuideSelection::validate_source(invalid_index));
    CHECK_THROWS(CutGuideSelection::validate_source(invalid_index, small_limit));

    TriangleMesh bad_vertex;
    bad_vertex.its = valid->its;
    bad_vertex.its.vertices[0].x() = std::numeric_limits<float>::quiet_NaN();
    CHECK_THROWS(CutGuideSelection::validate_source(bad_vertex));
    bad_vertex.its.vertices[0].x() = std::numeric_limits<float>::infinity();
    CHECK_THROWS(CutGuideSelection::validate_source(bad_vertex));

    TriangleMesh degenerate;
    degenerate.its.vertices = {{0,0,0}, {1,0,0}, {2,0,0}};
    degenerate.its.indices = {{0,1,2}};
    CHECK_THROWS(CutGuideSelection::validate_source(degenerate));
    for (float extent : {1e-25f, 1e15f, 1e30f}) {
        TriangleMesh unstable;
        unstable.its.vertices = {{0,0,0}, {extent,0,0}, {0,extent,0}};
        unstable.its.indices = {{0,1,2}};
        const Vec3d cross = unstable.its.vertices[1].cast<double>().cross(unstable.its.vertices[2].cast<double>());
        REQUIRE(cross.allFinite());
        REQUIRE(cross.norm() > 0.);
        // A positive finite double cross is insufficient: native float normals under/overflow.
        CHECK_THROWS(CutGuideSelection::validate_source(unstable));
    }
}

TEST_CASE("Cut guide invalid queries and guide-budget failures leave existing selections masks and mesh unchanged", "[cut_guide_selection]")
{
    const MeshPtr mesh = strip({0., 0., 0.});
    const auto before = mesh->its;
    const auto planes = measured_planes(*mesh);
    CutGuideSelectionLimits limits;
    limits.max_guides = 1;
    CutGuideSelection selection(mesh, planes, limits);
    const Vec3d hit = hit_sample(*mesh, 0);
    const auto cached = selection.angular_region(0, hit, 30.f);
    const auto cached_before = *cached;
    for (float angle : {-1.f, 91.f, std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN()})
        CHECK_THROWS(selection.angular_region(0, hit, angle));
    CHECK_THROWS(selection.angular_region(mesh->its.indices.size(), hit, 30.f));
    CHECK_THROWS(selection.angular_region(std::numeric_limits<size_t>::max(), hit, 30.f));
    for (double x : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                     2. * double(std::numeric_limits<float>::max())}) {
        const Vec3d invalid_hit(x, 0., 0.);
        CHECK_THROWS(selection.angular_region(0, invalid_hit, 30.f));
        CHECK_THROWS(selection.toggle({}, {0}, 0, invalid_hit));
    }
    CHECK(*cached == cached_before);
    CHECK(selection.angular_region(0, hit, 30.f).get() == cached.get());

    const Patches current = selection.toggle({}, {0, 1}, 0, hit);
    const Patches current_before = current;
    CHECK_THROWS(selection.toggle(current, {4, 5}, 4, hit_sample(*mesh, 4)));
    check_identical_patches(current, current_before);
    CHECK_THROWS(selection.toggle(current, {}, 0, hit));
    CHECK_THROWS(selection.toggle(current, {0, 0}, 0, hit));
    CHECK_THROWS(selection.toggle(current, {1}, 0, hit));
    CHECK_THROWS(selection.toggle({current.front(), current.front()}, {0}, 0, hit));
    CHECK_THROWS(selection.toggle(current, {-1}, 0, hit));
    CHECK_THROWS(selection.toggle(current, {int(mesh->its.indices.size())}, 0, hit));
    CHECK_THROWS(selection.toggle(current, {0}, mesh->its.indices.size(), hit));
    Patches invalid_current = current;
    invalid_current.front().seed_facet = 4;
    CHECK_THROWS(selection.toggle(invalid_current, {2}, 2, hit_sample(*mesh, 2)));
    check_identical_patches(current, current_before);
    check_unchanged_mesh(before, mesh->its);

    // Exercise the real default budget without allocating an unbounded source mesh.
    const MeshPtr many = islands(4097);
    const auto many_before = many->its;
    CutGuideSelection default_limit(many, measured_planes(*many));
    const Patches one = default_limit.toggle({}, {0}, 0, hit_sample(*many, 0));
    const Patches one_before = one;
    CHECK_THROWS(default_limit.toggle(one, all_facets(*many), 1, hit_sample(*many, 1)));
    check_identical_patches(one, one_before);
    check_unchanged_mesh(many_before, many->its);
}

TEST_CASE("Cut guide selection never reads or writes painted Model annotations IDs mesh or instance and volume poses", "[cut_guide_selection]")
{
    Model model;
    ModelObject *object = model.add_object();
    ModelVolume *volume = object->add_volume(TriangleMesh(its_make_cube(2., 3., 4.)), ModelVolumeType::MODEL_PART, false);
    ModelInstance *instance = object->add_instance();
    Geometry::Transformation volume_pose;
    volume_pose.set_offset(Vec3d(3., -2., 4.));
    volume_pose.set_rotation(Vec3d(0.2, -0.1, 0.3));
    volume_pose.set_scaling_factor(Vec3d(1.3, 0.7, 2.));
    volume_pose.set_mirror(Vec3d(-1., 1., 1.));
    volume->set_transformation(volume_pose);
    instance->set_offset(Vec3d(450., 100., 7.));
    instance->set_rotation(Vec3d(-0.3, 0.1, 0.7));

    TriangleSelector painted(volume->mesh());
    painted.set_facet(0, EnforcerBlockerType::ENFORCER);
    painted.set_facet(1, EnforcerBlockerType::BLOCKER);
    std::array<FacetsAnnotation *, 5> layers{&volume->supported_facets, &volume->seam_facets,
        &volume->mmu_segmentation_facets, &volume->fuzzy_skin_facets, &volume->exterior_facets};
    std::array<PaintData, 5> annotation_before;
    std::array<ObjectID, 5> annotation_ids;
    std::array<ObjectBase::Timestamp, 5> annotation_timestamps;
    for (size_t i = 0; i < layers.size(); ++i) {
        REQUIRE(layers[i]->set(painted));
        REQUIRE_FALSE(layers[i]->empty());
        annotation_before[i] = layers[i]->get_data();
        annotation_ids[i] = layers[i]->id();
        annotation_timestamps[i] = layers[i]->timestamp();
    }
    const auto mesh = volume->get_mesh_shared_ptr();
    const auto geometry_before = mesh->its;
    const Transform3d volume_before = volume->get_matrix(), instance_before = instance->get_matrix();
    const ObjectID object_id = object->id(), volume_id = volume->id(), instance_id = instance->id();
    const auto planes = measured_planes(*mesh);
    CutGuideSelection selection(mesh, planes);
    const Vec3d hit = hit_sample(*mesh, 0);
    const auto region = selection.angular_region(0, hit, 90.f);
    CHECK(*region == native_region(*mesh, 0, hit, 90.f));
    const Patches selected = selection.toggle({}, *region, 0, hit);
    REQUIRE(selected.size() == 6);
    CHECK(selection.toggle(selected, *region, 0, hit).empty());

    CHECK(object->id() == object_id);
    CHECK(volume->id() == volume_id);
    CHECK(instance->id() == instance_id);
    CHECK(volume->get_mesh_shared_ptr() == mesh);
    CHECK((volume->get_matrix().matrix().array() == volume_before.matrix().array()).all());
    CHECK((instance->get_matrix().matrix().array() == instance_before.matrix().array()).all());
    check_unchanged_mesh(geometry_before, volume->mesh().its);
    for (size_t i = 0; i < layers.size(); ++i) {
        CHECK(layers[i]->get_data() == annotation_before[i]);
        CHECK(layers[i]->id() == annotation_ids[i]);
        CHECK(layers[i]->timestamp() == annotation_timestamps[i]);
    }
    REQUIRE(model.objects.size() == 1);
    REQUIRE(object->volumes.size() == 1);
    REQUIRE(object->instances.size() == 1);
}

TEST_CASE("Smart L guide provenance reaches the real side classifier and placement-aware cut consumers unchanged", "[cut_guide_selection]")
{
    const auto generated = extrude_planar_profile(
        {Vec3d(0,0,0), Vec3d(10,0,0), Vec3d(10,0,20),
         Vec3d(40,0,20), Vec3d(40,0,24), Vec3d(0,0,24)}, Vec3d(0,10,0));
    INFO(generated.message);
    REQUIRE(generated);
    const MeshPtr mesh = mesh_from(generated.mesh);
    const auto source_before = mesh->its;
    const auto planes = measured_planes(*mesh);
    size_t seed = mesh->its.indices.size();
    for (size_t i = 0; i < mesh->its.indices.size(); ++i) {
        const Vec3i &f = mesh->its.indices[i];
        if (mesh->its.vertices[size_t(f[0])].z() == 20.f &&
            mesh->its.vertices[size_t(f[1])].z() == 20.f &&
            mesh->its.vertices[size_t(f[2])].z() == 20.f && precise_normal(*mesh, i).z() < 0.) {
            seed = i;
            break;
        }
    }
    REQUIRE(seed < mesh->its.indices.size());
    const Vec3d original_hit = hit_sample(*mesh, seed);
    std::vector<int> manual_plane_mask;
    for (size_t i = 0; i < planes.size(); ++i)
        if (planes[i] == planes[seed]) manual_plane_mask.push_back(int(i));
    CutGuideSelection selection(mesh, planes);

    for (float angle : {0.f, 30.f}) {
        const auto region = selection.angular_region(seed, original_hit, angle);
        REQUIRE(*region == manual_plane_mask);
        const Patches patches = selection.toggle({}, *region, seed, original_hit);
        REQUIRE(patches.size() == 1);
        CHECK(patches.front().seed_facet == seed);
        CHECK(patches.front().hit_local == original_hit);
        const Patches patches_before = patches;

        for (size_t placement = 0; placement < 3; ++placement) {
            INFO("angle=" << angle << ", placement=" << placement);
            const Vec3d scale = placement == 0 ? Vec3d::Ones() : Vec3d(placement == 2 ? -1.25 : 1.25, 0.8, 1.5);
            Transform3d pose = Transform3d::Identity();
            pose.linear() = Eigen::AngleAxisd(0.17 * double(placement), Vec3d::UnitZ()).toRotationMatrix() * scale.asDiagonal();
            pose.translation() = placement == 0 ? Vec3d::Zero() : Vec3d(450., 100., 7.);
            const Transform3d pose_before = pose;
            const Matrix3d normal_matrix = pose.linear().inverse().transpose();
            const auto world_guide = [&](const Vec3d &hit, const Vec3d &normal) {
                return CutPlaneFaceGuide{pose * hit, (normal_matrix * normal).normalized()};
            };
            std::vector<CutPlaneFaceGuide> smart_guides;
            for (const auto &patch : patches)
                smart_guides.push_back(world_guide(patch.hit_local, patch.normal_local));
            const CutPlaneFaceGuide manual = world_guide(original_hit, precise_normal(*mesh, seed));
            const auto &chosen = smart_guides.front();
            const auto smart_side = select_cut_moving_side(smart_guides, chosen.point_world, chosen.outward_normal_world);
            const auto manual_side = select_cut_moving_side({manual}, manual.point_world, manual.outward_normal_world);
            INFO(smart_side.message);
            REQUIRE(smart_side);
            REQUIRE(manual_side);
            CHECK(smart_side.side == -1);
            CHECK(smart_side.side == manual_side.side);
            CHECK(smart_side.moving_normal.isApprox(manual_side.moving_normal, 1e-12));
            CHECK(smart_side.moving_normal.isApprox(Vec3d::UnitZ(), 1e-12));

            const auto plan = prepare_cut_and_stretch_from_face(mesh->its, pose, patches.front().seed_facet, smart_side.moving_normal);
            const auto manual_plan = prepare_cut_and_stretch_from_face(mesh->its, pose, seed, manual_side.moving_normal);
            INFO(plan.message);
            REQUIRE(plan);
            REQUIRE(manual_plan);
            REQUIRE(plan.source_frame);
            CHECK(plan.contact_area_mm2 == Approx(100. * std::abs(scale.x() * scale.y())).epsilon(2e-5));
            CHECK(plan.source_volume_mm3 == Approx(3600. * std::abs(pose.linear().determinant())).epsilon(1e-10));
            CHECK(plan.contact_area_mm2 == Approx(manual_plan.contact_area_mm2).margin(1e-8));
            CHECK(plan.cut_to_world.matrix().isApprox(manual_plan.cut_to_world.matrix(), 1e-12));
            CHECK(std::abs(smart_side.moving_normal.dot(plan.cut_to_world.translation() - chosen.point_world)) < 1e-9);
            CHECK(plan.cut_to_world.linear().col(2).isApprox(smart_side.moving_normal, 1e-12));

            // Exercise the other real GUI consumer: an independent interior plane,
            // not an approximate world-space reconstruction of the source face.
            const Vec3d interior_origin = pose * Vec3d(5., 5., 12.);
            const Vec3d interior_normal = (normal_matrix * Vec3d::UnitZ()).normalized();
            const auto interior_side = select_cut_moving_side(smart_guides, interior_origin, interior_normal);
            REQUIRE(interior_side);
            CHECK(interior_side.side == 1);
            const auto interior_plan = prepare_cut_and_stretch(mesh->its, pose, interior_origin, interior_side.moving_normal);
            INFO(interior_plan.message);
            REQUIRE(interior_plan);
            REQUIRE(interior_plan.source_frame);
            CHECK(interior_plan.contact_area_mm2 == Approx(plan.contact_area_mm2).epsilon(2e-5));
            CHECK(std::abs(interior_normal.dot(interior_plan.cut_to_world.translation() - interior_origin)) < 1e-9);

            const auto stationary_before = plan.stationary, moving_before = plan.moving, bridge_before = plan.bridge_unit;
            constexpr double distance = 6.;
            const auto result = apply_cut_and_stretch_to_source(plan, distance);
            INFO(result.message);
            REQUIRE(result);
            const auto diagnostics = its_mesh_diagnostics(result.mesh);
            CHECK(diagnostics.open_edges == 0);
            CHECK(diagnostics.non_manifold_edges == 0);
            CHECK(diagnostics.non_manifold_vertices == 0);
            CHECK(diagnostics.same_direction_edges == 0);
            CHECK_FALSE(diagnostics.has_reversed_faces);
            for (const Vec3d &corner : {Vec3d(40,0,20), Vec3d(40,10,20), Vec3d(40,0,24), Vec3d(40,10,24)}) {
                const Vec3d expected_local = corner + Vec3d(0., 0., distance / scale.z());
                CHECK(std::any_of(result.mesh.vertices.begin(), result.mesh.vertices.end(), [&](const Vec3f &v) {
                    return (v.cast<double>() - expected_local).norm() < 1e-4;
                }));
                const Vec3d expected_world = pose * corner + distance * smart_side.moving_normal;
                CHECK(std::any_of(plan.moving.vertices.begin(), plan.moving.vertices.end(), [&](const Vec3f &v) {
                    return (plan.cut_to_world * (v.cast<double>() + Vec3d(0., 0., distance)) - expected_world).norm() < 1e-4;
                }));
            }
            check_unchanged_mesh(stationary_before, plan.stationary);
            check_unchanged_mesh(moving_before, plan.moving);
            check_unchanged_mesh(bridge_before, plan.bridge_unit);
            check_unchanged_mesh(source_before, mesh->its);
            check_identical_patches(patches, patches_before);
            CHECK((pose.matrix().array() == pose_before.matrix().array()).all());
        }
    }
}

TEST_CASE("Whole-cube smart flood is rejected by the actual side consumer when its derived guides contradict", "[cut_guide_selection]")
{
    Model model;
    ModelObject *object = model.add_object();
    ModelVolume *volume = object->add_volume(TriangleMesh(its_make_cube(2., 3., 4.)), ModelVolumeType::MODEL_PART, false);
    ModelInstance *instance = object->add_instance();
    volume->set_offset(Vec3d(3., -2., 4.));
    instance->set_offset(Vec3d(450., 100., 7.));
    instance->set_rotation(Vec3d(0., 0., 0.37));
    instance->set_scaling_factor(Vec3d(1.25, 0.8, 1.5));
    instance->set_mirror(Vec3d(-1., 1., 1.));
    TriangleSelector painted(volume->mesh());
    painted.set_facet(0, EnforcerBlockerType::ENFORCER);
    REQUIRE(volume->supported_facets.set(painted));
    const auto paint_before = volume->supported_facets.get_data();
    const auto paint_id = volume->supported_facets.id();
    const auto paint_timestamp = volume->supported_facets.timestamp();
    const std::array<ObjectID, 3> ids_before{object->id(), volume->id(), instance->id()};
    const Transform3d instance_before = instance->get_matrix(), volume_before = volume->get_matrix();
    const Transform3d pose = instance_before * volume_before;
    const auto mesh = volume->get_mesh_shared_ptr();
    const auto source_before = mesh->its;
    const auto planes = measured_planes(*mesh);
    CutGuideSelection selection(mesh, planes);
    size_t seed = mesh->its.indices.size();
    for (size_t i = 0; i < mesh->its.indices.size(); ++i)
        if (precise_normal(*mesh, i).z() > 0.99) { seed = i; break; }
    REQUIRE(seed < mesh->its.indices.size());
    const Vec3d hit = hit_sample(*mesh, seed);
    const auto region = selection.angular_region(seed, hit, 90.f);
    CHECK(*region == all_facets(*mesh));
    const Patches patches = selection.toggle({}, *region, seed, hit);
    REQUIRE(patches.size() == 6);
    const Patches patches_before = patches;
    std::vector<CutPlaneFaceGuide> guides;
    const Matrix3d normal_matrix = pose.linear().inverse().transpose();
    for (const auto &patch : patches)
        guides.push_back({pose * patch.hit_local, (normal_matrix * patch.normal_local).normalized()});
    const auto guides_before = guides;
    const Vec3d plane_origin = pose * Vec3d(1., 1.5, 1.25);
    const Vec3d plane_normal = (normal_matrix * Vec3d::UnitZ()).normalized();
    const auto signed_distance = [&](const CutPlaneFaceGuide &guide) {
        return plane_normal.dot(guide.point_world - plane_origin);
    };
    const auto above = std::find_if(guides.begin(), guides.end(), [&](const auto &guide) { return signed_distance(guide) > 1.; });
    const auto below = std::find_if(guides.begin(), guides.end(), [&](const auto &guide) { return signed_distance(guide) < -1.; });
    REQUIRE(above != guides.end());
    REQUIRE(below != guides.end());
    const auto require_contradiction = [&](const std::vector<CutPlaneFaceGuide> &input, const Vec3d &normal) {
        const auto side = select_cut_moving_side(input, plane_origin, normal);
        CHECK_FALSE(side);
        CHECK(side.side == 0);
        CHECK(side.moving_normal == Vec3d::Zero());
        CHECK(side.message.find("opposite sides") != std::string::npos);
        CHECK(side.message.find("no majority") != std::string::npos);
    };
    require_contradiction(guides, plane_normal);
    require_contradiction(guides, -plane_normal);
    std::vector<CutPlaneFaceGuide> reversed(guides.rbegin(), guides.rend());
    require_contradiction(reversed, plane_normal);
    // Every witness still comes from a real smart-derived patch; repeating one
    // side must not turn contradictory evidence into an accepted majority.
    std::vector<CutPlaneFaceGuide> majority(12, *above);
    majority.push_back(*below);
    require_contradiction(majority, plane_normal);

    check_identical_patches(patches, patches_before);
    REQUIRE(guides.size() == guides_before.size());
    for (size_t i = 0; i < guides.size(); ++i) {
        CHECK(guides[i].point_world == guides_before[i].point_world);
        CHECK(guides[i].outward_normal_world == guides_before[i].outward_normal_world);
    }
    CHECK(selection.angular_region(seed, hit, 90.f).get() == region.get());
    check_unchanged_mesh(source_before, volume->mesh().its);
    CHECK(volume->get_mesh_shared_ptr() == mesh);
    CHECK(object->id() == ids_before[0]);
    CHECK(volume->id() == ids_before[1]);
    CHECK(instance->id() == ids_before[2]);
    CHECK((volume->get_matrix().matrix().array() == volume_before.matrix().array()).all());
    CHECK((instance->get_matrix().matrix().array() == instance_before.matrix().array()).all());
    CHECK(volume->supported_facets.get_data() == paint_before);
    CHECK(volume->supported_facets.id() == paint_id);
    CHECK(volume->supported_facets.timestamp() == paint_timestamp);
}
