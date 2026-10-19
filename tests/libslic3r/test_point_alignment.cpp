// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include <catch2/catch.hpp>
#include "libslic3r/PointAlignment.hpp"
#include "libslic3r/FaceExtrusionAnnotations.hpp"
#include "libslic3r/Measure.hpp"
#include <array>
#include <limits>

using namespace Slic3r;

namespace {
struct AlignmentFixture {
    Model model;
    ModelObject *object;
    ModelInstance *instance;
    ModelVolume *fixed, *moving, *neighbor, *unrelated;
    AlignmentFixture()
    {
        object = model.add_object();
        fixed = object->add_volume(make_cube(10, 10, 10), false);
        moving = object->add_volume(make_cube(2, 3, 4), false);
        neighbor = object->add_volume(make_cube(3, 3, 3), false);
        moving->set_offset(Vec3d(20, 5, 3));
        neighbor->set_offset(Vec3d(-10, 0, 0));
        fixed->name = "FIXED"; moving->name = "MOVING"; neighbor->name = "NEIGHBOR";
        instance = object->add_instance();
        auto *other = model.add_object();
        unrelated = other->add_volume(make_cube(4, 4, 4), false);
        unrelated->set_offset(Vec3d(70, 20, 0));
        other->add_instance();
        for (auto *v : {fixed, moving, neighbor, unrelated})
            v->set_convex_hull_shared_ptr(v->get_mesh_shared_ptr()); // cubes are their own convex hull
    }
    PointAlignmentAnchor first(const Vec3d &local = Vec3d(10, 10, 10)) const
    { return point_alignment_anchor(*object, *instance, *fixed, local); }
    PointAlignmentAnchor second(const Vec3d &local = Vec3d::Zero()) const
    { return point_alignment_anchor(*object, *instance, *moving, local); }
};

struct SavedVolume {
    ModelVolume *volume;
    ObjectID id;
    std::shared_ptr<const TriangleMesh> mesh;
    Transform3d matrix;
    FaceExtrusionAnnotationStamp annotations;
    std::string name, material, guid;
    int extruder;
};
SavedVolume saved(ModelVolume *v)
{
    return {v, v->id(), v->get_mesh_shared_ptr(), v->get_matrix(), capture_face_extrusion_annotation_stamp(*v),
            v->name, v->material_id(), v->part_guid(), v->extruder_id()};
}
void unchanged(const SavedVolume &s, bool pose = true)
{
    CHECK(s.volume->id() == s.id);
    CHECK(s.volume->get_mesh_shared_ptr() == s.mesh);
    if (pose) CHECK(point_alignment_same_matrix(s.volume->get_matrix(), s.matrix));
    CHECK(capture_face_extrusion_annotation_stamp(*s.volume) == s.annotations);
    CHECK(s.volume->name == s.name);
    CHECK(s.volume->material_id() == s.material);
    CHECK(s.volume->part_guid() == s.guid);
    CHECK(s.volume->extruder_id() == s.extruder);
}
}

TEST_CASE("Point alignment preserves fixed and unrelated parts with one vector commit", "[point_alignment]")
{
    AlignmentFixture f;
    REQUIRE_FALSE(f.model.is_need_backup());
    f.moving->supported_facets.set_triangle_from_string(0, "4");
    f.neighbor->mmu_segmentation_facets.set_triangle_from_string(0, "8");
    f.moving->config.set_key_value("extruder", new ConfigOptionInt(3));
    f.moving->set_material_id("point-alignment-material");
    const std::array<SavedVolume, 4> before = {saved(f.fixed), saved(f.moving), saved(f.neighbor), saved(f.unrelated)};
    const Transform3d instance_before = f.instance->get_matrix();
    const auto fixed = f.first(), moving = f.second();
    auto plan = prepare_point_alignment(f.model, fixed, moving);
    REQUIRE(plan);
    REQUIRE_FALSE(plan.no_op);
    CHECK(plan.world_displacement == Vec3d(-10, 5, 7)); // REAL roles: first fixed, second moving
    for (const auto &s : before) unchanged(s); // planning has no mutations
    CHECK(point_alignment_same_matrix(f.instance->get_matrix(), instance_before));

    // Real Model mutation, deliberately NOT a claim to exercise GUI Undo/notifications.
    int snapshot_count = 0, vector_commits = 0;
    if (plan && !plan.no_op) {
        ++snapshot_count;
        f.moving->set_offset(plan.moving_matrix.translation());
        ++vector_commits;
    }
    CHECK(snapshot_count == 1);
    CHECK(vector_commits == 1);
    CHECK((f.second().world_point - fixed.world_point).norm() < 1e-8);
    CHECK((f.moving->get_matrix().linear().array() == before[1].matrix.linear().array()).all());
    unchanged(before[0]); unchanged(before[1], false); unchanged(before[2]); unchanged(before[3]);
    CHECK(point_alignment_same_matrix(f.instance->get_matrix(), instance_before));
    auto again = prepare_point_alignment(f.model, f.first(), f.second());
    REQUIRE(again);
    REQUIRE(again.no_op);
    if (again && !again.no_op) { ++snapshot_count; ++vector_commits; }
    CHECK(snapshot_count == 1);
    CHECK(vector_commits == 1);
}

TEST_CASE("Point alignment handles rotated nonuniform mirrored frames without changing the linear block", "[point_alignment]")
{
    for (int sample = 0; sample < 96; ++sample) {
        AlignmentFixture f;
        Transform3d instance = Transform3d::Identity(), moving = Transform3d::Identity();
        const Matrix3d ri = Eigen::AngleAxisd(0.071 * sample, Vec3d(1, 2, 3).normalized()).toRotationMatrix();
        const Matrix3d rv = Eigen::AngleAxisd(-0.113 * sample, Vec3d(3, 1, 2).normalized()).toRotationMatrix();
        instance.linear() = ri * Vec3d(sample & 1 ? -2.5 : 2.5, 0.6, sample & 2 ? -1.7 : 1.7).asDiagonal();
        instance.translation() = Vec3d(120.25, -53.5, 0);
        moving.linear() = rv * Vec3d(0.4, sample & 4 ? -3.0 : 3.0, 1.2).asDiagonal();
        moving.translation() = Vec3d(20, 8, 30);
        f.moving->set_transformation(moving);
        // Ground stationary geometry; moving geometry may legitimately sink, but
        // neither this preflight nor native do_move may shift the fixed geometry.
        double min_z = std::numeric_limits<double>::max();
        for (const auto *v : {f.fixed, f.neighbor})
            for (const Vec3f &p : v->mesh().its.vertices)
                min_z = std::min(min_z, (instance * v->get_matrix() * p.cast<double>()).z());
        instance.translation().z() -= min_z;
        f.instance->set_transformation(Geometry::Transformation(instance));
        const auto first = f.first(Vec3d(6, 4, 5)), second = f.second(Vec3d(0.3, 1.1, 2.7));
        const auto before = saved(f.moving);
        const Vec3d gap = sample % 3 == 0 ? Vec3d(25.4, -12.7, 2.54) : Vec3d::Zero();
        INFO("sample=" << sample);
        const auto plan = prepare_point_alignment(f.model, first, second, gap);
        REQUIRE(plan);
        CHECK((plan.moving_world_point - first.world_point - gap).norm() < 1e-8);
        CHECK((plan.moving_matrix.linear().array() == moving.linear().array()).all());
        const Vec3d legacy_local = (instance * moving).linear().inverse() * plan.world_displacement;
        Transform3d legacy = moving;
        legacy.translate(legacy_local);
        CHECK((plan.moving_matrix.matrix() - legacy.matrix()).norm() < 1e-10);
        unchanged(before); // pure prepare, including all mirrored/nonuniform cases
        f.moving->set_offset(plan.moving_matrix.translation());
        CHECK((f.second(Vec3d(0.3, 1.1, 2.7)).world_point - first.world_point - gap).norm() < 1e-8);
        const double min_after = f.object->get_instance_min_z(0);
        CHECK_FALSE((min_after > SINKING_Z_THRESHOLD && std::abs(min_after) > EPSILON));
    }
}

TEST_CASE("Point alignment keeps sub-centimetre offsets actionable and validates no-ops and nonfinite data", "[point_alignment]")
{
    for (double gap : {0.01, 0.005, 0.000001, 0.0}) {
        AlignmentFixture f;
        f.moving->set_offset(Vec3d(10 + gap, 10, 10));
        const auto plan = prepare_point_alignment(f.model, f.first(), f.second());
        REQUIRE(plan);
        CHECK(plan.no_op == (gap == 0));
        CHECK(std::abs(plan.world_displacement.x() + gap) < 1e-12);
    }
    AlignmentFixture f;
    auto first = f.first(), second = f.second();
    SECTION("nonfinite requested offset") {
        CHECK_FALSE(prepare_point_alignment(f.model, first, second, Vec3d(std::numeric_limits<double>::infinity(), 0, 0)));
    }
    SECTION("finite offset outside render/storage range") {
        CHECK_FALSE(prepare_point_alignment(f.model, first, second, Vec3d(1e100, 0, 0)));
    }
    SECTION("nonfinite local point") {
        second.local_point.x() = std::numeric_limits<double>::quiet_NaN();
        CHECK_FALSE(prepare_point_alignment(f.model, first, second));
    }
    SECTION("world-local mismatch") {
        second.world_point.z() += 0.02;
        CHECK_FALSE(prepare_point_alignment(f.model, first, second));
    }
    SECTION("singular transform") {
        Transform3d singular = f.moving->get_matrix();
        singular.linear().col(1).setZero();
        f.moving->set_transformation(singular);
        CHECK_FALSE(prepare_point_alignment(f.model, first, f.second()));
    }
}

TEST_CASE("Point alignment rejects stale identity mesh pose and cross-instance selections without touching the model", "[point_alignment]")
{
    AlignmentFixture f;
    auto first = f.first(), second = f.second();
    SECTION("new part ID") { f.moving->set_new_unique_id(); }
    SECTION("stale object ID") { first.object_id = second.object_id = f.unrelated->get_object()->id(); }
    SECTION("stale instance ID") { first.instance_id = second.instance_id = f.unrelated->get_object()->instances.front()->id(); }
    SECTION("same geometry different mesh pointer") { f.moving->set_mesh(f.moving->mesh()); }
    SECTION("moving pose changed") { f.moving->set_offset(Vec3d(21, 5, 3)); }
    SECTION("fixed pose changed") { f.fixed->set_offset(Vec3d(0.1, 0, 0)); }
    SECTION("instance pose changed") { f.instance->set_offset(Vec3d(0, 2, 0)); }
    SECTION("two points on same volume") { second = f.first(Vec3d::Zero()); }
    SECTION("different objects") {
        second = point_alignment_anchor(*f.unrelated->get_object(), *f.unrelated->get_object()->instances.front(), *f.unrelated, Vec3d::Zero());
    }
    SECTION("same object multiple instances even if picked instance matches") { f.object->add_instance(); }
    SECTION("different picked instance") {
        auto *copy = f.object->add_instance();
        second = point_alignment_anchor(*f.object, *copy, *f.moving, Vec3d::Zero());
    }
    SECTION("modifier instead of positive moving part") { f.moving->set_type(ModelVolumeType::PARAMETER_MODIFIER); }
    const auto before = saved(f.moving);
    REQUIRE_FALSE(prepare_point_alignment(f.model, first, second));
    unchanged(before);
}

TEST_CASE("Point alignment protects global exterior and native bed side effects but allows serialized legacy tags", "[point_alignment]")
{
    AlignmentFixture f;
    SECTION("legacy exterior is supported") {
        auto mesh = f.moving->mesh().its;
        mesh.properties.assign(mesh.indices.size(), FaceProperty{eExteriorAppearance, 1.0});
        f.moving->set_mesh(mesh);
        const auto before = saved(f.moving);
        CHECK_FALSE(point_alignment_has_exterior(f.model));
        CHECK(prepare_point_alignment(f.model, f.first(), f.second()));
        unchanged(before);
    }
    SECTION("first exterior payload only on unrelated modifier") {
        f.unrelated->set_type(ModelVolumeType::PARAMETER_MODIFIER);
        auto &data = const_cast<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>> &>(f.unrelated->exterior_facets.get_data());
        data.first.emplace_back(0, 0); // synthetic legacy/orphan payload, test only
        const auto before = saved(f.unrelated);
        CHECK(point_alignment_has_exterior(f.model));
        CHECK_FALSE(prepare_point_alignment(f.model, f.first(), f.second()));
        unchanged(before);
    }
    SECTION("second exterior payload only on unrelated object") {
        auto &data = const_cast<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>> &>(f.unrelated->exterior_facets.get_data());
        data.second.push_back(true);
        CHECK(point_alignment_has_exterior(f.model));
        CHECK_FALSE(prepare_point_alignment(f.model, f.first(), f.second()));
    }
    SECTION("floating unrelated object is not silently dropped") {
        f.unrelated->set_offset(Vec3d(70, 20, 3));
        const auto before = saved(f.unrelated);
        auto plan = prepare_point_alignment(f.model, f.first(), f.second());
        CHECK_FALSE(plan);
        CHECK(plan.message.find("bed") != std::string::npos);
        unchanged(before);
    }
    SECTION("moving the only bed-supporting part is rejected before mutation") {
        f.fixed->set_offset(Vec3d(0, 0, 10));
        f.neighbor->set_offset(Vec3d(-10, 0, 10));
        f.moving->set_offset(Vec3d(20, 5, 0));
        CHECK_FALSE(prepare_point_alignment(f.model, f.first(), f.second()));
    }
}

TEST_CASE("Point alignment plane gate rejects wrong plane IDs and off-plane picked points", "[point_alignment]")
{
    const Vec3d origin(7, -2, 3), normal(0, 0, 2);
    CHECK(point_alignment_plane_matches(4, 4, normal, origin, Vec3d(8, 6, 3)));
    CHECK_FALSE(point_alignment_plane_matches(4, 5, normal, origin, Vec3d(8, 6, 3)));
    CHECK_FALSE(point_alignment_plane_matches(-1, -1, normal, origin, Vec3d(8, 6, 3)));
    CHECK_FALSE(point_alignment_plane_matches(4, 4, Vec3d::Zero(), origin, Vec3d(8, 6, 3)));
    CHECK_FALSE(point_alignment_plane_matches(4, 4, normal, origin, Vec3d(8, 6, 3.001)));
    CHECK_FALSE(point_alignment_plane_matches(4, 4, normal, origin, Vec3d(std::numeric_limits<double>::quiet_NaN(), 0, 3)));
    CHECK(point_alignment_plane_matches(4, 4, normal, origin, Vec3d(300, 6, 3.000004))); // float ray rounding
    CHECK_FALSE(point_alignment_plane_matches(4, 4, normal, origin, Vec3d(300, 6, 3.05))); // raised GL highlight
}

TEST_CASE("Point alignment uses the real measured source plane rather than its raised highlight", "[point_alignment]")
{
    AlignmentFixture f;
    const auto &mesh = f.fixed->mesh().its;
    Measure::Measuring measuring(mesh);
    for (size_t index = 0; index < mesh.indices.size(); ++index) {
        const Vec3i triangle = mesh.indices[index];
        const Vec3d center = (mesh.vertices[triangle[0]].cast<double>() + mesh.vertices[triangle[1]].cast<double>() +
                              mesh.vertices[triangle[2]].cast<double>()) / 3.0;
        auto feature = measuring.get_feature(index, center, Transform3d::Identity(), true);
        REQUIRE(feature);
        REQUIRE(feature->get_type() == Measure::SurfaceFeatureType::Plane);
        const auto [plane_id, normal, point] = feature->get_plane();
        CHECK(point_alignment_plane_matches(plane_id, plane_id, normal, point, center));
        CHECK_FALSE(point_alignment_plane_matches(plane_id, plane_id, normal, point, center + 0.05 * normal.normalized()));
        CHECK_FALSE(point_alignment_plane_matches(plane_id, plane_id + 1, normal, point, center));
    }
}
