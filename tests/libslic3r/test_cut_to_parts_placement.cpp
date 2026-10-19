// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include <catch2/catch.hpp>

#include "libslic3r/CutUtils.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "slic3r/GUI/ModelObjectLoadPlacement.hpp"
#include "slic3r/GUI/Gizmos/NativeCutOutputPolicy.hpp"

#include <cereal/archives/binary.hpp>

#include <algorithm>
#include <array>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace Slic3r;

namespace {

constexpr double placement_tolerance_mm = 1e-5;

using WorldPoints = std::vector<Vec3d>;
using PaintData = std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>;

WorldPoints world_points(const ModelObject &object, const ModelVolume &volume)
{
    REQUIRE(object.instances.size() == 1);
    const Transform3d world = object.instances.front()->get_matrix() * volume.get_matrix();
    WorldPoints points;
    points.reserve(volume.mesh().its.vertices.size());
    for (const Vec3f &vertex : volume.mesh().its.vertices)
        points.emplace_back(world * vertex.cast<double>());
    return points;
}

std::array<Vec3d, 2> world_extents(const WorldPoints &points)
{
    REQUIRE_FALSE(points.empty());
    std::array<Vec3d, 2> extents{points.front(), points.front()};
    for (const Vec3d &point : points) {
        REQUIRE(point.allFinite());
        extents[0] = extents[0].cwiseMin(point);
        extents[1] = extents[1].cwiseMax(point);
    }
    return extents;
}

void check_extents(const WorldPoints &points, const Vec3d &minimum, const Vec3d &maximum)
{
    const auto extents = world_extents(points);
    CHECK((extents[0] - minimum).norm() <= placement_tolerance_mm);
    CHECK((extents[1] - maximum).norm() <= placement_tolerance_mm);
}

void check_translation(const WorldPoints &before, const WorldPoints &after, const Vec3d &translation)
{
    REQUIRE(after.size() == before.size());
    for (size_t i = 0; i < before.size(); ++i)
        CHECK((after[i] - before[i] - translation).norm() <= placement_tolerance_mm);
}

std::array<PaintData, 5> annotation_data(const ModelVolume &volume)
{
    return {volume.supported_facets.get_data(), volume.seam_facets.get_data(),
            volume.mmu_segmentation_facets.get_data(), volume.fuzzy_skin_facets.get_data(),
            volume.exterior_facets.get_data()};
}

struct SourceFixture {
    Model model;
    ModelObject *object;
    ModelVolume *volume;
    ModelInstance *instance;
    std::shared_ptr<const TriangleMesh> original_mesh;
    indexed_triangle_set original_geometry;
    Transform3d original_volume_matrix;
    Transform3d original_instance_matrix;
    std::array<PaintData, 5> original_annotations;
    ObjectID original_object_id;
    ObjectID original_volume_id;
    ObjectID original_instance_id;
    Vec3d source_world_offset;

    explicit SourceFixture(const Vec3d &offset, const Vec3d &volume_offset = Vec3d::Zero())
        : source_world_offset(offset)
    {
        object = model.add_object();
        object->name = "Native cut placement source";
        volume = object->add_volume(TriangleMesh(its_make_cube(20., 10., 30.)),
                                    ModelVolumeType::MODEL_PART, false);
        volume->name = "Source part";
        volume->set_offset(volume_offset);
        instance = object->add_instance();
        instance->set_offset(offset - volume_offset);
        object->invalidate_bounding_box();

        original_mesh = volume->get_mesh_shared_ptr();
        original_geometry = volume->mesh().its;
        original_volume_matrix = volume->get_matrix();
        original_instance_matrix = instance->get_matrix();
        original_annotations = annotation_data(*volume);
        original_object_id = object->id();
        original_volume_id = volume->id();
        original_instance_id = instance->id();
        check_extents(world_points(*object, *volume), offset, offset + Vec3d(20., 10., 30.));
    }

    Transform3d cut_matrix() const
    {
        // Same FFF world-to-instance-offset conversion as AdvancedCut::get_cut_matrix().
        const Vec3d plane_world = source_world_offset + Vec3d(10., 5., 12.);
        Transform3d matrix = Transform3d::Identity();
        matrix.translation() = plane_world - instance->get_offset();
        return matrix;
    }

    void check_unchanged() const
    {
        REQUIRE(model.objects.size() == 1);
        REQUIRE(object->volumes.size() == 1);
        REQUIRE(object->instances.size() == 1);
        CHECK(object->id() == original_object_id);
        CHECK(volume->id() == original_volume_id);
        CHECK(instance->id() == original_instance_id);
        CHECK(volume->get_mesh_shared_ptr() == original_mesh);
        CHECK(volume->mesh().its.vertices == original_geometry.vertices);
        CHECK(volume->mesh().its.indices == original_geometry.indices);
        REQUIRE(volume->mesh().its.properties.size() == original_geometry.properties.size());
        for (size_t i = 0; i < original_geometry.properties.size(); ++i) {
            CHECK(volume->mesh().its.properties[i].type == original_geometry.properties[i].type);
            CHECK(volume->mesh().its.properties[i].area == original_geometry.properties[i].area);
        }
        CHECK((volume->get_matrix().matrix().array() == original_volume_matrix.matrix().array()).all());
        CHECK((instance->get_matrix().matrix().array() == original_instance_matrix.matrix().array()).all());
        CHECK(annotation_data(*volume) == original_annotations);
        CHECK(object->name == "Native cut placement source");
        CHECK(volume->name == "Source part");
    }
};

const ModelVolume &named_part(const ModelObject &object, const std::string &name)
{
    const auto found = std::find_if(object.volumes.begin(), object.volumes.end(),
        [&name](const ModelVolume *volume) { return volume->name == name; });
    REQUIRE(found != object.volumes.end());
    REQUIRE((*found)->is_model_part());
    return **found;
}

ModelObject &load_with_placement_policy(Model &model, const ModelObject &cut_object,
    GUI::ModelObjectLoadPlacement placement = GUI::ModelObjectLoadPlacement::EnsureOnBed)
{
    // Real Model operations and the production placement helper used by the loader.
    // This deliberately does not claim to execute the wx/GL scene-loading path.
    ModelObject *loaded = model.add_object(cut_object);
    loaded->sort_volumes(true);
    GUI::apply_model_object_load_placement(*loaded, false, placement);
    return *loaded;
}

void check_source_corners_retained(const SourceFixture &source, const ModelObject &result)
{
    WorldPoints result_points;
    for (const ModelVolume *volume : result.volumes) {
        const WorldPoints points = world_points(result, *volume);
        result_points.insert(result_points.end(), points.begin(), points.end());
    }
    for (const Vec3d &corner : world_points(*source.object, *source.volume)) {
        CHECK(std::any_of(result_points.begin(), result_points.end(), [&corner](const Vec3d &point) {
            return (point - corner).norm() <= placement_tolerance_mm;
        }));
    }
}

void apply_source_pose(SourceFixture &source, int pose)
{
    if (pose != 0) {
        source.volume->set_rotation(Vec3d(0.31, -0.22, 0.17));
        source.volume->set_scaling_factor(Vec3d(1.1, 0.8, 1.3));
        source.instance->set_rotation(Vec3d(-0.13, 0.27, 0.41));
        source.instance->set_scaling_factor(Vec3d(0.9, 1.2, 1.05));
        if (pose == 2)
            source.volume->set_mirror(Vec3d(-1., 1., 1.));
        if (pose == 3)
            source.instance->set_mirror(Vec3d(1., -1., 1.));
        source.object->invalidate_bounding_box();
        source.original_volume_matrix = source.volume->get_matrix();
        source.original_instance_matrix = source.instance->get_matrix();
    }
}

Transform3d posed_cut_matrix(const SourceFixture &source)
{
    const Transform3d world = source.original_instance_matrix * source.original_volume_matrix;
    const Vec3d normal = (world.linear().inverse().transpose() * Vec3d::UnitZ()).normalized();
    Transform3d cut = Transform3d::Identity();
    cut.linear() = Eigen::Quaterniond::FromTwoVectors(Vec3d::UnitZ(), normal).toRotationMatrix();
    // The same world-center minus instance-offset convention used by AdvancedCut.
    cut.translation() = world * Vec3d(10., 5., 12.) - source.instance->get_offset();
    return cut;
}

std::array<WorldPoints, 2> cut_world_sides(const ModelObjectPtrs &objects)
{
    std::array<WorldPoints, 2> sides;
    REQUIRE((objects.size() == 1 || objects.size() == 2));
    if (objects.size() == 1) {
        REQUIRE(objects.front()->parts_count() == 2);
        sides[0] = world_points(*objects.front(), named_part(*objects.front(), "Source part_A"));
        sides[1] = world_points(*objects.front(), named_part(*objects.front(), "Source part_B"));
    } else {
        for (size_t side = 0; side != 2; ++side) {
            REQUIRE(objects[side]->volumes.size() == 1);
            sides[side] = world_points(*objects[side], *objects[side]->volumes.front());
        }
    }
    return sides;
}

void check_world_half(const SourceFixture &source, const WorldPoints &actual, bool upper)
{
    INFO("Half: " << (upper ? "A" : "B"));
    REQUIRE_FALSE(actual.empty());
    const Transform3d world = source.original_instance_matrix * source.original_volume_matrix;
    const Transform3d inverse = world.inverse();
    const Vec3d minimum(0., 0., upper ? 12. : 0.);
    const Vec3d maximum(20., 10., upper ? 30. : 12.);
    double max_outside_mm = 0.;
    double max_surface_error_mm = 0.;
    for (const Vec3d &point : actual) {
        REQUIRE(point.allFinite());
        const Vec3d local = inverse * point;
        double surface_distance_mm = std::numeric_limits<double>::infinity();
        for (int axis = 0; axis != 3; ++axis) {
            // Inverse-transpose plane distances are physical world millimeters,
            // also under nonuniform scaling or composed volume/instance poses.
            const double normal_length = inverse.linear().row(axis).norm();
            const double lower_distance = (local[axis] - minimum[axis]) / normal_length;
            const double upper_distance = (maximum[axis] - local[axis]) / normal_length;
            max_outside_mm = std::max({max_outside_mm, -lower_distance, -upper_distance});
            surface_distance_mm = std::min({surface_distance_mm, std::abs(lower_distance), std::abs(upper_distance)});
        }
        max_surface_error_mm = std::max(max_surface_error_mm, surface_distance_mm);
    }
    CHECK(max_outside_mm <= placement_tolerance_mm);
    CHECK(max_surface_error_mm <= placement_tolerance_mm);
    // All original and cut-plane corners must remain at their exact intended
    // world support. This checks each half, not just the union bounding box.
    for (int corner = 0; corner != 8; ++corner) {
        Vec3d local;
        for (int axis = 0; axis != 3; ++axis)
            local[axis] = (corner & (1 << axis)) ? maximum[axis] : minimum[axis];
        const Vec3d expected = world * local;
        double nearest_mm = std::numeric_limits<double>::infinity();
        for (const Vec3d &point : actual)
            nearest_mm = std::min(nearest_mm, (point - expected).norm());
        CHECK(nearest_mm <= placement_tolerance_mm);
    }
}

void check_production_cut_placement(SourceFixture &source, const GUI::NativeCutOutputDecision &decision)
{
    REQUIRE(decision.allowed());
    REQUIRE(decision.route == GUI::NativeCutRoute::Plane);
    Cut cut(source.object, 0, posed_cut_matrix(source), decision.attributes);
    cut.set_offset_for_two_part = decision.separate_half_offset;
    const ModelObjectPtrs &results = cut.perform_with_plane();
    REQUIRE_FALSE(cut.was_canceled());
    CHECK(results.size() == 1);
    auto before_load = cut_world_sides(results);
    for (size_t side = 0; side != 2; ++side) {
        INFO("Before loader");
        check_world_half(source, before_load[side], side == 0);
    }

    Model loaded;
    for (const ModelObject *object : results)
        load_with_placement_policy(loaded, *object,
            GUI::cut_result_load_placement(decision.cut_to_parts()));
    const auto after_load = cut_world_sides(loaded.objects);
    for (size_t side = 0; side != 2; ++side) {
        INFO("After production loader placement");
        check_world_half(source, after_load[side], side == 0);
        check_translation(before_load[side], after_load[side], Vec3d::Zero());
    }
    source.check_unchanged();
}

std::array<bool, 7> output_values(const GUI::NativeCutOutputSettings &settings)
{
    return {settings.keep_upper, settings.keep_lower, settings.cut_to_parts,
            settings.place_on_cut_upper, settings.place_on_cut_lower, settings.flip_upper, settings.flip_lower};
}

} // namespace

TEST_CASE("Native Cut to parts preserves world placement only with the explicit loader policy", "[native_cut_parts]")
{
    const ModelObjectCutAttributes keep_both = ModelObjectCutAttribute::KeepUpper | ModelObjectCutAttribute::KeepLower;

    SECTION("On-bed and later-plate cuts retain both parts in their original world positions") {
        for (const Vec3d &offset : {Vec3d(0., 0., 0.), Vec3d(450., 100., 0.)}) {
            INFO("Source offset: " << offset.transpose());
            // Nonzero volume and instance offsets must compose to the original world pose.
            SourceFixture source(offset, Vec3d(3., -2., 4.));
            Cut cut(source.object, 0, source.cut_matrix(), ModelObjectCutAttribute::CutToParts | keep_both);
            // AdvancedCut sets this even when CutToParts is enabled. The parts branch must bypass it.
            cut.set_offset_for_two_part = true;
            const ModelObjectPtrs &results = cut.perform_with_plane();
            REQUIRE_FALSE(cut.was_canceled());
            REQUIRE(results.size() == 1);
            const ModelObject &result = *results.front();
            REQUIRE(result.volumes.size() == 2);
            REQUIRE(result.parts_count() == 2);
            const WorldPoints upper = world_points(result, named_part(result, "Source part_A"));
            const WorldPoints lower = world_points(result, named_part(result, "Source part_B"));
            check_extents(upper, offset + Vec3d(0., 0., 12.), offset + Vec3d(20., 10., 30.));
            check_extents(lower, offset, offset + Vec3d(20., 10., 12.));
            check_source_corners_retained(source, result);

            Model loaded_model;
            const ModelObject &loaded = load_with_placement_policy(loaded_model, result,
                GUI::ModelObjectLoadPlacement::PreserveExisting);
            REQUIRE(loaded_model.objects.size() == 1);
            REQUIRE(loaded.parts_count() == 2);
            check_translation(upper, world_points(loaded, named_part(loaded, "Source part_A")), Vec3d::Zero());
            check_translation(lower, world_points(loaded, named_part(loaded, "Source part_B")), Vec3d::Zero());
            check_source_corners_retained(source, loaded);
            source.check_unchanged();
        }
    }

    SECTION("Keep orientation alone still separates objects and beds each half independently") {
        const Vec3d offset(450., 100., 0.);
        SourceFixture source(offset);
        // Keep orientation means neither PlaceOnCut nor Flip; it is not CutToParts.
        Cut cut(source.object, 0, source.cut_matrix(), keep_both);
        cut.set_offset_for_two_part = true;
        const ModelObjectPtrs &results = cut.perform_with_plane();
        REQUIRE_FALSE(cut.was_canceled());
        REQUIRE(results.size() == 2);
        const ModelObject &upper = *results[0];
        const ModelObject &lower = *results[1];
        REQUIRE(upper.volumes.size() == 1);
        REQUIRE(lower.volumes.size() == 1);
        const WorldPoints upper_before = world_points(upper, *upper.volumes.front());
        const WorldPoints lower_before = world_points(lower, *lower.volumes.front());
        // The native separate-object path adds 1.1 times the upper mesh's 20 mm X extent.
        check_extents(upper_before, offset + Vec3d(22., 0., 12.), offset + Vec3d(42., 10., 30.));
        check_extents(lower_before, offset, offset + Vec3d(20., 10., 12.));

        Model loaded_model;
        const ModelObject &loaded_upper = load_with_placement_policy(loaded_model, upper);
        const ModelObject &loaded_lower = load_with_placement_policy(loaded_model, lower);
        REQUIRE(loaded_model.objects.size() == 2);
        const WorldPoints upper_after = world_points(loaded_upper, *loaded_upper.volumes.front());
        const WorldPoints lower_after = world_points(loaded_lower, *loaded_lower.volumes.front());
        check_translation(upper_before, upper_after, Vec3d(0., 0., -12.));
        check_translation(lower_before, lower_after, Vec3d::Zero());
        check_extents(upper_after, offset + Vec3d(22., 0., 0.), offset + Vec3d(42., 10., 18.));
        check_extents(lower_after, offset, offset + Vec3d(20., 10., 12.));
        source.check_unchanged();
    }

    SECTION("PreserveExisting retains floating and sunk parts while EnsureOnBed still moves the whole object") {
        for (double bottom_z : {7., -7.}) {
            INFO("Source bottom Z: " << bottom_z);
            const Vec3d offset(450., 100., bottom_z);
            // Nonzero volume and instance offsets must compose to the original world pose.
            SourceFixture source(offset, Vec3d(3., -2., 4.));
            Cut cut(source.object, 0, source.cut_matrix(), ModelObjectCutAttribute::CutToParts | keep_both);
            cut.set_offset_for_two_part = true;
            const ModelObjectPtrs &results = cut.perform_with_plane();
            REQUIRE_FALSE(cut.was_canceled());
            REQUIRE(results.size() == 1);
            const ModelObject &result = *results.front();
            REQUIRE(result.parts_count() == 2);
            const WorldPoints upper = world_points(result, named_part(result, "Source part_A"));
            const WorldPoints lower = world_points(result, named_part(result, "Source part_B"));
            check_extents(upper, offset + Vec3d(0., 0., 12.), offset + Vec3d(20., 10., 30.));
            check_extents(lower, offset, offset + Vec3d(20., 10., 12.));
            check_source_corners_retained(source, result);

            Model preserved_model;
            const ModelObject &preserved = load_with_placement_policy(preserved_model, result,
                GUI::ModelObjectLoadPlacement::PreserveExisting);
            REQUIRE(preserved_model.objects.size() == 1);
            REQUIRE(preserved.parts_count() == 2);
            check_translation(upper, world_points(preserved, named_part(preserved, "Source part_A")), Vec3d::Zero());
            check_translation(lower, world_points(preserved, named_part(preserved, "Source part_B")), Vec3d::Zero());
            check_source_corners_retained(source, preserved);

            Model loaded_model;
            const ModelObject &loaded = load_with_placement_policy(loaded_model, result);
            REQUIRE(loaded_model.objects.size() == 1);
            REQUIRE(loaded.parts_count() == 2);
            const WorldPoints upper_after = world_points(loaded, named_part(loaded, "Source part_A"));
            const WorldPoints lower_after = world_points(loaded, named_part(loaded, "Source part_B"));
            const Vec3d whole_object_shift(0., 0., -bottom_z);
            check_translation(upper, upper_after, whole_object_shift);
            check_translation(lower, lower_after, whole_object_shift);
            check_extents(upper_after, Vec3d(450., 100., 12.), Vec3d(470., 110., 30.));
            check_extents(lower_after, Vec3d(450., 100., 0.), Vec3d(470., 110., 12.));
            source.check_unchanged();
        }
    }
}

TEST_CASE("Native Cut has preserving defaults and a complete output archive", "[native_cut_parts][native_cut_policy_characterization]")
{
    GUI::NativeCutOutputSettings settings;
    CHECK(settings.keep_upper);
    CHECK(settings.keep_lower);
    CHECK(settings.cut_to_parts);
    CHECK_FALSE(settings.place_on_cut_upper);
    CHECK_FALSE(settings.place_on_cut_lower);
    CHECK_FALSE(settings.flip_upper);
    CHECK_FALSE(settings.flip_lower);
    auto decision = GUI::resolve_native_cut_output(settings, {});
    CHECK(decision.allowed());
    CHECK(decision.route == GUI::NativeCutRoute::Plane);
    CHECK_FALSE(decision.separate_half_offset);
    CHECK(decision.cut_to_parts());
    CHECK(decision.attributes.has(ModelObjectCutAttribute::KeepUpper));
    CHECK(decision.attributes.has(ModelObjectCutAttribute::KeepLower));
    CHECK_FALSE(decision.attributes.has(ModelObjectCutAttribute::PlaceOnCutUpper));
    CHECK(decision.load_placement() == GUI::ModelObjectLoadPlacement::PreserveExisting);

    settings.keep_upper = false;
    settings.cut_to_parts = true;
    settings.flip_lower = true;
    const unsigned sentinel = 0x12345678;
    std::ostringstream extracted, expected, legacy_prefix;
    {
        cereal::BinaryOutputArchive archive(extracted);
        archive(settings, sentinel);
    }
    {
        cereal::BinaryOutputArchive archive(expected);
        archive(settings.keep_upper, settings.keep_lower, settings.flip_lower, settings.flip_upper,
                settings.cut_to_parts, settings.place_on_cut_upper, settings.place_on_cut_lower, sentinel);
    }
    {
        cereal::BinaryOutputArchive archive(legacy_prefix);
        archive(settings.keep_upper, settings.keep_lower, settings.flip_lower, settings.flip_upper);
    }
    CHECK(extracted.str() == expected.str());
    CHECK(extracted.str().substr(0, legacy_prefix.str().size()) == legacy_prefix.str());
}

TEST_CASE("Native Cut production defaults and Reset preserve both world halves", "[native_cut_parts][native_cut_policy_desired]")
{
    for (bool reset : {false, true}) {
        GUI::NativeCutOutputSettings settings;
        if (reset) {
            settings.keep_upper = settings.keep_lower = false;
            settings.cut_to_parts = true;
            settings.place_on_cut_lower = settings.flip_upper = settings.flip_lower = true;
            settings.reset();
        }
        INFO("Reset: " << reset);
        const auto decision = GUI::resolve_native_cut_output(settings, {});
        CHECK(decision.cut_to_parts());
        CHECK(decision.load_placement() == GUI::ModelObjectLoadPlacement::PreserveExisting);
        for (int pose = 0; pose != 4; ++pose) {
            INFO("Pose: " << pose);
            SourceFixture source(Vec3d(450., 100., pose == 0 ? 0. : 7.), Vec3d(3., -2., 4.));
            apply_source_pose(source, pose);
            check_production_cut_placement(source, decision);
        }
    }
}

TEST_CASE("Native Cut explicitly selected Parts supports composed poses without moving either half", "[native_cut_parts][native_cut_policy_control]")
{
    GUI::NativeCutOutputSettings settings;
    settings.cut_to_parts = true;
    settings.normalize_for_render();
    const auto decision = GUI::resolve_native_cut_output(settings, {});
    REQUIRE(decision.cut_to_parts());
    for (int pose = 0; pose != 4; ++pose) {
        INFO("Pose: " << pose);
        SourceFixture source(Vec3d(450., 100., pose == 0 ? 0. : 7.), Vec3d(3., -2., 4.));
        apply_source_pose(source, pose);
        check_production_cut_placement(source, decision);
    }
}

TEST_CASE("Native Cut output archive restores all placement decisions", "[native_cut_parts][native_cut_policy_desired]")
{
    for (unsigned pattern : {0u, 1u, 2u, 4u, 8u, 16u, 32u, 64u, 127u}) {
        INFO("Pattern: " << pattern);
        GUI::NativeCutOutputSettings saved;
        saved.keep_upper = (pattern & 1u) != 0;
        saved.keep_lower = (pattern & 2u) != 0;
        saved.cut_to_parts = (pattern & 4u) != 0;
        saved.place_on_cut_upper = (pattern & 8u) != 0;
        saved.place_on_cut_lower = (pattern & 16u) != 0;
        saved.flip_upper = (pattern & 32u) != 0;
        saved.flip_lower = (pattern & 64u) != 0;
        GUI::NativeCutOutputSettings restored;
        restored.keep_upper = !saved.keep_upper;
        restored.keep_lower = !saved.keep_lower;
        restored.cut_to_parts = !saved.cut_to_parts;
        restored.place_on_cut_upper = !saved.place_on_cut_upper;
        restored.place_on_cut_lower = !saved.place_on_cut_lower;
        restored.flip_upper = !saved.flip_upper;
        restored.flip_lower = !saved.flip_lower;
        std::stringstream buffer;
        {
            cereal::BinaryOutputArchive archive(buffer);
            archive(saved);
        }
        {
            cereal::BinaryInputArchive archive(buffer);
            archive(restored);
        }
        const auto wanted = output_values(saved);
        const auto actual = output_values(restored);
        for (size_t field = 0; field != wanted.size(); ++field) {
            INFO("Output field: " << field);
            CHECK(actual[field] == wanted[field]);
        }
    }
}

TEST_CASE("Native Cut temporary part reassignment retains the requested in-place output", "[native_cut_parts][native_cut_policy_desired]")
{
    GUI::NativeCutOutputSettings settings;
    settings.cut_to_parts = true;
    settings.normalize_for_render();
    for (bool modified : {false, true, false}) {
        INFO("Modified assignments: " << modified);
        // This is the actual state transition called before the GUI checkbox,
        // not an imitation of mouse events or of PartSelection's raycaster.
        settings.on_part_assignment_changed(modified);
        settings.normalize_for_render();
        GUI::NativeCutOutputContext context;
        context.has_valid_modified_parts = modified;
        const auto decision = GUI::resolve_native_cut_output(settings, context);
        CHECK(settings.cut_to_parts);
        CHECK(decision.allowed() == !modified);
    }
    const auto decision = GUI::resolve_native_cut_output(settings, {});
    CHECK(decision.cut_to_parts());
    SourceFixture source(Vec3d(450., 100., 0.), Vec3d(3., -2., 4.));
    check_production_cut_placement(source, decision);
}

TEST_CASE("Native Cut requested in-place output rejects incompatible consumers without clearing intent", "[native_cut_parts][native_cut_policy_desired]")
{
    for (int incompatible = 0; incompatible != 3; ++incompatible) {
        INFO("Unsupported context: " << incompatible);
        GUI::NativeCutOutputSettings settings;
        settings.cut_to_parts = true;
        settings.normalize_for_render();
        GUI::NativeCutOutputContext context;
        context.has_connectors = incompatible == 0;
        context.mode = incompatible == 1 ? CutMode::cutTongueAndGroove : CutMode::cutPlanar;
        context.has_valid_modified_parts = incompatible == 2;
        settings.on_mode_changed(context.mode);
        settings.on_part_assignment_changed(context.has_valid_modified_parts);
        settings.normalize_for_render();
        const auto decision = GUI::resolve_native_cut_output(settings, context);
        CHECK(settings.cut_to_parts);
        CHECK_FALSE(decision.allowed());
        const auto expected = incompatible == 0 ? GUI::NativeCutOutputBlock::InPlaceWithConnectors :
            incompatible == 1 ? GUI::NativeCutOutputBlock::InPlaceWithGroove : GUI::NativeCutOutputBlock::InPlaceWithReassignment;
        CHECK(decision.block == expected);

        // An explicit separate-output selection keeps the existing consumers.
        settings.cut_to_parts = false;
        const auto separate = GUI::resolve_native_cut_output(settings, context);
        CHECK(separate.allowed());
        CHECK_FALSE(separate.cut_to_parts());
        CHECK(separate.separate_half_offset);
        CHECK(separate.load_placement() == GUI::ModelObjectLoadPlacement::EnsureOnBed);
        CHECK(separate.route == (incompatible == 1 ? GUI::NativeCutRoute::Groove :
            incompatible == 2 ? GUI::NativeCutRoute::Contour : GUI::NativeCutRoute::Plane));
    }
}

TEST_CASE("Native Cut action normalizes requested in-place flags without a render", "[native_cut_parts][native_cut_policy_desired]")
{
    GUI::NativeCutOutputSettings settings;
    settings.cut_to_parts = true;
    settings.keep_upper = settings.keep_lower = false;
    settings.place_on_cut_upper = settings.place_on_cut_lower = true;
    settings.flip_upper = settings.flip_lower = true;
    auto decision = GUI::resolve_native_cut_output(settings, {});
    CHECK(decision.allowed());
    CHECK(decision.attributes.has(ModelObjectCutAttribute::KeepUpper));
    CHECK(decision.attributes.has(ModelObjectCutAttribute::KeepLower));
    CHECK_FALSE(decision.attributes.has(ModelObjectCutAttribute::PlaceOnCutUpper));
    CHECK_FALSE(decision.attributes.has(ModelObjectCutAttribute::PlaceOnCutLower));
    CHECK_FALSE(decision.attributes.has(ModelObjectCutAttribute::FlipUpper));
    CHECK_FALSE(decision.attributes.has(ModelObjectCutAttribute::FlipLower));
    CHECK_FALSE(decision.separate_half_offset);
}

TEST_CASE("Native Cut Reset restores in-place output after an explicit separate-output choice", "[native_cut_parts][native_cut_policy_extra]")
{
    GUI::NativeCutOutputSettings settings;
    settings.cut_to_parts = false;
    settings.keep_upper = settings.keep_lower = false;
    settings.place_on_cut_upper = settings.place_on_cut_lower = true;
    settings.flip_upper = settings.flip_lower = true;
    settings.reset();
    CHECK(settings.cut_to_parts);
    CHECK(settings.keep_upper);
    CHECK(settings.keep_lower);
    CHECK_FALSE(settings.place_on_cut_upper);
    CHECK_FALSE(settings.place_on_cut_lower);
    CHECK_FALSE(settings.flip_upper);
    CHECK_FALSE(settings.flip_lower);
    const auto decision = GUI::resolve_native_cut_output(settings, {});
    CHECK(decision.cut_to_parts());
    CHECK_FALSE(decision.separate_half_offset);
    CHECK(decision.load_placement() == GUI::ModelObjectLoadPlacement::PreserveExisting);
    SourceFixture source(Vec3d(450., 100., 7.), Vec3d(3., -2., 4.));
    check_production_cut_placement(source, decision);
}

TEST_CASE("Native Cut explicit separate output retains single-half behavior and rejects keep-neither", "[native_cut_parts][native_cut_policy_extra]")
{
    for (bool keep_upper : {true, false}) {
        INFO("Keep only " << (keep_upper ? "A" : "B"));
        GUI::NativeCutOutputSettings settings;
        settings.cut_to_parts = false;
        settings.keep_upper = keep_upper;
        settings.keep_lower = !keep_upper;
        auto decision = GUI::resolve_native_cut_output(settings, {});
        REQUIRE(decision.allowed());
        REQUIRE(decision.route == GUI::NativeCutRoute::Plane);
        CHECK(decision.attributes.has(ModelObjectCutAttribute::KeepUpper) == keep_upper);
        CHECK(decision.attributes.has(ModelObjectCutAttribute::KeepLower) == !keep_upper);
        CHECK_FALSE(decision.cut_to_parts());
        CHECK_FALSE(decision.attributes.has(ModelObjectCutAttribute::PlaceOnCutUpper));
        CHECK_FALSE(decision.attributes.has(ModelObjectCutAttribute::PlaceOnCutLower));
        CHECK_FALSE(decision.attributes.has(ModelObjectCutAttribute::FlipUpper));
        CHECK_FALSE(decision.attributes.has(ModelObjectCutAttribute::FlipLower));
        CHECK(decision.separate_half_offset);
        CHECK(decision.load_placement() == GUI::ModelObjectLoadPlacement::EnsureOnBed);

        SourceFixture source(Vec3d(450., 100., 7.), Vec3d(3., -2., 4.));
        Cut cut(source.object, 0, source.cut_matrix(), decision.attributes);
        cut.set_offset_for_two_part = decision.separate_half_offset;
        const ModelObjectPtrs &results = cut.perform_with_plane();
        REQUIRE_FALSE(cut.was_canceled());
        REQUIRE(results.size() == 1);
        REQUIRE(results.front()->parts_count() == 1);
        REQUIRE(results.front()->volumes.size() == 1);
        const WorldPoints before = world_points(*results.front(), *results.front()->volumes.front());
        // Native single-half output does not add the two-half XY separation.
        check_world_half(source, before, keep_upper);

        Model loaded;
        const ModelObject &object = load_with_placement_policy(loaded, *results.front(),
            GUI::cut_result_load_placement(decision.cut_to_parts()));
        REQUIRE(object.parts_count() == 1);
        const WorldPoints after = world_points(object, *object.volumes.front());
        const double original_bottom_z = 7. + (keep_upper ? 12. : 0.);
        check_translation(before, after, Vec3d(0., 0., -original_bottom_z));
        check_extents(after, Vec3d(450., 100., 0.), Vec3d(470., 110., keep_upper ? 18. : 12.));
        source.check_unchanged();
    }

    GUI::NativeCutOutputSettings neither;
    neither.cut_to_parts = false;
    neither.keep_upper = neither.keep_lower = false;
    auto blocked = GUI::resolve_native_cut_output(neither, {});
    CHECK_FALSE(blocked.allowed());
    CHECK(blocked.block == GUI::NativeCutOutputBlock::KeepNeither);
    CHECK_FALSE(blocked.attributes.has(ModelObjectCutAttribute::KeepUpper));
    CHECK_FALSE(blocked.attributes.has(ModelObjectCutAttribute::KeepLower));
    CHECK_FALSE(blocked.cut_to_parts());
    // No Cut/Model mutation is dispatched for a denied request.
}
