// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#ifndef slic3r_FaceExtrusionModelPolicy_hpp_
#define slic3r_FaceExtrusionModelPolicy_hpp_

#include "libslic3r/Model.hpp"
#include "libslic3r/FaceExtrusionAnnotations.hpp"
#include "libslic3r/MeshSolidValidation.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Slic3r::GUI {

// The real gizmo and native regressions share these Model decisions. Geometry,
// wx/GL interaction, annotation consent and the Undo transaction stay with their
// existing owners. Only Faces may edit a negative source; the GUI additionally
// requires that negative volume to be the single explicitly selected target.
enum class FaceExtrusionToolMode { Faces, PointProfile, CutAndStretch };

constexpr double face_extrusion_bed_tolerance_mm = 1e-4;

inline bool face_extrusion_source_role_allowed(FaceExtrusionToolMode mode, ModelVolumeType type)
{
    return type == ModelVolumeType::MODEL_PART ||
        (mode == FaceExtrusionToolMode::Faces && type == ModelVolumeType::NEGATIVE_VOLUME);
}

inline bool face_extrusion_selection_roles_allowed(
    FaceExtrusionToolMode mode, const std::vector<ModelVolumeType> &types)
{
    if (types.size() == 1 && mode == FaceExtrusionToolMode::Faces && types.front() == ModelVolumeType::NEGATIVE_VOLUME)
        return true;
    return !types.empty() && std::all_of(types.begin(), types.end(),
        [](ModelVolumeType type) { return type == ModelVolumeType::MODEL_PART; });
}

inline bool face_extrusion_source_type_matches(
    FaceExtrusionToolMode mode, ModelVolumeType expected, const ModelVolume &live)
{
    return live.type() == expected && face_extrusion_source_role_allowed(mode, live.type());
}

namespace FaceExtrusionModelDetail {

inline bool usable_transform(const Transform3d &matrix)
{
    if (!matrix.matrix().allFinite())
        return false;
    const Matrix3d linear = matrix.linear();
    const double scale = linear.cwiseAbs().maxCoeff();
    // A relative determinant test also rejects numerically singular transforms.
    return scale > 0.0 && std::isfinite(scale) &&
           std::abs((linear / scale).determinant()) > 1e-12 &&
           linear.inverse().allFinite();
}

inline bool usable_mesh(const indexed_triangle_set &mesh)
{
    if (mesh.vertices.empty() || mesh.indices.empty() ||
        mesh.indices.size() > size_t(std::numeric_limits<int>::max()))
        return false;
    for (const Vec3f &vertex : mesh.vertices)
        if (!vertex.allFinite())
            return false;
    for (const Vec3i &face : mesh.indices) {
        for (int i = 0; i < 3; ++i)
            if (face[i] < 0 || size_t(face[i]) >= mesh.vertices.size())
                return false;
        const Vec3d a = mesh.vertices[face[0]].cast<double>();
        const Vec3d b = mesh.vertices[face[1]].cast<double>();
        const Vec3d c = mesh.vertices[face[2]].cast<double>();
        if ((b - a).cross(c - a).squaredNorm() <= 0.0)
            return false;
    }
    return true;
}

inline bool finite_bounds(const BoundingBoxf3 &bounds)
{
    return bounds.defined && bounds.min.allFinite() && bounds.max.allFinite();
}

} // namespace FaceExtrusionModelDetail

struct FaceExtrusionPartPreparation {
    TriangleMesh mesh;
    bool valid = false;
    std::string message;

    explicit operator bool() const noexcept { return valid; }
};

// Shared, non-mutating NewPart preparation. The caller retains its transform,
// stale-state, source-role, Model and Undo checks. All geometry stays local until
// world, stored and materialized float-roundtrip solids have passed validation;
// any rejected preparation returns an empty mesh and a diagnostic.
inline FaceExtrusionPartPreparation prepare_face_extrusion_part_storage(
    const TriangleMesh &world_mesh, const Transform3d &instance_matrix,
    double max_world_error = face_extrusion_bed_tolerance_mm)
{
    FaceExtrusionPartPreparation result;
    if (!FaceExtrusionModelDetail::usable_transform(instance_matrix)) {
        result.message = "The new part has an invalid instance transform.";
        return result;
    }
    if (!std::isfinite(max_world_error) || max_world_error < 0.0) {
        result.message = "The new part has an invalid positional tolerance.";
        return result;
    }

    TriangleMesh mesh;
    const char *stage = "New part world mesh";
    try {
        validate_solid_mesh(world_mesh.its);
        stage = "New part object-space mesh";
        mesh = world_mesh;
        mesh.transform(instance_matrix.inverse(), true);
        validate_solid_mesh(mesh.its);

        // Retain the original double per-vertex positional budget. A small
        // position error alone cannot certify topology after float narrowing.
        for (size_t i = 0; i < mesh.its.vertices.size(); ++i) {
            const Vec3d round_trip = instance_matrix * mesh.its.vertices[i].cast<double>();
            if (!round_trip.allFinite() ||
                (round_trip - world_mesh.its.vertices[i].cast<double>()).norm() > max_world_error) {
                result.message = "A new part loses precision in object coordinates.";
                return result;
            }
        }

        stage = "New part float world-coordinate round trip";
        TriangleMesh round_trip_mesh = mesh;
        round_trip_mesh.transform(instance_matrix, true);
        validate_solid_mesh(round_trip_mesh.its);
        for (size_t i = 0; i < round_trip_mesh.its.vertices.size(); ++i) {
            const Vec3d round_trip = round_trip_mesh.its.vertices[i].cast<double>();
            if ((round_trip - world_mesh.its.vertices[i].cast<double>()).norm() > max_world_error) {
                result.message = "A new part loses precision after the world-coordinate round trip.";
                return result;
            }
        }
    } catch (const std::exception &error) {
        result.message = std::string(stage) + ": " + error.what();
        return result;
    }
    result.mesh = std::move(mesh);
    result.valid = true;
    return result;
}

struct FaceExtrusionVolumeState {
    ObjectID id;
    ModelVolumeType type = ModelVolumeType::MODEL_PART;
    std::shared_ptr<const TriangleMesh> identity;
    std::shared_ptr<const TriangleMesh> hull_identity;
    Transform3d matrix = Transform3d::Identity();
    TriangleMesh bed_mesh; // COPY of other-part hull, or min-Z fallback mesh
    bool has_hull = false;
};

// Only immutable mesh identities, copied geometry and value metadata cross into
// the worker. The caller keeps the existing all-volume stale-state checks.
struct FaceExtrusionPlacementContext {
    ObjectID source_id;
    ModelVolumeType source_type = ModelVolumeType::MODEL_PART;
    Transform3d world = Transform3d::Identity();
    Transform3d instance_matrix = Transform3d::Identity();
    Transform3d instance_without_translation = Transform3d::Identity();
    Transform3d volume_matrix = Transform3d::Identity();
    double instance_offset_z = 0.0;
    size_t model_parts = 0;
    std::vector<FaceExtrusionVolumeState> volumes;
};

inline FaceExtrusionPlacementContext capture_face_extrusion_placement(
    const ModelObject &object, const ModelInstance &instance, const ModelVolume &source,
    const Transform3d &world)
{
    FaceExtrusionPlacementContext context;
    context.source_id = source.id();
    context.source_type = source.type();
    context.world = world;
    context.instance_matrix = instance.get_matrix();
    context.instance_without_translation = instance.get_matrix(true);
    context.instance_offset_z = instance.get_offset(Z);
    context.volume_matrix = source.get_matrix();
    context.model_parts = object.parts_count();
    for (const ModelVolume *volume : object.volumes) {
        FaceExtrusionVolumeState stamp;
        stamp.id = volume->id();
        stamp.type = volume->type();
        stamp.identity = volume->get_mesh_shared_ptr();
        stamp.matrix = volume->get_matrix();
        if (volume->is_model_part()) {
            const TriangleMesh &hull = volume->get_convex_hull();
            stamp.hull_identity = volume->get_convex_hull_shared_ptr();
            stamp.has_hull = !hull.empty();
            if (volume != &source)
                stamp.bed_mesh = stamp.has_hull ? hull : volume->mesh();
        }
        context.volumes.push_back(std::move(stamp));
    }
    return context;
}

struct FaceExtrusionPlacementBounds {
    BoundingBoxf3 printable_bounds;
    BoundingBoxf3 visual_bounds; // result stays visible even when it is a negative cutter
    double min_z = 0.0;
    double max_z = 0.0;
};

inline FaceExtrusionPlacementBounds face_extrusion_replacement_placement(
    const FaceExtrusionPlacementContext &context, const TriangleMesh &result_local,
    const TriangleMesh &result_hull, const std::function<void()> &check_cancel = {})
{
    using namespace FaceExtrusionModelDetail;
    FaceExtrusionPlacementBounds placement;
    placement.visual_bounds = result_local.transformed_bounding_box(context.world);
    const bool negative = context.source_type == ModelVolumeType::NEGATIVE_VOLUME;
    if (negative) {
        // Ignoring a cutter for bed placement must NOT ignore invalid cutter
        // geometry, a bad hull/pose, or an object without a positive carrier.
        if (context.model_parts == 0 || !usable_transform(context.world) ||
            !usable_mesh(result_local.its) || !usable_mesh(result_hull.its) || !finite_bounds(placement.visual_bounds))
            throw std::runtime_error("The negative result has invalid geometry, bounds or carrier.");
    } else if (context.source_type != ModelVolumeType::MODEL_PART) {
        throw std::runtime_error("Unsupported source role for face-extrusion placement.");
    } else {
        placement.printable_bounds = placement.visual_bounds;
    }
    double min_z = std::numeric_limits<double>::max();
    double max_z = -std::numeric_limits<double>::max();
    const auto merge_z = [&](const TriangleMesh &mesh, const Transform3d &matrix, bool include_max) {
        for (size_t i = 0; i < mesh.its.indices.size(); ++i) {
            if ((i & 1023) == 0 && check_cancel) check_cancel();
            const Vec3i &face = mesh.its.indices[i];
            for (int index : {face[0], face[1], face[2]}) {
                const double z = (matrix * mesh.its.vertices[index].cast<double>()).z();
                if (!std::isfinite(z)) throw std::runtime_error("Invalid post-operation height.");
                min_z = std::min(min_z, z);
                if (include_max) max_z = std::max(max_z, z);
            }
        }
    };
    // EXACT arithmetic/order used by ModelObject::get_instance_min/max_z:
    // instance without translation, volume matrix, hull, then instance Z offset.
    // A negative result remains a cutter, never printable material. Its original
    // positive carrier(s), not either the old or new cutter, determine placement.
    if (!negative)
        merge_z(result_hull, context.instance_without_translation * context.volume_matrix, true);
    for (const auto &volume : context.volumes) {
        if (volume.id == context.source_id || volume.type != ModelVolumeType::MODEL_PART)
            continue; // never let the OLD source or modifiers fake safety
        if (volume.bed_mesh.empty()) continue;
        if (!usable_mesh(volume.bed_mesh.its) || !usable_transform(volume.matrix))
            throw std::runtime_error("Another part has unsupported geometry or transformation for a safe refresh.");
        merge_z(volume.bed_mesh, context.instance_without_translation * volume.matrix, volume.has_hull);
        const BoundingBoxf3 bounds = volume.bed_mesh.transformed_bounding_box(context.instance_matrix * volume.matrix);
        if (!finite_bounds(bounds)) throw std::runtime_error("Another part has invalid bounds.");
        placement.printable_bounds.merge(bounds);
        placement.visual_bounds.merge(bounds);
    }
    placement.min_z = min_z + context.instance_offset_z;
    placement.max_z = max_z + context.instance_offset_z;
    if (!finite_bounds(placement.printable_bounds) || !std::isfinite(placement.min_z) || !std::isfinite(placement.max_z))
        throw std::runtime_error("The post-operation object has invalid bounds.");
    return placement;
}

inline bool face_extrusion_bed_safe(size_t model_parts, double min_z, double max_z)
{
    // EXACT ensure_on_bed(true) branches, including the one-part thin-lamina
    // case. No epsilon: even a tiny non-zero implicit translation is refused.
    if (model_parts == 0 || !std::isfinite(min_z) || !std::isfinite(max_z))
        return false;
    return model_parts == 1 ? !(min_z != 0.0 && (min_z >= SINKING_Z_THRESHOLD || max_z < 0.0)) :
                             max_z >= SINKING_MIN_Z_THRESHOLD;
}

inline bool face_extrusion_new_part_above_bed(ModelVolumeType source_type, double min_z)
{
    if (source_type == ModelVolumeType::NEGATIVE_VOLUME)
        return std::isfinite(min_z); // cutters may cross the bed; the carrier must not move
    return source_type == ModelVolumeType::MODEL_PART && !(min_z < -face_extrusion_bed_tolerance_mm);
}

inline bool face_extrusion_new_parts_bed_safe(
    ModelVolumeType source_type, const ModelObject &object,
    const ModelInstance &instance, double combined_max_z)
{
    if (source_type == ModelVolumeType::MODEL_PART)
        // Positive NewParts retain the existing multiple-part refresh branch.
        return std::isfinite(combined_max_z) && combined_max_z >= SINKING_MIN_Z_THRESHOLD;
    if (source_type != ModelVolumeType::NEGATIVE_VOLUME || object.parts_count() == 0 ||
        object.instances.size() != 1 || object.instances.front() != &instance ||
        !FaceExtrusionModelDetail::usable_transform(instance.get_matrix()))
        return false;
    // A new negative changes neither printable count nor printable bounds. Read
    // the current positive carrier, never cached preview bounds or a cutter's Z.
    // Validate the same hull/min-Z fallback that the native getters will use.
    for (const ModelVolume *volume : object.volumes) {
        if (!volume->is_model_part()) continue;
        const TriangleMesh &hull = volume->get_convex_hull();
        const TriangleMesh &bed_mesh = hull.empty() ? volume->mesh() : hull;
        if (!FaceExtrusionModelDetail::usable_transform(volume->get_matrix()) ||
            !FaceExtrusionModelDetail::usable_mesh(bed_mesh.its) ||
            !FaceExtrusionModelDetail::finite_bounds(bed_mesh.transformed_bounding_box(instance.get_matrix() * volume->get_matrix())))
            return false;
    }
    return face_extrusion_bed_safe(object.parts_count(), object.get_instance_min_z(0), object.get_instance_max_z(0));
}

inline ModelVolume *insert_face_extrusion_part(
    ModelObject &object, TriangleMesh &&object_mesh, const std::string &name,
    FaceExtrusionToolMode mode, ModelVolumeType source_type, int effective_extruder)
{
    if (mode == FaceExtrusionToolMode::CutAndStretch || !face_extrusion_source_role_allowed(mode, source_type))
        throw std::invalid_argument("This source role does not support new face-extrusion parts.");
    ModelVolume *created = object.add_volume(std::move(object_mesh), source_type, false);
    created->calculate_convex_hull();
    created->name = name;
    created->source.is_from_builtin_objects = true;
    // Negative cutters get no extruder/config/paint/material clone. Positive
    // outputs retain ONLY the effective/inherited extruder, exactly as before.
    if (source_type == ModelVolumeType::MODEL_PART)
        created->config.set_key_value("extruder", new ConfigOptionInt(effective_extruder));
    return created;
}

inline void replace_face_extrusion_source(
    ModelVolume &source, FaceExtrusionToolMode mode, ModelVolumeType expected_type,
    std::shared_ptr<const TriangleMesh> result,
    const std::shared_ptr<const TriangleMesh> &hull, bool discard_annotations)
{
    if (!face_extrusion_source_type_matches(mode, expected_type, source))
        throw std::invalid_argument("The face-extrusion source role changed before replacement.");
    source.set_mesh(result);
    if (discard_annotations)
        // The caller staged legacy-property COW on a private candidate and
        // authorized these source-only resets inside its successful Undo step.
        clear_face_extrusion_annotations(source);
    source.set_convex_hull_shared_ptr(hull);
    source.invalidate_convex_hull_2d();
    source.set_new_unique_id(); // cache identity, NOT a new logical part GUID
    source.source = ModelVolume::Source();
    source.source.is_from_builtin_objects = true;
    source.get_object()->invalidate_bounding_box();
}

} // namespace Slic3r::GUI

#endif // slic3r_FaceExtrusionModelPolicy_hpp_
