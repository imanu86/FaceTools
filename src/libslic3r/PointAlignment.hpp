// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#ifndef slic3r_PointAlignment_hpp_
#define slic3r_PointAlignment_hpp_

#include "Model.hpp"
#include <Eigen/LU>
#include <cmath>
#include <limits>
#include <string>

namespace Slic3r {

// A picked point belongs to a mesh AND its original pose, not a transient GUI index.
struct PointAlignmentAnchor {
    ObjectID object_id, instance_id, volume_id;
    std::shared_ptr<const TriangleMesh> mesh;
    Transform3d instance_matrix = Transform3d::Identity();
    Transform3d volume_matrix = Transform3d::Identity();
    Vec3d local_point = Vec3d::Zero();
    Vec3d world_point = Vec3d::Zero();
};

struct PointAlignmentPlan {
    Transform3d moving_matrix = Transform3d::Identity();
    Vec3d world_displacement = Vec3d::Zero();
    Vec3d moving_world_point = Vec3d::Zero();
    bool valid = false;
    bool no_op = false;
    std::string message;
    explicit operator bool() const noexcept { return valid; }
};

inline bool point_alignment_same_matrix(const Transform3d &a, const Transform3d &b)
{
    return (a.matrix().array() == b.matrix().array()).all();
}

inline bool point_alignment_valid_transform(const Transform3d &matrix)
{
    if (!matrix.matrix().allFinite() || !matrix.matrix().cast<float>().allFinite() ||
        matrix.matrix().row(3) != Eigen::RowVector4d(0, 0, 0, 1))
        return false;
    const Matrix3d linear = matrix.linear();
    if (!linear.fullPivLu().isInvertible())
        return false;
    const Matrix3d inverse = linear.inverse();
    return inverse.allFinite() && linear.norm() * inverse.norm() <= 1e12;
}

inline bool point_alignment_plane_matches(int original_id, int selected_id, const Vec3d &normal,
                                          const Vec3d &plane_point, const Vec3d &local_point)
{
    const double norm = normal.norm();
    if (original_id < 0 || original_id != selected_id || !normal.allFinite() || !plane_point.allFinite() ||
        !local_point.allFinite() || !std::isfinite(norm) || norm == 0.0)
        return false;
    // Native ray hits are float. Allow a bounded float-rounding budget, not the
    // 0.05 mm raised highlight or an arbitrary projection onto another plane.
    const double scale = std::max(1.0, std::max(local_point.cwiseAbs().maxCoeff(), plane_point.cwiseAbs().maxCoeff()));
    const double tolerance = std::min(1e-4, std::max(1e-6, 8.0 * std::numeric_limits<float>::epsilon() * scale));
    return std::abs((normal / norm).dot(local_point - plane_point)) <= tolerance;
}

inline bool point_alignment_has_exterior(const Model &model)
{
    // This dedicated field is NOT in native Undo serialization. Legacy mesh tags
    // ARE serialized and must not trigger this guard. Check modifiers/neighbors too.
    for (const ModelObject *object : model.objects)
        if (object)
            for (const ModelVolume *volume : object->volumes)
                if (volume) {
                    const auto &data = volume->exterior_facets.get_data();
                    if (!data.first.empty() || !data.second.empty())
                        return true;
                }
    return false;
}

inline PointAlignmentAnchor point_alignment_anchor(const ModelObject &object, const ModelInstance &instance,
                                                    const ModelVolume &volume, const Vec3d &local_point)
{
    return {object.id(), instance.id(), volume.id(), volume.get_mesh_shared_ptr(),
            instance.get_matrix(), volume.get_matrix(), local_point,
            instance.get_matrix() * volume.get_matrix() * local_point};
}

// Read-only plan for the Prepare Volume/do_move path (FFF only at the GUI gate).
// desired_gap is Point 2 - Point 1 in WORLD millimetres. Zero means coincidence.
// No metadata, mesh, instance, selection, or Undo state is changed here.
inline PointAlignmentPlan prepare_point_alignment(const Model &model, const PointAlignmentAnchor &fixed,
                                                   const PointAlignmentAnchor &moving,
                                                   const Vec3d &desired_gap = Vec3d::Zero())
{
    PointAlignmentPlan result;
    const auto fail = [&](const char *message) { result.message = message; return result; };
    if (!desired_gap.allFinite() || !fixed.local_point.allFinite() || !moving.local_point.allFinite() ||
        !fixed.world_point.allFinite() || !moving.world_point.allFinite())
        return fail("Point coordinates and offsets must be finite.");
    if (fixed.object_id != moving.object_id || fixed.instance_id != moving.instance_id ||
        fixed.volume_id == moving.volume_id)
        return fail("Choose two different parts of the same object and instance.");
    const ModelObject *object = nullptr;
    for (const ModelObject *candidate : model.objects)
        if (candidate && candidate->id() == fixed.object_id) { object = candidate; break; }
    if (!object || object->instances.size() != 1 || !object->instances.front() || object->instances.front()->id() != fixed.instance_id)
        return fail("Use one instance of this object; moving a shared part would also move its copies.");
    const ModelVolume *fixed_volume = nullptr, *moving_volume = nullptr;
    for (const ModelVolume *volume : object->volumes)
        if (volume) {
            if (volume->id() == fixed.volume_id) fixed_volume = volume;
            if (volume->id() == moving.volume_id) moving_volume = volume;
        }
    if (!fixed_volume || !moving_volume || !fixed_volume->is_model_part() || !moving_volume->is_model_part())
        return fail("The picked parts changed. Pick the two points again.");
    const ModelInstance &instance = *object->instances.front();
    const auto current = [&](const PointAlignmentAnchor &anchor, const ModelVolume &volume) {
        return anchor.mesh && volume.get_mesh_shared_ptr() == anchor.mesh &&
               point_alignment_same_matrix(instance.get_matrix(), anchor.instance_matrix) &&
               point_alignment_same_matrix(volume.get_matrix(), anchor.volume_matrix) &&
               point_alignment_valid_transform(anchor.instance_matrix) && point_alignment_valid_transform(anchor.volume_matrix) &&
               (anchor.instance_matrix * anchor.volume_matrix * anchor.local_point - anchor.world_point).norm() <= 1e-8;
    };
    if (!current(fixed, *fixed_volume) || !current(moving, *moving_volume))
        return fail("The mesh, pose or picked point changed. Pick the two points again.");
    if (point_alignment_has_exterior(model))
        return fail("This project contains exterior-facet data that Undo cannot restore. Point moves are disabled; existing data is unchanged.");

    result.world_displacement = fixed.world_point + desired_gap - moving.world_point;
    result.moving_matrix = moving.volume_matrix;
    // Equivalent to V * T((I*V).linear^-1 * world_delta), without multiplying
    // the linear block at all: rotation, scale, shear and mirror remain bit-exact.
    result.moving_matrix.translation() += instance.get_matrix().linear().inverse() * result.world_displacement;
    result.moving_world_point = instance.get_matrix() * result.moving_matrix * moving.local_point;
    if (!result.moving_matrix.matrix().allFinite() || !result.moving_matrix.matrix().cast<float>().allFinite() ||
        !result.moving_world_point.allFinite() || !result.moving_world_point.cast<float>().allFinite() ||
        (result.moving_world_point - fixed.world_point - desired_gap).norm() > 1e-8)
        return fail("This translation loses coordinate precision; nothing was moved.");
    if (result.world_displacement.norm() <= 1e-9 || point_alignment_same_matrix(result.moving_matrix, moving.volume_matrix)) {
        result.valid = result.no_op = true;
        result.moving_matrix = moving.volume_matrix;
        result.moving_world_point = moving.world_point;
        return result;
    }

    // do_move traverses the entire canvas and can drop unrelated flying instances
    // to the bed. Match ModelObject::get_instance_min_z's calculation BEFORE Undo,
    // substituting only the candidate moving volume matrix. Never silently rebase.
    for (const ModelObject *other : model.objects) {
        if (!other) return fail("The project contains an unavailable object.");
        for (const ModelInstance *other_instance : other->instances) {
            if (!other_instance || !point_alignment_valid_transform(other_instance->get_matrix()))
                return fail("The project contains an unsupported instance transformation.");
            double min_z = std::numeric_limits<double>::max();
            for (const ModelVolume *volume : other->volumes) {
                if (!volume) return fail("The project contains an unavailable part.");
                if (!volume->is_model_part()) continue;
                const auto &hull = volume->get_convex_hull_shared_ptr();
                if (!hull || !volume->mesh_ptr())
                    return fail("A part has no validated mesh/hull for the bed preflight.");
                const auto &mesh = hull->its.indices.empty() ? volume->mesh().its : hull->its;
                const Transform3d &volume_matrix = volume == moving_volume ? result.moving_matrix : volume->get_matrix();
                if (!point_alignment_valid_transform(volume_matrix))
                    return fail("The project contains an unsupported part transformation.");
                const Transform3d transform = other_instance->get_matrix(true) * volume_matrix;
                for (const Vec3i &facet : mesh.indices)
                    for (int corner = 0; corner < 3; ++corner) {
                        const int index = facet[corner];
                        if (index < 0 || size_t(index) >= mesh.vertices.size() || !mesh.vertices[index].allFinite())
                            return fail("A part has invalid geometry for the bed preflight.");
                        const double z = (transform * mesh.vertices[index].cast<double>()).z();
                        if (!std::isfinite(z)) return fail("A part exceeds coordinate precision limits.");
                        min_z = std::min(min_z, z);
                    }
            }
            if (min_z == std::numeric_limits<double>::max()) min_z = 0.0;
            min_z += other_instance->get_offset(Z);
            if (!std::isfinite(min_z) || (min_z > SINKING_Z_THRESHOLD && std::abs(min_z) > EPSILON))
                return fail("The native move would also drop an object to the bed. Keep a stationary part on the bed before aligning points.");
        }
    }
    result.valid = true;
    return result;
}

} // namespace Slic3r
#endif
