// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include "FFFStlExportSafety.hpp"
#include "libslic3r/MeshBoolean.hpp"
#include "libslic3r/MeshSolidValidation.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
#include <utility>

namespace Slic3r {
namespace {
constexpr double stl_transform_error_mm = 1e-4;

struct PreparedOperand {
    indexed_triangle_set mesh; // Actual object-space float storage consumed by CSG.
    std::vector<Vec3d> original_vertices;
    Transform3d volume_transform = Transform3d::Identity();
};

// The OMC adapter catches std::exceptions internally. Never let it swallow the
// caller's cancellation: record the original exception, request cancellation via
// its bool callback, and rethrow the original exception after the adapter returns.
class CancelContext {
public:
    explicit CancelContext(const std::function<void()> &fn) : m_fn(fn) {}
    void check()
    {
        rethrow_if_cancelled();
        if (m_fn) {
            try { m_fn(); }
            catch (...) { m_exception = std::current_exception(); throw; }
        }
    }
    bool backend_cancel() noexcept
    {
        try { check(); return false; }
        catch (...) { return true; }
    }
    void rethrow_if_cancelled() const
    {
        if (m_exception) std::rethrow_exception(m_exception);
    }
private:
    const std::function<void()> &m_fn;
    std::exception_ptr m_exception;
};

template<class Result, class Fn>
Result guarded(const std::function<void()> &check_cancel, Fn &&fn)
{
    CancelContext cancel(check_cancel);
    try {
        cancel.check();
        Result result = fn(cancel);
        cancel.check();
        return result;
    } catch (const std::exception &error) {
        cancel.rethrow_if_cancelled();
        Result result;
        result.message = error.what();
        return result;
    } catch (...) {
        cancel.rethrow_if_cancelled();
        Result result;
        result.message = "STL geometry preparation failed.";
        return result;
    }
}

void validate(const indexed_triangle_set &mesh, CancelContext &cancel)
{
    validate_solid_mesh(mesh, [&cancel] { cancel.check(); });
}

void append_checked(indexed_triangle_set &destination, const indexed_triangle_set &source)
{
    constexpr size_t limit = size_t(std::numeric_limits<int>::max());
    if (destination.vertices.size() > limit || source.vertices.size() > limit - destination.vertices.size() ||
        destination.indices.size() > limit || source.indices.size() > limit - destination.indices.size())
        throw std::runtime_error("STL composition exceeds indexed mesh limits.");
    its_merge(destination, source);
}

std::vector<Vec3d> original_vertices(const indexed_triangle_set &mesh, CancelContext &cancel)
{
    std::vector<Vec3d> result;
    result.reserve(mesh.vertices.size());
    for (size_t i = 0; i < mesh.vertices.size(); ++i) {
        if ((i & 1023) == 0) cancel.check();
        result.emplace_back(mesh.vertices[i].cast<double>());
    }
    return result;
}

std::vector<Vec3d> transformed_reference(const std::vector<Vec3d> &original,
    const Transform3d &cumulative, CancelContext &cancel)
{
    if (!cumulative.matrix().allFinite())
        throw std::runtime_error("STL cumulative transformation is not finite.");
    std::vector<Vec3d> reference;
    reference.reserve(original.size());
    for (size_t i = 0; i < original.size(); ++i) {
        if ((i & 1023) == 0) cancel.check();
        const Vec3d expected = cumulative * original[i];
        if (!expected.allFinite()) throw std::runtime_error("STL physical reference is not finite.");
        reference.push_back(expected);
    }
    return reference;
}

double transform_checked(indexed_triangle_set &mesh, const Transform3d &transform, CancelContext &cancel,
                         const std::vector<Vec3d> *physical_reference = nullptr)
{
    if (!transform.matrix().allFinite())
        throw std::runtime_error("STL export has a non-finite transformation.");
    const double determinant = transform.linear().determinant();
    if (!std::isfinite(determinant) || determinant == 0.0)
        throw std::runtime_error("STL export has a singular or unrepresentable transformation.");
    std::vector<Vec3d> local_reference;
    if (!physical_reference) {
        local_reference = transformed_reference(original_vertices(mesh, cancel), transform, cancel);
        physical_reference = &local_reference;
    }
    if (physical_reference->size() != mesh.vertices.size())
        throw std::runtime_error("STL transform reference no longer matches vertex order.");
    cancel.check();
    // Vertex order stays unchanged (reflection reverses face winding only).
    // Physical error is measured against DOUBLE references, not cached bounds,
    // topology alone or the already-rounded input of a previous transform.
    its_transform(mesh, transform, true);
    double maximum_error = 0.0;
    for (size_t i = 0; i < mesh.vertices.size(); ++i) {
        if ((i & 1023) == 0) cancel.check();
        const Vec3d &expected = (*physical_reference)[i];
        const double error = (mesh.vertices[i].cast<double>() - expected).norm();
        if (!expected.allFinite() || !mesh.vertices[i].allFinite() ||
            !std::isfinite(error) || error > stl_transform_error_mm)
            throw std::runtime_error("STL transformation loses more than 0.0001 mm in float storage.");
        maximum_error = std::max(maximum_error, error);
    }
    validate(mesh, cancel);
    return maximum_error;
}

// CSG creates new/reordered vertices, so no source-vertex correspondence is
// invented for its result. Conservatively account for amplified operand storage
// error separately from the final result-transform rounding. This is NOT a proof
// of global Boolean geometric fidelity (the real checked adapter remains required).
double operand_world_error(const PreparedOperand &operand, const Transform3d &instance, CancelContext &cancel)
{
    const auto reference = transformed_reference(operand.original_vertices,
                                                  instance * operand.volume_transform, cancel);
    double maximum_error = 0.0;
    for (size_t i = 0; i < operand.mesh.vertices.size(); ++i) {
        if ((i & 1023) == 0) cancel.check();
        const Vec3d actual = instance * operand.mesh.vertices[i].cast<double>();
        const double error = (actual - reference[i]).norm();
        if (!actual.allFinite() || !std::isfinite(error) || error > stl_transform_error_mm)
            throw std::runtime_error("STL instance amplifies operand float-storage error beyond 0.0001 mm.");
        maximum_error = std::max(maximum_error, error);
    }
    return maximum_error;
}

StlGeometryResult result_of(indexed_triangle_set mesh, CancelContext &cancel)
{
    StlGeometryResult result;
    if (mesh.indices.empty() && mesh.vertices.empty()) {
        result.status = StlGeometryStatus::Empty;
        return result;
    }
    validate(mesh, cancel); // BEFORE TriangleMesh computes any cached statistics.
    result.mesh = TriangleMesh(std::move(mesh));
    result.status = StlGeometryStatus::Success;
    return result;
}

void validate_result(const StlGeometryResult &result, CancelContext &cancel)
{
    if (result.failed())
        throw std::runtime_error(result.message.empty() ? "An STL contributor failed." : result.message);
    if (result.status == StlGeometryStatus::Empty) {
        if (!result.mesh.its.indices.empty() || !result.mesh.its.vertices.empty())
            throw std::runtime_error("An empty STL contributor carries unexpected geometry.");
    } else if (result.status == StlGeometryStatus::Success) {
        validate(result.mesh.its, cancel);
    } else {
        throw std::runtime_error("Invalid STL contributor status.");
    }
}

indexed_triangle_set checked_boolean(const indexed_triangle_set &a, const indexed_triangle_set &b,
                                    const char *operation, CancelContext &cancel)
{
    if (!MeshBoolean::openmeshcraft::is_enabled())
        throw std::runtime_error("Strict Boolean STL export requires the checked OpenMeshCraft adapter.");
    if (!a.indices.empty() || !a.vertices.empty()) validate(a, cancel);
    if (!b.indices.empty() || !b.vertices.empty()) validate(b, cancel);
    // Construct only after validation. The empty set is an adapter identity, not
    // an input solid sent to B20's non-empty validator.
    TriangleMesh left = a.indices.empty() ? TriangleMesh{} : TriangleMesh(a);
    TriangleMesh right = b.indices.empty() ? TriangleMesh{} : TriangleMesh(b);
    std::vector<TriangleMesh> output;
    bool failed = false;
    bool completed = false;
    MeshBoolean::mcut::make_boolean(left, right, output, operation,
        [&cancel] { return cancel.backend_cancel(); },
        [&completed](float progress) { if (progress == 100.0f) completed = true; },
        [&failed] { failed = true; });
    cancel.rethrow_if_cancelled();
    cancel.check();
    if (failed || !completed)
        throw std::runtime_error("Boolean STL export did not complete successfully.");
    // Current enabled OMC adapter publishes 0 meshes for success-empty, 1 for a
    // non-empty result. Any reported failure discards ALL output, even if nonempty.
    indexed_triangle_set mesh;
    for (const auto &part : output) {
        validate(part.its, cancel);
        append_checked(mesh, part.its);
    }
    if (!mesh.indices.empty()) validate(mesh, cancel);
    return mesh;
}
} // namespace

StlGeometryResult prepare_fff_stl_object(const ModelObject &object, int instance_scope,
    FffStlExportMode mode, const std::function<void()> &check_cancel)
{
    return guarded<StlGeometryResult>(check_cancel, [&](CancelContext &cancel) {
        if (instance_scope < -1 || (instance_scope >= 0 && size_t(instance_scope) >= object.instances.size()))
            throw std::runtime_error("Invalid STL instance scope.");
        if (mode != FffStlExportMode::RawPositives && mode != FffStlExportMode::BooleanPositivesMinusNegatives)
            throw std::runtime_error("Invalid FFF STL export mode.");
        const bool boolean_mode = mode == FffStlExportMode::BooleanPositivesMinusNegatives;
        if (boolean_mode && !MeshBoolean::openmeshcraft::is_enabled())
            throw std::runtime_error("Strict Boolean STL export requires OpenMeshCraft.");
        std::vector<PreparedOperand> positives, negatives;
        for (const ModelVolume *volume : object.volumes) {
            cancel.check();
            if (!volume) throw std::runtime_error("Null STL source volume.");
            if (!volume->is_model_part() && !(boolean_mode && volume->is_negative_volume())) continue;
            if (!volume->mesh_ptr()) throw std::runtime_error("STL source volume has no mesh.");
            PreparedOperand operand;
            operand.mesh = volume->mesh().its;
            validate(operand.mesh, cancel);
            operand.original_vertices = original_vertices(operand.mesh, cancel);
            operand.volume_transform = volume->get_matrix();
            const auto reference = transformed_reference(operand.original_vertices, operand.volume_transform, cancel);
            transform_checked(operand.mesh, operand.volume_transform, cancel, &reference);
            (volume->is_model_part() ? positives : negatives).push_back(std::move(operand));
        }
        // Validate ALL contributing operands before the first Boolean call. Model
        // volume ordering is not changed; this explicit export policy uses copies.
        indexed_triangle_set object_mesh;
        for (const auto &part : positives) {
            cancel.check();
            if (boolean_mode) object_mesh = checked_boolean(object_mesh, part.mesh, "UNION", cancel);
            else append_checked(object_mesh, part.mesh);
        }
        for (const auto &cutter : negatives) {
            cancel.check();
            object_mesh = checked_boolean(object_mesh, cutter.mesh, "A_NOT_B", cancel);
        }
        if (!object_mesh.indices.empty()) validate(object_mesh, cancel);

        indexed_triangle_set composed;
        const auto add_instance = [&](const ModelInstance *instance) {
            if (!instance) throw std::runtime_error("Null STL source instance.");
            // Empty mathematical objects contribute nothing, but bad instance
            // metadata must not be hidden by that identity.
            const Transform3d transform = instance->get_matrix();
            if (!transform.matrix().allFinite())
                throw std::runtime_error("Non-finite STL instance transformation.");
            const double determinant = transform.linear().determinant();
            if (!std::isfinite(determinant) || determinant == 0.0)
                throw std::runtime_error("Invalid STL instance transformation.");
            double input_world_error = 0.0;
            if (boolean_mode) {
                for (const auto &operand : positives)
                    input_world_error = std::max(input_world_error, operand_world_error(operand, transform, cancel));
                for (const auto &operand : negatives)
                    input_world_error = std::max(input_world_error, operand_world_error(operand, transform, cancel));
            }
            if (object_mesh.indices.empty()) return;
            indexed_triangle_set copy = object_mesh;
            std::vector<Vec3d> physical_reference;
            if (!boolean_mode) {
                physical_reference.reserve(copy.vertices.size());
                for (const auto &operand : positives) {
                    const auto reference = transformed_reference(operand.original_vertices,
                                                                  transform * operand.volume_transform, cancel);
                    physical_reference.insert(physical_reference.end(), reference.begin(), reference.end());
                }
            }
            const double result_rounding = transform_checked(copy, transform, cancel,
                boolean_mode ? nullptr : &physical_reference);
            if (boolean_mode && input_world_error + result_rounding > stl_transform_error_mm)
                throw std::runtime_error("STL operand and final-transform storage errors exceed 0.0001 mm.");
            append_checked(composed, copy);
        };
        if (instance_scope == -1) {
            for (const ModelInstance *instance : object.instances) {
                cancel.check();
                add_instance(instance);
            }
        } else {
            add_instance(object.instances[size_t(instance_scope)]);
        }
        // Instances are deliberately UNFUSED. Contact/overlap here fails even in
        // Boolean mode, which authorizes volume-level operations only.
        return result_of(std::move(composed), cancel);
    });
}

StlGeometryResult prepare_stl_volume(const indexed_triangle_set &source,
    const std::vector<Transform3d> &transforms, const std::function<void()> &check_cancel)
{
    return guarded<StlGeometryResult>(check_cancel, [&](CancelContext &cancel) {
        indexed_triangle_set copy = source;
        validate(copy, cancel);
        const auto original = original_vertices(copy, cancel);
        Transform3d cumulative = Transform3d::Identity();
        for (const auto &transform : transforms) {
            cumulative = transform * cumulative;
            const auto reference = transformed_reference(original, cumulative, cancel);
            transform_checked(copy, transform, cancel, &reference);
        }
        return result_of(std::move(copy), cancel);
    });
}

StlGeometryResult combine_stl_export_results(const std::vector<StlGeometryResult> &contributors,
    const std::function<void()> &check_cancel)
{
    return guarded<StlGeometryResult>(check_cancel, [&](CancelContext &cancel) {
        indexed_triangle_set combined;
        for (const auto &contributor : contributors) {
            cancel.check();
            validate_result(contributor, cancel);
            if (contributor.status == StlGeometryStatus::Success)
                append_checked(combined, contributor.mesh.its);
        }
        return result_of(std::move(combined), cancel);
    });
}

StlExportBatchResult finalize_stl_export_batch(std::vector<StlGeometryResult> final_results,
    const std::function<void()> &check_cancel)
{
    return guarded<StlExportBatchResult>(check_cancel, [&](CancelContext &cancel) {
        StlExportBatchResult result;
        for (size_t i = 0; i < final_results.size(); ++i) {
            cancel.check();
            validate_result(final_results[i], cancel);
            if (final_results[i].status == StlGeometryStatus::Success)
                result.outputs.push_back({i, std::move(final_results[i].mesh)});
        }
        if (result.outputs.empty()) throw std::runtime_error("There are no non-empty solid STL outputs.");
        result.ready = true;
        return result;
    });
}
} // namespace Slic3r
