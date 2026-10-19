// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#ifndef slic3r_FFFStlExportSafety_hpp_
#define slic3r_FFFStlExportSafety_hpp_

#include "libslic3r/Model.hpp"
#include "libslic3r/TriangleMesh.hpp"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace Slic3r {

// Export-only contract. Does not change Model, Undo, selection or legacy Boolean APIs.
// Boolean mode means UNION(all MODEL_PART volumes) MINUS(all NEGATIVE_VOLUME volumes),
// regardless of their order in ModelObject::volumes. No implicit cross-instance/object union.
enum class FffStlExportMode { RawPositives, BooleanPositivesMinusNegatives };
enum class StlGeometryStatus { Success, Empty, Failure };

struct StlGeometryResult {
    StlGeometryStatus status = StlGeometryStatus::Failure;
    TriangleMesh mesh;
    std::string message;
    bool failed() const noexcept { return status == StlGeometryStatus::Failure; }
};

// instance_scope == -1 selects all instances; otherwise it must be a valid index.
// A mathematically empty Boolean result is Empty, NOT Failure and NOT raw positives.
// Any exception from check_cancel propagates unchanged. Other errors become Failure.
StlGeometryResult prepare_fff_stl_object(
    const ModelObject &object, int instance_scope, FffStlExportMode mode,
    const std::function<void()> &check_cancel = {});

// Existing explicit volume-selection route only. Apply the exact sequence of
// consumed Transform3d -> float-storage transforms, validating after each one.
StlGeometryResult prepare_stl_volume(
    const indexed_triangle_set &mesh, const std::vector<Transform3d> &transforms,
    const std::function<void()> &check_cancel = {});

// Empty contributors add nothing; ANY Failure aborts the complete result.
// Components are concatenated, NEVER implicitly united or repaired.
StlGeometryResult combine_stl_export_results(
    const std::vector<StlGeometryResult> &contributors,
    const std::function<void()> &check_cancel = {});

struct PreparedStlOutput {
    size_t request_index = 0; // Keep caller's name/selection mapping when Empty is skipped.
    TriangleMesh mesh;
};
struct StlExportBatchResult {
    bool ready = false;
    std::vector<PreparedStlOutput> outputs;
    std::string message;
};

// Call AFTER the caller's final align_to_origin/origin_translation and ALL merges.
// Geometry only: does not reserve filenames or write files. Empty outputs are skipped;
// an all-empty batch is not writable. Failure publishes no partial output vector.
StlExportBatchResult finalize_stl_export_batch(
    std::vector<StlGeometryResult> final_results,
    const std::function<void()> &check_cancel = {});

} // namespace Slic3r
#endif
