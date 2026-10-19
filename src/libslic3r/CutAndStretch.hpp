// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#ifndef slic3r_CutAndStretch_hpp_
#define slic3r_CutAndStretch_hpp_

#include "TriangleMesh.hpp"

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace Slic3r {

struct CutAndStretchSourceFrame; // opaque, immutable source-local provenance
struct CutAndStretchGuidedState; // opaque, immutable excluded upper operands

enum class CutAndStretchScope { WholeHalfSpace, GuidedComponent };

struct CutAndStretchSelectionSummary {
    size_t upper_component_count = 0;
    size_t selected_upper_component = 0; // internal zero-based index; NOT stable across preparations
    size_t guide_facet_count = 0;
    size_t selected_source_facet_count = 0; // surviving original provenance IDs; excludes generated caps
};

// Immutable prepared operands, in one physical-mm cut frame. The cut plane is
// z=0; only the moving operand is translated. bridge_unit occupies z=0..1
// within the validated float rounding budget (no cap vertices are flattened).
// The preview may scale the bridge along Z, but must NEVER scale moving.
// In GuidedComponent scope stationary is the whole fixed remainder, and may
// extend above z=0. The displayed moving mesh is the ACTUAL selected component,
// not merely the guide mask. The GUI must separately confirm this immutable plan.
struct CutAndStretchPlan {
    indexed_triangle_set stationary;
    indexed_triangle_set moving;
    indexed_triangle_set bridge_unit;
    Transform3d cut_to_world = Transform3d::Identity();
    double contact_area_mm2 = 0.0;
    double source_volume_mm3 = 0.0;
    std::string message;
    bool valid = false;
    std::shared_ptr<const CutAndStretchSourceFrame> source_frame; // set only by placement-aware preparation
    CutAndStretchScope scope = CutAndStretchScope::WholeHalfSpace;
    CutAndStretchSelectionSummary selection;
    std::shared_ptr<const CutAndStretchGuidedState> guided_state;
    explicit operator bool() const noexcept { return valid; }
};

struct CutAndStretchResult {
    indexed_triangle_set mesh; // WORLD mm for old apply; ORIGINAL SOURCE LOCAL for apply_to_source.
    std::string message;
    bool valid = false;
    explicit operator bool() const noexcept { return valid; }
};

// LEGACY WHOLE-HALF-SPACE operation. Source is not modified. A conservative
// contact-section policy rejects coplanar shelves which are not covered by the
// moving cap, disconnected cuts, malformed solids, and geometry inside the
// native cutter's uncertainty band. These APIs NEVER infer a local guided target.
// A failed/canceled result is invalid and contains no partial operands.
CutAndStretchPlan prepare_cut_and_stretch(
    const indexed_triangle_set &source_world,
    const Vec3d &plane_origin, const Vec3d &moving_normal,
    const std::function<bool()> &cancel = {},
    const std::function<void(int)> &progress = {});

// Placement-aware legacy path: convert directly from the unchanged local mesh
// into a near-source physical-mm cut frame. Never first narrow world coordinates
// to float. cut_to_world remains the actual double world pose for preview.
CutAndStretchPlan prepare_cut_and_stretch(
    const indexed_triangle_set &source_local, const Transform3d &source_to_world,
    const Vec3d &plane_origin_world, const Vec3d &moving_normal_world,
    const std::function<bool()> &cancel = {},
    const std::function<void(int)> &progress = {});

// Explicit face provenance: the ORIGINAL source triangle defines the plane
// exactly. The supplied world normal selects its moving side and must be parallel
// to that transformed facet. No near-plane inference or vertex projection.
CutAndStretchPlan prepare_cut_and_stretch_from_face(
    const indexed_triangle_set &source_local, const Transform3d &source_to_world,
    size_t plane_facet, const Vec3d &moving_normal_world,
    const std::function<bool()> &cancel = {},
    const std::function<void(int)> &progress = {});

// GUIDED COMPONENT operation, distinct from the legacy whole-half-space APIs.
// Pass ALL original triangle IDs in the selected guide masks, not just their
// clicked seeds. IDs must be nonempty, unique and in range. Every ID must map to
// the SAME closed upper component; a straddling facet is allowed only when all
// its upper fragments have that one owner. A wholly fixed/unassociated guide
// rejects. An exactly coplanar original guide may resolve through its original
// coplanar patch's adjacency; no nearest-point or majority fallback is used.
// The entire identified component moves, NOT an arbitrary CAD sub-feature.
// Source multishell/cavity surfaces, disconnected stationary halves/contact
// sections and partial lower-cap triangles remain explicitly unsupported.
// The fixed remainder is reconstructed by real CSG union, so source geometry is
// preserved within the existing physical-mm budget, not by identical triangles.
// No Model state, sibling volume or placement is modified by these CPU APIs.
CutAndStretchPlan prepare_cut_and_stretch_guided(
    const indexed_triangle_set &source_local, const Transform3d &source_to_world,
    const Vec3d &plane_origin_world, const Vec3d &moving_normal_world,
    const std::vector<int> &guide_facets,
    const std::function<bool()> &cancel = {},
    const std::function<void(int)> &progress = {});

CutAndStretchPlan prepare_cut_and_stretch_guided_from_face(
    const indexed_triangle_set &source_local, const Transform3d &source_to_world,
    size_t plane_facet, const Vec3d &moving_normal_world,
    const std::vector<int> &guide_facets,
    const std::function<bool()> &cancel = {},
    const std::function<void(int)> &progress = {});

// Positive extension only. Zero is a non-applicable preview/no-op, not a model
// edit; negative distances require a different remove-band operation and reject.
// Uses actual solid Boolean union, never concatenation or overlap-epsilon fixes.
// In guided scope any collision/touch with an excluded upper component rejects.
// Failure/cancellation leaves result.mesh empty and never modifies the plan.
CutAndStretchResult apply_cut_and_stretch(
    const CutAndStretchPlan &plan, double distance_mm,
    const std::function<bool()> &cancel = {},
    const std::function<void(int)> &progress = {});

// Requires a plan from a placement-aware API. Returns vertices directly in the
// original source-local frame, validating float storage error in physical
// millimeters. Do NOT inverse-transform this result through world coordinates.
CutAndStretchResult apply_cut_and_stretch_to_source(
    const CutAndStretchPlan &plan, double distance_mm,
    const std::function<bool()> &cancel = {},
    const std::function<void(int)> &progress = {});

} // namespace Slic3r

#endif // slic3r_CutAndStretch_hpp_
