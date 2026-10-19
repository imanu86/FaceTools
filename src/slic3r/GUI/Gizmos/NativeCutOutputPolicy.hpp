// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#ifndef slic3r_NativeCutOutputPolicy_hpp_
#define slic3r_NativeCutOutputPolicy_hpp_

#include "libslic3r/CutUtils.hpp"
#include "slic3r/GUI/ModelObjectLoadPlacement.hpp"

namespace Slic3r::GUI {

// Native Cut output intent, independent of the wx/GL lifetime. Context changes
// affect whether the request can run, never silently select a different output.
struct NativeCutOutputSettings {
    bool keep_upper{true};
    bool keep_lower{true};
    bool cut_to_parts{true};
    bool place_on_cut_upper{false};
    bool place_on_cut_lower{false};
    bool flip_upper{false};
    bool flip_lower{false};

    void reset() { *this = NativeCutOutputSettings{}; }
    bool keeps_any_piece() const { return keep_upper || keep_lower; }

    // These are the same transition points used by the gizmo and native state
    // regressions. The request survives entering/leaving an unsupported mode or
    // temporarily assigning a part to the other side. Only explicit output UI
    // changes, Reset, or archive restoration may change that request.
    void on_mode_changed(CutMode) const {}
    void on_part_assignment_changed(bool) const {}

    // Normalize a presentation copy. The action resolver independently applies
    // this to its own copy, so a missing render frame cannot change the cut.
    void normalize_for_render()
    {
        if (cut_to_parts) {
            keep_upper = keep_lower = true;
            place_on_cut_upper = place_on_cut_lower = false;
            flip_upper = flip_lower = false;
        }
    }

    template<class Archive> void serialize(Archive &ar)
    {
        // Retain the native four-field prefix and include every output choice.
        // This is the gizmo's in-memory Undo state, not a project-file format.
        ar(keep_upper, keep_lower, flip_lower, flip_upper,
           cut_to_parts, place_on_cut_upper, place_on_cut_lower);
    }
};

enum class NativeCutRoute { Plane, Contour, Groove };
enum class NativeCutOutputBlock {
    None,
    KeepNeither,
    InPlaceWithConnectors,
    InPlaceWithGroove,
    InPlaceWithReassignment
};

struct NativeCutOutputContext {
    CutMode mode{CutMode::cutPlanar};
    bool has_connectors{false};
    bool has_valid_modified_parts{false};
    int dowels_count{0};
    bool cut_id_invalid{false};
};

struct NativeCutOutputDecision {
    NativeCutRoute route{NativeCutRoute::Plane};
    ModelObjectCutAttributes attributes;
    NativeCutOutputBlock block{NativeCutOutputBlock::None};
    bool separate_half_offset{true};

    bool allowed() const { return block == NativeCutOutputBlock::None; }
    bool cut_to_parts() const
    {
        // enum_bitmask::has is not const in the native API.
        auto copy = attributes;
        return copy.has(ModelObjectCutAttribute::CutToParts);
    }
    ModelObjectLoadPlacement load_placement() const
    {
        return cut_result_load_placement(cut_to_parts());
    }
};

inline NativeCutOutputDecision resolve_native_cut_output(
    const NativeCutOutputSettings &settings, const NativeCutOutputContext &context)
{
    NativeCutOutputDecision decision;
    const bool groove = context.mode == CutMode::cutTongueAndGroove;
    decision.route = groove ? NativeCutRoute::Groove :
        context.has_valid_modified_parts ? NativeCutRoute::Contour : NativeCutRoute::Plane;

    if (settings.cut_to_parts) {
        decision.separate_half_offset = false;
        if (context.has_connectors)
            decision.block = NativeCutOutputBlock::InPlaceWithConnectors;
        else if (groove)
            decision.block = NativeCutOutputBlock::InPlaceWithGroove;
        else if (context.has_valid_modified_parts)
            decision.block = NativeCutOutputBlock::InPlaceWithReassignment;
        // Those native consumers do not implement same-object placement. Do not
        // silently dispatch to one with a weakened or misleading Parts flag.
        if (!decision.allowed())
            return decision;
    }

    NativeCutOutputSettings effective = settings;
    effective.normalize_for_render();
    if (!effective.keeps_any_piece()) {
        decision.block = NativeCutOutputBlock::KeepNeither;
        return decision;
    }

    // Explicit separate-object output retains the existing native flags and
    // separation behavior. In-place output is canonical regardless of raw UI
    // keep/orientation values or whether the render normalization ever ran.
    decision.attributes =
        only_if(context.has_connectors ? true : effective.keep_upper, ModelObjectCutAttribute::KeepUpper) |
        only_if(context.has_connectors ? true : effective.keep_lower, ModelObjectCutAttribute::KeepLower) |
        only_if(effective.cut_to_parts, ModelObjectCutAttribute::CutToParts) |
        only_if(effective.place_on_cut_upper, ModelObjectCutAttribute::PlaceOnCutUpper) |
        only_if(effective.place_on_cut_lower, ModelObjectCutAttribute::PlaceOnCutLower) |
        only_if(effective.flip_upper, ModelObjectCutAttribute::FlipUpper) |
        only_if(effective.flip_lower, ModelObjectCutAttribute::FlipLower) |
        only_if(context.dowels_count > 0, ModelObjectCutAttribute::CreateDowels) |
        only_if(!context.has_connectors && !groove && context.cut_id_invalid, ModelObjectCutAttribute::InvalidateCutInfo);
    return decision;
}

} // namespace Slic3r::GUI

#endif // slic3r_NativeCutOutputPolicy_hpp_
