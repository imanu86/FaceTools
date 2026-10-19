// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#ifndef slic3r_ModelObjectLoadPlacement_hpp_
#define slic3r_ModelObjectLoadPlacement_hpp_

#include "libslic3r/Model.hpp"

namespace Slic3r::GUI {

enum class ModelObjectLoadPlacement { EnsureOnBed, PreserveExisting };

inline ModelObjectLoadPlacement cut_result_load_placement(bool cut_to_parts)
{
    return cut_to_parts ? ModelObjectLoadPlacement::PreserveExisting : ModelObjectLoadPlacement::EnsureOnBed;
}

// Imports retain their existing bed-placement policy. Native Cut to parts
// supplies already-positioned instances: do not drop or lift that object after
// cutting it. This changes only automatic Z placement, not geometry validation,
// scaling limits, part transforms, or the behavior of other import callers.
inline void apply_model_object_load_placement(ModelObject &object, bool allow_negative_z,
                                              ModelObjectLoadPlacement placement)
{
    if (placement != ModelObjectLoadPlacement::PreserveExisting)
        object.ensure_on_bed(allow_negative_z);
}

} // namespace Slic3r::GUI

#endif // slic3r_ModelObjectLoadPlacement_hpp_
