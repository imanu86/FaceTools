// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#ifndef slic3r_FaceToolsPreviewBounds_hpp_
#define slic3r_FaceToolsPreviewBounds_hpp_

#include "libslic3r/BoundingBox.hpp"

namespace Slic3r::GUI {

// Translation only: the cut-plane vertices already contain their orientation.
// Preserve GLModel's defined planar box. BoundingBoxf3::transformed reconstructs
// it with a positive-volume constructor, losing defined for an axis-aligned plane.
// Do not pad geometry or repair invalid input; callers retain their finite/float
// precision guards. Undefined boxes stay undefined and nonfinite values propagate.
inline BoundingBoxf3 translated_face_tools_preview_bounds(const BoundingBoxf3 &local_bounds,
                                                          const Vec3d &translation)
{
    BoundingBoxf3 bounds = local_bounds;
    bounds.translate(translation);
    return bounds;
}

} // namespace Slic3r::GUI

#endif // slic3r_FaceToolsPreviewBounds_hpp_
