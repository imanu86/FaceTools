// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#ifndef slic3r_PlanarProfile_hpp_
#define slic3r_PlanarProfile_hpp_

#include "FaceExtrusion.hpp"

namespace Slic3r {

struct PlanarProfileResult {
    indexed_triangle_set cap; // Open cap; original 3D positions, checked for mesh precision.
    Vec3d origin = Vec3d::Zero();
    Vec3d normal = Vec3d::Zero(); // Unit normal following the outer contour ordering.
    double area_mm2 = 0.0; // Projected outer area minus holes, in world square millimeters.
    size_t boundary_edges = 0;
    FaceExtrusionError error = FaceExtrusionError::None;
    std::string message;

    explicit operator bool() const noexcept { return error == FaceExtrusionError::None; }
};

// Ordered world-mm contours, implicitly closed (do not repeat the first point).
// The outer contour and each optional hole require at least three distinct,
// finite, non-collinear points. Contours must be coplanar within options, simple,
// and mutually disjoint; holes must be strictly inside the outer contour and
// cannot overlap, nest or touch. Hole orientation may be either direction.
// Projection is used only for topology/triangulation: no hull, tolerance weld,
// arbitrary reordering or silent flattening of the original 3D coordinates.
// Float mesh storage is checked against max(length, planarity) in mm, and its
// topology/planarity are revalidated. A maximum of 4096 total points bounds the
// synchronous work; larger profiles return TooLarge. Input is immutable. On
// failure cap is empty and message explains the rejection.
PlanarProfileResult triangulate_planar_profile(
    const std::vector<Vec3d> &outer,
    const std::vector<std::vector<Vec3d>> &holes = {},
    const FaceExtrusionOptions &options = {});

// Reuses extrude_planar_face on the validated cap. Displacement is in the same
// world-mm frame, normal to the profile; either sign is supported. Failure never
// returns a partial mesh or changes input contours.
FaceExtrusionResult extrude_planar_profile(
    const std::vector<Vec3d> &outer, const Vec3d &displacement,
    const std::vector<std::vector<Vec3d>> &holes = {},
    const FaceExtrusionOptions &options = {});

} // namespace Slic3r

#endif // slic3r_PlanarProfile_hpp_
