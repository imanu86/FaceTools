// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#ifndef slic3r_CutPlaneSelection_hpp_
#define slic3r_CutPlaneSelection_hpp_

#include "Point.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace Slic3r {

struct CutPlaneFaceGuide {
    Vec3d point_world = Vec3d::Zero();
    Vec3d outward_normal_world = Vec3d::Zero();
};

struct CutPlaneSelectionResult {
    Vec3d moving_normal = Vec3d::Zero();
    int side = 0;
    std::string message;
    explicit operator bool() const noexcept { return side != 0; }
};

namespace CutPlaneSelectionDetail {
inline bool unit_normal(const Vec3d &normal, Vec3d &unit)
{
    if (!normal.allFinite()) return false;
    const double scale = normal.cwiseAbs().maxCoeff();
    if (!(scale > 0.0)) return false;
    // Avoid overflow/underflow of norm() for large or subnormal finite normals.
    const Vec3d scaled = normal / scale;
    unit = scaled / scaled.norm();
    return unit.allFinite();
}

inline CutPlaneSelectionResult failure(const char *message)
{
    CutPlaneSelectionResult result; // invalid results never expose a partial side/normal
    result.message = message;
    return result;
}
} // namespace CutPlaneSelectionDetail

// Guides identify the WHOLE moving half-space, not individually deformable faces.
// Coordinates and outward normals are already in WORLD space. Never mutate them,
// reproject a guide, average points/normals, or decide by a majority vote.
inline CutPlaneSelectionResult select_cut_moving_side(
    const std::vector<CutPlaneFaceGuide> &guides,
    const Vec3d &plane_origin,
    const Vec3d &plane_normal,
    double distance_tolerance_mm = 1e-4)
{
    using CutPlaneSelectionDetail::failure;
    if (guides.empty())
        return failure("Select at least one guide face before confirming the cut plane.");
    if (guides.size() > 4096)
        return failure("Too many guide faces (maximum 4096). Select fewer representative faces.");
    if (!plane_origin.allFinite() || !std::isfinite(distance_tolerance_mm) || !(distance_tolerance_mm > 0.0))
        return failure("The plane origin and positive distance tolerance must be finite.");
    Vec3d normal;
    if (!CutPlaneSelectionDetail::unit_normal(plane_normal, normal))
        return failure("The cut plane needs a finite, nonzero normal.");

    constexpr double roundoff_factor = 64.0 * std::numeric_limits<double>::epsilon();
    constexpr double parallel_roundoff = 128.0 * std::numeric_limits<double>::epsilon();
    // Arithmetic uncertainty may NEVER grow into the geometric near-plane band.
    // A far-away origin/ill-conditioned cancellation must reject, not guess a side.
    const double max_coplanar_roundoff = std::min(1e-10, distance_tolerance_mm * 1e-6);
    int selected_side = 0;
    for (const CutPlaneFaceGuide &guide : guides) {
        Vec3d outward;
        if (!guide.point_world.allFinite() || !CutPlaneSelectionDetail::unit_normal(guide.outward_normal_world, outward))
            return failure("Every guide needs a finite point and a finite, nonzero outward normal.");
        const Vec3d delta = guide.point_world - plane_origin;
        if (!delta.allFinite())
            return failure("Guide coordinates exceed the supported numeric range. Place the plane origin near the model.");
        const Vec3d terms = delta.cwiseProduct(normal);
        const double magnitude_sum = terms.cwiseAbs().sum();
        const double distance = terms.sum();
        if (!std::isfinite(magnitude_sum) || !std::isfinite(distance))
            return failure("The signed distance exceeds the supported numeric range.");
        // Relative rounding budget for subtraction, normalization and the dot
        // product. No absolute geometric floor: an axis-aligned offset of 8e-5
        // (or even 1e-12) is NOT coplanar merely because it is inside the 1e-4 band.
        const double arithmetic_error = roundoff_factor * magnitude_sum;
        const double abs_distance = std::abs(distance);
        int side = 0;
        if (abs_distance > distance_tolerance_mm && abs_distance - distance_tolerance_mm > arithmetic_error) {
            side = distance > 0.0 ? 1 : -1;
        } else if (abs_distance <= arithmetic_error && arithmetic_error <= max_coplanar_roundoff) {
            // Use angular error itself, not 1-|dot| (which hides small real tilts).
            // A point on the plane with a nonparallel face is ambiguous.
            if (outward.cross(normal).norm() > parallel_roundoff)
                return failure("A guide lies on the cut plane but its face is not parallel. Move the plane or choose another guide.");
            const double outward_dot_plane = outward.dot(normal);
            if (!std::isfinite(outward_dot_plane) || std::abs(outward_dot_plane) < 0.5)
                return failure("The coplanar guide normal is ambiguous.");
            // The material is INSIDE the face: e.g. underside normal -Z selects +Z.
            side = outward_dot_plane > 0.0 ? -1 : 1;
        } else {
            return failure("A guide is too close to the cut plane to determine its side reliably. Move the plane away from the guide or pick an exactly coplanar parallel face.");
        }
        if (selected_side != 0 && selected_side != side)
            return failure("The guide faces select opposite sides of the cut plane. Keep guides on one moving side; no majority vote is applied.");
        selected_side = side;
    }
    CutPlaneSelectionResult result;
    result.side = selected_side;
    result.moving_normal = double(selected_side) * normal;
    return result;
}

} // namespace Slic3r
#endif
