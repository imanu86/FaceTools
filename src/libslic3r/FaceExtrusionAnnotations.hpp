// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#ifndef slic3r_FaceExtrusionAnnotations_hpp_
#define slic3r_FaceExtrusionAnnotations_hpp_

#include "Model.hpp"

#include <array>
#include <cstddef>
#include <cstring>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

namespace Slic3r {

// Presence of stored per-facet annotation data, not a claim about who created
// it. Imported/generated annotations and manual painting share these fields.
// This summary is for UI reasons only; do NOT use flags/counts as a stale-job
// or consent check. Use the complete content stamp below for that purpose.
struct FaceExtrusionAnnotationInfo {
    bool support = false;
    bool seam = false;
    bool mmu = false;
    bool fuzzy_skin = false;
    bool exterior = false;
    size_t legacy = 0;
    size_t legacy_small_overhang = 0;
    size_t legacy_small_hole = 0;
    size_t legacy_exterior = 0;
    size_t legacy_other = 0;

    bool any() const noexcept
    {
        return support || seam || mmu || fuzzy_skin || exterior || legacy != 0;
    }
    bool operator==(const FaceExtrusionAnnotationInfo &other) const noexcept
    {
        return std::tie(support, seam, mmu, fuzzy_skin, exterior, legacy,
                        legacy_small_overhang, legacy_small_hole, legacy_exterior, legacy_other) ==
               std::tie(other.support, other.seam, other.mmu, other.fuzzy_skin, other.exterior, other.legacy,
                        other.legacy_small_overhang, other.legacy_small_hole, other.legacy_exterior, other.legacy_other);
    }
    bool operator!=(const FaceExtrusionAnnotationInfo &other) const noexcept { return !(*this == other); }
};

inline FaceExtrusionAnnotationInfo inspect_face_extrusion_annotations(const ModelVolume &volume)
{
    auto present = [](const FacetsAnnotation &annotation) {
        const auto &data = annotation.get_data();
        return !data.first.empty() || !data.second.empty();
    };
    FaceExtrusionAnnotationInfo info;
    info.support = present(volume.supported_facets);
    info.seam = present(volume.seam_facets);
    info.mmu = present(volume.mmu_segmentation_facets);
    info.fuzzy_skin = present(volume.fuzzy_skin_facets);
    info.exterior = present(volume.exterior_facets);
    if (const TriangleMesh *mesh = volume.mesh_ptr()) {
        for (const FaceProperty &property : mesh->its.properties) {
            if (property.type == eNormal)
                continue; // Ordinary initialized mesh properties are NOT annotations.
            ++info.legacy;
            switch (property.type) {
            case eSmallOverhang:       ++info.legacy_small_overhang; break;
            case eSmallHole:           ++info.legacy_small_hole; break;
            case eExteriorAppearance: ++info.legacy_exterior; break;
            default:                  ++info.legacy_other; break;
            }
        }
    }
    return info;
}

// Exact, collision-free content snapshot. Covers every triangle ID, stream
// offset and encoded bit in all five layers, and the type/area of EVERY legacy
// property in order, including normal entries. Double areas are compared by
// representation so even unchanged NaNs compare equal. No struct padding is
// compared. IDs, timestamps and derived caches deliberately are not payload.
//
// Capture at consent/job enqueue and compare with the live source immediately
// before commit. This may copy substantial data: it is not a per-frame hash.
// The caller still checks volume identity, geometry, pose/config and selection
// independently, and must not read a volume concurrently with a mutation.
struct FaceExtrusionAnnotationStamp {
    FaceExtrusionAnnotationStamp() = default;
    bool operator==(const FaceExtrusionAnnotationStamp &other) const
    {
        return m_layers == other.m_layers && m_legacy == other.m_legacy;
    }
    bool operator!=(const FaceExtrusionAnnotationStamp &other) const { return !(*this == other); }

private:
    using FacetPayload = std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>;
    struct LegacyProperty {
        int type = 0;
        std::array<unsigned char, sizeof(double)> area{};
        bool operator==(const LegacyProperty &other) const { return type == other.type && area == other.area; }
    };
    std::array<FacetPayload, 5> m_layers;
    std::vector<LegacyProperty> m_legacy;
    friend FaceExtrusionAnnotationStamp capture_face_extrusion_annotation_stamp(const ModelVolume &);
};

inline FaceExtrusionAnnotationStamp capture_face_extrusion_annotation_stamp(const ModelVolume &volume)
{
    FaceExtrusionAnnotationStamp stamp;
    stamp.m_layers = {{volume.supported_facets.get_data(), volume.seam_facets.get_data(),
                       volume.mmu_segmentation_facets.get_data(), volume.fuzzy_skin_facets.get_data(),
                       volume.exterior_facets.get_data()}};
    if (const TriangleMesh *mesh = volume.mesh_ptr()) {
        stamp.m_legacy.reserve(mesh->its.properties.size());
        for (const FaceProperty &property : mesh->its.properties) {
            FaceExtrusionAnnotationStamp::LegacyProperty entry;
            entry.type = static_cast<int>(property.type);
            std::memcpy(entry.area.data(), &property.area, sizeof(property.area));
            stamp.m_legacy.push_back(entry);
        }
    }
    return stamp;
}

// Destructive ONLY to the target volume's facet annotations. Call only on a
// prepared successful boolean result after explicit consent, with a matching
// source stamp, or in an equivalently guarded successful commit. Never clear
// the live source while computing, after an error/cancel, or without consent.
// This helper does not provide GUI/Undo/transaction synchronization itself.
//
// All volume/config/layer IDs and the part GUID are preserved. Nonempty layers
// use FacetsAnnotation::reset(), which updates their undo timestamps; empty
// layers are left untouched. No geometry, transform, material, name, config,
// text/emboss, cut/provenance, other volume or modifier is changed. Legacy
// properties use copy-on-write, so any model copy sharing the old mesh retains
// its annotations. Allocation is completed before the first mutation.
inline void clear_face_extrusion_annotations(ModelVolume &volume)
{
    const FaceExtrusionAnnotationInfo info = inspect_face_extrusion_annotations(volume);
    if (!info.any())
        return;

    std::shared_ptr<const TriangleMesh> prepared_mesh;
    if (info.legacy != 0) {
        auto detached = std::make_shared<TriangleMesh>(volume.mesh());
        for (FaceProperty &property : detached->its.properties)
            if (property.type != eNormal) {
                property.type = eNormal;
                property.area = 0.0;
            }
        prepared_mesh = std::move(detached);
    }

    // shared_ptr assignment and the reset operations below do not allocate.
    if (prepared_mesh)
        volume.set_mesh(prepared_mesh);
    auto reset_if_present = [](FacetsAnnotation &annotation) {
        const auto &data = annotation.get_data();
        if (!data.first.empty() || !data.second.empty())
            annotation.reset();
    };
    reset_if_present(volume.supported_facets);
    reset_if_present(volume.seam_facets);
    reset_if_present(volume.mmu_segmentation_facets);
    reset_if_present(volume.fuzzy_skin_facets);
    reset_if_present(volume.exterior_facets);
    if (info.mmu) {
        // Derived MMU cache: force normal get_extruders() recomputation, including
        // its private default-extruder flag, on the next query.
        volume.mmuseg_extruders.clear();
        volume.mmuseg_ts = 0;
    }
}

} // namespace Slic3r

#endif // slic3r_FaceExtrusionAnnotations_hpp_
