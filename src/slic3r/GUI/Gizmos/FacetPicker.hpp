#ifndef slic3r_FacetPicker_hpp_
#define slic3r_FacetPicker_hpp_

#include "slic3r/GUI/GLModel.hpp"
#include "libslic3r/Point.hpp"

#include <memory>

namespace Slic3r {
class ModelVolume;

namespace GUI {
class CommonGizmosDataPool;
class Selection;
struct Camera; // declared as a struct in Camera.hpp; matching it avoids MSVC ABI linker errors

// Picks a single triangular facet of the selected object under the mouse and draws it
// highlighted.
//
// It deliberately knows nothing about what a picked facet is *for* - the caller decides
// that - so it can serve any gizmo that needs "which triangle is under the cursor".
// GLGizmoAdvancedCut uses it to place the cut plane; GLGizmoFlatten does the same job
// inline today and could be moved onto this later.
class FacetPicker
{
public:
    struct Hit
    {
        int                facet             = -1;
        const ModelVolume *volume            = nullptr;
        std::shared_ptr<const TriangleMesh> mesh; // immutable source identity and lifetime
        Transform3d        trafo             = Transform3d::Identity(); // volume -> world
        Vec3d              world_pos         = Vec3d::Zero(); // approximate ray hit, for display
        Vec3d              world_plane_point = Vec3d::Zero(); // original triangle vertex, transformed in double
        Vec3d              world_normal      = Vec3d::Zero();

        bool valid() const { return facet >= 0 && volume != nullptr && mesh != nullptr; }
        bool same_geometry(const Hit &other) const {
            return valid() && other.valid() && volume == other.volume && mesh == other.mesh &&
                facet == other.facet && (trafo.matrix().array() == other.trafo.matrix().array()).all();
        }
    };

    bool is_active() const { return m_active; }
    void set_active(bool active);

    // Drops the cached hit and the highlight geometry.
    void reset();

    // Recomputes the cached hit from the mouse position. Returns true when the ray hits a
    // model part; on a miss the cache is left invalid.
    bool update(const Vec2d &mouse_position, const CommonGizmosDataPool *pool, const Selection &selection, const Camera &camera);

    // Draws the cached facet. No-op when there is no valid hit.
    void render(const Camera &camera);

    const Hit &hit() const { return m_hit; }

    // Explicit original-facet intent only: never infer a plane from a rounded ray hit.
    // Invalid indices, vertices, facets or nonfinite/non-affine/singular transforms
    // return an invalid hit without modifying the source. No GL context is needed.
    static Hit from_source_facet(const ModelVolume &volume, size_t facet, const Transform3d &world);

    // Inverse-transpose preserves outward orientation under non-uniform scale and
    // mirrors. Returns zero for an invalid transform or nonfinite/zero normal.
    static Vec3d world_normal(const Transform3d &trafo, const Vec3d &mesh_normal);

private:
    bool    m_active{false};
    Hit     m_hit;
    Hit     m_rendered_hit; // volume + immutable mesh + facet + exact pose, not facet index alone
    GLModel m_highlight;
};

} // namespace GUI
} // namespace Slic3r

#endif // slic3r_FacetPicker_hpp_
