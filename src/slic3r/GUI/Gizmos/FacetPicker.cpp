#include "FacetPicker.hpp"

#include "GLGizmoBase.hpp"
#include "GLGizmosCommon.hpp"
#include "slic3r/GUI/Camera.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/Selection.hpp"
#include "libslic3r/Model.hpp"

#include <GL/glew.h>

#include <Eigen/LU>
#include <array>
#include <cmath>
#include <limits>
#include <optional>

namespace Slic3r {
namespace GUI {

namespace {

std::optional<Matrix3d> affine_inverse(const Transform3d &trafo)
{
    if (!trafo.matrix().allFinite() ||
        (trafo.matrix().row(3).array() != Eigen::RowVector4d(0, 0, 0, 1).array()).any())
        return std::nullopt;
    Eigen::FullPivLU<Matrix3d> lu(trafo.linear());
    lu.setThreshold(0.0); // Reject singularity, not an arbitrary geometric scale.
    if (!lu.isInvertible())
        return std::nullopt;
    const Matrix3d inverse = lu.inverse();
    return inverse.allFinite() ? std::optional<Matrix3d>(inverse) : std::nullopt;
}

Vec3d finite_unit_vector(const Vec3d &value)
{
    if (!value.allFinite()) return Vec3d::Zero();
    const double scale = value.cwiseAbs().maxCoeff();
    if (!(scale > 0.0)) return Vec3d::Zero();
    const Vec3d scaled = value / scale;
    return Vec3d(scaled / scaled.norm());
}

} // namespace

Vec3d FacetPicker::world_normal(const Transform3d &trafo, const Vec3d &mesh_normal)
{
    const auto inverse = affine_inverse(trafo);
    const Vec3d local = finite_unit_vector(mesh_normal);
    return inverse && !local.isZero() ? finite_unit_vector(inverse->transpose() * local) : Vec3d::Zero();
}

FacetPicker::Hit FacetPicker::from_source_facet(const ModelVolume &volume, size_t facet, const Transform3d &world)
{
    if (!volume.is_model_part() || !volume.mesh_ptr() || facet > size_t(std::numeric_limits<int>::max()))
        return {};
    const auto mesh = volume.get_mesh_shared_ptr();
    const indexed_triangle_set &its = mesh->its;
    if (facet >= its.indices.size()) return {};
    const Vec3i &triangle = its.indices[facet];
    std::array<Vec3d, 3> vertices;
    for (int i = 0; i < 3; ++i) {
        if (triangle[i] < 0 || size_t(triangle[i]) >= its.vertices.size()) return {};
        vertices[i] = its.vertices[triangle[i]].cast<double>();
        if (!vertices[i].allFinite()) return {};
    }
    // Float mesh coordinates are promoted BEFORE subtraction and cross product.
    const Vec3d normal = world_normal(world, (vertices[1] - vertices[0]).cross(vertices[2] - vertices[0]));
    if (normal.isZero()) return {};
    std::array<Vec3d, 3> world_vertices;
    for (int i = 0; i < 3; ++i) {
        world_vertices[i] = world * vertices[i];
        if (!world_vertices[i].allFinite() || !world_vertices[i].cast<float>().allFinite()) return {};
    }
    Hit hit;
    hit.facet = int(facet);
    hit.volume = &volume;
    hit.mesh = mesh;
    hit.trafo = world;
    hit.world_pos = hit.world_plane_point = world_vertices[0];
    hit.world_normal = normal;
    return hit;
}

void FacetPicker::set_active(bool active)
{
    if (m_active == active)
        return;
    m_active = active;
    if (!m_active)
        reset();
}

void FacetPicker::reset()
{
    m_hit = Hit();
    m_rendered_hit = Hit();
    m_highlight.reset();
}

bool FacetPicker::update(const Vec2d &mouse_position, const CommonGizmosDataPool *pool, const Selection &selection, const Camera &camera)
{
    m_hit = Hit();

    if (!mouse_position.allFinite() || !camera.get_position().allFinite() ||
        pool == nullptr || pool->raycaster() == nullptr || pool->selection_info() == nullptr)
        return false;

    const ModelObject *mo           = pool->selection_info()->model_object();
    const int          instance_idx = selection.get_instance_idx();
    if (mo == nullptr || instance_idx < 0 || instance_idx >= int(mo->instances.size()))
        return false;

    // Raycaster builds one MeshRaycaster per model-part volume, in mo->volumes order, so
    // build the matching transform list the same way.
    const ModelInstance *mi        = mo->instances[instance_idx];
    const double         sla_shift = pool->selection_info()->get_sla_shift();

    std::vector<Transform3d>         trafo_matrices;
    std::vector<const ModelVolume *> part_volumes;
    for (const ModelVolume *mv : mo->volumes) {
        if (mv->is_model_part()) {
            trafo_matrices.emplace_back(Geometry::translation_transform(sla_shift * Vec3d::UnitZ()) *
                                        mi->get_transformation().get_matrix() * mv->get_matrix());
            part_volumes.emplace_back(mv);
        }
    }

    const std::vector<const MeshRaycaster *> raycasters = pool->raycaster()->raycasters();
    if (raycasters.size() != trafo_matrices.size())
        return false; // hollowed-mesh or stale-cache case; skip rather than mis-index

    // A same-size list may still contain a hollowed or stale mesh. Original facet
    // indices are meaningful only for the exact immutable source mesh referenced by
    // the raycaster; IndexedMesh exposes references to its original vertex/index arrays.
    for (size_t i = 0; i < raycasters.size(); ++i) {
        if (!raycasters[i] || !part_volumes[i]->mesh_ptr()) return false;
        const auto &ray_mesh = raycasters[i]->get_aabb_mesh();
        const auto &source = part_volumes[i]->mesh().its;
        if (&ray_mesh.vertices() != &source.vertices || &ray_mesh.indices() != &source.indices)
            return false;
    }

    double closest_sq = std::numeric_limits<double>::max();
    for (size_t mesh_id = 0; mesh_id < trafo_matrices.size(); ++mesh_id) {
        const Transform3d &world = trafo_matrices[mesh_id];
        if (!affine_inverse(world)) continue; // Never pass a singular matrix to the raycaster.
        Vec3f ray_hit = Vec3f::Zero(), ray_normal = Vec3f::Zero();
        size_t facet = 0;
        if (raycasters[mesh_id]->unproject_on_mesh(mouse_position, world, camera, ray_hit, ray_normal, nullptr, &facet)) {
            const Vec3d world_hit = world * ray_hit.cast<double>();
            if (!world_hit.allFinite() || !world_hit.cast<float>().allFinite()) continue;
            const double sq = (camera.get_position() - world_hit).squaredNorm();
            if (!std::isfinite(sq) || sq >= closest_sq) continue;
            Hit candidate = from_source_facet(*part_volumes[mesh_id], facet, world);
            closest_sq = sq;
            if (!candidate.valid()) { m_hit = Hit(); continue; } // Do not pick through an invalid nearest facet.
            // Keep volume, original facet and approximate ray position together at
            // the nearest-hit branch. The plane support/normal never use ray floats.
            candidate.world_pos = world_hit;
            m_hit = std::move(candidate);
        }
    }
    return m_hit.valid();
}

void FacetPicker::render(const Camera &camera)
{
    if (!m_hit.valid())
        return;

    if (!m_rendered_hit.same_geometry(m_hit)) {
        m_rendered_hit = Hit();
        m_highlight.reset();

        const indexed_triangle_set &its = m_hit.mesh->its;
        const auto &tri = its.indices[m_hit.facet];
        // Lift only the display overlay, never the original-facet plane support.
        const Vec3d bias = m_hit.world_normal * 0.05;
        std::array<Vec3f, 3> vertices;
        for (int i = 0; i < 3; ++i) {
            const Vec3d v = m_hit.trafo * its.vertices[tri[i]].cast<double>() + bias;
            if (!v.allFinite() || !v.cast<float>().allFinite()) return;
            vertices[i] = v.cast<float>();
        }
        if (finite_unit_vector((vertices[1].cast<double>() - vertices[0].cast<double>()).cross(
                vertices[2].cast<double>() - vertices[0].cast<double>())).isZero())
            return; // A collapsed float overlay must never upload an invalid normal.
        GLModel::Geometry geometry;
        geometry.format = {GLModel::PrimitiveType::Triangles, GLModel::Geometry::EVertexLayout::P3N3};
        for (const Vec3f &vertex : vertices)
            geometry.add_vertex(vertex, Vec3f(m_hit.world_normal.cast<float>()));
        geometry.add_triangle(0, 1, 2);
        m_highlight.init_from(std::move(geometry));
        m_rendered_hit = m_hit;
    }

    if (!m_highlight.is_initialized())
        return;

    const auto &p_flat_shader = wxGetApp().get_shader("flat");
    if (!p_flat_shader)
        return;

    glsafe(::glEnable(GL_DEPTH_TEST));
    glsafe(::glDisable(GL_CULL_FACE));
    wxGetApp().bind_shader(p_flat_shader);
    p_flat_shader->set_uniform("projection_matrix", camera.get_projection_matrix());
    // Vertices are already in world space, so the model part of the matrix is identity.
    p_flat_shader->set_uniform("view_model_matrix", camera.get_view_matrix());
    m_highlight.set_color(GLGizmoBase::FLATTEN_HOVER_COLOR);
    m_highlight.render_geometry();
    wxGetApp().unbind_shader();
    glsafe(::glEnable(GL_CULL_FACE));
}

} // namespace GUI
} // namespace Slic3r
