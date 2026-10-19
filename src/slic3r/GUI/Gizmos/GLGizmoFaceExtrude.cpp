#include "GLGizmoFaceExtrude.hpp"
#include "FaceToolsPreviewBounds.hpp"

#include "libslic3r/FaceExtrusion.hpp"
#include "libslic3r/FaceExtrusionAnnotations.hpp"
#include "libslic3r/PlanarProfile.hpp"
#include "libslic3r/CutAndStretch.hpp"
#include "libslic3r/CutPlaneSelection.hpp"
#include "libslic3r/CutGuideSelection.hpp"
#include "libslic3r/MeshBoolean.hpp"
#include "libslic3r/MeshSolidValidation.hpp"
#include "slic3r/GUI/Jobs/Worker.hpp"
#include "slic3r/GUI/Camera.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/GUI_ObjectList.hpp"
#include "slic3r/GUI/ImGuiWrapper.hpp"
#include "slic3r/GUI/MeshUtils.hpp"
#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/Selection.hpp"

#include <GL/glew.h>
#include <imgui/imgui.h>
#include <wx/glcanvas.h>

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <utility>
#include <cmath>
#include <exception>
#include <limits>
#include <numeric>
#include <utility>

namespace Slic3r {
namespace GUI {
namespace {

constexpr double bed_tolerance_mm = face_extrusion_bed_tolerance_mm;
using FaceExtrusionModelDetail::usable_transform;
using FaceExtrusionModelDetail::usable_mesh;
using FaceExtrusionModelDetail::finite_bounds;

CutAndStretchScope cut_scope(bool guided)
{
    return guided ? CutAndStretchScope::GuidedComponent : CutAndStretchScope::WholeHalfSpace;
}

std::optional<Vec3d> world_normal(const Vec3d &local, const Transform3d &world)
{
    if (!local.allFinite() || !usable_transform(world))
        return std::nullopt;
    Vec3d normal = world.linear().inverse().transpose() * local;
    const double length = normal.norm();
    if (!normal.allFinite() || !std::isfinite(length) || length <= 0.0)
        return std::nullopt;
    return Vec3d(normal / length);
}

std::optional<Vec3d> positioned_cut_plane(const Vec3d &reference, const Vec3d &normal, double coordinate)
{
    if (!reference.allFinite() || !normal.allFinite() || !std::isfinite(coordinate) ||
        std::abs(normal.squaredNorm() - 1.0) > 1e-8) return std::nullopt;
    const double displacement = coordinate - normal.dot(reference);
    const Vec3d origin = reference + displacement * normal;
    // Geometry is stored in a source-relative physical frame, NOT absolute WORLD
    // floats. A later plate must not spend the mesh precision budget on its offset.
    // Still reject huge relative moves / ill-conditioned double world arithmetic
    // BEFORE publishing the origin to GL/camera bounds; never clamp or move source.
    const double float_error = 2.0 * double(std::numeric_limits<float>::epsilon()) * std::abs(displacement);
    const double world_scale = std::max({1.0, reference.cwiseAbs().maxCoeff(), origin.cwiseAbs().maxCoeff(), std::abs(coordinate)});
    const double double_error = 8.0 * std::numeric_limits<double>::epsilon() * world_scale;
    return origin.allFinite() && std::isfinite(float_error) && std::isfinite(double_error) &&
        float_error + double_error <= 1e-4 ? std::optional<Vec3d>(origin) : std::nullopt;
}

std::optional<size_t> exact_cut_plane_facet(bool face_basis, const std::optional<size_t> &facet,
    size_t triangle_count, const Vec3d &normal, const Vec3d &reference, double coordinate)
{
    // Only explicit original-face provenance at EXACTLY zero offset authorizes
    // the facet-plane kernel. Never promote a nearby/typed plane by tolerance.
    return face_basis && facet && *facet < triangle_count && normal.allFinite() && reference.allFinite() &&
        std::isfinite(coordinate) && coordinate == normal.dot(reference) ? facet : std::nullopt;
}

std::optional<Vec2d> profile_screen_point(const Vec3d &point, const Camera &camera)
{
    const Vec4d clip = camera.get_projection_matrix().matrix() * camera.get_view_matrix().matrix() *
        Vec4d(point.x(), point.y(), point.z(), 1.0);
    const auto &viewport = camera.get_viewport();
    if (!clip.allFinite() || clip.w() <= 1e-12 || viewport[2] <= 0 || viewport[3] <= 0)
        return std::nullopt;
    const Vec3d ndc = clip.head<3>() / clip.w();
    if ((ndc.array() < -1.0).any() || (ndc.array() > 1.0).any())
        return std::nullopt;
    return Vec2d(viewport[0] + 0.5 * (ndc.x() + 1.0) * viewport[2],
                 viewport[1] + 0.5 * (1.0 - ndc.y()) * viewport[3]);
}

std::shared_ptr<GLModel> profile_edge_model(const std::vector<std::vector<Vec3d>> &contours)
{
    GLModel::Geometry geometry;
    geometry.format = {GLModel::PrimitiveType::Lines, GLModel::Geometry::EVertexLayout::P3N3};
    for (const auto &contour : contours) {
        if (contour.size() < 2) continue;
        const size_t edges = contour.size() > 2 ? contour.size() : contour.size() - 1;
        for (size_t i = 0; i < edges; ++i) {
            const Vec3f a = contour[i].cast<float>(), b = contour[(i + 1) % contour.size()].cast<float>();
            if (!a.allFinite() || !b.allFinite()) return {};
            const unsigned int start = unsigned(geometry.vertices_count());
            const Vec3f normal = Vec3f::UnitZ();
            geometry.add_vertex(a, normal);
            geometry.add_vertex(b, normal);
            geometry.add_line(start, start + 1);
        }
    }
    if (geometry.is_empty()) return {};
    auto model = std::make_shared<GLModel>();
    model->init_from(std::move(geometry));
    return model;
}

std::vector<float> bridge_base_heights(const indexed_triangle_set &unit)
{
    // Exact XY pairing, not tolerance welding: base Z may contain the cutter's
    // validated float residual and MUST NOT be multiplied by the distance.
    std::map<std::pair<float, float>, float> lower;
    for (const Vec3f &vertex : unit.vertices) {
        if (!vertex.allFinite()) throw std::runtime_error("Invalid unit bridge vertex.");
        if (std::abs(vertex.z()) <= 1e-4f) {
            const auto inserted = lower.emplace(std::make_pair(vertex.x(), vertex.y()), vertex.z());
            if (!inserted.second && inserted.first->second != vertex.z())
                throw std::runtime_error("Ambiguous exact XY pairing in the contact bridge.");
        } else if (std::abs(vertex.z() - 1.0f) > 1e-4f) {
            throw std::runtime_error("The prepared bridge is not a unit extrusion.");
        }
    }
    std::vector<float> result;
    result.reserve(unit.vertices.size());
    size_t upper = 0;
    for (const Vec3f &vertex : unit.vertices) {
        const auto found = lower.find(std::make_pair(vertex.x(), vertex.y()));
        if (found == lower.end()) throw std::runtime_error("Missing exact lower-cap bridge vertex.");
        if (vertex.z() > 0.5f) {
            ++upper;
            if (vertex.z() != float(double(found->second) + 1.0))
                throw std::runtime_error("The unit bridge does not preserve its exact base heights.");
        }
        result.push_back(found->second);
    }
    if (lower.empty() || upper == 0) throw std::runtime_error("The prepared bridge has no paired caps.");
    return result;
}

void update_bridge_preview(indexed_triangle_set &preview, const indexed_triangle_set &unit,
                           const std::vector<float> &base_z, double distance)
{
    if (!std::isfinite(distance) || distance <= 0.0 || preview.vertices.size() != unit.vertices.size() ||
        base_z.size() != unit.vertices.size()) throw std::runtime_error("Invalid bridge preview distance or cache.");
    for (size_t i = 0; i < unit.vertices.size(); ++i) {
        if (unit.vertices[i].z() > 0.5f) {
            const float z = float(double(base_z[i]) + distance);
            if (!std::isfinite(z)) throw std::runtime_error("Bridge preview exceeds float precision.");
            preview.vertices[i].z() = z; // bottom and XY stay bit-exact; never accumulate a delta
        }
    }
}

std::unique_ptr<GLModel> bridge_preview_model(const indexed_triangle_set &mesh)
{
    // This is the GL upload's normal generation, not a topology/CSG pass. Use
    // double intermediates so an extreme numeric distance cannot send NaNs to GL.
    if (mesh.indices.size() > size_t(std::numeric_limits<unsigned int>::max()) / 3)
        throw std::runtime_error("Bridge preview exceeds GL index limits.");
    GLModel::Geometry geometry;
    geometry.format = {GLModel::PrimitiveType::Triangles, GLModel::Geometry::EVertexLayout::P3N3};
    geometry.reserve_vertices(3 * mesh.indices.size());
    geometry.reserve_indices(3 * mesh.indices.size());
    for (const Vec3i &face : mesh.indices) {
        const Vec3d a = mesh.vertices[face[0]].cast<double>();
        const Vec3d cross = (mesh.vertices[face[1]].cast<double>() - a).cross(mesh.vertices[face[2]].cast<double>() - a);
        const double length = cross.stableNorm();
        if (!cross.allFinite() || !std::isfinite(length) || length <= 0.0)
            throw std::runtime_error("The extension is too small or large for a finite bridge preview.");
        const Vec3f normal = (cross / length).cast<float>();
        const unsigned int start = unsigned(geometry.vertices_count());
        for (int index : {face[0], face[1], face[2]}) geometry.add_vertex(mesh.vertices[index], normal);
        geometry.add_triangle(start, start + 1, start + 2);
    }
    auto model = std::make_unique<GLModel>();
    model->init_from(std::move(geometry));
    return model;
}

std::optional<double> ray_axis_parameter(const Linef3 &ray, const Vec3d &anchor, const Vec3d &axis)
{
    if (!ray.a.allFinite() || !ray.b.allFinite() || !anchor.allFinite() || !axis.allFinite())
        return std::nullopt;
    const Vec3d delta = ray.b - ray.a;
    const double ray_length = delta.norm();
    if (!std::isfinite(ray_length) || ray_length <= 1e-12 || std::abs(axis.squaredNorm() - 1.0) > 1e-6)
        return std::nullopt;
    const Vec3d direction = delta / ray_length;
    const Vec3d offset = ray.a - anchor;
    const double alignment = axis.dot(direction);
    const double denominator = 1.0 - alignment * alignment;
    // An axis viewed end-on has no stable mouse-to-distance solution. Keep the
    // last value instead of dividing by almost zero or inventing depth motion.
    if (!std::isfinite(denominator) || denominator <= 1e-4)
        return std::nullopt;
    const double parameter = (axis.dot(offset) - alignment * direction.dot(offset)) / denominator;
    const double along_ray = alignment * parameter - direction.dot(offset);
    if (!std::isfinite(parameter) || !std::isfinite(along_ray) || along_ray < 0.0)
        return std::nullopt;
    return parameter;
}

bool same_screen_camera(const Transform3d &view, const Transform3d &projection,
                        const Transform3d &previous_view, const Transform3d &previous_projection)
{
    const Matrix4d current = projection.matrix() * view.matrix();
    const Matrix4d previous = previous_projection.matrix() * previous_view.matrix();
    // Ignore clip-depth-only changes: updating preview bounds changes near/far
    // every frame. Rebasing on those changes would make the arrow drag stall.
    // The x/y AND homogeneous-w rows fully determine screen positions, including
    // the distinction between perspective and orthographic projection.
    return current.allFinite() && previous.allFinite() &&
        current.topRows<2>().isApprox(previous.topRows<2>(), 1e-10) &&
        current.row(3).isApprox(previous.row(3), 1e-10);
}

bool concordant_normals(const Vec3d &a, const Vec3d &b)
{
    // Same signed WORLD direction; opposite-facing parallel surfaces are not
    // part of a batch. Each patch is still revalidated by the unchanged kernel.
    return a.dot(b) >= std::cos(FaceExtrusionOptions{}.angular_tolerance_rad);
}

std::shared_ptr<GLModel> make_hover_model(const indexed_triangle_set &source,
                                         const std::vector<int> &triangles,
                                         const Transform3d &world)
{
    if (triangles.empty() || triangles.size() > size_t(std::numeric_limits<int>::max()) / 3 ||
        !usable_transform(world))
        return nullptr;
    // Copy only the picked connected patch, not the entire source mesh on each
    // hovered plane transition. Original triangles preserve all holes and caps.
    indexed_triangle_set patch;
    patch.vertices.reserve(3 * triangles.size());
    patch.indices.reserve(triangles.size());
    for (int triangle : triangles) {
        if (triangle < 0 || size_t(triangle) >= source.indices.size())
            return nullptr;
        const int first = int(patch.vertices.size());
        for (int corner = 0; corner < 3; ++corner) {
            const int vertex = source.indices[triangle][corner];
            if (vertex < 0 || size_t(vertex) >= source.vertices.size())
                return nullptr;
            patch.vertices.push_back(source.vertices[vertex]);
        }
        patch.indices.emplace_back(first, first + 1, first + 2);
    }
    TriangleMesh mesh(std::move(patch));
    mesh.transform(world, true);
    if (!usable_mesh(mesh.its))
        return nullptr;
    std::vector<int> patch_triangles(triangles.size());
    std::iota(patch_triangles.begin(), patch_triangles.end(), 0);
    auto model = init_plane_data(mesh.its, patch_triangles, 0.05f);
    if (!model || !finite_bounds(model->get_bounding_box()) ||
        !model->mesh || !usable_mesh(model->mesh->its))
        return nullptr;
    return model;
}

// Preserve all the GL states changed by the transparent world-space overlay.
class OverlayState
{
public:
    OverlayState()
        : blend(glIsEnabled(GL_BLEND)), cull(glIsEnabled(GL_CULL_FACE)),
          depth_test(glIsEnabled(GL_DEPTH_TEST))
    {
        glGetBooleanv(GL_DEPTH_WRITEMASK, &depth_mask);
        glGetIntegerv(GL_FRONT_FACE, &front_face);
        glGetIntegerv(GL_CULL_FACE_MODE, &cull_face);
        glGetIntegerv(GL_BLEND_SRC_RGB, &src_rgb);
        glGetIntegerv(GL_BLEND_DST_RGB, &dst_rgb);
        glGetIntegerv(GL_BLEND_SRC_ALPHA, &src_alpha);
        glGetIntegerv(GL_BLEND_DST_ALPHA, &dst_alpha);
        glGetIntegerv(GL_BLEND_EQUATION_RGB, &equation_rgb);
        glGetIntegerv(GL_BLEND_EQUATION_ALPHA, &equation_alpha);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendEquationSeparate(GL_FUNC_ADD, GL_FUNC_ADD);
        glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        glFrontFace(GL_CCW);
    }

    ~OverlayState()
    {
        glBlendFuncSeparate(src_rgb, dst_rgb, src_alpha, dst_alpha);
        glBlendEquationSeparate(equation_rgb, equation_alpha);
        glFrontFace(front_face);
        glCullFace(cull_face);
        glDepthMask(depth_mask);
        if (blend) glEnable(GL_BLEND); else glDisable(GL_BLEND);
        if (cull) glEnable(GL_CULL_FACE); else glDisable(GL_CULL_FACE);
        if (depth_test) glEnable(GL_DEPTH_TEST); else glDisable(GL_DEPTH_TEST);
    }

private:
    GLboolean blend, cull, depth_test, depth_mask;
    GLint front_face, cull_face, src_rgb, dst_rgb, src_alpha, dst_alpha;
    GLint equation_rgb, equation_alpha;
};

// An inward solid is otherwise entirely hidden by the opaque original. Show
// its edges through the source without modifying that source's visibility.
class WireframeState
{
public:
    WireframeState()
        : depth_test(glIsEnabled(GL_DEPTH_TEST)), cull(glIsEnabled(GL_CULL_FACE))
    {
        glGetIntegerv(GL_POLYGON_MODE, polygon_mode);
        glGetFloatv(GL_LINE_WIDTH, &line_width);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glLineWidth(1.0f);
    }

    ~WireframeState()
    {
        if (polygon_mode[0] == polygon_mode[1])
            glPolygonMode(GL_FRONT_AND_BACK, polygon_mode[0]);
        else {
            // Separate face modes can only occur in a compatibility context.
            glPolygonMode(GL_FRONT, polygon_mode[0]);
            glPolygonMode(GL_BACK, polygon_mode[1]);
        }
        glLineWidth(line_width);
        if (depth_test) glEnable(GL_DEPTH_TEST);
        if (cull) glEnable(GL_CULL_FACE);
    }

private:
    GLboolean depth_test, cull;
    GLint polygon_mode[2];
    GLfloat line_width;
};

class CommitGuard
{
public:
    explicit CommitGuard(bool &flag) : flag(flag) { flag = true; }
    ~CommitGuard() { flag = false; }
private:
    bool &flag;
};

std::optional<double> snapped_arrow_distance(double start_distance, double raw_distance)
{
    if (!std::isfinite(start_distance) || !std::isfinite(raw_distance))
        return std::nullopt;
    // Quantize only the drag DELTA. A custom typed start value must not jump on
    // Down, and the exact unsnapped accumulator must survive camera rebases.
    const double delta = raw_distance - start_distance;
    const double value = start_distance + std::round(delta); // physical 1 mm steps
    if (!std::isfinite(delta) || !std::isfinite(value))
        return std::nullopt;
    return value;
}

std::string arrow_distance_label(double distance)
{
    return Slic3r::string_printf("%.9g mm", distance);
}

std::string boolean_source_limit(const ModelVolume &source)
{
    // Procedural metadata is NEVER covered by annotation-discard consent.
    const ModelObject *object = source.get_object();
    if (object && (object->is_cut() || !object->cut_connectors.empty() ||
                   object->cut_id.check_sum() != 1 || object->cut_id.connectors_cnt() != 0 ||
                   std::any_of(object->volumes.begin(), object->volumes.end(),
                       [](const ModelVolume *part) { return part && part->cut_info.is_connector; })))
        return _u8L("Union/Subtract cannot preserve this object's cut relationships or connectors. Use New part; the original cut object will remain unchanged.");
    if (source.is_text() || source.emboss_shape || source.cut_info.is_connector)
        return _u8L("Union/Subtract is disabled for procedural text, SVG/emboss or connectors. Use New part, or explicitly convert the source to ordinary mesh geometry first.");
    return {};
}

std::string project_exterior_limit(const Model &model)
{
    // Undo captures the ENTIRE Model but omits the dedicated fifth layer. Check
    // every volume, not just the selected source / MODEL_PARTs. Legacy mesh tags
    // are serialized with the mesh and are deliberately NOT blocked here.
    for (const ModelObject *object : model.objects) {
        if (!object) continue;
        for (const ModelVolume *volume : object->volumes) {
            if (!volume) continue;
            const auto &exterior = volume->exterior_facets.get_data();
            if (!exterior.first.empty() || !exterior.second.empty())
                return _u8L("This project contains exterior-facet annotation data that Undo cannot restore. All face-extrusion Apply actions, including New part, are disabled to protect every part. Existing data and selected faces are unchanged. Legacy exterior-appearance mesh tags are separate and supported.");
        }
    }
    return {};
}

std::string annotation_reasons(const FaceExtrusionAnnotationInfo &info)
{
    std::string reasons;
    const auto add = [&](const std::string &reason) {
        if (!reasons.empty()) reasons += "\n";
        reasons += "- " + reason;
    };
    if (info.support) add(_u8L("Support painting data"));
    if (info.seam) add(_u8L("Seam painting data"));
    if (info.mmu) add(_u8L("Color / multi-material painting data"));
    if (info.fuzzy_skin) add(_u8L("Fuzzy-skin painting data"));
    if (info.exterior) add(_u8L("Exterior-facet annotations (Undo unsupported; cannot discard)"));
    if (info.legacy_small_overhang) add(_u8L("Legacy small-overhang tags") + ": " + std::to_string(info.legacy_small_overhang));
    if (info.legacy_small_hole) add(_u8L("Legacy small-hole tags") + ": " + std::to_string(info.legacy_small_hole));
    if (info.legacy_exterior) add(_u8L("Legacy exterior-appearance tags") + ": " + std::to_string(info.legacy_exterior));
    if (info.legacy_other) add(_u8L("Other legacy facet tags") + ": " + std::to_string(info.legacy_other));
    return reasons;
}

struct BooleanCanceled {};

void validate_boolean_mesh(const indexed_triangle_set &mesh, const std::function<bool()> &canceled)
{
    // The same non-mutating solid gate protects Boolean results and NewPart
    // storage. Keep the worker's cancellation type and transactional caller.
    validate_solid_mesh(mesh, [&] { if (canceled()) throw BooleanCanceled{}; });
}

// No gizmo pointer is retained by this callback. Re-resolve the committed IDs
// after scene callbacks, including the one-volume object-tree fallback.
void refresh_boolean_source(ObjectID object_id, ObjectID source_id)
{
    Plater *plater = wxGetApp().plater();
    ObjectList *list = wxGetApp().obj_list();
    if (!plater || !list)
        return;
    int object_index = -1, volume_index = -1;
    const auto locate = [object_id, source_id](int &oi, int &vi) {
        Plater *current = wxGetApp().plater();
        if (!current) return false;
        const auto &objects = current->model().objects;
        for (size_t i = 0; i < objects.size(); ++i) {
            if (objects[i]->id() != object_id) continue;
            for (size_t j = 0; j < objects[i]->volumes.size(); ++j) {
                if (objects[i]->volumes[j]->id() == source_id) {
                    oi = int(i); vi = int(j); return true;
                }
            }
        }
        return false;
    };
    if (!locate(object_index, volume_index))
        return;
    list->update_item_error_icon(object_index, volume_index);
    if (!locate(object_index, volume_index) || !(list = wxGetApp().obj_list()))
        return;
    list->update_info_items(object_index);
    if (!locate(object_index, volume_index) || !(plater = wxGetApp().plater()))
        return;
    plater->changed_object(object_index); // standard FFF refresh; preflight proved no auto-shift
    if (!locate(object_index, volume_index))
        return;
    list = wxGetApp().obj_list();
    plater = wxGetApp().plater();
    if (!list || !plater)
        return;
    if (GLCanvas3D *canvas = plater->get_view3D_canvas3D()) {
        canvas->update_instance_printable_state_for_object(object_index);
        canvas->requires_check_outside_state();
        canvas->request_extra_frame();
    }
    list->select_item([locate] {
        int oi = -1, vi = -1;
        ObjectList *current_list = wxGetApp().obj_list();
        if (!current_list || !locate(oi, vi)) return wxDataViewItem{};
        auto *model = current_list->GetModel();
        const wxDataViewItem part = model->GetItemByVolumeId(oi, vi);
        return part.IsOk() ? part : model->GetItemById(oi);
    });
    if (ObjectList *current_list = wxGetApp().obj_list())
        current_list->selection_changed();
}

// Same callback lifetime boundary for the frozen NewParts path: all created
// IDs are values, and no gizmo/member/ModelVolume pointer survives the reload.
void refresh_created_parts(ObjectID object_id, std::vector<ObjectID> created_ids)
{
    Plater *plater = wxGetApp().plater();
    ObjectList *list = wxGetApp().obj_list();
    if (!plater || !list) return;
    const auto object_index = [object_id] {
        Plater *current = wxGetApp().plater();
        if (current) {
            const auto &objects = current->model().objects;
            for (size_t i = 0; i < objects.size(); ++i)
                if (objects[i]->id() == object_id) return int(i);
        }
        return -1;
    };
    const int index = object_index();
    if (index < 0) return;
    const wxDataViewItemArray selected = list->reorder_volumes_and_get_selection(
        index, [ids = std::move(created_ids)](const ModelVolume *volume) {
            return std::find(ids.begin(), ids.end(), volume->id()) != ids.end();
        });
    list = wxGetApp().obj_list();
    if (!list) return;
    if (!selected.IsEmpty())
        list->select_items(selected);
    if (ObjectList *current_list = wxGetApp().obj_list())
        current_list->selection_changed();
    const int current_index = object_index();
    if ((plater = wxGetApp().plater()) && current_index >= 0) {
        if (GLCanvas3D *canvas = plater->get_view3D_canvas3D()) {
            canvas->update_instance_printable_state_for_object(current_index);
            canvas->request_extra_frame();
        }
    }
    // reorder_volumes_and_get_selection already schedules changed_object.
    // Never explicitly ensure_on_bed(): that would move the original instance.
}

} // namespace

struct GLGizmoFaceExtrude::BooleanLifetime {
    // Accessed only on the GUI thread. Worker callbacks never dereference this.
    bool alive = true;
    uint64_t epoch = 0;
};

struct GLGizmoFaceExtrude::AnnotationState {
    ObjectID source_id;
    std::shared_ptr<const TriangleMesh> mesh;
    FaceExtrusionAnnotationInfo info;
    std::array<uint64_t, 5> timestamps{};
    std::array<std::pair<size_t, size_t>, 5> sizes{};
    std::optional<FaceExtrusionAnnotationStamp> consent; // absent = default false
};

struct GLGizmoFaceExtrude::BooleanRequest : FaceExtrusionPlacementContext {
    uint64_t epoch = 0;
    OutputMode mode = OutputMode::NewParts;
    ToolMode tool_mode = ToolMode::Faces;
    bool prepare_stretch = false;
    std::vector<std::vector<int>> profile_contours;
    double plane_position = 0.0;
    CutPhase cut_phase = CutPhase::Guides;
    PlaneBasis plane_basis = PlaneBasis::XY;
    Vec3d plane_normal = Vec3d::UnitZ(), plane_reference = Vec3d::Zero();
    std::optional<Vec3d> plane_reference_local;
    std::optional<size_t> plane_source_facet;
    std::vector<Vec3d> guide_hits, guide_normals;
    std::vector<std::vector<int>> guide_facets; // exact immutable masks, not only seed facets
    bool cut_guided = true;
    bool moving_opposite = false;
    bool moving_side_chosen = false;
    Vec3d plane_origin = Vec3d::Zero(), moving_normal = Vec3d::Zero();
    std::shared_ptr<const CutAndStretchPlan> stretch_plan;
    std::vector<float> bridge_base_z;
    ObjectID object_id, instance_id;
    std::shared_ptr<const TriangleMesh> source_identity;
    double distance = 0.0;
    std::vector<std::pair<size_t, int>> faces;
    TriangleMesh source_world;
    std::vector<TriangleMesh> cutters;
    FaceExtrusionAnnotationStamp annotation_stamp;
    bool discard_annotations = false;
    std::atomic<bool> cancel{false};
    std::string status, error;
    bool canceled = false;
    bool empty = false;
    std::shared_ptr<const TriangleMesh> result;
    std::shared_ptr<const TriangleMesh> hull;
    BoundingBoxf3 post_bounds;
    double post_min_z = 0.0, post_max_z = 0.0;
};

GLGizmoFaceExtrude::GLGizmoFaceExtrude(GLCanvas3D &parent, unsigned int sprite_id, ToolMode mode)
    : GLGizmoBase(parent, sprite_id), m_tool_mode(mode), m_boolean_lifetime(std::make_shared<BooleanLifetime>())
{
    m_distance = mode == ToolMode::Faces ? 5.0 : 10.0;
    m_output_mode = mode == ToolMode::CutAndStretch ? OutputMode::Union : OutputMode::NewParts;
}

GLGizmoFaceExtrude::~GLGizmoFaceExtrude()
{
    m_boolean_lifetime->alive = false;
    invalidate_boolean();
    m_boolean_lifetime.reset();
    // Do not touch the parent during destruction; the wx window may already be
    // gone. Bindings exist only while this gizmo owns an arrow drag.
    if (m_capture_window) {
        m_capture_window->Unbind(wxEVT_MOUSE_CAPTURE_LOST, &GLGizmoFaceExtrude::on_arrow_capture_lost, this);
        m_capture_window->Unbind(wxEVT_KILL_FOCUS, &GLGizmoFaceExtrude::on_arrow_focus_lost, this);
        if (m_capture_window->HasCapture())
            m_capture_window->ReleaseMouse();
    }
}

bool GLGizmoFaceExtrude::on_init()
{
    m_grabbers.emplace_back();
    m_grabbers.front().enabled = false;
    return true;
}

std::string GLGizmoFaceExtrude::get_icon_filename(bool is_dark_mode) const
{
    switch (m_tool_mode) {
    case ToolMode::PointProfile: return is_dark_mode ? "point_profile_dark.svg" : "point_profile.svg";
    case ToolMode::CutAndStretch: return is_dark_mode ? "cut_stretch_dark.svg" : "cut_stretch.svg";
    default: return is_dark_mode ? "face_extrude_dark.svg" : "face_extrude.svg";
    }
}

std::string GLGizmoFaceExtrude::on_get_name() const
{
    const std::string name = m_tool_mode == ToolMode::PointProfile ? _u8L("Extrude point profile") :
        m_tool_mode == ToolMode::CutAndStretch ? _u8L("Cut and extend") : _u8L("Extrude faces");
    return m_state == Off && !on_is_activable()
        ? name + "\n" + (m_tool_mode == ToolMode::Faces ?
            _u8L("Select a positive part, or one negative part in the Objects list, in Prepare. A negative needs a positive carrier. The object must have exactly one instance.") :
            _u8L("Select a model part in Prepare. The object must have exactly one instance."))
        : name;
}

std::string GLGizmoFaceExtrude::get_tooltip() const
{
    if (boolean_busy())
        return _u8L("Calculating geometry. Drag to orbit; Esc requests cancellation.");
    if (m_arrow_drag)
        return {}; // immediate foreground badge replaces the delayed/capture-suppressed native tooltip
    if (m_arrow_hovered)
        return arrow_parameter_label() + "\n" + (editing_cut_plane() ?
            _u8L("Drag the yellow arrow to move ONLY the cut plane. Esc restores its position.") :
            _u8L("Drag: 1 mm steps. Enter a custom distance in the panel. Esc cancels the drag."));
    if (m_tool_mode == ToolMode::PointProfile)
        return _u8L("Click mesh vertices in boundary order. The first point locks the source part. Drag elsewhere to orbit.");
    if (m_tool_mode == ToolMode::CutAndStretch && m_cut_phase == CutPhase::Guides && m_smart_guide_fill)
        return _u8L("Smart fill: hover to preview the native angle region. Click adds missing guide triangles, or removes the region when all are selected. No painting or geometry change.");
    if (m_tool_mode == ToolMode::CutAndStretch)
        return m_cut_phase == CutPhase::Guides ? _u8L("Click to add/remove guide faces of the portion to move, even with different normals. Then continue to the independent plane.") :
            m_cut_phase == CutPhase::Plane ? _u8L("Position the yellow plane by clicking the source, typing its absolute coordinate or dragging the yellow arrow. Confirm before moving geometry.") :
            _u8L("The entire connected portion indicated by the guides moves normally to the cut plane. Preview first; only Apply changes the source.");
    return _u8L("Click a planar face; Shift + click adds or removes faces. Drag elsewhere to orbit.");
}

BoundingBoxf3 GLGizmoFaceExtrude::get_bounding_box() const
{
    BoundingBoxf3 bounds;
    if (m_state == On) {
        // The canvas uses this WORLD-space box to choose camera near/far planes.
        // Without it a long extrusion can be clipped outside the source object.
        if (m_preview_valid) {
            for (const PreviewPart &part : m_previews) {
                const BoundingBoxf3 preview_bounds = part.render_bounds.transformed(part.render_world);
                if (finite_bounds(preview_bounds))
                    bounds.merge(preview_bounds);
            }
        }
        if (m_face_model) {
            const BoundingBoxf3 face_bounds = m_tool_mode == ToolMode::CutAndStretch ?
                translated_face_tools_preview_bounds(m_face_model->get_bounding_box(), m_plane_render_matrix.translation()) : m_face_model->get_bounding_box();
            if (finite_bounds(face_bounds))
                bounds.merge(face_bounds);
        }
        if (m_hover_model) {
            const BoundingBoxf3 &hover_bounds = m_hover_model->get_bounding_box();
            if (finite_bounds(hover_bounds))
                bounds.merge(hover_bounds);
        }
        if (!m_hidden_source)
            for (const auto &guide : m_guide_models)
                if (guide && finite_bounds(guide->get_bounding_box())) bounds.merge(guide->get_bounding_box());
        if (m_profile_lines && finite_bounds(m_profile_lines->get_bounding_box()))
            bounds.merge(m_profile_lines->get_bounding_box());
        if (m_tool_mode == ToolMode::PointProfile) {
            const double radius = 0.28 * get_grabber_size();
            if (std::isfinite(radius) && radius > 0.0) {
                auto include_marker = [&](const Vec3d &point) {
                    if (point.allFinite()) {
                        bounds.merge(Vec3d(point - Vec3d::Constant(radius)));
                        bounds.merge(Vec3d(point + Vec3d::Constant(radius)));
                    }
                };
                for (const auto &contour : m_profile_world)
                    for (const Vec3d &point : contour) include_marker(point);
                if (m_hovered && m_hovered->vertex_id >= 0)
                    include_marker(m_hovered->world_at_pick * m_hovered->hit_local);
            }
        }
        if (m_arrow_available && finite_bounds(m_arrow_bounds))
            bounds.merge(m_arrow_bounds);
    }
    return bounds;
}

bool GLGizmoFaceExtrude::on_is_selectable() const
{
    const Plater *plater = wxGetApp().plater();
    return plater && plater->printer_technology() == ptFFF &&
           m_parent.get_canvas_type() == GLCanvas3D::CanvasView3D;
}

bool GLGizmoFaceExtrude::get_target(ModelObject *&object, ModelInstance *&instance,
                                   int *object_index, ModelVolume **negative_source) const
{
    object = nullptr;
    instance = nullptr;
    if (negative_source) *negative_source = nullptr;
    const Plater *plater = wxGetApp().plater();
    const Selection &selection = m_parent.get_selection();
    if (!plater || !on_is_selectable() || !plater->is_view3D_shown() ||
        !std::isfinite(m_parent.get_explosion_ratio()) ||
        std::abs(m_parent.get_explosion_ratio() - 1.0f) > 1e-6f ||
        selection.is_empty() || selection.is_wipe_tower() ||
        !selection.is_from_single_object() || !selection.is_from_single_instance())
        return false;

    const Model *model = selection.get_model();
    const int index = selection.get_object_idx();
    const int instance_index = selection.get_instance_idx();
    if (!model || model != m_parent.get_model() ||
        model != &plater->model() || index < 0 || size_t(index) >= model->objects.size())
        return false;
    ModelObject *candidate = model->objects[index];
    if (!candidate || candidate->instances.size() != 1 || instance_index != 0)
        return false;

    // SingleModifier includes negative volumes AND unsupported modifiers. Resolve
    // real Model roles, not the broad GL flag. A negative must be an explicit
    // singleton, never an implicit child of a selected parent/mixed selection.
    // These GLVolume pointers are used only during this call, never retained.
    ModelVolume *negative = nullptr;
    std::vector<ModelVolumeType> selected_roles;
    selected_roles.reserve(selection.volumes_count());
    for (unsigned int selected_index : selection.get_volume_idxs()) {
        const GLVolume *selected = selection.get_volume(selected_index);
        if (!selected || selected->object_idx() != index || selected->volume_idx() < 0 ||
            size_t(selected->volume_idx()) >= candidate->volumes.size())
            return false;
        ModelVolume *volume = candidate->volumes[selected->volume_idx()];
        selected_roles.push_back(volume->type());
        if (volume->is_negative_volume()) {
            // A reload/reorder must not reinterpret an old GL index as another
            // negative source before its scene identities have caught up.
            if (selected->geometry_id.first != volume->id().id ||
                selected->geometry_id.second != candidate->instances.front()->id().id)
                return false;
            negative = volume;
        }
    }
    if (!face_extrusion_selection_roles_allowed(m_tool_mode, selected_roles) ||
        !usable_transform(candidate->instances.front()->get_matrix()) ||
        (negative && (m_tool_mode != ToolMode::Faces || selected_roles.size() != 1 || candidate->parts_count() == 0)))
        return false;
    if (negative_source) *negative_source = negative;
    object = candidate;
    instance = candidate->instances.front();
    if (object_index)
        *object_index = index;
    return true;
}

bool GLGizmoFaceExtrude::on_is_activable() const
{
    ModelObject *object;
    ModelInstance *instance;
    return get_target(object, instance);
}

CommonGizmosDataID GLGizmoFaceExtrude::on_get_requirements() const
{
    // Positive tools retain the shared MODEL_PART-only pool. The negative path
    // owns one private raycaster; invalid/zero-carrier targets must not reach the
    // common pool's single-MODEL_PART assertion during a scene/type change.
    ModelObject *object;
    ModelInstance *instance;
    ModelVolume *negative = nullptr;
    int requirements = int(CommonGizmosDataID::SelectionInfo) | int(CommonGizmosDataID::InstancesHider);
    if (get_target(object, instance, nullptr, &negative) && !negative)
        requirements |= int(CommonGizmosDataID::Raycaster);
    return CommonGizmosDataID(requirements);
}

void GLGizmoFaceExtrude::restore_negative_visibility()
{
    if (!m_negative_visibility) return;
    const NegativeVisibility saved = std::move(*m_negative_visibility);
    m_negative_visibility.reset();
    m_parent.set_clipping_plane(0, saved.lower_clip);
    const Model *model = m_parent.get_model();
    if (model) {
        for (const ModelObject *object : model->objects) {
            if (object->id() != saved.object_id) continue;
            // A common-pool refresh owns visibility of a different active object.
            // Never re-show an old object over its newly established isolation.
            const auto *info = m_c ? m_c->selection_info() : nullptr;
            if (info && info->model_object() != object) break;
            for (GLVolume *volume : m_parent.get_volumes().volumes) {
                if (volume->volume_idx() < 0 || volume->geometry_id.second != saved.instance_id.id) continue;
                const auto previous = std::find_if(saved.previous_visibility.begin(), saved.previous_visibility.end(),
                    [volume](const auto &entry) { return entry.first.id == volume->geometry_id.first; });
                if (previous != saved.previous_visibility.end()) volume->is_active = previous->second;
            }
            break;
        }
    }
    if (m_state == On && !m_committing) m_parent.request_extra_frame();
}

void GLGizmoFaceExtrude::update_negative_visibility()
{
    ModelObject *object;
    ModelInstance *instance;
    ModelVolume *negative = nullptr;
    int object_index = -1;
    const auto *info = m_c ? m_c->selection_info() : nullptr;
    if (m_state != On || m_committing || m_tool_mode != ToolMode::Faces ||
        !get_target(object, instance, &object_index, &negative) || !negative ||
        !info || info->model_object() != object || info->get_active_instance() != 0) {
        restore_negative_visibility();
        return;
    }
    if (m_negative_visibility && m_negative_visibility->object_id == object->id() &&
        m_negative_visibility->instance_id == instance->id() && m_negative_visibility->volume_id == negative->id())
        return;
    restore_negative_visibility();
    NegativeVisibility saved;
    saved.object_id = object->id();
    saved.instance_id = instance->id();
    saved.volume_id = negative->id();
    saved.lower_clip = m_parent.get_clipping_planes()[0];
    for (const GLVolume *volume : m_parent.get_volumes().volumes) {
        if (volume->object_idx() == object_index && volume->volume_idx() >= 0 && volume->geometry_id.second == instance->id().id)
            saved.previous_visibility.emplace_back(ObjectID(volume->geometry_id.first), volume->is_active);
    }
    m_negative_visibility = std::move(saved);
    // Match stable scene geometry IDs, not GL pointers or indices retained over
    // reloads. Change only visibility; the common hider owns color/view policy.
    for (GLVolume *volume : m_parent.get_volumes().volumes)
        if (volume->object_idx() == object_index && volume->volume_idx() >= 0 && volume->geometry_id.second == instance->id().id)
            volume->is_active = volume->geometry_id.first == negative->id().id;
    // A cutter may cross the bed. Remove only the native lower-bed clip while
    // it is isolated; retain camera clipping and restore the exact saved plane.
    // Canvas uses z_min = -plane[0].offset, not the plane's signed distance.
    m_parent.set_clipping_plane(0, ClippingPlane(Vec3d::UnitZ(), double(std::numeric_limits<float>::max())));
    m_parent.request_extra_frame();
}

void GLGizmoFaceExtrude::clear_hover(bool invalidate_query)
{
    const bool was_visible = bool(m_hover_model);
    m_hovered.reset();
    m_hover_model.reset();
    m_hover_region.reset();
    m_guide_hover_error.clear();
    if (invalidate_query)
        m_hover_query.reset();
    if (was_visible && m_state == On && !m_committing)
        m_parent.request_extra_frame();
}

void GLGizmoFaceExtrude::update_annotation_state(const ModelVolume &source, bool refresh)
{
    if (!m_annotations || m_annotations->source_id != source.id() || m_annotations->mesh != source.get_mesh_shared_ptr()) {
        m_annotations = std::make_unique<AnnotationState>();
        m_annotations->source_id = source.id();
        m_annotations->mesh = source.get_mesh_shared_ptr();
        refresh = true;
    }
    const std::array<const FacetsAnnotation *, 5> layers{{&source.supported_facets, &source.seam_facets,
        &source.mmu_segmentation_facets, &source.fuzzy_skin_facets, &source.exterior_facets}};
    std::array<uint64_t, 5> timestamps;
    std::array<std::pair<size_t, size_t>, 5> sizes;
    for (size_t i = 0; i < layers.size(); ++i) {
        timestamps[i] = layers[i]->timestamp();
        const auto &data = layers[i]->get_data();
        sizes[i] = {data.first.size(), data.second.size()};
    }
    const bool changed = timestamps != m_annotations->timestamps || sizes != m_annotations->sizes;
    if (changed) {
        // Cheap cache invalidation only, NEVER a substitute for exact payload
        // comparison at consent/start/finish (some setters do not touch time).
        m_annotations->consent.reset();
        invalidate_boolean();
    }
    if (refresh || changed) {
        m_annotations->info = inspect_face_extrusion_annotations(source);
        m_annotations->timestamps = timestamps;
        m_annotations->sizes = sizes;
    }
}

void GLGizmoFaceExtrude::invalidate_boolean()
{
    if (m_boolean_request)
        m_boolean_request->cancel.store(true, std::memory_order_relaxed);
    m_boolean_request.reset();
    if (m_boolean_lifetime)
        ++m_boolean_lifetime->epoch;
    m_boolean_error.clear();
}

void GLGizmoFaceExtrude::restore_source_visibility()
{
    if (!m_hidden_source) return;
    if (const Model *model = m_parent.get_model())
        for (const ModelObject *object : model->objects)
            if (object && object->id() == m_hidden_source->first)
                for (const ModelVolume *volume : object->volumes)
                    if (volume && volume->id() == m_hidden_source->second)
                        m_parent.toggle_model_objects_visibility(true, object, 0, volume);
    m_hidden_source.reset(); // no GLVolume/model pointers survive scene changes
}

void GLGizmoFaceExtrude::update_stretch_visibility()
{
    const bool show = m_tool_mode == ToolMode::CutAndStretch && m_cut_phase == CutPhase::Distance && m_stretch_plan && m_previews.size() == 3 &&
        std::isfinite(m_distance) && m_distance >= 0.0 && (m_preview_valid || m_distance == 0.0);
    ModelObject *object; ModelInstance *instance; ModelVolume *source;
    if (!show || !resolve_source(object, instance, source)) {
        restore_source_visibility();
        return;
    }
    m_parent.toggle_model_objects_visibility(false, object, 0, source);
    m_hidden_source = std::make_pair(object->id(), source->id());
}

void GLGizmoFaceExtrude::invalidate_preview()
{
    restore_source_visibility();
    invalidate_boolean();
    m_profile_valid = false;
    m_profile_world.clear();
    m_profile_lines.reset();
    m_stretch_plan.reset();
    m_cut_scope_confirmation.clear();
    m_rendered_stretch_plan.reset();
    m_preview_bridge = indexed_triangle_set();
    m_bridge_base_z.clear();
    m_bridge_distance.reset();
    m_prepare_requested = false; // only explicit plane confirmation/retry can prepare a cut
    m_guide_models.clear();
    m_world_source = TriangleMesh();
    m_previews.clear();
    m_face_model.reset();
    m_arrow_available = false;
    m_arrow_shaft_available = false;
    m_arrow_hovered = false;
    m_arrow_bounds.reset();
    m_handle_distance = 0.0;
    if (!m_grabbers.empty())
        m_grabbers.front().enabled = false;
    m_world_dirty = true;
    m_preview_dirty = true;
    m_preview_valid = false;
    m_error.clear();
}

void GLGizmoFaceExtrude::clear_selection()
{
    finish_arrow_drag(false);
    clear_hover();
    restore_source_visibility();
    m_selected.clear();
    m_exact_guide_selection = false;
    m_profile_contours.assign(1, {});
    m_cut_phase = CutPhase::Guides;
    m_cut_guided = true;
    m_cut_scope_confirmation.clear();
    m_plane_basis = PlaneBasis::XY;
    m_plane_initialized = false;
    m_face_plane_ready = true;
    m_plane_position = 0.0;
    m_plane_normal_world = Vec3d::UnitZ();
    m_plane_reference_world = Vec3d::Zero();
    m_plane_reference_local.reset();
    m_plane_source_facet.reset();
    m_cut_guide_index = 0;
    m_plane_render_matrix = Transform3d::Identity();
    m_moving_opposite = false;
    m_moving_side_chosen = false;
    m_prepare_requested = false;
    m_annotations.reset();
    m_selection_error.clear();
    invalidate_preview();
}

void GLGizmoFaceExtrude::on_set_state()
{
    restore_negative_visibility();
    clear_selection();
    m_negative_raycaster.reset();
    m_measurements.clear();
    m_raycaster_meshes.clear();
    m_left_gesture = LeftGesture::None;
    m_consume_left_up = false;
    if (m_state == On) {
        m_smart_guide_fill = false;
        m_smart_guide_angle = 30.0f;
        m_output_mode = m_tool_mode == ToolMode::CutAndStretch ? OutputMode::Union : OutputMode::NewParts;
        if (m_c)
            m_c->raycaster_ptr()->set_only_support_model_part_flag(true);
    }
    GLGizmoBase::on_set_state();
}

void GLGizmoFaceExtrude::on_load(cereal::BinaryInputArchive &)
{
    // Preview geometry is intentionally not serialized into undo snapshots.
    data_changed(true);
}

void GLGizmoFaceExtrude::data_changed(bool /*is_serializing*/)
{
    if (m_committing)
        return;
    const bool orbiting = m_left_gesture == LeftGesture::Orbit;
    const bool owned_press = m_left_gesture != LeftGesture::None || m_consume_left_up;
    restore_negative_visibility();
    clear_selection();
    m_negative_raycaster.reset();
    m_measurements.clear();
    m_raycaster_meshes.clear();
    // Async completion may reload the scene WHILE the user is orbiting. Keep
    // its release owned; never turn that Up into plate/object selection.
    m_left_gesture = orbiting ? LeftGesture::Orbit : owned_press ? LeftGesture::Cancelled : LeftGesture::None;
    m_consume_left_up = owned_press;
    // Common data has normally just re-established instance visibility. Apply
    // the narrower negative isolation afterwards, including Undo/scene reloads.
    update_negative_visibility();
    if (m_state == On)
        m_parent.request_extra_frame();
}

Measure::Measuring *GLGizmoFaceExtrude::get_measuring(const ModelVolume &volume, std::string &error)
{
    MeasurementCache &cache = m_measurements[volume.id()];
    if (cache.mesh != volume.get_mesh_shared_ptr()) {
        cache.guide_selection.reset();
        cache.measuring.reset();
        cache.mesh = volume.get_mesh_shared_ptr();
        cache.error.clear();
        cache.guide_error.clear();
    }
    if (!cache.measuring && cache.error.empty()) {
        // Cache failed validation too: hovering an invalid source must not
        // repeatedly scan the complete mesh. The immutable mesh is its revision.
        try {
            // Both native measurement and angular fill consume float normals.
            // Bound Cut work and reject unusable arithmetic BEFORE either native constructor.
            if (m_tool_mode == ToolMode::CutAndStretch)
                CutGuideSelection::validate_source(*cache.mesh);
            if (!usable_mesh(cache.mesh->its))
                cache.error = _u8L("The source mesh contains invalid or degenerate triangles.");
            else
                cache.measuring = std::make_unique<Measure::Measuring>(cache.mesh->its);
        } catch (const std::exception &e) {
            cache.error = m_tool_mode == ToolMode::CutAndStretch ? e.what() :
                _u8L("Unable to identify a planar face in this mesh.");
        }
    }
    if (!cache.measuring)
        error = cache.error;
    return cache.measuring.get();
}

std::optional<GLGizmoFaceExtrude::SelectedFace> GLGizmoFaceExtrude::pick_face(
    const Vec2d &mouse_position, std::string &error, bool surface_only)
{
    // Shared by hover and click, but does NOT write the locked selection,
    // extrusion preview or UI error. The caller decides what to do with a hit.
    try {
        ModelObject *object;
        ModelInstance *instance;
        ModelVolume *negative = nullptr;
        if (!get_target(object, instance, nullptr, &negative) || !m_c || !mouse_position.allFinite())
            return std::nullopt;
        auto *info = m_c->selection_info();
        if (!info || info->model_object() != object || info->get_active_instance() != 0)
            return std::nullopt;
        update_negative_visibility();

        std::vector<ModelVolume *> parts;
        std::vector<const MeshRaycaster *> raycasters;
        if (negative) {
            // Never raycast the positive carrier in front of the selected cutter,
            // nor enable modifiers/other negatives in the shared positive pool.
            if (!m_negative_raycaster || m_negative_raycaster->object_id != object->id() ||
                m_negative_raycaster->instance_id != instance->id() || m_negative_raycaster->volume_id != negative->id() ||
                m_negative_raycaster->mesh != negative->get_mesh_shared_ptr()) {
                m_negative_raycaster.reset();
                if (!get_measuring(*negative, error)) return std::nullopt; // validate BEFORE native construction
                NegativeRaycasterCache cache;
                cache.object_id = object->id();
                cache.instance_id = instance->id();
                cache.volume_id = negative->id();
                cache.mesh = negative->get_mesh_shared_ptr();
                cache.raycaster = std::make_unique<MeshRaycaster>(*cache.mesh);
                m_negative_raycaster = std::move(cache);
            }
            parts.push_back(negative);
            raycasters.push_back(m_negative_raycaster->raycaster.get());
        } else {
            m_negative_raycaster.reset();
            auto *raycaster = m_c->raycaster();
            if (!raycaster) return std::nullopt;
            // Unchanged positive path: indexed by MODEL_PART order, not volumes.
            std::vector<std::shared_ptr<const TriangleMesh>> meshes;
            for (ModelVolume *volume : object->volumes) {
                if (volume->is_model_part()) {
                    parts.push_back(volume);
                    meshes.push_back(volume->get_mesh_shared_ptr());
                }
            }
            // Camera/pose changes need a new query, not new mesh raycasters.
            if (meshes != m_raycaster_meshes) {
                raycaster->set_only_support_model_part_flag(true);
                raycaster->update();
                m_raycaster_meshes = std::move(meshes);
            }
            raycasters = raycaster->raycasters();
        }
        if (raycasters.size() != parts.size()) {
            error = _u8L("The object changed. Select the face again.");
            m_raycaster_meshes.clear();
            return std::nullopt;
        }

        const Camera &camera = m_parent.get_active_camera();
        ModelVolume *closest_volume = nullptr;
        Vec3d closest_hit = Vec3d::Zero();
        size_t closest_facet = 0;
        Transform3d closest_world = Transform3d::Identity();
        double closest_distance = std::numeric_limits<double>::max();
        for (size_t i = 0; i < parts.size(); ++i) {
            const Transform3d world = instance->get_matrix() * parts[i]->get_matrix();
            if (!usable_transform(world)) {
                error = _u8L("The object has an invalid or singular transformation.");
                return std::nullopt;
            }
            Vec3f hit = Vec3f::Zero();
            Vec3f normal = Vec3f::Zero();
            size_t facet = 0;
            if (!raycasters[i]->unproject_on_mesh(mouse_position, world, camera, hit, normal, nullptr, &facet, negative == nullptr) ||
                !hit.allFinite() || !normal.allFinite() || facet >= parts[i]->mesh().its.indices.size())
                continue;
            const Vec3d hit_world = world * hit.cast<double>();
            const double distance = (camera.get_position() - hit_world).squaredNorm();
            if (hit_world.allFinite() && std::isfinite(distance) && distance < closest_distance) {
                closest_distance = distance;
                closest_volume = parts[i];
                closest_facet = facet;
                closest_hit = hit.cast<double>();
                closest_world = world;
            }
        }
        if (!closest_volume)
            return std::nullopt;

        if (m_tool_mode == ToolMode::PointProfile) {
            const auto &mesh = closest_volume->mesh().its;
            const Vec3i &triangle = mesh.indices[closest_facet];
            int vertex = -1;
            double best = std::pow(14.0 * m_parent.get_scale(), 2);
            for (int i = 0; i < 3; ++i) {
                if (triangle[i] < 0 || size_t(triangle[i]) >= mesh.vertices.size())
                    return std::nullopt;
                const auto screen = profile_screen_point(closest_world * mesh.vertices[triangle[i]].cast<double>(), camera);
                if (screen && (*screen - mouse_position).squaredNorm() <= best) {
                    best = (*screen - mouse_position).squaredNorm();
                    vertex = triangle[i];
                }
            }
            if (vertex < 0) {
                error = _u8L("Click near a visible mesh vertex; the profile uses exact source vertices, not arbitrary surface points.");
                return std::nullopt;
            }
            SelectedFace point;
            point.object_id = object->id();
            point.instance_id = instance->id();
            point.volume_id = closest_volume->id();
            point.type_at_pick = closest_volume->type();
            point.mesh = closest_volume->get_mesh_shared_ptr();
            point.world_at_pick = closest_world;
            point.hit_local = mesh.vertices[vertex].cast<double>();
            point.vertex_id = vertex;
            point.seed_facet = closest_facet;
            return point; // profile normal is computed only from the validated ordered contour
        }
        if (surface_only) {
            SelectedFace point;
            point.object_id = object->id(); point.instance_id = instance->id();
            point.volume_id = closest_volume->id(); point.type_at_pick = closest_volume->type();
            point.mesh = closest_volume->get_mesh_shared_ptr();
            point.world_at_pick = closest_world; point.hit_local = closest_hit; point.seed_facet = closest_facet;
            return point; // independent plane placement may click any source surface, not only planar patches
        }
        Measure::Measuring *measuring = get_measuring(*closest_volume, error);
        if (!measuring)
            return std::nullopt;
        // Keep the feature LOCAL: SurfaceFeature::translate is not suitable for
        // normals under non-uniform scale. Also never use the plane's centroid
        // as a seed, since that point may lie inside a hole in the planar patch.
        const auto feature = measuring->get_feature(closest_facet, closest_hit, Transform3d::Identity(), true);
        if (!feature || feature->get_type() != Measure::SurfaceFeatureType::Plane) {
            error = _u8L("Select a planar face.");
            return std::nullopt;
        }
        const auto plane = feature->get_plane();
        const int plane_id = std::get<0>(plane);
        if (plane_id < 0 || plane_id >= measuring->get_num_of_planes())
            return std::nullopt;
        SelectedFace face;
        face.object_id = object->id();
        face.instance_id = instance->id();
        face.volume_id = closest_volume->id();
        face.type_at_pick = closest_volume->type();
        face.mesh = closest_volume->get_mesh_shared_ptr();
        face.plane_id = plane_id;
        face.world_at_pick = closest_world;
        face.normal_local = std::get<1>(plane);
        if (m_tool_mode == ToolMode::CutAndStretch) {
            // Measure's plane normal is float-derived. The explicit facet-plane
            // kernel needs the original seed triangle's DOUBLE geometric normal.
            const auto &mesh = face.mesh->its;
            const Vec3i &triangle = mesh.indices.at(closest_facet);
            const Vec3d a = mesh.vertices.at(size_t(triangle[0])).cast<double>();
            const Vec3d b = mesh.vertices.at(size_t(triangle[1])).cast<double>();
            const Vec3d c = mesh.vertices.at(size_t(triangle[2])).cast<double>();
            face.normal_local = (b - a).cross(c - a);
        }
        if (!face.normal_local.allFinite() || face.normal_local.norm() <= 0.0) {
            error = _u8L("The selected face has an invalid normal.");
            return std::nullopt;
        }
        face.normal_local.normalize();
        face.hit_local = closest_hit;
        face.seed_facet = closest_facet;
        // Triangle lists are copied only on a plane transition or a click, not
        // for every pixel movement within the same plane.
        return face;
    } catch (const std::exception &) {
        error = _u8L("Unable to identify a planar face in this mesh.");
        return std::nullopt;
    }
}

bool GLGizmoFaceExtrude::load_face_triangles(SelectedFace &face, std::string &error) const
{
    const auto entry = m_measurements.find(face.volume_id);
    if (entry == m_measurements.end() || entry->second.mesh != face.mesh ||
        !entry->second.measuring || face.plane_id < 0 ||
        face.plane_id >= entry->second.measuring->get_num_of_planes()) {
        error = _u8L("The object changed. Select the face again.");
        return false;
    }
    face.triangles = entry->second.measuring->get_plane_triangle_indices(face.plane_id);
    std::sort(face.triangles.begin(), face.triangles.end());
    face.triangles.erase(std::unique(face.triangles.begin(), face.triangles.end()), face.triangles.end());
    if (face.triangles.empty() || face.triangles.front() < 0 ||
        size_t(face.triangles.back()) >= face.mesh->its.indices.size() ||
        !std::binary_search(face.triangles.begin(), face.triangles.end(), int(face.seed_facet))) {
        error = _u8L("Select a valid planar face.");
        return false;
    }
    return true;
}

bool GLGizmoFaceExtrude::same_source(const SelectedFace &a, const SelectedFace &b)
{
    return a.object_id == b.object_id && a.instance_id == b.instance_id &&
           a.volume_id == b.volume_id && a.type_at_pick == b.type_at_pick && a.mesh == b.mesh;
}

bool GLGizmoFaceExtrude::same_face(const SelectedFace &a, const SelectedFace &b)
{
    return same_source(a, b) && a.plane_id == b.plane_id && a.vertex_id == b.vertex_id;
}

CutGuideSelection *GLGizmoFaceExtrude::get_cut_guide_selection(const SelectedFace &face, std::string &error)
{
    const auto found = m_measurements.find(face.volume_id);
    if (found == m_measurements.end() || found->second.mesh != face.mesh || !found->second.measuring) {
        error = _u8L("The source changed. Select the guide region again.");
        return nullptr;
    }
    auto &cache = found->second;
    if (!cache.guide_selection && cache.guide_error.empty()) {
        // Before a guide locks the source, hovering many parts must not retain
        // a heavy native TriangleSelector for every previously visited volume.
        for (auto &entry : m_measurements)
            if (entry.first != face.volume_id) entry.second.guide_selection.reset();
        try {
            CutGuideSelection::validate_source(*face.mesh);
            const size_t count = face.mesh->its.indices.size();
            const int planes = cache.measuring->get_num_of_planes();
            if (planes <= 0 || size_t(planes) > count)
                throw std::runtime_error(_u8L("The source has an invalid measurement-plane map."));
            std::vector<int> plane_ids(count, -1);
            for (int plane = 0; plane < planes; ++plane) {
                for (int triangle : cache.measuring->get_plane_triangle_indices(plane)) {
                    if (triangle < 0 || size_t(triangle) >= count ||
                        (plane_ids[triangle] != -1 && plane_ids[triangle] != plane))
                        throw std::runtime_error(_u8L("The source has an inconsistent measurement-plane map."));
                    plane_ids[triangle] = plane; // duplicate entries in the SAME native plane do not expand the mask
                }
            }
            cache.guide_selection = std::make_unique<CutGuideSelection>(face.mesh, std::move(plane_ids));
        } catch (const std::exception &e) {
            cache.guide_error = e.what(); // cache failed source validation too, not a per-pixel retry
        }
    }
    if (!cache.guide_selection) error = cache.guide_error;
    return cache.guide_selection.get();
}

std::vector<CutGuidePatch> GLGizmoFaceExtrude::selected_cut_patches() const
{
    std::vector<CutGuidePatch> patches;
    patches.reserve(m_selected.size());
    for (const auto &face : m_selected) {
        CutGuidePatch patch;
        patch.plane_id = face.plane_id;
        patch.triangles = face.triangles; // click/region-transition only; canonical legacy masks are already unique
        patch.seed_facet = face.seed_facet;
        patch.hit_local = face.hit_local;
        patch.normal_local = face.normal_local;
        patches.push_back(std::move(patch));
    }
    return patches;
}

void GLGizmoFaceExtrude::select_cut_region(SelectedFace face)
{
    if (m_tool_mode != ToolMode::CutAndStretch || m_cut_phase != CutPhase::Guides || boolean_busy()) return;
    if (std::any_of(m_selected.begin(), m_selected.end(),
                    [&face](const SelectedFace &entry) { return !same_source(entry, face); })) {
        m_selection_error = _u8L("Select guide regions from the same locked source part.");
        return;
    }
    std::string error;
    auto *selector = get_cut_guide_selection(face, error);
    if (!selector) { m_selection_error = std::move(error); return; }
    try {
        std::shared_ptr<const std::vector<int>> angular;
        if (m_smart_guide_fill)
            angular = selector->angular_region(face.seed_facet, face.hit_local, m_smart_guide_angle);
        else if (!load_face_triangles(face, error)) {
            m_selection_error = std::move(error);
            return;
        }
        auto patches = selector->toggle(selected_cut_patches(), angular ? *angular : face.triangles,
                                        face.seed_facet, face.hit_local);
        std::vector<SelectedFace> next;
        next.reserve(patches.size());
        for (auto &patch : patches) {
            if (!world_normal(patch.normal_local, face.world_at_pick) ||
                !(face.world_at_pick * patch.hit_local).allFinite())
                throw std::runtime_error(_u8L("A guide region has an invalid normal or world position."));
            SelectedFace selected;
            selected.object_id = face.object_id; selected.instance_id = face.instance_id;
            selected.volume_id = face.volume_id; selected.mesh = face.mesh;
            selected.world_at_pick = face.world_at_pick;
            selected.plane_id = patch.plane_id; selected.seed_facet = patch.seed_facet;
            selected.hit_local = patch.hit_local; selected.normal_local = patch.normal_local;
            selected.triangles = std::move(patch.triangles);
            next.push_back(std::move(selected));
        }
        // Keep the explicit guide-combo identity by its real seed, not its old
        // vector index. For the first angular batch, default to the CLICKED guide.
        const size_t preferred_seed = m_selected.empty() ? face.seed_facet :
            m_selected[std::min(m_cut_guide_index, m_selected.size() - 1)].seed_facet;
        auto preferred = std::find_if(next.begin(), next.end(), [preferred_seed](const SelectedFace &guide) {
            return std::binary_search(guide.triangles.begin(), guide.triangles.end(), int(preferred_seed));
        });
        if (preferred == next.end())
            preferred = std::find_if(next.begin(), next.end(), [&face](const SelectedFace &guide) {
                return std::binary_search(guide.triangles.begin(), guide.triangles.end(), int(face.seed_facet));
            });
        const size_t next_guide = preferred == next.end() ? 0 : size_t(preferred - next.begin());
        // Publish ONE validated batch. Invalid masks, resource caps or new metadata
        // never erase a subset of the existing guides or touch Model annotations.
        if (next.empty()) clear_selection();
        else {
            m_selected.swap(next);
            m_exact_guide_selection = true;
            m_cut_guide_index = next_guide;
            clear_hover();
            m_selection_error.clear();
            invalidate_preview();
        }
        request_preview_redraw();
    } catch (const std::exception &e) {
        m_selection_error = e.what();
    }
}

void GLGizmoFaceExtrude::update_cut_region_hover(SelectedFace face)
{
    if (std::any_of(m_selected.begin(), m_selected.end(),
                    [&face](const SelectedFace &entry) { return !same_source(entry, face); })) {
        clear_hover(false);
        return;
    }
    if (!m_smart_guide_fill && m_hovered && same_face(*m_hovered, face) &&
        (m_hovered->world_at_pick.matrix().array() == face.world_at_pick.matrix().array()).all())
        return;
    std::string error;
    auto *selector = get_cut_guide_selection(face, error);
    if (!selector) {
        clear_hover(false);
        m_guide_hover_error = std::move(error);
        return;
    }
    std::shared_ptr<const std::vector<int>> region;
    try {
        if (m_smart_guide_fill)
            region = selector->angular_region(face.seed_facet, face.hit_local, m_smart_guide_angle);
        else {
            if (!load_face_triangles(face, error)) throw std::runtime_error(error);
            region = std::make_shared<const std::vector<int>>(std::move(face.triangles));
        }
        if (m_hover_region == region && m_hovered && same_source(*m_hovered, face) &&
            (m_hovered->world_at_pick.matrix().array() == face.world_at_pick.matrix().array()).all())
            return; // same native region: no mask copy, repartition or GL rebuild on every pixel
        clear_hover(false);
        m_hovered = face;
        m_hover_region = region;
        // Preflight the exact same batch/caps as a click. In particular, a smooth
        // flood that would exceed the guide limit must not build a huge misleading overlay.
        (void)selector->toggle(selected_cut_patches(), *region, face.seed_facet, face.hit_local);
        std::vector<unsigned char> selected(face.mesh->its.indices.size(), 0);
        for (const auto &guide : m_selected)
            for (int triangle : guide.triangles)
                selected.at(size_t(triangle)) = 1;
        std::vector<int> visible;
        visible.reserve(region->size());
        for (int triangle : *region)
            if (!selected.at(size_t(triangle))) visible.push_back(triangle);
        // Existing selected triangles remain amber; green marks only the missing
        // portion, without drawing a second coincident overlay over selected material.
        if (!visible.empty()) {
            m_hover_model = make_hover_model(face.mesh->its, visible, face.world_at_pick);
            if (!m_hover_model) throw std::runtime_error(_u8L("Unable to highlight this guide region."));
            m_parent.request_extra_frame();
        }
    } catch (const std::exception &e) {
        clear_hover(false);
        m_hovered = std::move(face);
        m_hover_region = std::move(region); // memoize a rejected region until its input or selection changes
        m_guide_hover_error = e.what();
    }
}

void GLGizmoFaceExtrude::select_face(const Vec2d &mouse_position, bool additive)
{
    if (boolean_busy())
        return;
    if (m_tool_mode == ToolMode::PointProfile) {
        select_profile_vertex(mouse_position);
        return;
    }
    if (m_tool_mode == ToolMode::CutAndStretch) {
        if (m_cut_phase == CutPhase::Plane) {
            place_cut_plane(mouse_position);
            return;
        }
        if (m_cut_phase == CutPhase::Distance) {
            m_selection_error = _u8L("Return to the plane or guide-face step before changing the selection.");
            return;
        }
        additive = true; // guide faces click-toggle; unlike extrusion, their normals may differ
    }
    if (!additive)
        clear_selection();
    std::string error;
    auto face = pick_face(mouse_position, error); // ALWAYS a fresh click-time raycast.
    if (!face) {
        m_selection_error = std::move(error);
        return;
    }
    if (m_tool_mode == ToolMode::CutAndStretch && (m_smart_guide_fill || m_exact_guide_selection)) {
        select_cut_region(std::move(*face));
        return;
    }
    if (additive) {
        const auto selected = std::find_if(m_selected.begin(), m_selected.end(),
            [&face](const SelectedFace &entry) { return same_face(entry, *face); });
        if (selected != m_selected.end()) {
            clear_hover();
            m_selected.erase(selected);
            m_selection_error.clear();
            if (m_tool_mode == ToolMode::CutAndStretch && m_selected.empty()) clear_selection();
            else invalidate_preview();
            return;
        }
    }
    const auto normal = world_normal(face->normal_local, face->world_at_pick);
    if (!normal) {
        m_selection_error = _u8L("The selected face has an invalid normal or transformation.");
        return;
    }
    if (!m_selected.empty()) {
        const SelectedFace &first = m_selected.front();
        if (!same_source(*face, first)) {
            m_selection_error = _u8L("Select additional faces from the same source part.");
            return;
        }
        const auto reference = world_normal(first.normal_local, face->world_at_pick);
        if (m_tool_mode == ToolMode::Faces && (!reference || !concordant_normals(*normal, *reference))) {
            m_selection_error = _u8L("Additional faces must face the same world direction. The existing selection is unchanged.");
            return;
        }
    }
    if (m_tool_mode == ToolMode::CutAndStretch && m_selected.size() >= CutGuideSelectionLimits{}.max_guides) {
        m_selection_error = _u8L("Cut guide selection is limited to 4096 guides. The existing selection is unchanged.");
        return;
    }
    try {
        if (!load_face_triangles(*face, error)) {
            m_selection_error = std::move(error);
            return;
        }
        m_selected.push_back(std::move(*face));
        clear_hover();
        m_selection_error.clear();
        invalidate_preview();
    } catch (const std::exception &) {
        m_selection_error = _u8L("Unable to select this face.");
    }
    // Building is deferred to render/UI; a click never mutates the model.
}

void GLGizmoFaceExtrude::select_profile_vertex(const Vec2d &mouse_position)
{
    std::string error;
    auto point = pick_face(mouse_position, error);
    if (!point || point->vertex_id < 0) {
        m_selection_error = std::move(error);
        return;
    }
    if (!m_selected.empty()) {
        const auto &first = m_selected.front();
        if (point->object_id != first.object_id || point->instance_id != first.instance_id ||
            point->volume_id != first.volume_id || point->mesh != first.mesh) {
            m_selection_error = _u8L("The first point locks the source part and instance. Clear the profile before choosing another part.");
            return;
        }
        for (const auto &contour : m_profile_contours)
            for (int vertex : contour)
                if (vertex == point->vertex_id ||
                    (first.mesh->its.vertices[vertex].array() == first.mesh->its.vertices[point->vertex_id].array()).all()) {
                    m_selection_error = _u8L("This vertex is already in the profile. Contours close automatically; do not repeat the first point.");
                    return;
                }
    }
    if (m_selected.empty())
        m_selected.push_back(*point); // stable target token; never repick Boolean source from later points
    if (m_profile_contours.empty()) m_profile_contours.emplace_back();
    m_profile_contours.back().push_back(point->vertex_id);
    clear_hover();
    m_selection_error.clear();
    invalidate_preview();
    request_preview_redraw();
}

void GLGizmoFaceExtrude::undo_profile_vertex()
{
    if (m_profile_contours.empty()) return;
    if (!m_profile_contours.back().empty())
        m_profile_contours.back().pop_back();
    else if (m_profile_contours.size() > 1)
        m_profile_contours.pop_back();
    if (m_profile_contours.size() == 1 && m_profile_contours.front().empty()) {
        clear_selection();
    } else {
        m_selection_error.clear();
        invalidate_preview();
    }
    request_preview_redraw();
}

void GLGizmoFaceExtrude::update_hover()
{
    if (boolean_busy() || m_arrow_hovered || m_arrow_drag || m_left_gesture != LeftGesture::None ||
        (m_tool_mode == ToolMode::CutAndStretch &&
            (m_cut_phase == CutPhase::Distance || (m_cut_phase == CutPhase::Plane &&
             (m_plane_basis != PlaneBasis::Face || m_face_plane_ready))))) {
        clear_hover();
        return;
    }
    ModelObject *object;
    ModelInstance *instance;
    ModelVolume *negative = nullptr;
    if (!get_target(object, instance, nullptr, &negative) || !m_c || !m_c->selection_info() ||
        (!negative && !m_c->raycaster()) ||
        m_c->selection_info()->model_object() != object || m_c->selection_info()->get_active_instance() != 0) {
        clear_hover();
        return;
    }
    if (ImGui::GetCurrentContext()) {
        const ImGuiIO &io = ImGui::GetIO();
        // Like painting, never highlight through an ImGui window or while the
        // pointer is navigating/dragging. This also handles popup capture.
        if (io.WantCaptureMouse || io.MouseDown[0] || io.MouseDown[1] || io.MouseDown[2]) {
            clear_hover();
            return;
        }
    }
    const Camera &camera = m_parent.get_active_camera();
    const Vec2d mouse = m_parent.get_local_mouse_position();
    const auto &viewport = camera.get_viewport();
    if (!mouse.allFinite() || viewport[2] <= 0 || viewport[3] <= 0 ||
        mouse.x() < 0.0 || mouse.y() < 0.0 || mouse.x() >= viewport[2] || mouse.y() >= viewport[3] ||
        !camera.get_view_matrix().matrix().allFinite() || !camera.get_projection_matrix().matrix().allFinite()) {
        clear_hover();
        return;
    }

    HoverQuery query;
    query.object_id = object->id();
    query.instance_id = instance->id();
    query.mouse = mouse;
    query.view = camera.get_view_matrix();
    query.projection = camera.get_projection_matrix();
    query.viewport = viewport;
    for (const ModelVolume *volume : object->volumes) {
        if (negative ? volume != negative : !volume->is_model_part())
            continue;
        const Transform3d world = instance->get_matrix() * volume->get_matrix();
        if (!usable_transform(world)) {
            clear_hover();
            return;
        }
        query.volumes.push_back({volume->id(), volume->get_mesh_shared_ptr(), world, volume->type()});
    }
    // Painting's temporary seed selection is the model for this interaction,
    // but its mouse-only cache is insufficient for a stationary cursor while
    // orbiting, zooming, switching ortho/perspective or resizing the viewport.
    if (m_hover_query && m_hover_query->object_id == query.object_id &&
        m_hover_query->instance_id == query.instance_id && m_hover_query->mouse == query.mouse &&
        m_hover_query->viewport == query.viewport && m_hover_query->volumes == query.volumes &&
        (m_hover_query->view.matrix().array() == query.view.matrix().array()).all() &&
        (m_hover_query->projection.matrix().array() == query.projection.matrix().array()).all())
        return;
    m_hover_query = std::move(query);

    std::string hover_error; // Never overwrite an error concerning the locked face.
    auto face = pick_face(mouse, hover_error);
    if (!face) {
        clear_hover(false); // Cache misses too; no continuous raycasting while idle.
        if (m_tool_mode == ToolMode::CutAndStretch && m_cut_phase == CutPhase::Guides)
            m_guide_hover_error = std::move(hover_error);
        return;
    }
    if (m_tool_mode == ToolMode::CutAndStretch && m_cut_phase == CutPhase::Guides &&
        (m_smart_guide_fill || m_exact_guide_selection)) {
        update_cut_region_hover(std::move(*face));
        return;
    }
    if (m_hovered && same_face(*m_hovered, *face) &&
        (m_hovered->world_at_pick.matrix().array() == face->world_at_pick.matrix().array()).all())
        return; // Same plane: neither copy triangles nor regenerate a GL model.

    clear_hover(false);
    m_hovered = std::move(face);
    if (m_tool_mode == ToolMode::PointProfile) {
        m_parent.request_extra_frame(); // include the exact vertex marker in the frustum
        return; // rendered below, without planar-face measurement
    }
    // A locked face already has its amber highlight. Avoid z-fighting and a
    // second overlay there; other faces can be hovered without replacing it.
    if (std::any_of(m_selected.begin(), m_selected.end(),
                    [this](const SelectedFace &entry) { return same_face(entry, *m_hovered); }))
        return;
    try {
        if (!load_face_triangles(*m_hovered, hover_error))
            return;
        m_hover_model = make_hover_model(m_hovered->mesh->its, m_hovered->triangles,
                                        m_hovered->world_at_pick);
        if (m_hover_model)
            m_parent.request_extra_frame(); // Include the visual offset in the next frustum.
    } catch (const std::exception &) {
        m_hover_model.reset();
        // Retain the plane identity so a failed GL build is not retried on every
        // movement inside that same plane. Reload/transform/plane changes retry.
    }
}

bool GLGizmoFaceExtrude::resolve_source(ModelObject *&object, ModelInstance *&instance,
                                       ModelVolume *&volume, int *object_index) const
{
    volume = nullptr;
    ModelVolume *negative = nullptr;
    if (m_selected.empty() || !get_target(object, instance, object_index, &negative))
        return false;
    const SelectedFace &first = m_selected.front();
    if (object->id() != first.object_id || instance->id() != first.instance_id ||
        (negative ? (first.type_at_pick != ModelVolumeType::NEGATIVE_VOLUME || negative->id() != first.volume_id) :
                    first.type_at_pick != ModelVolumeType::MODEL_PART))
        return false;
    for (const SelectedFace &face : m_selected) {
        if (!same_source(face, first))
            return false;
    }
    for (ModelVolume *candidate : object->volumes) {
        if (candidate->id() == first.volume_id) {
            if (!face_extrusion_source_type_matches(m_tool_mode, first.type_at_pick, *candidate) ||
                candidate->get_mesh_shared_ptr() != first.mesh)
                return false;
            volume = candidate;
            return true;
        }
    }
    return false;
}

void GLGizmoFaceExtrude::update_preview()
{
    if (m_tool_mode != ToolMode::Faces) {
        update_variant_preview();
        return;
    }
    if (m_selected.empty() || m_committing)
        return;
    ModelObject *object;
    ModelInstance *instance;
    ModelVolume *volume;
    if (!resolve_source(object, instance, volume)) {
        clear_selection();
        m_error = _u8L("The object changed. Select the face again.");
        return;
    }
    const Transform3d world = instance->get_matrix() * volume->get_matrix();
    if (!(world.matrix().array() == m_world_matrix.matrix().array()).all() ||
        !(instance->get_matrix().matrix().array() == m_instance_matrix.matrix().array()).all()) {
        invalidate_boolean();
        m_world_dirty = true;
    }
    if (!m_world_dirty && !m_preview_dirty)
        return;

    m_preview_valid = false;
    m_previews.clear();
    m_error.clear();
    try {
        if (m_world_dirty) {
            m_world_dirty = false;
            m_world_matrix = world;
            m_instance_matrix = instance->get_matrix();
            m_face_model.reset();
            m_world_source = TriangleMesh();
            m_normal_world = Vec3d::Zero();
            m_preview_dirty = true;
            if (!usable_transform(world) || !usable_transform(m_instance_matrix)) {
                m_error = _u8L("The object has an invalid or singular transformation.");
                m_preview_dirty = false;
                return;
            }
            const auto reference_normal = world_normal(m_selected.front().normal_local, world);
            if (!reference_normal) {
                m_error = _u8L("The selected face has an invalid normal.");
                m_preview_dirty = false;
                return;
            }
            m_normal_world = *reference_normal;
            std::vector<int> highlight_triangles;
            for (const SelectedFace &face : m_selected) {
                const auto normal = world_normal(face.normal_local, world);
                if (!normal || !concordant_normals(*normal, m_normal_world)) {
                    m_error = _u8L("The selected faces no longer face the same world direction.");
                    m_preview_dirty = false;
                    return;
                }
                highlight_triangles.insert(highlight_triangles.end(), face.triangles.begin(), face.triangles.end());
            }
            std::sort(highlight_triangles.begin(), highlight_triangles.end());
            highlight_triangles.erase(std::unique(highlight_triangles.begin(), highlight_triangles.end()), highlight_triangles.end());
            m_world_source = *m_selected.front().mesh;
            // Handedness correction changes winding, not triangle IDs.
            m_world_source.transform(world, true);
            if (!usable_mesh(m_world_source.its)) {
                m_error = _u8L("The transformed mesh loses precision or contains invalid triangles.");
                m_preview_dirty = false;
                return;
            }
            // A combined highlight is visual only; these disconnected patches
            // are NEVER combined into one kernel input or one output mesh.
            m_face_model = init_plane_data(m_world_source.its, highlight_triangles, 0.05f);
            // The upstream highlight builder computes normals in FLOAT. Even
            // finite source coordinates may overflow there; never feed invalid
            // generated geometry into rendering or the camera's frustum box.
            if (!m_face_model || !finite_bounds(m_face_model->get_bounding_box()) ||
                !m_face_model->mesh || !usable_mesh(m_face_model->mesh->its)) {
                m_face_model.reset();
                m_preview_dirty = false;
                m_error = _u8L("The selected face exceeds the precision available for preview.");
                return;
            }
        }
        m_preview_dirty = false;
        // An unsuccessful transform rebuild must never leave usable stale data
        // behind when the distance is edited on a subsequent frame.
        if (!m_face_model || m_world_source.empty()) {
            m_error = _u8L("Select a valid planar face with a non-singular transformation.");
            return;
        }
        if (!std::isfinite(m_distance) || std::abs(m_distance) <= 1e-6) {
            m_error = _u8L("Enter a finite, non-zero distance in millimeters.");
            return;
        }
        // Validate EACH connected patch separately, including patches at
        // different elevations. Never relax the kernel's connectedness or
        // planarity constraints, or concatenate outputs into intersecting shells.
        std::vector<PreviewPart> previews;
        previews.reserve(m_selected.size());
        double combined_max_z = object->get_max_z();
        for (size_t face_index = 0; face_index < m_selected.size(); ++face_index) {
            FaceExtrusionResult result = extrude_planar_face(m_world_source.its,
                m_selected[face_index].triangles, m_distance * m_normal_world);
            if (!result || result.mesh.vertices.empty() || result.mesh.indices.empty()) {
                m_error = _u8L("Cannot extrude selected face") + " " + std::to_string(face_index + 1) + ":\n" + result.message;
                return;
            }
            PreviewPart part;
            part.world_mesh = TriangleMesh(std::move(result.mesh));
            const BoundingBoxf3 bounds = part.render_bounds = part.world_mesh.bounding_box();
            if (!finite_bounds(bounds)) {
                m_error = _u8L("The extrusion has invalid bounds.");
                return;
            }
            if (m_output_mode == OutputMode::NewParts && !face_extrusion_new_part_above_bed(volume->type(), bounds.min.z())) {
                m_error = _u8L("A new part would extend below the build plate. Reverse or reduce the distance.");
                return;
            }
            if (volume->is_model_part())
                combined_max_z = std::max(combined_max_z, bounds.max.z());
            part.model = std::make_unique<GLModel>();
            part.model->init_from(part.world_mesh.its);
            previews.push_back(std::move(part));
        }
        // changed_object calls ensure_on_bed(true) in FFF. Adding any number of
        // parts must not lift the original object due to this implicit refresh.
        if (m_output_mode == OutputMode::NewParts &&
            !face_extrusion_new_parts_bed_safe(volume->type(), *object, *instance, combined_max_z)) {
            m_error = _u8L("Creating these parts would move the original object onto the build plate.");
            return;
        }
        m_previews = std::move(previews);
        m_preview_valid = true;
        // Projection may already have been chosen before this render-time
        // rebuild. The next frame must account for the new world-space box.
        m_parent.request_extra_frame();
    } catch (const std::exception &) {
        m_preview_dirty = false;
        m_error = _u8L("Unable to build the extrusion preview.");
    }
}

void GLGizmoFaceExtrude::set_cut_phase(CutPhase phase)
{
    if (m_tool_mode != ToolMode::CutAndStretch || phase == CutPhase::Distance || m_selected.empty()) return;
    finish_arrow_drag(true);
    m_cut_phase = phase;
    m_moving_side_chosen = false;
    m_selection_error.clear();
    clear_hover();
    invalidate_preview(); // cancel OWN request, restore original; keep all guide faces and plane settings
    update_preview();
    request_preview_redraw();
}

void GLGizmoFaceExtrude::choose_cut_plane_basis(PlaneBasis basis)
{
    if (!editing_cut_plane() || m_arrow_drag) return;
    const Vec3d previous_origin = m_plane_origin_world;
    m_plane_basis = basis;
    m_face_plane_ready = basis != PlaneBasis::Face; // a later click, never a guide face, chooses this normal
    if (basis != PlaneBasis::Face) {
        m_plane_reference_local.reset();
        m_plane_source_facet.reset();
        m_plane_normal_world = basis == PlaneBasis::XY ? Vec3d::UnitZ() :
            basis == PlaneBasis::XZ ? Vec3d::UnitY() : Vec3d::UnitX();
        m_plane_position = m_plane_normal_world.dot(previous_origin);
    }
    m_moving_side_chosen = false;
    m_selection_error.clear();
    clear_hover();
    invalidate_preview();
    update_preview();
    request_preview_redraw();
}

bool GLGizmoFaceExtrude::use_cut_face_plane(const SelectedFace &face)
{
    if (!editing_cut_plane() || m_arrow_drag || boolean_busy() || m_selected.empty()) return false;
    ModelObject *object; ModelInstance *instance; ModelVolume *source;
    const SelectedFace &first = m_selected.front();
    if (!resolve_source(object, instance, source) || face.object_id != first.object_id ||
        face.instance_id != first.instance_id || face.volume_id != first.volume_id || face.mesh != first.mesh) {
        m_selection_error = _u8L("Choose a face of the locked source part.");
        return false;
    }
    const Transform3d world = instance->get_matrix() * source->get_matrix();
    const auto normal = world_normal(face.normal_local, world);
    const Vec3d reference = world * face.hit_local;
    if (!normal || !reference.allFinite() || !positioned_cut_plane(reference, *normal, normal->dot(reference))) {
        m_selection_error = _u8L("This face cannot define a plane within safe physical and world-coordinate precision.");
        return false;
    }
    m_plane_basis = PlaneBasis::Face;
    m_plane_normal_world = *normal;
    m_plane_reference_local = face.hit_local;
    m_plane_source_facet = face.seed_facet;
    m_plane_reference_world = reference;
    m_plane_position = normal->dot(reference);
    m_plane_initialized = true;
    m_face_plane_ready = true;
    m_moving_side_chosen = false;
    m_selection_error.clear();
    clear_hover();
    // In particular, a coplanar guide must use THIS exact guide point, not a
    // bbox-center projection or a second raycast sample with a different residual.
    invalidate_preview();
    update_preview();
    request_preview_redraw();
    return true;
}

void GLGizmoFaceExtrude::use_cut_guide_plane(size_t index)
{
    if (index < m_selected.size() && use_cut_face_plane(m_selected[index]))
        m_cut_guide_index = index;
}

void GLGizmoFaceExtrude::place_cut_plane(const Vec2d &mouse_position)
{
    if (!editing_cut_plane() || m_arrow_drag || m_selected.empty()) return;
    const bool pick_orientation = m_plane_basis == PlaneBasis::Face && !m_face_plane_ready;
    std::string error;
    auto hit = pick_face(mouse_position, error, !pick_orientation);
    if (!hit) { m_selection_error = error.empty() ? _u8L("Click a surface of the locked source part.") : error; return; }
    const auto &first = m_selected.front();
    if (hit->object_id != first.object_id || hit->instance_id != first.instance_id ||
        hit->volume_id != first.volume_id || hit->mesh != first.mesh) {
        m_selection_error = _u8L("Place the plane on the same locked source part; guide faces are unchanged.");
        return;
    }
    if (m_plane_basis == PlaneBasis::Face) {
        for (size_t i = 0; i < m_selected.size(); ++i) {
            const SelectedFace &guide = m_selected[i];
            // Reuse only a guide whose EXACT selected mask contains the clicked
            // source triangle. A shared measurement-plane ID is insufficient for partial regions.
            if (std::any_of(guide.triangles.begin(), guide.triangles.end(),
                    [&hit](int triangle) { return triangle >= 0 && size_t(triangle) == hit->seed_facet; })) {
                use_cut_guide_plane(i);
                return;
            }
        }
    }
    if (pick_orientation) {
        use_cut_face_plane(*hit); // other planar faces remain selectable too
        return;
    }
    m_plane_position = m_plane_normal_world.dot(hit->world_at_pick * hit->hit_local);
    invalidate_boolean();
    m_preview_dirty = true;
    m_selection_error.clear();
    clear_hover();
    update_preview(); // only the small graphical plane/matrix; never a cut job
    request_preview_redraw();
}

void GLGizmoFaceExtrude::update_cut_plane()
{
    if (!m_plane_initialized || !m_face_plane_ready) return;
    const auto origin = positioned_cut_plane(m_plane_reference_world, m_plane_normal_world, m_plane_position);
    if (!origin) throw std::runtime_error(_u8L("The plane offset or world-coordinate arithmetic exceeds safe precision. Keep the plane near the source."));
    m_plane_origin_world = *origin;
    m_plane_render_matrix = Transform3d::Identity();
    m_plane_render_matrix.translation() = *origin;
    if (!m_plane_render_matrix.matrix().cast<float>().allFinite())
        throw std::runtime_error(_u8L("The plane position exceeds preview precision."));
    if (!m_face_model) {
        const Vec3d u = m_plane_normal_world.unitOrthogonal(), v = m_plane_normal_world.cross(u);
        GLModel::Geometry plane;
        plane.format = {GLModel::PrimitiveType::Triangles, GLModel::Geometry::EVertexLayout::P3N3};
        for (const Vec2d &corner : {Vec2d(-1,-1), Vec2d(1,-1), Vec2d(1,1), Vec2d(-1,1)}) {
            const Vec3f point = (m_plane_radius * (corner.x() * u + corner.y() * v)).cast<float>();
            if (!point.allFinite()) throw std::runtime_error(_u8L("The cut plane exceeds preview precision."));
            plane.add_vertex(point, Vec3f(m_plane_normal_world.cast<float>()));
        }
        plane.add_triangle(0,1,2); plane.add_triangle(0,2,3);
        m_face_model = std::make_shared<GLModel>();
        m_face_model->init_from(std::move(plane));
    }
    const BoundingBoxf3 bounds = translated_face_tools_preview_bounds(m_face_model->get_bounding_box(), m_plane_render_matrix.translation());
    if (!finite_bounds(bounds) || !bounds.min.cast<float>().allFinite() || !bounds.max.cast<float>().allFinite())
        throw std::runtime_error(_u8L("The cut plane has invalid bounds."));
    m_normal_world = m_cut_phase == CutPhase::Distance && m_moving_opposite ? -m_plane_normal_world : m_plane_normal_world;
}

void GLGizmoFaceExtrude::confirm_cut_plane()
{
    if (!editing_cut_plane() || m_arrow_drag || boolean_busy() || !m_face_plane_ready) return;
    update_preview();
    ModelObject *object; ModelInstance *instance; ModelVolume *source;
    if (!resolve_source(object, instance, source) || !m_face_model || !m_error.empty()) return;
    const Transform3d world = instance->get_matrix() * source->get_matrix();
    std::vector<CutPlaneFaceGuide> guides;
    for (const auto &face : m_selected) {
        const auto normal = world_normal(face.normal_local, world);
        if (!normal) { m_selection_error = _u8L("A guide face has an invalid normal."); return; }
        guides.push_back({world * face.hit_local, *normal});
    }
    const auto side = select_cut_moving_side(guides, m_plane_origin_world, m_plane_normal_world);
    if (!side) { m_selection_error = side.message; request_preview_redraw(); return; }
    invalidate_boolean();
    m_cut_phase = CutPhase::Distance;
    m_cut_scope_confirmation.clear();
    m_moving_opposite = side.side < 0;
    m_moving_side_chosen = true;
    m_normal_world = side.moving_normal;
    m_prepare_requested = true; // explicit confirmation only, never plane dragging or orientation edits
    m_preview_dirty = true;
    m_selection_error.clear();
    clear_hover();
    request_preview_redraw();
    start_boolean(true);
}

void GLGizmoFaceExtrude::update_variant_preview()
{
    if (m_selected.empty() || m_committing) return;
    ModelObject *object; ModelInstance *instance; ModelVolume *source;
    if (!resolve_source(object, instance, source)) {
        clear_selection();
        m_error = _u8L("The source changed. Choose the vertices or cut plane again.");
        return;
    }
    const Transform3d world = instance->get_matrix() * source->get_matrix();
    if (!(world.matrix().array() == m_world_matrix.matrix().array()).all() ||
        !(instance->get_matrix().matrix().array() == m_instance_matrix.matrix().array()).all()) {
        invalidate_boolean();
        if (m_tool_mode == ToolMode::CutAndStretch && m_cut_phase == CutPhase::Distance) {
            m_cut_phase = CutPhase::Plane;
            m_moving_side_chosen = false;
            m_prepare_requested = false;
            m_selection_error = _u8L("The source pose changed. Confirm the independent cut plane again.");
        }
        m_world_dirty = true;
    }
    if (!m_world_dirty && !m_preview_dirty) {
        update_stretch_visibility();
        return;
    }
    m_preview_valid = false;
    m_error.clear();
    try {
        if (m_world_dirty) {
            restore_source_visibility();
            m_world_dirty = false;
            m_preview_dirty = true;
            m_previews.clear();
            m_stretch_plan.reset();
            m_cut_scope_confirmation.clear();
            m_profile_valid = false;
            m_profile_world.clear();
            m_profile_lines.reset();
            m_face_model.reset();
            m_normal_world = Vec3d::Zero();
            m_world_source = TriangleMesh();
            m_world_matrix = world;
            m_instance_matrix = instance->get_matrix();
            if (!usable_transform(world) || !usable_transform(m_instance_matrix))
                throw std::runtime_error(_u8L("The object has an invalid or singular transformation."));
            if (m_tool_mode == ToolMode::PointProfile) {
                m_world_source = *m_selected.front().mesh;
                m_world_source.transform(world, true);
                if (!usable_mesh(m_world_source.its))
                    throw std::runtime_error(_u8L("The source loses precision in world coordinates."));
                for (const auto &contour : m_profile_contours) {
                    std::vector<Vec3d> points;
                    for (int vertex : contour) {
                        if (vertex < 0 || size_t(vertex) >= source->mesh().its.vertices.size())
                            throw std::runtime_error(_u8L("The profile vertex identity changed."));
                        points.push_back(world * source->mesh().its.vertices[vertex].cast<double>());
                    }
                    m_profile_world.push_back(std::move(points));
                }
                m_profile_lines = profile_edge_model(m_profile_world);
                const std::vector<std::vector<Vec3d>> holes(m_profile_world.begin() + 1, m_profile_world.end());
                PlanarProfileResult profile = triangulate_planar_profile(m_profile_world.front(), holes);
                if (!profile)
                    throw std::runtime_error(profile.message);
                m_normal_world = profile.normal;
                m_profile_valid = true;
                m_face_model = std::make_shared<GLModel>();
                m_face_model->init_from(profile.cap);
            } else {
                m_guide_models.clear();
                for (const auto &face : m_selected) {
                    auto model = make_hover_model(face.mesh->its, face.triangles, world);
                    if (!model) throw std::runtime_error(_u8L("Unable to highlight a selected guide face."));
                    m_guide_models.push_back(std::move(model));
                }
                const BoundingBoxf3 bounds = source->mesh().bounding_box().transformed(world);
                if (!finite_bounds(bounds)) throw std::runtime_error(_u8L("The source has invalid bounds."));
                m_plane_reference_world = m_plane_reference_local ? Vec3d(world * *m_plane_reference_local) : bounds.center();
                m_plane_radius = std::max(1.0, 0.65 * bounds.size().norm() +
                    (m_plane_reference_world - bounds.center()).norm());
                if (m_cut_phase != CutPhase::Guides && !m_plane_initialized) {
                    m_plane_basis = PlaneBasis::XY;
                    m_plane_normal_world = Vec3d::UnitZ();
                    m_plane_position = m_plane_reference_world.z(); // absolute world Z, independent of every guide normal
                    m_face_plane_ready = true;
                    m_plane_initialized = true;
                }
            }
            if (m_tool_mode == ToolMode::PointProfile && (!m_face_model || !finite_bounds(m_face_model->get_bounding_box())))
                throw std::runtime_error(_u8L("The selected profile has invalid bounds."));
        }
        m_preview_dirty = false;
        if (m_tool_mode == ToolMode::PointProfile && m_world_source.empty()) return;
        if (m_tool_mode == ToolMode::CutAndStretch) {
            if (m_cut_phase == CutPhase::Guides) return; // selection first: no inferred plane/arrow/side
            if (!m_face_plane_ready) {
                m_error = _u8L("Click a planar face now to choose the independent plane orientation.");
                return;
            }
            update_cut_plane(); // cached small quad; position edits update only its render matrix
            if (m_cut_phase == CutPhase::Plane) return; // no preparation from positioning or dragging
            if (!m_moving_side_chosen) {
                m_error = _u8L("Confirm a plane that identifies one consistent side from the guide faces.");
                return;
            }
            if (!std::isfinite(m_distance) || m_distance < 0.0) {
                restore_source_visibility();
                m_error = _u8L("Cut and extend supports positive extension only. Negative shortening is not implemented.");
                return;
            }
            if (!m_stretch_plan) {
                m_error = _u8L("Prepare the cut preview before applying. The original is unchanged.");
                return;
            }
            if (m_previews.size() != 3 || m_rendered_stretch_plan != m_stretch_plan) {
                m_previews.clear();
                m_rendered_stretch_plan = m_stretch_plan;
                m_preview_bridge = m_stretch_plan->bridge_unit;
                m_bridge_distance.reset();
                for (const indexed_triangle_set *mesh : {&m_stretch_plan->stationary, &m_stretch_plan->moving, &m_stretch_plan->bridge_unit}) {
                    PreviewPart part;
                    for (const Vec3f &vertex : mesh->vertices) part.render_bounds.merge(vertex.cast<double>());
                    part.model = std::make_unique<GLModel>();
                    if (mesh != &m_stretch_plan->bridge_unit) part.model->init_from(*mesh);
                    m_previews.push_back(std::move(part));
                }
            }
            if (m_distance > 0.0 && (!m_bridge_distance || *m_bridge_distance != m_distance)) {
                update_bridge_preview(m_preview_bridge, m_stretch_plan->bridge_unit, m_bridge_base_z, m_distance);
                auto &bridge = m_previews[2];
                bridge.render_bounds.reset();
                for (const Vec3f &vertex : m_preview_bridge.vertices) bridge.render_bounds.merge(vertex.cast<double>());
                bridge.model = bridge_preview_model(m_preview_bridge); // bridge only; no CSG, weld or retriangulation
                m_bridge_distance = m_distance;
            }
            if (m_previews.size() != 3)
                throw std::runtime_error(_u8L("The prepared preview is incomplete."));
            m_previews[0].render_world = m_stretch_plan->cut_to_world;
            m_previews[1].render_world = m_stretch_plan->cut_to_world * Eigen::Translation3d(0.0, 0.0, m_distance);
            m_previews[2].render_world = m_stretch_plan->cut_to_world; // exact bottom + translated top, never Z-scale residuals
            for (size_t i = 0; i < m_previews.size(); ++i) {
                if (i == 2 && m_distance == 0.0) continue; // no singular bridge render at zero
                const auto &part = m_previews[i];
                const BoundingBoxf3 bounds = part.render_bounds.transformed(part.render_world);
                if (!usable_transform(part.render_world) || !part.render_world.matrix().cast<float>().allFinite() ||
                    !finite_bounds(bounds) || !bounds.min.cast<float>().allFinite() || !bounds.max.cast<float>().allFinite())
                    throw std::runtime_error(_u8L("The prepared preview exceeds available precision."));
            }
            m_preview_valid = m_distance > 0.0;
            if (!m_preview_valid) m_error = _u8L("Zero extension is a preview-only no-op; no Undo step is created.");
            update_stretch_visibility();
        } else {
            m_previews.clear();
            if (!m_profile_valid) {
                m_error = _u8L("Complete a valid ordered coplanar profile with at least three vertices per contour.");
                return;
            }
            if (!std::isfinite(m_distance) || std::abs(m_distance) <= 1e-6) {
                m_error = _u8L("Enter a finite, non-zero distance in millimeters.");
                return;
            }
            const std::vector<std::vector<Vec3d>> holes(m_profile_world.begin() + 1, m_profile_world.end());
            FaceExtrusionResult result = extrude_planar_profile(m_profile_world.front(), m_distance * m_normal_world, holes);
            if (!result) throw std::runtime_error(result.message);
            PreviewPart part;
            part.world_mesh = TriangleMesh(std::move(result.mesh));
            const BoundingBoxf3 bounds = part.render_bounds = part.world_mesh.bounding_box();
            if (!usable_mesh(part.world_mesh.its) || !finite_bounds(bounds))
                throw std::runtime_error(_u8L("The profile extrusion has invalid geometry or bounds."));
            if (m_output_mode == OutputMode::NewParts && (bounds.min.z() < -bed_tolerance_mm ||
                std::max(object->get_max_z(), bounds.max.z()) < SINKING_MIN_Z_THRESHOLD))
                throw std::runtime_error(_u8L("Creating this profile would move the original onto the build plate or extend below it."));
            part.model = std::make_unique<GLModel>();
            part.model->init_from(part.world_mesh.its);
            m_previews.push_back(std::move(part));
            m_preview_valid = true;
        }
        request_preview_redraw();
    } catch (const std::exception &exception) {
        restore_source_visibility();
        m_preview_dirty = false;
        m_error = exception.what();
    }
}

bool GLGizmoFaceExtrude::arrow_frame(const Transform3d &world, Vec3d &anchor, Vec3d &normal) const
{
    if (m_selected.empty()) return false;
    anchor = world * m_selected.front().hit_local;
    if (m_tool_mode == ToolMode::Faces) {
        const auto n = world_normal(m_selected.front().normal_local, world);
        if (!n) return false;
        normal = *n;
    } else {
        if (m_world_dirty || !(world.matrix().array() == m_world_matrix.matrix().array()).all()) return false;
        if (m_tool_mode == ToolMode::PointProfile && !m_profile_valid) return false;
        if (m_tool_mode == ToolMode::CutAndStretch && (m_cut_phase == CutPhase::Guides ||
            !m_plane_initialized || !m_face_plane_ready ||
            !positioned_cut_plane(m_plane_reference_world, m_plane_normal_world, m_plane_position) ||
            (m_cut_phase == CutPhase::Distance && !m_moving_side_chosen))) return false;
        normal = editing_cut_plane() ? m_plane_normal_world : m_normal_world;
        if (m_tool_mode == ToolMode::CutAndStretch) anchor = m_plane_origin_world;
    }
    return anchor.allFinite() && normal.allFinite() && std::abs(normal.squaredNorm() - 1.0) <= 1e-6;
}

void GLGizmoFaceExtrude::render_profile_markers()
{
    if (m_tool_mode != ToolMode::PointProfile) return;
    OverlayState state;
    glDisable(GL_DEPTH_TEST); // chosen vertices/contours stay legible while orbiting
    const Camera &camera = m_parent.get_active_camera();
    const double radius = 0.28 * get_grabber_size();
    if (!std::isfinite(radius) || radius <= 0.0) return;
    auto marker = [&](const Vec3d &point, const std::array<float,4> &color) {
        if (!point.allFinite()) return;
        Transform3d matrix = Transform3d::Identity();
        matrix.translate(point); matrix.scale(radius);
        if (matrix.matrix().cast<float>().allFinite())
            render_glmodel(m_sphere, color, camera.get_view_matrix() * matrix, camera.get_projection_matrix(), false, 0.5f);
    };
    if (m_profile_lines)
        render_glmodel(*m_profile_lines, {{1.0f,0.7f,0.1f,1.0f}}, camera.get_view_matrix(), camera.get_projection_matrix(), false, 0.7f);
    for (size_t c = 0; c < m_profile_world.size(); ++c)
        for (const Vec3d &point : m_profile_world[c])
            marker(point, c == 0 ? std::array<float,4>{{1.0f,0.65f,0.1f,1.0f}} : std::array<float,4>{{0.95f,0.25f,0.65f,1.0f}});
    if (m_hovered && m_hovered->vertex_id >= 0)
        marker(m_hovered->world_at_pick * m_hovered->hit_local, {{0.15f,1.0f,0.5f,1.0f}});
}

void GLGizmoFaceExtrude::request_preview_redraw()
{
    // request_extra_frame alone cannot wake an idle, clean canvas: its idle
    // handler tests m_dirty before consulting the extra-frame flag. This custom
    // gesture bypasses the standard gizmo update path that normally sets dirty.
    set_dirty();
    m_parent.set_as_dirty();
    m_parent.request_extra_frame();
}

std::string GLGizmoFaceExtrude::arrow_parameter_label() const
{
    const bool plane = m_arrow_drag ? m_arrow_drag->plane_drag : editing_cut_plane();
    if (!plane) return arrow_distance_label(m_distance);
    const std::string coordinate = m_plane_basis == PlaneBasis::XY ? "Z" :
        m_plane_basis == PlaneBasis::XZ ? "Y" : m_plane_basis == PlaneBasis::YZ ? "X" : "n.XYZ";
    return _u8L("Plane") + " " + coordinate + ": " + arrow_distance_label(m_plane_position);
}

void GLGizmoFaceExtrude::render_drag_tooltip()
{
    if (!m_arrow_drag || !ImGui::GetCurrentContext() || !ImGui::IsMousePosValid())
        return;
    const ImGuiIO &io = ImGui::GetIO();
    if (!std::isfinite(io.MousePos.x) || !std::isfinite(io.MousePos.y) ||
        io.DisplaySize.x <= 0.f || io.DisplaySize.y <= 0.f)
        return;
    const std::string text = arrow_parameter_label() + "\n" + _u8L("1 mm steps; Esc cancels");
    const ImVec2 size = ImGui::CalcTextSize(text.c_str());
    const float pad = std::max(4.f, ImGui::GetFontSize() * 0.4f);
    const ImVec2 box(size.x + 2.f * pad, size.y + 2.f * pad);
    const ImVec2 position(
        std::clamp(io.MousePos.x + 2.f * pad, 0.f, std::max(0.f, io.DisplaySize.x - box.x)),
        std::clamp(io.MousePos.y + 2.f * pad, 0.f, std::max(0.f, io.DisplaySize.y - box.y)));
    // No window, input item, focus or capture: this remains visible even when
    // the pointer crosses the ImGui panel while the arrow owns the drag.
    ImDrawList *draw = ImGui::GetForegroundDrawList();
    draw->AddRectFilled(position, ImVec2(position.x + box.x, position.y + box.y), IM_COL32(24, 28, 32, 240), pad);
    draw->AddText(ImVec2(position.x + pad, position.y + pad), IM_COL32(255, 255, 255, 255), text.c_str());
}

void GLGizmoFaceExtrude::update_arrow()
{
    const BoundingBoxf3 previous_bounds = m_arrow_bounds;
    m_arrow_available = false;
    m_arrow_shaft_available = false;
    m_arrow_bounds.reset();
    if (!m_grabbers.empty())
        m_grabbers.front().enabled = false;
    ModelObject *object;
    ModelInstance *instance;
    ModelVolume *source;
    if (m_state != On || m_committing || boolean_busy() || !resolve_source(object, instance, source) || m_grabbers.empty())
        return;
    const Transform3d world = instance->get_matrix() * source->get_matrix();
    Vec3d anchor, normal_value;
    if (!arrow_frame(world, anchor, normal_value)) return;
    const std::optional<Vec3d> normal = normal_value;
    const double radius = 0.75 * get_grabber_size();
    if (!normal || !anchor.allFinite() || !std::isfinite(radius) || radius <= 0.0)
        return;

    // Zero distance has no solid preview, but MUST retain a usable handle.
    // Invalid/extreme numeric input must not inject huge bounds into the camera:
    // in that case retain the last valid handle position (or the initial face).
    if (editing_cut_plane())
        m_handle_distance = 0.0; // plane widget sits ON the plane; default extension10 never offsets this arrow
    else if (std::isfinite(m_distance) && ((m_preview_valid && !m_preview_dirty && !m_world_dirty) || std::abs(m_distance) <= 1e-6))
        m_handle_distance = m_distance;
    const double shaft_length = m_handle_distance + 2.0 * radius;
    const Matrix3d orientation = Eigen::Quaterniond::FromTwoVectors(Vec3d::UnitZ(), *normal).toRotationMatrix();
    m_arrow_matrix = Transform3d::Identity();
    m_arrow_matrix.linear() = orientation * Vec3d(radius, radius, 2.0 * radius).asDiagonal();
    m_arrow_matrix.translation() = anchor + shaft_length * *normal;
    if (!usable_transform(m_arrow_matrix) || !m_arrow_matrix.matrix().cast<float>().allFinite() ||
        !m_arrow_matrix.linear().inverse().transpose().cast<float>().allFinite())
        return;
    m_arrow_bounds = m_cone.get_bounding_box().transformed(m_arrow_matrix);
    if (!finite_bounds(m_arrow_bounds) || !m_arrow_bounds.min.cast<float>().allFinite() || !m_arrow_bounds.max.cast<float>().allFinite())
        return;

    if (std::abs(shaft_length) > 1e-9) {
        const Vec3d shaft_direction = shaft_length < 0.0 ? -*normal : *normal;
        m_arrow_shaft_matrix = Transform3d::Identity();
        m_arrow_shaft_matrix.linear() = Eigen::Quaterniond::FromTwoVectors(Vec3d::UnitZ(), shaft_direction).toRotationMatrix() *
            Vec3d(0.12 * radius, 0.12 * radius, std::abs(shaft_length)).asDiagonal();
        m_arrow_shaft_matrix.translation() = anchor;
        if (usable_transform(m_arrow_shaft_matrix) && m_arrow_shaft_matrix.matrix().cast<float>().allFinite() &&
            m_arrow_shaft_matrix.linear().inverse().transpose().cast<float>().allFinite()) {
            const BoundingBoxf3 shaft_bounds = m_cylinder.get_bounding_box().transformed(m_arrow_shaft_matrix);
            if (finite_bounds(shaft_bounds) && shaft_bounds.min.cast<float>().allFinite() && shaft_bounds.max.cast<float>().allFinite()) {
                m_arrow_shaft_available = true;
                m_arrow_bounds.merge(shaft_bounds);
            }
        }
    }
    if (!m_arrow_raycaster) {
        m_arrow_pick_mesh = TriangleMesh(its_make_cone(1.0, 1.0, 2.0 * PI / 24.0));
        m_arrow_raycaster = std::make_unique<MeshRaycaster>(m_arrow_pick_mesh);
    }
    m_arrow_available = true;
    m_grabbers.front().enabled = true;
    m_grabbers.front().center = m_arrow_matrix.translation();
    m_grabbers.front().set_model_matrix(m_arrow_matrix);
    if (!finite_bounds(previous_bounds) || previous_bounds.min != m_arrow_bounds.min || previous_bounds.max != m_arrow_bounds.max)
        m_parent.request_extra_frame();
}

bool GLGizmoFaceExtrude::arrow_contains_mouse(const Vec2d &mouse) const
{
    if (!m_arrow_available || !m_arrow_raycaster || !mouse.allFinite())
        return false;
    const Camera &camera = m_parent.get_active_camera();
    const auto &viewport = camera.get_viewport();
    if (mouse.x() < 0.0 || mouse.y() < 0.0 || mouse.x() >= viewport[2] || mouse.y() >= viewport[3] ||
        !camera.get_view_matrix().matrix().allFinite() || !camera.get_projection_matrix().matrix().allFinite())
        return false;
    Vec3f hit, normal;
    // The cone is an overlay, including for inward extrusions. Pick the actual
    // cone mesh, not the source surface behind it or a stale framebuffer ID.
    return m_arrow_raycaster->closest_hit(mouse, m_arrow_matrix, camera, hit, normal) && hit.allFinite();
}

bool GLGizmoFaceExtrude::begin_arrow_drag(const Vec2d &mouse)
{
    ModelObject *object;
    ModelInstance *instance;
    ModelVolume *source;
    wxWindow *canvas = m_parent.get_wxglcanvas();
    const bool plane_drag = editing_cut_plane();
    const double value = plane_drag ? m_plane_position : m_distance;
    if (boolean_busy() || !m_arrow_available || !canvas || !std::isfinite(value) || !resolve_source(object, instance, source))
        return false;
    if (wxWindow::GetCapture() && wxWindow::GetCapture() != canvas)
        return false;
    ArrowDrag drag;
    drag.source_world = instance->get_matrix() * source->get_matrix();
    if (!arrow_frame(drag.source_world, drag.anchor, drag.normal))
        return false; // valid profile normal / chosen moving side, never an unrelated facet normal
    const auto parameter = ray_axis_parameter(m_parent.mouse_ray(Point(mouse.x(), mouse.y())), drag.anchor, drag.normal);
    if (!parameter) {
        m_selection_error = _u8L("Rotate the view to drag this direction, or enter the distance numerically.");
        return false;
    }
    const Camera &camera = m_parent.get_active_camera();
    drag.view = camera.get_view_matrix();
    drag.projection = camera.get_projection_matrix();
    drag.viewport = camera.get_viewport();
    drag.initial_parameter = *parameter;
    drag.initial_distance = value;
    drag.raw_distance = value;
    drag.restore_distance = value;
    drag.plane_drag = plane_drag; // target stays fixed throughout camera rebases and cancellation
    drag.last_mouse = mouse;
    m_arrow_drag = drag;
    m_selection_error.clear();
    clear_hover();
    m_left_gesture = LeftGesture::Arrow;
    m_hover_id = 0;
    start_dragging();
    // NoSnapPreview: do not call Selection::start_dragging, use_grabbers or
    // do_stop_dragging. Their model refresh would invalidate our locked faces.
    m_parent.set_mouse_as_dragging();
    m_capture_window = canvas;
    canvas->Bind(wxEVT_MOUSE_CAPTURE_LOST, &GLGizmoFaceExtrude::on_arrow_capture_lost, this);
    canvas->Bind(wxEVT_KILL_FOCUS, &GLGizmoFaceExtrude::on_arrow_focus_lost, this);
    if (!canvas->HasCapture())
        canvas->CaptureMouse();
    m_parent.post_event(SimpleEvent(EVT_GLCANVAS_MOUSE_DRAGGING_STARTED));
    request_preview_redraw();
    return true;
}

void GLGizmoFaceExtrude::update_arrow_drag(const Vec2d &mouse)
{
    if (!m_arrow_drag || !mouse.allFinite())
        return;
    ModelObject *object;
    ModelInstance *instance;
    ModelVolume *source;
    if (!resolve_source(object, instance, source) ||
        !((instance->get_matrix() * source->get_matrix()).matrix().array() == m_arrow_drag->source_world.matrix().array()).all()) {
        clear_selection();
        return;
    }
    const Camera &camera = m_parent.get_active_camera();
    const bool same_camera = same_screen_camera(camera.get_view_matrix(), camera.get_projection_matrix(),
                                                m_arrow_drag->view, m_arrow_drag->projection) &&
                             camera.get_viewport() == m_arrow_drag->viewport;
    if (same_camera && mouse == m_arrow_drag->last_mouse)
        return; // duplicate event/render sample; near/far-only changes cannot add drift
    const auto parameter = ray_axis_parameter(m_parent.mouse_ray(Point(mouse.x(), mouse.y())),
                                              m_arrow_drag->anchor, m_arrow_drag->normal);
    if (!parameter)
        return;
    m_arrow_drag->last_mouse = mouse;
    request_preview_redraw(); // also move the live badge between snapped steps
    if (!same_camera) {
        // Preserve the RAW fractional residue, not the displayed snapped value.
        // Repeated zoom/rebase must never accumulate or discard sub-step motion.
        m_arrow_drag->view = camera.get_view_matrix();
        m_arrow_drag->projection = camera.get_projection_matrix();
        m_arrow_drag->viewport = camera.get_viewport();
        m_arrow_drag->initial_parameter = *parameter;
        m_arrow_drag->initial_distance = m_arrow_drag->raw_distance;
        return;
    }
    const double raw_distance = m_arrow_drag->initial_distance + (*parameter - m_arrow_drag->initial_parameter);
    const auto distance = snapped_arrow_distance(m_arrow_drag->restore_distance, raw_distance);
    if (!distance)
        return;
    m_arrow_drag->raw_distance = raw_distance;
    double &value = m_arrow_drag->plane_drag ? m_plane_position : m_distance;
    if (*distance == value)
        return;
    invalidate_boolean();
    value = *distance;
    m_preview_dirty = true;
    m_parent.set_mouse_as_dragging();
}

void GLGizmoFaceExtrude::finish_arrow_drag(bool cancel)
{
    if (!m_arrow_drag)
        return;
    if (cancel) {
        double &value = m_arrow_drag->plane_drag ? m_plane_position : m_distance;
        value = m_arrow_drag->restore_distance;
        m_preview_dirty = true;
    }
    m_arrow_drag.reset();
    stop_dragging();
    m_left_gesture = LeftGesture::Cancelled;
    m_consume_left_up = true;
    const wxWeakRef<wxWindow> capture_window = m_capture_window;
    if (m_capture_window) {
        m_capture_window->Unbind(wxEVT_MOUSE_CAPTURE_LOST, &GLGizmoFaceExtrude::on_arrow_capture_lost, this);
        m_capture_window->Unbind(wxEVT_KILL_FOCUS, &GLGizmoFaceExtrude::on_arrow_focus_lost, this);
        if (m_capture_window->HasCapture())
            m_capture_window->ReleaseMouse();
        m_capture_window = nullptr;
    }
    if (capture_window && m_parent.get_wxglcanvas() == capture_window.get()) {
        m_parent.mouse_up_cleanup();
        m_parent.post_event(SimpleEvent(EVT_GLCANVAS_MOUSE_DRAGGING_FINISHED));
    }
    request_preview_redraw();
}

void GLGizmoFaceExtrude::on_arrow_capture_lost(wxMouseCaptureLostEvent &)
{
    finish_arrow_drag(true);
}

void GLGizmoFaceExtrude::on_arrow_focus_lost(wxFocusEvent &event)
{
    finish_arrow_drag(true);
    event.Skip();
}

void GLGizmoFaceExtrude::render_arrow(bool picking)
{
    if (!m_arrow_available)
        return;
    OverlayState state;
    glDisable(GL_DEPTH_TEST); // standard gizmo overlay; also usable inside the source
    glDisable(GL_BLEND);
    const Camera &camera = picking ? wxGetApp().plater()->get_picking_camera() : m_parent.get_active_camera();
    const std::array<float, 4> color = picking ? picking_color_component(0) :
        (editing_cut_plane() ? std::array<float,4>{{1.0f,0.95f,0.05f,1.0f}} :
         m_arrow_hovered || m_arrow_drag ? GRABBER_HOVER_COL : GRABBER_NORMAL_COL);
    if (!picking && m_arrow_shaft_available)
        render_glmodel(m_cylinder, color, camera.get_view_matrix() * m_arrow_shaft_matrix,
                       camera.get_projection_matrix(), false, 0.2f);
    render_glmodel(m_cone, color, camera.get_view_matrix() * m_arrow_matrix,
                   camera.get_projection_matrix(), picking, 0.2f);
}

void GLGizmoFaceExtrude::on_render_for_picking()
{
    if (m_state != On || m_committing)
        return;
    update_arrow();
    render_arrow(true);
}

bool GLGizmoFaceExtrude::click_threshold_exceeded(const Vec2d &mouse) const
{
    const double canvas_scale = m_parent.get_scale();
    const double scale = std::isfinite(canvas_scale) && canvas_scale > 0.0 ? canvas_scale : 1.0;
    return !mouse.allFinite() || (mouse - m_press_mouse).squaredNorm() > std::pow(5.0 * scale, 2);
}

bool GLGizmoFaceExtrude::click_camera_unchanged() const
{
    const Camera &camera = m_parent.get_active_camera();
    // Ignore near/far-only changes caused by clearing the hover offset. Zoom,
    // orbit and viewport changes still cancel a pending face click.
    return same_screen_camera(camera.get_view_matrix(), camera.get_projection_matrix(), m_press_view, m_press_projection) &&
        camera.get_viewport() == m_press_viewport;
}

bool GLGizmoFaceExtrude::on_mouse(const wxMouseEvent &event)
{
    if (m_state != On)
        return false;
    const Vec2d mouse(event.GetX(), event.GetY()); // already Retina-scaled by Canvas
    const bool in_ui = ImGui::GetCurrentContext() && ImGui::GetIO().WantCaptureMouse;
    if (event.Moving()) {
        finish_arrow_drag(false); // recover a missing release; never create a part
        m_left_gesture = LeftGesture::None;
        m_consume_left_up = false;
        m_parent.request_extra_frame();
        return false;
    }
    if (event.Leaving()) {
        finish_arrow_drag(false);
        if (m_left_gesture == LeftGesture::PendingFace)
            m_left_gesture = LeftGesture::Cancelled;
        clear_hover();
        m_arrow_hovered = false;
        return false;
    }
    if (event.RightDown() || event.MiddleDown()) {
        if (m_arrow_drag || m_left_gesture == LeftGesture::PendingFace) {
            finish_arrow_drag(false);
            m_left_gesture = LeftGesture::Cancelled;
        }
        return false;
    }
    if (event.RightUp()) {
        // Native RightDown/Dragging already pan. The canvas RightUp handler also
        // selects the volume under the pointer: suppress ONLY that selection.
        // Canvas performs mouse_up_cleanup even for a consumed button release.
        return true;
    }
    if (event.LeftDown() || event.LeftDClick()) {
        // Idempotent if manager and canvas dispatch the same press twice.
        if (m_left_gesture == LeftGesture::PendingFace || m_left_gesture == LeftGesture::Orbit || m_left_gesture == LeftGesture::Arrow)
            return true;
        m_consume_left_up = true;
        m_press_mouse = mouse;
        m_press_shift = event.ShiftDown();
        const Camera &camera = m_parent.get_active_camera();
        m_press_view = camera.get_view_matrix();
        m_press_projection = camera.get_projection_matrix();
        m_press_viewport = camera.get_viewport();
        clear_hover();
        if (event.CmdDown() || boolean_busy()) {
            m_left_gesture = LeftGesture::Orbit;
        } else if (m_committing || in_ui) {
            m_left_gesture = LeftGesture::Cancelled;
        } else {
            update_preview();
            update_arrow();
            if (arrow_contains_mouse(mouse)) {
                if (!begin_arrow_drag(mouse))
                    m_left_gesture = LeftGesture::Cancelled;
            } else {
                m_left_gesture = LeftGesture::PendingFace;
            }
        }
        m_parent.request_extra_frame();
        // Always consume Down, including Cmd: otherwise the canvas can arm
        // object translation, toggle object selection or start a Shift rectangle.
        return true;
    }
    if (event.Dragging()) {
        if (!event.LeftIsDown()) {
            finish_arrow_drag(false);
            if (m_left_gesture == LeftGesture::PendingFace)
                m_left_gesture = LeftGesture::Cancelled;
            return false; // right/middle navigation
        }
        if (m_left_gesture == LeftGesture::Arrow) {
            if (event.CmdDown()) {
                finish_arrow_drag(false);
                m_left_gesture = LeftGesture::Orbit;
            } else {
                update_arrow_drag(mouse);
                return true;
            }
        }
        if (m_left_gesture == LeftGesture::PendingFace) {
            if (event.CmdDown() || click_threshold_exceeded(mouse) || !click_camera_unchanged())
                m_left_gesture = LeftGesture::Orbit;
            else
                return true;
        }
        if (m_left_gesture == LeftGesture::Orbit) {
            // Down never armed move_volume_idx or a selection rectangle. Canvas
            // uses the first passed drag to seed its orbit origin, then follows
            // its normal camera code. Shift never toggles a face during orbit.
            m_consume_left_up = true;
            return false;
        }
        return m_left_gesture == LeftGesture::Cancelled;
    }
    if (event.LeftUp()) {
        const LeftGesture gesture = m_left_gesture;
        if (gesture == LeftGesture::Arrow) {
            update_arrow_drag(mouse);
            finish_arrow_drag(false);
        } else if (gesture == LeftGesture::PendingFace && !m_committing && !in_ui && !event.CmdDown() &&
                   !click_threshold_exceeded(mouse) && click_camera_unchanged()) {
            // Select only after a real click, with a fresh CURRENT-camera hit.
            select_face(mouse, m_press_shift);
        }
        m_left_gesture = LeftGesture::None;
        if (gesture != LeftGesture::None)
            m_consume_left_up = true;
        // Canvas's consumed-Up path runs mouse_up_cleanup (including rotation
        // center reset), without entering its plate/real-object selection path.
        m_parent.request_extra_frame();
        // Retain this result until Moving/next press for duplicate Up dispatch.
        return m_consume_left_up;
    }
    return false; // middle release / wheel / unrelated navigation events
}

bool GLGizmoFaceExtrude::gizmo_event(SLAGizmoEventType action, const Vec2d &,
                                    bool, bool, bool)
{
    if (m_state != On)
        return false;
    if (boolean_busy() && (action == SLAGizmoEventType::Escape || action == SLAGizmoEventType::Delete)) {
        m_boolean_request->cancel.store(true, std::memory_order_relaxed);
        m_parent.request_extra_frame();
        return true;
    }
    if (action == SLAGizmoEventType::Escape && m_arrow_drag) {
        finish_arrow_drag(true); // keep faces, restore pre-drag distance; no Undo
        return true;
    }
    if (action == SLAGizmoEventType::Delete || action == SLAGizmoEventType::Escape) {
        const bool had_selection = !m_selected.empty();
        clear_selection();
        m_left_gesture = LeftGesture::Cancelled;
        m_consume_left_up = true;
        m_parent.request_extra_frame();
        // Delete NEVER deletes the original object. The next Escape, when no
        // face remains, may be handled by the manager to exit this gizmo.
        return action == SLAGizmoEventType::Delete || had_selection;
    }
    // Deliberately no Enter/ApplyChanges path: an active distance input must not
    // accidentally commit. Creation is only through the explicit UI button.
    return false;
}

void GLGizmoFaceExtrude::on_render()
{
    if (m_state != On || m_committing)
        return;
    update_negative_visibility();
    if (m_arrow_drag && (!m_capture_window || !m_capture_window->HasCapture()))
        finish_arrow_drag(true);
    else if (m_arrow_drag && ImGui::GetCurrentContext() && !ImGui::GetIO().MouseDown[0])
        finish_arrow_drag(false); // release consumed by ImGui/toolbar before gizmo dispatch
    else if (m_arrow_drag && ImGui::GetCurrentContext() && ImGui::IsMousePosValid()) {
        // Canvas may render immediately on ImGui capture BEFORE dispatching the
        // drag to us. Sample its current pointer now, then rebuild ALL prisms and
        // the arrow in this same frame. Duplicate samples are idempotent.
        const ImVec2 mouse = ImGui::GetIO().MousePos;
        update_arrow_drag(Vec2d(mouse.x, mouse.y));
    }
    update_preview();
    update_arrow();
    m_arrow_hovered = !m_arrow_drag && m_left_gesture == LeftGesture::None &&
        !(ImGui::GetCurrentContext() && ImGui::GetIO().WantCaptureMouse) &&
        arrow_contains_mouse(m_parent.get_local_mouse_position());
    m_hover_id = m_arrow_hovered || m_arrow_drag ? 0 : -1;
    update_hover();
    if (!m_face_model && m_guide_models.empty() && !m_preview_valid && !m_hover_model && !m_arrow_available &&
        !(m_tool_mode == ToolMode::PointProfile && (!m_selected.empty() || m_hovered)))
        return;
    const Camera &camera = m_parent.get_active_camera();
    OverlayState state;
    // Stretch operands remain in one cached cut frame. Bridge caps use exact
    // base+distance heights; the whole moving portion gets rigid +Z translation.
    const bool stretch_zero = m_tool_mode == ToolMode::CutAndStretch && m_stretch_plan &&
        m_previews.size() == 3 && m_distance == 0.0;
    if (m_preview_valid || stretch_zero) {
        for (size_t i = 0; i < m_previews.size(); ++i) {
            if (stretch_zero && i == 2) continue;
            const auto &part = m_previews[i];
            const std::array<float,4> color = m_tool_mode == ToolMode::CutAndStretch ?
                (i == 1 ? std::array<float,4>{{1.0f,0.55f,0.08f,0.9f}} :
                 i == 2 ? std::array<float,4>{{0.1f,0.8f,0.95f,0.75f}} : std::array<float,4>{{0.45f,0.55f,0.65f,0.85f}}) :
                m_output_mode == OutputMode::Subtract ? std::array<float,4>{{0.95f,0.22f,0.18f,0.42f}} :
                                                       std::array<float,4>{{0.15f,0.72f,0.95f,0.42f}};
            render_glmodel(*part.model, color, camera.get_view_matrix() * part.render_world, camera.get_projection_matrix());
        }
    }
    if (!m_hidden_source)
        for (const auto &guide : m_guide_models)
            render_glmodel(*guide, {{1.0f, 0.45f, 0.06f, 0.8f}}, camera.get_view_matrix(), camera.get_projection_matrix());
    if (m_hover_model)
        render_glmodel(*m_hover_model, {{0.15f, 0.9f, 0.5f, 0.55f}},
                       camera.get_view_matrix(), camera.get_projection_matrix());
    if (m_face_model) {
        const bool culling = glIsEnabled(GL_CULL_FACE);
        if (m_tool_mode == ToolMode::CutAndStretch) glDisable(GL_CULL_FACE); // explicit plane visible from either side
        render_glmodel(*m_face_model, m_tool_mode == ToolMode::CutAndStretch ?
                       std::array<float,4>{{1.0f, 0.95f, 0.05f, 0.35f}} : std::array<float,4>{{1.0f, 0.65f, 0.12f, 0.65f}},
                       camera.get_view_matrix() * (m_tool_mode == ToolMode::CutAndStretch ? m_plane_render_matrix : Transform3d::Identity()),
                       camera.get_projection_matrix());
        if (culling) glEnable(GL_CULL_FACE);
    }
    if (m_preview_valid && (m_distance < 0.0 || m_output_mode == OutputMode::Subtract)) {
        WireframeState wireframe;
        for (const PreviewPart &part : m_previews)
            render_glmodel(*part.model, m_output_mode == OutputMode::Subtract ?
                           std::array<float, 4>{{1.0f, 0.25f, 0.15f, 0.8f}} : std::array<float, 4>{{0.95f, 0.3f, 0.65f, 0.75f}},
                           camera.get_view_matrix() * part.render_world, camera.get_projection_matrix(), false, 0.4f);
    }
    render_profile_markers();
    render_arrow(false);
}

void GLGizmoFaceExtrude::on_render_input_window(float x, float y, float bottom_limit)
{
    if (m_state != On || m_committing) return;
    update_preview();
    update_negative_visibility();
    const bool profile = m_tool_mode == ToolMode::PointProfile;
    const bool stretch = m_tool_mode == ToolMode::CutAndStretch;
    ModelObject *target_object;
    ModelInstance *target_instance;
    ModelVolume *negative_source = nullptr;
    const bool negative = m_tool_mode == ToolMode::Faces &&
        get_target(target_object, target_instance, nullptr, &negative_source) && negative_source;
    bool to_plane = false, confirm_plane = false, back_to_plane = false, back_to_guides = false;
    y = std::max(0.0f, std::min(y, bottom_limit - m_window_height));
    GizmoImguiSetNextWIndowPos(x, y, ImGuiCond_Always);
    ImGuiWrapper::push_toolbar_style(m_parent.get_scale());
    GizmoImguiBegin(on_get_name(), ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse);
    const float width = m_imgui->scaled(27.0f);
    if (profile) {
        m_imgui->text_wrapped(_u8L("Click visible mesh vertices in boundary order: at least three coplanar, non-collinear points. The first point locks the source part and instance."), width);
        m_imgui->text_wrapped(_u8L("Contours close automatically. No convex hull, flattening or inherited holes: encode each hole explicitly. The arrow follows the ordered profile normal."), width);
    } else if (stretch) {
        m_imgui->text(m_cut_phase == CutPhase::Guides ? _u8L("1 / 3 - Select moving guide faces") :
            m_cut_phase == CutPhase::Plane ? _u8L("2 / 3 - Position the independent plane") : _u8L("3 / 3 - Inspect the target and extend"));
        m_imgui->text_wrapped(m_cut_phase == CutPhase::Guides ?
            _u8L("Choose guide faces on one moving piece. Guides identify a component separated by the plane, not an arbitrary branch within a connected solid.") :
            m_cut_phase == CutPhase::Plane ? _u8L("Yellow plane: click the source, enter its absolute coordinate or drag the yellow arrow. The model stays unchanged.") :
            _u8L("Orange is the WHOLE piece that will move, not just the guide faces. Gray stays fixed; blue extends its contact. Inspect and confirm this actual region before applying."), width);
    } else {
        if (negative)
            m_imgui->text_wrapped(_u8L("The selected negative part is isolated while this tool is open. Pick its own planar faces, including below the bed. Other parts are restored on exit."), width);
        m_imgui->text_wrapped(_u8L("Hover a planar face to highlight it; click to select."), width);
        m_imgui->text_wrapped(_u8L("Shift + click: add or remove faces facing the same direction, on the same source part."), width);
    }
    if (!stretch)
        m_imgui->text_wrapped(_u8L("Arrow: live 1 mm steps; type any custom distance. Drag elsewhere to orbit. Esc restores the distance before a drag."), width);
    ImGui::Separator();
    m_imgui->disabled_begin(boolean_busy() || m_arrow_drag.has_value());
    if (profile) {
        if (ImGui::Button(_u8L("Undo last point").c_str())) undo_profile_vertex();
        ImGui::SameLine();
        if (ImGui::Button(_u8L("Clear profile").c_str())) { clear_selection(); request_preview_redraw(); }
        m_imgui->disabled_begin(!m_profile_valid);
        if (ImGui::Button(_u8L("Add hole contour").c_str())) {
            m_profile_contours.emplace_back();
            invalidate_preview();
            request_preview_redraw();
        }
        m_imgui->disabled_end();
        for (size_t i = 0; i < m_profile_contours.size(); ++i)
            m_imgui->text((i == 0 ? _u8L("Outer contour") : _u8L("Hole") + " " + std::to_string(i)) +
                ": " + std::to_string(m_profile_contours[i].size()) + " " + _u8L("vertices"));
    } else if (stretch) {
        if (m_cut_phase == CutPhase::Guides) {
            int scope = m_cut_guided ? 0 : 1;
            bool scope_changed = ImGui::RadioButton(_u8L("Guide-associated piece").c_str(), &scope, 0);
            scope_changed |= ImGui::RadioButton(_u8L("Whole moving half-space (explicit)").c_str(), &scope, 1);
            if (scope_changed) {
                m_cut_guided = scope == 0;
                invalidate_preview();
                request_preview_redraw();
            }
            m_imgui->text_wrapped(m_cut_guided ?
                _u8L("Every selected guide triangle must identify the same cut component. Ambiguous or fixed-side guides are refused; no automatic whole-model fallback.") :
                _u8L("This explicitly moves the entire supported half-space of the locked part. Guide faces choose the side, not a smaller region."), width);
            ImGui::Separator();
            if (ImGui::Checkbox(_u8L("Smart fill by angle").c_str(), &m_smart_guide_fill)) {
                clear_hover(); // changing the selector never changes the existing locked mask
                request_preview_redraw();
            }
            m_imgui->tooltip(_u8L("Uses the native neighbor-angle fill in local mesh coordinates, without painting. A smooth chain can select an entire connected surface. 0 degrees still uses native tolerance; 90 degrees can cross opposing normals. Check the highlight before clicking."), width);
            if (m_smart_guide_fill) {
                float angle = m_smart_guide_angle;
                ImGui::SetNextItemWidth(m_imgui->scaled(13.0f));
                if (ImGui::SliderFloat(_u8L("Guide fill angle").c_str(), &angle, 0.0f, 90.0f, "%.1f deg")) {
                    m_smart_guide_angle = std::isfinite(angle) ? std::clamp(angle, 0.0f, 90.0f) : 30.0f;
                    clear_hover();
                    request_preview_redraw();
                }
            }
            if (m_smart_guide_fill || m_exact_guide_selection)
                m_imgui->text_wrapped(_u8L("Green previews missing triangles; selected triangles stay orange. Click adds missing triangles to the guide selection, or removes the whole region if it is already selected. Turning this option off keeps the selected regions."), width);
            m_imgui->text(_u8L("Selected guide faces") + ": " + std::to_string(m_selected.size()));
            if (m_exact_guide_selection) {
                size_t triangles = 0;
                for (const auto &guide : m_selected) triangles += guide.triangles.size();
                m_imgui->text(_u8L("Selected source triangles") + ": " + std::to_string(triangles));
            }
        }
        else if (m_cut_phase == CutPhase::Plane) {
            if (!m_selected.empty()) {
                m_cut_guide_index = std::min(m_cut_guide_index, m_selected.size() - 1);
                if (m_selected.size() > 1) {
                    const std::string chosen = "G" + std::to_string(m_cut_guide_index + 1);
                    ImGui::SetNextItemWidth(m_imgui->scaled(8.0f));
                    if (ImGui::BeginCombo(_u8L("Guide face").c_str(), chosen.c_str())) {
                        for (size_t i = 0; i < m_selected.size(); ++i) {
                            const std::string label = "G" + std::to_string(i + 1);
                            if (ImGui::Selectable(label.c_str(), i == m_cut_guide_index)) m_cut_guide_index = i;
                        }
                        ImGui::EndCombo();
                    }
                }
                if (ImGui::Button(_u8L("Use guide face as plane").c_str())) use_cut_guide_plane(m_cut_guide_index);
            }
            int basis = int(m_plane_basis);
            bool orientation = ImGui::RadioButton("XY / Z", &basis, int(PlaneBasis::XY));
            ImGui::SameLine(); orientation |= ImGui::RadioButton("XZ / Y", &basis, int(PlaneBasis::XZ));
            ImGui::SameLine(); orientation |= ImGui::RadioButton("YZ / X", &basis, int(PlaneBasis::YZ));
            orientation |= ImGui::RadioButton(_u8L("From a face").c_str(), &basis, int(PlaneBasis::Face));
            if (orientation) choose_cut_plane_basis(PlaneBasis(basis));
            m_imgui->disabled_begin(!m_face_plane_ready);
            const std::string label = m_plane_basis == PlaneBasis::XY ? _u8L("Plane Z (world mm)") :
                m_plane_basis == PlaneBasis::XZ ? _u8L("Plane Y (world mm)") :
                m_plane_basis == PlaneBasis::YZ ? _u8L("Plane X (world mm)") : _u8L("Plane coordinate (world mm)");
            ImGui::SetNextItemWidth(m_imgui->scaled(10.0f));
            if (ImGui::InputDouble(label.c_str(), &m_plane_position, 0.1, 1.0, "%.3f")) {
                invalidate_boolean();
                m_preview_dirty = true;
                m_selection_error.clear();
                update_preview();
                request_preview_redraw();
            }
            m_imgui->disabled_end();
            if (m_plane_basis == PlaneBasis::Face)
                m_imgui->text_wrapped(m_face_plane_ready ? _u8L("Coordinate = normal dot world XYZ. Click a selected guide to reuse its exact plane, or another surface to position. From a face again chooses another orientation.") :
                    _u8L("Now click a planar face to choose its orientation. Guide faces remain selected."), width);
        }
    }
    if (!stretch) {
        int mode = int(m_output_mode);
        bool mode_changed = ImGui::RadioButton((negative ? _u8L("New negative part") : _u8L("New part")).c_str(), &mode, int(OutputMode::NewParts));
        if (!negative) ImGui::SameLine();
        mode_changed |= ImGui::RadioButton((negative ? _u8L("Union with negative part") : _u8L("Union")).c_str(), &mode, int(OutputMode::Union));
        if (!negative) ImGui::SameLine();
        mode_changed |= ImGui::RadioButton((negative ? _u8L("Subtract from negative part") : _u8L("Subtract")).c_str(), &mode, int(OutputMode::Subtract));
        if (mode_changed) {
            invalidate_boolean();
            m_output_mode = OutputMode(mode);
            if (m_output_mode == OutputMode::NewParts && m_annotations) m_annotations->consent.reset();
            m_preview_dirty = true;
            update_preview();
            request_preview_redraw();
        }
    }
    m_imgui->disabled_end();
    if (!stretch) {
        if (m_output_mode == OutputMode::NewParts)
            m_imgui->text_wrapped(profile ? _u8L("One closed profile extrusion becomes a new part. The source remains unchanged; no union is performed.") :
                negative ? _u8L("One separate negative cutter is created per face. The original negative part and all positive parts stay unchanged; no union is performed.") :
                _u8L("One separate closed part is created per face. The original part is not changed. No union is performed."), width);
        else {
            m_imgui->text_wrapped(_u8L("The colored extrusion is only an operand preview, NOT the Boolean result. Apply validates the result before replacing the locked source part."), width);
            if (negative)
                m_imgui->text_wrapped(m_output_mode == OutputMode::Union ?
                    _u8L("Union enlarges the negative cutter and its void. Only the selected negative part is replaced; positive parts are not Boolean operands.") :
                    _u8L("Subtraction removes overlap from the negative cutter, reducing its void. Choose a distance that enters the negative mesh. The result remains negative."), width);
            else if (m_output_mode == OutputMode::Subtract)
                m_imgui->text_wrapped(_u8L("Subtraction removes overlapping source material. Choose a distance that enters the source."), width);
        }
    }
    const Model *project = m_parent.get_model();
    const std::string project_limit = project ? project_exterior_limit(*project) : _u8L("The project is unavailable.");
    ModelObject *object = nullptr; ModelInstance *instance = nullptr; ModelVolume *source = nullptr;
    const bool has_source = resolve_source(object, instance, source);
    if (has_source) m_imgui->text_wrapped((negative ? _u8L("Locked negative source") : _u8L("Locked source part")) + ": " + source->name, width);
    const bool annotation_source = (!stretch || m_cut_phase == CutPhase::Distance) &&
        (stretch || m_output_mode != OutputMode::NewParts) && has_source;
    std::string source_limit;
    if (annotation_source) {
        source_limit = boolean_source_limit(*source);
        update_annotation_state(*source);
    }
    if (!stretch || m_cut_phase == CutPhase::Distance) {
        ImGui::Separator();
        m_imgui->disabled_begin(boolean_busy() || m_arrow_drag.has_value());
        ImGui::SetNextItemWidth(m_imgui->scaled(10.0f));
        bool changed = ImGui::InputDouble(_u8L("Distance (mm)").c_str(), &m_distance, 0.1, 1.0, "%.3f");
        if (!stretch && ImGui::Button(_u8L("Reverse direction").c_str())) { m_distance = -m_distance; changed = true; }
        m_imgui->disabled_end();
        if (changed) {
            invalidate_boolean();
            m_preview_dirty = true;
            update_preview();
            request_preview_redraw();
        }
        m_imgui->text_wrapped(stretch ? _u8L("Arrow: 1 mm steps along the plane normal. Zero is no-op; negative shortening is unsupported. Esc cancels the drag.") :
            profile ? _u8L("Positive: along the ordered profile normal. Negative: opposite direction. Only explicitly encoded holes are preserved.") :
            negative ? _u8L("Positive distance: outward from the negative mesh. Negative distance: inward. Normals and distance signs are not inverted for cutters.") :
            _u8L("Positive distance: outward. Negative distance: inward (parts may overlap)."), width);
    }
    if (!profile && !stretch) m_imgui->text(_u8L("Selected faces") + ": " + std::to_string(m_selected.size()));
    if (stretch && m_cut_phase == CutPhase::Distance && m_stretch_plan) {
        ImGui::Separator();
        if (m_stretch_plan->scope == CutAndStretchScope::GuidedComponent)
            m_imgui->text(_u8L("Moving cut component") + ": " +
                std::to_string(m_stretch_plan->selection.selected_upper_component + 1) + " / " +
                std::to_string(m_stretch_plan->selection.upper_component_count));
        else
            m_imgui->text_wrapped(_u8L("Explicit whole-half scope: every branch in this highlighted half moves."), width);
        m_imgui->text_wrapped(_u8L("A T-shaped bar may remain one connected piece: both arms then move. If orange includes an unwanted branch, go back and choose a separating plane; guides do not create a CAD boundary."), width);
        const bool region_ready = has_source && m_rendered_stretch_plan == m_stretch_plan &&
            m_previews.size() == 3 && std::isfinite(m_distance) && m_distance >= 0.0 &&
            (m_preview_valid || m_distance == 0.0);
        bool confirmed = m_cut_scope_confirmation.matches(m_stretch_plan, cut_scope(m_cut_guided));
        m_imgui->disabled_begin(!region_ready || boolean_busy() || m_arrow_drag.has_value());
        if (ImGui::Checkbox(_u8L("I confirm the entire orange region as the moving target").c_str(), &confirmed)) {
            invalidate_boolean();
            if (confirmed) m_cut_scope_confirmation.confirm(m_stretch_plan, cut_scope(m_cut_guided));
            else m_cut_scope_confirmation.clear();
            request_preview_redraw();
        }
        m_imgui->disabled_end();
    }
    bool annotation_blocked = false;
    if (annotation_source && m_annotations && m_annotations->info.any()) {
        ImGui::Separator();
        m_imgui->text_wrapped(_u8L("Stored annotations detected on the source part. These may come from the imported file, not manual painting:"), width);
        m_imgui->text_wrapped(annotation_reasons(m_annotations->info), width);
        if (!source_limit.empty() || !project_limit.empty()) {
            annotation_blocked = true;
            m_imgui->text_wrapped(_u8L("The restriction below cannot be bypassed by annotation consent."), width);
        } else {
            bool discard = m_annotations->consent.has_value();
            m_imgui->disabled_begin(boolean_busy() || m_arrow_drag.has_value());
            if (ImGui::Checkbox(_u8L("Discard source annotations on success").c_str(), &discard)) {
                invalidate_boolean();
                if (discard) {
                    try {
                        update_annotation_state(*source, true);
                        m_annotations->consent = capture_face_extrusion_annotation_stamp(*source);
                    } catch (const std::exception &) {
                        m_annotations->consent.reset();
                        m_boolean_error = _u8L("Unable to record annotation consent. Nothing was changed.");
                    }
                } else {
                    m_annotations->consent.reset();
                }
                // Only consent changes. Keep all picked faces, distance and preview.
                request_preview_redraw();
            }
            m_imgui->disabled_end();
            annotation_blocked = !m_annotations->consent.has_value();
            m_imgui->text_wrapped(_u8L("Nothing is removed now. Only a successful source replacement discards these supported source annotations, in the same Undo step. Other parts stay untouched. Procedural and geometry checks still apply."), width);
        }
    }
    if (!m_selection_error.empty()) m_imgui->error_text_wrapped(m_selection_error.c_str(), width);
    else if (stretch && m_cut_phase == CutPhase::Guides && !m_guide_hover_error.empty())
        m_imgui->error_text_wrapped(m_guide_hover_error.c_str(), width);
    if (!m_error.empty()) m_imgui->error_text_wrapped(m_error.c_str(), width);
    if (!project_limit.empty()) m_imgui->error_text_wrapped(project_limit.c_str(), width);
    else if (!source_limit.empty()) m_imgui->error_text_wrapped(source_limit.c_str(), width);
    if (!m_boolean_error.empty() && m_boolean_error != project_limit) m_imgui->error_text_wrapped(m_boolean_error.c_str(), width);
    if (boolean_busy())
        m_imgui->text_wrapped(m_boolean_request->cancel.load(std::memory_order_relaxed) ?
            _u8L("Cancellation requested. The original part is still unchanged.") :
            m_boolean_request->prepare_stretch ? _u8L("Preparing cut preview in the background. No Model or Undo change is made.") :
            _u8L("Calculating... The source is unchanged until a valid result is applied. Camera navigation remains available."), width);
    ImGui::Separator();
    bool prepare = false, create = false;
    if (stretch && m_cut_phase == CutPhase::Guides) {
        m_imgui->disabled_begin(m_selected.empty() || !m_error.empty() || boolean_busy() || m_arrow_drag.has_value());
        to_plane = ImGui::Button(_u8L("Next: choose plane").c_str());
        m_imgui->disabled_end();
    } else if (stretch && m_cut_phase == CutPhase::Plane) {
        m_imgui->disabled_begin(!m_face_plane_ready || !m_face_model || !m_error.empty() || boolean_busy() || m_arrow_drag.has_value());
        confirm_plane = ImGui::Button(_u8L("Confirm plane and preview").c_str());
        m_imgui->disabled_end();
    } else {
        if (stretch && !m_stretch_plan) {
            m_imgui->disabled_begin(boolean_busy() || m_arrow_drag.has_value() || !source_limit.empty() || !project_limit.empty());
            prepare = ImGui::Button(_u8L("Retry preview").c_str());
            m_imgui->disabled_end();
        }
        m_imgui->disabled_begin(!m_preview_valid || m_arrow_drag.has_value() || boolean_busy() || !source_limit.empty() || !project_limit.empty() || annotation_blocked ||
            (stretch && !m_cut_scope_confirmation.matches(m_stretch_plan, cut_scope(m_cut_guided))));
        const std::string create_label = stretch ? _u8L("Apply cut and extend") :
            m_output_mode == OutputMode::Union ? (negative ? _u8L("Apply union to negative part") : _u8L("Apply union")) :
            m_output_mode == OutputMode::Subtract ? (negative ? _u8L("Apply subtraction to negative part") : _u8L("Apply subtraction")) :
            negative ? (m_selected.size() == 1 ? _u8L("Create new negative part") : _u8L("Create negative parts") + " (" + std::to_string(m_selected.size()) + ")") :
            m_selected.size() == 1 ? _u8L("Create new part") : _u8L("Create parts") + " (" + std::to_string(m_selected.size()) + ")";
        create = ImGui::Button(create_label.c_str());
        m_imgui->disabled_end();
    }
    if (stretch && m_cut_phase != CutPhase::Guides) {
        m_imgui->disabled_begin(m_arrow_drag.has_value());
        if (m_cut_phase == CutPhase::Distance) {
            back_to_plane = ImGui::Button(_u8L("Back: plane").c_str());
            ImGui::SameLine();
        }
        back_to_guides = ImGui::Button(_u8L("Back: guide faces").c_str());
        m_imgui->disabled_end();
    }
    ImGui::SameLine();
    const bool cancel = ImGui::Button((boolean_busy() ? _u8L("Cancel calculation") : _u8L("Clear selection")).c_str());
    if (!stretch) m_imgui->text_wrapped(_u8L("Esc: cancel, clear selection, exit. Delete: clear selection."), width);
    const float height = ImGui::GetWindowHeight();
    if (height != m_window_height) { m_window_height = height; m_imgui->set_requires_extra_frame(); }
    GizmoImguiEnd();
    ImGuiWrapper::pop_toolbar_style();
    render_drag_tooltip();
    if (profile && ImGui::GetCurrentContext()) {
        ImDrawList *draw = ImGui::GetForegroundDrawList();
        size_t number = 0;
        for (const auto &contour : m_profile_world)
            for (const Vec3d &point : contour) {
                ++number;
                if (const auto screen = profile_screen_point(point, m_parent.get_active_camera())) {
                    const std::string label = std::to_string(number);
                    draw->AddText(ImVec2(float(screen->x() + 6.0), float(screen->y() + 6.0)), IM_COL32(255,200,45,255), label.c_str());
                }
            }
    }
    if (stretch && m_cut_phase != CutPhase::Distance && !m_world_dirty && ImGui::GetCurrentContext()) {
        ImDrawList *draw = ImGui::GetForegroundDrawList();
        for (size_t i = 0; i < m_selected.size(); ++i)
            if (const auto screen = profile_screen_point(m_world_matrix * m_selected[i].hit_local, m_parent.get_active_camera())) {
                const std::string label = "G" + std::to_string(i + 1);
                draw->AddText(ImVec2(float(screen->x() + 6.0), float(screen->y() + 6.0)), IM_COL32(255,150,40,255), label.c_str());
            }
    }
    // End ImGui before any action that can refresh the model/selection. No member
    // access follows Apply, because its callbacks may destroy this gizmo.
    if (to_plane || back_to_plane) {
        set_cut_phase(CutPhase::Plane);
    } else if (back_to_guides) {
        set_cut_phase(CutPhase::Guides);
    } else if (confirm_plane) {
        confirm_cut_plane();
    } else if (cancel) {
        if (boolean_busy()) m_boolean_request->cancel.store(true, std::memory_order_relaxed);
        else clear_selection();
        request_preview_redraw();
    } else if (create) {
        if (!stretch && m_output_mode == OutputMode::NewParts) create_part();
        else start_boolean();
    } else if (stretch && m_cut_phase == CutPhase::Distance && (prepare || m_prepare_requested) && !boolean_busy() && !m_arrow_drag &&
               m_left_gesture == LeftGesture::None && !ImGui::GetIO().MouseDown[0] &&
               !ImGui::GetIO().MouseDown[1] && !ImGui::GetIO().MouseDown[2] &&
               !ImGui::IsAnyItemActive() && source_limit.empty() && project_limit.empty()) {
        start_boolean(true);
    }
}

void GLGizmoFaceExtrude::start_boolean(bool prepare_stretch)
{
    if (m_state != On || m_committing || m_arrow_drag || boolean_busy() || m_output_mode == OutputMode::NewParts)
        return;
    update_preview();
    ModelObject *object;
    ModelInstance *instance;
    ModelVolume *source;
    Plater *plater = wxGetApp().plater();
    if (prepare_stretch) m_prepare_requested = false;
    const bool stretch = m_tool_mode == ToolMode::CutAndStretch;
    if (!plater || !resolve_source(object, instance, source) || (prepare_stretch && !stretch)) return;
    if (stretch) {
        if (m_cut_phase != CutPhase::Distance || !m_face_plane_ready || !m_moving_side_chosen || !m_face_model || !m_plane_origin_world.allFinite() ||
            !m_normal_world.allFinite() || std::abs(m_normal_world.squaredNorm() - 1.0) > 1e-6) return;
        if (!prepare_stretch && (!m_preview_valid || !m_stretch_plan || m_distance <= 0.0 ||
            !m_cut_scope_confirmation.matches(m_stretch_plan, cut_scope(m_cut_guided)))) return;
    } else if (!m_preview_valid || m_previews.size() != m_selected.size()) return;
    m_boolean_error = project_exterior_limit(plater->model());
    if (m_boolean_error.empty())
        m_boolean_error = boolean_source_limit(*source);
    if (!m_boolean_error.empty()) {
        request_preview_redraw();
        return;
    }
    Worker &worker = plater->get_ui_job_worker();
    if (!worker.is_idle()) {
        m_boolean_error = _u8L("Another operation is still running. Wait before calculating geometry.");
        if (prepare_stretch) m_prepare_requested = true;
        return;
    }
    if (!usable_transform(m_world_matrix) || !usable_transform(m_world_matrix.inverse())) {
        m_boolean_error = _u8L("The source transformation is unsupported for a precise Boolean result.");
        return;
    }
    invalidate_boolean();
    try {
        update_annotation_state(*source, true);
        auto annotation_stamp = capture_face_extrusion_annotation_stamp(*source);
        const bool discard_annotations = !prepare_stretch && m_annotations->info.any();
        if (discard_annotations && (!m_annotations->consent || *m_annotations->consent != annotation_stamp)) {
            m_annotations->consent.reset();
            m_boolean_error = _u8L("Review the detected source annotations and explicitly consent to discard them. If they changed, consent is required again. Selected faces are kept.");
            request_preview_redraw();
            return;
        }
        auto request = std::make_shared<BooleanRequest>();
        request->annotation_stamp = std::move(annotation_stamp);
        request->discard_annotations = discard_annotations;
        request->epoch = m_boolean_lifetime->epoch;
        request->mode = m_output_mode;
        request->tool_mode = m_tool_mode;
        request->prepare_stretch = prepare_stretch;
        request->cut_guided = m_cut_guided;
        request->profile_contours = m_profile_contours;
        request->plane_position = m_plane_position;
        request->cut_phase = m_cut_phase;
        request->plane_basis = m_plane_basis;
        request->plane_normal = m_plane_normal_world;
        request->plane_reference = m_plane_reference_world;
        request->plane_reference_local = m_plane_reference_local;
        request->plane_source_facet = m_plane_source_facet;
        request->moving_opposite = m_moving_opposite;
        request->moving_side_chosen = m_moving_side_chosen;
        request->plane_origin = m_plane_origin_world;
        request->moving_normal = m_normal_world;
        request->stretch_plan = prepare_stretch ? nullptr : m_stretch_plan;
        request->distance = m_distance;
        request->object_id = object->id();
        request->instance_id = instance->id();
        request->source_identity = source->get_mesh_shared_ptr();
        if (!stretch) request->source_world = m_world_source;
        request->status = prepare_stretch ? (m_cut_guided ? _u8L("Identifying the guide-associated moving piece...") :
            _u8L("Preparing the explicitly selected whole half-space...")) :
            stretch ? _u8L("Applying cut and extension...") :
            m_output_mode == OutputMode::Union ? _u8L("Computing union...") : _u8L("Computing subtraction...");
        for (const SelectedFace &face : m_selected) {
            request->faces.emplace_back(face.seed_facet, face.plane_id);
            request->guide_hits.push_back(face.hit_local);
            request->guide_normals.push_back(face.normal_local);
            request->guide_facets.push_back(face.triangles);
        }
        if (!stretch)
            for (const PreviewPart &part : m_previews)
                request->cutters.push_back(part.world_mesh);
        static_cast<FaceExtrusionPlacementContext &>(*request) =
            capture_face_extrusion_placement(*object, *instance, *source, m_world_matrix);
        m_boolean_request = request;
        clear_hover();
        m_arrow_hovered = false;
        m_left_gesture = LeftGesture::None;
        const std::weak_ptr<BooleanLifetime> lifetime = m_boolean_lifetime;
        const bool queued = queue_job(worker, [request](JobNew::Ctl &ctl) {
            // ONLY immutable copied geometry / value metadata and the request's
            // atomic cancel flag cross into this worker. No Model/GUI pointers.
            const auto canceled = [&] { return request->cancel.load(std::memory_order_relaxed) || ctl.was_canceled(); };
            const auto check_cancel = [&] { if (canceled()) throw BooleanCanceled{}; };
            // Every exit (success, rejected plan, cancellation or exception) must
            // close the progress notification. Early returns after request->error
            // used to leave it frozen at the last reported percentage (e.g. 36%)
            // with the Stop button still shown although the job had finished.
            struct CloseProgress {
                JobNew::Ctl &ctl; const decltype(request->status) &status;
                ~CloseProgress() { try { ctl.update_status(100, status); } catch (...) {} }
            } close_progress{ctl, request->status};
            try {
                check_cancel();
                TriangleMesh current;
                if (request->tool_mode == ToolMode::CutAndStretch) {
                    const auto progress = [&](int value) {
                        ctl.update_status(std::clamp(value, 0, 100) * 80 / 100, request->status);
                    };
                    if (request->prepare_stretch) {
                        // Immutable LOCAL mesh + double placement go directly to the
                        // near-source physical frame. Never narrow later-plate WORLD
                        // coordinates to float before cutting; no live Model mutation.
                        const auto facet = exact_cut_plane_facet(request->plane_basis == PlaneBasis::Face,
                            request->plane_source_facet, request->source_identity->its.indices.size(),
                            request->plane_normal, request->plane_reference, request->plane_position);
                        CutAndStretchPlan plan;
                        if (request->cut_guided) {
                            std::vector<int> guides;
                            for (const auto &mask : request->guide_facets) {
                                check_cancel();
                                if (mask.empty()) throw std::runtime_error("A selected guide has no original facet mask.");
                                guides.insert(guides.end(), mask.begin(), mask.end());
                            }
                            std::sort(guides.begin(), guides.end());
                            guides.erase(std::unique(guides.begin(), guides.end()), guides.end());
                            check_cancel();
                            plan = facet ?
                                prepare_cut_and_stretch_guided_from_face(request->source_identity->its, request->world,
                                    *facet, request->moving_normal, guides, canceled, progress) :
                                prepare_cut_and_stretch_guided(request->source_identity->its, request->world,
                                    request->plane_origin, request->moving_normal, guides, canceled, progress);
                        } else {
                            plan = facet ?
                                prepare_cut_and_stretch_from_face(request->source_identity->its, request->world,
                                    *facet, request->moving_normal, canceled, progress) :
                                prepare_cut_and_stretch(request->source_identity->its, request->world,
                                    request->plane_origin, request->moving_normal, canceled, progress);
                        }
                        check_cancel();
                        if (!plan) { request->error = plan.message; return; }
                        if (plan.scope != cut_scope(request->cut_guided))
                            throw std::runtime_error("The prepared region does not match the requested cut scope. Nothing was changed.");
                        if (!usable_transform(plan.cut_to_world) ||
                            !(plan.cut_to_world.linear().transpose() * plan.cut_to_world.linear()).isApprox(Matrix3d::Identity(), 1e-8) ||
                            std::abs(plan.cut_to_world.linear().determinant() - 1.0) > 1e-8 ||
                            !std::isfinite(plan.contact_area_mm2) || plan.contact_area_mm2 <= 0.0 ||
                            !std::isfinite(plan.source_volume_mm3) || plan.source_volume_mm3 <= 0.0)
                            throw std::runtime_error("The prepared cut frame or contact section is invalid.");
                        for (const indexed_triangle_set *mesh : {&plan.stationary, &plan.moving, &plan.bridge_unit})
                            validate_boolean_mesh(*mesh, canceled);
                        check_cancel();
                        request->bridge_base_z = bridge_base_heights(plan.bridge_unit);
                        request->stretch_plan = std::make_shared<const CutAndStretchPlan>(std::move(plan));
                        ctl.update_status(100, request->status);
                        return; // CPU cache only; no result, Snapshot or live Model change
                    }
                    if (!request->stretch_plan || request->distance <= 0.0)
                        throw std::runtime_error("A prepared plane and positive extension are required.");
                    CutAndStretchResult result = apply_cut_and_stretch_to_source(*request->stretch_plan,
                        request->distance, canceled, progress);
                    check_cancel();
                    if (!result) { request->error = result.message; return; }
                    current = TriangleMesh(std::move(result.mesh)); // SOURCE LOCAL replacement, never union with the old arm
                } else {
                    current = std::move(request->source_world);
                }
                const std::string operation = request->mode == OutputMode::Union ? "UNION" : "A_NOT_B";
                for (size_t i = 0; i < request->cutters.size(); ++i) {
                    check_cancel();
                    bool failed = false;
                    std::vector<TriangleMesh> outputs;
                    MeshBoolean::mcut::make_boolean(current, request->cutters[i], outputs, operation, canceled,
                        [&](float progress) {
                            if (!std::isfinite(progress)) return;
                            const double fraction = (double(i) + std::clamp(double(progress), 0.0, 100.0) / 100.0) / request->cutters.size();
                            ctl.update_status(int(80.0 * fraction), request->status);
                        }, [&] { failed = true; });
                    check_cancel();
                    if (failed) {
                        request->error = "The Boolean backend could not validate this operation.";
                        return;
                    }
                    indexed_triangle_set combined;
                    for (const TriangleMesh &output : outputs) {
                        check_cancel();
                        if (output.empty()) continue;
                        if (!usable_mesh(output.its) ||
                            output.its.vertices.size() > size_t(std::numeric_limits<int>::max()) - combined.vertices.size() ||
                            output.its.indices.size() > size_t(std::numeric_limits<int>::max()) - combined.indices.size())
                            throw std::runtime_error("Boolean output exceeds safe indexed mesh limits.");
                        // Preserve EVERY output and its winding, including cavity
                        // shells. Never TriangleMesh::split or per-shell flip.
                        its_merge(combined, output.its);
                    }
                    if (combined.indices.empty()) {
                        request->empty = true;
                        return; // valid empty difference, NOT a backend failure
                    }
                    current = TriangleMesh(std::move(combined));
                }
                ctl.update_status(82, request->status);
                validate_boolean_mesh(current.its, canceled);
                TriangleMesh local = current;
                if (request->tool_mode != ToolMode::CutAndStretch) {
                    local.transform(request->world.inverse(), true);
                    validate_boolean_mesh(local.its, canceled);
                    for (size_t i = 0; i < local.its.vertices.size(); ++i) {
                        if ((i & 1023) == 0) check_cancel();
                        const Vec3d world = request->world * local.its.vertices[i].cast<double>();
                        if (!world.allFinite() || (world - current.its.vertices[i].cast<double>()).norm() > bed_tolerance_mm)
                            throw std::runtime_error("Boolean result loses precision in the source's local coordinates.");
                    }
                    // Faces/profiles retain their existing WORLD-float round trip.
                    TriangleMesh round_trip = local;
                    round_trip.transform(request->world, true);
                    validate_boolean_mesh(round_trip.its, canceled);
                }
                // Cut already returns source-local storage. Its kernel validates
                // physical-mm float error and local/cut round trips (including
                // mirror winding), so do NOT reintroduce an absolute WORLD float.
                ctl.update_status(94, request->status);
                check_cancel();
                TriangleMesh hull = local.convex_hull_3d();
                if (!usable_mesh(hull.its))
                    throw std::runtime_error("Unable to prepare a valid convex hull for the Boolean result.");
                const FaceExtrusionPlacementBounds placement =
                    face_extrusion_replacement_placement(*request, local, hull, check_cancel);
                request->post_bounds = placement.visual_bounds;
                request->post_min_z = placement.min_z;
                request->post_max_z = placement.max_z;
                if (!face_extrusion_bed_safe(request->model_parts, request->post_min_z, request->post_max_z)) {
                    request->error = "The normal build-plate refresh would move this object. This operation is refused, even for a tiny non-zero shift; no coordinates or transforms were adjusted.";
                    return;
                }
                check_cancel();
                request->result = std::make_shared<const TriangleMesh>(std::move(local.its));
                request->hull = std::make_shared<const TriangleMesh>(std::move(hull));
                ctl.update_status(100, request->status);
            } catch (const BooleanCanceled &) {
                request->canceled = true;
            }
        }, [lifetime, request, self = this](bool canceled, std::exception_ptr &eptr) {
            // MUST consume exceptions even when this gizmo was closed/destroyed
            // or superseded: the worker rethrows unhandled eptr on the UI thread.
            const std::exception_ptr error = std::exchange(eptr, nullptr);
            const auto alive = lifetime.lock();
            if (!alive || !alive->alive || alive->epoch != request->epoch)
                return; // do not dereference self before this lifetime/epoch guard
            try {
                std::string message;
                if (error) {
                    try { std::rethrow_exception(error); }
                    catch (const std::exception &e) { message = e.what(); }
                    catch (...) { message = "Unknown Boolean computation error."; }
                }
                self->finish_boolean(request, canceled, message);
            } catch (...) {
                // A successful commit advances epoch BEFORE callbacks. Never
                // access self after callbacks, including from an exception path.
                if (alive->alive && alive->epoch == request->epoch && self->m_boolean_request == request) {
                    self->m_boolean_request.reset();
                    self->m_boolean_error = _u8L("Unable to apply the Boolean result. The original part is unchanged.");
                    self->m_parent.request_extra_frame();
                }
            }
        });
        if (!queued) {
            invalidate_boolean();
            m_boolean_error = _u8L("The operation could not be queued. The original part is unchanged.");
        }
    } catch (const std::exception &) {
        invalidate_boolean();
        m_boolean_error = _u8L("Unable to prepare the Boolean operation. The original part is unchanged.");
    }
    m_parent.request_extra_frame();
}

void GLGizmoFaceExtrude::finish_boolean(const std::shared_ptr<BooleanRequest> &request, bool canceled, const std::string &exception)
{
    if (m_boolean_request != request || m_state != On || m_boolean_lifetime->epoch != request->epoch)
        return;
    const auto reject = [&](const std::string &message) {
        m_boolean_request.reset();
        m_boolean_error = message;
        request_preview_redraw();
    };
    if (canceled || request->canceled || request->cancel.load(std::memory_order_relaxed)) {
        reject(_u8L("Calculation canceled. The original part is unchanged."));
        return;
    }
    if (!exception.empty() || !request->error.empty()) {
        reject(_u8L("Boolean operation not applied. The original part is unchanged.") + "\n" +
               (exception.empty() ? request->error : exception));
        return;
    }
    if (request->empty) {
        reject(_u8L("No solid would remain. The original part was left unchanged."));
        return;
    }
    ModelObject *object;
    ModelInstance *instance;
    ModelVolume *source;
    int object_index = -1;
    const bool result_ready = request->prepare_stretch ? bool(request->stretch_plan) : bool(request->result && request->hull);
    bool current = result_ready && !m_committing && !m_arrow_drag &&
        m_tool_mode == request->tool_mode && m_output_mode == request->mode && m_distance == request->distance &&
        m_profile_contours == request->profile_contours && m_plane_position == request->plane_position &&
        m_cut_phase == request->cut_phase && m_plane_basis == request->plane_basis &&
        m_cut_guided == request->cut_guided &&
        (m_plane_normal_world.array() == request->plane_normal.array()).all() &&
        (m_plane_reference_world.array() == request->plane_reference.array()).all() &&
        m_plane_source_facet == request->plane_source_facet &&
        bool(m_plane_reference_local) == bool(request->plane_reference_local) &&
        (!m_plane_reference_local || (m_plane_reference_local->array() == request->plane_reference_local->array()).all()) &&
        m_moving_opposite == request->moving_opposite && m_moving_side_chosen == request->moving_side_chosen &&
        (m_tool_mode != ToolMode::CutAndStretch ||
            ((m_plane_origin_world.array() == request->plane_origin.array()).all() &&
             (m_normal_world.array() == request->moving_normal.array()).all() &&
             (request->prepare_stretch || (m_stretch_plan == request->stretch_plan &&
              m_cut_scope_confirmation.matches(request->stretch_plan, cut_scope(request->cut_guided)))))) &&
        resolve_source(object, instance, source, &object_index);
    if (current) {
        current = object->id() == request->object_id && instance->id() == request->instance_id && source->id() == request->source_id &&
            source->get_mesh_shared_ptr() == request->source_identity &&
            face_extrusion_source_type_matches(m_tool_mode, request->source_type, *source) && boolean_source_limit(*source).empty() &&
            (instance->get_matrix().matrix().array() == request->instance_matrix.matrix().array()).all() &&
            (instance->get_matrix(true).matrix().array() == request->instance_without_translation.matrix().array()).all() &&
            instance->get_offset(Z) == request->instance_offset_z &&
            (source->get_matrix().matrix().array() == request->volume_matrix.matrix().array()).all() &&
            object->volumes.size() == request->volumes.size() && object->parts_count() == request->model_parts &&
            m_selected.size() == request->faces.size() && m_selected.size() == request->guide_hits.size() &&
            m_selected.size() == request->guide_normals.size() && m_selected.size() == request->guide_facets.size();
        for (size_t i = 0; current && i < object->volumes.size(); ++i) {
            const ModelVolume *volume = object->volumes[i];
            const auto &stamp = request->volumes[i];
            current = volume->id() == stamp.id && volume->type() == stamp.type && volume->get_mesh_shared_ptr() == stamp.identity &&
                (!volume->is_model_part() || volume->get_convex_hull_shared_ptr() == stamp.hull_identity) &&
                (volume->get_matrix().matrix().array() == stamp.matrix.matrix().array()).all();
        }
        for (size_t i = 0; current && i < m_selected.size(); ++i)
            current = request->faces[i] == std::make_pair(m_selected[i].seed_facet, m_selected[i].plane_id) &&
                (request->guide_hits[i].array() == m_selected[i].hit_local.array()).all() &&
                (request->guide_normals[i].array() == m_selected[i].normal_local.array()).all() &&
                request->guide_facets[i] == m_selected[i].triangles;
    }
    if (!current) {
        reject(_u8L("The target changed during calculation. Nothing was applied; select the faces again."));
        return;
    }
    Plater *plater = wxGetApp().plater();
    if (!plater || !wxGetApp().obj_list() || plater->printer_technology() != ptFFF) {
        reject(_u8L("The editing context changed. Nothing was applied."));
        return;
    }
    if (const std::string limit = project_exterior_limit(plater->model()); !limit.empty()) {
        reject(limit); // protect annotations on any other object/volume as well
        return;
    }
    std::shared_ptr<const TriangleMesh> prepared_result = request->result;
    try {
        const bool consent_current = !request->discard_annotations ||
            (m_annotations && m_annotations->source_id == source->id() &&
             m_annotations->mesh == source->get_mesh_shared_ptr() && m_annotations->consent &&
             *m_annotations->consent == request->annotation_stamp);
        if (!consent_current || capture_face_extrusion_annotation_stamp(*source) != request->annotation_stamp) {
            update_annotation_state(*source, true);
            m_annotations->consent.reset();
            reject(_u8L("Source annotations or consent changed during calculation. Nothing was applied. Review the detected annotations; selected faces are kept."));
            return;
        }
        if (request->prepare_stretch) {
            m_stretch_plan = request->stretch_plan; // CPU data only; GL upload deferred to render/UI
            m_cut_scope_confirmation.clear(); // a new region always needs fresh consent
            m_bridge_base_z = std::move(request->bridge_base_z);
            m_boolean_request.reset();
            m_prepare_requested = false;
            m_preview_dirty = true;
            m_boolean_error.clear();
            request_preview_redraw();
            return; // preparation NEVER takes a Snapshot or changes any live volume
        }
        if (request->discard_annotations) {
            // Stage any legacy-property COW allocation on a PRIVATE candidate.
            // No instances => no backup task; no live Model, paint or Undo change.
            // Candidate metadata is never assigned back to the live volume.
            Model scratch;
            ModelVolume *candidate = scratch.add_object()->add_volume(*source);
            candidate->set_mesh(prepared_result);
            clear_face_extrusion_annotations(*candidate);
            prepared_result = candidate->get_mesh_shared_ptr();
        }
    } catch (const std::exception &) {
        reject(_u8L("Unable to prepare annotation-safe replacement. The original part is unchanged."));
        return;
    }
    if (request->source_type == ModelVolumeType::NEGATIVE_VOLUME &&
        (!face_extrusion_source_type_matches(m_tool_mode, request->source_type, *source) ||
         !face_extrusion_new_parts_bed_safe(request->source_type, *object, *instance, 0.0))) {
        reject(_u8L("Refreshing the negative part would move its positive carrier onto the build plate. Nothing was applied."));
        return;
    }
    // Heavy geometry, hull, float round-trip and actual ensure_on_bed preflight
    // are complete. Live IDs/meshes/transforms and source annotations match the
    // snapshot request, so its POST-operation bounds remain valid here.
    std::optional<Plater::TakeSnapshot> snapshot;
    const ObjectID object_id = object->id();
    ObjectID result_id;
    {
        CommitGuard guard(m_committing);
        snapshot.emplace(plater, request->tool_mode == ToolMode::CutAndStretch ? _u8L("Cut and extend source part") :
            request->tool_mode == ToolMode::PointProfile ?
                (request->mode == OutputMode::Union ? _u8L("Union point profile with source") : _u8L("Subtract point profile from source")) :
            request->source_type == ModelVolumeType::NEGATIVE_VOLUME ?
                (request->mode == OutputMode::Union ? _u8L("Union face extrusions with negative source") : _u8L("Subtract face extrusions from negative source")) :
            request->mode == OutputMode::Union ? _u8L("Union face extrusions with source") : _u8L("Subtract face extrusions from source"));
        restore_negative_visibility();
        m_negative_raycaster.reset();
        // Invalidate epoch and release all cached source/preview references
        // BEFORE any callback. Only local owning result pointers survive.
        clear_selection();
        m_measurements.clear();
        m_raycaster_meshes.clear();
        // Preserve an in-flight camera gesture across this asynchronous commit.
        replace_face_extrusion_source(*source, request->tool_mode, request->source_type,
                                      prepared_result, request->hull, request->discard_annotations);
        result_id = source->id();
    } // guard writes m_committing BEFORE callbacks, not after them
    refresh_boolean_source(object_id, result_id);
    // No member access here or in any post-callback guard/destructor. Snapshot
    // remains alive through reselection, suppressing nested selection snapshots.
}

void GLGizmoFaceExtrude::create_part()
{
    if (m_tool_mode == ToolMode::CutAndStretch || m_committing || m_arrow_drag || boolean_busy() || m_output_mode != OutputMode::NewParts)
        return;
    update_preview();
    if (!m_preview_valid)
        return;
    ModelObject *object;
    ModelInstance *instance;
    ModelVolume *source;
    int object_index = -1;
    if (!resolve_source(object, instance, source, &object_index)) {
        clear_selection();
        m_error = _u8L("The object changed. Select the face again.");
        return;
    }
    Plater *plater = wxGetApp().plater();
    ObjectList *object_list = wxGetApp().obj_list();
    if (!plater || !object_list || !usable_transform(instance->get_matrix()))
        return;
    if (const std::string limit = project_exterior_limit(plater->model()); !limit.empty()) {
        m_boolean_error = limit;
        request_preview_redraw();
        return; // no Snapshot, data reset, or face-selection reset
    }

    // Preparation phase: all N closed solids and storage-frame round trips
    // must be valid BEFORE any ModelObject mutation or undo snapshot. This
    // NewParts keeps one independent volume per patch, with the source role;
    // explicit Boolean modes have their separate async preflight above.
    if (m_previews.empty() || m_previews.size() != m_selected.size())
        return;
    const Transform3d instance_matrix = instance->get_matrix();
    const ModelVolumeType source_type = source->type();
    const bool negative = source_type == ModelVolumeType::NEGATIVE_VOLUME;
    const std::string base_name = source->name + " - " +
        (m_tool_mode == ToolMode::PointProfile ? _u8L("Point profile extrusion") :
         negative ? _u8L("Negative face extrusion") : _u8L("Face extrusion"));
    const int extruder = negative ? 0 : source->extruder_id();
    std::vector<TriangleMesh> prepared;
    std::vector<std::string> names;
    std::vector<ObjectID> new_volume_ids;
    try {
        prepared.reserve(m_previews.size());
        names.reserve(m_previews.size());
        new_volume_ids.reserve(m_previews.size());
        for (const PreviewPart &part : m_previews) {
            auto preparation = prepare_face_extrusion_part_storage(part.world_mesh, instance_matrix, bed_tolerance_mm);
            if (!preparation) {
                m_preview_valid = false;
                m_error = _u8L(preparation.message.c_str());
                return;
            }
            prepared.push_back(std::move(preparation.mesh));
            names.push_back(base_name + (m_previews.size() > 1 ? " " + std::to_string(names.size() + 1) : ""));
        }
    } catch (const std::exception &) {
        m_preview_valid = false;
        m_error = _u8L("Unable to prepare the new parts. The original object is unchanged.");
        return;
    }

    // A preview may be cached while an unrelated carrier changes. Re-resolve
    // the explicit source role/mesh/pose and native carrier placement NOW, before
    // any snapshot: new negatives do not increase the printable part count.
    ModelObject *live_object;
    ModelInstance *live_instance;
    ModelVolume *live_source;
    if (!resolve_source(live_object, live_instance, live_source) || live_object != object ||
        live_instance != instance || live_source != source ||
        !face_extrusion_source_type_matches(m_tool_mode, source_type, *live_source) ||
        !(live_instance->get_matrix().matrix().array() == instance_matrix.matrix().array()).all() ||
        !((live_instance->get_matrix() * live_source->get_matrix()).matrix().array() == m_world_matrix.matrix().array()).all()) {
        clear_selection();
        m_error = _u8L("The object changed. Select the face again.");
        return;
    }
    if (negative && !face_extrusion_new_parts_bed_safe(source_type, *live_object, *live_instance, 0.0)) {
        m_preview_valid = false;
        m_preview_dirty = true;
        m_error = _u8L("Creating these negative parts would move their positive carrier onto the build plate.");
        return;
    }
    // Mutation phase: one snapshot for the complete batch. No selected-face or
    // preview references are retained across the subsequent scene reload.
    std::optional<Plater::TakeSnapshot> snapshot;
    const ObjectID object_id = object->id();
    {
        CommitGuard guard(m_committing);
        snapshot.emplace(plater, m_tool_mode == ToolMode::PointProfile ? _u8L("Extrude point profile as new part") :
            negative ? (prepared.size() == 1 ? _u8L("Extrude face as new negative part") : _u8L("Extrude faces as new negative parts")) :
            prepared.size() == 1 ? _u8L("Extrude face as new part") : _u8L("Extrude faces as new parts"));
        restore_negative_visibility();
        m_negative_raycaster.reset();
        for (size_t i = 0; i < prepared.size(); ++i) {
            ModelVolume *created = insert_face_extrusion_part(*object, std::move(prepared[i]), names[i],
                                                            m_tool_mode, source_type, extruder);
            new_volume_ids.push_back(created->id());
        }
        object->invalidate_bounding_box();

        // Sorting/reload happens only after member cleanup and guard teardown.
        clear_selection();
        m_measurements.clear();
        m_raycaster_meshes.clear();
        m_left_gesture = LeftGesture::None;
        m_consume_left_up = false;
    }
    refresh_created_parts(object_id, std::move(new_volume_ids));
    // No member access after callbacks, including from an RAII guard. Snapshot
    // still spans reselection, preserving one Undo for the complete N-part batch.
}

} // namespace GUI
} // namespace Slic3r
