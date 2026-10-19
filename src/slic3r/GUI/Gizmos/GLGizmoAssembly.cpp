#include "GLGizmoAssembly.hpp"
#include "slic3r/GUI/GLCanvas3D.hpp"
#include "slic3r/GUI/3DScene.hpp"
#include "slic3r/GUI/MeshUtils.hpp"
#include "slic3r/GUI/GUI_App.hpp"
#include "slic3r/GUI/Plater.hpp"
#include "slic3r/GUI/Gizmos/GizmoObjectManipulation.hpp"
#include "slic3r/Utils/UndoRedo.hpp"

#include "libslic3r/PresetBundle.hpp"
#include "libslic3r/MeasureUtils.hpp"

#include <imgui/imgui_internal.h>

#include <numeric>

#include <GL/glew.h>

#include <tbb/parallel_for.h>
#include <future>
#include <wx/clipbrd.h>

namespace Slic3r {
namespace GUI {

GLGizmoAssembly::GLGizmoAssembly(GLCanvas3D& parent, unsigned int sprite_id) :
    GLGizmoMeasure(parent, sprite_id)
{
    m_measure_mode       = EMeasureMode::ONLY_ASSEMBLY;
}

std::string GLGizmoAssembly::get_icon_filename(bool is_dark_mode) const
{
    return is_dark_mode ? "toolbar_assembly_dark.svg" : "toolbar_assembly.svg";
}

bool GLGizmoAssembly::on_mouse(const wxMouseEvent &mouse_event)
{
    if (is_prepare_point_mode() && mouse_event.LeftDown() && !mouse_event.CmdDown() &&
        m_mode == EMode::PointSelection && m_hover_id == PLANE_ID && m_curr_feature &&
        m_curr_feature->get_type() == Measure::SurfaceFeatureType::Plane) {
        // Native Shift/Plane picking hits the raised highlight (0.05 local mm),
        // not the source surface. Recast the REAL mesh and require the same plane;
        // never project a point from a different plane or silently use the offset.
        std::string reason;
        const auto reject = [&] {
            m_point_message = reason.empty() ? _u8L("The real surface no longer matches the highlighted plane. Pick the point again.") : reason;
            m_mouse_left_down = true;
            m_mouse_left_down_mesh_deal = false;
            m_parent.set_as_dirty();
            m_parent.request_extra_frame();
            return true;
        };
        auto *volume = static_cast<GLVolume *>(m_curr_feature->volume);
        const auto &live = m_parent.get_volumes().volumes;
        if (!volume || std::find(live.begin(), live.end(), volume) == live.end() || !point_scene_is_current(reason))
            return reject();
        const auto raycaster = m_mesh_raycaster_map.find(volume);
        const auto measuring = m_mesh_measure_map.find(volume);
        if (raycaster == m_mesh_raycaster_map.end() || !raycaster->second || !raycaster->second->mesh_raycaster ||
            measuring == m_mesh_measure_map.end() || !measuring->second ||
            !point_alignment_same_matrix(volume->world_matrix(), m_curr_feature->world_tran))
            return reject();
        Vec3f hit, normal;
        size_t facet = 0;
        if (!raycaster->second->mesh_raycaster->unproject_on_mesh(m_parent.get_local_mouse_position(), volume->world_matrix(),
                wxGetApp().plater()->get_camera(), hit, normal, nullptr, &facet, false))
            return reject();
        const auto plane = measuring->second->get_feature(facet, hit.cast<double>(), Transform3d::Identity(), true);
        if (!plane || plane->get_type() != Measure::SurfaceFeatureType::Plane ||
            std::get<0>(plane->get_plane()) != std::get<0>(m_curr_feature->get_plane()))
            return reject();
        m_curr_point_on_feature_position = volume->world_matrix() * hit.cast<double>();
    }
    return GLGizmoMeasure::on_mouse(mouse_event);
}

bool GLGizmoAssembly::on_init() {
    GLGizmoMeasure::on_init();
    m_shortcut_key = WXK_CONTROL_Y;
    return true;
}

std::string GLGizmoAssembly::on_get_name() const
{
    if (!on_is_activable() && m_state == EState::Off) {
        if (wxGetApp().plater()->canvas3D()->get_canvas_type() == GLCanvas3D::ECanvasType::CanvasAssembleView) {
            return _u8L("Assemble") + ":\n" + _u8L("Please confirm explosion ratio = 1 and select at least two volumes.");
        }
        else {
            return _u8L("Assemble") + ":\n" + _u8L("Please select at least two volumes.");
        }
    } else {
        return _u8L("Assemble");
    }
}

bool GLGizmoAssembly::on_is_activable() const
{
    const Selection& selection = m_parent.get_selection();
    if (selection.is_wipe_tower()) {
        return false;
    }
    const int    selection_volumes_count = 2;
    if (wxGetApp().plater()->canvas3D()->get_canvas_type() == GLCanvas3D::ECanvasType::CanvasAssembleView) {
        if (abs(m_parent.get_explosion_ratio() - 1.0f) < 1e-2 && selection.volumes_count() >= selection_volumes_count) {
            return true;
        }
        return false;
    } else {
        return selection.volumes_count() >= selection_volumes_count;
    }
}

void GLGizmoAssembly::on_render_input_window(float x, float y, float bottom_limit)
{
    static std::optional<Measure::SurfaceFeature> last_feature;
    static EMode last_mode = EMode::FeatureSelection;
    static SelectedFeatures last_selected_features;

    static float last_y = 0.0f;
    static float last_h = 0.0f;

    if (is_prepare_point_mode()) {
        // Never dereference an old GLVolume from the native picking vectors.
        const auto &live = m_parent.get_volumes().volumes;
        const auto gone = [&](const GLVolume *v) { return !v || std::find(live.begin(), live.end(), v) == live.end(); };
        if (std::any_of(m_hit_different_volumes.begin(), m_hit_different_volumes.end(), gone) ||
            std::any_of(m_hit_order_volumes.begin(), m_hit_order_volumes.end(), gone)) {
            reset_all_feature();
            m_point_message = _u8L("The picked scene changed. Pick the two points again.");
        }
    }
    if (m_editing_distance)
        return;
    m_current_active_imgui_id      = ImGui::GetActiveID();
    // adjust window position to avoid overlap the view toolbar
    const float win_h = ImGui::GetWindowHeight();
    y = std::min(y, bottom_limit - win_h);
    GizmoImguiSetNextWIndowPos(x, y, ImGuiCond_Always, 0.0f, 0.0f);
    if (last_h != win_h || last_y != y) {
        // ask canvas for another frame to render the window in the correct position
        m_imgui->set_requires_extra_frame();
        if (last_h != win_h)
            last_h = win_h;
        if (last_y != y)
            last_y = y;
    }
    // Orca
    ImGuiWrapper::push_toolbar_style(m_parent.get_scale());
    GizmoImguiBegin(get_name(), ImGuiWindowFlags_NoMove | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoTitleBar);
    init_render_input_window();

    float moving_size = m_imgui->calc_text_size(_L("(Moving)")).x;
    float combo_text_width = std::max(m_imgui->calc_text_size(_L("Face and face assembly")).x,
                                      m_imgui->calc_text_size(_L("Point and point assembly")).x);
    float combox_content_size = combo_text_width * 1.2f + m_imgui->calc_text_size("xxx").x + m_imgui->scaled(3.0f);
    float caption_size = moving_size + 2 * m_space_size;
    if (render_assembly_mode_combo(caption_size + 0.5 * m_space_size,  combox_content_size)) {
        ;
    }
    show_selection_ui();
    std::optional<Vec3d> point_apply;
    const bool native_circle_pair = m_selected_features.first.feature && m_selected_features.second.feature &&
        m_selected_features.first.feature->get_type() == Measure::SurfaceFeatureType::Circle &&
        m_selected_features.second.feature->get_type() == Measure::SurfaceFeatureType::Circle;
    if (is_prepare_point_mode() && !native_circle_pair) {
        point_apply = show_prepare_point_alignment();
    } else {
        // Native Assemble-tab poses, face alignment and raw Circle-Circle XYZ stay unchanged.
        show_face_face_assembly_common();
        ImGui::Separator();
        show_face_face_assembly_senior();
        show_distance_xyz_ui();
        render_input_window_warning(m_same_model_object);
        if (is_prepare_point_mode() && native_circle_pair) {
            m_imgui->disabled_begin(true);
            m_imgui->button(_L("Coincide points"));
            m_imgui->disabled_end();
            m_imgui->text_wrapped(_u8L("Select each circle's center as a Point to use point coincidence. Native circle-distance controls are unchanged."), m_imgui->scaled(28.0f));
        }
    }
    ImGui::Separator();

    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(6.0f, 10.0f));
    float get_cur_y = ImGui::GetContentRegionMax().y + ImGui::GetFrameHeight() + y;
    float caption_max    = 0.f;
    float total_text_max = 0.f;
    for (const auto &t : std::array<std::string, 3>{"point_selection", "reset", "unselect"}) {
        caption_max    = std::max(caption_max, m_imgui->calc_text_size(m_desc[t + "_caption"]).x);
        total_text_max = std::max(total_text_max, m_imgui->calc_text_size(m_desc[t]).x);
    }
    show_tooltip_information(caption_max, x, get_cur_y);

    float f_scale = m_parent.get_main_toolbar_scale();
    ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(6.0f, 4.0f * f_scale));

    ImGui::PopStyleVar(2);

    if (last_feature != m_curr_feature || last_mode != m_mode || last_selected_features != m_selected_features) {
        // the dialog may have changed its size, ask for an extra frame to render it properly
        last_feature = m_curr_feature;
        last_mode = m_mode;
        last_selected_features = m_selected_features;
        m_imgui->set_requires_extra_frame();
    }
    m_last_active_item_imgui = m_current_active_imgui_id;
    GizmoImguiEnd();
    // Orca
    ImGuiWrapper::pop_toolbar_style();
    // A scene commit can trigger callbacks/reselection: close ImGui first and do
    // not access gizmo members after the final native notification.
    if (point_apply)
        apply_point_alignment(*point_apply);
}

void GLGizmoAssembly::render_input_window_warning(bool same_model_object)
{
    if (wxGetApp().plater()->canvas3D()->get_canvas_type() == GLCanvas3D::ECanvasType::CanvasView3D) {
        if (m_hit_different_volumes.size() == 2) {
            if (same_model_object == false) {
                m_imgui->warning_text(_L("Warning") + ": " +
               _L("It is recommended to assemble the objects first,\nbecause the objects is restriced to bed \nand only parts can be lifted."));
            }
        }
    }
}

bool GLGizmoAssembly::render_assembly_mode_combo(double label_width, float item_width)
{
    ImGui::AlignTextToFramePadding();
    size_t                   selection_idx = int(m_assembly_mode);
    std::vector<std::string> modes         = {_u8L("Face and face assembly"), _u8L("Point and point assembly")};
    bool                     is_changed    = false;

    ImGuiWrapper::push_combo_style(m_parent.get_scale());
    if (render_combo(_u8L("Mode"), modes, selection_idx, label_width, item_width)) {
        is_changed = true;
        switch_to_mode((AssemblyMode) selection_idx);
    }
    ImGuiWrapper::pop_combo_style();
    return is_changed;
}

void GLGizmoAssembly::switch_to_mode(AssemblyMode new_mode)
{
    m_assembly_mode = new_mode;
    reset_all_feature();
    m_alignment_picks.reset();
    m_point_anchors = {};
    m_point_volumes = {{nullptr, nullptr}};
    m_point_buffer_initialized = false;
    m_point_message.clear();
}

bool GLGizmoAssembly::is_prepare_point_mode() const
{
    return m_parent.get_canvas_type() == GLCanvas3D::ECanvasType::CanvasView3D &&
           m_assembly_mode == AssemblyMode::POINT_POINT;
}

bool GLGizmoAssembly::same_point_pick(const SelectedFeatures::Item &a, const SelectedFeatures::Item &b)
{
    const auto same = [](const std::optional<Measure::SurfaceFeature> &x, const std::optional<Measure::SurfaceFeature> &y) {
        if (x.has_value() != y.has_value()) return false;
        if (!x) return true;
        return x->get_type() == y->get_type() && x->volume == y->volume &&
               x->origin_surface_feature == y->origin_surface_feature && x->get_pt1() == y->get_pt1() &&
               x->get_pt2() == y->get_pt2() && x->get_pt3() == y->get_pt3() && x->get_value() == y->get_value() &&
               (!x->volume || point_alignment_same_matrix(x->world_tran, y->world_tran));
    };
    return a.is_center == b.is_center && same(a.feature, b.feature) && same(a.source, b.source);
}

bool GLGizmoAssembly::resolve_point_pick(const SelectedFeatures::Item &item, PointAlignmentAnchor &anchor,
                                         GLVolume *&volume, std::string &reason) const
{
    const auto fail = [&](const char *text) { reason = _u8L(text); return false; };
    if (!item.feature || item.feature->get_type() != Measure::SurfaceFeatureType::Point)
        return fail("Select two Point features. Use Shift to pick a point on a face, edge or circle.");
    const auto &point = *item.feature;
    const auto *owner = point.volume ? &point : (item.source ? &*item.source : nullptr);
    if (!owner || !owner->volume || (item.source && item.source->volume && item.source->volume != owner->volume))
        return fail("The point and its source feature do not belong to the same part. Pick it again.");
    volume = static_cast<GLVolume *>(owner->volume);
    const auto &live = m_parent.get_volumes().volumes;
    const auto found = std::find(live.begin(), live.end(), volume);
    if (found == live.end() || !m_parent.get_selection().contains_volume(unsigned(found - live.begin())))
        return fail("The picked part is no longer selected. Pick the points again.");
    const Model *model = m_parent.get_model();
    const int oi = volume->object_idx(), vi = volume->volume_idx(), ii = volume->instance_idx();
    if (!model || oi < 0 || size_t(oi) >= model->objects.size() || !model->objects[oi])
        return fail("The picked object is no longer available.");
    const ModelObject &object = *model->objects[oi];
    if (vi < 0 || ii < 0 || size_t(vi) >= object.volumes.size() || size_t(ii) >= object.instances.size() ||
        !object.volumes[vi] || !object.instances[ii])
        return fail("The picked part or instance is no longer available.");
    const ModelVolume &part = *object.volumes[vi];
    const ModelInstance &instance = *object.instances[ii];
    const Transform3d world = instance.get_matrix() * part.get_matrix();
    if (!part.is_model_part() || !part.mesh_ptr() || volume->ori_mesh != part.mesh_ptr() ||
        volume->get_sla_shift_z() != 0 || !point_alignment_valid_transform(world) ||
        !point_alignment_same_matrix(part.get_matrix(), volume->get_volume_transformation().get_matrix()) ||
        !point_alignment_same_matrix(instance.get_matrix(), volume->get_instance_transformation().get_matrix()) ||
        !point_alignment_same_matrix(world, owner->world_tran))
        return fail("The picked mesh or displayed pose changed. Pick the points again.");
    Vec3d local;
    if (point.origin_surface_feature && point.origin_surface_feature->get_type() == Measure::SurfaceFeatureType::Point)
        local = point.origin_surface_feature->get_point();
    else if (item.is_center && item.source)
        local = world.inverse() * point.get_point();
    else
        return fail("The picked point has no valid local point reference. Pick it again.");
    if (item.source && item.source->get_type() == Measure::SurfaceFeatureType::Plane) {
        const auto &origin = item.source->origin_surface_feature;
        if (!origin || origin->get_type() != Measure::SurfaceFeatureType::Plane ||
            std::get<0>(origin->get_plane()) != std::get<0>(item.source->get_plane()))
            return fail("The point's source plane changed. Pick it again.");
        const auto [plane_id, normal, plane_point] = origin->get_plane();
        if (!point_alignment_plane_matches(plane_id, std::get<0>(item.source->get_plane()), normal, plane_point, local))
            return fail("The picked point no longer lies on its source plane.");
    }
    anchor = point_alignment_anchor(object, instance, part, local);
    anchor.world_point = point.get_point();
    if (!local.allFinite() || !anchor.world_point.allFinite() || (world * local - anchor.world_point).norm() > 1e-8)
        return fail("The picked local/world coordinates disagree. Pick the point again.");
    return true;
}

bool GLGizmoAssembly::sync_point_anchors(std::string &reason)
{
    const SelectedFeatures::Item *items[] = {&m_selected_features.first, &m_selected_features.second};
    SelectedFeatures::Item *previous[] = {&m_alignment_picks.first, &m_alignment_picks.second};
    bool ready = true;
    for (size_t i = 0; i < 2; ++i) {
        if (!same_point_pick(*items[i], *previous[i])) {
            *previous[i] = *items[i];
            m_point_anchors[i].reset();
            m_point_volumes[i] = nullptr;
            m_point_buffer_initialized = false;
            m_point_message.clear();
        }
        PointAlignmentAnchor current;
        GLVolume *volume = nullptr;
        std::string error;
        if (!resolve_point_pick(*items[i], current, volume, error)) {
            if (reason.empty()) reason = error;
            ready = false;
            continue;
        }
        if (!m_point_anchors[i]) {
            m_point_anchors[i] = current;
            m_point_volumes[i] = volume;
        } else {
            const auto &saved = *m_point_anchors[i];
            if (volume != m_point_volumes[i] || current.object_id != saved.object_id || current.instance_id != saved.instance_id ||
                current.volume_id != saved.volume_id || current.mesh != saved.mesh ||
                !point_alignment_same_matrix(current.instance_matrix, saved.instance_matrix) ||
                !point_alignment_same_matrix(current.volume_matrix, saved.volume_matrix)) {
                reason = _u8L("The selected mesh, ID or pose changed. Pick the points again.");
                ready = false;
            }
        }
    }
    if (!ready) return false;
    const auto &fixed = *m_point_anchors[0];
    const auto &moving = *m_point_anchors[1];
    if (fixed.object_id != moving.object_id || fixed.instance_id != moving.instance_id || fixed.volume_id == moving.volume_id) {
        reason = _u8L("Select two different parts of the SAME object and instance. Point 1 is fixed; only Point 2's part moves.");
        return false;
    }
    const ModelObject *object = m_parent.get_model()->objects[m_point_volumes[0]->object_idx()];
    if (object->instances.size() != 1) {
        reason = _u8L("This object has multiple instances. Use a single instance so its other copies cannot move.");
        return false;
    }
    if (m_hit_different_volumes.size() != 2 || m_hit_order_volumes.size() != 2 ||
        m_hit_different_volumes[0] != m_point_volumes[0] || m_hit_different_volumes[1] != m_point_volumes[1] ||
        m_hit_order_volumes[0] != m_point_volumes[0] || m_hit_order_volumes[1] != m_point_volumes[1]) {
        reason = _u8L("The point selection order changed. Pick Point 1 (Fixed), then Point 2 (Moving).");
        return false;
    }
    return true;
}

bool GLGizmoAssembly::point_scene_is_current(std::string &reason) const
{
    const Model *model = m_parent.get_model();
    if (!model || m_parent.get_selection().get_model() != model) {
        reason = _u8L("The selection belongs to a different scene. Pick the points again.");
        return false;
    }
    for (const GLVolume *volume : m_parent.get_volumes().volumes) {
        if (!volume) { reason = _u8L("The scene contains an unavailable volume."); return false; }
        const int oi = volume->object_idx(), vi = volume->volume_idx(), ii = volume->instance_idx();
        if (oi < 0 || size_t(oi) >= model->objects.size()) continue; // native wipe-tower proxies
        const ModelObject *object = model->objects[oi];
        if (!object || vi < 0 || ii < 0 || size_t(vi) >= object->volumes.size() || size_t(ii) >= object->instances.size() ||
            !object->volumes[vi] || !object->instances[ii] || volume->ori_mesh != object->volumes[vi]->mesh_ptr() ||
            !point_alignment_same_matrix(volume->get_volume_transformation().get_matrix(), object->volumes[vi]->get_matrix()) ||
            !point_alignment_same_matrix(volume->get_instance_transformation().get_matrix(), object->instances[ii]->get_matrix())) {
            reason = _u8L("Another displayed part has an uncommitted pose. Finish that edit before aligning points.");
            return false;
        }
        const auto cached = m_mesh_measure_map.find(const_cast<GLVolume *>(volume));
        if (cached != m_mesh_measure_map.end() && (!cached->second || !object->volumes[vi]->mesh_ptr() ||
            cached->second->get_its().vertices != object->volumes[vi]->mesh().its.vertices ||
            cached->second->get_its().indices != object->volumes[vi]->mesh().its.indices)) {
            reason = _u8L("The picked mesh cache changed. Reopen the tool and pick the points again.");
            return false;
        }
    }
    return true;
}

std::optional<Vec3d> GLGizmoAssembly::show_prepare_point_alignment()
{
    std::string reason;
    bool ready = sync_point_anchors(reason);
    const Model *model = m_parent.get_model();
    if (model && point_alignment_has_exterior(*model)) {
        reason = _u8L("This project contains exterior-facet data that Undo cannot restore. Point moves are disabled; existing data is unchanged.");
        ready = false;
    }
    if (wxGetApp().plater()->printer_technology() != ptFFF) {
        reason = _u8L("Point coincidence in Prepare currently requires an FFF project.");
        ready = false;
    }
    const float width = m_imgui->scaled(28.0f);
    ImGui::Separator();
    m_imgui->text_wrapped(_u8L("Point 1 stays FIXED. Only the part containing Point 2 MOVES, by translation only (no rotation or scaling)."), width);
    m_imgui->text_wrapped(_u8L("Target offset = Point 2 minus Point 1, in world XYZ. Enter or Apply commits all three components together."), width);
    if (ready && (!m_point_buffer_initialized || m_point_buffer_inches != m_use_inches)) {
        m_point_gap_buffer = m_point_anchors[1]->world_point - m_point_anchors[0]->world_point;
        if (m_use_inches) m_point_gap_buffer *= GizmoObjectManipulation::mm_to_in;
        m_point_buffer_initialized = true;
        m_point_buffer_inches = m_use_inches;
    }
    std::optional<Vec3d> apply;
    bool enter = false;
    m_imgui->disabled_begin(!ready);
    const char *axis_labels[] = {"X", "Y", "Z"};
    for (int axis = 0; axis < 3; ++axis) {
        ImGui::SetNextItemWidth(m_imgui->scaled(12.0f));
        const std::string label = std::string(axis_labels[axis]) + m_units + "##point_alignment_" + axis_labels[axis];
        enter |= ImGui::InputDouble(label.c_str(), &m_point_gap_buffer[axis], 0.0, 0.0, "%.9g", ImGuiInputTextFlags_EnterReturnsTrue);
    }
    const bool apply_xyz = m_imgui->button(_L("Apply XYZ offset"));
    if ((enter || apply_xyz) && ready)
        apply = m_use_inches ? Vec3d(m_point_gap_buffer * GizmoObjectManipulation::in_to_mm) : m_point_gap_buffer;
    ImGui::SameLine();
    if (m_imgui->button(_L("Coincide points")) && ready)
        apply = Vec3d::Zero(); // never commits an intermediate scalar edit
    m_imgui->disabled_end();
    if (!reason.empty()) m_imgui->warning_text(reason);
    if (!m_point_message.empty()) m_imgui->text_wrapped(m_point_message, width);
    return apply;
}

void GLGizmoAssembly::apply_point_alignment(const Vec3d &gap_mm)
{
    if (!is_prepare_point_mode() || m_state != On) return;
    m_parent.set_as_dirty();
    m_parent.request_extra_frame(); // show failure/no-op feedback without waiting for mouse movement
    Plater *plater = wxGetApp().plater();
    std::string reason;
    if (!plater || plater->printer_technology() != ptFFF || m_parent.get_model() != &plater->model() ||
        !sync_point_anchors(reason) || !point_scene_is_current(reason)) {
        m_point_message = reason.empty() ? _u8L("The editing context changed; nothing was moved.") : reason;
        return;
    }
    const PointAlignmentPlan plan = prepare_point_alignment(plater->model(), *m_point_anchors[0], *m_point_anchors[1], gap_mm);
    if (!plan) { m_point_message = _u8L(plan.message); return; }
    if (plan.no_op) {
        m_point_message = _u8L("The points already have this offset. No move or Undo entry was created.");
        return;
    }
    GLVolume *moving_volume = m_point_volumes[1];
    // Allocate/copy point state before Snapshot. Center picks also acquire a proper
    // local point reference, so later native measurement updates remain valid.
    SelectedFeatures::Item moved = m_selected_features.second;
    moved.feature->clone(Measure::SurfaceFeature(plan.moving_world_point));
    moved.feature->volume = moving_volume;
    moved.feature->world_tran = m_point_anchors[1]->instance_matrix * plan.moving_matrix;
    if (!moved.feature->origin_surface_feature)
        moved.feature->origin_surface_feature = std::make_shared<Measure::SurfaceFeature>(m_point_anchors[1]->local_point);
    if (moved.source && moved.source->origin_surface_feature) {
        moved.source->clone(*moved.source->origin_surface_feature);
        moved.source->world_tran = moved.feature->world_tran;
        moved.source->translate(moved.source->world_tran);
    }
    Selection &selection = m_parent.get_selection();
    selection.setup_cache();
    Plater::TakeSnapshot snapshot(plater, _u8L("Align picked points"), UndoRedo::SnapshotType::GizmoAction);
    selection.set_mode(Selection::Volume);
    ++m_pending_scale;
    moving_volume->set_volume_offset(plan.moving_matrix.translation());
    m_selected_features.second = std::move(moved);
    m_alignment_picks = m_selected_features;
    m_point_anchors[1]->volume_matrix = plan.moving_matrix;
    m_point_anchors[1]->world_point = plan.moving_world_point;
    m_point_gap_buffer = m_use_inches ? Vec3d(gap_mm * GizmoObjectManipulation::mm_to_in) : gap_mm;
    m_point_buffer_initialized = true;
    m_point_message.clear();
    update_measurement_result();
    register_single_mesh_pick();
    m_parent.set_as_dirty();
    m_parent.request_extra_frame();
    m_parent.do_move(""); // ONE vector commit, no second snapshot, no member access afterwards
}

} // namespace GUI
} // namespace Slic3r
