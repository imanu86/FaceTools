// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#ifndef slic3r_GLGizmoFaceExtrude_hpp_
#define slic3r_GLGizmoFaceExtrude_hpp_

#include "GLGizmoBase.hpp"
#include "FaceExtrusionModelPolicy.hpp"
#include "CutAndStretchScopeConfirmation.hpp"
#include "libslic3r/Measure.hpp"

#include <array>
#include <map>
#include <memory>
#include <optional>
#include <utility>
#include <vector>
#include <wx/weakref.h>
#include <wx/window.h>

namespace Slic3r {
struct CutAndStretchPlan;
class CutGuideSelection;
struct CutGuidePatch;
namespace GUI {

class MeshRaycaster;

// Three separately registered variants share interaction, exact annotation
// consent, async validation and one source-only replacement/Undo pipeline.
// Faces/PointProfile default to new parts; CutAndStretch always replaces source.
class GLGizmoFaceExtrude final : public GLGizmoBase
{
public:
    using ToolMode = FaceExtrusionToolMode;
    GLGizmoFaceExtrude(GLCanvas3D &parent, unsigned int sprite_id, ToolMode mode = ToolMode::Faces);
    ~GLGizmoFaceExtrude() override;

    std::string get_icon_filename(bool is_dark_mode) const override;
    std::string get_tooltip() const override;
    BoundingBoxf3 get_bounding_box() const override;
    bool on_mouse(const wxMouseEvent &event) override;
    bool gizmo_event(SLAGizmoEventType action, const Vec2d &mouse_position,
                     bool shift_down, bool alt_down, bool control_down) override;
    void data_changed(bool is_serializing) override;

protected:
    bool on_init() override;
    std::string on_get_name() const override;
    bool on_is_activable() const override;
    bool on_is_selectable() const override;
    CommonGizmosDataID on_get_requirements() const override;
    void on_set_state() override;
    void on_load(cereal::BinaryInputArchive &) override;
    void on_render() override;
    void on_render_for_picking() override;
    void on_render_input_window(float x, float y, float bottom_limit) override;

private:
    enum class OutputMode { NewParts, Union, Subtract };
    enum class CutPhase { Guides, Plane, Distance };
    enum class PlaneBasis { XY, XZ, YZ, Face };
    struct BooleanRequest;
    struct BooleanLifetime;
    struct AnnotationState;
    enum class LeftGesture { None, PendingFace, Orbit, Arrow, Cancelled };

    struct ArrowDrag {
        Vec3d anchor = Vec3d::Zero();
        Vec3d normal = Vec3d::UnitZ();
        Transform3d source_world = Transform3d::Identity();
        Transform3d view = Transform3d::Identity();
        Transform3d projection = Transform3d::Identity();
        std::array<int, 4> viewport{};
        double initial_parameter = 0.0;
        double initial_distance = 0.0; // unsnapped rebase baseline
        double raw_distance = 0.0;
        double restore_distance = 0.0; // fixed snap origin / Escape value
        Vec2d last_mouse = Vec2d::Zero();
        bool plane_drag = false; // captured at Down; Escape/rebase restore this parameter, not distance
    };

    struct MeasurementCache {
        // Keep the immutable mesh alive, and use its identity together with the
        // ModelVolume ID. GLVolume pointers do not survive scene reloads.
        std::shared_ptr<const TriangleMesh> mesh;
        std::unique_ptr<Measure::Measuring> measuring;
        std::unique_ptr<CutGuideSelection> guide_selection; // private, unpainted CPU selector
        std::string error;
        std::string guide_error;
    };

    struct SelectedFace {
        ObjectID object_id;
        ObjectID instance_id;
        ObjectID volume_id;
        ModelVolumeType type_at_pick = ModelVolumeType::MODEL_PART;
        std::shared_ptr<const TriangleMesh> mesh;
        std::vector<int> triangles;
        Vec3d normal_local = Vec3d::Zero();
        Vec3d hit_local = Vec3d::Zero();
        size_t seed_facet = 0;
        int plane_id = -1;
        int vertex_id = -1; // PointProfile: exact source-mesh vertex, never a free/flattened point
        Transform3d world_at_pick = Transform3d::Identity();
    };

    struct HoverVolume {
        ObjectID id;
        std::shared_ptr<const TriangleMesh> mesh;
        Transform3d world = Transform3d::Identity();
        ModelVolumeType type = ModelVolumeType::MODEL_PART;
        bool operator==(const HoverVolume &other) const {
            return id == other.id && type == other.type && mesh == other.mesh &&
                   (world.matrix().array() == other.world.matrix().array()).all();
        }
    };

    struct HoverQuery {
        ObjectID object_id;
        ObjectID instance_id;
        Vec2d mouse = Vec2d::Zero();
        Transform3d view = Transform3d::Identity();
        Transform3d projection = Transform3d::Identity();
        std::array<int, 4> viewport{};
        std::vector<HoverVolume> volumes;
    };

    struct NegativeRaycasterCache {
        ObjectID object_id, instance_id, volume_id;
        // The mesh outlives its private raycaster, including scene reloads.
        std::shared_ptr<const TriangleMesh> mesh;
        std::unique_ptr<MeshRaycaster> raycaster;
    };

    struct NegativeVisibility {
        ObjectID object_id, instance_id, volume_id;
        std::vector<std::pair<ObjectID, bool>> previous_visibility;
        ClippingPlane lower_clip;
    };

    bool get_target(ModelObject *&object, ModelInstance *&instance,
                    int *object_index = nullptr, ModelVolume **negative_source = nullptr) const;
    void update_negative_visibility();
    void restore_negative_visibility();
    bool resolve_source(ModelObject *&object, ModelInstance *&instance,
                        ModelVolume *&volume, int *object_index = nullptr) const;
    Measure::Measuring *get_measuring(const ModelVolume &volume, std::string &error);
    std::optional<SelectedFace> pick_face(const Vec2d &mouse_position, std::string &error, bool surface_only = false);
    bool load_face_triangles(SelectedFace &face, std::string &error) const;
    static bool same_source(const SelectedFace &a, const SelectedFace &b);
    static bool same_face(const SelectedFace &a, const SelectedFace &b);
    CutGuideSelection *get_cut_guide_selection(const SelectedFace &face, std::string &error);
    std::vector<CutGuidePatch> selected_cut_patches() const;
    void select_cut_region(SelectedFace face);
    void update_cut_region_hover(SelectedFace face);
    void select_face(const Vec2d &mouse_position, bool additive);
    void clear_selection();
    void invalidate_preview();
    void clear_hover(bool invalidate_query = true);
    void update_hover();
    void update_preview();
    void update_variant_preview();
    void set_cut_phase(CutPhase phase);
    void choose_cut_plane_basis(PlaneBasis basis);
    void place_cut_plane(const Vec2d &mouse_position);
    bool use_cut_face_plane(const SelectedFace &face);
    void use_cut_guide_plane(size_t index);
    void update_cut_plane();
    void confirm_cut_plane();
    bool editing_cut_plane() const { return m_tool_mode == ToolMode::CutAndStretch && m_cut_phase == CutPhase::Plane; }
    std::string arrow_parameter_label() const;
    void render_profile_markers();
    void restore_source_visibility();
    void update_stretch_visibility();
    bool arrow_frame(const Transform3d &world, Vec3d &anchor, Vec3d &normal) const;
    void select_profile_vertex(const Vec2d &mouse_position);
    void undo_profile_vertex();
    void update_arrow();
    bool arrow_contains_mouse(const Vec2d &mouse) const;
    bool begin_arrow_drag(const Vec2d &mouse);
    void update_arrow_drag(const Vec2d &mouse);
    void finish_arrow_drag(bool cancel);
    void on_arrow_capture_lost(wxMouseCaptureLostEvent &event);
    void on_arrow_focus_lost(wxFocusEvent &event);
    void render_arrow(bool picking);
    void render_drag_tooltip();
    void request_preview_redraw();
    bool click_threshold_exceeded(const Vec2d &mouse) const;
    bool click_camera_unchanged() const;
    void create_part();
    void start_boolean(bool prepare_stretch = false);
    void finish_boolean(const std::shared_ptr<BooleanRequest> &request, bool canceled, const std::string &exception);
    void invalidate_boolean();
    void update_annotation_state(const ModelVolume &source, bool refresh = false);
    bool boolean_busy() const { return bool(m_boolean_request); }

    struct PreviewPart {
        TriangleMesh world_mesh;
        std::unique_ptr<GLModel> model;
        Transform3d render_world = Transform3d::Identity(); // identity except cached stretch cut-frame operands
        BoundingBoxf3 render_bounds;
    };

    const ToolMode m_tool_mode;
    std::vector<std::vector<int>> m_profile_contours{{}};
    std::vector<std::vector<Vec3d>> m_profile_world;
    std::shared_ptr<GLModel> m_profile_lines;
    bool m_profile_valid = false;
    CutPhase m_cut_phase = CutPhase::Guides;
    bool m_cut_guided = true; // whole-half behavior requires an explicit choice
    CutAndStretchScopeConfirmation m_cut_scope_confirmation;
    bool m_smart_guide_fill = false;
    float m_smart_guide_angle = 30.0f;
    bool m_exact_guide_selection = false; // retain partial masks when smart fill is subsequently disabled
    PlaneBasis m_plane_basis = PlaneBasis::XY;
    bool m_plane_initialized = false;
    bool m_face_plane_ready = true;
    double m_plane_position = 0.0; // absolute world coordinate: Z/Y/X, or n.dot(world point)
    Vec3d m_plane_normal_world = Vec3d::UnitZ();
    Vec3d m_plane_reference_world = Vec3d::Zero(); // near source, never an origin/mesh translation
    std::optional<Vec3d> m_plane_reference_local; // exact picked face/guide point; absent = bbox-center reference
    std::optional<size_t> m_plane_source_facet; // original defining triangle; used ONLY at exactly zero Face offset
    size_t m_cut_guide_index = 0;
    double m_plane_radius = 1.0;
    Transform3d m_plane_render_matrix = Transform3d::Identity();
    std::vector<std::shared_ptr<GLModel>> m_guide_models;
    bool m_moving_opposite = false;
    bool m_moving_side_chosen = false;
    Vec3d m_plane_origin_world = Vec3d::Zero();
    std::shared_ptr<const CutAndStretchPlan> m_stretch_plan;
    std::shared_ptr<const CutAndStretchPlan> m_rendered_stretch_plan;
    indexed_triangle_set m_preview_bridge;
    std::vector<float> m_bridge_base_z;
    std::optional<double> m_bridge_distance;
    bool m_prepare_requested = false;
    std::optional<std::pair<ObjectID, ObjectID>> m_hidden_source;
    std::map<ObjectID, MeasurementCache> m_measurements;
    std::optional<NegativeRaycasterCache> m_negative_raycaster;
    std::optional<NegativeVisibility> m_negative_visibility;
    std::vector<SelectedFace> m_selected;
    // Transient hover is deliberately independent of the locked selection and
    // extrusion preview. Moving the mouse never invokes the extrusion kernel.
    std::optional<SelectedFace> m_hovered;
    std::optional<HoverQuery> m_hover_query;
    std::shared_ptr<GLModel> m_hover_model;
    std::shared_ptr<const std::vector<int>> m_hover_region; // immutable exact candidate mask, not a plane ID
    std::string m_guide_hover_error;
    std::vector<std::shared_ptr<const TriangleMesh>> m_raycaster_meshes;
    TriangleMesh m_world_source;
    std::vector<PreviewPart> m_previews;
    std::shared_ptr<GLModel> m_face_model;
    Transform3d m_world_matrix = Transform3d::Identity();
    Transform3d m_instance_matrix = Transform3d::Identity();
    Vec3d m_normal_world = Vec3d::Zero();
    double m_distance = 5.0;
    OutputMode m_output_mode = OutputMode::NewParts;
    std::shared_ptr<BooleanLifetime> m_boolean_lifetime;
    std::shared_ptr<BooleanRequest> m_boolean_request;
    std::string m_boolean_error;
    std::unique_ptr<AnnotationState> m_annotations;
    bool m_world_dirty = true;
    bool m_preview_dirty = true;
    bool m_preview_valid = false;
    LeftGesture m_left_gesture = LeftGesture::None;
    Vec2d m_press_mouse = Vec2d::Zero();
    Transform3d m_press_view = Transform3d::Identity();
    Transform3d m_press_projection = Transform3d::Identity();
    std::array<int, 4> m_press_viewport{};
    bool m_press_shift = false;
    bool m_consume_left_up = false;
    std::optional<ArrowDrag> m_arrow_drag;
    wxWeakRef<wxWindow> m_capture_window;
    TriangleMesh m_arrow_pick_mesh;
    std::unique_ptr<MeshRaycaster> m_arrow_raycaster;
    Transform3d m_arrow_matrix = Transform3d::Identity();
    Transform3d m_arrow_shaft_matrix = Transform3d::Identity();
    BoundingBoxf3 m_arrow_bounds;
    double m_handle_distance = 0.0;
    bool m_arrow_available = false;
    bool m_arrow_shaft_available = false;
    bool m_arrow_hovered = false;
    bool m_committing = false;
    float m_window_height = 0.0f;
    std::string m_error;
    std::string m_selection_error;
};

} // namespace GUI
} // namespace Slic3r

#endif // slic3r_GLGizmoFaceExtrude_hpp_
