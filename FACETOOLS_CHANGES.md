# FaceTools changes relative to Bambu Studio (base commit 79df7ae, "2.8.3 Beta Build Fixes")

Modifications by Emmanuele Semprini, 2026 (see `git log` for dates). Summary:
- New (libslic3r): FaceExtrusion, PlanarProfile, CutAndStretch, CutGuideSelection, CutPlaneSelection, PointAlignment, MeshSolidValidation, CheckedStlExport, FFFStlExportSafety, annotations/diagnostics headers.
- New (GUI): GLGizmoFaceExtrude gizmo, FacetPicker, scope-confirmation and preview-bounds helpers, native-cut output policy.
- Modified: MeshBoolean*, GLGizmoAdvancedCut (cut from face / cut and extend), GLGizmoAssembly, GLGizmoMeasure, GLGizmoMeshBoolean, GLGizmosManager, GLCanvas3D, Plater (undo/redo and coincidence placement), Jobs/Job.cpp, PrintJob.cpp and SelectMachine.cpp (worker-to-UI error dispatch for the Send dialog; no change to plug-in, signature checks or protocol).
- Tests: tests/libslic3r/test_face_extrusion*.cpp, test_cut_*.cpp, test_planar_profile.cpp, test_point_alignment.cpp, etc.; tests/slic3r_gui/test_job_ui_dispatch.cpp.
- Tools: tools/facetools (fixture generators, boolean probe).
