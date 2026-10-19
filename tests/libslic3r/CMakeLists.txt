get_filename_component(_TEST_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

add_executable(${_TEST_NAME}_tests 
	${_TEST_NAME}_tests.cpp
	test_3mf.cpp
	test_aabbindirect.cpp
	test_clipper_offset.cpp
	test_clipper_utils.cpp
	test_config.cpp
	test_elephant_foot_compensation.cpp
	test_geometry.cpp
	test_placeholder_parser.cpp
	test_polygon.cpp
	test_mutable_polygon.cpp
	test_mutable_priority_queue.cpp
	test_stl.cpp
	test_meshboolean.cpp
	test_marchingsquares.cpp
	test_timeutils.cpp
	test_voronoi.cpp
    test_optimizers.cpp
    test_png_io.cpp
    test_timeutils.cpp
    test_indexed_triangle_set.cpp
    test_face_extrusion.cpp
    test_face_extrusion_model.cpp
    test_mesh_quality.cpp
    test_cut_guided_scope.cpp
    test_fff_stl_export_safety.cpp
    test_checked_stl_export.cpp
    test_negative_face_extrusion.cpp
    test_planar_profile.cpp
    test_cut_and_stretch.cpp
    test_cut_and_stretch_placement.cpp
    test_point_alignment.cpp
    test_cut_plane_selection.cpp
    test_cut_plane_selection_workflow.cpp
    test_triangle_selector.cpp
    test_debounce.cpp
    test_emboss.cpp
    ../libnest2d/printer_parts.cpp
	)

if (TARGET OpenVDB::openvdb)
    target_sources(${_TEST_NAME}_tests PRIVATE test_hollowing.cpp)
endif()

target_link_libraries(${_TEST_NAME}_tests test_common libslic3r)
# MeshDiagnostics.hpp currently includes Point.hpp without its directory.
target_include_directories(${_TEST_NAME}_tests PRIVATE ${PROJECT_SOURCE_DIR}/src/libslic3r)
set_property(TARGET ${_TEST_NAME}_tests PROPERTY FOLDER "tests")

if (WIN32)
    # bambustudio_copy_dlls now requires 4 args (target, config, postfix,
    # output_dlls); legacy 1-arg call sites in tests/* predate that signature
    # change.  Test executables link libslic3r and do not need the heavyweight
    # OCCT / GMP / MPFR / FFmpeg DLL bundle, so skip the copy step here.
    # Run tests with PATH pointing at deps build dir if extra DLLs are needed.
    # bambustudio_copy_dlls(${_TEST_NAME}_tests <config> <postfix> <out>)
endif()

# catch_discover_tests(${_TEST_NAME}_tests TEST_PREFIX "${_TEST_NAME}: ")
add_test(${_TEST_NAME}_tests ${_TEST_NAME}_tests ${CATCH_EXTRA_ARGS})

# Keep the feature regressions buildable independently of historical test sources.
# Native BBS 3MF pulls in SVG/logging helpers currently implemented in the GUI
# library. Link the production library and its crypto/HTTP dependencies, not stubs.
# The Catch2 runner does not initialize wxWidgets or enter its event loop.
if (TARGET libslic3r_gui)
    add_executable(face_extrusion_tests
        libslic3r_tests.cpp
        test_face_extrusion.cpp
        test_face_extrusion_model.cpp
        test_mesh_quality.cpp
        test_cut_guided_scope.cpp
        test_fff_stl_export_safety.cpp
        test_checked_stl_export.cpp
        test_negative_face_extrusion.cpp
        test_planar_profile.cpp
        test_cut_and_stretch.cpp
        test_cut_and_stretch_placement.cpp
        test_point_alignment.cpp
        test_cut_plane_selection.cpp
        test_cut_plane_selection_workflow.cpp
        test_cut_plane_preview_bounds.cpp
        test_facet_picker_geometry.cpp
        test_cut_to_parts_placement.cpp
        test_cut_guide_selection.cpp
    )
    target_link_libraries(face_extrusion_tests PRIVATE test_common libslic3r libslic3r_gui)
    target_include_directories(face_extrusion_tests PRIVATE ${PROJECT_SOURCE_DIR}/src/libslic3r)
    set_target_properties(face_extrusion_tests PROPERTIES FOLDER "tests" MACOSX_BUNDLE FALSE)
    add_test(NAME face_extrusion_tests COMMAND face_extrusion_tests
        "[face_extrusion],[face_extrusion_model],[mesh_quality],[planar_profile],[cut_and_stretch],[cut_and_stretch_placement],[point_alignment],[cut_plane_selection],[cut_plane_preview_bounds],[facet_picker_geometry],[native_cut_parts],[cut_guide_selection],[cut_guided_scope],[strict_fff_export],[strict_stl_writer]" ${CATCH_EXTRA_ARGS})
endif()
