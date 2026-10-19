# Native wx event-queue regressions, without a printer, plugin or GUI windows.
# Keep the console-app lifetime separate from the geometry/Model test runner.
add_executable(job_ui_dispatch_tests
    ../libslic3r/libslic3r_tests.cpp
    test_job_ui_dispatch.cpp
)
target_link_libraries(job_ui_dispatch_tests PRIVATE test_common libslic3r libslic3r_gui)
# wxWidgets' include paths and platform definitions are directory-scoped under
# src/, not public usage requirements of the production GUI library.
get_directory_property(_job_ui_includes DIRECTORY "${PROJECT_SOURCE_DIR}/src" INCLUDE_DIRECTORIES)
get_directory_property(_job_ui_definitions DIRECTORY "${PROJECT_SOURCE_DIR}/src" COMPILE_DEFINITIONS)
get_directory_property(_job_ui_options DIRECTORY "${PROJECT_SOURCE_DIR}/src" COMPILE_OPTIONS)
target_include_directories(job_ui_dispatch_tests PRIVATE ${_job_ui_includes} "${PROJECT_SOURCE_DIR}/src/libslic3r")
target_compile_definitions(job_ui_dispatch_tests PRIVATE ${_job_ui_definitions})
target_compile_options(job_ui_dispatch_tests PRIVATE ${_job_ui_options})
set_target_properties(job_ui_dispatch_tests PROPERTIES FOLDER "tests" MACOSX_BUNDLE FALSE)
add_test(NAME job_ui_dispatch_tests COMMAND job_ui_dispatch_tests "[job_ui_dispatch]" ${CATCH_EXTRA_ARGS})
