if(TARGET fltkged AND BUILD_TESTING)
  include(FLTKCADTest)
  if(BRLCAD_ENABLE_TCL)
    add_test(NAME fltkged_batch COMMAND ${CMAKE_COMMAND}
      -DFLTKGED=$<TARGET_FILE:fltkged> -DWORK=${CMAKE_CURRENT_BINARY_DIR}/batch
      -P ${CMAKE_CURRENT_SOURCE_DIR}/batch.cmake)
    set_tests_properties(fltkged_batch PROPERTIES TIMEOUT 30 ENVIRONMENT "DISPLAY=")
  endif()

  brlcad_addexec(fltkged_layout_test
    "layout_test.cpp;../FltkGedWindow.cpp;../Fl_Ged_Localization.cpp" "libfltkcad" TEST)
  target_compile_features(fltkged_layout_test PRIVATE cxx_std_17)
  target_include_directories(fltkged_layout_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
  fltkcad_add_gui_test(fltkged_layout $<TARGET_FILE:fltkged_layout_test>)
  fltkcad_add_gui_test(fltkged_combined_smoke
    $<TARGET_FILE:fltkged> --combined --smoke --software)

  if(EXISTS "${FLTKCAD_TEST_DATABASE}")
    brlcad_addexec(fltkged_application_test
      "application_test.cpp;../FltkGedWindow.cpp;../Fl_Ged_Localization.cpp" "libfltkcad" TEST)
    target_compile_features(fltkged_application_test PRIVATE cxx_std_17)
    target_include_directories(fltkged_application_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
    fltkcad_add_gui_test(fltkged_application
      $<TARGET_FILE:fltkged_application_test> "${FLTKCAD_TEST_DATABASE}")

    fltkcad_add_gui_test(fltkged_software_smoke
      $<TARGET_FILE:fltkged> --smoke --software
      "${FLTKCAD_TEST_DATABASE}")
    fltkcad_add_gui_test(fltkged_opengl_smoke
      $<TARGET_FILE:fltkged> --smoke --opengl
      "${FLTKCAD_TEST_DATABASE}")

    if(TARGET rt)
      fltkcad_add_gui_test(fltkged_render_software
        $<TARGET_FILE:fltkged> --render-smoke --software
        "${FLTKCAD_TEST_DATABASE}")
      set_tests_properties(fltkged_render_software PROPERTIES
        ENVIRONMENT "FLTKGED_FORCE_PROCESS_POLL=1")
      fltkcad_add_gui_test(fltkged_render_opengl
        $<TARGET_FILE:fltkged> --render-smoke --opengl
        "${FLTKCAD_TEST_DATABASE}")
    endif()
  else()
    message(STATUS
      "FLTKCAD_TEST_DATABASE does not exist; fltkged database tests are disabled")
  endif()
  fltkcad_add_gui_test(fltkged_missing_database
    $<TARGET_FILE:fltkged> --smoke --software
    "${CMAKE_CURRENT_BINARY_DIR}/missing.g")
  set_tests_properties(fltkged_missing_database PROPERTIES WILL_FAIL TRUE)
endif()

cmakefiles(CMakeLists.txt layout_test.cpp application_test.cpp batch.cmake regression_audit.py)
