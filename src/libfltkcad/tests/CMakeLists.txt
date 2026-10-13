if(TARGET libfltkcad AND BUILD_TESTING)
  include(FLTKCADTest)

  brlcad_addexec(fltkcad_uri_test
    "uri_test.cpp;../Uri.cpp" "libbu" TEST)
  target_compile_features(fltkcad_uri_test PRIVATE cxx_std_17)
  target_include_directories(fltkcad_uri_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
  add_test(NAME fltkcad_uri_test COMMAND fltkcad_uri_test)

  brlcad_addexec(fltkcad_html_test "html_smoke.cpp" "libfltkcad" TEST)
  fltkcad_add_gui_test(fltkcad_html_test $<TARGET_FILE:fltkcad_html_test>)

  brlcad_addexec(fltkcad_help_test "help_test.cpp" "libfltkcad" TEST)
  fltkcad_add_gui_test(fltkcad_help_test $<TARGET_FILE:fltkcad_help_test>
    "${CMAKE_CURRENT_SOURCE_DIR}/fixtures/integration.html"
    "${BRLCAD_SOURCE_DIR}/doc/html/manuals/mged/ged.html")

  brlcad_addexec(fltkcad_history_test "terminal_history_test.cpp" "libfltkcad" TEST)
  add_test(NAME fltkcad_history_test COMMAND fltkcad_history_test)

  brlcad_addexec(fltkcad_terminal_test "terminal_widget_test.cpp" "libfltkcad" TEST)
  fltkcad_add_gui_test(fltkcad_terminal_test $<TARGET_FILE:fltkcad_terminal_test>)

  brlcad_addexec(fltkcad_dm_test "dm_backend_test.cpp" "libfltkcad" TEST)
  fltkcad_add_gui_test(fltkcad_dm_osmesa_test
    $<TARGET_FILE:fltkcad_dm_test> osmesa)
  fltkcad_add_gui_test(fltkcad_dm_opengl_test
    $<TARGET_FILE:fltkcad_dm_test> opengl)

  if(EXISTS "${FLTKCAD_TEST_DATABASE}")
    brlcad_addexec(fltkcad_widget_suite "widget_suite.cpp" "libfltkcad" TEST)
    fltkcad_add_gui_test(fltkcad_widgets_software
      $<TARGET_FILE:fltkcad_widget_suite> software
      "${FLTKCAD_TEST_DATABASE}")
    fltkcad_add_gui_test(fltkcad_widgets_opengl
      $<TARGET_FILE:fltkcad_widget_suite> opengl
      "${FLTKCAD_TEST_DATABASE}")
  else()
    message(STATUS
      "FLTKCAD_TEST_DATABASE does not exist; database widget tests are disabled")
  endif()

  if(BRLCAD_ENABLE_TCL)
    brlcad_addexec(fltkcad_tcl_events_test "tcl_events_test.cpp" "libfltkcad" TEST)
    fltkcad_add_gui_test(fltkcad_tcl_events $<TARGET_FILE:fltkcad_tcl_events_test>)
  endif()

  if(BRLCAD_ENABLE_TCL AND BRLCAD_ENABLE_TK)
    brlcad_addexec(fltkcad_tcl_tk_test "tcl_tk_test.cpp" "libfltkcad" TEST)
    fltkcad_add_gui_test(fltkcad_tcl_tk $<TARGET_FILE:fltkcad_tcl_tk_test>)
  endif()

  if(BRLCAD_X11_PROVIDER_RESOLVED STREQUAL "XMIN"
     AND FLTKCAD_XMIN_RUN AND FLTKCAD_XMINCTL
     AND TARGET fltkcad_html_demo)
    add_test(NAME fltkcad_xmin_integration
      COMMAND sh "${CMAKE_CURRENT_SOURCE_DIR}/xmin_integration.sh"
        "${FLTKCAD_XMIN_RUN}"
        "${FLTKCAD_XMINCTL}"
        $<TARGET_FILE:fltkcad_html_demo>
        "${CMAKE_CURRENT_SOURCE_DIR}/fixtures/integration.html"
        "${CMAKE_CURRENT_BINARY_DIR}"
    )
    set_tests_properties(fltkcad_xmin_integration PROPERTIES TIMEOUT 20)

    add_test(NAME fltkcad_xmin_dependency_audit
      COMMAND sh "${CMAKE_CURRENT_SOURCE_DIR}/xmin_dependency_audit.sh"
        $<TARGET_FILE:fltkcad_html_demo>
    )
    set_tests_properties(fltkcad_xmin_dependency_audit PROPERTIES
      SKIP_RETURN_CODE 77)
  endif()
endif()

cmakefiles(CMakeLists.txt help_test.cpp uri_test.cpp html_smoke.cpp html_test_utils.h terminal_history_test.cpp
  terminal_widget_test.cpp dm_backend_test.cpp widget_suite.cpp tcl_tk_test.cpp tcl_events_test.cpp
  xmin_integration.sh xmin_dependency_audit.sh fixtures/integration.html
  fixtures/integration.css fixtures/integration.svg
  package/CMakeLists.txt package/main.cpp)
