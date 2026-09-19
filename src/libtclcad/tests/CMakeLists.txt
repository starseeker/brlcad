if(BRLCAD_ENABLE_TCL)
  brlcad_addexec(tclcad_check_ged_cmds check_tclcad_cmds.cpp "libtclcad;libged" TEST)
  if(TARGET tclcad_check_ged_cmds)
    add_dependencies(tclcad_check_ged_cmds libtclcad ged_plugins)
  endif(TARGET tclcad_check_ged_cmds)
  brlcad_add_test(NAME tclcad_test_ged_cmds_exist COMMAND tclcad_check_ged_cmds)
  brlcad_addexec(tclcad_check_itcl check_itcl.cpp "libtclcad" TEST)
  brlcad_addexec(tclcad_check_multi_interp check_multi_interp.cpp "libtclcad" TEST)
  brlcad_add_test(NAME tclcad_test_itcl_version COMMAND tclcad_check_itcl)
  brlcad_add_test(NAME tclcad_test_multi_interp COMMAND tclcad_check_multi_interp)

  if(BRLCAD_ENABLE_TK AND
     BRLCAD_X11_PROVIDER_RESOLVED STREQUAL "XMIN" AND
     Xmin_RUN_EXECUTABLE)
    brlcad_add_gui_test(
      NAME tclcad_test_multi_interp_gui
      SCREEN 1280x960x24
      TIMEOUT 120
      COMMAND $<TARGET_FILE:tclcad_check_multi_interp> --gui
    )
  endif()
  add_test(
    NAME tclcad_test_new_view_cleanup
    COMMAND $<TARGET_FILE:btclsh> ${CMAKE_CURRENT_SOURCE_DIR}/new_view_cleanup.tcl
            ${BRLCAD_SOURCE_DIR}/src/libged/tests/draw/moss.g
  )
  set_property(TEST tclcad_test_new_view_cleanup APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  if(TCLCAD_TKOBOL_HOST AND UNIX)
    find_program(XVFB_RUN_EXECUTABLE xvfb-run)
    if(XVFB_RUN_EXECUTABLE)
      set(TCLCAD_XVFB_TEST_WRAPPER
        ${BRLCAD_SOURCE_DIR}/misc/CMake/scripts/xvfb_test.sh)
      brlcad_addexec(tclcad_tkobol_thread_affinity
        tkobol_thread_affinity.cpp "libtclcad;libBObol;Threads::Threads" TEST)
      if(TARGET tclcad_tkobol_thread_affinity)
        add_dependencies(tclcad_tkobol_thread_affinity libtclcad)
      endif(TARGET tclcad_tkobol_thread_affinity)
      add_test(
        NAME tclcad_test_tkobol_thread_affinity
        COMMAND ${SH_EXEC} ${TCLCAD_XVFB_TEST_WRAPPER} ${XVFB_RUN_EXECUTABLE}
                $<TARGET_FILE:tclcad_tkobol_thread_affinity>
      )
      set_property(TEST tclcad_test_tkobol_thread_affinity APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
      set_property(TEST tclcad_test_tkobol_thread_affinity PROPERTY SKIP_RETURN_CODE 125)
      add_test(
        NAME tclcad_test_tkobol_widget_smoke
        COMMAND ${SH_EXEC} ${TCLCAD_XVFB_TEST_WRAPPER} ${XVFB_RUN_EXECUTABLE}
                $<TARGET_FILE:btclsh>
                ${CMAKE_CURRENT_SOURCE_DIR}/tkobol_widget_smoke.tcl
                ${BRLCAD_SOURCE_DIR}/src/libged/tests/draw/moss.g
                ${CMAKE_CURRENT_BINARY_DIR}/tkobol_widget_smoke.png
      )
      set_property(TEST tclcad_test_tkobol_widget_smoke APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
      set_property(TEST tclcad_test_tkobol_widget_smoke PROPERTY SKIP_RETURN_CODE 125)
      add_test(
        NAME tclcad_test_archer_obol_smoke
        COMMAND ${SH_EXEC} ${TCLCAD_XVFB_TEST_WRAPPER} ${XVFB_RUN_EXECUTABLE}
                $<TARGET_FILE:btclsh>
                ${CMAKE_CURRENT_SOURCE_DIR}/archer_obol_smoke.tcl
                ${BRLCAD_SOURCE_DIR}/src/libged/tests/draw/moss.g
                ${CMAKE_CURRENT_BINARY_DIR}/archer_obol_smoke.g
      )
      set_property(TEST tclcad_test_archer_obol_smoke APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
      set_property(TEST tclcad_test_archer_obol_smoke PROPERTY SKIP_RETURN_CODE 125)
    endif(XVFB_RUN_EXECUTABLE)
  endif(TCLCAD_TKOBOL_HOST AND UNIX)
endif(BRLCAD_ENABLE_TCL)

cmakefiles(
  CMakeLists.txt
  archer_obol_smoke.tcl
  check_tclcad_cmds.cpp
  check_itcl.cpp
  check_multi_interp.cpp
  itcl4_gui.tcl
  new_view_cleanup.tcl
  tkobol_thread_affinity.cpp
  tkobol_widget_smoke.tcl
  )

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
