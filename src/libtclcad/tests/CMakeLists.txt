if(BRLCAD_ENABLE_TCL)
  brlcad_addexec(tclcad_check_ged_cmds check_tclcad_cmds.cpp "libtclcad;libged" TEST)
  if(TARGET tclcad_check_ged_cmds)
    add_dependencies(tclcad_check_ged_cmds libtclcad ged_plugins)
  endif(TARGET tclcad_check_ged_cmds)
  brlcad_add_test(NAME tclcad_test_ged_cmds_exist COMMAND tclcad_check_ged_cmds)
  add_test(
    NAME tclcad_test_new_view_cleanup
    COMMAND $<TARGET_FILE:btclsh> ${CMAKE_CURRENT_SOURCE_DIR}/new_view_cleanup.tcl
            ${BRLCAD_SOURCE_DIR}/src/libged/tests/draw/moss.g
  )
  set_property(TEST tclcad_test_new_view_cleanup APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  if(TCLCAD_TKOBOL_HOST AND UNIX)
    find_program(XVFB_RUN_EXECUTABLE xvfb-run)
    if(XVFB_RUN_EXECUTABLE)
      brlcad_addexec(tclcad_tkobol_thread_affinity
        tkobol_thread_affinity.cpp "libtclcad;libbrlobol;Threads::Threads" TEST)
      if(TARGET tclcad_tkobol_thread_affinity)
        add_dependencies(tclcad_tkobol_thread_affinity libtclcad)
      endif(TARGET tclcad_tkobol_thread_affinity)
      add_test(
        NAME tclcad_test_tkobol_thread_affinity
        COMMAND ${XVFB_RUN_EXECUTABLE} -a $<TARGET_FILE:tclcad_tkobol_thread_affinity>
      )
      set_property(TEST tclcad_test_tkobol_thread_affinity APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
      add_test(
        NAME tclcad_test_tkobol_widget_smoke
        COMMAND ${XVFB_RUN_EXECUTABLE} -a $<TARGET_FILE:btclsh>
                ${CMAKE_CURRENT_SOURCE_DIR}/tkobol_widget_smoke.tcl
                ${BRLCAD_SOURCE_DIR}/src/libged/tests/draw/moss.g
                ${CMAKE_CURRENT_BINARY_DIR}/tkobol_widget_smoke.png
      )
      set_property(TEST tclcad_test_tkobol_widget_smoke APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
      add_test(
        NAME tclcad_test_archer_obol_smoke
        COMMAND ${XVFB_RUN_EXECUTABLE} -a $<TARGET_FILE:btclsh>
                ${CMAKE_CURRENT_SOURCE_DIR}/archer_obol_smoke.tcl
                ${BRLCAD_SOURCE_DIR}/src/libged/tests/draw/moss.g
      )
      set_property(TEST tclcad_test_archer_obol_smoke APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
    endif(XVFB_RUN_EXECUTABLE)
  endif(TCLCAD_TKOBOL_HOST AND UNIX)
endif(BRLCAD_ENABLE_TCL)

cmakefiles(
  CMakeLists.txt
  archer_obol_smoke.tcl
  check_tclcad_cmds.cpp
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
