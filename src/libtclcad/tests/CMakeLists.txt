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
endif(BRLCAD_ENABLE_TCL)

cmakefiles(
  CMakeLists.txt
  check_tclcad_cmds.cpp
  new_view_cleanup.tcl
  )

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
