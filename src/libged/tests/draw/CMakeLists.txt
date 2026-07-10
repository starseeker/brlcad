set(
  test_deps
  dm-obol
  ged-close
  ged-draw
  ged-dm
  ged-brep
  ged-facetize
  ged-open
  ged-select
  ged-tol
  ged-view
  ged-zap
)

set(
  obol_test_deps
  ged-close
  ged-draw
  ged-brep
  ged-facetize
  ged-open
  ged-select
  ged-tol
  ged-view
  ged-zap
)

set(basic_view_test_ctrls basic.apng)
brlcad_addexec(ged_test_draw "basic.cpp;util.cpp" "libged;libbrlobol;PNG::PNG;ZLIB::ZLIB;Threads::Threads" TEST)
if(TARGET ged_test_draw)
  foreach(tdp ${obol_test_deps})
    if (TARGET ${tdp})
      add_dependencies(ged_test_draw ${tdp})
    endif ()
  endforeach()
endif(TARGET ged_test_draw)
brlcad_add_test(NAME ged_test_drawing_basic COMMAND ged_test_draw -c "${CMAKE_CURRENT_SOURCE_DIR}")
if(TARGET ged_test_draw)
  set_property(TEST ged_test_drawing_basic APPEND PROPERTY ENVIRONMENT GED_TEST_DRAW_OBOL_CAPTURE=1)
endif(TARGET ged_test_draw)

set(faceplate_test_ctrls faceplate.apng)
brlcad_addexec(ged_test_faceplate "faceplate.cpp;util.cpp" "libged;libbrlobol;PNG::PNG;ZLIB::ZLIB;Threads::Threads" TEST)
if(TARGET ged_test_faceplate)
  foreach(tdp ${obol_test_deps})
    if (TARGET ${tdp})
      add_dependencies(ged_test_faceplate ${tdp})
    endif ()
  endforeach()
endif(TARGET ged_test_faceplate)
brlcad_add_test(NAME ged_test_drawing_faceplate COMMAND ged_test_faceplate -c "${CMAKE_CURRENT_SOURCE_DIR}")
if(TARGET ged_test_faceplate)
  set_property(TEST ged_test_drawing_faceplate APPEND PROPERTY ENVIRONMENT GED_TEST_DRAW_OBOL_CAPTURE=1)
endif(TARGET ged_test_faceplate)

set(lod_test_ctrls lod.apng)
brlcad_addexec(ged_test_lod "lod.cpp;util.cpp" "libged;libbrlobol;PNG::PNG;ZLIB::ZLIB;Threads::Threads" TEST)
if(TARGET ged_test_lod)
  foreach(tdp ${obol_test_deps})
    if (TARGET ${tdp})
      add_dependencies(ged_test_lod ${tdp})
    endif ()
  endforeach()
endif(TARGET ged_test_lod)
brlcad_add_test(NAME ged_test_drawing_lod COMMAND ged_test_lod -c "${CMAKE_CURRENT_SOURCE_DIR}")
if(TARGET ged_test_lod)
  set_property(TEST ged_test_drawing_lod APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST ged_test_drawing_lod APPEND PROPERTY ENVIRONMENT GED_TEST_DRAW_OBOL_CAPTURE=1)
endif(TARGET ged_test_lod)
distclean(${CMAKE_CURRENT_BINARY_DIR}/moss_lod_tmp.g)

set(select_test_ctrls select.apng)
brlcad_addexec(ged_test_select_draw "select.cpp;util.cpp" "libged;libbrlobol;PNG::PNG;ZLIB::ZLIB;Threads::Threads" TEST)
if(TARGET ged_test_select_draw)
  foreach(tdp ${obol_test_deps})
    if (TARGET ${tdp})
      add_dependencies(ged_test_select_draw ${tdp})
    endif ()
  endforeach()
endif(TARGET ged_test_select_draw)
brlcad_add_test(NAME ged_test_drawing_select COMMAND ged_test_select_draw -c "${CMAKE_CURRENT_SOURCE_DIR}")
if(TARGET ged_test_select_draw)
  set_property(TEST ged_test_drawing_select APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST ged_test_drawing_select APPEND PROPERTY ENVIRONMENT GED_TEST_DRAW_OBOL_CAPTURE=1)
endif(TARGET ged_test_select_draw)
distclean(${CMAKE_CURRENT_BINARY_DIR}/moss_select_tmp.g)

set(quad_test_ctrls quad_00.apng quad_01.apng quad_02.apng quad_03.apng)
brlcad_addexec(ged_test_quad "quad.cpp;util.cpp" "libged;libbrlobol;PNG::PNG;ZLIB::ZLIB;Threads::Threads" TEST)
if(TARGET ged_test_quad)
  foreach(tdp ${obol_test_deps})
    if (TARGET ${tdp})
      add_dependencies(ged_test_quad ${tdp})
    endif ()
  endforeach()
endif(TARGET ged_test_quad)
brlcad_add_test(NAME ged_test_drawing_quad COMMAND ged_test_quad -c "${CMAKE_CURRENT_SOURCE_DIR}")
if(TARGET ged_test_quad)
  set_property(TEST ged_test_drawing_quad APPEND PROPERTY ENVIRONMENT GED_TEST_DRAW_OBOL_CAPTURE=1)
endif(TARGET ged_test_quad)
distclean(${CMAKE_CURRENT_BINARY_DIR}/moss_quad_tmp.g)

# The idea here is to test view stability across resize events
set(aet_test_ctrls aet_00.apng aet_01.apng aet_02.apng aet_03.apng)
brlcad_addexec(ged_test_aet "aet.cpp;util.cpp" "libged;libdm;libbrlobol;PNG::PNG;ZLIB::ZLIB;Threads::Threads" TEST)
if(TARGET ged_test_aet)
  foreach(tdp ${test_deps})
    if (TARGET ${tdp})
      add_dependencies(ged_test_aet ${tdp})
    endif ()
  endforeach()
endif(TARGET ged_test_aet)

# Typed selection, snap, and measure semantic records
brlcad_addexec(ged_test_selection_semantics "selection_semantics.cpp" "libged;libbu" TEST)
if(TARGET ged_test_selection_semantics)
  foreach(tdp ${test_deps})
    if (TARGET ${tdp})
      add_dependencies(ged_test_selection_semantics ${tdp})
    endif ()
  endforeach()
endif(TARGET ged_test_selection_semantics)
brlcad_add_test(NAME ged_test_selection_semantics
  COMMAND ged_test_selection_semantics "${CMAKE_CURRENT_SOURCE_DIR}")
if(TARGET ged_test_selection_semantics)
  set_property(TEST ged_test_selection_semantics
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
endif()
distclean(${CMAKE_CURRENT_BINARY_DIR}/moss_selection_semantics_tmp.g)
distclean(${CMAKE_CURRENT_BINARY_DIR}/ged_selection_semantics_cache)

brlcad_addexec(ged_test_snap_semantics "snap_semantics.cpp" "libbv;libbu" TEST)
if(TARGET ged_test_snap_semantics)
  foreach(tdp ${test_deps})
    if (TARGET ${tdp})
      add_dependencies(ged_test_snap_semantics ${tdp})
    endif ()
  endforeach()
endif(TARGET ged_test_snap_semantics)
brlcad_add_test(NAME ged_test_snap_semantics COMMAND ged_test_snap_semantics)
if(TARGET ged_test_snap_semantics)
  set_property(TEST ged_test_snap_semantics
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
endif()

# LoD cross-run cache stability (populate cache → close → reopen → same pixels)
brlcad_addexec(ged_test_lod_crossrun "lod_crossrun.cpp;util.cpp" "libged;libdm;libbrlobol" TEST)
if(TARGET ged_test_lod_crossrun)
  foreach(tdp ${test_deps})
    if (TARGET ${tdp})
      add_dependencies(ged_test_lod_crossrun ${tdp})
    endif ()
  endforeach()
endif(TARGET ged_test_lod_crossrun)
brlcad_add_test(NAME ged_test_lod_crossrun
  COMMAND ged_test_lod_crossrun -c "${CMAKE_CURRENT_SOURCE_DIR}")
if(TARGET ged_test_lod_crossrun)
  set_property(TEST ged_test_lod_crossrun
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
endif()
distclean(${CMAKE_CURRENT_BINARY_DIR}/moss_lod_cr_tmp.g)
distclean(${CMAKE_CURRENT_BINARY_DIR}/lod_crossrun_cache)

# MGED view-state regression tests not covered by image baselines.
brlcad_addexec(ged_test_mged_view_state "mged_view_state.cpp;util.cpp" "libged;libdm;libbrlobol" TEST)
if(TARGET ged_test_mged_view_state)
  foreach(tdp ${test_deps})
    if (TARGET ${tdp})
      add_dependencies(ged_test_mged_view_state ${tdp})
    endif ()
  endforeach()
endif(TARGET ged_test_mged_view_state)
brlcad_add_test(NAME ged_test_mged_view_state COMMAND ged_test_mged_view_state -c "${CMAKE_CURRENT_SOURCE_DIR}")
if(TARGET ged_test_mged_view_state)
  set_property(TEST ged_test_mged_view_state APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
endif()
distclean(${CMAKE_CURRENT_BINARY_DIR}/mged_view_state_t1.g)
distclean(${CMAKE_CURRENT_BINARY_DIR}/mged_view_state_t2.g)
distclean(${CMAKE_CURRENT_BINARY_DIR}/mged_view_state_t3.g)

# rtwizard view pipeline regression tests.
brlcad_addexec(ged_test_rtwizard_view "rtwizard_view.cpp;util.cpp" "libged;libdm;libbrlobol;libicv" TEST)
if(TARGET ged_test_rtwizard_view)
  foreach(tdp ${test_deps})
    if (TARGET ${tdp})
      add_dependencies(ged_test_rtwizard_view ${tdp})
    endif ()
  endforeach()
endif(TARGET ged_test_rtwizard_view)
brlcad_add_test(NAME ged_test_rtwizard_view COMMAND ged_test_rtwizard_view -c "${CMAKE_CURRENT_SOURCE_DIR}")
if(TARGET ged_test_rtwizard_view)
  set_property(TEST ged_test_rtwizard_view APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
endif()
distclean(${CMAKE_CURRENT_BINARY_DIR}/rtw_view_t1.g)
distclean(${CMAKE_CURRENT_BINARY_DIR}/rtw_view_t2.g)
distclean(${CMAKE_CURRENT_BINARY_DIR}/rtw_view_t3.g)
distclean(${CMAKE_CURRENT_BINARY_DIR}/rtw_view_t4.g)
distclean(${CMAKE_CURRENT_BINARY_DIR}/rtw_view_t5.g)

# GED draw transaction to Obol controller bridge regression.
brlcad_addexec(ged_test_obol_draw_sync "obol_sync.cpp" "libged;libbrlobol;libwdb" TEST)
if(TARGET ged_test_obol_draw_sync)
  foreach(tdp ${test_deps})
    if (TARGET ${tdp})
      add_dependencies(ged_test_obol_draw_sync ${tdp})
    endif ()
  endforeach()
endif(TARGET ged_test_obol_draw_sync)
brlcad_add_test(NAME ged_test_obol_draw_sync COMMAND ged_test_obol_draw_sync)
if(TARGET ged_test_obol_draw_sync)
  set_property(TEST ged_test_obol_draw_sync APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
endif()
distclean(${CMAKE_CURRENT_BINARY_DIR}/ged_obol_draw_sync_tmp.g)

# Independent-view state is structural scope, not a mode bit.
brlcad_addexec(ged_test_view_independent "view_independent.cpp" "libged;libbu" TEST)
if(TARGET ged_test_view_independent)
  foreach(tdp ${test_deps})
    if (TARGET ${tdp})
      add_dependencies(ged_test_view_independent ${tdp})
    endif ()
  endforeach()
endif(TARGET ged_test_view_independent)
brlcad_add_test(NAME ged_test_view_independent
  COMMAND ged_test_view_independent "${CMAKE_CURRENT_SOURCE_DIR}")
if(TARGET ged_test_view_independent)
  set_property(TEST ged_test_view_independent APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
endif()

# Phase V7: consolidated view command grammar and aliases
brlcad_addexec(ged_test_view_command "view_command.cpp" "libged;libbu" TEST)
if(TARGET ged_test_view_command)
  foreach(tdp ${test_deps})
    if (TARGET ${tdp})
      add_dependencies(ged_test_view_command ${tdp})
    endif ()
  endforeach()
endif(TARGET ged_test_view_command)
brlcad_add_test(NAME ged_test_view_command
  COMMAND ged_test_view_command "${CMAKE_CURRENT_SOURCE_DIR}")
if(TARGET ged_test_view_command)
  set_property(TEST ged_test_view_command APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
endif()

# Lighting/zbuffer/zclip regression test (mged_shaded_mode_helper path)
brlcad_addexec(ged_test_mged_shaded_mode "mged_shaded_mode.cpp;util.cpp" "libged;libdm;libbrlobol;libicv" TEST)
if(TARGET ged_test_mged_shaded_mode)
  foreach(tdp ${test_deps})
    if (TARGET ${tdp})
      add_dependencies(ged_test_mged_shaded_mode ${tdp})
    endif ()
  endforeach()
endif(TARGET ged_test_mged_shaded_mode)
brlcad_add_test(NAME ged_test_mged_shaded_mode
  COMMAND ged_test_mged_shaded_mode "${CMAKE_CURRENT_SOURCE_DIR}")
if(TARGET ged_test_mged_shaded_mode)
  set_property(TEST ged_test_mged_shaded_mode
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
endif()
distclean(${CMAKE_CURRENT_BINARY_DIR}/smb_t1.g)
distclean(${CMAKE_CURRENT_BINARY_DIR}/smb_t2.g)
distclean(${CMAKE_CURRENT_BINARY_DIR}/smb_t3.g)
distclean(${CMAKE_CURRENT_BINARY_DIR}/smb_t4.g)

foreach(_ged_drawing_stage0_test
    ged_test_drawing_basic
    ged_test_drawing_faceplate
    ged_test_drawing_lod
    ged_test_drawing_select
    ged_test_drawing_quad
    ged_test_selection_semantics
    ged_test_snap_semantics
    ged_test_lod_crossrun
    ged_test_mged_view_state
    ged_test_rtwizard_view
    ged_test_view_independent
    ged_test_view_command
    ged_test_mged_shaded_mode)
  if(TEST ${_ged_drawing_stage0_test})
    set_property(TEST ${_ged_drawing_stage0_test} APPEND PROPERTY LABELS drawing_baseline)
  endif()
endforeach()

cmakefiles(
  README
  empty.png
  moss.g
  moss.png
  rook.g
  selection_semantics.cpp
  snap_semantics.cpp
  view_independent.cpp
  view_command.cpp
  lod_crossrun.cpp
  mged_view_state.cpp
  mged_shaded_mode.cpp
  rtwizard_view.cpp
  util.cpp
  scripts/verify_main_parity.sh
  ${basic_view_test_ctrls}
  ${faceplate_test_ctrls}
  ${lod_test_ctrls}
  ${select_test_ctrls}
  ${quad_test_ctrls}
  ${aet_test_ctrls}
)

cmakefiles(CMakeLists.txt)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
