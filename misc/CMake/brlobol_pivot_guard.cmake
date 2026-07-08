#                 B R L O B O L _ P I V O T _ G U A R D . C M A K E
# BRL-CAD
#
# Copyright (c) 2026 United States Government as represented by
# the U.S. Army Research Laboratory.
#
# Guardrails for the Obol-canonical drawing migration.
#
# Keep this file small.  Its job is to prevent known-retired boundaries from
# returning, not to encode historical implementation details token by token.

if(NOT DEFINED BRLCAD_SOURCE_DIR)
  message(FATAL_ERROR "BRLCAD_SOURCE_DIR is required")
endif()

if(NOT DEFINED BRLOBOL_PIVOT_GUARD_MODE)
  set(BRLOBOL_PIVOT_GUARD_MODE "transitional")
endif()

set(_brlobol_guard_extensions
  [[(CMakeLists\.txt|\.(c|cc|cpp|cxx|h|hh|hpp|cmake)$)]])

set(_brlobol_guard_active_boundary_patterns
  [[#[ \t]*include[ \t]*[<"]bsg]]
  [[#[ \t]*include[ \t]*[<"]dm]]
  [[(^|[^A-Za-z0-9_])bsg_view([^A-Za-z0-9_]|$)]]
  [[(^|[^A-Za-z0-9_])bsg_render_item([^A-Za-z0-9_]|$)]]
  [[(^|[^A-Za-z0-9_])bsg_backend_scene([^A-Za-z0-9_]|$)]]
  [[(^|[^A-Za-z0-9_])bsg_vlist([^A-Za-z0-9_]|$)]]
  [[(^|[^A-Za-z0-9_])bsg_vlblock([^A-Za-z0-9_]|$)]]
  [[(^|[^A-Za-z0-9_])BSG_[A-Za-z0-9_]*VLIST[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
  [[(^|[^A-Za-z0-9_])dm_draw_]]
  [[(^|[^A-Za-z0-9_])dm_plugins([^A-Za-z0-9_]|$)]]
  [[(^|[^A-Za-z0-9_])libbsg([^A-Za-z0-9_]|$)]]
  [[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]]
)

set(_brlobol_guard_retired_ged_draw_symbols
  [[(^|[^A-Za-z0-9_])ged_draw_shape_draft_[A-Za-z0-9_]*[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_create_evaluated_path_shape_ref[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_scene_ref_create_draft_pair[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_group_ref_append_scene_ref[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_scene_ref_append_source_owner_to_group[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_scene_ref_attach_source_state_record[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_source_data_free[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_scene_ref_apply_display_settings[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_obol_database_source_publish_primitive_wireframe_for_path_internal[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_obol_database_source_publish_annotation_record_from_internal_for_path[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_bsg_appearance_from_neutral[ \t\r\n]*\(]]
)

set(_brlobol_guard_retired_test_files
  src/libged/tests/draw/bsg_quad_stability.cpp
  src/libged/tests/draw/bsg_render_stability.cpp
  src/libged/tests/draw/mged_bsg.cpp
  src/libged/tests/draw/mged_shaded_mode_bsg.cpp
  src/libged/tests/draw/rtwizard_bsg.cpp
  src/libged/tests/draw/tcl_overlay_bsg.cpp
)

set(_brlobol_guard_retired_test_tokens
  [[ged_test_bsg_quad_stability]]
  [[ged_test_bsg_render_stability]]
  [[ged_test_mged_bsg]]
  [[ged_test_mged_shaded_mode_bsg]]
  [[ged_test_rtwizard_bsg]]
  [[ged_test_tcl_overlay_bsg]]
  [[bsg_quad_stability\.cpp]]
  [[bsg_render_stability\.cpp]]
  [[mged_bsg\.cpp]]
  [[mged_shaded_mode_bsg\.cpp]]
  [[rtwizard_bsg\.cpp]]
  [[tcl_overlay_bsg\.cpp]]
)

function(_brlobol_guard_fail _msg)
  set_property(GLOBAL APPEND PROPERTY BRLOBOL_PIVOT_GUARD_FAILURES "${_msg}")
endfunction()

function(_brlobol_guard_collect _outvar)
  set(_ret)
  foreach(_entry IN LISTS ARGN)
    if(IS_DIRECTORY "${BRLCAD_SOURCE_DIR}/${_entry}")
      file(GLOB_RECURSE _dir_files LIST_DIRECTORIES false
	"${BRLCAD_SOURCE_DIR}/${_entry}/*")
      list(APPEND _ret ${_dir_files})
    elseif(EXISTS "${BRLCAD_SOURCE_DIR}/${_entry}")
      list(APPEND _ret "${BRLCAD_SOURCE_DIR}/${_entry}")
    endif()
  endforeach()
  list(REMOVE_DUPLICATES _ret)
  set(${_outvar} "${_ret}" PARENT_SCOPE)
endfunction()

function(_brlobol_guard_read_rel _outvar _rel)
  set(_path "${BRLCAD_SOURCE_DIR}/${_rel}")
  if(NOT EXISTS "${_path}")
    _brlobol_guard_fail("${_rel} is required for Obol pivot guard checks")
    set(${_outvar} "" PARENT_SCOPE)
    return()
  endif()
  file(READ "${_path}" _contents)
  set(${_outvar} "${_contents}" PARENT_SCOPE)
endfunction()

function(_brlobol_guard_forbid_regexes _rel _contents)
  foreach(_pat IN LISTS ARGN)
    string(REGEX MATCH "${_pat}" _hit "${_contents}")
    if(_hit)
      _brlobol_guard_fail("${_rel}: ${_hit}")
    endif()
  endforeach()
endfunction()

function(_brlobol_guard_check_file _file)
  file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
  if("${_rel}" STREQUAL "misc/CMake/brlobol_pivot_guard.cmake")
    return()
  endif()
  if(NOT "${_rel}" MATCHES "${_brlobol_guard_extensions}")
    return()
  endif()

  file(READ "${_file}" _contents)
  _brlobol_guard_forbid_regexes("${_rel}" "${_contents}"
    ${_brlobol_guard_active_boundary_patterns})
endfunction()

function(_brlobol_guard_check_dependency_inventory)
  _brlobol_guard_read_rel(_inventory
    "doc/notes/obol_legacy_dependency_inventory.txt")
  _brlobol_guard_read_rel(_source_dirs "src/source_dirs.cmake")

  set(_expected_deps_libbrlobol [[libwdb;librt;libimgstream;libbg;libbu]])
  set(_expected_deps_libbv [[libbg;libbn;libbu]])
  set(_expected_deps_libnmg [[libbg;libbn;libbu]])
  set(_expected_deps_libbrep [[libbg;libbu]])
  set(_expected_deps_librt [[libbrep;libnmg;libbsg;libbv;libbg;libbn;libbu]])
  set(_expected_deps_libanalyze [[librt;libbg;libbn;libbu]])
  set(_expected_deps_libimgstream [[libicv;libbn;libbu]])
  set(_expected_deps_libdm [[libbsg;librt;libicv;libbn;libpkg;libbu]])
  set(_expected_deps_libged [[libbrlobol;libicv;libanalyze;libwdb;liboptical;libdm;libbsg;libbu]])
  set(_expected_deps_libqtcad [[libbrlobol;libged;libdm;libbg;libbn;libbu]])
  set(_expected_deps_libtclcad [[libged;libdm;libbn;libbu]])
  set(_expected_deps_qged [[libqtcad]])
  set(_expected_deps_mged [[libtclcad]])

  foreach(_target
      libbrlobol
      libbv
      libnmg
      libbrep
      librt
      libanalyze
      libimgstream
      libdm
      libged
      libqtcad
      libtclcad
      qged
      mged)
    set(_expected "${_expected_deps_${_target}}")
    string(REGEX MATCH "set_deps\\(${_target}[ \t]+\"([^\"]*)\"\\)"
      _row "${_source_dirs}")
    set(_actual "${CMAKE_MATCH_1}")
    if(NOT _row OR NOT "${_actual}" STREQUAL "${_expected}")
      string(REPLACE ";" "," _expected_msg "${_expected}")
      _brlobol_guard_fail(
	"src/source_dirs.cmake direct dependency row changed for ${_target}; expected ${_expected_msg}")
    endif()
    string(FIND "${_inventory}" "* ${_target}: ${_expected}" _inventory_idx)
    if(_inventory_idx EQUAL -1)
      string(REPLACE ";" "," _expected_msg "${_expected}")
      _brlobol_guard_fail(
	"doc/notes/obol_legacy_dependency_inventory.txt missing ${_target}: ${_expected_msg}")
    endif()
  endforeach()

  _brlobol_guard_read_rel(_qged_cmake "src/qged/CMakeLists.txt")
  string(REGEX MATCH [[dm_plugins]] _qged_dm_plugins "${_qged_cmake}")
  if(_qged_dm_plugins)
    _brlobol_guard_fail("qged reintroduced dm_plugins as a build prerequisite")
  endif()
endfunction()

function(_brlobol_guard_check_retired_tests)
  foreach(_rel IN LISTS _brlobol_guard_retired_test_files)
    if(EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _brlobol_guard_fail(
	"${_rel} was retired; keep migrated coverage in Obol/native GED draw tests")
    endif()
  endforeach()

  _brlobol_guard_read_rel(_cmake "src/libged/tests/draw/CMakeLists.txt")
  foreach(_token IN LISTS _brlobol_guard_retired_test_tokens)
    string(REGEX MATCH "${_token}" _hit "${_cmake}")
    if(_hit)
      _brlobol_guard_fail(
	"src/libged/tests/draw/CMakeLists.txt reintroduced retired BSG test token ${_token}")
    endif()
  endforeach()

  _brlobol_guard_read_rel(_librt_tests "src/librt/tests/CMakeLists.txt")
  foreach(_token
      "brlcad_addexec(rt_view_legacy_bsg"
      "brlcad_add_test(NAME rt_view_legacy_bsg")
    string(FIND "${_librt_tests}" "${_token}" _idx)
    if(NOT _idx EQUAL -1)
      _brlobol_guard_fail(
	"src/librt/tests/CMakeLists.txt reintroduced retired RT view legacy test ${_token}")
    endif()
  endforeach()
endfunction()

function(_brlobol_guard_check_public_headers)
  foreach(_rel include/ged.h include/ged/defines.h include/ged/view.h)
    _brlobol_guard_read_rel(_contents "${_rel}")
    _brlobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[#[ \t]*include[ \t]*[<"]bsg]]
      [[#[ \t]*include[ \t]*[<"]dm/fbserv\.h]]
      [[#[ \t]*include[ \t]*[<"]dm]]
      [[(^|[^A-Za-z0-9_])bsg_view([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])bsg_scene_ref([^A-Za-z0-9_]|$)]])
  endforeach()
endfunction()

function(_brlobol_guard_check_legacy_include_allowlist)
  set(_allowed
    src/libged/ged_view_state.cpp
    src/librt/view_legacy_bsg.c)

  _brlobol_guard_collect(_files
    include
    src/libged
    src/librt)
  foreach(_file IN LISTS _files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    list(FIND _allowed "${_rel}" _allowed_idx)
    if(NOT _allowed_idx EQUAL -1)
      continue()
    endif()
    if(NOT "${_rel}" MATCHES "${_brlobol_guard_extensions}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    string(REGEX MATCH [[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
      _hit "${_contents}")
    if(_hit)
      _brlobol_guard_fail(
	"${_rel} directly includes rt/view_legacy_bsg.h outside the retained adapter allowlist")
    endif()
  endforeach()
endfunction()

function(_brlobol_guard_check_retired_ged_draw_symbols)
  foreach(_rel
      src/libged/ged_draw_source.c
      src/libged/ged_draw_transactions.c
      src/libged/draw_obol.cpp
      include/ged/draw.h)
    if(EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      file(READ "${BRLCAD_SOURCE_DIR}/${_rel}" _contents)
      _brlobol_guard_forbid_regexes("${_rel}" "${_contents}"
	${_brlobol_guard_retired_ged_draw_symbols})
    endif()
  endforeach()
endfunction()

function(_brlobol_guard_collect_active_scan_files _outvar)
  if("${BRLOBOL_PIVOT_GUARD_MODE}" STREQUAL "strict")
    file(GLOB_RECURSE _files LIST_DIRECTORIES false
      "${BRLCAD_SOURCE_DIR}/include/*"
      "${BRLCAD_SOURCE_DIR}/src/*"
      "${BRLCAD_SOURCE_DIR}/misc/CMake/*")
  else()
    _brlobol_guard_collect(_files
      include/brlobol
      include/brlobol.h
      include/qtcad/QgObolExport.h
      include/qtcad/QgObolMeasure.h
      include/qtcad/QgObolPick.h
      include/qtcad/QgObolSnap.h
      include/qtcad/QgObolWindowHost.h
      src/libbrlobol
      src/libqtcad/QgObolDatabaseSync.cpp
      src/libqtcad/QgObolDatabaseSyncPrivate.h
      src/libqtcad/QgObolEditPreview.cpp
      src/libqtcad/QgObolEditPreviewPrivate.h
      src/libqtcad/QgObolExport.cpp
      src/libqtcad/QgObolMeasure.cpp
      src/libqtcad/QgObolOverlaySync.cpp
      src/libqtcad/QgObolPick.cpp
      src/libqtcad/QgObolSelectionSync.cpp
      src/libqtcad/QgObolSnap.cpp
      src/libqtcad/QgObolWindowHost.cpp)
  endif()
  set(${_outvar} "${_files}" PARENT_SCOPE)
endfunction()

_brlobol_guard_check_dependency_inventory()
_brlobol_guard_check_retired_tests()
_brlobol_guard_check_public_headers()
_brlobol_guard_check_legacy_include_allowlist()
_brlobol_guard_check_retired_ged_draw_symbols()

_brlobol_guard_collect_active_scan_files(_brlobol_guard_scan_files)
foreach(_file IN LISTS _brlobol_guard_scan_files)
  _brlobol_guard_check_file("${_file}")
endforeach()

get_property(_failures GLOBAL PROPERTY BRLOBOL_PIVOT_GUARD_FAILURES)
if(_failures)
  list(SORT _failures)
  string(REPLACE ";" "\n  " _failure_text "${_failures}")
  message(FATAL_ERROR
    "Obol pivot guard failed in ${BRLOBOL_PIVOT_GUARD_MODE} mode:\n"
    "  ${_failure_text}\n"
    "New Obol-canonical code must not depend on retired BSG/DM APIs.")
endif()

message(STATUS "Obol pivot guard passed in ${BRLOBOL_PIVOT_GUARD_MODE} mode")
