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
  _brlobol_guard_read_rel(_libdm_cmake "src/libdm/CMakeLists.txt")

  set(_expected_deps_libbrlobol [[libwdb;librt;libimgstream;libbg;libbu]])
  set(_expected_deps_libbv [[libbg;libbn;libbu]])
  set(_expected_deps_libnmg [[libbg;libbn;libbu]])
  set(_expected_deps_libbrep [[libbg;libbu]])
  set(_expected_deps_librt [[libbrep;libnmg;libbv;libbg;libbn;libbu]])
  set(_expected_deps_libanalyze [[librt;libbg;libbn;libbu]])
  set(_expected_deps_libimgstream [[libicv;libbn;libpkg;libbu]])
  set(_expected_deps_libdm [[librt;libicv;libbn;libpkg;libbu]])
  set(_expected_deps_libged [[libbrlobol;libimgstream;libicv;libanalyze;libwdb;liboptical;libbu]])
  set(_expected_deps_libqtcad [[libbrlobol;libged;libbg;libbn;libbu]])
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

  string(REGEX MATCH [[PRIVATE_LIBS[^\n]*libimgstream]] _libdm_imgstream_private
    "${_libdm_cmake}")
  if(NOT _libdm_imgstream_private)
    _brlobol_guard_fail(
      "src/libdm/CMakeLists.txt no longer links libimgstream privately for display-host framebuffer compatibility")
  endif()
  string(FIND "${_inventory}" [[libdm retains only the narrow `struct fb` fbserv backend adapter]] _libdm_adapter_inventory_idx)
  if(_libdm_adapter_inventory_idx EQUAL -1)
    _brlobol_guard_fail(
      "doc/notes/obol_legacy_dependency_inventory.txt no longer documents the narrow libdm fbserv compatibility adapter")
  endif()
endfunction()

function(_brlobol_guard_check_libbsg_retired)
  foreach(_rel include/bsg/CMakeLists.txt include/bsg.h src/libbsg/CMakeLists.txt)
    if(EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _brlobol_guard_fail("${_rel} reintroduced the retired libbsg API or implementation")
    endif()
  endforeach()

  _brlobol_guard_collect(_files include src misc/CMake)
  foreach(_file IN LISTS _files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if("${_rel}" STREQUAL "misc/CMake/brlobol_pivot_guard.cmake" OR
       NOT "${_rel}" MATCHES "${_brlobol_guard_extensions}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    _brlobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[#[ \t]*include[ \t]*[<"]bsg]]
      [[(^|[^A-Za-z0-9_])libbsg([^A-Za-z0-9_]|$)]])
  endforeach()
endfunction()

function(_brlobol_guard_check_legacy_graphical_dm_retired)
  foreach(_rel
      src/libdm/dm-gl.c
      src/libdm/dm-gl.h
      src/libdm/X/CMakeLists.txt
      src/libdm/glx/CMakeLists.txt
      src/libdm/qtgl/CMakeLists.txt
      src/libdm/swrast/CMakeLists.txt
      src/libdm/tkswrast/CMakeLists.txt
      src/libdm/wgl/CMakeLists.txt
      src/libdm/fb_imgstream.c
      src/libdm/tests/fb_imgstream.cpp)
    if(EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _brlobol_guard_fail("${_rel} reintroduced a retired graphical display-manager path")
    endif()
  endforeach()

  _brlobol_guard_read_rel(_libdm_cmake "src/libdm/CMakeLists.txt")
  _brlobol_guard_forbid_regexes("src/libdm/CMakeLists.txt" "${_libdm_cmake}"
    [[add_subdirectory\((X|glx|qtgl|swrast|tkswrast|wgl)\)]]
    [[(^|[^A-Za-z0-9_])libdmgl([^A-Za-z0-9_]|$)]])
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

function(_brlobol_guard_check_removed_rt_view_aliases)
  _brlobol_guard_read_rel(_rt_view "include/rt/view.h")
  foreach(_token
      rt_view_info
      rt_view_lod_policy
      rt_view_lod_settings
      RT_VIEW_INFO_INIT
      RT_VIEW_LOD_POLICY_INIT
      RT_VIEW_LOD_SETTINGS_INIT)
    string(FIND "${_rt_view}" "${_token}" _idx)
    if(NOT _idx EQUAL -1)
      _brlobol_guard_fail(
	"include/rt/view.h reintroduced removed libbv compatibility alias ${_token}")
    endif()
  endforeach()
  foreach(_rel src/librt/view.c src/librt/tests/view_info.c)
    if(EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _brlobol_guard_fail("${_rel} reintroduced removed libbv wrapper code")
    endif()
  endforeach()
endfunction()

function(_brlobol_guard_check_native_drawing_tests)
  foreach(_rel
      src/libged/tests/draw/basic.cpp
      src/libged/tests/draw/faceplate.cpp
      src/libged/tests/draw/lod.cpp
      src/libged/tests/draw/select.cpp
      src/libged/tests/draw/quad.cpp
      src/libged/tests/draw/util.cpp)
    _brlobol_guard_read_rel(_contents "${_rel}")
    _brlobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[#[ \t]*include[ \t]*[<"]dm]]
      [[(^|[^A-Za-z0-9_])DM_SWRAST([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])ged_exec_dm[ \t\r\n]*\(]]
      [[(^|[^A-Za-z0-9_])ged_view_context_display_manager_get[ \t\r\n]*\(]])
  endforeach()

  _brlobol_guard_read_rel(_cmake "src/libged/tests/draw/CMakeLists.txt")
  foreach(_target
      ged_test_draw
      ged_test_faceplate
      ged_test_lod
      ged_test_select_draw
      ged_test_quad)
    string(REGEX MATCH
      "brlcad_addexec\\(${_target}[ \t][^\n]*libdm"
      _legacy_link "${_cmake}")
    if(_legacy_link)
      _brlobol_guard_fail(
	"src/libged/tests/draw/CMakeLists.txt links native Obol test ${_target} to libdm")
    endif()
  endforeach()
endfunction()

function(_brlobol_guard_check_obol_controller_quarantine)
  set(_allowed
    include/dm/obol.h
    src/libdm/dm-generic.c
    src/libged/dm/dm.c
    src/mged/attach.c)

  _brlobol_guard_collect(_files
    include
    src/libbrlobol
    src/libdm
    src/libged
    src/libqtcad
    src/qged
    src/gtools/gsh
    src/mged
    src/libtclcad)
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
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])dm_obol_controller[ \t\r\n]*\(]]
      _hit "${_contents}")
    if(_hit)
      _brlobol_guard_fail(
	"${_rel} calls dm_obol_controller outside the libdm/GED Obol bridge quarantine")
    endif()
  endforeach()
endfunction()

function(_brlobol_guard_check_qtcad_legacy_dm_open_quarantine)
  _brlobol_guard_collect(_files
    include/qtcad
    src/libdm
    src/libqtcad
    src/qged)
  foreach(_file IN LISTS _files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if(NOT "${_rel}" MATCHES "${_brlobol_guard_extensions}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])qg_legacy_view_dm_open_(qtgl|swrast)[ \t\r\n]*\(]]
      _hit "${_contents}")
    if(_hit)
      _brlobol_guard_fail(
	"${_rel} opens a qtcad legacy DM outside the standalone framebuffer fallback quarantine")
    endif()
  endforeach()
endfunction()

function(_brlobol_guard_check_qged_edit_preview_policy)
  _brlobol_guard_collect(_files
    src/qged)
  foreach(_file IN LISTS _files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if(NOT "${_rel}" MATCHES "${_brlobol_guard_extensions}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    string(REGEX MATCH [[QgObolEditPreviewPrivate\.h|(^|[^A-Za-z0-9_])qg_obol_edit_preview_[A-Za-z0-9_]*[ \t\r\n]*\(]]
      _hit "${_contents}")
    if(_hit)
      _brlobol_guard_fail(
	"${_rel} uses the retired qtcad direct-widget edit-preview helper; route qged edit previews through GED view features")
    endif()
  endforeach()
endfunction()

function(_brlobol_guard_check_tclcad_obol_readback_bridge)
  _brlobol_guard_read_rel(_commands "src/libtclcad/commands.c")
  foreach(_needle
      [[#include "ged/draw_obol.h"]]
      [[ged_draw_obol_controller_opaque_for_view]]
      [[ged_draw_obol_framebuffer_present]]
      [[ged_draw_obol_view_display_image]])
    string(FIND "${_commands}" "${_needle}" _idx)
    if(_idx EQUAL -1)
      _brlobol_guard_fail(
	"src/libtclcad/commands.c no longer routes image export readback through the GED/Obol bridge before raw DM fallback")
    endif()
  endforeach()
endfunction()

function(_brlobol_guard_check_fbserv_legacy_framebuffer_quarantine)
  _brlobol_guard_collect(_files
    src/libged
    src/qged
    src/libtclcad)
  foreach(_file IN LISTS _files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if(NOT "${_rel}" MATCHES "${_brlobol_guard_extensions}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])fbs_fbp([^A-Za-z0-9_]|$)]]
      _hit "${_contents}")
    if(_hit)
      _brlobol_guard_fail(
	"${_rel} reaches into fbserv_obj::fbs_fbp; use fbs_set_legacy_framebuffer/fbs_legacy_framebuffer instead")
    endif()
  endforeach()
endfunction()

function(_brlobol_guard_check_fbserv_transport_quarantine)
  set(_transport_fields
    [[(^|[^A-Za-z0-9_])fbs_listener([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])fbs_clients([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])MAX_CLIENTS([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])fbserv_client([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])fbserv_listener([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])fbs_is_listening([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])fbs_listen_on_port([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])fbs_open_server_handler([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])fbs_close_server_handler([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])fbs_open_client_handler([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])fbs_close_client_handler([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])fbs_open_ipc_client_handler([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])fbs_close_ipc_client_handler([^A-Za-z0-9_]|$)]])

  _brlobol_guard_collect(_files
    src/libged
    src/gtools/gsh
    src/qged/fbserv.cpp
    src/libtclcad/fbserv.c)
  foreach(_file IN LISTS _files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if(NOT "${_rel}" MATCHES "${_brlobol_guard_extensions}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    foreach(_pat IN LISTS _transport_fields)
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_guard_fail(
	  "${_rel} reaches into fbserv transport internals; use fbs transport/accessor helpers instead")
      endif()
    endforeach()
  endforeach()
endfunction()

function(_brlobol_guard_check_fbserv_transport_setup_helpers)
  _brlobol_guard_collect(_files
    src/qged/fbserv.cpp
    src/libtclcad/fb.c
    src/libtclcad/fbserv.c
    src/libtclcad/commands.c)

  set(_callback_assignments
    [[(^|[^A-Za-z0-9_])fbs_is_listening[ \t\r\n]*=]]
    [[(^|[^A-Za-z0-9_])fbs_listen_on_port[ \t\r\n]*=]]
    [[(^|[^A-Za-z0-9_])fbs_open_server_handler[ \t\r\n]*=]]
    [[(^|[^A-Za-z0-9_])fbs_close_server_handler[ \t\r\n]*=]]
    [[(^|[^A-Za-z0-9_])fbs_open_client_handler[ \t\r\n]*=]]
    [[(^|[^A-Za-z0-9_])fbs_close_client_handler[ \t\r\n]*=]]
    [[(^|[^A-Za-z0-9_])fbs_open_ipc_client_handler[ \t\r\n]*=]]
    [[(^|[^A-Za-z0-9_])fbs_close_ipc_client_handler[ \t\r\n]*=]])

  foreach(_file IN LISTS _files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    file(READ "${_file}" _contents)
    foreach(_pat IN LISTS _callback_assignments)
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_guard_fail(
	  "${_rel} assigns fbserv transport callbacks directly; use fbs_set_transport or tclcad_fbserv_set_transport")
      endif()
    endforeach()
  endforeach()

  _brlobol_guard_read_rel(_qged "src/qged/fbserv.cpp")
  string(FIND "${_qged}" [[fbs_set_transport]] _qged_transport_idx)
  if(_qged_transport_idx EQUAL -1)
    _brlobol_guard_fail(
      "src/qged/fbserv.cpp no longer installs fbserv transport callbacks through fbs_set_transport")
  endif()

  _brlobol_guard_read_rel(_tcl_fb "src/libtclcad/fb.c")
  _brlobol_guard_read_rel(_tcl_fbserv "src/libtclcad/fbserv.c")
  string(FIND "${_tcl_fb}" [[tclcad_fbserv_set_transport]] _tcl_fb_transport_idx)
  string(FIND "${_tcl_fbserv}" [[tclcad_fbserv_set_transport]] _tcl_fbserv_transport_idx)
  if(_tcl_fb_transport_idx EQUAL -1 OR _tcl_fbserv_transport_idx EQUAL -1)
    _brlobol_guard_fail(
      "libtclcad no longer centralizes fbserv transport callback setup through tclcad_fbserv_set_transport")
  endif()
endfunction()

function(_brlobol_guard_check_qged_fbserv_backend_access)
  _brlobol_guard_collect(_files
    src/qged)
  foreach(_file IN LISTS _files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if(NOT "${_rel}" MATCHES "${_brlobol_guard_extensions}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])fbs_fb_(ops|ctx)([^A-Za-z0-9_]|$)]]
      _hit "${_contents}")
    if(_hit)
      _brlobol_guard_fail(
	"${_rel} reaches into fbserv backend internals; use fbs_framebuffer_* helpers instead")
    endif()
  endforeach()
endfunction()

function(_brlobol_guard_check_fbserv_backend_contract)
  _brlobol_guard_read_rel(_imgstream_fbserv "include/imgstream/fbserv.h")
  _brlobol_guard_read_rel(_imgstream_auth "src/libimgstream/fbserv_auth.c")
  _brlobol_guard_read_rel(_imgstream_framebuffer "src/libimgstream/fbserv_framebuffer.c")
  _brlobol_guard_read_rel(_imgstream_pkg "src/libimgstream/fbserv_pkg.c")
  _brlobol_guard_read_rel(_imgstream_server "src/libimgstream/fbserv_server.c")
  _brlobol_guard_read_rel(_imgstream_state "src/libimgstream/fbserv_state.c")
  _brlobol_guard_read_rel(_libdm_cmake "src/libdm/CMakeLists.txt")
  _brlobol_guard_read_rel(_dm_h "include/dm.h")
  _brlobol_guard_read_rel(_dm_fbserv "include/dm/fbserv.h")
  _brlobol_guard_read_rel(_protocol_note
    "doc/notes/obol_framebuffer_protocol.txt")

  string(FIND "${_imgstream_fbserv}" [[struct fbserv_fb_ops]] _ops_idx)
  if(_ops_idx EQUAL -1)
    _brlobol_guard_fail(
      "include/imgstream/fbserv.h no longer owns the fbserv framebuffer backend operation table")
  endif()
  string(FIND "${_imgstream_fbserv}" [[FBSERV_MSG_FBOPEN]] _msg_idx)
  if(_msg_idx EQUAL -1)
    _brlobol_guard_fail(
      "include/imgstream/fbserv.h no longer owns the fbserv wire-protocol message IDs")
  endif()
  string(FIND "${_imgstream_fbserv}" [[struct fbserv_obj]] _obj_idx)
  string(FIND "${_imgstream_fbserv}" [[struct fbserv_transport_ops]] _transport_ops_idx)
  if(_obj_idx EQUAL -1 OR _transport_ops_idx EQUAL -1)
    _brlobol_guard_fail(
      "include/imgstream/fbserv.h no longer owns the transitional fbserv object/transport type contract")
  endif()
  string(FIND "${_imgstream_state}" [[fbserv_obj_init]] _state_impl_idx)
  string(FIND "${_imgstream_server}" [[fbs_init]] _server_state_wrapper_idx)
  if(_state_impl_idx EQUAL -1 OR _server_state_wrapper_idx EQUAL -1)
    _brlobol_guard_fail(
      "libimgstream no longer owns fbserv state helpers and complete fbs_* lifecycle wiring")
  endif()
  string(FIND "${_imgstream_pkg}" [[fbserv_pkg_switch_init]] _pkg_switch_impl_idx)
  string(FIND "${_imgstream_server}" [[fbserv_pkg_switch_init]] _server_pkg_switch_wrapper_idx)
  if(_pkg_switch_impl_idx EQUAL -1 OR _server_pkg_switch_wrapper_idx EQUAL -1)
    _brlobol_guard_fail(
      "libimgstream no longer owns the fbserv PKG message switch and fbs_* wiring")
  endif()
  string(FIND "${_imgstream_framebuffer}" [[fbserv_backend_info]] _backend_dispatch_idx)
  string(FIND "${_imgstream_server}" [[fbserv_backend_info]] _server_backend_dispatch_idx)
  if(_backend_dispatch_idx EQUAL -1 OR _server_backend_dispatch_idx EQUAL -1)
    _brlobol_guard_fail(
      "libimgstream no longer owns fbserv backend-operation dispatch and server wiring")
  endif()
  foreach(_entry
      [[fbs_open(]]
      [[fbs_close(]]
      [[fbs_new_client(]]
      [[fbs_open_ipc(]]
      [[fbs_process_client_bytes(]]
      [[fbs_drain_client(]])
    string(FIND "${_imgstream_server}" "${_entry}" _server_lifecycle_idx)
    if(_server_lifecycle_idx EQUAL -1)
      _brlobol_guard_fail(
	"src/libimgstream/fbserv_server.c no longer owns required lifecycle entry ${_entry}")
    endif()
  endforeach()
  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libdm/fbserv.c")
    _brlobol_guard_fail(
      "src/libdm/fbserv.c reintroduced libdm ownership of the framebuffer protocol server")
  endif()
  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/fbserv/server.c")
    _brlobol_guard_fail(
      "src/fbserv/server.c reintroduced a duplicate standalone framebuffer protocol server")
  endif()
  _brlobol_guard_read_rel(_mged_fbserv "src/mged/fbserv.c")
  foreach(_token [[pkg_process(]] [[pkg_switch]] [[MSG_FBOPEN]])
    string(FIND "${_mged_fbserv}" "${_token}" _mged_protocol_idx)
    if(NOT _mged_protocol_idx EQUAL -1)
      _brlobol_guard_fail(
	"src/mged/fbserv.c reintroduced private framebuffer protocol processing (${_token})")
    endif()
  endforeach()
  string(FIND "${_libdm_cmake}" "  fbserv.c" _libdm_fbserv_source_idx)
  if(NOT _libdm_fbserv_source_idx EQUAL -1)
    _brlobol_guard_fail(
      "src/libdm/CMakeLists.txt reintroduced the framebuffer protocol server core")
  endif()
  string(FIND "${_imgstream_fbserv}" [[fbserv_generate_token]] _auth_decl_idx)
  string(FIND "${_imgstream_auth}" [[fbserv_verify_token]] _auth_impl_idx)
  if(_auth_decl_idx EQUAL -1 OR _auth_impl_idx EQUAL -1)
    _brlobol_guard_fail(
      "libimgstream no longer owns fbserv session-token authentication helpers")
  endif()
  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/fbserv/auth.h")
    _brlobol_guard_fail(
      "src/fbserv/auth.h reintroduced header-local fbserv auth helpers; use include/imgstream/fbserv.h")
  endif()
  string(REGEX MATCH [[#[ \t]*define[ \t]+MSG_FB[A-Z0-9_]*]]
    _dm_msg_def "${_dm_h}")
  if(_dm_msg_def)
    _brlobol_guard_fail(
      "include/dm.h redefined fbserv wire-protocol message IDs instead of using include/imgstream/fbserv.h")
  endif()
  string(FIND "${_dm_fbserv}" [[#include "imgstream/fbserv.h"]] _include_idx)
  if(_include_idx EQUAL -1)
    _brlobol_guard_fail(
      "include/dm/fbserv.h no longer imports the imgstream-owned fbserv backend contract")
  endif()
  string(REGEX MATCH [[struct[ \t\r\n]+fbserv_fb_ops[ \t\r\n]*\{]]
    _dm_ops_def "${_dm_fbserv}")
  if(_dm_ops_def)
    _brlobol_guard_fail(
      "include/dm/fbserv.h redefined struct fbserv_fb_ops instead of using include/imgstream/fbserv.h")
  endif()
  string(REGEX MATCH [[struct[ \t\r\n]+fbserv_(obj|listener|client|transport_ops)[ \t\r\n]*\{]]
    _dm_obj_def "${_dm_fbserv}")
  if(_dm_obj_def)
    _brlobol_guard_fail(
      "include/dm/fbserv.h redefined fbserv object/transport types instead of using include/imgstream/fbserv.h")
  endif()
  string(FIND "${_protocol_note}" [[include/imgstream/fbserv.h]] _note_idx)
  if(_note_idx EQUAL -1)
    _brlobol_guard_fail(
      "doc/notes/obol_framebuffer_protocol.txt no longer documents imgstream fbserv backend/protocol contract ownership")
  endif()
endfunction()

function(_brlobol_guard_check_tkobol_host_ownership)
  _brlobol_guard_read_rel(_tclcad_cmake "src/libtclcad/CMakeLists.txt")
  _brlobol_guard_read_rel(_tclcad_dm "src/libtclcad/dm.c")
  _brlobol_guard_read_rel(_dm_plugins "src/libdm/dm_plugins.cpp")
  foreach(_needle
      [[tkobol/dm-tkobol.cpp]]
      [[tkobol/vendor/togl/togl.c]])
    string(FIND "${_tclcad_cmake}" "${_needle}" _host_source_idx)
    if(_host_source_idx EQUAL -1)
      _brlobol_guard_fail(
	"src/libtclcad/CMakeLists.txt no longer owns required Tk Obol host source ${_needle}")
    endif()
  endforeach()
  foreach(_needle [[dm_register_backend]] [[BrlcadTkObolHost_Init]])
    string(FIND "${_tclcad_dm}" "${_needle}" _host_init_idx)
    if(_host_init_idx EQUAL -1)
      _brlobol_guard_fail(
	"src/libtclcad/dm.c no longer registers the built-in Tk Obol host (${_needle})")
    endif()
  endforeach()
  _brlobol_guard_forbid_regexes("src/libtclcad/CMakeLists.txt" "${_tclcad_cmake}"
    [[add_subdirectory\(tkobol\)]]
    [[dm_plugin_library\(dm-tkobol]])
  _brlobol_guard_forbid_regexes("src/libdm/dm_plugins.cpp" "${_dm_plugins}"
    [[DM_SWRAST]]
    [[(^|[^A-Za-z0-9_])(ogl|wgl|swrast|qtgl|tkswrast)([^A-Za-z0-9_]|$)]])

  _brlobol_guard_read_rel(_gsh_cmake "src/gtools/gsh/CMakeLists.txt")
  _brlobol_guard_forbid_regexes("src/gtools/gsh/CMakeLists.txt" "${_gsh_cmake}"
    [[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]])
endfunction()

function(_brlobol_guard_collect_active_scan_files _outvar)
  _brlobol_guard_collect(_files
    include/brlobol
    include/brlobol.h
    include/qtcad
    src/libbrlobol
    src/libqtcad
    src/qged
    src/gtools/gsh)
  set(${_outvar} "${_files}" PARENT_SCOPE)
endfunction()

_brlobol_guard_check_dependency_inventory()
_brlobol_guard_check_libbsg_retired()
_brlobol_guard_check_legacy_graphical_dm_retired()
_brlobol_guard_check_retired_tests()
_brlobol_guard_check_public_headers()
_brlobol_guard_check_retired_ged_draw_symbols()
_brlobol_guard_check_removed_rt_view_aliases()
_brlobol_guard_check_native_drawing_tests()
_brlobol_guard_check_obol_controller_quarantine()
_brlobol_guard_check_qtcad_legacy_dm_open_quarantine()
_brlobol_guard_check_qged_edit_preview_policy()
_brlobol_guard_check_tclcad_obol_readback_bridge()
_brlobol_guard_check_fbserv_legacy_framebuffer_quarantine()
_brlobol_guard_check_fbserv_transport_quarantine()
_brlobol_guard_check_fbserv_transport_setup_helpers()
_brlobol_guard_check_qged_fbserv_backend_access()
_brlobol_guard_check_fbserv_backend_contract()
_brlobol_guard_check_tkobol_host_ownership()

_brlobol_guard_collect_active_scan_files(_brlobol_guard_scan_files)
foreach(_file IN LISTS _brlobol_guard_scan_files)
  _brlobol_guard_check_file("${_file}")
endforeach()

get_property(_failures GLOBAL PROPERTY BRLOBOL_PIVOT_GUARD_FAILURES)
if(_failures)
  list(SORT _failures)
  string(REPLACE ";" "\n  " _failure_text "${_failures}")
  message(FATAL_ERROR
    "Obol pivot guard failed:\n"
    "  ${_failure_text}\n"
    "New Obol-canonical code must not depend on retired BSG/DM APIs.")
endif()

message(STATUS "Obol pivot guard passed")
