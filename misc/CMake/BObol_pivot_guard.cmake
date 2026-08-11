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

include("${CMAKE_CURRENT_LIST_DIR}/BObol_pivot_guard_rules.cmake")

set(_bobol_guard_extensions "${BOBOL_GUARD_SOURCE_EXTENSIONS}")
set(_bobol_guard_active_boundary_patterns
  ${BOBOL_GUARD_FORBIDDEN_DEPENDENCY_PATTERNS}
  ${BOBOL_GUARD_FORBIDDEN_OWNERSHIP_PATTERNS})
set(_bobol_guard_retired_ged_draw_symbols
  ${BOBOL_GUARD_RETIRED_GED_SYMBOL_PATTERNS})
set(_bobol_guard_retired_test_files ${BOBOL_GUARD_TEMPORARY_TEST_FILES})
set(_bobol_guard_retired_test_tokens ${BOBOL_GUARD_TEMPORARY_TEST_TOKENS})

function(_bobol_guard_fail _msg)
  set_property(GLOBAL APPEND PROPERTY BOBOL_PIVOT_GUARD_FAILURES "${_msg}")
endfunction()

function(_bobol_guard_collect _outvar)
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

function(_bobol_guard_read_rel _outvar _rel)
  set(_path "${BRLCAD_SOURCE_DIR}/${_rel}")
  if(NOT EXISTS "${_path}")
    _bobol_guard_fail("${_rel} is required for Obol pivot guard checks")
    set(${_outvar} "" PARENT_SCOPE)
    return()
  endif()
  file(READ "${_path}" _contents)
  set(${_outvar} "${_contents}" PARENT_SCOPE)
endfunction()

function(_bobol_guard_forbid_regexes _rel _contents)
  foreach(_pat IN LISTS ARGN)
    string(REGEX MATCH "${_pat}" _hit "${_contents}")
    if(_hit)
      _bobol_guard_fail(
        "${_rel}: forbidden dependency/ownership invariant matched '${_hit}'; route drawing through the typed display endpoint and direct libBObol service boundary")
    endif()
  endforeach()
endfunction()

function(_bobol_guard_check_file _file)
  file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
  if("${_rel}" STREQUAL "misc/CMake/BObol_pivot_guard.cmake" OR
     "${_rel}" STREQUAL "misc/CMake/BObol_pivot_guard_rules.cmake")
    return()
  endif()
  if(NOT "${_rel}" MATCHES "${_bobol_guard_extensions}")
    return()
  endif()

  file(READ "${_file}" _contents)
  _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
    ${_bobol_guard_active_boundary_patterns})
endfunction()

function(_bobol_guard_check_dependency_inventory)
  _bobol_guard_read_rel(_inventory
    "doc/notes/obol_legacy_dependency_inventory.txt")
  _bobol_guard_read_rel(_source_dirs "src/source_dirs.cmake")

  set(_expected_deps_libBObol [[libwdb;librt;libimgstream;libbg;libbu]])
  set(_expected_deps_libbv [[libbg;libbn;libbu]])
  set(_expected_deps_libnmg [[libbg;libbn;libbu]])
  set(_expected_deps_libbrep [[libbg;libbu]])
  set(_expected_deps_librt [[libbrep;libnmg;libbv;libbg;libbn;libbu]])
  set(_expected_deps_libanalyze [[librt;libbg;libbn;libbu]])
  set(_expected_deps_libimgstream [[libicv;libbn;libpkg;libbu]])
  set(_expected_deps_libged [[libBObol;libimgstream;libicv;libanalyze;libwdb;liboptical;libbu]])
  set(_expected_deps_libqtcad [[libBObol;libged;libbg;libbn;libbu]])
  set(_expected_deps_libtclcad [[libged;libimgstream;libpkg;libbn;libbu]])
  set(_expected_deps_qged [[libqtcad]])
  set(_expected_deps_mged [[libtclcad]])
  set(_expected_deps_fbserv [[libimgstream;libpkg;libbu;libqtcad;libtclcad]])
  set(_expected_deps_rt [[librt;liboptical;libimgstream;libicv;libqtcad;libtclcad]])

  foreach(_target
      libBObol
      libbv
      libnmg
      libbrep
      librt
      libanalyze
      libimgstream
      libged
      libqtcad
      libtclcad
      qged
      mged
      fbserv
      rt)
    set(_expected "${_expected_deps_${_target}}")
    string(REGEX MATCH "set_deps\\(${_target}[ \t]+\"([^\"]*)\"\\)"
      _row "${_source_dirs}")
    set(_actual "${CMAKE_MATCH_1}")
    if(NOT _row OR NOT "${_actual}" STREQUAL "${_expected}")
      string(REPLACE ";" "," _expected_msg "${_expected}")
      _bobol_guard_fail(
	"src/source_dirs.cmake direct dependency row changed for ${_target}; expected ${_expected_msg}")
    endif()
    string(FIND "${_inventory}" "* ${_target}: ${_expected}" _inventory_idx)
    if(_inventory_idx EQUAL -1)
      string(REPLACE ";" "," _expected_msg "${_expected}")
      _bobol_guard_fail(
	"doc/notes/obol_legacy_dependency_inventory.txt missing ${_target}: ${_expected_msg}")
    endif()
  endforeach()

  _bobol_guard_read_rel(_qged_cmake "src/qged/CMakeLists.txt")
  foreach(_retired
      include/dm.h
      include/dm/CMakeLists.txt
      src/libdm/CMakeLists.txt
      misc/pkgconfig/libdm.pc.in)
    if(EXISTS "${BRLCAD_SOURCE_DIR}/${_retired}")
      _bobol_guard_fail("${_retired} restored the retired libdm compatibility layer")
    endif()
  endforeach()
endfunction()

function(_bobol_guard_check_libbsg_retired)
  foreach(_rel include/bsg/CMakeLists.txt include/bsg.h src/libbsg/CMakeLists.txt)
    if(EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _bobol_guard_fail("${_rel} reintroduced the retired libbsg API or implementation")
    endif()
  endforeach()

  _bobol_guard_collect(_files include src misc/CMake)
  foreach(_file IN LISTS _files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if("${_rel}" STREQUAL "misc/CMake/BObol_pivot_guard.cmake" OR
       "${_rel}" STREQUAL "misc/CMake/BObol_pivot_guard_rules.cmake" OR
       NOT "${_rel}" MATCHES "${_bobol_guard_extensions}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[#[ \t]*include[ \t]*[<"]bsg]]
      [[(^|[^A-Za-z0-9_])libbsg([^A-Za-z0-9_]|$)]])
  endforeach()
endfunction()

function(_bobol_guard_check_legacy_graphical_dm_retired)
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
      _bobol_guard_fail("${_rel} reintroduced a retired graphical display-manager path")
    endif()
  endforeach()
  # ISST is an ADRT image producer.  Its Tk photo presenter is intentionally
  # pixel-only; restoring a private OpenGL/display-manager window would create
  # a second graphical path outside the endpoint architecture.
  _bobol_guard_read_rel(_isst_source "src/adrt/isst.c")
  _bobol_guard_forbid_regexes("src/adrt/isst.c" "${_isst_source}"
    [[#[ \t]*include[ \t]*[<\"]GL/]]
    [[#[ \t]*include[ \t]*[<\"]dm\.h]]
    [[(^|[^A-Za-z0-9_])dm_(open|make_current|draw_end|set_bg)[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])gl[A-Z][A-Za-z0-9_]*[ \t\r\n]*\(]])
  _bobol_guard_read_rel(_isst_cmake "src/adrt/CMakeLists.txt")
  _bobol_guard_forbid_regexes("src/adrt/CMakeLists.txt" "${_isst_cmake}"
    [[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])dm_plugins([^A-Za-z0-9_]|$)]]
    [[BRLCAD_ENABLE_OPENGL[ \t]*AND[ \t]*BRLCAD_ENABLE_TK]])

  # QISST is also a CPU image producer.  Its renderer runs in a worker thread
  # and publishes RGB QImage frames to a QWidget; it must not recover the
  # former OpenGL context-transfer and texture-presenter path.
  foreach(_rel src/isst/isstgl.h src/isst/isstgl.cpp src/isst/main_window.h)
    _bobol_guard_read_rel(_qisst_source "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_qisst_source}"
      [[QOpenGL]]
      [[QGLWidget]]
      [[(^|[^A-Za-z0-9_])gl[A-Z][A-Za-z0-9_]*[ \t\r\n]*\(]])
  endforeach()
  _bobol_guard_read_rel(_qisst_cmake "src/isst/CMakeLists.txt")
  _bobol_guard_forbid_regexes("src/isst/CMakeLists.txt" "${_qisst_cmake}"
    [[BRLCAD_ENABLE_QT[ \t]*AND[ \t]*BRLCAD_ENABLE_OPENGL]]
    [[Qt[56]::OpenGL]]
    [[OPENGL_(opengl|gl)_LIBRARY]])

  foreach(_rel src/liboptical/sh_points.c src/liboptical/sh_spm.c)
    _bobol_guard_read_rel(_optical_source "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_optical_source}"
      [[#[ \t]*include[ \t]*[<\"]dm\.h]])
  endforeach()
endfunction()

function(_bobol_guard_check_retired_tests)
  foreach(_rel IN LISTS _bobol_guard_retired_test_files)
    if(EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _bobol_guard_fail(
	"${_rel} was retired; keep migrated coverage in Obol/native GED draw tests")
    endif()
  endforeach()

  _bobol_guard_read_rel(_cmake "src/libged/tests/draw/CMakeLists.txt")
  foreach(_token IN LISTS _bobol_guard_retired_test_tokens)
    string(REGEX MATCH "${_token}" _hit "${_cmake}")
    if(_hit)
      _bobol_guard_fail(
	"src/libged/tests/draw/CMakeLists.txt reintroduced retired BSG test token ${_token}")
    endif()
  endforeach()

endfunction()

function(_bobol_guard_check_public_headers)
  foreach(_rel include/ged.h include/ged/defines.h include/ged/view.h)
    _bobol_guard_read_rel(_contents "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[#[ \t]*include[ \t]*[<"]bsg]]
      [[#[ \t]*include[ \t]*[<"]dm/fbserv\.h]]
      [[#[ \t]*include[ \t]*[<"]dm]]
      [[(^|[^A-Za-z0-9_])bsg_view([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])bsg_scene_ref([^A-Za-z0-9_]|$)]])
  endforeach()
  foreach(_retired_qtcad_namespace_header
      include/qtcad/QgNamespace.h
      include/qtcad/QgNamespaceCompat.h)
    if(EXISTS "${BRLCAD_SOURCE_DIR}/${_retired_qtcad_namespace_header}")
      _bobol_guard_fail(
        "${_retired_qtcad_namespace_header} restored a retired Qt namespace compatibility surface")
    endif()
  endforeach()
  _bobol_guard_read_rel(_qg_types "include/qtcad/QgTypes.h")
  foreach(_alias
      QgView_AUTO
      QgView_SW
      QgView_GL
      UPPER_RIGHT_QUADRANT
      UPPER_LEFT_QUADRANT
      LOWER_LEFT_QUADRANT
      LOWER_RIGHT_QUADRANT
      QgViewMode
      QgInstanceEditMode
      QgPrimitiveEditMode
      G_STANDARD_OBJ
      G_REGION
      G_AIR
      G_AIR_REGION
      G_ASSEMBLY)
    string(FIND "${_qg_types}" "${_alias}" _qg_type_alias_idx)
    if(NOT _qg_type_alias_idx EQUAL -1)
      _bobol_guard_fail(
        "include/qtcad/QgTypes.h reintroduced retired integer compatibility alias ${_alias}")
    endif()
  endforeach()
endfunction()

function(_bobol_guard_check_retired_ged_draw_symbols)
  foreach(_rel
      src/libged/ged_draw_source.cpp
      src/libged/ged_draw_transactions.c
      src/libged/draw_obol.cpp
      include/ged/draw.h)
    if(EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      file(READ "${BRLCAD_SOURCE_DIR}/${_rel}" _contents)
      _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
	${_bobol_guard_retired_ged_draw_symbols})
    endif()
  endforeach()
endfunction()

function(_bobol_guard_check_removed_rt_view_aliases)
  _bobol_guard_read_rel(_rt_view "include/rt/view.h")
  foreach(_token
      rt_view_info
      rt_view_lod_policy
      rt_view_lod_settings
      RT_VIEW_INFO_INIT
      RT_VIEW_LOD_POLICY_INIT
      RT_VIEW_LOD_SETTINGS_INIT)
    string(FIND "${_rt_view}" "${_token}" _idx)
    if(NOT _idx EQUAL -1)
      _bobol_guard_fail(
	"include/rt/view.h reintroduced removed libbv compatibility alias ${_token}")
    endif()
  endforeach()
  foreach(_rel src/librt/view.c src/librt/tests/view_info.c)
    if(EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _bobol_guard_fail("${_rel} reintroduced removed libbv wrapper code")
    endif()
  endforeach()
endfunction()

function(_bobol_guard_check_native_drawing_tests)
  foreach(_rel
      src/libged/tests/draw/basic.cpp
      src/libged/tests/draw/faceplate.cpp
      src/libged/tests/draw/lod.cpp
      src/libged/tests/draw/select.cpp
      src/libged/tests/draw/quad.cpp
      src/libged/tests/draw/util.cpp)
    _bobol_guard_read_rel(_contents "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[#[ \t]*include[ \t]*[<"]dm]]
      [[(^|[^A-Za-z0-9_])DM_SWRAST([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])ged_exec_dm[ \t\r\n]*\(]]
      [[(^|[^A-Za-z0-9_])ged_view_context_display_manager_get[ \t\r\n]*\(]])
  endforeach()

  _bobol_guard_read_rel(_cmake "src/libged/tests/draw/CMakeLists.txt")
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
      _bobol_guard_fail(
	"src/libged/tests/draw/CMakeLists.txt links native Obol test ${_target} to libdm")
    endif()
  endforeach()
endfunction()

function(_bobol_guard_check_obol_controller_quarantine)
  foreach(_retired
      include/dm/obol.h
      src/libdm/obol/CMakeLists.txt
      src/libdm/obol/dm-obol.h
      src/libdm/obol/dm-obol.cpp)
    if(EXISTS "${BRLCAD_SOURCE_DIR}/${_retired}")
      _bobol_guard_fail("${_retired} restored the retired dm-obol adapter")
    endif()
  endforeach()
  _bobol_guard_collect(_files
    include
    src/libBObol
    src/libdm
    src/libged
    src/libqtcad
    src/qged
    src/gtools/gsh
    src/mged
    src/libtclcad)
  foreach(_file IN LISTS _files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if(NOT "${_rel}" MATCHES "${_bobol_guard_extensions}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])dm_obol_controller[ \t\r\n]*\(]]
      _hit "${_contents}")
    if(_hit)
      _bobol_guard_fail(
	"${_rel} calls dm_obol_controller outside the libdm/GED Obol bridge quarantine")
    endif()
  endforeach()
endfunction()

function(_bobol_guard_check_display_endpoint_boundary)
  foreach(_rel
      include/BObol/BDisplayEndpoint.h
      include/BObol/BHostFactory.h
      src/libBObol/display_endpoint.cpp
      src/libBObol/host_factory.cpp)
    if(NOT EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _bobol_guard_fail("${_rel} is required for the typed display endpoint boundary")
    endif()
  endforeach()

  _bobol_guard_read_rel(_endpoint_h "include/BObol/BDisplayEndpoint.h")
  _bobol_guard_read_rel(_endpoint_cpp "src/libBObol/display_endpoint.cpp")
  foreach(_needle
      [[bobol_display_endpoint_property_count]]
      [[bobol_display_endpoint_property_descriptor]]
      [[bobol_display_endpoint_property_get]]
      [[bobol_display_endpoint_property_set]]
      [[bobol_display_endpoint_property_provider_set]]
      [[bobol_display_endpoint_capture_plane]]
      [[bobol_display_endpoint_framebuffer_capture_provider_set]]
      [[bobol_display_endpoint_host_capabilities]])
    string(FIND "${_endpoint_h}" "${_needle}" _endpoint_property_decl_idx)
    string(FIND "${_endpoint_cpp}" "${_needle}" _endpoint_property_impl_idx)
    if(_endpoint_property_decl_idx EQUAL -1 OR
       _endpoint_property_impl_idx EQUAL -1)
      _bobol_guard_fail(
	"typed display endpoint property contract is incomplete (${_needle})")
    endif()
  endforeach()

  _bobol_guard_read_rel(_ged_view_state "src/libged/ged_view_state.cpp")
  foreach(_needle
      [[bobol_display_endpoint_t *display_endpoint]]
      [[ged_view_context_display_endpoint_get]]
      [[ged_view_context_display_endpoint_set]])
    string(FIND "${_ged_view_state}" "${_needle}" _endpoint_idx)
    if(_endpoint_idx EQUAL -1)
      _bobol_guard_fail(
	"src/libged/ged_view_state.cpp no longer maintains ${_needle}")
    endif()
  endforeach()
  _bobol_guard_forbid_regexes("src/libged/ged_view_state.cpp"
    "${_ged_view_state}"
    [[ged_draw_obol_scene_controller_owned]])

  foreach(_rel
      include/ged/view.h
      src/libged/ged_view_state.cpp
      src/mged/cmd.cpp
      src/mged/overlay.c
      src/mged/rtif.c)
    _bobol_guard_read_rel(_contents "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[(^|[^A-Za-z0-9_])ged_view_context_display_manager_(get|set)[ \t\r\n]*\(]])
  endforeach()

  foreach(_rel
      include/ged/display.h
      src/libged/draw_obol.cpp)
    _bobol_guard_read_rel(_contents "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[GED_DRAW_OBOL_SOFTWARE_WIRE_(AUTO|QUALITY|FAST)]]
      [[ged_draw_obol_software_wire_mode_(get|set)_for_view]])
  endforeach()

  _bobol_guard_read_rel(_bobol_input_h "include/BObol/BInput.h")
  _bobol_guard_forbid_regexes("include/BObol/BInput.h"
    "${_bobol_input_h}"
    [[(^|[^A-Za-z0-9_])BOBOL_INPUT_KEY([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])BOBOL_INPUT_POINTER([^A-Za-z0-9_]|$)]])

  foreach(_rel
      src/libged/ged.cpp
      src/libged/ged_private.h
      src/libged/dm/dm.cpp
      src/libged/rt/rt.c)
    _bobol_guard_read_rel(_contents "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[ged_rt_fb_(set|get|refresh)]]
      [[rt_fb_dev]])
  endforeach()
  foreach(_rel
      include/ged/defines.h
      src/libged/ged.cpp
      src/libged/ged_private.h
      src/libged/dm/dm.cpp)
    _bobol_guard_read_rel(_contents "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[ged_dm_ctx_]]
      [[dm_map]])
  endforeach()
  _bobol_guard_read_rel(_ged_rt "src/libged/rt/rt.c")
  _bobol_guard_forbid_regexes("src/libged/rt/rt.c" "${_ged_rt}"
    [[/dev/ogl]])

  _bobol_guard_read_rel(_ged_dm "src/libged/dm/dm.cpp")
  _bobol_guard_forbid_regexes("src/libged/dm/dm.cpp" "${_ged_dm}"
    [[dm_get_bg]]
    [[dm_set_bg]]
    [[dm_get_zclip]]
    [[dm_set_zclip]]
    [[dm_get_zbuffer]]
    [[dm_set_zbuffer]]
    [[dm_get_light]]
    [[dm_set_light]]
    [[dm_get_vparse]]
    [[dm_get_mvars]]
    [[bu_struct_parse]]
    [[dm_open[ \t\r\n]*\(]]
    [[dm_list_types[ \t\r\n]*\(]]
    [[dm_init_msgs[ \t\r\n]*\(]]
    [[dm_obol_controller[ \t\r\n]*\(]]
    [[_dm_cmd_debug]]
    [[_dm_cmd_width]]
    [[_dm_cmd_height]]
    [[_dm_cmd_attach]]
    [[_dm_cmd_type]]
    [[_dm_cmd_types]]
    [[_dm_cmd_initmsg]]
    [[#[ \t]*include[ \t]*[<"]dm\.h]]
    [[struct[ \t]+dm]])
  foreach(_needle
      [[_dm_cmd_status]]
      [[_dm_cmd_renderer]]
      [[_dm_cmd_open]]
      [[_dm_cmd_close]]
      [[_dm_cmd_host]]
      [[_dm_cmd_diagnostics]]
      [[bobol_display_endpoint_property_get]]
      [[bobol_display_endpoint_property_set]]
      [[bobol_host_factory_registry_count]])
    string(FIND "${_ged_dm}" "${_needle}" _ged_dm_endpoint_idx)
    if(_ged_dm_endpoint_idx EQUAL -1)
      _bobol_guard_fail(
	"GED dm command no longer exposes endpoint-native operations (${_needle})")
    endif()
  endforeach()

  _bobol_guard_read_rel(_ged_screengrab "src/libged/dm/screengrab.cpp")
  string(FIND "${_ged_screengrab}"
    [[bobol_display_endpoint_capture]] _endpoint_capture_idx)
  if(_endpoint_capture_idx EQUAL -1)
    _bobol_guard_fail(
      "GED screengrab no longer captures through the display endpoint")
  endif()
  _bobol_guard_forbid_regexes("src/libged/dm/screengrab.cpp"
    "${_ged_screengrab}"
    [[#[ \t]*include[ \t]*[<"]dm\.h]]
    [[struct[ \t]+dm]]
    [[dm_get_fb]]
    [[dm_get_display_image]]
    [[ged_draw_obol_view_display_image]]
    [[fb_write_icv]]
    [[framebuffer-only capture is not supported]])
  foreach(_needle
      [[bobol_display_endpoint_capture_plane]]
      [[BOBOL_CAPTURE_FRAMEBUFFER]])
    string(FIND "${_ged_screengrab}" "${_needle}" _capture_plane_idx)
    if(_capture_plane_idx EQUAL -1)
      _bobol_guard_fail(
	"GED screengrab no longer selects endpoint framebuffer composition (${_needle})")
    endif()
  endforeach()

  _bobol_guard_read_rel(_ged_fbserv "src/libged/obol_fbserv.cpp")
  foreach(_needle
      [[ged_obol_capture_framebuffer]]
      [[bobol_display_endpoint_framebuffer_capture_provider_set]])
    string(FIND "${_ged_fbserv}" "${_needle}" _fb_provider_idx)
    if(_fb_provider_idx EQUAL -1)
      _bobol_guard_fail(
	"libged framebuffer bridge no longer owns endpoint capture provider lifecycle (${_needle})")
    endif()
  endforeach()

  _bobol_guard_read_rel(_ged_ert "src/libged/dm/ert.cpp")
  _bobol_guard_forbid_regexes("src/libged/dm/ert.cpp" "${_ged_ert}"
    [[#[ \t]*include[ \t]*[<"]dm\.h]]
    [[#[ \t]*include[ \t]*[<"]dm/fbserv_legacy\.h]]
    [[struct[ \t]+dm]]
    [[dm_get_fb]]
    [[dm_fbserv_set_framebuffer]])

  _bobol_guard_read_rel(_ged_dm_cmake "src/libged/dm/CMakeLists.txt")
  _bobol_guard_forbid_regexes("src/libged/dm/CMakeLists.txt"
    "${_ged_dm_cmake}"
    [[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]])

  _bobol_guard_read_rel(_ged_get_autoview_cmake
    "src/libged/get_autoview/CMakeLists.txt")
  _bobol_guard_forbid_regexes("src/libged/get_autoview/CMakeLists.txt"
    "${_ged_get_autoview_cmake}"
    [[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]])

  _bobol_guard_collect(_libged_files src/libged)
  foreach(_file IN LISTS _libged_files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if(NOT "${_rel}" MATCHES "${_bobol_guard_extensions}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[#[ \t]*include[ \t]*[<"]dm\.h]])
  endforeach()

  _bobol_guard_read_rel(_tclcad_commands "src/libtclcad/commands.c")
  string(FIND "${_tclcad_commands}" [[{"dm",]] _tclcad_dm_name_idx)
  string(FIND "${_tclcad_commands}"
    [[to_pass_through_func, ged_exec_dm}]] _tclcad_dm_exec_idx)
  if(_tclcad_dm_name_idx EQUAL -1 OR _tclcad_dm_exec_idx EQUAL -1)
    _bobol_guard_fail(
      "TclCAD no longer exposes the endpoint-native GED dm command")
  endif()

  _bobol_guard_read_rel(_canvas_state "src/libqtcad/QgCanvasState.h")
  foreach(_needle
      [[bool owns_obol = false]]
      [[if (s.owns_obol)]]
      [[if (!s.obol && create_controller)]]
      [[qgcanvas_bind_obol_controller]])
    string(FIND "${_canvas_state}" "${_needle}" _canvas_idx)
    if(_canvas_idx EQUAL -1)
      _bobol_guard_fail(
	"src/libqtcad/QgCanvasState.h no longer preserves endpoint-owned controller borrowing (${_needle})")
    endif()
  endforeach()

  _bobol_guard_read_rel(_qg_view "src/libqtcad/QgView.cpp")
  foreach(_needle
      [[bobol_display_endpoint_create]]
      [[bobol_display_endpoint_render_engine_set]]
      [[bobol_display_endpoint_host_open]]
      [[BOBOL_HOST_MODE_EMBEDDED]]
      [[new QgGL(parent, nullptr, false)]]
      [[new QgSW(parent, nullptr, false)]])
    string(FIND "${_qg_view}" "${_needle}" _qg_endpoint_idx)
    if(_qg_endpoint_idx EQUAL -1)
      _bobol_guard_fail(
	"QgView no longer owns its visible canvas through an embedded endpoint (${_needle})")
    endif()
  endforeach()

  _bobol_guard_read_rel(_qg_session "src/libqtcad/QgSession.cpp")
  _bobol_guard_forbid_regexes("src/libqtcad/QgSession.cpp" "${_qg_session}"
    [[#[ \t]*include[ \t]*[<\"]dm]]
    [[(^|[^A-Za-z0-9_])DM_SWRAST([^A-Za-z0-9_]|$)]])

  _bobol_guard_read_rel(_qg_host "src/libqtcad/QgObolWindowHost.cpp")
  foreach(_needle
      [[new QgGL(NULL, nullptr, false)]]
      [[new QgSW(NULL, nullptr, false)]]
      [[factory.set_title = qg_factory_set_title]]
      [[factory.set_visible = qg_factory_set_visible]]
      [[factory.set_vsync = qg_factory_set_vsync]]
      [[BOBOL_HOST_CAP_PRESENT_VSYNC]])
    string(FIND "${_qg_host}" "${_needle}" _qg_host_controller_idx)
    if(_qg_host_controller_idx EQUAL -1)
      _bobol_guard_fail(
	"Qt endpoint host factories reintroduced a temporary canvas-owned controller (${_needle})")
    endif()
  endforeach()

  _bobol_guard_read_rel(_endpoint "src/libBObol/display_endpoint.cpp")
  foreach(_needle
      [["endpoint.title"]]
      [["endpoint.visible"]]
      [["endpoint.vsync"]]
      [[bobol_host_factory_instance_set_title]]
      [[bobol_host_factory_instance_set_visible]]
      [[bobol_host_factory_instance_set_vsync]])
    string(FIND "${_endpoint}" "${_needle}" _endpoint_property_idx)
    if(_endpoint_property_idx EQUAL -1)
      _bobol_guard_fail(
	"display endpoint lost typed toplevel host policy (${_needle})")
    endif()
  endforeach()

  _bobol_guard_read_rel(_tk_host "src/libtclcad/tkobol/tk-obol-host.cpp")
  foreach(_needle
      [[factory.set_title = tk_factory_set_title]]
      [[factory.set_visible = tk_factory_set_visible]]
      [[factory.set_vsync = tk_factory_set_vsync]]
      [[BOBOL_HOST_CAP_PRESENT_VSYNC]])
    string(FIND "${_tk_host}" "${_needle}" _tk_host_property_idx)
    if(_tk_host_property_idx EQUAL -1)
      _bobol_guard_fail(
	"Tk endpoint host lost typed toplevel policy (${_needle})")
    endif()
  endforeach()

  _bobol_guard_read_rel(_qg_scene_sync "src/libqtcad/QgSceneSync.cpp")
  _bobol_guard_forbid_regexes("src/libqtcad/QgSceneSync.cpp"
    "${_qg_scene_sync}"
    [[ged_draw_obol_scene_controller_owned]])

  _bobol_guard_read_rel(_qg_quad "src/libqtcad/QgQuadView.cpp")
  string(FIND "${_qg_quad}" [[ged_view_context_obol_endpoint_set]]
    _qg_quad_endpoint_idx)
  if(_qg_quad_endpoint_idx EQUAL -1)
    _bobol_guard_fail(
      "QgQuadView no longer associates visible endpoints with GED view records")
  endif()

  _bobol_guard_read_rel(_qged_fbserv "src/qged/fbserv.cpp")
  foreach(_needle
      [[displayEndpoint()]]
      [[ged_view_context_obol_endpoint_set]]
      [[ged_view_framebuffer_backend_ensure]]
      [[qged_fbserv_release_ged_handlers]])
    string(FIND "${_qged_fbserv}" "${_needle}" _qged_fb_endpoint_idx)
    if(_qged_fb_endpoint_idx EQUAL -1)
      _bobol_guard_fail(
	"qged fbserv no longer binds its visible endpoint through libged (${_needle})")
    endif()
  endforeach()
  _bobol_guard_forbid_regexes("src/qged/fbserv.cpp" "${_qged_fbserv}"
    [[(^|[^A-Za-z0-9_])qdm_]]
    [[QgObolWindowHost[ \t\r\n]+host]]
    [[bobol_display_endpoint_host]]
    [[ged_draw_obol_framebuffer_backend_install_for_view]])

  _bobol_guard_read_rel(_ged_fbserv "src/libged/obol_fbserv.cpp")
  foreach(_needle
      [[ged_view_context_obol_endpoint_get]]
      [[present_on_flush = 1]])
    string(FIND "${_ged_fbserv}" "${_needle}" _ged_fb_endpoint_host_idx)
    if(_ged_fb_endpoint_host_idx EQUAL -1)
      _bobol_guard_fail(
	"libged no longer binds framebuffer publication to its active endpoint (${_needle})")
    endif()
  endforeach()

  _bobol_guard_read_rel(_qged_app "src/qged/QgEdApp.cpp")
  string(FIND "${_qged_app}" [[qtcad_obol_host_factories_register]] _qged_factory_idx)
  if(_qged_factory_idx EQUAL -1)
    _bobol_guard_fail("qged no longer registers libqtcad Obol host factories")
  endif()
  string(FIND "${_qged_app}" [[qged_fbserv_release_ged_handlers]]
    _qged_fbserv_release_idx)
  if(_qged_fbserv_release_idx EQUAL -1)
    _bobol_guard_fail(
      "qged no longer clears its Qt fbserv transport target before GED teardown")
  endif()

  _bobol_guard_read_rel(_draw_obol "src/libged/draw_obol.cpp")
  string(REGEX MATCHALL [[ged_draw_obol_scene_controller_owned\(gedp\)]]
    _ownership_checks "${_draw_obol}")
  list(LENGTH _ownership_checks _ownership_check_count)
  if(NOT _ownership_check_count EQUAL 0)
    _bobol_guard_fail(
      "src/libged/draw_obol.cpp uses controller ownership as a scene capability gate")
  endif()

  _bobol_guard_read_rel(_scene_api "src/libged/ged_scene_api.cpp")
  foreach(_needle
      [[ged_scene_backend_apply_private]]
      [[scene_backend_ops->apply]]
      [[ged_scene_backend_snapshot_private]]
      [[scene_backend_ops->snapshot]]
      [[ged_scene_backend_set_private]])
    string(FIND "${_scene_api}" "${_needle}" _scene_backend_idx)
    if(_scene_backend_idx EQUAL -1)
      _bobol_guard_fail(
	"libged lost its single semantic scene-backend seam (${_needle})")
    endif()
  endforeach()
endfunction()

function(_bobol_guard_check_imgstream_display_host_ownership)
  set(_display_host_files
    include/imgstream/fb_compat.h
    src/libimgstream/fb_compat.c
    include/BObol/BWindowHost.h
    src/libBObol/window_host.cpp)
  set(_combined)
  foreach(_rel IN LISTS _display_host_files)
    _bobol_guard_read_rel(_contents "${_rel}")
    string(APPEND _combined "${_contents}")
    foreach(_retired
        imgstream_fb_display_host_set
        imgstream_fb_display_host_clear
        bobol_window_host_register_display_host
        bobol_window_host_unregister_display_host)
      string(FIND "${_contents}" "${_retired}" _retired_idx)
      if(NOT _retired_idx EQUAL -1)
        _bobol_guard_fail(
          "display framebuffer code reintroduced process-global host API ${_retired}")
      endif()
    endforeach()
  endforeach()

  foreach(_needle
      imgstream_fb_open_display
      bobol_window_host_open_display_framebuffer)
    string(FIND "${_combined}" "${_needle}" _explicit_idx)
    if(_explicit_idx EQUAL -1)
      _bobol_guard_fail(
        "display framebuffer code no longer exposes explicit per-open host API ${_needle}")
    endif()
  endforeach()
endfunction()

function(_bobol_guard_check_imgstream_remote_client_ownership)
  _bobol_guard_read_rel(_fb_header "include/imgstream/fb_compat.h")
  _bobol_guard_read_rel(_fb_compat "src/libimgstream/fb_compat.c")
  _bobol_guard_read_rel(_fb_remote "src/libimgstream/fb_remote.c")
  _bobol_guard_read_rel(_stress "regress/fbserv/fbserv_stress.cpp")
  _bobol_guard_read_rel(_stress_cmake "regress/fbserv/CMakeLists.txt")

  foreach(_needle
      [[struct imgstream_fb_remote_options]]
      [[imgstream_fb_open_remote]])
    string(FIND "${_fb_header}" "${_needle}" _header_idx)
    if(_header_idx EQUAL -1)
      _bobol_guard_fail(
	"libimgstream no longer exposes explicit per-open remote framebuffer state (${_needle})")
    endif()
  endforeach()
  foreach(_needle
      [[imgstream_fb_remote_open]]
      [[FBSERV_MSG_FBAUTH]]
      [[FBSERV_MSG_FBOPEN]]
      [[pkg_connect_addr]]
      [[pkg_open]])
    string(FIND "${_fb_remote}" "${_needle}" _remote_idx)
    if(_remote_idx EQUAL -1)
      _bobol_guard_fail(
	"libimgstream remote framebuffer client lost required protocol ownership (${_needle})")
    endif()
  endforeach()
  foreach(_needle
      [[imgstream_fb_open_remote]]
      [[--ipc]]
      [[imgstream_fb_writerect]]
      [[imgstream_fb_readrect]])
    string(FIND "${_stress}" "${_needle}" _stress_idx)
    if(_stress_idx EQUAL -1)
      _bobol_guard_fail(
	"fbserv stress coverage no longer exercises imgstream remote operation ${_needle}")
    endif()
  endforeach()
  string(FIND "${_stress_cmake}" [[libimgstream]] _stress_link_idx)
  if(_stress_link_idx EQUAL -1)
    _bobol_guard_fail(
      "fbserv stress coverage no longer links the owning libimgstream client")
  endif()
  _bobol_guard_forbid_regexes("src/libimgstream/fb_remote.c"
    "${_fb_remote}"
    [[#[ 	]*include[ 	]*[<"]dm]])
  _bobol_guard_forbid_regexes("regress/fbserv/fbserv_stress.cpp"
    "${_stress}"
    [[(^|[^A-Za-z0-9_])MSG_FB(OPEN|AUTH)([^A-Za-z0-9_]|$)]])
endfunction()

function(_bobol_guard_check_standalone_fb_ownership)
  _bobol_guard_read_rel(_fb_cmake "src/fb/CMakeLists.txt")
  _bobol_guard_forbid_regexes("src/fb/CMakeLists.txt" "${_fb_cmake}"
    [[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])dm_plugins([^A-Za-z0-9_]|$)]])

  _bobol_guard_collect(_fb_files src/fb)
  foreach(_file IN LISTS _fb_files)
    if(NOT "${_file}" MATCHES [[\.c$]])
      continue()
    endif()
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    file(READ "${_file}" _contents)
    _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[#[ \t]*include[ \t]*[<"]dm(\.h|/)]]
      [[struct[ \t\r\n]+fb([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])FB_NULL([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])RGBPIXEL_NULL([^A-Za-z0-9_]|$)]])
  endforeach()

  foreach(_rel
      src/sig/ddisp.c
      src/util/bwhist.c
      src/util/pixhist.c
      src/util/pixhist3d.c)
    _bobol_guard_read_rel(_contents "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[#[ \t]*include[ \t]*[<"]dm(\.h|/)]]
      [[struct[ \t\r\n]+fb([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])FB_NULL([^A-Za-z0-9_]|$)]])
    string(FIND "${_contents}" [[imgstream_fb_open]] _open_idx)
    if(_open_idx EQUAL -1)
      _bobol_guard_fail("${_rel} no longer opens its output through imgstream")
    endif()
  endforeach()

  _bobol_guard_read_rel(_sig_cmake "src/sig/CMakeLists.txt")
  _bobol_guard_read_rel(_util_cmake "src/util/CMakeLists.txt")
  foreach(_row
      [[brlcad_addexec(ddisp ddisp.c "libimgstream;libbu"]])
    string(FIND "${_sig_cmake}" "${_row}" _row_idx)
    if(_row_idx EQUAL -1)
      _bobol_guard_fail("signal display target no longer has its imgstream-only dependency row")
    endif()
  endforeach()
  foreach(_target bwhist pixhist pixhist3d)
    string(REGEX MATCH "brlcad_addexec\\(${_target}[^\n]*libimgstream;libbu"
      _row "${_util_cmake}")
    if(NOT _row)
      _bobol_guard_fail("${_target} no longer has its imgstream-only dependency row")
    endif()
  endforeach()

  set(_image_utility_sources
    bw-a.c
    bw-png.c
    double-asc.c
    halftone.c
    imgdims.c
    pix-png.c
    pix-ppm.c
    pix-spm.c
    pixbgstrip.c
    pixborder.c
    pixelswap.c
    pixhalve.c
    pixmorph.c
    wavelet.c)
  foreach(_source IN LISTS _image_utility_sources)
    set(_rel "src/util/${_source}")
    _bobol_guard_read_rel(_contents "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[#[ \t]*include[ \t]*[<"]dm(\.h|/)]])
  endforeach()

  set(_image_utility_targets
    bw-a
    bw-png
    double-asc
    halftone
    imgdims
    pix-png
    pix-ppm
    pixbgstrip
    pixborder
    pixelswap
    pixhalve
    pixmorph
    wavelet
    pix-spm)
  foreach(_target IN LISTS _image_utility_targets)
    string(REGEX MATCH "brlcad_addexec\\(${_target}[^\n]*libdm"
      _legacy_link "${_util_cmake}")
    string(REGEX MATCH "add_target_deps\\(${_target}[^\n]*dm_plugins"
      _legacy_plugins "${_util_cmake}")
    if(_legacy_link OR _legacy_plugins)
      _bobol_guard_fail("${_target} regained a libdm or DM plugin dependency")
    endif()
  endforeach()

  _bobol_guard_collect(_util_files src/util)
  foreach(_file IN LISTS _util_files)
    if(NOT "${_file}" MATCHES [[\.c$]])
      continue()
    endif()
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    file(READ "${_file}" _contents)
    _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[#[ \t]*include[ \t]*[<"]dm(\.h|/)]])
  endforeach()

  _bobol_guard_forbid_regexes("src/util/CMakeLists.txt"
    "${_util_cmake}"
    [[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])dm_plugins([^A-Za-z0-9_]|$)]])
  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/util/plot3-dm.c")
    _bobol_guard_fail("retired plot3-dm libdm sample was restored")
  endif()
endfunction()

function(_bobol_guard_check_rt_imgstream_ownership)
  _bobol_guard_collect(_rt_files src/rt src/remrt src/art)
  foreach(_file IN LISTS _rt_files)
    if(NOT "${_file}" MATCHES [[(CMakeLists\.txt|\.(c|cpp|h)$)]])
      continue()
    endif()
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    file(READ "${_file}" _contents)
    _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[#[ \t]*include[ \t]*[<"]dm(\.h|/)]]
      [[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])dm_plugins([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])FB_NULL([^A-Za-z0-9_]|$)]]
      [[struct[ \t\r\n]+fb([^A-Za-z0-9_]|$)]])
  endforeach()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/rt/libfb-dummy.c")
    _bobol_guard_fail("src/rt/libfb-dummy.c restored the retired RT framebuffer shim")
  endif()

  _bobol_guard_read_rel(_rt_ext "src/rt/ext.h")
  _bobol_guard_read_rel(_rt_main "src/rt/main.c")
  _bobol_guard_read_rel(_remrt "src/remrt/remrt.c")
  _bobol_guard_read_rel(_art_tile "src/art/tile.cpp")
  string(FIND "${_rt_ext}" [[imgstream_fb_t *fbp]] _rt_type_idx)
  string(FIND "${_rt_main}" [[imgstream_fb_open]] _rt_open_idx)
  string(FIND "${_remrt}" [[imgstream_fb_open]] _remrt_open_idx)
  string(FIND "${_art_tile}" [[imgstream_fb_writerect]] _art_rect_idx)
  if(_rt_type_idx EQUAL -1 OR _rt_open_idx EQUAL -1 OR
      _remrt_open_idx EQUAL -1 OR _art_rect_idx EQUAL -1)
    _bobol_guard_fail("RT/imgstream ownership contract is incomplete")
  endif()
endfunction()

function(_bobol_guard_check_standalone_display_provider_boundary)
  _bobol_guard_read_rel(_rt_main "src/rt/main.c")
  _bobol_guard_read_rel(_rt_cmake "src/rt/CMakeLists.txt")
  _bobol_guard_read_rel(_fbserv_main "src/fbserv/fbserv.c")
  _bobol_guard_read_rel(_fbserv_cmake "src/fbserv/CMakeLists.txt")

  foreach(_source_name "src/rt/main.c" "src/fbserv/fbserv.c")
    if(_source_name STREQUAL "src/rt/main.c")
      set(_source_contents "${_rt_main}")
    else()
      set(_source_contents "${_fbserv_main}")
    endif()
    _bobol_guard_forbid_regexes("${_source_name}" "${_source_contents}"
      [[#[ \t]*include[ \t]*[<"](tcl|tk)(\.h|/)]]
      [[(^|[^A-Za-z0-9_])(Tcl|Tk)_[A-Za-z0-9_]+]])
    string(FIND "${_source_contents}" [[BObol/BDisplaySession.h]]
      _header_idx)
    string(FIND "${_source_contents}" [[qtcad/display_provider.h]]
      _qtcad_header_idx)
    string(FIND "${_source_contents}" [[tclcad/setup.h]] _tclcad_header_idx)
    string(FIND "${_source_contents}" [[qtcad_obol_display_provider_register]]
      _qtcad_register_idx)
    string(FIND "${_source_contents}" [[tclcad_obol_display_provider_register]]
      _provider_register_idx)
    string(FIND "${_source_contents}" [[bobol_display_session_open]]
      _open_idx)
    if(_header_idx EQUAL -1 OR _qtcad_header_idx EQUAL -1 OR
       _tclcad_header_idx EQUAL -1 OR _qtcad_register_idx EQUAL -1 OR
       _provider_register_idx EQUAL -1 OR _open_idx EQUAL -1)
      _bobol_guard_fail(
        "${_source_name} no longer selects the Qt-first/TclCAD-fallback Obol provider boundary")
    endif()
  endforeach()

  foreach(_cmake_name "src/rt/CMakeLists.txt" "src/fbserv/CMakeLists.txt")
    if(_cmake_name STREQUAL "src/rt/CMakeLists.txt")
      set(_cmake_contents "${_rt_cmake}")
    else()
      set(_cmake_contents "${_fbserv_cmake}")
    endif()
    _bobol_guard_forbid_regexes("${_cmake_name}" "${_cmake_contents}"
      [[(^|[^A-Za-z0-9_])libtkobolprovider([^A-Za-z0-9_]|$)]])
    string(FIND "${_cmake_contents}" [[libqtcad]] _qtcad_target_idx)
    string(FIND "${_cmake_contents}" [[libtclcad]] _tclcad_target_idx)
    string(FIND "${_cmake_contents}" [[HAVE_QTCAD_OBOL_DISPLAY_PROVIDER]]
      _qtcad_definition_idx)
    string(FIND "${_cmake_contents}" [[HAVE_TCLCAD_OBOL_DISPLAY_PROVIDER]]
      _tclcad_definition_idx)
    if(_qtcad_target_idx EQUAL -1 OR _tclcad_target_idx EQUAL -1 OR
       _qtcad_definition_idx EQUAL -1 OR _tclcad_definition_idx EQUAL -1)
      _bobol_guard_fail(
        "${_cmake_name} bypasses the Qt-first/TclCAD-fallback Obol provider targets")
    endif()
  endforeach()
endfunction()

function(_bobol_guard_check_rtwizard_host_boundary)
  _bobol_guard_read_rel(_rtwizard_cmake "src/rtwizard/CMakeLists.txt")
  _bobol_guard_forbid_regexes("src/rtwizard/CMakeLists.txt"
    "${_rtwizard_cmake}"
    [[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])dm_plugins([^A-Za-z0-9_]|$)]])

  _bobol_guard_read_rel(_rtwizard_fb_page
    "src/tclscripts/rtwizard/lib/FbPage.itk")
  _bobol_guard_forbid_regexes("src/tclscripts/rtwizard/lib/FbPage.itk"
    "${_rtwizard_fb_page}"
    [[/dev/(oglsp|wglsp)]])
  string(REGEX MATCH [[exec.*fbserv.*-A.*-F[ \t]+\$fbType]]
    _rtwizard_obol_fbserv "${_rtwizard_fb_page}")
  string(REGEX MATCH [[set[ \t]+fbType[ \t]+"/dev/ogl"]]
    _rtwizard_system_gl "${_rtwizard_fb_page}")
  string(REGEX MATCH [[set[ \t]+fbType[ \t]+"/dev/wgl"]]
    _rtwizard_windows_gl "${_rtwizard_fb_page}")
  string(FIND "${_rtwizard_fb_page}" [[FBSERV_TOKEN]]
    _rtwizard_token_idx)
  if(NOT _rtwizard_obol_fbserv OR NOT _rtwizard_system_gl OR
     NOT _rtwizard_windows_gl OR _rtwizard_token_idx EQUAL -1)
    _bobol_guard_fail(
      "rtwizard fixed-size preview must use a token-isolated Obol fbserv surface")
  endif()
endfunction()

function(_bobol_guard_check_ged_framebuffer_image_ownership)
  foreach(_dir fb2pix pix2fb png2fb fbclear overlay)
    _bobol_guard_read_rel(_cmake "src/libged/${_dir}/CMakeLists.txt")
    _bobol_guard_forbid_regexes("src/libged/${_dir}/CMakeLists.txt"
      "${_cmake}"
      [[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])dm_plugins([^A-Za-z0-9_]|$)]])

    _bobol_guard_collect(_files "src/libged/${_dir}")
    foreach(_file IN LISTS _files)
      if(NOT "${_file}" MATCHES [[\.(c|cpp|h)$]])
	continue()
      endif()
      file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
      file(READ "${_file}" _contents)
      _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
	[[#[ \t]*include[ \t]*[<"]dm(\.h|/)]]
	[[struct[ \t\r\n]+dm([^A-Za-z0-9_]|$)]]
	[[struct[ \t\r\n]+fb([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])dm_(draw|open|close|get|set|put|refresh)[A-Za-z0-9_]*[ \t\r\n]*\(]]
	[[(^|[^A-Za-z0-9_])fb_(read|write|clear|refresh|get)[A-Za-z0-9_]*[ \t\r\n]*\(]])
    endforeach()
  endforeach()

  _bobol_guard_read_rel(_bridge "src/libged/obol_fbserv.cpp")
  string(FIND "${_bridge}" [[ged_view_framebuffer_apply]]
    _apply_idx)
  if(_apply_idx EQUAL -1)
    _bobol_guard_fail(
      "libged no longer exposes serialized Obol framebuffer command access")
  endif()
endfunction()

function(_bobol_guard_check_qged_edit_preview_policy)
  _bobol_guard_collect(_files
    src/qged)
  foreach(_file IN LISTS _files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if(NOT "${_rel}" MATCHES "${_bobol_guard_extensions}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    string(REGEX MATCH [[QgObolEditPreviewPrivate\.h|(^|[^A-Za-z0-9_])qg_obol_edit_preview_[A-Za-z0-9_]*[ \t\r\n]*\(]]
      _hit "${_contents}")
    if(_hit)
      _bobol_guard_fail(
	"${_rel} uses the retired qtcad direct-widget edit-preview helper; route qged edit previews through GED view features")
    endif()
  endforeach()
endfunction()

function(_bobol_guard_check_tclcad_obol_readback_bridge)
  _bobol_guard_read_rel(_commands "src/libtclcad/commands.c")
  foreach(_needle
      [[#include "ged/display.h"]]
      [[ged_view_display_image_capture]]
      [[tclcad_commands_capture_rgb]]
      [[ged_view_context_bv_const]])
    string(FIND "${_commands}" "${_needle}" _idx)
    if(_idx EQUAL -1)
      _bobol_guard_fail(
	"src/libtclcad/commands.c no longer routes image export through endpoint-native capture (${_needle})")
    endif()
  endforeach()
  _bobol_guard_forbid_regexes("src/libtclcad/commands.c" "${_commands}"
    [[(^|[^A-Za-z0-9_])ged_draw_obol_view_display_image[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])dm_get_display_image[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])bobol_display_endpoint_capture[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])dm_(get|set)_(bg|light|transparency|zbuffer|zclip)[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])to_bounds[ \t\r\n]*\(]])
  foreach(_needle
      [[tclcad_commands_endpoint_bool_get]]
      [[tclcad_commands_endpoint_bool_set]]
      [[controller.background.bottom]]
      [[renderer.lighting]]
      [[renderer.transparency]]
      [[renderer.depth_test]]
      [[view.zclip]])
    string(FIND "${_commands}" "${_needle}" _property_idx)
    if(_property_idx EQUAL -1)
      _bobol_guard_fail(
	"TclCAD commands lost endpoint-native display policy (${_needle})")
    endif()
  endforeach()
  _bobol_guard_read_rel(_ged_tcl "src/tclscripts/lib/Ged.tcl")
  _bobol_guard_read_rel(_archer_core "src/tclscripts/archer/ArcherCore.tcl")
  _bobol_guard_forbid_regexes("src/tclscripts/lib/Ged.tcl" "${_ged_tcl}"
    [[(^|[^A-Za-z0-9_])(bounds|bounds_all)[ \t]*\{]])
  _bobol_guard_forbid_regexes("src/tclscripts/archer/ArcherCore.tcl" "${_archer_core}"
    [[bounds_all]]
    [[mZClip(Front|Back)]]
    [[updateZClipPlanes]])
  _bobol_guard_read_rel(_archer "src/tclscripts/archer/Archer.tcl")
  _bobol_guard_read_rel(_bot_edit "src/tclscripts/archer/BotEditFrame.tcl")
  _bobol_guard_forbid_regexes("src/tclscripts/archer/Archer.tcl" "${_archer}"
    [[mZClip(Front|Back)]]
    [[getZClipState]]
    [[zclip(Front|Back)]])
  _bobol_guard_forbid_regexes("src/tclscripts/archer/BotEditFrame.tcl" "${_bot_edit}"
    [[getZClipState]])

  _bobol_guard_read_rel(_controller "src/libBObol/view_controller.cpp")
  foreach(_needle
      [[SoClipPlane]]
      [[BObolClipMinimum]]
      [[BObolClipMaximum]])
    string(FIND "${_controller}" "${_needle}" _clip_idx)
    if(_clip_idx EQUAL -1)
      _bobol_guard_fail(
	"libBObol lost retained clipping-plane policy (${_needle})")
    endif()
  endforeach()
  foreach(_needle
      [[$mGed mouse_ray]]
      [[$mGed mouse_pick_detail]]
      [[set mLastMouseRayStart]]
      [[set mLastMouseRayTarget]])
    string(FIND "${_ged_tcl}" "${_needle}" _ray_idx)
    if(_ray_idx EQUAL -1)
      _bobol_guard_fail(
	"TclCAD mouse picking lost display-independent model rays (${_needle})")
    endif()
  endforeach()
  _bobol_guard_forbid_regexes("src/tclscripts/lib/Ged.tcl" "${_ged_tcl}"
    [[shoot_ray[ \t]+obj_ray]])
endfunction()

function(_bobol_guard_check_fbserv_legacy_framebuffer_retired)
  _bobol_guard_collect(_files
    src/libged
    src/qged
    src/libtclcad)
  foreach(_file IN LISTS _files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if(NOT "${_rel}" MATCHES "${_bobol_guard_extensions}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    string(REGEX MATCH
      [[(^|[^A-Za-z0-9_])(fbs_fbp|fbs_set_legacy_framebuffer|fbs_legacy_framebuffer|fbserv_set_legacy_framebuffer|fbserv_legacy_framebuffer)([^A-Za-z0-9_]|$)]]
      _hit "${_contents}")
    if(_hit)
      _bobol_guard_fail(
	"${_rel} restored the retired legacy framebuffer pointer path")
    endif()
  endforeach()
endfunction()

function(_bobol_guard_check_fbserv_transport_quarantine)
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

  _bobol_guard_collect(_files
    src/libged
    src/gtools/gsh
    src/qged/fbserv.cpp
    src/libtclcad/fbserv.c)
  foreach(_file IN LISTS _files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if(NOT "${_rel}" MATCHES "${_bobol_guard_extensions}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    foreach(_pat IN LISTS _transport_fields)
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_bobol_guard_fail(
	  "${_rel} reaches into fbserv transport internals; use fbs transport/accessor helpers instead")
      endif()
    endforeach()
  endforeach()
endfunction()

function(_bobol_guard_check_fbserv_transport_setup_helpers)
  _bobol_guard_collect(_files
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
	_bobol_guard_fail(
	  "${_rel} assigns fbserv transport callbacks directly; use fbs_set_transport or tclcad_fbserv_set_transport")
      endif()
    endforeach()
  endforeach()

  _bobol_guard_read_rel(_qged "src/qged/fbserv.cpp")
  string(FIND "${_qged}" [[fbs_set_transport]] _qged_transport_idx)
  if(_qged_transport_idx EQUAL -1)
    _bobol_guard_fail(
      "src/qged/fbserv.cpp no longer installs fbserv transport callbacks through fbs_set_transport")
  endif()

  _bobol_guard_read_rel(_tcl_fbserv "src/libtclcad/fbserv.c")
  string(FIND "${_tcl_fbserv}" [[tclcad_fbserv_set_transport]] _tcl_fbserv_transport_idx)
  if(_tcl_fbserv_transport_idx EQUAL -1)
    _bobol_guard_fail(
      "src/libtclcad/fbserv.c no longer centralizes the fallback fbserv transport through tclcad_fbserv_set_transport")
  endif()
endfunction()

function(_bobol_guard_check_qged_fbserv_backend_access)
  _bobol_guard_collect(_files
    src/qged)
  foreach(_file IN LISTS _files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if(NOT "${_rel}" MATCHES "${_bobol_guard_extensions}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])fbs_fb_(ops|ctx)([^A-Za-z0-9_]|$)]]
      _hit "${_contents}")
    if(_hit)
      _bobol_guard_fail(
	"${_rel} reaches into fbserv backend internals; use fbs_framebuffer_* helpers instead")
    endif()
  endforeach()
endfunction()

function(_bobol_guard_check_fbserv_backend_contract)
  _bobol_guard_read_rel(_imgstream_fbserv "include/imgstream/fbserv.h")
  _bobol_guard_read_rel(_imgstream_auth "src/libimgstream/fbserv_auth.c")
  _bobol_guard_read_rel(_imgstream_framebuffer "src/libimgstream/fbserv_framebuffer.c")
  _bobol_guard_read_rel(_imgstream_pkg "src/libimgstream/fbserv_pkg.c")
  _bobol_guard_read_rel(_imgstream_server "src/libimgstream/fbserv_server.c")
  _bobol_guard_read_rel(_imgstream_state "src/libimgstream/fbserv_state.c")
  _bobol_guard_read_rel(_standalone_fbserv "src/fbserv/fbserv.c")
  _bobol_guard_read_rel(_standalone_fbserv_cmake "src/fbserv/CMakeLists.txt")
  _bobol_guard_read_rel(_protocol_note
    "doc/notes/obol_framebuffer_protocol.txt")

  string(FIND "${_imgstream_fbserv}" [[struct fbserv_fb_ops]] _ops_idx)
  if(_ops_idx EQUAL -1)
    _bobol_guard_fail(
      "include/imgstream/fbserv.h no longer owns the fbserv framebuffer backend operation table")
  endif()
  string(FIND "${_imgstream_fbserv}" [[FBSERV_MSG_FBOPEN]] _msg_idx)
  if(_msg_idx EQUAL -1)
    _bobol_guard_fail(
      "include/imgstream/fbserv.h no longer owns the fbserv wire-protocol message IDs")
  endif()
  string(FIND "${_imgstream_fbserv}" [[struct fbserv_obj]] _obj_idx)
  string(FIND "${_imgstream_fbserv}" [[struct fbserv_transport_ops]] _transport_ops_idx)
  if(_obj_idx EQUAL -1 OR _transport_ops_idx EQUAL -1)
    _bobol_guard_fail(
      "include/imgstream/fbserv.h no longer owns the transitional fbserv object/transport type contract")
  endif()
  string(FIND "${_imgstream_state}" [[fbserv_obj_init]] _state_impl_idx)
  string(FIND "${_imgstream_server}" [[fbs_init]] _server_state_wrapper_idx)
  if(_state_impl_idx EQUAL -1 OR _server_state_wrapper_idx EQUAL -1)
    _bobol_guard_fail(
      "libimgstream no longer owns fbserv state helpers and complete fbs_* lifecycle wiring")
  endif()
  string(FIND "${_imgstream_pkg}" [[fbserv_pkg_switch_init]] _pkg_switch_impl_idx)
  string(FIND "${_imgstream_server}" [[fbserv_pkg_switch_init]] _server_pkg_switch_wrapper_idx)
  if(_pkg_switch_impl_idx EQUAL -1 OR _server_pkg_switch_wrapper_idx EQUAL -1)
    _bobol_guard_fail(
      "libimgstream no longer owns the fbserv PKG message switch and fbs_* wiring")
  endif()
  string(FIND "${_imgstream_framebuffer}" [[fbserv_backend_info]] _backend_dispatch_idx)
  string(FIND "${_imgstream_server}" [[fbserv_backend_info]] _server_backend_dispatch_idx)
  if(_backend_dispatch_idx EQUAL -1 OR _server_backend_dispatch_idx EQUAL -1)
    _bobol_guard_fail(
      "libimgstream no longer owns fbserv backend-operation dispatch and server wiring")
  endif()

  foreach(_needle
      [[imgstream_fbserv_set_framebuffer]]
      [[imgstream_fb_open]])
    string(FIND "${_standalone_fbserv}" "${_needle}" _standalone_native_idx)
    if(_standalone_native_idx EQUAL -1)
      _bobol_guard_fail(
	"standalone fbserv no longer uses its imgstream-native framebuffer (${_needle})")
    endif()
  endforeach()
  _bobol_guard_forbid_regexes("src/fbserv/fbserv.c"
    "${_standalone_fbserv}"
    [[#[ 	]*include[ 	]*[<"]dm]]
    [[(^|[^A-Za-z0-9_])struct[ 	]+fb([^A-Za-z0-9_]|$)]]
    [[dm_fbserv_set_framebuffer]])
  _bobol_guard_forbid_regexes("src/fbserv/CMakeLists.txt"
    "${_standalone_fbserv_cmake}"
    [[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])dm_plugins([^A-Za-z0-9_]|$)]])
  foreach(_entry
      [[fbs_open(]]
      [[fbs_close(]]
      [[fbs_new_client(]]
      [[fbs_open_ipc(]]
      [[fbs_process_client_bytes(]]
      [[fbs_drain_client(]])
    string(FIND "${_imgstream_server}" "${_entry}" _server_lifecycle_idx)
    if(_server_lifecycle_idx EQUAL -1)
      _bobol_guard_fail(
	"src/libimgstream/fbserv_server.c no longer owns required lifecycle entry ${_entry}")
    endif()
  endforeach()
  foreach(_retired
      include/dm.h
      include/dm/fbserv.h
      src/libdm/CMakeLists.txt
      src/libdm/fbserv.c
      src/libdm/fbserv_legacy.c)
    if(EXISTS "${BRLCAD_SOURCE_DIR}/${_retired}")
      _bobol_guard_fail("${_retired} restored retired framebuffer compatibility ownership")
    endif()
  endforeach()
  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/fbserv/server.c")
    _bobol_guard_fail(
      "src/fbserv/server.c reintroduced a duplicate standalone framebuffer protocol server")
  endif()
  _bobol_guard_read_rel(_mged_fbserv "src/mged/fbserv.c")
  foreach(_token [[pkg_process(]] [[pkg_switch]] [[MSG_FBOPEN]])
    string(FIND "${_mged_fbserv}" "${_token}" _mged_protocol_idx)
    if(NOT _mged_protocol_idx EQUAL -1)
      _bobol_guard_fail(
	"src/mged/fbserv.c reintroduced private framebuffer protocol processing (${_token})")
    endif()
  endforeach()
  string(FIND "${_imgstream_fbserv}" [[fbserv_generate_token]] _auth_decl_idx)
  string(FIND "${_imgstream_auth}" [[fbserv_verify_token]] _auth_impl_idx)
  if(_auth_decl_idx EQUAL -1 OR _auth_impl_idx EQUAL -1)
    _bobol_guard_fail(
      "libimgstream no longer owns fbserv session-token authentication helpers")
  endif()
  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/fbserv/auth.h")
    _bobol_guard_fail(
      "src/fbserv/auth.h reintroduced header-local fbserv auth helpers; use include/imgstream/fbserv.h")
  endif()
  string(FIND "${_protocol_note}" [[include/imgstream/fbserv.h]] _note_idx)
  if(_note_idx EQUAL -1)
    _bobol_guard_fail(
      "doc/notes/obol_framebuffer_protocol.txt no longer documents imgstream fbserv backend/protocol contract ownership")
  endif()
endfunction()

function(_bobol_guard_check_tkobol_host_ownership)
  _bobol_guard_read_rel(_tclcad_cmake "src/libtclcad/CMakeLists.txt")
  _bobol_guard_read_rel(_mged_cmake "src/mged/CMakeLists.txt")
  _bobol_guard_read_rel(_mged_main "src/mged/mged.c")
  _bobol_guard_read_rel(_mged_attach "src/mged/attach.c")
  _bobol_guard_read_rel(_tclcad_tkobol_init "src/libtclcad/tkobol_init.c")
  _bobol_guard_read_rel(_tkobol_host "src/libtclcad/tkobol/tk-obol-host.cpp")
  foreach(_needle
      [[tkobol/tk-obol-host.cpp]]
      [[tkobol/vendor/togl/togl.c]])
    string(FIND "${_tclcad_cmake}" "${_needle}" _host_source_idx)
    if(_host_source_idx EQUAL -1)
      _bobol_guard_fail(
	"src/libtclcad/CMakeLists.txt no longer owns required Tk Obol host source ${_needle}")
    endif()
  endforeach()
  foreach(_needle
      [[BrlcadTkObolHost_Init]]
      [[tclcad_obol_host_factories_register]])
    string(FIND "${_tclcad_tkobol_init}" "${_needle}" _host_init_idx)
    if(_host_init_idx EQUAL -1)
      _bobol_guard_fail(
	"src/libtclcad/tkobol_init.c no longer initializes the built-in Tk Obol host (${_needle})")
    endif()
  endforeach()
  _bobol_guard_forbid_regexes("src/libtclcad/CMakeLists.txt" "${_tclcad_cmake}"
    [[add_subdirectory\(tkobol\)]]
    [[tkobol/dm-tkobol\.cpp]]
    [[dm_plugin_library\(dm-tkobol]])
  _bobol_guard_forbid_regexes("src/libtclcad/tkobol_init.c" "${_tclcad_tkobol_init}"
    [[(^|[^A-Za-z0-9_])dm_register_backend[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])dm_bestXType_tcl([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])dm_validXType_tcl([^A-Za-z0-9_]|$)]])
  _bobol_guard_forbid_regexes("src/mged/CMakeLists.txt" "${_mged_cmake}"
    [[(^|[^A-Za-z0-9_])dm_plugins([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]])
  foreach(_helper [[ged_plugins]] [[ged_subprocesses]] [[ged_exec]])
    string(FIND "${_mged_cmake}" "${_helper}" _mged_helper_idx)
    if(_mged_helper_idx EQUAL -1)
      _bobol_guard_fail(
	"src/mged/CMakeLists.txt no longer builds runtime GED helper ${_helper}")
    endif()
  endforeach()
  _bobol_guard_forbid_regexes("src/mged/attach.c" "${_mged_attach}"
    [[(^|[^A-Za-z0-9_])dm_(make_current|set_win_bounds|processOptions)[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])dm_get_(dname|dm_name|dm_lname)[ \t\r\n]*\(]])
  _bobol_guard_forbid_regexes("src/mged/mged.c" "${_mged_main}"
    [[dm_init_msgs[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])dm_version[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])fb_version[ \t\r\n]*\(]])
  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libtclcad/tkobol/dm-tkobol.cpp")
    _bobol_guard_fail("dm-tkobol compatibility adapter was restored")
  endif()
  foreach(_needle [["tk-gl"]] [["tk-photo"]] [[Tk_PhotoPutBlock]]
      [[::brlcad::tkObol::host]])
    string(FIND "${_tkobol_host}" "${_needle}" _native_host_idx)
    if(_native_host_idx EQUAL -1)
      _bobol_guard_fail(
	"Tk-native Obol host lost required factory/presentation logic (${_needle})")
    endif()
  endforeach()
  _bobol_guard_read_rel(_tclcad_commands "src/libtclcad/commands.c")
  _bobol_guard_forbid_regexes("src/libtclcad/commands.c" "${_tclcad_commands}"
    [[(^|[^A-Za-z0-9_])dm_list_tcl([^A-Za-z0-9_]|$)]])
  foreach(_rel
      src/tclscripts/lib/Ged.tcl
      src/tclscripts/archer/ArcherCore.tcl
      src/tclscripts/mged/rt.tcl
      src/tclscripts/mged/openw.tcl)
    _bobol_guard_read_rel(_tclcad_script "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_tclcad_script}"
      [[(^|[^A-Za-z0-9_])dm_list[ \t\r\n]*\(]]
      [[(^|[^A-Za-z0-9_])dm_bestXType[ \t\r\n]*\(]]
      [[(^|[^A-Za-z0-9_])dm_validXType[ \t\r\n]*\(]])
  endforeach()
  _bobol_guard_read_rel(_tclcad_ged_script "src/tclscripts/lib/Ged.tcl")
  _bobol_guard_forbid_regexes("src/tclscripts/lib/Ged.tcl"
    "${_tclcad_ged_script}"
    [[itk_option[ \t]+define[ \t]+-type[ \t]+type[ \t]+Type]])
  _bobol_guard_read_rel(_archer_core_script
    "src/tclscripts/archer/ArcherCore.tcl")
  _bobol_guard_forbid_regexes("src/tclscripts/archer/ArcherCore.tcl"
    "${_archer_core_script}"
    [[(^|[^A-Za-z0-9_])mDisplayType([^A-Za-z0-9_]|$)]])
  foreach(_direct_source _mged_attach _tclcad_commands)
    foreach(_needle
        [[bobol_display_endpoint_create]]
        [[bobol_display_endpoint_render_engine_set]]
        [[ged_view_context_obol_endpoint_set]]
        [[bobol_display_endpoint_host_open]])
      string(FIND "${${_direct_source}}" "${_needle}" _direct_host_idx)
      if(_direct_host_idx EQUAL -1)
        _bobol_guard_fail(
          "direct Tk endpoint ownership lost ${_needle} in ${_direct_source}")
      endif()
    endforeach()
  endforeach()
  foreach(_retired
      [[tclcad_tkobol_controller]]
      [[tclcad_tkobol_endpoint]]
      [[ged_draw_obol_controller_attach_opaque(]])
    string(FIND "${_tclcad_commands}" "${_retired}" _retired_attach_idx)
    if(NOT _retired_attach_idx EQUAL -1)
      _bobol_guard_fail(
	"src/libtclcad/commands.c reintroduced direct controller attachment (${_retired})")
    endif()
  endforeach()
  _bobol_guard_read_rel(_gsh_cmake "src/gtools/gsh/CMakeLists.txt")
  _bobol_guard_forbid_regexes("src/gtools/gsh/CMakeLists.txt" "${_gsh_cmake}"
    [[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]])
  foreach(_helper [[ged_plugins]] [[ged_subprocesses]] [[ged_exec]])
    string(FIND "${_gsh_cmake}" "${_helper}" _gsh_helper_idx)
    if(_gsh_helper_idx EQUAL -1)
      _bobol_guard_fail(
	"src/gtools/gsh/CMakeLists.txt no longer builds runtime GED helper ${_helper}")
    endif()
  endforeach()
  _bobol_guard_read_rel(_gsh "src/gtools/gsh/gsh.cpp")
  string(FIND "${_gsh}" [[gsh_headless_endpoint_ensure]]
    _gsh_endpoint_idx)
  if(_gsh_endpoint_idx EQUAL -1)
    _bobol_guard_fail(
      "src/gtools/gsh/gsh.cpp no longer obtains drawing through its GED view endpoint")
  endif()
  _bobol_guard_forbid_regexes("src/gtools/gsh/gsh.cpp" "${_gsh}"
    [[new[ \t\r\n]+BObolViewController]]
    [[ged_draw_obol_controller_(attach|detach)]])
endfunction()

function(_bobol_guard_check_retained_export_ownership)
  _bobol_guard_read_rel(_mged_cmd "src/mged/cmd.cpp")
  foreach(_retired
      src/libdm/plot/CMakeLists.txt
      src/libdm/plot/dm-plot.c
      src/libdm/plot/dm-plot.h
      src/libdm/postscript/CMakeLists.txt
      src/libdm/postscript/dm-ps.c
      src/libdm/postscript/dm-ps.h)
    if(EXISTS "${BRLCAD_SOURCE_DIR}/${_retired}")
      _bobol_guard_fail("${_retired} restored a retired export display manager")
    endif()
  endforeach()
  _bobol_guard_forbid_regexes("src/mged/cmd.cpp" "${_mged_cmd}"
    [[mged_attach[ \t\r\n]*\([^\n]*postscript]])
  foreach(_needle [[ged_exec(s->gedp]] [[bu_argv_dup]])
    string(FIND "${_mged_cmd}" "${_needle}" _mged_ps_idx)
    if(_mged_ps_idx EQUAL -1)
      _bobol_guard_fail(
	"MGED postscript command no longer delegates to the retained GED exporter (${_needle})")
    endif()
  endforeach()
endfunction()

function(_bobol_guard_check_retired_text_dm)
  foreach(_retired
      src/libdm/txt/CMakeLists.txt
      src/libdm/txt/dm-txt.c)
    if(EXISTS "${BRLCAD_SOURCE_DIR}/${_retired}")
      _bobol_guard_fail("${_retired} restored the unreachable text display manager")
    endif()
  endforeach()
endfunction()

function(_bobol_guard_check_unreachable_drawing_shims)
  foreach(_retired
      include/qtcad/QgLegacyView.h
      src/libqtcad/QgLegacyView.cpp
      src/libqtcad/QgLegacyViewContext.h
      src/libqtcad/QgObolDatabaseSync.cpp
      src/libqtcad/QgObolDatabaseSyncPrivate.h)
    if(EXISTS "${BRLCAD_SOURCE_DIR}/${_retired}")
      _bobol_guard_fail("${_retired} restored the retired direct Qt database-source synchronization path")
    endif()
  endforeach()

  _bobol_guard_collect(_scene_files
    src/libqtcad
    src/qged)
  foreach(_file IN LISTS _scene_files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if(NOT "${_rel}" MATCHES "${_bobol_guard_extensions}" OR
       "${_rel}" MATCHES "/tests/")
      continue()
    endif()
    file(READ "${_file}" _contents)
    _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[QgObolDatabaseSync]]
      [[qg_obol_(sync|remove|clear)_database_sources]]
      [[(^|[^A-Za-z0-9_])(replace|remove|clear)DatabaseSource(Instance|s)[A-Za-z0-9_]*[ \t\r\n]*\(]])
  endforeach()

  foreach(_rel
      src/mged/scene_refresh.c
      src/libtclcad/view/draw.c
      src/libtclcad/view/refresh.c)
    _bobol_guard_read_rel(_contents "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[(^|[^A-Za-z0-9_])dm_draw_(objs|viewobjs)[ \t\r\n]*\(]])
  endforeach()
  _bobol_guard_read_rel(_mged_scene_refresh "src/mged/scene_refresh.c")
  _bobol_guard_forbid_regexes("src/mged/scene_refresh.c"
    "${_mged_scene_refresh}"
    [[(^|[^A-Za-z0-9_])dm_[A-Za-z0-9_]*[ \t\r\n]*\(]]
    [[ged_view_context_display_manager_set]])
  foreach(_needle
      [[mged_obol_scene_refresh]]
      [[ged_scene_redraw]])
    string(FIND "${_mged_scene_refresh}" "${_needle}" _scene_refresh_idx)
    if(_scene_refresh_idx EQUAL -1)
      _bobol_guard_fail(
	"MGED scene refresh lost retained transaction ownership (${_needle})")
    endif()
  endforeach()

  _bobol_guard_read_rel(_tclcad_refresh "src/libtclcad/view/refresh.c")
  _bobol_guard_forbid_regexes("src/libtclcad/view/refresh.c" "${_tclcad_refresh}"
    [[(^|[^A-Za-z0-9_])dm_draw_[A-Za-z0-9_]*[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])dm_set_(depth_mask|zbuffer)[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])fb_refresh[ \t\r\n]*\(]])
  foreach(_needle
      [[go_draw(draw_view_ctx)]]
      [[ged_view_framebuffer_present]]
      [[bobol_display_endpoint_view_sync]]
      [[bobol_display_endpoint_request_frame]])
    string(FIND "${_tclcad_refresh}" "${_needle}" _tclcad_refresh_idx)
    if(_tclcad_refresh_idx EQUAL -1)
      _bobol_guard_fail(
	"TclCAD refresh lost retained endpoint rendering (${_needle})")
    endif()
  endforeach()

  foreach(_rel
      src/libtclcad/commands.c
      src/libtclcad/mouse.c
      src/libtclcad/polygons.c
      src/libtclcad/wrapper.c
      src/libtclcad/view/axes.c)
    _bobol_guard_read_rel(_tclcad_dimensions "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_tclcad_dimensions}"
      [[(^|[^A-Za-z0-9_])dm_get_(width|height|aspect)[ \t\r\n]*\(]]
      [[(^|[^A-Za-z0-9_])dm_geometry_request[ \t\r\n]*\(]])
  endforeach()
  _bobol_guard_read_rel(_tclcad_fontsize "src/libtclcad/commands.c")
  _bobol_guard_forbid_regexes("src/libtclcad/commands.c"
    "${_tclcad_fontsize}"
    [[(^|[^A-Za-z0-9_])dm_(get|set)_fontsize[ \t\r\n]*\(]])

  _bobol_guard_read_rel(_mged_refresh "src/mged/mged.c")
  _bobol_guard_forbid_regexes("src/mged/mged.c" "${_mged_refresh}"
    [[(^|[^A-Za-z0-9_])(draw_grid|adcursor|draw_m_axes|draw_v_axes)[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])dm_hud_begin[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])dm_draw_point_2d[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])dm_(configure_win|draw_begin|draw_end|make_current|get_fontsize|get_zbuffer|set_zbuffer|get_stereo|get_native_repaint_pending|set_native_repaint_pending)[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])fb_refresh[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])paint_rect_area[ \t\r\n]*\(]])
  foreach(_needle
      [[mged_obol_framebuffer_composition_sync]]
      [[composition.framebuffer.mode]]
      [[mged_obol_faceplate_state_sync]]
      [[ged_view_faceplate_sync]]
      [[bobol_display_endpoint_view_sync]]
      [[bobol_display_endpoint_request_frame]])
    string(FIND "${_mged_refresh}" "${_needle}" _mged_refresh_idx)
    if(_mged_refresh_idx EQUAL -1)
      _bobol_guard_fail(
	"MGED refresh lost retained endpoint rendering (${_needle})")
    endif()
  endforeach()
  _bobol_guard_read_rel(_mged_repaint "src/mged/mged_display.h")
  _bobol_guard_forbid_regexes("src/mged/mged_display.h" "${_mged_repaint}"
    [[(^|[^A-Za-z0-9_])dm_set_native_repaint_pending[ \t\r\n]*\(]])
  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/mged/mged_dm.h")
    _bobol_guard_fail(
      "src/mged/mged_dm.h restored a retired display-manager state record")
  endif()
  foreach(_rel
      src/mged/mged_display.h
      src/mged/mged.h
      src/mged/attach.c
      src/mged/mged.c)
    _bobol_guard_read_rel(_mged_display_identity "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_mged_display_identity}"
      [[(^|[^A-Za-z0-9_])mged_dm([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])mged_curr_dm([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])active_dm_set([^A-Za-z0-9_]|$)]])
  endforeach()
  _bobol_guard_read_rel(_ged_dm_endpoint "src/libged/dm/dm.cpp")
  foreach(_needle
      [[host_width = bview && bv_width_get(bview) > 0]]
      [[host_height = bview && bv_height_get(bview) > 0]]
      [[bv_dimensions_set(bview]])
    string(FIND "${_ged_dm_endpoint}" "${_needle}" _endpoint_size_idx)
    if(_endpoint_size_idx EQUAL -1)
      _bobol_guard_fail(
	"GED endpoint open lost its canonical host-size bootstrap (${_needle})")
    endif()
  endforeach()
  foreach(_rel
      src/mged/mged_obol_smoke.sh
      src/mged/mged_obol_ert_smoke.sh
      src/mged/mged_obol_progressive_lod.sh)
    _bobol_guard_read_rel(_mged_smoke "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_mged_smoke}"
      [[-a[ \t]+obol]])
  endforeach()
  foreach(_rel src/mged/adc.c src/mged/grid.c)
    _bobol_guard_read_rel(_retained_faceplate_source "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_retained_faceplate_source}"
      [[(^|[^A-Za-z0-9_])dm_draw_(line|point)_2d[ \t\r\n]*\(]]
      [[(^|[^A-Za-z0-9_])dm_set_(fg|line_attr)[ \t\r\n]*\(]])
  endforeach()
  _bobol_guard_read_rel(_mged_axes "src/mged/axes.c")
  _bobol_guard_forbid_regexes("src/mged/axes.c" "${_mged_axes}"
    [[(^|[^A-Za-z0-9_])dm_draw_hud_axes[ \t\r\n]*\(]])
  foreach(_needle
      [[ged_view_feature_batch_hud_axes_replace]]
      [[_faceplate/edit_axes/initial]]
      [[_faceplate/edit_axes/current]])
    string(FIND "${_mged_axes}" "${_needle}" _retained_edit_axes_idx)
    if(_retained_edit_axes_idx EQUAL -1)
      _bobol_guard_fail(
	"MGED edit axes lost retained publication (${_needle})")
    endif()
  endforeach()
  _bobol_guard_read_rel(_mged_rect "src/mged/rect.c")
  _bobol_guard_forbid_regexes("src/mged/rect.c" "${_mged_rect}"
    [[(^|[^A-Za-z0-9_])dm_draw_line_2d[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])dm_set_(fg|line_attr)[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])dm_(Xx2Normal|Xy2Normal|Normal2Xx|Normal2Xy|get_width|get_height|get_aspect)[ \t\r\n]*\(]])
  foreach(_retired_renderer
      src/libdm/adc.c
      src/libdm/axes.c
      src/libdm/clip.c
      src/libdm/grid.c
      src/libdm/labels.c
      src/libdm/rect.c
      src/libdm/scale.c)
    if(EXISTS "${BRLCAD_SOURCE_DIR}/${_retired_renderer}")
      _bobol_guard_fail(
	"${_retired_renderer} restored a retired direct faceplate renderer")
    endif()
  endforeach()
  foreach(_retired_surface
      include/dm/view.h
      include/dm/vlist.h
      src/libdm/view.c
      src/libdm/null/dm-Null.h)
    if(EXISTS "${BRLCAD_SOURCE_DIR}/${_retired_surface}")
      _bobol_guard_fail(
	"${_retired_surface} restored a retired libdm drawing compatibility surface")
    endif()
  endforeach()
  _bobol_guard_read_rel(_mged_events "src/mged/doevent.c")
  _bobol_guard_forbid_regexes("src/mged/doevent.c" "${_mged_events}"
    [[(^|[^A-Za-z0-9_])dm_(Xx2Normal|Xy2Normal|get_width|get_height|get_aspect|configure_win|doevent)[ \t\r\n]*\(]])
  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/mged/dm-generic.c")
    _bobol_guard_fail(
      "src/mged/dm-generic.c restored a retired display-manager command implementation")
  endif()
  _bobol_guard_read_rel(_mged_display_commands "src/mged/display-command.c")
  _bobol_guard_forbid_regexes("src/mged/display-command.c" "${_mged_display_commands}"
    [[(^|[^A-Za-z0-9_])dm_(Xx2Normal|Xy2Normal|get_width|get_height|set_width|set_height)[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])dm_(get_bg|set_bg|make_current|internal_var)[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])(view_state_flag_hook|dirty_hook|zclip_hook|set_hook_data)[ \t\r\n]*\(]]
    [[mged_view_hook_state]])
  _bobol_guard_read_rel(_mged_draw_command "src/mged/cmd.cpp")
  _bobol_guard_forbid_regexes("src/mged/cmd.cpp" "${_mged_draw_command}"
    [[(^|[^A-Za-z0-9_])dm_get_(width|height)[ \t\r\n]*\(]])
  _bobol_guard_read_rel(_mged_settings "src/mged/set.c")
  _bobol_guard_forbid_regexes("src/mged/set.c" "${_mged_settings}"
    [[(^|[^A-Za-z0-9_])dm_set_perspective[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])bv_perspective_set[ \t\r\n]*\(]])
  string(FIND "${_mged_settings}" [[ged_view_context_display_property_set]]
    _mged_perspective_policy_idx)
  if(_mged_perspective_policy_idx EQUAL -1)
    _bobol_guard_fail(
      "src/mged/set.c no longer routes perspective settings through endpoint policy")
  endif()
  _bobol_guard_read_rel(_mged_fbserv "src/mged/fbserv.c")
  _bobol_guard_forbid_regexes("src/mged/fbserv.c" "${_mged_fbserv}"
    [[#[ \t]*include[ \t]*[<"]dm/fbserv_legacy\.h]]
    [[(^|[^A-Za-z0-9_])dm_interp[ \t\r\n]*\(]]
    [[(^|[^A-Za-z0-9_])dm_fbserv_set_framebuffer[ \t\r\n]*\(]])
  foreach(_needle
      [[s->gedp->ged_fbs]]
      [[mged_obol_framebuffer_ensure]]
      [[fbs_framebuffer_backend_installed]])
    string(FIND "${_mged_fbserv}" "${_needle}" _mged_fbserv_endpoint_idx)
    if(_mged_fbserv_endpoint_idx EQUAL -1)
      _bobol_guard_fail(
	"MGED fbserv no longer uses the session-owned Obol/imgstream backend (${_needle})")
    endif()
  endforeach()
  foreach(_rel
      src/mged/mged_display.h
      src/mged/attach.c
      src/mged/mged.c)
    _bobol_guard_read_rel(_mged_null_dm_shell "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_mged_null_dm_shell}"
      [[(^|[^A-Za-z0-9_])dm_dmp([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])dm_fbp([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])dm_fbserv([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])DMP([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])dm_(open|close|get_fb)[ \t\r\n]*\(]])
  endforeach()
  _bobol_guard_read_rel(_mged_attach "src/mged/attach.c")
  foreach(_needle
      [[mged_obol_framebuffer_ensure]]
      [[ged_view_framebuffer_backend_ensure]])
    string(FIND "${_mged_attach}" "${_needle}" _mged_fb_backend_idx)
    if(_mged_fb_backend_idx EQUAL -1)
      _bobol_guard_fail(
	"MGED attachment lost endpoint-owned framebuffer setup (${_needle})")
    endif()
  endforeach()
  foreach(_needle
      [[mged_obol_framebuffer_ensure]]
      [[fbs_framebuffer_backend_installed]])
    string(FIND "${_mged_rect}" "${_needle}" _mged_rect_fb_idx)
    if(_mged_rect_fb_idx EQUAL -1)
      _bobol_guard_fail(
	"MGED rectangle rendering lost the endpoint framebuffer contract (${_needle})")
    endif()
  endforeach()
  foreach(_rel
      src/mged/mged_display.h
      src/mged/attach.c
      src/mged/mged.c
      src/mged/set.c
      src/mged/share.c
      src/mged/cmd.cpp)
    _bobol_guard_read_rel(_mged_backend_cache "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_mged_backend_cache}"
      [[(^|[^A-Za-z0-9_])mv_backend_cache([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])_backend_cache_state([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])dm_backend_cache_state([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])mged_default_backend_cache([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])share_backend_cache([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])dm_get_backend_cache[ \t\r\n]*\(]]
      [[--cache]]
      [[no-cache]])
  endforeach()
  foreach(_rel src/mged/chgview.c src/mged/f_cmd.h src/mged/setup.c)
    _bobol_guard_read_rel(_mged_register_debug "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_mged_register_debug}"
      [[(^|[^A-Za-z0-9_])f_regdebug([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])dm_set_debug[ \t\r\n]*\(]])
  endforeach()
  _bobol_guard_read_rel(_mged_help "src/tclscripts/mged/help.tcl")
  _bobol_guard_forbid_regexes("src/tclscripts/mged/help.tcl" "${_mged_help}"
    [[mged_help_data\(regdebug\)]])
  foreach(_rel src/tclscripts/mged/openw.tcl src/tclscripts/mged/mgedrc.tcl)
    _bobol_guard_read_rel(_mged_backend_cache_tcl "${_rel}")
    foreach(_needle [[mged_gui($id,cache)]] [[mged_default(cache)]] [[Backend Cache]])
      string(FIND "${_mged_backend_cache_tcl}" "${_needle}" _mged_cache_idx)
      if(NOT _mged_cache_idx EQUAL -1)
        _bobol_guard_fail("${_rel} restored retired backend-cache policy (${_needle})")
      endif()
    endforeach()
  endforeach()
  foreach(_rel
      src/mged/mged_display.h
      src/mged/mged.c
      src/mged/cmd.cpp
      src/mged/edsol.c
      src/mged/share.c)
    _bobol_guard_read_rel(_mged_identity "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_mged_identity}"
      [[(^|[^A-Za-z0-9_])dm_get_(pathname|id|tkname)[ \t\r\n]*\(]])
  endforeach()
  _bobol_guard_read_rel(_mged_bindings "src/tclscripts/mged/bindings.tcl")
  _bobol_guard_forbid_regexes("src/tclscripts/mged/bindings.tcl"
    "${_mged_bindings}"
    [[(^|[^A-Za-z0-9_])tkswrast([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])mged_bind_dm([^A-Za-z0-9_]|$)]])
  foreach(_rel
      src/mged/attach.c
      src/mged/f_cmd.h
      src/mged/cmd.h
      src/mged/cmd.cpp
      src/mged/setup.c
      src/mged/mged.c
      src/qged/QgEdMainWindow.cpp
      src/qged/QgEdMainWindow.h
      src/libged/dm/dm.cpp
      src/tclscripts/mged/apply.tcl
      src/tclscripts/mged/helpdevel.tcl)
    _bobol_guard_read_rel(_mged_stale_endpoint_names "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_mged_stale_endpoint_names}"
      [[(^|[^A-Za-z0-9_])get_dm_list([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])f_get_dm_list([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])cmd_ged_dm_wrapper([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])attach_display_manager([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])qged_dm_during_clbk([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])do_dm_init([^A-Za-z0-9_]|$)]]
      [[\.dm_tkobol]]
      [[GED_DM_DURING_DEBUG]])
  endforeach()
  foreach(_rel
      src/tclscripts/mged/adc.tcl
      src/tclscripts/mged/apply.tcl
      src/tclscripts/mged/bindings.tcl
      src/tclscripts/mged/collaborate.tcl
      src/tclscripts/mged/color_scheme.tcl
      src/tclscripts/mged/comb.tcl
      src/tclscripts/mged/grid.tcl
      src/tclscripts/mged/grouper.tcl
      src/tclscripts/mged/mged.tcl
      src/tclscripts/mged/mgedrc.tcl
      src/tclscripts/mged/mouse.tcl
      src/tclscripts/mged/mview.tcl
      src/tclscripts/mged/openw.tcl
      src/tclscripts/mged/qray.tcl
      src/tclscripts/mged/rt.tcl)
    _bobol_guard_read_rel(_mged_stale_pane_state "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_mged_stale_pane_state}"
      [[(^|[^A-Za-z0-9_])active_dm([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])show_dm([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])set_dm_win([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])dm_win_hide([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])dm_loc([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])dm_key_bindings([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])dm_id([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])new_dm([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])set_active_dm([^A-Za-z0-9_]|$)]])
  endforeach()
  foreach(_rel
      src/mged/mged.c
      src/tclscripts/mged/openw.tcl
      src/tclscripts/mged/mview.tcl
      src/tclscripts/mged/mgedrc.tcl
      src/tclscripts/mged/help.tcl
      doc/asciidoc/system/man1/mged.adoc
      doc/asciidoc/system/mann/gui.adoc
      doc/asciidoc/articles/mged.adoc)
    _bobol_guard_read_rel(_mged_host_term "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_mged_host_term}"
      [[--dm-type]]
      [[mged_default\(dm_type\)]]
      [[gui[ \t]+-dt]]
      [[(^|[^A-Za-z0-9_])-dt([^A-Za-z0-9_]|$)]])
  endforeach()
  foreach(_rel include/tclcad/misc.h include/tclcad/setup.h)
    _bobol_guard_read_rel(_tclcad_public_header "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_tclcad_public_header}"
      [[#[ \t]*include[ \t]*[<"]dm\.h]])
  endforeach()
  _bobol_guard_read_rel(_tclcad_wrapper "src/libtclcad/wrapper.c")
  _bobol_guard_forbid_regexes("src/libtclcad/wrapper.c" "${_tclcad_wrapper}"
    [[(^|[^A-Za-z0-9_])to_dm_func([^A-Za-z0-9_]|$)]])
  _bobol_guard_read_rel(_rt_edit_header "include/rt/edit.h")
  _bobol_guard_forbid_regexes("include/rt/edit.h" "${_rt_edit_header}"
    [[(^|[^A-Za-z0-9_])dm_loadmatrix([^A-Za-z0-9_]|$)]])
  foreach(_rel
      src/libtclcad/commands.c
      src/tclscripts/lib/Ged.tcl
      src/tclscripts/archer/Archer.tcl
      src/tclscripts/archer/ArcherCore.tcl)
    _bobol_guard_read_rel(_retired_backend_cache "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_retired_backend_cache}"
      [[(^|[^A-Za-z0-9_])(dm_view_data|cache_on|mBackendCacheMode)([^A-Za-z0-9_]|$)]])
  endforeach()
  _bobol_guard_read_rel(_mged_edsol "src/mged/edsol.c")
  _bobol_guard_forbid_regexes("src/mged/edsol.c" "${_mged_edsol}"
    [[(^|[^A-Za-z0-9_])dm_get_dname[ \t\r\n]*\(]])
  foreach(_needle
      [[mged_rubber_band_state_sync]]
      [[ged_view_feature_batch_hud_line_set_replace]]
      [[_faceplate/rubber_band]])
    string(FIND "${_mged_rect}" "${_needle}" _retained_rubber_band_idx)
    if(_retained_rubber_band_idx EQUAL -1)
      _bobol_guard_fail(
	"MGED rubber band lost retained publication (${_needle})")
    endif()
  endforeach()

  foreach(_rel src/mged/titles.c src/mged/menu.c src/mged/scroll.c)
    _bobol_guard_read_rel(_retained_mged_hud "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_retained_mged_hud}"
      [[(^|[^A-Za-z0-9_])dm_draw_(line_2d|string_2d)[ \t\r\n]*\(]]
      [[(^|[^A-Za-z0-9_])dm_set_(fg|line_attr)[ \t\r\n]*\(]])
  endforeach()
  _bobol_guard_read_rel(_mged_titles "src/mged/titles.c")
  _bobol_guard_forbid_regexes("src/mged/titles.c" "${_mged_titles}"
    [[(^|[^A-Za-z0-9_])dm_get_(fontsize|aspect)[ \t\r\n]*\(]])
  foreach(_needle
      [[mged_hud_builder_init]]
      [[mged_hud_builder_publish]]
      [[scroll_display(s, &hud]]
      [[mmenu_display(s, &hud]])
    string(FIND "${_mged_titles}" "${_needle}" _retained_mged_hud_idx)
    if(_retained_mged_hud_idx EQUAL -1)
      _bobol_guard_fail(
	"MGED titles lost retained menu/scroll HUD batching (${_needle})")
    endif()
  endforeach()

  foreach(_rel
      include/tclcad/setup.h
      src/libtclcad/tclcad_private.h
      src/libtclcad/tkobol_init.c
      src/libtclcad/fb.c)
    _bobol_guard_read_rel(_contents "${_rel}")
    _bobol_guard_forbid_regexes("${_rel}" "${_contents}"
      [[(^|[^A-Za-z0-9_])Fbo_Init([^A-Za-z0-9_]|$)]]
      [[Tcl_CreateCommand\([^\n]*"fb_open"]])
  endforeach()

endfunction()

function(_bobol_guard_collect_active_scan_files _outvar)
  _bobol_guard_collect(_files
    include/BObol
    include/BObol.h
    include/qtcad
    src/libBObol
    src/libqtcad
    src/qged
    src/gtools/gsh)
  set(${_outvar} "${_files}" PARENT_SCOPE)
endfunction()

_bobol_guard_check_dependency_inventory()
_bobol_guard_check_libbsg_retired()
_bobol_guard_check_legacy_graphical_dm_retired()
_bobol_guard_check_retired_tests()
_bobol_guard_check_public_headers()
_bobol_guard_check_retired_ged_draw_symbols()
_bobol_guard_check_removed_rt_view_aliases()
_bobol_guard_check_native_drawing_tests()
_bobol_guard_check_obol_controller_quarantine()
_bobol_guard_check_display_endpoint_boundary()
_bobol_guard_check_imgstream_display_host_ownership()
_bobol_guard_check_imgstream_remote_client_ownership()
_bobol_guard_check_standalone_fb_ownership()
_bobol_guard_check_rt_imgstream_ownership()
_bobol_guard_check_standalone_display_provider_boundary()
_bobol_guard_check_rtwizard_host_boundary()
_bobol_guard_check_ged_framebuffer_image_ownership()
_bobol_guard_check_qged_edit_preview_policy()
_bobol_guard_check_tclcad_obol_readback_bridge()
_bobol_guard_check_fbserv_legacy_framebuffer_retired()
_bobol_guard_check_fbserv_transport_quarantine()
_bobol_guard_check_fbserv_transport_setup_helpers()
_bobol_guard_check_qged_fbserv_backend_access()
_bobol_guard_check_fbserv_backend_contract()
_bobol_guard_check_tkobol_host_ownership()
_bobol_guard_check_retained_export_ownership()
_bobol_guard_check_retired_text_dm()
_bobol_guard_check_unreachable_drawing_shims()

_bobol_guard_collect_active_scan_files(_bobol_guard_scan_files)
foreach(_file IN LISTS _bobol_guard_scan_files)
  _bobol_guard_check_file("${_file}")
endforeach()

get_property(_failures GLOBAL PROPERTY BOBOL_PIVOT_GUARD_FAILURES)
if(_failures)
  list(SORT _failures)
  string(REPLACE ";" "\n  " _failure_text "${_failures}")
  message(FATAL_ERROR
    "Obol pivot guard failed:\n"
    "  ${_failure_text}\n"
    "New Obol-canonical code must not depend on retired BSG/DM APIs.")
endif()

message(STATUS "Obol pivot guard passed")
