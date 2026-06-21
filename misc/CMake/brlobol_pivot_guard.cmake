#                 B R L O B O L _ P I V O T _ G U A R D . C M A K E
# BRL-CAD
#
# Copyright (c) 2026 United States Government as represented by
# the U.S. Army Research Laboratory.
#
# Guardrails for the Obol-canonical drawing migration.

if(NOT DEFINED BRLCAD_SOURCE_DIR)
  message(FATAL_ERROR "BRLCAD_SOURCE_DIR is required")
endif()

if(NOT DEFINED BRLOBOL_PIVOT_GUARD_MODE)
  set(BRLOBOL_PIVOT_GUARD_MODE "transitional")
endif()

set(_patterns
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

set(_extensions [[(CMakeLists\.txt|\.(c|cc|cpp|cxx|h|hh|hpp|cmake)$)]])

function(_brlobol_pivot_guard_fail _msg)
  set_property(GLOBAL APPEND PROPERTY BRLOBOL_PIVOT_GUARD_FAILURES "${_msg}")
endfunction()

function(_brlobol_pivot_guard_collect _outvar)
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

function(_brlobol_pivot_guard_check_file _file)
  file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
  if("${_rel}" STREQUAL "misc/CMake/brlobol_pivot_guard.cmake")
    return()
  endif()
  if(NOT "${_rel}" MATCHES "${_extensions}")
    return()
  endif()

  file(READ "${_file}" _contents)
  foreach(_pat IN LISTS _patterns)
    string(REGEX MATCH "${_pat}" _hit "${_contents}")
    if(_hit)
      _brlobol_pivot_guard_fail("${_rel}: ${_hit}")
    endif()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_require_inventory_token _contents _token)
  string(FIND "${_contents}" "${_token}" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail(
      "doc/notes/obol_legacy_dependency_inventory.txt missing: ${_token}")
  endif()
endfunction()

function(_brlobol_pivot_guard_require_dependency_row _source_dirs _inventory _target _deps)
  string(REGEX MATCH "set_deps\\(${_target}[ \t]+\"${_deps}\"\\)" _source_hit
    "${_source_dirs}")
  if(NOT _source_hit)
    _brlobol_pivot_guard_fail(
      "src/source_dirs.cmake direct dependency row changed for ${_target}; update the Obol legacy dependency inventory and guard allowlist")
  endif()
  _brlobol_pivot_guard_require_inventory_token("${_inventory}"
    "${_target}: ${_deps}")
endfunction()

function(_brlobol_pivot_guard_check_dependency_inventory)
  set(_inventory_file
    "${BRLCAD_SOURCE_DIR}/doc/notes/obol_legacy_dependency_inventory.txt")
  if(NOT EXISTS "${_inventory_file}")
    _brlobol_pivot_guard_fail(
      "doc/notes/obol_legacy_dependency_inventory.txt is required")
    return()
  endif()

  file(READ "${_inventory_file}" _inventory)
  file(READ "${BRLCAD_SOURCE_DIR}/src/source_dirs.cmake" _source_dirs)

  _brlobol_pivot_guard_require_dependency_row("${_source_dirs}" "${_inventory}"
    libbrlobol [[libwdb;librt;libimgstream;libbg;libbu]])
  _brlobol_pivot_guard_require_dependency_row("${_source_dirs}" "${_inventory}"
    libnmg [[libbg;libbn;libbu]])
  _brlobol_pivot_guard_require_dependency_row("${_source_dirs}" "${_inventory}"
    libbrep [[libbg;libbu]])
  _brlobol_pivot_guard_require_dependency_row("${_source_dirs}" "${_inventory}"
    libanalyze [[librt;libbg;libbn;libbu]])
  _brlobol_pivot_guard_require_dependency_row("${_source_dirs}" "${_inventory}"
    libimgstream [[libicv;libbn;libbu]])
  _brlobol_pivot_guard_require_dependency_row("${_source_dirs}" "${_inventory}"
    libdm [[libbsg;librt;libicv;libbn;libpkg;libbu]])
  _brlobol_pivot_guard_require_dependency_row("${_source_dirs}" "${_inventory}"
    libged [[libicv;libanalyze;libwdb;liboptical;libdm;libbsg;libbu]])
  _brlobol_pivot_guard_require_dependency_row("${_source_dirs}" "${_inventory}"
    libqtcad [[libbrlobol;libged;libdm;libbsg;libbg;libbn;libbu]])
  _brlobol_pivot_guard_require_dependency_row("${_source_dirs}" "${_inventory}"
    libtclcad [[libged;libdm;libbsg;libbn;libbu]])
  _brlobol_pivot_guard_require_dependency_row("${_source_dirs}" "${_inventory}"
    qged [[libqtcad]])
  _brlobol_pivot_guard_require_dependency_row("${_source_dirs}" "${_inventory}"
    mged [[libtclcad]])

  file(READ "${BRLCAD_SOURCE_DIR}/src/qged/CMakeLists.txt" _qged_cmake)
  string(FIND "${_qged_cmake}" "add_dependencies(qged dm_plugins ged_plugins)"
    _qged_plugins_idx)
  if(_qged_plugins_idx EQUAL -1)
    _brlobol_pivot_guard_fail(
      "qged dm_plugins dependency changed; update the Obol legacy dependency inventory and guard allowlist")
  endif()
  _brlobol_pivot_guard_require_inventory_token("${_inventory}"
    "qged plugin dependency: dm_plugins")
endfunction()

function(_brlobol_pivot_guard_check_obol_realization_coverage)
  set(_coverage_file
    "${BRLCAD_SOURCE_DIR}/doc/notes/obol_realization_coverage.txt")
  if(NOT EXISTS "${_coverage_file}")
    _brlobol_pivot_guard_fail(
      "doc/notes/obol_realization_coverage.txt is required")
    return()
  endif()

  file(READ "${_coverage_file}" _coverage)
  foreach(_token
      "share/db/pinewood.g"
      "share/db/havoc.g"
      "m35.g"
      "share/db/faa/Generic_Twin.g"
      "largest in-tree example"
      "valid and invalid BoTs"
      "implicit solids"
      "FAA maturity"
      "BRLOBOL_QTCAD_GENERIC_TWIN=1"
      "controlled Generic_Twin.g opt-in workflow gate")
    string(FIND "${_coverage}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail(
	"doc/notes/obol_realization_coverage.txt missing coverage token ${_token}")
    endif()
  endforeach()

  set(_qtcad_real_models_test
    "${BRLCAD_SOURCE_DIR}/src/libqtcad/tests/test_qtcad_obol_real_models.cpp")
  if(EXISTS "${_qtcad_real_models_test}")
    file(READ "${_qtcad_real_models_test}" _qtcad_real_models_contents)
    foreach(_token
	"BRLOBOL_QTCAD_GENERIC_TWIN"
	"faa/Generic_Twin.g"
	"qtcad Obol Generic_Twin maturity workflow should pass")
      string(FIND "${_qtcad_real_models_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/tests/test_qtcad_obol_real_models.cpp missing Generic_Twin opt-in token ${_token}")
      endif()
    endforeach()
  endif()

  set(_qtcad_tests_cmake
    "${BRLCAD_SOURCE_DIR}/src/libqtcad/tests/CMakeLists.txt")
  if(EXISTS "${_qtcad_tests_cmake}")
    file(READ "${_qtcad_tests_cmake}" _qtcad_tests_cmake_contents)
    string(FIND "${_qtcad_tests_cmake_contents}" "Generic_Twin.g" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail(
	"src/libqtcad/tests/CMakeLists.txt must keep the Generic_Twin.g opt-in model dependency")
    endif()
  endif()

  set(_prototype_test
    "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_prototype.cpp")
  if(EXISTS "${_prototype_test}")
    file(READ "${_prototype_test}" _prototype_test_contents)
    foreach(_token
	"exercise_required_hierarchy_model"
	"required pinewood hierarchy coverage should pass"
	"required havoc hierarchy coverage should pass"
	"havoc.g stays wire-only"
	"hierarchyExport.getLineCount() < min_wire_segments"
	"hierarchyWireMeasure.setAngleComputationEnabled(FALSE)"
	"hierarchyWireMeasure.getShapeCount() < min_wire_shapes"
	"hierarchyMeshExport.getTriangleCount() < min_mesh_triangles"
	"hierarchyMeshMeasure.getShapeCount() < min_mesh_shapes"
	"view controller should expose cache-backed RT exact ray picks outside qtcad"
	"sparse bucketed angle measure should find the connected corner")
      string(FIND "${_prototype_test_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/libbrlobol/tests/test_prototype.cpp missing required hierarchy coverage token ${_token}")
      endif()
    endforeach()
    string(FIND "${_prototype_test_contents}" "total_triangle_count" _idx)
    if(_idx GREATER -1)
      _brlobol_pivot_guard_fail(
	"src/libbrlobol/tests/test_prototype.cpp must not reintroduce slow required-hierarchy total_triangle_count traversal")
    endif()
  endif()

  set(_measure_action
    "${BRLCAD_SOURCE_DIR}/src/libbrlobol/measure_action.cpp")
  if(EXISTS "${_measure_action}")
    file(READ "${_measure_action}" _measure_action_contents)
    foreach(_token
	"measure_collect_angle_endpoint_candidates"
	"measure_add_angle_endpoint"
	"measure_source_local_query_distance_limit"
	"setQueryDistanceLimit"
	"queryDistanceAllows"
	"shape->point.getNum()"
	"shape->command[i] != SoBRLVListShape::DRAW"
	"const int currentSegment = segmentIndex++"
	"connectedPairs")
      string(FIND "${_measure_action_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/libbrlobol/measure_action.cpp missing bucketed connected-angle token ${_token}")
      endif()
    endforeach()
    string(FIND "${_measure_action_contents}" "shape->getSegment(" _idx)
    if(_idx GREATER -1)
      _brlobol_pivot_guard_fail(
	"src/libbrlobol/measure_action.cpp must keep vlist measurement on a linear command-stream traversal")
    endif()
  endif()

  set(_snap_action
    "${BRLCAD_SOURCE_DIR}/src/libbrlobol/snap_action.cpp")
  if(EXISTS "${_snap_action}")
    file(READ "${_snap_action}" _snap_action_contents)
    foreach(_token
	"shape->point.getNum()"
	"shape->command[i] == SoBRLVListShape::DRAW"
	"const int currentSegment = segmentIndex++")
      string(FIND "${_snap_action_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/libbrlobol/snap_action.cpp missing linear vlist traversal token ${_token}")
      endif()
    endforeach()
    string(FIND "${_snap_action_contents}" "shape->getSegment(" _idx)
    if(_idx GREATER -1)
      _brlobol_pivot_guard_fail(
	"src/libbrlobol/snap_action.cpp must keep vlist snapping on a linear command-stream traversal")
    endif()
  endif()
endfunction()

function(_brlobol_pivot_guard_check_libimgstream_boundary)
  set(_required_files
    include/imgstream.h
    include/imgstream/fb_compat.h
    include/imgstream/stream.h
    include/imgstream/CMakeLists.txt
    src/libimgstream/CMakeLists.txt
    src/libimgstream/fb_compat.c
    src/libimgstream/stream.c
    src/libimgstream/tests/CMakeLists.txt
    src/libimgstream/tests/fb_compat.c
    src/libimgstream/tests/stream.c)

  foreach(_rel IN LISTS _required_files)
    if(NOT EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _brlobol_pivot_guard_fail("${_rel} is required for the initial libimgstream migration batch")
    endif()
  endforeach()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libimgstream/CMakeLists.txt")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libimgstream/CMakeLists.txt" _imgstream_cmake)
    foreach(_token
	"brlcad_addlib(libimgstream"
	"fb_compat.c"
	"PRIVATE_LIBS \${libimgstream_deps}"
	"add_subdirectory(tests)"
      )
      string(FIND "${_imgstream_cmake}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libimgstream/CMakeLists.txt missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libimgstream/tests/CMakeLists.txt")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libimgstream/tests/CMakeLists.txt" _imgstream_tests_cmake)
    foreach(_token
	"brlcad_addexec(imgstream_test"
	"brlcad_add_test(NAME imgstream_test"
	"brlcad_addexec(imgstream_fb_compat_test"
	"brlcad_add_test(NAME imgstream_fb_compat_test"
      )
      string(FIND "${_imgstream_tests_cmake}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libimgstream/tests/CMakeLists.txt missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/include/imgstream/fb_compat.h")
    file(READ "${BRLCAD_SOURCE_DIR}/include/imgstream/fb_compat.h" _imgstream_fb_header)
    foreach(_token
	"struct imgstream_fb"
	"IMGSTREAM_FB_SPEC_MEMORY"
	"IMGSTREAM_FB_SPEC_FILE"
	"IMGSTREAM_FB_SPEC_FANOUT"
	"IMGSTREAM_FB_SPEC_DIAGNOSTIC"
	"IMGSTREAM_FB_SPEC_REMOTE"
	"IMGSTREAM_FB_SPEC_DISPLAY"
	"enum imgstream_fb_remote_kind"
	"IMGSTREAM_FB_REMOTE_TCP"
	"IMGSTREAM_FB_REMOTE_IPC"
	"enum imgstream_fb_display_kind"
	"IMGSTREAM_FB_DISPLAY_QTGL"
	"IMGSTREAM_FB_DISPLAY_SWRAST"
	"enum imgstream_fb_diagnostic_kind"
	"struct imgstream_fb_spec_info"
	"struct imgstream_fb_display_host"
	"struct imgstream_fb_colormap"
	"imgstream_fb_spec_kind"
	"imgstream_fb_spec_info"
	"imgstream_fb_display_host_set"
	"imgstream_fb_display_host_clear"
	"imgstream_fb_open"
	"imgstream_fb_colormap_linear"
	"imgstream_fb_rmap"
	"imgstream_fb_wmap"
	"imgstream_fb_readrect"
	"imgstream_fb_writerect"
	"imgstream_fb_bwreadrect"
	"imgstream_fb_bwwriterect"
	"imgstream_fb_flush"
	"imgstream_fb_reset"
	"imgstream_fb_viewport"
	"imgstream_fb_view"
	"imgstream_fb_window"
	"imgstream_fb_zoom"
	"imgstream_fb_cursor"
	"imgstream_fb_scursor"
	"imgstream_fb_setcursor"
	"imgstream_fb_poll_rate")
      string(FIND "${_imgstream_fb_header}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("include/imgstream/fb_compat.h missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libimgstream/fb_compat.c")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libimgstream/fb_compat.c" _imgstream_fb_impl)
    foreach(_token
	"IMGSTREAM_FB_DEFAULT_SIZE"
	"imgstream_fb_spec_supported"
	"fb_remote_info_from_spec"
	"fb_remote_info_from_legacy_spec"
	"fb_display_host_registered"
	"imgstream_fb_display_host_set"
	"imgstream_fb_display_host_clear"
	"imgstream_fb_open"
	"fb_file_path_from_spec"
	"fb_stream_flush_file"
	"IMGSTREAM_FB_FANOUT_MAX"
	"fb_fanout_parse"
	"fb_open_fanout"
	"IMGSTREAM_FB_DIAGNOSTIC_DEBUG"
	"fb_diagnostic_kind_from_spec"
	"fb_diagnostic_log_pixels"
	"imgstream_fb_colormap_linear"
	"imgstream_fb_rmap"
	"imgstream_fb_wmap"
	"imgstream_fb_read"
	"imgstream_fb_write"
	"imgstream_fb_bwreadrect"
	"imgstream_fb_bwwriterect"
	"imgstream_fb_flush"
	"imgstream_fb_reset"
	"imgstream_fb_viewport"
	"imgstream_fb_window"
	"imgstream_fb_zoom"
	"imgstream_fb_scursor"
	"imgstream_fb_setcursor"
	"imgstream_fb_poll"
	"imgstream_fb_poll_rate")
      string(FIND "${_imgstream_fb_impl}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libimgstream/fb_compat.c missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libimgstream/tests/fb_compat.c")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libimgstream/tests/fb_compat.c" _imgstream_fb_test)
    foreach(_token
	"test_file_compat_stream"
	"test_colormap_compat_state"
	"test_fanout_compat_stream"
	"test_diagnostic_compat_stream"
	"imgstream_fb_spec_info"
	"IMGSTREAM_FB_REMOTE_TCP"
	"IMGSTREAM_FB_REMOTE_IPC"
	"IMGSTREAM_FB_DISPLAY_QTGL"
	"test_display_host_compat_stream"
	"imgstream_fb_display_host_set"
	"imgstream_fb_display_host_clear"
	"imgstream_fb_flush"
	"imgstream_fb_bwreadrect"
	"imgstream_fb_bwwriterect"
	"imgstream_fb_reset"
	"imgstream_fb_viewport"
	"imgstream_fb_window"
	"imgstream_fb_zoom"
	"imgstream_fb_scursor"
	"imgstream_fb_setcursor"
	"imgstream_fb_poll_rate"
	"/dev/disk:"
	"/dev/stack"
	"/dev/debug"
	"/dev/null")
      string(FIND "${_imgstream_fb_test}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libimgstream/tests/fb_compat.c missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/include/imgstream/stream.h")
    file(READ "${BRLCAD_SOURCE_DIR}/include/imgstream/stream.h" _imgstream_header)
    foreach(_token
	"struct imgstream"
	"IMGSTREAM_PIXEL_RGB8"
	"IMGSTREAM_PIXEL_RGBA8"
	"imgstream_create"
	"imgstream_get_info"
	"imgstream_resize"
	"imgstream_write_rect"
	"imgstream_dirty_rect"
	"imgstream_subscribe"
	"imgstream_producer_begin"
	"imgstream_create_from_icv"
	"imgstream_snapshot_icv")
      string(FIND "${_imgstream_header}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("include/imgstream/stream.h missing ${_token}")
      endif()
    endforeach()
  endif()

  _brlobol_pivot_guard_collect(_imgstream_files
    include/imgstream
    include/imgstream.h
    src/libimgstream)

  foreach(_file IN LISTS _imgstream_files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if(NOT "${_rel}" MATCHES "${_extensions}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    set(_forbidden
      [[#[ \t]*include[ \t]*[<"]bsg]]
      [[#[ \t]*include[ \t]*[<"]dm]]
      [[(^|[^A-Za-z0-9_])libbsg([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])libbrlobol([^A-Za-z0-9_]|$)]]
      [[#[ \t]*include[ \t]*[<"]brlobol]]
      [[(^|[^A-Za-z0-9_])SoBRL[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]])
    foreach(_pat IN LISTS _forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail("${_rel} breaks libimgstream non-display boundary: ${_hit}")
      endif()
    endforeach()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_brlobol_image_source)
  set(_required_files
    include/brlobol/image_source.h
    src/libbrlobol/image_source.cpp
    src/libbrlobol/tests/test_image_source.cpp)

  foreach(_rel IN LISTS _required_files)
    if(NOT EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _brlobol_pivot_guard_fail("${_rel} is required for Obol image source ownership")
    endif()
  endforeach()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/include/brlobol/image_source.h")
    file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/image_source.h" _image_source_header)
    foreach(_token
	"SoBRLImageSource"
	"SourceKind"
	"SOURCE_STATIC_IMAGE"
	"SOURCE_IMAGE_STREAM"
	"Status"
	"PixelFormat"
	"SoSFString imageId"
	"SoSFEnum sourceKind"
	"SoSFUInt32 dirtyRevision"
	"SoSFBool streamConnected"
	"setImage"
	"setStream"
	"imgstream_t"
	"icv_image_t")
      string(FIND "${_image_source_header}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("include/brlobol/image_source.h missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libbrlobol/image_source.cpp")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/image_source.cpp" _image_source_impl)
    foreach(_token
	"SO_NODE_SOURCE(SoBRLImageSource)"
	"imgstream_subscribe"
	"imgstream_unsubscribe"
	"imgstream_create_from_icv"
	"imgstream_get_info"
	"dirtyCB")
      string(FIND "${_image_source_impl}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libbrlobol/image_source.cpp missing ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg]]
	[[#[ \t]*include[ \t]*[<"]dm]]
	[[(^|[^A-Za-z0-9_])libbsg([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _hit "${_image_source_impl}")
      if(_hit)
	_brlobol_pivot_guard_fail("src/libbrlobol/image_source.cpp reintroduced a legacy display dependency: ${_hit}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libbrlobol/init.cpp")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/init.cpp" _brlobol_init)
    string(FIND "${_brlobol_init}" "SoBRLImageSource::initClass()" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/init.cpp must initialize SoBRLImageSource")
    endif()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libbrlobol/CMakeLists.txt")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/CMakeLists.txt" _brlobol_cmake)
    foreach(_token
	"image_source.cpp"
	"libimgstream")
      string(FIND "${_brlobol_cmake}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libbrlobol/CMakeLists.txt missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_image_source.cpp")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_image_source.cpp" _image_source_test)
    foreach(_token
	"imgstream_create"
	"imgstream_write_rect"
	"imgstream_producer_begin"
	"icv_create"
	"setImage"
	"setStream")
      string(FIND "${_image_source_test}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libbrlobol/tests/test_image_source.cpp missing ${_token}")
      endif()
    endforeach()
  endif()
endfunction()

function(_brlobol_pivot_guard_check_brlobol_image_display)
  set(_required_files
    include/brlobol/viewport_image.h
    include/brlobol/image_plane.h
    src/libbrlobol/viewport_image.cpp
    src/libbrlobol/image_plane.cpp
    src/libbrlobol/image_display_util.h
    src/libbrlobol/image_display_util.cpp
    src/libbrlobol/tests/test_image_display.cpp)

  foreach(_rel IN LISTS _required_files)
    if(NOT EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _brlobol_pivot_guard_fail("${_rel} is required for Obol image display ownership")
    endif()
  endforeach()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/include/brlobol/viewport_image.h")
    file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/viewport_image.h" _viewport_header)
    foreach(_token
	"SoBRLViewportImage"
	"SoSFNode imageSource"
	"SoSFEnum layer"
	"SoSFVec2f sourceCenter"
	"SoSFFloat sourceZoom"
	"cursorVisible"
	"realizedDirtyRevision"
	"syncFromSource")
      string(FIND "${_viewport_header}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("include/brlobol/viewport_image.h missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/include/brlobol/image_plane.h")
    file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/image_plane.h" _plane_header)
    foreach(_token
	"SoBRLImagePlane"
	"SoSFNode imageSource"
	"SoSFString sourcePath"
	"SoSFBool selectable"
	"SoSFBool depthTest"
	"SoSFBool depthWrite"
	"realizedDirtyRevision"
	"syncFromSource")
      string(FIND "${_plane_header}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("include/brlobol/image_plane.h missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libbrlobol/viewport_image.cpp")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/viewport_image.cpp" _viewport_impl)
    foreach(_token
	"SO_NODE_SOURCE(SoBRLViewportImage)"
	"SoHUDKit"
	"brlobol_image_make_textured_quad"
	"realizedDirtyRevision"
	"syncFromSource"
	"SoTexture2")
      string(FIND "${_viewport_impl}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libbrlobol/viewport_image.cpp missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libbrlobol/image_plane.cpp")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/image_plane.cpp" _plane_impl)
    foreach(_token
	"SO_NODE_SOURCE(SoBRLImagePlane)"
	"brlobol_image_make_textured_quad"
	"depthTest"
	"depthWrite"
	"SoTexture2"
	"SoFaceSet")
      string(FIND "${_plane_impl}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libbrlobol/image_plane.cpp missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libbrlobol/image_display_util.cpp")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/image_display_util.cpp" _display_util_impl)
    foreach(_token
	"imgstream_read"
	"SoTexture2"
	"SoTextureCoordinate2"
	"SoCoordinate3"
	"SoFaceSet"
	"SoDepthBuffer"
	"SoPickStyle"
	"SoShapeHints")
      string(FIND "${_display_util_impl}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libbrlobol/image_display_util.cpp missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libbrlobol/init.cpp")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/init.cpp" _brlobol_init)
    foreach(_token
	"SoBRLViewportImage::initClass()"
	"SoBRLImagePlane::initClass()")
      string(FIND "${_brlobol_init}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libbrlobol/init.cpp must initialize ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libbrlobol/CMakeLists.txt")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/CMakeLists.txt" _brlobol_cmake)
    foreach(_token
	"viewport_image.cpp"
	"image_plane.cpp"
	"image_display_util.cpp")
      string(FIND "${_brlobol_cmake}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libbrlobol/CMakeLists.txt missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/CMakeLists.txt")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/CMakeLists.txt" _brlobol_tests_cmake)
    string(FIND "${_brlobol_tests_cmake}" "test_brlobol_image_display" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/tests/CMakeLists.txt missing test_brlobol_image_display")
    endif()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_image_display.cpp")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_image_display.cpp" _display_test)
    foreach(_token
	"SoBRLViewportImage"
	"SoBRLImagePlane"
	"imgstream_write_rect"
	"syncFromSource"
	"getTextureNode"
	"getImageFaceSet")
      string(FIND "${_display_test}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libbrlobol/tests/test_image_display.cpp missing ${_token}")
      endif()
    endforeach()
  endif()

  foreach(_rel IN LISTS _required_files)
    if(NOT EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      continue()
    endif()
    file(READ "${BRLCAD_SOURCE_DIR}/${_rel}" _contents)
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg]]
	[[#[ \t]*include[ \t]*[<"]dm]]
	[[(^|[^A-Za-z0-9_])libbsg([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail("${_rel} reintroduced a legacy display dependency: ${_hit}")
      endif()
    endforeach()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_brlobol_window_host)
  set(_required_files
    include/brlobol/headless_window_host.h
    include/brlobol/window_host.h
    src/libbrlobol/headless_window_host.cpp
    src/libbrlobol/window_host.cpp
    src/libbrlobol/tests/test_headless_window_host.cpp
    src/libbrlobol/tests/test_window_host.cpp)

  foreach(_rel IN LISTS _required_files)
    if(NOT EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _brlobol_pivot_guard_fail("${_rel} is required for Obol window/display-host ownership")
    endif()
  endforeach()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/include/brlobol/window_host.h")
    file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/window_host.h" _window_header)
    foreach(_token
	"BRLObolWindowHost"
	"BRLObolWindowDesc"
	"BRLObolInputProfile"
	"BRLObolInputEvent"
	"BRLObolInputAction"
	"attachController"
	"openFramebuffer"
	"brlobol_window_host_register_display_host"
	"brlobol_window_host_unregister_display_host")
      string(FIND "${_window_header}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("include/brlobol/window_host.h missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/include/brlobol/headless_window_host.h")
    file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/headless_window_host.h" _headless_header)
    foreach(_token
	"BRLObolHeadlessWindowHost"
	"renderPending"
	"setBackgroundColor"
	"setOutputComponents"
	"getLastFrameBuffer"
	"getLastRenderReason")
      string(FIND "${_headless_header}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("include/brlobol/headless_window_host.h missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/include/brlobol/view_controller.h")
    file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/view_controller.h" _view_controller_header)
    foreach(_token
	"BRLObolViewController"
	"getRtViewInfo"
	"setLodService"
	"setLodAutoSubmit"
	"isLodAutoSubmitEnabled"
	"setLodForcedLevel"
	"clearLodForcedLevel"
	"hasLodForcedLevel"
	"getLodForcedLevel"
	"setExactFullDetailBudget"
	"getMaxExactFullDetailFaceCount"
	"getMaxExactFullDetailPointCount"
	"consumeExportSourceFullDetail"
	"consumeMeasureSourceFullDetail"
	"consumeSnapSourceFullDetail"
	"prepareRtPickCaches"
	"getRtPickCache"
	"getRtPickCacheSourceRevision"
	"pickSourceMeshExactRay"
	"pickRtExactRay"
	"clearRtPickCaches"
	"setMeshResidencyBudget"
	"clearMeshResidencyBudget"
	"hasMeshResidencyBudget"
	"getMaxResidentMeshBytes"
	"isMeshResidencyDisplayEvictionEnabled"
	"evictMeshPayloadsToBudget"
	"getLastMeshBudgetInitialResidentBytes"
	"getLastMeshBudgetFinalResidentBytes"
	"getLastMeshBudgetEvictedDisplayMeshCount"
	"hasPendingLodResults"
	"processPendingLodResults"
	"submitLodRequestsIfNeeded"
	"submitLodRequests"
	"applyLodResults"
	"getLodViewRevision"
	"setLodPolicyRevision"
	"getLastLodSubmittedTaskCount"
	"getLastLodAppliedResultCount"
	"getLastLodDiagnostics"
	"struct rt_view_info")
      string(FIND "${_view_controller_header}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("include/brlobol/view_controller.h missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libbrlobol/view_controller.cpp")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/view_controller.cpp" _view_controller_impl)
    foreach(_token
	"BRLObolViewController::getRtViewInfo"
	"source->configureDatabaseSource"
	"BRLObolViewController::setLodService"
	"BRLObolViewController::setLodAutoSubmit"
	"BRLObolViewController::setExactFullDetailBudget"
	"BRLObolViewController::getMaxExactFullDetailFaceCount"
	"BRLObolViewController::getMaxExactFullDetailPointCount"
	"BRLObolViewController::consumeExportSourceFullDetail"
	"BRLObolViewController::consumeMeasureSourceFullDetail"
	"BRLObolViewController::consumeSnapSourceFullDetail"
	"controller_database_source_for_request"
	"controller_source_request_template"
	"SoBRLExportAction"
	"SoBRLMeasureAction"
	"SoBRLSnapAction"
	"BRLObolViewController::prepareRtPickCaches"
	"BRLObolViewController::getRtPickCache"
	"BRLObolViewController::getRtPickCacheSourceRevision"
	"BRLObolViewController::pickSourceMeshExactRay"
	"SoBRLSourceMeshPickAction"
	"drainMatchingResults"
	"requestMatched"
	"submitRequestIndices"
	"brlobol_lod_submit_rt_source_full_detail_request"
	"BRLObolViewController::pickRtExactRay"
	"BRLObolViewController::clearRtPickCaches"
	"BRLObolViewController::setMeshResidencyBudget"
	"BRLObolViewController::clearMeshResidencyBudget"
	"BRLObolViewController::enforceMeshResidencyBudget"
	"BRLObolViewController::evictMeshPayloadsToBudget"
	"SoBRLMeshResidencyAction"
	"lod-memory-budget"
	"BRLObolViewController::submitLodRequestsIfNeeded"
	"BRLObolViewController::processPendingLodResults"
	"BRLObolViewController::submitLodRequests"
	"BRLObolViewController::applyLodResults"
	"lodResultReadyCB"
	"lodResultsPending"
	"lodAutoSubmit"
	"lodViewSignature"
	"lodUseForcedLevel"
	"maxExactFullDetailFaceCount"
	"maxExactFullDetailPointCount"
	"rtPickCaches"
	"rtPickCacheSourceRevisions"
	"meshResidencyBudgetEnabled"
	"meshResidencyEvictDisplayPayloads"
	"lastMeshBudgetInitialResidentBytes"
	"lastMeshBudgetEvictedDisplayMeshCount"
	"advanceLodPolicyRevision"
	"controller_lod_view_signature"
	"controller_lod_source_signature"
	"source->lodBotThreshold.getValue()"
	"stale LoD result revision rejected"
	"SoBRLMeshLodSubmitAction"
	"SoBRLLodUpdateAction"
	"advanceLodViewRevision"
	"lodViewRevision"
	"rt_view_info_init"
	"rt_view_info_sanitize"
	"SoOrthographicCamera"
	"SoPerspectiveCamera")
      string(FIND "${_view_controller_impl}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libbrlobol/view_controller.cpp missing ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[source->[ \t\r\n]*setDatabase]]
	[[source->[ \t\r\n]*path[ \t\r\n]*=]]
	[[source->[ \t\r\n]*drawMode[ \t\r\n]*=]]
	[[source->[ \t\r\n]*sourceRevision[ \t\r\n]*=]]
	[[source->[ \t\r\n]*markStale]])
      string(REGEX MATCH "${_pat}" _hit "${_view_controller_impl}")
      if(_hit)
	_brlobol_pivot_guard_fail("src/libbrlobol/view_controller.cpp must configure controller-owned database sources atomically: ${_hit}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libbrlobol/window_host.cpp")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/window_host.cpp" _window_impl)
    foreach(_token
	"imgstream_fb_display_host_set"
	"imgstream_fb_display_host_clear"
	"BRLObolViewController"
	"SoBRLImageSource"
	"SoBRLViewportImage"
	"syncFromSource"
	"sourceZoom"
	"cursorVisible")
      string(FIND "${_window_impl}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libbrlobol/window_host.cpp missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libbrlobol/headless_window_host.cpp")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/headless_window_host.cpp" _headless_impl)
    foreach(_token
	"BRLObolHeadlessWindowHost"
	"BRLObolViewController"
	"SoOffscreenRenderer"
	"SoOrthographicCamera"
	"SoDB::getContextManager"
	"controller->consumeRenderRequest")
      string(FIND "${_headless_impl}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libbrlobol/headless_window_host.cpp missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_window_host.cpp")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_window_host.cpp" _window_test)
    foreach(_token
	"BRLObolWindowHost"
	"BRLObolInputProfile"
	"brlobol_window_host_register_display_host"
	"imgstream_fb_open"
	"\"/dev/qtgl\""
	"getFramebufferViewportImage"
	"imgstream_fb_view"
	"imgstream_fb_cursor"
	"imgstream_fb_flush")
      string(FIND "${_window_test}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libbrlobol/tests/test_window_host.cpp missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_headless_window_host.cpp")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_headless_window_host.cpp" _headless_test)
    foreach(_token
	"BRLObolHeadlessWindowHost"
	"HeadlessTestContextManager"
	"renderScene"
	"imgstream_fb_open"
	"\"/dev/swrast\""
	"imgstream_fb_poll"
	"getLastFrameBuffer"
	"getLastRenderReason")
      string(FIND "${_headless_test}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libbrlobol/tests/test_headless_window_host.cpp missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libbrlobol/CMakeLists.txt")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/CMakeLists.txt" _brlobol_cmake)
    foreach(_token
	"headless_window_host.cpp"
	"window_host.cpp")
      string(FIND "${_brlobol_cmake}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libbrlobol/CMakeLists.txt missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/CMakeLists.txt")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/CMakeLists.txt" _brlobol_tests_cmake)
    foreach(_token
	"test_brlobol_headless_window_host"
	"test_brlobol_window_host")
      string(FIND "${_brlobol_tests_cmake}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libbrlobol/tests/CMakeLists.txt missing ${_token}")
      endif()
    endforeach()
  endif()

  foreach(_rel IN LISTS _required_files)
    if(NOT EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      continue()
    endif()
    file(READ "${BRLCAD_SOURCE_DIR}/${_rel}" _contents)
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg]]
	[[#[ \t]*include[ \t]*[<"]dm]]
	[[(^|[^A-Za-z0-9_])libbsg([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail("${_rel} reintroduced a legacy display dependency: ${_hit}")
      endif()
    endforeach()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_plot3_ownership)
  set(_libbg_plot3 "${BRLCAD_SOURCE_DIR}/src/libbg/plot3.c")
  if(NOT EXISTS "${_libbg_plot3}")
    _brlobol_pivot_guard_fail(
      "src/libbg/plot3.c is required as the default plot3 implementation unit")
  else()
    file(READ "${_libbg_plot3}" _libbg_plot3_contents)
    string(FIND "${_libbg_plot3_contents}" "#define PLOT3_IMPLEMENTATION"
      _plot3_impl_idx)
    if(_plot3_impl_idx EQUAL -1)
      _brlobol_pivot_guard_fail(
	"src/libbg/plot3.c must define PLOT3_IMPLEMENTATION")
    endif()
  endif()

  set(_libbg_cmake "${BRLCAD_SOURCE_DIR}/src/libbg/CMakeLists.txt")
  if(NOT EXISTS "${_libbg_cmake}")
    _brlobol_pivot_guard_fail(
      "src/libbg/CMakeLists.txt is required for plot3 ownership checks")
  else()
    file(READ "${_libbg_cmake}" _libbg_cmake_contents)
    string(REGEX MATCH
      "set\\([ \t\r\n]*LIBBG_SOURCES[^)]*[ \t\r\n]plot3\\.c([ \t\r\n]|\\))"
      _plot3_source_hit "${_libbg_cmake_contents}")
    if(NOT _plot3_source_hit)
      _brlobol_pivot_guard_fail(
	"src/libbg/CMakeLists.txt must list plot3.c in LIBBG_SOURCES")
    endif()
    string(REGEX MATCH
      "UNITY_BUILD_SKIP[^)]*[ \t\r\n]plot3\\.c([ \t\r\n]|\\))"
      _plot3_unity_hit "${_libbg_cmake_contents}")
    if(NOT _plot3_unity_hit)
      _brlobol_pivot_guard_fail(
	"src/libbg/CMakeLists.txt must keep plot3.c out of unity builds")
    endif()
    string(FIND "${_libbg_cmake_contents}" "PLOT3_DLL_EXPORTS"
      _plot3_exports_idx)
    if(_plot3_exports_idx EQUAL -1)
      _brlobol_pivot_guard_fail(
	"src/libbg/CMakeLists.txt must export default plot3 symbols from libbg")
    endif()
  endif()

  set(_plot3_compat_files
    src/libbsg/vlist.c
    src/libbg/tests/chull.c
    src/libbg/tests/obr.c
  )
  foreach(_rel IN LISTS _plot3_compat_files)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail("${_rel} is required for plot3 ownership checks")
      continue()
    endif()
    file(READ "${_file}" _contents)
    string(FIND "${_contents}" "PLOT3_IMPLEMENTATION" _compat_impl_idx)
    if(NOT _compat_impl_idx EQUAL -1)
      _brlobol_pivot_guard_fail(
	"${_rel} must not define the default plot3 implementation")
    endif()
  endforeach()

  set(_libbsg_cmake "${BRLCAD_SOURCE_DIR}/src/libbsg/CMakeLists.txt")
  if(NOT EXISTS "${_libbsg_cmake}")
    _brlobol_pivot_guard_fail(
      "src/libbsg/CMakeLists.txt is required for plot3 ownership checks")
  else()
    file(READ "${_libbsg_cmake}" _libbsg_cmake_contents)
    string(FIND "${_libbsg_cmake_contents}" "PLOT3_DLL_EXPORTS"
      _bsg_plot3_exports_idx)
    if(NOT _bsg_plot3_exports_idx EQUAL -1)
      _brlobol_pivot_guard_fail(
	"src/libbsg/CMakeLists.txt must not export default plot3 symbols")
    endif()
  endif()

  set(_line_layer_header "${BRLCAD_SOURCE_DIR}/include/bg/line_layer.h")
  set(_line_layer_impl "${BRLCAD_SOURCE_DIR}/src/libbg/line_layer.c")
  if(NOT EXISTS "${_line_layer_header}")
    _brlobol_pivot_guard_fail(
      "include/bg/line_layer.h is required for bg_line_layer plot3 export checks")
  else()
    file(READ "${_line_layer_header}" _line_layer_header_contents)
    foreach(_token
	"bg_line_layer_write_plot3"
	"bg_line_layer_builder_write_plot3")
      string(FIND "${_line_layer_header_contents}" "${_token}" _ll_header_idx)
      if(_ll_header_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "include/bg/line_layer.h must declare ${_token}")
      endif()
    endforeach()
  endif()
  if(NOT EXISTS "${_line_layer_impl}")
    _brlobol_pivot_guard_fail(
      "src/libbg/line_layer.c is required for bg_line_layer plot3 export checks")
  else()
    file(READ "${_line_layer_impl}" _line_layer_impl_contents)
    foreach(_token
	"pl_color"
	"pdv_3move"
	"pdv_3cont"
	"pdv_3point")
      string(FIND "${_line_layer_impl_contents}" "${_token}" _ll_impl_idx)
      if(_ll_impl_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/libbg/line_layer.c bg_line_layer plot3 export must use ${_token}")
      endif()
    endforeach()
  endif()

  set(_plot3_util_cmake "${BRLCAD_SOURCE_DIR}/src/util/CMakeLists.txt")
  if(NOT EXISTS "${_plot3_util_cmake}")
    _brlobol_pivot_guard_fail(
      "src/util/CMakeLists.txt is required for plot3 utility dependency checks")
  else()
    file(READ "${_plot3_util_cmake}" _plot3_util_cmake_contents)
    set(_plot3_util_targets
      asc-plot3
      pixhist3d-plot3
      plot3-asc
      plot3-plot3
      plot3-ps
      plot3color
      plot3getframe
      plot3line2
      plot3rot
      plot3stat
    )
    foreach(_target IN LISTS _plot3_util_targets)
      string(REGEX MATCH "brlcad_addexec\\(${_target}[^)]*\\)"
	_plot3_util_row "${_plot3_util_cmake_contents}")
      if(NOT _plot3_util_row)
	_brlobol_pivot_guard_fail(
	  "src/util/CMakeLists.txt must define plot3 utility target ${_target}")
	continue()
      endif()
      string(FIND "${_plot3_util_row}" "libbsg" _plot3_util_libbsg_idx)
      if(NOT _plot3_util_libbsg_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "plot3 utility target ${_target} must not link libbsg")
      endif()
    endforeach()
  endif()
endfunction()

function(_brlobol_pivot_guard_check_libbg_bsg_neutralization)
  _brlobol_pivot_guard_collect(_libbg_neutral_files
    include/bg
    src/libbg
  )

  foreach(_file IN LISTS _libbg_neutral_files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if(NOT "${_rel}" MATCHES "${_extensions}")
      continue()
    endif()

    file(READ "${_file}" _contents)
    string(REGEX MATCH [[#[ \t]*include[ \t]*[<"]bsg/]]
      _hit "${_contents}")
    if(_hit)
      _brlobol_pivot_guard_fail(
	"${_rel} reintroduced a BSG include into BG-owned code after BG neutralization: ${_hit}")
    endif()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_utility_vlist_neutralization)
  set(_utility_vlist_neutral_files
    src/conv/dxf/dxf-g.c
    src/conv/iges/conv_drawings.c
    src/conv/jack/g-jack.c
    src/conv/off/g-off.c
    src/gtools/ganalyze.cpp
    src/proc-db/tea_nmg.c
  )

  set(_utility_vlist_forbidden
    [[#[ \t]*include[ \t]*[<"]bsg/vlist\.h]]
    [[(^|[^A-Za-z0-9_])bsg_vlist[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])bsg_vlblock[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])bsg_plot_vlblock([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])BSG_(GET_VLIST|FREE_VLIST|ADD_VLIST|CK_VLIST|VLIST_[A-Za-z0-9_]*)([^A-Za-z0-9_]|$)]]
  )

  foreach(_rel IN LISTS _utility_vlist_neutral_files)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for utility vlist neutralization checks")
      continue()
    endif()

    file(READ "${_file}" _contents)
    string(FIND "${_contents}" "bg/vlist.h" _bg_vlist_idx)
    if(_bg_vlist_idx EQUAL -1)
      _brlobol_pivot_guard_fail(
	"${_rel} must include bg/vlist.h after utility vlist neutralization")
    endif()
    foreach(_pat IN LISTS _utility_vlist_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced BSG vlist names after utility vlist neutralization: ${_hit}")
      endif()
    endforeach()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_libanalyze_vlist_neutralization)
  set(_libanalyze_vlist_neutral_files
    src/libanalyze/nirt/nirt.h
    src/libanalyze/nirt/nirt.cpp
  )

  set(_libanalyze_vlist_forbidden
    [[#[ \t]*include[ \t]*[<"]bsg/vlist\.h]]
    [[(^|[^A-Za-z0-9_])bsg_vlist[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])bsg_vlblock[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])bsg_plot_vlblock([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])BSG_(GET_VLIST|FREE_VLIST|ADD_VLIST|CK_VLIST|VLIST_[A-Za-z0-9_]*)([^A-Za-z0-9_]|$)]]
  )

  foreach(_rel IN LISTS _libanalyze_vlist_neutral_files)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for libanalyze vlist neutralization checks")
      continue()
    endif()

    file(READ "${_file}" _contents)
    foreach(_pat IN LISTS _libanalyze_vlist_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced BSG vlist names after BG vlist neutralization: ${_hit}")
      endif()
    endforeach()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_libbrep_bsg_neutralization)
  _brlobol_pivot_guard_collect(_libbrep_neutral_files
    include/brep
    src/libbrep
  )

  set(_libbrep_forbidden
    [[#[ \t]*include[ \t]*[<"]bsg/]]
    [[(^|[^A-Za-z0-9_])bsg_[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])BSG_[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])libbsg([^A-Za-z0-9_]|$)]]
  )

  foreach(_file IN LISTS _libbrep_neutral_files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if(NOT "${_rel}" MATCHES "${_extensions}")
      continue()
    endif()

    file(READ "${_file}" _contents)
    foreach(_pat IN LISTS _libbrep_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced BSG names after libbrep BSG neutralization: ${_hit}")
      endif()
    endforeach()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_libnmg_vlist_neutralization)
  set(_libnmg_vlist_neutral_files
    src/libnmg/nmg_private.h
    src/libnmg/eval.c
    src/libnmg/fcut.c
    src/libnmg/plot.c
    src/libnmg/misc.c
  )

  set(_libnmg_vlist_forbidden
    [[#[ \t]*include[ \t]*[<"]bsg/vlist\.h]]
    [[(^|[^A-Za-z0-9_])bsg_vlist[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])bsg_vlblock[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])bsg_plot_vlblock([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])BSG_(GET_VLIST|FREE_VLIST|ADD_VLIST|CK_VLIST|VLIST_[A-Za-z0-9_]*)([^A-Za-z0-9_]|$)]]
  )

  foreach(_rel IN LISTS _libnmg_vlist_neutral_files)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for libnmg vlist neutralization checks")
      continue()
    endif()

    file(READ "${_file}" _contents)
    foreach(_pat IN LISTS _libnmg_vlist_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced BSG vlist names after BG vlist neutralization: ${_hit}")
      endif()
    endforeach()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_librt_vlist_neutralization)
  set(_librt_vlist_neutral_files
    include/rt/nmg_conv.h
    src/librt/prep.cpp
    src/librt/vlist.c
    src/librt/primitives/annot/annot.c
    src/librt/primitives/arb8/arb8.c
    src/librt/primitives/arbn/arbn.c
    src/librt/primitives/ars/ars.c
    src/librt/primitives/bot/bot_plot.cpp
    src/librt/primitives/brep/brep.cpp
    src/librt/primitives/brep/brep_debug.h
    src/librt/primitives/bspline/bspline.cpp
    src/librt/primitives/cline/cline.c
    src/librt/primitives/datum/datum.c
    src/librt/primitives/dsp/dsp.c
    src/librt/primitives/ebm/ebm.c
    src/librt/primitives/ehy/ehy.c
    src/librt/primitives/ell/ell.c
    src/librt/primitives/epa/epa.c
    src/librt/primitives/eto/eto.c
    src/librt/primitives/extrude/extrude.c
    src/librt/primitives/grip/grip.c
    src/librt/primitives/half/half.c
    src/librt/primitives/hf/hf.c
    src/librt/primitives/hrt/hrt.c
    src/librt/primitives/hyp/hyp.c
    src/librt/primitives/joint/joint.c
    src/librt/primitives/metaball/metaball.c
    src/librt/primitives/part/part.c
    src/librt/primitives/pipe/pipe.c
    src/librt/primitives/pnts/pnts.c
    src/librt/primitives/poly/poly.c
    src/librt/primitives/primitive_util.c
    src/librt/primitives/revolve/revolve.c
    src/librt/primitives/rhc/rhc.c
    src/librt/primitives/rpc/rpc.c
    src/librt/primitives/sketch/sketch.c
    src/librt/primitives/superell/superell.c
    src/librt/primitives/tgc/tgc.c
    src/librt/primitives/tor/tor.c
    src/librt/primitives/vol/vol.c
    src/librt/tests/edit/nmg.cpp
  )

  set(_librt_vlist_forbidden
    [[#[ \t]*include[ \t]*[<"]bsg/vlist\.h]]
    [[(^|[^A-Za-z0-9_])bsg_vlist[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])bsg_vlblock[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])BSG_(GET_VLIST|FREE_VLIST|ADD_VLIST|CK_VLIST|VLIST_[A-Za-z0-9_]*)([^A-Za-z0-9_]|$)]]
  )

  foreach(_rel IN LISTS _librt_vlist_neutral_files)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for librt vlist neutralization checks")
      continue()
    endif()

    file(READ "${_file}" _contents)
    foreach(_pat IN LISTS _librt_vlist_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced BSG vlist names after RT vlist neutralization: ${_hit}")
      endif()
    endforeach()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_librt_primitive_line_command_neutralization)
  file(GLOB_RECURSE _librt_primitive_line_command_files LIST_DIRECTORIES false
    "${BRLCAD_SOURCE_DIR}/src/librt/primitives/*.c"
    "${BRLCAD_SOURCE_DIR}/src/librt/primitives/*.cpp"
    "${BRLCAD_SOURCE_DIR}/src/librt/primitives/*.h"
  )
  list(APPEND _librt_primitive_line_command_files
    "${BRLCAD_SOURCE_DIR}/include/rt/functab.h")

  set(_librt_primitive_line_command_forbidden
    [[#[ \t]*include[ \t]*[<"]bsg/geometry\.h]]
    [[(^|[^A-Za-z0-9_])BSG_GEOMETRY_[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
  )

  foreach(_file IN LISTS _librt_primitive_line_command_files)
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_file} is required for librt primitive line-command neutralization checks")
      continue()
    endif()

    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    file(READ "${_file}" _contents)
    foreach(_pat IN LISTS _librt_primitive_line_command_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced BSG primitive line-command names after RT neutralization: ${_hit}")
      endif()
    endforeach()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_librt_bsg_umbrella_neutralization)
  file(GLOB_RECURSE _librt_bsg_umbrella_files LIST_DIRECTORIES false
    "${BRLCAD_SOURCE_DIR}/include/rt/*.h"
    "${BRLCAD_SOURCE_DIR}/src/librt/*.c"
    "${BRLCAD_SOURCE_DIR}/src/librt/*.cpp"
    "${BRLCAD_SOURCE_DIR}/src/librt/*.h"
  )

  foreach(_file IN LISTS _librt_bsg_umbrella_files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if("${_rel}" MATCHES "^src/librt/tests/")
      continue()
    endif()

    file(READ "${_file}" _contents)
    string(REGEX MATCH "#[ \t]*include[ \t]*[<\"]bsg\\.h[>\"]"
      _hit "${_contents}")
    if(_hit)
      _brlobol_pivot_guard_fail(
	"${_rel} reintroduced the BSG umbrella header after librt include neutralization: ${_hit}")
    endif()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_librt_edit_knob_neutralization)
  set(_edit_header "${BRLCAD_SOURCE_DIR}/include/rt/edit.h")
  if(NOT EXISTS "${_edit_header}")
    _brlobol_pivot_guard_fail(
      "include/rt/edit.h is required for librt edit knob neutralization checks")
  else()
    file(READ "${_edit_header}" _edit_header_contents)
    foreach(_token
	"struct rt_edit_knobs"
	"rt_edit_knobs_reset"
	"rt_edit_knobs_hash"
	"struct rt_edit_view"
	"rt_edit_set_view")
      string(FIND "${_edit_header_contents}" "${_token}" _token_idx)
      if(_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "include/rt/edit.h must expose RT-owned edit state via ${_token}")
      endif()
    endforeach()

    set(_edit_header_forbidden
      [[#[ \t]*include[ \t]*[<"]bsg/defines\.h]]
      [[(^|[^A-Za-z0-9_])bsg_view_knobs([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])bsg_knobs_reset([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])bsg_knobs_hash([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])BSG_KNOBS_[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])bsg_view([^A-Za-z0-9_]|$)]]
      [[struct[ \t\r\n]+bsg_view[ \t\r\n]*\*[ \t\r\n]*vp]]
    )
    foreach(_pat IN LISTS _edit_header_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_edit_header_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "include/rt/edit.h reintroduced BSG-owned edit state: ${_hit}")
      endif()
    endforeach()
  endif()

  set(_edit_legacy_bsg_header
    "${BRLCAD_SOURCE_DIR}/include/rt/edit_legacy_bsg.h")
  if(NOT EXISTS "${_edit_legacy_bsg_header}")
    _brlobol_pivot_guard_fail(
      "include/rt/edit_legacy_bsg.h is required for transitional BSG edit compatibility")
  else()
    file(READ "${_edit_legacy_bsg_header}" _edit_legacy_bsg_contents)
    foreach(_token
	"rt_edit_view_from_bsg"
	"rt_edit_create_bsg"
	"rt_edit_reinit_bsg"
	"rt_edit_knob_cmd_process_bsg")
      string(FIND "${_edit_legacy_bsg_contents}" "${_token}" _legacy_token_idx)
      if(_legacy_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "include/rt/edit_legacy_bsg.h must expose transitional edit adapter ${_token}")
      endif()
    endforeach()
  endif()

  set(_edit_impl "${BRLCAD_SOURCE_DIR}/src/librt/edit.cpp")
  if(NOT EXISTS "${_edit_impl}")
    _brlobol_pivot_guard_fail(
      "src/librt/edit.cpp is required for librt edit knob neutralization checks")
  else()
    file(READ "${_edit_impl}" _edit_impl_contents)
    set(_edit_impl_forbidden
      [[#[ \t]*include[ \t]*[<"]bsg/defines\.h]]
      [[#[ \t]*include[ \t]*[<"]bsg/util\.h]]
      [[#[ \t]*include[ \t]*[<"]bsg/view_state\.h]]
      [[#[ \t]*include[ \t]*[<"]rt/edit_legacy_bsg\.h]]
      [[(^|[^A-Za-z0-9_])bsg_knobs_reset([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])bsg_knobs_hash([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])BSG_KNOBS_[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])bsg_view([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])rt_edit_view_from_bsg([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])rt_edit_create_bsg([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])rt_edit_reinit_bsg([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])rt_edit_knob_cmd_process_bsg([^A-Za-z0-9_]|$)]]
    )
    foreach(_pat IN LISTS _edit_impl_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_edit_impl_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "src/librt/edit.cpp reintroduced BSG edit adapter coupling: ${_hit}")
      endif()
    endforeach()
  endif()

  set(_edit_test_view_helper
    "${BRLCAD_SOURCE_DIR}/src/librt/tests/edit/edit_test_view.h")
  if(NOT EXISTS "${_edit_test_view_helper}")
    _brlobol_pivot_guard_fail(
      "src/librt/tests/edit/edit_test_view.h is required for neutral edit view test coverage")
  else()
    file(READ "${_edit_test_view_helper}" _edit_test_view_helper_contents)
    foreach(_token
	rt_edit_view
	rt_edit_test_view_update
	rt_edit_test_view_set_aet
	rt_edit_test_view_set_size
	rt_edit_test_view_init_size
	rt_edit_test_view_init_identity_size
	rt_edit_test_view_init
	"rt/view.h")
      string(FIND "${_edit_test_view_helper_contents}" "${_token}"
	_edit_test_view_helper_token_idx)
      if(_edit_test_view_helper_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/librt/tests/edit/edit_test_view.h must expose neutral edit view test helper ${_token}")
      endif()
    endforeach()
  endif()

  foreach(_rel
      src/librt/tests/edit/arbn.cpp
      src/librt/tests/edit/annot.cpp
      src/librt/tests/edit/datum.cpp
      src/librt/tests/edit/revolve.cpp
      src/librt/tests/edit/ell.cpp
      src/librt/tests/edit/tor.cpp
      src/librt/tests/edit/tgc.cpp
      src/librt/tests/edit/epa.cpp
      src/librt/tests/edit/ehy.cpp
      src/librt/tests/edit/eto.cpp
      src/librt/tests/edit/hyp.cpp
      src/librt/tests/edit/rpc.cpp
      src/librt/tests/edit/rhc.cpp
      src/librt/tests/edit/part.cpp
      src/librt/tests/edit/superell.cpp
      src/librt/tests/edit/cline.cpp
      src/librt/tests/edit/extrude.cpp
      src/librt/tests/edit/metaball.cpp
      src/librt/tests/edit/hrt.cpp
      src/librt/tests/edit/ars.cpp
      src/librt/tests/edit/pipe.cpp
      src/librt/tests/edit/dsp.cpp
      src/librt/tests/edit/vol.cpp
      src/librt/tests/edit/ebm.cpp
      src/librt/tests/edit/sph.cpp
      src/librt/tests/edit/arb8.cpp
      src/librt/tests/edit/bot.cpp
      src/librt/tests/edit/sketch.cpp
      src/librt/tests/edit/bspline.cpp
      src/librt/tests/edit/brep.cpp
      src/librt/tests/edit/nmg.cpp
      src/librt/tests/edit/comb.cpp)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for neutral edit view test coverage")
      continue()
    endif()
    file(READ "${_file}" _contents)
    foreach(_token
	"edit_test_view.h"
	"rt_edit_test_view_init"
	"rt_edit_create(&fp")
      string(FIND "${_contents}" "${_token}" _edit_test_token_idx)
      if(_edit_test_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "${_rel} must exercise RT-owned edit views through ${_token}")
      endif()
    endforeach()
    set(_edit_test_forbidden
      [[#[ \t]*include[ \t]*[<"]bsg/util\.h]]
      [[#[ \t]*include[ \t]*[<"]rt/edit_legacy_bsg\.h]]
      [[(^|[^A-Za-z0-9_])rt_edit_create_bsg([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])struct[ \t\r\n]+bsg_view([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])bsg_init([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])bsg_mat_aet([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])bsg_update([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])INV_BV([^A-Za-z0-9_]|$)]]
    )
    foreach(_pat IN LISTS _edit_test_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced BSG edit view test setup: ${_hit}")
      endif()
    endforeach()
  endforeach()

  set(_edit_legacy_bsg_impl
    "${BRLCAD_SOURCE_DIR}/src/librt/edit_legacy_bsg.cpp")
  if(NOT EXISTS "${_edit_legacy_bsg_impl}")
    _brlobol_pivot_guard_fail(
      "src/librt/edit_legacy_bsg.cpp is required for transitional BSG edit adapters")
  else()
    file(READ "${_edit_legacy_bsg_impl}" _edit_legacy_bsg_impl_contents)
    foreach(_token
	"rt/edit_legacy_bsg.h"
	"rt/view_legacy_bsg.h"
	"rt_edit_view_from_bsg"
	"rt_edit_create_bsg"
	"rt_edit_reinit_bsg"
	"rt_edit_knob_cmd_process_bsg"
	"rt_view_scale_from_bsg"
	"rt_view_base2local_from_bsg"
	"rt_view_local2base_from_bsg"
	"rt_view_coord_from_bsg"
	"rt_view_rotate_about_from_bsg"
	"rt_view_rotation_from_bsg"
	"rt_view_center_from_bsg"
	"rt_view_model2view_from_bsg"
	"rt_view_view2model_from_bsg")
      string(FIND "${_edit_legacy_bsg_impl_contents}" "${_token}" _legacy_impl_idx)
      if(_legacy_impl_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/librt/edit_legacy_bsg.cpp must own transitional edit adapter ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg/view_state\.h]]
	[[(^|[^A-Za-z0-9_])v->[ \t\r\n]*gv_(scale|base2local|local2base|coord|rotate_about|rotation|center|model2view|view2model)]])
      string(REGEX MATCH "${_pat}" _edit_legacy_bsg_direct_hit
	"${_edit_legacy_bsg_impl_contents}")
      if(_edit_legacy_bsg_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/librt/edit_legacy_bsg.cpp must route BSG edit-view snapshots through rt/view_legacy_bsg.h: ${_edit_legacy_bsg_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_librt_cmake "${BRLCAD_SOURCE_DIR}/src/librt/CMakeLists.txt")
  if(NOT EXISTS "${_librt_cmake}")
    _brlobol_pivot_guard_fail(
      "src/librt/CMakeLists.txt is required for librt edit adapter registration checks")
  else()
    file(READ "${_librt_cmake}" _librt_cmake_contents)
    string(REGEX MATCH
      "set\\([ \t\r\n]*LIBRT_SOURCES[^)]*[ \t\r\n]edit_legacy_bsg\\.cpp([ \t\r\n]|\\))"
      _edit_legacy_bsg_source_hit "${_librt_cmake_contents}")
    if(NOT _edit_legacy_bsg_source_hit)
      _brlobol_pivot_guard_fail(
	"src/librt/CMakeLists.txt must list edit_legacy_bsg.cpp in LIBRT_SOURCES")
    endif()
  endif()

  set(_mged_chgview "${BRLCAD_SOURCE_DIR}/src/mged/chgview.c")
  if(EXISTS "${_mged_chgview}")
    file(READ "${_mged_chgview}" _mged_chgview_contents)
    string(REGEX MATCH
      "bsg_knobs_reset[ \t\r\n]*\\([ \t\r\n]*&[ \t\r\n]*MEDIT"
      _mged_edit_reset_hit "${_mged_chgview_contents}")
    if(_mged_edit_reset_hit)
      _brlobol_pivot_guard_fail(
	"src/mged/chgview.c must use rt_edit_knobs_reset for MEDIT(s)->k")
    endif()
  endif()

  set(_mged_cmd "${BRLCAD_SOURCE_DIR}/src/mged/cmd.cpp")
  if(EXISTS "${_mged_cmd}")
    file(READ "${_mged_cmd}" _mged_cmd_contents)
    string(REGEX MATCH
      "bsg_knobs_hash[ \t\r\n]*\\([ \t\r\n]*&[ \t\r\n]*e->[ \t\r\n]*k"
      _mged_edit_hash_hit "${_mged_cmd_contents}")
    if(_mged_edit_hash_hit)
      _brlobol_pivot_guard_fail(
	"src/mged/cmd.cpp must use rt_edit_knobs_hash for rt_edit::k")
    endif()
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_orientation_quat_from_bsg]]
	[[rt_view_eye_pos_from_bsg]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_size_from_bsg]])
      string(REGEX MATCH "${_token}" _mged_cmd_view_arg_token_hit
	"${_mged_cmd_contents}")
      if(NOT _mged_cmd_view_arg_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/cmd.cpp must route simulation view argument injection through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[bv->[ \t\r\n]*gv_rotation]]
	[[bv->[ \t\r\n]*gv_eye_pos]]
	[[bv->[ \t\r\n]*gv_size]]
	[[vrp->[ \t\r\n]*vr_scale[ \t\r\n]*=[^\n;]*gv_scale]])
      string(REGEX MATCH "${_pat}" _mged_cmd_view_arg_direct_hit
	"${_mged_cmd_contents}")
      if(_mged_cmd_view_arg_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/cmd.cpp reintroduced direct BSG simulation view argument reads: ${_mged_cmd_view_arg_direct_hit}")
      endif()
    endforeach()
    foreach(_token
	[[rt_view_width_from_bsg]]
	[[rt_view_height_from_bsg]]
	[[rt_view_dimensions_set_bsg]]
	[[rt_view_scale_state_set_bsg]]
	[[rt_view_initial_scale_from_bsg]]
	[[rt_view_absolute_scale_from_bsg]]
	[[rt_view_local2base_from_bsg]]
	[[rt_view_base2local_from_bsg]]
	[[rt_view_unit_conversion_set_bsg]]
	[[rt_view_inverse_size_from_bsg]]
	[[rt_view_perspective_from_bsg]]
	[[rt_view_perspective_set_bsg]]
	[[rt_view_coord_from_bsg]]
	[[rt_view_coord_set_bsg]]
	[[rt_view_rotate_about_from_bsg]]
	[[rt_view_rotate_about_set_bsg]]
	[[rt_view_keypoint_from_bsg]]
	[[rt_view_keypoint_set_bsg]]
	[[rt_view_aet_from_bsg]]
	[[rt_view_aet_state_set_bsg]]
	[[rt_view_rotation_from_bsg]]
	[[rt_view_rotation_set_bsg]]
	[[rt_view_center_from_bsg]]
	[[rt_view_center_vec_set_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_model2view_set_bsg]]
	[[rt_view_view2model_from_bsg]]
	[[rt_view_view2model_set_bsg]]
	[[rt_view_pmodel2view_from_bsg]]
	[[rt_view_pmodel2view_set_bsg]]
	[[rt_view_pmat_from_bsg]]
	[[rt_view_pmat_set_bsg]]
	[[rt_view_eye_pos_set_bsg]])
      string(REGEX MATCH "${_token}" _mged_cmd_cache_token_hit
	"${_mged_cmd_contents}")
      if(NOT _mged_cmd_cache_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/cmd.cpp must route adapter-covered staging/cache/hash view reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(c|dst)->gv_(scale|i_scale|a_scale|local2base|base2local|size|isize|perspective|coord|rotate_about)[ \t]*=[^\n;]*(v|src)->gv_]]
	[[staging->gv_(width|height)[ \t]*=[^\n;]*mged_view->gv_]]
	[[gvp->[ \t\r\n]*gv_(width|height)[ \t\r\n]*=]]
	[[staging->[ \t\r\n]*gv_(width|height)[ \t\r\n]*=]]
	[[staging->[ \t\r\n]*gv_(local2base|base2local)[ \t\r\n]*=]]
	[[v->[ \t\r\n]*gv_(perspective|coord|rotate_about)[ \t\r\n]*=]]
	[[dst->[ \t\r\n]*gv_(perspective|coord|rotate_about)[ \t\r\n]*=]]
	[[v->[ \t\r\n]*gv_(scale|i_scale|a_scale|size|isize|local2base|base2local)[ \t\r\n]*=]]
	[[dst->[ \t\r\n]*gv_(scale|i_scale|a_scale|size|isize|local2base|base2local)[ \t\r\n]*=]]
	[[VMOVE[^\n;]*v->gv_(eye_pos|keypoint)[^\n;]*c->gv_]]
	[[VMOVE[^\n;]*v->gv_aet[^\n;]*c->gv_]]
	[[MAT_COPY[^\n;]*v->gv_(rotation|center|pmat)[^\n;]*c->gv_]]
	[[MAT_COPY[^\n;]*v->gv_(model2view|view2model|pmodel2view)[^\n;]*c->gv_]]
	[[rt_view_(eye_pos|keypoint)_from_bsg[^\n;]*dst->gv_]]
	[[rt_view_aet_from_bsg[^\n;]*dst->gv_]]
	[[rt_view_(rotation|center|pmat)_from_bsg[^\n;]*dst->gv_]]
	[[rt_view_(model2view|view2model|pmodel2view)_from_bsg[^\n;]*dst->gv_]]
	[[VMOVE[^\n;]*(c|dst)->gv_(eye_pos|keypoint|aet)[^\n;]*(v|src)->gv_]]
	[[MAT_COPY[^\n;]*(c|dst)->gv_(rotation|center|model2view|view2model|pmodel2view|pmat)[^\n;]*(v|src)->gv_]]
	[[bu_data_hash_update[^\n;]*&v->gv_(scale|i_scale|a_scale|size|isize|perspective|center|rotation|model2view|view2model|pmodel2view|pmat|eye_pos|keypoint|aet|coord|rotate_about)]])
      string(REGEX MATCH "${_pat}" _mged_cmd_cache_direct_hit
	"${_mged_cmd_contents}")
      if(_mged_cmd_cache_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/cmd.cpp reintroduced direct BSG staging/cache/hash view reads: ${_mged_cmd_cache_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_edit "${BRLCAD_SOURCE_DIR}/src/libged/edit/edit.cpp")
  if(EXISTS "${_libged_edit}")
    file(READ "${_libged_edit}" _libged_edit_contents)
    string(REGEX MATCH "s->[ \t\r\n]*vp[ \t\r\n]*=" _libged_edit_vp_hit
      "${_libged_edit_contents}")
    if(_libged_edit_vp_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/edit/edit.cpp must install edit views via rt_edit_set_view, not direct s->vp assignment")
    endif()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_scale[ \t\r\n]*=]]
      _libged_edit_scale_hit "${_libged_edit_contents}")
    if(_libged_edit_scale_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/edit/edit.cpp reintroduced direct edit-view scale writes: ${_libged_edit_scale_hit}")
    endif()
  endif()
endfunction()

function(_brlobol_pivot_guard_check_librt_view_info_neutralization)
  foreach(_rel
      include/rt/view.h
      include/rt/functab.h
      include/rt/primitives/bot.h
      include/rt/primitives/brep.h
      include/rt/primitives/pg.h)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for RT view-info neutralization checks")
      continue()
    endif()

    file(READ "${_file}" _contents)
    if("${_rel}" STREQUAL "include/rt/view.h")
      foreach(_token
	  "struct rt_view_info"
	  "struct rt_view_lod_settings"
	  "enum rt_view_lod_policy_mode"
	  "struct rt_view_lod_policy"
	  "RT_VIEW_LOD_POLICY_INIT"
	  "RT_VIEW_MIN_SIZE"
	  "RT_VIEW_MIN_SCALE"
	  "rt_view_info_init"
	  "rt_view_info_sanitize"
	  "rt_view_lod_policy_init"
	  "rt_view_lod_policy_sanitize"
	  "rt_view_lod_curve_scale"
	  "rt_view_lod_bot_threshold"
	  "rt_view_avg_sample_spacing"
	  "rt_view_solid_point_spacing"
	  "struct rt_mesh_lod"
	  "struct rt_mesh_lod_data"
	  "struct rt_mesh_lod_detail"
	  "rt_mesh_lod_detail_setup_callback"
	  "struct rt_mesh_lod_info"
	  "struct rt_mesh_lod_cache_status"
	  "rt_mesh_lod_cache_clear_all"
	  "rt_mesh_lod_cache_status_init"
	  "db_mesh_lod_status"
	  "db_mesh_lod_refresh"
	  "db_mesh_lod_invalidate"
	  "db_mesh_lod_store_mesh"
	  "faces order, i.e. face_count * 3 normals"
	  "rt_mesh_lod_load_level"
	  "rt_mesh_lod_load_view"
	  "rt_mesh_lod_current_level"
	  "rt_mesh_lod_has_active_data"
	  "rt_mesh_lod_data_get"
	  "rt_mesh_lod_info_get"
	  "rt_mesh_lod_detail_init"
	  "rt_mesh_lod_detail_callbacks_set"
	  "rt_mesh_lod_detail_callbacks_clear"
	  "rt_mesh_lod_memshrink"
	  "rt_mesh_lod_destroy")
	string(FIND "${_contents}" "${_token}" _token_idx)
	if(_token_idx EQUAL -1)
	  _brlobol_pivot_guard_fail(
	    "include/rt/view.h must expose neutral RT view state via ${_token}")
	endif()
      endforeach()
    endif()
    if("${_rel}" STREQUAL "include/rt/functab.h")
      string(FIND "${_contents}" "const struct rt_view_info *" _view_info_idx)
      if(_view_info_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "include/rt/functab.h must use rt_view_info for primitive plot/LoD/face-set callbacks")
      endif()
    endif()

    set(_view_info_forbidden
      [[#[ \t]*include[ \t]*[<"]bsg/]]
      [[(^|[^A-Za-z0-9_])bsg_view([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])bsg_mesh_lod([^A-Za-z0-9_]|$)]]
      [[struct[ \t\r\n]+bsg_view]]
      [[struct[ \t\r\n]+bsg_mesh_lod]]
      [[const[ \t\r\n]+struct[ \t\r\n]+bsg_view[ \t\r\n]*\*]]
    )
    foreach(_pat IN LISTS _view_info_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced BSG view state into the neutral RT view-info API: ${_hit}")
      endif()
    endforeach()
  endforeach()

  set(_view_legacy_bsg_header
    "${BRLCAD_SOURCE_DIR}/include/rt/view_legacy_bsg.h")
  if(NOT EXISTS "${_view_legacy_bsg_header}")
    _brlobol_pivot_guard_fail(
      "include/rt/view_legacy_bsg.h is required for transitional BSG view compatibility")
  else()
    file(READ "${_view_legacy_bsg_header}" _view_legacy_bsg_contents)
    foreach(_token
	rt_view_info_from_bsg
	rt_view_width_from_bsg
	rt_view_height_from_bsg
	rt_view_dimensions_set_bsg
	rt_view_orientation_quat_from_bsg
	rt_view_aet_from_bsg
	rt_view_aet_set_bsg
	rt_view_aet_state_set_bsg
	rt_view_perspective_from_bsg
	rt_view_perspective_set_bsg
	rt_view_model2view_from_bsg
	rt_view_model2view_set_bsg
	rt_view_view2model_from_bsg
	rt_view_view2model_set_bsg
	rt_view_pmodel2view_from_bsg
	rt_view_pmodel2view_set_bsg
	rt_view_pmat_from_bsg
	rt_view_pmat_set_bsg
	rt_view_rotation_from_bsg
	rt_view_rotation_set_bsg
	rt_view_center_from_bsg
	rt_view_center_vec_set_bsg
	rt_view_plane_from_bsg
	rt_view_lod_policy_from_bsg
	rt_view_lod_policy_apply_bsg
	rt_view_lod_policy_copy_bsg
	rt_view_lod_bounds_update_bsg
	rt_view_lod_bounds_callback_set_bsg
	rt_view_lod_bounds_callback_is_bsg
	rt_view_bounds_update_callback_bsg_t
	rt_view_bounds_update_callback_from_bsg
	rt_view_bounds_update_callback_set_bsg
	rt_view_bounds_update_callback_call_bsg
	rt_view_scale_from_bsg
	rt_view_scale_set_bsg
	rt_view_scale_storage_from_bsg
	rt_view_scale_state_set_bsg
	rt_view_initial_scale_from_bsg
	rt_view_initial_scale_set_bsg
	rt_view_absolute_scale_from_bsg
	rt_view_absolute_scale_set_bsg
	rt_view_local2base_from_bsg
	rt_view_base2local_from_bsg
	rt_view_unit_conversion_set_bsg
	rt_view_width_from_bsg
	rt_view_height_from_bsg
	rt_view_radius_from_bsg
	rt_view_dimensions_set_bsg
	rt_view_screen_to_view_from_bsg
	rt_view_screen_point_from_bsg
	rt_view_interactive_rect_from_bsg
	rt_view_interactive_rect_set_bsg
	rt_view_adc_from_bsg
	rt_view_adc_set_bsg
	rt_view_grid_from_bsg
	rt_view_grid_set_bsg
	rt_view_model_axes_from_bsg
	rt_view_model_axes_set_bsg
	rt_view_view_axes_from_bsg
	rt_view_view_axes_set_bsg
	rt_view_center_dot_from_bsg
	rt_view_center_dot_set_bsg
	rt_view_scale_overlay_from_bsg
	rt_view_scale_overlay_set_bsg
	rt_view_params_from_bsg
	rt_view_params_set_bsg
	rt_view_refresh_request_bsg
	rt_view_refresh_dirty_from_bsg
	rt_view_refresh_consume_bsg
	rt_view_refresh_complete_bsg
	rt_view_refresh_enabled_from_bsg
	rt_view_refresh_enabled_set_bsg
	rt_view_refresh_suppressed_from_bsg
	rt_view_refresh_suppress_begin_bsg
	rt_view_refresh_suppress_end_bsg
	rt_view_refresh_drawn_count_from_bsg
	rt_view_refresh_drawn_count_set_bsg
	rt_view_size_from_bsg
	rt_view_size_set_bsg
	rt_view_inverse_size_from_bsg
	rt_view_eye_pos_from_bsg
	rt_view_eye_pos_set_bsg
	rt_view_keypoint_from_bsg
	rt_view_keypoint_set_bsg
	rt_view_rotate_about_from_bsg
	rt_view_rotate_about_set_bsg
	rt_view_coord_from_bsg
	rt_view_coord_set_bsg
	rt_view_snap_lines_from_bsg
	rt_view_snap_lines_set_bsg
	rt_view_snap_source_flags_from_bsg
	rt_view_snap_source_flags_set_bsg
	rt_view_snap_kind_mask_from_bsg
	rt_view_prepare_tcl_snap_bsg
	rt_view_center_linesnap_bsg
	rt_view_zclip_from_bsg
	rt_view_zclip_set_bsg
	rt_view_framebuffer_mode_from_bsg
	rt_view_framebuffer_mode_set_bsg
	rt_view_cleared_from_bsg
	rt_view_cleared_set_bsg
	rt_view_settings_shared_bsg
	rt_view_snap_tolerance_factor_from_bsg
	rt_view_snap_tolerance_factor_set_bsg
	rt_mesh_lod_load_view_scene_ref_bsg
	rt_mesh_lod_free_scene_ref_bsg)
      string(FIND "${_view_legacy_bsg_contents}" "${_token}"
	_view_legacy_adapter_idx)
      if(_view_legacy_adapter_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "include/rt/view_legacy_bsg.h must expose transitional ${_token} adapter")
      endif()
    endforeach()
  endif()

  set(_librt_view_impl "${BRLCAD_SOURCE_DIR}/src/librt/view.c")
  if(EXISTS "${_librt_view_impl}")
    file(READ "${_librt_view_impl}" _librt_view_impl_contents)
    foreach(_token
	rt_view_info_init
	rt_view_info_sanitize
	rt_view_lod_policy_init
	rt_view_lod_policy_sanitize
	rt_view_lod_curve_scale
	rt_view_lod_bot_threshold
	rt_view_avg_sample_spacing
	rt_view_solid_point_spacing)
      string(FIND "${_librt_view_impl_contents}" "${_token}"
	_librt_view_impl_token_idx)
      if(_librt_view_impl_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/librt/view.c must keep neutral RT view helper ${_token}")
      endif()
    endforeach()
    set(_librt_view_impl_forbidden
      [[#[ \t]*include[ \t]*[<"]bsg/view_state\.h]]
      [[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
      [[(^|[^A-Za-z0-9_])rt_view_info_from_bsg([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])bsg_view([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])primitive_lod_curve_scale([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])primitive_lod_bot_threshold([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])view_avg_sample_spacing([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])solid_point_spacing([^A-Za-z0-9_]|$)]]
    )
    foreach(_pat IN LISTS _librt_view_impl_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_librt_view_impl_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "src/librt/view.c reintroduced BSG view adapter logic into neutral RT view helpers: ${_hit}")
      endif()
    endforeach()
  endif()

  set(_librt_cache_lod_impl "${BRLCAD_SOURCE_DIR}/src/librt/cache_lod.c")
  if(EXISTS "${_librt_cache_lod_impl}")
    file(READ "${_librt_cache_lod_impl}" _librt_cache_lod_impl_contents)
    foreach(_token
	rt_mesh_lod_context_ensure
	rt_mesh_lod_cache_clear_all
	rt_mesh_lod_cache_status_init
	db_mesh_lod_status
	db_mesh_lod_refresh
	db_mesh_lod_invalidate
	db_mesh_lod_store_mesh
	db_mesh_lod_get
	rt_mesh_lod_load_level
	rt_mesh_lod_load_view
	rt_mesh_lod_current_level
	rt_mesh_lod_has_active_data
	rt_mesh_lod_data_get
	rt_mesh_lod_info_get
	rt_mesh_lod_mesh_arrays_validate
	rt_mesh_lod_detail_has_payload
	rt_mesh_lod_detail_init
	rt_mesh_lod_detail_callbacks_set
	rt_mesh_lod_detail_callbacks_clear
	rt_mesh_lod_detail_setup_bsg)
      string(FIND "${_librt_cache_lod_impl_contents}" "${_token}"
	_librt_cache_lod_impl_token_idx)
      if(_librt_cache_lod_impl_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	"src/librt/cache_lod.c must keep neutral RT mesh-LoD helper ${_token}")
      endif()
    endforeach()
    foreach(_token
	"bsg->points_orig && bsg->porig_cnt > 0"
	"data->points_orig && data->point_orig_count"
	"info->has_faces && info->has_points && info->has_original_points")
      string(FIND "${_librt_cache_lod_impl_contents}" "${_token}"
	_librt_cache_lod_active_data_token_idx)
      if(_librt_cache_lod_active_data_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/librt/cache_lod.c must require original points for active RT mesh-LoD data token ${_token}")
      endif()
    endforeach()
    foreach(_token
	"detail.normal_count != index_count"
	"!detail.normals && detail.normal_count != 0"
	"bot->faces, bot->num_faces"
	"invalid BoT mesh arrays"
	"!faces || fcnt <= 0"
	"!points || pcnt <= 0"
	"!points_orig || porig_cnt <= 0"
	"faces[i] < 0 || faces[i] >= pcnt"
	"faces[i] >= porig_cnt"
	"data->normals = bsg->normals"
	"data->normal_count = (bsg->normals && data->face_count)"
	"info->normal_count = (bsg->normals && info->face_count)"
	"info->has_normals = (bsg->normals && info->normal_count)")
      string(FIND "${_librt_cache_lod_impl_contents}" "${_token}"
	_librt_cache_lod_normal_token_idx)
      if(_librt_cache_lod_normal_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/librt/cache_lod.c must keep RT mesh-LoD normal metadata token ${_token}")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^_A-Za-z0-9])rt_mesh_lod_bsg[ \t\r\n]*\(]]
      _librt_cache_lod_public_bsg_hit "${_librt_cache_lod_impl_contents}")
    if(_librt_cache_lod_public_bsg_hit)
      _brlobol_pivot_guard_fail(
	"src/librt/cache_lod.c must keep public BSG extraction behind view_legacy_bsg.c")
    endif()
  endif()

  set(_libbsg_mesh_lod_impl "${BRLCAD_SOURCE_DIR}/src/libbsg/mesh_lod.cpp")
  if(EXISTS "${_libbsg_mesh_lod_impl}")
    file(READ "${_libbsg_mesh_lod_impl}" _libbsg_mesh_lod_impl_contents)
    foreach(_token
	"level > max_pop_threshold_level && !full_detail_setup_clbk"
	"if (!sp->set_level(level))"
	"bsg_mesh_lod_active_data_clear"
	"bsg_mesh_lod_active_pop_data_publish"
	"s->curr_level = -1"
	"return false"
	"else if (s->full_detail_clear_clbk)"
	"s->full_detail_setup_clbk = NULL"
	"s->detail_clbk_data = NULL"
	"sp->lod_tri_norms.size() >= sp->lod_tris.size() * 3"
	"l->normals = (const vect_t *)sp->lod_tri_norms.data()"
	"l->normals = NULL")
      string(FIND "${_libbsg_mesh_lod_impl_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/libbsg/mesh_lod.cpp must keep bounded POP normal publication token ${_token}")
      endif()
    endforeach()
  endif()

  set(_librt_view_legacy_impl
    "${BRLCAD_SOURCE_DIR}/src/librt/view_legacy_bsg.c")
  if(NOT EXISTS "${_librt_view_legacy_impl}")
    _brlobol_pivot_guard_fail(
      "src/librt/view_legacy_bsg.c is required for transitional BSG view adapters")
  else()
    file(READ "${_librt_view_legacy_impl}" _librt_view_legacy_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]bsg/lod\.h]]
	[[#[ \t]*include[ \t]*[<"]bsg/faceplate\.h]]
	[[#[ \t]*include[ \t]*[<"]bsg/view_state\.h]]
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	rt_view_info_from_bsg
	rt_view_orientation_quat_from_bsg
	rt_view_aet_from_bsg
	rt_view_aet_set_bsg
	rt_view_aet_state_set_bsg
	rt_view_perspective_from_bsg
	rt_view_perspective_set_bsg
	rt_view_model2view_from_bsg
	rt_view_model2view_set_bsg
	rt_view_view2model_from_bsg
	rt_view_view2model_set_bsg
	rt_view_pmodel2view_from_bsg
	rt_view_pmodel2view_set_bsg
	rt_view_pmat_from_bsg
	rt_view_pmat_set_bsg
	rt_view_rotation_from_bsg
	rt_view_rotation_set_bsg
	rt_view_center_from_bsg
	rt_view_center_vec_set_bsg
	rt_view_plane_from_bsg
	rt_view_lod_policy_from_bsg
	rt_view_lod_policy_apply_bsg
	rt_view_lod_policy_copy_bsg
	rt_view_lod_bounds_update_bsg
	rt_view_lod_bounds_callback_set_bsg
	rt_view_lod_bounds_callback_is_bsg
	rt_view_bounds_update_callback_from_bsg
	rt_view_bounds_update_callback_set_bsg
	rt_view_bounds_update_callback_call_bsg
	rt_view_scale_from_bsg
	rt_view_scale_set_bsg
	rt_view_scale_storage_from_bsg
	rt_view_scale_state_set_bsg
	rt_view_initial_scale_from_bsg
	rt_view_initial_scale_set_bsg
	rt_view_absolute_scale_from_bsg
	rt_view_absolute_scale_set_bsg
	rt_view_local2base_from_bsg
	rt_view_base2local_from_bsg
	rt_view_unit_conversion_set_bsg
	rt_view_width_from_bsg
	rt_view_height_from_bsg
	rt_view_radius_from_bsg
	rt_view_dimensions_set_bsg
	rt_view_screen_to_view_from_bsg
	rt_view_screen_point_from_bsg
	rt_view_interactive_rect_from_bsg
	rt_view_interactive_rect_set_bsg
	rt_view_adc_from_bsg
	rt_view_adc_set_bsg
	rt_view_grid_from_bsg
	rt_view_grid_set_bsg
	rt_view_model_axes_from_bsg
	rt_view_model_axes_set_bsg
	rt_view_view_axes_from_bsg
	rt_view_view_axes_set_bsg
	rt_view_center_dot_from_bsg
	rt_view_center_dot_set_bsg
	rt_view_scale_overlay_from_bsg
	rt_view_scale_overlay_set_bsg
	rt_view_params_from_bsg
	rt_view_params_set_bsg
	rt_view_refresh_request_bsg
	rt_view_refresh_dirty_from_bsg
	rt_view_refresh_consume_bsg
	rt_view_refresh_complete_bsg
	rt_view_refresh_enabled_from_bsg
	rt_view_refresh_enabled_set_bsg
	rt_view_refresh_suppressed_from_bsg
	rt_view_refresh_suppress_begin_bsg
	rt_view_refresh_suppress_end_bsg
	rt_view_refresh_drawn_count_from_bsg
	rt_view_refresh_drawn_count_set_bsg
	rt_view_size_from_bsg
	rt_view_size_set_bsg
	rt_view_inverse_size_from_bsg
	rt_view_eye_pos_from_bsg
	rt_view_eye_pos_set_bsg
	rt_view_keypoint_from_bsg
	rt_view_keypoint_set_bsg
	rt_view_rotate_about_from_bsg
	rt_view_rotate_about_set_bsg
	rt_view_coord_from_bsg
	rt_view_coord_set_bsg
	rt_view_snap_lines_from_bsg
	rt_view_snap_lines_set_bsg
	rt_view_snap_source_flags_from_bsg
	rt_view_snap_source_flags_set_bsg
	rt_view_snap_kind_mask_from_bsg
	rt_view_prepare_tcl_snap_bsg
	rt_view_center_linesnap_bsg
	rt_view_zclip_from_bsg
	rt_view_zclip_set_bsg
	rt_view_framebuffer_mode_from_bsg
	rt_view_framebuffer_mode_set_bsg
	rt_view_cleared_from_bsg
	rt_view_cleared_set_bsg
	rt_view_settings_shared_bsg
	rt_view_snap_tolerance_factor_from_bsg
	rt_view_snap_tolerance_factor_set_bsg
	rt_mesh_lod_load_view_scene_ref_bsg
	rt_mesh_lod_free_scene_ref_bsg
	bsg_view_adc_get
	bsg_view_adc_set
	bsg_view_grid_get
	bsg_view_grid_set
	bsg_view_model_axes_get
	bsg_view_model_axes_set
	bsg_view_view_axes_get
	bsg_view_view_axes_set
	bsg_view_center_dot_get
	bsg_view_center_dot_set
	bsg_view_scale_state_get
	bsg_view_scale_state_set
	bsg_view_params_get
	bsg_view_params_set
	bsg_view_lod_source_policy_get
	bsg_view_lod_source_policy_set
	bsg_view_bounds
	quat_mat2quat
	bsg_mesh_lod_load_view_scene_ref
	bsg_mesh_lod_free_scene_ref
	_rt_mesh_lod_bsg)
      string(REGEX MATCH "${_token}" _legacy_token_hit
	"${_librt_view_legacy_contents}")
      if(NOT _legacy_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/librt/view_legacy_bsg.c must own transitional BSG view adapter ${_token}")
      endif()
    endforeach()
  endif()

  set(_librt_cmake "${BRLCAD_SOURCE_DIR}/src/librt/CMakeLists.txt")
  if(EXISTS "${_librt_cmake}")
    file(READ "${_librt_cmake}" _librt_cmake_contents)
    string(FIND "${_librt_cmake_contents}" "view_legacy_bsg.c"
      _librt_view_legacy_cmake_idx)
    if(_librt_view_legacy_cmake_idx EQUAL -1)
      _brlobol_pivot_guard_fail(
	"src/librt/CMakeLists.txt must build src/librt/view_legacy_bsg.c")
    endif()
  endif()

  set(_librt_view_info_test
    "${BRLCAD_SOURCE_DIR}/src/librt/tests/view_info.c")
  if(NOT EXISTS "${_librt_view_info_test}")
    _brlobol_pivot_guard_fail(
      "src/librt/tests/view_info.c is required for neutral RT view-info coverage")
  else()
    file(READ "${_librt_view_info_test}" _librt_view_info_test_contents)
    foreach(_token
	rt_view_info_init
	rt_view_info_sanitize
	rt_view_lod_policy_init
	rt_view_lod_policy_sanitize
	test_lod_policy_defaults
	rt_view_lod_curve_scale
	rt_view_lod_bot_threshold
	rt_view_avg_sample_spacing
	rt_view_solid_point_spacing
	rt_mesh_lod_load_view
	rt_mesh_lod_current_level
	rt_mesh_lod_has_active_data
	rt_mesh_lod_data_get
	rt_mesh_lod_info_get
	rt_mesh_lod_detail_init
	rt_mesh_lod_detail_callbacks_set
	rt_mesh_lod_detail_callbacks_clear
	rt_mesh_lod_memshrink
	rt_mesh_lod_cache_status_init
	db_mesh_lod_status
	db_mesh_lod_refresh
	db_mesh_lod_invalidate
	db_mesh_lod_store_mesh
	db_mesh_lod_update
	"mesh lod full-detail callback did not publish valid normals"
	"mesh lod full-detail callback replacement did not preserve active POP data"
	"mesh lod full-detail callback replacement did not invalidate stale borrowed data"
	"mesh lod full-detail callback clear did not release callback ownership"
	"mesh lod full-detail callback setup failure did not clear stale active data"
	"mesh lod full-detail callback missing arrays did not clear stale active data"
	"mesh lod full-detail callback bad indices did not clear stale active data"
	"mesh lod full-detail callback bad original indices did not clear stale active data"
	"mesh lod full-detail callback producer failure did not clear partial active data"
	"mesh lod invalid BoT cache rejection"
	"mesh lod invalid generated mesh store should preserve existing cache status"
	"mesh lod generated mesh normals"
	"(const vect_t *)detail_normals"
	"mesh lod get after reopen")
      string(FIND "${_librt_view_info_test_contents}" "${_token}"
	_librt_view_info_test_token_idx)
      if(_librt_view_info_test_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/librt/tests/view_info.c must cover neutral RT view-info token ${_token}")
      endif()
    endforeach()
    set(_librt_view_info_test_forbidden
      [[#[ \t]*include[ \t]*[<"]bsg/]]
      [[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
      [[(^|[^A-Za-z0-9_])rt_view_info_from_bsg([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])bsg_view([^A-Za-z0-9_]|$)]]
      [[struct[ \t\r\n]+bsg_view]]
    )
    foreach(_pat IN LISTS _librt_view_info_test_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_librt_view_info_test_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "src/librt/tests/view_info.c reintroduced BSG adapter coverage into the neutral RT view-info test: ${_hit}")
      endif()
    endforeach()
  endif()

  set(_librt_view_legacy_bsg_test
    "${BRLCAD_SOURCE_DIR}/src/librt/tests/view_legacy_bsg.c")
  if(NOT EXISTS "${_librt_view_legacy_bsg_test}")
    _brlobol_pivot_guard_fail(
      "src/librt/tests/view_legacy_bsg.c is required for transitional RT view-info adapter coverage")
  else()
    file(READ "${_librt_view_legacy_bsg_test}" _librt_view_legacy_bsg_test_contents)
    foreach(_token
	"bsg/util.h"
	"bsg/view_state.h"
	"rt/edit_legacy_bsg.h"
	"rt/view_legacy_bsg.h"
	rt_edit_view_from_bsg
	rt_view_info_from_bsg
	rt_view_orientation_quat_from_bsg
	rt_view_aet_from_bsg
	rt_view_aet_set_bsg
	rt_view_aet_state_set_bsg
	rt_view_perspective_from_bsg
	rt_view_perspective_set_bsg
	rt_view_model2view_from_bsg
	rt_view_model2view_set_bsg
	rt_view_view2model_from_bsg
	rt_view_view2model_set_bsg
	rt_view_pmodel2view_from_bsg
	rt_view_pmodel2view_set_bsg
	rt_view_pmat_from_bsg
	rt_view_pmat_set_bsg
	rt_view_rotation_from_bsg
	rt_view_rotation_set_bsg
	rt_view_center_from_bsg
	rt_view_center_vec_set_bsg
	rt_view_plane_from_bsg
	rt_view_lod_policy_from_bsg
	rt_view_lod_policy_apply_bsg
	rt_view_lod_policy_copy_bsg
	rt_view_lod_bounds_update_bsg
	rt_view_lod_bounds_callback_set_bsg
	rt_view_lod_bounds_callback_is_bsg
	rt_view_bounds_update_callback_from_bsg
	rt_view_bounds_update_callback_set_bsg
	rt_view_bounds_update_callback_call_bsg
	rt_view_scale_from_bsg
	rt_view_scale_set_bsg
	rt_view_scale_storage_from_bsg
	rt_view_scale_state_set_bsg
	rt_view_initial_scale_from_bsg
	rt_view_initial_scale_set_bsg
	rt_view_absolute_scale_from_bsg
	rt_view_absolute_scale_set_bsg
	rt_view_local2base_from_bsg
	rt_view_base2local_from_bsg
	rt_view_unit_conversion_set_bsg
	rt_view_width_from_bsg
	rt_view_height_from_bsg
	rt_view_radius_from_bsg
	rt_view_dimensions_set_bsg
	rt_view_screen_to_view_from_bsg
	rt_view_screen_point_from_bsg
	rt_view_interactive_rect_from_bsg
	rt_view_interactive_rect_set_bsg
	rt_view_adc_from_bsg
	rt_view_adc_set_bsg
	rt_view_grid_from_bsg
	rt_view_grid_set_bsg
	rt_view_model_axes_from_bsg
	rt_view_model_axes_set_bsg
	rt_view_view_axes_from_bsg
	rt_view_view_axes_set_bsg
	rt_view_center_dot_from_bsg
	rt_view_center_dot_set_bsg
	rt_view_scale_overlay_from_bsg
	rt_view_scale_overlay_set_bsg
	rt_view_params_from_bsg
	rt_view_params_set_bsg
	rt_view_refresh_request_bsg
	rt_view_refresh_dirty_from_bsg
	rt_view_refresh_consume_bsg
	rt_view_refresh_complete_bsg
	rt_view_refresh_enabled_from_bsg
	rt_view_refresh_enabled_set_bsg
	rt_view_refresh_suppressed_from_bsg
	rt_view_refresh_suppress_begin_bsg
	rt_view_refresh_suppress_end_bsg
	rt_view_refresh_drawn_count_from_bsg
	rt_view_refresh_drawn_count_set_bsg
	rt_view_size_from_bsg
	rt_view_size_set_bsg
	rt_view_inverse_size_from_bsg
	rt_view_eye_pos_from_bsg
	rt_view_eye_pos_set_bsg
	rt_view_keypoint_from_bsg
	rt_view_keypoint_set_bsg
	rt_view_rotate_about_from_bsg
	rt_view_rotate_about_set_bsg
	rt_view_coord_from_bsg
	rt_view_coord_set_bsg
	rt_view_snap_lines_from_bsg
	rt_view_snap_lines_set_bsg
	rt_view_snap_source_flags_from_bsg
	rt_view_snap_source_flags_set_bsg
	rt_view_snap_kind_mask_from_bsg
	rt_view_prepare_tcl_snap_bsg
	rt_view_center_linesnap_bsg
	rt_view_zclip_from_bsg
	rt_view_zclip_set_bsg
	rt_view_framebuffer_mode_from_bsg
	rt_view_framebuffer_mode_set_bsg
	rt_view_cleared_from_bsg
	rt_view_cleared_set_bsg
	rt_view_settings_shared_bsg
	rt_view_snap_tolerance_factor_from_bsg
	rt_view_snap_tolerance_factor_set_bsg
	rt_mesh_lod_load_view_scene_ref_bsg
	rt_mesh_lod_free_scene_ref_bsg
	test_bsg_lod_policy_adapter
	bsg_view_lod_source_policy_set
	test_bsg_null_adapter
	test_bsg_adapter
	test_bsg_adapter_sanitizes
	test_bsg_orientation_adapter
	test_bsg_aet_adapter
	test_bsg_perspective_adapter
	test_bsg_camera_adapter
	test_bsg_mesh_lod_adapter_boundary
	test_bsg_faceplate_state_adapter
	bsg_view_adc_get
	bsg_view_grid_get
	bsg_view_model_axes_get
	bsg_view_view_axes_get
	bsg_view_center_dot_get
	bsg_view_scale_state_get
	bsg_view_params_get
	test_bsg_edit_view_adapter)
      string(FIND "${_librt_view_legacy_bsg_test_contents}" "${_token}"
	_librt_view_legacy_bsg_test_token_idx)
      if(_librt_view_legacy_bsg_test_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/librt/tests/view_legacy_bsg.c must cover RT view-info adapter token ${_token}")
      endif()
    endforeach()
  endif()

  foreach(_rel
      src/libged/check/check.c
      src/libged/get_eyemodel/get_eyemodel.c
      src/libged/rtwizard/rtwizard.c
      src/libged/ged_util.cpp)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for GED RT view-info adapter checks")
      continue()
    endif()

    file(READ "${_file}" _ged_view_snapshot_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_info_from_bsg]]
	[[rt_view_orientation_quat_from_bsg]])
      string(REGEX MATCH "${_token}" _ged_view_snapshot_token_hit
	"${_ged_view_snapshot_contents}")
      if(NOT _ged_view_snapshot_token_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must route view size/orientation extraction through rt/view_legacy_bsg.h")
      endif()
    endforeach()

    foreach(_pat
	[[quat_mat2quat[ \t\r\n]*\([^;]*gv_rotation]]
	[=[gv_size[ \t\r\n]*[,;)]]=])
      string(REGEX MATCH "${_pat}" _ged_view_snapshot_direct_hit
	"${_ged_view_snapshot_contents}")
      if(_ged_view_snapshot_direct_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced direct BSG view size/orientation reads: ${_ged_view_snapshot_direct_hit}")
      endif()
    endforeach()
  endforeach()

  foreach(_rel
      src/libged/ged_util.cpp
      src/libged/rt/rt.c
      src/libged/rtwizard/rtwizard.c
      src/libged/view/saveview.c
      src/libged/dm/ert.cpp)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for RT command perspective adapter checks")
      continue()
    endif()

    file(READ "${_file}" _ged_perspective_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_perspective_from_bsg]])
      string(REGEX MATCH "${_token}" _ged_perspective_token_hit
	"${_ged_perspective_contents}")
      if(NOT _ged_perspective_token_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must route RT command perspective reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()

    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_perspective([^A-Za-z0-9_]|$)]]
      _ged_direct_perspective_hit "${_ged_perspective_contents}")
    if(_ged_direct_perspective_hit)
      _brlobol_pivot_guard_fail(
	"${_rel} reintroduced direct BSG view perspective reads: ${_ged_direct_perspective_hit}")
    endif()

    if(_rel STREQUAL "src/libged/dm/ert.cpp")
      foreach(_token
	  [[rt_view_framebuffer_mode_from_bsg]]
	  [[rt_view_framebuffer_mode_set_bsg]])
	string(REGEX MATCH "${_token}" _ged_ert_fb_token_hit
	  "${_ged_perspective_contents}")
	if(NOT _ged_ert_fb_token_hit)
	  _brlobol_pivot_guard_fail(
	    "src/libged/dm/ert.cpp must route framebuffer mode view-state access through rt/view_legacy_bsg.h")
	endif()
      endforeach()
      string(REGEX MATCH [[(^|[^A-Za-z0-9_])bsg_view_(set_)?framebuffer_mode([^A-Za-z0-9_]|$)]]
	_ged_ert_fb_direct_hit "${_ged_perspective_contents}")
      if(_ged_ert_fb_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/dm/ert.cpp reintroduced direct BSG framebuffer mode view-state access: ${_ged_ert_fb_direct_hit}")
      endif()
    endif()
  endforeach()

  set(_ged_eye_model_helper "${BRLCAD_SOURCE_DIR}/src/libged/ged_util.cpp")
  if(EXISTS "${_ged_eye_model_helper}")
    file(READ "${_ged_eye_model_helper}" _ged_eye_model_helper_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_view2model_from_bsg]]
	[[rt_view_center_from_bsg]]
	[[rt_view_rotation_from_bsg]]
	[[rt_view_zclip_from_bsg]])
      string(REGEX MATCH "${_token}" _ged_eye_model_helper_token_hit
	"${_ged_eye_model_helper_contents}")
      if(NOT _ged_eye_model_helper_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/ged_util.cpp must route eye-model view snapshot reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_center([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_rotation([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_zclip([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _ged_eye_model_helper_direct_hit
	"${_ged_eye_model_helper_contents}")
      if(_ged_eye_model_helper_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/ged_util.cpp reintroduced direct BSG eye-model view snapshot reads: ${_ged_eye_model_helper_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_nirt_impl "${BRLCAD_SOURCE_DIR}/src/libged/nirt/nirt.cpp")
  if(EXISTS "${_libged_nirt_impl}")
    file(READ "${_libged_nirt_impl}" _libged_nirt_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_center_from_bsg]]
	[[rt_view_rotation_from_bsg]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_view2model_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_nirt_token_hit
	"${_libged_nirt_contents}")
      if(NOT _libged_nirt_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/nirt/nirt.cpp must route NIRT view snapshot reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_center([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_rotation([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libged_nirt_direct_hit
	"${_libged_nirt_contents}")
      if(_libged_nirt_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/nirt/nirt.cpp reintroduced direct BSG NIRT view snapshot reads: ${_libged_nirt_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_bot_pick_impl "${BRLCAD_SOURCE_DIR}/src/libged/bot/pick.cpp")
  if(EXISTS "${_libged_bot_pick_impl}")
    file(READ "${_libged_bot_pick_impl}" _libged_bot_pick_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_center_from_bsg]]
	[[rt_view_rotation_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_bot_pick_token_hit
	"${_libged_bot_pick_contents}")
      if(NOT _libged_bot_pick_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/bot/pick.cpp must route viewport pick ray reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_center([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_rotation([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libged_bot_pick_direct_hit
	"${_libged_bot_pick_contents}")
      if(_libged_bot_pick_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/bot/pick.cpp reintroduced direct BSG viewport pick ray reads: ${_libged_bot_pick_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_brep_pick_impl "${BRLCAD_SOURCE_DIR}/src/libged/brep/pick.cpp")
  if(EXISTS "${_libged_brep_pick_impl}")
    file(READ "${_libged_brep_pick_impl}" _libged_brep_pick_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_center_from_bsg]]
	[[rt_view_rotation_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_brep_pick_token_hit
	"${_libged_brep_pick_contents}")
      if(NOT _libged_brep_pick_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/brep/pick.cpp must route viewport pick ray reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_center([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_rotation([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libged_brep_pick_direct_hit
	"${_libged_brep_pick_contents}")
      if(_libged_brep_pick_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/brep/pick.cpp reintroduced direct BSG viewport pick ray reads: ${_libged_brep_pick_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_ged_view_coord_helper "${BRLCAD_SOURCE_DIR}/src/libged/ged_util.cpp")
  if(EXISTS "${_ged_view_coord_helper}")
    file(READ "${_ged_view_coord_helper}" _ged_view_coord_helper_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_coord_from_bsg]])
      string(REGEX MATCH "${_token}" _ged_view_coord_helper_token_hit
	"${_ged_view_coord_helper_contents}")
      if(NOT _ged_view_coord_helper_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/ged_util.cpp must route default coord reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[\*[ \t]*coord[ \t]*=[^\n;]*gv_coord]]
      _ged_view_coord_helper_direct_hit "${_ged_view_coord_helper_contents}")
    if(_ged_view_coord_helper_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/ged_util.cpp reintroduced direct BSG default coord reads: ${_ged_view_coord_helper_direct_hit}")
    endif()
  endif()

  set(_libged_arot_coord_impl "${BRLCAD_SOURCE_DIR}/src/libged/arot/arot.c")
  if(EXISTS "${_libged_arot_coord_impl}")
    file(READ "${_libged_arot_coord_impl}" _libged_arot_coord_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_coord_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_arot_coord_adapter_hit
	"${_libged_arot_coord_contents}")
      if(NOT _libged_arot_coord_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/arot/arot.c must route default coord reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[_ged_do_rot[^\n;]*gv_coord]]
      _libged_arot_coord_direct_hit "${_libged_arot_coord_contents}")
    if(_libged_arot_coord_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/arot/arot.c reintroduced direct BSG default coord reads: ${_libged_arot_coord_direct_hit}")
    endif()
  endif()

  set(_ged_savekey_file "${BRLCAD_SOURCE_DIR}/src/libged/savekey/savekey.c")
  if(NOT EXISTS "${_ged_savekey_file}")
    _brlobol_pivot_guard_fail(
      "src/libged/savekey/savekey.c is required for RT -M savekey adapter checks")
  else()
    file(READ "${_ged_savekey_file}" _ged_savekey_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_info_from_bsg]]
	[[rt_view_rotation_from_bsg]]
	[[rt_view_view2model_from_bsg]])
      string(REGEX MATCH "${_token}" _ged_savekey_token_hit
	"${_ged_savekey_contents}")
      if(NOT _ged_savekey_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/savekey/savekey.c must route RT -M view reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()

    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_size([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_rotation([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _ged_savekey_direct_hit
	"${_ged_savekey_contents}")
      if(_ged_savekey_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/savekey/savekey.c reintroduced direct BSG RT -M view reads: ${_ged_savekey_direct_hit}")
      endif()
    endforeach()
  endif()

  foreach(_rel
      src/libged/plot/plot.c
      src/libged/ps/ps.c
      src/libged/png/png.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for export camera adapter checks")
      continue()
    endif()

    file(READ "${_file}" _ged_export_camera_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_model2view_from_bsg]])
      string(REGEX MATCH "${_token}" _ged_export_camera_token_hit
	"${_ged_export_camera_contents}")
      if(NOT _ged_export_camera_token_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must route export camera reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()

    if("${_rel}" STREQUAL "src/libged/plot/plot.c")
      foreach(_token
	  [[rt_view_center_from_bsg]]
	  [[rt_view_scale_from_bsg]])
	string(REGEX MATCH "${_token}" _ged_plot_camera_token_hit
	  "${_ged_export_camera_contents}")
	if(NOT _ged_plot_camera_token_hit)
	  _brlobol_pivot_guard_fail(
	    "${_rel} must route plot camera center/scale reads through rt/view_legacy_bsg.h")
	endif()
      endforeach()
    else()
      foreach(_token
	  [[rt_view_perspective_from_bsg]]
	  [[rt_view_eye_pos_from_bsg]])
	string(REGEX MATCH "${_token}" _ged_bitmap_camera_token_hit
	  "${_ged_export_camera_contents}")
	if(NOT _ged_bitmap_camera_token_hit)
	  _brlobol_pivot_guard_fail(
	    "${_rel} must route bitmap/vector perspective camera reads through rt/view_legacy_bsg.h")
	endif()
      endforeach()
    endif()

    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_center([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_perspective([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_eye_pos([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _ged_export_camera_direct_hit
	"${_ged_export_camera_contents}")
      if(_ged_export_camera_direct_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced direct BSG export camera reads: ${_ged_export_camera_direct_hit}")
      endif()
    endforeach()
  endforeach()

  foreach(_rel
      src/libged/m2v_point/m2v_point.c
      src/libged/model2view/model2view.c
      src/libged/model2view_lu/model2view_lu.c
      src/libged/model2grid_lu/model2grid_lu.c
      src/libged/view2grid_lu/view2grid_lu.c
      src/libged/grid2model_lu/grid2model_lu.c
      src/libged/grid2view_lu/grid2view_lu.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for GED model2view conversion adapter checks")
      continue()
    endif()
    file(READ "${_file}" _ged_model2view_conv_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_model2view_from_bsg]])
      string(REGEX MATCH "${_token}" _ged_model2view_conv_token_hit
	"${_ged_model2view_conv_contents}")
      if(NOT _ged_model2view_conv_token_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must route model2view conversion reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
      _ged_model2view_conv_direct_hit "${_ged_model2view_conv_contents}")
    if(_ged_model2view_conv_direct_hit)
      _brlobol_pivot_guard_fail(
	"${_rel} reintroduced direct BSG model2view conversion reads: ${_ged_model2view_conv_direct_hit}")
    endif()
  endforeach()

  foreach(_rel
      src/libged/v2m_point/v2m_point.c
      src/libged/view2model/view2model.c
      src/libged/view2model_lu/view2model_lu.c
      src/libged/grid2model_lu/grid2model_lu.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for GED view2model conversion adapter checks")
      continue()
    endif()
    file(READ "${_file}" _ged_view2model_conv_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_view2model_from_bsg]])
      string(REGEX MATCH "${_token}" _ged_view2model_conv_token_hit
	"${_ged_view2model_conv_contents}")
      if(NOT _ged_view2model_conv_token_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must route view2model conversion reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]]
      _ged_view2model_conv_direct_hit "${_ged_view2model_conv_contents}")
    if(_ged_view2model_conv_direct_hit)
      _brlobol_pivot_guard_fail(
	"${_rel} reintroduced direct BSG view2model conversion reads: ${_ged_view2model_conv_direct_hit}")
    endif()
  endforeach()

  foreach(_rel
      src/libged/view2model_vec/view2model_vec.c
      src/libged/view/viewdir.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for GED view rotation adapter checks")
      continue()
    endif()

    file(READ "${_file}" _ged_view_rotation_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_rotation_from_bsg]])
      string(REGEX MATCH "${_token}" _ged_view_rotation_token_hit
	"${_ged_view_rotation_contents}")
      if(NOT _ged_view_rotation_token_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must route view rotation reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_rotation([^A-Za-z0-9_]|$)]]
      _ged_view_rotation_direct_hit "${_ged_view_rotation_contents}")
    if(_ged_view_rotation_direct_hit)
      _brlobol_pivot_guard_fail(
	"${_rel} reintroduced direct BSG view rotation reads: ${_ged_view_rotation_direct_hit}")
    endif()
  endforeach()

  foreach(_rel
      src/libged/model2view_lu/model2view_lu.c
      src/libged/view2model_lu/view2model_lu.c
      src/libged/model2grid_lu/model2grid_lu.c
      src/libged/view2grid_lu/view2grid_lu.c
      src/libged/grid2model_lu/grid2model_lu.c
      src/libged/grid2view_lu/grid2view_lu.c
      src/libged/overlay/overlay.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for GED view-scale conversion adapter checks")
      continue()
    endif()
    file(READ "${_file}" _ged_view_scale_conv_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_from_bsg]])
      string(REGEX MATCH "${_token}" _ged_view_scale_conv_token_hit
	"${_ged_view_scale_conv_contents}")
      if(NOT _ged_view_scale_conv_token_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must route view-scale conversion reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]]
      _ged_view_scale_conv_direct_hit "${_ged_view_scale_conv_contents}")
    if(_ged_view_scale_conv_direct_hit)
      _brlobol_pivot_guard_fail(
	"${_rel} reintroduced direct BSG view-scale conversion reads: ${_ged_view_scale_conv_direct_hit}")
    endif()
  endforeach()

  foreach(_rel
      src/libqtcad/QgQuadView.cpp
      src/mged/chgview.c
      src/mged/dm-generic.c
      src/mged/attach.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      continue()
    endif()
    file(READ "${_file}" _lod_policy_boundary_contents)
    string(REGEX MATCH [[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
      _lod_policy_boundary_header_hit "${_lod_policy_boundary_contents}")
    if(NOT _lod_policy_boundary_header_hit)
      _brlobol_pivot_guard_fail(
	"${_rel} must use rt/view_legacy_bsg.h for transitional BSG LoD policy adaptation")
    endif()
    if("${_rel}" STREQUAL "src/libqtcad/QgQuadView.cpp")
      foreach(_token
	  rt_view_lod_policy_copy_bsg
	  rt_view_lod_bounds_update_bsg)
	string(FIND "${_lod_policy_boundary_contents}" "${_token}"
	  _lod_policy_boundary_qtcad_idx)
	if(_lod_policy_boundary_qtcad_idx EQUAL -1)
	  _brlobol_pivot_guard_fail(
	    "${_rel} must route quad-view LoD policy/bounds behavior through ${_token}")
	endif()
      endforeach()
    elseif("${_rel}" STREQUAL "src/mged/attach.c")
      foreach(_token
	  rt_view_lod_policy_from_bsg
	  rt_view_lod_policy_apply_bsg)
	string(FIND "${_lod_policy_boundary_contents}" "${_token}"
	  _lod_policy_boundary_attach_idx)
	if(_lod_policy_boundary_attach_idx EQUAL -1)
	  _brlobol_pivot_guard_fail(
	    "${_rel} must route copied display LoD policy through ${_token}")
	endif()
      endforeach()
    else()
      string(FIND "${_lod_policy_boundary_contents}" "rt_view_lod_policy_from_bsg"
	_lod_policy_boundary_from_idx)
      if(_lod_policy_boundary_from_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "${_rel} must inspect zoom-refresh LoD policy through rt_view_lod_policy_from_bsg")
      endif()
    endif()
    foreach(_pat
	[[struct[ \t\r\n]+bsg_lod_source_policy_settings]]
	[[bsg_view_lod_source_policy_get]]
	[[bsg_view_lod_source_policy_set]]
	[[bsg_view_bounds]])
      string(REGEX MATCH "${_pat}" _lod_policy_boundary_bsg_hit
	"${_lod_policy_boundary_contents}")
      if(_lod_policy_boundary_bsg_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must not call BSG LoD policy/bounds APIs directly: ${_lod_policy_boundary_bsg_hit}")
      endif()
    endforeach()
  endforeach()

  set(_librt_lod_test
    "${BRLCAD_SOURCE_DIR}/src/librt/tests/lod.c")
  if(NOT EXISTS "${_librt_lod_test}")
    _brlobol_pivot_guard_fail(
      "src/librt/tests/lod.c is required for RT mesh-LoD API coverage")
  else()
    file(READ "${_librt_lod_test}" _librt_lod_test_contents)
    foreach(_token
	db_mesh_lod_update
	db_mesh_lod_status
	db_mesh_lod_refresh
	db_mesh_lod_invalidate
	db_mesh_lod_get
	rt_mesh_lod_load_level
	rt_mesh_lod_current_level
	rt_mesh_lod_has_active_data
	rt_mesh_lod_data_get
	rt_mesh_lod_info_get
	rt_mesh_lod_memshrink
	rt_mesh_lod_destroy)
      string(FIND "${_librt_lod_test_contents}" "${_token}"
	_librt_lod_test_token_idx)
      if(_librt_lod_test_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/librt/tests/lod.c must exercise neutral RT mesh-LoD helper ${_token}")
      endif()
    endforeach()
    set(_librt_lod_test_forbidden
      [[#[ \t]*include[ \t]*[<"]bsg/]]
      [[(^|[^A-Za-z0-9_])bsg_[A-Za-z0-9_]*]]
      [[struct[ \t\r\n]+bsg_[A-Za-z0-9_]*]]
      [[(^|[^A-Za-z0-9_])BSG_[A-Za-z0-9_]*]]
    )
    foreach(_pat IN LISTS _librt_lod_test_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_librt_lod_test_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "src/librt/tests/lod.c reintroduced direct BSG mesh-LoD coverage: ${_hit}")
      endif()
    endforeach()
  endif()

  foreach(_rel
      src/libged/solids_on_ray/solids_on_ray.c
      src/libged/plot/plot.c
      src/libged/nirt/nirt.cpp
      src/libged/adc/adc.c
      src/libged/dm/dm.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for RT view-coordinate constant neutralization checks")
      continue()
    endif()
    file(READ "${_file}" _libged_rt_view_constant_contents)
    foreach(_token
	"rt/view.h"
	"RT_VIEW_")
      string(FIND "${_libged_rt_view_constant_contents}" "${_token}"
	_libged_rt_view_constant_token_idx)
      if(_libged_rt_view_constant_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "${_rel} must use RT-owned view-coordinate constants via ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])BSG_VIEW_(MIN|MAX|RANGE)([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])BSG_INV_(VIEW|4096)([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])INV_BV([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libged_rt_view_constant_hit
	"${_libged_rt_view_constant_contents}")
      if(_libged_rt_view_constant_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must not use BSG-owned spellings for normalized RT view coordinates: ${_libged_rt_view_constant_hit}")
      endif()
    endforeach()
  endforeach()

  set(_libged_dm_view_dims_impl "${BRLCAD_SOURCE_DIR}/src/libged/dm/dm.c")
  if(EXISTS "${_libged_dm_view_dims_impl}")
    file(READ "${_libged_dm_view_dims_impl}" _libged_dm_view_dims_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_width_from_bsg]]
	[[rt_view_height_from_bsg]]
	[[rt_view_dimensions_set_bsg]]
	[[rt_view_scale_storage_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_dm_view_dims_token_hit
	"${_libged_dm_view_dims_contents}")
      if(NOT _libged_dm_view_dims_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/dm/dm.c must route DM fallback view reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[if[ \t\r\n]*\([ \t\r\n]*![ \t\r\n]*target_view->[ \t\r\n]*gv_width[ \t\r\n]*\)]]
	[[if[ \t\r\n]*\([ \t\r\n]*![ \t\r\n]*target_view->[ \t\r\n]*gv_height[ \t\r\n]*\)]]
	[[target_view->[ \t\r\n]*gv_(width|height)[ \t\r\n]*=]]
	[[dm_set_vp[ \t\r\n]*\([^;]*gv_scale]])
      string(REGEX MATCH "${_pat}" _libged_dm_view_dims_direct_hit
	"${_libged_dm_view_dims_contents}")
      if(_libged_dm_view_dims_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/dm/dm.c reintroduced direct BSG DM fallback view reads: ${_libged_dm_view_dims_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_solids_on_ray_impl
    "${BRLCAD_SOURCE_DIR}/src/libged/solids_on_ray/solids_on_ray.c")
  if(EXISTS "${_libged_solids_on_ray_impl}")
    file(READ "${_libged_solids_on_ray_impl}" _libged_solids_on_ray_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_center_from_bsg]]
	[[rt_view_rotation_from_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_scale_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_solids_on_ray_token_hit
	"${_libged_solids_on_ray_contents}")
      if(NOT _libged_solids_on_ray_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/solids_on_ray/solids_on_ray.c must route ray view reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_center([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_rotation([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libged_solids_on_ray_direct_hit
	"${_libged_solids_on_ray_contents}")
      if(_libged_solids_on_ray_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/solids_on_ray/solids_on_ray.c reintroduced direct BSG ray view reads: ${_libged_solids_on_ray_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_adc_impl "${BRLCAD_SOURCE_DIR}/src/libged/adc/adc.c")
  if(EXISTS "${_libged_adc_impl}")
    file(READ "${_libged_adc_impl}" _libged_adc_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_view2model_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_adc_token_hit
	"${_libged_adc_contents}")
      if(NOT _libged_adc_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/adc/adc.c must route ADC view scale/matrix reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libged_adc_direct_hit
	"${_libged_adc_contents}")
      if(_libged_adc_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/adc/adc.c reintroduced direct BSG ADC view reads: ${_libged_adc_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_grid_impl "${BRLCAD_SOURCE_DIR}/src/libged/grid/grid.c")
  if(EXISTS "${_libged_grid_impl}")
    file(READ "${_libged_grid_impl}" _libged_grid_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_center_from_bsg]]
	[[rt_view_center_vec_set_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_view2model_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_grid_token_hit
	"${_libged_grid_contents}")
      if(NOT _libged_grid_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/grid/grid.c must route vsnap view reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]]
	[[MAT_DELTAS_GET_NEG[ \t\r\n]*\([^;]*gv_center]]
	[[MAT_DELTAS_VEC_NEG[ \t\r\n]*\([^;]*gv_center]])
      string(REGEX MATCH "${_pat}" _libged_grid_direct_hit
	"${_libged_grid_contents}")
      if(_libged_grid_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/grid/grid.c reintroduced direct BSG vsnap view reads: ${_libged_grid_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_eye_cmd "${BRLCAD_SOURCE_DIR}/src/libged/view/eye.c")
  if(EXISTS "${_libged_eye_cmd}")
    file(READ "${_libged_eye_cmd}" _libged_eye_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_view2model_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_eye_cmd_token_hit
	"${_libged_eye_cmd_contents}")
      if(NOT _libged_eye_cmd_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/eye.c must route eye view-to-model reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]]
      _libged_eye_cmd_direct_hit "${_libged_eye_cmd_contents}")
    if(_libged_eye_cmd_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/eye.c reintroduced direct BSG eye view-to-model reads: ${_libged_eye_cmd_direct_hit}")
    endif()
    string(REGEX MATCH [[MAT_DELTAS_VEC_NEG[ \t\r\n]*\([^;]*gv_center]]
      _libged_eye_cmd_center_direct_hit "${_libged_eye_cmd_contents}")
    if(_libged_eye_cmd_center_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/eye.c reintroduced direct BSG eye center writes: ${_libged_eye_cmd_center_direct_hit}")
    endif()
  endif()

  set(_libged_rect_impl "${BRLCAD_SOURCE_DIR}/src/libged/rect/rect.c")
  if(EXISTS "${_libged_rect_impl}")
    file(READ "${_libged_rect_impl}" _libged_rect_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_center_from_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_view2model_from_bsg]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_center_vec_set_bsg]]
	[[rt_view_scale_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_rect_token_hit
	"${_libged_rect_contents}")
      if(NOT _libged_rect_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/rect/rect.c must route rectangle-zoom view reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]]
	[[MAT_DELTAS_GET_NEG[ \t\r\n]*\([^;]*gv_center]]
	[[MAT_DELTAS_VEC_NEG[ \t\r\n]*\([^;]*gv_center]]
	[[(^|[^A-Za-z0-9_])gv_scale[ \t]*[*/+-]?=]])
      string(REGEX MATCH "${_pat}" _libged_rect_direct_hit
	"${_libged_rect_contents}")
      if(_libged_rect_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/rect/rect.c reintroduced direct BSG rectangle-zoom view reads: ${_libged_rect_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_lookat_cmd "${BRLCAD_SOURCE_DIR}/src/libged/view/lookat.c")
  if(EXISTS "${_libged_lookat_cmd}")
    file(READ "${_libged_lookat_cmd}" _libged_lookat_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_view2model_from_bsg]]
	[[rt_view_aet_from_bsg]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_aet_set_bsg]]
	[[rt_view_center_vec_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_lookat_cmd_token_hit
	"${_libged_lookat_cmd_contents}")
      if(NOT _libged_lookat_cmd_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/lookat.c must route lookat view reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]]
	"gv_aet[ \t\r\n]*\\[[ \t\r\n]*Z[ \t\r\n]*\\]"
	[[MAT_DELTAS_VEC_NEG[ \t\r\n]*\([^;]*gv_center]]
	[[VSET[ \t\r\n]*\([^;]*gv_aet]]
	[[bsg_mat_aet[ \t\r\n]*\(]])
      string(REGEX MATCH "${_pat}" _libged_lookat_cmd_direct_hit
	"${_libged_lookat_cmd_contents}")
      if(_libged_lookat_cmd_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/lookat.c reintroduced direct BSG lookat view reads: ${_libged_lookat_cmd_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_qvrot_cmd "${BRLCAD_SOURCE_DIR}/src/libged/view/qvrot.c")
  if(EXISTS "${_libged_qvrot_cmd}")
    file(READ "${_libged_qvrot_cmd}" _libged_qvrot_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_rotation_from_bsg]]
	[[rt_view_rotation_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_qvrot_cmd_token_hit
	"${_libged_qvrot_cmd_contents}")
      if(NOT _libged_qvrot_cmd_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/qvrot.c must route incremental rotation reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[bn_mat_mul2[^\n;]*gv_rotation]]
	[[bn_mat_angles[^\n;]*gv_rotation]]
	[[MAT_COPY[^\n;]*gv_rotation]])
      string(REGEX MATCH "${_pat}" _libged_qvrot_cmd_direct_hit
	"${_libged_qvrot_cmd_contents}")
      if(_libged_qvrot_cmd_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/qvrot.c reintroduced direct BSG incremental rotation access: ${_libged_qvrot_cmd_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_rtcheck2_impl
    "${BRLCAD_SOURCE_DIR}/src/libged/rtcheck/rtcheck2.cpp")
  if(EXISTS "${_libged_rtcheck2_impl}")
    file(READ "${_libged_rtcheck2_impl}" _libged_rtcheck2_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_rtcheck2_token_hit
	"${_libged_rtcheck2_contents}")
      if(NOT _libged_rtcheck2_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/rtcheck/rtcheck2.cpp must route overlap glyph scale through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]]
      _libged_rtcheck2_direct_hit "${_libged_rtcheck2_contents}")
    if(_libged_rtcheck2_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/rtcheck/rtcheck2.cpp reintroduced direct BSG overlap glyph scale reads: ${_libged_rtcheck2_direct_hit}")
    endif()
  endif()

  set(_libged_zoom_impl "${BRLCAD_SOURCE_DIR}/src/libged/zoom/zoom.c")
  if(EXISTS "${_libged_zoom_impl}")
    file(READ "${_libged_zoom_impl}" _libged_zoom_contents)
    foreach(_token
	"rt/view.h"
	"RT_VIEW_MIN_SCALE")
      string(FIND "${_libged_zoom_contents}" "${_token}" _libged_zoom_token_idx)
      if(_libged_zoom_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/libged/zoom/zoom.c must use RT-owned view scale constant ${_token}")
      endif()
    endforeach()
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_zoom_adapter_hit
	"${_libged_zoom_contents}")
      if(NOT _libged_zoom_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/zoom/zoom.c must route zoom view-scale reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])BSG_MINVIEWSCALE([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])BSG_MINVIEWSIZE([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libged_zoom_bsg_scale_hit
	"${_libged_zoom_contents}")
      if(_libged_zoom_bsg_scale_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/zoom/zoom.c must not use BSG-owned spellings for view scale limits: ${_libged_zoom_bsg_scale_hit}")
      endif()
    endforeach()
    foreach(_pat
	[[gv_scale[ \t]*=]]
	[[gv_scale[ \t]*/=]]
	[[gv_size[ \t]*=[^\n;]*gv_scale]]
	[[gv_isize[ \t]*=[^\n;]*gv_size]])
      string(REGEX MATCH "${_pat}" _libged_zoom_bsg_readback_hit
	"${_libged_zoom_contents}")
      if(_libged_zoom_bsg_readback_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/zoom/zoom.c reintroduced direct BSG zoom readbacks: ${_libged_zoom_bsg_readback_hit}")
      endif()
    endforeach()
  endif()

  foreach(_rel
      src/libged/view/size.c
      src/libged/scale/scale.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for RT view size limit neutralization checks")
      continue()
    endif()
    file(READ "${_file}" _libged_view_size_contents)
    foreach(_token
	"rt/view.h"
	"RT_VIEW_MIN_SIZE")
      string(FIND "${_libged_view_size_contents}" "${_token}"
	_libged_view_size_token_idx)
      if(_libged_view_size_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "${_rel} must use RT-owned view size limit ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])BSG_MINVIEWSIZE([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])BSG_MINVIEWSCALE([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libged_view_size_bsg_hit
	"${_libged_view_size_contents}")
      if(_libged_view_size_bsg_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must not use BSG-owned spellings for view size/scale limits: ${_libged_view_size_bsg_hit}")
      endif()
    endforeach()
  endforeach()

  set(_libged_scale_impl "${BRLCAD_SOURCE_DIR}/src/libged/scale/scale.c")
  if(EXISTS "${_libged_scale_impl}")
    file(READ "${_libged_scale_impl}" _libged_scale_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_scale_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_scale_adapter_hit
	"${_libged_scale_contents}")
      if(NOT _libged_scale_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/scale/scale.c must route scale-command view-scale reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[gv_scale[ \t]*=]]
	[[gv_scale[ \t]*\*=]]
	[[gv_size[ \t]*=[^\n;]*gv_scale]]
	[[gv_isize[ \t]*=[^\n;]*gv_size]])
      string(REGEX MATCH "${_pat}" _libged_scale_readback_hit
	"${_libged_scale_contents}")
      if(_libged_scale_readback_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/scale/scale.c reintroduced direct BSG scale-command readbacks: ${_libged_scale_readback_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_view_size_cmd "${BRLCAD_SOURCE_DIR}/src/libged/view/size.c")
  if(EXISTS "${_libged_view_size_cmd}")
    file(READ "${_libged_view_size_cmd}" _libged_view_size_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_info_from_bsg]]
	[[rt_view_size_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_view_size_cmd_token_hit
	"${_libged_view_size_cmd_contents}")
      if(NOT _libged_view_size_cmd_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/size.c must route query output through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[bu_vls_printf[^\n;]*gv_size]]
      _libged_view_size_cmd_direct_hit "${_libged_view_size_cmd_contents}")
    if(_libged_view_size_cmd_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/size.c reintroduced direct BSG query output: ${_libged_view_size_cmd_direct_hit}")
    endif()
    foreach(_pat
	[[gv_size[ \t]*=]]
	[[gv_isize[ \t]*=[^\n;]*gv_size]]
	[[gv_scale[ \t]*=[^\n;]*gv_size]])
      string(REGEX MATCH "${_pat}" _libged_view_size_cmd_readback_hit
	"${_libged_view_size_cmd_contents}")
      if(_libged_view_size_cmd_readback_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/size.c reintroduced direct BSG size-command readbacks: ${_libged_view_size_cmd_readback_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_view_center_cmd "${BRLCAD_SOURCE_DIR}/src/libged/view/center.cpp")
  if(EXISTS "${_libged_view_center_cmd}")
    file(READ "${_libged_view_center_cmd}" _libged_view_center_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_center_from_bsg]]
	[[rt_view_center_vec_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_view_center_cmd_token_hit
	"${_libged_view_center_cmd_contents}")
      if(NOT _libged_view_center_cmd_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/center.cpp must route query output through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[MAT_DELTAS_GET_NEG[ \t\r\n]*\([^;]*gv_center]]
      _libged_view_center_cmd_direct_hit
      "${_libged_view_center_cmd_contents}")
    if(_libged_view_center_cmd_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/center.cpp reintroduced direct BSG query output: ${_libged_view_center_cmd_direct_hit}")
    endif()
    string(REGEX MATCH [[MAT_DELTAS_VEC_NEG[ \t\r\n]*\([^;]*gv_center]]
      _libged_view_center_cmd_set_direct_hit
      "${_libged_view_center_cmd_contents}")
    if(_libged_view_center_cmd_set_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/center.cpp reintroduced direct BSG center writes: ${_libged_view_center_cmd_set_direct_hit}")
    endif()
  endif()

  set(_libged_view_aet_cmd "${BRLCAD_SOURCE_DIR}/src/libged/view/aet.c")
  if(EXISTS "${_libged_view_aet_cmd}")
    file(READ "${_libged_view_aet_cmd}" _libged_view_aet_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_aet_from_bsg]]
	[[rt_view_aet_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_view_aet_cmd_token_hit
	"${_libged_view_aet_cmd_contents}")
      if(NOT _libged_view_aet_cmd_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/aet.c must route query output through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[bn_encode_vect[^\n;]*gv_aet]]
	[[VADD2[^\n;]*gv_aet[^\n;]*gv_aet]]
	[[VMOVE[^\n;]*gv_aet]]
	[[VSET[^\n;]*gv_aet]]
	[[bsg_mat_aet[ \t\r\n]*\(]])
      string(REGEX MATCH "${_pat}" _libged_view_aet_cmd_direct_hit
	"${_libged_view_aet_cmd_contents}")
      if(_libged_view_aet_cmd_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/aet.c reintroduced direct BSG AET reads: ${_libged_view_aet_cmd_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_brep_tikz_cmd "${BRLCAD_SOURCE_DIR}/src/libged/brep/tikz.cpp")
  if(EXISTS "${_libged_brep_tikz_cmd}")
    file(READ "${_libged_brep_tikz_cmd}" _libged_brep_tikz_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_aet_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_brep_tikz_token_hit
	"${_libged_brep_tikz_cmd_contents}")
      if(NOT _libged_brep_tikz_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/brep/tikz.cpp must route TikZ AET output through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_aet([^A-Za-z0-9_]|$)]]
      _libged_brep_tikz_direct_hit "${_libged_brep_tikz_cmd_contents}")
    if(_libged_brep_tikz_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/brep/tikz.cpp reintroduced direct BSG TikZ AET output: ${_libged_brep_tikz_direct_hit}")
    endif()
  endif()

  set(_libged_select_cmd "${BRLCAD_SOURCE_DIR}/src/libged/select/select.c")
  if(EXISTS "${_libged_select_cmd}")
    file(READ "${_libged_select_cmd}" _libged_select_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_model2view_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_select_cmd_token_hit
	"${_libged_select_cmd_contents}")
      if(NOT _libged_select_cmd_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/select/select.c must route selection model-to-view reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
      _libged_select_cmd_direct_hit "${_libged_select_cmd_contents}")
    if(_libged_select_cmd_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/select/select.c reintroduced direct BSG selection model-to-view reads: ${_libged_select_cmd_direct_hit}")
    endif()
  endif()

  set(_libged_polyclip_impl "${BRLCAD_SOURCE_DIR}/src/libged/polyclip.cpp")
  if(EXISTS "${_libged_polyclip_impl}")
    file(READ "${_libged_polyclip_impl}" _libged_polyclip_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_view2model_from_bsg]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_plane_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_polyclip_token_hit
	"${_libged_polyclip_contents}")
      if(NOT _libged_polyclip_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/polyclip.cpp must route polygon view transform reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_plane([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libged_polyclip_direct_hit
	"${_libged_polyclip_contents}")
      if(_libged_polyclip_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/polyclip.cpp reintroduced direct BSG polygon view transform reads: ${_libged_polyclip_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_view_snap_cmd "${BRLCAD_SOURCE_DIR}/src/libged/view/snap.c")
  if(EXISTS "${_libged_view_snap_cmd}")
    file(READ "${_libged_view_snap_cmd}" _libged_view_snap_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_view2model_from_bsg]]
	[[rt_view_snap_tolerance_factor_from_bsg]]
	[[rt_view_snap_tolerance_factor_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_view_snap_cmd_token_hit
	"${_libged_view_snap_cmd_contents}")
      if(NOT _libged_view_snap_cmd_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/snap.c must route snap view transform reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_(set_)?snap_tolerance_factor([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libged_view_snap_cmd_direct_hit
	"${_libged_view_snap_cmd_contents}")
      if(_libged_view_snap_cmd_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/snap.c reintroduced direct BSG snap view transform reads: ${_libged_view_snap_cmd_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_view_data_lines_cmd "${BRLCAD_SOURCE_DIR}/src/libged/view/data_lines.c")
  if(EXISTS "${_libged_view_data_lines_cmd}")
    file(READ "${_libged_view_data_lines_cmd}" _libged_view_data_lines_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_snap_lines_from_bsg]]
	[[rt_view_snap_lines_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_view_data_lines_cmd_token_hit
	"${_libged_view_data_lines_cmd_contents}")
      if(NOT _libged_view_data_lines_cmd_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/data_lines.c must route snap-line view-state access through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg/view_state\.h]]
	[[(^|[^A-Za-z0-9_])bsg_view_(set_)?snap_lines([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libged_view_data_lines_cmd_direct_hit
	"${_libged_view_data_lines_cmd_contents}")
      if(_libged_view_data_lines_cmd_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/data_lines.c reintroduced direct BSG snap-line view-state access: ${_libged_view_data_lines_cmd_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_snap_semantics_test "${BRLCAD_SOURCE_DIR}/src/libged/tests/draw/snap_semantics.cpp")
  if(EXISTS "${_libged_snap_semantics_test}")
    file(READ "${_libged_snap_semantics_test}" _libged_snap_semantics_test_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_snap_source_flags_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_snap_semantics_token_hit
	"${_libged_snap_semantics_test_contents}")
      if(NOT _libged_snap_semantics_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/tests/draw/snap_semantics.cpp must route snap-state setup through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])bsg_view_set_snap_source_flags([^A-Za-z0-9_]|$)]]
      _libged_snap_semantics_direct_hit "${_libged_snap_semantics_test_contents}")
    if(_libged_snap_semantics_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/tests/draw/snap_semantics.cpp reintroduced direct BSG snap-state setup: ${_libged_snap_semantics_direct_hit}")
    endif()
  endif()

  set(_libged_draw_tests_cmake "${BRLCAD_SOURCE_DIR}/src/libged/tests/draw/CMakeLists.txt")
  if(EXISTS "${_libged_draw_tests_cmake}")
    file(READ "${_libged_draw_tests_cmake}" _libged_draw_tests_cmake_contents)
    string(REGEX MATCH [[brlcad_addexec[ \t\r\n]*\([ \t\r\n]*ged_test_snap_semantics[^)]*librt]]
      _libged_snap_semantics_librt_hit "${_libged_draw_tests_cmake_contents}")
    if(NOT _libged_snap_semantics_librt_hit)
      _brlobol_pivot_guard_fail(
	"ged_test_snap_semantics must link librt for rt/view_legacy_bsg.h snap-state setup")
    endif()
  endif()

  set(_libged_faceplate_cmd "${BRLCAD_SOURCE_DIR}/src/libged/view/faceplate/faceplate.c")
  if(EXISTS "${_libged_faceplate_cmd}")
    file(READ "${_libged_faceplate_cmd}" _libged_faceplate_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_framebuffer_mode_from_bsg]]
	[[rt_view_framebuffer_mode_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_faceplate_fb_token_hit
	"${_libged_faceplate_cmd_contents}")
      if(NOT _libged_faceplate_fb_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/faceplate/faceplate.c must route framebuffer mode view-state access through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])bsg_view_(set_)?framebuffer_mode([^A-Za-z0-9_]|$)]]
      _libged_faceplate_fb_direct_hit "${_libged_faceplate_cmd_contents}")
    if(_libged_faceplate_fb_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/faceplate/faceplate.c reintroduced direct BSG framebuffer mode view-state access: ${_libged_faceplate_fb_direct_hit}")
    endif()
  endif()

  foreach(_rel
      src/libged/view/faceplate/faceplate.c
      src/libtclcad/view/faceplate.c
      src/qged/plugins/view/settings/CADViewSettings.cpp)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(EXISTS "${_file}")
      file(READ "${_file}" _contents)
      foreach(_token
	  [[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	  [[rt_view_center_dot_from_bsg]]
	  [[rt_view_center_dot_set_bsg]]
	  [[rt_view_scale_overlay_from_bsg]]
	  [[rt_view_scale_overlay_set_bsg]]
	  [[rt_view_params_from_bsg]]
	  [[rt_view_params_set_bsg]])
	string(REGEX MATCH "${_token}" _faceplate_state_token_hit
	  "${_contents}")
	if(NOT _faceplate_state_token_hit)
	  _brlobol_pivot_guard_fail(
	    "${_rel} must route faceplate scalar state through rt/view_legacy_bsg.h token ${_token}")
	endif()
      endforeach()
      string(REGEX MATCH [[(^|[^A-Za-z0-9_])bsg_view_(center_dot|scale_state|params)_(get|set)([^A-Za-z0-9_]|$)]]
	_faceplate_state_direct_hit "${_contents}")
      if(_faceplate_state_direct_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced direct BSG faceplate scalar-state access: ${_faceplate_state_direct_hit}")
      endif()
    endif()
  endforeach()

  set(_qged_view_settings "${BRLCAD_SOURCE_DIR}/src/qged/plugins/view/settings/CADViewSettings.cpp")
  if(EXISTS "${_qged_view_settings}")
    file(READ "${_qged_view_settings}" _qged_view_settings_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_adc_from_bsg]]
	[[rt_view_adc_set_bsg]]
	[[rt_view_grid_from_bsg]]
	[[rt_view_grid_set_bsg]]
	[[rt_view_model_axes_from_bsg]]
	[[rt_view_model_axes_set_bsg]]
	[[rt_view_view_axes_from_bsg]]
	[[rt_view_view_axes_set_bsg]])
      string(REGEX MATCH "${_token}" _qged_view_settings_token_hit
	"${_qged_view_settings_contents}")
      if(NOT _qged_view_settings_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/qged/plugins/view/settings/CADViewSettings.cpp must route top-level faceplate state through rt/view_legacy_bsg.h token ${_token}")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])bsg_view_(adc|grid|model_axes|view_axes)_(get|set)([^A-Za-z0-9_]|$)]]
      _qged_view_settings_direct_hit "${_qged_view_settings_contents}")
    if(_qged_view_settings_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/qged/plugins/view/settings/CADViewSettings.cpp reintroduced direct BSG top-level faceplate state access: ${_qged_view_settings_direct_hit}")
    endif()
  endif()

  set(_libged_fbclear_cmd "${BRLCAD_SOURCE_DIR}/src/libged/fbclear/fbclear.c")
  if(EXISTS "${_libged_fbclear_cmd}")
    file(READ "${_libged_fbclear_cmd}" _libged_fbclear_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_framebuffer_mode_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_fbclear_fb_token_hit
	"${_libged_fbclear_cmd_contents}")
      if(NOT _libged_fbclear_fb_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/fbclear/fbclear.c must route framebuffer mode view-state access through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg/view_state\.h]]
	[[(^|[^A-Za-z0-9_])bsg_view_(set_)?framebuffer_mode([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libged_fbclear_fb_direct_hit
	"${_libged_fbclear_cmd_contents}")
      if(_libged_fbclear_fb_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/fbclear/fbclear.c reintroduced direct BSG framebuffer mode view-state access: ${_libged_fbclear_fb_direct_hit}")
      endif()
    endforeach()
  endif()

  foreach(_rel
      src/libdm/view.c
      src/libdm/swrast/fb-swrast.cpp
      src/libdm/qtgl/fb-qtgl.cpp)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(EXISTS "${_file}")
      file(READ "${_file}" _contents)
      foreach(_token
	  [[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	  [[rt_view_framebuffer_mode]])
	string(REGEX MATCH "${_token}" _libdm_fb_mode_token_hit "${_contents}")
	if(NOT _libdm_fb_mode_token_hit)
	  _brlobol_pivot_guard_fail(
	    "${_rel} must route framebuffer mode view-state access through rt/view_legacy_bsg.h")
	endif()
      endforeach()
      string(REGEX MATCH [[(^|[^A-Za-z0-9_])bsg_view_(set_)?framebuffer_mode([^A-Za-z0-9_]|$)]]
	_libdm_fb_mode_direct_hit "${_contents}")
      if(_libdm_fb_mode_direct_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced direct BSG framebuffer mode view-state access: ${_libdm_fb_mode_direct_hit}")
      endif()
    endif()
  endforeach()

  set(_libged_zap_cmd "${BRLCAD_SOURCE_DIR}/src/libged/zap/zap2.cpp")
  if(EXISTS "${_libged_zap_cmd}")
    file(READ "${_libged_zap_cmd}" _libged_zap_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_cleared_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_zap_cleared_token_hit
	"${_libged_zap_cmd_contents}")
      if(NOT _libged_zap_cleared_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/zap/zap2.cpp must route cleared-state view-state access through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg/view_state\.h]]
	[[(^|[^A-Za-z0-9_])bsg_view_(set_)?cleared([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libged_zap_cleared_direct_hit
	"${_libged_zap_cmd_contents}")
      if(_libged_zap_cleared_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/zap/zap2.cpp reintroduced direct BSG cleared-state view-state access: ${_libged_zap_cleared_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_gtools_gsh_impl "${BRLCAD_SOURCE_DIR}/src/gtools/gsh/gsh.cpp")
  if(EXISTS "${_gtools_gsh_impl}")
    file(READ "${_gtools_gsh_impl}" _gtools_gsh_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_cleared_from_bsg]]
	[[rt_view_cleared_set_bsg]])
      string(REGEX MATCH "${_token}" _gtools_gsh_cleared_token_hit
	"${_gtools_gsh_contents}")
      if(NOT _gtools_gsh_cleared_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/gtools/gsh/gsh.cpp must route cleared-state view-state access through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])bsg_view_(set_)?cleared([^A-Za-z0-9_]|$)]]
      _gtools_gsh_cleared_direct_hit "${_gtools_gsh_contents}")
    if(_gtools_gsh_cleared_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/gtools/gsh/gsh.cpp reintroduced direct BSG cleared-state view-state access: ${_gtools_gsh_cleared_direct_hit}")
    endif()
  endif()

  set(_libged_view_polygons_cmd
    "${BRLCAD_SOURCE_DIR}/src/libged/view/polygons.c")
  if(EXISTS "${_libged_view_polygons_cmd}")
    file(READ "${_libged_view_polygons_cmd}" _libged_view_polygons_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_view_polygons_cmd_token_hit
	"${_libged_view_polygons_cmd_contents}")
      if(NOT _libged_view_polygons_cmd_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/polygons.c must route polygon scale reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]]
      _libged_view_polygons_cmd_direct_hit
      "${_libged_view_polygons_cmd_contents}")
    if(_libged_view_polygons_cmd_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/polygons.c reintroduced direct BSG polygon scale reads: ${_libged_view_polygons_cmd_direct_hit}")
    endif()
  endif()

  set(_libged_view_vz_cmd "${BRLCAD_SOURCE_DIR}/src/libged/view/view.c")
  if(EXISTS "${_libged_view_vz_cmd}")
    file(READ "${_libged_view_vz_cmd}" _libged_view_vz_cmd_contents)
    string(REGEX MATCH [[rt_view_model2view_from_bsg]]
      _libged_view_vz_cmd_token_hit "${_libged_view_vz_cmd_contents}")
    if(NOT _libged_view_vz_cmd_token_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/view.c must route vZ model-to-view reads through rt/view_legacy_bsg.h")
    endif()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
      _libged_view_vz_cmd_direct_hit "${_libged_view_vz_cmd_contents}")
    if(_libged_view_vz_cmd_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/view.c reintroduced direct BSG vZ model-to-view reads: ${_libged_view_vz_cmd_direct_hit}")
    endif()
  endif()

  foreach(_rel
      src/libged/pipe/pipe.c
      src/libged/metaball/metaball.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for primitive edit transform adapter checks")
      continue()
    endif()
    file(READ "${_file}" _libged_prim_edit_transform_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_view2model_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_prim_edit_transform_token_hit
	"${_libged_prim_edit_transform_contents}")
      if(NOT _libged_prim_edit_transform_token_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must route primitive edit transform reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libged_prim_edit_transform_direct_hit
	"${_libged_prim_edit_transform_contents}")
      if(_libged_prim_edit_transform_direct_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced direct BSG primitive edit transform reads: ${_libged_prim_edit_transform_direct_hit}")
      endif()
    endforeach()
  endforeach()

  foreach(_rel
      src/libged/bot/edbot.c
      src/libged/move_arb_edge/move_arb_edge.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for subentity selection transform adapter checks")
      continue()
    endif()
    file(READ "${_file}" _libged_subentity_transform_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_model2view_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_subentity_transform_token_hit
	"${_libged_subentity_transform_contents}")
      if(NOT _libged_subentity_transform_token_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must route subentity selection model-to-view reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
      _libged_subentity_transform_direct_hit
      "${_libged_subentity_transform_contents}")
    if(_libged_subentity_transform_direct_hit)
      _brlobol_pivot_guard_fail(
	"${_rel} reintroduced direct BSG subentity selection model-to-view reads: ${_libged_subentity_transform_direct_hit}")
    endif()
  endforeach()

  set(_libged_rot_point_cmd "${BRLCAD_SOURCE_DIR}/src/libged/rot_point/rot_point.c")
  if(EXISTS "${_libged_rot_point_cmd}")
    file(READ "${_libged_rot_point_cmd}" _libged_rot_point_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_rotation_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_rot_point_cmd_token_hit
	"${_libged_rot_point_cmd_contents}")
      if(NOT _libged_rot_point_cmd_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/rot_point/rot_point.c must route rotation reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_rotation([^A-Za-z0-9_]|$)]]
      _libged_rot_point_cmd_direct_hit "${_libged_rot_point_cmd_contents}")
    if(_libged_rot_point_cmd_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/rot_point/rot_point.c reintroduced direct BSG rotation reads: ${_libged_rot_point_cmd_direct_hit}")
    endif()
  endif()

  set(_libged_view_labels_cmd "${BRLCAD_SOURCE_DIR}/src/libged/view/labels.c")
  if(EXISTS "${_libged_view_labels_cmd}")
    file(READ "${_libged_view_labels_cmd}" _libged_view_labels_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_screen_to_view_from_bsg]]
	[[rt_view_view2model_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_view_labels_cmd_token_hit
	"${_libged_view_labels_cmd_contents}")
      if(NOT _libged_view_labels_cmd_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/labels.c must route label view-to-model reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[((^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)|(^|[^A-Za-z0-9_])bsg_screen_to_view([^A-Za-z0-9_]|$))]]
      _libged_view_labels_cmd_direct_hit "${_libged_view_labels_cmd_contents}")
    if(_libged_view_labels_cmd_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/labels.c reintroduced direct BSG label view-to-model reads: ${_libged_view_labels_cmd_direct_hit}")
    endif()
  endif()

  foreach(_rel
      src/libged/view/ypr.c
      src/libged/rmat/rmat.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for view rotation query adapter checks")
      continue()
    endif()
    file(READ "${_file}" _libged_view_rotation_query_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_rotation_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_view_rotation_query_token_hit
	"${_libged_view_rotation_query_contents}")
      if(NOT _libged_view_rotation_query_token_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must route rotation query output through rt/view_legacy_bsg.h")
      endif()
    endforeach()
  endforeach()

  set(_libged_view_ypr_cmd "${BRLCAD_SOURCE_DIR}/src/libged/view/ypr.c")
  if(EXISTS "${_libged_view_ypr_cmd}")
    file(READ "${_libged_view_ypr_cmd}" _libged_view_ypr_cmd_contents)
    string(REGEX MATCH [[bn_mat_trn[ \t\r\n]*\([^,;]*,[^;]*gv_rotation]]
      _libged_view_ypr_cmd_direct_hit "${_libged_view_ypr_cmd_contents}")
    if(_libged_view_ypr_cmd_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/ypr.c reintroduced direct BSG query output: ${_libged_view_ypr_cmd_direct_hit}")
    endif()
    string(REGEX MATCH [[bn_mat_trn[ \t\r\n]*\([^;]*gv_rotation]]
      _libged_view_ypr_cmd_set_direct_hit "${_libged_view_ypr_cmd_contents}")
    if(_libged_view_ypr_cmd_set_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/ypr.c reintroduced direct BSG rotation writes: ${_libged_view_ypr_cmd_set_direct_hit}")
    endif()
  endif()

  set(_libged_rmat_cmd "${BRLCAD_SOURCE_DIR}/src/libged/rmat/rmat.c")
  if(EXISTS "${_libged_rmat_cmd}")
    file(READ "${_libged_rmat_cmd}" _libged_rmat_cmd_contents)
    string(REGEX MATCH [[bn_encode_mat[^\n;]*gv_rotation]]
      _libged_rmat_cmd_direct_hit "${_libged_rmat_cmd_contents}")
    if(_libged_rmat_cmd_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/rmat/rmat.c reintroduced direct BSG query output: ${_libged_rmat_cmd_direct_hit}")
    endif()
    string(REGEX MATCH [[MAT_COPY[^\n;]*gv_rotation]]
      _libged_rmat_cmd_set_direct_hit "${_libged_rmat_cmd_contents}")
    if(_libged_rmat_cmd_set_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/rmat/rmat.c reintroduced direct BSG rotation writes: ${_libged_rmat_cmd_set_direct_hit}")
    endif()
  endif()

  foreach(_rel
      src/libged/orient/orient.c
      src/libged/setview/setview.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(EXISTS "${_file}")
      file(READ "${_file}" _libged_rotation_setter_contents)
      foreach(_token
	  [[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	  [[rt_view_rotation_set_bsg]])
	string(REGEX MATCH "${_token}" _libged_rotation_setter_token_hit
	  "${_libged_rotation_setter_contents}")
	if(NOT _libged_rotation_setter_token_hit)
	  _brlobol_pivot_guard_fail(
	    "${_rel} must route rotation writes through rt/view_legacy_bsg.h")
	endif()
      endforeach()
      foreach(_pat
	  [[quat_quat2mat[^\n;]*gv_rotation]]
	  [[bn_mat_angles[^\n;]*gv_rotation]]
	  [[MAT_COPY[^\n;]*gv_rotation]])
	string(REGEX MATCH "${_pat}" _libged_rotation_setter_direct_hit
	  "${_libged_rotation_setter_contents}")
	if(_libged_rotation_setter_direct_hit)
	  _brlobol_pivot_guard_fail(
	    "${_rel} reintroduced direct BSG rotation writes: ${_libged_rotation_setter_direct_hit}")
	endif()
      endforeach()
    endif()
  endforeach()

  set(_libged_pmat_cmd "${BRLCAD_SOURCE_DIR}/src/libged/pmat/pmat.c")
  if(EXISTS "${_libged_pmat_cmd}")
    file(READ "${_libged_pmat_cmd}" _libged_pmat_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_pmat_from_bsg]]
	[[rt_view_pmat_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_pmat_cmd_token_hit
	"${_libged_pmat_cmd_contents}")
      if(NOT _libged_pmat_cmd_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/pmat/pmat.c must route projection-matrix access through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[bn_encode_mat[^\n;]*gv_pmat|MAT_COPY[^\n;]*gv_pmat]]
      _libged_pmat_cmd_direct_hit "${_libged_pmat_cmd_contents}")
    if(_libged_pmat_cmd_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/pmat/pmat.c reintroduced direct BSG projection-matrix access: ${_libged_pmat_cmd_direct_hit}")
    endif()
  endif()

  set(_libged_pmodel2view_cmd "${BRLCAD_SOURCE_DIR}/src/libged/pmodel2view/pmodel2view.c")
  if(EXISTS "${_libged_pmodel2view_cmd}")
    file(READ "${_libged_pmodel2view_cmd}" _libged_pmodel2view_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_pmodel2view_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_pmodel2view_cmd_token_hit
	"${_libged_pmodel2view_cmd_contents}")
      if(NOT _libged_pmodel2view_cmd_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/pmodel2view/pmodel2view.c must route projected-model-to-view query output through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[bn_encode_mat[^\n;]*gv_pmodel2view]]
      _libged_pmodel2view_cmd_direct_hit "${_libged_pmodel2view_cmd_contents}")
    if(_libged_pmodel2view_cmd_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/pmodel2view/pmodel2view.c reintroduced direct BSG projected-model-to-view query output: ${_libged_pmodel2view_cmd_direct_hit}")
    endif()
  endif()

  set(_libged_view_quat_cmd "${BRLCAD_SOURCE_DIR}/src/libged/view/quat.c")
  if(EXISTS "${_libged_view_quat_cmd}")
    file(READ "${_libged_view_quat_cmd}" _libged_view_quat_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_orientation_quat_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_view_quat_cmd_token_hit
	"${_libged_view_quat_cmd_contents}")
      if(NOT _libged_view_quat_cmd_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/quat.c must route quaternion query output through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[quat_mat2quat[ \t\r\n]*\([^;]*gv_rotation]]
      _libged_view_quat_cmd_direct_hit "${_libged_view_quat_cmd_contents}")
    if(_libged_view_quat_cmd_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/quat.c reintroduced direct BSG query output: ${_libged_view_quat_cmd_direct_hit}")
    endif()
    foreach(_token
	[[rt_view_rotation_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_view_quat_cmd_set_token_hit
	"${_libged_view_quat_cmd_contents}")
      if(NOT _libged_view_quat_cmd_set_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/quat.c must route quaternion rotation writes through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[quat_quat2mat[^\n;]*gv_rotation]]
      _libged_view_quat_cmd_set_direct_hit "${_libged_view_quat_cmd_contents}")
    if(_libged_view_quat_cmd_set_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/quat.c reintroduced direct BSG rotation writes: ${_libged_view_quat_cmd_set_direct_hit}")
    endif()
  endif()

  set(_libged_perspective_cmd
    "${BRLCAD_SOURCE_DIR}/src/libged/perspective/perspective.c")
  if(EXISTS "${_libged_perspective_cmd}")
    file(READ "${_libged_perspective_cmd}" _libged_perspective_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_perspective_from_bsg]]
	[[rt_view_perspective_set_bsg]]
	[[rt_view_pmat_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_perspective_cmd_token_hit
	"${_libged_perspective_cmd_contents}")
      if(NOT _libged_perspective_cmd_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/perspective/perspective.c must route query output plus perspective/projection-matrix writes through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[gedp->[ \t\r\n]*ged_gvp->[ \t\r\n]*gv_perspective[ \t\r\n]*=]]
	[[bu_vls_printf[^\n;]*gv_perspective]]
	[[if[ \t]*\([^\n;]*gv_perspective]]
	[[persp_mat[^\n;]*gv_perspective]]
	[[persp_mat[^\n;]*gv_pmat]]
	[[MAT_COPY[^\n;]*gv_pmat]])
      string(REGEX MATCH "${_pat}" _libged_perspective_cmd_direct_hit
	"${_libged_perspective_cmd_contents}")
      if(_libged_perspective_cmd_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/perspective/perspective.c reintroduced direct BSG perspective reads/writes: ${_libged_perspective_cmd_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_eye_pos_cmd "${BRLCAD_SOURCE_DIR}/src/libged/eye_pos/eye_pos.c")
  if(EXISTS "${_libged_eye_pos_cmd}")
    file(READ "${_libged_eye_pos_cmd}" _libged_eye_pos_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_eye_pos_from_bsg]]
	[[rt_view_eye_pos_set_bsg]]
	[[rt_view_pmat_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_eye_pos_cmd_token_hit
	"${_libged_eye_pos_cmd_contents}")
      if(NOT _libged_eye_pos_cmd_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/eye_pos/eye_pos.c must route eye-position and projection-matrix access through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[VSCALE[ \t\r\n]*\([^,;]*,[^;]*gv_eye_pos]]
	[[VMOVE[^\n;]*gv_eye_pos]]
	[[mike_persp_mat[^\n;]*gv_eye_pos]]
	[[mike_persp_mat[^\n;]*gv_pmat]])
      string(REGEX MATCH "${_pat}" _libged_eye_pos_cmd_direct_hit
	"${_libged_eye_pos_cmd_contents}")
      if(_libged_eye_pos_cmd_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/eye_pos/eye_pos.c reintroduced direct BSG eye-position/projection-matrix access: ${_libged_eye_pos_cmd_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_view_cmd_impl "${BRLCAD_SOURCE_DIR}/src/libged/view/view.c")
  if(EXISTS "${_libged_view_cmd_impl}")
    file(READ "${_libged_view_cmd_impl}" _libged_view_cmd_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_info_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_view_cmd_token_hit
	"${_libged_view_cmd_contents}")
      if(NOT _libged_view_cmd_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/view.c must route width/height query output through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[bu_vls_printf[^\n;]*gv_width]]
	[[bu_vls_printf[^\n;]*gv_height]])
      string(REGEX MATCH "${_pat}" _libged_view_cmd_direct_hit
	"${_libged_view_cmd_contents}")
      if(_libged_view_cmd_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/view.c reintroduced direct BSG width/height query output: ${_libged_view_cmd_direct_hit}")
      endif()
    endforeach()
  endif()

  foreach(_rel
      src/libtclcad/commands.c
      src/mged/adc.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for TclCAD/MGED RT view numeric policy checks")
      continue()
    endif()
    file(READ "${_file}" _tcl_mged_rt_view_contents)
    foreach(_token
	"rt/view.h"
	"RT_VIEW_")
      string(FIND "${_tcl_mged_rt_view_contents}" "${_token}"
	_tcl_mged_rt_view_token_idx)
      if(_tcl_mged_rt_view_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "${_rel} must use RT-owned view numeric policy via ${_token}")
      endif()
    endforeach()
    if("${_rel}" STREQUAL "src/mged/adc.c")
      string(FIND "${_tcl_mged_rt_view_contents}" "RT_INV_VIEW"
	_mged_adc_rt_inv_view_idx)
      if(_mged_adc_rt_inv_view_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/mged/adc.c must use RT_INV_VIEW for ADC view-unit conversion")
      endif()
    endif()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])BSG_VIEW_(MIN|MAX|RANGE)([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])BSG_INV_(VIEW|4096)([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])INV_BV([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _tcl_mged_rt_view_hit
	"${_tcl_mged_rt_view_contents}")
      if(_tcl_mged_rt_view_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must not use BSG-owned spellings for RT view numeric policy: ${_tcl_mged_rt_view_hit}")
      endif()
    endforeach()
  endforeach()

  foreach(_rel
      src/mged/mged_dm.h
      src/mged/axes.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for shared MGED RT inverse-view constant checks")
      continue()
    endif()
    file(READ "${_file}" _mged_rt_inverse_view_contents)
    foreach(_token
	"rt/view.h"
	"RT_INV_VIEW")
      string(FIND "${_mged_rt_inverse_view_contents}" "${_token}"
	_mged_rt_inverse_view_token_idx)
      if(_mged_rt_inverse_view_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "${_rel} must use RT_INV_VIEW via rt/view.h for view-unit conversion")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])BSG_INV_(VIEW|4096)([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])INV_BV([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _mged_rt_inverse_view_hit
	"${_mged_rt_inverse_view_contents}")
      if(_mged_rt_inverse_view_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must not use BSG-owned spellings for inverse RT view-unit conversion: ${_mged_rt_inverse_view_hit}")
      endif()
    endforeach()
    if("${_rel}" STREQUAL "src/mged/mged_dm.h")
      foreach(_token
	  [[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	  [[rt_view_size_from_bsg]]
	  [[rt_view_scale_from_bsg]]
	  [[rt_view_settings_shared_bsg]])
	string(REGEX MATCH "${_token}" _mged_dm_view_macro_token_hit
	  "${_mged_rt_inverse_view_contents}")
	if(NOT _mged_dm_view_macro_token_hit)
	  _brlobol_pivot_guard_fail(
	    "src/mged/mged_dm.h must route VIEWSIZE/VIEWFACTOR macros through rt/view_legacy_bsg.h")
	endif()
      endforeach()
      foreach(_pat
	  [[#[ \t]*define[ \t]+VIEWSIZE[^\n]*gv_size]]
	  [[#[ \t]*define[ \t]+VIEWFACTOR[^\n]*gv_scale]]
	  [[(^|[^A-Za-z0-9_])bsg_view_settings_shared([^A-Za-z0-9_]|$)]])
	string(REGEX MATCH "${_pat}" _mged_dm_view_macro_direct_hit
	  "${_mged_rt_inverse_view_contents}")
	if(_mged_dm_view_macro_direct_hit)
	  _brlobol_pivot_guard_fail(
	    "src/mged/mged_dm.h reintroduced direct BSG view-size/scale macro reads: ${_mged_dm_view_macro_direct_hit}")
	endif()
      endforeach()
    endif()
  endforeach()

  set(_mged_attach_impl "${BRLCAD_SOURCE_DIR}/src/mged/attach.c")
  if(EXISTS "${_mged_attach_impl}")
    file(READ "${_mged_attach_impl}" _mged_attach_contents)
    foreach(_token
	"rt/view.h"
	"rt/view_legacy_bsg.h"
	"rt_view_scale_storage_from_bsg"
	"RT_VIEW_MIN"
	"RT_VIEW_MAX")
      string(FIND "${_mged_attach_contents}" "${_token}" _mged_attach_token_idx)
      if(_mged_attach_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/mged/attach.c must use RT-owned view adapters/constants ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])BSG_VIEW_(MIN|MAX)([^A-Za-z0-9_]|$)]]
	[[dm_set_vp[ \t\r\n]*\([^;]*gv_scale]])
      string(REGEX MATCH "${_pat}" _mged_attach_bsg_view_hit
	"${_mged_attach_contents}")
      if(_mged_attach_bsg_view_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/attach.c must not use BSG-owned spellings for DM view state: ${_mged_attach_bsg_view_hit}")
      endif()
    endforeach()
  endif()

  set(_mged_chgview_impl "${BRLCAD_SOURCE_DIR}/src/mged/chgview.c")
  if(EXISTS "${_mged_chgview_impl}")
    file(READ "${_mged_chgview_impl}" _mged_chgview_contents)
    foreach(_token
	"rt/view.h"
	"RT_VIEW_MIN_SIZE")
      string(FIND "${_mged_chgview_contents}" "${_token}" _mged_chgview_token_idx)
      if(_mged_chgview_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/mged/chgview.c must use RT-owned view size clamp ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])BSG_MINVIEWSIZE([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _mged_chgview_bsg_size_hit
	"${_mged_chgview_contents}")
      if(_mged_chgview_bsg_size_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/chgview.c must not use BSG-owned spellings for view size limits: ${_mged_chgview_bsg_size_hit}")
      endif()
    endforeach()
  endif()

  foreach(_rel
      src/mged/menu.c
      src/mged/titles.c
      src/mged/scroll.c
      src/mged/usepen.c
      src/mged/doevent.c
      src/mged/dm-generic.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail(
	"${_rel} is required for MGED HUD/menu/input RT view-coordinate checks")
      continue()
    endif()
    file(READ "${_file}" _mged_hud_menu_rt_view_contents)
    foreach(_token
	"rt/view.h"
	"RT_VIEW_")
      string(FIND "${_mged_hud_menu_rt_view_contents}" "${_token}"
	_mged_hud_menu_rt_view_token_idx)
      if(_mged_hud_menu_rt_view_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "${_rel} must use RT-owned view-coordinate constants via ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])BSG_VIEW_(MIN|MAX|RANGE)([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])BSG_INV_(VIEW|4096)([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])INV_BV([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _mged_hud_menu_rt_view_hit
	"${_mged_hud_menu_rt_view_contents}")
      if(_mged_hud_menu_rt_view_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must not use BSG-owned spellings for normalized RT view coordinates: ${_mged_hud_menu_rt_view_hit}")
      endif()
    endforeach()
  endforeach()

  set(_mged_titles_impl "${BRLCAD_SOURCE_DIR}/src/mged/titles.c")
  if(EXISTS "${_mged_titles_impl}")
    file(READ "${_mged_titles_impl}" _mged_titles_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_info_from_bsg]]
	[[rt_view_center_from_bsg]]
	[[rt_view_aet_from_bsg]]
	[[rt_view_perspective_from_bsg]]
	[[rt_view_scale_from_bsg]])
      string(REGEX MATCH "${_token}" _mged_titles_adapter_hit
	"${_mged_titles_contents}")
      if(NOT _mged_titles_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/titles.c must route title/status view reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_center([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_size([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_aet([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_perspective([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _mged_titles_direct_view_hit
	"${_mged_titles_contents}")
      if(_mged_titles_direct_view_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/titles.c reintroduced direct BSG title/status view reads: ${_mged_titles_direct_view_hit}")
      endif()
    endforeach()
  endif()

  set(_mged_chgmodel_impl "${BRLCAD_SOURCE_DIR}/src/mged/chgmodel.c")
  if(EXISTS "${_mged_chgmodel_impl}")
    file(READ "${_mged_chgmodel_impl}" _mged_chgmodel_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_center_from_bsg]]
	[[rt_view_scale_from_bsg]])
      string(REGEX MATCH "${_token}" _mged_chgmodel_adapter_hit
	"${_mged_chgmodel_contents}")
      if(NOT _mged_chgmodel_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/chgmodel.c must route make default view center/scale reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_center([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _mged_chgmodel_direct_hit
	"${_mged_chgmodel_contents}")
      if(_mged_chgmodel_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/chgmodel.c reintroduced direct BSG make default view reads: ${_mged_chgmodel_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_mged_plot_impl "${BRLCAD_SOURCE_DIR}/src/mged/plot.c")
  if(EXISTS "${_mged_plot_impl}")
    file(READ "${_mged_plot_impl}" _mged_plot_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_rotation_from_bsg]])
      string(REGEX MATCH "${_token}" _mged_plot_adapter_hit
	"${_mged_plot_contents}")
      if(NOT _mged_plot_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/plot.c must route area view rotation reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_rotation([^A-Za-z0-9_]|$)]]
      _mged_plot_direct_hit "${_mged_plot_contents}")
    if(_mged_plot_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/mged/plot.c reintroduced direct BSG area view rotation reads: ${_mged_plot_direct_hit}")
    endif()
  endif()

  set(_mged_scroll_impl "${BRLCAD_SOURCE_DIR}/src/mged/scroll.c")
  if(EXISTS "${_mged_scroll_impl}")
    file(READ "${_mged_scroll_impl}" _mged_scroll_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_from_bsg]])
      string(REGEX MATCH "${_token}" _mged_scroll_adapter_hit
	"${_mged_scroll_contents}")
      if(NOT _mged_scroll_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/scroll.c must route slider view-scale reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]]
      _mged_scroll_direct_hit "${_mged_scroll_contents}")
    if(_mged_scroll_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/mged/scroll.c reintroduced direct BSG slider view-scale reads: ${_mged_scroll_direct_hit}")
    endif()
  endif()

  set(_mged_axes_impl "${BRLCAD_SOURCE_DIR}/src/mged/axes.c")
  if(EXISTS "${_mged_axes_impl}")
    file(READ "${_mged_axes_impl}" _mged_axes_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_info_from_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_rotation_from_bsg]])
      string(REGEX MATCH "${_token}" _mged_axes_adapter_hit
	"${_mged_axes_contents}")
      if(NOT _mged_axes_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/axes.c must route axes overlay view reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_size([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_rotation([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _mged_axes_direct_hit
	"${_mged_axes_contents}")
      if(_mged_axes_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/axes.c reintroduced direct BSG axes overlay view reads: ${_mged_axes_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_mged_rect_impl "${BRLCAD_SOURCE_DIR}/src/mged/rect.c")
  if(EXISTS "${_mged_rect_impl}")
    file(READ "${_mged_rect_impl}" _mged_rect_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_center_from_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_view2model_from_bsg]])
      string(REGEX MATCH "${_token}" _mged_rect_adapter_hit
	"${_mged_rect_contents}")
      if(NOT _mged_rect_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/rect.c must route rectangle zoom view reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_center([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _mged_rect_direct_hit
	"${_mged_rect_contents}")
      if(_mged_rect_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/rect.c reintroduced direct BSG rectangle zoom view reads: ${_mged_rect_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_mged_edpipe_impl "${BRLCAD_SOURCE_DIR}/src/mged/edpipe.c")
  if(EXISTS "${_mged_edpipe_impl}")
    file(READ "${_mged_edpipe_impl}" _mged_edpipe_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_view2model_from_bsg]])
      string(REGEX MATCH "${_token}" _mged_edpipe_adapter_hit
	"${_mged_edpipe_contents}")
      if(NOT _mged_edpipe_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/edpipe.c must route pipe pick direction reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]]
      _mged_edpipe_direct_hit "${_mged_edpipe_contents}")
    if(_mged_edpipe_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/mged/edpipe.c reintroduced direct BSG pipe pick direction reads: ${_mged_edpipe_direct_hit}")
    endif()
  endif()

  set(_mged_usepen_impl "${BRLCAD_SOURCE_DIR}/src/mged/usepen.c")
  if(EXISTS "${_mged_usepen_impl}")
    file(READ "${_mged_usepen_impl}" _mged_usepen_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_center_from_bsg]])
      string(REGEX MATCH "${_token}" _mged_usepen_adapter_hit
	"${_mged_usepen_contents}")
      if(NOT _mged_usepen_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/usepen.c must route view-centered transform reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_center([^A-Za-z0-9_]|$)]]
      _mged_usepen_direct_hit "${_mged_usepen_contents}")
    if(_mged_usepen_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/mged/usepen.c reintroduced direct BSG view-centered transform reads: ${_mged_usepen_direct_hit}")
    endif()
  endif()

  set(_mged_setup_impl "${BRLCAD_SOURCE_DIR}/src/mged/setup.c")
  if(EXISTS "${_mged_setup_impl}")
    file(READ "${_mged_setup_impl}" _mged_setup_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_center_from_bsg]])
      string(REGEX MATCH "${_token}" _mged_setup_adapter_hit
	"${_mged_setup_contents}")
      if(NOT _mged_setup_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/setup.c must route initial view-center reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_center([^A-Za-z0-9_]|$)]]
      _mged_setup_direct_hit "${_mged_setup_contents}")
    if(_mged_setup_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/mged/setup.c reintroduced direct BSG initial view-center reads: ${_mged_setup_direct_hit}")
    endif()
  endif()

  set(_mged_chgview_status_impl "${BRLCAD_SOURCE_DIR}/src/mged/chgview.c")
  if(EXISTS "${_mged_chgview_status_impl}")
    file(READ "${_mged_chgview_status_impl}" _mged_chgview_status_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_initial_scale_from_bsg]]
	[[rt_view_absolute_scale_from_bsg]]
	[[rt_view_center_from_bsg]]
	[[rt_view_rotation_from_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_view2model_from_bsg]])
      string(REGEX MATCH "${_token}" _mged_chgview_status_adapter_hit
	"${_mged_chgview_status_contents}")
      if(NOT _mged_chgview_status_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/chgview.c must route status query view reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[bu_vls_printf[^;]*gv_scale]]
	[[mged_bn_mat_print[^\n;]*gv_(center|rotation|model2view|view2model)]])
      string(REGEX MATCH "${_pat}" _mged_chgview_status_direct_hit
	"${_mged_chgview_status_contents}")
      if(_mged_chgview_status_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/chgview.c reintroduced direct BSG status query view reads: ${_mged_chgview_status_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_mged_grid_impl "${BRLCAD_SOURCE_DIR}/src/mged/grid.c")
  if(EXISTS "${_mged_grid_impl}")
    file(READ "${_mged_grid_impl}" _mged_grid_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_info_from_bsg]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_center_from_bsg]]
	[[rt_view_center_vec_set_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_view2model_from_bsg]])
      string(REGEX MATCH "${_token}" _mged_grid_adapter_hit
	"${_mged_grid_contents}")
      if(NOT _mged_grid_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/grid.c must route grid view read snapshots through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_size([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]]
	[[MAT_DELTAS_GET_NEG[ \t\r\n]*\([^;]*gv_center]]
	[[MAT_DELTAS_VEC_NEG[ \t\r\n]*\([^;]*gv_center]])
      string(REGEX MATCH "${_pat}" _mged_grid_direct_hit
	"${_mged_grid_contents}")
      if(_mged_grid_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/grid.c reintroduced direct BSG grid view reads: ${_mged_grid_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_mged_adc_impl "${BRLCAD_SOURCE_DIR}/src/mged/adc.c")
  if(EXISTS "${_mged_adc_impl}")
    file(READ "${_mged_adc_impl}" _mged_adc_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_view2model_from_bsg]])
      string(REGEX MATCH "${_token}" _mged_adc_adapter_hit
	"${_mged_adc_contents}")
      if(NOT _mged_adc_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/adc.c must route ADC view scale/matrix reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _mged_adc_direct_hit
	"${_mged_adc_contents}")
      if(_mged_adc_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/adc.c reintroduced direct BSG ADC view scale/matrix reads: ${_mged_adc_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_mged_set_impl "${BRLCAD_SOURCE_DIR}/src/mged/set.c")
  if(EXISTS "${_mged_set_impl}")
    file(READ "${_mged_set_impl}" _mged_set_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_center_from_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_perspective_set_bsg]]
	[[rt_view_coord_set_bsg]]
	[[rt_view_rotate_about_set_bsg]])
      string(REGEX MATCH "${_token}" _mged_set_adapter_hit
	"${_mged_set_contents}")
      if(NOT _mged_set_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/set.c must route view reads/simple setters through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_center([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]]
	[[vs_gvp->[ \t\r\n]*gv_perspective[ \t\r\n]*=]]
	[[vs_gvp->[ \t\r\n]*gv_coord[ \t\r\n]*=]]
	[[vs_gvp->[ \t\r\n]*gv_rotate_about[ \t\r\n]*=]])
      string(REGEX MATCH "${_pat}" _mged_set_direct_hit
	"${_mged_set_contents}")
      if(_mged_set_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/set.c reintroduced direct BSG view reads/simple setters: ${_mged_set_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_mged_chgview_edit_reads_impl "${BRLCAD_SOURCE_DIR}/src/mged/chgview.c")
  if(EXISTS "${_mged_chgview_edit_reads_impl}")
    file(READ "${_mged_chgview_edit_reads_impl}" _mged_chgview_edit_reads_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_center_from_bsg]]
	[[rt_view_rotation_from_bsg]]
	[[rt_view_scale_set_bsg]]
	[[rt_view_initial_scale_set_bsg]]
	[[rt_view_absolute_scale_set_bsg]]
	[[rt_view_center_vec_set_bsg]]
	[[rt_view_rotation_set_bsg]]
	[[rt_view_view2model_from_bsg]])
      string(REGEX MATCH "${_token}" _mged_chgview_edit_reads_adapter_hit
	"${_mged_chgview_edit_reads_contents}")
      if(NOT _mged_chgview_edit_reads_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/chgview.c must route edit/slew view read snapshots through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_token
	[[rt_view_coord_from_bsg]]
	[[rt_view_coord_set_bsg]]
	[[rt_view_rotate_about_from_bsg]])
      string(REGEX MATCH "${_token}" _mged_chgview_coord_adapter_hit
	"${_mged_chgview_edit_reads_contents}")
      if(NOT _mged_chgview_coord_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/chgview.c must route edit rotation coord/rotate-about reads and temporary coord writes through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[fval[^\n;]*/[^\n;]*gv_scale]]
	[[last_arr[^\n;]*gv_scale]]
	[[bn_mat_inv[^\n;]*gv_rotation]]
	[[bn_mat_mul[^\n;]*gv_rotation]]
	[[MAT4X3PNT[^\n;]*gv_view2model]]
	[[VSCALE[^\n;]*/[^\n;]*gv_scale]]
	[[MAT_DELTAS_GET_NEG[ \t\r\n]*\([^;]*gv_center]]
	[[gv_a_scale[^\n;]*gv_scale]]
	[[f[ \t]*=[^\n;]*gv_scale[ \t]*/[^\n;]*gv_i_scale]]
	[[bu_vls_printf[^\n;]*gv_a_scale]]
	[[if[ \t]*\([^\n;]*gv_scale[ \t]*<]]
	[[gv_i_scale[ \t]*=[^\n;]*gv_scale]]
	[[vrp->[ \t]*vr_scale[ \t]*=[^\n;]*gv_scale]]
	[[vs_current_view->vr_scale[ \t]*=[^\n;]*gv_scale]]
	[[save_coord[ \t]*=[^\n;]*gv_coord]]
	[[v->[ \t\r\n]*gv_coord[ \t\r\n]*=]]
	[[mged_erot[ \t\r\n]*\([^;]*gv_coord]]
	[[mged_erot[ \t\r\n]*\([^;]*gv_rotate_about]]
	[[view_state->[ \t\r\n]*vs_gvp->[ \t\r\n]*gv_scale[ \t\r\n]*=[^\n;]*view_scale]]
	[[view_state->[ \t\r\n]*vs_gvp->[ \t\r\n]*gv_i_scale[ \t\r\n]*=]]
	[[view_state->[ \t\r\n]*vs_gvp->[ \t\r\n]*gv_a_scale[ \t\r\n]*=]]
	[[MAT_COPY[ \t\r\n]*\([^\n,]*view_state->[ \t\r\n]*vs_gvp->[ \t\r\n]*gv_rotation[^;]*vs_current_view->[ \t\r\n]*vr_rot_mat]]
	[[MAT_COPY[ \t\r\n]*\([^\n,]*view_state->[ \t\r\n]*vs_gvp->[ \t\r\n]*gv_center[^;]*vs_current_view->[ \t\r\n]*vr_tvc_mat]]
	[[view_state->[ \t\r\n]*vs_gvp->[ \t\r\n]*gv_scale[ \t\r\n]*=[^\n;]*vs_current_view->[ \t\r\n]*vr_scale]])
      string(REGEX MATCH "${_pat}" _mged_chgview_edit_reads_direct_hit
	"${_mged_chgview_edit_reads_contents}")
      if(_mged_chgview_edit_reads_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/chgview.c reintroduced direct BSG edit/slew view reads: ${_mged_chgview_edit_reads_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_mged_dm_generic_event_impl "${BRLCAD_SOURCE_DIR}/src/mged/dm-generic.c")
  if(EXISTS "${_mged_dm_generic_event_impl}")
    file(READ "${_mged_dm_generic_event_impl}" _mged_dm_generic_event_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_center_from_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_view2model_from_bsg]]
	[[rt_view_zclip_set_bsg]])
      string(REGEX MATCH "${_token}" _mged_dm_generic_event_adapter_hit
	"${_mged_dm_generic_event_contents}")
      if(NOT _mged_dm_generic_event_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/dm-generic.c must route mouse/ADC input view reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_center([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_set_zclip([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _mged_dm_generic_event_direct_hit
	"${_mged_dm_generic_event_contents}")
      if(_mged_dm_generic_event_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/dm-generic.c reintroduced direct BSG mouse/ADC input view reads: ${_mged_dm_generic_event_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_mged_doevent_impl "${BRLCAD_SOURCE_DIR}/src/mged/doevent.c")
  if(EXISTS "${_mged_doevent_impl}")
    file(READ "${_mged_doevent_impl}" _mged_doevent_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_center_from_bsg]]
	[[rt_view_view2model_from_bsg]])
      string(REGEX MATCH "${_token}" _mged_doevent_adapter_hit
	"${_mged_doevent_contents}")
      if(NOT _mged_doevent_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/doevent.c must route motion-event view reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_center([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _mged_doevent_direct_hit
	"${_mged_doevent_contents}")
      if(_mged_doevent_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/doevent.c reintroduced direct BSG motion-event view reads: ${_mged_doevent_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_mged_main_view_reads_impl "${BRLCAD_SOURCE_DIR}/src/mged/mged.c")
  if(EXISTS "${_mged_main_view_reads_impl}")
    file(READ "${_mged_main_view_reads_impl}" _mged_main_view_reads_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_zclip_set_bsg]])
      string(REGEX MATCH "${_token}" _mged_main_view_reads_adapter_hit
	"${_mged_main_view_reads_contents}")
      if(NOT _mged_main_view_reads_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/mged.c must route edit matrix and rate-knob scale reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_set_zclip([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _mged_main_view_reads_direct_hit
	"${_mged_main_view_reads_contents}")
      if(_mged_main_view_reads_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/mged.c reintroduced direct BSG edit matrix/rate-knob view reads: ${_mged_main_view_reads_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_mged_dozoom_view_reads_impl "${BRLCAD_SOURCE_DIR}/src/mged/dozoom.c")
  if(EXISTS "${_mged_dozoom_view_reads_impl}")
    file(READ "${_mged_dozoom_view_reads_impl}" _mged_dozoom_view_reads_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_perspective_from_bsg]]
	[[rt_view_eye_pos_from_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_pmat_from_bsg]]
	[[rt_view_pmat_set_bsg]])
      string(REGEX MATCH "${_token}" _mged_dozoom_view_reads_adapter_hit
	"${_mged_dozoom_view_reads_contents}")
      if(NOT _mged_dozoom_view_reads_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/dozoom.c must route perspective, eye, model-to-view, and projection-matrix access through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[if[ \t]*\([^\n;]*v->[ \t]*gv_perspective[ \t]*>=]]
	[[fastf_t[ \t\r\n]+persp[ \t\r\n]*=[^\n;]*v->[ \t]*gv_perspective]]
	[[if[ \t]*\([^\n;]*v->[ \t]*gv_perspective[ \t]*<]]
	[[v->[ \t]*gv_eye_pos]]
	[[dm_loadmatrix[^\n;]*v->[ \t]*gv_model2view]]
	[[MAT_COPY[ \t]*\([ \t]*(saved_pmat|perspective_mat)[^,]*,[^\n;]*v->[ \t]*gv_pmat]]
	[[MAT_COPY[ \t]*\([ \t]*v->[ \t]*gv_pmat]])
      string(REGEX MATCH "${_pat}" _mged_dozoom_view_reads_direct_hit
	"${_mged_dozoom_view_reads_contents}")
      if(_mged_dozoom_view_reads_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/dozoom.c reintroduced direct BSG perspective/eye/model-to-view/projection-matrix access: ${_mged_dozoom_view_reads_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_mged_buttons_view_save_impl "${BRLCAD_SOURCE_DIR}/src/mged/buttons.c")
  if(EXISTS "${_mged_buttons_view_save_impl}")
    file(READ "${_mged_buttons_view_save_impl}" _mged_buttons_view_save_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_rotation_from_bsg]]
	[[rt_view_center_from_bsg]]
	[[rt_view_scale_set_bsg]]
	[[rt_view_rotation_set_bsg]]
	[[rt_view_center_vec_set_bsg]])
      string(REGEX MATCH "${_token}" _mged_buttons_view_save_adapter_hit
	"${_mged_buttons_view_save_contents}")
      if(NOT _mged_buttons_view_save_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/buttons.c must route saved-view snapshots through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[sav_vscale[ \t]*=[^;\n]*gv_scale]]
	[[MAT_COPY[ \t\r\n]*\([ \t\r\n]*sav_viewrot[^;]*gv_rotation]]
	[[MAT_COPY[ \t\r\n]*\([ \t\r\n]*sav_toviewcenter[^;]*gv_center]]
	[[view_state->[ \t\r\n]*vs_gvp->[ \t\r\n]*gv_scale[ \t\r\n]*=[^\n;]*sav_vscale]]
	[[MAT_COPY[ \t\r\n]*\([^\n,]*view_state->[ \t\r\n]*vs_gvp->[ \t\r\n]*gv_rotation[^;]*sav_viewrot]]
	[[MAT_COPY[ \t\r\n]*\([^\n,]*view_state->[ \t\r\n]*vs_gvp->[ \t\r\n]*gv_center[^;]*sav_toviewcenter]])
      string(REGEX MATCH "${_pat}" _mged_buttons_view_save_direct_hit
	"${_mged_buttons_view_save_contents}")
      if(_mged_buttons_view_save_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/buttons.c reintroduced direct BSG saved-view reads: ${_mged_buttons_view_save_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_mged_rtif_view_reads_impl "${BRLCAD_SOURCE_DIR}/src/mged/rtif.c")
  if(EXISTS "${_mged_rtif_view_reads_impl}")
    file(READ "${_mged_rtif_view_reads_impl}" _mged_rtif_view_reads_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_set_bsg]]
	[[rt_view_rotation_set_bsg]]
	[[rt_view_center_vec_set_bsg]]
	[[rt_view_view2model_from_bsg]])
      string(REGEX MATCH "${_token}" _mged_rtif_view_reads_adapter_hit
	"${_mged_rtif_view_reads_contents}")
      if(NOT _mged_rtif_view_reads_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/rtif.c must route RT command view access through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[view_state->[ \t\r\n]*vs_gvp->[ \t\r\n]*gv_scale[ \t\r\n]*=]]
	[[MAT_COPY[ \t\r\n]*\([^\n,]*view_state->[ \t\r\n]*vs_gvp->[ \t\r\n]*gv_rotation]]
	[[MAT_IDN[ \t\r\n]*\([^\n;]*view_state->[ \t\r\n]*vs_gvp->[ \t\r\n]*gv_rotation]]
	[[MAT_DELTAS_VEC_NEG[ \t\r\n]*\([^;]*view_state->[ \t\r\n]*vs_gvp->[ \t\r\n]*gv_center]]
	[[MAT4X3PNT[^\n;]*gv_view2model]])
      string(REGEX MATCH "${_pat}" _mged_rtif_view_reads_direct_hit
	"${_mged_rtif_view_reads_contents}")
      if(_mged_rtif_view_reads_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/rtif.c reintroduced direct BSG RT command view access: ${_mged_rtif_view_reads_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_mged_edsol_abs_tran_impl "${BRLCAD_SOURCE_DIR}/src/mged/edsol.c")
  if(EXISTS "${_mged_edsol_abs_tran_impl}")
    file(READ "${_mged_edsol_abs_tran_impl}" _mged_edsol_abs_tran_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_view2model_from_bsg]])
      string(REGEX MATCH "${_token}" _mged_edsol_abs_tran_adapter_hit
	"${_mged_edsol_abs_tran_contents}")
      if(NOT _mged_edsol_abs_tran_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/edsol.c must route absolute translation bookkeeping reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[inv_Viewscale[ \t]*=[^\n;]*/[^\n;]*gv_scale]]
	[[MAT4X3PNT[ \t\r\n]*\([ \t\r\n]*model_pos[^;]*gv_view2model]]
	[[MAT4X3PNT[ \t\r\n]*\([ \t\r\n]*ea_view_pos[^;]*gv_model2view]]
	[[MAT4X3PNT[ \t\r\n]*\([ \t\r\n]*rot_point[^;]*gv_view2model]]
	[[MAT4X3VEC[ \t\r\n]*\([ \t\r\n]*view_dir[^;]*gv_view2model[^;]*(view_z_dir|z_dir)]]
	[[MAT4X3VEC[ \t\r\n]*\([ \t\r\n]*view_pl[^;]*gv_view2model[^;]*view_dir]]
	[[MAT4X3VEC[ \t\r\n]*\([ \t\r\n]*dir[^;]*gv_view2model[^;]*work]]
	[[rt_bot_find_[ve]_nearest_pt2[^;]*gv_model2view]]
	[[nmg_find_e_nearest_pt2[^;]*gv_model2view]]
	[[MAT4X3PNT[ \t\r\n]*\([ \t\r\n]*start_pt[^;]*gv_view2model]]
	[[MAT4X3VEC[ \t\r\n]*\([ \t\r\n]*dir[^;]*gv_view2model[^;]*tmp]]
	[[view_state->[ \t\r\n]*vs_gvp->[ \t\r\n]*gv_model2view]]
	[[view_state->[ \t\r\n]*vs_gvp->[ \t\r\n]*gv_view2model]])
      string(REGEX MATCH "${_pat}" _mged_edsol_abs_tran_direct_hit
	"${_mged_edsol_abs_tran_contents}")
      if(_mged_edsol_abs_tran_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/edsol.c reintroduced direct BSG edit view snapshot reads: ${_mged_edsol_abs_tran_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_align_view_reads_impl "${BRLCAD_SOURCE_DIR}/src/libged/view/align.c")
  if(EXISTS "${_libged_align_view_reads_impl}")
    file(READ "${_libged_align_view_reads_impl}" _libged_align_view_reads_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_center_from_bsg]]
	[[rt_view_view2model_from_bsg]]
	[[rt_view_aet_from_bsg]]
	[[rt_view_aet_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_align_view_reads_adapter_hit
	"${_libged_align_view_reads_contents}")
      if(NOT _libged_align_view_reads_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/align.c must route center/eye/AET reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[MAT4X3PNT[^\n;]*gv_center]]
	[[MAT4X3PNT[^\n;]*gv_view2model]]
	[[gv_aet[ \t]*\[]]
	[[VSET[^\n;]*gv_aet]]
	[[bsg_mat_aet[ \t\r\n]*\(]])
      string(REGEX MATCH "${_pat}" _libged_align_view_reads_direct_hit
	"${_libged_align_view_reads_contents}")
      if(_libged_align_view_reads_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/align.c reintroduced direct BSG align view reads: ${_libged_align_view_reads_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_vutil_view_reads_impl "${BRLCAD_SOURCE_DIR}/src/libged/vutil.c")
  if(EXISTS "${_libged_vutil_view_reads_impl}")
    file(READ "${_libged_vutil_view_reads_impl}" _libged_vutil_view_reads_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_rotation_from_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_view2model_from_bsg]]
	[[rt_view_center_from_bsg]]
	[[rt_view_inverse_size_from_bsg]]
	[[rt_view_keypoint_from_bsg]]
	[[rt_view_rotate_about_from_bsg]]
	[[rt_view_center_vec_set_bsg]]
	[[rt_view_rotation_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_vutil_view_reads_adapter_hit
	"${_libged_vutil_view_reads_contents}")
      if(NOT _libged_vutil_view_reads_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/vutil.c must route rotation/translation helper view reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[bn_mat_inv[ \t\r\n]*\([^;]*gv_rotation]]
	[[bn_mat_mul[ \t\r\n]*\([^;]*gv_rotation]]
	[[bn_mat_mul2[ \t\r\n]*\([^;]*gv_rotation]]
	[[MAT4X3PNT[^\n;]*gv_model2view]]
	[[MAT4X3PNT[^\n;]*gv_view2model]]
	[[MAT_DELTAS_GET_NEG[ \t\r\n]*\([^;]*gv_center]]
	[[MAT_DELTAS_VEC_NEG[ \t\r\n]*\([^;]*gv_center]]
	[[MAT_COPY[^\n;]*gv_rotation]]
	[[(^|[^A-Za-z0-9_])gv_isize([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_keypoint([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_rotate_about([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libged_vutil_view_reads_direct_hit
	"${_libged_vutil_view_reads_contents}")
      if(_libged_vutil_view_reads_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/vutil.c reintroduced direct BSG utility view reads: ${_libged_vutil_view_reads_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_isize_impl "${BRLCAD_SOURCE_DIR}/src/libged/isize/isize.c")
  if(EXISTS "${_libged_isize_impl}")
    file(READ "${_libged_isize_impl}" _libged_isize_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_inverse_size_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_isize_adapter_hit
	"${_libged_isize_contents}")
      if(NOT _libged_isize_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/isize/isize.c must route inverse view-size reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_isize([^A-Za-z0-9_]|$)]]
      _libged_isize_direct_hit "${_libged_isize_contents}")
    if(_libged_isize_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/isize/isize.c reintroduced direct BSG inverse view-size reads: ${_libged_isize_direct_hit}")
    endif()
  endif()

  set(_libged_keypoint_impl "${BRLCAD_SOURCE_DIR}/src/libged/keypoint/keypoint.c")
  if(EXISTS "${_libged_keypoint_impl}")
    file(READ "${_libged_keypoint_impl}" _libged_keypoint_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_keypoint_from_bsg]]
	[[rt_view_keypoint_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_keypoint_adapter_hit
	"${_libged_keypoint_contents}")
      if(NOT _libged_keypoint_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/keypoint/keypoint.c must route keypoint get/set through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_keypoint([^A-Za-z0-9_]|$)]]
      _libged_keypoint_direct_hit "${_libged_keypoint_contents}")
    if(_libged_keypoint_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/keypoint/keypoint.c reintroduced direct BSG keypoint access: ${_libged_keypoint_direct_hit}")
    endif()
  endif()

  set(_libged_rotate_about_impl "${BRLCAD_SOURCE_DIR}/src/libged/rot/rotate_about.c")
  if(EXISTS "${_libged_rotate_about_impl}")
    file(READ "${_libged_rotate_about_impl}" _libged_rotate_about_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_rotate_about_from_bsg]]
	[[rt_view_rotate_about_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_rotate_about_adapter_hit
	"${_libged_rotate_about_contents}")
      if(NOT _libged_rotate_about_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/rot/rotate_about.c must route rotate-about get/set through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[gedp->[ \t\r\n]*ged_gvp->[ \t\r\n]*gv_rotate_about]]
      _libged_rotate_about_direct_hit "${_libged_rotate_about_contents}")
    if(_libged_rotate_about_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/rot/rotate_about.c reintroduced direct BSG rotate-about access: ${_libged_rotate_about_direct_hit}")
    endif()
  endif()

  set(_libged_view_knob_impl "${BRLCAD_SOURCE_DIR}/src/libged/view/knob.c")
  if(EXISTS "${_libged_view_knob_impl}")
    file(READ "${_libged_view_knob_impl}" _libged_view_knob_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_keypoint_from_bsg]]
	[[rt_view_absolute_scale_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_view_knob_adapter_hit
	"${_libged_view_knob_contents}")
      if(NOT _libged_view_knob_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/knob.c must route keypoint-origin rotation reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[gv_coord[ \t]*==[ \t]*'m']]
      _libged_view_knob_coord_direct_hit "${_libged_view_knob_contents}")
    if(_libged_view_knob_coord_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/knob.c reintroduced direct BSG coord default reads: ${_libged_view_knob_coord_direct_hit}")
    endif()
    string(REGEX MATCH [[bsg_knobs_rot[^\n;]*gv_keypoint]]
      _libged_view_knob_direct_hit "${_libged_view_knob_contents}")
    if(_libged_view_knob_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/knob.c reintroduced direct BSG keypoint-origin rotation reads: ${_libged_view_knob_direct_hit}")
    endif()
    string(REGEX MATCH [[bu_vls_printf[^\n;]*gv_a_scale]]
      _libged_view_knob_a_scale_direct_hit "${_libged_view_knob_contents}")
    if(_libged_view_knob_a_scale_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/knob.c reintroduced direct BSG absolute-scale query output: ${_libged_view_knob_a_scale_direct_hit}")
    endif()
  endif()

  set(_libged_preview_view_reads_impl "${BRLCAD_SOURCE_DIR}/src/libged/draw/preview.cpp")
  if(EXISTS "${_libged_preview_view_reads_impl}")
    file(READ "${_libged_preview_view_reads_impl}" _libged_preview_view_reads_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_rotation_from_bsg]]
	[[rt_view_rotation_set_bsg]]
	[[rt_view_center_vec_set_bsg]]
	[[rt_view_view2model_from_bsg]])
      string(REGEX MATCH "${_token}" _libged_preview_view_reads_adapter_hit
	"${_libged_preview_view_reads_contents}")
      if(NOT _libged_preview_view_reads_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/draw/preview.cpp must route preview camera view access through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[MAT_COPY[ \t\r\n]*\([^\n,]*gedp->[ \t\r\n]*ged_gvp->[ \t\r\n]*gv_rotation]]
	[[MAT_DELTAS_VEC_NEG[ \t\r\n]*\([^\n,]*gedp->[ \t\r\n]*ged_gvp->[ \t\r\n]*gv_center]]
	[[MAT4X3PNT[^\n;]*gv_view2model]])
      string(REGEX MATCH "${_pat}" _libged_preview_view_reads_direct_hit
	"${_libged_preview_view_reads_contents}")
      if(_libged_preview_view_reads_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/draw/preview.cpp reintroduced direct BSG preview camera view access: ${_libged_preview_view_reads_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_loadview_impl "${BRLCAD_SOURCE_DIR}/src/libged/draw/loadview.cpp")
  if(EXISTS "${_libged_loadview_impl}")
    file(READ "${_libged_loadview_impl}" _libged_loadview_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_size_set_bsg]]
	[[rt_view_rotation_set_bsg]]
	[[rt_view_center_vec_set_bsg]]
	[[rt_view_perspective_set_bsg]])
      string(REGEX MATCH "${_token}" _libged_loadview_adapter_hit
	"${_libged_loadview_contents}")
      if(NOT _libged_loadview_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/draw/loadview.cpp must route view restore writes through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])bsg_view_set_size([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_set_rotation([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_set_center_vec([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_set_perspective([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libged_loadview_direct_hit
	"${_libged_loadview_contents}")
      if(_libged_loadview_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/draw/loadview.cpp reintroduced direct BSG view restore writes: ${_libged_loadview_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libtclcad_view_draw_impl "${BRLCAD_SOURCE_DIR}/src/libtclcad/view/draw.c")
  if(EXISTS "${_libtclcad_view_draw_impl}")
    file(READ "${_libtclcad_view_draw_impl}" _libtclcad_view_draw_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_perspective_from_bsg]]
	[[rt_view_pmat_from_bsg]])
      string(REGEX MATCH "${_token}" _libtclcad_view_draw_adapter_hit
	"${_libtclcad_view_draw_contents}")
      if(NOT _libtclcad_view_draw_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/view/draw.c must route draw matrix/perspective reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_perspective([^A-Za-z0-9_]|$)]]
	[[dm_loadpmatrix[^\n;]*gv_pmat]])
      string(REGEX MATCH "${_pat}" _libtclcad_view_draw_direct_hit
	"${_libtclcad_view_draw_contents}")
      if(_libtclcad_view_draw_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/view/draw.c reintroduced direct BSG draw matrix/perspective reads: ${_libtclcad_view_draw_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libtclcad_view_refresh_impl "${BRLCAD_SOURCE_DIR}/src/libtclcad/view/refresh.c")
  if(EXISTS "${_libtclcad_view_refresh_impl}")
    file(READ "${_libtclcad_view_refresh_impl}" _libtclcad_view_refresh_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_local2base_from_bsg]]
	[[rt_view_base2local_from_bsg]])
      string(REGEX MATCH "${_token}" _libtclcad_view_refresh_adapter_hit
	"${_libtclcad_view_refresh_contents}")
      if(NOT _libtclcad_view_refresh_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/view/refresh.c must route unit-conversion stash reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[double[^\n;]*=[^\n;]*gdvp->[ \t\r\n]*gv_(local2base|base2local)]]
      _libtclcad_view_refresh_direct_hit "${_libtclcad_view_refresh_contents}")
    if(_libtclcad_view_refresh_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libtclcad/view/refresh.c reintroduced direct BSG unit-conversion stash reads: ${_libtclcad_view_refresh_direct_hit}")
    endif()
  endif()

  set(_libtclcad_wrapper_impl "${BRLCAD_SOURCE_DIR}/src/libtclcad/wrapper.c")
  if(EXISTS "${_libtclcad_wrapper_impl}")
    file(READ "${_libtclcad_wrapper_impl}" _libtclcad_wrapper_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_dimensions_set_bsg]]
	[[rt_view_perspective_from_bsg]]
	[[rt_view_lod_policy_from_bsg]])
      string(REGEX MATCH "${_token}" _libtclcad_wrapper_adapter_hit
	"${_libtclcad_wrapper_contents}")
      if(NOT _libtclcad_wrapper_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/wrapper.c must route perspective, dimensions, and zoom-refresh LoD policy access through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_perspective([^A-Za-z0-9_]|$)]]
	[[gedp->[ \t\r\n]*ged_gvp->[ \t\r\n]*gv_(width|height)[ \t\r\n]*=]]
	[[gdvp->[ \t\r\n]*gv_(width|height)[ \t\r\n]*=]]
	[[struct[ \t\r\n]+bsg_lod_source_policy_settings]]
	[[bsg_view_lod_source_policy_get]])
      string(REGEX MATCH "${_pat}" _libtclcad_wrapper_direct_hit
	"${_libtclcad_wrapper_contents}")
      if(_libtclcad_wrapper_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/wrapper.c reintroduced direct BSG perspective/LoD policy/dimension access: ${_libtclcad_wrapper_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libtclcad_view_axes_impl "${BRLCAD_SOURCE_DIR}/src/libtclcad/view/axes.c")
  if(EXISTS "${_libtclcad_view_axes_impl}")
    file(READ "${_libtclcad_view_axes_impl}" _libtclcad_view_axes_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_info_from_bsg]])
      string(REGEX MATCH "${_token}" _libtclcad_view_axes_adapter_hit
	"${_libtclcad_view_axes_contents}")
      if(NOT _libtclcad_view_axes_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/view/axes.c must route axes size scaling through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])gv_size([^A-Za-z0-9_]|$)]]
      _libtclcad_view_axes_direct_hit "${_libtclcad_view_axes_contents}")
    if(_libtclcad_view_axes_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libtclcad/view/axes.c reintroduced direct BSG axes size scaling reads: ${_libtclcad_view_axes_direct_hit}")
    endif()
  endif()

  set(_libtclcad_commands_impl "${BRLCAD_SOURCE_DIR}/src/libtclcad/commands.c")
  if(EXISTS "${_libtclcad_commands_impl}")
    file(READ "${_libtclcad_commands_impl}" _libtclcad_commands_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_lod_policy_from_bsg]]
	[[rt_view_lod_policy_apply_bsg]]
	[[rt_view_dimensions_set_bsg]]
	[[rt_view_screen_to_view_from_bsg]]
	[[rt_view_center_from_bsg]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_view2model_from_bsg]]
	[[rt_view_snap_lines_from_bsg]]
	[[rt_view_prepare_tcl_snap_bsg]]
	[[rt_view_center_linesnap_bsg]]
	[[rt_view_zclip_from_bsg]]
	[[rt_view_zclip_set_bsg]])
      string(REGEX MATCH "${_token}" _libtclcad_commands_adapter_hit
	"${_libtclcad_commands_contents}")
      if(NOT _libtclcad_commands_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/commands.c must route TclCAD LoD policy and command view/dimension access through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_token
	[[rt_view_coord_from_bsg]]
	[[rt_view_coord_set_bsg]])
      string(REGEX MATCH "${_token}" _libtclcad_commands_coord_adapter_hit
	"${_libtclcad_commands_contents}")
      if(NOT _libtclcad_commands_coord_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/commands.c must route set_coord get/set through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[struct[ \t\r\n]+bsg_lod_source_policy_settings]]
	[[bsg_view_lod_source_policy_get]]
	[[bsg_view_lod_source_policy_set]]
	[[gedp->[ \t\r\n]*ged_gvp->[ \t\r\n]*gv_(width|height)[ \t\r\n]*=]]
	[[gdvp->[ \t\r\n]*gv_(width|height)[ \t\r\n]*=]]
	[[(^|[^A-Za-z0-9_])bsg_view_snap_lines([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_prepare_tcl_snap([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_center_linesnap([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_(set_)?zclip([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_screen_to_view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_center([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]]
	[[gdvp->[ \t\r\n]*gv_coord[ \t\r\n]*=]]
	[[bu_vls_printf[^\n;]*gv_coord]]
	[[MAT4X3PNT[^\n;]*gdvp->gv_model2view]]
	[[MAT4X3PNT[^\n;]*gdvp->gv_view2model]])
      string(REGEX MATCH "${_pat}" _libtclcad_commands_bsg_policy_hit
	"${_libtclcad_commands_contents}")
      if(_libtclcad_commands_bsg_policy_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/commands.c reintroduced direct BSG LoD policy/view reads or dimension writes: ${_libtclcad_commands_bsg_policy_hit}")
      endif()
    endforeach()
  endif()

  set(_libtclcad_mouse_first_impl "${BRLCAD_SOURCE_DIR}/src/libtclcad/mouse.c")
  if(EXISTS "${_libtclcad_mouse_first_impl}")
    file(READ "${_libtclcad_mouse_first_impl}" _libtclcad_mouse_first_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_size_from_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_view2model_from_bsg]]
	[[rt_view_rotation_from_bsg]]
	[[rt_view_dimensions_set_bsg]]
	[[rt_view_screen_to_view_from_bsg]]
	[[rt_view_prepare_tcl_snap_bsg]])
      string(REGEX MATCH "${_token}" _libtclcad_mouse_first_adapter_hit
	"${_libtclcad_mouse_first_contents}")
      if(NOT _libtclcad_mouse_first_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/mouse.c must route simple selection/point input view reads and dimension writes through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_size([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_rotation([^A-Za-z0-9_]|$)]]
	[[gdvp->[ \t\r\n]*gv_(width|height)[ \t\r\n]*=]]
	[[(^|[^A-Za-z0-9_])bsg_screen_to_view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_prepare_tcl_snap([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libtclcad_mouse_first_direct_hit
	"${_libtclcad_mouse_first_contents}")
      if(_libtclcad_mouse_first_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/mouse.c reintroduced direct BSG mouse view snapshot reads or dimension writes: ${_libtclcad_mouse_first_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libtclcad_polygons_impl "${BRLCAD_SOURCE_DIR}/src/libtclcad/polygons.c")
  if(EXISTS "${_libtclcad_polygons_impl}")
    file(READ "${_libtclcad_polygons_impl}" _libtclcad_polygons_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_from_bsg]]
	[[rt_view_center_from_bsg]]
	[[rt_view_rotation_from_bsg]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_view2model_from_bsg]]
	[[rt_view_plane_from_bsg]]
	[[rt_view_dimensions_set_bsg]]
	[[rt_view_screen_to_view_from_bsg]]
	[[rt_view_prepare_tcl_snap_bsg]])
      string(REGEX MATCH "${_token}" _libtclcad_polygons_adapter_hit
	"${_libtclcad_polygons_contents}")
      if(NOT _libtclcad_polygons_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/polygons.c must route polygon view snapshots and dimension writes through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_center([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_rotation([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_view2model([^A-Za-z0-9_]|$)]]
	[[gdvp->[ \t\r\n]*gv_(width|height)[ \t\r\n]*=]]
	[[(^|[^A-Za-z0-9_])bsg_screen_to_view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_plane([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_prepare_tcl_snap([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libtclcad_polygons_direct_hit
	"${_libtclcad_polygons_contents}")
      if(_libtclcad_polygons_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/polygons.c reintroduced direct BSG polygon view snapshot reads or dimension writes: ${_libtclcad_polygons_direct_hit}")
      endif()
    endforeach()
  endif()

  file(GLOB_RECURSE _mged_rt_numeric_files LIST_DIRECTORIES false
    "${BRLCAD_SOURCE_DIR}/src/mged/*.c"
    "${BRLCAD_SOURCE_DIR}/src/mged/*.cc"
    "${BRLCAD_SOURCE_DIR}/src/mged/*.cpp"
    "${BRLCAD_SOURCE_DIR}/src/mged/*.cxx"
    "${BRLCAD_SOURCE_DIR}/src/mged/*.h"
    "${BRLCAD_SOURCE_DIR}/src/mged/*.hh"
    "${BRLCAD_SOURCE_DIR}/src/mged/*.hpp"
  )
  foreach(_file IN LISTS _mged_rt_numeric_files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    file(READ "${_file}" _mged_rt_numeric_contents)
    foreach(_pat
	[[(^|[^A-Za-z0-9_])BSG_VIEW_(MIN|MAX|RANGE)([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])BSG_INV_(VIEW|4096)([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])INV_BV([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])BSG_MINVIEWSIZE([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])BSG_MINVIEWSCALE([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _mged_rt_numeric_hit
	"${_mged_rt_numeric_contents}")
      if(_mged_rt_numeric_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must use RT-owned view numeric policy instead of BSG numeric aliases: ${_mged_rt_numeric_hit}")
      endif()
    endforeach()
  endforeach()

  set(_librt_edit_tor_test
    "${BRLCAD_SOURCE_DIR}/src/librt/tests/edit/tor.cpp")
  if(EXISTS "${_librt_edit_tor_test}")
    file(READ "${_librt_edit_tor_test}" _librt_edit_tor_test_contents)
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])BSG_VIEW_(MIN|MAX|RANGE)([^A-Za-z0-9_]|$)]]
      _librt_edit_tor_bsg_view_hit "${_librt_edit_tor_test_contents}")
    if(_librt_edit_tor_bsg_view_hit)
      _brlobol_pivot_guard_fail(
	"src/librt/tests/edit/tor.cpp must describe edit view units with RT_VIEW_* names")
    endif()
  endif()

  set(_librt_tests_cmake
    "${BRLCAD_SOURCE_DIR}/src/librt/tests/CMakeLists.txt")
  if(EXISTS "${_librt_tests_cmake}")
    file(READ "${_librt_tests_cmake}" _librt_tests_cmake_contents)
    foreach(_token
	"brlcad_addexec(rt_view_info"
	"brlcad_add_test(NAME rt_view_info"
	"brlcad_addexec(rt_view_legacy_bsg"
	"brlcad_add_test(NAME rt_view_legacy_bsg")
      string(FIND "${_librt_tests_cmake_contents}" "${_token}"
	_librt_view_info_cmake_token_idx)
      if(_librt_view_info_cmake_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/librt/tests/CMakeLists.txt must register RT view-info test ${_token}")
      endif()
    endforeach()
  endif()

  file(GLOB_RECURSE _primitive_files LIST_DIRECTORIES false
    "${BRLCAD_SOURCE_DIR}/src/librt/primitives/*.c"
    "${BRLCAD_SOURCE_DIR}/src/librt/primitives/*.cpp"
    "${BRLCAD_SOURCE_DIR}/src/librt/primitives/*.h"
  )
  foreach(_file IN LISTS _primitive_files)
    file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
    if("${_rel}" STREQUAL "src/librt/primitives/sketch/polygons.c")
      continue()
    endif()
    if("${_rel}" STREQUAL "src/librt/primitives/sketch/polygons_legacy_bsg.c")
      continue()
    endif()

    file(READ "${_file}" _contents)
    set(_primitive_view_forbidden
      [[#[ \t]*include[ \t]*[<"]bsg/view_state\.h]]
      [[const[ \t\r\n]+struct[ \t\r\n]+bsg_view[ \t\r\n]*\*]]
      [[(^|[^A-Za-z0-9_])primitive_lod_curve_scale([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])primitive_lod_bot_threshold([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])view_avg_sample_spacing([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])solid_point_spacing([^A-Za-z0-9_]|$)]]
    )
    foreach(_pat IN LISTS _primitive_view_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced BSG view state into RT primitive callbacks: ${_hit}")
      endif()
    endforeach()
  endforeach()

  set(_ged_bsg_draw_source
    "${BRLCAD_SOURCE_DIR}/src/libged/bsg_ged_draw_source.c")
  if(EXISTS "${_ged_bsg_draw_source}")
    file(READ "${_ged_bsg_draw_source}" _ged_bsg_draw_source_contents)
    foreach(_token
	"const struct rt_view_info *view_info"
	"ft_lod_realize(&realization, ip, tol, view_info"
	"ft_indexed_face_set(&face_set, ip, ttol, tol")
      string(FIND "${_ged_bsg_draw_source_contents}" "${_token}"
	_ged_view_info_token_idx)
      if(_ged_view_info_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/libged/bsg_ged_draw_source.c must pass RT view-info snapshots into RT primitive callbacks via ${_token}")
      endif()
    endforeach()
    set(_ged_bsg_draw_source_forbidden
      [[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
      [[(^|[^A-Za-z0-9_])rt_view_info_from_bsg([^A-Za-z0-9_]|$)]]
    )
    foreach(_pat IN LISTS _ged_bsg_draw_source_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_ged_bsg_draw_source_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/bsg_ged_draw_source.c reintroduced BSG-to-RT view adaptation below the publish boundary: ${_hit}")
      endif()
    endforeach()
    string(REGEX MATCH
      "ft_(lod_realize|indexed_face_set)[^;]*,[ \t\r\n]*v[ \t\r\n]*[),]"
      _ged_direct_bsg_view_hit "${_ged_bsg_draw_source_contents}")
    if(_ged_direct_bsg_view_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/bsg_ged_draw_source.c must not pass bsg_view directly to RT primitive callbacks: ${_ged_direct_bsg_view_hit}")
    endif()
  endif()

  set(_ged_bsg_draw_view
    "${BRLCAD_SOURCE_DIR}/src/libged/bsg_ged_draw_view.c")
  if(NOT EXISTS "${_ged_bsg_draw_view}")
    _brlobol_pivot_guard_fail(
      "src/libged/bsg_ged_draw_view.c is required to isolate GED BSG-to-RT view snapshot adaptation")
  else()
    file(READ "${_ged_bsg_draw_view}" _ged_bsg_draw_view_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	ged_draw_view_info_from_bsg
	ged_draw_view_perspective_from_bsg
	ged_draw_view_scale_from_bsg
	ged_draw_view_lod_policy_from_bsg
	ged_draw_view_lod_policy_apply_bsg
	ged_draw_view_lod_policy_apply_bsg_bot_threshold
	ged_draw_view_set_lod_bounds_update
	ged_draw_view_has_lod_bounds_update
	ged_draw_scene_ref_realization_set_bsg_view_policy
	ged_draw_mesh_lod_load_view_scene_ref
	ged_draw_mesh_lod_free_scene_ref
	rt_view_lod_bounds_callback_set_bsg
	rt_view_lod_bounds_callback_is_bsg
	rt_view_perspective_from_bsg
	rt_view_scale_from_bsg
	rt_mesh_lod_load_view_scene_ref_bsg
	rt_mesh_lod_free_scene_ref_bsg
	rt_view_info_from_bsg
	rt_view_lod_policy_from_bsg
	rt_view_lod_policy_apply_bsg)
      string(REGEX MATCH "${_token}" _ged_bsg_draw_view_token_hit
	"${_ged_bsg_draw_view_contents}")
      if(NOT _ged_bsg_draw_view_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/bsg_ged_draw_view.c must own GED BSG-to-RT view adapter token ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[struct[ \t\r\n]+bsg_lod_source_policy_settings]]
	[[#[ \t]*include[ \t]*[<"]bsg/lod\.h]]
	[[bsg_view_lod_source_policy_get]]
	[[bsg_view_lod_source_policy_set]]
	[[bsg_view_bounds]]
	[[bsg_mesh_lod_load_view_scene_ref]]
	[[bsg_mesh_lod_free_scene_ref]]
	[[rt_mesh_lod_bsg]]
	[[->[ \t\r\n]*gv_scale]]
	[[->[ \t\r\n]*gv_perspective]])
      string(REGEX MATCH "${_pat}" _ged_bsg_draw_view_bsg_policy_hit
	"${_ged_bsg_draw_view_contents}")
      if(_ged_bsg_draw_view_bsg_policy_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/bsg_ged_draw_view.c must delegate BSG LoD view hooks through rt/view_legacy_bsg.h: ${_ged_bsg_draw_view_bsg_policy_hit}")
      endif()
    endforeach()
  endif()

  set(_ged_bsg_draw_private
    "${BRLCAD_SOURCE_DIR}/src/libged/bsg_ged_draw_private.h")
  if(EXISTS "${_ged_bsg_draw_private}")
    file(READ "${_ged_bsg_draw_private}" _ged_bsg_draw_private_contents)
    string(FIND "${_ged_bsg_draw_private_contents}" "ged_draw_view_info_from_bsg"
      _ged_private_view_adapter_idx)
    if(_ged_private_view_adapter_idx EQUAL -1)
      _brlobol_pivot_guard_fail(
	"src/libged/bsg_ged_draw_private.h must expose ged_draw_view_info_from_bsg for BSG-owning wrapper files")
    endif()
    string(FIND "${_ged_bsg_draw_private_contents}" "ged_draw_view_lod_policy"
      _ged_private_lod_policy_idx)
    if(_ged_private_lod_policy_idx EQUAL -1)
      _brlobol_pivot_guard_fail(
	"src/libged/bsg_ged_draw_private.h must expose the private GED view LoD policy snapshot")
    endif()
    foreach(_func
	ged_draw_view_perspective_from_bsg
	ged_draw_view_scale_from_bsg
	ged_draw_view_lod_policy_from_bsg
	ged_draw_view_lod_policy_apply_bsg
	ged_draw_view_lod_policy_apply_bsg_bot_threshold
	ged_draw_view_set_lod_bounds_update
	ged_draw_view_has_lod_bounds_update
	ged_draw_scene_ref_realization_set_bsg_view_policy
	ged_draw_mesh_lod_load_view_scene_ref
	ged_draw_mesh_lod_free_scene_ref
	ged_draw_brep_mesh_lod_detail_setup)
      string(FIND "${_ged_bsg_draw_private_contents}" "${_func}"
	_ged_private_lod_adapter_idx)
      if(_ged_private_lod_adapter_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/libged/bsg_ged_draw_private.h must expose ${_func} for BSG-owning wrapper files")
      endif()
    endforeach()
    string(FIND "${_ged_bsg_draw_private_contents}" "struct rt_mesh_lod *rt_mesh_lod"
      _ged_private_lod_owner_idx)
    if(_ged_private_lod_owner_idx EQUAL -1)
      _brlobol_pivot_guard_fail(
	"src/libged/bsg_ged_draw_private.h must keep RT-owned mesh-LoD state visible through struct rt_mesh_lod *rt_mesh_lod")
    endif()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])ged_draw_rt_mesh_lod_bsg([^A-Za-z0-9_]|$)]]
	[[struct[ \t\r\n]+bsg_mesh_lod]]
	[[struct[ \t\r\n]+bsg_mesh_lod[ \t\r\n]*[*][ \t\r\n]*mesh_lod]])
      string(REGEX MATCH "${_pat}" _ged_private_bsg_lod_state_hit
	"${_ged_bsg_draw_private_contents}")
      if(_ged_private_bsg_lod_state_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/bsg_ged_draw_private.h must not expose borrowed BSG mesh-LoD state: ${_ged_private_bsg_lod_state_hit}")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])mesh_c([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_mesh_lod_context([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _ged_private_draw_mesh_c_hit
	"${_ged_bsg_draw_private_contents}")
      if(_ged_private_draw_mesh_c_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/bsg_ged_draw_private.h must not carry BSG mesh-LoD context through draw source state: ${_ged_private_draw_mesh_c_hit}")
      endif()
    endforeach()
    foreach(_func
	ged_draw_scene_ref_publish_primitive_face_set
	ged_draw_scene_ref_publish_primitive_wireframe)
      string(REGEX MATCH
	"${_func}[^;]*const[ \t\r\n]+struct[ \t\r\n]+rt_view_info[ \t\r\n]*\\*[ \t\r\n]*view_info"
	_ged_private_view_info_hit "${_ged_bsg_draw_private_contents}")
      if(NOT _ged_private_view_info_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/bsg_ged_draw_private.h must route ${_func} through rt_view_info")
      endif()
    endforeach()
    string(REGEX MATCH
      "ged_draw_scene_ref_publish_primitive_face_set[^;]*const[ \t\r\n]+struct[ \t\r\n]+bsg_view[ \t\r\n]*\\*"
      _ged_private_face_set_bsg_view_hit "${_ged_bsg_draw_private_contents}")
    if(_ged_private_face_set_bsg_view_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/bsg_ged_draw_private.h reintroduced bsg_view for primitive face-set publication: ${_ged_private_face_set_bsg_view_hit}")
    endif()
  endif()

  set(_ged_bsg_draw_brep
    "${BRLCAD_SOURCE_DIR}/src/libged/bsg_ged_draw_brep.cpp")
  if(EXISTS "${_ged_bsg_draw_brep}")
    file(READ "${_ged_bsg_draw_brep}" _ged_bsg_draw_brep_contents)
    foreach(_token
	ged_draw_brep_mesh_lod_detail_setup
	rt_mesh_lod_detail_callbacks_set
	struct rt_mesh_lod_detail
	_ged_draw_brep_mesh_info_clbk)
      string(REGEX MATCH "${_token}" _ged_brep_lod_token_hit
	"${_ged_bsg_draw_brep_contents}")
      if(NOT _ged_brep_lod_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/bsg_ged_draw_brep.cpp must own BRep mesh-LoD detail callback token ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg/lod\.h]]
	[[struct[ \t\r\n]+bsg_mesh_lod]]
	[[bsg_mesh_lod_detail_setup_clbk]]
	[[bsg_mesh_lod_detail_clear_clbk]]
	[[bsg_mesh_lod_detail_free_clbk]]
	[[rt_mesh_lod_bsg]])
      string(REGEX MATCH "${_pat}" _ged_brep_bsg_lod_hit
	"${_ged_bsg_draw_brep_contents}")
      if(_ged_brep_bsg_lod_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/bsg_ged_draw_brep.cpp must register BRep full-detail callbacks through RT mesh-LoD detail helpers: ${_ged_brep_bsg_lod_hit}")
      endif()
    endforeach()
  endif()

  set(_ged_private_header "${BRLCAD_SOURCE_DIR}/src/libged/ged_private.h")
  if(EXISTS "${_ged_private_header}")
    file(READ "${_ged_private_header}" _ged_private_header_contents)
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg/lod\.h]]
	[[(^|[^A-Za-z0-9_])mesh_c([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_mesh_lod_context([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _ged_private_header_mesh_c_hit
	"${_ged_private_header_contents}")
      if(_ged_private_header_mesh_c_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/ged_private.h must not carry BSG mesh-LoD context or include bsg/lod.h for draw_data_t: ${_ged_private_header_mesh_c_hit}")
      endif()
    endforeach()
  endif()

  set(_ged_bsg_draw_cmake "${BRLCAD_SOURCE_DIR}/src/libged/CMakeLists.txt")
  if(EXISTS "${_ged_bsg_draw_cmake}")
    file(READ "${_ged_bsg_draw_cmake}" _ged_bsg_draw_cmake_contents)
    string(FIND "${_ged_bsg_draw_cmake_contents}" "bsg_ged_draw_view.c"
      _ged_bsg_draw_view_cmake_idx)
    if(_ged_bsg_draw_view_cmake_idx EQUAL -1)
      _brlobol_pivot_guard_fail(
	"src/libged/CMakeLists.txt must build src/libged/bsg_ged_draw_view.c")
    endif()
  endif()

  foreach(_rel
      include/ged/defines.h
      include/ged/bsg_ged_draw.h)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg/lod\.h]]
	[[(^|[^A-Za-z0-9_])bsg_mesh_lod_context([^A-Za-z0-9_]|$)]]
	[[struct[ \t\r\n]+bsg_mesh_lod_context[ \t\r\n]*[*][ \t\r\n]*ged_lod]])
      string(REGEX MATCH "${_pat}" _ged_public_lod_context_hit
	"${_contents}")
      if(_ged_public_lod_context_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must not expose GED-owned BSG mesh-LoD context state: ${_ged_public_lod_context_hit}")
      endif()
    endforeach()
  endforeach()

  foreach(_rel
      src/libged/ged.cpp
      src/libged/open/open.cpp
      src/libged/close/close.cpp)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg/lod\.h]]
	[[(^|[^A-Za-z0-9_])bsg_mesh_lod_context_create([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_mesh_lod_context_destroy([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])ged_lod([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _ged_lifecycle_lod_context_hit
	"${_contents}")
      if(_ged_lifecycle_lod_context_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must not create, destroy, or store a GED-owned BSG mesh-LoD context: ${_ged_lifecycle_lod_context_hit}")
      endif()
    endforeach()
  endforeach()

  foreach(_rel
      src/libged/bsg_ged_draw_draft.c
      src/libged/bsg_ged_draw_refs.c
      src/libged/bsg_ged_draw_transactions.c
      src/libged/draw.cpp
      src/libged/draw/draw.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    string(REGEX MATCH
      [[#[ \t]*include[ \t]*[<"]([.]/|[.][.]/)?bsg_ged_draw_private[.]h]]
      _ged_wrapper_private_include_hit "${_contents}")
    if(NOT _ged_wrapper_private_include_hit)
      _brlobol_pivot_guard_fail(
	"${_rel} must include bsg_ged_draw_private.h for GED BSG-to-RT view adapter declarations")
    endif()
    string(REGEX MATCH
      [[extern[ \t\r\n]+"C"[ \t\r\n]+void[ \t\r\n]+ged_draw_view_info_from_bsg]]
      _ged_wrapper_ad_hoc_cxx_adapter_hit "${_contents}")
    if(_ged_wrapper_ad_hoc_cxx_adapter_hit)
      _brlobol_pivot_guard_fail(
	"${_rel} must use bsg_ged_draw_private.h instead of redeclaring ged_draw_view_info_from_bsg: ${_ged_wrapper_ad_hoc_cxx_adapter_hit}")
    endif()
    string(FIND "${_contents}" "ged_draw_view_info_from_bsg" _ged_wrapper_adapter_idx)
    if(_ged_wrapper_adapter_idx EQUAL -1)
      _brlobol_pivot_guard_fail(
	"${_rel} must use ged_draw_view_info_from_bsg at the BSG-owning wrapper boundary")
    endif()
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[(^|[^A-Za-z0-9_])rt_view_info_from_bsg([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must not call the RT legacy BSG adapter directly: ${_hit}")
      endif()
    endforeach()
  endforeach()

  foreach(_rel
      src/libged/bsg_ged_draw_transactions.c
      src/libged/draw.cpp
      src/libged/draw/draw.c
      src/libged/view/gobjs.cpp)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])ged_lod([^A-Za-z0-9_]|$)]]
      _ged_draw_ged_lod_hit "${_contents}")
    if(_ged_draw_ged_lod_hit)
      _brlobol_pivot_guard_fail(
	"${_rel} must not thread ged_lod/BSG mesh-LoD context into draw source state: ${_ged_draw_ged_lod_hit}")
    endif()
  endforeach()

  set(_ged_view_independent_test "${BRLCAD_SOURCE_DIR}/src/libged/tests/draw/view_independent.cpp")
  if(EXISTS "${_ged_view_independent_test}")
    file(READ "${_ged_view_independent_test}" _ged_view_independent_test_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_scale_state_set_bsg]])
      string(REGEX MATCH "${_token}" _ged_view_independent_test_token_hit
	"${_ged_view_independent_test_contents}")
      if(NOT _ged_view_independent_test_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/tests/draw/view_independent.cpp must route oversized view setup through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[=[views\[[^]]+\]->[ \t\r\n]*gv_(scale|size|isize)[ \t\r\n]*=]=])
      string(REGEX MATCH "${_pat}" _ged_view_independent_test_direct_hit
	"${_ged_view_independent_test_contents}")
      if(_ged_view_independent_test_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/tests/draw/view_independent.cpp reintroduced direct BSG oversized view setup writes: ${_ged_view_independent_test_direct_hit}")
      endif()
    endforeach()
  endif()

  foreach(_rel
      src/libged/tests/draw/ged_draw_scene.cpp
      src/libged/tests/draw/lod.cpp
      src/libged/tests/draw/lod_crossrun.cpp
      src/libged/tests/draw/quad.cpp
      src/libged/tests/draw/util.cpp
      src/libged/tests/draw/view_command.cpp
      src/libged/tests/draw/view_independent.cpp)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg/lod\.h]]
	[[(^|[^A-Za-z0-9_])ged_lod([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_mesh_lod_key_get([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_mesh_lod_key_put([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_mesh_lod_clear_cache([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_mesh_lod_context_destroy([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_bounds([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _ged_draw_test_lod_context_hit
	"${_contents}")
      if(_ged_draw_test_lod_context_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must use RT mesh-LoD helpers instead of direct BSG LoD cache/context APIs: ${_ged_draw_test_lod_context_hit}")
      endif()
    endforeach()
  endforeach()

  set(_libged_selection_semantics_test "${BRLCAD_SOURCE_DIR}/src/libged/tests/draw/selection_semantics.cpp")
  if(EXISTS "${_libged_selection_semantics_test}")
    file(READ "${_libged_selection_semantics_test}" _selection_semantics_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_width_from_bsg]]
	[[rt_view_height_from_bsg]])
      string(REGEX MATCH "${_token}" _selection_semantics_adapter_hit
	"${_selection_semantics_contents}")
      if(NOT _selection_semantics_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/tests/draw/selection_semantics.cpp must route pick-point size reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[v->[ \t\r\n]*gv_width]]
	[[v->[ \t\r\n]*gv_height]])
      string(REGEX MATCH "${_pat}" _selection_semantics_direct_hit
	"${_selection_semantics_contents}")
      if(_selection_semantics_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/tests/draw/selection_semantics.cpp reintroduced direct pick-point BSG size reads: ${_selection_semantics_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libged_rtwizard_bsg_test "${BRLCAD_SOURCE_DIR}/src/libged/tests/draw/rtwizard_bsg.cpp")
  if(EXISTS "${_libged_rtwizard_bsg_test}")
    file(READ "${_libged_rtwizard_bsg_test}" _rtwizard_bsg_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_model2view_from_bsg]])
      string(REGEX MATCH "${_token}" _rtwizard_bsg_adapter_hit
	"${_rtwizard_bsg_contents}")
      if(NOT _rtwizard_bsg_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/tests/draw/rtwizard_bsg.cpp must route view matrix comparisons through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    string(REGEX MATCH [[v->[ \t\r\n]*gv_model2view]]
      _rtwizard_bsg_direct_hit "${_rtwizard_bsg_contents}")
    if(_rtwizard_bsg_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/tests/draw/rtwizard_bsg.cpp reintroduced direct view matrix comparison reads: ${_rtwizard_bsg_direct_hit}")
    endif()
  endif()

  foreach(_rel
      src/libged/bsg_ged_draw_binding.c
      src/libged/bsg_ged_draw_draft.c
      src/libged/bsg_ged_draw_refs.c
      src/libged/bsg_ged_draw_highlight.c
      src/libged/bsg_ged_draw_material.c
      src/libged/bsg_ged_draw_scene_root.c
      src/libged/bsg_ged_draw_state.c
      src/libged/bsg_ged_draw_transactions.c
      src/libged/bsg_ged_draw_tree.c
      src/libged/bsg_ged_draw_util.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    string(REGEX MATCH [[#[ \t]*include[ \t]*[<"]bsg/lod\.h]]
      _ged_bsg_draw_lod_header_hit "${_contents}")
    if(_ged_bsg_draw_lod_header_hit)
      _brlobol_pivot_guard_fail(
	"${_rel} must not include bsg/lod.h unless it owns transitional mesh-LoD payload operations: ${_ged_bsg_draw_lod_header_hit}")
    endif()
  endforeach()

  set(_ged_bsg_draw_binding
    "${BRLCAD_SOURCE_DIR}/src/libged/bsg_ged_draw_binding.c")
  if(EXISTS "${_ged_bsg_draw_binding}")
    file(READ "${_ged_bsg_draw_binding}" _ged_bsg_draw_binding_contents)
    string(FIND "${_ged_bsg_draw_binding_contents}" "rt_mesh_lod_destroy"
      _ged_binding_rt_lod_destroy_idx)
    if(_ged_binding_rt_lod_destroy_idx EQUAL -1)
      _brlobol_pivot_guard_fail(
	"src/libged/bsg_ged_draw_binding.c must destroy draw source mesh-LoD handles through rt_mesh_lod_destroy")
    endif()
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg/lod\.h]]
	[[bsg_mesh_lod_destroy]]
	[[->[ \t\r\n]*mesh_lod]])
      string(REGEX MATCH "${_pat}" _ged_binding_bsg_lod_destroy_hit
	"${_ged_bsg_draw_binding_contents}")
      if(_ged_binding_bsg_lod_destroy_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/bsg_ged_draw_binding.c must not destroy or store borrowed BSG mesh-LoD state: ${_ged_binding_bsg_lod_destroy_hit}")
      endif()
    endforeach()
  endif()

  foreach(_rel
      src/libged/bsg_ged_draw_draft.c
      src/libged/bsg_ged_draw_refs.c
      src/libged/bsg_ged_draw_transactions.c
      src/libged/draw.cpp
      src/libged/draw/draw.c
      src/libged/lod/lod.c
      src/libged/view/lod.cpp)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      continue()
    endif()
    file(READ "${_file}" _contents)
    foreach(_pat
	[[struct[ \t\r\n]+bsg_lod_source_policy_settings]]
	[[bsg_view_lod_source_policy_get]]
	[[bsg_view_lod_source_policy_set]]
	[[bsg_view_bounds]])
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} must route BSG LoD policy/bounds hooks through bsg_ged_draw_view.c private adapters: ${_hit}")
      endif()
    endforeach()
  endforeach()

  set(_librt_private_header "${BRLCAD_SOURCE_DIR}/src/librt/librt_private.h")
  if(EXISTS "${_librt_private_header}")
    file(READ "${_librt_private_header}" _librt_private_contents)
    string(FIND "${_librt_private_contents}" "struct rt_mesh_lod_context"
      _rt_mesh_lod_context_idx)
    if(_rt_mesh_lod_context_idx EQUAL -1)
      _brlobol_pivot_guard_fail(
	"src/librt/librt_private.h must keep mesh LoD cache state behind rt_mesh_lod_context")
    endif()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])bsg_mesh_lod_context([^A-Za-z0-9_]|$)]]
      _private_bsg_lod_context_hit "${_librt_private_contents}")
    if(_private_bsg_lod_context_hit)
      _brlobol_pivot_guard_fail(
	"src/librt/librt_private.h reintroduced BSG mesh LoD context state: ${_private_bsg_lod_context_hit}")
    endif()
  endif()

  set(_db_open_impl "${BRLCAD_SOURCE_DIR}/src/librt/db_open.c")
  if(EXISTS "${_db_open_impl}")
    file(READ "${_db_open_impl}" _db_open_contents)
    string(FIND "${_db_open_contents}" "_rt_mesh_lod_context_destroy"
      _rt_mesh_lod_context_destroy_idx)
    if(_rt_mesh_lod_context_destroy_idx EQUAL -1)
      _brlobol_pivot_guard_fail(
	"src/librt/db_open.c must destroy private mesh LoD cache state through _rt_mesh_lod_context_destroy")
    endif()
    set(_db_open_lod_forbidden
      [[#[ \t]*include[ \t]*[<"]bsg/lod\.h]]
      [[(^|[^A-Za-z0-9_])bsg_mesh_lod_context_destroy([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])bsg_mesh_lod_context([^A-Za-z0-9_]|$)]]
    )
    foreach(_pat IN LISTS _db_open_lod_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_db_open_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "src/librt/db_open.c reintroduced direct BSG mesh LoD context ownership: ${_hit}")
      endif()
    endforeach()
  endif()

  set(_ged_db_index_impl "${BRLCAD_SOURCE_DIR}/src/libged/db_index.cpp")
  if(EXISTS "${_ged_db_index_impl}")
    file(READ "${_ged_db_index_impl}" _ged_db_index_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view\.h]]
	db_mesh_lod_invalidate)
      string(REGEX MATCH "${_token}" _ged_db_index_token_hit
	"${_ged_db_index_contents}")
      if(NOT _ged_db_index_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/db_index.cpp must invalidate LoD cache through RT-owned token ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg/lod\.h]]
	[[ged_lod]]
	[[bsg_mesh_lod_key_get]]
	[[bsg_mesh_lod_clear_cache]]
	[[bsg_mesh_lod_key_put]])
      string(REGEX MATCH "${_pat}" _ged_db_index_lod_hit
	"${_ged_db_index_contents}")
      if(_ged_db_index_lod_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/db_index.cpp must not invalidate LoD cache through direct BSG APIs: ${_ged_db_index_lod_hit}")
      endif()
    endforeach()
  endif()

  set(_ged_draw_impl "${BRLCAD_SOURCE_DIR}/src/libged/draw.cpp")
  if(EXISTS "${_ged_draw_impl}")
    file(READ "${_ged_draw_impl}" _ged_draw_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view\.h]]
	db_mesh_lod_status
	db_mesh_lod_refresh
	db_mesh_lod_store_mesh
	db_mesh_lod_get
	rt_mesh_lod_info_get
	rt_mesh_lod_data_get
	ged_draw_mesh_lod_load_view_scene_ref
	ged_draw_brep_mesh_lod_detail_setup)
      string(REGEX MATCH "${_token}" _ged_draw_lod_token_hit
	"${_ged_draw_contents}")
      if(NOT _ged_draw_lod_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/draw.cpp must route BoT draw LoD cache ownership through RT/private GED token ${_token}")
      endif()
    endforeach()
    foreach(_pat
	bsg_mesh_lod_key_get
	bsg_mesh_lod_key_put
	bsg_mesh_lod_clear_cache
	bsg_mesh_lod_cache
	bsg_mesh_lod_create)
      string(REGEX MATCH "${_pat}" _ged_draw_bsg_cache_hit
	"${_ged_draw_contents}")
      if(_ged_draw_bsg_cache_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/draw.cpp must route mesh-LoD cache lookup/write/create through RT helpers: ${_ged_draw_bsg_cache_hit}")
      endif()
    endforeach()
    string(REGEX MATCH [[bsg_mesh_lod_load_view_scene_ref]]
      _ged_draw_bsg_lod_view_load_hit "${_ged_draw_contents}")
    if(_ged_draw_bsg_lod_view_load_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/draw.cpp must route BSG scene-ref mesh-LoD view loading through ged_draw_mesh_lod_load_view_scene_ref: ${_ged_draw_bsg_lod_view_load_hit}")
    endif()
    foreach(_token
	ged_draw_view_perspective_from_bsg
	ged_draw_view_scale_from_bsg)
      string(REGEX MATCH "${_token}" _ged_draw_view_snapshot_token_hit
	"${_ged_draw_contents}")
      if(NOT _ged_draw_view_snapshot_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/draw.cpp must route draw view snapshot reads through private GED draw view helpers")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_perspective([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _ged_draw_view_snapshot_direct_hit
	"${_ged_draw_contents}")
      if(_ged_draw_view_snapshot_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/draw.cpp reintroduced direct BSG draw view snapshot reads: ${_ged_draw_view_snapshot_direct_hit}")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])ged_draw_rt_mesh_lod_bsg([^A-Za-z0-9_]|$)]]
	[[struct[ \t\r\n]+bsg_mesh_lod]]
	)
      string(REGEX MATCH "${_pat}" _ged_draw_bsg_lod_state_hit
	"${_ged_draw_contents}")
      if(_ged_draw_bsg_lod_state_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/draw.cpp must keep mesh-LoD realization on RT handles and private adapters: ${_ged_draw_bsg_lod_state_hit}")
      endif()
    endforeach()
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg/lod\.h]]
	bsg_mesh_lod_detail_setup_clbk
	bsg_mesh_lod_detail_clear_clbk
	bsg_mesh_lod_detail_free_clbk
	brep_mesh_info_clbk
	brep_mesh_info_clear_clbk
	brep_mesh_info_free_clbk)
      string(REGEX MATCH "${_pat}" _ged_draw_brep_lod_callback_hit
	"${_ged_draw_contents}")
      if(_ged_draw_brep_lod_callback_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/draw.cpp must route BRep BSG mesh-LoD detail callbacks through bsg_ged_draw_brep.cpp: ${_ged_draw_brep_lod_callback_hit}")
      endif()
    endforeach()
    foreach(_pat
	bot_mesh_info_clbk
	bot_mesh_info_clear_clbk
	bot_mesh_info_free_clbk)
      string(REGEX MATCH "${_pat}" _ged_draw_bot_lod_callback_hit
	"${_ged_draw_contents}")
      if(_ged_draw_bot_lod_callback_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/draw.cpp must use RT-owned BoT LoD detail callbacks instead of local ${_pat}")
      endif()
    endforeach()
  endif()

  set(_ged_view_lod_impl "${BRLCAD_SOURCE_DIR}/src/libged/view/lod.cpp")
  if(EXISTS "${_ged_view_lod_impl}")
    file(READ "${_ged_view_lod_impl}" _ged_view_lod_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view\.h]]
	db_mesh_lod_clear
	db_mesh_lod_refresh
	db_mesh_lod_status
	db_mesh_lod_store_mesh
	rt_mesh_lod_cache_clear_all)
      string(REGEX MATCH "${_token}" _ged_view_lod_token_hit
	"${_ged_view_lod_contents}")
      if(NOT _ged_view_lod_token_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/lod.cpp must clear LoD caches through RT-owned token ${_token}")
      endif()
    endforeach()
    string(REGEX MATCH [[bsg_mesh_lod_clear_cache]]
      _ged_view_lod_bsg_clear_hit "${_ged_view_lod_contents}")
    if(_ged_view_lod_bsg_clear_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/lod.cpp must not clear LoD caches through direct BSG APIs: ${_ged_view_lod_bsg_clear_hit}")
    endif()
    string(REGEX MATCH [[bsg_mesh_lod_key_get]]
      _ged_view_lod_bsg_key_get_hit "${_ged_view_lod_contents}")
    if(_ged_view_lod_bsg_key_get_hit)
      _brlobol_pivot_guard_fail(
	"src/libged/view/lod.cpp must inspect LoD cache keys through RT status helpers: ${_ged_view_lod_bsg_key_get_hit}")
    endif()
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg/lod\.h]]
	[[bsg_mesh_lod_cache]]
	[[bsg_mesh_lod_key_put]])
      string(REGEX MATCH "${_pat}" _ged_view_lod_bsg_cache_hit
	"${_ged_view_lod_contents}")
      if(_ged_view_lod_bsg_cache_hit)
	_brlobol_pivot_guard_fail(
	  "src/libged/view/lod.cpp must store generated LoD meshes through RT helpers: ${_ged_view_lod_bsg_cache_hit}")
      endif()
    endforeach()
  endif()

  set(_libbsg_lod_header "${BRLCAD_SOURCE_DIR}/include/bsg/lod.h")
  if(EXISTS "${_libbsg_lod_header}")
    file(READ "${_libbsg_lod_header}" _libbsg_lod_header_contents)
    foreach(_token
	"bsg_mesh_lod_cache"
	"one normal per triangle corner"
	"f order, i.e. fcnt * 3 normals")
      string(FIND "${_libbsg_lod_header_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "include/bsg/lod.h must document mesh-LoD normal cache contract token ${_token}")
      endif()
    endforeach()
  endif()
endfunction()

function(_brlobol_pivot_guard_check_librt_sketch_polygon_neutralization)
  set(_sketch_header "${BRLCAD_SOURCE_DIR}/include/rt/primitives/sketch.h")
  if(NOT EXISTS "${_sketch_header}")
    _brlobol_pivot_guard_fail(
      "include/rt/primitives/sketch.h is required for RT sketch polygon checks")
  else()
    file(READ "${_sketch_header}" _sketch_contents)
    foreach(_token
	rt_sketch_polygon
	db_sketch_to_polygon
	db_sketch_polygon_to_sketch
	rt_sketch_polygon_destroy)
      string(FIND "${_sketch_contents}" "${_token}" _token_idx)
      if(_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "include/rt/primitives/sketch.h must expose neutral sketch polygon API ${_token}")
      endif()
    endforeach()
    set(_sketch_forbidden
      [[#[ \t]*include[ \t]*[<"]bsg/]]
      [[(^|[^A-Za-z0-9_])bsg_polygon_ref([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])bsg_view([^A-Za-z0-9_]|$)]]
      [[struct[ \t\r\n]+bsg_view]]
    )
    foreach(_pat IN LISTS _sketch_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_sketch_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "include/rt/primitives/sketch.h reintroduced BSG polygon/view state: ${_hit}")
      endif()
    endforeach()
  endif()

  set(_sketch_impl
    "${BRLCAD_SOURCE_DIR}/src/librt/primitives/sketch/polygons.c")
  if(NOT EXISTS "${_sketch_impl}")
    _brlobol_pivot_guard_fail(
      "src/librt/primitives/sketch/polygons.c is required for neutral RT sketch polygon conversion")
  else()
    file(READ "${_sketch_impl}" _sketch_impl_contents)
    foreach(_token
	"RT_SKETCH_POLYGON_GENERAL"
	"db_sketch_to_polygon"
	"db_sketch_polygon_to_sketch"
	"rt_sketch_polygon_bg_polygon")
      string(FIND "${_sketch_impl_contents}" "${_token}" _sketch_impl_token_idx)
      if(_sketch_impl_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/librt/primitives/sketch/polygons.c must keep neutral sketch polygon implementation token ${_token}")
      endif()
    endforeach()
    set(_sketch_impl_forbidden
      [[#[ \t]*include[ \t]*[<"]bsg/]]
      [[#[ \t]*include[ \t]*[<"]rt/primitives/sketch_legacy_bsg\.h]]
      [[(^|[^A-Za-z0-9_])bsg_polygon_ref([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])bsg_polygon([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])bsg_view([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])BSG_POLYGON_[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])db_sketch_to_view_polygon_ref([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])db_sketch_to_view_polygon_scoped_ref([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])db_view_polygon_ref_to_sketch([^A-Za-z0-9_]|$)]]
    )
    foreach(_pat IN LISTS _sketch_impl_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_sketch_impl_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "src/librt/primitives/sketch/polygons.c reintroduced BSG sketch polygon adapter logic: ${_hit}")
      endif()
    endforeach()
  endif()

  set(_sketch_private_header
    "${BRLCAD_SOURCE_DIR}/src/librt/primitives/sketch/polygons_private.h")
  if(NOT EXISTS "${_sketch_private_header}")
    _brlobol_pivot_guard_fail(
      "src/librt/primitives/sketch/polygons_private.h is required for RT-owned sketch polygon state")
  endif()

  set(_sketch_legacy_header
    "${BRLCAD_SOURCE_DIR}/include/rt/primitives/sketch_legacy_bsg.h")
  if(NOT EXISTS "${_sketch_legacy_header}")
    _brlobol_pivot_guard_fail(
      "include/rt/primitives/sketch_legacy_bsg.h is required for transitional BSG sketch polygon compatibility")
  else()
    file(READ "${_sketch_legacy_header}" _sketch_legacy_contents)
    foreach(_token
	db_sketch_to_view_polygon_ref
	db_sketch_to_view_polygon_scoped_ref
	db_view_polygon_ref_to_sketch)
      string(FIND "${_sketch_legacy_contents}" "${_token}"
	_legacy_token_idx)
      if(_legacy_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "include/rt/primitives/sketch_legacy_bsg.h must expose transitional ${_token} adapter")
      endif()
    endforeach()
  endif()

  set(_sketch_legacy_impl
    "${BRLCAD_SOURCE_DIR}/src/librt/primitives/sketch/polygons_legacy_bsg.c")
  if(NOT EXISTS "${_sketch_legacy_impl}")
    _brlobol_pivot_guard_fail(
      "src/librt/primitives/sketch/polygons_legacy_bsg.c is required for transitional BSG sketch polygon adapters")
  else()
    file(READ "${_sketch_legacy_impl}" _sketch_legacy_impl_contents)
    foreach(_token
	"bsg/polygon.h"
	"bsg/view_state.h"
	"db_sketch_to_view_polygon_ref"
	"db_sketch_to_view_polygon_scoped_ref"
	"db_view_polygon_ref_to_sketch"
	"rt_sketch_polygon_to_bsg"
	"rt_sketch_polygon_from_bsg")
      string(FIND "${_sketch_legacy_impl_contents}" "${_token}"
	_sketch_legacy_impl_token_idx)
      if(_sketch_legacy_impl_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/librt/primitives/sketch/polygons_legacy_bsg.c must own transitional sketch adapter ${_token}")
      endif()
    endforeach()
  endif()

  set(_librt_cmake "${BRLCAD_SOURCE_DIR}/src/librt/CMakeLists.txt")
  if(EXISTS "${_librt_cmake}")
    file(READ "${_librt_cmake}" _librt_cmake_contents)
    string(FIND "${_librt_cmake_contents}" "primitives/sketch/polygons_legacy_bsg.c"
      _sketch_legacy_cmake_idx)
    if(_sketch_legacy_cmake_idx EQUAL -1)
      _brlobol_pivot_guard_fail(
	"src/librt/CMakeLists.txt must build src/librt/primitives/sketch/polygons_legacy_bsg.c")
    endif()
  endif()

  foreach(_rel
      src/libged/view/polygons.c
      src/librt/tests/bsg_poly_sketch.c
      src/qged/plugins/polygon/QPolyCreate.cpp
      src/qged/plugins/polygon/QPolyMod.cpp)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(EXISTS "${_file}")
      file(READ "${_file}" _contents)
      string(FIND "${_contents}" "rt/primitives/sketch_legacy_bsg.h"
	_legacy_include_idx)
      if(_legacy_include_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "${_rel} must include rt/primitives/sketch_legacy_bsg.h for BSG polygon conversion adapters")
      endif()
      if(_rel STREQUAL "src/qged/plugins/polygon/QPolyCreate.cpp" OR
	  _rel STREQUAL "src/qged/plugins/polygon/QPolyMod.cpp")
	foreach(_token
	    "rt/view_legacy_bsg.h"
	    "rt_view_snap_source_flags_set_bsg"
	    "rt_view_snap_lines_set_bsg"
	    "rt_view_snap_lines_from_bsg")
	  string(FIND "${_contents}" "${_token}" _qged_poly_snap_token_idx)
	  if(_qged_poly_snap_token_idx EQUAL -1)
	    _brlobol_pivot_guard_fail(
	      "${_rel} must route polygon snap-state access through rt/view_legacy_bsg.h token ${_token}")
	  endif()
	endforeach()
	foreach(_pat
	    [[(^|[^A-Za-z0-9_])bsg_view_set_snap_source_flags([^A-Za-z0-9_]|$)]]
	    [[(^|[^A-Za-z0-9_])bsg_view_(set_)?snap_lines([^A-Za-z0-9_]|$)]])
	  string(REGEX MATCH "${_pat}" _qged_poly_snap_direct_hit "${_contents}")
	  if(_qged_poly_snap_direct_hit)
	    _brlobol_pivot_guard_fail(
	      "${_rel} reintroduced direct BSG polygon snap-state access: ${_qged_poly_snap_direct_hit}")
	  endif()
	endforeach()
      endif()
      if(_rel STREQUAL "src/qged/plugins/polygon/QPolyMod.cpp")
	foreach(_token
	    "rt_view_center_vec_set_bsg"
	    "rt_view_aet_set_bsg")
	  string(FIND "${_contents}" "${_token}" _qged_poly_view_setter_token_idx)
	  if(_qged_poly_view_setter_token_idx EQUAL -1)
	    _brlobol_pivot_guard_fail(
	      "${_rel} must route polygon alignment view writes through rt/view_legacy_bsg.h token ${_token}")
	  endif()
	endforeach()
	foreach(_pat
	    [[MAT_DELTAS_VEC_NEG[ \t\r\n]*\([^;]*gv_center]]
	    [[VSET[^\n;]*gv_aet]]
	    [[VMOVE[^\n;]*gv_aet]]
	    [[bsg_mat_aet[ \t\r\n]*\(]])
	  string(REGEX MATCH "${_pat}" _qged_poly_view_setter_direct_hit
	    "${_contents}")
	  if(_qged_poly_view_setter_direct_hit)
	    _brlobol_pivot_guard_fail(
	      "${_rel} reintroduced direct BSG polygon alignment view writes: ${_qged_poly_view_setter_direct_hit}")
	  endif()
	endforeach()
      endif()
    endif()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_qtcad_obol_test_links)
  set(_qtcad_tests_cmake
    "${BRLCAD_SOURCE_DIR}/src/libqtcad/tests/CMakeLists.txt")
  if(NOT EXISTS "${_qtcad_tests_cmake}")
    _brlobol_pivot_guard_fail(
      "src/libqtcad/tests/CMakeLists.txt is required for qtcad Obol dependency checks")
    return()
  endif()

  file(READ "${_qtcad_tests_cmake}" _qtcad_tests)
  set(_obol_test_targets
    test_qtcad_obol_controller
    test_qtcad_obol_draw_sync
    test_qtcad_obol_edit_preview
    test_qtcad_obol_export
    test_qtcad_obol_faceplate_sync
    test_qtcad_obol_selection_sync
    test_qtcad_obol_pick
    test_qtcad_obol_snap
    test_qtcad_obol_measure
    test_qtcad_obol_overlay_sync
    test_qtcad_obol_real_models
    test_qtcad_obol_window_host
  )
  set(_legacy_link_targets libdm libbsg dm_plugins)

  foreach(_target IN LISTS _obol_test_targets)
    string(REGEX MATCHALL
      "target_link_libraries\\([ \t\r\n]*${_target}[^)]*\\)"
      _link_blocks "${_qtcad_tests}")
    if(NOT _link_blocks)
      _brlobol_pivot_guard_fail(
	"qtcad Obol test target ${_target} has no guarded target_link_libraries block")
      continue()
    endif()

    foreach(_block IN LISTS _link_blocks)
      foreach(_legacy IN LISTS _legacy_link_targets)
	string(REGEX MATCH "(^|[ \t\r\n])${_legacy}([ \t\r\n]|$)"
	  _legacy_hit "${_block}")
	if(_legacy_hit)
	  _brlobol_pivot_guard_fail(
	    "qtcad Obol test target ${_target} directly links ${_legacy}")
	endif()
      endforeach()
    endforeach()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_qtcad_obol_edit_preview_intent)
  foreach(_rel
      include/qtcad/QgObolEditPreview.h
      src/libqtcad/QgObolEditPreview.cpp
      src/libqtcad/tests/test_qtcad_obol_edit_preview.cpp)
    if(NOT EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _brlobol_pivot_guard_fail("${_rel} is required for qtcad Obol edit-preview intent coverage")
      continue()
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/qtcad/QgObolEditPreview.h" _qtcad_edit_preview_header)
  string(FIND "${_qtcad_edit_preview_header}" "qg_obol_edit_preview_update_with_intent" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("include/qtcad/QgObolEditPreview.h missing explicit edit-preview intent helper")
  endif()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libqtcad/tests/test_qtcad_obol_edit_preview.cpp" _qtcad_edit_preview_test)
  foreach(_token
      "SoBRLExportAction"
      "SoBRLSnapAction"
      "SoBRLMeasureAction"
      "editIntentId"
      "editIntentRole"
      "qg_obol_edit_preview_update_with_intent"
      "custom edit preview export should preserve explicit intent metadata"
      "custom edit preview snap should preserve explicit intent metadata"
      "custom edit preview measure should preserve explicit intent metadata"
      "qtcad edit preview export should preserve edit-intent metadata"
      "qtcad edit preview snap should preserve edit-intent metadata"
      "qtcad edit preview measure should preserve edit-intent metadata"
      "replacement edit preview export should keep edit-intent metadata current")
    string(FIND "${_qtcad_edit_preview_test}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libqtcad/tests/test_qtcad_obol_edit_preview.cpp missing edit-preview intent coverage token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libqtcad/QgObolEditPreview.cpp" _qtcad_edit_preview_impl)
  foreach(_token
      "qg_obol_edit_preview_update_with_intent"
      "obol->replaceEditPreviewWithIntent"
      "display->need_update(QG_VIEW_DRAWN)"
      "obol->removeEditPreview")
    string(FIND "${_qtcad_edit_preview_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libqtcad/QgObolEditPreview.cpp missing Obol edit-preview bridge token ${_token}")
    endif()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_qtcad_obol_export_source_exact)
  set(_qtcad_export_header
    "${BRLCAD_SOURCE_DIR}/include/qtcad/QgObolExport.h")
  set(_qtcad_export_impl
    "${BRLCAD_SOURCE_DIR}/src/libqtcad/QgObolExport.cpp")
  set(_qtcad_export_test
    "${BRLCAD_SOURCE_DIR}/src/libqtcad/tests/test_qtcad_obol_export.cpp")

  foreach(_file IN ITEMS
      "${_qtcad_export_header}"
      "${_qtcad_export_impl}"
      "${_qtcad_export_test}")
    if(NOT EXISTS "${_file}")
      file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
      _brlobol_pivot_guard_fail("${_rel} is required for qtcad Obol source-backed exact export coverage")
    endif()
  endforeach()

  if(EXISTS "${_qtcad_export_header}")
    file(READ "${_qtcad_export_header}" _qtcad_export_header_contents)
    foreach(_token
	"QgObolExportGeometryRecord"
	"QgObolExportTriangleRecord"
	"editIntentId"
	"editIntentRole"
	"submittedSourceRequestCount"
	"sourceFullDetailPending"
	"qg_obol_export_geometry_full_detail"
	"source-backed LoD")
      string(FIND "${_qtcad_export_header_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("include/qtcad/QgObolExport.h missing source-backed exact export token ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${_qtcad_export_impl}")
    file(READ "${_qtcad_export_impl}" _qtcad_export_impl_contents)
    foreach(_token
	"brlobol/export_action.h"
	"SoBRLExportAction::FULL_DETAIL"
	"qg_obol_export_consume_source_full_detail"
	"consumeExportSourceFullDetail"
	"submittedSourceRequestCount"
	"sourceFullDetailPending"
	"triangle.editIntentId"
	"triangle.editIntentRole")
      string(FIND "${_qtcad_export_impl_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libqtcad/QgObolExport.cpp missing source-backed exact export token ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]brlobol/lod_service\.h]]
	[[#[ \t]*include[ \t]*[<"]brlobol/database_source\.h]]
	[[drainMatchingResults]]
	[[submitSourceBackedFullDetailRequests]]
	[[getMaxExactFullDetail(Face|Point)Count]])
      string(REGEX MATCH "${_pat}" _hit "${_qtcad_export_impl_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail("src/libqtcad/QgObolExport.cpp reintroduced controller-owned exact export plumbing: ${_hit}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${_qtcad_export_test}")
    file(READ "${_qtcad_export_test}" _qtcad_export_test_contents)
    foreach(_token
	"SoBRLLodMeshShape"
	"BRLObolLodService"
	"BRLOBOL_LOD_RESULT_FULL_DETAIL"
	"SoBRLExportAction"
	"makeSourceBackedFullDetailLodRequest"
	"qg_obol_export_geometry_full_detail"
	"qtcad exact Obol export should consume ready source-backed full detail"
	"qtcad exact Obol export should preserve source-backed triangle identity and edit intent"
	"qtcad exact Obol export should pass controller full-detail budget to source provider"
	"qtcad exact Obol export should report pending submitted source detail"
	"qtcad source-backed exact Obol export should not initialize the legacy display manager")
      string(FIND "${_qtcad_export_test_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libqtcad/tests/test_qtcad_obol_export.cpp missing source-backed exact export token ${_token}")
      endif()
    endforeach()
  endif()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libqtcad/CMakeLists.txt" _qtcad_cmake)
  string(FIND "${_qtcad_cmake}" "QgObolExport.cpp" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("src/libqtcad/CMakeLists.txt missing QgObolExport.cpp")
  endif()

  file(READ "${BRLCAD_SOURCE_DIR}/include/qtcad/CMakeLists.txt" _qtcad_header_cmake)
  string(FIND "${_qtcad_header_cmake}" "QgObolExport.h" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("include/qtcad/CMakeLists.txt missing QgObolExport.h")
  endif()
endfunction()

function(_brlobol_pivot_guard_check_qtcad_obol_pick_source_exact)
  set(_qtcad_pick_header
    "${BRLCAD_SOURCE_DIR}/include/qtcad/QgObolPick.h")
  set(_qtcad_pick_impl
    "${BRLCAD_SOURCE_DIR}/src/libqtcad/QgObolPick.cpp")
  set(_qtcad_select_filter
    "${BRLCAD_SOURCE_DIR}/src/libqtcad/QgSelectFilter.cpp")
  set(_qtcad_pick_test
    "${BRLCAD_SOURCE_DIR}/src/libqtcad/tests/test_qtcad_obol_pick.cpp")

  foreach(_file IN ITEMS
      "${_qtcad_pick_header}"
      "${_qtcad_pick_impl}"
      "${_qtcad_select_filter}"
      "${_qtcad_pick_test}")
    if(NOT EXISTS "${_file}")
      file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
      _brlobol_pivot_guard_fail("${_rel} is required for qtcad Obol source-backed exact pick coverage")
    endif()
  endforeach()

  if(EXISTS "${_qtcad_pick_header}")
    file(READ "${_qtcad_pick_header}" _qtcad_pick_header_contents)
    foreach(_token
	"QgObolPickRecord"
	"IMPLICIT_SOLID"
	"editIntentId"
	"editIntentRole"
	"faceVertexIndexA"
	"nearestFaceEdgeSlot"
	"nearestFaceEdgeVertexIndexA"
	"nearestFaceVertexSlot"
	"float distance"
	"submittedSourceRequestCount"
	"qg_obol_pick_ray")
      string(FIND "${_qtcad_pick_header_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("include/qtcad/QgObolPick.h missing source-backed exact pick token ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${_qtcad_pick_impl}")
    file(READ "${_qtcad_pick_impl}" _qtcad_pick_impl_contents)
    foreach(_token
		"BRLObolSourceMeshPickResult"
		"BRLObolRtPickResult"
		"qg_obol_normalized_pick_path"
		"qg_obol_pick_path_already_recorded"
		"qg_obol_pick_source_full_detail"
	"qg_obol_pick_camera_line"
	"qg_obol_pick_rt_exact"
	"qg_obol_pick_point_internal"
	"qg_obol_pick_ray"
	"submittedSourceRequestCount"
	"pickAction.setRay(rayOrigin, rayDirection)"
	"qg_obol_insert_pick_record"
	"getEditIntentId"
	"getEditIntentRole"
	"getFaceVertexIndexA"
	"getNearestFaceEdgeSlot"
	"getNearestFaceEdgeVertexIndexA"
	"getNearestFaceVertexSlot"
		"pickSourceMeshExactRay"
		"pickRtExactRay"
	"getLine"
	"distance")
      string(FIND "${_qtcad_pick_impl_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libqtcad/QgObolPick.cpp missing source-backed exact pick token ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${_qtcad_select_filter}")
    file(READ "${_qtcad_select_filter}" _qtcad_select_filter_contents)
    foreach(_token
	"qg_obol_pick_point(view_widget()"
	"qg_obol_pick_rect(view_widget()"
	"qg_obol_pick_ray(view_widget()"
	"submittedSourceRequests"
	"submittedSourceRequests > 0"
	"set_selected_paths(v, paths);")
      string(FIND "${_qtcad_select_filter_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libqtcad/QgSelectFilter.cpp missing source-backed exact pick selection token ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${_qtcad_pick_test}")
    file(READ "${_qtcad_pick_test}" _qtcad_pick_test_contents)
    foreach(_token
	"SoBRLLodMeshShape"
	"BRLObolLodService"
	"BRLOBOL_LOD_RESULT_FULL_DETAIL"
	"SoRayPickAction"
	"SoBRLSourceMeshPickAction"
	"mk_bot(wdbp, \"lod-pick.bot\""
	"makeSourceBackedFullDetailLodRequest"
	"qtcad LoD pick fixture should build a bounded source full-detail request"
	"queuedResultCountForDiagnostics"
	"qtcad Obol point pick should consume ready source-backed full detail"
	"qtcad Obol source-backed pick should preserve exact mesh, sub-entity identity, and edit intent"
	"nearestFaceEdgeVertexIndexA"
	"qtcad source-backed exact Obol pick should pass controller full-detail budget to source provider"
	"qtcad source-backed exact Obol pick should report pending submitted source detail"
	"qtcad point select filter should consume pending source-backed exact pick"
	"qtcad point select filter should not fall back to legacy selection while exact source pick is pending"
	"qtcad point select filter should submit the pending source-backed exact pick request"
	"qtcad box select filter should consume pending source-backed exact pick"
	"qtcad box select filter should not fall back to legacy selection while exact source pick is pending"
	"qtcad box select filter should submit pending source-backed exact pick requests"
	"qtcad ray select filter should consume pending source-backed exact pick"
	"qtcad ray select filter should not fall back to legacy selection while exact source pick is pending"
	"qtcad ray select filter should submit the pending source-backed exact pick request"
	"qtcad Obol point pick should use librt exact implicit pick"
	"qtcad Obol librt exact implicit pick should preserve RT hit identity"
	"qtcad Obol librt exact implicit pick should retain a controller RT pick cache"
	"qtcad Obol rectangle pick should reuse controller-cached librt exact implicit picks"
	"qtcad Obol rectangle pick should keep reusing the controller RT pick cache"
	"qtcad Obol cached rectangle librt pick should preserve RT hit identity"
	"qtcad Obol ray pick should reuse controller-cached librt exact implicit picks"
	"qtcad Obol cached ray librt pick should preserve RT hit identity"
	"qtcad select ray filter should expose Obol explicit-ray implicit paths"
	"qtcad select ray filter should reuse the controller RT pick cache"
	"qtcad select ray filter should run Obol explicit-ray selection without legacy dbip"
	"qtcad select ray filter without legacy dbip should expose Obol explicit-ray implicit paths"
	"qtcad select ray filter without legacy dbip should reuse the controller RT pick cache"
	"qtcad mixed Obol/RT point pick should choose the nearer librt hit"
	"qtcad mixed Obol/RT pick-all should order merged hits by distance"
	"qtcad mixed Obol/RT ray pick-all should order merged hits by distance"
	"QgObolPickRecord::IMPLICIT_SOLID")
      string(FIND "${_qtcad_pick_test_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libqtcad/tests/test_qtcad_obol_pick.cpp missing source-backed exact pick token ${_token}")
      endif()
    endforeach()
  endif()
endfunction()

function(_brlobol_pivot_guard_check_qtcad_obol_snap_source_exact)
  set(_qtcad_snap_header
    "${BRLCAD_SOURCE_DIR}/include/qtcad/QgObolSnap.h")
  set(_qtcad_snap_impl
    "${BRLCAD_SOURCE_DIR}/src/libqtcad/QgObolSnap.cpp")
  set(_qtcad_snap_test
    "${BRLCAD_SOURCE_DIR}/src/libqtcad/tests/test_qtcad_obol_snap.cpp")

  foreach(_file IN ITEMS
      "${_qtcad_snap_header}"
      "${_qtcad_snap_impl}"
      "${_qtcad_snap_test}")
    if(NOT EXISTS "${_file}")
      file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
      _brlobol_pivot_guard_fail("${_rel} is required for qtcad Obol source-backed exact snap coverage")
    endif()
  endforeach()

  if(EXISTS "${_qtcad_snap_header}")
    file(READ "${_qtcad_snap_header}" _qtcad_snap_header_contents)
    foreach(_token
	"VERTEX"
	"EDGE_NEAREST"
	"vertexIndex"
	"edgeSlot"
	"edgeVertexIndexA"
	"edgeVertexIndexB"
	"submittedSourceRequestCount"
	"sourceFullDetailPending"
	"editIntentId"
	"editIntentRole"
	"qg_obol_snap_point_full_detail"
	"source-backed LoD")
      string(FIND "${_qtcad_snap_header_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("include/qtcad/QgObolSnap.h missing source-backed exact snap token ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${_qtcad_snap_impl}")
    file(READ "${_qtcad_snap_impl}" _qtcad_snap_impl_contents)
    foreach(_token
	"qg_obol_snap_point_with_policy"
	"qg_obol_snap_point_full_detail"
	"SoBRLSnapAction::FULL_DETAIL"
	"getVertexIndex"
	"getEdgeSlot"
	"getEdgeVertexIndexA"
	"getEdgeVertexIndexB"
	"getEditIntentId"
	"getEditIntentRole"
	"qg_obol_snap_consume_source_full_detail"
	"submittedSourceRequestCount"
	"sourceFullDetailPending"
	"consumeSnapSourceFullDetail")
      string(FIND "${_qtcad_snap_impl_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libqtcad/QgObolSnap.cpp missing source-backed exact snap token ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]brlobol/lod_service\.h]]
	[[#[ \t]*include[ \t]*[<"]brlobol/database_source\.h]]
	[[drainMatchingResults]]
	[[submitSourceBackedFullDetailRequests]]
	[[getMaxExactFullDetail(Face|Point)Count]])
      string(REGEX MATCH "${_pat}" _hit "${_qtcad_snap_impl_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail("src/libqtcad/QgObolSnap.cpp reintroduced controller-owned exact snap plumbing: ${_hit}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${_qtcad_snap_test}")
    file(READ "${_qtcad_snap_test}" _qtcad_snap_test_contents)
    foreach(_token
	"SoBRLLodMeshShape"
	"BRLObolLodService"
	"rt/view_legacy_bsg.h"
	"rt_view_dimensions_set_bsg"
	"rt_view_size_set_bsg"
	"rt_view_view2model_set_bsg"
	"rt_view_snap_source_flags_set_bsg"
	"rt_view_snap_lines_set_bsg"
	"BRLOBOL_LOD_RESULT_FULL_DETAIL"
	"SoBRLSnapAction"
	"mk_bot(wdbp, \"lod-snap.bot\""
	"makeSourceBackedFullDetailLodRequest"
	"qtcad LoD snap fixture should build a bounded source full-detail request"
	"qg_obol_snap_point_full_detail"
	"qtcad exact Obol snap should consume ready source-backed full detail"
	"qtcad exact Obol snap should preserve source-backed face identity and edit intent"
	"qtcad exact Obol snap should preserve source-backed vertex identity"
	"qtcad exact Obol snap should preserve source-backed edge identity"
	"qtcad exact Obol snap should report pending submitted source detail"
	"qtcad exact Obol snap should pass controller full-detail budget to source provider")
      string(FIND "${_qtcad_snap_test_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libqtcad/tests/test_qtcad_obol_snap.cpp missing source-backed exact snap token ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[v->[ \t\r\n]*gv_(width|height|size)[ \t\r\n]*=]]
	[[MAT_(IDN|DELTAS)[ \t\r\n]*\([^;]*v->[ \t\r\n]*gv_view2model]]
	[[(^|[^A-Za-z0-9_])bsg_view_set_snap_source_flags([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_set_snap_lines([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _hit "${_qtcad_snap_test_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/tests/test_qtcad_obol_snap.cpp reintroduced direct BSG snap test setup writes: ${_hit}")
      endif()
    endforeach()
  endif()
endfunction()

function(_brlobol_pivot_guard_check_qtcad_obol_measure_source_exact)
  set(_qtcad_measure_header
    "${BRLCAD_SOURCE_DIR}/include/qtcad/QgObolMeasure.h")
  set(_qtcad_measure_impl
    "${BRLCAD_SOURCE_DIR}/src/libqtcad/QgObolMeasure.cpp")
  set(_qtcad_measure_test
    "${BRLCAD_SOURCE_DIR}/src/libqtcad/tests/test_qtcad_obol_measure.cpp")

  foreach(_file IN ITEMS
      "${_qtcad_measure_header}"
      "${_qtcad_measure_impl}"
      "${_qtcad_measure_test}")
    if(NOT EXISTS "${_file}")
      file(RELATIVE_PATH _rel "${BRLCAD_SOURCE_DIR}" "${_file}")
      _brlobol_pivot_guard_fail("${_rel} is required for qtcad Obol source-backed exact measure coverage")
    endif()
  endforeach()

  if(EXISTS "${_qtcad_measure_header}")
    file(READ "${_qtcad_measure_header}" _qtcad_measure_header_contents)
    foreach(_token
	"QgObolMeasureGeometryRecord"
	"nearestFaceVertexIndexA"
	"nearestFaceEdgeSlot"
	"nearestFaceEdgeVertexIndexA"
	"nearestFaceVertexSlot"
	"nearestEditIntentId"
	"nearestEditIntentRole"
	"submittedSourceRequestCount"
	"sourceFullDetailPending"
	"qg_obol_measure_geometry_full_detail")
      string(FIND "${_qtcad_measure_header_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("include/qtcad/QgObolMeasure.h missing source-backed exact measure token ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${_qtcad_measure_impl}")
    file(READ "${_qtcad_measure_impl}" _qtcad_measure_impl_contents)
    foreach(_token
	"brlobol/measure_action.h"
	"SoBRLMeasureAction::FULL_DETAIL"
	"qg_obol_measure_consume_source_full_detail"
	"consumeMeasureSourceFullDetail"
	"submittedSourceRequestCount"
	"sourceFullDetailPending"
	"getNearestFaceVertexIndexA"
	"getNearestFaceEdgeSlot"
	"getNearestFaceEdgeVertexIndexA"
	"getNearestFaceVertexSlot"
	"getNearestEditIntentId"
	"getNearestEditIntentRole")
      string(FIND "${_qtcad_measure_impl_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libqtcad/QgObolMeasure.cpp missing source-backed exact measure token ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]brlobol/lod_service\.h]]
	[[#[ \t]*include[ \t]*[<"]brlobol/database_source\.h]]
	[[drainMatchingResults]]
	[[submitSourceBackedFullDetailRequests]]
	[[getMaxExactFullDetail(Face|Point)Count]])
      string(REGEX MATCH "${_pat}" _hit "${_qtcad_measure_impl_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail("src/libqtcad/QgObolMeasure.cpp reintroduced controller-owned exact measure plumbing: ${_hit}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${_qtcad_measure_test}")
    file(READ "${_qtcad_measure_test}" _qtcad_measure_test_contents)
    foreach(_token
	"SoBRLLodMeshShape"
	"BRLObolLodService"
	"BRLOBOL_LOD_RESULT_FULL_DETAIL"
	"SoBRLMeasureAction"
	"makeSourceBackedFullDetailLodRequest"
	"query-scoped source full-detail request"
	"qg_obol_measure_geometry_full_detail"
	"qtcad exact Obol measure should consume ready source-backed full detail"
	"qtcad exact Obol measure should preserve source-backed face, sub-entity identity, and edit intent"
	"nearestFaceEdgeVertexIndexA"
	"qtcad exact Obol measure should report pending submitted source detail"
	"qtcad exact Obol measure should pass controller full-detail budget to source provider")
      string(FIND "${_qtcad_measure_test_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libqtcad/tests/test_qtcad_obol_measure.cpp missing source-backed exact measure token ${_token}")
      endif()
    endforeach()
  endif()
endfunction()

function(_brlobol_pivot_guard_check_qtcad_window_host_adapter)
  set(_required_files
    include/qtcad/QgObolWindowHost.h
    src/libqtcad/QgObolWindowHost.cpp
    src/libqtcad/tests/test_qtcad_obol_window_host.cpp)

  foreach(_rel IN LISTS _required_files)
    if(NOT EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _brlobol_pivot_guard_fail("${_rel} is required for qtcad Obol window-host adapter coverage")
    endif()
  endforeach()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/include/qtcad/QgObolWindowHost.h")
    file(READ "${BRLCAD_SOURCE_DIR}/include/qtcad/QgObolWindowHost.h" _qtcad_window_header)
    foreach(_token
	"QgObolWindowHost"
	"BRLObolWindowHost"
	"setCanvas"
	"lastFrameImage"
	"lastRenderReason")
      string(FIND "${_qtcad_window_header}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("include/qtcad/QgObolWindowHost.h missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libqtcad/QgObolWindowHost.cpp")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libqtcad/QgObolWindowHost.cpp" _qtcad_window_impl)
    foreach(_token
	"QgCanvasBase"
	"QgSW"
	"ownsCanvas"
	"qg_obol_window_host_create_canvas"
	"widget->show()"
	"obolViewController"
	"get_obol_viewport_image"
	"consumeRenderRequest"
	"BRLOBOL_WINDOW_BACKEND_QT")
      string(FIND "${_qtcad_window_impl}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libqtcad/QgObolWindowHost.cpp missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libqtcad/tests/test_qtcad_obol_window_host.cpp")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libqtcad/tests/test_qtcad_obol_window_host.cpp" _qtcad_window_test)
    foreach(_token
	"QgObolWindowHost"
	"QgSW"
	"brlobol_window_host_register_display_host"
	"imgstream_fb_open"
	"\"/dev/qtgl\""
	"imgstream_fb_poll"
	"lastFrameImage"
	"test_qtcad_owned_window_host"
	"test_qtcad_owned_display_host_bridge"
	"isWindow"
	"host.canvas() == NULL")
      string(FIND "${_qtcad_window_test}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("src/libqtcad/tests/test_qtcad_obol_window_host.cpp missing ${_token}")
      endif()
    endforeach()
  endif()

  if(EXISTS "${BRLCAD_SOURCE_DIR}/src/libqtcad/CMakeLists.txt")
    file(READ "${BRLCAD_SOURCE_DIR}/src/libqtcad/CMakeLists.txt" _qtcad_cmake)
    string(FIND "${_qtcad_cmake}" "QgObolWindowHost.cpp" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libqtcad/CMakeLists.txt missing QgObolWindowHost.cpp")
    endif()
  endif()

  foreach(_rel IN LISTS _required_files)
    if(NOT EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      continue()
    endif()
    file(READ "${BRLCAD_SOURCE_DIR}/${_rel}" _contents)
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg]]
	[[#[ \t]*include[ \t]*[<"]dm]]
	[[(^|[^A-Za-z0-9_])libbsg([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail("${_rel} reintroduced a legacy display dependency: ${_hit}")
      endif()
    endforeach()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_qtcad_measure_filter)
  set(_measure_filter_files
    include/qtcad/QgMeasureFilter.h
    src/libqtcad/QgMeasureFilter.cpp
  )
  set(_measure_forbidden
    [[#[ \t]*include[ \t]*[<"]bsg/interaction\.h]]
    [[#[ \t]*include[ \t]*[<"]bsg/feature\.h]]
    [[#[ \t]*include[ \t]*[<"]bsg/geometry\.h]]
    [[#[ \t]*include[ \t]*[<"]raytrace\.h]]
    [[(^|[^A-Za-z0-9_])bsg_interaction_[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])bsg_pick_point([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])bsg_feature_ref([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])bsg_feature_[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])bsg_measure_result([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])bsg_measure_candidates([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])rt_shootray([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])rt_gettrees_and_attrs([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])struct[ \t]+application([^A-Za-z0-9_]|$)]]
    [[(^|[^A-Za-z0-9_])struct[ \t]+rt_i([^A-Za-z0-9_]|$)]]
  )

  foreach(_rel IN LISTS _measure_filter_files)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail("${_rel} is required for qtcad measurement migration checks")
      continue()
    endif()

    file(READ "${_file}" _contents)
    foreach(_pat IN LISTS _measure_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel}: qtcad measurement filter reintroduced legacy BSG measurement fallback symbol: ${_hit}")
      endif()
    endforeach()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_qtcad_selection_api)
  set(_select_header "${BRLCAD_SOURCE_DIR}/include/qtcad/QgSelectFilter.h")
  if(NOT EXISTS "${_select_header}")
    _brlobol_pivot_guard_fail(
      "include/qtcad/QgSelectFilter.h is required for qtcad selection API migration checks")
  else()
    file(READ "${_select_header}" _select_header_contents)
    set(_select_header_forbidden
      [[#[ \t]*include[ \t]*[<"]bsg]]
      [[#[ \t]*include[ \t]*[<"]raytrace\.h]]
      [[(^|[^A-Za-z0-9_])bsg_pick_result([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])bsg_interaction_[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])pick_result[ \t]*\([^A-Za-z0-9_]*\)]]
      [[(^|[^A-Za-z0-9_])interaction_result[ \t]*\([^A-Za-z0-9_]*\)]]
    )
    foreach(_pat IN LISTS _select_header_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_select_header_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "include/qtcad/QgSelectFilter.h reintroduced public legacy BSG/raytrace selection API: ${_hit}")
      endif()
    endforeach()
  endif()

  set(_selector_plugin
    "${BRLCAD_SOURCE_DIR}/src/qged/plugins/view/select/CADViewSelector.cpp")
  if(NOT EXISTS "${_selector_plugin}")
    _brlobol_pivot_guard_fail(
      "src/qged/plugins/view/select/CADViewSelector.cpp is required for qged selector migration checks")
  else()
    file(READ "${_selector_plugin}" _selector_contents)
    set(_selector_forbidden
      [[#[ \t]*include[ \t]*[<"]bsg/interaction\.h]]
      [[#[ \t]*include[ \t]*[<"]bsg/node\.h]]
      [[(^|[^A-Za-z0-9_])bsg_interaction_[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
    )
    foreach(_pat IN LISTS _selector_forbidden)
      string(REGEX MATCH "${_pat}" _hit "${_selector_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "qged selector plugin reintroduced direct BSG interaction consumption: ${_hit}")
      endif()
    endforeach()
  endif()
endfunction()

function(_brlobol_pivot_guard_check_qtcad_view_snapshot_adapters)
  set(_qtcad_canvas_state "${BRLCAD_SOURCE_DIR}/src/libqtcad/QgCanvasState.h")
  if(EXISTS "${_qtcad_canvas_state}")
    file(READ "${_qtcad_canvas_state}" _qtcad_canvas_state_contents)
    foreach(_token
	"rt/view_legacy_bsg.h"
	"rt_view_rotation_from_bsg"
	"rt_view_center_from_bsg"
	"rt_view_scale_from_bsg"
	"rt_view_perspective_from_bsg"
	"rt_view_dimensions_set_bsg"
	"rt_view_aet_set_bsg")
      string(FIND "${_qtcad_canvas_state_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/QgCanvasState.h must route Obol camera view snapshots through rt/view_legacy_bsg.h token ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[s[.]v->[ \t\r\n]*gv_(width|height)[ \t\r\n]*=]]
	[[s[.]v->[ \t\r\n]*gv_(rotation|scale|perspective)]]
	[[bsg_view_get_center_vec]]
	[[(^|[^A-Za-z0-9_])bsg_view_set_aet([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _hit "${_qtcad_canvas_state_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/QgCanvasState.h reintroduced direct qtcad Obol camera BSG view reads: ${_hit}")
      endif()
    endforeach()
  endif()

  set(_qtcad_draw_sync_test "${BRLCAD_SOURCE_DIR}/src/libqtcad/tests/test_qtcad_obol_draw_sync.cpp")
  if(EXISTS "${_qtcad_draw_sync_test}")
    file(READ "${_qtcad_draw_sync_test}" _qtcad_draw_sync_test_contents)
    foreach(_token
	"rt/view_legacy_bsg.h"
	"rt_view_center_vec_set_bsg"
	"rt_view_scale_set_bsg")
      string(FIND "${_qtcad_draw_sync_test_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/tests/test_qtcad_obol_draw_sync.cpp must route view sync setup writes through rt/view_legacy_bsg.h token ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])bsg_view_set_center_vec([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_set_scale([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _hit "${_qtcad_draw_sync_test_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/tests/test_qtcad_obol_draw_sync.cpp reintroduced direct BSG utility view setters: ${_hit}")
      endif()
    endforeach()
  endif()

  set(_qtcad_qsketch_test "${BRLCAD_SOURCE_DIR}/src/libqtcad/tests/qsketch.cpp")
  if(EXISTS "${_qtcad_qsketch_test}")
    file(READ "${_qtcad_qsketch_test}" _qtcad_qsketch_test_contents)
    foreach(_token
	"rt/view_legacy_bsg.h"
	"rt_view_aet_set_bsg"
	"rt_view_scale_set_bsg"
	"rt_view_dimensions_set_bsg"
	"rt_view_center_vec_set_bsg"
	"rt_view_size_set_bsg")
      string(FIND "${_qtcad_qsketch_test_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/tests/qsketch.cpp must route view setup writes through rt/view_legacy_bsg.h token ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[m_bv->[ \t\r\n]*gv_(width|height|scale|size|isize|aet)[ \t\r\n]*=]]
	[[MAT_(IDN|DELTAS_VEC_NEG)[ \t\r\n]*\([^;]*m_bv->[ \t\r\n]*gv_center]]
	[[(^|[^A-Za-z0-9_])bsg_mat_aet[ \t\r\n]*\([^;]*m_bv]])
      string(REGEX MATCH "${_pat}" _hit "${_qtcad_qsketch_test_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/tests/qsketch.cpp reintroduced direct BSG qsketch view setup writes: ${_hit}")
      endif()
    endforeach()
  endif()

  set(_qtcad_view_ctrl "${BRLCAD_SOURCE_DIR}/src/libqtcad/QgViewCtrl.cpp")
  if(EXISTS "${_qtcad_view_ctrl}")
    file(READ "${_qtcad_view_ctrl}" _qtcad_view_ctrl_contents)
    foreach(_token
	"rt/view_legacy_bsg.h"
	"rt_view_framebuffer_mode_from_bsg"
	"rt_view_framebuffer_mode_set_bsg")
      string(FIND "${_qtcad_view_ctrl_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/QgViewCtrl.cpp must route framebuffer mode view-state access through rt/view_legacy_bsg.h token ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg/view_state\.h]]
	[[(^|[^A-Za-z0-9_])bsg_view_(set_)?framebuffer_mode([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _qtcad_view_ctrl_fb_direct
	"${_qtcad_view_ctrl_contents}")
      if(_qtcad_view_ctrl_fb_direct)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/QgViewCtrl.cpp reintroduced direct BSG framebuffer mode view-state access: ${_qtcad_view_ctrl_fb_direct}")
      endif()
    endforeach()
  endif()

  set(_qtcad_view_filter "${BRLCAD_SOURCE_DIR}/src/libqtcad/QgViewFilter.cpp")
  if(EXISTS "${_qtcad_view_filter}")
    file(READ "${_qtcad_view_filter}" _qtcad_view_filter_contents)
    foreach(_token
	"rt/view_legacy_bsg.h"
	"rt_view_width_from_bsg"
	"rt_view_height_from_bsg"
	"rt_view_size_from_bsg"
	"rt_view_snap_lines_from_bsg"
	"rt_view_snap_source_flags_from_bsg"
	"rt_view_snap_kind_mask_from_bsg"
	"rt_view_screen_point_from_bsg"
	"rt_view_snap_tolerance_factor_from_bsg")
      string(FIND "${_qtcad_view_filter_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/QgViewFilter.cpp must route snap tolerance view reads through rt/view_legacy_bsg.h token ${_token}")
      endif()
    endforeach()
    string(REGEX MATCH [[v->[ \t\r\n]*gv_(width|height|size)]]
      _qtcad_view_filter_direct "${_qtcad_view_filter_contents}")
    if(_qtcad_view_filter_direct)
      _brlobol_pivot_guard_fail(
	"src/libqtcad/QgViewFilter.cpp reintroduced direct snap tolerance BSG view reads: ${_qtcad_view_filter_direct}")
    endif()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])bsg_view_snap_tolerance_factor([^A-Za-z0-9_]|$)]]
      _qtcad_view_filter_snap_tol_direct "${_qtcad_view_filter_contents}")
    if(_qtcad_view_filter_snap_tol_direct)
      _brlobol_pivot_guard_fail(
	"src/libqtcad/QgViewFilter.cpp reintroduced direct BSG snap tolerance factor reads: ${_qtcad_view_filter_snap_tol_direct}")
    endif()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])bsg_view_snap_lines([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_snap_source_flags([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_view_snap_kind_mask([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _qtcad_view_filter_snap_policy_direct "${_qtcad_view_filter_contents}")
      if(_qtcad_view_filter_snap_policy_direct)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/QgViewFilter.cpp reintroduced direct BSG snap policy reads: ${_qtcad_view_filter_snap_policy_direct}")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^A-Za-z0-9_])bsg_screen_pt([^A-Za-z0-9_]|$)]]
      _qtcad_view_filter_screen_pt_direct "${_qtcad_view_filter_contents}")
    if(_qtcad_view_filter_screen_pt_direct)
      _brlobol_pivot_guard_fail(
	"src/libqtcad/QgViewFilter.cpp reintroduced direct BSG screen-point conversion: ${_qtcad_view_filter_screen_pt_direct}")
    endif()
  endif()

  foreach(_rel
      src/libqtcad/QgCanvasInput.cpp
      src/libqtcad/bindings.cpp)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail("${_rel} is required for qtcad input view adapter checks")
      continue()
    endif()
    file(READ "${_file}" _contents)
    foreach(_token
	"rt/view_legacy_bsg.h"
	"rt_view_height_from_bsg"
	"rt_view_center_from_bsg"
	"rt_view_bounds_update_callback_from_bsg"
	"rt_view_bounds_update_callback_set_bsg"
	"rt_view_bounds_update_callback_call_bsg"
	"rt_view_aet_set_bsg")
      string(FIND "${_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "${_rel} must route input view reads/callbacks through rt/view_legacy_bsg.h token ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[v->[ \t\r\n]*gv_height]]
	[[gv_bounds_update]]
	[[MAT_DELTAS_GET_NEG[ \t\r\n]*\([^;]*v->[ \t\r\n]*gv_center]]
	[[(^|[^A-Za-z0-9_])bsg_view_set_aet([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _hit "${_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced direct qtcad input BSG view reads: ${_hit}")
      endif()
    endforeach()
  endforeach()

  set(_qtcad_sketch_filter "${BRLCAD_SOURCE_DIR}/src/libqtcad/QgSketchFilter.cpp")
  if(EXISTS "${_qtcad_sketch_filter}")
    file(READ "${_qtcad_sketch_filter}" _qtcad_sketch_filter_contents)
    foreach(_token
	"rt/view_legacy_bsg.h"
	"rt_view_model2view_from_bsg"
	"rt_view_screen_to_view_from_bsg"
	"rt_view_screen_point_from_bsg"
	"rt_view_width_from_bsg")
      string(FIND "${_qtcad_sketch_filter_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/QgSketchFilter.cpp must route sketch view snapshots through rt/view_legacy_bsg.h token ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[v->[ \t\r\n]*gv_model2view]]
	[[v->[ \t\r\n]*gv_width]]
	[[bsg_screen_to_view]]
	[[bsg_screen_pt]])
      string(REGEX MATCH "${_pat}" _hit "${_qtcad_sketch_filter_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/QgSketchFilter.cpp reintroduced direct sketch BSG view reads: ${_hit}")
      endif()
    endforeach()
  endif()

  foreach(_rel
      src/libqtcad/QgPolyFilter.cpp
      src/libged/view/objs.cpp
      src/libged/view/polygons.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(EXISTS "${_file}")
      file(READ "${_file}" _contents)
      foreach(_token
	  [[rt/view_legacy_bsg.h]]
	  [[rt_view_screen_point_from_bsg]])
	string(FIND "${_contents}" "${_token}" _idx)
	if(_idx EQUAL -1)
	  _brlobol_pivot_guard_fail(
	    "${_rel} must route screen-point conversion through rt/view_legacy_bsg.h token ${_token}")
	endif()
      endforeach()
      string(REGEX MATCH [[(^|[^A-Za-z0-9_])bsg_screen_pt([^A-Za-z0-9_]|$)]]
	_screen_point_direct "${_contents}")
      if(_screen_point_direct)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced direct BSG screen-point conversion: ${_screen_point_direct}")
      endif()
    endif()
  endforeach()

  set(_qtcad_select_filter "${BRLCAD_SOURCE_DIR}/src/libqtcad/QgSelectFilter.cpp")
  if(EXISTS "${_qtcad_select_filter}")
    file(READ "${_qtcad_select_filter}" _qtcad_select_filter_contents)
    foreach(_token
	"rt/view_legacy_bsg.h"
	"rt_view_width_from_bsg"
	"rt_view_height_from_bsg"
	"rt_view_radius_from_bsg"
	"rt_view_screen_to_view_from_bsg"
	"rt_view_view2model_from_bsg"
	"rt_view_rotation_from_bsg"
	"_qg_select_ray_from_view"
	"qg_obol_pick_ray(view_widget()")
      string(FIND "${_qtcad_select_filter_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/QgSelectFilter.cpp must route selection view snapshots through rt/view_legacy_bsg.h token ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[v->[ \t\r\n]*gv_(width|height|view2model|rotation)]]
	[[v->[ \t\r\n]*radius]]
	[[bsg_screen_to_view]])
      string(REGEX MATCH "${_pat}" _hit "${_qtcad_select_filter_contents}")
      if(_hit)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/QgSelectFilter.cpp reintroduced direct selection BSG view reads: ${_hit}")
      endif()
    endforeach()
  endif()

  foreach(_rel
      src/libqtcad/QgSelectFilter.cpp
      src/libged/rect/rect.c
      src/libged/select/select.c
      src/libged/view/faceplate/interactive_rect.c
      src/libtclcad/commands.c
      src/libtclcad/view/refresh.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(EXISTS "${_file}")
      file(READ "${_file}" _contents)
      foreach(_token
	  [[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	  [[rt_view_interactive_rect_from_bsg]])
	string(REGEX MATCH "${_token}" _interactive_rect_token_hit
	  "${_contents}")
	if(NOT _interactive_rect_token_hit)
	  _brlobol_pivot_guard_fail(
	    "${_rel} must route interactive rectangle state through rt/view_legacy_bsg.h token ${_token}")
	endif()
      endforeach()
      string(REGEX MATCH [[(^|[^A-Za-z0-9_])bsg_view_interactive_rect_(get|set)([^A-Za-z0-9_]|$)]]
	_interactive_rect_direct_hit "${_contents}")
      if(_interactive_rect_direct_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced direct BSG interactive rectangle state access: ${_interactive_rect_direct_hit}")
      endif()
    endif()
  endforeach()

  foreach(_rel
      src/libqtcad/QgCanvasState.h
      src/libqtcad/QgView.cpp
      src/libqtcad/QgGL.cpp
      src/libqtcad/QgSW.cpp
      src/libtclcad/view/refresh.c
      src/mged/usepen.c
      src/mged/mged.c
      src/mged/dm-generic.c
      src/mged/dozoom.c)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(EXISTS "${_file}")
      file(READ "${_file}" _contents)
      foreach(_token
	  [[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	  [[rt_view_refresh_]])
	string(REGEX MATCH "${_token}" _refresh_token_hit
	  "${_contents}")
	if(NOT _refresh_token_hit)
	  _brlobol_pivot_guard_fail(
	    "${_rel} must route refresh state through rt/view_legacy_bsg.h token ${_token}")
	endif()
      endforeach()
      string(REGEX MATCH [[(^|[^A-Za-z0-9_])bsg_view_refresh_(enabled|set_enabled|suppressed|suppress_begin|suppress_end|request|consume|complete|dirty|drawn_count|set_drawn_count)([^A-Za-z0-9_]|$)]]
	_refresh_direct_hit "${_contents}")
      if(_refresh_direct_hit)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced direct BSG refresh-state access: ${_refresh_direct_hit}")
      endif()
    endif()
  endforeach()

  set(_qtcad_measure_filter "${BRLCAD_SOURCE_DIR}/src/libqtcad/QgMeasureFilter.cpp")
  if(EXISTS "${_qtcad_measure_filter}")
    file(READ "${_qtcad_measure_filter}" _qtcad_measure_filter_contents)
    foreach(_token
	"rt/view_legacy_bsg.h"
	"rt_view_screen_to_view_from_bsg"
	"rt_view_view2model_from_bsg")
      string(FIND "${_qtcad_measure_filter_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/QgMeasureFilter.cpp must route 2D measurement view-to-model reads through rt/view_legacy_bsg.h token ${_token}")
      endif()
    endforeach()
    string(REGEX MATCH [[(MAT4X3PNT[ \t\r\n]*\([^;]*v->[ \t\r\n]*gv_view2model|bsg_screen_to_view)]]
      _qtcad_measure_direct "${_qtcad_measure_filter_contents}")
    if(_qtcad_measure_direct)
      _brlobol_pivot_guard_fail(
	"src/libqtcad/QgMeasureFilter.cpp reintroduced direct 2D measurement BSG view reads: ${_qtcad_measure_direct}")
    endif()
  endif()

  set(_qtcad_quad_view "${BRLCAD_SOURCE_DIR}/src/libqtcad/QgQuadView.cpp")
  if(EXISTS "${_qtcad_quad_view}")
    file(READ "${_qtcad_quad_view}" _qtcad_quad_view_contents)
    foreach(_token
	"rt/view_legacy_bsg.h"
	"rt_view_base2local_from_bsg"
	"rt_view_local2base_from_bsg"
	"rt_view_width_from_bsg"
	"rt_view_height_from_bsg"
	"rt_view_dimensions_set_bsg"
	"brlobol/export_action.h"
	"SoBRLExportAction"
	"SoBRLExportAction::DISPLAY_LEVEL"
	"qg_quad_obol_visible_paths")
      string(FIND "${_qtcad_quad_view_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/QgQuadView.cpp missing qtcad quad-view migration token ${_token}")
      endif()
    endforeach()
    string(REGEX MATCH [[=[^;]*->[ \t\r\n]*gv_(base2local|local2base|width|height)]]
      _qtcad_quad_direct "${_qtcad_quad_view_contents}")
    if(_qtcad_quad_direct)
      _brlobol_pivot_guard_fail(
	"src/libqtcad/QgQuadView.cpp reintroduced direct source BSG view copies: ${_qtcad_quad_direct}")
    endif()
    foreach(_pat
	[[#[ \t]*include[ \t]*[<"]bsg/export\.h]]
	[[#[ \t]*include[ \t]*[<"]bsg/render\.h]]
	[[(^|[^A-Za-z0-9_])bsg_export_query([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])bsg_export_request([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])BSG_EXPORT_[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])BSG_RENDER_FLAG_[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _qtcad_quad_export_direct
	"${_qtcad_quad_view_contents}")
      if(_qtcad_quad_export_direct)
	_brlobol_pivot_guard_fail(
	  "src/libqtcad/QgQuadView.cpp reintroduced direct BSG export path discovery: ${_qtcad_quad_export_direct}")
      endif()
    endforeach()
  endif()

  foreach(_rel
      src/libqtcad/QgSW.cpp
      src/libqtcad/QgGL.cpp)
    set(_file "${BRLCAD_SOURCE_DIR}/${_rel}")
    if(NOT EXISTS "${_file}")
      _brlobol_pivot_guard_fail("${_rel} is required for qtcad framebuffer view adapter checks")
      continue()
    endif()
    file(READ "${_file}" _contents)
    foreach(_token
	"rt/view_legacy_bsg.h"
	"rt_view_width_from_bsg"
	"rt_view_height_from_bsg"
	"rt_view_dimensions_set_bsg"
	"rt_view_scale_storage_from_bsg")
      string(FIND "${_contents}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "${_rel} must route framebuffer size view reads through rt/view_legacy_bsg.h token ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[->[ \t\r\n]*gv_(width|height)[ \t\r\n]*=]]
	[[fb_configure_window[ \t\r\n]*\([^;]*->[ \t\r\n]*gv_(width|height)]]
	[[dm_set_vp[ \t\r\n]*\([^;]*gv_scale]])
      string(REGEX MATCH "${_pat}" _qtcad_fb_direct "${_contents}")
      if(_qtcad_fb_direct)
	_brlobol_pivot_guard_fail(
	  "${_rel} reintroduced direct framebuffer BSG view reads: ${_qtcad_fb_direct}")
      endif()
    endforeach()
  endforeach()

  set(_qtcad_sw "${BRLCAD_SOURCE_DIR}/src/libqtcad/QgSW.cpp")
  if(EXISTS "${_qtcad_sw}")
    file(READ "${_qtcad_sw}" _qtcad_sw_contents)
    string(FIND "${_qtcad_sw_contents}" "rt_view_model2view_from_bsg" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail(
	"src/libqtcad/QgSW.cpp must route legacy SW matrix loads through rt_view_model2view_from_bsg")
    endif()
    string(REGEX MATCH [[dm_loadmatrix[ \t\r\n]*\([^;]*->[ \t\r\n]*gv_model2view]]
      _qtcad_sw_direct "${_qtcad_sw_contents}")
    if(_qtcad_sw_direct)
      _brlobol_pivot_guard_fail(
	"src/libqtcad/QgSW.cpp reintroduced direct legacy SW model2view reads: ${_qtcad_sw_direct}")
    endif()
  endif()
endfunction()

function(_brlobol_pivot_guard_check_qged_view_info_rt_adapter)
  set(_view_info_model
    "${BRLCAD_SOURCE_DIR}/src/qged/plugins/view/info/CADViewModel.cpp")
  if(NOT EXISTS "${_view_info_model}")
    _brlobol_pivot_guard_fail(
      "src/qged/plugins/view/info/CADViewModel.cpp is required for qged view-info migration checks")
    return()
  endif()

  file(READ "${_view_info_model}" _view_info_contents)
  foreach(_token
      [[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
      [[rt_view_info_from_bsg]]
      [[rt_view_aet_from_bsg]]
      [[rt_view_center_from_bsg]])
    string(REGEX MATCH "${_token}" _view_info_adapter_hit
      "${_view_info_contents}")
    if(NOT _view_info_adapter_hit)
      _brlobol_pivot_guard_fail(
	"src/qged/plugins/view/info/CADViewModel.cpp must route view-info numeric reads through rt/view_legacy_bsg.h")
    endif()
  endforeach()

  foreach(_pat
      [[(^|[^A-Za-z0-9_])gv_size([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])gv_width([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])gv_height([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])gv_aet([^A-Za-z0-9_]|$)]]
      [[(^|[^A-Za-z0-9_])gv_center([^A-Za-z0-9_]|$)]])
    string(REGEX MATCH "${_pat}" _view_info_direct_hit
      "${_view_info_contents}")
    if(_view_info_direct_hit)
      _brlobol_pivot_guard_fail(
	"src/qged/plugins/view/info/CADViewModel.cpp reintroduced direct BSG view-info numeric reads: ${_view_info_direct_hit}")
    endif()
  endforeach()

  set(_qged_app "${BRLCAD_SOURCE_DIR}/src/qged/QgEdApp.cpp")
  if(EXISTS "${_qged_app}")
    file(READ "${_qged_app}" _qged_app_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_dimensions_set_bsg]])
      string(REGEX MATCH "${_token}" _qged_app_dim_adapter_hit
	"${_qged_app_contents}")
      if(NOT _qged_app_dim_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/qged/QgEdApp.cpp must route current-view dimension writes through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[view[ \t\r\n]*\([ \t\r\n]*\)->[ \t\r\n]*gv_width[ \t\r\n]*=]]
	[[view[ \t\r\n]*\([ \t\r\n]*\)->[ \t\r\n]*gv_height[ \t\r\n]*=]])
      string(REGEX MATCH "${_pat}" _qged_app_dim_direct_hit
	"${_qged_app_contents}")
      if(_qged_app_dim_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/qged/QgEdApp.cpp reintroduced direct current-view dimension writes: ${_qged_app_dim_direct_hit}")
      endif()
    endforeach()
  endif()
endfunction()

function(_brlobol_pivot_guard_check_brlobol_mesh_identity)
  foreach(_rel
      include/brlobol/pick_detail.h
      include/brlobol/export_action.h
      include/brlobol/measure_action.h
      include/brlobol/mesh_shape.h
      include/brlobol/snap_action.h
      src/libbrlobol/pick_detail.cpp
      src/libbrlobol/mesh_shape.cpp
      src/libbrlobol/export_action.cpp
      src/libbrlobol/measure_action.cpp
      src/libbrlobol/snap_action.cpp
      src/libbrlobol/tests/test_prototype.cpp)
    if(NOT EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _brlobol_pivot_guard_fail("${_rel} is required for Obol mesh vertex identity coverage")
      continue()
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/pick_detail.h" _pick_header)
  foreach(_token
      "setFaceVertexIndices"
      "getFaceVertexIndexA"
      "getFaceVertexIndexB"
      "getFaceVertexIndexC"
      "setNearestFaceEdge"
      "getNearestFaceEdgeSlot"
      "getNearestFaceEdgeVertexIndexA"
      "getNearestFaceEdgeVertexIndexB"
      "setNearestFaceVertex"
      "getNearestFaceVertexSlot"
      "getNearestFaceVertexIndex"
      "faceVertexIndex"
      "BRLObolSourceMeshPickResult"
      "BRLObolRtPickResult"
      "BRLObolRtPickCache"
      "prepare"
      "pickRay"
      "getObjectPathCount"
      "brlobol_pick_source_full_detail_result"
      "brlobol_pick_rt_ray"
      "IMPLICIT_SOLID"
      "SoBRLSourceMeshPickAction"
      "setRay"
      "getSourceBackedFullDetailRequestCount"
      "getSourceBackedFullDetailRequest"
      "makeSourceBackedFullDetailLodRequest"
      "submitSourceBackedFullDetailRequests"
      "consumeSourceBackedFullDetailResults"
      "rayIntersectsBounds"
      "clear")
    string(FIND "${_pick_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/pick_detail.h missing mesh vertex identity token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/pick_detail.cpp" _pick_impl)
  foreach(_token
      "faceVertexIndex[0] = -1"
      "other.faceVertexIndex"
      "SoBRLPickDetail::setFaceVertexIndices"
      "SoBRLPickDetail::getFaceVertexIndex"
      "nearestFaceEdgeSlot"
      "SoBRLPickDetail::setNearestFaceEdge"
      "SoBRLPickDetail::setNearestFaceVertex"
      "BRLObolSourceMeshPickResult::clear"
      "BRLObolRtPickResult::clear"
      "BRLObolRtPickCache::prepare"
      "BRLObolRtPickCache::pickRay"
      "BRLObolRtPickCache::clear"
      "pick_rt_same_object_paths"
      "brlobol_pick_rt_ray"
      "brlobol_pick_rt_hit"
      "rt_shootray"
      "RT_HIT_NORMAL"
      "pick_source_ray_triangle"
      "pick_source_fill_detail"
      "pick_source_mesh_face_index"
      "pick_source_mesh_vertex_index"
      "brlobol_pick_source_full_detail_result"
      "pick_source_full_detail_result_valid"
      "sourceRequest.queryRayValid"
      "SO_ACTION_SOURCE(SoBRLSourceMeshPickAction)"
      "SoBRLSourceMeshPickAction::initClass"
      "SoBRLSourceMeshPickAction::meshShapeAction"
      "SoBRLSourceMeshPickAction::appendSourceBackedFullDetailRequest"
      "SoBRLSourceMeshPickAction::rayIntersectsBounds"
      "pick_action_ray_intersects_box"
      "SoModelMatrixElement::get"
      "shape->getPickGeometryPolicy()"
      "shape->needsSourceBackedFullDetail()"
      "request.queryRayValid = 1"
      "request.queryRayOrigin"
      "request.queryRayDirection"
      "brlobol_lod_submit_rt_source_full_detail_request")
    string(FIND "${_pick_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/pick_detail.cpp missing mesh vertex identity token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/mesh_shape.h" _mesh_header)
  string(FIND "${_mesh_header}" "getTriangleVertexIndices" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("include/brlobol/mesh_shape.h missing getTriangleVertexIndices")
  endif()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/mesh_shape.cpp" _mesh_impl)
  foreach(_token
      "SoBRLMeshShape::getTriangleVertexIndices"
      "detail->setFaceVertexIndices"
      "mesh_nearest_face_edge_slot"
      "mesh_nearest_face_vertex_slot"
      "detail->setNearestFaceEdge"
      "detail->setNearestFaceVertex"
      "pp->getObjectPoint(this)"
      "faceDetail->getPoint(0)"
      "p0->getCoordinateIndex()")
    string(FIND "${_mesh_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/mesh_shape.cpp missing mesh vertex identity token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/export_action.h" _export_header)
  foreach(_token
      "vertexIndexA"
      "vertexIndexB"
      "vertexIndexC")
    string(FIND "${_export_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/export_action.h missing mesh vertex identity field ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/export_action.cpp" _export_impl)
  foreach(_token
      "shape->getTriangleVertexIndices"
      "record.vertexIndexA"
      "record.vertexIndexB"
      "record.vertexIndexC"
      "export_source_mesh_vertex_index")
    string(FIND "${_export_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/export_action.cpp missing mesh vertex identity token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/snap_action.h" _snap_header)
  foreach(_token
      "VERTEX"
      "EDGE_NEAREST"
      "getVertexIndex"
      "getEdgeSlot"
      "getEdgeVertexIndexA"
      "getEdgeVertexIndexB")
    string(FIND "${_snap_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/snap_action.h missing mesh sub-entity snap token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/snap_action.cpp" _snap_impl)
  foreach(_token
      "SoBRLSnapAction::getVertexIndex"
      "SoBRLSnapAction::getEdgeSlot"
      "SoBRLSnapAction::getEdgeVertexIndexA"
      "SoBRLSnapAction::getEdgeVertexIndexB"
      "snapAction->consider(VERTEX"
      "snapAction->consider(EDGE_NEAREST"
      "this->consider(VERTEX"
      "this->consider(EDGE_NEAREST"
      "shape->getTriangleVertexIndices"
      "shape->getFullDetailTriangleVertexIndices"
      "snap_source_mesh_vertex_index")
    string(FIND "${_snap_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/snap_action.cpp missing mesh sub-entity snap token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/measure_action.h" _measure_header)
  foreach(_token
      "getNearestFaceVertexIndexA"
      "getNearestFaceVertexIndexB"
      "getNearestFaceVertexIndexC"
      "getNearestFaceEdgeSlot"
      "getNearestFaceEdgeVertexIndexA"
      "getNearestFaceEdgeVertexIndexB"
      "getNearestFaceVertexSlot"
      "getNearestFaceVertexIndex")
    string(FIND "${_measure_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/measure_action.h missing mesh sub-entity measure token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/measure_action.cpp" _measure_impl)
  foreach(_token
      "SoBRLMeasureAction::getNearestFaceVertexIndexA"
      "SoBRLMeasureAction::getNearestFaceEdgeSlot"
      "nearest_face_edge_slot"
      "nearest_face_vertex_slot"
      "nearestFaceEdgeVertexIndex"
      "shape->getTriangleVertexIndices"
      "shape->getFullDetailTriangleVertexIndices"
      "this->measureTriangle(sourceRequest.path"
      "vertexIndices[edges[edgeSlot][0]]"
      "measure_source_mesh_vertex_index")
    string(FIND "${_measure_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/measure_action.cpp missing mesh sub-entity measure token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_prototype.cpp" _prototype_test)
  foreach(_token
      "getFaceVertexIndexA() != 0"
      "getFaceVertexIndexB() != 1"
      "getFaceVertexIndexC() != 2"
      "getNearestFaceEdgeSlot() != 0"
      "getNearestFaceEdgeVertexIndexA() != 0"
      "getNearestFaceEdgeVertexIndexB() != 1"
      "getNearestFaceVertexSlot() != 0"
      "getNearestFaceVertexIndex() != 0"
      "tri.vertexIndexA != 0"
      "tri.vertexIndexB != 1"
      "tri.vertexIndexC != 2"
      "mesh measure should report nearest transformed mesh face and sub-entity identity"
      "meshVertexSnap"
      "SoBRLSnapAction::VERTEX"
      "getVertexIndex() != 0"
      "meshEdgeSnap"
      "SoBRLSnapAction::EDGE_NEAREST"
      "getEdgeVertexIndexA() != 0"
      "getEdgeVertexIndexB() != 1")
    string(FIND "${_prototype_test}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/tests/test_prototype.cpp missing mesh vertex identity coverage token ${_token}")
    endif()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_brlobol_edit_intent_identity)
  foreach(_rel
      include/brlobol/export_action.h
      include/brlobol/measure_action.h
      include/brlobol/mesh_shape.h
      include/brlobol/pick_detail.h
      include/brlobol/snap_action.h
      include/brlobol/source_mesh_request.h
      include/brlobol/edit_preview.h
      include/brlobol/view_controller.h
      include/brlobol/vlist_shape.h
      src/libbrlobol/edit_preview.cpp
      src/libbrlobol/export_action.cpp
      src/libbrlobol/measure_action.cpp
      src/libbrlobol/mesh_shape.cpp
      src/libbrlobol/pick_detail.cpp
      src/libbrlobol/snap_action.cpp
      src/libbrlobol/vlist_shape.cpp
      src/libbrlobol/tests/test_prototype.cpp)
    if(NOT EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _brlobol_pivot_guard_fail("${_rel} is required for Obol edit-intent identity coverage")
      continue()
    endif()
  endforeach()

  foreach(_rel
      include/brlobol/mesh_shape.h
      include/brlobol/vlist_shape.h)
    file(READ "${BRLCAD_SOURCE_DIR}/${_rel}" _shape_header)
    foreach(_token
	"editIntentId"
	"editIntentRole")
      string(FIND "${_shape_header}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("${_rel} missing edit-intent field ${_token}")
      endif()
    endforeach()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/source_mesh_request.h" _source_request_header)
  foreach(_token
      "SbString editIntentId"
      "SbString editIntentRole"
      "editIntentId = \"\""
      "editIntentRole = \"\"")
    string(FIND "${_source_request_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/source_mesh_request.h missing edit-intent request token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/pick_detail.h" _pick_header)
  foreach(_token
      "setEditIntent"
      "getEditIntentId"
      "getEditIntentRole")
    string(FIND "${_pick_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/pick_detail.h missing edit-intent detail token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/export_action.h" _export_header)
  foreach(_token
      "editIntentId"
      "editIntentRole")
    string(FIND "${_export_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/export_action.h missing edit-intent export token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/snap_action.h" _snap_header)
  foreach(_token
      "getEditIntentId"
      "getEditIntentRole"
      "candidateEditIntentId"
      "candidateEditIntentRole")
    string(FIND "${_snap_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/snap_action.h missing edit-intent snap token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/measure_action.h" _measure_header)
  foreach(_token
      "getNearestEditIntentId"
      "getNearestEditIntentRole"
      "nearestEditIntentId"
      "nearestEditIntentRole")
    string(FIND "${_measure_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/measure_action.h missing edit-intent measure token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/edit_preview.h" _edit_preview_header)
  foreach(_token
      "SoSFString editIntentId"
      "SoSFString editIntentRole"
      "setEditIntent")
    string(FIND "${_edit_preview_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/edit_preview.h missing edit-preview intent token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/view_controller.h" _view_controller_header)
  string(FIND "${_view_controller_header}" "replaceEditPreviewWithIntent" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("include/brlobol/view_controller.h missing explicit edit-preview intent API")
  endif()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/pick_detail.cpp" _pick_impl)
  foreach(_token
      "SoBRLPickDetail::setEditIntent"
      "SoBRLPickDetail::getEditIntentId"
      "SoBRLPickDetail::getEditIntentRole"
      "sourceRequest.editIntentId")
    string(FIND "${_pick_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/pick_detail.cpp missing edit-intent detail implementation token ${_token}")
    endif()
  endforeach()

  foreach(_rel
      src/libbrlobol/vlist_shape.cpp
      src/libbrlobol/mesh_shape.cpp)
    file(READ "${BRLCAD_SOURCE_DIR}/${_rel}" _shape_impl)
    foreach(_token
	"SO_NODE_ADD_FIELD(editIntentId"
	"SO_NODE_ADD_FIELD(editIntentRole"
	"setEditIntent(this->editIntentId")
      string(FIND "${_shape_impl}" "${_token}" _idx)
      if(_idx EQUAL -1)
	_brlobol_pivot_guard_fail("${_rel} missing edit-intent shape implementation token ${_token}")
      endif()
    endforeach()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/export_action.cpp" _export_impl)
  foreach(_token
      "record.editIntentId = editIntentId"
      "record.editIntentRole = editIntentRole"
      "shape->editIntentId.getValue()"
      "sourceRequest.editIntentId")
    string(FIND "${_export_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/export_action.cpp missing edit-intent export token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/snap_action.cpp" _snap_impl)
  foreach(_token
      "SoBRLSnapAction::getEditIntentId"
      "candidateEditIntentId = editIntentId"
      "shape->editIntentId.getValue()"
      "sourceRequest.editIntentId")
    string(FIND "${_snap_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/snap_action.cpp missing edit-intent snap token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/measure_action.cpp" _measure_impl)
  foreach(_token
      "SoBRLMeasureAction::getNearestEditIntentId"
      "nearestEditIntentId = editIntentId"
      "shape->editIntentId.getValue()"
      "sourceRequest.editIntentId")
    string(FIND "${_measure_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/measure_action.cpp missing edit-intent measure token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/edit_preview.cpp" _edit_preview_impl)
  foreach(_token
      "SO_NODE_ADD_FIELD(editIntentId"
      "SO_NODE_ADD_FIELD(editIntentRole"
      "SoBRLEditPreview::setEditIntent"
      "shape->editIntentId = intentId.getLength()"
      "shape->editIntentRole = intentRole.getLength()")
    string(FIND "${_edit_preview_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/edit_preview.cpp missing edit-preview intent token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/view_controller.cpp" _view_controller_impl)
  foreach(_token
      "BRLObolViewController::replaceEditPreviewWithIntent"
      "preview->setEditIntent")
    string(FIND "${_view_controller_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/view_controller.cpp missing explicit edit-preview intent token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_prototype.cpp" _prototype_test)
  foreach(_token
      "edit preview should publish typed edit-intent metadata"
      "edit preview should publish explicit live edit intent fields"
      "edit preview snap should expose explicit live edit intent"
      "edit preview measure should expose explicit live edit intent"
      "edit preview export should expose explicit live edit intent"
      "getEditIntentId()"
      "edit preview pick should preserve preview and edit-intent identity"
      "edit preview snap should preserve edit-intent metadata"
      "getNearestEditIntentId()"
      "edit preview measure should preserve edit-intent metadata"
      "edit preview export should preserve edit-intent metadata")
    string(FIND "${_prototype_test}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/tests/test_prototype.cpp missing edit-intent coverage token ${_token}")
    endif()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_brlobol_lod_metadata)
  foreach(_rel
      include/brlobol.h
      include/brlobol/CMakeLists.txt
      include/brlobol/database_source.h
      include/brlobol/lod_mesh_shape.h
      include/brlobol/lod_realization.h
      include/brlobol/lod_service.h
      include/brlobol/lod_update_action.h
      include/brlobol/mesh_lod_submit_action.h
      include/brlobol/mesh_residency_action.h
      include/brlobol/export_action.h
      include/brlobol/measure_action.h
      include/brlobol/mesh_shape.h
      include/brlobol/snap_action.h
      src/libbrlobol/CMakeLists.txt
      src/libbrlobol/init.cpp
      src/libbrlobol/lod_mesh_shape.cpp
      src/libbrlobol/lod_realization.cpp
      src/libbrlobol/lod_service.cpp
      src/libbrlobol/lod_update_action.cpp
      src/libbrlobol/mesh_lod_submit_action.cpp
      src/libbrlobol/mesh_residency_action.cpp
      src/libbrlobol/export_action.cpp
      src/libbrlobol/measure_action.cpp
      src/libbrlobol/mesh_shape.cpp
      src/libbrlobol/snap_action.cpp
      src/libbrlobol/database_source.cpp
      src/libbrlobol/tests/CMakeLists.txt
      src/libbrlobol/tests/test_lod_realization.cpp
      src/libbrlobol/tests/test_lod_service.cpp
      src/libbrlobol/tests/test_lod_update_action.cpp
      src/libbrlobol/tests/test_prototype.cpp)
    if(NOT EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _brlobol_pivot_guard_fail("${_rel} is required for Obol mesh LoD metadata coverage")
      continue()
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol.h" _brlobol_header)
  string(FIND "${_brlobol_header}" "brlobol/lod_realization.h" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("include/brlobol.h must include brlobol/lod_realization.h")
  endif()
  string(FIND "${_brlobol_header}" "brlobol/lod_service.h" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("include/brlobol.h must include brlobol/lod_service.h")
  endif()
  string(FIND "${_brlobol_header}" "brlobol/lod_mesh_shape.h" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("include/brlobol.h must include brlobol/lod_mesh_shape.h")
  endif()
  string(FIND "${_brlobol_header}" "brlobol/lod_update_action.h" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("include/brlobol.h must include brlobol/lod_update_action.h")
  endif()
  string(FIND "${_brlobol_header}" "brlobol/mesh_lod_submit_action.h" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("include/brlobol.h must include brlobol/mesh_lod_submit_action.h")
  endif()
  string(FIND "${_brlobol_header}" "brlobol/mesh_residency_action.h" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("include/brlobol.h must include brlobol/mesh_residency_action.h")
  endif()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/CMakeLists.txt" _brlobol_include_cmake)
  string(FIND "${_brlobol_include_cmake}" "lod_realization.h" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("include/brlobol/CMakeLists.txt must install lod_realization.h")
  endif()
  string(FIND "${_brlobol_include_cmake}" "lod_service.h" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("include/brlobol/CMakeLists.txt must install lod_service.h")
  endif()
  string(FIND "${_brlobol_include_cmake}" "lod_mesh_shape.h" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("include/brlobol/CMakeLists.txt must install lod_mesh_shape.h")
  endif()
  string(FIND "${_brlobol_include_cmake}" "lod_update_action.h" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("include/brlobol/CMakeLists.txt must install lod_update_action.h")
  endif()
  string(FIND "${_brlobol_include_cmake}" "mesh_lod_submit_action.h" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("include/brlobol/CMakeLists.txt must install mesh_lod_submit_action.h")
  endif()
  string(FIND "${_brlobol_include_cmake}" "mesh_residency_action.h" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("include/brlobol/CMakeLists.txt must install mesh_residency_action.h")
  endif()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/CMakeLists.txt" _brlobol_cmake)
  string(FIND "${_brlobol_cmake}" "lod_realization.cpp" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("src/libbrlobol/CMakeLists.txt must build lod_realization.cpp")
  endif()
  foreach(_token
      "lod_service.cpp"
      "lod_mesh_shape.cpp"
      "lod_update_action.cpp"
      "mesh_lod_submit_action.cpp"
      "mesh_residency_action.cpp"
      "Threads::Threads")
    string(FIND "${_brlobol_cmake}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/CMakeLists.txt must build/link staged LoD service token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/init.cpp" _brlobol_init_impl)
  foreach(_token
      "brlobol/lod_mesh_shape.h"
      "brlobol/lod_update_action.h"
      "brlobol/mesh_lod_submit_action.h"
      "brlobol/mesh_residency_action.h"
      "SoBRLLodMeshShape::initClass"
      "SoBRLSourceMeshPickAction::initClass"
      "SoBRLLodUpdateAction::initClass"
      "SoBRLMeshLodSubmitAction::initClass"
      "SoBRLMeshResidencyAction::initClass")
    string(FIND "${_brlobol_init_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/init.cpp must initialize staged LoD update action token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/lod_mesh_shape.h" _lod_mesh_header)
  foreach(_token
      "SoBRLLodMeshShape"
      "SoBRLMeshShape"
      "SO_NODE_HEADER(SoBRLLodMeshShape)"
      "initClass")
    string(FIND "${_lod_mesh_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/lod_mesh_shape.h missing LoD-backed mesh token ${_token}")
    endif()
  endforeach()
  foreach(_pat
      [[#[ \t]*include[ \t]*[<"]bsg/]]
      [[(^|[^A-Za-z0-9_])bsg_mesh_lod([^A-Za-z0-9_]|$)]])
    string(REGEX MATCH "${_pat}" _hit "${_lod_mesh_header}")
    if(_hit)
      _brlobol_pivot_guard_fail("include/brlobol/lod_mesh_shape.h must not expose BSG LoD APIs: ${_hit}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/lod_mesh_shape.cpp" _lod_mesh_impl)
  foreach(_token
      "SO_NODE_SOURCE(SoBRLLodMeshShape)"
      "SO_NODE_CONSTRUCTOR(SoBRLLodMeshShape)"
      "setLodBackedMesh(TRUE)"
      "setLodPreserveFullDetail(FALSE)"
      "SO_NODE_INIT_CLASS(SoBRLLodMeshShape, SoBRLMeshShape")
    string(FIND "${_lod_mesh_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/lod_mesh_shape.cpp missing LoD-backed mesh token ${_token}")
    endif()
  endforeach()
  foreach(_pat
      [[#[ \t]*include[ \t]*[<"]bsg/]]
      [[(^|[^A-Za-z0-9_])bsg_mesh_lod([^A-Za-z0-9_]|$)]])
    string(REGEX MATCH "${_pat}" _hit "${_lod_mesh_impl}")
    if(_hit)
      _brlobol_pivot_guard_fail("src/libbrlobol/lod_mesh_shape.cpp must not use BSG LoD APIs: ${_hit}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/lod_realization.h" _lod_header)
  foreach(_token
      "BRLObolLodRequest"
      "BRLObolLodResult"
      "BRLObolLodGeometryHandle"
      "BRLObolLodMeshPayload"
      "BRLObolLodCacheKey"
      "BRLObolLodProviderParam"
      "BRLObolLodDependency"
      "BRLObolLodAttribute"
      "BRLObolLodProxy"
      "BRLOBOL_LOD_QUALITY_FAST_DISPLAY"
      "BRLOBOL_LOD_GEOMETRY_RT_MESH_CACHE"
      "BRLOBOL_LOD_PROXY_OBB"
      "databaseRevision"
      "sourceRevision"
      "sourceContentHash"
      "objectPath"
      "viewRevision"
      "policyRevision"
      "providerId"
      "providerVersion"
      "qualityTier"
      "estimatedError"
      "dependencies"
      "attributes"
      "proxy"
      "mesh"
      "normals"
      "faceIndex"
      "vertexIndex"
      "terminal"
      "fallback"
      "diagnostic"
      "brlobol_lod_cache_key"
      "brlobol_lod_mesh_payload_from_rt_mesh_data"
      "brlobol_lod_result_matches_request"
      "brlobol_lod_result_from_rt_mesh_info"
      "brlobol_lod_directory_result"
      "brlobol_lod_attributes_result"
      "brlobol_lod_aabb_result"
      "brlobol_lod_proxy_result")
    string(FIND "${_lod_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/lod_realization.h missing staged LoD contract token ${_token}")
    endif()
  endforeach()
  foreach(_pat
      [[#[ \t]*include[ \t]*[<"]bsg/]]
      [[(^|[^A-Za-z0-9_])bsg_mesh_lod([^A-Za-z0-9_]|$)]])
    string(REGEX MATCH "${_pat}" _hit "${_lod_header}")
    if(_hit)
      _brlobol_pivot_guard_fail("include/brlobol/lod_realization.h must not expose BSG LoD APIs: ${_hit}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/lod_realization.cpp" _lod_impl)
  foreach(_token
      "brlobol-lod-v1"
      "database_revision"
      "source_revision"
      "source_content_hash"
      "object_path"
      "view_revision"
      "policy_revision"
      "provider_id"
      "provider_version"
      "provider_param_count"
      "BRLOBOL_LOD_PROVIDER_STALE"
      "BRLOBOL_LOD_PROVIDER_CACHE_MISS"
      "BRLObolLodMeshPayload::isValid"
      "normals.clear"
      "normals.empty"
      "normals.size() == coordIndex.size()"
      "data.normal_count != index_count"
      "faceIndex.clear"
      "faceIndex.empty"
      "vertexIndex.clear"
      "vertexIndex.empty"
      "brlobol_lod_mesh_payload_from_rt_mesh_data"
      "BRLOBOL_LOD_RESULT_DIRECTORY"
      "BRLOBOL_LOD_RESULT_ATTRIBUTES"
      "BRLOBOL_LOD_RESULT_AABB"
      "BRLOBOL_LOD_RESULT_PROXY"
      "BRLOBOL_LOD_PROXY_AABB"
      "BRLOBOL_LOD_PROXY_OBB"
      "brlobol_lod_directory_result"
      "brlobol_lod_attributes_result"
      "brlobol_lod_aabb_result"
      "brlobol_lod_proxy_result")
    string(FIND "${_lod_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/lod_realization.cpp missing staged LoD implementation token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/lod_service.h" _lod_service_header)
  foreach(_token
      "BRLObolLodTaskProc"
      "BRLObolLodTaskDataFreeProc"
      "BRLObolLodCacheWriteProc"
      "BRLObolLodSubscriberId"
      "BRLObolLodResultReadyCB"
      "BRLObolRtMeshLodProvider"
      "brlobol_rt_mesh_lod_provider_task"
      "brlobol_rt_mesh_lod_provider_free"
      "BRLObolRtSourceFullDetailProvider"
      "brlobol_rt_source_full_detail_provider_task"
      "brlobol_rt_source_full_detail_provider_free"
      "brlobol_lod_rt_source_full_detail_request_from_source_mesh_request"
      "brlobol_lod_submit_rt_source_full_detail_request"
      "maxFullDetailFaceCount"
      "BRLObolLodTask"
      "BRLObolLodService"
      "start"
      "stop"
      "beginGeneration"
      "cancelGeneration"
      "isGenerationCancelled"
      "submit"
      "submitIfNotActive"
      "hasActiveRequest"
      "drainResults"
      "drainMatchingResults"
      "subscribeResultReady"
      "unsubscribeResultReady"
      "inFlightCount"
      "realizeDataFree"
      "useForcedLevel"
      "shrinkAfterCopy"
      "delayedTaskCountForDiagnostics"
      "queuedCacheWriteCountForDiagnostics")
    string(FIND "${_lod_service_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/lod_service.h missing async LoD service token ${_token}")
    endif()
  endforeach()
  foreach(_pat
      [[#[ \t]*include[ \t]*[<"]bsg/]]
      [[(^|[^A-Za-z0-9_])bsg_mesh_lod([^A-Za-z0-9_]|$)]])
    string(REGEX MATCH "${_pat}" _hit "${_lod_service_header}")
    if(_hit)
      _brlobol_pivot_guard_fail("include/brlobol/lod_service.h must not expose BSG LoD APIs: ${_hit}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/lod_service.cpp" _lod_service_impl)
  foreach(_token
      "std::condition_variable"
      "std::thread"
      "lod_worker_loop"
      "lod_cache_writer_loop"
      "lod_execute_task"
      "lod_task_free_realize_data"
      "BRLObolLodSubscriber"
      "lod_notify_result_ready"
      "subscriberCv"
      "lod_normalize_result"
      "BRLObolLodService::drainMatchingResults"
      "BRLObolLodService::submitIfNotActive"
      "BRLObolLodService::hasActiveRequest"
      "activeRequestKeys"
      "lod_active_request_key_recorded_unlocked"
      "BRLOBOL_LOD_PROVIDER_CANCELLED"
      "BRLOBOL_LOD_PROVIDER_STALE"
      "brlobol_rt_mesh_lod_provider_task"
      "brlobol_rt_mesh_lod_provider_free"
      "brlobol_rt_source_full_detail_provider_task"
      "brlobol_rt_source_full_detail_provider_free"
      "brlobol_lod_rt_source_full_detail_request_from_source_mesh_request"
      "brlobol_lod_submit_rt_source_full_detail_request"
      "db_mesh_lod_status"
      "db_mesh_lod_refresh"
      "rt_db_get_internal"
      "ID_BOT"
      "RT_BOT_CK_MAGIC"
      "BRLOBOL_LOD_RESULT_FULL_DETAIL"
      "source_query.bounds"
      "source_query.ray.origin"
      "source_query.ray.direction"
      "source_query.tolerance"
      "lod_remove_source_query_provider_params"
      "lod_provider_param_has_no_trailing_tokens"
      "lod_request_query_space_is_source_local"
      "lod_request_snap_query_bounds"
      "lod_request_pick_query_ray"
      "lod_request_has_scoped_subset_query"
      "lod_parse_bounds_provider_param"
      "lod_ray_intersects_triangle"
      "lod_bounds_intersect"
      "selectedFaces"
      "RT source full-detail provider scoped query matched no faces"
      "result.mesh.faceIndex"
      "result.mesh.vertexIndex"
      "sourceToLocal"
      "rt_mesh_lod_has_active_data"
      "rt_mesh_lod_load_level"
      "rt_mesh_lod_info_get"
      "rt_mesh_lod_data_get"
      "rt_mesh_lod_memshrink"
      "brlobol_lod_mesh_payload_from_rt_mesh_data"
      "cacheWriteInFlight"
      "debugDelayMilliseconds"
      "completed.insert")
    string(FIND "${_lod_service_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/lod_service.cpp missing async LoD service token ${_token}")
    endif()
  endforeach()
  foreach(_pat
      [[#[ \t]*include[ \t]*[<"]bsg/]]
      [[(^|[^A-Za-z0-9_])bsg_mesh_lod([^A-Za-z0-9_]|$)]])
    string(REGEX MATCH "${_pat}" _hit "${_lod_service_impl}")
    if(_hit)
      _brlobol_pivot_guard_fail("src/libbrlobol/lod_service.cpp must not use BSG LoD APIs: ${_hit}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/lod_update_action.h" _lod_update_header)
  foreach(_token
      "SoBRLLodUpdateAction"
      "drainService"
      "setResults"
      "addResult"
      "getMatchedResultCount"
      "getAppliedResultCount"
      "getRejectedResultCount"
      "getUnmatchedResultCount"
      "getDiagnostics")
    string(FIND "${_lod_update_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/lod_update_action.h missing staged LoD update action token ${_token}")
    endif()
  endforeach()
  foreach(_pat
      [[#[ \t]*include[ \t]*[<"]bsg/]]
      [[(^|[^A-Za-z0-9_])bsg_mesh_lod([^A-Za-z0-9_]|$)]])
    string(REGEX MATCH "${_pat}" _hit "${_lod_update_header}")
    if(_hit)
      _brlobol_pivot_guard_fail("include/brlobol/lod_update_action.h must not expose BSG LoD APIs: ${_hit}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/lod_update_action.cpp" _lod_update_impl)
  foreach(_token
      "SO_ACTION_SOURCE(SoBRLLodUpdateAction)"
      "BRLObolLodService"
      "service.drainResults"
      "SoBRLMeshShape"
      "applyStagedLodResult"
      "lod_update_path_matches"
      "matchedResultCount"
      "appliedResultCount"
      "rejectedResultCount"
      "unmatchedResultCount"
      "staged LoD result did not match any mesh")
    string(FIND "${_lod_update_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/lod_update_action.cpp missing staged LoD update action token ${_token}")
    endif()
  endforeach()
  foreach(_pat
      [[#[ \t]*include[ \t]*[<"]bsg/]]
      [[(^|[^A-Za-z0-9_])bsg_mesh_lod([^A-Za-z0-9_]|$)]])
    string(REGEX MATCH "${_pat}" _hit "${_lod_update_impl}")
    if(_hit)
      _brlobol_pivot_guard_fail("src/libbrlobol/lod_update_action.cpp must not use BSG LoD APIs: ${_hit}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/mesh_residency_action.h" _mesh_residency_header)
  foreach(_token
      "SoBRLMeshResidencyAction"
      "setMaxResidentMeshBytes"
      "setEvictDisplayPayloads"
      "getInitialResidentMeshBytes"
      "getFinalResidentMeshBytes"
      "getFreedFullDetailBytes"
      "getFreedDisplayBytes"
      "getEvictedFullDetailMeshCount"
      "getEvictedDisplayMeshCount")
    string(FIND "${_mesh_residency_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/mesh_residency_action.h missing mesh residency budget token ${_token}")
    endif()
  endforeach()
  foreach(_pat
      [[#[ \t]*include[ \t]*[<"]bsg/]]
      [[(^|[^A-Za-z0-9_])bsg_mesh_lod([^A-Za-z0-9_]|$)]])
    string(REGEX MATCH "${_pat}" _hit "${_mesh_residency_header}")
    if(_hit)
      _brlobol_pivot_guard_fail("include/brlobol/mesh_residency_action.h must not expose BSG LoD APIs: ${_hit}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/mesh_residency_action.cpp" _mesh_residency_impl)
  foreach(_token
      "SO_ACTION_SOURCE(SoBRLMeshResidencyAction)"
      "SoBRLMeshShape"
      "estimateResidentMeshBytes"
      "estimateFullDetailMeshBytes"
      "estimateDisplayMeshBytes"
      "evictFullDetailMesh"
      "evictActiveDisplayMesh"
      "needsSourceBackedFullDetail"
      "maxResidentMeshBytes"
      "freedDisplayBytes"
      "skippedDisplayMeshCount")
    string(FIND "${_mesh_residency_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/mesh_residency_action.cpp missing mesh residency budget token ${_token}")
    endif()
  endforeach()
  foreach(_pat
      [[#[ \t]*include[ \t]*[<"]bsg/]]
      [[(^|[^A-Za-z0-9_])bsg_mesh_lod([^A-Za-z0-9_]|$)]])
    string(REGEX MATCH "${_pat}" _hit "${_mesh_residency_impl}")
    if(_hit)
      _brlobol_pivot_guard_fail("src/libbrlobol/mesh_residency_action.cpp must not use BSG LoD APIs: ${_hit}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/mesh_lod_submit_action.h" _lod_submit_header)
  foreach(_token
      "SoBRLMeshLodSubmitAction"
      "setService"
      "setDatabase"
      "setViewInfo"
      "setGeneration"
      "setRevisions"
      "getSubmittedTaskCount"
      "getSkippedMeshCount"
      "setRequireLodBacked"
      "getRequireLodBacked"
      "setForcedLevel"
      "clearForcedLevel"
      "hasForcedLevel"
      "getForcedLevel"
      "getDiagnostics")
    string(FIND "${_lod_submit_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/mesh_lod_submit_action.h missing view-driven LoD submit action token ${_token}")
    endif()
  endforeach()
  foreach(_pat
      [[#[ \t]*include[ \t]*[<"]bsg/]]
      [[(^|[^A-Za-z0-9_])bsg_mesh_lod([^A-Za-z0-9_]|$)]])
    string(REGEX MATCH "${_pat}" _hit "${_lod_submit_header}")
    if(_hit)
      _brlobol_pivot_guard_fail("include/brlobol/mesh_lod_submit_action.h must not expose BSG LoD APIs: ${_hit}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/mesh_lod_submit_action.cpp" _lod_submit_impl)
  foreach(_token
      "SO_ACTION_SOURCE(SoBRLMeshLodSubmitAction)"
      "BRLObolLodService"
      "SoBRLMeshShape"
      "makeLodRequest"
      "brlobol_lod_cache_key"
      "current LoD request is already resident"
      "BRLObolRtMeshLodProvider"
      "brlobol_rt_mesh_lod_provider_task"
      "brlobol_rt_mesh_lod_provider_free"
      "realizeDataFree"
      "service->submit"
      "service->hasActiveRequest"
      "service->submitIfNotActive"
      "current LoD request is already active"
      "isLodBackedMesh"
      "provider->useForcedLevel"
      "provider->forcedLevel"
      "mesh is not LoD-backed"
      "BRLOBOL_LOD_DRAW_SHADED")
    string(FIND "${_lod_submit_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/mesh_lod_submit_action.cpp missing view-driven LoD submit action token ${_token}")
    endif()
  endforeach()
  foreach(_pat
      [[#[ \t]*include[ \t]*[<"]bsg/]]
      [[(^|[^A-Za-z0-9_])bsg_mesh_lod([^A-Za-z0-9_]|$)]])
    string(REGEX MATCH "${_pat}" _hit "${_lod_submit_impl}")
    if(_hit)
      _brlobol_pivot_guard_fail("src/libbrlobol/mesh_lod_submit_action.cpp must not use BSG LoD APIs: ${_hit}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/export_action.h" _export_header)
  foreach(_token
      "enum GeometryPolicy"
      "FULL_DETAIL"
      "DISPLAY_LEVEL"
      "setGeometryPolicy"
      "getGeometryPolicy"
      "getSkippedLodDisplayMeshCount"
      "getSourceBackedFullDetailRequestCount"
      "makeSourceBackedFullDetailLodRequest"
      "appendSourceBackedFullDetailResult"
      "submitSourceBackedFullDetailRequests"
      "consumeSourceBackedFullDetailResults"
      "struct TriangleRecord"
      "lodAvailable"
      "lodActiveLevel"
      "lodFaceCount"
      "lodPointCount"
      "lodOriginalPointCount"
      "lodBoundsMin"
      "lodBoundsMax")
    string(FIND "${_export_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/export_action.h missing LoD metadata field ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/mesh_shape.h" _mesh_header)
  foreach(_token
      "lodAvailable"
      "lodActiveLevel"
      "lodFaceCount"
      "lodPointCount"
      "lodOriginalPointCount"
      "lodBoundsMin"
      "lodBoundsMax"
      "lodStagedAvailable"
      "lodResultKind"
      "lodQualityTier"
      "lodProviderStatus"
      "lodCacheKey"
      "lodDependencyPath"
      "lodAttributeName"
      "lodProxyKind"
      "lodProxyHalfExtents"
      "lodBacked"
      "lodDisplayActive"
      "lodPreserveFullDetail"
      "PickGeometryPolicy"
      "PICK_DISPLAY_LEVEL"
      "PICK_FULL_DETAIL"
      "pickGeometryPolicy"
      "setLodBackedMesh"
      "isLodBackedMesh"
      "setLodPreserveFullDetail"
      "isLodPreserveFullDetailEnabled"
      "setPickGeometryPolicy"
      "getPickGeometryPolicy"
      "isLodDisplayActive"
      "hasFullDetailMesh"
      "getFullDetailTriangleCount"
      "getFullDetailTriangle"
      "makeLodRequest"
      "makeSourceMeshRequest"
      "needsSourceBackedFullDetail"
      "estimateDisplayMeshBytes"
      "estimateResidentMeshBytes"
      "evictFullDetailMesh"
      "evictActiveDisplayMesh"
      "applyStagedLodResult"
      "clearStagedLodResult")
    string(FIND "${_mesh_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/mesh_shape.h missing LoD metadata field ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/mesh_shape.cpp" _mesh_impl)
  foreach(_token
      "SoBRLMeshShape::applyStagedLodResult"
      "SoBRLMeshShape::makeLodRequest"
      "brlobol_lod_result_matches_request"
      "BRLOBOL_LOD_PROVIDER_STALE"
      "BRLOBOL_LOD_PROVIDER_CANCELLED"
      "BRLOBOL_LOD_PROVIDER_CACHE_MISS"
      "BRLOBOL_LOD_PROVIDER_ERROR"
      "lodDependencySourceRevision"
      "lodAttributeValue"
      "lodProxyAxisX"
      "lodProxyRadius"
      "SoBRLMeshShape::setLodBackedMesh"
      "SoBRLMeshShape::isLodBackedMesh"
      "SoBRLMeshShape::setLodPreserveFullDetail"
      "SoBRLMeshShape::isLodPreserveFullDetailEnabled"
      "SoBRLMeshShape::setPickGeometryPolicy"
      "SoBRLMeshShape::getPickGeometryPolicy"
      "SoBRLMeshShape::isLodDisplayActive"
      "SoBRLMeshShape::makeSourceMeshRequest"
      "SoBRLMeshShape::needsSourceBackedFullDetail"
      "SoBRLMeshShape::estimateDisplayMeshBytes"
      "SoBRLMeshShape::evictActiveDisplayMesh"
      "SoBRLMeshShape::updateSourceMeshMetricsFromFields"
      "SoBRLMeshShape::captureFullDetailMesh"
      "SoBRLMeshShape::restoreFullDetailMesh"
      "SoBRLMeshShape::getFullDetailTriangle"
      "getFullDetailTriangleVertexIndices"
      "fullDetailPoint"
      "fullDetailCoordIndex"
      "sourceMeshMetricsValid"
      "sourceMeshBounds"
      "BRLObolSourceMeshRequest"
      "lodDisplayActive"
      "lodPreserveFullDetail"
      "pickGeometryPolicy"
      "BRLOBOL_LOD_RESULT_MESH"
      "result.mesh.isValid"
      "setIndexedTriangles"
      "mesh_lod_uint64_string")
    string(FIND "${_mesh_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/mesh_shape.cpp missing staged LoD mesh-consumption token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/source_mesh_request.h" _source_mesh_request_header)
  foreach(_token
      "BRLObolSourceMeshRequest"
      "sourceName"
      "sourceType"
      "sourceId"
      "faceCount"
      "pointCount"
      "bounds"
      "localToWorld"
      "queryBoundsValid"
      "queryBounds"
      "queryRayValid"
      "queryRayOrigin"
      "queryRayDirection"
      "queryToleranceValid"
      "queryTolerance"
      "materialColor"
      "selected"
      "lodActiveLevel"
      "colorOverride")
    string(FIND "${_source_mesh_request_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/source_mesh_request.h missing source-backed exact request token ${_token}")
    endif()
  endforeach()
  foreach(_pat
      [[#[ \t]*include[ \t]*[<"]bsg/]]
      [[(^|[^A-Za-z0-9_])bsg_mesh_lod([^A-Za-z0-9_]|$)]])
    string(REGEX MATCH "${_pat}" _hit "${_mesh_impl}")
    if(_hit)
      _brlobol_pivot_guard_fail("src/libbrlobol/mesh_shape.cpp must not use BSG LoD APIs: ${_hit}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/database_source.cpp" _db_source_impl)
  foreach(_token
      "publish_lod_metadata_if_cached"
      "lodBotThreshold"
      "BRLObolLodRequest"
      "BRLObolLodResult"
      "brlobol_lod_result_from_rt_mesh_info"
      "db_mesh_lod_status"
      "db_mesh_lod_get"
      "rt_mesh_lod_load_view"
      "rt_mesh_lod_info_get"
      "SoBRLLodMeshShape"
      "applyStagedLodResult")
    string(FIND "${_db_source_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/database_source.cpp missing LoD metadata token ${_token}")
    endif()
  endforeach()
  foreach(_pat
      [[#[ \t]*include[ \t]*[<"]bsg/]]
      [[(^|[^A-Za-z0-9_])bsg_mesh_lod([^A-Za-z0-9_]|$)]])
    string(REGEX MATCH "${_pat}" _hit "${_db_source_impl}")
    if(_hit)
      _brlobol_pivot_guard_fail("src/libbrlobol/database_source.cpp must use RT LoD metadata, not BSG LoD APIs: ${_hit}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/database_source.h" _db_source_header)
  foreach(_token
      "configureDatabaseSource"
      "lodBotThreshold"
      "lodBotThresholdSensor")
    string(FIND "${_db_source_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/database_source.h missing LoD threshold token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/export_action.cpp" _export_impl)
  foreach(_token
      "SoBRLExportAction::FULL_DETAIL"
      "SoBRLExportAction::DISPLAY_LEVEL"
      "shape->isLodDisplayActive()"
      "shape->hasFullDetailMesh()"
      "shape->needsSourceBackedFullDetail()"
      "shape->getFullDetailTriangle"
      "appendSourceBackedFullDetailRequest"
      "getSourceBackedFullDetailRequestCount"
      "makeSourceBackedFullDetailLodRequest"
      "appendSourceBackedFullDetailResult"
      "submitSourceBackedFullDetailRequests"
      "consumeSourceBackedFullDetailResults"
      "skippedLodDisplayMeshCount"
      "export_source_full_detail_result_valid"
      "export_source_mesh_face_index"
      "brlobol_lod_rt_source_full_detail_request_from_source_mesh_request"
      "brlobol_lod_submit_rt_source_full_detail_request"
      "request.localToWorld = localToWorld"
      "shape->lodAvailable.getValue()"
      "shape->lodActiveLevel.getValue()"
      "record.lodFaceCount"
      "record.lodBoundsMin"
      "record.lodBoundsMax")
    string(FIND "${_export_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/export_action.cpp missing LoD metadata token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/measure_action.h" _measure_header)
  foreach(_token
      "enum GeometryPolicy"
      "FULL_DETAIL"
      "DISPLAY_LEVEL"
      "setGeometryPolicy"
      "getGeometryPolicy"
      "getSkippedLodDisplayMeshCount"
      "getSourceBackedFullDetailRequestCount"
      "makeSourceBackedFullDetailLodRequest"
      "consumeSourceBackedFullDetailResult"
      "submitSourceBackedFullDetailRequests"
      "consumeSourceBackedFullDetailResults")
    string(FIND "${_measure_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/measure_action.h missing LoD geometry policy token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/measure_action.cpp" _measure_impl)
  foreach(_token
      "SoBRLMeasureAction::FULL_DETAIL"
      "SoBRLMeasureAction::DISPLAY_LEVEL"
      "shape->isLodDisplayActive()"
      "shape->hasFullDetailMesh()"
      "shape->needsSourceBackedFullDetail()"
      "shape->getFullDetailTriangle"
      "appendSourceBackedFullDetailRequest"
      "getSourceBackedFullDetailRequestCount"
      "makeSourceBackedFullDetailLodRequest"
      "consumeSourceBackedFullDetailResult"
      "submitSourceBackedFullDetailRequests"
      "consumeSourceBackedFullDetailResults"
      "measure_source_full_detail_result_valid"
      "measure_source_mesh_face_index"
      "brlobol_lod_rt_source_full_detail_request_from_source_mesh_request"
      "brlobol_lod_submit_rt_source_full_detail_request"
      "request.localToWorld = localToWorld"
      "this->haveQueryPoint"
      "request.queryBoundsValid = 1"
      "request.queryBounds.extendBy(localQuery)"
      "skippedLodDisplayMeshCount")
    string(FIND "${_measure_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/measure_action.cpp missing LoD geometry policy token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/snap_action.h" _snap_header)
  foreach(_token
      "enum GeometryPolicy"
      "FULL_DETAIL"
      "DISPLAY_LEVEL"
      "setGeometryPolicy"
      "getGeometryPolicy"
      "getSkippedLodDisplayMeshCount"
      "getSourceBackedFullDetailRequestCount"
      "makeSourceBackedFullDetailLodRequest"
      "consumeSourceBackedFullDetailResult"
      "submitSourceBackedFullDetailRequests"
      "consumeSourceBackedFullDetailResults")
    string(FIND "${_snap_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/snap_action.h missing LoD geometry policy token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/snap_action.cpp" _snap_impl)
  foreach(_token
      "SoBRLSnapAction::FULL_DETAIL"
      "SoBRLSnapAction::DISPLAY_LEVEL"
      "snap_source_local_tolerance"
      "shape->isLodDisplayActive()"
      "shape->hasFullDetailMesh()"
      "shape->needsSourceBackedFullDetail()"
      "shape->getFullDetailTriangle"
      "appendSourceBackedFullDetailRequest"
      "getSourceBackedFullDetailRequestCount"
      "makeSourceBackedFullDetailLodRequest"
      "consumeSourceBackedFullDetailResult"
      "submitSourceBackedFullDetailRequests"
      "consumeSourceBackedFullDetailResults"
      "snap_source_full_detail_result_valid"
      "snap_source_mesh_face_index"
      "brlobol_lod_rt_source_full_detail_request_from_source_mesh_request"
      "brlobol_lod_submit_rt_source_full_detail_request"
      "request.localToWorld = localToWorld"
      "request.queryBoundsValid = 1"
      "request.queryToleranceValid = 1"
      "skippedLodDisplayMeshCount")
    string(FIND "${_snap_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/snap_action.cpp missing LoD geometry policy token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_prototype.cpp" _prototype_test)
  foreach(_token
      "db_mesh_lod_update(dbip, \"tet.bot\")"
      "lodBotThreshold"
      "isLodBackedMesh"
      "SoBRLLodMeshShape"
      "lodAvailable"
      "lodStagedAvailable"
      "applyStagedLodResult"
      "BRLOBOL_LOD_RESULT_DIRECTORY"
      "BRLOBOL_LOD_RESULT_MESH"
      "BRLOBOL_LOD_PROXY_OBB"
      "lodDependencyPath"
      "lodAttributeName"
      "lodProxyKind"
      "lodProviderId"
      "lodFaceCount"
      "lodBoundsMin"
      "lodBoundsMax")
    string(FIND "${_prototype_test}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/tests/test_prototype.cpp missing LoD metadata coverage token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/CMakeLists.txt" _brlobol_tests_cmake)
  string(FIND "${_brlobol_tests_cmake}" "test_brlobol_lod_realization" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("src/libbrlobol/tests/CMakeLists.txt missing test_brlobol_lod_realization")
  endif()
  foreach(_token
      "test_brlobol_lod_service"
      "libbrlobol_lod_service"
      "libwdb"
      "test_brlobol_lod_update_action"
      "libbrlobol_lod_update_action")
    string(FIND "${_brlobol_tests_cmake}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/tests/CMakeLists.txt missing async LoD service test token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_lod_update_action.cpp" _lod_update_test)
  foreach(_token
      "SoBRLLodUpdateAction"
      "drainService"
      "BRLObolLodService"
      "getMatchedResultCount"
      "getAppliedResultCount"
      "getRejectedResultCount"
      "getUnmatchedResultCount"
      "BRLOBOL_LOD_PROVIDER_STALE"
      "BRLOBOL_LOD_PROVIDER_CACHE_MISS"
      "BRLOBOL_LOD_GEOMETRY_OBOL_MESH"
      "mesh_payload_result"
      "hasFullDetailMesh"
      "getFullDetailTriangleCount"
      "full-detail snap did not use preserved full-detail mesh"
      "getSkippedLodDisplayMeshCount"
      "SoBRLExportAction::DISPLAY_LEVEL"
      "SoBRLMeasureAction::DISPLAY_LEVEL"
      "SoBRLSnapAction::FULL_DETAIL"
      "SoGetBoundingBoxAction"
      "SoRayPickAction"
      "SoBRLMeshResidencyAction"
      "PICK_DISPLAY_LEVEL"
      "PICK_FULL_DETAIL"
      "setPickGeometryPolicy"
      "full-detail pick did not use preserved full-detail mesh"
      "camera->height = 90.0f"
      "setLodForcedLevel"
      "clearLodForcedLevel"
      "setForcedLevel"
      "makeLodRequest"
      "getRtViewInfo"
      "BRLOBOL_LOD_RESULT_AABB"
      "SoBRLMeshLodSubmitAction"
      "SoBRLLodMeshShape"
      "isLodBackedMesh"
      "isLodPreserveFullDetailEnabled"
      "isLodDisplayActive"
      "expected_view_lod_level"
      "rt_mesh_lod_load_view"
      "LoD submit action did not use view-policy active level"
      "LoD submit action did not skip already-resident view/policy request"
      "LoD view controller did not skip resident changed-scene LoD request"
      "LoD view controller did not skip resident threshold policy request"
      "LoD view controller queued duplicate resident request result"
      "LoD view controller queued duplicate threshold request result"
      "LoD-backed mesh retained full-detail payload after display LoD update"
      "LoD-backed mesh request did not keep source metrics without full-detail payload"
      "evictActiveDisplayMesh"
      "active display eviction did not preserve source-backed LoD identity"
      "evicted active display mesh did not keep source-backed exact export"
      "evicted active display mesh did not keep source-backed exact measure"
      "evicted active display mesh did not keep source-backed exact snap"
      "evicted active display mesh did not keep source-backed exact pick"
      "test_mesh_residency_budget_action"
      "mesh residency budget did not evict preserved full detail first"
      "controller mesh residency budget did not preserve source-backed exact identity"
      "test_view_controller_mesh_residency_budget_auto"
      "controller mesh residency budget did not auto-evict applied LoD payload"
      "setMeshResidencyBudget"
      "clearMeshResidencyBudget"
      "hasMeshResidencyBudget"
      "evictMeshPayloadsToBudget"
      "getLastMeshBudgetEvictedDisplayMeshCount"
      "exact export did not request source-backed full-detail LoD mesh"
      "exact export source-backed request did not convert to RT full-detail LoD request"
      "exact export did not consume source-backed full-detail LoD result"
      "controller source-backed exact export helper did not consume matching LoD result"
      "exact export did not submit source-backed full-detail LoD request"
      "exact export source-backed submit helper did not publish stale source result"
      "controller multi-source source-backed exact submit did not use matching database source"
      "controller multi-source source-backed exact submit did not consume matching database-scoped result"
      "test_view_controller_source_backed_partial_ready_submit"
      "controller partial-ready exact helper did not consume ready result and submit missing request"
      "controller partial-ready exact pick helper did not consume ready result and submit missing request"
      "exact measure did not request source-backed full-detail LoD mesh"
      "exact measure source-backed request did not carry bounded query metadata"
      "exact measure source-backed request did not convert to RT full-detail LoD request"
      "exact measure did not consume source-backed full-detail LoD result"
      "controller source-backed exact measure helper did not consume matching LoD result"
      "exact measure query distance limit did not filter resident nearest primitives"
      "bounded exact measure source-backed request did not carry explicit query tolerance"
      "bounded exact measure source-backed request did not convert explicit query tolerance"
      "exact snap did not request source-backed full-detail LoD mesh"
      "exact snap source-backed request did not carry bounded query metadata"
      "exact snap source-backed request did not convert to RT full-detail LoD request"
      "exact snap did not consume source-backed full-detail LoD result"
      "controller source-backed exact snap helper did not consume matching LoD result"
      "exact pick did not consume source-backed full-detail LoD result"
	      "SoBRLSourceMeshPickAction"
	      "exact pick action did not collect and consume source-backed full-detail LoD result"
	      "controller source-backed exact pick helper did not consume matching LoD result"
	      "mappedSourceResult.mesh.faceIndex"
      "mappedSourceResult.mesh.vertexIndex"
      "exact export did not preserve source face and vertex index mapping"
      "exact measure did not preserve source face and vertex index mapping"
      "exact measure did not accept compact bounded source subset with identity mapping"
      "exact measure accepted compact bounded source subset without face index mapping"
      "exact measure accepted compact bounded source subset without vertex index mapping"
      "exact snap did not preserve source face index mapping"
      "exact pick did not preserve source face and vertex index mapping"
      "exact snap did not preserve source vertex index mapping"
      "exact snap did not accept compact source vertex subset with vertex index mapping"
      "exact snap accepted compact source vertex subset without vertex index mapping"
      "exact pick did not accept ray-scoped source face and vertex subset"
      "exact pick accepted compact ray-scoped source point subset without vertex index mapping"
      "exact pick accepted non-query source face subset"
      "test_rt_exact_pick_provider"
      "RT exact pick provider did not return implicit comb hit identity"
      "BRLObolRtPickCache pickCache"
      "RT exact pick provider did not prepare reusable pick cache"
      "RT exact pick provider one-shot wrapper did not use cache-backed ray path"
      "SoBRLPickDetail::IMPLICIT_SOLID"
      "source_query.bounds"
      "source_query.tolerance"
      "source_query.ray.origin"
      "source_query.ray.direction"
      "test_mesh_lod_submit_action"
      "make_submit_test_db"
      "BU_DIR_CACHE"
      "getSubmittedTaskCount"
      "getSkippedMeshCount"
      "activeDuplicateSubmit"
      "LoD submit action did not skip active duplicate view/policy request"
      "lod-submit.bot"
      "test_view_controller_lod_submit_and_apply"
      "submitLodRequests"
      "getLastLodSubmittedTaskCount"
      "getLastLodAppliedResultCount"
      "getLodViewRevision"
      "setLodPolicyRevision"
      "setLodService"
      "setLodAutoSubmit"
      "isLodAutoSubmitEnabled"
      "hasPendingLodResults"
      "processPendingLodResults"
      "submitLodRequestsIfNeeded"
      "lod-result")
    string(FIND "${_lod_update_test}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/tests/test_lod_update_action.cpp missing staged LoD update action coverage token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_lod_realization.cpp" _lod_test)
  foreach(_token
      "brlobol_lod_cache_key"
      "brlobol_lod_mesh_payload_from_rt_mesh_data"
      "brlobol_lod_result_matches_request"
      "brlobol_lod_result_from_rt_mesh_info"
      "viewRevision++"
      "stale_cache_entry"
      "BRLOBOL_LOD_GEOMETRY_RT_MESH_CACHE"
      "brlobol_lod_directory_result"
      "brlobol_lod_attributes_result"
      "brlobol_lod_aabb_result"
      "brlobol_lod_proxy_result"
      "BRLOBOL_LOD_PROXY_OBB"
      "addDependency"
      "addAttribute")
    string(FIND "${_lod_test}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/tests/test_lod_realization.cpp missing staged LoD coverage token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_lod_service.cpp" _lod_service_test)
  foreach(_token
      "BRLObolLodService"
      "addDependency"
      "cancelGeneration"
      "BRLOBOL_LOD_PROVIDER_CANCELLED"
      "BRLOBOL_LOD_PROVIDER_STALE"
      "queuedCacheWriteCountForDiagnostics"
      "delayedTaskCountForDiagnostics"
      "submitIfNotActive"
      "hasActiveRequest"
      "test_active_request_duplicate_suppression"
      "LoD service accepted duplicate active request"
      "subscribeResultReady"
      "unsubscribeResultReady"
      "result_ready_cb"
      "test_result_ready_subscription"
      "realizeDataFree"
      "test_task_realize_data_cleanup"
      "cacheWriteOrder"
      "drainResults"
      "drainMatchingResults"
      "test_filtered_result_drain"
      "LoD filtered result drain did not isolate requested result"
      "staged_payload_task"
      "brlobol_rt_mesh_lod_provider_task"
      "BRLObolRtMeshLodProvider"
      "BU_DIR_CACHE"
      "useForcedLevel"
      "forcedLevel"
      "shrinkAfterCopy"
      "cachedNormals"
      "LoD RT provider did not store cached mesh normals"
      "cachedNormalProvider"
      "cachedNormalResult.counts.normalCount == 0"
      "cachedNormalResult.mesh.normals.size() !="
      "LoD RT provider did not return cached mesh normals"
      "db_mesh_lod_invalidate"
      "cleared_cache_entry"
      "refreshMissing = FALSE"
      "RT mesh LoD provider has no cache payload"
      "brlobol_rt_source_full_detail_provider_task"
      "brlobol_lod_rt_source_full_detail_request_from_source_mesh_request"
      "brlobol_lod_submit_rt_source_full_detail_request"
      "BRLObolRtSourceFullDetailProvider"
      "BRLOBOL_LOD_RESULT_FULL_DETAIL"
      "RT source full-detail provider source metrics changed"
      "RT source full-detail provider request exceeds full-detail limits"
      "LoD RT source full-detail provider query bounds did not reduce returned face payload"
      "LoD RT source full-detail provider scoped bounds limit should apply after payload reduction"
      "LoD RT source full-detail provider scoped bounds miss should not expand to whole-object payload"
      "LoD RT source full-detail provider should ignore non-source-local bounds when reducing payloads"
      "LoD RT source full-detail provider should not bypass whole-object limits for non-source-local bounds"
      "LoD RT source full-detail provider should ignore malformed source-local bounds when reducing payloads"
      "LoD RT source full-detail provider should not bypass whole-object limits for malformed source-local tolerance"
      "LoD RT source full-detail provider should ignore duplicate query-space params when reducing payloads"
      "LoD RT source full-detail provider should not bypass whole-object limits for duplicate bounds params"
      "LoD RT source full-detail provider should ignore mixed scoped query kinds when reducing payloads"
      "LoD RT source full-detail provider should not bypass whole-object limits for mixed scoped query kinds"
      "LoD RT source full-detail provider query ray did not reduce returned face payload"
      "scopedResult.mesh.faceIndex"
      "scopedResult.mesh.vertexIndex"
      "compactRayResult.mesh.vertexIndex"
      "LoD RT source full-detail provider query ray did not compact source vertex payload"
      "LoD RT source full-detail provider scoped ray limit should apply after payload reduction"
      "LoD RT source full-detail provider scoped ray miss should not expand to whole-object payload"
      "LoD RT source full-detail provider should ignore non-source-local rays when reducing payloads"
      "LoD RT source full-detail provider should ignore malformed source-local rays when reducing payloads"
      "LoD RT source full-detail provider should ignore duplicate ray params when reducing payloads"
      "LoD RT source full-detail provider should keep measure query hints whole-object without tolerance"
      "LoD RT source full-detail provider should reduce explicit bounded measure query payloads"
      "LoD RT source full-detail helper should replace stale template query params"
      "LoD RT source full-detail helper should submit current source query params after template cleanup"
      "LoD RT source full-detail helper did not submit source request"
      "test_rt_source_full_detail_provider_task"
      "mesh.isValid"
      "BRLOBOL_LOD_RESULT_ATTRIBUTES")
    string(FIND "${_lod_service_test}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/tests/test_lod_service.cpp missing async LoD service coverage token ${_token}")
    endif()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_brlobol_material_object)
  foreach(_rel
      include/brlobol/material_object.h
      include/brlobol/database_source.h
      include/brlobol/CMakeLists.txt
      src/libbrlobol/material_object.cpp
      src/libbrlobol/database_source.cpp
      src/libbrlobol/init.cpp
      src/libbrlobol/CMakeLists.txt
      src/libbrlobol/tests/test_prototype.cpp)
    if(NOT EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _brlobol_pivot_guard_fail("${_rel} is required for Obol material object coverage")
      continue()
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/material_object.h" _material_header)
  foreach(_token
      "class BRLOBOL_EXPORT SoBRLMaterialObject"
      "materialName"
      "parentName"
      "materialSource"
      "propertyGroup"
      "findProperty")
    string(FIND "${_material_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/material_object.h missing material metadata token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/database_source.h" _db_source_header)
  foreach(_token
      "SoBRLMaterialObject"
      "getRealizedMaterialObject"
      "getRealizedMaterialObjectCount")
    string(FIND "${_db_source_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/database_source.h missing material object accessor ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/database_source.cpp" _db_source_impl)
  foreach(_token
      "material_object_from_internal"
      "RT_CHECK_MATERIAL"
      "internalType == ID_MATERIAL"
      "assign_material_identity"
      "count_material_objects_in_node")
    string(FIND "${_db_source_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/database_source.cpp missing material object realization token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/init.cpp" _init_impl)
  string(FIND "${_init_impl}" "SoBRLMaterialObject::initClass()" _idx)
  if(_idx EQUAL -1)
    _brlobol_pivot_guard_fail("src/libbrlobol/init.cpp must register SoBRLMaterialObject")
  endif()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_prototype.cpp" _prototype_test)
  foreach(_token
      "exercise_generated_material_object"
      "database-backed material wire object coverage"
      "database-backed material shaded object coverage"
      "findProperty(\"physical\", \"density\""
      "material object contributed export geometry")
    string(FIND "${_prototype_test}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/tests/test_prototype.cpp missing material object coverage token ${_token}")
    endif()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_brlobol_curve_shaded_vlist)
  foreach(_rel
      src/libbrlobol/database_source.cpp
      src/libbrlobol/tests/test_prototype.cpp)
    if(NOT EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _brlobol_pivot_guard_fail("${_rel} is required for Obol curve shaded-vlist coverage")
      continue()
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/database_source.cpp" _db_source_impl)
  foreach(_token
      "vlist_from_plot_internal"
      "rt_obj_plot(&vhead, intern"
      "internalType == ID_SKETCH || internalType == ID_ANNOT")
    string(FIND "${_db_source_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/database_source.cpp missing curve shaded-vlist realization token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_prototype.cpp" _prototype_test)
  foreach(_token
      "exercise_generated_primitive_shaded_vlist"
      "shaded vlist pick did not hit line geometry"
      "shaded vlist snap did not preserve line/path identity"
      "database-backed sketch shaded vlist coverage"
      "database-backed annotation shaded vlist coverage")
    string(FIND "${_prototype_test}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/tests/test_prototype.cpp missing curve shaded-vlist coverage token ${_token}")
    endif()
  endforeach()
endfunction()

function(_brlobol_pivot_guard_check_brlobol_point_identity)
  foreach(_rel
      include/brlobol/export_action.h
      include/brlobol/vlist_shape.h
      src/libbrlobol/database_source.cpp
      src/libbrlobol/export_action.cpp
      src/libbrlobol/snap_action.cpp
      src/libbrlobol/vlist_shape.cpp
      src/libbrlobol/tests/test_prototype.cpp)
    if(NOT EXISTS "${BRLCAD_SOURCE_DIR}/${_rel}")
      _brlobol_pivot_guard_fail("${_rel} is required for Obol point identity coverage")
      continue()
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/vlist_shape.h" _vlist_header)
  foreach(_token
      "getPointPrimitiveCount"
      "getPointPrimitive"
      "pointColorValid"
      "pointScaleValid"
      "pointNormalValid"
      "getPointColor"
      "getPointScale"
      "getPointNormal")
    string(FIND "${_vlist_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/vlist_shape.h missing point primitive API ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/include/brlobol/export_action.h" _export_header)
  foreach(_token
      "struct PointRecord"
      "getPointCount"
      "appendPoint"
      "pointColorValid"
      "pointScaleValid"
      "pointNormalValid"
      "std::vector<PointRecord> points")
    string(FIND "${_export_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/export_action.h missing point export token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/database_source.cpp" _db_source_impl)
  foreach(_token
      "vlist_from_pnts"
      "RT_PNT_TYPE_PNT"
      "RT_PNT_TYPE_COL"
      "RT_PNT_TYPE_SCA"
      "RT_PNT_TYPE_NRM"
      "RT_PNT_TYPE_COL_SCA"
      "RT_PNT_TYPE_COL_NRM"
      "RT_PNT_TYPE_SCA_NRM"
      "RT_PNT_TYPE_COL_SCA_NRM"
      "defaultScale"
      "SoBRLVListShape::POINT"
      "setPointAttributes"
      "internalType == ID_PNTS")
    string(FIND "${_db_source_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/database_source.cpp missing PNTS point realization token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/export_action.cpp" _export_impl)
  foreach(_token
      "SoBRLExportAction::getPointCount"
      "shape->getPointPrimitiveCount()"
      "export_transform_point_normal"
      "export_transform_point_scale"
      "localToWorld.inverse().transpose()"
      "shape->getPointColor"
      "shape->getPointScale"
      "shape->getPointNormal"
      "SoBRLExportAction::appendPoint"
      "record.pointColorValid"
      "record.pointScaleValid"
      "record.pointNormalValid"
      "record.point = point")
    string(FIND "${_export_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/export_action.cpp missing point export implementation token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/snap_action.cpp" _snap_impl)
  foreach(_token
      "shape->command[i] != SoBRLVListShape::POINT"
      "snapAction->consider(ENDPOINT"
      "editIntentRole, i, query, worldPoint")
    string(FIND "${_snap_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/snap_action.cpp missing point snap token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_prototype.cpp" _prototype_test)
  foreach(_token
      "exercise_generated_pnts_shaded_points"
      "exercise_generated_pnts_shaded_attributes"
      "exercise_generated_pnts_attribute_variant"
      "shape->getPointPrimitiveCount() != 2"
      "RT_PNT_TYPE_COL_SCA_NRM"
      "pnts_col.s"
      "pnts_sca.s"
      "pnts_nrm.s"
      "pnts_col_sca.s"
      "pnts_col_nrm.s"
      "pnts_sca_nrm.s"
      "pnts_global_scale.s"
      "pnts_sca_precedence.s"
      "database-backed PNTS shaded point attribute realization"
      "database-backed PNTS color-scale attribute variant"
      "database-backed PNTS scale-normal attribute variant"
      "database-backed PNTS global scale fallback"
      "database-backed PNTS per-point scale precedence"
      "setScale(SbVec3f(2.0f, 3.0f, 4.0f))"
      "transformed shaded point-attribute export records are wrong"
      "pointExport.getPointCount() != 2"
      "pointSnap.getPrimitiveIndex() != 0"
      "database-backed PNTS shaded point realization")
    string(FIND "${_prototype_test}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/tests/test_prototype.cpp missing point realization coverage token ${_token}")
    endif()
  endforeach()
endfunction()

_brlobol_pivot_guard_check_dependency_inventory()
_brlobol_pivot_guard_check_obol_realization_coverage()
_brlobol_pivot_guard_check_libimgstream_boundary()
_brlobol_pivot_guard_check_brlobol_image_source()
_brlobol_pivot_guard_check_brlobol_image_display()
_brlobol_pivot_guard_check_brlobol_window_host()
_brlobol_pivot_guard_check_brlobol_mesh_identity()
_brlobol_pivot_guard_check_brlobol_edit_intent_identity()
_brlobol_pivot_guard_check_brlobol_material_object()
_brlobol_pivot_guard_check_brlobol_curve_shaded_vlist()
_brlobol_pivot_guard_check_brlobol_point_identity()
_brlobol_pivot_guard_check_brlobol_lod_metadata()
_brlobol_pivot_guard_check_plot3_ownership()
_brlobol_pivot_guard_check_libbg_bsg_neutralization()
_brlobol_pivot_guard_check_libbrep_bsg_neutralization()
_brlobol_pivot_guard_check_utility_vlist_neutralization()
_brlobol_pivot_guard_check_libanalyze_vlist_neutralization()
_brlobol_pivot_guard_check_libnmg_vlist_neutralization()
_brlobol_pivot_guard_check_librt_vlist_neutralization()
_brlobol_pivot_guard_check_librt_primitive_line_command_neutralization()
_brlobol_pivot_guard_check_librt_bsg_umbrella_neutralization()
_brlobol_pivot_guard_check_librt_edit_knob_neutralization()
_brlobol_pivot_guard_check_librt_view_info_neutralization()
_brlobol_pivot_guard_check_librt_sketch_polygon_neutralization()
_brlobol_pivot_guard_check_qtcad_obol_test_links()
_brlobol_pivot_guard_check_qtcad_obol_edit_preview_intent()
_brlobol_pivot_guard_check_qtcad_obol_export_source_exact()
_brlobol_pivot_guard_check_qtcad_obol_pick_source_exact()
_brlobol_pivot_guard_check_qtcad_obol_snap_source_exact()
_brlobol_pivot_guard_check_qtcad_obol_measure_source_exact()
_brlobol_pivot_guard_check_qtcad_window_host_adapter()
_brlobol_pivot_guard_check_qtcad_measure_filter()
_brlobol_pivot_guard_check_qtcad_selection_api()
_brlobol_pivot_guard_check_qtcad_view_snapshot_adapters()
_brlobol_pivot_guard_check_qged_view_info_rt_adapter()

if("${BRLOBOL_PIVOT_GUARD_MODE}" STREQUAL "strict")
  file(GLOB_RECURSE _files LIST_DIRECTORIES false
    "${BRLCAD_SOURCE_DIR}/include/*"
    "${BRLCAD_SOURCE_DIR}/src/*"
    "${BRLCAD_SOURCE_DIR}/misc/CMake/*"
  )
else()
  _brlobol_pivot_guard_collect(_files
    include/brlobol
    include/brlobol.h
    src/libbrlobol
    # QgObolDrawSync remains the BSG/GED transaction adapter; the helpers
    # listed below are the qtcad Obol-canonical API surface.
    include/qtcad/QgObolDatabaseSync.h
    include/qtcad/QgObolEditPreview.h
    include/qtcad/QgObolExport.h
    include/qtcad/QgObolMeasure.h
    include/qtcad/QgObolOverlaySync.h
    include/qtcad/QgObolPick.h
    include/qtcad/QgObolSelectionSync.h
    include/qtcad/QgObolSnap.h
    include/qtcad/QgObolWindowHost.h
    src/libqtcad/QgObolDatabaseSync.cpp
    src/libqtcad/QgObolEditPreview.cpp
    src/libqtcad/QgObolExport.cpp
    src/libqtcad/QgObolMeasure.cpp
    src/libqtcad/QgObolOverlaySync.cpp
    src/libqtcad/QgObolPick.cpp
    src/libqtcad/QgObolSelectionSync.cpp
    src/libqtcad/QgObolSnap.cpp
    src/libqtcad/QgObolWindowHost.cpp
  )
endif()

foreach(_file IN LISTS _files)
  _brlobol_pivot_guard_check_file("${_file}")
endforeach()

get_property(_failures GLOBAL PROPERTY BRLOBOL_PIVOT_GUARD_FAILURES)
if(_failures)
  list(SORT _failures)
  string(REPLACE ";" "\n  " _failure_text "${_failures}")
  message(FATAL_ERROR
    "Obol pivot guard failed in ${BRLOBOL_PIVOT_GUARD_MODE} mode:\n"
    "  ${_failure_text}\n"
    "New Obol-canonical code must not depend on BSG/DM APIs.")
endif()

message(STATUS "Obol pivot guard passed in ${BRLOBOL_PIVOT_GUARD_MODE} mode")
