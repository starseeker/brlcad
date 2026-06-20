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
	"BRLObolViewController::setLodService"
	"BRLObolViewController::setLodAutoSubmit"
	"BRLObolViewController::submitLodRequestsIfNeeded"
	"BRLObolViewController::processPendingLodResults"
	"BRLObolViewController::submitLodRequests"
	"BRLObolViewController::applyLodResults"
	"lodResultReadyCB"
	"lodResultsPending"
	"lodAutoSubmit"
	"lodViewSignature"
	"lodUseForcedLevel"
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
	"bsg/view_state.h"
	"rt/edit_legacy_bsg.h"
	"rt_edit_view_from_bsg"
	"rt_edit_create_bsg"
	"rt_edit_reinit_bsg"
	"rt_edit_knob_cmd_process_bsg")
      string(FIND "${_edit_legacy_bsg_impl_contents}" "${_token}" _legacy_impl_idx)
      if(_legacy_impl_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/librt/edit_legacy_bsg.cpp must own transitional edit adapter ${_token}")
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
	  "rt_mesh_lod_load_level"
	  "rt_mesh_lod_load_view"
	  "rt_mesh_lod_current_level"
	  "rt_mesh_lod_has_active_data"
	  "rt_mesh_lod_data_get"
	  "rt_mesh_lod_info_get"
	  "rt_mesh_lod_detail_init"
	  "rt_mesh_lod_detail_callbacks_set"
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
	rt_view_orientation_quat_from_bsg
	rt_view_aet_from_bsg
	rt_view_perspective_from_bsg
	rt_view_model2view_from_bsg
	rt_view_view2model_from_bsg
	rt_view_rotation_from_bsg
	rt_view_center_from_bsg
	rt_view_lod_policy_from_bsg
	rt_view_lod_policy_apply_bsg
	rt_view_lod_policy_copy_bsg
	rt_view_lod_bounds_update_bsg
	rt_view_lod_bounds_callback_set_bsg
	rt_view_lod_bounds_callback_is_bsg
	rt_view_scale_from_bsg
	rt_view_eye_pos_from_bsg
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
	rt_mesh_lod_detail_init
	rt_mesh_lod_detail_callbacks_set
	rt_mesh_lod_detail_setup_bsg)
      string(FIND "${_librt_cache_lod_impl_contents}" "${_token}"
	_librt_cache_lod_impl_token_idx)
      if(_librt_cache_lod_impl_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/librt/cache_lod.c must keep neutral RT mesh-LoD helper ${_token}")
      endif()
    endforeach()
    string(REGEX MATCH [[(^|[^_A-Za-z0-9])rt_mesh_lod_bsg[ \t\r\n]*\(]]
      _librt_cache_lod_public_bsg_hit "${_librt_cache_lod_impl_contents}")
    if(_librt_cache_lod_public_bsg_hit)
      _brlobol_pivot_guard_fail(
	"src/librt/cache_lod.c must keep public BSG extraction behind view_legacy_bsg.c")
    endif()
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
	[[#[ \t]*include[ \t]*[<"]bsg/view_state\.h]]
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	rt_view_info_from_bsg
	rt_view_orientation_quat_from_bsg
	rt_view_aet_from_bsg
	rt_view_perspective_from_bsg
	rt_view_model2view_from_bsg
	rt_view_view2model_from_bsg
	rt_view_rotation_from_bsg
	rt_view_center_from_bsg
	rt_view_lod_policy_from_bsg
	rt_view_lod_policy_apply_bsg
	rt_view_lod_policy_copy_bsg
	rt_view_lod_bounds_update_bsg
	rt_view_lod_bounds_callback_set_bsg
	rt_view_lod_bounds_callback_is_bsg
	rt_view_scale_from_bsg
	rt_view_eye_pos_from_bsg
	rt_mesh_lod_load_view_scene_ref_bsg
	rt_mesh_lod_free_scene_ref_bsg
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
	rt_mesh_lod_memshrink
	rt_mesh_lod_cache_status_init
	db_mesh_lod_status
	db_mesh_lod_refresh
	db_mesh_lod_invalidate
	db_mesh_lod_store_mesh
	db_mesh_lod_update
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
	"rt/view_legacy_bsg.h"
	rt_view_info_from_bsg
	rt_view_orientation_quat_from_bsg
	rt_view_aet_from_bsg
	rt_view_perspective_from_bsg
	rt_view_model2view_from_bsg
	rt_view_view2model_from_bsg
	rt_view_rotation_from_bsg
	rt_view_center_from_bsg
	rt_view_lod_policy_from_bsg
	rt_view_lod_policy_apply_bsg
	rt_view_lod_policy_copy_bsg
	rt_view_lod_bounds_update_bsg
	rt_view_lod_bounds_callback_set_bsg
	rt_view_lod_bounds_callback_is_bsg
	rt_view_scale_from_bsg
	rt_view_eye_pos_from_bsg
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
	test_bsg_mesh_lod_adapter_boundary)
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
  endforeach()

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
  endforeach()

  set(_mged_attach_impl "${BRLCAD_SOURCE_DIR}/src/mged/attach.c")
  if(EXISTS "${_mged_attach_impl}")
    file(READ "${_mged_attach_impl}" _mged_attach_contents)
    foreach(_token
	"rt/view.h"
	"RT_VIEW_MIN"
	"RT_VIEW_MAX")
      string(FIND "${_mged_attach_contents}" "${_token}" _mged_attach_token_idx)
      if(_mged_attach_token_idx EQUAL -1)
	_brlobol_pivot_guard_fail(
	  "src/mged/attach.c must use RT-owned view window bounds ${_token}")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])BSG_VIEW_(MIN|MAX)([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _mged_attach_bsg_view_hit
	"${_mged_attach_contents}")
      if(_mged_attach_bsg_view_hit)
	_brlobol_pivot_guard_fail(
	  "src/mged/attach.c must not use BSG-owned spellings for view window bounds: ${_mged_attach_bsg_view_hit}")
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

  set(_libtclcad_view_draw_impl "${BRLCAD_SOURCE_DIR}/src/libtclcad/view/draw.c")
  if(EXISTS "${_libtclcad_view_draw_impl}")
    file(READ "${_libtclcad_view_draw_impl}" _libtclcad_view_draw_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_model2view_from_bsg]]
	[[rt_view_perspective_from_bsg]])
      string(REGEX MATCH "${_token}" _libtclcad_view_draw_adapter_hit
	"${_libtclcad_view_draw_contents}")
      if(NOT _libtclcad_view_draw_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/view/draw.c must route draw matrix/perspective reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_model2view([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_perspective([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libtclcad_view_draw_direct_hit
	"${_libtclcad_view_draw_contents}")
      if(_libtclcad_view_draw_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/view/draw.c reintroduced direct BSG draw matrix/perspective reads: ${_libtclcad_view_draw_direct_hit}")
      endif()
    endforeach()
  endif()

  set(_libtclcad_wrapper_impl "${BRLCAD_SOURCE_DIR}/src/libtclcad/wrapper.c")
  if(EXISTS "${_libtclcad_wrapper_impl}")
    file(READ "${_libtclcad_wrapper_impl}" _libtclcad_wrapper_contents)
    foreach(_token
	[[#[ \t]*include[ \t]*[<"]rt/view_legacy_bsg\.h]]
	[[rt_view_perspective_from_bsg]]
	[[rt_view_lod_policy_from_bsg]])
      string(REGEX MATCH "${_token}" _libtclcad_wrapper_adapter_hit
	"${_libtclcad_wrapper_contents}")
      if(NOT _libtclcad_wrapper_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/wrapper.c must route perspective and zoom-refresh LoD policy reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[(^|[^A-Za-z0-9_])gv_perspective([^A-Za-z0-9_]|$)]]
	[[struct[ \t\r\n]+bsg_lod_source_policy_settings]]
	[[bsg_view_lod_source_policy_get]])
      string(REGEX MATCH "${_pat}" _libtclcad_wrapper_direct_hit
	"${_libtclcad_wrapper_contents}")
      if(_libtclcad_wrapper_direct_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/wrapper.c reintroduced direct BSG perspective/LoD policy reads: ${_libtclcad_wrapper_direct_hit}")
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
	[[rt_view_center_from_bsg]]
	[[rt_view_scale_from_bsg]])
      string(REGEX MATCH "${_token}" _libtclcad_commands_adapter_hit
	"${_libtclcad_commands_contents}")
      if(NOT _libtclcad_commands_adapter_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/commands.c must route TclCAD LoD policy and command view reads through rt/view_legacy_bsg.h")
      endif()
    endforeach()
    foreach(_pat
	[[struct[ \t\r\n]+bsg_lod_source_policy_settings]]
	[[bsg_view_lod_source_policy_get]]
	[[bsg_view_lod_source_policy_set]]
	[[(^|[^A-Za-z0-9_])gv_center([^A-Za-z0-9_]|$)]]
	[[(^|[^A-Za-z0-9_])gv_scale([^A-Za-z0-9_]|$)]])
      string(REGEX MATCH "${_pat}" _libtclcad_commands_bsg_policy_hit
	"${_libtclcad_commands_contents}")
      if(_libtclcad_commands_bsg_policy_hit)
	_brlobol_pivot_guard_fail(
	  "src/libtclcad/commands.c reintroduced direct BSG LoD policy/view reads: ${_libtclcad_commands_bsg_policy_hit}")
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
	[[->[ \t\r\n]*gv_scale]])
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

function(_brlobol_pivot_guard_check_brlobol_mesh_identity)
  foreach(_rel
      include/brlobol/pick_detail.h
      include/brlobol/export_action.h
      include/brlobol/mesh_shape.h
      src/libbrlobol/pick_detail.cpp
      src/libbrlobol/mesh_shape.cpp
      src/libbrlobol/export_action.cpp
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
      "faceVertexIndex")
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
      "SoBRLPickDetail::setNearestFaceVertex")
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
      "record.vertexIndexC")
    string(FIND "${_export_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/export_action.cpp missing mesh vertex identity token ${_token}")
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
      "tri.vertexIndexC != 2")
    string(FIND "${_prototype_test}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/tests/test_prototype.cpp missing mesh vertex identity coverage token ${_token}")
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
      "SoBRLLodMeshShape::initClass"
      "SoBRLLodUpdateAction::initClass"
      "SoBRLMeshLodSubmitAction::initClass")
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
      "BRLObolLodTask"
      "BRLObolLodService"
      "start"
      "stop"
      "beginGeneration"
      "cancelGeneration"
      "isGenerationCancelled"
      "submit"
      "drainResults"
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
      "BRLOBOL_LOD_PROVIDER_CANCELLED"
      "BRLOBOL_LOD_PROVIDER_STALE"
      "brlobol_rt_mesh_lod_provider_task"
      "brlobol_rt_mesh_lod_provider_free"
      "db_mesh_lod_status"
      "db_mesh_lod_refresh"
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
      "BRLObolRtMeshLodProvider"
      "brlobol_rt_mesh_lod_provider_task"
      "brlobol_rt_mesh_lod_provider_free"
      "realizeDataFree"
      "service->submit"
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
      "SoBRLMeshShape::updateSourceMeshMetricsFromFields"
      "SoBRLMeshShape::captureFullDetailMesh"
      "SoBRLMeshShape::restoreFullDetailMesh"
      "SoBRLMeshShape::getFullDetailTriangle"
      "getFullDetailTriangleVertexIndices"
      "fullDetailPoint"
      "fullDetailCoordIndex"
      "sourceMeshMetricsValid"
      "sourceMeshBounds"
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
      "shape->getFullDetailTriangle"
      "skippedLodDisplayMeshCount"
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
      "getSkippedLodDisplayMeshCount")
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
      "shape->getFullDetailTriangle"
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
      "getSkippedLodDisplayMeshCount")
    string(FIND "${_snap_header}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("include/brlobol/snap_action.h missing LoD geometry policy token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/snap_action.cpp" _snap_impl)
  foreach(_token
      "SoBRLSnapAction::FULL_DETAIL"
      "SoBRLSnapAction::DISPLAY_LEVEL"
      "shape->isLodDisplayActive()"
      "shape->hasFullDetailMesh()"
      "shape->getFullDetailTriangle"
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
      "SoRayPickAction"
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
      "LoD-backed mesh retained full-detail payload after display LoD update"
      "LoD-backed mesh request did not keep source metrics without full-detail payload"
      "test_mesh_lod_submit_action"
      "make_submit_test_db"
      "BU_DIR_CACHE"
      "getSubmittedTaskCount"
      "getSkippedMeshCount"
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
      "subscribeResultReady"
      "unsubscribeResultReady"
      "result_ready_cb"
      "test_result_ready_subscription"
      "realizeDataFree"
      "test_task_realize_data_cleanup"
      "cacheWriteOrder"
      "drainResults"
      "staged_payload_task"
      "brlobol_rt_mesh_lod_provider_task"
      "BRLObolRtMeshLodProvider"
      "BU_DIR_CACHE"
      "useForcedLevel"
      "forcedLevel"
      "shrinkAfterCopy"
      "db_mesh_lod_invalidate"
      "cleared_cache_entry"
      "refreshMissing = FALSE"
      "RT mesh LoD provider has no cache payload"
      "mesh.isValid"
      "BRLOBOL_LOD_RESULT_ATTRIBUTES")
    string(FIND "${_lod_service_test}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/tests/test_lod_service.cpp missing async LoD service coverage token ${_token}")
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
      "getPointPrimitive")
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
      "SoBRLVListShape::POINT"
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
      "SoBRLExportAction::appendPoint"
      "record.point = point")
    string(FIND "${_export_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/export_action.cpp missing point export implementation token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/snap_action.cpp" _snap_impl)
  foreach(_token
      "shape->getPointPrimitiveCount()"
      "snapAction->consider(ENDPOINT"
      "primitiveIndex, query, worldPoint")
    string(FIND "${_snap_impl}" "${_token}" _idx)
    if(_idx EQUAL -1)
      _brlobol_pivot_guard_fail("src/libbrlobol/snap_action.cpp missing point snap token ${_token}")
    endif()
  endforeach()

  file(READ "${BRLCAD_SOURCE_DIR}/src/libbrlobol/tests/test_prototype.cpp" _prototype_test)
  foreach(_token
      "exercise_generated_pnts_shaded_points"
      "shape->getPointPrimitiveCount() != 2"
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
_brlobol_pivot_guard_check_qtcad_window_host_adapter()
_brlobol_pivot_guard_check_qtcad_measure_filter()
_brlobol_pivot_guard_check_qtcad_selection_api()

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
    include/qtcad/QgObolMeasure.h
    include/qtcad/QgObolOverlaySync.h
    include/qtcad/QgObolPick.h
    include/qtcad/QgObolSelectionSync.h
    include/qtcad/QgObolSnap.h
    include/qtcad/QgObolWindowHost.h
    src/libqtcad/QgObolDatabaseSync.cpp
    src/libqtcad/QgObolEditPreview.cpp
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
