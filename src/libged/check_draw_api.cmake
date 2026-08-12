# Verify the deliberately stable GED draw/view API inventory.  This check is
# intentionally declaration based: implementation-only helpers and dynamic
# plugin command wrappers do not belong in this manifest.

cmake_policy(SET CMP0007 NEW)

if(NOT DEFINED SOURCE_ROOT OR NOT DEFINED MANIFEST)
  message(FATAL_ERROR "check_draw_api.cmake requires SOURCE_ROOT and MANIFEST")
endif()

set(_headers
  include/ged/view_types.h
  include/ged/view.h
  include/ged/scene.h
  include/ged/scene_types.h
  include/ged/view_feature.h
  include/ged/view_feature_types.h
  include/ged/view_edit.h
  include/ged/view_edit_types.h
  include/ged/view_export.h
  include/ged/view_export_types.h
  include/ged/view_polygon.h
  include/ged/view_polygon_types.h
  include/ged/selection.h
  include/ged/selection_types.h
  include/ged/view_feature_batch.h
  include/ged/display.h
  include/ged/event.h
)

set(_actual)
foreach(_header IN LISTS _headers)
  file(READ "${SOURCE_ROOT}/${_header}" _contents)
  string(REGEX MATCH
    "#[ \t]*include[ \t]*[<\"]common\\.h[>\"]"
    _common_header_include "${_contents}")
  if(_common_header_include)
    message(FATAL_ERROR
      "GED draw API contract: installed header ${_header} includes common.h")
  endif()
  string(REGEX MATCH
    "#[ \t]*include[ \t]*[<\"]BObol/"
    _concrete_renderer_include "${_contents}")
  if(_concrete_renderer_include)
    message(FATAL_ERROR
      "GED draw API contract: ${_header} includes a concrete BObol renderer header")
  endif()
  string(REGEX MATCH
    "(struct[ \t]+bobol_|bobol_display_endpoint|bobol_endpoint_property)"
    _concrete_renderer_type "${_contents}")
  if(_concrete_renderer_type)
    message(FATAL_ERROR
      "GED draw API contract: ${_header} exposes a concrete BObol renderer type: ${_concrete_renderer_type}")
  endif()
  string(REGEX MATCH
    "#[ \t]*include[ \t]*[<\"][^>\"]*(private|_private)\\.h[>\"]"
    _private_header_include "${_contents}")
  if(_private_header_include)
    message(FATAL_ERROR
      "GED draw API contract: ${_header} includes a private implementation header")
  endif()
  string(REGEX MATCH "(Obol|Coin3D|OpenGL|OSMesa)"
    _renderer_vocabulary "${_contents}")
  if(_renderer_vocabulary)
    message(FATAL_ERROR
      "GED draw API contract: ${_header} references renderer vocabulary: ${_renderer_vocabulary}")
  endif()
  foreach(_obsolete_init
      GED_VIEW_FEATURE_STYLE_INIT
      GED_VIEW_FEATURE_SUMMARY_INIT
      GED_VIEW_FEATURE_BATCH_DESC_INIT
      GED_VIEW_FEATURE_BATCH_EVENT_INIT
      GED_VIEW_FEATURE_LINE_LAYER_INIT
      GED_VIEW_FEATURE_METADATA_INIT
      GED_VIEW_EDIT_TRANSACTION_INIT
      GED_PICK_DETAIL_INIT)
    string(REGEX MATCH "#[ \t]*define[ \t]+${_obsolete_init}([ \t]|$)"
      _obsolete_init_match "${_contents}")
    if(_obsolete_init_match)
      message(FATAL_ERROR
        "GED draw API contract: obsolete initializer alias ${_obsolete_init}; use its typed init/default function")
    endif()
  endforeach()
  string(REGEX MATCH
    "(const[ \t\r\n]+)?void[ \t\r\n]*\\*[ \t\r\n]*(view|view_ctx)[ \t\r\n]*[,)]"
    _untyped_view_parameter "${_contents}")
  if(_untyped_view_parameter)
    message(FATAL_ERROR
      "GED draw API contract: ${_header} exposes an untyped view parameter: ${_untyped_view_parameter}")
  endif()
  string(REGEX MATCHALL
    "(^|\n)[ \t]*GED_EXPORT[ \t\r\n]+(extern[ \t\r\n]+)?[^;]*;"
    _declarations "${_contents}")
  foreach(_declaration IN LISTS _declarations)
    string(REGEX MATCH "ged_[A-Za-z0-9_]+" _declaration_name
      "${_declaration}")
    if(NOT _declaration_name)
      continue()
    endif()
    string(FIND "${_contents}" "${_declaration}" _declaration_offset)
    string(SUBSTRING "${_contents}" 0 ${_declaration_offset}
      _declaration_prefix)
    string(REGEX MATCH
      "/\\*\\*([^*]|\\*+[^*/])*\\*+/[ \t\r\n]*$"
      _declaration_documentation "${_declaration_prefix}")
    if(NOT _declaration_documentation)
      message(FATAL_ERROR
        "GED draw API contract: ${_header} declaration ${_declaration_name} has no immediately preceding Doxygen documentation")
    endif()
    string(REGEX MATCH "ged_[A-Za-z0-9_]+[ \t\r\n]*\\(" _name "${_declaration}")
    if(_name)
      string(REGEX REPLACE "[ \t\r\n]*\\($" "" _name "${_name}")
      list(APPEND _actual "${_name}")
    endif()
  endforeach()
endforeach()
list(REMOVE_DUPLICATES _actual)
list(SORT _actual)

foreach(_obsolete_prefix
    "^ged_draw_view_context_"
    "^ged_draw_feature_batch_"
    "^ged_draw_obol_database_source_"
    "^ged_draw_obol_(controller|scene_controller)"
    "^ged_view_context_(tclcad|user_data)_")
  foreach(_name IN LISTS _actual)
    if(_name MATCHES "${_obsolete_prefix}")
      message(FATAL_ERROR
        "GED draw API contract: obsolete public vocabulary ${_name}; use the typed scene/source/view-feature/result or endpoint owner")
    endif()
  endforeach()
endforeach()

if(PRINT_ONLY)
  foreach(_name IN LISTS _actual)
    message(STATUS "GED_DRAW_API ${_name}")
  endforeach()
  return()
endif()

file(STRINGS "${MANIFEST}" _expected)
list(FILTER _expected EXCLUDE REGEX "^[ \t]*(#|$)")
list(REMOVE_DUPLICATES _expected)
list(SORT _expected)

if(NOT "${_actual}" STREQUAL "${_expected}")
  set(_added ${_actual})
  set(_removed ${_expected})
  foreach(_name IN LISTS _expected)
    list(REMOVE_ITEM _added "${_name}")
  endforeach()
  foreach(_name IN LISTS _actual)
    list(REMOVE_ITEM _removed "${_name}")
  endforeach()
  string(REPLACE ";" "\n  " _added_text "${_added}")
  string(REPLACE ";" "\n  " _removed_text "${_removed}")
  message(FATAL_ERROR
    "GED draw API manifest mismatch.\n"
    "New declarations (review and add deliberately):\n  ${_added_text}\n"
    "Missing declarations (review and remove deliberately):\n  ${_removed_text}")
endif()
