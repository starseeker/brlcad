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
  include/ged/draw_types.h
  include/ged/scene.h
  include/ged/scene_types.h
  include/ged/view_feature.h
  include/ged/view_feature_types.h
  include/ged/polygon.h
  include/ged/polygon_types.h
  include/ged/selection.h
  include/ged/selection_types.h
  include/ged/view_feature_batch.h
  include/ged/display.h
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
    "(const[ \t\r\n]+)?void[ \t\r\n]*\\*[ \t\r\n]*(view|view_ctx)[ \t\r\n]*[,)]"
    _untyped_view_parameter "${_contents}")
  if(_untyped_view_parameter)
    message(FATAL_ERROR
      "GED draw API contract: ${_header} exposes an untyped view parameter: ${_untyped_view_parameter}")
  endif()
  string(REGEX MATCHALL "GED_EXPORT[^;]*;" _declarations "${_contents}")
  foreach(_declaration IN LISTS _declarations)
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
    "^ged_draw_obol_(controller|scene_controller)")
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
