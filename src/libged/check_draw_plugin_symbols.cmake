# Audit the deliberately exported, non-installed drawing ABI used by in-tree
# GED command modules and backend adapters.

cmake_policy(SET CMP0057 NEW)

if(NOT DEFINED NM OR NOT EXISTS "${NM}")
  message(FATAL_ERROR "libged drawing plugin symbol contract: nm is unavailable")
endif()
if(NOT DEFINED LIBRARY OR NOT EXISTS "${LIBRARY}")
  message(FATAL_ERROR
    "libged drawing plugin symbol contract: library is unavailable: ${LIBRARY}")
endif()
if(NOT DEFINED MANIFEST OR NOT EXISTS "${MANIFEST}")
  message(FATAL_ERROR
    "libged drawing plugin symbol contract: manifest is unavailable: ${MANIFEST}")
endif()

execute_process(
  COMMAND "${NM}" -D --defined-only "${LIBRARY}"
  RESULT_VARIABLE _nm_result
  OUTPUT_VARIABLE _nm_output
  ERROR_VARIABLE _nm_error)
if(NOT _nm_result EQUAL 0)
  message(FATAL_ERROR
    "libged drawing plugin symbol contract: nm failed: ${_nm_error}")
endif()

string(REPLACE "\n" ";" _nm_lines "${_nm_output}")
set(_exported)
foreach(_line IN LISTS _nm_lines)
  if(_line MATCHES "[ \t]([A-Za-z_][A-Za-z0-9_]*)$")
    list(APPEND _exported "${CMAKE_MATCH_1}")
  endif()
endforeach()

file(STRINGS "${MANIFEST}" _manifest_lines)
set(_approved)
foreach(_line IN LISTS _manifest_lines)
  string(STRIP "${_line}" _symbol)
  if(_symbol STREQUAL "" OR _symbol MATCHES "^#")
    continue()
  endif()
  list(APPEND _approved "${_symbol}")
  if(NOT _symbol IN_LIST _exported)
    message(FATAL_ERROR
      "libged drawing plugin symbol contract: approved symbol ${_symbol} is not exported")
  endif()
endforeach()

# Families in this private ABI must be enumerated, so a new GED_EXPORT cannot
# silently expand the surface.
foreach(_symbol IN LISTS _exported)
  if(_symbol MATCHES "^_ged_(uplot|view_feature_batch_)" OR
      _symbol MATCHES "^ged_edit_buf_" OR
      _symbol STREQUAL "ged_scene_backend_set_private" OR
      _symbol STREQUAL "ged_view_data_lines" OR
      _symbol STREQUAL "ged_view_feature_gobject_create")
    if(NOT _symbol IN_LIST _approved)
      message(FATAL_ERROR
        "libged drawing plugin symbol contract: unreviewed export ${_symbol}")
    endif()
  endif()
endforeach()
