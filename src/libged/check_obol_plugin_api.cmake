# Verify the deliberately opt-in Obol extension separately from GED's
# renderer-neutral draw/view API.

if(NOT DEFINED SOURCE_ROOT OR NOT DEFINED MANIFEST)
  message(FATAL_ERROR
    "check_obol_plugin_api.cmake requires SOURCE_ROOT and MANIFEST")
endif()

set(_header "${SOURCE_ROOT}/include/ged/plugin/obol.h")
file(READ "${_header}" _contents)
string(REGEX MATCHALL "GED_EXPORT[^;]*;" _declarations "${_contents}")
set(_actual)
foreach(_declaration IN LISTS _declarations)
  string(REGEX MATCH "ged_plugin_obol_[A-Za-z0-9_]+[ \t\r\n]*\\(" _name
    "${_declaration}")
  if(_name)
    string(REGEX REPLACE "[ \t\r\n]*\\($" "" _name "${_name}")
    list(APPEND _actual "${_name}")
  endif()
endforeach()
list(REMOVE_DUPLICATES _actual)
list(SORT _actual)

file(STRINGS "${MANIFEST}" _expected)
list(FILTER _expected EXCLUDE REGEX "^[ \t]*(#|$)")
list(REMOVE_DUPLICATES _expected)
list(SORT _expected)

if(NOT "${_actual}" STREQUAL "${_expected}")
  message(FATAL_ERROR
    "GED Obol plugin API manifest mismatch: actual=${_actual}; expected=${_expected}")
endif()
