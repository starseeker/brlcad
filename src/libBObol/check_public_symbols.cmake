# Verify the installed libBObol boundary does not re-export bridge helpers or
# adapters removed by the direct libBObol pivot.

if(NOT DEFINED LIBRARY OR NOT EXISTS "${LIBRARY}")
  message(FATAL_ERROR "libBObol symbol contract: library is unavailable: ${LIBRARY}")
endif()

if(DEFINED NM AND EXISTS "${NM}")
  execute_process(
    COMMAND "${NM}" -D -C --defined-only "${LIBRARY}"
    RESULT_VARIABLE _symbol_result
    OUTPUT_VARIABLE _symbols
    ERROR_VARIABLE _symbol_error)
  set(_symbol_tool "nm")
elseif(DEFINED DUMPBIN AND EXISTS "${DUMPBIN}")
  execute_process(
    COMMAND "${DUMPBIN}" /nologo /exports "${LIBRARY}"
    RESULT_VARIABLE _symbol_result
    OUTPUT_VARIABLE _symbols
    ERROR_VARIABLE _symbol_error)
  set(_symbol_tool "dumpbin")
else()
  message(FATAL_ERROR
    "libBObol symbol contract: neither nm nor dumpbin is available")
endif()
if(NOT _symbol_result EQUAL 0)
  message(FATAL_ERROR
    "libBObol symbol contract: ${_symbol_tool} failed: ${_symbol_error}")
endif()

foreach(_required
    "bobol_display_endpoint_create"
    "bobol_display_endpoint_destroy"
    "bobol_display_endpoint_property_count")
  string(FIND "${_symbols}" "${_required}" _required_position)
  if(_required_position EQUAL -1)
    message(FATAL_ERROR
      "libBObol symbol contract: required stable symbol ${_required} is absent")
  endif()
endforeach()

foreach(_forbidden
    "ged_bobol_"
    "bobol_evaluated_wire_"
    "BObolViewController::Impl::"
    "BObolSceneController::Impl::"
    "SoBRLDatabaseSource::Impl::")
  string(FIND "${_symbols}" "${_forbidden}" _forbidden_position)
  if(NOT _forbidden_position EQUAL -1)
    message(FATAL_ERROR
      "libBObol symbol contract: private or removed symbol matches ${_forbidden}")
  endif()
endforeach()
