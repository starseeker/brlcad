# Verify the installed libBObol boundary does not re-export bridge helpers or
# adapters removed by the direct libBObol pivot.

if(NOT DEFINED NM OR NOT EXISTS "${NM}")
  message(FATAL_ERROR "libBObol symbol contract: nm executable is unavailable")
endif()
if(NOT DEFINED LIBRARY OR NOT EXISTS "${LIBRARY}")
  message(FATAL_ERROR "libBObol symbol contract: library is unavailable: ${LIBRARY}")
endif()

execute_process(
  COMMAND "${NM}" -D -C --defined-only "${LIBRARY}"
  RESULT_VARIABLE _nm_result
  OUTPUT_VARIABLE _symbols
  ERROR_VARIABLE _nm_error)
if(NOT _nm_result EQUAL 0)
  message(FATAL_ERROR "libBObol symbol contract: nm failed: ${_nm_error}")
endif()

foreach(_required
    "bobol_display_endpoint_create"
    "bobol_display_endpoint_destroy"
    "bobol_display_endpoint_property_count")
  string(FIND "${_symbols}" " ${_required}" _required_position)
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
