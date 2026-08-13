if(NOT DEFINED JAVA_EXECUTABLE OR NOT EXISTS "${JAVA_EXECUTABLE}")
  message(FATAL_ERROR "JAVA_EXECUTABLE does not name an available Java runtime")
endif()

if(NOT DEFINED TLA2TOOLS_JAR OR NOT EXISTS "${TLA2TOOLS_JAR}")
  message(FATAL_ERROR "TLA2TOOLS_JAR does not name an available tla2tools.jar")
endif()

if(NOT DEFINED SPEC_DIR OR NOT IS_DIRECTORY "${SPEC_DIR}")
  message(FATAL_ERROR "SPEC_DIR does not name the TLA+ specification directory")
endif()

if(NOT DEFINED MODULE OR MODULE STREQUAL "")
  message(FATAL_ERROR "MODULE is required")
endif()

if(NOT DEFINED CONFIG OR CONFIG STREQUAL "")
  message(FATAL_ERROR "CONFIG is required")
endif()

if(NOT EXISTS "${SPEC_DIR}/${MODULE}.tla")
  message(FATAL_ERROR "TLA+ module not found: ${SPEC_DIR}/${MODULE}.tla")
endif()

if(NOT EXISTS "${SPEC_DIR}/${CONFIG}")
  message(FATAL_ERROR "TLC configuration not found: ${SPEC_DIR}/${CONFIG}")
endif()

file(SHA256 "${TLA2TOOLS_JAR}" _tla2tools_actual_sha256)
if(DEFINED TLA2TOOLS_SHA256
   AND NOT TLA2TOOLS_SHA256 STREQUAL ""
   AND NOT _tla2tools_actual_sha256 STREQUAL TLA2TOOLS_SHA256)
  message(
    FATAL_ERROR
    "tla2tools.jar SHA-256 mismatch: expected ${TLA2TOOLS_SHA256}, got ${_tla2tools_actual_sha256}")
endif()

if(NOT DEFINED TLC_WORKERS OR TLC_WORKERS STREQUAL "")
  set(TLC_WORKERS 1)
endif()
if(NOT DEFINED TLC_HEAP OR TLC_HEAP STREQUAL "")
  set(TLC_HEAP 2g)
endif()

execute_process(
  COMMAND
    "${JAVA_EXECUTABLE}" "-Xmx${TLC_HEAP}" -cp "${TLA2TOOLS_JAR}"
    tlc2.TLC -cleanup -workers "${TLC_WORKERS}" -config "${CONFIG}" "${MODULE}"
  WORKING_DIRECTORY "${SPEC_DIR}"
  RESULT_VARIABLE _tlc_result
  OUTPUT_VARIABLE _tlc_stdout
  ERROR_VARIABLE _tlc_stderr
  TIMEOUT 300
)

set(_tlc_output "${_tlc_stdout}\n${_tlc_stderr}")

if(DEFINED EXPECT_FAILURE AND EXPECT_FAILURE)
  if(_tlc_result EQUAL 0)
    message(FATAL_ERROR "TLC mutation test unexpectedly passed:\n${_tlc_output}")
  endif()
  if(DEFINED EXPECTED_VIOLATION AND NOT EXPECTED_VIOLATION STREQUAL "")
    string(FIND "${_tlc_output}" "${EXPECTED_VIOLATION}" _violation_index)
    if(_violation_index EQUAL -1)
      message(
        FATAL_ERROR
        "TLC failed, but not for expected property ${EXPECTED_VIOLATION}:\n${_tlc_output}")
    endif()
  endif()
  message(
    STATUS
    "TLC rejected ${MODULE}/${CONFIG} as expected (${EXPECTED_VIOLATION}); jar SHA-256 ${_tla2tools_actual_sha256}")
  return()
endif()

if(NOT _tlc_result EQUAL 0)
  message(FATAL_ERROR "TLC failed for ${MODULE}/${CONFIG}:\n${_tlc_output}")
endif()

message(
  STATUS
  "TLC passed ${MODULE}/${CONFIG}; jar SHA-256 ${_tla2tools_actual_sha256}\n${_tlc_output}")
