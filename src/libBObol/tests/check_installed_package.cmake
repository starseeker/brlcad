if(NOT DEFINED BUILD_DIR OR NOT DEFINED SOURCE_DIR OR NOT DEFINED STAGE_DIR)
  message(FATAL_ERROR "installed package check requires build, source, and stage directories")
endif()

file(REMOVE_RECURSE "${STAGE_DIR}")
set(_install_command "${CMAKE_COMMAND}" --install "${BUILD_DIR}"
  --prefix "${STAGE_DIR}")
if(DEFINED CONFIG AND NOT "${CONFIG}" STREQUAL "")
  list(APPEND _install_command --config "${CONFIG}")
endif()
execute_process(
  COMMAND ${_install_command}
  RESULT_VARIABLE _install_result
  OUTPUT_VARIABLE _install_output
  ERROR_VARIABLE _install_error)
if(NOT _install_result EQUAL 0)
  message(FATAL_ERROR
    "libBObol install staging failed:\n${_install_output}\n${_install_error}")
endif()

set(_consumer_build "${STAGE_DIR}/consumer-build")
execute_process(
  COMMAND "${CMAKE_COMMAND}" -S "${SOURCE_DIR}" -B "${_consumer_build}"
    "-DCMAKE_PREFIX_PATH=${STAGE_DIR}"
    "-DCMAKE_BUILD_TYPE=Release"
  RESULT_VARIABLE _configure_result
  OUTPUT_VARIABLE _configure_output
  ERROR_VARIABLE _configure_error)
if(NOT _configure_result EQUAL 0)
  message(FATAL_ERROR
    "installed libBObol consumer configure failed:\n${_configure_output}\n${_configure_error}")
endif()

execute_process(
  COMMAND "${CMAKE_COMMAND}" --build "${_consumer_build}"
  RESULT_VARIABLE _build_result
  OUTPUT_VARIABLE _build_output
  ERROR_VARIABLE _build_error)
if(NOT _build_result EQUAL 0)
  message(FATAL_ERROR
    "installed libBObol consumer build failed:\n${_build_output}\n${_build_error}")
endif()
