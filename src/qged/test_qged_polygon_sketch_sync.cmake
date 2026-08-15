if(NOT DEFINED BRLCAD_ROOT OR NOT DEFINED MGED OR NOT DEFINED QGED OR
    NOT DEFINED SCRIPT)
  message(FATAL_ERROR "BRLCAD_ROOT, MGED, QGED, and SCRIPT are required")
endif()

string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef _suffix)
set(_test_dir "${BRLCAD_ROOT}/Testing/Temporary/qged_polygon_${_suffix}")
set(_database "${_test_dir}/polygon_sync.g")
set(_setup "${_test_dir}/setup.mged")
file(MAKE_DIRECTORY "${_test_dir}")
file(WRITE "${_setup}"
  "in anchor.s sph 0 0 0 1\n"
  "q\n")

execute_process(
  COMMAND "${MGED}" -c "${_database}"
  INPUT_FILE "${_setup}"
  RESULT_VARIABLE _setup_result
  OUTPUT_QUIET
  ERROR_VARIABLE _setup_error)
if(NOT _setup_result EQUAL 0)
  file(REMOVE_RECURSE "${_test_dir}")
  message(FATAL_ERROR "polygon sync fixture setup failed: ${_setup_error}")
endif()

set(_qged_renderer_args --swrast)
if(QUAD)
  list(APPEND _qged_renderer_args --quad)
endif()
set(_qged_environment "BRLCAD_ROOT=${BRLCAD_ROOT}" "QT_QPA_PLATFORM=offscreen")
if(USE_SYSTEM_GL)
  set(_qged_renderer_args)
  if(QUAD)
    list(APPEND _qged_renderer_args --quad)
  endif()
  set(_qged_environment "BRLCAD_ROOT=${BRLCAD_ROOT}")
endif()

execute_process(
  COMMAND "${CMAKE_COMMAND}" -E env
    ${_qged_environment}
    "${QGED}" ${_qged_renderer_args} --test-script "${SCRIPT}" "${_database}"
  RESULT_VARIABLE _qged_result
  OUTPUT_VARIABLE _qged_output
  ERROR_VARIABLE _qged_error)

file(REMOVE_RECURSE "${_test_dir}")
if(NOT _qged_result EQUAL 0)
  message(FATAL_ERROR
    "qged polygon sketch sync replay failed (${_qged_result})\n"
    "stdout:\n${_qged_output}\n"
    "stderr:\n${_qged_error}")
endif()
