if(NOT DEFINED BRLCAD_ROOT OR NOT DEFINED MGED OR NOT DEFINED QGED OR
    NOT DEFINED SCRIPT)
  message(FATAL_ERROR "BRLCAD_ROOT, MGED, QGED, and SCRIPT are required")
endif()

string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef _suffix)
set(_test_dir "${BRLCAD_ROOT}/Testing/Temporary/qged_edit_${_suffix}")
set(_database "${_test_dir}/primitive_edit.g")
set(_setup "${_test_dir}/setup.mged")
file(MAKE_DIRECTORY "${_test_dir}")
file(WRITE "${_setup}"
  "in ell.s ell 0 0 0 2 0 0 0 1 0 0 0 1\n"
  "in rename_ell.s ell 0 0 0 3 0 0 0 2 0 0 0 1\n"
  "in remove_ell.s ell 0 0 0 4 0 0 0 2 0 0 0 1\n"
  "put sketch.s sketch V {0 0 0} A {1 0 0} B {0 1 0} VL {{1 0} {1 1} {0 1} {0 0}} SL {{line S 0 E 1} {line S 1 E 2} {line S 2 E 3} {line S 3 E 0}}\n"
  "in extrude.s extrude 0 0 0 0 0 2 1 0 0 0 1 0 sketch.s\n"
  "in revolve.s revolve 0 0 0 0 0 1 1 0 0 180 sketch.s\n"
  "in arb.s arb8 0 0 0 2 0 0 2 2 0 0 2 0 0 0 2 2 0 2 2 2 2 0 2 2\n"
  "put bot.s bot mode volume orient rh flags {} V {{0 0 0} {2 0 0} {0 2 0} {0 0 2}} F {{0 2 1} {0 1 3} {0 3 2} {1 2 3}}\n"
  "g remove_group.c remove_ell.s\n"
  "g group.c ell.s extrude.s revolve.s arb.s sketch.s bot.s\n"
  "q\n")

execute_process(
  COMMAND "${MGED}" -c "${_database}"
  INPUT_FILE "${_setup}"
  RESULT_VARIABLE _setup_result
  OUTPUT_QUIET
  ERROR_VARIABLE _setup_error)
if(NOT _setup_result EQUAL 0)
  file(REMOVE_RECURSE "${_test_dir}")
  message(FATAL_ERROR "primitive edit fixture setup failed: ${_setup_error}")
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
    "qged primitive edit replay failed (${_qged_result})\n"
    "stdout:\n${_qged_output}\n"
    "stderr:\n${_qged_error}")
endif()
if(_qged_error MATCHES "db_string_to_path\\(\\) failed")
  message(FATAL_ERROR
    "qged primitive edit replay emitted stale database-path diagnostics\n"
    "stderr:\n${_qged_error}")
endif()
