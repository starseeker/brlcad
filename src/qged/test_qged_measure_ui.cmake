if(NOT DEFINED BRLCAD_ROOT OR NOT DEFINED MGED OR NOT DEFINED QGED OR
    NOT DEFINED SCRIPT)
  message(FATAL_ERROR "BRLCAD_ROOT, MGED, QGED, and SCRIPT are required")
endif()

string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef _suffix)
set(_test_dir "${BRLCAD_ROOT}/Testing/Temporary/qged_measure_${_suffix}")
set(_database "${_test_dir}/measure.g")
set(_setup "${_test_dir}/setup.mged")
file(MAKE_DIRECTORY "${_test_dir}")
file(WRITE "${_setup}"
  "put measure_box.s bot mode volume orient rh flags {} V {{-2 -2 -1} {2 -2 -1} {2 2 -1} {-2 2 -1} {-2 -2 1} {2 -2 1} {2 2 1} {-2 2 1}} F {{0 2 1} {0 3 2} {4 5 6} {4 6 7} {0 1 5} {0 5 4} {1 2 6} {1 6 5} {2 3 7} {2 7 6} {3 0 4} {3 4 7}}\n"
  "r measure.r u measure_box.s\n"
  "q\n")

execute_process(
  COMMAND "${MGED}" -c "${_database}"
  INPUT_FILE "${_setup}"
  RESULT_VARIABLE _setup_result
  OUTPUT_QUIET
  ERROR_VARIABLE _setup_error)
if(NOT _setup_result EQUAL 0)
  file(REMOVE_RECURSE "${_test_dir}")
  message(FATAL_ERROR "measure fixture setup failed: ${_setup_error}")
endif()

set(_qged_args)
set(_qged_environment "BRLCAD_ROOT=${BRLCAD_ROOT}")
if(NOT USE_SYSTEM_GL)
  list(APPEND _qged_args --swrast)
  list(APPEND _qged_environment "QT_QPA_PLATFORM=offscreen")
endif()
if(QUAD)
  list(APPEND _qged_args --quad)
endif()
set(_run_dir "${_test_dir}")
if(DEFINED ARTIFACT_DIR AND NOT ARTIFACT_DIR STREQUAL "")
  get_filename_component(_run_dir "${ARTIFACT_DIR}" ABSOLUTE)
  file(MAKE_DIRECTORY "${_run_dir}")
endif()
execute_process(
  COMMAND "${CMAKE_COMMAND}" -E env
    ${_qged_environment}
    "${QGED}" ${_qged_args} --test-script "${SCRIPT}" "${_database}"
  WORKING_DIRECTORY "${_run_dir}"
  RESULT_VARIABLE _qged_result
  OUTPUT_VARIABLE _qged_output
  ERROR_VARIABLE _qged_error)

if(NOT _qged_result EQUAL 0)
  file(REMOVE_RECURSE "${_test_dir}")
  message(FATAL_ERROR
    "qged measure replay failed (${_qged_result})\n"
    "stdout:\n${_qged_output}\n"
    "stderr:\n${_qged_error}")
endif()

foreach(_image
    measurement-base.png
    measurement-2d.png
    measurement-cleared.png
    measurement-resized-base.png
    measurement-3d.png)
  if(NOT EXISTS "${_run_dir}/${_image}")
    file(REMOVE_RECURSE "${_test_dir}")
    message(FATAL_ERROR "qged measure replay did not create ${_image}")
  endif()
endforeach()
file(REMOVE_RECURSE "${_test_dir}")
