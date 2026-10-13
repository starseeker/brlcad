include_guard(GLOBAL)

set(FLTKCAD_TEST_DATABASE "${BRLCAD_SOURCE_DIR}/src/libged/tests/draw/moss.g"
  CACHE FILEPATH "BRL-CAD database used by FLTK integration tests")

if(BRLCAD_X11_PROVIDER_RESOLVED STREQUAL "XMIN")
  get_filename_component(_fltkcad_xmin_prefix "${Xmin_DIR}/../../.." ABSOLUTE)
  find_program(FLTKCAD_XMIN_RUN NAMES xmin-run HINTS "${_fltkcad_xmin_prefix}/bin")
  find_program(FLTKCAD_XMINCTL NAMES xminctl HINTS "${_fltkcad_xmin_prefix}/bin")
  unset(_fltkcad_xmin_prefix)
endif()
mark_as_advanced(FLTKCAD_TEST_DATABASE FLTKCAD_XMIN_RUN FLTKCAD_XMINCTL)

function(fltkcad_add_gui_test name)
  if(BRLCAD_X11_PROVIDER_RESOLVED STREQUAL "XMIN" AND FLTKCAD_XMIN_RUN)
    add_test(NAME "${name}" COMMAND "${FLTKCAD_XMIN_RUN}" -- ${ARGN})
  else()
    add_test(NAME "${name}" COMMAND ${ARGN})
  endif()
  set_tests_properties("${name}" PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 30)
endfunction()
