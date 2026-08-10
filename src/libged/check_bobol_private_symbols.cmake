# The command-local C++ integration helpers are deliberately not part of the
# libged ABI.  HIDE_INTERNAL_SYMBOLS gives unannotated definitions hidden ELF
# visibility; this check prevents a future GED_EXPORT from silently undoing
# that boundary.

if(NOT DEFINED NM OR NOT EXISTS "${NM}")
  message(FATAL_ERROR "libged BObol symbol contract: nm executable is unavailable")
endif()
if(NOT DEFINED LIBRARY OR NOT EXISTS "${LIBRARY}")
  message(FATAL_ERROR "libged BObol symbol contract: library is unavailable: ${LIBRARY}")
endif()

execute_process(
  COMMAND "${NM}" -D -C --defined-only "${LIBRARY}"
  RESULT_VARIABLE _nm_result
  OUTPUT_VARIABLE _symbols
  ERROR_VARIABLE _nm_error)
if(NOT _nm_result EQUAL 0)
  message(FATAL_ERROR "libged BObol symbol contract: nm failed: ${_nm_error}")
endif()

foreach(_private_prefix
    "ged_bobol_"
    "ged_draw_obol_source_"
    "ged_draw_obol_controller"
    "ged_draw_obol_scene_controller"
    "ged_draw_obol_database_source_"
    "ged_draw_obol_database_sources_"
    "ged_draw_obol_view_context_"
    "ged_draw_obol_local_shape_")
  string(FIND "${_symbols}" "${_private_prefix}" _private_position)
  if(NOT _private_position EQUAL -1)
    message(FATAL_ERROR
      "libged BObol symbol contract: private helper prefix ${_private_prefix} was exported")
  endif()
endforeach()
