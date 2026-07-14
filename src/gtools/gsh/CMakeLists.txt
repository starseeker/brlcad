brlcad_addexec(gsh "gsh.cpp" "libged;libbrlobol;libimgstream;libbu;Threads::Threads")
target_include_directories(gsh BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../libtermio)
add_dependencies(gsh ged_plugins)
add_dependencies(gsh ged_subprocesses ged_exec)
if(TARGET rt)
  add_dependencies(gsh rt)
endif(TARGET rt)

if(SH_EXEC AND TARGET gsh)
  brlcad_add_test(NAME gsh_obol_smoke
    COMMAND ${SH_EXEC} "${CMAKE_CURRENT_SOURCE_DIR}/gsh_obol_smoke.sh"
    "$<TARGET_FILE:gsh>" "${CMAKE_SOURCE_DIR}/src/libged/tests/draw/moss.g"
    "${CMAKE_CURRENT_BINARY_DIR}")
  set_property(TEST gsh_obol_smoke APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})

  brlcad_add_test(NAME gsh_obol_ert_smoke
    COMMAND ${SH_EXEC} "${CMAKE_CURRENT_SOURCE_DIR}/gsh_obol_ert_smoke.sh"
    "$<TARGET_FILE:gsh>" "${CMAKE_SOURCE_DIR}/src/libged/tests/draw/moss.g"
    "${CMAKE_CURRENT_BINARY_DIR}")
  set_tests_properties(gsh_obol_ert_smoke PROPERTIES TIMEOUT 90)
  set_property(TEST gsh_obol_ert_smoke APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})

  brlcad_add_test(NAME gsh_obol_progressive_lod
    COMMAND ${SH_EXEC} "${CMAKE_CURRENT_SOURCE_DIR}/gsh_obol_progressive_lod.sh"
    "$<TARGET_FILE:gsh>" "${CMAKE_SOURCE_DIR}/src/libged/tests/draw/moss.g"
    "${CMAKE_CURRENT_BINARY_DIR}")
  set_tests_properties(gsh_obol_progressive_lod PROPERTIES TIMEOUT 90)
  set_property(TEST gsh_obol_progressive_lod APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
endif(SH_EXEC AND TARGET gsh)

cmakefiles(
  CMakeLists.txt
  gsh_obol_ert_smoke.sh
  gsh_obol_progressive_lod.sh
  gsh_obol_smoke.sh
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
