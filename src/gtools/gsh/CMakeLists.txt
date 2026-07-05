brlcad_addexec(gsh "gsh.cpp" "libged;libdm;libbu;Threads::Threads")
target_include_directories(gsh BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../libtermio)
add_dependencies(gsh dm_plugins ged_plugins)

if(SH_EXEC AND TARGET gsh)
  brlcad_add_test(NAME gsh_obol_smoke
    COMMAND ${SH_EXEC} "${CMAKE_CURRENT_SOURCE_DIR}/gsh_obol_smoke.sh"
    "$<TARGET_FILE:gsh>" "${CMAKE_SOURCE_DIR}/src/libged/tests/draw/moss.g"
    "${CMAKE_CURRENT_BINARY_DIR}")
  set_property(TEST gsh_obol_smoke APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})

  brlcad_add_test(NAME gsh_obol_progressive_lod
    COMMAND ${SH_EXEC} "${CMAKE_CURRENT_SOURCE_DIR}/gsh_obol_progressive_lod.sh"
    "$<TARGET_FILE:gsh>" "${CMAKE_SOURCE_DIR}/src/libged/tests/draw/moss.g"
    "${CMAKE_CURRENT_BINARY_DIR}")
  set_tests_properties(gsh_obol_progressive_lod PROPERTIES TIMEOUT 90)
  set_property(TEST gsh_obol_progressive_lod APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
endif(SH_EXEC AND TARGET gsh)

cmakefiles(
  CMakeLists.txt
  gsh_obol_progressive_lod.sh
  gsh_obol_smoke.sh
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
