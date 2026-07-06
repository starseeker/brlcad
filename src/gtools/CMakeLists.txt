brlcad_adddata(gtransfer.c sample_applications)

brlcad_addexec(ganalyze ganalyze.cpp "libanalyze;librt;libbu;libpkg;${LINENOISE_LIBRARIES}" NO_INSTALL)
target_include_directories(ganalyze BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../libtermio)

brlcad_addexec(gchecker gchecker.cpp "libged;librt;libbu")
if (TARGET ged-check)
  add_dependencies(gchecker ged-check)
endif()
brlcad_addexec(ged_draw_perf "ged_draw_perf.cpp" "libged;libdm;libbrlobol")
if(TARGET ged_draw_perf)
  foreach(tdp dm-obol ged-draw ged-dm ged-view)
    if (TARGET ${tdp})
      add_dependencies(ged_draw_perf ${tdp})
    endif ()
  endforeach()
endif(TARGET ged_draw_perf)
brlcad_addexec(gex "gex.cpp" "libbu;librt")
brlcad_addexec(glint glint.cpp "librt;${M_LIBRARY}")
brlcad_addexec(gqa gqa.c "libged;librt")

brlcad_addexec(gtransfer gtransfer.c "librt;libpkg" NO_INSTALL)

add_subdirectory(beset)
add_subdirectory(gdiff)
add_subdirectory(gist)
add_subdirectory(gsh)
add_subdirectory(tests)

cmakefiles(
  CMakeLists.txt
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
