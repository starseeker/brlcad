if(NOT WIN32)
  brlcad_addexec(remrt "../rt/opt.c;ihost.c;remrt.cpp" "liboptical;libdm")
  target_include_directories(remrt BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_features(remrt PRIVATE cxx_std_17)
  add_target_deps(remrt dm_plugins)

  brlcad_addexec(rtsrv "../rt/usage.cpp;../rt/view.c;../rt/do.c;../rt/grid.c;../rt/heatgraph.c;../rt/opt.c;../rt/scanline.c;../rt/worker.c;rtsrv.cpp" "libdm;liboptical;libpkg;libicv")
  target_include_directories(rtsrv BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_features(rtsrv PRIVATE cxx_std_17)
  set_property(TARGET rtsrv APPEND PROPERTY COMPILE_DEFINITIONS "RTSRV")
  add_target_deps(rtsrv dm_plugins)
endif(NOT WIN32)

cmakefiles(
  ihost.c
  ihost.h
  protocol.h
  remrt.c
  remrt.cpp
  rtsrv.c
  rtsrv.cpp
)
cmakefiles(CMakeLists.txt)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
