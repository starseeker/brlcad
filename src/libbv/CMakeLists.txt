set(
  LIBBV_SRCS
  view.c
)

set(BV_LOCAL_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})

# Note - libbv_deps is defined by ${BRLCAD_SOURCE_DIR}/src/source_dirs.cmake
brlcad_addlib(libbv "${LIBBV_SRCS}" "" "${BV_LOCAL_INCLUDE_DIRS}"
  PUBLIC_LIBS ${libbv_deps}
  PRIVATE_LIBS ${M_LIBRARY}
)
set_target_properties(libbv PROPERTIES VERSION 1.0.0 SOVERSION 1)
if(HIDE_INTERNAL_SYMBOLS)
  set_property(TARGET libbv APPEND PROPERTY COMPILE_DEFINITIONS "BV_DLL_EXPORTS")
  if(TARGET libbv-obj)
    set_property(TARGET libbv-obj APPEND PROPERTY COMPILE_DEFINITIONS "BV_DLL_EXPORTS")
  endif(TARGET libbv-obj)
endif(HIDE_INTERNAL_SYMBOLS)

add_subdirectory(tests)

cmakefiles(
  CMakeLists.txt
  view.c
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
