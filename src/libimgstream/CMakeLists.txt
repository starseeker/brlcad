set(
  LIBIMGSTREAM_SOURCES
  fb_compat.c
  stream.c
)

# Note - libimgstream_deps is defined by ${BRLCAD_SOURCE_DIR}/src/source_dirs.cmake
brlcad_addlib(libimgstream "${LIBIMGSTREAM_SOURCES}" "" ""
  PRIVATE_LIBS ${libimgstream_deps}
)
set_target_properties(libimgstream PROPERTIES VERSION 20.0.1 SOVERSION 20)

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
