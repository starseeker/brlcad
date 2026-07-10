set(
  LIBIMGSTREAM_SOURCES
  fb_compat.c
  fbserv_auth.c
  fbserv_framebuffer.c
  fbserv_pkg.c
  fbserv_server.c
  fbserv_state.c
  stream.c
)

find_package(OpenSSL QUIET)
set(LIBIMGSTREAM_PRIVATE_LIBS ${libimgstream_deps})
if(OpenSSL_FOUND OR OPENSSL_FOUND)
  list(APPEND LIBIMGSTREAM_PRIVATE_LIBS OpenSSL::SSL OpenSSL::Crypto)
endif()

# Note - libimgstream_deps is defined by ${BRLCAD_SOURCE_DIR}/src/source_dirs.cmake
brlcad_addlib(libimgstream "${LIBIMGSTREAM_SOURCES}" "" ""
  PRIVATE_LIBS ${LIBIMGSTREAM_PRIVATE_LIBS}
)
set_target_properties(libimgstream PROPERTIES VERSION 20.0.1 SOVERSION 20)

if(OpenSSL_FOUND OR OPENSSL_FOUND)
  target_compile_definitions(libimgstream PRIVATE HAVE_OPENSSL_SSL_H=1)
  if(TARGET libimgstream-obj)
    target_compile_definitions(libimgstream-obj PRIVATE HAVE_OPENSSL_SSL_H=1)
  endif()
  if(TARGET libimgstream-brlcad-obj)
    target_compile_definitions(libimgstream-brlcad-obj PRIVATE HAVE_OPENSSL_SSL_H=1)
    brlcad_register_aggregate_link_libs(OpenSSL::SSL OpenSSL::Crypto)
  endif()
endif()

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
