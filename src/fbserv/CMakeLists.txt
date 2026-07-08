include_directories(${DM_INCLUDE_DIRS})

# Optional: link against OpenSSL for TLS support.  Session-token auth helpers
# are provided by libimgstream's fbserv protocol contract.
find_package(OpenSSL QUIET)

brlcad_addexec(fbserv "fbserv.c;server.c" "libdm;libimgstream")
target_link_libraries(fbserv ${WINSOCK_LIB})
add_target_deps(fbserv dm_plugins)
if(OpenSSL_FOUND OR OPENSSL_FOUND)
  target_link_libraries(fbserv OpenSSL::SSL OpenSSL::Crypto)
  target_compile_definitions(fbserv PRIVATE HAVE_OPENSSL_SSL_H=1)
endif()

cmakefiles(
  tls_wrap.h
  CMakeLists.txt
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
