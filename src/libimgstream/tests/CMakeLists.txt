brlcad_addexec(imgstream_test stream.c "libimgstream;libicv;libbu" TEST)
brlcad_add_test(NAME imgstream_test COMMAND imgstream_test)

brlcad_addexec(imgstream_fb_compat_test fb_compat.c "libimgstream;libicv;libbu" TEST)
brlcad_add_test(NAME imgstream_fb_compat_test COMMAND imgstream_fb_compat_test)

brlcad_addexec(imgstream_fbserv_auth_test fbserv_auth.c "libimgstream;libpkg;libbu;${WINSOCK_LIB}" TEST)
brlcad_add_test(NAME imgstream_fbserv_auth_test COMMAND imgstream_fbserv_auth_test)

brlcad_addexec(imgstream_fbserv_lifecycle_test fbserv_lifecycle.c "libimgstream;libpkg;libbu" TEST)
brlcad_add_test(NAME imgstream_fbserv_lifecycle_test COMMAND imgstream_fbserv_lifecycle_test)

if(OpenSSL_FOUND OR OPENSSL_FOUND)
  brlcad_addexec(imgstream_fbserv_tls_test fbserv_tls.cpp "libimgstream;libpkg;libbu;Threads::Threads" TEST ALL_CXX)
  brlcad_add_test(NAME imgstream_fbserv_tls_test COMMAND imgstream_fbserv_tls_test)
  set_tests_properties(imgstream_fbserv_tls_test PROPERTIES TIMEOUT 20)
endif()

cmakefiles(CMakeLists.txt fbserv_tls.cpp)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
