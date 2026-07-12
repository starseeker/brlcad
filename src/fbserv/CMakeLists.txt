brlcad_addexec(fbserv "fbserv.c" "libimgstream;libpkg;libbu")
target_link_libraries(fbserv ${WINSOCK_LIB})

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
