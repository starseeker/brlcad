include_directories(${DM_INCLUDE_DIRS})

brlcad_addexec(fbserv "fbserv.c" "libdm;libimgstream;libpkg")
target_link_libraries(fbserv ${WINSOCK_LIB})
add_target_deps(fbserv dm_plugins)

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
