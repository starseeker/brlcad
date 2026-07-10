set(
  imgstream_headers
  fb_compat.h
  fbserv.h
  stream.h
)
brlcad_manage_files(imgstream_headers ${INCLUDE_DIR}/brlcad/imgstream REQUIRED libimgstream)

cmakefiles(CMakeLists.txt)
cmakefiles(
  fbserv.h
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
