set(
  bv_headers
  display.h
  view.h
)

brlcad_manage_files(bv_headers ${INCLUDE_DIR}/brlcad/bv REQUIRED libbv)

cmakefiles(
  CMakeLists.txt
  ${bv_headers}
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
