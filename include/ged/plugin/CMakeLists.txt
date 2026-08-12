set(
  ged_plugin_headers
  obol.h
)
brlcad_manage_files(ged_plugin_headers ${INCLUDE_DIR}/brlcad/ged/plugin REQUIRED libged)

cmakefiles(CMakeLists.txt)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
