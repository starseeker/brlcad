set(PLUGIN_LIBS libged libBObol libbv libbu)

# Keep this command as a real dynamically loaded plugin.  Besides providing
# the optional faceplate control, it is the build-tested example of the
# installed ged/plugin/obol.h extension boundary available to third parties.
ged_plugin_library(ged-navgizmo
  navgizmo.cpp
  LIBS ${PLUGIN_LIBS}
  PLUGIN_ONLY
)

cmakefiles(
  CMakeLists.txt
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
