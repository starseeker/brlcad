set(PLUGIN_LIBS libged libimgstream libbu)

ged_plugin_library(ged-fbclear
  PLUGIN_ONLY
  fbclear.cpp
  LIBS ${PLUGIN_LIBS}
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
