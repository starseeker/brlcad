set(PLUGIN_LIBS libged libimgstream libicv libbu)

ged_plugin_library(ged-png2fb
  PLUGIN_ONLY
  png2fb.cpp
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
