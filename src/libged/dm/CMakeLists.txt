set(PLUGIN_LIBS libged libBObol libbu)

set(DM_SRCS
  dm.cpp
  ert.cpp
  screengrab.cpp
)

ged_plugin_library(ged-dm
  PLUGIN_ONLY
  ${DM_SRCS}
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
