set(PLUGIN_LIBS libged libbu)

set(WHO_SRCS
  who.cpp
  who2.cpp
  who_solids.cpp
)
set_source_files_properties(who.cpp TARGET_DIRECTORY libged
  PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)

ged_plugin_library(ged-who
  ${WHO_SRCS}
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
