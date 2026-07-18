set(PLUGIN_LIBS libged libbu)
set_source_files_properties(lod.cpp TARGET_DIRECTORY libged
  PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)

set(LOD_SRCS
  lod.cpp
  )

ged_plugin_library(ged-lod
  ${LOD_SRCS}
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
