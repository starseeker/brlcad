set(PLUGIN_LIBS libged libbu)
set_source_files_properties(shaded_mode.cpp TARGET_DIRECTORY libged
  PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)

ged_plugin_library(ged-shaded_mode
  shaded_mode.cpp
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
