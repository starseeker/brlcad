set(PLUGIN_LIBS libged libbu)
set_source_files_properties(bigE.cpp draw.cpp TARGET_DIRECTORY libged
  PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)

set(DRAW_SRCS
  bigE.cpp
  draw.cpp
  draw2.cpp
  loadview.cpp
  preview.cpp
)

ged_plugin_library(ged-draw
  ${DRAW_SRCS}
  INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
  LIBS ${PLUGIN_LIBS}
)

cmakefiles(
  CMakeLists.txt
  ged_draw.h
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
