set(PLUGIN_LIBS libged libbu)

set(
  VIEW_SRCS
  aet.c
  align.c
  autoview.c
  autoview2.cpp
  axes.cpp
  center.cpp
  eye.c
  faceplate/faceplate.cpp
  faceplate/faceplate_axes.cpp
  faceplate/faceplate_grid.cpp
  faceplate/interactive_rect.cpp
  feature_resolver.cpp
  knob.c
  labels.cpp
  lines.cpp
  lod.cpp
  lookat.c
  objs.cpp
  polygons.cpp
  quat.c
  qvrot.c
  saveview.cpp
  size.c
  snap.c
  view.cpp
  viewdir.c
  ypr.c
)
set_source_files_properties(saveview.cpp TARGET_DIRECTORY libged
  PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)

ged_plugin_library(ged-view
  ${VIEW_SRCS}
  INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
  LIBS ${PLUGIN_LIBS}
)

cmakefiles(
  CMakeLists.txt
  feature_resolver.hpp
  ged_view.h
  faceplate/faceplate.h
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
