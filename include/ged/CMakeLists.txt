set(
  ged_headers
  commands.h
  database.h
  db_index.h
  debug.h
  defines.h
  display.h
  draw.h
  edit.h
  event.h
  objects.h
  view_polygon.h
  view_polygon_types.h
  view_edit.h
  view_edit_types.h
  view_export.h
  view_export_types.h
  view_feature_batch.h
  scene.h
  scene_types.h
  selection.h
  selection_types.h
  version.h
  view.h
  view_feature.h
  view_feature_types.h
  view_types.h
)
brlcad_manage_files(ged_headers ${INCLUDE_DIR}/brlcad/ged REQUIRED libged)
add_subdirectory(plugin)

# Deliberately not installed: concrete retained-render host adapter.
cmakefiles(
  CMakeLists.txt
  display_obol_private.h
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
