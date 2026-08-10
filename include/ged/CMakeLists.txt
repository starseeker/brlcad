set(
  ged_headers
  commands.h
  database.h
  db_index.h
  debug.h
  defines.h
  display.h
  draw.h
  draw_types.h
  event_txn.h
  objects.h
  polygon.h
  polygon_types.h
  view_feature_batch.h
  scene.h
  scene_record_types.h
  scene_types.h
  selection.h
  selection_types.h
  selection_state.h
  version.h
  view.h
  view_feature.h
  view_feature_types.h
  view_types.h
)
brlcad_manage_files(ged_headers ${INCLUDE_DIR}/brlcad/ged REQUIRED libged)

# Deliberately not installed: concrete retained-render host adapter.
cmakefiles(
  CMakeLists.txt
  display_obol_private.h
  draw_source.h
  draw_scene.h
  scene_internal.h
  scene_record_api_internal.h
  view_feature_internal.h
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
