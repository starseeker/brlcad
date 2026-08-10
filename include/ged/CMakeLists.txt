set(
  ged_headers
  commands.h
  database.h
  db_index.h
  debug.h
  defines.h
  display.h
  draw.h
  draw_source.h
  draw_types.h
  draw_scene.h
  event_txn.h
  objects.h
  polygon.h
  result.h
  scene.h
  selection.h
  selection_state.h
  version.h
  view.h
  view_feature.h
  view_types.h
)
brlcad_manage_files(ged_headers ${INCLUDE_DIR}/brlcad/ged REQUIRED libged)

# Deliberately not installed: concrete retained-render host adapter.
cmakefiles(CMakeLists.txt display_obol_private.h)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
