set(
  brlobol_headers
  adc.h
  axes.h
  database_source.h
  defines.h
  draw_cache.h
  edit_preview.h
  evaluated_wire.h
  export_action.h
  framebuffer.h
  grid.h
  headless_window_host.h
  hud_label_overlay.h
  image_plane.h
  image_source.h
  init.h
  line_layer_overlay.h
  lod_realization.h
  lod_service.h
  lod_mesh_shape.h
  lod_update_action.h
  material_object.h
  measure_action.h
  mesh_lod_cache.h
  mesh_lod_submit_action.h
  mesh_residency_action.h
  mesh_shape.h
  pick_detail.h
  performance.h
  realize_action.h
  scene_controller.h
  scene_group.h
  snap_action.h
  source_mesh_request.h
  view_controller.h
  view_lod.h
  view_store.h
  viewport_image.h
  vlist_shape.h
  window_host.h
)
brlcad_manage_files(brlobol_headers ${INCLUDE_DIR}/brlcad/brlobol REQUIRED libbrlobol)

cmakefiles(CMakeLists.txt)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
