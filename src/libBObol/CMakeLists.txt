set(
  LIBBOBOL_SRCS
  adc.cpp
  axes.cpp
  cad_assembly.cpp
  cad_normals.cpp
  cad_publication.cpp
  compact_occurrence_stream.cpp
  compact_presentation_staging.cpp
  database_source_cache.cpp
  database_source_compact_access.cpp
  database_source_interaction.cpp
  database_source_mesh_geometry.cpp
  database_source.cpp
  database_source_presentation.cpp
  database_source_summary.cpp
  database_source_types.cpp
  display_endpoint.cpp
  display_session.cpp
  draw_cache.cpp
  draw_cache_private.h
  edit_manipulator.cpp
  edit_preview.cpp
  evaluated_points.cpp
  exact_sample_identity_private.h
  export_action.cpp
  framebuffer.cpp
  grid.cpp
  headless_window_host.cpp
  host_factory.cpp
  hud_label_overlay.cpp
  image_display_util.cpp
  image_plane.cpp
  image_source.cpp
  identity_counter_private.h
  input.cpp
  init.cpp
  line_layer_overlay.cpp
  lod_admission_policy_private.h
  lod_control_private.h
  lod_capacity_policy_private.h
  lod_coordinator_private.h
  lod_delivery_policy_private.h
  lod_presentation_policy_private.h
  lod_population_identity_private.h
  lod_progress_estimator_private.h
  lod_quiet_successor_private.h
  lod_result_authentication_private.h
  lod_scene_evidence_private.h
  lod_terminal_quality_private.h
  lod_visual_importance_private.h
  lod_view_policy_private.h
  lod_view_snapshot_private.h
  retained_allocation_private.cpp
  retained_allocation_private.h
  lod_realization.cpp
  lod_service.cpp
  lod_mesh_shape.cpp
  lod_presentation_policy.cpp
  lod_admission_policy.cpp
  lod_update_action.cpp
  lod_view_quality_history.cpp
  material_object.cpp
  measure_action.cpp
  mesh_lod_cache.cpp
  mesh_lod_submit_action.cpp
  mesh_residency_action.cpp
  mesh_shape.cpp
  navigation_gizmo.cpp
  pick_detail.cpp
  performance.cpp
  realize_action.cpp
  rt_render.cpp
  scene_controller.cpp
  scene_group.cpp
  serialized_bot_source.cpp
  snap_action.cpp
  source_realization.cpp
  transaction_fault_private.h
  view_attachment.cpp
  view_controller.cpp
  view_controller_host.cpp
  view_controller_lighting.cpp
  view_controller_picking.cpp
  view_controller_progressive.cpp
  view_controller_render.cpp
  view_controller_residency.cpp
  view_controller_scene.cpp
  view_controller_private.h
  view_lod_coordinator_state_private.h
  view_lod.cpp
  view_query.cpp
  view_store.cpp
  viewport_image.cpp
  vlist_shape.cpp
  window_host.cpp
)

# Keep responsibility-split database source units independently compilable.
# Unity inclusion can otherwise conceal an accidental private dependency
# between realization/mutation and read-only inspection.
set_source_files_properties(
  cad_normals.cpp
  compact_presentation_staging.cpp
  database_source_cache.cpp
  database_source_compact_access.cpp
  database_source_interaction.cpp
  database_source_mesh_geometry.cpp
  database_source.cpp
  database_source_presentation.cpp
  database_source_summary.cpp
  database_source_types.cpp
  lod_presentation_policy.cpp
  lod_admission_policy.cpp
  lod_view_quality_history.cpp
  serialized_bot_source.cpp
  view_controller.cpp
  view_controller_host.cpp
  view_controller_lighting.cpp
  view_controller_picking.cpp
  view_controller_progressive.cpp
  view_controller_render.cpp
  view_controller_residency.cpp
  view_controller_scene.cpp
  PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
)

find_package(Obol 2.0 CONFIG REQUIRED)
set(BOBOL_OPENGL_INCLUDE_DIRS)
if(OPENGL_INCLUDE_DIR_GL)
  list(APPEND BOBOL_OPENGL_INCLUDE_DIRS ${OPENGL_INCLUDE_DIR_GL})
elseif(OPENGL_INCLUDE_DIR)
  list(APPEND BOBOL_OPENGL_INCLUDE_DIRS ${OPENGL_INCLUDE_DIR})
endif()
set(BOBOL_LOCAL_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})

brlcad_addlib(libBObol "${LIBBOBOL_SRCS}" "" "${BOBOL_LOCAL_INCLUDE_DIRS}"
  PUBLIC_LIBS "librt;libimgstream;libbv;libbg;libbu;Obol::Obol;Threads::Threads"
  ALL_CXX
)

# libBObol's public headers mark the supported ABI explicitly with
# BOBOL_EXPORT.  Keep implementation-only C++ entities (in particular PImpl
# constructors, destructors, and inline helpers) out of the ELF/Mach-O export
# table.  Besides tightening the ABI, this prevents consumers from binding to
# compiler-generated implementation symbols whose shape may change freely.
set_target_properties(libBObol PROPERTIES
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN YES
)
if(TARGET libBObol-obj)
  set_target_properties(libBObol-obj PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN YES
  )
endif(TARGET libBObol-obj)
if(BOBOL_OPENGL_INCLUDE_DIRS)
  target_include_directories(libBObol SYSTEM PRIVATE ${BOBOL_OPENGL_INCLUDE_DIRS})
  if(TARGET libBObol-obj)
    target_include_directories(libBObol-obj SYSTEM PRIVATE ${BOBOL_OPENGL_INCLUDE_DIRS})
  endif(TARGET libBObol-obj)
endif(BOBOL_OPENGL_INCLUDE_DIRS)
set_target_properties(libBObol PROPERTIES VERSION 1.0.0 SOVERSION 1)
if(HIDE_INTERNAL_SYMBOLS)
  set_property(TARGET libBObol APPEND PROPERTY COMPILE_DEFINITIONS "BOBOL_DLL_EXPORTS")
  if(TARGET libBObol-obj)
    set_property(TARGET libBObol-obj APPEND PROPERTY COMPILE_DEFINITIONS "BOBOL_DLL_EXPORTS")
  endif(TARGET libBObol-obj)
endif(HIDE_INTERNAL_SYMBOLS)

include(${CMAKE_SOURCE_DIR}/include/BObol/api_tiers.cmake)
add_subdirectory(tests)

cmakefiles(
  CMakeLists.txt
  check_public_symbols.cmake
  adc.cpp
  axes.cpp
  cad_assembly.cpp
  cad_assembly_private.h
  cad_normals_private.h
  cad_normals.cpp
  cad_publication.cpp
  cad_publication_private.h
  compact_occurrence_stream.cpp
  compact_occurrence_registry_private.h
  compact_presentation_staging.cpp
  compact_presentation_staging_private.h
  database_source_cache.cpp
  database_source_compact_access.cpp
  database_source_interaction.cpp
  database_source_mesh_geometry.cpp
  database_source_mesh_geometry_private.h
  database_source.cpp
  database_source_presentation.cpp
  database_source_presentation_private.h
  database_source_private.h
  database_source_realization.h
  database_source_summary.cpp
  database_source_types.cpp
  display_endpoint.cpp
  display_session.cpp
  draw_cache.cpp
  draw_cache_private.h
  edit_manipulator.cpp
  edit_preview.cpp
  evaluated_points.cpp
  exact_sample_identity_private.h
  export_action.cpp
  framebuffer.cpp
  grid.cpp
  headless_window_host.cpp
  host_factory.cpp
  hud_label_overlay.cpp
  image_display_util.cpp
  image_display_util.h
  image_plane.cpp
  image_source.cpp
  identity_counter_private.h
  input.cpp
  init.cpp
  line_layer_overlay.cpp
  lod_admission_policy_private.h
  lod_control_private.h
  lod_capacity_policy_private.h
  lod_coordinator_private.h
  lod_delivery_policy_private.h
  lod_presentation_policy_private.h
  lod_population_identity_private.h
  lod_progress_estimator_private.h
  lod_presentation_policy.cpp
  lod_admission_policy.cpp
  lod_quiet_successor_private.h
  lod_result_authentication_private.h
  lod_scene_evidence_private.h
  lod_terminal_quality_private.h
  lod_visual_importance_private.h
  lod_view_policy_private.h
  lod_view_snapshot_private.h
  lod_realization.cpp
  lod_service.cpp
  lod_mesh_shape.cpp
  lod_update_action.cpp
  lod_view_quality_history.cpp
  material_object.cpp
  measure_action.cpp
  mesh_lod_cache.cpp
  mesh_lod_submit_action.cpp
  mesh_residency_action.cpp
  mesh_shape.cpp
  navigation_gizmo.cpp
  pick_detail.cpp
  performance.cpp
  performance_private.h
  retained_allocation_private.cpp
  retained_allocation_private.h
  realize_action.cpp
  scene_controller.cpp
  scene_group.cpp
  snap_action.cpp
  source_realization.cpp
  transaction_fault_private.h
  view_attachment.cpp
  view_controller.cpp
  view_controller_host.cpp
  view_controller_lighting.cpp
  view_controller_picking.cpp
  view_controller_private.h
  view_controller_progressive.cpp
  view_controller_render.cpp
  view_controller_residency.cpp
  view_controller_scene.cpp
  view_lod_coordinator_state_private.h
  view_lod.cpp
  view_query.cpp
  view_store.cpp
  viewport_image.cpp
  vlist_shape.cpp
  window_host.cpp
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
