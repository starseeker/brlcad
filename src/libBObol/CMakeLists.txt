set(
  LIBBOBOL_SRCS
  adc.cpp
  axes.cpp
  cad_assembly.cpp
  compact_occurrence_stream.cpp
  database_source_cache.cpp
  database_source_compact_access.cpp
  database_source_interaction.cpp
  database_source.cpp
  database_source_summary.cpp
  database_source_types.cpp
  display_endpoint.cpp
  display_session.cpp
  draw_cache.cpp
  edit_manipulator.cpp
  edit_preview.cpp
  evaluated_points.cpp
  export_action.cpp
  framebuffer.cpp
  grid.cpp
  headless_window_host.cpp
  host_factory.cpp
  hud_label_overlay.cpp
  image_display_util.cpp
  image_plane.cpp
  image_source.cpp
  input.cpp
  init.cpp
  line_layer_overlay.cpp
  lod_coordinator_private.h
  retained_allocation_private.cpp
  retained_allocation_private.h
  lod_realization.cpp
  lod_service.cpp
  lod_mesh_shape.cpp
  lod_update_action.cpp
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
  snap_action.cpp
  source_realization.cpp
  view_attachment.cpp
  view_controller.cpp
  view_controller_lighting.cpp
  view_controller_progressive.cpp
  view_controller_render.cpp
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
  database_source_cache.cpp
  database_source_compact_access.cpp
  database_source_interaction.cpp
  database_source.cpp
  database_source_summary.cpp
  database_source_types.cpp
  view_controller.cpp
  view_controller_lighting.cpp
  view_controller_progressive.cpp
  view_controller_render.cpp
  view_controller_scene.cpp
  PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
)

find_package(Obol CONFIG REQUIRED
  PATHS
  ${PROJECT_BINARY_DIR}/lib/cmake/Obol-1.0.0
  ${PROJECT_BINARY_DIR}/lib64/cmake/Obol-1.0.0
)
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
  compact_occurrence_stream.cpp
  compact_occurrence_registry_private.h
  database_source_cache.cpp
  database_source_compact_access.cpp
  database_source_interaction.cpp
  database_source.cpp
  database_source_private.h
  database_source_realization.h
  database_source_summary.cpp
  database_source_types.cpp
  display_endpoint.cpp
  display_session.cpp
  draw_cache.cpp
  edit_manipulator.cpp
  edit_preview.cpp
  evaluated_points.cpp
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
  input.cpp
  init.cpp
  line_layer_overlay.cpp
  lod_realization.cpp
  lod_service.cpp
  lod_mesh_shape.cpp
  lod_update_action.cpp
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
  view_attachment.cpp
  view_controller.cpp
  view_controller_lighting.cpp
  view_controller_private.h
  view_controller_progressive.cpp
  view_controller_render.cpp
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
