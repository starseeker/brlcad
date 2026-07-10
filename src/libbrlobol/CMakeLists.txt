set(
  LIBBRLOBOL_SRCS
  adc.cpp
  axes.cpp
  cad_assembly.cpp
  database_source.cpp
  draw_cache.cpp
  edit_preview.cpp
  evaluated_points.cpp
  evaluated_wire.cpp
  export_action.cpp
  framebuffer.cpp
  grid.cpp
  headless_window_host.cpp
  hud_label_overlay.cpp
  image_display_util.cpp
  image_plane.cpp
  image_source.cpp
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
  pick_detail.cpp
  performance.cpp
  realize_action.cpp
  scene_controller.cpp
  scene_group.cpp
  snap_action.cpp
  view_attachment.cpp
  view_controller.cpp
  view_lod.cpp
  view_store.cpp
  viewport_image.cpp
  vlist_shape.cpp
  window_host.cpp
)

if(NOT TARGET Obol::Obol)
  find_package(Obol CONFIG QUIET)
endif()
if(NOT TARGET Obol::Obol AND EXISTS "${CMAKE_BINARY_DIR}/lib/cmake/Obol-1.0.0/obol-config.cmake")
  set(Obol_DIR "${CMAKE_BINARY_DIR}/lib/cmake/Obol-1.0.0" CACHE PATH "Obol package directory" FORCE)
  find_package(Obol CONFIG REQUIRED)
endif()
if(NOT TARGET Obol::Obol)
  message(FATAL_ERROR "libbrlobol requires the Obol::Obol CMake package target")
endif()

set(BRLOBOL_OPENGL_INCLUDE_DIRS)
if(OPENGL_INCLUDE_DIR_GL)
  list(APPEND BRLOBOL_OPENGL_INCLUDE_DIRS ${OPENGL_INCLUDE_DIR_GL})
elseif(OPENGL_INCLUDE_DIR)
  list(APPEND BRLOBOL_OPENGL_INCLUDE_DIRS ${OPENGL_INCLUDE_DIR})
endif()
set(BRLOBOL_LOCAL_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})

brlcad_addlib(libbrlobol "${LIBBRLOBOL_SRCS}" "" "${BRLOBOL_LOCAL_INCLUDE_DIRS}"
  PUBLIC_LIBS "librt;libimgstream;libbv;libbg;libbu;Obol::Obol;Threads::Threads"
  ALL_CXX
)
if(BRLOBOL_OPENGL_INCLUDE_DIRS)
  target_include_directories(libbrlobol SYSTEM PRIVATE ${BRLOBOL_OPENGL_INCLUDE_DIRS})
  if(TARGET libbrlobol-obj)
    target_include_directories(libbrlobol-obj SYSTEM PRIVATE ${BRLOBOL_OPENGL_INCLUDE_DIRS})
  endif(TARGET libbrlobol-obj)
endif(BRLOBOL_OPENGL_INCLUDE_DIRS)
set_target_properties(libbrlobol PROPERTIES VERSION 1.0.0 SOVERSION 1)
if(HIDE_INTERNAL_SYMBOLS)
  set_property(TARGET libbrlobol APPEND PROPERTY COMPILE_DEFINITIONS "BRLOBOL_DLL_EXPORTS")
  if(TARGET libbrlobol-obj)
    set_property(TARGET libbrlobol-obj APPEND PROPERTY COMPILE_DEFINITIONS "BRLOBOL_DLL_EXPORTS")
  endif(TARGET libbrlobol-obj)
endif(HIDE_INTERNAL_SYMBOLS)

add_subdirectory(tests)

cmakefiles(
  CMakeLists.txt
  adc.cpp
  axes.cpp
  cad_assembly.cpp
  cad_assembly_private.h
  database_source.cpp
  database_source_realization.h
  draw_cache.cpp
  edit_preview.cpp
  evaluated_points.cpp
  evaluated_wire.cpp
  export_action.cpp
  framebuffer.cpp
  grid.cpp
  headless_window_host.cpp
  hud_label_overlay.cpp
  image_display_util.cpp
  image_display_util.h
  image_plane.cpp
  image_source.cpp
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
  pick_detail.cpp
  performance.cpp
  performance_private.h
  realize_action.cpp
  scene_controller.cpp
  scene_group.cpp
  snap_action.cpp
  view_attachment.cpp
  view_controller.cpp
  view_lod.cpp
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
