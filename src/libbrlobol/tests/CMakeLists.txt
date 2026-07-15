brlcad_addexec(test_brlobol_prototype test_prototype.cpp "libbrlobol;libwdb;libbg;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libbrlobol_prototype COMMAND test_brlobol_prototype)
set_property(TEST libbrlobol_prototype APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})

brlcad_find_package(OSMESA REQUIRED)
brlcad_addexec(test_brlobol_render test_render.cpp "libbrlobol;libbu;Obol::Obol;${OSMESA_LIBRARY}" TEST ALL_CXX)
target_include_directories(test_brlobol_render SYSTEM PRIVATE ${OSMESA_INCLUDE_DIR} ${BRLOBOL_OPENGL_INCLUDE_DIRS})
brlcad_add_test(NAME libbrlobol_render COMMAND test_brlobol_render)

brlcad_addexec(test_brlobol_image_source test_image_source.cpp "libbrlobol;libimgstream;libicv;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libbrlobol_image_source COMMAND test_brlobol_image_source)

brlcad_addexec(test_brlobol_image_display test_image_display.cpp "libbrlobol;libimgstream;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libbrlobol_image_display COMMAND test_brlobol_image_display)

brlcad_addexec(test_brlobol_lod_realization test_lod_realization.cpp "libbrlobol;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libbrlobol_lod_realization COMMAND test_brlobol_lod_realization)

brlcad_addexec(test_brlobol_lod_service test_lod_service.cpp "libbrlobol;libwdb;libbu;Obol::Obol;Threads::Threads" TEST ALL_CXX)
brlcad_add_test(NAME libbrlobol_lod_service COMMAND test_brlobol_lod_service)
set_tests_properties(libbrlobol_lod_service PROPERTIES TIMEOUT 10)

brlcad_addexec(test_brlobol_lod_update_action test_lod_update_action.cpp "libbrlobol;libwdb;libbu;Obol::Obol;Threads::Threads" TEST ALL_CXX)
brlcad_add_test(NAME libbrlobol_lod_update_action COMMAND test_brlobol_lod_update_action)

brlcad_addexec(test_brlobol_draw_cache test_draw_cache.cpp "libbrlobol;libwdb;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libbrlobol_draw_cache COMMAND test_brlobol_draw_cache)

brlcad_addexec(test_brlobol_evaluated_points test_evaluated_points.cpp "libbrlobol;libwdb;librt;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libbrlobol_evaluated_points COMMAND test_brlobol_evaluated_points)

brlcad_addexec(test_brlobol_mesh_lod_cache test_mesh_lod_cache.cpp "libbrlobol;libwdb;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libbrlobol_mesh_lod_cache COMMAND test_brlobol_mesh_lod_cache)

brlcad_addexec(test_brlobol_view_store test_view_store.cpp "libbrlobol;libwdb;librt;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libbrlobol_view_store COMMAND test_brlobol_view_store)

brlcad_addexec(test_brlobol_rt_render test_rt_render.cpp "libbrlobol;libwdb;librt;libbu;Obol::Obol;Threads::Threads" TEST ALL_CXX)
brlcad_add_test(NAME libbrlobol_rt_render COMMAND test_brlobol_rt_render)

brlcad_addexec(test_brlobol_window_host test_window_host.cpp "libbrlobol;libimgstream;libbu;Obol::Obol;Threads::Threads" TEST ALL_CXX)
brlcad_add_test(NAME libbrlobol_window_host COMMAND test_brlobol_window_host)

brlcad_addexec(test_brlobol_headless_window_host test_headless_window_host.cpp "libbrlobol;libimgstream;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libbrlobol_headless_window_host COMMAND test_brlobol_headless_window_host)

if(Obol_NanoRTContextManager_FOUND)
  brlcad_addexec(test_brlobol_nanort_render test_nanort_render.cpp "libbrlobol;libwdb;libbu;Obol::Obol;Threads::Threads;${OSMESA_LIBRARY}" TEST ALL_CXX)
  target_include_directories(test_brlobol_nanort_render SYSTEM PRIVATE
    ${OSMESA_INCLUDE_DIR}
    ${BRLOBOL_OPENGL_INCLUDE_DIRS}
    ${Obol_NanoRTContextManager_INCLUDE_DIR}
  )
  brlcad_add_test(NAME libbrlobol_nanort_render COMMAND test_brlobol_nanort_render)
endif()

cmakefiles(
  CMakeLists.txt
  test_headless_window_host.cpp
  test_draw_cache.cpp
  test_evaluated_points.cpp
  test_image_display.cpp
  test_image_source.cpp
  test_lod_realization.cpp
  test_lod_service.cpp
  test_lod_update_action.cpp
  test_mesh_lod_cache.cpp
  test_nanort_render.cpp
  test_prototype.cpp
  test_rt_render.cpp
  test_render.cpp
  test_view_store.cpp
  test_window_host.cpp
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
