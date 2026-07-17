brlcad_addexec(test_bobol_prototype test_prototype.cpp "libBObol;libwdb;libbg;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libBObol_prototype COMMAND test_bobol_prototype)
set_property(TEST libBObol_prototype APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})

brlcad_find_package(OSMESA REQUIRED)
brlcad_addexec(test_bobol_render test_render.cpp "libBObol;libbu;Obol::Obol;${OSMESA_LIBRARY}" TEST ALL_CXX)
target_include_directories(test_bobol_render SYSTEM PRIVATE ${OSMESA_INCLUDE_DIR} ${BOBOL_OPENGL_INCLUDE_DIRS})
brlcad_add_test(NAME libBObol_render COMMAND test_bobol_render)

brlcad_addexec(test_bobol_image_source test_image_source.cpp "libBObol;libimgstream;libicv;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libBObol_image_source COMMAND test_bobol_image_source)

brlcad_addexec(test_bobol_image_display test_image_display.cpp "libBObol;libimgstream;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libBObol_image_display COMMAND test_bobol_image_display)

brlcad_addexec(test_bobol_lod_realization test_lod_realization.cpp "libBObol;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libBObol_lod_realization COMMAND test_bobol_lod_realization)

brlcad_addexec(test_bobol_lod_service test_lod_service.cpp "libBObol;libwdb;libbu;Obol::Obol;Threads::Threads" TEST ALL_CXX)
brlcad_add_test(NAME libBObol_lod_service COMMAND test_bobol_lod_service)
set_tests_properties(libBObol_lod_service PROPERTIES TIMEOUT 10)

brlcad_addexec(test_bobol_lod_update_action test_lod_update_action.cpp "libBObol;libwdb;libbu;Obol::Obol;Threads::Threads" TEST ALL_CXX)
brlcad_add_test(NAME libBObol_lod_update_action COMMAND test_bobol_lod_update_action)

brlcad_addexec(test_bobol_draw_cache test_draw_cache.cpp "libBObol;libwdb;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libBObol_draw_cache COMMAND test_bobol_draw_cache)

brlcad_addexec(test_bobol_evaluated_points test_evaluated_points.cpp "libBObol;libwdb;librt;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libBObol_evaluated_points COMMAND test_bobol_evaluated_points)

brlcad_addexec(test_bobol_mesh_lod_cache test_mesh_lod_cache.cpp "libBObol;libwdb;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libBObol_mesh_lod_cache COMMAND test_bobol_mesh_lod_cache)

brlcad_addexec(test_bobol_view_store test_view_store.cpp "libBObol;libwdb;librt;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libBObol_view_store COMMAND test_bobol_view_store)

brlcad_addexec(test_bobol_rt_render test_rt_render.cpp "libBObol;libwdb;librt;libbu;Obol::Obol;Threads::Threads" TEST ALL_CXX)
brlcad_add_test(NAME libBObol_rt_render COMMAND test_bobol_rt_render)

brlcad_addexec(test_bobol_window_host test_window_host.cpp "libBObol;libimgstream;libbu;Obol::Obol;Threads::Threads" TEST ALL_CXX)
brlcad_add_test(NAME libBObol_window_host COMMAND test_bobol_window_host)

brlcad_addexec(test_bobol_headless_window_host test_headless_window_host.cpp "libBObol;libimgstream;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libBObol_headless_window_host COMMAND test_bobol_headless_window_host)

if(Obol_NanoRTContextManager_FOUND)
  brlcad_addexec(test_bobol_nanort_render test_nanort_render.cpp "libBObol;libwdb;libbu;Obol::Obol;Threads::Threads;${OSMESA_LIBRARY}" TEST ALL_CXX)
  target_include_directories(test_bobol_nanort_render SYSTEM PRIVATE
    ${OSMESA_INCLUDE_DIR}
    ${BOBOL_OPENGL_INCLUDE_DIRS}
    ${Obol_NanoRTContextManager_INCLUDE_DIR}
  )
  brlcad_add_test(NAME libBObol_nanort_render COMMAND test_bobol_nanort_render)
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
