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
if(BRLCAD_ENABLE_ADDRESS_SANITIZER OR BRLCAD_ENABLE_UNDEFINED_SANITIZER OR BRLCAD_SANITIZE_THREAD)
  set_tests_properties(libBObol_lod_service PROPERTIES TIMEOUT 60)
else()
  set_tests_properties(libBObol_lod_service PROPERTIES TIMEOUT 10)
endif()

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

brlcad_addexec(test_bobol_compact_edit test_compact_edit.cpp "libBObol;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libBObol_compact_edit COMMAND test_bobol_compact_edit)

brlcad_addexec(test_bobol_rt_render test_rt_render.cpp "libBObol;libwdb;librt;libbu;Obol::Obol;Threads::Threads" TEST ALL_CXX)
brlcad_add_test(NAME libBObol_rt_render COMMAND test_bobol_rt_render)

brlcad_addexec(test_bobol_window_host test_window_host.cpp "libBObol;libimgstream;libbu;Obol::Obol;Threads::Threads" TEST ALL_CXX)
brlcad_add_test(NAME libBObol_window_host COMMAND test_bobol_window_host)

brlcad_addexec(test_bobol_headless_window_host test_headless_window_host.cpp "libBObol;libimgstream;libbu;Obol::Obol" TEST ALL_CXX)
brlcad_add_test(NAME libBObol_headless_window_host COMMAND test_bobol_headless_window_host)

# Every installed header is a supported include boundary, including the
# advanced source tier.  Compile each header in an otherwise empty translation
# unit so accidental reliance on include order is caught at build time.
set(_bobol_header_self_containment_sources)
foreach(_bobol_header IN LISTS LIBBOBOL_INSTALLED_HEADERS)
  string(REPLACE ".h" "" _bobol_header_stem "${_bobol_header}")
  set(_bobol_header_source
    "${CMAKE_CURRENT_BINARY_DIR}/header-self-containment/${_bobol_header_stem}.cpp")
  file(GENERATE
    OUTPUT "${_bobol_header_source}"
    CONTENT "#include \"BObol/${_bobol_header}\"\n")
  list(APPEND _bobol_header_self_containment_sources
    "${_bobol_header_source}")
endforeach()
add_library(test_bobol_header_self_containment OBJECT
  ${_bobol_header_self_containment_sources})
target_link_libraries(test_bobol_header_self_containment PRIVATE libBObol)
set_target_properties(test_bobol_header_self_containment PROPERTIES
  UNITY_BUILD OFF
  FOLDER "BRL-CAD Tests/libBObol")

brlcad_addexec(test_bobol_public_umbrella test_public_umbrella.cpp
  "libBObol" TEST ALL_CXX)
brlcad_add_test(NAME libBObol_public_umbrella
  COMMAND test_bobol_public_umbrella)
set_property(TEST libBObol_public_umbrella APPEND PROPERTY LABELS bobol_api)

brlcad_add_test(NAME libBObol_public_symbol_contract
  COMMAND "${CMAKE_COMMAND}"
    "-DNM=${CMAKE_NM}"
    "-DLIBRARY=$<TARGET_FILE:libBObol>"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/../check_public_symbols.cmake")
set_tests_properties(libBObol_public_symbol_contract PROPERTIES
  LABELS "bobol_api;bobol_headless")

# Exercise the actual install tree, including BRLCADConfig's Obol dependency
# resolution.  The script stages the configured installation under the build
# tree and then configures and links a fresh out-of-tree consumer.
brlcad_add_test(NAME libBObol_installed_package_consumer
  COMMAND "${CMAKE_COMMAND}"
    "-DBUILD_DIR=${CMAKE_BINARY_DIR}"
    "-DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/package_consumer"
    "-DSTAGE_DIR=${CMAKE_CURRENT_BINARY_DIR}/installed-package"
    "-DCONFIG=$<CONFIG>"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/check_installed_package.cmake")
set_tests_properties(libBObol_installed_package_consumer PROPERTIES
  LABELS "bobol_api;bobol_package"
  TIMEOUT 300)

if(Obol_NanoRTContextManager_FOUND)
  brlcad_addexec(test_bobol_nanort_render test_nanort_render.cpp "libBObol;libwdb;libbu;Obol::Obol;Threads::Threads;${OSMESA_LIBRARY}" TEST ALL_CXX)
  target_include_directories(test_bobol_nanort_render SYSTEM PRIVATE
    ${OSMESA_INCLUDE_DIR}
    ${BOBOL_OPENGL_INCLUDE_DIRS}
    ${Obol_NanoRTContextManager_INCLUDE_DIR}
  )
  brlcad_add_test(NAME libBObol_nanort_render COMMAND test_bobol_nanort_render)
endif()

# This label is a contract: every listed test must pass without an X11 or
# Wayland display and must not open a native display provider.  OSMesa and the
# explicit headless host are permitted.
set(_bobol_headless_tests
  libBObol_prototype
  libBObol_render
  libBObol_image_source
  libBObol_image_display
  libBObol_lod_realization
  libBObol_lod_service
  libBObol_lod_update_action
  libBObol_draw_cache
  libBObol_evaluated_points
  libBObol_mesh_lod_cache
  libBObol_view_store
  libBObol_compact_edit
  libBObol_rt_render
  libBObol_window_host
  libBObol_headless_window_host
  libBObol_public_umbrella
  libBObol_public_symbol_contract)
if(Obol_NanoRTContextManager_FOUND)
  list(APPEND _bobol_headless_tests libBObol_nanort_render)
endif()
foreach(_bobol_test IN LISTS _bobol_headless_tests)
  if(TEST ${_bobol_test})
    set_property(TEST ${_bobol_test} APPEND PROPERTY LABELS bobol_headless)
  endif()
endforeach()

# Focused contracts for the ASan/LSan/UBSan and TSan build jobs.  The label is
# intentionally smaller than the full headless suite so sanitizer failures map
# to ownership, stale-handle, endpoint, or worker/publication boundaries.
set(_bobol_sanitizer_tests
  libBObol_lod_service
  libBObol_view_store
  libBObol_compact_edit
  libBObol_rt_render
  libBObol_headless_window_host)
foreach(_bobol_test IN LISTS _bobol_sanitizer_tests)
  if(TEST ${_bobol_test})
    set_property(TEST ${_bobol_test} APPEND PROPERTY LABELS bobol_sanitizer)
  endif()
endforeach()

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
  test_public_umbrella.cpp
  check_installed_package.cmake
  package_consumer/CMakeLists.txt
  package_consumer/main.cpp
  test_prototype.cpp
  test_rt_render.cpp
  test_render.cpp
  test_view_store.cpp
  test_compact_edit.cpp
  test_window_host.cpp
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
