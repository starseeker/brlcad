if(BRLCAD_ENABLE_QT)
  find_package_qt(REQUIRED)
endif(BRLCAD_ENABLE_QT)

# TODO - need this now for code in lower level directories,
# need to do some cleanup.
include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${BRLCAD_SOURCE_DIR}/src/libqtcad)

# NOTE: We can't use the CMake global values for automoc et. al.
# because most of BRL-CAD's targets are not Qt targets.  We don't
# want to run the automatic moc logic for all of them.

set(
  qged_srcs
  QgEdMainWindow.cpp
  QgEdApp.cpp
  QgEdFilter.cpp
  QgGuiTestDriver.cpp
  QgProgressiveDiagnostics.cpp
  fbserv.cpp
  QgEdPalette.cpp
)

set(
  qth_names
  QgEdApp
  QgEdFilter
  fbserv
  QgEdMainWindow
  QgEdPalette
)
foreach(qh ${qth_names})
  set(moc_headers ${moc_headers} ${qh}.h)
endforeach(qh ${qth_names})

if(BRLCAD_ENABLE_QT)
  if(Qt6Widgets_FOUND)
    qt6_wrap_cpp(moc_srcs ${moc_headers})
  else()
    qt5_wrap_cpp(moc_srcs ${moc_headers})
  endif(Qt6Widgets_FOUND)
  set(fbserv_moc_src "${CMAKE_CURRENT_BINARY_DIR}/moc_fbserv.cpp")

  # https://github.com/Alexhuszagh/BreezeStyleSheets
  if(Qt6Widgets_FOUND)
    qt6_add_resources(theme_qrc theme/breeze.qrc)
  else()
    qt5_add_resources(theme_qrc theme/breeze.qrc)
  endif(Qt6Widgets_FOUND)

  add_executable(qged WIN32 main.cpp ${qged_srcs} ${moc_srcs} ${theme_qrc})
  # Note - Network is needed for async embedded raytracing display
  if(Qt6Widgets_FOUND)
    target_link_libraries(
      qged
      libqtcad
      libged
      libimgstream
      libpkg
      librt
      libbu
      Qt6::Widgets
      Qt6::Network
      Qt6::OpenGL
      Qt6::OpenGLWidgets
    )
  else()
    target_link_libraries(
      qged
      libqtcad
      libged
      libimgstream
      libpkg
      librt
      libbu
      Qt5::Widgets
      Qt5::Network
    )
  endif(Qt6Widgets_FOUND)
  target_compile_definitions(qged PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  install(TARGETS qged RUNTIME DESTINATION ${BIN_DIR} LIBRARY DESTINATION ${LIB_DIR} ARCHIVE DESTINATION ${LIB_DIR})

  brlcad_add_test(NAME qged_event_script_replay
    COMMAND qged --swrast --test-script
    "${CMAKE_CURRENT_SOURCE_DIR}/test_qged_event_script.json")
  set_property(TEST qged_event_script_replay APPEND PROPERTY ENVIRONMENT
    BRLCAD_ROOT=${CMAKE_BINARY_DIR} QT_QPA_PLATFORM=offscreen)
  set_property(TEST qged_event_script_replay APPEND PROPERTY LABELS drawing_baseline)

  brlcad_add_test(NAME qged_polygon_ui_replay
    COMMAND qged --swrast --test-script
    "${CMAKE_CURRENT_SOURCE_DIR}/test_qged_polygon_ui.json")
  set_property(TEST qged_polygon_ui_replay APPEND PROPERTY ENVIRONMENT
    BRLCAD_ROOT=${CMAKE_BINARY_DIR} QT_QPA_PLATFORM=offscreen)
  set_property(TEST qged_polygon_ui_replay APPEND PROPERTY LABELS drawing_baseline)

  brlcad_add_test(NAME qged_polygon_ui_quad_replay
    COMMAND qged --swrast --quad --test-script
    "${CMAKE_CURRENT_SOURCE_DIR}/test_qged_polygon_ui.json")
  set_property(TEST qged_polygon_ui_quad_replay APPEND PROPERTY ENVIRONMENT
    BRLCAD_ROOT=${CMAKE_BINARY_DIR} QT_QPA_PLATFORM=offscreen)
  set_property(TEST qged_polygon_ui_quad_replay APPEND PROPERTY LABELS drawing_baseline)

  brlcad_add_test(NAME qged_polygon_sketch_sync_replay
    COMMAND ${CMAKE_COMMAND}
    -DBRLCAD_ROOT=${CMAKE_BINARY_DIR}
    -DMGED=$<TARGET_FILE:mged>
    -DQGED=$<TARGET_FILE:qged>
    -DSCRIPT=${CMAKE_CURRENT_SOURCE_DIR}/test_qged_polygon_sketch_sync.json
    -P ${CMAKE_CURRENT_SOURCE_DIR}/test_qged_polygon_sketch_sync.cmake)
  set_property(TEST qged_polygon_sketch_sync_replay APPEND PROPERTY LABELS
    drawing_baseline)

  brlcad_add_test(NAME qged_polygon_sketch_sync_quad_replay
    COMMAND ${CMAKE_COMMAND}
    -DBRLCAD_ROOT=${CMAKE_BINARY_DIR}
    -DMGED=$<TARGET_FILE:mged>
    -DQGED=$<TARGET_FILE:qged>
    -DSCRIPT=${CMAKE_CURRENT_SOURCE_DIR}/test_qged_polygon_sketch_sync.json
    -DQUAD=ON
    -P ${CMAKE_CURRENT_SOURCE_DIR}/test_qged_polygon_sketch_sync.cmake)
  set_property(TEST qged_polygon_sketch_sync_quad_replay APPEND PROPERTY LABELS
    drawing_baseline)

  brlcad_add_test(NAME qged_primitive_edit_ui_replay
    COMMAND ${CMAKE_COMMAND}
    -DBRLCAD_ROOT=${CMAKE_BINARY_DIR}
    -DMGED=$<TARGET_FILE:mged>
    -DQGED=$<TARGET_FILE:qged>
    -DSCRIPT=${CMAKE_CURRENT_SOURCE_DIR}/test_qged_primitive_edit_ui.json
    -P ${CMAKE_CURRENT_SOURCE_DIR}/test_qged_primitive_edit_ui.cmake)
  set_property(TEST qged_primitive_edit_ui_replay APPEND PROPERTY LABELS
    drawing_baseline)

  brlcad_add_test(NAME qged_primitive_edit_ui_quad_replay
    COMMAND ${CMAKE_COMMAND}
    -DBRLCAD_ROOT=${CMAKE_BINARY_DIR}
    -DMGED=$<TARGET_FILE:mged>
    -DQGED=$<TARGET_FILE:qged>
    -DSCRIPT=${CMAKE_CURRENT_SOURCE_DIR}/test_qged_primitive_edit_ui.json
    -DQUAD=ON
    -P ${CMAKE_CURRENT_SOURCE_DIR}/test_qged_primitive_edit_ui.cmake)
  set_property(TEST qged_primitive_edit_ui_quad_replay APPEND PROPERTY LABELS
    drawing_baseline)

  # qged's view and embedded framebuffer paths are Obol/libBObol based.
  # GED plugins are still loaded at runtime; legacy display-manager plugins are
  # not a qged build prerequisite.  qged's own Qt plugins are added below and
  # wired in after src/qged/plugins defines the aggregate target.
  add_dependencies(qged ged_plugins)

  add_executable(test_qged_obol_fbserv test_qged_obol_fbserv.cpp fbserv.cpp ${fbserv_moc_src})
  if(Qt6Widgets_FOUND)
    target_link_libraries(
      test_qged_obol_fbserv
      libqtcad
      libged
      libimgstream
      libpkg
      librt
      libbu
      Qt6::Widgets
      Qt6::Network
      Qt6::OpenGL
      Qt6::OpenGLWidgets
    )
  else()
    target_link_libraries(
      test_qged_obol_fbserv
      libqtcad
      libged
      libimgstream
      libpkg
      librt
      libbu
      Qt5::Widgets
      Qt5::Network
    )
  endif(Qt6Widgets_FOUND)
  target_compile_definitions(test_qged_obol_fbserv PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  brlcad_add_test(NAME test_qged_obol_fbserv COMMAND test_qged_obol_fbserv)
  set_property(TEST test_qged_obol_fbserv
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qged_obol_fbserv
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qged_obol_fbserv
    APPEND PROPERTY LABELS drawing_baseline)

  foreach(qh ${qth_names})
    distclean(moc_${qh}.cpp_parameters)
  endforeach(qh ${qth_names})
  distclean(${CMAKE_CURRENT_BINARY_DIR}/theme/breeze.qrc.depends)

  # On some platforms distclean has additional cleanup to do
  distclean("${CMAKE_BINARY_DIR}/${BIN_DIR}/qged.pdb")

endif(BRLCAD_ENABLE_QT)

# The GUI matrix uses jq for its offline, observational progressive-control
# refinement checker.  Keep jq optional for ordinary qged builds, but exercise
# the checker itself whenever the qualification dependency is available.
find_program(QGED_JQ_EXECUTABLE NAMES jq)
if(QGED_JQ_EXECUTABLE)
  brlcad_add_test(
    NAME qged_lod_control_trace_contract
    COMMAND
      ${SH_EXEC}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_lod_control_trace.sh
      ${QGED_JQ_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/lod_control_trace.jq
  )
  set_property(TEST qged_lod_control_trace_contract APPEND PROPERTY LABELS bobol_contract)
endif(QGED_JQ_EXECUTABLE)

# Core application is defined, now handle plugins
add_subdirectory(plugins)

if(BRLCAD_ENABLE_QT AND TARGET qged AND TARGET qged_plugins)
  add_dependencies(qged qged_plugins)
endif()

brlcad_addexec(
  qged_unique_mesh_stress_generator
  tests/generate_unique_mesh_stress.cpp
  "libwdb;librt;libbu"
  NO_INSTALL
  ALL_CXX
)

brlcad_addexec(
  qged_apng_encode
  tests/qged_apng_encode.cpp
  "libicv;libbu"
  NO_INSTALL
  ALL_CXX
)

# Theme subdirectory handles file bookkeeping for CMake
add_subdirectory(theme)

set(
  qged_ignore_srcs
  ${qged_srcs}
  ${moc_headers}
  BUILDING_QT
  CMakeLists.txt
  QgEdCategories.h
  QgGuiTestDriver.h
  TODO
  fonts/Inconsolata-LICENSE.txt
  fonts/Inconsolata.otf
  main.cpp
  test_qged_obol_fbserv.cpp
  test_qged_event_script.json
  test_qged_polygon_ui.json
  test_qged_polygon_sketch_sync.cmake
  test_qged_polygon_sketch_sync.json
  test_qged_primitive_edit_ui.cmake
  test_qged_primitive_edit_ui.json
  tests/qged_framebuffer_matrix.sh
  tests/qged_polygon_visual_matrix.sh
  tests/README.md
  tests/generate_many_lucy_stress.tcl
  tests/generate_unique_mesh_stress.cpp
  tests/lod_control_trace.jq
  tests/qged_apng_encode.cpp
  tests/qged_gui_matrix.sh
  tests/qged_resize_matrix.sh
  tests/qged_selection_ui_matrix.sh
  tests/test_lod_control_trace.sh
  qged_signals.eps
  qged_signals.dia
  view_update_cycle.odg
  view_update_cycle.pdf
)

cmakefiles(${qged_ignore_srcs})

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
