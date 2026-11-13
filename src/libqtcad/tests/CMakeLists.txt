include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)

set(TEST_QTCAD_SESSION_SRCS test_qtcad_session.cpp)
set(TEST_QTCAD_MODEL_INVARIANTS_SRCS test_qtcad_model_invariants.cpp)
set(TEST_QTCAD_PRIMITIVE_EDIT_SRCS test_qtcad_primitive_edit.cpp)
set(TEST_QTCAD_CONSOLE_SRCS test_qtcad_console.cpp)
set(TEST_QTCAD_WIDGETS_SRCS test_qtcad_widgets.cpp)
set(TEST_QTCAD_PLUGIN_MANAGER_SRCS test_qtcad_plugin_manager.cpp)
set(TEST_QTCAD_PLUGIN_CONTEXT_SRCS test_qtcad_plugin_context.cpp)
set(TEST_QTCAD_PLUGIN_BASE_HDRS test_qtcad_plugin_base.h)
set(TEST_QTCAD_PLUGIN_VALID_SRCS test_qtcad_plugin_valid.cpp)
set(TEST_QTCAD_PLUGIN_MISSING_API_SRCS test_qtcad_plugin_missing_api.cpp)
set(TEST_QTCAD_PLUGIN_PANEL_SRCS test_qtcad_plugin_panel.cpp)
set(TEST_QTCAD_PLUGIN_DIALOG_SRCS test_qtcad_plugin_dialog.cpp)
set(TEST_QTCAD_EVENT_PLAYER_SRCS test_qtcad_event_player.cpp)
set(TEST_QTCAD_PLUGIN_METADATA
  test_qtcad_plugin_valid.json
  test_qtcad_plugin_missing_api.json
  test_qtcad_plugin_panel.json
  test_qtcad_plugin_dialog.json
)
set(TEST_QTCAD_PLUGIN_REUSE_SRCS test_qtcad_plugin_reuse.cpp)
set(TEST_QTCAD_OBOL_CONTROLLER_SRCS test_qtcad_obol_controller.cpp)
set(TEST_QTCAD_OBOL_DRAW_SYNC_SRCS test_qtcad_obol_draw_sync.cpp ../QgCanvasInput.cpp)
set(TEST_QTCAD_OBOL_EDIT_PREVIEW_SRCS test_qtcad_obol_edit_preview.cpp)
set(TEST_QTCAD_OBOL_EXPORT_SRCS test_qtcad_obol_export.cpp)
set(TEST_QTCAD_OBOL_FACEPLATE_SYNC_SRCS test_qtcad_obol_faceplate_sync.cpp)
set(TEST_QTCAD_OBOL_MEASURE_SRCS test_qtcad_obol_measure.cpp)
set(TEST_QTCAD_OBOL_PICK_SRCS test_qtcad_obol_pick.cpp)
set(TEST_QTCAD_OBOL_PRESENTATION_HDRS qtcad_obol_test_presentation.h)
set(TEST_QTCAD_OBOL_PROGRESSIVE_LOD_SRCS
  test_qtcad_obol_progressive_lod.cpp
  ${TEST_QTCAD_OBOL_PRESENTATION_HDRS})
set(TEST_QTCAD_OBOL_REAL_MODELS_SRCS
  test_qtcad_obol_real_models.cpp
  ${TEST_QTCAD_OBOL_PRESENTATION_HDRS})
set(TEST_QTCAD_OBOL_SELECTION_SYNC_SRCS test_qtcad_obol_selection_sync.cpp)
set(TEST_QTCAD_OBOL_SNAP_SRCS test_qtcad_obol_snap.cpp)
set(TEST_QTCAD_OBOL_WINDOW_HOST_SRCS test_qtcad_obol_window_host.cpp)
set(TEST_QGMODEL_MODEL_TESTER_SRCS test_qgmodel_model_tester.cpp)
set(TEST_QGMODEL_SCALE_SRCS test_qgmodel_scale.cpp)

set(QTCAD_TEST_PLUGIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/plugins)
if(TARGET Qt6::Test)
  set(QTCAD_TEST_QT_CORE Qt6::Core)
  set(QTCAD_TEST_QT_TEST Qt6::Test)
  set(QTCAD_TEST_QT_WIDGETS Qt6::Widgets)
elseif(TARGET Qt5::Test)
  set(QTCAD_TEST_QT_CORE Qt5::Core)
  set(QTCAD_TEST_QT_TEST Qt5::Test)
  set(QTCAD_TEST_QT_WIDGETS Qt5::Widgets)
endif()

if(QTCAD_TEST_QT_TEST)
  add_executable(test_qtcad_model_invariants ${TEST_QTCAD_MODEL_INVARIANTS_SRCS})
  set_target_properties(test_qtcad_model_invariants PROPERTIES AUTOMOC ON)
  target_compile_definitions(test_qtcad_model_invariants PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  target_link_libraries(test_qtcad_model_invariants
    libqtcad
    libged
    libwdb
    librt
    libbu
    ${QTCAD_TEST_QT_TEST}
    ${QTCAD_TEST_QT_WIDGETS}
  )
  add_dependencies(test_qtcad_model_invariants ged_plugins)
  brlcad_add_test(NAME test_qtcad_model_invariants COMMAND test_qtcad_model_invariants)
  set_property(TEST test_qtcad_model_invariants APPEND PROPERTY ENVIRONMENT
    BRLCAD_ROOT=${CMAKE_BINARY_DIR}
    QT_QPA_PLATFORM=offscreen
  )

  add_executable(test_qtcad_plugin_context ${TEST_QTCAD_PLUGIN_CONTEXT_SRCS})
  set_target_properties(test_qtcad_plugin_context PROPERTIES AUTOMOC ON)
  target_compile_definitions(test_qtcad_plugin_context PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  target_link_libraries(test_qtcad_plugin_context
    libqtcad
    ${QTCAD_TEST_QT_CORE}
    ${QTCAD_TEST_QT_TEST}
  )
  brlcad_add_test(NAME test_qtcad_plugin_context COMMAND test_qtcad_plugin_context)

  add_library(qtcad_test_plugin_valid MODULE
    ${TEST_QTCAD_PLUGIN_VALID_SRCS}
    ${TEST_QTCAD_PLUGIN_BASE_HDRS}
  )
  add_library(qtcad_test_plugin_missing_api MODULE
    ${TEST_QTCAD_PLUGIN_MISSING_API_SRCS}
    ${TEST_QTCAD_PLUGIN_BASE_HDRS}
  )
  add_library(qtcad_test_plugin_panel MODULE
    ${TEST_QTCAD_PLUGIN_PANEL_SRCS}
  )
  add_library(qtcad_test_plugin_dialog MODULE
    ${TEST_QTCAD_PLUGIN_DIALOG_SRCS}
  )
  foreach(plugin_target
      qtcad_test_plugin_valid
      qtcad_test_plugin_missing_api
      qtcad_test_plugin_panel
      qtcad_test_plugin_dialog)
    set_target_properties(${plugin_target} PROPERTIES
      AUTOMOC ON
      PREFIX ""
      LIBRARY_OUTPUT_DIRECTORY ${QTCAD_TEST_PLUGIN_DIR}
      RUNTIME_OUTPUT_DIRECTORY ${QTCAD_TEST_PLUGIN_DIR}
    )
    target_compile_definitions(${plugin_target} PRIVATE BRLCADBUILD HAVE_CONFIG_H)
    target_link_libraries(${plugin_target}
      libqtcad
      ${QTCAD_TEST_QT_CORE}
      ${QTCAD_TEST_QT_WIDGETS}
    )
  endforeach()

  add_executable(test_qtcad_plugin_manager ${TEST_QTCAD_PLUGIN_MANAGER_SRCS})
  set_target_properties(test_qtcad_plugin_manager PROPERTIES AUTOMOC ON)
  target_compile_definitions(test_qtcad_plugin_manager PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  target_link_libraries(test_qtcad_plugin_manager
    libqtcad
    ${QTCAD_TEST_QT_CORE}
    ${QTCAD_TEST_QT_TEST}
    ${QTCAD_TEST_QT_WIDGETS}
  )
  add_dependencies(test_qtcad_plugin_manager
    qtcad_test_plugin_valid
    qtcad_test_plugin_missing_api
    qtcad_test_plugin_panel
    qtcad_test_plugin_dialog
  )
  brlcad_add_test(NAME test_qtcad_plugin_manager COMMAND test_qtcad_plugin_manager)
  set_property(TEST test_qtcad_plugin_manager APPEND PROPERTY ENVIRONMENT
    QTCAD_TEST_PLUGIN_DIR=${QTCAD_TEST_PLUGIN_DIR}
    BRLCAD_ROOT=${CMAKE_BINARY_DIR}
    QT_QPA_PLATFORM=offscreen
  )
endif()

if(TARGET Qt6::Test)
  add_executable(test_qtcad_console ${TEST_QTCAD_CONSOLE_SRCS})
  target_link_libraries(test_qtcad_console libqtcad Qt6::Test Qt6::Widgets)
  add_executable(test_qtcad_widgets ${TEST_QTCAD_WIDGETS_SRCS})
  target_link_libraries(test_qtcad_widgets libqtcad libwdb Qt6::Test Qt6::Widgets)
  add_executable(test_qtcad_session ${TEST_QTCAD_SESSION_SRCS})
  target_link_libraries(test_qtcad_session libqtcad Qt6::Test)
elseif(TARGET Qt5::Test)
  add_executable(test_qtcad_console ${TEST_QTCAD_CONSOLE_SRCS})
  target_link_libraries(test_qtcad_console libqtcad Qt5::Test Qt5::Widgets)
  add_executable(test_qtcad_widgets ${TEST_QTCAD_WIDGETS_SRCS})
  target_link_libraries(test_qtcad_widgets libqtcad libwdb Qt5::Test Qt5::Widgets)
  add_executable(test_qtcad_session ${TEST_QTCAD_SESSION_SRCS})
  target_link_libraries(test_qtcad_session libqtcad Qt5::Test)
endif()
foreach(qtcad_test_target test_qtcad_console test_qtcad_widgets test_qtcad_session)
  if(TARGET ${qtcad_test_target})
    set_target_properties(${qtcad_test_target} PROPERTIES AUTOMOC ON)
    target_compile_definitions(${qtcad_test_target} PRIVATE BRLCADBUILD HAVE_CONFIG_H)
    add_dependencies(${qtcad_test_target} ged_plugins)
    brlcad_add_test(NAME ${qtcad_test_target} COMMAND ${qtcad_test_target})
  endif()
endforeach()
if(TARGET test_qtcad_widgets)
  set_property(TEST test_qtcad_widgets APPEND PROPERTY ENVIRONMENT
    QT_QPA_PLATFORM=offscreen)
endif()
if(TARGET test_qtcad_session)
  set_property(TEST test_qtcad_session APPEND PROPERTY ENVIRONMENT
    BRLCAD_ROOT=${CMAKE_BINARY_DIR} QT_QPA_PLATFORM=offscreen)
endif()

if(TARGET Qt6::Widgets)
  add_executable(test_qtcad_primitive_edit ${TEST_QTCAD_PRIMITIVE_EDIT_SRCS})
  target_link_libraries(test_qtcad_primitive_edit libqtcad libged libwdb librt libbu Qt6::Widgets)
elseif(TARGET Qt5::Widgets)
  add_executable(test_qtcad_primitive_edit ${TEST_QTCAD_PRIMITIVE_EDIT_SRCS})
  target_link_libraries(test_qtcad_primitive_edit libqtcad libged libwdb librt libbu Qt5::Widgets)
endif()
if(TARGET test_qtcad_primitive_edit)
  target_compile_definitions(test_qtcad_primitive_edit PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  add_dependencies(test_qtcad_primitive_edit ged_plugins)
  brlcad_add_test(NAME test_qtcad_primitive_edit COMMAND test_qtcad_primitive_edit)
  set_property(TEST test_qtcad_primitive_edit APPEND PROPERTY ENVIRONMENT
    BRLCAD_ROOT=${CMAKE_BINARY_DIR} QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_primitive_edit APPEND PROPERTY LABELS drawing_baseline)
endif()

if(TARGET Qt6::Widgets)
  add_executable(test_qtcad_event_player ${TEST_QTCAD_EVENT_PLAYER_SRCS})
  target_link_libraries(test_qtcad_event_player libqtcad Qt6::Widgets)
elseif(TARGET Qt5::Widgets)
  add_executable(test_qtcad_event_player ${TEST_QTCAD_EVENT_PLAYER_SRCS})
  target_link_libraries(test_qtcad_event_player libqtcad Qt5::Widgets)
endif()
if(TARGET test_qtcad_event_player)
  target_compile_definitions(test_qtcad_event_player PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  brlcad_add_test(NAME test_qtcad_event_player COMMAND test_qtcad_event_player)
  set_property(TEST test_qtcad_event_player APPEND PROPERTY ENVIRONMENT
    BRLCAD_ROOT=${CMAKE_BINARY_DIR}
    QT_QPA_PLATFORM=offscreen
  )
endif()

if(TARGET Qt6::Test)
  add_executable(test_qgmodel_model_tester ${TEST_QGMODEL_MODEL_TESTER_SRCS})
  target_link_libraries(
    test_qgmodel_model_tester
    libqtcad
    libged
    librt
    libbu
    Qt6::Widgets
    Qt6::Test
    )

  add_executable(test_qgmodel_scale ${TEST_QGMODEL_SCALE_SRCS})
  target_link_libraries(
    test_qgmodel_scale
    libqtcad
    libged
    librt
    libbu
    Qt6::Widgets
    Qt6::Test
    )
endif(TARGET Qt6::Test)
if(TARGET Qt5::Test)
  add_executable(test_qgmodel_model_tester ${TEST_QGMODEL_MODEL_TESTER_SRCS})
  target_link_libraries(
    test_qgmodel_model_tester
    libqtcad
    libged
    librt
    libbu
    Qt5::Widgets
    Qt5::Test
    )

  add_executable(test_qgmodel_scale ${TEST_QGMODEL_SCALE_SRCS})
  target_link_libraries(
    test_qgmodel_scale
    libqtcad
    libged
    librt
    libbu
    Qt5::Widgets
    Qt5::Test
    )
endif(TARGET Qt5::Test)
if(TARGET test_qgmodel_model_tester)
  target_compile_definitions(test_qgmodel_model_tester PRIVATE BRLCADBUILD HAVE_CONFIG_H)
endif(TARGET test_qgmodel_model_tester)
if(TARGET test_qgmodel_scale)
  target_compile_definitions(test_qgmodel_scale PRIVATE BRLCADBUILD HAVE_CONFIG_H)
endif(TARGET test_qgmodel_scale)

if(TARGET Qt6::Widgets)
  add_executable(test_qtcad_obol_controller ${TEST_QTCAD_OBOL_CONTROLLER_SRCS})
  target_link_libraries(
    test_qtcad_obol_controller
    libqtcad
    libbu
    Qt6::Widgets
    )
endif(TARGET Qt6::Widgets)
if(NOT TARGET test_qtcad_obol_controller AND TARGET Qt5::Widgets)
  add_executable(test_qtcad_obol_controller ${TEST_QTCAD_OBOL_CONTROLLER_SRCS})
  target_link_libraries(
    test_qtcad_obol_controller
    libqtcad
    libbu
    Qt5::Widgets
    )
endif(NOT TARGET test_qtcad_obol_controller AND TARGET Qt5::Widgets)
if(TARGET test_qtcad_obol_controller)
  target_compile_definitions(test_qtcad_obol_controller PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  brlcad_add_test(NAME test_qtcad_obol_controller COMMAND test_qtcad_obol_controller)
  set_property(TEST test_qtcad_obol_controller
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_controller
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  # Exercise retained software-frame presentation at a fractional device
  # scale.  A DPR-1 retained image can have the right physical dimensions yet
  # make the whole viewport (including the faceplate) grow on an idle repaint.
  set_property(TEST test_qtcad_obol_controller
    APPEND PROPERTY ENVIRONMENT QT_SCALE_FACTOR=1.5)
  set_property(TEST test_qtcad_obol_controller
    APPEND PROPERTY LABELS drawing_baseline)
endif(TARGET test_qtcad_obol_controller)

if(TARGET Qt6::Widgets)
  add_executable(test_qtcad_obol_window_host ${TEST_QTCAD_OBOL_WINDOW_HOST_SRCS})
  target_link_libraries(
    test_qtcad_obol_window_host
    libqtcad
    libimgstream
    libged
    libwdb
    libbu
    Qt6::Widgets
    )
endif(TARGET Qt6::Widgets)
if(NOT TARGET test_qtcad_obol_window_host AND TARGET Qt5::Widgets)
  add_executable(test_qtcad_obol_window_host ${TEST_QTCAD_OBOL_WINDOW_HOST_SRCS})
  target_link_libraries(
    test_qtcad_obol_window_host
    libqtcad
    libimgstream
    libged
    libwdb
    libbu
    Qt5::Widgets
    )
endif(NOT TARGET test_qtcad_obol_window_host AND TARGET Qt5::Widgets)
if(TARGET test_qtcad_obol_window_host)
  add_dependencies(test_qtcad_obol_window_host ged_plugins)
  target_compile_definitions(test_qtcad_obol_window_host PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  brlcad_add_test(NAME test_qtcad_obol_window_host COMMAND test_qtcad_obol_window_host)
  set_property(TEST test_qtcad_obol_window_host
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_window_host
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_window_host
    APPEND PROPERTY LABELS
    "drawing_baseline;bobol_graphics_qualification;bobol_graphics_offscreen")
  if(BRLCAD_ENABLE_OPENGL AND UNIX)
    find_program(XVFB_RUN_EXECUTABLE xvfb-run)
    if(XVFB_RUN_EXECUTABLE)
      if(TARGET test_qtcad_obol_controller)
        add_test(
          NAME test_qtcad_obol_controller_x11_gl
          COMMAND ${SH_EXEC} ${BRLCAD_SOURCE_DIR}/misc/CMake/scripts/xvfb_test.sh
                  ${XVFB_RUN_EXECUTABLE} ${CMAKE_COMMAND} -E env
                  QT_QPA_PLATFORM=xcb QT_SCALE_FACTOR=1.5
                  $<TARGET_FILE:test_qtcad_obol_controller>)
        set_tests_properties(test_qtcad_obol_controller_x11_gl PROPERTIES
          ENVIRONMENT "BRLCAD_ROOT=${CMAKE_BINARY_DIR}"
          SKIP_RETURN_CODE 125 TIMEOUT 120
          LABELS "drawing_baseline;bobol_graphics_qualification;bobol_graphics_system_gl")
      endif()
      add_test(
        NAME test_qtcad_obol_window_host_x11_gl
        COMMAND ${SH_EXEC} ${BRLCAD_SOURCE_DIR}/misc/CMake/scripts/xvfb_test.sh
                ${XVFB_RUN_EXECUTABLE} ${CMAKE_COMMAND} -E env QT_QPA_PLATFORM=xcb
                $<TARGET_FILE:test_qtcad_obol_window_host> --system-gl)
      set_property(TEST test_qtcad_obol_window_host_x11_gl
        APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
      set_property(TEST test_qtcad_obol_window_host_x11_gl
        PROPERTY SKIP_RETURN_CODE 125)
      set_property(TEST test_qtcad_obol_window_host_x11_gl
        APPEND PROPERTY LABELS
        "drawing_baseline;bobol_graphics_qualification;bobol_graphics_system_gl")
    endif(XVFB_RUN_EXECUTABLE)
  endif(BRLCAD_ENABLE_OPENGL AND UNIX)
endif(TARGET test_qtcad_obol_window_host)

if(TARGET Qt6::Widgets)
  add_executable(test_qtcad_obol_draw_sync ${TEST_QTCAD_OBOL_DRAW_SYNC_SRCS})
  target_link_libraries(
    test_qtcad_obol_draw_sync
    libqtcad
    libged
    libwdb
    librt
    libbu
    Qt6::Widgets
    )
endif(TARGET Qt6::Widgets)
if(NOT TARGET test_qtcad_obol_draw_sync AND TARGET Qt5::Widgets)
  add_executable(test_qtcad_obol_draw_sync ${TEST_QTCAD_OBOL_DRAW_SYNC_SRCS})
  target_link_libraries(
    test_qtcad_obol_draw_sync
    libqtcad
    libged
    libwdb
    librt
    libbu
    Qt5::Widgets
    )
endif(NOT TARGET test_qtcad_obol_draw_sync AND TARGET Qt5::Widgets)
if(TARGET test_qtcad_obol_draw_sync)
  target_compile_definitions(test_qtcad_obol_draw_sync PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  brlcad_add_test(NAME test_qtcad_obol_draw_sync COMMAND test_qtcad_obol_draw_sync)
  brlcad_add_test(
    NAME test_qtcad_obol_production_flow_lifecycle
    COMMAND test_qtcad_obol_draw_sync production-flow-lifecycle)
  set_tests_properties(
    test_qtcad_obol_draw_sync
    test_qtcad_obol_production_flow_lifecycle
    PROPERTIES
      ENVIRONMENT "BRLCAD_ROOT=${CMAKE_BINARY_DIR};QT_QPA_PLATFORM=offscreen"
      RESOURCE_LOCK qtcad_obol_draw_sync_fixture)
  set_property(TEST test_qtcad_obol_draw_sync
    APPEND PROPERTY LABELS drawing_baseline)
  set_tests_properties(test_qtcad_obol_production_flow_lifecycle PROPERTIES
    LABELS "bobol_flow;bobol_graphical"
    TIMEOUT 30)
  distclean(${CMAKE_CURRENT_BINARY_DIR}/qtcad_obol_draw_sync_tmp.g)
endif(TARGET test_qtcad_obol_draw_sync)

if(TARGET Qt6::Widgets)
  add_executable(test_qtcad_obol_edit_preview ${TEST_QTCAD_OBOL_EDIT_PREVIEW_SRCS})
  target_link_libraries(
    test_qtcad_obol_edit_preview
    libqtcad
    libbu
    Qt6::Widgets
    )
endif(TARGET Qt6::Widgets)
if(NOT TARGET test_qtcad_obol_edit_preview AND TARGET Qt5::Widgets)
  add_executable(test_qtcad_obol_edit_preview ${TEST_QTCAD_OBOL_EDIT_PREVIEW_SRCS})
  target_link_libraries(
    test_qtcad_obol_edit_preview
    libqtcad
    libbu
    Qt5::Widgets
    )
endif(NOT TARGET test_qtcad_obol_edit_preview AND TARGET Qt5::Widgets)
if(TARGET test_qtcad_obol_edit_preview)
  target_compile_definitions(test_qtcad_obol_edit_preview PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  brlcad_add_test(NAME test_qtcad_obol_edit_preview COMMAND test_qtcad_obol_edit_preview)
  set_property(TEST test_qtcad_obol_edit_preview
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_edit_preview
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_edit_preview
    APPEND PROPERTY LABELS drawing_baseline)
endif(TARGET test_qtcad_obol_edit_preview)

if(TARGET Qt6::Widgets)
  add_executable(test_qtcad_obol_faceplate_sync ${TEST_QTCAD_OBOL_FACEPLATE_SYNC_SRCS})
  target_link_libraries(
    test_qtcad_obol_faceplate_sync
    libqtcad
    libbu
    Qt6::Widgets
    )
endif(TARGET Qt6::Widgets)
if(NOT TARGET test_qtcad_obol_faceplate_sync AND TARGET Qt5::Widgets)
  add_executable(test_qtcad_obol_faceplate_sync ${TEST_QTCAD_OBOL_FACEPLATE_SYNC_SRCS})
  target_link_libraries(
    test_qtcad_obol_faceplate_sync
    libqtcad
    libbu
    Qt5::Widgets
    )
endif(NOT TARGET test_qtcad_obol_faceplate_sync AND TARGET Qt5::Widgets)
if(TARGET test_qtcad_obol_faceplate_sync)
  target_compile_definitions(test_qtcad_obol_faceplate_sync PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  brlcad_add_test(NAME test_qtcad_obol_faceplate_sync COMMAND test_qtcad_obol_faceplate_sync)
  set_property(TEST test_qtcad_obol_faceplate_sync
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_faceplate_sync
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_faceplate_sync
    APPEND PROPERTY LABELS drawing_baseline)
endif(TARGET test_qtcad_obol_faceplate_sync)

if(TARGET Qt6::Widgets)
  add_executable(test_qtcad_obol_selection_sync ${TEST_QTCAD_OBOL_SELECTION_SYNC_SRCS})
  target_link_libraries(
    test_qtcad_obol_selection_sync
    libqtcad
    libged
    libwdb
    librt
    libbu
    Qt6::Widgets
    )
endif(TARGET Qt6::Widgets)
if(NOT TARGET test_qtcad_obol_selection_sync AND TARGET Qt5::Widgets)
  add_executable(test_qtcad_obol_selection_sync ${TEST_QTCAD_OBOL_SELECTION_SYNC_SRCS})
  target_link_libraries(
    test_qtcad_obol_selection_sync
    libqtcad
    libged
    libwdb
    librt
    libbu
    Qt5::Widgets
    )
endif(NOT TARGET test_qtcad_obol_selection_sync AND TARGET Qt5::Widgets)
if(TARGET test_qtcad_obol_selection_sync)
  target_compile_definitions(test_qtcad_obol_selection_sync PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  brlcad_add_test(NAME test_qtcad_obol_selection_sync COMMAND test_qtcad_obol_selection_sync)
  set_property(TEST test_qtcad_obol_selection_sync
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_selection_sync
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_selection_sync
    APPEND PROPERTY LABELS drawing_baseline)
  distclean(${CMAKE_CURRENT_BINARY_DIR}/qtcad_obol_selection_sync_tmp.g)
endif(TARGET test_qtcad_obol_selection_sync)

if(TARGET Qt6::Widgets)
  add_executable(test_qtcad_obol_pick ${TEST_QTCAD_OBOL_PICK_SRCS})
  target_link_libraries(
    test_qtcad_obol_pick
    libqtcad
    libged
    libwdb
    librt
    libbu
    Qt6::Widgets
    )
endif(TARGET Qt6::Widgets)
if(NOT TARGET test_qtcad_obol_pick AND TARGET Qt5::Widgets)
  add_executable(test_qtcad_obol_pick ${TEST_QTCAD_OBOL_PICK_SRCS})
  target_link_libraries(
    test_qtcad_obol_pick
    libqtcad
    libged
    libwdb
    librt
    libbu
    Qt5::Widgets
    )
endif(NOT TARGET test_qtcad_obol_pick AND TARGET Qt5::Widgets)
if(TARGET test_qtcad_obol_pick)
  target_compile_definitions(test_qtcad_obol_pick PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  brlcad_add_test(NAME test_qtcad_obol_pick COMMAND test_qtcad_obol_pick)
  set_property(TEST test_qtcad_obol_pick
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_pick
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_pick
    APPEND PROPERTY LABELS drawing_baseline)
  distclean(${CMAKE_CURRENT_BINARY_DIR}/qtcad_obol_pick_tmp.g)
endif(TARGET test_qtcad_obol_pick)

if(TARGET Qt6::Widgets)
  add_executable(test_qtcad_obol_export ${TEST_QTCAD_OBOL_EXPORT_SRCS})
  target_link_libraries(
    test_qtcad_obol_export
    libqtcad
    libged
    libwdb
    librt
    libbu
    Qt6::Widgets
    )
endif(TARGET Qt6::Widgets)
if(NOT TARGET test_qtcad_obol_export AND TARGET Qt5::Widgets)
  add_executable(test_qtcad_obol_export ${TEST_QTCAD_OBOL_EXPORT_SRCS})
  target_link_libraries(
    test_qtcad_obol_export
    libqtcad
    libged
    libwdb
    librt
    libbu
    Qt5::Widgets
    )
endif(NOT TARGET test_qtcad_obol_export AND TARGET Qt5::Widgets)
if(TARGET test_qtcad_obol_export)
  target_compile_definitions(test_qtcad_obol_export PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  brlcad_add_test(NAME test_qtcad_obol_export COMMAND test_qtcad_obol_export)
  set_property(TEST test_qtcad_obol_export
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_export
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_export
    APPEND PROPERTY LABELS drawing_baseline)
  distclean(${CMAKE_CURRENT_BINARY_DIR}/qtcad_obol_export_tmp.g)
endif(TARGET test_qtcad_obol_export)

if(TARGET Qt6::Widgets)
  add_executable(test_qtcad_obol_snap ${TEST_QTCAD_OBOL_SNAP_SRCS})
  target_link_libraries(
    test_qtcad_obol_snap
    libqtcad
    libged
    libwdb
    librt
    libbu
    Qt6::Widgets
    )
endif(TARGET Qt6::Widgets)
if(NOT TARGET test_qtcad_obol_snap AND TARGET Qt5::Widgets)
  add_executable(test_qtcad_obol_snap ${TEST_QTCAD_OBOL_SNAP_SRCS})
  target_link_libraries(
    test_qtcad_obol_snap
    libqtcad
    libged
    libwdb
    librt
    libbu
    Qt5::Widgets
    )
endif(NOT TARGET test_qtcad_obol_snap AND TARGET Qt5::Widgets)
if(TARGET test_qtcad_obol_snap)
  target_compile_definitions(test_qtcad_obol_snap PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  brlcad_add_test(NAME test_qtcad_obol_snap COMMAND test_qtcad_obol_snap)
  set_property(TEST test_qtcad_obol_snap
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_snap
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_snap
    APPEND PROPERTY LABELS drawing_baseline)
  distclean(${CMAKE_CURRENT_BINARY_DIR}/qtcad_obol_snap_tmp.g)
endif(TARGET test_qtcad_obol_snap)

if(TARGET Qt6::Widgets)
  add_executable(test_qtcad_obol_measure ${TEST_QTCAD_OBOL_MEASURE_SRCS})
  target_link_libraries(
    test_qtcad_obol_measure
    libqtcad
    libged
    libwdb
    librt
    libbg
    libbu
    Qt6::Widgets
    )
endif(TARGET Qt6::Widgets)
if(NOT TARGET test_qtcad_obol_measure AND TARGET Qt5::Widgets)
  add_executable(test_qtcad_obol_measure ${TEST_QTCAD_OBOL_MEASURE_SRCS})
  target_link_libraries(
    test_qtcad_obol_measure
    libqtcad
    libged
    libwdb
    librt
    libbg
    libbu
    Qt5::Widgets
    )
endif(NOT TARGET test_qtcad_obol_measure AND TARGET Qt5::Widgets)
if(TARGET test_qtcad_obol_measure)
  target_compile_definitions(test_qtcad_obol_measure PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  brlcad_add_test(NAME test_qtcad_obol_measure COMMAND test_qtcad_obol_measure)
  set_property(TEST test_qtcad_obol_measure
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_measure
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_measure
    APPEND PROPERTY LABELS drawing_baseline)
  distclean(${CMAKE_CURRENT_BINARY_DIR}/qtcad_obol_measure_tmp.g)
endif(TARGET test_qtcad_obol_measure)

if(TARGET Qt6::Widgets)
  add_executable(test_qtcad_obol_real_models ${TEST_QTCAD_OBOL_REAL_MODELS_SRCS})
  target_link_libraries(
    test_qtcad_obol_real_models
    libqtcad
    libged
    libicv
    librt
    libbu
    Qt6::Widgets
    )
endif(TARGET Qt6::Widgets)
if(NOT TARGET test_qtcad_obol_real_models AND TARGET Qt5::Widgets)
  add_executable(test_qtcad_obol_real_models ${TEST_QTCAD_OBOL_REAL_MODELS_SRCS})
  target_link_libraries(
    test_qtcad_obol_real_models
    libqtcad
    libged
    libicv
    librt
    libbu
    Qt5::Widgets
    )
endif(NOT TARGET test_qtcad_obol_real_models AND TARGET Qt5::Widgets)
if(TARGET test_qtcad_obol_real_models)
  target_compile_definitions(test_qtcad_obol_real_models PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  if(TARGET pinewood.g)
    add_dependencies(test_qtcad_obol_real_models pinewood.g)
  endif(TARGET pinewood.g)
  if(TARGET havoc.g)
    add_dependencies(test_qtcad_obol_real_models havoc.g)
  endif(TARGET havoc.g)
  if(TARGET m35.g)
    add_dependencies(test_qtcad_obol_real_models m35.g)
  endif(TARGET m35.g)
  if(TARGET Generic_Twin.g)
    add_dependencies(test_qtcad_obol_real_models Generic_Twin.g)
  endif(TARGET Generic_Twin.g)
  foreach(_qtcad_real_model_case
      pinewood_wire
      pinewood_shaded
      havoc_wire
      m35_wire_interactions
      m35_wire_pick_all_stress
      m35_deferred_wire_color
      m35_deferred_shaded_color
	  nist_pmi1_shaded
	  nist_pmi2_shaded
	  nist_pmi3_shaded
	  nist_pmi4_shaded
	  nist_pmi5_shaded
	  nist_pmi6_shaded
	  nist_pmi11_shaded
	  nist_pmi7_10_shaded
	  nist_pmi7_10_wire
      m35_color_table_mutation)
    brlcad_add_test(
      NAME test_qtcad_obol_real_models_${_qtcad_real_model_case}
      COMMAND test_qtcad_obol_real_models ${_qtcad_real_model_case})
    set_tests_properties(test_qtcad_obol_real_models_${_qtcad_real_model_case}
      PROPERTIES TIMEOUT 180)
    set_property(TEST test_qtcad_obol_real_models_${_qtcad_real_model_case}
      APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
    set_property(TEST test_qtcad_obol_real_models_${_qtcad_real_model_case}
      APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
    # These rows intentionally validate the software/OSMesa endpoint.  Do not
    # let a developer's graphical DISPLAY leak into Qt/Obol backend discovery;
    # System-GL behavior has its own explicit graphical matrix.
    set_property(TEST test_qtcad_obol_real_models_${_qtcad_real_model_case}
      APPEND PROPERTY ENVIRONMENT DISPLAY=)
    set_property(TEST test_qtcad_obol_real_models_${_qtcad_real_model_case}
      APPEND PROPERTY ENVIRONMENT XAUTHORITY=)
    set_property(TEST test_qtcad_obol_real_models_${_qtcad_real_model_case}
      APPEND PROPERTY LABELS drawing_baseline)
  endforeach()
endif(TARGET test_qtcad_obol_real_models)

if(TARGET Qt6::Widgets)
  add_executable(test_qtcad_obol_progressive_lod ${TEST_QTCAD_OBOL_PROGRESSIVE_LOD_SRCS})
  target_link_libraries(
    test_qtcad_obol_progressive_lod
    libqtcad
    libged
    librt
    libbu
    Qt6::Widgets
    )
endif(TARGET Qt6::Widgets)
if(NOT TARGET test_qtcad_obol_progressive_lod AND TARGET Qt5::Widgets)
  add_executable(test_qtcad_obol_progressive_lod ${TEST_QTCAD_OBOL_PROGRESSIVE_LOD_SRCS})
  target_link_libraries(
    test_qtcad_obol_progressive_lod
    libqtcad
    libged
    librt
    libbu
    Qt5::Widgets
    )
endif(NOT TARGET test_qtcad_obol_progressive_lod AND TARGET Qt5::Widgets)
if(TARGET test_qtcad_obol_progressive_lod)
  target_compile_definitions(test_qtcad_obol_progressive_lod PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  # Earlier staged/startup-refine fixtures encoded superseded controller-level
  # AABB/OBB and source-node expansion designs.  The production contract is
  # one compact occurrence registry with per-leaf box-to-view-mesh refinement,
  # exercised end-to-end by the automatic cold-start case below.
  brlcad_add_test(
    NAME test_qtcad_obol_progressive_lod_generic_twin_startup_auto_expand
    COMMAND test_qtcad_obol_progressive_lod --startup-auto-expand "${CMAKE_BINARY_DIR}/share/db/faa/Generic_Twin.g" - all direct 4 1)
  set_tests_properties(test_qtcad_obol_progressive_lod_generic_twin_startup_auto_expand PROPERTIES TIMEOUT 120)
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_startup_auto_expand
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_startup_auto_expand
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_startup_auto_expand
    APPEND PROPERTY ENVIRONMENT QTCAD_OBOL_FORCE_OSMESA=1)
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_startup_auto_expand
    APPEND PROPERTY LABELS drawing_baseline)
  brlcad_add_test(
    NAME test_qtcad_obol_progressive_lod_generic_twin_wire
    COMMAND test_qtcad_obol_progressive_lod --startup-wire-lod "${CMAKE_BINARY_DIR}/share/db/faa/Generic_Twin.g" - all direct 4 1)
  set_tests_properties(test_qtcad_obol_progressive_lod_generic_twin_wire PROPERTIES TIMEOUT 120)
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_wire
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_wire
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_wire
    APPEND PROPERTY ENVIRONMENT QTCAD_OBOL_FORCE_OSMESA=1)
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_wire
    APPEND PROPERTY LABELS drawing_baseline)
  distclean(${CMAKE_CURRENT_BINARY_DIR}/qtcad_obol_progressive_lod_tmp.g)
  distclean(${CMAKE_CURRENT_BINARY_DIR}/qtcad_obol_progressive_lod_cache)
  set_tests_properties(
    test_qtcad_obol_progressive_lod_generic_twin_startup_auto_expand
    test_qtcad_obol_progressive_lod_generic_twin_wire
    PROPERTIES RESOURCE_LOCK qtcad_obol_progressive_lod_fixture)
endif(TARGET test_qtcad_obol_progressive_lod)

foreach(_qtcad_obol_offscreen_test
    test_qtcad_obol_controller
    test_qtcad_obol_window_host
    test_qtcad_obol_draw_sync
    test_qtcad_obol_production_flow_lifecycle
    test_qtcad_obol_edit_preview
    test_qtcad_obol_faceplate_sync
    test_qtcad_obol_selection_sync
    test_qtcad_obol_pick
    test_qtcad_obol_export
    test_qtcad_obol_snap
    test_qtcad_obol_measure
    test_qtcad_obol_progressive_lod_generic_twin_startup_auto_expand
    test_qtcad_obol_progressive_lod_generic_twin_wire
    test_qtcad_obol_real_models_pinewood_wire
    test_qtcad_obol_real_models_pinewood_shaded
    test_qtcad_obol_real_models_havoc_wire
    test_qtcad_obol_real_models_m35_wire_interactions
    test_qtcad_obol_real_models_m35_wire_pick_all_stress
    test_qtcad_obol_real_models_m35_deferred_wire_color
    test_qtcad_obol_real_models_m35_deferred_shaded_color
    test_qtcad_obol_real_models_m35_color_table_mutation)
  if(TEST ${_qtcad_obol_offscreen_test})
    set_property(TEST ${_qtcad_obol_offscreen_test}
      APPEND PROPERTY ENVIRONMENT QTCAD_OBOL_FORCE_OSMESA=1)
  endif()
endforeach()

if(TARGET Qt6::Widgets)
  set(QGMODEL_TEST_DATA_G "${CMAKE_SOURCE_DIR}/src/libged/tests/draw/moss.g")
  # Phase 8: standalone libqtcad plugin reuse test.
  # Deliberately links only libqtcad + Qt6::Widgets (NOT libged/qged) to
  # prove the plugin framework is host-agnostic.
  add_executable(test_qtcad_plugin_reuse ${TEST_QTCAD_PLUGIN_REUSE_SRCS})
  target_link_libraries(
    test_qtcad_plugin_reuse
    libqtcad
    libbu
    Qt6::Widgets
    )
  target_compile_definitions(test_qtcad_plugin_reuse PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  # The test scans the qged plugin directory at runtime; add dependencies on
  # the Phase 7 extension plugins so they are built before the test runs.
  if(TARGET qged_host_status_command)
    add_dependencies(test_qtcad_plugin_reuse
      qged_host_status_command
      qged_host_status_panel
      qged_host_status_dialog
    )
  endif()
  if(TARGET qged_ell)
    add_dependencies(test_qtcad_plugin_reuse qged_ell)
  endif()
  brlcad_add_test(NAME test_qtcad_plugin_reuse
    COMMAND test_qtcad_plugin_reuse "${CMAKE_BINARY_DIR}/libexec/qged")
  set_property(TEST test_qtcad_plugin_reuse
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_plugin_reuse
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_plugin_reuse
    APPEND PROPERTY LABELS drawing_baseline)

  if(TARGET test_qgmodel_model_tester)
    brlcad_add_test(NAME test_qgmodel_model_tester
      COMMAND test_qgmodel_model_tester "${QGMODEL_TEST_DATA_G}")
    set_property(TEST test_qgmodel_model_tester
      APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
    set_property(TEST test_qgmodel_model_tester
      APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
    set_property(TEST test_qgmodel_model_tester
      APPEND PROPERTY LABELS drawing_baseline)
  endif()
  if(TARGET test_qgmodel_scale)
    brlcad_add_test(NAME test_qgmodel_scale
      COMMAND test_qgmodel_scale "${QGMODEL_TEST_DATA_G}")
    set_property(TEST test_qgmodel_scale
      APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
    set_property(TEST test_qgmodel_scale
      APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
    set_property(TEST test_qgmodel_scale
      APPEND PROPERTY LABELS drawing_baseline)
    brlcad_add_test(NAME slow-test_qgmodel_scale_16384
      COMMAND test_qgmodel_scale "${QGMODEL_TEST_DATA_G}")
    set_property(TEST slow-test_qgmodel_scale_16384
      PROPERTY SKIP_RETURN_CODE 123)
    set_property(TEST slow-test_qgmodel_scale_16384
      APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
    set_property(TEST slow-test_qgmodel_scale_16384
      APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
    set_property(TEST slow-test_qgmodel_scale_16384
      APPEND PROPERTY ENVIRONMENT BRLCAD_RUN_SLOW_TESTS=1)
    set_property(TEST slow-test_qgmodel_scale_16384
      APPEND PROPERTY ENVIRONMENT BRLCAD_QGMODEL_SCALE_CHILDREN=16384)
    set_property(TEST slow-test_qgmodel_scale_16384
      APPEND PROPERTY ENVIRONMENT BRLCAD_QGMODEL_DRAW_CHILDREN=64)
    set_property(TEST slow-test_qgmodel_scale_16384
      APPEND PROPERTY LABELS drawing_baseline)
    set_property(TEST slow-test_qgmodel_scale_16384
      APPEND PROPERTY LABELS slow)
  endif()
endif(TARGET Qt6::Widgets)

# On some platforms distclean has additional cleanup to do

cmakefiles(
  CMakeLists.txt
  ${TEST_QTCAD_SESSION_SRCS}
  ${TEST_QTCAD_MODEL_INVARIANTS_SRCS}
  ${TEST_QTCAD_PRIMITIVE_EDIT_SRCS}
  ${TEST_QTCAD_CONSOLE_SRCS}
  ${TEST_QTCAD_PLUGIN_MANAGER_SRCS}
  ${TEST_QTCAD_PLUGIN_CONTEXT_SRCS}
  ${TEST_QTCAD_PLUGIN_BASE_HDRS}
  ${TEST_QTCAD_PLUGIN_VALID_SRCS}
  ${TEST_QTCAD_PLUGIN_MISSING_API_SRCS}
  ${TEST_QTCAD_PLUGIN_PANEL_SRCS}
  ${TEST_QTCAD_PLUGIN_DIALOG_SRCS}
  ${TEST_QTCAD_EVENT_PLAYER_SRCS}
  ${TEST_QTCAD_PLUGIN_METADATA}
  qtcad_obol_nist_matrix.sh
  ${TEST_QTCAD_OBOL_CONTROLLER_SRCS}
  ${TEST_QTCAD_OBOL_DRAW_SYNC_SRCS}
  ${TEST_QTCAD_OBOL_EDIT_PREVIEW_SRCS}
  ${TEST_QTCAD_OBOL_FACEPLATE_SYNC_SRCS}
  ${TEST_QTCAD_OBOL_MEASURE_SRCS}
  ${TEST_QTCAD_OBOL_PICK_SRCS}
  ${TEST_QTCAD_OBOL_PROGRESSIVE_LOD_SRCS}
  ${TEST_QTCAD_OBOL_REAL_MODELS_SRCS}
  ${TEST_QTCAD_OBOL_SELECTION_SYNC_SRCS}
  ${TEST_QTCAD_OBOL_SNAP_SRCS}
  ${TEST_QTCAD_OBOL_WINDOW_HOST_SRCS}
  ${TEST_QTCAD_PLUGIN_REUSE_SRCS}
  ${TEST_QGMODEL_MODEL_TESTER_SRCS}
  ${TEST_QGMODEL_SCALE_SRCS}
  )

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
