include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)

set(QGMODEL_SRCS qgmodel.cpp)
set(QGVIEW_SRCS qgview.cpp)
set(QSKETCH_SRCS qsketch.cpp)
set(TEST_QTCAD_PLUGIN_REUSE_SRCS test_qtcad_plugin_reuse.cpp)
set(TEST_QTCAD_OBOL_CONTROLLER_SRCS test_qtcad_obol_controller.cpp)
set(TEST_QTCAD_OBOL_DRAW_SYNC_SRCS test_qtcad_obol_draw_sync.cpp)
set(TEST_QTCAD_OBOL_EDIT_PREVIEW_SRCS test_qtcad_obol_edit_preview.cpp)
set(TEST_QTCAD_OBOL_EXPORT_SRCS test_qtcad_obol_export.cpp)
set(TEST_QTCAD_OBOL_FACEPLATE_SYNC_SRCS test_qtcad_obol_faceplate_sync.cpp)
set(TEST_QTCAD_OBOL_MEASURE_SRCS test_qtcad_obol_measure.cpp)
set(TEST_QTCAD_OBOL_OVERLAY_SYNC_SRCS test_qtcad_obol_overlay_sync.cpp)
set(TEST_QTCAD_OBOL_PICK_SRCS test_qtcad_obol_pick.cpp)
set(TEST_QTCAD_OBOL_PROGRESSIVE_LOD_SRCS test_qtcad_obol_progressive_lod.cpp)
set(TEST_QTCAD_OBOL_REAL_MODELS_SRCS test_qtcad_obol_real_models.cpp)
set(TEST_QTCAD_OBOL_SELECTION_SYNC_SRCS test_qtcad_obol_selection_sync.cpp)
set(TEST_QTCAD_OBOL_SNAP_SRCS test_qtcad_obol_snap.cpp)
set(TEST_QTCAD_OBOL_WINDOW_HOST_SRCS test_qtcad_obol_window_host.cpp)
set(TEST_QGMODEL_MODEL_TESTER_SRCS test_qgmodel_model_tester.cpp)
set(TEST_QGMODEL_SCALE_SRCS test_qgmodel_scale.cpp)

if(TARGET Qt6::Test)
  add_executable(qgmodel qgmodel.cpp)
  target_link_libraries(
    qgmodel
    libqtcad
    libged
    librt
    libbu
    Qt6::Widgets
    Qt6::Test
    )

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
  add_executable(qgmodel ${QGMODEL_SRCS})
  target_link_libraries(
    qgmodel
    libqtcad
    libged
    librt
    libbu
    Qt5::Widgets
    Qt5::Test
    )

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
if(TARGET qgmodel)
  target_compile_definitions(qgmodel PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  install(TARGETS qgmodel RUNTIME DESTINATION ${BIN_DIR} LIBRARY DESTINATION ${LIB_DIR} ARCHIVE DESTINATION ${LIB_DIR})
endif(TARGET qgmodel)
if(TARGET test_qgmodel_model_tester)
  target_compile_definitions(test_qgmodel_model_tester PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  install(TARGETS test_qgmodel_model_tester RUNTIME DESTINATION ${BIN_DIR} LIBRARY DESTINATION ${LIB_DIR} ARCHIVE DESTINATION ${LIB_DIR})
endif(TARGET test_qgmodel_model_tester)
if(TARGET test_qgmodel_scale)
  target_compile_definitions(test_qgmodel_scale PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  install(TARGETS test_qgmodel_scale RUNTIME DESTINATION ${BIN_DIR} LIBRARY DESTINATION ${LIB_DIR} ARCHIVE DESTINATION ${LIB_DIR})
endif(TARGET test_qgmodel_scale)

if(TARGET Qt6::Widgets)
  add_executable(qgview ${QGVIEW_SRCS})
  target_link_libraries(
    qgview
    libqtcad
    libged
    librt
    libbu
    Qt6::Widgets
    )
endif(TARGET Qt6::Widgets)
if(TARGET Qt5::Widgets)
  add_executable(qgview ${QGVIEW_SRCS})
  target_link_libraries(
    qgview
    libqtcad
    libged
    librt
    libbu
    Qt5::Widgets
    )
endif(TARGET Qt5::Widgets)
if(TARGET qgview)
  target_compile_definitions(qgview PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  install(TARGETS qgview RUNTIME DESTINATION ${BIN_DIR} LIBRARY DESTINATION ${LIB_DIR} ARCHIVE DESTINATION ${LIB_DIR})
endif(TARGET qgview)

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
  set_property(TEST test_qtcad_obol_controller
    APPEND PROPERTY LABELS drawing_baseline)
endif(TARGET test_qtcad_obol_controller)

if(TARGET Qt6::Widgets)
  add_executable(test_qtcad_obol_window_host ${TEST_QTCAD_OBOL_WINDOW_HOST_SRCS})
  target_link_libraries(
    test_qtcad_obol_window_host
    libqtcad
    libimgstream
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
    libbu
    Qt5::Widgets
    )
endif(NOT TARGET test_qtcad_obol_window_host AND TARGET Qt5::Widgets)
if(TARGET test_qtcad_obol_window_host)
  target_compile_definitions(test_qtcad_obol_window_host PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  brlcad_add_test(NAME test_qtcad_obol_window_host COMMAND test_qtcad_obol_window_host)
  set_property(TEST test_qtcad_obol_window_host
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_window_host
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_window_host
    APPEND PROPERTY LABELS drawing_baseline)
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
  set_property(TEST test_qtcad_obol_draw_sync
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_draw_sync
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_draw_sync
    APPEND PROPERTY LABELS drawing_baseline)
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
  add_executable(test_qtcad_obol_overlay_sync ${TEST_QTCAD_OBOL_OVERLAY_SYNC_SRCS})
  target_link_libraries(
    test_qtcad_obol_overlay_sync
    libqtcad
    libged
    libwdb
    librt
    libbg
    libbu
    Qt6::Widgets
    )
endif(TARGET Qt6::Widgets)
if(NOT TARGET test_qtcad_obol_overlay_sync AND TARGET Qt5::Widgets)
  add_executable(test_qtcad_obol_overlay_sync ${TEST_QTCAD_OBOL_OVERLAY_SYNC_SRCS})
  target_link_libraries(
    test_qtcad_obol_overlay_sync
    libqtcad
    libged
    libwdb
    librt
    libbg
    libbu
    Qt5::Widgets
    )
endif(NOT TARGET test_qtcad_obol_overlay_sync AND TARGET Qt5::Widgets)
if(TARGET test_qtcad_obol_overlay_sync)
  target_compile_definitions(test_qtcad_obol_overlay_sync PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  brlcad_add_test(NAME test_qtcad_obol_overlay_sync COMMAND test_qtcad_obol_overlay_sync)
  set_property(TEST test_qtcad_obol_overlay_sync
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_overlay_sync
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_overlay_sync
    APPEND PROPERTY LABELS drawing_baseline)
  distclean(${CMAKE_CURRENT_BINARY_DIR}/qtcad_obol_overlay_sync_tmp.g)
endif(TARGET test_qtcad_obol_overlay_sync)

if(TARGET Qt6::Widgets)
  add_executable(test_qtcad_obol_real_models ${TEST_QTCAD_OBOL_REAL_MODELS_SRCS})
  target_link_libraries(
    test_qtcad_obol_real_models
    libqtcad
    libged
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
      m35_wire_pick_all_stress)
    brlcad_add_test(
      NAME test_qtcad_obol_real_models_${_qtcad_real_model_case}
      COMMAND test_qtcad_obol_real_models ${_qtcad_real_model_case})
    set_tests_properties(test_qtcad_obol_real_models_${_qtcad_real_model_case}
      PROPERTIES TIMEOUT 180)
    set_property(TEST test_qtcad_obol_real_models_${_qtcad_real_model_case}
      APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
    set_property(TEST test_qtcad_obol_real_models_${_qtcad_real_model_case}
      APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
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
  brlcad_add_test(
    NAME test_qtcad_obol_progressive_lod
    COMMAND test_qtcad_obol_progressive_lod "${CMAKE_SOURCE_DIR}/src/libged/tests/draw/moss.g")
  set_tests_properties(test_qtcad_obol_progressive_lod PROPERTIES TIMEOUT 120)
  set_property(TEST test_qtcad_obol_progressive_lod
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_progressive_lod
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_progressive_lod
    APPEND PROPERTY ENVIRONMENT QTCAD_OBOL_FORCE_OSMESA=1)
  set_property(TEST test_qtcad_obol_progressive_lod
    APPEND PROPERTY LABELS drawing_baseline)
  brlcad_add_test(
    NAME test_qtcad_obol_progressive_lod_startup_refine
    COMMAND test_qtcad_obol_progressive_lod --startup-refine "${CMAKE_SOURCE_DIR}/src/libged/tests/draw/moss.g" all.g all.bot facetize 2 2)
  set_tests_properties(test_qtcad_obol_progressive_lod_startup_refine PROPERTIES TIMEOUT 120)
  set_property(TEST test_qtcad_obol_progressive_lod_startup_refine
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_progressive_lod_startup_refine
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_progressive_lod_startup_refine
    APPEND PROPERTY ENVIRONMENT QTCAD_OBOL_FORCE_OSMESA=1)
  set_property(TEST test_qtcad_obol_progressive_lod_startup_refine
    APPEND PROPERTY LABELS drawing_baseline)
  brlcad_add_test(
    NAME test_qtcad_obol_progressive_lod_generic_twin_startup
    COMMAND test_qtcad_obol_progressive_lod --startup-only "${CMAKE_BINARY_DIR}/share/db/faa/Generic_Twin.g" - all direct 1 1)
  set_tests_properties(test_qtcad_obol_progressive_lod_generic_twin_startup PROPERTIES TIMEOUT 45)
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_startup
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_startup
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_startup
    APPEND PROPERTY ENVIRONMENT QTCAD_OBOL_FORCE_OSMESA=1)
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_startup
    APPEND PROPERTY LABELS drawing_baseline)
  brlcad_add_test(
    NAME test_qtcad_obol_progressive_lod_generic_twin_startup_expand
    COMMAND test_qtcad_obol_progressive_lod --startup-expand --startup-expand-passes 4 "${CMAKE_BINARY_DIR}/share/db/faa/Generic_Twin.g" - all direct 4 1)
  set_tests_properties(test_qtcad_obol_progressive_lod_generic_twin_startup_expand PROPERTIES TIMEOUT 120)
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_startup_expand
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_startup_expand
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_startup_expand
    APPEND PROPERTY ENVIRONMENT QTCAD_OBOL_FORCE_OSMESA=1)
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_startup_expand
    APPEND PROPERTY LABELS drawing_baseline)
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
    NAME test_qtcad_obol_progressive_lod_generic_twin_wire_lod
    COMMAND test_qtcad_obol_progressive_lod --startup-wire-lod --startup-expand-passes 4 "${CMAKE_BINARY_DIR}/share/db/faa/Generic_Twin.g" - all direct 4 1)
  set_tests_properties(test_qtcad_obol_progressive_lod_generic_twin_wire_lod PROPERTIES TIMEOUT 180)
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_wire_lod
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_wire_lod
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_wire_lod
    APPEND PROPERTY ENVIRONMENT QTCAD_OBOL_FORCE_OSMESA=1)
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_wire_lod
    APPEND PROPERTY LABELS drawing_baseline)
  brlcad_add_test(
    NAME test_qtcad_obol_progressive_lod_generic_twin_shaded_lod
    COMMAND test_qtcad_obol_progressive_lod --shaded-lod "${CMAKE_BINARY_DIR}/share/db/faa/Generic_Twin.g" - all direct 4 1)
  set_tests_properties(test_qtcad_obol_progressive_lod_generic_twin_shaded_lod PROPERTIES TIMEOUT 180)
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_shaded_lod
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_shaded_lod
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_shaded_lod
    APPEND PROPERTY ENVIRONMENT QTCAD_OBOL_FORCE_OSMESA=1)
  set_property(TEST test_qtcad_obol_progressive_lod_generic_twin_shaded_lod
    APPEND PROPERTY LABELS drawing_baseline)
  distclean(${CMAKE_CURRENT_BINARY_DIR}/qtcad_obol_progressive_lod_tmp.g)
  distclean(${CMAKE_CURRENT_BINARY_DIR}/qtcad_obol_progressive_lod_cache)
  brlcad_add_test(
    NAME slow-test_qtcad_obol_progressive_lod_repeated_moss
    COMMAND test_qtcad_obol_progressive_lod --slow "${CMAKE_SOURCE_DIR}/src/libged/tests/draw/moss.g" all.g all.bot facetize 24 2 4)
  set_tests_properties(slow-test_qtcad_obol_progressive_lod_repeated_moss
    PROPERTIES TIMEOUT 300)
  set_property(TEST slow-test_qtcad_obol_progressive_lod_repeated_moss
    PROPERTY SKIP_RETURN_CODE 123)
  set_property(TEST slow-test_qtcad_obol_progressive_lod_repeated_moss
    APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
  set_property(TEST slow-test_qtcad_obol_progressive_lod_repeated_moss
    APPEND PROPERTY ENVIRONMENT QT_QPA_PLATFORM=offscreen)
  set_property(TEST slow-test_qtcad_obol_progressive_lod_repeated_moss
    APPEND PROPERTY ENVIRONMENT QTCAD_OBOL_FORCE_OSMESA=1)
  set_property(TEST slow-test_qtcad_obol_progressive_lod_repeated_moss
    APPEND PROPERTY LABELS drawing_baseline)
endif(TARGET test_qtcad_obol_progressive_lod)

foreach(_qtcad_obol_offscreen_test
    test_qtcad_obol_controller
    test_qtcad_obol_window_host
    test_qtcad_obol_draw_sync
    test_qtcad_obol_edit_preview
    test_qtcad_obol_faceplate_sync
    test_qtcad_obol_selection_sync
    test_qtcad_obol_pick
    test_qtcad_obol_export
    test_qtcad_obol_snap
    test_qtcad_obol_measure
    test_qtcad_obol_overlay_sync
    test_qtcad_obol_progressive_lod
    test_qtcad_obol_progressive_lod_startup_refine
    test_qtcad_obol_progressive_lod_generic_twin_startup
    test_qtcad_obol_progressive_lod_generic_twin_startup_expand
    test_qtcad_obol_progressive_lod_generic_twin_wire_lod
    test_qtcad_obol_progressive_lod_generic_twin_shaded_lod
    slow-test_qtcad_obol_progressive_lod_repeated_moss
    test_qtcad_obol_real_models_pinewood_wire
    test_qtcad_obol_real_models_pinewood_shaded
    test_qtcad_obol_real_models_havoc_wire
    test_qtcad_obol_real_models_m35_wire_interactions
    test_qtcad_obol_real_models_m35_wire_pick_all_stress)
  if(TEST ${_qtcad_obol_offscreen_test})
    set_property(TEST ${_qtcad_obol_offscreen_test}
      APPEND PROPERTY ENVIRONMENT QTCAD_OBOL_FORCE_OSMESA=1)
  endif()
endforeach()

if(TARGET Qt6::Widgets)
  set(QGMODEL_TEST_DATA_G "${CMAKE_SOURCE_DIR}/src/libged/tests/draw/moss.g")
  # qsketch - standalone .g sketch editing demo
  # qsketch.cpp defines QSketchEditWindow inline and uses #include "qsketch.moc"
  # at the bottom, so we need AUTOMOC ON rather than qt6_wrap_cpp.
  add_executable(qsketch ${QSKETCH_SRCS})
  set_target_properties(qsketch PROPERTIES AUTOMOC ON)
  target_link_libraries(
    qsketch
    libqtcad
    librt
    libbu
    Qt6::Widgets
    )

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
  endif()
endif(TARGET Qt6::Widgets)

if(TARGET Qt5::Widgets AND NOT TARGET Qt6::Widgets)
  add_executable(qsketch ${QSKETCH_SRCS})
  set_target_properties(qsketch PROPERTIES AUTOMOC ON)
  target_link_libraries(
    qsketch
    libqtcad
    librt
    libbu
    Qt5::Widgets
    )
endif(TARGET Qt5::Widgets AND NOT TARGET Qt6::Widgets)
if(TARGET qsketch)
  target_compile_definitions(qsketch PRIVATE BRLCADBUILD HAVE_CONFIG_H)
  install(TARGETS qsketch RUNTIME DESTINATION ${BIN_DIR} LIBRARY DESTINATION ${LIB_DIR} ARCHIVE DESTINATION ${LIB_DIR})
endif(TARGET qsketch)

cmakefiles(
  CMakeLists.txt
  ${QGMODEL_SRCS}
  ${QGVIEW_SRCS}
  ${QSKETCH_SRCS}
  ${TEST_QTCAD_OBOL_CONTROLLER_SRCS}
  ${TEST_QTCAD_OBOL_DRAW_SYNC_SRCS}
  ${TEST_QTCAD_OBOL_EDIT_PREVIEW_SRCS}
  ${TEST_QTCAD_OBOL_FACEPLATE_SYNC_SRCS}
  ${TEST_QTCAD_OBOL_MEASURE_SRCS}
  ${TEST_QTCAD_OBOL_OVERLAY_SYNC_SRCS}
  ${TEST_QTCAD_OBOL_PICK_SRCS}
  ${TEST_QTCAD_OBOL_PROGRESSIVE_LOD_SRCS}
  ${TEST_QTCAD_OBOL_REAL_MODELS_SRCS}
  ${TEST_QTCAD_OBOL_SELECTION_SYNC_SRCS}
  ${TEST_QTCAD_OBOL_SNAP_SRCS}
  ${TEST_QTCAD_PLUGIN_REUSE_SRCS}
  ${TEST_QGMODEL_MODEL_TESTER_SRCS}
  )

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
