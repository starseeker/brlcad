if(BRLCAD_ENABLE_QT)
  include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${BRLCAD_BINARY_DIR}/${INCLUDE_DIR}
    ${BRLCAD_SOURCE_DIR}/include
    ${QTCAD_INCLUDE_DIRS}
    ${GED_INCLUDE_DIRS}
    ${RT_INCLUDE_DIRS}
    ${BU_INCLUDE_DIRS}
  )
endif(BRLCAD_ENABLE_QT)

set(sketch_srcs QSketch.cpp sketch.cpp ../qged_edit_preview_util.cpp)
set(sketch_moc_hdrs QSketch.h EditSketchFactory.h)

if(BRLCAD_ENABLE_QT)
  if(Qt6Widgets_FOUND)
    qt6_wrap_cpp(sketch_moc_srcs ${sketch_moc_hdrs})
    qt6_add_resources(icon_qrc icon.qrc)
  else()
    qt5_wrap_cpp(sketch_moc_srcs ${sketch_moc_hdrs})
    qt5_add_resources(icon_qrc icon.qrc)
  endif(Qt6Widgets_FOUND)

  qged_plugin(qged_sketch SHARED ${sketch_srcs} ${sketch_moc_srcs} ${icon_qrc})
  if(Qt6Widgets_FOUND)
    target_link_libraries(qged_sketch libqtcad libBObol libged librt libbu Qt6::Widgets)
  else()
    target_link_libraries(qged_sketch libqtcad libBObol libged librt libbu Qt5::Widgets)
  endif(Qt6Widgets_FOUND)
  set_property(TARGET qged_sketch APPEND PROPERTY COMPILE_DEFINITIONS BRLCADBUILD HAVE_CONFIG_H)
  validate_style(qged_sketch sketch.cpp)
  plugin_setup(qged_sketch qged)

  distclean(
    icon.qrc.depends
    moc_QSketch.cpp_parameters
    moc_EditSketchFactory.cpp_parameters
  )
endif(BRLCAD_ENABLE_QT)

set(
  ignore_srcs
  CMakeLists.txt
  ${sketch_srcs}
  ${sketch_moc_hdrs}
  icon.qrc
  sketch_edit.svg
  sketch.json
)

cmakefiles(${ignore_srcs})

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
