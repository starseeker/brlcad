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

set(arb_srcs QArb.cpp arb.cpp ../qged_edit_preview_util.cpp)
set(arb_moc_hdrs QArb.h EditArbFactory.h)

if(BRLCAD_ENABLE_QT)
  if(Qt6Widgets_FOUND)
    qt6_wrap_cpp(arb_moc_srcs ${arb_moc_hdrs})
    qt6_add_resources(icon_qrc icon.qrc)
  else()
    qt5_wrap_cpp(arb_moc_srcs ${arb_moc_hdrs})
    qt5_add_resources(icon_qrc icon.qrc)
  endif(Qt6Widgets_FOUND)

  qged_plugin(qged_arb SHARED ${arb_srcs} ${arb_moc_srcs} ${icon_qrc})
  if(Qt6Widgets_FOUND)
    target_link_libraries(
      qged_arb
      libqtcad
      libBObol
      libged
      librt
      libbu
      Qt6::Widgets
    )
  else()
    target_link_libraries(
      qged_arb
      libqtcad
      libBObol
      libged
      librt
      libbu
      Qt5::Widgets
    )
  endif(Qt6Widgets_FOUND)
  set_property(TARGET qged_arb APPEND PROPERTY COMPILE_DEFINITIONS BRLCADBUILD HAVE_CONFIG_H)
  validate_style(qged_arb arb.cpp)
  plugin_setup(qged_arb qged)

  distclean(
    icon.qrc.depends
    moc_QArb.cpp_parameters
    moc_EditArbFactory.cpp_parameters
  )
endif(BRLCAD_ENABLE_QT)

set(
  ignore_srcs
  CMakeLists.txt
  ${arb_srcs}
  ${arb_moc_hdrs}
  icon.qrc
  arb.svg
  arb.json
)

cmakefiles(${ignore_srcs})

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
