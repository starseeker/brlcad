if(BRLCAD_ENABLE_QT)
  include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${BRLCAD_BINARY_DIR}/${INCLUDE_DIR}
    ${BRLCAD_SOURCE_DIR}/include
    ${QTCAD_INCLUDE_DIRS}
    ${GED_INCLUDE_DIRS}
    ${RT_INCLUDE_DIRS}
    ${BREP_INCLUDE_DIRS}
    ${BU_INCLUDE_DIRS}
  )
endif(BRLCAD_ENABLE_QT)

set(brep_srcs QBrep.cpp brep.cpp ../qged_edit_preview_util.cpp)
set(brep_moc_hdrs QBrep.h EditBrepFactory.h)

if(BRLCAD_ENABLE_QT)
  if(Qt6Widgets_FOUND)
    qt6_wrap_cpp(brep_moc_srcs ${brep_moc_hdrs})
  else()
    qt5_wrap_cpp(brep_moc_srcs ${brep_moc_hdrs})
  endif(Qt6Widgets_FOUND)

  qged_plugin(qged_brep SHARED ${brep_srcs} ${brep_moc_srcs})
  if(Qt6Widgets_FOUND)
    target_link_libraries(
      qged_brep
      libqtcad
      libged
      librt
      libbrep
      libbu
      Qt6::Widgets
    )
  else()
    target_link_libraries(
      qged_brep
      libqtcad
      libged
      librt
      libbrep
      libbu
      Qt5::Widgets
    )
  endif(Qt6Widgets_FOUND)
  set_property(TARGET qged_brep APPEND PROPERTY COMPILE_DEFINITIONS BRLCADBUILD HAVE_CONFIG_H)
  validate_style(qged_brep brep.cpp)
  plugin_setup(qged_brep qged)

  distclean(
    moc_QBrep.cpp_parameters
    moc_EditBrepFactory.cpp_parameters
  )
endif(BRLCAD_ENABLE_QT)

set(
  ignore_srcs
  CMakeLists.txt
  ${brep_srcs}
  ${brep_moc_hdrs}
  brep.json
)

cmakefiles(${ignore_srcs})

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
