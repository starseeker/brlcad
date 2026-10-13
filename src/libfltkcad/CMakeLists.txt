set(FLTKCAD_SRCS
  CommandHistory.cpp
  Fl_Cad_Quad_View.cpp
  Fl_Cad_Terminal.cpp
  Fl_Cad_View.cpp
  Fl_Ged_Accessibility.cpp
  Fl_Ged_Attributes.cpp
  Fl_Ged_Combination_Editor.cpp
  Fl_Ged_Create_Primitive.cpp
  Fl_Ged_Help.cpp
  Fl_Ged_Hierarchy.cpp
  Fl_Ged_Interaction.cpp
  Fl_Ged_Preferences.cpp
  Fl_Ged_Primitive_Editor.cpp
  Fl_Ged_Render_Tools.cpp
  Fl_Ged_Search.cpp
  Fl_Ged_Sketch_Editor.cpp
  Fl_Ged_Specialist_Tools.cpp
  Fl_Ged_Tcl.cpp
  Fl_Ged_Transform_Editor.cpp
  Fl_Ged_View_Controls.cpp
  Fl_Litehtml_View.cpp
  Uri.cpp
)

if(BRLCAD_ENABLE_FLTK)
  option(FLTKCAD_BUILD_EXAMPLES "Build libfltkcad example applications" OFF)
  mark_as_advanced(FLTKCAD_BUILD_EXAMPLES)

  set(FLTKCAD_PRIVATE_LIBS fltk::gl fltk::images litehtml)
  set(FLTKCAD_PUBLIC_LIBS ${libfltkcad_deps})
  list(REMOVE_ITEM FLTKCAD_PUBLIC_LIBS libtclcad)
  set(FLTKCAD_PRIVATE_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
  if(BRLCAD_ENABLE_TCL)
    find_package_tcl(REQUIRED)
    list(APPEND FLTKCAD_PRIVATE_LIBS libtclcad ${TCL_LIBRARY})
    list(APPEND FLTKCAD_PRIVATE_INCLUDES ${TCL_INCLUDE_PATH})
  endif()

  # Widget implementations have independent file-local helper names.
  brlcad_addlib(libfltkcad "${FLTKCAD_SRCS}" "" "${FLTKCAD_PRIVATE_INCLUDES}" SHARED NO_UNITY
    PUBLIC_LIBS ${FLTKCAD_PUBLIC_LIBS} fltk::fltk
    PRIVATE_LIBS ${FLTKCAD_PRIVATE_LIBS}
  )
  set_target_properties(libfltkcad PROPERTIES VERSION 20.0.1 SOVERSION 20)

  foreach(fltkcad_target libfltkcad libfltkcad-obj libfltkcad-brlcad-obj)
    if(TARGET ${fltkcad_target})
      target_compile_features(${fltkcad_target} PUBLIC cxx_std_17)
      target_compile_definitions(${fltkcad_target} PRIVATE
        FLTKCAD_BRLCAD_VERSION="${BRLCAD_VERSION}")
      if(BRLCAD_ENABLE_TCL)
        target_compile_definitions(${fltkcad_target} PRIVATE FLTKCAD_HAVE_TCL=1)
      endif()
    endif()
  endforeach()

  add_dependencies(libfltkcad dm_plugins ged_plugins)

  install(FILES README.md docs/html.md docs/terminal.md docs/display.md
    DESTINATION ${DOC_DIR}/libfltkcad)
endif()

add_subdirectory(examples)
add_subdirectory(tests)

cmakefiles(CMakeLists.txt ${FLTKCAD_SRCS} Uri.h README.md
  docs/html.md docs/terminal.md docs/display.md)
