if(BRLCAD_ENABLE_TCL AND BRLCAD_ENABLE_TK AND BRLCAD_ENABLE_OPENGL)
  set(TCL_ENABLE_TK ON)
  find_package_tcl(TCL)
  find_package(OpenGL)
  set(TKOBOL_PLATFORM_FOUND OFF)
  set(TKOBOL_PLATFORM_LIBS)
  set(TKOBOL_PLATFORM_DEFINITION)
  if("${TK_WINDOWING_SYSTEM}" STREQUAL "x11")
    find_package(X11)
    if(X11_FOUND)
      set(TKOBOL_PLATFORM_FOUND ON)
      set(TKOBOL_PLATFORM_LIBS ${X11_LIBRARIES})
      set(TKOBOL_PLATFORM_DEFINITION TKOBOL_X11)
    endif()
  elseif(WIN32)
    set(TKOBOL_PLATFORM_FOUND ON)
    set(TKOBOL_PLATFORM_LIBS ${OPENGL32_LIBRARY})
    set(TKOBOL_PLATFORM_DEFINITION TKOBOL_WGL)
  endif()

  if(TCL_FOUND AND TK_FOUND AND OPENGL_FOUND AND TKOBOL_PLATFORM_FOUND AND
      TARGET libtclcad)
    set(TKOBOL_SRCS dm-tkobol.cpp)
    include_directories(
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_SOURCE_DIR}/../..
      ${BRLCAD_BINARY_DIR}/include
      ${BRLCAD_SOURCE_DIR}/include
      ${BU_INCLUDE_DIRS}
      ${TCL_INCLUDE_PATH}
      ${TK_INCLUDE_PATH}
      ${X11_INCLUDE_DIR}
      ${OPENGL_INCLUDE_DIR}
    )

    add_definitions(-DDM_PLUGIN)
    dm_plugin_library(dm-tkobol SHARED ${TKOBOL_SRCS})
    target_link_libraries(
      dm-tkobol
      libtclcad
      libbrlobol
      libdm
      libbu
      ${TCL_LIBRARY}
      ${TK_LIBRARY}
      ${TKOBOL_PLATFORM_LIBS}
      ${OPENGL_LIBRARIES}
    )
    set_property(
      TARGET dm-tkobol
      APPEND
      PROPERTY COMPILE_DEFINITIONS BRLCADBUILD HAVE_CONFIG_H HAVE_TK ${TKOBOL_PLATFORM_DEFINITION}
    )
    validate_style(dm-tkobol "${TKOBOL_SRCS}")
    plugin_setup(dm-tkobol dm)
  endif()
endif()

cmakefiles(CMakeLists.txt dm-tkobol.cpp)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
