if(BRLCAD_ENABLE_TCL)
  if(BRLCAD_ENABLE_TK)
    set(TCL_ENABLE_TK ON)
  endif(BRLCAD_ENABLE_TK)
  find_package_tcl(REQUIRED)
  if(EXISTS ${CMAKE_BINARY_DIR}/${LIB_DIR}/itcl3.4)
    set(ITCL_VERSION "3.4")
  endif(EXISTS ${CMAKE_BINARY_DIR}/${LIB_DIR}/itcl3.4)
  if(EXISTS ${CMAKE_BINARY_DIR}/${LIB_DIR}/itk3.4)
    set(ITK_VERSION "3.4")
  endif(EXISTS ${CMAKE_BINARY_DIR}/${LIB_DIR}/itk3.4)
  if(EXISTS ${CMAKE_BINARY_DIR}/${LIB_DIR}/Iwidgets4.1.1)
    set(IWIDGETS_VERSION "4.1.1")
  endif(EXISTS ${CMAKE_BINARY_DIR}/${LIB_DIR}/Iwidgets4.1.1)
endif()

set(
  LIBTCLCAD_SRCS
  auto_path.c
  bn.c
  bu.c
  command_io.cpp
  commands.c
  cmdhist.c
  dm.c
  eval.c
  fb.c
  fbserv.c
  global.c
  init.c
  mouse.c
  polygons.c
  rt.c
  tkImgFmtPIX.c
  view/arrows.c
  view/autoview.c
  view/axes.c
  view/draw.c
  view/faceplate.c
  view/labels.c
  view/lines.c
  view/refresh.c
  view/util.c
  view_data.c
  wrapper.c
)

set(TCLCAD_TKOBOL_HOST OFF)
set(TCLCAD_TKOBOL_DEFINITIONS)
set(TCLCAD_TKOBOL_DM_DEFINITION)
if(BRLCAD_ENABLE_TK AND BRLCAD_ENABLE_OPENGL)
  find_package(OpenGL)
  if("${TK_WINDOWING_SYSTEM}" STREQUAL "x11")
    find_package(X11)
    if(OPENGL_FOUND AND X11_FOUND)
      set(TCLCAD_TKOBOL_HOST ON)
      set(TCLCAD_TKOBOL_DEFINITIONS TOGL_X11 TOGL_USE_OVERLAY=0)
      set(TCLCAD_TKOBOL_DM_DEFINITION TKOBOL_X11)
    endif()
  elseif(WIN32)
    if(OPENGL_FOUND)
      set(TCLCAD_TKOBOL_HOST ON)
      set(TCLCAD_TKOBOL_DEFINITIONS TOGL_WGL TOGL_USE_OVERLAY=0)
      set(TCLCAD_TKOBOL_DM_DEFINITION TKOBOL_WGL)
    endif()
  endif()
  if(TCLCAD_TKOBOL_HOST)
    set(TCLCAD_TKOBOL_HOST ON)
    list(APPEND LIBTCLCAD_SRCS
      tkobol/tk-obol-host.cpp
      tkobol/vendor/togl/togl.c
    )
    set_property(
      SOURCE tkobol/tk-obol-host.cpp
      APPEND
      PROPERTY COMPILE_DEFINITIONS ${TCLCAD_TKOBOL_DM_DEFINITION}
    )
    if(WIN32)
      list(APPEND LIBTCLCAD_SRCS tkobol/vendor/togl/toglFont.c)
      set_property(
        SOURCE tkobol/vendor/togl/toglFont.c
        APPEND
        PROPERTY COMPILE_DEFINITIONS TOGL_WGL
      )
    endif()
    set_property(
      SOURCE tkobol/vendor/togl/togl.c
      APPEND
      PROPERTY COMPILE_DEFINITIONS ${TCLCAD_TKOBOL_DEFINITIONS}
    )
  endif(TCLCAD_TKOBOL_HOST)
endif(BRLCAD_ENABLE_TK AND BRLCAD_ENABLE_OPENGL)

if(BRLCAD_ENABLE_TCL)
  # Include directories needed by libtclcad users
  set(TCLCAD_INCLUDE_DIRS ${TCL_INCLUDE_PATH})
  if(TK_INCLUDE_PATH)
    set(TCLCAD_INCLUDE_DIRS ${TCLCAD_INCLUDE_DIRS} ${TK_INCLUDE_PATH} ${X11_INCLUDE_DIR})
  endif(TK_INCLUDE_PATH)

  # Include directories only needed by the implementation
  set(TCLCAD_LOCAL_INCLUDE_DIRS ${PNG_INCLUDE_DIRS})
  if(BRLCAD_ENABLE_X11)
    set(TCLCAD_XLIBS ${TCLCAD_XLIBS} ${X11_LIBRARIES})
  endif(BRLCAD_ENABLE_X11)

  # Note - libtclcad_deps is defined by ${BRLCAD_SOURCE_DIR}/src/source_dirs.cmake
  set(TCLCAD_PUBLIC_LIBS ${libtclcad_deps} ${TCL_LIBRARY})
  set(TCLCAD_PRIVATE_LIBS ${TCLCAD_XLIBS} PNG::PNG)
  if(TCLCAD_TKOBOL_HOST)
    list(APPEND TCLCAD_PRIVATE_LIBS ${OPENGL_LIBRARIES})
  endif()
  if(TK_LIBRARY)
    set(TCLCAD_PUBLIC_LIBS ${TCLCAD_PUBLIC_LIBS} ${TK_LIBRARY})
  endif(TK_LIBRARY)

  brlcad_addlib(libtclcad "${LIBTCLCAD_SRCS}" "${TCLCAD_INCLUDE_DIRS}" "${TCLCAD_LOCAL_INCLUDE_DIRS}"
    PUBLIC_LIBS ${TCLCAD_PUBLIC_LIBS}
    PRIVATE_LIBS ${TCLCAD_PRIVATE_LIBS}
    # refresh.c defines DM_WITH_RT before including dm/view.h to expose
    # dm_draw_viewobjs; when batched after draw.c the DM_VIEW_H guard fires
    # without DM_WITH_RT set and the declaration is lost.
    # bu.c and dm.c both define a non-static register_cmds function.
    UNITY_BUILD_SKIP view/refresh.c bu.c dm.c
    # The full static library link is too large for some testing environments -
    # disable this test
    NO_STATIC_LINK_TEST
  )
  set(ptargets libtclcad)
  if(TARGET libtclcad-obj)
    list(APPEND ptargets libtclcad-obj)
  endif(TARGET libtclcad-obj)
  foreach(pt ${ptargets})
    target_compile_definitions(${pt} PRIVATE ITCL_VERSION=\"${ITCL_VERSION}\")
    target_compile_definitions(${pt} PRIVATE ITK_VERSION=\"${ITK_VERSION}\")
    target_compile_definitions(${pt} PRIVATE IWIDGETS_VERSION=\"${IWIDGETS_VERSION}\")
    if(TK_LIBRARY)
      target_compile_definitions(${pt} PRIVATE HAVE_TK)
    endif(TK_LIBRARY)
    if(TCLCAD_TKOBOL_HOST)
      target_compile_definitions(${pt} PRIVATE HAVE_TKOBOL_HOST)
    endif()
  endforeach(pt ${ptargets})

  set_target_properties(libtclcad PROPERTIES VERSION 20.0.1 SOVERSION 20)
  if(TARGET itcl_pkgIndex)
    add_dependencies(libtclcad itcl_pkgIndex)
  endif(TARGET itcl_pkgIndex)
  if(TARGET itk_pkgIndex)
    add_dependencies(libtclcad itk_pkgIndex)
  endif(TARGET itk_pkgIndex)

  # For anyone to do anything with tclcad, we need the libged plugins
  add_dependencies(libtclcad ged_plugins)
endif(BRLCAD_ENABLE_TCL)

add_subdirectory(tests)

cmakefiles(
  ${LIBTCLCAD_SRCS}
  draw_view_move_helpers.h
  tclcad_private.h
  view/view.h
  tkobol/tk-obol-host.cpp
  tkobol/vendor/togl/LICENSE
  tkobol/vendor/togl/README.brlcad
  tkobol/vendor/togl/Tcl3D_License.txt
  tkobol/vendor/togl/togl.c
  tkobol/vendor/togl/togl.h
  tkobol/vendor/togl/toglFont.c
  tkobol/vendor/togl/toglFont.h
  CMakeLists.txt
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
