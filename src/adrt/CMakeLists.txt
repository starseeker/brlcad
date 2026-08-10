find_package_zlib(REQUIRED)

set(_zlib_location "")
if(TARGET ZLIB::ZLIB)
  get_target_property(_zlib_location ZLIB::ZLIB IMPORTED_LOCATION)
  if(NOT _zlib_location)
    get_target_property(_zlib_location ZLIB::ZLIB IMPORTED_LOCATION_RELEASE)
  endif()
  if(NOT _zlib_location)
    get_target_property(_zlib_location ZLIB::ZLIB IMPORTED_LOCATION_DEBUG)
  endif()
endif()
is_subpath("${CMAKE_BINARY_DIR}" "${_zlib_location}" LOCAL_TEST)
if(LOCAL_TEST OR TARGET ZLIB_BLD)
  add_definitions(-DZ_PREFIX)
endif(LOCAL_TEST OR TARGET ZLIB_BLD)

set(
  LIBRENDER_SRCS
  load.c
  load_g.c
  librender/camera.c
  librender/component.c
  librender/cut.c
  librender/depth.c
  librender/flat.c
  librender/flos.c
  librender/grid.c
  librender/hit.c
  librender/normal.c
  librender/path.c
  librender/phong.c
  librender/render_util.c
  librender/spall.c
  librender/surfel.c
  librender/texture_blend.c
  librender/texture_bump.c
  librender/texture_camo.c
  librender/texture_checker.c
  librender/texture_clouds.c
  librender/texture_gradient.c
  librender/texture_image.c
  librender/texture_mix.c
  librender/texture_perlin.c
  librender/texture_stack.c
)

set(
  ADRT_LOCAL_INCLUDE_DIRS
  ${BRLCAD_SOURCE_DIR}/src/adrt
  ${BRLCAD_SOURCE_DIR}/src/adrt/libcommon
  ${BRLCAD_SOURCE_DIR}/src/adrt/librender
)

brlcad_addlib(librender "${LIBRENDER_SRCS}" "" "${ADRT_LOCAL_INCLUDE_DIRS}"
  PRIVATE_LIBS libgcv ZLIB::ZLIB ${M_LIBRARY}
)
set_target_properties(librender PROPERTIES VERSION 20.0.1 SOVERSION 20)
add_dependencies(librender gcv_plugins)

set(
  tie_HDRS
  adrt.h
  adrt_struct.h
  librender/camera.h
  librender/render.h
  librender/render_internal.h
  librender/render_util.h
  librender/texture.h
  librender/texture_internal.h
)
# TODO - should these be in brlcad/tie instead?
brlcad_manage_files(tie_HDRS ${INCLUDE_DIR}/tie)

# Tcl/Tk based ISST.  ADRT produces an RGB image; TkPhoto presents those
# pixels directly, so this tool does not need a graphical display manager or
# an OpenGL context.

if(BRLCAD_ENABLE_TK AND TK_LIBRARY)
  set(local_isst_dirs ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/librender)

  set(
    isst_libs
    libtclcad
    librender
    ${TCL_LIBRARY}
    ${TK_LIBRARY}
    ${M_LIBRARY}
    )

  brlcad_adddata(isst.tcl tclscripts/isst)
  brlcad_addexec(isst isst.c "${isst_libs}" GUI)
  target_include_directories(isst BEFORE PRIVATE ${local_isst_dirs})
  add_target_deps(isst isst_tcl_cp TclIndexBld)
  if(UNIX AND SH_EXEC)
    find_program(ISST_XVFB_RUN_EXECUTABLE xvfb-run)
    if(ISST_XVFB_RUN_EXECUTABLE)
      brlcad_add_test(NAME adrt_isst_tk_photo_smoke
	COMMAND ${SH_EXEC} ${BRLCAD_SOURCE_DIR}/misc/CMake/scripts/xvfb_test.sh
	${ISST_XVFB_RUN_EXECUTABLE} $<TARGET_FILE:isst>)
      set_property(TEST adrt_isst_tk_photo_smoke APPEND PROPERTY ENVIRONMENT
	"BRLCAD_ROOT=${CMAKE_BINARY_DIR};BRLCAD_ISST_SMOKE=1")
      set_property(TEST adrt_isst_tk_photo_smoke PROPERTY SKIP_RETURN_CODE 125)
    endif(ISST_XVFB_RUN_EXECUTABLE)
  endif(UNIX AND SH_EXEC)
endif(BRLCAD_ENABLE_TK AND TK_LIBRARY)

cmakefiles(isst.tcl)
cmakefiles(isst.c)

add_subdirectory(master)
add_subdirectory(slave)

set(adrt_ignore_files load.h tienet.h)
cmakefiles(${adrt_ignore_files})
cmakefiles(CMakeLists.txt)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
