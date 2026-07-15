if(BRLCAD_ENABLE_OPENCL)
  add_definitions(-DUSE_OPENCL)
endif(BRLCAD_ENABLE_OPENCL)

set(
  LIBRTUIF_SOURCES
  do.c
  grid.c
  heatgraph.c
  opt.c
  progressive.c
  scanline.c
  usage.cpp
  worker.c
)

brlcad_addlib(librtuif "${LIBRTUIF_SOURCES}" "" "${CMAKE_CURRENT_SOURCE_DIR}"
  STATIC NO_INSTALL NO_STATIC_LINK_TEST
  PRIVATE_LIBS librt libimgstream ${M_LIBRARY}
)

list(APPEND LIBRTUIF_SOURCES main.c)

if(BRLCAD_ENABLE_OPENCL)
  set(OPENCL_LIBS ${OPENCL_LIBRARIES})
endif(BRLCAD_ENABLE_OPENCL)

set(RT_STD_LIBS Threads::Threads)
if(MPI_C_LIBRARIES AND BRLCAD_ENABLE_MPI)
  set(RT_STD_LIBS ${RT_STD_LIBS} ${MPI_C_LIBRARIES})
endif(MPI_C_LIBRARIES AND BRLCAD_ENABLE_MPI)
set(
  RT_STD_LIBS
  ${RT_STD_LIBS}
  librt
  libnmg
  libimgstream
  liboptical
  libicv
  ${M_LIBRARY}
)

brlcad_adddata(rtexample.c sample_applications)

brlcad_addexec(rtexample "rtexample.c" "librt;liboptical;libicv" NO_INSTALL)

brlcad_addexec(reshoot "reshoot.c" "Threads::Threads;librt;liboptical;libicv" FOLDER Rt)
target_include_directories(reshoot BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

brlcad_addexec(rtshot "rtshot.c" "${RT_STD_LIBS}" FOLDER Rt)
target_include_directories(rtshot BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

brlcad_addexec(rtscale "read-rtlog.c;rtscale.c" "Threads::Threads;librt;liboptical;libicv;${M_LIBRARY}" FOLDER Rt)
target_include_directories(rtscale BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

brlcad_addexec(rt "${LIBRTUIF_SOURCES};view.c" "${OPENCL_LIBS};${RT_STD_LIBS}" FOLDER Rt)
target_include_directories(rt BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# /dev/null is the no-presentation raytrace benchmark sink.  It must remain
# usable without any toolkit display or registered Obol provider.
if(UNIX)
  add_test(
    NAME rt_null_sink
    COMMAND ${SH_EXEC} ${CMAKE_CURRENT_SOURCE_DIR}/rt_null_sink_smoke.sh
            $<TARGET_FILE:rt> ${CMAKE_BINARY_DIR}/share/db/moss.g
  )
  set_property(TEST rt_null_sink APPEND PROPERTY ENVIRONMENT
    BRLCAD_ROOT=${CMAKE_BINARY_DIR})
endif(UNIX)
set(RT_OBOL_DISPLAY_PROVIDER FALSE)
set(RT_OBOL_SYSTEM_GL_PROVIDER FALSE)
# Prefer Qt when it is enabled.  TclCAD supplies the same small provider API
# when Qt is unavailable; neither branch changes rt's session code.
if(TARGET libqtcad)
  get_property(RT_QTCAD_OBOL_PROVIDER TARGET libqtcad
    PROPERTY QTCAD_OBOL_DISPLAY_PROVIDER_AVAILABLE)
  if(RT_QTCAD_OBOL_PROVIDER)
    target_link_libraries(rt libqtcad)
    target_compile_definitions(rt PRIVATE HAVE_QTCAD_OBOL_DISPLAY_PROVIDER)
    set(RT_OBOL_DISPLAY_PROVIDER TRUE)
    get_property(RT_QTCAD_OBOL_SYSTEM_GL_PROVIDER TARGET libqtcad
      PROPERTY QTCAD_OBOL_SYSTEM_GL_PROVIDER_AVAILABLE)
    if(RT_QTCAD_OBOL_SYSTEM_GL_PROVIDER)
      set(RT_OBOL_SYSTEM_GL_PROVIDER TRUE)
    endif(RT_QTCAD_OBOL_SYSTEM_GL_PROVIDER)
  endif(RT_QTCAD_OBOL_PROVIDER)
endif(TARGET libqtcad)
if(NOT RT_OBOL_DISPLAY_PROVIDER AND TARGET libtclcad)
  get_property(RT_TCLCAD_OBOL_PROVIDER TARGET libtclcad
    PROPERTY TCLCAD_OBOL_DISPLAY_PROVIDER_AVAILABLE)
  if(RT_TCLCAD_OBOL_PROVIDER)
    target_link_libraries(rt libtclcad)
    target_compile_definitions(rt PRIVATE HAVE_TCLCAD_OBOL_DISPLAY_PROVIDER)
    set(RT_OBOL_DISPLAY_PROVIDER TRUE)
    set(RT_OBOL_SYSTEM_GL_PROVIDER TRUE)
  endif(RT_TCLCAD_OBOL_PROVIDER)
endif(NOT RT_OBOL_DISPLAY_PROVIDER AND TARGET libtclcad)
if(RT_OBOL_DISPLAY_PROVIDER AND UNIX)
  find_program(XVFB_RUN_EXECUTABLE xvfb-run)
  if(XVFB_RUN_EXECUTABLE)
    set(RT_XVFB_TEST_WRAPPER
      ${BRLCAD_SOURCE_DIR}/misc/CMake/scripts/xvfb_test.sh)
    set(RT_DISPLAY_TARGETS /dev/swrast)
    if(RT_OBOL_SYSTEM_GL_PROVIDER)
      list(APPEND RT_DISPLAY_TARGETS /dev/ogl)
    endif(RT_OBOL_SYSTEM_GL_PROVIDER)
    foreach(RT_DISPLAY_TARGET ${RT_DISPLAY_TARGETS})
      string(REPLACE "/dev/" "" RT_DISPLAY_NAME ${RT_DISPLAY_TARGET})
      add_test(
        NAME rt_display_${RT_DISPLAY_NAME}
        COMMAND ${SH_EXEC} ${RT_XVFB_TEST_WRAPPER} ${XVFB_RUN_EXECUTABLE}
                ${SH_EXEC} ${CMAKE_CURRENT_SOURCE_DIR}/rt_display_smoke.sh
                $<TARGET_FILE:rt> ${CMAKE_BINARY_DIR}/share/db/moss.g
                ${RT_DISPLAY_TARGET}
      )
      set_property(TEST rt_display_${RT_DISPLAY_NAME}
        APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
      set_property(TEST rt_display_${RT_DISPLAY_NAME}
        PROPERTY SKIP_RETURN_CODE 125)
    endforeach(RT_DISPLAY_TARGET)
  endif(XVFB_RUN_EXECUTABLE)
endif(RT_OBOL_DISPLAY_PROVIDER AND UNIX)

brlcad_addexec(rtarea "${LIBRTUIF_SOURCES};viewarea.c" "${RT_STD_LIBS}" FOLDER Rt)
target_include_directories(rtarea BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TARGET rtarea APPEND PROPERTY COMPILE_DEFINITIONS "RT_TXT_OUTPUT")

brlcad_addexec(rtcheck "${LIBRTUIF_SOURCES};viewcheck.c" "${RT_STD_LIBS}" FOLDER Rt)
target_include_directories(rtcheck BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TARGET rtcheck APPEND PROPERTY COMPILE_DEFINITIONS "RT_TXT_OUTPUT")

brlcad_addexec(rtedge "${LIBRTUIF_SOURCES};viewedge.c" "${RT_STD_LIBS}" FOLDER Rt)
target_include_directories(rtedge BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

brlcad_addexec(rtg3 "${LIBRTUIF_SOURCES};viewg3.c" "${RT_STD_LIBS}" FOLDER Rt)
target_include_directories(rtg3 BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TARGET rtg3 APPEND PROPERTY COMPILE_DEFINITIONS "RT_TXT_OUTPUT")

brlcad_addexec(rthide "${LIBRTUIF_SOURCES};viewhide.c" "${RT_STD_LIBS}" FOLDER Rt)
target_include_directories(rthide BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TARGET rthide APPEND PROPERTY COMPILE_DEFINITIONS "RT_TXT_OUTPUT")

brlcad_addexec(rtsil "${LIBRTUIF_SOURCES};viewsil.c" "${RT_STD_LIBS}" FOLDER Rt)
target_include_directories(rtsil BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

brlcad_addexec(rtsurf "rtsurf.c;rtsurf_hits.cpp" "libanalyze;librt;liboptical;libicv")
target_include_directories(rtsurf BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

brlcad_addexec(rtweight "${LIBRTUIF_SOURCES};viewweight.c" "libanalyze;${RT_STD_LIBS}" FOLDER Rt)
target_include_directories(rtweight BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TARGET rtweight APPEND PROPERTY COMPILE_DEFINITIONS "RT_TXT_OUTPUT")

brlcad_addexec(rtxray "${LIBRTUIF_SOURCES};viewxray.c" "${RT_STD_LIBS}" FOLDER Rt)
target_include_directories(rtxray BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

brlcad_addexec(rtdepth "${LIBRTUIF_SOURCES};viewdepth.c" "${RT_STD_LIBS}" FOLDER Rt)
target_include_directories(rtdepth BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

set(
  rt_ignore_files
  CMakeLists.txt
  ext.h
  main.c
  rt_display_smoke.sh
  rt_null_sink_smoke.sh
  rtuif.h
  rtsurf_hits.h
  scanline.h
  viewdir.c
  viewdummy.c
)
cmakefiles(${rt_ignore_files})

# -------------------------------------------------------------------
# Regression tests for opt.c / get_args()
#
# We compile all LIBRTUIF_SOURCES except main.c, plus viewcheck.c as
# a minimal view stub, plus the test driver.  This avoids duplicating
# main() while still exercising get_args() directly so that gcov/lcov
# can measure branch coverage of opt.c.
# -------------------------------------------------------------------

set(
  RT_OPT_TEST_SRCS
  do.c
  grid.c
  heatgraph.c
  opt.c
  progressive.c
  scanline.c
  usage.cpp
  worker.c
  viewcheck.c
  tests/test_rt_opt.c
)

brlcad_addexec(test_rt_opt "${RT_OPT_TEST_SRCS}" "${RT_STD_LIBS}" TEST
  UNITY_BUILD_SKIP tests/test_rt_opt.c
)
target_include_directories(test_rt_opt BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set_property(
  TARGET test_rt_opt
  APPEND PROPERTY COMPILE_DEFINITIONS "RT_TXT_OUTPUT"
)

add_subdirectory(tests)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
