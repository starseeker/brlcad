brlcad_addexec(fbserv "fbserv.c" "libimgstream;libpkg;libbu")
target_link_libraries(fbserv ${WINSOCK_LIB})
set(FBSERV_OBOL_DISPLAY_PROVIDER FALSE)
set(FBSERV_OBOL_SYSTEM_GL_PROVIDER FALSE)
# Prefer Qt when it is enabled.  TclCAD supplies the same small provider API
# when Qt is unavailable; neither branch changes fbserv's session code.
if(TARGET libqtcad)
  get_property(FBSERV_QTCAD_OBOL_PROVIDER TARGET libqtcad
    PROPERTY QTCAD_OBOL_DISPLAY_PROVIDER_AVAILABLE)
  if(FBSERV_QTCAD_OBOL_PROVIDER)
    target_link_libraries(fbserv libqtcad)
    target_compile_definitions(fbserv PRIVATE HAVE_QTCAD_OBOL_DISPLAY_PROVIDER)
    set(FBSERV_OBOL_DISPLAY_PROVIDER TRUE)
    get_property(FBSERV_QTCAD_OBOL_SYSTEM_GL_PROVIDER TARGET libqtcad
      PROPERTY QTCAD_OBOL_SYSTEM_GL_PROVIDER_AVAILABLE)
    if(FBSERV_QTCAD_OBOL_SYSTEM_GL_PROVIDER)
      set(FBSERV_OBOL_SYSTEM_GL_PROVIDER TRUE)
    endif(FBSERV_QTCAD_OBOL_SYSTEM_GL_PROVIDER)
  endif(FBSERV_QTCAD_OBOL_PROVIDER)
endif(TARGET libqtcad)
if(NOT FBSERV_OBOL_DISPLAY_PROVIDER AND TARGET libtclcad)
  get_property(FBSERV_TCLCAD_OBOL_PROVIDER TARGET libtclcad
    PROPERTY TCLCAD_OBOL_DISPLAY_PROVIDER_AVAILABLE)
  if(FBSERV_TCLCAD_OBOL_PROVIDER)
    target_link_libraries(fbserv libtclcad)
    target_compile_definitions(fbserv PRIVATE HAVE_TCLCAD_OBOL_DISPLAY_PROVIDER)
    set(FBSERV_OBOL_DISPLAY_PROVIDER TRUE)
    set(FBSERV_OBOL_SYSTEM_GL_PROVIDER TRUE)
  endif(FBSERV_TCLCAD_OBOL_PROVIDER)
endif(NOT FBSERV_OBOL_DISPLAY_PROVIDER AND TARGET libtclcad)
if(FBSERV_OBOL_DISPLAY_PROVIDER AND UNIX)
  find_program(XVFB_RUN_EXECUTABLE xvfb-run)
  if(XVFB_RUN_EXECUTABLE)
    set(FBSERV_XVFB_TEST_WRAPPER
      ${BRLCAD_SOURCE_DIR}/misc/CMake/scripts/xvfb_test.sh)
    set(FBSERV_DISPLAY_TARGETS /dev/swrast)
    if(FBSERV_OBOL_SYSTEM_GL_PROVIDER)
      list(APPEND FBSERV_DISPLAY_TARGETS /dev/ogl)
    endif(FBSERV_OBOL_SYSTEM_GL_PROVIDER)
    foreach(FBSERV_DISPLAY_TARGET ${FBSERV_DISPLAY_TARGETS})
      string(REPLACE "/dev/" "" FBSERV_DISPLAY_NAME ${FBSERV_DISPLAY_TARGET})
      add_test(
        NAME fbserv_display_${FBSERV_DISPLAY_NAME}
        COMMAND ${SH_EXEC} ${FBSERV_XVFB_TEST_WRAPPER} ${XVFB_RUN_EXECUTABLE}
                ${SH_EXEC} ${CMAKE_CURRENT_SOURCE_DIR}/fbserv_display_smoke.sh
                $<TARGET_FILE:fbserv> ${FBSERV_DISPLAY_TARGET}
      )
      set_property(TEST fbserv_display_${FBSERV_DISPLAY_NAME}
        APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
      set_property(TEST fbserv_display_${FBSERV_DISPLAY_NAME}
        PROPERTY SKIP_RETURN_CODE 125)
    endforeach(FBSERV_DISPLAY_TARGET)
  endif(XVFB_RUN_EXECUTABLE)
endif(FBSERV_OBOL_DISPLAY_PROVIDER AND UNIX)

cmakefiles(
  tls_wrap.h
  fbserv_display_smoke.sh
  CMakeLists.txt
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
