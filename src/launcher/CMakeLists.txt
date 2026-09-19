set(
  LAUNCHER_SRCS
  launcher.c
  registry.c
  launch.c
  fbtext.c
  ui_fb.c
  ui_text.c
)

brlcad_addexec(brlcad-launcher "${LAUNCHER_SRCS}" "libBObol;libimgstream;libicv;libbu")
target_include_directories(brlcad-launcher BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if(TARGET libqtcad)
  get_property(LAUNCHER_QT_PROVIDER TARGET libqtcad
    PROPERTY QTCAD_OBOL_DISPLAY_PROVIDER_AVAILABLE)
endif()
if(LAUNCHER_QT_PROVIDER)
  target_link_libraries(brlcad-launcher libqtcad)
  target_compile_definitions(brlcad-launcher PRIVATE HAVE_QTCAD_OBOL_DISPLAY_PROVIDER)
elseif(TARGET libtclcad)
  get_property(LAUNCHER_TCL_PROVIDER TARGET libtclcad
    PROPERTY TCLCAD_OBOL_DISPLAY_PROVIDER_AVAILABLE)
  if(LAUNCHER_TCL_PROVIDER)
    target_link_libraries(brlcad-launcher libtclcad)
    target_compile_definitions(brlcad-launcher PRIVATE HAVE_TCLCAD_OBOL_DISPLAY_PROVIDER)
  endif()
endif()

# Application manifests, discovered at run time by the launcher (nothing is
# hard-coded).  New GUI programs appear in the launcher simply by installing an
# additional manifest here.
brlcad_adddata(apps/mged.app launcher/apps)
brlcad_adddata(apps/archer.app launcher/apps)
brlcad_adddata(apps/rtwizard.app launcher/apps)
brlcad_adddata(apps/isst.app launcher/apps)
brlcad_adddata(apps/shell.app launcher/apps)

# Rendered splash image (optional).  If present it is installed and shown by the
# windowed UI; otherwise the launcher falls back to the bundled BRL-CAD logo.
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/assets/splash.png")
  brlcad_adddata(assets/splash.png launcher)
endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/assets/splash.png")

cmakefiles(
  CMakeLists.txt
  README
  launcher.h
  fbtext.h
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
