if(BRLCAD_ENABLE_TK)
  set(TCL_ENABLE_TK ON)
endif(BRLCAD_ENABLE_TK)

find_package_tcl()

set(
  MGED_SOURCES
  adc.c
  arbs.c
  attach.c
  axes.c
  buttons.c
  chgmodel.c
  chgtree.c
  chgview.c
  clone.c
  cmd.cpp
  color_scheme.c
  columns.c
  dm-generic.c
  dodraw.c
  doevent.c
  scene_refresh.c
  edsol.c
  f_db.c
  facedef.c
  fbserv.c
  grid.c
  history.c
  hud.c
  mater.c
  menu.c
  mged.c
  mged_impl.cpp
  mover.c
  overlay.c
  plot.c
  rect.c
  rtif.c
  scroll.c
  set.c
  setup.c
  share.c
  tclsync.c
  tedit.c
  titles.c
  track.c
  update.c
  usepen.c
  utility1.c
  vparse.c
  wdb_obj.c
)

cmakefiles(hud.h)

if(HAVE_WINDOWS_H)
  # To associate an icon with the application for Windows (needed
  # for the WiX installer) we must use an rc file.
  enable_language(RC)
  set(MGED_SOURCES ${MGED_SOURCES} mged.rc)
endif(HAVE_WINDOWS_H)

set(
  mged_libs
  libtclcad
  libged
  libwdb
  libdm
  libimgstream
  libnmg
  Threads::Threads
  ${WS2_32_LIBRARY}
  ${OPENGL32_LIBRARY}
)

# Ideally we'd build as GUI for Windows to avoid the command window popping up,
# but I'm not sure if that would be a problem for scripting with mged.exe on
# Windows so leave it for now...
if(BRLCAD_ENABLE_TCL)
  set(MGED_INCLUDE_DIRS
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../libtermio
  )
  brlcad_addexec(mged "${MGED_SOURCES}" "${mged_libs}")
  if(TK_LIBRARY)
    target_compile_definitions(mged PRIVATE HAVE_TK=1)
  endif(TK_LIBRARY)
  target_include_directories(mged BEFORE PRIVATE ${MGED_INCLUDE_DIRS})
  add_target_deps(mged Tkhtml Tktable itcl_pkgIndex itk_pkgIndex TclIndexBld)
  if(TARGET rt)
    add_dependencies(mged rt)
  endif(TARGET rt)

  # MGED needs the tclIndex and pkgIndex files set up correctly
	  add_dependencies(mged TclIndexBld)

	  if(SH_EXEC)
	    brlcad_add_test(NAME mged_obol_smoke
	      COMMAND ${SH_EXEC} "${CMAKE_CURRENT_SOURCE_DIR}/mged_obol_smoke.sh"
	      "$<TARGET_FILE:mged>" "${CMAKE_SOURCE_DIR}/src/libged/tests/draw/moss.g"
	      "${CMAKE_CURRENT_BINARY_DIR}")
	    set_property(TEST mged_obol_smoke APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})

	    if(UNIX)
	      find_program(MGED_XVFB_RUN_EXECUTABLE xvfb-run)
	      if(MGED_XVFB_RUN_EXECUTABLE)
		brlcad_add_test(NAME mged_tkobol_attach_smoke
		  COMMAND ${MGED_XVFB_RUN_EXECUTABLE} -a ${SH_EXEC}
		  "${CMAKE_CURRENT_SOURCE_DIR}/mged_tkobol_attach_smoke.sh"
		  "$<TARGET_FILE:mged>" "${CMAKE_SOURCE_DIR}/src/libged/tests/draw/moss.g"
		  "${CMAKE_CURRENT_BINARY_DIR}")
		set_property(TEST mged_tkobol_attach_smoke APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
	      endif(MGED_XVFB_RUN_EXECUTABLE)
	    endif(UNIX)

	    brlcad_add_test(NAME mged_obol_ert_smoke
	      COMMAND ${SH_EXEC} "${CMAKE_CURRENT_SOURCE_DIR}/mged_obol_ert_smoke.sh"
	      "$<TARGET_FILE:mged>" "${CMAKE_SOURCE_DIR}/src/libged/tests/draw/moss.g"
	      "${CMAKE_CURRENT_BINARY_DIR}")
	    set_tests_properties(mged_obol_ert_smoke PROPERTIES TIMEOUT 90)
	    set_property(TEST mged_obol_ert_smoke APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})

	    brlcad_add_test(NAME mged_obol_progressive_lod
	      COMMAND ${SH_EXEC} "${CMAKE_CURRENT_SOURCE_DIR}/mged_obol_progressive_lod.sh"
	      "$<TARGET_FILE:mged>" "${CMAKE_SOURCE_DIR}/src/libged/tests/draw/moss.g"
	      "${CMAKE_CURRENT_BINARY_DIR}")
	    set_tests_properties(mged_obol_progressive_lod PROPERTIES TIMEOUT 120)
	    set_property(TEST mged_obol_progressive_lod APPEND PROPERTY ENVIRONMENT BRLCAD_ROOT=${CMAKE_BINARY_DIR})
	  endif(SH_EXEC)

	  # MGED is one of the programs that gets a start menu entry
	  set_property(
	    INSTALL "${BIN_DIR}/$<TARGET_FILE_NAME:mged>"
    PROPERTY CPACK_START_MENU_SHORTCUTS "MGED ${BRLCAD_VERSION}"
  )
  set_property(INSTALL "${BIN_DIR}/$<TARGET_FILE_NAME:mged>" PROPERTY CPACK_DESKTOP_SHORTCUTS "MGED ${BRLCAD_VERSION}")
else(BRLCAD_ENABLE_TCL)
  cmakefiles(${MGED_SOURCES})
endif(BRLCAD_ENABLE_TCL)

brlcad_addexec(cad_boundp cad_boundp.c "libbu;${M_LIBRARY}")

brlcad_addexec(cad_parea cad_parea.c libbu)

set(
  mged_ignore_files
  CMakeLists.txt
  cmd.h
  f_cmd.h
  menu.h
  mged_obol_ert_smoke.sh
  mged_obol_progressive_lod.sh
  mged_obol_smoke.sh
  mged_tkobol_attach_smoke.sh
  mged.h
  mged.ico
  mged.rc
  mged_dm.h
  mged_impl.h
  mged_wdb.h
  sedit.h
)
cmakefiles(${mged_ignore_files})

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
