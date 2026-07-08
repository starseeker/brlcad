#  ************ libdm tests *************

brlcad_addexec(dm_test dm_test.c "libdm;libbu" TEST)
brlcad_addexec(test_dm_fb_imgstream fb_imgstream.cpp "libdm;libbrlobol;libimgstream;libbu" TEST ALL_CXX)
brlcad_add_test(NAME libdm_fb_imgstream COMMAND test_dm_fb_imgstream)

#TODO - these should be portable without X11, but we need to set up the Tk Xlib
# and provide an appropriate include first...
if(BRLCAD_ENABLE_TK AND BRLCAD_ENABLE_X11)
  set(TCL_ENABLE_TK ON)
  find_package_tcl()

  if(TCL_FOUND AND TK_FOUND)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${TCL_INCLUDE_PATH} ${TK_INCLUDE_PATH} ${X11_INCLUDE_DIR})

    if(BUILD_TESTING)
      brlcad_addexec(tcl_img tcl_img.cpp "libdm;libbu;${TCL_LIBRARY};${TK_LIBRARY}" TEST_USEDATA)
      brlcad_addexec(tcl_tevent tcl_tevent.cpp "libdm;libbu;${TCL_LIBRARY};${TK_LIBRARY}" TEST_USEDATA)
    endif(BUILD_TESTING)
  endif(TCL_FOUND AND TK_FOUND)
endif()

cmakefiles(CMakeLists.txt)
cmakefiles(
  dm_test.c
  fb_imgstream.cpp
  tcl_img.cpp
  tcl_tevent.cpp
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
