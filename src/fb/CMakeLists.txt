brlcad_find_package(PNG REQUIRED)

if(MSVC)
  add_definitions(-DIF_WGL)
endif(MSVC)
brlcad_addexec(bw-fb bw-fb.c "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(fb-bw fb-bw.c "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(fb-fb fb-fb.c "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(fb-pix "fb-pix.c;cmap-crunch.c" "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(fb-png "fb-png.c;cmap-crunch.c" "libimgstream;libbu;PNG::PNG" FOLDER Fb)
target_include_directories(fb-png BEFORE PRIVATE ${PNG_INCLUDE_DIRS})
brlcad_addexec(fbcbars fbcbars.c "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(fbclear fbclear.c "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(fbcmap fbcmap.c "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(fbcolor fbcolor.c "libimgstream;libbu" FOLDER Fb)
target_include_directories(fbcolor BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../libtermio)
brlcad_add_test(
  NAME fbcolor_full_ramp
  COMMAND
    "${CMAKE_COMMAND}" "-DFBCOLOR=$<TARGET_FILE:fbcolor>" "-DTEST_DIR=${CMAKE_CURRENT_BINARY_DIR}"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/test_fbcolor_full_ramp.cmake"
)
brlcad_addexec(fbfade "fbfade.c;ioutil.c" "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(fbframe fbframe.c "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(fbfree fbfree.c "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(fbgamma fbgamma.c "libimgstream;libbu;${M_LIBRARY}" FOLDER Fb)
brlcad_addexec(fbgrid fbgrid.c "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(fbhelp fbhelp.c "libimgstream;libbu" FOLDER Fb)
brlcad_add_test(
  NAME fbhelp_output_stream
  COMMAND "${CMAKE_COMMAND}" "-DFBHELP=$<TARGET_FILE:fbhelp>" -P "${CMAKE_CURRENT_SOURCE_DIR}/test_fbhelp_output.cmake"
)
brlcad_addexec(fblabel fblabel.c "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(fbline fbline.c "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(fbpoint fbpoint.c "libimgstream;libbu" FOLDER Fb)
target_include_directories(fbpoint BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../libtermio)
brlcad_addexec(fbstretch "fbstretch.c;ioutil.c" "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(fbzoom fbzoom.c "libimgstream;libbu" FOLDER Fb)
target_include_directories(fbzoom BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../libtermio)
brlcad_addexec(gif-fb "gif-fb.c;ioutil.c" "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(gif2fb gif2fb.c "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(pix-fb pix-fb.c "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(pixautosize pixautosize.c "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(plot3-fb plot3-fb.c "libimgstream;libbu" FOLDER Fb)
brlcad_addexec(png-fb png-fb.c "libimgstream;libbu;PNG::PNG" FOLDER Fb)
target_include_directories(png-fb BEFORE PRIVATE ${PNG_INCLUDE_DIRS})

brlcad_addexec(spm-fb spm-fb.c "libimgstream;libbn;libbu" FOLDER Fb)

cmakefiles(
  CMakeLists.txt
  fbcolor.c
  fbpoint.c
  fbzoom.c
  test_fbcolor_full_ramp.cmake
  test_fbhelp_output.cmake
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
