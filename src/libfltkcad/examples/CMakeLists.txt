if(TARGET libfltkcad AND FLTKCAD_BUILD_EXAMPLES)
  brlcad_addexec(fltkcad_html_demo html_viewer.cpp libfltkcad NO_INSTALL)
  brlcad_addexec(fltkcad_terminal_demo terminal_demo.cpp libfltkcad NO_INSTALL)
  brlcad_addexec(fltkcad_display_demo display_demo.cpp libfltkcad NO_INSTALL)
endif()

cmakefiles(CMakeLists.txt html_viewer.cpp terminal_demo.cpp display_demo.cpp)
