set(FLTKGED_SRCS
  Batch.cpp
  Fl_Ged_Localization.cpp
  FltkGedWindow.cpp
  main.cpp
)

if(BRLCAD_ENABLE_FLTK)
  find_package(Intl QUIET)
  brlcad_addexec(fltkged "${FLTKGED_SRCS}" "libfltkcad" NO_MAN)
  target_compile_features(fltkged PRIVATE cxx_std_17)
  target_include_directories(fltkged PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  if(Intl_FOUND)
    target_compile_definitions(fltkged PRIVATE
      FLTKGED_HAVE_GETTEXT=1
      FLTKGED_LOCALE_DIRECTORY="${CMAKE_INSTALL_PREFIX}/${DATA_DIR}/locale")
    target_link_libraries(fltkged Intl::Intl)
  endif()
  add_dependencies(fltkged dm_plugins ged_plugins)
  if(BRLCAD_ENABLE_TCL)
    add_dependencies(fltkged TclIndexBld)
  endif()

  install(FILES README.md COMPATIBILITY.md REGRESSION_STATUS.md RELEASE_CHECKLIST.md
    DESTINATION ${DOC_DIR}/fltkged)
  install(FILES po/fltkged.pot po/README.md DESTINATION ${DATA_DIR}/fltkged/po)
endif()

add_subdirectory(tests)

cmakefiles(CMakeLists.txt ${FLTKGED_SRCS} Batch.H FltkGedWindow.H Fl_Ged_Localization.H
  README.md COMPATIBILITY.md REGRESSION_STATUS.md RELEASE_CHECKLIST.md po/fltkged.pot po/README.md)
