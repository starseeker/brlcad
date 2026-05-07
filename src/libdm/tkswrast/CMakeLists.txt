brlcad_find_package(OSMESA REQUIRED)

if(BRLCAD_ENABLE_TK)
  set(TCL_ENABLE_TK ON)
  find_package_tcl(TCL)

  if(TCL_FOUND AND TK_FOUND)
    include_directories(
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_SOURCE_DIR}/..
      ${BRLCAD_BINARY_DIR}/include
      ${BRLCAD_SOURCE_DIR}/include
      ${BU_INCLUDE_DIRS}
      ${OSMESA_INCLUDE_DIR}
      ${TCL_INCLUDE_PATH}
      ${TK_INCLUDE_PATH}
    )

    add_definitions(-DDM_PLUGIN)
    add_definitions(-DUSE_MGL_NAMESPACE)

    set(TKSWRAST_SRCS dm-tkswrast.cpp)

    dm_plugin_library(dm-tkswrast SHARED ${TKSWRAST_SRCS})
    target_link_libraries(dm-tkswrast libdm libbu ${OSMESA_LIBRARY} ${TCL_LIBRARY} ${TK_LIBRARY})
    set_property(TARGET dm-tkswrast APPEND PROPERTY COMPILE_DEFINITIONS BRLCADBUILD HAVE_CONFIG_H HAVE_TK OSMESA)
    validate_style(dm-tkswrast "${TKSWRAST_SRCS}")

    plugin_setup(dm-tkswrast dm)
  endif(TCL_FOUND AND TK_FOUND)
endif(BRLCAD_ENABLE_TK)

cmakefiles(
  CMakeLists.txt
  dm-tkswrast.cpp
)
