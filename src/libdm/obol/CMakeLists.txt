set(OBOL_SRCS dm-obol.cpp)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${BRLCAD_BINARY_DIR}/include
  ${BRLCAD_SOURCE_DIR}/include
  ${BU_INCLUDE_DIRS}
)

add_definitions(-DDM_PLUGIN)

dm_plugin_library(dm-obol SHARED ${OBOL_SRCS})
target_link_libraries(dm-obol libbrlobol libdm libbu)
set_property(TARGET dm-obol APPEND PROPERTY COMPILE_DEFINITIONS BRLCADBUILD HAVE_CONFIG_H)
validate_style(dm-obol "${OBOL_SRCS}")

plugin_setup(dm-obol dm)

cmakefiles(
  CMakeLists.txt
  ${OBOL_SRCS}
  dm-obol.h
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
