include(${CMAKE_CURRENT_SOURCE_DIR}/api_tiers.cmake)
brlcad_manage_files(LIBBOBOL_INSTALLED_HEADERS
  ${INCLUDE_DIR}/brlcad/BObol REQUIRED libBObol)

cmakefiles(CMakeLists.txt api_tiers.cmake)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
