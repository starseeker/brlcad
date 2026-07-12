# Main libdm library
set(
  LIBDM_SRCS
  null/dm-Null.c
  null/if_null.c
  txt/if_debug.c
  adc.c
  asize.c
  axes.c
  clip.c
  dm-generic.c
  dm_init.cpp
  dm_plugins.cpp
  fbserv_legacy.c
  dm_util.c
  fb_generic.c
  fb_log.c
  fb_paged_io.c
  fb_rect.c
  fb_util.c
  grid.c
  if_disk.c
  if_mem.c
  if_remote.c
  if_stack.c
  labels.c
  options.c
  rect.c
  scale.c
  view.c
  vers.c
)
set_property(
  SOURCE dm_init.cpp
  APPEND
  PROPERTY COMPILE_DEFINITIONS "DM_PLUGIN_SUFFIX=\"${CMAKE_SHARED_LIBRARY_SUFFIX}\""
)

brlcad_find_package(PNG REQUIRED)

if(BRLCAD_ENABLE_OPENGL)
  set(DM_INCLUDE_DIRS ${DM_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIR})
endif(BRLCAD_ENABLE_OPENGL)

# local includes
set(DM_LOCAL_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${PNG_INCLUDE_DIRS} ${BRLCAD_SOURCE_DIR}/src/libbg)

# Note - libdm_deps is defined by ${BRLCAD_SOURCE_DIR}/src/source_dirs.cmake
brlcad_addlib(libdm "${LIBDM_SRCS}" "${DM_INCLUDE_DIRS}" "${DM_LOCAL_INCLUDE_DIRS}"
  PUBLIC_LIBS ${libdm_deps}
  PRIVATE_LIBS libimgstream PNG::PNG ${WINSOCK_LIB}
  UNITY_BUILD_SKIP dm_init.cpp
  )
set_target_properties(libdm PROPERTIES VERSION 20.0.1 SOVERSION 20)

if(TARGET profont_ProFont_ttf_cp)
  add_dependencies(libdm profont_ProFont_ttf_cp)
endif(TARGET profont_ProFont_ttf_cp)

add_custom_command(
  OUTPUT "${CMAKE_BINARY_DIR}/dm_plugins.sentinel"
  COMMAND ${CMAKE_COMMAND} -E touch "${CMAKE_BINARY_DIR}/dm_plugins.sentinel"
)
add_custom_target(dm_plugins_sentinel DEPENDS ${CMAKE_BINARY_DIR}/dm_plugins.sentinel)
set_target_properties(dm_plugins_sentinel PROPERTIES FOLDER "BRL-CAD Plugins/dm")

function(dm_plugin_library name)
  add_library(${name} ${ARGN})
  add_dependencies(dm_plugins_sentinel ${name})
endfunction(dm_plugin_library name)

# This target is supplied so applications wanting to run libdm
# commands at build time can depend on the plugins being built as well as
# libdm itself (otherwise LIBGED commands might fail due to their plugins not
# yet having been built.)
add_custom_target(dm_plugins ALL DEPENDS dm_plugins_sentinel)
set_target_properties(dm_plugins PROPERTIES FOLDER "BRL-CAD Plugins")

# Retired plugin targets can survive indefinitely in an existing build tree,
# where libdm will still discover and attempt to load them.  Remove those
# orphaned outputs when configuring the consolidated Obol backend.
set(DM_RETIRED_PLUGINS dm-X dm-obol dm-ogl dm-plot dm-ps dm-qtgl dm-swrast dm-tkswrast dm-tkobol dm-txt)
foreach(retired_plugin IN LISTS DM_RETIRED_PLUGINS)
  file(
    REMOVE
    "${CMAKE_BINARY_DIR}/${LIBEXEC_DIR}/dm/${CMAKE_SHARED_LIBRARY_PREFIX}${retired_plugin}${CMAKE_SHARED_LIBRARY_SUFFIX}"
  )
endforeach()

add_subdirectory(tests)

set(
  libdm_ignore_files
  CMakeLists.txt
  README
  TODO
  if_TEMPLATE.c
  null/dm-Null.h
  include/private.h
  include/calltable.h
  fontstash/LICENSE.fontstash
  fontstash/README.md
  fontstash/fontstash.h
  fontstash/glfontstash.h
  fontstash/stb_truetype.h
)
cmakefiles(${libdm_ignore_files})

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
