# Copyright (c) 2010-2024 United States Government as represented by
# the U.S. Army Research Laboratory.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above
# copyright notice, this list of conditions and the following
# disclaimer in the documentation and/or other materials provided
# with the distribution.
#
# 3. The name of the author may not be used to endorse or promote
# products derived from this software without specific prior written
# permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
# OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

include(CMakeParseArguments)

# Build a local copy of bext if we were unable to locate one

function(brlcad_ext_setup)

  set(BRLCAD_EXT_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/bext_build)
  set(BRLCAD_EXT_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR})

  # If we don't have
  if (NOT DEFINED BRLCAD_EXT_SOURCE_DIR)
    set(BRLCAD_EXT_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/bext)
    if (NOT EXISTS ${BRLCAD_EXT_SOURCE_DIR})
      find_program(GIT_EXEC git REQUIRED)
      execute_process(COMMAND ${GIT_EXEC} clone https://github.com/BRL-CAD/bext.git)
    endif (NOT EXISTS ${BRLCAD_EXT_SOURCE_DIR})
  endif (NOT DEFINED BRLCAD_EXT_SOURCE_DIR)
  if (NOT EXISTS ${BRLCAD_EXT_SOURCE_DIR})
    message(FATAL_ERROR "bext directory ${BRLCAD_EXT_SOURCE_DIR} is not present")
  endif (NOT EXISTS ${BRLCAD_EXT_SOURCE_DIR})

  if (NOT EXISTS ${BRLCAD_EXT_BUILD_DIR})
    file(MAKE_DIRECTORY ${BRLCAD_EXT_BUILD_DIR})
  endif (NOT EXISTS ${BRLCAD_EXT_BUILD_DIR})

  # Need to control options for this based on BRL-CAD configure settings.
  # Unlike an independent bext build, we know for this one what we can turn on
  # and off
  set(BEXT_ENABLE_ALL OFF)
  if ("${BRLCAD_BUNDLED_LIBS}" STREQUAL "BUNDLED")
    set(BEXT_ENABLE_ALL ON)
  endif ("${BRLCAD_BUNDLED_LIBS}" STREQUAL "BUNDLED")
  set(BEXT_USE_GDAL ON)
  if (NOT BRLCAD_ENABLE_GDAL)
    set(BEXT_USE_GDAL OFF)
  endif (NOT BRLCAD_ENABLE_GDAL)
  set(BEXT_USE_QT ON)
  if (NOT BRLCAD_ENABLE_QT)
    set(BEXT_USE_QT OFF)
  endif (NOT BRLCAD_ENABLE_QT)
  set(BEXT_USE_TCL ON)
  if (NOT BRLCAD_ENABLE_TCL)
    set(BEXT_USE_TCL OFF)
  endif (NOT BRLCAD_ENABLE_TCL)
  set(BEXT_USE_APPLESEED OFF)
  if (BRLCAD_ENABLE_APPLESEED)
    set(BEXT_USE_APPLESEED ON)
  endif (BRLCAD_ENABLE_APPLESEED)

  set(EXT_CONFIG_STATUS 0)
  if (BRLCAD_COMPONENTS)
    set(active_dirs ${BRLCAD_COMPONENTS})
    foreach(wc ${BRLCAD_COMPONENTS})
      deps_expand(${wc} active_dirs)
    endforeach(wc ${BRLCAD_COMPONENTS})
    string(REPLACE ";" "\\;" active_dirs "${active_dirs}")
    message("${CMAKE_COMMAND} ${BRLCAD_EXT_SOURCE_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${BRLCAD_EXT_INSTALL_DIR} -DBRLCAD_COMPONENTS=${active_dirs}")
    execute_process(COMMAND ${CMAKE_COMMAND} ${BRLCAD_EXT_SOURCE_DIR}
      -DGIT_SHALLOW_CLONE=ON
      -DENABLE_ALL=${BEXT_ENABLE_ALL}
      -DUSE_GDAL=${BEXT_USE_GDAL}
      -DUSE_QT=${BEXT_USE_QT}
      -DUSE_TCL=${BEXT_USE_TCL}
      -DUSE_APPLESEED=${BEXT_USE_APPLESEED}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_INSTALL_PREFIX=${BRLCAD_EXT_INSTALL_DIR}
      -DBRLCAD_COMPONENTS=${active_dirs}
      WORKING_DIRECTORY ${BRLCAD_EXT_BUILD_DIR}
      RESULT_VARIABLE EXT_CONFIG_STATUS
      )
  else (BRLCAD_COMPONENTS)
    execute_process(COMMAND ${CMAKE_COMMAND} ${BRLCAD_EXT_SOURCE_DIR}
      -DGIT_SHALLOW_CLONE=ON
      -DENABLE_ALL=${BEXT_ENABLE_ALL}
      -DUSE_GDAL=${BEXT_USE_GDAL}
      -DUSE_QT=${BEXT_USE_QT}
      -DUSE_TCL=${BEXT_USE_TCL}
      -DUSE_APPLESEED=${BEXT_USE_APPLESEED}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_INSTALL_PREFIX=${BRLCAD_EXT_INSTALL_DIR}
      WORKING_DIRECTORY ${BRLCAD_EXT_BUILD_DIR}
      RESULT_VARIABLE EXT_CONFIG_STATUS
      )
    message("${CMAKE_COMMAND} ${BRLCAD_EXT_SOURCE_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${BRLCAD_EXT_INSTALL_DIR}")
  endif (BRLCAD_COMPONENTS)
  if (EXT_CONFIG_STATUS)
    message(FATAL_ERROR "Unable to successfully configure bext dependency repository for building")
  endif (EXT_CONFIG_STATUS)

  set(EXT_BUILD_STATUS 0)
  if (NOT DEFINED BRLCAD_EXT_PARALLEL)
    set(BRLCAD_EXT_PARALLEL 8)
  endif (NOT DEFINED BRLCAD_EXT_PARALLEL)
  if (NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
  endif (NOT DEFINED CMAKE_BUILD_TYPE)
  message("BRLCAD_EXT_PARALLEL: ${BRLCAD_EXT_PARALLEL}")
  if (BRLCAD_EXT_PARALLEL EQUAL 1)
    if (CMAKE_CONFIGURATION_TYPES)
      execute_process(COMMAND ${CMAKE_COMMAND} --build ${BRLCAD_EXT_BUILD_DIR} --parallel 1 --config ${CMAKE_BUILD_TYPE} WORKING_DIRECTORY ${BRLCAD_EXT_BUILD_DIR} RESULT_VARIABLE EXT_BUILD_STATUS)
    else (CMAKE_CONFIGURATION_TYPES)
      execute_process(COMMAND ${CMAKE_COMMAND} --build ${BRLCAD_EXT_BUILD_DIR} --parallel 1 WORKING_DIRECTORY ${BRLCAD_EXT_BUILD_DIR} RESULT_VARIABLE EXT_BUILD_STATUS)
    endif (CMAKE_CONFIGURATION_TYPES)
  else (BRLCAD_EXT_PARALLEL EQUAL 1)
    if (CMAKE_CONFIGURATION_TYPES)
      execute_process(COMMAND ${CMAKE_COMMAND} --build ${BRLCAD_EXT_BUILD_DIR} --parallel ${BRLCAD_EXT_PARALLEL} --config ${CMAKE_BUILD_TYPE} WORKING_DIRECTORY ${BRLCAD_EXT_BUILD_DIR} RESULT_VARIABLE EXT_BUILD_STATUS)
    else (CMAKE_CONFIGURATION_TYPES)
      execute_process(COMMAND ${CMAKE_COMMAND} --build ${BRLCAD_EXT_BUILD_DIR} --parallel ${BRLCAD_EXT_PARALLEL} WORKING_DIRECTORY ${BRLCAD_EXT_BUILD_DIR} RESULT_VARIABLE EXT_BUILD_STATUS)
    endif (CMAKE_CONFIGURATION_TYPES)
  endif (BRLCAD_EXT_PARALLEL EQUAL 1)
  if (EXT_BUILD_STATUS)
    message(FATAL_ERROR "Unable to successfully build bext dependency repository")
  endif (EXT_BUILD_STATUS)

  set(BRLCAD_EXT_DIR ${BRLCAD_EXT_INSTALL_DIR}/bext_output CACHE PATH "Local bext install" FORCE)
  set(BRLCAD_EXT_INSTALL_DIR ${BRLCAD_EXT_DIR}/install CACHE PATH "Local bext install" FORCE)
  set(BRLCAD_EXT_NOINSTALL_DIR ${BRLCAD_EXT_DIR}/noinstall CACHE PATH "Local bext install" FORCE)
endfunction(brlcad_ext_setup)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
