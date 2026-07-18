# Copyright (c) 2010-2026 United States Government as represented by
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

# *******************************************************************
# ***                   env2c CMakeLists.txt                      ***
# *******************************************************************

# Minimum required version of CMake
cmake_minimum_required(VERSION 3.20)

# set CMake project name
project(ENV2C)

include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

find_package(Threads)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-O3 O3_COMPILER_FLAG)

add_executable(env2c env2c.cpp)
target_link_libraries(env2c Threads::Threads)
if(BRLCAD_ENABLE_ADDRESS_SANITIZER OR BRLCAD_ENABLE_UNDEFINED_SANITIZER OR BRLCAD_SANITIZE_THREAD)
  # This executable runs during the build to generate C sources.  Keep host
  # build tools native; sanitizer coverage belongs on the generated targets.
  target_compile_options(env2c PRIVATE "-fno-sanitize=all")
  target_link_options(env2c PRIVATE "-fno-sanitize=all")
endif()

set(STACK_SIZE_FLAGS "-Wl,-z,stack-size=4194304")
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.17")
  include(CheckLinkerFlag)
  check_linker_flag("CXX" "${STACK_SIZE_FLAGS}" HAVE_STACK_SIZE_FLAG)
else(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.17")
  include(CheckCXXCompilerFlag)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${STACK_SIZE_FLAGS}")
  check_cxx_compiler_flag("${STACK_SIZE_FLAGS}" HAVE_STACK_SIZE_FLAG)
endif(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.17")
if(HAVE_STACK_SIZE_FLAG)
  target_link_options(env2c PRIVATE "${STACK_SIZE_FLAGS}")
endif(HAVE_STACK_SIZE_FLAG)

if(O3_COMPILER_FLAG AND BRLCAD_OPTIMIZED)
  # If we have the O3 flag, use it
  target_compile_options(env2c PRIVATE "-O3")
endif()

if(COMMAND CMAKEFILES)
  cmakefiles(
    CMakeLists.txt
    env2c.cpp
  )
endif()

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
