# Copyright (c) 2026 United States Government as represented by
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

include_guard(GLOBAL)

set(_BRLCAD_REPOSITORY_DEFAULT_ROOT "https://github.com/BRL-CAD")
set(BRLCAD_REPOSITORY_ROOT "${_BRLCAD_REPOSITORY_DEFAULT_ROOT}" CACHE STRING
  "Repository URL or local mirror directory containing bext and its dependency repositories")

function(brlcad_repository_root out_root out_local)
  set(root "${BRLCAD_REPOSITORY_ROOT}")
  set(local FALSE)
  if(root STREQUAL "")
    message(FATAL_ERROR "BRLCAD_REPOSITORY_ROOT must be a repository URL or local mirror directory")
  endif()

  # Match bext's URL handling, but resolve relative paths against BRL-CAD.
  # Drive paths must be checked before scp-style SSH locations.
  if(root MATCHES "^file://")
    set(local TRUE)
    string(REPLACE "=" "%3D" root "${root}")
  elseif(root MATCHES "^[A-Za-z]:[/\\\\]" OR
      NOT (root MATCHES "^[A-Za-z][A-Za-z0-9+.-]*://" OR root MATCHES "^[^/]+:.+"))
    set(local TRUE)
    string(REPLACE "\\" "/" root "${root}")
    if(WIN32 OR NOT root MATCHES "^[A-Za-z]:/")
      get_filename_component(root "${root}" ABSOLUTE BASE_DIR "${CMAKE_SOURCE_DIR}")
    endif()
    # Escape literal percent signs first so native paths cannot introduce
    # URL escapes.  An equals sign also has meaning in Git's -c assignments.
    string(REPLACE "%" "%25" root "${root}")
    string(REPLACE " " "%20" root "${root}")
    string(REPLACE "#" "%23" root "${root}")
    string(REPLACE "?" "%3F" root "${root}")
    string(REPLACE "=" "%3D" root "${root}")
    if(root MATCHES "^[A-Za-z]:/")
      set(root "file:///${root}")
    elseif(root MATCHES "^//")
      set(root "file:${root}")
    else()
      set(root "file://${root}")
    endif()
  endif()

  if(NOT root MATCHES "^file:///*$")
    string(REGEX REPLACE "/+$" "" root "${root}")
    string(APPEND root "/")
  else()
    set(root "file:///")
  endif()
  set(${out_root} "${root}" PARENT_SCOPE)
  set(${out_local} "${local}" PARENT_SCOPE)
endfunction()

function(brlcad_repository_git_command out_command)
  set(git_command "${GIT_EXEC}")
  if(WIN32)
    list(APPEND git_command -c core.longpaths=true)
  endif()
  if(_BRLCAD_REPOSITORY_LOCAL)
    list(APPEND git_command -c protocol.file.allow=always)
    # This also constrains recursive Git processes and user URL rewrites.
    list(PREPEND git_command "${CMAKE_COMMAND}" -E env GIT_ALLOW_PROTOCOL=file --)
  endif()
  set(${out_command} "${git_command}" PARENT_SCOPE)
endfunction()

brlcad_repository_root(_BRLCAD_REPOSITORY_ROOT _BRLCAD_REPOSITORY_LOCAL)
set(_BRLCAD_BEXT_REPOSITORY "${_BRLCAD_REPOSITORY_ROOT}bext.git")
message(STATUS "BRL-CAD repository root: ${_BRLCAD_REPOSITORY_ROOT}")
if(_BRLCAD_REPOSITORY_LOCAL)
  message(STATUS "BRL-CAD Git acquisition: file transport only")
endif()

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
