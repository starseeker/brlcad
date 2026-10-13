# BRL-CAD
#
# Copyright (c) 2026 United States Government as represented by
# the U.S. Army Research Laboratory.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# version 2.1 as published by the Free Software Foundation.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this file; see the file named COPYING for more
# information.
#
file(MAKE_DIRECTORY "${WORK}")
set(database "${WORK}/batch.g")
set(input "${WORK}/input.tcl")
file(REMOVE "${database}")
file(WRITE "${input}" [=[
in first.s sph 0 0 0 1
in second.s sph 5 0 0 1
g all *.s
set definition [get all]
if {[string first first.s $definition] < 0 || [string first second.s $definition] < 0} {
    error "Geometry globbing did not create the expected combination"
}
set lines {first
second}
puts $lines
proc visit {path} {attr set $path visited yes; return 1}
search first.s -exec visit "{}" ";"
if {[attr get first.s visited] ne "yes"} {error "Search did not invoke a Tcl procedure"}
puts BATCH_PASS
]=])
execute_process(COMMAND "${FLTKGED}" -c "${database}" "source {${input}}"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE errors)
if(NOT status EQUAL 0 OR NOT output MATCHES "BATCH_PASS")
  message(FATAL_ERROR "Batch Tcl script failed: ${status}
${output}
${errors}")
endif()
execute_process(COMMAND "${FLTKGED}" -c "${database}" "exit 7; error {exit did not stop evaluation}"
  RESULT_VARIABLE status)
if(NOT status EQUAL 7)
  message(FATAL_ERROR "Tcl exit did not preserve its exit status: ${status}")
endif()
set(marker "${WORK}/quit-reached")
file(REMOVE "${marker}")
execute_process(COMMAND "${FLTKGED}" -c "${database}" "catch {quit}; close [open {${marker}} w]"
  RESULT_VARIABLE status)
if(NOT status EQUAL 0 OR EXISTS "${marker}")
  message(FATAL_ERROR "Tcl quit did not stop a caught script: ${status}")
endif()
file(WRITE "${input}" "set unfinished {\n")
execute_process(COMMAND "${FLTKGED}" -c "${database}" INPUT_FILE "${input}"
  RESULT_VARIABLE status)
if(NOT status EQUAL 1)
  message(FATAL_ERROR "An incomplete Tcl command was accepted")
endif()
execute_process(COMMAND "${FLTKGED}" -c --shell=ged "${database}" size 42
  RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Raw GED batch command failed")
endif()
execute_process(COMMAND "${FLTKGED}" -c --shell=ged "${database}" size
  RESULT_VARIABLE status OUTPUT_VARIABLE output)
if(NOT status EQUAL 0 OR NOT output MATCHES "1000")
  message(FATAL_ERROR "Batch view was not initialized")
endif()
file(REMOVE "${database}" "${input}")
