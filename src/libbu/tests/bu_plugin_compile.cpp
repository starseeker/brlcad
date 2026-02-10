/*                B U _ P L U G I N _ C O M P I L E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file bu_plugin_compile.cpp
 *
 * Compile-time test for bu_plugin.h to verify it compiles correctly
 * with BU_PLUGIN_IMPLEMENTATION defined, particularly on Windows where
 * the C2491 error ("definition of dllimport function not allowed") would
 * occur if BU_PLUGIN_API incorrectly defaulted to dllimport.
 *
 * This test simulates the problematic scenario from ged_init.cpp where:
 * - BU_PLUGIN_IMPLEMENTATION is defined (to compile the implementation)
 * - Windows DLL macros might be present in the environment
 * - We need to ensure BU_PLUGIN_API does NOT result in dllimport decoration
 *   when compiling the implementation
 *
 * Success: This file compiles without errors
 * Failure: Compilation errors like C2491 (MSVC) or similar errors on other compilers
 */

#include "common.h"

/*
 * Simulate the scenario from ged_init.cpp:
 * Define BU_PLUGIN_IMPLEMENTATION to include the implementation.
 * Do NOT define BU_PLUGIN_DLL_EXPORTS or BU_PLUGIN_DLL_IMPORTS.
 * This should result in BU_PLUGIN_API being empty on Windows (no dllimport/dllexport).
 */
#define BU_PLUGIN_IMPLEMENTATION
#include "../bu_plugin.h"

/*
 * Test that we can call bu_plugin functions without linker/compiler errors.
 * If BU_PLUGIN_API was incorrectly set to dllimport, this would fail on Windows
 * with error C2491.
 */
static int test_cmd(void) {
    return 42;
}

int main(int argc, const char **argv) {
    (void)argc;
    (void)argv;

    /* Initialize the plugin system */
    int init_result = bu_plugin_init();
    if (init_result != 0) {
        return 1;
    }

    /* Register a test command */
    int reg_result = bu_plugin_cmd_register("test_cmd", test_cmd);
    if (reg_result != 0) {
        return 1;
    }

    /* Verify command exists */
    if (!bu_plugin_cmd_exists("test_cmd")) {
        return 1;
    }

    /* Get command count */
    size_t count = bu_plugin_cmd_count();
    if (count == 0) {
        return 1;
    }

    /* Test passed - compilation succeeded and basic functionality works */
    return 0;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
