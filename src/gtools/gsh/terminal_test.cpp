/* BRL-CAD
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

#include "common.h"
#include "bu/app.h"
#include "linenoise.hpp"

int
main(int, char **argv)
{
    bu_setprogname(argv[0]);
    linenoise::linenoiseState terminal("test> ", 0, 1);
    std::string line;
    if (terminal.Readline(line)) return EXIT_FAILURE;
    std::cout << "RESULT:" << line << std::endl;
    return EXIT_SUCCESS;
}
