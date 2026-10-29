/*                     T E S T _ A R G V . C
 * BRL-CAD
 *
 * Copyright (c) 2014-2026 United States Government as represented by
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "bu.h"

#include "vmath.h"


static size_t
argv_test_limited(const char *s, size_t expected, size_t lmax, const char * const expected_argv[])
{
    char *ts = bu_strdup(s);
    char **av = (char **)bu_calloc(lmax+1, sizeof(char *), "av");
    /* Require the parser to write the terminator rather than rely on calloc. */
    av[expected] = ts;
    size_t ac = bu_argv_from_string(av, lmax, ts);
    if (ac != expected)
	bu_exit(EXIT_FAILURE, "Test failed: %s\n", s);
    bu_log("Test results (%zd elements): %s\n", ac,  s);
    for (size_t i = 0; i < ac; i++) {
	bu_log("%zd\t: %s\n", i, av[i]);
	if (expected_argv && !BU_STR_EQUAL(av[i], expected_argv[i]))
	    bu_exit(EXIT_FAILURE, "Test failed: %s, argument %zu: expected '%s', got '%s'\n",
		    s, i, expected_argv[i], av[i]);
    }
    if (av[ac] != NULL)
	bu_exit(EXIT_FAILURE, "Test failed: %s, missing argv terminator\n", s);
    bu_free(av, "array");
    bu_free(ts, "string copy");
    return ac;
}


static size_t
argv_test(const char *s, size_t expected)
{
    return argv_test_limited(s, expected, strlen(s), NULL);
}


int
main(int UNUSED(argc), char *argv[])
{
    // Normally this file is part of bu_test, so only set this if it
    // looks like the program name is still unset.
    if (bu_getprogname()[0] == '\0')
	bu_setprogname(argv[0]);

    argv_test("a0", 1);
    argv_test("a0 a1", 2);
    argv_test("a0 a1 a2", 3);
    argv_test("a0 \"a 1\"", 2);
    argv_test("\"a 0\" a1", 2);
    argv_test("\"a 0\" \"a 1\"", 2);
    argv_test("a0 \"a 1\" a2", 3);
    argv_test("a0 a\\ 1", 2);
    argv_test("a\\ 0 a1", 2);
    argv_test("a\\ 0 a\\ 1", 2);
    argv_test("a\\ 0 \"a 1\"", 2);
    argv_test("\"a 0\" a\\ 1", 2);
    argv_test("\\\"a 0\\\" a\\ 1", 3);
    argv_test("a\\ 0 \\\"a 1\\\"", 3);
    argv_test("a\\\\ 0 a\\ 1", 3);
    argv_test("a\\ 0 a\\\\ 1", 3);
    argv_test("\"a\\ 0\" a\\ 1", 2);
    argv_test("\"a\\\\ 0\" a\\ 1", 2);


    argv_test_limited("a0 a1 a2", 0, 0, NULL);
    argv_test_limited("a0 a1 a2", 1, 1, NULL);
    argv_test_limited("a0 a1 a2", 2, 2, NULL);
    argv_test_limited("a0 a1 a2", 3, 3, NULL);
    argv_test_limited("a0 a1 a2", 3, 4, NULL);

    /* Empty quoted arguments must not absorb following arguments.  Check
     * values as well as counts so a lost quote cannot go unnoticed. */
    {
	const struct {
	    const char *input;
	    const char *expected[5];
	} tests[] = {
	    {"", {NULL}},
	    {" \t\r\n", {NULL}},
	    {"a0 a1 ", {"a0", "a1", NULL}},
	    {"a0 a1\t", {"a0", "a1", NULL}},
	    {" r1 r2 ", {"r1", "r2", NULL}},
	    {"\"\"", {"", NULL}},
	    {"\"\" \"\"", {"", "", NULL}},
	    {"\"\" a1", {"", "a1", NULL}},
	    {"a0 \"\"", {"a0", "", NULL}},
	    {"a0 \"\"  ", {"a0", "", NULL}},
	    {"a0 \"\" a2", {"a0", "", "a2", NULL}},
	    {"\"a0\" \"\" \"a2\"", {"a0", "", "a2", NULL}},
	    {"a0 \"\" \"\" a3", {"a0", "", "", "a3", NULL}},
	    {"check_cycle \"\" \"{}\" \";\"", {"check_cycle", "", "{}", ";", NULL}},
	    {"a0 \\\"a1\\\"", {"a0", "\"a1\"", NULL}},
	    {"a0 \\ a1", {"a0", " a1", NULL}},
	    {"a0 a1\\", {"a0", "a1", NULL}},
	    {"a0 caf\xc3\xa9", {"a0", "caf\xc3\xa9", NULL}},
	    {"a0 \"\\\"a1\\\"\"", {"a0", "\"a1\"", NULL}}
	};
	for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); i++) {
	    size_t expected = 0;
	    while (tests[i].expected[expected])
		expected++;
	    for (size_t limit = 0; limit <= expected + 1; limit++)
		argv_test_limited(tests[i].input, limit < expected ? limit : expected,
				  limit, tests[i].expected);
	}
    }

    return 0;
}


/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8 cino=N-s
 */
