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
#include <cstdlib>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "bu/app.h"
#include "bu/log.h"
#include "../ged_private.h"

static void
check(bool condition, const char *message)
{
    if (!condition) {
	bu_log("terminal columns: %s\n", message);
	exit(EXIT_FAILURE);
    }
}

static size_t
line_count(const std::string &text)
{
    size_t count = 0;
    for (char c : text) if (c == '\n') count++;
    return count;
}

int
main(int, char **argv)
{
    bu_setprogname(argv[0]);
    struct ged *gedp = ged_open("inmem", "0", 0);
    struct ged *other = ged_create();
    check(gedp && other, "session initialization failed");
    check(ged_terminal_columns(gedp) == 80, "default width changed");
    std::set<std::string> names;
    for (int i = 0; i < 9; i++) {
	std::string name = "object_" + std::to_string(i);
	const char *make[] = {"make", name.c_str(), "sph"};
	check(ged_exec(gedp, 3, make) == BRLCAD_OK, "could not create test geometry");
	names.insert(name);
    }
    size_t previous_lines = names.size();
    for (size_t width : {size_t(1), size_t(40), size_t(80), size_t(160)}) {
	ged_set_terminal_columns(gedp, width);
	check(ged_terminal_columns(other) == 80, "width leaked between sessions");
	const char *ls[] = {"ls"};
	check(ged_exec(gedp, 1, ls) == BRLCAD_OK, "ls failed");
	std::string text = bu_vls_cstr(gedp->ged_result_str);
	std::istringstream input(text);
	std::set<std::string> found;
	std::string word;
	while (input >> word) check(found.insert(word).second, "listing duplicated a name");
	check(found == names, "listing lost a name");
	check(line_count(text) <= previous_lines, "wider listing used more rows");
	previous_lines = line_count(text);
	std::istringstream rows(text);
	while (std::getline(rows, word))
	    check(word.size() <= width || width < names.begin()->size(), "listing exceeds width");
    }
    for (const char *option : {"-p", "-l"}) {
	const char *ls[] = {"ls", option};
	ged_set_terminal_columns(gedp, 1);
	check(ged_exec(gedp, 2, ls) == BRLCAD_OK, "ls option failed");
	std::string narrow = bu_vls_cstr(gedp->ged_result_str);
	ged_set_terminal_columns(gedp, 160);
	check(ged_exec(gedp, 2, ls) == BRLCAD_OK, "ls option failed");
	check(narrow == bu_vls_cstr(gedp->ged_result_str), "explicit listing format changed with width");
    }
    ged_set_terminal_columns(gedp, 0);
    check(ged_terminal_columns(gedp) == 80, "zero did not restore default width");
    ged_set_terminal_columns(NULL, 20);
    check(ged_terminal_columns(NULL) == 80, "null session has no default width");

    std::vector<std::string> labels = {"short", std::string(18, 'r'), std::string(100, 'x')};
    struct directory entries[3] = {};
    struct directory *pointers[3];
    for (size_t i = 0; i < labels.size(); i++) {
	entries[i].d_namep = const_cast<char *>(labels[i].c_str());
	entries[i].d_flags = RT_DIR_COMB | RT_DIR_REGION;
	pointers[i] = &entries[i];
    }
    struct bu_vls text = BU_VLS_INIT_ZERO;
    for (size_t width : {size_t(1), size_t(80), std::numeric_limits<size_t>::max()}) {
	bu_vls_trunc(&text, 0);
	_ged_vls_col_pr4v(&text, pointers, 3, 0, 0, width);
	for (const auto &label : labels)
	    check(std::string(bu_vls_cstr(&text)).find(label + "/R") != std::string::npos,
		"decorated or long name was truncated");
    }
    bu_vls_trunc(&text, 0);
    _ged_vls_col_pr4v(&text, pointers, 0, 0, 0, 1);
    check(bu_vls_strlen(&text) == 0, "empty listing produced output");
    bu_vls_free(&text);
    ged_close(gedp);
    ged_destroy(other);
    return EXIT_SUCCESS;
}
