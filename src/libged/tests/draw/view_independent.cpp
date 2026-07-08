/*         V I E W _ I N D E P E N D E N T . C P P
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
/** @file view_independent.cpp
 *
 * Phase V6 regression test: `view independent` is backed by a private
 * structural BSG scope rather than the legacy mode bit alone.
 */

#include "common.h"

#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

#include <bu.h>
#include <ged.h>
#include <rt/view.h>
#include "view_test_util.h"
#include "ged/draw.h"
#include "../../ged_draw_view_private.h"

#define ASSERT(cond) do { \
    nchecks++; \
    if (!(cond)) { \
	bu_log("FAIL [%s:%d] %s\n", __FILE__, __LINE__, #cond); \
	nfails++; \
    } \
} while (0)

static int nchecks = 0;
static int nfails = 0;

static std::vector<std::string>
drawn_paths(struct ged *gedp, void *v)
{
    std::vector<std::string> ret;
    struct bu_vls paths = BU_VLS_INIT_ZERO;
    ged_draw_list_paths(gedp, v, -1, 0, &paths);

    const char *s = bu_vls_cstr(&paths);
    const char *start = s;
    for (const char *p = s; p && *p; p++) {
	if (*p != '\n')
	    continue;
	ret.push_back(std::string(start, p - start));
	start = p + 1;
    }
    if (start && *start)
	ret.push_back(std::string(start));

    bu_vls_free(&paths);
    return ret;
}

static int
has_path(const std::vector<std::string> &paths, const char *path)
{
    return (std::find(paths.begin(), paths.end(), std::string(path)) != paths.end()) ? 1 : 0;
}

static std::string
who_view(struct ged *gedp, const char *view_name)
{
    const char *av[4] = {"who", "-V", view_name, NULL};
    bu_vls_trunc(gedp->ged_result_str, 0);
    ASSERT(ged_exec_who(gedp, 3, av) == BRLCAD_OK);
    return std::string(bu_vls_cstr(gedp->ged_result_str));
}

static int
who_has_path(const std::string &who, const char *path)
{
    size_t start = 0;
    while (start < who.size()) {
	size_t end = who.find('\n', start);
	if (end == std::string::npos)
	    end = who.size();
	if (who.compare(start, end - start, path) == 0)
	    return 1;
	start = end + 1;
    }
    return 0;
}

static int
set_view_independent(struct ged *gedp, const char *view_name, int independent)
{
    const char *av[5] = {"view", "independent", view_name, independent ? "1" : "0", NULL};
    return ged_exec_view(gedp, 4, av);
}

static int
draw_shared(struct ged *gedp, const char *path)
{
    const char *av[4] = {"draw", "-R", path, NULL};
    return ged_exec_draw(gedp, 3, av);
}

static int
draw_shared_autoview(struct ged *gedp, const char *path)
{
    const char *av[3] = {"draw", path, NULL};
    return ged_exec_draw(gedp, 2, av);
}

static int
draw_view(struct ged *gedp, const char *view_name, const char *path)
{
    const char *av[6] = {"draw", "-R", "-V", view_name, path, NULL};
    return ged_exec_draw(gedp, 5, av);
}

static int
redraw_view(struct ged *gedp, const char *view_name)
{
    const char *av[4] = {"redraw", "-V", view_name, NULL};
    return ged_exec_redraw(gedp, 3, av);
}

static int
erase_view(struct ged *gedp, const char *view_name, const char *path)
{
    const char *av[5] = {"erase", "-V", view_name, path, NULL};
    return ged_exec_erase(gedp, 4, av);
}

static int
erase_current(struct ged *gedp, const char *path)
{
    const char *av[3] = {"erase", path, NULL};
    return ged_exec_erase(gedp, 2, av);
}

static int
zap_view_db(struct ged *gedp, const char *view_name)
{
    const char *av[5] = {"zap", "-V", view_name, "-g", NULL};
    return ged_exec_zap(gedp, 4, av);
}

static int
zap_current(struct ged *gedp)
{
    const char *av[2] = {"zap", NULL};
    return ged_exec_zap(gedp, 1, av);
}

int
main(int argc, const char **argv)
{
    bu_setprogname(argv[0]);

    if (argc != 2)
	bu_exit(EXIT_FAILURE, "Usage: ged_test_view_independent <directory-containing-moss.g>\n");

    struct bu_vls gpath = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&gpath, "%s/moss.g", argv[1]);

    struct ged *gedp = ged_open("db", bu_vls_cstr(&gpath), 1);
    ASSERT(gedp != NULL);
    if (!gedp)
	return EXIT_FAILURE;

    void *view_set_ctx = ged_view_set_ctx(gedp);
    ASSERT(ged_view_set_context_remove(view_set_ctx, NULL));
    void *views[2] = {NULL, NULL};
    for (int i = 0; i < 2; i++) {
	char view_name[16];
	snprintf(view_name, sizeof(view_name), "V%d", i);
	views[i] = ged_view_context_create();
	ASSERT(views[i] != NULL);
	ASSERT(bv_name_set(DRAW_TEST_BV(views[i]), view_name));
	ASSERT(ged_view_set_context_add(view_set_ctx, views[i]));
	ged_view_context_owned_add(gedp, views[i]);
	if (!i)
	    ged_view_active_ctx_set(gedp, views[i]);
    }

    for (int i = 0; i < 2; i++) {
	ASSERT(bv_scale_state_set(DRAW_TEST_BV(views[i]), 1.0e9, 1.0, 0.0, 1.0e9, 1.0 / 1.0e9));
    }
    ASSERT(draw_shared_autoview(gedp, "all.g") == BRLCAD_OK);
    ASSERT(ged_draw_view_context_lod_bounds_callback_is(views[0]));
    ASSERT(ged_draw_view_context_lod_bounds_callback_is(views[1]));
    ASSERT(bv_size_get(DRAW_TEST_BV_CONST(views[0])) < 1.0e8);
    ASSERT(bv_size_get(DRAW_TEST_BV_CONST(views[1])) < 1.0e8);
    ASSERT(zap_current(gedp) == BRLCAD_OK);
    ASSERT(drawn_paths(gedp, views[0]).size() == 0);
    ASSERT(drawn_paths(gedp, views[1]).size() == 0);

    ASSERT(draw_shared(gedp, "all.g") == BRLCAD_OK);
    ASSERT(!ged_view_context_is_independent(views[0]));
    ASSERT(drawn_paths(gedp, views[0]).size() == 1);
    ASSERT(has_path(drawn_paths(gedp, views[0]), "all.g"));
    ASSERT(drawn_paths(gedp, views[1]).size() == 1);

    ASSERT(set_view_independent(gedp, "V0", 1) == BRLCAD_OK);
    ASSERT(ged_view_context_is_independent(views[0]));
    ASSERT(!ged_view_context_independent_scope_is_null(views[0], 0));
    ASSERT(drawn_paths(gedp, views[0]).size() == 1);
    ASSERT(has_path(drawn_paths(gedp, views[0]), "all.g"));

    ASSERT(draw_shared(gedp, "box.r") == BRLCAD_OK);
    ASSERT(drawn_paths(gedp, views[1]).size() == 2);
    ASSERT(has_path(drawn_paths(gedp, views[1]), "box.r"));
    ASSERT(drawn_paths(gedp, views[0]).size() == 1);
    ASSERT(!has_path(drawn_paths(gedp, views[0]), "box.r"));
    ASSERT(who_has_path(who_view(gedp, "V1"), "box.r"));
    ASSERT(!who_has_path(who_view(gedp, "V0"), "box.r"));

    ASSERT(draw_view(gedp, "V0", "tor.r") == BRLCAD_OK);
    ASSERT(drawn_paths(gedp, views[0]).size() == 2);
    ASSERT(has_path(drawn_paths(gedp, views[0]), "tor.r"));
    ASSERT(!has_path(drawn_paths(gedp, views[1]), "tor.r"));
    ASSERT(who_has_path(who_view(gedp, "V0"), "tor.r"));
    ASSERT(!who_has_path(who_view(gedp, "V1"), "tor.r"));
    ASSERT(redraw_view(gedp, "V0") == BRLCAD_OK);
    ASSERT(has_path(drawn_paths(gedp, views[0]), "tor.r"));
    ASSERT(!has_path(drawn_paths(gedp, views[1]), "tor.r"));
    ASSERT(erase_view(gedp, "V0", "tor.r") == BRLCAD_OK);
    ASSERT(!has_path(drawn_paths(gedp, views[0]), "tor.r"));
    ASSERT(has_path(drawn_paths(gedp, views[0]), "all.g"));
    ASSERT(has_path(drawn_paths(gedp, views[1]), "box.r"));

    ASSERT(draw_view(gedp, "V0", "tor.r") == BRLCAD_OK);
    ASSERT(has_path(drawn_paths(gedp, views[0]), "tor.r"));
    ASSERT(zap_view_db(gedp, "V0") == BRLCAD_OK);
    ASSERT(drawn_paths(gedp, views[0]).size() == 0);
    ASSERT(has_path(drawn_paths(gedp, views[1]), "all.g"));
    ASSERT(has_path(drawn_paths(gedp, views[1]), "box.r"));
    ASSERT(!has_path(drawn_paths(gedp, views[1]), "tor.r"));

    ASSERT(set_view_independent(gedp, "V0", 0) == BRLCAD_OK);
    ASSERT(!ged_view_context_is_independent(views[0]));
    ASSERT(ged_view_context_independent_scope_is_null(views[0], 0));
    ASSERT(drawn_paths(gedp, views[0]).size() == 2);
    ASSERT(has_path(drawn_paths(gedp, views[0]), "all.g"));
    ASSERT(has_path(drawn_paths(gedp, views[0]), "box.r"));
    ASSERT(!has_path(drawn_paths(gedp, views[0]), "tor.r"));
    ASSERT(erase_current(gedp, "box.r") == BRLCAD_OK);
    ASSERT(!has_path(drawn_paths(gedp, views[0]), "box.r"));
    ASSERT(!has_path(drawn_paths(gedp, views[1]), "box.r"));
    ASSERT(has_path(drawn_paths(gedp, views[0]), "all.g"));
    ASSERT(zap_current(gedp) == BRLCAD_OK);
    ASSERT(drawn_paths(gedp, views[0]).size() == 0);
    ASSERT(drawn_paths(gedp, views[1]).size() == 0);

    bu_vls_free(&gpath);
    ged_close(gedp);

    bu_log("view_independent: %d checks, %d failures\n", nchecks, nfails);
    return nfails ? EXIT_FAILURE : EXIT_SUCCESS;
}
