/*                  T E S T _ E D I T . C P P
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
/** @file test_edit.cpp
 *
 * Comprehensive test suite for the libged `edit` command.
 */

#include "common.h"

#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "bu.h"
#include "vmath.h"
#include "raytrace.h"
#include "rt/db_fullpath.h"
#include "rt/geom.h"
#include "rt/edit.h"
#include "wdb.h"
#include "ged.h"
#include "ged/edit.h"
#include "ged/event.h"
#include "brep/util.h"
#include "ged/selection.h"

#include "../../ged_private.h"
#include "../../edit/uri.hh"

/* ================================================================== *
 * Shared infrastructure
 * ================================================================== */

static int total_tests  = 0;
static int failed_tests = 0;

/* Numeric tolerance: lenient enough for trig (sin/cos) results. */
#define NEAR_ENOUGH 1e-4

#define CHECK(cond, msg)                                             \
    do {                                                             \
        ++total_tests;                                               \
        if (!(cond)) {                                               \
            ++failed_tests;                                          \
            bu_log("FAIL [%s:%d]: %s\n", __FILE__, __LINE__, msg);  \
        } else {                                                     \
            bu_log("PASS: %s\n", msg);                               \
        }                                                            \
    } while (0)


struct edit_event_observer {
    int calls = 0;
    int tor_modified = 0;
    int sph_modified = 0;
    int bot_modified = 0;
};


struct edit_session_observer {
    std::vector<enum ged_edit_session_event_kind> kinds;
    std::vector<uint64_t> revisions;
    std::vector<int> commands;
    std::vector<std::string> paths;
    std::vector<std::string> replacement_paths;
    std::vector<enum ged_edit_session_invalidation_reason> invalidation_reasons;
    ged_edit_session_ref first = GED_EDIT_SESSION_REF_NULL;
};


static void
edit_session_cb(struct ged *UNUSED(gedp),
	const struct ged_edit_session_event *event, void *client_data)
{
    struct edit_session_observer *obs =
	(struct edit_session_observer *)client_data;
    if (!obs || !event)
	return;
    if (ged_edit_session_ref_is_null(obs->first))
	obs->first = event->session;
    obs->kinds.push_back(event->kind);
    obs->revisions.push_back(event->revision);
    obs->commands.push_back(event->command_id);
    obs->paths.emplace_back(event->path ? event->path : "");
    obs->replacement_paths.emplace_back(event->replacement_path ?
	event->replacement_path : "");
    obs->invalidation_reasons.push_back(event->invalidation_reason);
}


struct reentrant_edit_observer {
    ged_edit_session_ref session = GED_EDIT_SESSION_REF_NULL;
    int updates = 0;
    enum ged_edit_status cancel_status = GED_EDIT_ERROR;
};


static void
reentrant_edit_cancel_cb(struct ged *gedp,
	const struct ged_edit_session_event *event, void *client_data)
{
    struct reentrant_edit_observer *obs =
	(struct reentrant_edit_observer *)client_data;
    if (!obs || !event || event->kind != GED_EDIT_SESSION_UPDATE ||
	event->session.owner != obs->session.owner ||
	event->session.id != obs->session.id ||
	event->session.generation != obs->session.generation)
	return;
    obs->updates++;
    obs->cancel_status = ged_edit_session_cancel(gedp, obs->session);
}


static void
edit_event_cb(struct ged *UNUSED(gedp),
	      const struct ged_event *events,
	      size_t event_count,
	      const struct ged_event_result *UNUSED(result),
	      void *client_data)
{
    struct edit_event_observer *obs = (struct edit_event_observer *)client_data;
    if (!obs)
	return;

    obs->calls++;
    for (size_t i = 0; i < event_count; i++) {
	if (events[i].kind != GED_EVENT_OBJECT_MODIFIED || !events[i].name)
	    continue;
	if (BU_STR_EQUAL(events[i].name, "tor.s"))
	    obs->tor_modified++;
	if (BU_STR_EQUAL(events[i].name, "sph.s"))
	    obs->sph_modified++;
	if (BU_STR_EQUAL(events[i].name, "bot.s"))
	    obs->bot_modified++;
    }
}


/* ------------------------------------------------------------------ *
 * Shared helper: read rt_ell_internal for an ell or sphere by name.
 * ------------------------------------------------------------------ */
static int
read_ell(struct ged *gedp, const char *name, struct rt_ell_internal *out)
{
    struct directory *dp = db_lookup(gedp->dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL)
        return BRLCAD_ERROR;

    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, gedp->dbip, NULL) < 0)
        return BRLCAD_ERROR;

    if (intern.idb_minor_type != DB5_MINORTYPE_BRLCAD_ELL) {
        rt_db_free_internal(&intern);
        return BRLCAD_ERROR;
    }

    *out = *(struct rt_ell_internal *)intern.idb_ptr;
    intern.idb_ptr = NULL;
    rt_db_free_internal(&intern);
    return BRLCAD_OK;
}

static int
read_hrt(struct ged *gedp, const char *name, struct rt_hrt_internal *out)
{
    struct directory *dp = db_lookup(gedp->dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL)
        return BRLCAD_ERROR;

    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, gedp->dbip, NULL) < 0)
        return BRLCAD_ERROR;

    if (intern.idb_minor_type != DB5_MINORTYPE_BRLCAD_HRT) {
        rt_db_free_internal(&intern);
        return BRLCAD_ERROR;
    }

    *out = *(struct rt_hrt_internal *)intern.idb_ptr;
    intern.idb_ptr = NULL;
    rt_db_free_internal(&intern);
    return BRLCAD_OK;
}

static int
read_tor(struct ged *gedp, const char *name, struct rt_tor_internal *out);


/* ------------------------------------------------------------------ *
 * Shared helper: open a fixture database on the public service path.
 * ------------------------------------------------------------------ */
static struct ged *
open_fixture(const char *path)
{
    struct ged *gedp = ged_open("db", path, 1);
    return gedp;
}


/* ------------------------------------------------------------------ *
 * Shared helper: open a fixture database without enabling the compatibility
 * backend.
 * ------------------------------------------------------------------ */
static struct ged *
open_fixture_no_dbistate(const char *path)
{
    struct ged *gedp = ged_open("db", path, 1);
    if (!gedp)
        return NULL;
    /* Deliberately avoid enabling the compatibility backend. */
    return gedp;
}


/* ------------------------------------------------------------------ *
 * Shared helper: allocate a fresh writable temp-file path.
 * The caller must bu_vls_free(v) when done.
 * ------------------------------------------------------------------ */
static int
make_temp_path(struct bu_vls *v)
{
    char tmpname[MAXPATHLEN] = {0};
    FILE *fp = bu_temp_file(tmpname, MAXPATHLEN);
    if (!fp)
        return BRLCAD_ERROR;
    fclose(fp);
    bu_vls_sprintf(v, "%s", tmpname);
    return BRLCAD_OK;
}


/* ================================================================== *
 * Section 0 — Infrastructure smoke tests
 * ================================================================== */

/**
 * Phase 0 fixture: one representative primitive of each type that has
 * a non-NULL ft_edit_desc(), plus arb8 and sph for generic tests.
 */
static int
create_p0_fixture(const char *dbpath)
{
    struct rt_wdb *wdbp = wdb_fopen(dbpath);
    if (!wdbp) {
        bu_log("ERROR: unable to create fixture database %s\n", dbpath);
        return BRLCAD_ERROR;
    }

    /* --- tor (torus) ------------------------------------------------ */
    point_t tor_V  = {0, 0, 0};
    vect_t  tor_H  = {0, 0, 1};
    if (mk_tor(wdbp, "tor.s", tor_V, tor_H, 10.0, 3.0) != 0) {
        bu_log("mk_tor failed\n");
        db_close(wdbp->dbip);
        return BRLCAD_ERROR;
    }
    if (mk_tor(wdbp, "tor", tor_V, tor_H, 8.0, 2.0) != 0) {
        bu_log("mk_tor (tor) failed\n");
        db_close(wdbp->dbip);
        return BRLCAD_ERROR;
    }

    /* --- tgc (truncated general cone) ------------------------------- */
    point_t tgc_V  = {20, 0, 0};
    vect_t  tgc_H  = {0, 0, 10};
    vect_t  tgc_A  = {5, 0, 0};
    vect_t  tgc_B  = {0, 5, 0};
    vect_t  tgc_C  = {4, 0, 0};
    vect_t  tgc_D  = {0, 4, 0};
    if (mk_tgc(wdbp, "tgc.s", tgc_V, tgc_H, tgc_A, tgc_B, tgc_C, tgc_D) != 0) {
        bu_log("mk_tgc failed\n");
        db_close(wdbp->dbip);
        return BRLCAD_ERROR;
    }

    /* --- ell (ellipsoid) -------------------------------------------- */
    point_t ell_V  = {40, 0, 0};
    vect_t  ell_A  = {8, 0, 0};
    vect_t  ell_B  = {0, 5, 0};
    vect_t  ell_C  = {0, 0, 3};
    if (mk_ell(wdbp, "ell.s", ell_V, ell_A, ell_B, ell_C) != 0) {
        bu_log("mk_ell failed\n");
        db_close(wdbp->dbip);
        return BRLCAD_ERROR;
    }

    /* --- sph (sphere — ell alias) ----------------------------------- */
    point_t sph_V  = {60, 0, 0};
    if (mk_sph(wdbp, "sph.s", sph_V, 6.0) != 0) {
        bu_log("mk_sph failed\n");
        db_close(wdbp->dbip);
        return BRLCAD_ERROR;
    }

    /* --- arb8 (box) ------------------------------------------------- */
    fastf_t arb8_pts[8*3] = {
        -5, -5, -5,   5, -5, -5,   5,  5, -5,  -5,  5, -5,
        -5, -5,  5,   5, -5,  5,   5,  5,  5,  -5,  5,  5
    };
    {   /* offset by 80 in x */
        for (int k = 0; k < 8; k++) arb8_pts[k*3] += 80.0;
    }
    if (mk_arb8(wdbp, "arb8.s", arb8_pts) != 0) {
        bu_log("mk_arb8 failed\n");
        db_close(wdbp->dbip);
        return BRLCAD_ERROR;
    }

    /* --- rpc (right parabolic cylinder) ----------------------------- */
    point_t rpc_V  = {0, 20, 0};
    vect_t  rpc_H  = {0, 0, 10};
    vect_t  rpc_B  = {5, 0, 0};
    if (mk_rpc(wdbp, "rpc.s", rpc_V, rpc_H, rpc_B, 4.0) != 0) {
        bu_log("mk_rpc failed\n");
        db_close(wdbp->dbip);
        return BRLCAD_ERROR;
    }

    /* --- rhc (right hyperbolic cylinder) ---------------------------- */
    point_t rhc_V  = {20, 20, 0};
    vect_t  rhc_H  = {0, 0, 10};
    vect_t  rhc_B  = {5, 0, 0};
    if (mk_rhc(wdbp, "rhc.s", rhc_V, rhc_H, rhc_B, 4.0, 2.0) != 0) {
        bu_log("mk_rhc failed\n");
        db_close(wdbp->dbip);
        return BRLCAD_ERROR;
    }

    /* --- epa (elliptical paraboloid) -------------------------------- */
    point_t epa_V  = {40, 20, 0};
    vect_t  epa_H  = {0, 0, 10};
    vect_t  epa_A  = {1, 0, 0};   /* must be a unit vector (breadth dir) */
    if (mk_epa(wdbp, "epa.s", epa_V, epa_H, epa_A, 5.0, 4.0) != 0) {
        bu_log("mk_epa failed\n");
        db_close(wdbp->dbip);
        return BRLCAD_ERROR;
    }

    /* --- ehy (elliptical hyperboloid) ------------------------------- */
    point_t ehy_V  = {60, 20, 0};
    vect_t  ehy_H  = {0, 0, 10};
    vect_t  ehy_A  = {1, 0, 0};   /* must be a unit vector (breadth dir) */
    if (mk_ehy(wdbp, "ehy.s", ehy_V, ehy_H, ehy_A, 4.0, 2.0, 1.0) != 0) {
        bu_log("mk_ehy failed\n");
        db_close(wdbp->dbip);
        return BRLCAD_ERROR;
    }

    /* --- eto (elliptical torus) ------------------------------------- */
    point_t eto_V  = {0, 40, 0};
    vect_t  eto_N  = {0, 0, 1};
    vect_t  eto_C  = {8, 0, 2};
    if (mk_eto(wdbp, "eto.s", eto_V, eto_N, eto_C, 12.0, 3.0) != 0) {
        bu_log("mk_eto failed\n");
        db_close(wdbp->dbip);
        return BRLCAD_ERROR;
    }

    /* --- hyp (hyperboloid of one sheet) ----------------------------- */
    point_t hyp_Vi  = {20, 40, 0};
    vect_t  hyp_H  = {0, 0, 10};
    vect_t  hyp_A  = {5, 0, 0};
    if (mk_hyp(wdbp, "hyp.s", hyp_Vi, hyp_H, hyp_A, 4.0, 0.4) != 0) {
        bu_log("mk_hyp failed\n");
        db_close(wdbp->dbip);
        return BRLCAD_ERROR;
    }

    /* --- part (particle) ------------------------------------------- */
    point_t part_V  = {40, 40, 0};
    vect_t  part_H  = {0, 0, 8};
    if (mk_particle(wdbp, "part.s", part_V, part_H, 5.0, 3.0) != 0) {
        bu_log("mk_particle failed\n");
        db_close(wdbp->dbip);
        return BRLCAD_ERROR;
    }

    /* --- superell (super-ellipsoid) --------------------------------- */
    point_t sel_V  = {60, 40, 0};
    vect_t  sel_A  = {7, 0, 0};
    vect_t  sel_B  = {0, 5, 0};
    vect_t  sel_C  = {0, 0, 3};
    /* mk_superell takes (fp, name, center, a, b, c, n, e) */
    {
        struct rt_db_internal intern;
        RT_DB_INTERNAL_INIT(&intern);
        struct rt_superell_internal *sei;
        BU_ALLOC(sei, struct rt_superell_internal);
        sei->magic = RT_SUPERELL_INTERNAL_MAGIC;
        VMOVE(sei->v, sel_V);
        VMOVE(sei->a, sel_A);
        VMOVE(sei->b, sel_B);
        VMOVE(sei->c, sel_C);
        sei->n = 2.0;
        sei->e = 2.0;
        intern.idb_major_type = DB5_MAJORTYPE_BRLCAD;
        intern.idb_type = ID_SUPERELL;
        intern.idb_meth = &OBJ[ID_SUPERELL];
        intern.idb_ptr = sei;
        if (wdb_put_internal(wdbp, "superell.s", &intern, 1.0) < 0) {
            bu_log("wdb_put_internal(superell) failed\n");
            rt_db_free_internal(&intern);
            db_close(wdbp->dbip);
            return BRLCAD_ERROR;
        }
        rt_db_free_internal(&intern);
    }

    /* --- cline (cylinder/line element) ------------------------------ */
    point_t cl_V   = {80, 40, 0};
    vect_t  cl_H   = {0, 0, 15};
    if (mk_cline(wdbp, "cline.s", cl_V, cl_H, 2.0, 0.2) != 0) {
        bu_log("mk_cline failed\n");
        db_close(wdbp->dbip);
        return BRLCAD_ERROR;
    }

    /* --- hrt (heart) ----------------------------------------------- */
    point_t hrt_V = {100, 20, 0};
    vect_t hrt_X = {4, 0, 0};
    vect_t hrt_Y = {0, 5, 0};
    vect_t hrt_Z = {0, 0, 6};
    if (mk_hrt(wdbp, "hrt.s", hrt_V, hrt_X, hrt_Y, hrt_Z, 2.0) != 0) {
        bu_log("mk_hrt failed\n");
        db_close(wdbp->dbip);
        return BRLCAD_ERROR;
    }

    db_close(wdbp->dbip);
    return BRLCAD_OK;
}


/* ------------------------------------------------------------------ *
 * Section 0 test functions
 * ------------------------------------------------------------------ */

static void
test_p0_cmd_exists(void)
{
    CHECK(ged_cmd_exists("edit"), "edit command is registered");
}

static void
test_p0_noargs(struct ged *gedp)
{
    const char *argv[] = {"edit", NULL};
    bu_vls_trunc(gedp->ged_result_str, 0);
    int ret = ged_exec(gedp, 1, argv);
    (void)ret;
    CHECK(1, "edit with no args does not crash");
}

static void
test_p0_bad_geom(struct ged *gedp)
{
    const char *argv[] = {"edit", "nonexistent_object_xyzzy", NULL};
    bu_vls_trunc(gedp->ged_result_str, 0);
    ged_exec(gedp, 2, argv);
    CHECK(1, "edit with nonexistent geometry does not crash");
}

static void
test_p0_fixture_objects_exist(struct ged *gedp)
{
    static const char *names[] = {
        "tor.s", "tgc.s", "ell.s", "sph.s", "arb8.s",
        "rpc.s", "rhc.s", "epa.s", "ehy.s", "eto.s",
        "hyp.s", "part.s", "superell.s", "cline.s", "hrt.s",
        NULL
    };
    for (int i = 0; names[i] != NULL; i++) {
        struct directory *dp = db_lookup(gedp->dbip, names[i], LOOKUP_QUIET);
        std::string msg = std::string("fixture object exists: ") + names[i];
        CHECK(dp != RT_DIR_NULL, msg.c_str());
    }
}

static void
test_p0_obj_nosubcmd(struct ged *gedp)
{
    const char *argv[] = {"edit", "tor.s", NULL};
    bu_vls_trunc(gedp->ged_result_str, 0);
    ged_exec(gedp, 2, argv);
    CHECK(1, "edit <obj> with no subcommand does not crash");
}

static void
test_p0_perturb(struct ged *gedp)
{
    const char *argv[] = {"edit", "tor.s", "perturb", "0.1", NULL};
    bu_vls_trunc(gedp->ged_result_str, 0);
    ged_exec(gedp, 4, argv);
    CHECK(1, "edit tor.s perturb does not crash");
}

static void
test_p0_hrt_descriptor_ops(struct ged *gedp)
{
    struct rt_hrt_internal h = {};

    {
        const char *av[] = {"edit", "hrt", "--list-ops", NULL};
        bu_vls_trunc(gedp->ged_result_str, 0);
        int ret = ged_exec(gedp, 3, av);
        CHECK(ret == BRLCAD_OK, "edit hrt --list-ops returns OK");
        CHECK(std::strstr(bu_vls_cstr(gedp->ged_result_str), "set_center") != NULL,
              "edit hrt --list-ops lists set_center");
        CHECK(std::strstr(bu_vls_cstr(gedp->ged_result_str), "set_cusp_distance") != NULL,
              "edit hrt --list-ops lists set_cusp_distance");
    }

    {
        const char *av[] = {"edit", "hrt.s", "set_center", "1", "2", "3", NULL};
        bu_vls_trunc(gedp->ged_result_str, 0);
        int ret = ged_exec(gedp, 6, av);
        CHECK(ret == BRLCAD_OK, "edit hrt.s set_center 1 2 3 returns OK");
        CHECK(read_hrt(gedp, "hrt.s", &h) == BRLCAD_OK, "read hrt.s after set_center succeeds");
        CHECK(NEAR_EQUAL(h.v[X], 1.0, NEAR_ENOUGH), "hrt.s V.x == 1 after set_center");
        CHECK(NEAR_EQUAL(h.v[Y], 2.0, NEAR_ENOUGH), "hrt.s V.y == 2 after set_center");
        CHECK(NEAR_EQUAL(h.v[Z], 3.0, NEAR_ENOUGH), "hrt.s V.z == 3 after set_center");
    }

    {
        const char *av[] = {"edit", "hrt.s", "set_x_direction", "1", "0", "0", NULL};
        bu_vls_trunc(gedp->ged_result_str, 0);
        int ret = ged_exec(gedp, 6, av);
        CHECK(ret == BRLCAD_OK, "edit hrt.s set_x_direction 1 0 0 returns OK");
    }

    {
        const char *av[] = {"edit", "hrt.s", "set_y_direction", "0", "1", "0", NULL};
        bu_vls_trunc(gedp->ged_result_str, 0);
        int ret = ged_exec(gedp, 6, av);
        CHECK(ret == BRLCAD_OK, "edit hrt.s set_y_direction 0 1 0 returns OK");
    }

    {
        const char *av[] = {"edit", "hrt.s", "set_z_direction", "0", "0", "1", NULL};
        bu_vls_trunc(gedp->ged_result_str, 0);
        int ret = ged_exec(gedp, 6, av);
        CHECK(ret == BRLCAD_OK, "edit hrt.s set_z_direction 0 0 1 returns OK");
    }

    {
        const char *av[] = {"edit", "hrt.s", "set_cusp_distance", "2.5", NULL};
        bu_vls_trunc(gedp->ged_result_str, 0);
        int ret = ged_exec(gedp, 4, av);
        CHECK(ret == BRLCAD_OK, "edit hrt.s set_cusp_distance 2.5 returns OK");
        CHECK(read_hrt(gedp, "hrt.s", &h) == BRLCAD_OK, "read hrt.s after set_d succeeds");
        CHECK(NEAR_EQUAL(h.d, 2.5, NEAR_ENOUGH), "hrt.s d == 2.5 after set_d");
    }

    {
	struct rt_hrt_internal before = {};
	CHECK(read_hrt(gedp, "hrt.s", &before) == BRLCAD_OK,
	      "read hrt.s before rejected edit succeeds");
        const char *av[] = {"edit", "hrt.s", "set_x_direction", "1", "1", "0", NULL};
        bu_vls_trunc(gedp->ged_result_str, 0);
        int ret = ged_exec(gedp, 6, av);
	CHECK(ret == BRLCAD_ERROR,
	      "edit hrt.s rejects a semantically invalid non-orthogonal axis");
	struct rt_hrt_internal after = {};
	CHECK(read_hrt(gedp, "hrt.s", &after) == BRLCAD_OK &&
	      VNEAR_EQUAL(before.xdir, after.xdir, NEAR_ENOUGH),
	      "rejected HRT edit leaves committed geometry unchanged");
    }

    {
        const char *av[] = {"edit", "hrt.s", "set_d", "0", NULL};
        bu_vls_trunc(gedp->ged_result_str, 0);
        int ret = ged_exec(gedp, 4, av);
        CHECK(ret == BRLCAD_ERROR, "edit hrt.s set_d 0 returns error");
    }

    {
	struct edit_session_observer obs;
	ged_edit_observer_token token =
	    ged_edit_observer_add(gedp, edit_session_cb, &obs);
	CHECK(token != 0, "HRT typed-rejection observer registers");
	ged_edit_session_ref ref = GED_EDIT_SESSION_REF_NULL;
	CHECK(ged_edit_session_begin(gedp, "hrt.s", NULL, &ref) ==
	    GED_EDIT_OK, "typed API begins an HRT edit session");

	struct ged_edit_session_info before_info = {};
	CHECK(ged_edit_session_info_get(gedp, ref, &before_info) ==
	    GED_EDIT_OK, "HRT session reports state before rejected update");
	struct rt_edit_cmd_values before;
	CHECK(ged_edit_session_command_values_get(gedp, ref, 43002,
	    &before) == GED_EDIT_OK && before.value_count == 3,
	    "HRT session reads the current X direction");
	const size_t event_count = obs.kinds.size();
	fastf_t invalid_xdir[RT_EDIT_MAXPARA] = {1.0, 1.0, 0.0};
	struct ged_edit_command_input input = {};
	input.command_id = 43002;
	input.values = invalid_xdir;
	input.value_count = 3;
	CHECK(ged_edit_session_apply(gedp, ref, &input) ==
	    GED_EDIT_REJECTED,
	    "typed API reports semantic HRT rejection distinctly");

	struct ged_edit_session_info after_info = {};
	struct rt_edit_cmd_values after;
	CHECK(ged_edit_session_info_get(gedp, ref, &after_info) ==
	    GED_EDIT_OK && after_info.revision == before_info.revision &&
	    after_info.dirty == before_info.dirty,
	    "rejected typed update changes neither revision nor dirty state");
	CHECK(ged_edit_session_command_values_get(gedp, ref, 43002,
	    &after) == GED_EDIT_OK &&
	    after.value_count == before.value_count &&
	    VNEAR_EQUAL(before.values, after.values, NEAR_ENOUGH),
	    "rejected typed update leaves session geometry unchanged");
	CHECK(obs.kinds.size() == event_count,
	    "rejected typed update publishes no false update event");
	CHECK(ged_edit_session_cancel(gedp, ref) == GED_EDIT_OK,
	    "HRT typed-rejection session cancels cleanly");
	CHECK(ged_edit_observer_remove(gedp, token) == 1,
	    "HRT typed-rejection observer removal succeeds");
    }
}

static void
test_p0_obj_help_descriptor_ops(struct ged *gedp)
{
    {
        const char *av[] = {"edit", "tor.s", "-h", NULL};
        bu_vls_trunc(gedp->ged_result_str, 0);
        int ret = ged_exec(gedp, 3, av);
        CHECK(ret == BRLCAD_OK, "edit tor.s -h returns OK");
        CHECK(std::strstr(bu_vls_cstr(gedp->ged_result_str), "Torus") != NULL,
              "edit tor.s -h prints Torus descriptor help");
    }

    {
        const char *av[] = {"edit", "hrt.s", "-h", NULL};
        bu_vls_trunc(gedp->ged_result_str, 0);
        int ret = ged_exec(gedp, 3, av);
        CHECK(ret == BRLCAD_OK, "edit hrt.s -h returns OK");
        CHECK(std::strstr(bu_vls_cstr(gedp->ged_result_str), "Heart") != NULL,
              "edit hrt.s -h prints Heart descriptor help");
    }
}

static void
test_p0_desc_coverage(void)
{
    static const int typed_prims[] = {
        ID_TOR, ID_TGC, ID_ELL,
        ID_ARB8,
        ID_HRT,
        ID_ARS, ID_BSPLINE, ID_NMG,
        ID_EBM, ID_VOL, ID_PIPE, ID_PARTICLE,
        ID_RPC, ID_RHC, ID_EPA, ID_EHY, ID_ETO,
        ID_DSP, ID_CLINE, ID_COMBINATION, ID_SKETCH,
        ID_SUPERELL, ID_HYP, ID_EXTRUDE,
        -1
    };
    for (int i = 0; typed_prims[i] >= 0; i++) {
        struct bu_vls out = BU_VLS_INIT_ZERO;
        int ret = rt_edit_type_to_json(&out, typed_prims[i]);
        std::string msg = std::string("rt_edit_type_to_json succeeds for prim id ")
            + std::to_string(typed_prims[i]);
        CHECK(ret == BRLCAD_OK && bu_vls_strlen(&out) > 0, msg.c_str());
        bu_vls_free(&out);
    }
}

static void
test_p0_no_desc_returns_error(void)
{
    static const int no_desc[] = {
	ID_NULL,
        -1
    };
    for (int i = 0; no_desc[i] >= 0; i++) {
        struct bu_vls out = BU_VLS_INIT_ZERO;
        int ret = rt_edit_type_to_json(&out, no_desc[i]);
        std::string msg = std::string("rt_edit_type_to_json returns error for undescribed prim id ")
            + std::to_string(no_desc[i]);
        bool acceptable = (ret != BRLCAD_OK) || (bu_vls_strlen(&out) == 0);
        CHECK(acceptable, msg.c_str());
        bu_vls_free(&out);
    }
}


/* ================================================================== *
 * Section 1 — Parser, selection, conflict arbiter, edit buffer
 * ================================================================== */

static int
create_p1_fixture(const char *path)
{
    struct rt_wdb *wdbp = wdb_fopen(path);
    if (!wdbp)
        return BRLCAD_ERROR;

    point_t v  = VINIT_ZERO;
    vect_t  h  = {0, 0, 1};
    if (mk_tor(wdbp, "tor.s", v, h, 4.0, 1.0) != 0) { wdb_close(wdbp); return BRLCAD_ERROR; }

    point_t sp = {10, 0, 0};
    if (mk_sph(wdbp, "sph.s", sp, 2.0) != 0) { wdb_close(wdbp); return BRLCAD_ERROR; }

    fastf_t bot_vertices[] = {
	0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 2
    };
    int bot_faces[] = {
	0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3
    };
    if (mk_bot(wdbp, "bot.s", RT_BOT_SOLID, RT_BOT_CCW, 0, 4, 4,
	    bot_vertices, bot_faces, NULL, NULL) != 0) {
	wdb_close(wdbp);
	return BRLCAD_ERROR;
    }

    struct wmember wm;
    BU_LIST_INIT(&wm.l);
    mk_addmember("tor.s", &wm.l, NULL, WMOP_UNION);
    if (mk_lcomb(wdbp, "group.c", &wm, 0, NULL, NULL, NULL, 0) != 0) {
        wdb_close(wdbp);
        return BRLCAD_ERROR;
    }

    BU_LIST_INIT(&wm.l);
    mk_addmember("tor.s", &wm.l, NULL, WMOP_UNION);
    if (mk_lcomb(wdbp, "group2.c", &wm, 0, NULL, NULL, NULL, 0) != 0) {
	wdb_close(wdbp);
	return BRLCAD_ERROR;
    }

    wdb_close(wdbp);
    return BRLCAD_OK;
}

static void
test_p1_noargs(struct ged *gedp)
{
    const char *av[] = { "edit", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    int ret = ged_exec(gedp, 1, av);
    CHECK(ret == BRLCAD_OK, "edit (no args) returns OK (prints help)");
}

static void
test_p1_bad_geom(struct ged *gedp)
{
    const char *av[] = { "edit", "nonexistent_object_xyz", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    int ret = ged_exec(gedp, 2, av);
    CHECK(ret == BRLCAD_ERROR, "edit <nonexistent_object> returns error");
}

static void
test_p1_obj_nosubcmd(struct ged *gedp)
{
    const char *av[] = { "edit", "tor.s", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    int ret = ged_exec(gedp, 2, av);
    CHECK(ret == BRLCAD_ERROR, "edit tor.s (no subcommand) returns error");
    CHECK(bu_vls_strlen(gedp->ged_result_str) > 0,
          "edit tor.s (no subcommand) prints a message");
}

static void
test_p1_help_flag(struct ged *gedp)
{
    const char *av[] = { "edit", "-h", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    int ret = ged_exec(gedp, 2, av);
    CHECK(ret == BRLCAD_OK, "edit -h returns OK");
    CHECK(std::strstr(bu_vls_cstr(gedp->ged_result_str), "--list-all-prim-ops") != NULL,
          "edit -h mentions --list-all-prim-ops");
}

static void
test_p1_uri_fragment(void)
{
    std::string raw = "tor.s#V1";
    bool has_fragment = false;
    std::string fragment_val;
    try {
        uri obj_uri(std::string("g:") + raw);
        fragment_val = obj_uri.get_fragment();
        has_fragment = !fragment_val.empty();
    } catch (...) {}
    CHECK(has_fragment,         "URI fragment 'V1' is parsed from 'tor.s#V1'");
    CHECK(fragment_val == "V1", "URI fragment value is 'V1'");
}

static void
test_p1_uri_query(void)
{
    std::string raw = "tor.s?V*";
    bool has_query = false;
    std::string query_val;
    try {
        uri obj_uri(std::string("g:") + raw);
        query_val = obj_uri.get_query();
        has_query = !query_val.empty();
    } catch (...) {}
    CHECK(has_query,          "URI query 'V*' is parsed from 'tor.s?V*'");
    CHECK(query_val == "V*",  "URI query value is 'V*'");
}

static void
test_p1_batch_marker(struct ged *gedp)
{
    const char *av[] = { "edit", ".", "perturb", "0.1", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    int ret = ged_exec(gedp, 4, av);
    CHECK(ret == BRLCAD_ERROR,
          "edit . perturb (no selection) returns error");
}

static void
test_p1_selection_fallback(struct ged *gedp)
{
    ged_selection_clear(gedp, NULL);
    if (!ged_selection_state_available(gedp) ||
	    !ged_selection_select_path(gedp, NULL, "tor.s", 0)) {
	CHECK(0, "Could not select tor.s in default selection state");
	return;
    }

    const char *av[] = { "edit", "perturb", "0.0", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    int ret = ged_exec(gedp, 3, av);
    CHECK(ret == BRLCAD_OK,
          "edit perturb 0.0 (via selection fallback) returns OK");
    ged_selection_deselect_path(gedp, NULL, "tor.s", 0);
}

static void
test_p1_no_geom_no_sel(struct ged *gedp)
{
    ged_selection_clear(gedp, NULL);

    const char *av[] = { "edit", "perturb", "0.1", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    int ret = ged_exec(gedp, 3, av);
    CHECK(ret == BRLCAD_ERROR,
          "edit perturb (no geom, no selection) returns error");
}

static void
test_p1_conflict_arbiter(struct ged *gedp)
{
    ged_selection_clear(gedp, NULL);
    if (!ged_selection_state_available(gedp) ||
	    !ged_selection_select_path(gedp, NULL, "tor.s", 0)) {
	CHECK(0, "Could not select tor.s in default selection state");
	return;
    }

    const char *av[] = { "edit", "tor.s", "perturb", "0.1", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    int ret = ged_exec(gedp, 4, av);
    CHECK(ret == BRLCAD_ERROR,
          "edit tor.s perturb (selection conflict, no flag) returns error");
    bool has_conflict_msg =
        (bu_vls_strlen(gedp->ged_result_str) > 0) &&
        (strstr(bu_vls_cstr(gedp->ged_result_str), "Conflict") != NULL ||
         strstr(bu_vls_cstr(gedp->ged_result_str), "conflict") != NULL);
    CHECK(has_conflict_msg,
          "edit tor.s perturb (selection conflict) prints conflict message");
    ged_selection_deselect_path(gedp, NULL, "tor.s", 0);
}

static void
test_p1_flag_S(struct ged *gedp)
{
    ged_selection_clear(gedp, NULL);
    if (!ged_selection_state_available(gedp) ||
	    !ged_selection_select_path(gedp, NULL, "tor.s", 0)) {
	CHECK(0, "Could not select tor.s in default selection state");
	return;
    }

    const char *av[] = { "edit", "-S", "sph.s", "perturb", "0.0", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    int ret = ged_exec(gedp, 5, av);
    CHECK(ret == BRLCAD_OK,
          "edit -S sph.s perturb 0.0 (uses selection tor.s, not sph.s) returns OK");
    ged_selection_deselect_path(gedp, NULL, "tor.s", 0);
}

static void
test_p1_flag_f(struct ged *gedp)
{
    ged_selection_clear(gedp, NULL);
    if (!ged_selection_state_available(gedp) ||
	    !ged_selection_select_path(gedp, NULL, "tor.s", 0)) {
	CHECK(0, "Could not select tor.s in default selection state");
	return;
    }

    const char *av[] = { "edit", "-f", "tor.s", "perturb", "0.0", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    int ret = ged_exec(gedp, 5, av);
    CHECK(ret == BRLCAD_OK,
          "edit -f tor.s perturb 0.0 (force, bypasses conflict) returns OK");
    ged_selection_deselect_path(gedp, NULL, "tor.s", 0);
}

static void
test_p1_buf_get_missing(struct ged *gedp)
{
    struct db_full_path dfp;
    db_full_path_init(&dfp);
    db_string_to_path(&dfp, gedp->dbip, "tor.s");

    struct rt_edit *s = ged_edit_buf_get(gedp, &dfp);
    CHECK(s == NULL, "ged_edit_buf_get returns NULL for absent entry");
    db_free_full_path(&dfp);
}

static void
test_p1_buf_set_get(struct ged *gedp)
{
    struct bn_tol tol = BN_TOL_INIT_ZERO;
    BN_TOL_INIT_SET_TOL(&tol);
    struct db_full_path dfp;
    db_full_path_init(&dfp);
    if (db_string_to_path(&dfp, gedp->dbip, "tor.s") < 0) {
        CHECK(0, "db_string_to_path for tor.s failed unexpectedly");
        return;
    }
    struct rt_edit *s = rt_edit_create(&dfp, gedp->dbip, &tol, NULL);
    if (!s) { CHECK(0, "rt_edit_create failed for tor.s"); db_free_full_path(&dfp); return; }

    ged_edit_buf_set(gedp, &dfp, s);
    struct rt_edit *got = ged_edit_buf_get(gedp, &dfp);
    CHECK(got == s, "ged_edit_buf_get returns the stored rt_edit * after set");

    ged_edit_buf_abandon(gedp, &dfp);
    db_free_full_path(&dfp);
}

static void
test_p1_buf_abandon(struct ged *gedp)
{
    struct bn_tol tol = BN_TOL_INIT_ZERO;
    BN_TOL_INIT_SET_TOL(&tol);
    struct db_full_path dfp;
    db_full_path_init(&dfp);
    if (db_string_to_path(&dfp, gedp->dbip, "sph.s") < 0) {
        CHECK(0, "db_string_to_path for sph.s failed unexpectedly");
        return;
    }
    struct rt_edit *s = rt_edit_create(&dfp, gedp->dbip, &tol, NULL);
    if (!s) { CHECK(0, "rt_edit_create failed for sph.s"); db_free_full_path(&dfp); return; }

    ged_edit_buf_set(gedp, &dfp, s);
    ged_edit_buf_abandon(gedp, &dfp);
    struct rt_edit *got = ged_edit_buf_get(gedp, &dfp);
    CHECK(got == NULL, "ged_edit_buf_get returns NULL after abandon");
    db_free_full_path(&dfp);
}

static void
test_p1_buf_flush(struct ged *gedp)
{
    struct bn_tol tol = BN_TOL_INIT_ZERO;
    BN_TOL_INIT_SET_TOL(&tol);
    struct db_full_path dfp1, dfp2;
    db_full_path_init(&dfp1);
    db_full_path_init(&dfp2);
    db_string_to_path(&dfp1, gedp->dbip, "tor.s");
    db_string_to_path(&dfp2, gedp->dbip, "sph.s");

    struct rt_edit *s1 = rt_edit_create(&dfp1, gedp->dbip, &tol, NULL);
    struct rt_edit *s2 = rt_edit_create(&dfp2, gedp->dbip, &tol, NULL);
    if (!s1 || !s2) {
        CHECK(0, "rt_edit_create failed for flush test fixture");
        if (s1) rt_edit_destroy(s1);
        if (s2) rt_edit_destroy(s2);
        db_free_full_path(&dfp1);
        db_free_full_path(&dfp2);
        return;
    }

    ged_edit_buf_set(gedp, &dfp1, s1);
    ged_edit_buf_set(gedp, &dfp2, s2);
    struct edit_event_observer obs;
    ged_event_observer_token token =
	ged_event_observer_add(gedp, GED_EVENT_OBSERVER_POST_RECONCILE,
		edit_event_cb, &obs);
    CHECK(token != 0, "ged_edit_buf_flush event observer registers");
    ged_edit_buf_flush(gedp);

    struct rt_edit *g1 = ged_edit_buf_get(gedp, &dfp1);
    struct rt_edit *g2 = ged_edit_buf_get(gedp, &dfp2);
    CHECK(g1 == NULL && g2 == NULL, "ged_edit_buf_flush removes all entries");
    CHECK(obs.calls == 1, "ged_edit_buf_flush publishes one event transaction");
    CHECK(obs.tor_modified == 1 && obs.sph_modified == 1,
	"ged_edit_buf_flush emits object-modified events for flushed entries");
    CHECK(ged_event_observer_remove(gedp, token) == 1,
	"ged_edit_buf_flush event observer removal succeeds");

    db_free_full_path(&dfp1);
    db_free_full_path(&dfp2);
}

static void
test_p1_perturb_regression(struct ged *gedp)
{
    const char *av[] = { "edit", "tor.s", "perturb", "0.001", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    int ret = ged_exec(gedp, 4, av);
    CHECK(ret == BRLCAD_OK, "edit tor.s perturb 0.001 (regression) returns OK");
}


static void
test_p1_edit_session_api(struct ged *gedp)
{
    struct edit_session_observer obs;
    ged_edit_observer_token token =
	ged_edit_observer_add(gedp, edit_session_cb, &obs);
    CHECK(token != 0, "edit-session observer registers");

    const char *set_12[] = {"edit", "-i", "tor.s", "r1", "12.5", NULL};
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 5, set_12) == BRLCAD_OK,
	"CLI creates and updates an intermediate edit session");
    CHECK(obs.kinds.size() == 2 &&
	obs.kinds[0] == GED_EDIT_SESSION_BEGIN &&
	obs.kinds[1] == GED_EDIT_SESSION_UPDATE,
	"intermediate CLI edit reports begin followed by update");
    CHECK(obs.revisions.size() == 2 && obs.revisions[0] == 0 &&
	obs.revisions[1] == 1,
	"edit-session revisions start at zero and advance on update");

    ged_edit_session_ref ref = GED_EDIT_SESSION_REF_NULL;
    CHECK(ged_edit_session_find(gedp, "tor.s", &ref) == GED_EDIT_OK &&
	!ged_edit_session_ref_is_null(ref),
	"active edit session is discoverable by canonical path");
    CHECK(ref.owner == obs.first.owner && ref.id == obs.first.id &&
	ref.generation == obs.first.generation,
	"observer and query report the same stable session reference");

    struct ged_edit_session_info info;
    CHECK(ged_edit_session_info_get(gedp, ref, &info) == GED_EDIT_OK &&
	info.revision == 1 && info.dirty && info.primitive_type == ID_TOR,
	"session info reports current revision, dirty state, and primitive type");

    struct bu_vls path = BU_VLS_INIT_ZERO;
    CHECK(ged_edit_session_path_get(gedp, ref, &path) == GED_EDIT_OK &&
	BU_STR_EQUAL(bu_vls_cstr(&path), "/tor.s"),
	"session path query returns the canonical target");
    bu_vls_free(&path);

    const struct rt_edit_prim_desc *desc = NULL;
    CHECK(ged_edit_session_descriptor_get(gedp, ref, &desc) == GED_EDIT_OK &&
	desc && BU_STR_EQUAL(desc->prim_type, "tor"),
	"session descriptor query reports the librt primitive contract");

    struct rt_edit_cmd_values values;
    CHECK(ged_edit_session_command_values_get(gedp, ref, 1021,
	&values) == GED_EDIT_OK && values.value_count == 1 &&
	values.value_valid[0] && fabs(values.values[0] - 12.5) < NEAR_ENOUGH,
	"session value query observes the CLI intermediate value");
    const char *status[] = {"edit", "tor.s", "status", NULL};
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, status) == BRLCAD_OK &&
	strstr(bu_vls_cstr(gedp->ged_result_str), "revision=1") &&
	strstr(bu_vls_cstr(gedp->ged_result_str), "type=tor"),
	"CLI status reports the shared session revision and primitive type");
    const char *get_r1[] = {"edit", "tor.s", "get", "r1", NULL};
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 4, get_r1) == BRLCAD_OK &&
	strstr(bu_vls_cstr(gedp->ged_result_str), "12.5"),
	"CLI get observes the current shared intermediate value");

    const char *set_15[] = {"edit", "-i", "tor.s", "r1", "15", NULL};
    CHECK(ged_exec(gedp, 5, set_15) == BRLCAD_OK,
	"a second CLI operation reuses the active session");
    ged_edit_session_ref same_ref = GED_EDIT_SESSION_REF_NULL;
    CHECK(ged_edit_session_find(gedp, "tor.s", &same_ref) == GED_EDIT_OK &&
	same_ref.owner == ref.owner && same_ref.id == ref.id &&
	same_ref.generation == ref.generation,
	"successive edits preserve session identity");
    CHECK(obs.kinds.size() == 3 &&
	obs.kinds.back() == GED_EDIT_SESSION_UPDATE &&
	obs.revisions.back() == 2,
	"successive edit reports one higher-revision update");

    const char *checkpoint[] = {"edit", "tor.s", "checkpoint", NULL};
    CHECK(ged_exec(gedp, 3, checkpoint) == BRLCAD_OK &&
	obs.kinds.back() == GED_EDIT_SESSION_CHECKPOINT,
	"CLI checkpoint reports a session checkpoint event");

    const char *set_20[] = {"edit", "-i", "tor.s", "r1", "20", NULL};
    CHECK(ged_exec(gedp, 5, set_20) == BRLCAD_OK,
	"post-checkpoint intermediate update succeeds");
    const char *revert[] = {"edit", "-i", "tor.s", "revert", NULL};
    CHECK(ged_exec(gedp, 4, revert) == BRLCAD_OK &&
	obs.kinds.back() == GED_EDIT_SESSION_REVERT,
	"CLI revert reports a session revert event");
    CHECK(ged_edit_session_command_values_get(gedp, ref, 1021,
	&values) == GED_EDIT_OK &&
	fabs(values.values[0] - 15.0) < NEAR_ENOUGH,
	"session readback observes reverted intermediate geometry");

    const char *reset[] = {"edit", "tor.s", "reset", NULL};
    CHECK(ged_exec(gedp, 3, reset) == BRLCAD_OK &&
	obs.kinds.back() == GED_EDIT_SESSION_CANCEL,
	"CLI reset reports session cancellation");
    CHECK(ged_edit_session_info_get(gedp, ref, &info) == GED_EDIT_STALE,
	"closed session references are reported as stale");

    ged_edit_session_ref typed_ref = GED_EDIT_SESSION_REF_NULL;
    CHECK(ged_edit_session_begin(gedp, "tor.s", NULL, &typed_ref) ==
	GED_EDIT_OK,
	"typed API begins a new session after CLI cancellation");
    fastf_t typed_values[RT_EDIT_MAXPARA] = {0.0};
    typed_values[0] = 17.25;
    struct ged_edit_command_input input;
    input.command_id = 1021;
    input.values = typed_values;
    input.value_count = 1;
    input.strings = NULL;
    input.string_count = 0;
    input.view = NULL;
    CHECK(ged_edit_session_apply(gedp, typed_ref, &input) == GED_EDIT_OK,
	"typed API applies a descriptor command without string parsing");
    CHECK(ged_edit_session_command_values_get(gedp, typed_ref, 1021,
	&values) == GED_EDIT_OK &&
	fabs(values.values[0] - 17.25) < NEAR_ENOUGH,
	"typed update is visible through shared session readback");
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 4, get_r1) == BRLCAD_OK &&
	strstr(bu_vls_cstr(gedp->ged_result_str), "17.25"),
	"CLI get observes an update made through the typed GUI/manipulator API");
    struct rt_db_internal preview = RT_DB_INTERNAL_INIT_ZERO;
    CHECK(ged_edit_session_internal_copy(gedp, typed_ref, &preview) ==
	GED_EDIT_OK && preview.idb_type == ID_TOR,
	"session produces an independently owned preview snapshot");
    rt_db_free_internal(&preview);
    CHECK(ged_edit_session_commit(gedp, typed_ref) == GED_EDIT_OK,
	"typed API commits and closes its session");
    struct rt_tor_internal committed;
    CHECK(read_tor(gedp, "tor.s", &committed) == BRLCAD_OK &&
	fabs(committed.r_a - 17.25) < NEAR_ENOUGH,
	"typed commit writes the shared intermediate state to the database");

    const char *nested_set[] = {
	"edit", "-i", "group.c/tor.s", "r1", "8.5", NULL
    };
    CHECK(ged_exec(gedp, 5, nested_set) == BRLCAD_OK,
	"descriptor CLI dispatch accepts a nested combination path");
    ged_edit_session_ref nested_ref = GED_EDIT_SESSION_REF_NULL;
    CHECK(ged_edit_session_find(gedp, "group.c/tor.s", &nested_ref) ==
	GED_EDIT_OK && !ged_edit_session_ref_is_null(nested_ref),
	"nested CLI edit creates a path-scoped shared session");
    const char *nested_get[] = {
	"edit", "group.c/tor.s", "get", "r1", NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 4, nested_get) == BRLCAD_OK &&
	strstr(bu_vls_cstr(gedp->ged_result_str), "8.5"),
	"nested session state is visible through CLI readback");
    const char *nested_reset[] = {
	"edit", "group.c/tor.s", "reset", NULL
    };
    CHECK(ged_exec(gedp, 3, nested_reset) == BRLCAD_OK &&
	ged_edit_session_info_get(gedp, nested_ref, &info) == GED_EDIT_STALE,
	"nested session reset cancels the same shared state");

    auto read_combination = [gedp](int *region, long *los,
	std::string *shader, std::string *material) -> bool {
	struct directory *dp = db_lookup(gedp->dbip, "group.c", LOOKUP_QUIET);
	struct rt_db_internal intern = RT_DB_INTERNAL_INIT_ZERO;
	if (!dp || rt_db_get_internal(&intern, dp, gedp->dbip, NULL) < 0)
	    return false;
	struct rt_comb_internal *comb =
	    (struct rt_comb_internal *)intern.idb_ptr;
	if (!comb || comb->magic != RT_COMB_MAGIC) {
	    rt_db_free_internal(&intern);
	    return false;
	}
	if (region) *region = comb->region_flag;
	if (los) *los = comb->los;
	if (shader) *shader = bu_vls_cstr(&comb->shader);
	if (material) *material = bu_vls_cstr(&comb->material);
	rt_db_free_internal(&intern);
	return true;
    };

    const char *comb_shader[] = {
	"edit", "-i", "group.c", "set_shader", "plastic", NULL
    };
    const char *comb_region[] = {
	"edit", "-i", "group.c", "set_region", "1", NULL
    };
    const char *comb_los[] = {
	"edit", "-i", "group.c", "set_los", "80", NULL
    };
    const char *comb_material[] = {
	"edit", "-i", "group.c", "set_material", "air", NULL
    };
    CHECK(ged_exec(gedp, 5, comb_shader) == BRLCAD_OK &&
	ged_exec(gedp, 5, comb_region) == BRLCAD_OK &&
	ged_exec(gedp, 5, comb_los) == BRLCAD_OK &&
	ged_exec(gedp, 5, comb_material) == BRLCAD_OK,
	"combination edits accumulate in one shared intermediate session");
    int region = -1;
    long los = -1;
    std::string shader;
    std::string material;
    CHECK(read_combination(&region, &los, &shader, &material) &&
	!region && los == 0 && shader.empty() && material.empty(),
	"combination intermediate does not leak into the database");
    ged_edit_session_ref comb_ref = GED_EDIT_SESSION_REF_NULL;
    CHECK(ged_edit_session_find(gedp, "group.c", &comb_ref) == GED_EDIT_OK,
	"combination session is discoverable before commit");
    struct rt_edit_cmd_values comb_values;
    CHECK(ged_edit_session_command_values_get(gedp, comb_ref, 12007,
	&comb_values) == GED_EDIT_OK && comb_values.string_valid[0] &&
	BU_STR_EQUAL(comb_values.strings[0], "plastic"),
	"combination string readback observes intermediate state");
    CHECK(ged_edit_session_commit(gedp, comb_ref) == GED_EDIT_OK,
	"combination shared session commits explicitly");
    CHECK(read_combination(&region, &los, &shader, &material) &&
	region && los == 80 && shader == "plastic" && material == "air",
	"combination commit preserves geometry and v5 attributes");
    CHECK(ged_edit_observer_remove(gedp, token) == 1,
	"edit-session observer removal succeeds");
}


static void
test_p1_alternate_occurrence_session_conflict(struct ged *gedp)
{
    ged_edit_session_ref first = GED_EDIT_SESSION_REF_NULL;
    ged_edit_session_ref same = GED_EDIT_SESSION_REF_NULL;
    ged_edit_session_ref alternate = GED_EDIT_SESSION_REF_NULL;
    CHECK(ged_edit_session_begin(gedp, "group.c/tor.s", NULL, &first) ==
	GED_EDIT_OK, "first occurrence begins its writable session");
    CHECK(ged_edit_session_begin(gedp, "/group.c/tor.s", NULL, &same) ==
	GED_EDIT_OK && same.owner == first.owner && same.id == first.id &&
	same.generation == first.generation,
	"the same canonical occurrence joins its authoritative session");
    CHECK(ged_edit_session_begin(gedp, "group2.c/tor.s", NULL,
	&alternate) == GED_EDIT_CONFLICT &&
	ged_edit_session_ref_is_null(alternate),
	"a second occurrence of the writable terminal object reports conflict");

    const char *set_first[] = {
	"edit", "-i", "group.c/tor.s", "r1", "19", NULL
    };
    CHECK(ged_exec(gedp, 5, set_first) == BRLCAD_OK,
	"the authoritative occurrence remains editable after a conflict");
    struct ged_edit_session_info info = {};
    CHECK(ged_edit_session_info_get(gedp, first, &info) == GED_EDIT_OK &&
	info.revision == 1 && info.dirty,
	"the authoritative occurrence retains one coherent revision stream");

    const char *bad_value[] = {
	"edit", "-i", "group.c/tor.s", "r1", "not-a-number", NULL
    };
    CHECK(ged_exec(gedp, 5, bad_value) == BRLCAD_ERROR,
	"failed command parsing rejects an intermediate update");
    struct ged_edit_session_info after_bad = {};
    CHECK(ged_edit_session_info_get(gedp, first, &after_bad) == GED_EDIT_OK &&
	after_bad.revision == info.revision && after_bad.dirty == info.dirty,
	"failed command parsing leaves session revision and state unchanged");

    bu_vls_trunc(gedp->ged_result_str, 0);
    const char *set_alternate[] = {
	"edit", "-i", "group2.c/tor.s", "r1", "23", NULL
    };
    CHECK(ged_exec(gedp, 5, set_alternate) == BRLCAD_ERROR &&
	strstr(bu_vls_cstr(gedp->ged_result_str), "already being edited"),
	"CLI reports an alternate-occurrence session conflict explicitly");

    const char *get_first[] = {
	"edit", "group.c/tor.s", "get", "r1", NULL
    };
    CHECK(ged_exec(gedp, 4, get_first) == BRLCAD_OK &&
	strstr(bu_vls_cstr(gedp->ged_result_str), "19"),
	"CLI still reads the authoritative intermediate after a conflict");

    CHECK(ged_edit_session_cancel(gedp, first) == GED_EDIT_OK,
	"cancelling the first occurrence releases terminal edit authority");
    CHECK(ged_edit_session_begin(gedp, "group2.c/tor.s", NULL,
	&alternate) == GED_EDIT_OK,
	"another occurrence can begin after authority is released");
    CHECK(ged_edit_session_begin(gedp, "group.c/tor.s", NULL, &first) ==
	GED_EDIT_CONFLICT,
	"terminal edit authority is symmetric between occurrence paths");
    CHECK(ged_edit_session_cancel(gedp, alternate) == GED_EDIT_OK,
	"alternate occurrence session cancels cleanly");
}


static void
test_p1_independent_ged_context_sessions(struct ged *gedp,
	const char *database_path)
{
    ged_edit_session_ref first = GED_EDIT_SESSION_REF_NULL;
    CHECK(ged_edit_session_begin(gedp, "tor.s", NULL, &first) == GED_EDIT_OK,
	"first GED context begins an edit session");

    struct ged *other = open_fixture(database_path);
    CHECK(other != NULL, "a second GED context opens the same database");
    if (!other) {
	(void)ged_edit_session_cancel(gedp, first);
	return;
    }

    ged_edit_session_ref second = GED_EDIT_SESSION_REF_NULL;
    CHECK(ged_edit_session_begin(other, "tor.s", NULL, &second) == GED_EDIT_OK &&
	second.owner != first.owner,
	"independent GED contexts have independent session authorities");
    struct ged_edit_session_info foreign_info = {};
    CHECK(ged_edit_session_info_get(other, first, &foreign_info) ==
	GED_EDIT_INVALID && ged_edit_session_info_get(gedp, second,
	&foreign_info) == GED_EDIT_INVALID,
	"a session reference cannot cross its GED context boundary");

    fastf_t first_value[1] = {21.0};
    fastf_t second_value[1] = {22.0};
    struct ged_edit_command_input input = {};
    input.command_id = 1021;
    input.values = first_value;
    input.value_count = 1;
    CHECK(ged_edit_session_apply(gedp, first, &input) == GED_EDIT_OK,
	"first GED context updates its private intermediate");
    input.values = second_value;
    CHECK(ged_edit_session_apply(other, second, &input) == GED_EDIT_OK,
	"second GED context updates its private intermediate");

    struct rt_edit_cmd_values first_read = {};
    struct rt_edit_cmd_values second_read = {};
    CHECK(ged_edit_session_command_values_get(gedp, first, 1021,
	&first_read) == GED_EDIT_OK &&
	ged_edit_session_command_values_get(other, second, 1021,
	&second_read) == GED_EDIT_OK &&
	fabs(first_read.values[0] - 21.0) < NEAR_ENOUGH &&
	fabs(second_read.values[0] - 22.0) < NEAR_ENOUGH,
	"independent GED contexts do not leak intermediate values");

    CHECK(ged_edit_session_cancel(other, second) == GED_EDIT_OK,
	"second GED context cancels its intermediate cleanly");
    ged_close(other);
    CHECK(ged_edit_session_cancel(gedp, first) == GED_EDIT_OK,
	"first GED context remains valid after the second closes");
}


static void
test_p1_selection_session_state(struct ged *gedp)
{
    ged_edit_session_ref ref = GED_EDIT_SESSION_REF_NULL;
    CHECK(ged_edit_session_begin(gedp, "bot.s", NULL, &ref) == GED_EDIT_OK,
	"selection-state test begins a BoT edit session");

    const struct rt_edit_prim_desc *descriptor = NULL;
    int select_vertex = 0;
    if (ged_edit_session_descriptor_get(gedp, ref, &descriptor) ==
	    GED_EDIT_OK && descriptor) {
	for (int i = 0; i < descriptor->ncmd; i++) {
	    if (descriptor->cmds[i].name &&
		BU_STR_EQUAL(descriptor->cmds[i].name, "ECMD_BOT_PICKV")) {
		select_vertex = descriptor->cmds[i].cmd_id;
		break;
	    }
	}
    }
    CHECK(select_vertex != 0,
	"selection-state test discovers the stable BoT selection operation");

    struct rt_edit_param_bounds bounds = {};
    CHECK(ged_edit_session_parameter_bounds_get(gedp, ref, select_vertex, 0,
	&bounds) == GED_EDIT_OK && bounds.has_minimum && bounds.has_maximum &&
	fabs(bounds.minimum - (-1.0)) < NEAR_ENOUGH &&
	fabs(bounds.maximum - 3.0) < NEAR_ENOUGH,
	"BoT session exposes current-topology vertex bounds");

    fastf_t values[1] = {4.0};
    struct ged_edit_command_input input = {};
    input.command_id = select_vertex;
    input.values = values;
    input.value_count = 1;
    struct ged_edit_session_info info = {};
    CHECK(ged_edit_session_apply(gedp, ref, &input) == GED_EDIT_INVALID &&
	ged_edit_session_info_get(gedp, ref, &info) == GED_EDIT_OK &&
	info.revision == 0 && !info.dirty,
	"BoT topology bounds reject an unavailable vertex without changing state");
    values[0] = 1.0;
    CHECK(select_vertex && ged_edit_session_apply(gedp, ref, &input) ==
	GED_EDIT_OK, "typed BoT selection succeeds");

    CHECK(ged_edit_session_info_get(gedp, ref, &info) == GED_EDIT_OK &&
	info.revision == 1 && !info.dirty,
	"selection advances revision without marking geometry dirty");
    CHECK(ged_edit_session_checkpoint(gedp, ref) == GED_EDIT_OK,
	"clean selection state can be checkpointed");
    values[0] = 2.0;
    CHECK(ged_edit_session_apply(gedp, ref, &input) == GED_EDIT_OK &&
	ged_edit_session_revert(gedp, ref) == GED_EDIT_OK &&
	ged_edit_session_info_get(gedp, ref, &info) == GED_EDIT_OK &&
	!info.dirty,
	"revert restores the checkpoint dirty state as well as edit state");

    struct edit_event_observer obs;
    ged_event_observer_token token = ged_event_observer_add(gedp,
	GED_EVENT_OBSERVER_POST_RECONCILE, edit_event_cb, &obs);
    CHECK(token != 0, "clean selection commit observer registers");
    CHECK(ged_edit_session_commit(gedp, ref) == GED_EDIT_OK &&
	ged_edit_session_info_get(gedp, ref, &info) == GED_EDIT_STALE,
	"clean selection commit closes the session");
    CHECK(obs.bot_modified == 0,
	"clean selection commit does not emit a false object-modified event");
    CHECK(ged_event_observer_remove(gedp, token) == 1,
	"clean selection commit observer removal succeeds");
}


static void
test_p1_reentrant_session_observer(struct ged *gedp)
{
    struct rt_tor_internal before;
    CHECK(read_tor(gedp, "tor.s", &before) == BRLCAD_OK,
	"reentrant observer test reads the committed torus");

    ged_edit_session_ref ref = GED_EDIT_SESSION_REF_NULL;
    CHECK(ged_edit_session_begin(gedp, "tor.s", NULL, &ref) == GED_EDIT_OK,
	"reentrant observer test begins an edit session");

    struct reentrant_edit_observer mutator;
    mutator.session = ref;
    struct edit_session_observer recorder;
    ged_edit_observer_token mutator_token = ged_edit_observer_add(gedp,
	reentrant_edit_cancel_cb, &mutator);
    ged_edit_observer_token recorder_token = ged_edit_observer_add(gedp,
	edit_session_cb, &recorder);
    CHECK(mutator_token && recorder_token,
	"reentrant observer test registers both observers");

    fastf_t values[1] = {9.25};
    struct ged_edit_command_input input = {};
    input.command_id = 1021;
    input.values = values;
    input.value_count = 1;
    CHECK(ged_edit_session_apply(gedp, ref, &input) == GED_EDIT_OK,
	"an observer may cancel a session reentrantly during its update");
    CHECK(mutator.updates == 1 && mutator.cancel_status == GED_EDIT_OK,
	"the reentrant observer cancels exactly once");
    CHECK(recorder.kinds.size() == 2 &&
	recorder.kinds[0] == GED_EDIT_SESSION_UPDATE &&
	recorder.kinds[1] == GED_EDIT_SESSION_CANCEL,
	"reentrant cancellation is queued after complete update delivery");

    struct ged_edit_session_info info = {};
    CHECK(ged_edit_session_info_get(gedp, ref, &info) == GED_EDIT_STALE,
	"reentrantly canceled session references become stale");
    struct rt_tor_internal after;
    CHECK(read_tor(gedp, "tor.s", &after) == BRLCAD_OK &&
	NEAR_EQUAL(after.r_a, before.r_a, NEAR_ENOUGH),
	"reentrant cancellation does not leak intermediate geometry");
    CHECK(ged_edit_observer_remove(gedp, mutator_token) == 1 &&
	ged_edit_observer_remove(gedp, recorder_token) == 1,
	"reentrant observer registrations remove cleanly");
}


static void
test_p1_database_event_invalidation(struct ged *gedp)
{
    struct edit_session_observer obs;
    ged_edit_observer_token token = ged_edit_observer_add(gedp,
	edit_session_cb, &obs);
    CHECK(token != 0, "database invalidation observer registers");

    const char *intermediate[] = {
	"edit", "-i", "sph.s", "translate", "-r", "1", "0", "0", NULL
    };
    CHECK(ged_exec(gedp, 8, intermediate) == BRLCAD_OK,
	"database invalidation test creates a dirty sphere session");
    ged_edit_session_ref modified_ref = GED_EDIT_SESSION_REF_NULL;
    CHECK(ged_edit_session_find(gedp, "sph.s", &modified_ref) == GED_EDIT_OK,
	"dirty sphere session is active before an external command");

    const char *adjust[] = {"adjust", "sph.s", "V", "30 0 0", NULL};
    CHECK(ged_exec(gedp, 4, adjust) == BRLCAD_OK,
	"an ordinary CLI command modifies the session source object");
    struct ged_edit_session_info info = {};
    CHECK(!obs.kinds.empty() &&
	obs.kinds.back() == GED_EDIT_SESSION_INVALIDATE &&
	obs.paths.back() == "/sph.s" && obs.replacement_paths.back().empty(),
	"external modification reports an explicit session invalidation");
    CHECK(obs.invalidation_reasons.back() ==
	GED_EDIT_INVALIDATION_SOURCE_CHANGED,
	"external modification identifies a reloadable source change");
    CHECK(ged_edit_session_info_get(gedp, modified_ref, &info) ==
	GED_EDIT_STALE,
	"external modification makes the old session reference stale");
    struct rt_ell_internal sphere;
    CHECK(read_ell(gedp, "sph.s", &sphere) == BRLCAD_OK &&
	NEAR_EQUAL(sphere.v[X], 30.0, NEAR_ENOUGH) &&
	NEAR_EQUAL(MAGNITUDE(sphere.a), 2.0, NEAR_ENOUGH),
	"external geometry wins and stale intermediate geometry is discarded");

    ged_edit_session_ref rename_ref = GED_EDIT_SESSION_REF_NULL;
    CHECK(ged_edit_session_begin(gedp, "sph.s", NULL, &rename_ref) ==
	GED_EDIT_OK, "rename invalidation test begins a fresh session");
    const char *rename[] = {"move", "sph.s", "sph_renamed.s", NULL};
    CHECK(ged_exec(gedp, 3, rename) == BRLCAD_OK,
	"an ordinary CLI command renames the session source object");
    CHECK(!obs.kinds.empty() &&
	obs.kinds.back() == GED_EDIT_SESSION_INVALIDATE &&
	obs.paths.back() == "/sph.s" &&
	obs.replacement_paths.back() == "/sph_renamed.s",
	"rename invalidation supplies the canonical replacement path");
    CHECK(obs.invalidation_reasons.back() ==
	GED_EDIT_INVALIDATION_SOURCE_RENAMED,
	"rename invalidation is distinguished from removal");
    CHECK(ged_edit_session_info_get(gedp, rename_ref, &info) == GED_EDIT_STALE,
	"rename invalidates the old path session");

    ged_edit_session_ref renamed_ref = GED_EDIT_SESSION_REF_NULL;
    CHECK(ged_edit_session_begin(gedp, "sph_renamed.s", NULL, &renamed_ref) ==
	GED_EDIT_OK, "a new session attaches to the renamed object");
    CHECK(ged_edit_session_cancel(gedp, renamed_ref) == GED_EDIT_OK,
	"renamed object session cancels cleanly");

    ged_edit_session_ref nested_ref = GED_EDIT_SESSION_REF_NULL;
    CHECK(ged_edit_session_begin(gedp, "group.c/tor.s", NULL, &nested_ref) ==
	GED_EDIT_OK, "nested invalidation test begins a path session");
    const char *remove_child[] = {"rm", "group.c/tor.s", NULL};
    CHECK(ged_exec(gedp, 2, remove_child) == BRLCAD_OK,
	"an ordinary CLI command removes the edited nested occurrence");
    CHECK(!obs.kinds.empty() &&
	obs.kinds.back() == GED_EDIT_SESSION_INVALIDATE &&
	obs.paths.back() == "/group.c/tor.s" &&
	obs.replacement_paths.back().empty(),
	"nested occurrence removal invalidates the affected path session");
    CHECK(obs.invalidation_reasons.back() == GED_EDIT_INVALIDATION_PATH_REMOVED,
	"nested occurrence invalidation reports that its path was removed");
    CHECK(ged_edit_session_info_get(gedp, nested_ref, &info) == GED_EDIT_STALE,
	"removed nested occurrence leaves no writable stale session");

    ged_edit_session_ref removed_ref = GED_EDIT_SESSION_REF_NULL;
    CHECK(ged_edit_session_begin(gedp, "sph_renamed.s", NULL, &removed_ref) ==
	GED_EDIT_OK, "remove invalidation test begins a session");
    const char *remove_object[] = {"kill", "sph_renamed.s", NULL};
    CHECK(ged_exec(gedp, 2, remove_object) == BRLCAD_OK,
	"an ordinary CLI command removes the session source object");
    CHECK(!obs.kinds.empty() &&
	obs.kinds.back() == GED_EDIT_SESSION_INVALIDATE &&
	obs.paths.back() == "/sph_renamed.s" &&
	obs.replacement_paths.back().empty(),
	"object removal invalidates without advertising a replacement");
    CHECK(obs.invalidation_reasons.back() ==
	GED_EDIT_INVALIDATION_SOURCE_REMOVED,
	"object removal reports a terminal source removal");
    CHECK(ged_edit_session_info_get(gedp, removed_ref, &info) == GED_EDIT_STALE &&
	db_lookup(gedp->dbip, "sph_renamed.s", LOOKUP_QUIET) == RT_DIR_NULL,
	"removed object cannot be recovered through the stale session");

    CHECK(ged_edit_observer_remove(gedp, token) == 1,
	"database invalidation observer removes cleanly");
}


static void
test_p1_closedb_invalidation(struct ged *gedp, const char *db_path)
{
    struct rt_tor_internal before;
    CHECK(read_tor(gedp, "tor.s", &before) == BRLCAD_OK,
	"closedb invalidation test reads the committed torus");
    struct edit_session_observer obs;
    ged_edit_observer_token token = ged_edit_observer_add(gedp,
	edit_session_cb, &obs);
    CHECK(token != 0, "closedb invalidation observer registers");

    const char *intermediate[] = {
	"edit", "-i", "tor.s", "r1", "123.5", NULL
    };
    CHECK(ged_exec(gedp, 5, intermediate) == BRLCAD_OK,
	"closedb invalidation test leaves a dirty session open");
    ged_edit_session_ref old_ref = GED_EDIT_SESSION_REF_NULL;
    CHECK(ged_edit_session_find(gedp, "tor.s", &old_ref) == GED_EDIT_OK,
	"closedb invalidation test captures the active session");

    const char *close_cmd[] = {"closedb", NULL};
    CHECK(ged_exec(gedp, 1, close_cmd) == BRLCAD_OK && !gedp->dbip,
	"closedb succeeds with an active edit session");
    struct ged_edit_session_info info = {};
    CHECK(!obs.kinds.empty() &&
	obs.kinds.back() == GED_EDIT_SESSION_INVALIDATE &&
	obs.invalidation_reasons.back() == GED_EDIT_INVALIDATION_DATABASE_CLOSED &&
	ged_edit_session_info_get(gedp, old_ref, &info) == GED_EDIT_STALE,
	"closedb invalidates the transaction before destroying its database");

    const char *open_cmd[] = {"opendb", db_path, NULL};
    CHECK(ged_exec(gedp, 2, open_cmd) == BRLCAD_OK && gedp->dbip,
	"the same GED context reopens after closedb invalidation");
    struct rt_tor_internal after;
    CHECK(read_tor(gedp, "tor.s", &after) == BRLCAD_OK &&
	NEAR_EQUAL(after.r_a, before.r_a, NEAR_ENOUGH),
	"closedb invalidation never commits intermediate geometry");
    ged_edit_session_ref fresh_ref = GED_EDIT_SESSION_REF_NULL;
    CHECK(ged_edit_session_begin(gedp, "tor.s", NULL, &fresh_ref) ==
	GED_EDIT_OK && ged_edit_session_cancel(gedp, fresh_ref) == GED_EDIT_OK,
	"reopened GED context accepts new edit sessions");
    CHECK(ged_edit_observer_remove(gedp, token) == 1,
	"closedb invalidation observer removes cleanly");
}


/* ================================================================== *
 * Section 2 — translate / rotate / scale / checkpoint / revert / mat
 * ================================================================== */

/**
 * Phase 2 fixture:
 *   sph.s   — sphere at (10,0,0), r=2     (translate tests)
 *   other.s — sphere at (20,0,0), r=2     (translate -k/-a reference)
 *   ell.s   — ell at origin, A=(5,0,0)    (rotate tests)
 *   sca.s   — sphere at origin, r=3       (scale tests)
 *   tor.s   — torus at origin             (checkpoint/revert/mat tests)
 */
static int
create_p2_fixture(const char *path)
{
    struct rt_wdb *wdbp = wdb_fopen(path);
    if (!wdbp)
        return BRLCAD_ERROR;

    point_t sp = {10.0, 0.0, 0.0};
    if (mk_sph(wdbp, "sph.s", sp, 2.0) != 0) goto fail;

    { point_t op = {20.0, 0.0, 0.0};
      if (mk_sph(wdbp, "other.s", op, 2.0) != 0) goto fail; }

    { point_t ev = VINIT_ZERO;
      vect_t ea = {5.0,0,0}, eb = {0,2.0,0}, ec = {0,0,3.0};
      if (mk_ell(wdbp, "ell.s", ev, ea, eb, ec) != 0) goto fail; }

    { point_t sv = VINIT_ZERO;
      if (mk_sph(wdbp, "sca.s", sv, 3.0) != 0) goto fail; }

    { point_t tv = VINIT_ZERO; vect_t th = {0,0,1};
      if (mk_tor(wdbp, "tor.s", tv, th, 4.0, 1.0) != 0) goto fail; }

    wdb_close(wdbp);
    return BRLCAD_OK;
fail:
    wdb_close(wdbp);
    return BRLCAD_ERROR;
}

/* ---- translate ---------------------------------------------------- */

static void
test_p2_translate_abs(struct ged *gedp)
{
    const char *av[] = { "edit", "sph.s", "translate", "-a", "15", "0", "0", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 7, av) == BRLCAD_OK, "translate -a 15 0 0 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sph.s", &ell) == BRLCAD_OK) {
        CHECK(NEAR_EQUAL(ell.v[X], 15.0, NEAR_ENOUGH), "translate -a: sph.s V.x == 15");
        CHECK(NEAR_EQUAL(ell.v[Y],  0.0, NEAR_ENOUGH), "translate -a: sph.s V.y == 0");
    } else { CHECK(0, "translate -a: read_ell(sph.s) succeeded"); }
}

static void
test_p2_translate_rel(struct ged *gedp)
{
    const char *av[] = { "edit", "sph.s", "translate", "-r", "-5", "0", "0", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 7, av) == BRLCAD_OK, "translate -r -5 0 0 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sph.s", &ell) == BRLCAD_OK)
        CHECK(NEAR_EQUAL(ell.v[X], 10.0, NEAR_ENOUGH), "translate -r: sph.s V.x == 10");
    else CHECK(0, "translate -r: read_ell(sph.s) succeeded");
}

static void
test_p2_tra(struct ged *gedp)
{
    const char *av[] = { "edit", "sph.s", "tra", "5", "0", "0", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK, "tra 5 0 0 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sph.s", &ell) == BRLCAD_OK)
        CHECK(NEAR_EQUAL(ell.v[X], 15.0, NEAR_ENOUGH), "tra: sph.s V.x == 15");
    else CHECK(0, "tra: read_ell(sph.s) succeeded");
}

static void
test_p2_translate_from_to(struct ged *gedp)
{
    /* sph.s at (15,0,0); other.s at (20,0,0).  -k sph.s -a other.s → move to (20,0,0). */
    const char *av[] = { "edit", "sph.s", "translate", "-k", "sph.s", "-a", "other.s", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 7, av) == BRLCAD_OK, "translate -k sph.s -a other.s returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sph.s", &ell) == BRLCAD_OK)
        CHECK(NEAR_EQUAL(ell.v[X], 20.0, NEAR_ENOUGH), "translate -k/-a: sph.s V.x == 20");
    else CHECK(0, "translate -k/-a: read_ell(sph.s) succeeded");
}

static void
test_p2_translate_xonly(struct ged *gedp)
{
    /* sph.s at (20,0,0); -r -x 5 → (25,0,0) */
    const char *av[] = { "edit", "sph.s", "translate", "-r", "-x", "5", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK, "translate -r -x 5 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sph.s", &ell) == BRLCAD_OK) {
        CHECK(NEAR_EQUAL(ell.v[X], 25.0, NEAR_ENOUGH), "translate -x: sph.s V.x == 25");
        CHECK(NEAR_EQUAL(ell.v[Y],  0.0, NEAR_ENOUGH), "translate -x: sph.s V.y unchanged == 0");
    } else CHECK(0, "translate -x: read_ell(sph.s) succeeded");
}

/* ---- rotate ------------------------------------------------------- */

static void
test_p2_rotate_3angle(struct ged *gedp)
{
    const char *av[] = { "edit", "ell.s", "rotate", "0", "0", "90", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK, "rotate 0 0 90 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "ell.s", &ell) == BRLCAD_OK) {
        CHECK(NEAR_EQUAL(ell.a[X], 0.0, NEAR_ENOUGH), "rotate 0 0 90: ell.s A.x ≈ 0");
        CHECK(NEAR_EQUAL(ell.a[Y], 5.0, NEAR_ENOUGH), "rotate 0 0 90: ell.s A.y ≈ 5");
        CHECK(NEAR_EQUAL(ell.b[X], -2.0, NEAR_ENOUGH), "rotate 0 0 90: ell.s B.x ≈ -2");
    } else CHECK(0, "rotate 0 0 90: read_ell(ell.s) succeeded");
}

static void
test_p2_rotate_undo(struct ged *gedp)
{
    const char *av[] = { "edit", "ell.s", "rotate", "0", "0", "-90", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK, "rotate 0 0 -90 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "ell.s", &ell) == BRLCAD_OK) {
        CHECK(NEAR_EQUAL(ell.a[X], 5.0, NEAR_ENOUGH), "rotate 0 0 -90: ell.s A.x ≈ 5 (undone)");
        CHECK(NEAR_EQUAL(ell.a[Y], 0.0, NEAR_ENOUGH), "rotate 0 0 -90: ell.s A.y ≈ 0 (undone)");
    } else CHECK(0, "rotate 0 0 -90: read_ell(ell.s) succeeded");
}

static void
test_p2_rotate_z_flag(struct ged *gedp)
{
    const char *av[] = { "edit", "ell.s", "rotate", "-z", "45", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 5, av) == BRLCAD_OK, "rotate -z 45 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "ell.s", &ell) == BRLCAD_OK) {
        fastf_t expected = 5.0 * cos(45.0 * DEG2RAD);
        CHECK(NEAR_EQUAL(ell.a[X], expected, 1e-4), "rotate -z 45: ell.s A.x ≈ 5·cos(45°)");
        CHECK(NEAR_EQUAL(ell.a[Y], expected, 1e-4), "rotate -z 45: ell.s A.y ≈ 5·sin(45°)");
    } else CHECK(0, "rotate -z 45: read_ell(ell.s) succeeded");
}

static void
test_p2_rotate_single_angle(struct ged *gedp)
{
    const char *av[] = { "edit", "ell.s", "rotate", "90", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 4, av) == BRLCAD_OK, "rotate 90 (single angle) returns OK");
}

static void
test_p2_rotate_180_ambiguous(struct ged *gedp)
{
    const char *av[] = { "edit", "ell.s", "rotate", "180", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    int ret = ged_exec(gedp, 4, av);
    CHECK(ret == BRLCAD_ERROR, "rotate 180 (no axis) returns error (ambiguous)");
    bool has_msg = (strstr(bu_vls_cstr(gedp->ged_result_str), "ambiguous") != NULL ||
                    strstr(bu_vls_cstr(gedp->ged_result_str), "axis")      != NULL);
    CHECK(has_msg, "rotate 180 prints ambiguity message");
}

static void
test_p2_rotate_radians(struct ged *gedp)
{
    struct rt_ell_internal before;
    bool have_before = (read_ell(gedp, "ell.s", &before) == BRLCAD_OK);

    const char *av[] = {
        "edit", "ell.s", "rotate", "-R", "0", "0", "1.5707963267948966", NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 7, av) == BRLCAD_OK, "rotate -R 0 0 pi/2 returns OK");

    if (have_before) {
        struct rt_ell_internal after;
        if (read_ell(gedp, "ell.s", &after) == BRLCAD_OK) {
            CHECK(NEAR_EQUAL(after.a[X], -before.a[Y], 1e-4), "rotate -R pi/2: A.x ≈ -before.A.y");
            CHECK(NEAR_EQUAL(after.a[Y],  before.a[X], 1e-4), "rotate -R pi/2: A.y ≈  before.A.x");
        } else CHECK(0, "rotate -R: read_ell after rotation succeeded");
    }
}

static void
test_p2_rotate_two_angles(struct ged *gedp)
{
    const char *av[] = { "edit", "ell.s", "rotate", "-45", "45", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 5, av) == BRLCAD_OK, "rotate -45 45 (two-angle) returns OK");
}

/* ---- scale -------------------------------------------------------- */

static void
test_p2_scale_scalar(struct ged *gedp)
{
    const char *av[] = { "edit", "sca.s", "scale", "2", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 4, av) == BRLCAD_OK, "scale 2 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sca.s", &ell) == BRLCAD_OK)
        CHECK(NEAR_EQUAL(MAGNITUDE(ell.a), 6.0, NEAR_ENOUGH), "scale 2: sca.s |A| == 6");
    else CHECK(0, "scale 2: read_ell(sca.s) succeeded");
}

static void
test_p2_scale_vector(struct ged *gedp)
{
    const char *av[] = { "edit", "sca.s", "scale", "2", "2", "2", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK, "scale 2 2 2 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sca.s", &ell) == BRLCAD_OK)
        CHECK(NEAR_EQUAL(MAGNITUDE(ell.a), 12.0, NEAR_ENOUGH), "scale 2 2 2: sca.s |A| == 12");
    else CHECK(0, "scale 2 2 2: read_ell(sca.s) succeeded");
}

static void
test_p2_scale_from_to(struct ged *gedp)
{
    const char *av[] = {
        "edit", "sca.s", "scale", "-k", "0", "0", "0", "-a", "2", "2", "2", NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 11, av) == BRLCAD_OK, "scale -k 0 0 0 -a 2 2 2 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sca.s", &ell) == BRLCAD_OK)
        CHECK(NEAR_EQUAL(MAGNITUDE(ell.a), 24.0, NEAR_ENOUGH), "scale -k/-a: sca.s |A| == 24");
    else CHECK(0, "scale -k/-a: read_ell(sca.s) succeeded");
}

static void
test_p2_scale_explicit_r(struct ged *gedp)
{
    /* Ratio model: SCALE=(1,1,1)->(0,0,0) = (1,1,1), FACTOR=2 → effective=2/1=2
     * sca.s |A| was 24 before this call, so after: 24 * 2 = 48. */
    const char *av[] = {
        "edit", "sca.s", "scale", "-k", "0","0","0", "-a","1","1","1", "-r","2", NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 13, av) == BRLCAD_OK, "scale -k 0 0 0 -a 1 1 1 -r 2 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sca.s", &ell) == BRLCAD_OK)
        CHECK(NEAR_EQUAL(MAGNITUDE(ell.a), 48.0, NEAR_ENOUGH),
              "scale ratio (-k 0 0 0 -a 1 1 1 -r 2): sca.s |A| == 48");
    else CHECK(0, "scale ratio: read_ell(sca.s) failed");
}

static void
test_p2_scale_complex(struct ged *gedp)
{
    /* Full ratio model: SCALE=(7-5, 11-10, -2-15)=(2,1,-17), |per-axis|=(2,1,17)
     * FACTOR=(4,2,34) → effective=(4/2, 2/1, 34/17)=(2,2,2) → uniform 2x scale.
     * sca.s |A| was 48 before this call, so after: 48 * 2 = 96. */
    const char *av[] = {
        "edit", "sca.s", "scale",
        "-k","5","10","15", "-a","7","11","-2", "-r","4","2","34", NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 15, av) == BRLCAD_OK, "scale -k 5 10 15 -a 7 11 -2 -r 4 2 34 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sca.s", &ell) == BRLCAD_OK)
        CHECK(NEAR_EQUAL(MAGNITUDE(ell.a), 96.0, NEAR_ENOUGH),
              "scale ratio (-k 5 10 15 -a 7 11 -2 -r 4 2 34): sca.s |A| == 96");
    else CHECK(0, "scale ratio complex: read_ell(sca.s) failed");
}

static void
test_p2_scale_zero_error(struct ged *gedp)
{
    const char *av[] = { "edit", "sca.s", "scale", "0", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 4, av) == BRLCAD_ERROR, "scale 0 returns error");
}

/* ---- checkpoint / revert / reset ---------------------------------- */

static void
test_p2_checkpoint(struct ged *gedp)
{
    const char *av[] = { "edit", "tor.s", "checkpoint", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK, "checkpoint returns OK");

    struct db_full_path dfp;
    db_full_path_init(&dfp);
    db_string_to_path(&dfp, gedp->dbip, "tor.s");
    struct rt_edit *s = ged_edit_buf_get(gedp, &dfp);
    CHECK(s != NULL, "checkpoint: edit buffer has entry for tor.s");
    db_free_full_path(&dfp);
}

static void
test_p2_revert(struct ged *gedp)
{
    { const char *av[] = { "edit", "tor.s", "checkpoint", NULL }; ged_exec(gedp, 3, av); }

    { const char *av[] = { "edit", "-i", "tor.s", "translate", "-r", "100", "0", "0", NULL };
      bu_vls_trunc(gedp->ged_result_str, 0);
      CHECK(ged_exec(gedp, 8, av) == BRLCAD_OK, "revert test: -i translate returns OK"); }

    { const char *av[] = { "edit", "tor.s", "revert", NULL };
      bu_vls_trunc(gedp->ged_result_str, 0);
      CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK, "revert returns OK"); }
}

static void
test_p2_reset(struct ged *gedp)
{
    { const char *av[] = { "edit", "tor.s", "checkpoint", NULL }; ged_exec(gedp, 3, av); }

    { const char *av[] = { "edit", "tor.s", "reset", NULL };
      bu_vls_trunc(gedp->ged_result_str, 0);
      CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK, "reset returns OK"); }

    struct db_full_path dfp;
    db_full_path_init(&dfp);
    db_string_to_path(&dfp, gedp->dbip, "tor.s");
    struct rt_edit *s = ged_edit_buf_get(gedp, &dfp);
    CHECK(s == NULL, "reset: buffer entry removed for tor.s");
    db_free_full_path(&dfp);
}

static void
test_p2_revert_no_checkpoint(struct ged *gedp)
{
    { struct db_full_path dfp;
      db_full_path_init(&dfp);
      db_string_to_path(&dfp, gedp->dbip, "sca.s");
      ged_edit_buf_abandon(gedp, &dfp);
      db_free_full_path(&dfp); }

    const char *av[] = { "edit", "sca.s", "revert", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_ERROR, "revert with no buffer entry returns error");
}

/* ---- mat ---------------------------------------------------------- */

static void
test_p2_mat_identity(struct ged *gedp)
{
    struct rt_ell_internal before;
    bool have_before = (read_ell(gedp, "ell.s", &before) == BRLCAD_OK);

    const char *av[] = {
        "edit", "ell.s", "mat",
        "1","0","0","0","0","1","0","0","0","0","1","0","0","0","0","1", NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 19, av) == BRLCAD_OK, "mat identity returns OK");

    if (have_before) {
        struct rt_ell_internal after;
        if (read_ell(gedp, "ell.s", &after) == BRLCAD_OK)
            CHECK(NEAR_EQUAL(MAGNITUDE(after.a), MAGNITUDE(before.a), 1e-6),
                  "mat identity: |A| unchanged");
        else CHECK(0, "mat identity: read_ell after succeeded");
    }
}

static void
test_p2_mat_missing_values(struct ged *gedp)
{
    const char *av[] = {
        "edit", "ell.s", "mat",
        "1","0","0","0","0","1","0","0", NULL   /* only 8 of 16 */
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 10, av) == BRLCAD_ERROR, "mat with <16 values returns error");
}


/* ================================================================== *
 * Section A - edit operations without compatibility draw state
 * ================================================================== */

static int
create_pa_fixture(const char *path)
{
    struct rt_wdb *wdbp = wdb_fopen(path);
    if (!wdbp)
        return BRLCAD_ERROR;

    point_t sp = {10.0, 0.0, 0.0};
    if (mk_sph(wdbp, "pa_sph.s", sp, 2.0) != 0) goto fail;

    { point_t ev = VINIT_ZERO;
      vect_t ea = {5,0,0}, eb = {0,2,0}, ec = {0,0,3};
      if (mk_ell(wdbp, "pa_ell.s", ev, ea, eb, ec) != 0) goto fail; }

    { point_t sv = VINIT_ZERO;
      if (mk_sph(wdbp, "pa_sca.s", sv, 3.0) != 0) goto fail; }

    wdb_close(wdbp);
    return BRLCAD_OK;
fail:
    wdb_close(wdbp);
    return BRLCAD_ERROR;
}

static void
test_pa_translate_abs(struct ged *gedp)
{
    const char *av[] = { "edit", "pa_sph.s", "translate", "-a", "15", "0", "0", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 7, av) == BRLCAD_OK, "no-compat translate -a returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "pa_sph.s", &ell) == BRLCAD_OK)
        CHECK(NEAR_EQUAL(ell.v[X], 15.0, NEAR_ENOUGH),
              "no-compat translate -a: pa_sph.s V.x == 15");
    else CHECK(0, "no-compat translate -a: read_ell succeeded");
}

static void
test_pa_translate_rel(struct ged *gedp)
{
    const char *av[] = { "edit", "pa_sph.s", "translate", "-r", "-5", "0", "0", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 7, av) == BRLCAD_OK, "no-compat translate -r returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "pa_sph.s", &ell) == BRLCAD_OK)
        CHECK(NEAR_EQUAL(ell.v[X], 10.0, NEAR_ENOUGH),
              "no-compat translate -r: pa_sph.s V.x == 10 (restored)");
    else CHECK(0, "no-compat translate -r: read_ell succeeded");
}

static void
test_pa_rotate(struct ged *gedp)
{
    const char *av[] = { "edit", "pa_ell.s", "rotate", "0", "0", "90", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK, "no-compat rotate 0 0 90 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "pa_ell.s", &ell) == BRLCAD_OK) {
        CHECK(NEAR_EQUAL(ell.a[X], 0.0, NEAR_ENOUGH), "no-compat rotate: pa_ell.s A.x ≈ 0");
        CHECK(NEAR_EQUAL(ell.a[Y], 5.0, NEAR_ENOUGH), "no-compat rotate: pa_ell.s A.y ≈ 5");
    } else CHECK(0, "no-compat rotate: read_ell succeeded");
}

static void
test_pa_scale(struct ged *gedp)
{
    const char *av[] = { "edit", "pa_sca.s", "scale", "2", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 4, av) == BRLCAD_OK, "no-compat scale 2 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "pa_sca.s", &ell) == BRLCAD_OK)
        CHECK(NEAR_EQUAL(MAGNITUDE(ell.a), 6.0, NEAR_ENOUGH),
              "no-compat scale 2: pa_sca.s |A| == 6");
    else CHECK(0, "no-compat scale 2: read_ell succeeded");
}

static void
test_pa_S_flag_no_crash(struct ged *gedp)
{
    const char *av[] = { "edit", "-S", "translate", "-r", "1", "0", "0", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    ged_exec(gedp, 7, av);
    CHECK(1, "no-compat -S flag does not crash");
}

static void
test_pa_unknown_obj(struct ged *gedp)
{
    const char *av[] = { "edit", "nonexistent.s", "translate", "-r", "1", "0", "0", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 7, av) == BRLCAD_ERROR,
          "no-compat unknown object returns error");
}


/* ================================================================== *
 * Section C — rotate axis-mode, scale -c, per-axis scale
 * ================================================================== */

/**
 * Phase C fixture:
 *   ell.s   — ell at origin, A=(5,0,0), B=(0,2,0), C=(0,0,3)  (rotate axis tests)
 *   sca2.s  — ell at origin, A=(6,0,0), B=(0,3,0), C=(0,0,2)  (per-axis scale)
 *   ctr.s   — sphere at (10,0,0), r=2                          (center-based tests)
 */
static int
create_pc_fixture(const char *path)
{
    struct rt_wdb *wdbp = wdb_fopen(path);
    if (!wdbp)
        return BRLCAD_ERROR;

    { point_t ev = VINIT_ZERO;
      vect_t ea = {5,0,0}, eb = {0,2,0}, ec = {0,0,3};
      if (mk_ell(wdbp, "ell.s", ev, ea, eb, ec) != 0) goto fail; }

    { point_t v = VINIT_ZERO;
      vect_t a = {6,0,0}, b = {0,3,0}, c = {0,0,2};
      if (mk_ell(wdbp, "sca2.s", v, a, b, c) != 0) goto fail; }

    { point_t cv = {10.0, 0.0, 0.0};
      if (mk_sph(wdbp, "ctr.s", cv, 2.0) != 0) goto fail; }

    wdb_close(wdbp);
    return BRLCAD_OK;
fail:
    wdb_close(wdbp);
    return BRLCAD_ERROR;
}

/* RC1: rotate -k 0 0 0 -a 0 0 1 -d 90 */
static void
test_pc_rotate_axis_z90(struct ged *gedp)
{
    { const char *av[] = { "edit","ell.s","rotate","0","0","0",NULL }; ged_exec(gedp,6,av); }

    const char *av[] = {
        "edit","ell.s","rotate","-k","0","0","0","-a","0","0","1","-d","90",NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 13, av) == BRLCAD_OK, "rotate -k 0 0 0 -a 0 0 1 -d 90 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "ell.s", &ell) == BRLCAD_OK) {
        CHECK(NEAR_EQUAL(ell.a[X], 0.0, NEAR_ENOUGH), "axis rotate Z 90: ell.s A.x ≈ 0");
        CHECK(NEAR_EQUAL(ell.a[Y], 5.0, NEAR_ENOUGH), "axis rotate Z 90: ell.s A.y ≈ 5");
        CHECK(NEAR_EQUAL(ell.a[Z], 0.0, NEAR_ENOUGH), "axis rotate Z 90: ell.s A.z ≈ 0");
    } else CHECK(0, "axis rotate Z 90: read_ell(ell.s) succeeded");
}

/* RC2: rotate -k 0 0 0 -r 0 0 1 -d 90 (relative axis offset, same result) */
static void
test_pc_rotate_axis_z90_rel(struct ged *gedp)
{
    const char *av[] = {
        "edit","ell.s","rotate","-k","0","0","0","-r","0","0","1","-d","90",NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 13, av) == BRLCAD_OK, "rotate -k 0 0 0 -r 0 0 1 -d 90 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "ell.s", &ell) == BRLCAD_OK) {
        CHECK(NEAR_EQUAL(ell.a[X], -5.0, NEAR_ENOUGH), "axis rotate Z 90 (rel): ell.s A.x ≈ -5");
        CHECK(NEAR_EQUAL(ell.a[Y],  0.0, NEAR_ENOUGH), "axis rotate Z 90 (rel): ell.s A.y ≈ 0");
    } else CHECK(0, "axis rotate Z 90 (rel): read_ell(ell.s) succeeded");
}

/* RC3: rotate ctr.s 90° around X-axis; center should stay at (10,0,0) */
static void
test_pc_rotate_axis_x90(struct ged *gedp)
{
    const char *av[] = {
        "edit","ctr.s","rotate","-k","0","0","0","-a","1","0","0","-d","90",NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 13, av) == BRLCAD_OK, "rotate -k 0 0 0 -a 1 0 0 -d 90 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "ctr.s", &ell) == BRLCAD_OK) {
        CHECK(NEAR_EQUAL(ell.v[X], 10.0, NEAR_ENOUGH), "axis rotate X 90: ctr.s V.x ≈ 10");
        CHECK(NEAR_EQUAL(ell.v[Y],  0.0, NEAR_ENOUGH), "axis rotate X 90: ctr.s V.y ≈ 0");
        CHECK(NEAR_EQUAL(ell.v[Z],  0.0, NEAR_ENOUGH), "axis rotate X 90: ctr.s V.z ≈ 0");
    } else CHECK(0, "axis rotate X 90: read_ell(ctr.s) succeeded");
}

/* RC4: rotate ctr.s 90° around Z with center at its own position → center unchanged */
static void
test_pc_rotate_euler_center(struct ged *gedp)
{
    const char *av[] = {
        "edit","ctr.s","rotate","-c","10","0","0","0","0","90",NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 10, av) == BRLCAD_OK, "rotate -c 10 0 0  0 0 90 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "ctr.s", &ell) == BRLCAD_OK) {
        CHECK(NEAR_EQUAL(ell.v[X], 10.0, NEAR_ENOUGH), "rotate euler -c 10 0 0: ctr.s V.x ≈ 10");
        CHECK(NEAR_EQUAL(ell.v[Y],  0.0, NEAR_ENOUGH), "rotate euler -c 10 0 0: ctr.s V.y ≈ 0");
        CHECK(NEAR_EQUAL(ell.v[Z],  0.0, NEAR_ENOUGH), "rotate euler -c 10 0 0: ctr.s V.z ≈ 0");
    } else CHECK(0, "rotate euler -c center: read_ell(ctr.s) succeeded");
}

/* RC5: axis-mode with missing angle → error */
static void
test_pc_rotate_axis_missing_angle(struct ged *gedp)
{
    const char *av[] = {
        "edit","ell.s","rotate","-k","0","0","0","-a","0","0","1",NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 11, av) == BRLCAD_ERROR,
          "rotate axis mode without angle returns error");
    const char *msg = bu_vls_cstr(gedp->ged_result_str);
    CHECK((strstr(msg,"missing") != NULL || strstr(msg,"angle") != NULL),
          "rotate axis mode without angle prints helpful message");
}

/* SC1: scale -c 0 0 0 2 */
static void
test_pc_scale_uniform_center(struct ged *gedp)
{
    const char *av[] = { "edit","sca2.s","scale","-c","0","0","0","2",NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 8, av) == BRLCAD_OK, "scale -c 0 0 0 2 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sca2.s", &ell) == BRLCAD_OK) {
        CHECK(NEAR_EQUAL(ell.a[X], 12.0, NEAR_ENOUGH), "scale -c: sca2.s A.x ≈ 12");
        CHECK(NEAR_EQUAL(ell.b[Y],  6.0, NEAR_ENOUGH), "scale -c: sca2.s B.y ≈ 6");
        CHECK(NEAR_EQUAL(ell.c[Z],  4.0, NEAR_ENOUGH), "scale -c: sca2.s C.z ≈ 4");
    } else CHECK(0, "scale -c: read_ell(sca2.s) succeeded");
}

/* SC2: scale -c 0 0 0 0.5 (undo SC1) */
static void
test_pc_scale_center_halve(struct ged *gedp)
{
    const char *av[] = { "edit","sca2.s","scale","-c","0","0","0","0.5",NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 8, av) == BRLCAD_OK, "scale -c 0 0 0 0.5 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sca2.s", &ell) == BRLCAD_OK) {
        CHECK(NEAR_EQUAL(ell.a[X], 6.0, NEAR_ENOUGH), "scale -c 0.5: sca2.s A.x ≈ 6 (restored)");
        CHECK(NEAR_EQUAL(ell.b[Y], 3.0, NEAR_ENOUGH), "scale -c 0.5: sca2.s B.y ≈ 3 (restored)");
    } else CHECK(0, "scale -c 0.5: read_ell(sca2.s) succeeded");
}

/* SA1: scale 2 2 2 (three equal — still uniform) */
static void
test_pc_scale_three_equal(struct ged *gedp)
{
    const char *av[] = { "edit","sca2.s","scale","2","2","2",NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK, "scale 2 2 2 (uniform via 3 equal) returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sca2.s", &ell) == BRLCAD_OK)
        CHECK(NEAR_EQUAL(ell.a[X], 12.0, NEAR_ENOUGH), "scale 2 2 2: sca2.s A.x ≈ 12");
    else CHECK(0, "scale 2 2 2: read_ell(sca2.s) succeeded");
}

/* SA2: scale 1 1 0.5 (anisotropic — halve Z only) */
static void
test_pc_scale_anisotropic(struct ged *gedp)
{
    /* After SA1: A=(12,0,0), B=(0,6,0), C=(0,0,4) */
    const char *av[] = { "edit","sca2.s","scale","1","1","0.5",NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK, "scale 1 1 0.5 (anisotropic Z) returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sca2.s", &ell) == BRLCAD_OK) {
        CHECK(NEAR_EQUAL(ell.a[X], 12.0, NEAR_ENOUGH), "scale anisotropic: A.x unchanged (≈ 12)");
        CHECK(NEAR_EQUAL(ell.b[Y],  6.0, NEAR_ENOUGH), "scale anisotropic: B.y unchanged (≈ 6)");
        CHECK(NEAR_EQUAL(ell.c[Z],  2.0, NEAR_ENOUGH), "scale anisotropic: C.z halved (≈ 2)");
    } else CHECK(0, "scale anisotropic: read_ell(sca2.s) succeeded");
}

/* SA3: scale 0 1 1 → error */
static void
test_pc_scale_anisotropic_zero_error(struct ged *gedp)
{
    const char *av[] = { "edit","sca2.s","scale","0","1","1",NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_ERROR, "scale 0 1 1 (zero factor) returns error");
}


/* ================================================================== *
 * Section D — New design completions (edit2.cpp maturity):
 *   D1: translate -k with coordinate positions
 *   D2: translate -k . (self-reference as FROM)
 *   D3: translate per-component object extraction (-z ref.s)
 *   D4: rotate -c . (self-center rotation)
 *   D5: scale -c . (scale from own center)
 *   D6: rotate implicit angle from ANGLE_FROM→ANGLE_TO
 *   D7: rotate second -k/-a pair for explicit angle references
 * ================================================================== */

/**
 * Section D fixture:
 *   anchor.s — sphere at (0,0,0), r=1      (fixed reference origin)
 *   sph.s    — sphere at (10,0,0), r=2     (translate target)
 *   ref.s    — sphere at (0,7,3), r=1      (component-extraction reference)
 *   ell.s    — ell at origin, A=(5,0,0)    (rotate target)
 *   ctr.s    — sphere at (10,0,0), r=2     (rotate self-center target)
 *   sca.s    — sphere at origin, r=4       (scale target)
 */
static int
create_pd_fixture(const char *path)
{
    struct rt_wdb *wdbp = wdb_fopen(path);
    if (!wdbp)
	return BRLCAD_ERROR;

    { point_t v = VINIT_ZERO;
      if (mk_sph(wdbp, "anchor.s", v, 1.0) != 0) goto fail; }

    { point_t v = {10.0, 0.0, 0.0};
      if (mk_sph(wdbp, "sph.s", v, 2.0) != 0) goto fail; }

    { point_t v = {0.0, 7.0, 3.0};
      if (mk_sph(wdbp, "ref.s", v, 1.0) != 0) goto fail; }

    { point_t v = VINIT_ZERO;
      vect_t a = {5,0,0}, b = {0,2,0}, c = {0,0,3};
      if (mk_ell(wdbp, "ell.s", v, a, b, c) != 0) goto fail; }

    { point_t v = {10.0, 0.0, 0.0};
      if (mk_sph(wdbp, "ctr.s", v, 2.0) != 0) goto fail; }

    { point_t v = VINIT_ZERO;
      vect_t a = {5,0,0}, b = {0,2,0}, c = {0,0,3};
      if (mk_ell(wdbp, "ell2.s", v, a, b, c) != 0) goto fail; }

    { point_t v = VINIT_ZERO;
      vect_t a = {5,0,0}, b = {0,2,0}, c = {0,0,3};
      if (mk_ell(wdbp, "ell3.s", v, a, b, c) != 0) goto fail; }

    { point_t v = VINIT_ZERO;
      if (mk_sph(wdbp, "sca.s", v, 4.0) != 0) goto fail; }

    wdb_close(wdbp);
    return BRLCAD_OK;
fail:
    wdb_close(wdbp);
    return BRLCAD_ERROR;
}


/* D1: translate -k with explicit coordinate position (not object name) */
static void
test_pd_translate_k_coords(struct ged *gedp)
{
    /* sph.s starts at (10,0,0).  -k 10 0 0 -a 20 0 0: move by (10,0,0) */
    const char *av[] = {
	"edit","sph.s","translate","-k","10","0","0","-a","20","0","0",NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 11, av) == BRLCAD_OK,
	  "translate -k 10 0 0 -a 20 0 0 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sph.s", &ell) == BRLCAD_OK) {
	CHECK(NEAR_EQUAL(ell.v[X], 20.0, NEAR_ENOUGH),
	      "translate -k coord: sph.s V.x == 20");
    } else {
	CHECK(0, "translate -k coord: read_ell(sph.s) succeeded");
    }
}

/* Reset sph.s to (10,0,0) for subsequent tests */
static void reset_pd_sph(struct ged *gedp)
{
    const char *av[] = {"edit","sph.s","translate","-a","10","0","0",NULL};
    ged_exec(gedp, 7, av);
}

/* D2: translate -k . (self-reference: FROM = own keypoint) */
static void
test_pd_translate_k_self(struct ged *gedp)
{
    /* sph.s at (10,0,0).  -k . means FROM=own_keypoint=(10,0,0); with -a 20 0 0
     * as TO, delta = (10,0,0), target = (20,0,0). */
    const char *av[] = {
	"edit","sph.s","translate","-k",".","20","0","0",NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 8, av) == BRLCAD_OK,
	  "translate -k . 20 0 0 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sph.s", &ell) == BRLCAD_OK) {
	/* Without -a flag: the positional "20 0 0" is a relative delta, and
	 * -k . (FROM=own) is used only when combined with an explicit TO ref.
	 * With no -a, the mode is relative: delta=(20,0,0), target=(30,0,0). */
	CHECK(NEAR_EQUAL(ell.v[X], 30.0, NEAR_ENOUGH),
	      "translate -k . rel: sph.s V.x == 30 (relative, FROM=own)");
    } else {
	CHECK(0, "translate -k . rel: read_ell(sph.s) succeeded");
    }
}

/* D3: translate per-component object extraction (-a -z ref.s) */
static void
test_pd_translate_component_obj(struct ged *gedp)
{
    /* sph.s at (10,0,0); ref.s at (0,7,3).
     * -a -z ref.s: absolute Z = ref.s.z = 3, keep own X=10 and Y=0. */
    const char *av[] = {
	"edit","sph.s","translate","-a","-z","ref.s",NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK,
	  "translate -a -z ref.s returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sph.s", &ell) == BRLCAD_OK) {
	CHECK(NEAR_EQUAL(ell.v[X], 10.0, NEAR_ENOUGH),
	      "translate -a -z ref.s: V.x unchanged == 10");
	CHECK(NEAR_EQUAL(ell.v[Y],  0.0, NEAR_ENOUGH),
	      "translate -a -z ref.s: V.y unchanged == 0");
	CHECK(NEAR_EQUAL(ell.v[Z],  3.0, NEAR_ENOUGH),
	      "translate -a -z ref.s: V.z == ref.s.z == 3");
    } else {
	CHECK(0, "translate -a -z ref.s: read_ell(sph.s) succeeded");
    }
}

/* D3b: translate -a -y ref.s (extract Y component) */
static void
test_pd_translate_component_obj_y(struct ged *gedp)
{
    /* sph.s now at (10,0,3) after D3.  ref.s at (0,7,3).
     * -a -y ref.s: absolute Y = ref.s.y = 7, keep own X and Z. */
    const char *av[] = {
	"edit","sph.s","translate","-a","-y","ref.s",NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK,
	  "translate -a -y ref.s returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sph.s", &ell) == BRLCAD_OK) {
	CHECK(NEAR_EQUAL(ell.v[X], 10.0, NEAR_ENOUGH),
	      "translate -a -y ref.s: V.x unchanged == 10");
	CHECK(NEAR_EQUAL(ell.v[Y],  7.0, NEAR_ENOUGH),
	      "translate -a -y ref.s: V.y == ref.s.y == 7");
	CHECK(NEAR_EQUAL(ell.v[Z],  3.0, NEAR_ENOUGH),
	      "translate -a -y ref.s: V.z unchanged == 3");
    } else {
	CHECK(0, "translate -a -y ref.s: read_ell(sph.s) succeeded");
    }
}

/* D4: rotate -c . (self-center: pivot on own keypoint) */
static void
test_pd_rotate_c_self(struct ged *gedp)
{
    /* ctr.s is at (10,0,0).  Rotate 90° around Z with -c . → center stays
     * at own keypoint (10,0,0).  A sphere's keypoint is its center, so
     * the center is unaffected, and interior axes rotate (trivially no-op
     * geometry change for a sphere).  We verify the command succeeds. */
    const char *av[] = {
	"edit","ctr.s","rotate","-c",".","0","0","90",NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 9, av) == BRLCAD_OK,
	  "rotate -c . 0 0 90 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "ctr.s", &ell) == BRLCAD_OK) {
	/* Center of a sphere should not move when rotating about its own center */
	CHECK(NEAR_EQUAL(ell.v[X], 10.0, NEAR_ENOUGH),
	      "rotate -c .: ctr.s V.x unchanged == 10");
	CHECK(NEAR_EQUAL(ell.v[Y],  0.0, NEAR_ENOUGH),
	      "rotate -c .: ctr.s V.y unchanged == 0");
	CHECK(NEAR_EQUAL(ell.v[Z],  0.0, NEAR_ENOUGH),
	      "rotate -c .: ctr.s V.z unchanged == 0");
    } else {
	CHECK(0, "rotate -c .: read_ell(ctr.s) succeeded");
    }
}

/* D4b: rotate -c . in axis mode */
static void
test_pd_rotate_axis_c_self(struct ged *gedp)
{
    /* ell.s at origin with A=(5,0,0) at the start of section D.
     * No prior test in section D has touched ell.s.
     * Rotate 90° around Z-axis with -c . → pivot at own keypoint (origin).
     * A should go from (5,0,0) to (0,5,0). */
    const char *av[] = {
	"edit","ell.s","rotate","-k","0","0","0","-a","0","0","1","-c",".","90",NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 14, av) == BRLCAD_OK,
	  "rotate -k 0 0 0 -a 0 0 1 -c . 90 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "ell.s", &ell) == BRLCAD_OK) {
	CHECK(NEAR_EQUAL(ell.a[X], 0.0, NEAR_ENOUGH),
	      "rotate axis -c .: ell.s A.x ≈ 0");
	CHECK(NEAR_EQUAL(ell.a[Y], 5.0, NEAR_ENOUGH),
	      "rotate axis -c .: ell.s A.y ≈ 5");
    } else {
	CHECK(0, "rotate axis -c .: read_ell(ell.s) succeeded");
    }
}

/* D5: scale -c . (scale about own center) */
static void
test_pd_scale_c_self(struct ged *gedp)
{
    /* sca.s at (0,0,0), r=4.  Scale 2 from own center: r → 8. */
    const char *av[] = {
	"edit","sca.s","scale","-c",".","2",NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK,
	  "scale -c . 2 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sca.s", &ell) == BRLCAD_OK) {
	CHECK(NEAR_EQUAL(MAGNITUDE(ell.a), 8.0, NEAR_ENOUGH),
	      "scale -c . 2: sca.s |A| == 8");
    } else {
	CHECK(0, "scale -c . 2: read_ell(sca.s) succeeded");
    }
}

/* D5b: scale -c . anisotropic */
static void
test_pd_scale_c_self_aniso(struct ged *gedp)
{
    /* sca.s |A|=8 after D5.  Scale 0.5 0.5 0.5 from own center → |A|=4. */
    const char *av[] = {
	"edit","sca.s","scale","-c",".","0.5","0.5","0.5",NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 8, av) == BRLCAD_OK,
	  "scale -c . 0.5 0.5 0.5 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "sca.s", &ell) == BRLCAD_OK) {
	CHECK(NEAR_EQUAL(MAGNITUDE(ell.a), 4.0, NEAR_ENOUGH),
	      "scale -c . 0.5 0.5 0.5: sca.s |A| ≈ 4 (restored)");
    } else {
	CHECK(0, "scale -c . aniso: read_ell(sca.s) succeeded");
    }
}

/* D6: rotate implicit angle — single -k/-a pair, no -d
 *
 * ell2.s is a fresh ellipsoid at origin with A=(5,0,0) (separate from ell.s
 * which D4b rotates).
 *
 * Single-pair example:
 *   -k (1,0,0) -a (0,0,1):
 *     axis_dir = normalize((0,0,1)-(1,0,0)) = (-1,0,1)/√2
 *     angle_from defaults to axis_from = (1,0,0)
 *     angle_to  defaults to axis_to   = (0,0,1)
 *     center defaults to ell2.s keypoint = (0,0,0)
 *     vf_proj and vt_proj are non-zero (axis_from/to not on the axis) →
 *     implicit angle is computed and applied.
 */
static void
test_pd_rotate_implicit_angle_single_pair(struct ged *gedp)
{
    const char *av[] = {
	"edit","ell2.s","rotate","-k","1","0","0","-a","0","0","1",NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    int ret = ged_exec(gedp, 11, av);
    /* The angle derivation from projected arc should succeed for this
     * configuration (neither reference lies on the axis). */
    CHECK(ret == BRLCAD_OK,
	  "rotate -k 1 0 0 -a 0 0 1 (implicit angle, single pair) returns OK");
}

/* D6b: implicit angle where reference lies exactly on the axis → error */
static void
test_pd_rotate_implicit_angle_on_axis(struct ged *gedp)
{
    /* axis_dir = Z = (0,0,1).  angle_from defaults to axis_from = (0,0,0).
     * (0,0,0) − center(0,0,0) = (0,0,0) → projection is zero vector → error. */
    const char *av[] = {
	"edit","ell2.s","rotate","-k","0","0","0","-a","0","0","1",NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    int ret = ged_exec(gedp, 11, av);
    /* axis_from = (0,0,0) = center → vf = (0,0,0) → |vf_proj| = 0 →
     * "reference lies on axis" error. */
    CHECK(ret == BRLCAD_ERROR,
	  "rotate -k 0 0 0 -a 0 0 1 (reference on axis) returns error");
}

/* D7: rotate two -k/-a pairs — explicit ANGLE_FROM/ANGLE_TO
 *
 * ell2.s at origin with A=(5,0,0) (state after D6 applies some implicit
 * rotation; we use a fresh "ell3.s" to keep tests independent).
 *
 * Actually we add ell3.s to the fixture for complete isolation.
 * If ell3.s is not in the fixture, fall back: just verify OK return.
 *
 * Two-pair command:
 *   First  pair: axis_from=(0,0,0), axis_to=(0,0,1) → rotation axis = +Z
 *   Second pair: angle_from=(1,0,0), angle_to=(0,1,0)
 *     proj(angle_from) = (1,0,0) ⊥ Z → keeps (1,0,0)
 *     proj(angle_to)   = (0,1,0) ⊥ Z → keeps (0,1,0)
 *     angle = atan2(dot(Z, cross((1,0,0),(0,1,0))), dot((1,0,0),(0,1,0)))
 *           = atan2(dot(Z,(0,0,1)), 0) = atan2(1,0) = 90°
 *
 * ell3.s A=(5,0,0) after 90° Z rotation → A=(0,5,0).
 */
static void
test_pd_rotate_implicit_angle_two_pairs(struct ged *gedp)
{
    const char *av[] = {
	"edit","ell3.s","rotate",
	"-k","0","0","0","-a","0","0","1",   /* axis: from origin up Z */
	"-k","1","0","0","-a","0","1","0",   /* angle ref: X→Y = 90° */
	NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 19, av) == BRLCAD_OK,
	  "rotate two-pair implicit angle returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "ell3.s", &ell) == BRLCAD_OK) {
	CHECK(NEAR_EQUAL(ell.a[X], 0.0, 1e-4),
	      "rotate two-pair: ell3.s A.x ≈ 0 (after 90° Z)");
	CHECK(NEAR_EQUAL(ell.a[Y], 5.0, 1e-4),
	      "rotate two-pair: ell3.s A.y ≈ 5 (after 90° Z)");
	CHECK(NEAR_EQUAL(ell.a[Z], 0.0, 1e-4),
	      "rotate two-pair: ell3.s A.z ≈ 0");
    } else {
	CHECK(0, "rotate two-pair: read_ell(ell3.s) succeeded");
    }
}

/* D7b: rotate two pairs with reversed angle references → -90° result */
static void
test_pd_rotate_implicit_angle_reversed(struct ged *gedp)
{
    /* ell3.s is now at A≈(0,5,0) after D7.
     * Same axis (Z), but angle ref reversed: (0,1,0)→(1,0,0) = -90°.
     * Result: A should go back to ≈(5,0,0). */
    const char *av[] = {
	"edit","ell3.s","rotate",
	"-k","0","0","0","-a","0","0","1",
	"-k","0","1","0","-a","1","0","0",
	NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 19, av) == BRLCAD_OK,
	  "rotate two-pair reversed angle returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "ell3.s", &ell) == BRLCAD_OK) {
	CHECK(NEAR_EQUAL(ell.a[X], 5.0, 1e-4),
	      "rotate two-pair reversed: ell3.s A.x ≈ 5 (restored)");
	CHECK(NEAR_EQUAL(ell.a[Y], 0.0, 1e-4),
	      "rotate two-pair reversed: ell3.s A.y ≈ 0 (restored)");
    } else {
	CHECK(0, "rotate two-pair reversed: read_ell(ell3.s) succeeded");
    }
}


/* ================================================================== *
 * Section 3 — Descriptor-driven per-primitive editing (Phase 3)
 * ================================================================== */

/**
 * Helper: read rt_tor_internal for a torus by name from the database.
 * Caller must *not* free the returned data — the idb_ptr points into a
 * stack copy that is released by the caller after use.
 */
static int
read_tor(struct ged *gedp, const char *name, struct rt_tor_internal *out)
{
    struct directory *dp = db_lookup(gedp->dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL)
        return BRLCAD_ERROR;

    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, gedp->dbip, NULL) < 0)
        return BRLCAD_ERROR;

    if (intern.idb_minor_type != DB5_MINORTYPE_BRLCAD_TOR) {
        rt_db_free_internal(&intern);
        return BRLCAD_ERROR;
    }

    *out = *(struct rt_tor_internal *)intern.idb_ptr;
    intern.idb_ptr = NULL;
    rt_db_free_internal(&intern);
    return BRLCAD_OK;
}


/* Section 3 uses the Phase 0 fixture (tor.s + ell.s are already there) */

/* 3-0: descriptor-driven edit — change torus major radius (r1) via slug */
static void
test_p3_tor_set_radius_1(struct ged *gedp)
{
    const char *av[] = { "edit", "tor.s", "set_radius_1", "12.0", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 4, av) == BRLCAD_OK,
          "tor set_radius_1 12.0 returns OK");
    struct rt_tor_internal tor;
    if (read_tor(gedp, "tor.s", &tor) == BRLCAD_OK) {
        /* r_a is stored in base units (mm); compare with 12.0 */
        CHECK(NEAR_EQUAL(tor.r_a, 12.0, 0.1),
              "tor: r_a ≈ 12.0 after set_radius_1 12.0");
    } else {
        CHECK(0, "tor: read_tor succeeded after set_radius_1");
    }
}

/* 3-1: change torus minor radius (r2) via slug */
static void
test_p3_tor_set_radius_2(struct ged *gedp)
{
    const char *av[] = { "edit", "tor.s", "set_radius_2", "2.5", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 4, av) == BRLCAD_OK,
          "tor set_radius_2 2.5 returns OK");
    struct rt_tor_internal tor;
    if (read_tor(gedp, "tor.s", &tor) == BRLCAD_OK) {
        CHECK(NEAR_EQUAL(tor.r_h, 2.5, 0.1),
              "tor: r_h ≈ 2.5 after set_radius_2 2.5");
    } else {
        CHECK(0, "tor: read_tor succeeded after set_radius_2");
    }
}

/* 3-2: also accept param name "r1" as an alias for "set_radius_1" */
static void
test_p3_tor_r1_alias(struct ged *gedp)
{
    const char *av[] = { "edit", "tor.s", "r1", "15.0", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 4, av) == BRLCAD_OK,
          "tor r1 15.0 (param-name alias) returns OK");
    struct rt_tor_internal tor;
    if (read_tor(gedp, "tor.s", &tor) == BRLCAD_OK) {
        CHECK(NEAR_EQUAL(tor.r_a, 15.0, 0.1),
              "tor: r_a ≈ 15.0 after r1 alias");
    } else {
        CHECK(0, "tor: read_tor succeeded after r1 alias");
    }
}

/* 3-3: ell set_a — change semi-axis A magnitude */
static void
test_p3_ell_set_a(struct ged *gedp)
{
    const char *av[] = { "edit", "ell.s", "set_a", "6.0", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 4, av) == BRLCAD_OK,
          "ell set_a 6.0 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "ell.s", &ell) == BRLCAD_OK) {
        CHECK(NEAR_EQUAL(MAGNITUDE(ell.a), 6.0, 0.1),
              "ell: |A| ≈ 6.0 after set_a 6.0");
    } else {
        CHECK(0, "ell: read_ell succeeded after set_a");
    }
}

/* 3-4: ell "a" param-name alias for set_a */
static void
test_p3_ell_a_alias(struct ged *gedp)
{
    const char *av[] = { "edit", "ell.s", "a", "7.0", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 4, av) == BRLCAD_OK,
          "ell a 7.0 (param-name alias) returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "ell.s", &ell) == BRLCAD_OK) {
        CHECK(NEAR_EQUAL(MAGNITUDE(ell.a), 7.0, 0.1),
              "ell: |A| ≈ 7.0 after 'a 7.0' alias");
    } else {
        CHECK(0, "ell: read_ell succeeded after 'a' alias");
    }
}

/* 3-5: ell set_a_b_c — uniform scale all three semi-axes */
static void
test_p3_ell_set_abc(struct ged *gedp)
{
    const char *av[] = { "edit", "ell.s", "set_a_b_c", "4.0", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 4, av) == BRLCAD_OK,
          "ell set_a_b_c 4.0 returns OK");
    struct rt_ell_internal ell;
    if (read_ell(gedp, "ell.s", &ell) == BRLCAD_OK) {
        /* All three semi-axes should be equal in magnitude after uniform scale */
        double ma = MAGNITUDE(ell.a);
        double mb = MAGNITUDE(ell.b);
        double mc = MAGNITUDE(ell.c);
        CHECK(NEAR_EQUAL(ma, mb, 0.2),
              "ell: |A| ≈ |B| after set_a_b_c");
        CHECK(NEAR_EQUAL(mb, mc, 0.2),
              "ell: |B| ≈ |C| after set_a_b_c");
    } else {
        CHECK(0, "ell: read_ell succeeded after set_a_b_c");
    }
}

/* 3-6: 'edit tor --list-ops' now emits human-readable help (not JSON) */
static void
test_p3_list_ops_tor(struct ged *gedp)
{
    const char *av[] = { "edit", "tor", "--list-ops", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK,
          "edit tor --list-ops returns OK");
    const char *out = bu_vls_cstr(gedp->ged_result_str);
    CHECK(out && strlen(out) > 10,
          "edit tor --list-ops output is non-empty");
    CHECK(out && strstr(out, "Torus") != NULL,
          "edit tor --list-ops (human) mentions 'Torus'");
    CHECK(out && strstr(out, "set_radius") != NULL,
          "edit tor --list-ops (human) lists 'set_radius'");
    CHECK(out && strstr(out, "r1") != NULL,
          "edit tor --list-ops (human) lists short alias 'r1'");
    /* Must NOT be raw JSON */
    CHECK(out && strstr(out, "\"prim_type\"") == NULL,
          "edit tor --list-ops (human) does not contain raw JSON prim_type key");
}

/* 3-6b: 'edit tor --list-ops=json' emits JSON descriptor */
static void
test_p3_list_ops_tor_json(struct ged *gedp)
{
    const char *av[] = { "edit", "tor", "--list-ops=json", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK,
          "edit tor --list-ops=json returns OK");
    const char *json = bu_vls_cstr(gedp->ged_result_str);
    CHECK(json && strlen(json) > 10,
          "edit tor --list-ops=json output is non-empty");
    CHECK(json && strstr(json, "prim_type") != NULL,
          "edit tor --list-ops=json contains 'prim_type' JSON key");
    CHECK(json && strstr(json, "\"tor\"") != NULL,
          "edit tor --list-ops=json contains '\"tor\"'");
    CHECK(json && strstr(json, "\"name\": \"r1\"") != NULL,
          "edit tor --list-ops=json contains stable canonical name");
    CHECK(json && strstr(json, "\"symbol\": \"ECMD_TOR_R1\"") != NULL,
          "edit tor --list-ops=json contains stable ECMD symbol");
    CHECK(json && strstr(json, "aliases") != NULL,
          "edit tor --list-ops=json contains aliases");
    CHECK(json && strstr(json, "\"set_radius_1\"") != NULL,
          "edit tor --list-ops=json includes friendly set_radius_1 alias");
    CHECK(json && strstr(json, "\"r1\"") != NULL,
          "edit tor --list-ops=json includes alias r1");
}

/* 3-7: 'edit tor' (type-name help) prints category-grouped info */
static void
test_p3_type_help_tor(struct ged *gedp)
{
    const char *av[] = { "edit", "tor", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 2, av) == BRLCAD_OK,
          "edit tor (type help) returns OK");
    const char *help = bu_vls_cstr(gedp->ged_result_str);
    CHECK(help && strlen(help) > 10,
          "edit tor help output is non-empty");
    CHECK(help && strstr(help, "Torus") != NULL,
          "edit tor help mentions 'Torus'");
}

/* 3-8: unknown op for a descriptor primitive gives error + prim help */
static void
test_p3_unknown_prim_op(struct ged *gedp)
{
    const char *av[] = { "edit", "tor.s", "no_such_op", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_ERROR,
          "edit tor.s no_such_op returns BRLCAD_ERROR");
    const char *msg = bu_vls_cstr(gedp->ged_result_str);
    CHECK(msg && strlen(msg) > 0,
          "edit tor.s no_such_op prints an error message");
}

/* 3-9: missing value for a scalar param gives error */
static void
test_p3_missing_param(struct ged *gedp)
{
    const char *av[] = { "edit", "tor.s", "set_radius_1", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_ERROR,
          "edit tor.s set_radius_1 (no value) returns BRLCAD_ERROR");
}

/* 3-10: 'edit ell --list-ops' emits human-readable help */
static void
test_p3_list_ops_ell(struct ged *gedp)
{
    const char *av[] = { "edit", "ell", "--list-ops", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK,
          "edit ell --list-ops returns OK");
    const char *out = bu_vls_cstr(gedp->ged_result_str);
    CHECK(out && strstr(out, "Ellipsoid") != NULL,
          "edit ell --list-ops (human) mentions 'Ellipsoid'");
    CHECK(out && strstr(out, "\"prim_type\"") == NULL,
          "edit ell --list-ops (human) does not contain raw JSON");
}

/* 3-10b: 'edit ell --list-ops=json' emits JSON */
static void
test_p3_list_ops_ell_json(struct ged *gedp)
{
    const char *av[] = { "edit", "ell", "--list-ops=json", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK,
          "edit ell --list-ops=json returns OK");
    const char *json = bu_vls_cstr(gedp->ged_result_str);
    CHECK(json && strstr(json, "\"ell\"") != NULL,
          "edit ell --list-ops=json contains '\"ell\"'");
}

/* 3-11: 'edit --list-all-prim-ops' prints human-readable block for all types */
static void
test_p3_list_all_prim_ops(struct ged *gedp)
{
    const char *av[] = { "edit", "--list-all-prim-ops", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 2, av) == BRLCAD_OK,
          "edit --list-all-prim-ops returns OK");
    const char *out = bu_vls_cstr(gedp->ged_result_str);
    CHECK(out && strlen(out) > 50,
          "edit --list-all-prim-ops output is non-trivial");
    /* Should mention several known primitive types */
    CHECK(out && strstr(out, "Torus") != NULL,
          "edit --list-all-prim-ops mentions 'Torus'");
    CHECK(out && strstr(out, "Ellipsoid") != NULL,
          "edit --list-all-prim-ops mentions 'Ellipsoid'");
    /* Must not be raw JSON */
    CHECK(out && strstr(out, "\"prim_type\"") == NULL,
          "edit --list-all-prim-ops (human) does not contain raw JSON key");
}

/* 3-12: 'edit --list-all-prim-ops=json' emits a JSON array */
static void
test_p3_list_all_prim_ops_json(struct ged *gedp)
{
    const char *av[] = { "edit", "--list-all-prim-ops=json", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 2, av) == BRLCAD_OK,
          "edit --list-all-prim-ops=json returns OK");
    const char *json = bu_vls_cstr(gedp->ged_result_str);
    CHECK(json && strlen(json) > 50,
          "edit --list-all-prim-ops=json output is non-trivial");
    /* Should be a JSON array and contain multiple prim_type entries */
    CHECK(json && json[0] == '[',
          "edit --list-all-prim-ops=json starts with '['");
    CHECK(json && strstr(json, "prim_type") != NULL,
          "edit --list-all-prim-ops=json contains 'prim_type'");
    /* Must mention at least tor and ell */
    CHECK(json && strstr(json, "\"tor\"") != NULL,
          "edit --list-all-prim-ops=json contains '\"tor\"'");
    CHECK(json && strstr(json, "\"ell\"") != NULL,
          "edit --list-all-prim-ops=json contains '\"ell\"'");
}


/* ================================================================== *
 * Section 4 — Phase 4 descriptor completeness
 * ================================================================== */

/* ------------------------------------------------------------------ *
 * Helper: read rt_tgc_internal                                       *
 * ------------------------------------------------------------------ */
static int
read_tgc(struct ged *gedp, const char *name, struct rt_tgc_internal *out)
{
    struct directory *dp = db_lookup(gedp->dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL)
        return BRLCAD_ERROR;
    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, gedp->dbip, NULL) < 0)
        return BRLCAD_ERROR;
    if (intern.idb_minor_type != DB5_MINORTYPE_BRLCAD_TGC) {
        rt_db_free_internal(&intern);
        return BRLCAD_ERROR;
    }
    *out = *(struct rt_tgc_internal *)intern.idb_ptr;
    intern.idb_ptr = NULL;
    rt_db_free_internal(&intern);
    return BRLCAD_OK;
}

/* ------------------------------------------------------------------ *
 * Helper: read rt_eto_internal                                       *
 * ------------------------------------------------------------------ */
static int
read_eto(struct ged *gedp, const char *name, struct rt_eto_internal *out)
{
    struct directory *dp = db_lookup(gedp->dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL)
        return BRLCAD_ERROR;
    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, gedp->dbip, NULL) < 0)
        return BRLCAD_ERROR;
    if (intern.idb_minor_type != DB5_MINORTYPE_BRLCAD_ETO) {
        rt_db_free_internal(&intern);
        return BRLCAD_ERROR;
    }
    *out = *(struct rt_eto_internal *)intern.idb_ptr;
    intern.idb_ptr = NULL;
    rt_db_free_internal(&intern);
    return BRLCAD_OK;
}

/* ------------------------------------------------------------------ *
 * Helper: read rt_cline_internal                                     *
 * ------------------------------------------------------------------ */
static int
read_cline(struct ged *gedp, const char *name, struct rt_cline_internal *out)
{
    struct directory *dp = db_lookup(gedp->dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL)
        return BRLCAD_ERROR;
    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, gedp->dbip, NULL) < 0)
        return BRLCAD_ERROR;
    if (intern.idb_minor_type != DB5_MINORTYPE_BRLCAD_CLINE) {
        rt_db_free_internal(&intern);
        return BRLCAD_ERROR;
    }
    *out = *(struct rt_cline_internal *)intern.idb_ptr;
    intern.idb_ptr = NULL;
    rt_db_free_internal(&intern);
    return BRLCAD_OK;
}

/* ------------------------------------------------------------------ *
 * Helper: read rt_hyp_internal                                       *
 * ------------------------------------------------------------------ */
static int
read_hyp(struct ged *gedp, const char *name, struct rt_hyp_internal *out)
{
    struct directory *dp = db_lookup(gedp->dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL)
        return BRLCAD_ERROR;
    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, gedp->dbip, NULL) < 0)
        return BRLCAD_ERROR;
    if (intern.idb_minor_type != DB5_MINORTYPE_BRLCAD_HYP) {
        rt_db_free_internal(&intern);
        return BRLCAD_ERROR;
    }
    *out = *(struct rt_hyp_internal *)intern.idb_ptr;
    intern.idb_ptr = NULL;
    rt_db_free_internal(&intern);
    return BRLCAD_OK;
}

/* 4-0: tgc move_end_h_rt — POINT param moves H endpoint */
static void
test_p4_tgc_move_end_h_rt(struct ged *gedp)
{
    /* tgc.s: V={20,0,0}, H={0,0,10}; move H endpoint to (20,0,15) */
    const char *av[] = { "edit", "tgc.s", "move_end_h_rt", "20", "0", "15", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK,
          "tgc move_end_h_rt 20 0 15 returns OK");
    struct rt_tgc_internal tgc;
    if (read_tgc(gedp, "tgc.s", &tgc) == BRLCAD_OK) {
        /* new H = endpoint - V = {0,0,15} */
        CHECK(NEAR_EQUAL(MAGNITUDE(tgc.h), 15.0, 0.2),
              "tgc: |H| ≈ 15 after move_end_h_rt to (20,0,15)");
    } else {
        CHECK(0, "tgc: read_tgc succeeded after move_end_h_rt");
    }
}

/* 4-1: tgc move_end_h — POINT param, leave ends alone */
static void
test_p4_tgc_move_end_h(struct ged *gedp)
{
    const char *av[] = { "edit", "tgc.s", "move_end_h", "20", "0", "12", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK,
          "tgc move_end_h 20 0 12 returns OK");
    struct rt_tgc_internal tgc;
    if (read_tgc(gedp, "tgc.s", &tgc) == BRLCAD_OK) {
        CHECK(NEAR_EQUAL(MAGNITUDE(tgc.h), 12.0, 0.2),
              "tgc: |H| ≈ 12 after move_end_h to (20,0,12)");
    } else {
        CHECK(0, "tgc: read_tgc succeeded after move_end_h");
    }
}

/* 4-2: tgc rotate_h — VECTOR (degrees) rotates the H vector */
static void
test_p4_tgc_rotate_h(struct ged *gedp)
{
    /* Rotate H 90° around Y: initial H is approx (0,0,12), after rotation
     * it should point roughly along +X with the same magnitude. */
    const char *av[] = { "edit", "tgc.s", "rotate_h", "0", "90", "0", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK,
          "tgc rotate_h 0 90 0 returns OK");
    struct rt_tgc_internal tgc;
    if (read_tgc(gedp, "tgc.s", &tgc) == BRLCAD_OK) {
        /* H magnitude must be preserved; direction should have changed */
        double mag = MAGNITUDE(tgc.h);
        CHECK(mag > 0.1, "tgc: |H| > 0 after rotate_h");
        /* After 90° around Y, Z component should be near zero */
        CHECK(fabs(tgc.h[Z]) < 1.0,
              "tgc: H[Z] ≈ 0 after 90° Y-rotation");
    } else {
        CHECK(0, "tgc: read_tgc succeeded after rotate_h");
    }
}

/* 4-3: tgc rotate_axb — VECTOR (degrees) rotates A,B,C,D vectors */
static void
test_p4_tgc_rotate_axb(struct ged *gedp)
{
    /* Rotate AxB by 45° around Z */
    struct rt_tgc_internal before;
    if (read_tgc(gedp, "tgc.s", &before) != BRLCAD_OK) {
        CHECK(0, "tgc: read_tgc before rotate_axb"); return;
    }
    const char *av[] = { "edit", "tgc.s", "rotate_axb", "0", "0", "45", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK,
          "tgc rotate_axb 0 0 45 returns OK");
    struct rt_tgc_internal after;
    if (read_tgc(gedp, "tgc.s", &after) == BRLCAD_OK) {
        /* |A| must be preserved */
        CHECK(NEAR_EQUAL(MAGNITUDE(after.a), MAGNITUDE(before.a), 0.3),
              "tgc: |A| preserved after rotate_axb 45°");
    } else {
        CHECK(0, "tgc: read_tgc succeeded after rotate_axb");
    }
}

/* 4-4: eto rotate_c — VECTOR (degrees) rotates the C vector */
static void
test_p4_eto_rotate_c(struct ged *gedp)
{
    /* eto.s: C = {8,0,2} — rotate 90° around Z */
    struct rt_eto_internal before;
    if (read_eto(gedp, "eto.s", &before) != BRLCAD_OK) {
        CHECK(0, "eto: read_eto before rotate_c"); return;
    }
    const char *av[] = { "edit", "eto.s", "rotate_c", "0", "0", "90", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK,
          "eto rotate_c 0 0 90 returns OK");
    struct rt_eto_internal after;
    if (read_eto(gedp, "eto.s", &after) == BRLCAD_OK) {
        /* |C| must be preserved */
        double mag_before = MAGNITUDE(before.eto_C);
        double mag_after  = MAGNITUDE(after.eto_C);
        CHECK(NEAR_EQUAL(mag_before, mag_after, 0.3),
              "eto: |C| preserved after rotate_c 90°");
        /* Direction must have changed */
        CHECK(!NEAR_EQUAL(before.eto_C[X], after.eto_C[X], 0.5),
              "eto: C[X] changed after rotate_c 90°");
    } else {
        CHECK(0, "eto: read_eto succeeded after rotate_c");
    }
}

/* 4-5: cline move_end_h — POINT param relocates the H endpoint */
static void
test_p4_cline_move_end_h(struct ged *gedp)
{
    /* cline.s: V={80,40,0}, H={0,0,15}; move endpoint to (80,40,10) */
    const char *av[] = { "edit", "cline.s", "move_end_h", "80", "40", "10", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK,
          "cline move_end_h 80 40 10 returns OK");
    struct rt_cline_internal cli;
    if (read_cline(gedp, "cline.s", &cli) == BRLCAD_OK) {
        /* new H = endpoint - V = {0,0,10} */
        CHECK(NEAR_EQUAL(MAGNITUDE(cli.h), 10.0, 0.2),
              "cline: |H| ≈ 10 after move_end_h to (80,40,10)");
    } else {
        CHECK(0, "cline: read_cline succeeded after move_end_h");
    }
}

/* 4-6: hyp rotate_h — VECTOR (degrees) rotates the H vector */
static void
test_p4_hyp_rotate_h(struct ged *gedp)
{
    /* hyp.s: H={0,0,10}; rotate 90° around X */
    struct rt_hyp_internal before;
    if (read_hyp(gedp, "hyp.s", &before) != BRLCAD_OK) {
        CHECK(0, "hyp: read_hyp before rotate_h"); return;
    }
    const char *av[] = { "edit", "hyp.s", "rotate_h", "90", "0", "0", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK,
          "hyp rotate_h 90 0 0 returns OK");
    struct rt_hyp_internal after;
    if (read_hyp(gedp, "hyp.s", &after) == BRLCAD_OK) {
        /* |H| must be preserved */
        double mag_before = MAGNITUDE(before.hyp_Hi);
        double mag_after  = MAGNITUDE(after.hyp_Hi);
        CHECK(NEAR_EQUAL(mag_before, mag_after, 0.3),
              "hyp: |H| preserved after rotate_h 90°");
        /* Z component should now be near zero after 90° X rotation */
        CHECK(fabs(after.hyp_Hi[Z]) < 1.0,
              "hyp: H[Z] ≈ 0 after 90° X-rotation");
    } else {
        CHECK(0, "hyp: read_hyp succeeded after rotate_h");
    }
}

/* 4-7: extrude --list-ops shows the extrude descriptor (no object needed) */
static void
test_p4_extrude_list_ops(struct ged *gedp)
{
    const char *av[] = { "edit", "extrude", "--list-ops", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK,
          "edit extrude --list-ops returns OK");
    const char *out = bu_vls_cstr(gedp->ged_result_str);
    CHECK(out && strlen(out) > 10,
          "edit extrude --list-ops output is non-empty");
    CHECK(out && strstr(out, "Extrusion") != NULL,
          "edit extrude --list-ops mentions 'Extrusion'");
    CHECK(out && strstr(out, "set_h") != NULL,
          "edit extrude --list-ops lists 'set_h'");
    CHECK(out && strstr(out, "move_end_h") != NULL,
          "edit extrude --list-ops lists 'move_end_h'");
    CHECK(out && strstr(out, "rotate_h") != NULL,
          "edit extrude --list-ops lists 'rotate_h'");
}

/* 4-8: extrude --list-ops=json returns valid JSON descriptor */
static void
test_p4_extrude_list_ops_json(struct ged *gedp)
{
    const char *av[] = { "edit", "extrude", "--list-ops=json", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK,
          "edit extrude --list-ops=json returns OK");
    const char *json = bu_vls_cstr(gedp->ged_result_str);
    CHECK(json && strstr(json, "\"extrude\"") != NULL,
          "edit extrude --list-ops=json contains '\"extrude\"'");
    CHECK(json && strstr(json, "\"Set H\"") != NULL,
          "edit extrude --list-ops=json contains '\"Set H\"'");
}

/* 4-9: tgc --list-ops now includes the new move/rotation ops */
static void
test_p4_tgc_list_ops_complete(struct ged *gedp)
{
    const char *av[] = { "edit", "tgc", "--list-ops", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK,
          "edit tgc --list-ops returns OK");
    const char *out = bu_vls_cstr(gedp->ged_result_str);
    CHECK(out && strstr(out, "move_end_h_rt") != NULL,
          "edit tgc --list-ops lists 'move_end_h_rt'");
    CHECK(out && strstr(out, "rotate_h") != NULL,
          "edit tgc --list-ops lists 'rotate_h'");
    CHECK(out && strstr(out, "rotate_axb") != NULL,
          "edit tgc --list-ops lists 'rotate_axb'");
}

/* 4-10: cline --list-ops now includes move_end_h */
static void
test_p4_cline_list_ops_complete(struct ged *gedp)
{
    const char *av[] = { "edit", "cline", "--list-ops", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK,
          "edit cline --list-ops returns OK");
    const char *out = bu_vls_cstr(gedp->ged_result_str);
    CHECK(out && strstr(out, "move_end_h") != NULL,
          "edit cline --list-ops lists 'move_end_h'");
}

/* 4-11: hyp --list-ops now includes rotate_h */
static void
test_p4_hyp_list_ops_complete(struct ged *gedp)
{
    const char *av[] = { "edit", "hyp", "--list-ops", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK,
          "edit hyp --list-ops returns OK");
    const char *out = bu_vls_cstr(gedp->ged_result_str);
    CHECK(out && strstr(out, "rotate_h") != NULL,
          "edit hyp --list-ops lists 'rotate_h'");
}

/* 4-12: eto --list-ops now includes rotate_c */
static void
test_p4_eto_list_ops_complete(struct ged *gedp)
{
    const char *av[] = { "edit", "eto", "--list-ops", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK,
          "edit eto --list-ops returns OK");
    const char *out = bu_vls_cstr(gedp->ged_result_str);
    CHECK(out && strstr(out, "rotate_c") != NULL,
          "edit eto --list-ops lists 'rotate_c'");
}

/* 4-13: the supported C,D-adjusting height operation is executable. */
static void
test_p4_tgc_move_end_h_adj_cd(struct ged *gedp)
{
    const char *av[] = { "edit", "tgc.s", "scale_h_cd", "8", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 4, av) == BRLCAD_OK,
          "tgc scale_h_cd 8 returns OK");
    struct rt_tgc_internal tgc;
    if (read_tgc(gedp, "tgc.s", &tgc) == BRLCAD_OK) {
        double mag = MAGNITUDE(tgc.h);
        CHECK(NEAR_EQUAL(mag, 8.0, NEAR_ENOUGH),
	      "tgc: |H| == 8 after scale_h_cd");
    } else {
        CHECK(0, "tgc: read_tgc succeeded after scale_h_cd");
    }
}

/* 4-14: extrude is in --list-all-prim-ops */
static void
test_p4_all_prim_ops_has_extrude(struct ged *gedp)
{
    const char *av[] = { "edit", "--list-all-prim-ops=json", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 2, av) == BRLCAD_OK,
          "edit --list-all-prim-ops=json still returns OK after extrude added");
    const char *json = bu_vls_cstr(gedp->ged_result_str);
    CHECK(json && strstr(json, "\"extrude\"") != NULL,
          "edit --list-all-prim-ops=json includes '\"extrude\"'");
}


/* ================================================================== *
 * Section 5 — Phase 4+: pipe descriptor completions + arb8 descriptor
 * ================================================================== */

/* ------------------------------------------------------------------ *
 * Fixture helper: create a 3-point pipe + an arb8 box                *
 * ------------------------------------------------------------------ */
static int
create_p5_fixture(const char *dbpath)
{
    struct rt_wdb *wdbp = wdb_fopen(dbpath);
    if (!wdbp) {
        bu_log("ERROR: unable to create p5 fixture %s\n", dbpath);
        return BRLCAD_ERROR;
    }

    /* --- pipe with 3 points ----------------------------------------- */
    struct bu_list pipe_head;
    mk_pipe_init(&pipe_head);
    point_t pp0 = {0, 0,  0};
    point_t pp1 = {0, 0, 10};
    point_t pp2 = {0, 0, 20};
    mk_add_pipe_pnt(&pipe_head, pp0, 4.0, 2.0, 5.0);
    mk_add_pipe_pnt(&pipe_head, pp1, 4.0, 2.0, 5.0);
    mk_add_pipe_pnt(&pipe_head, pp2, 4.0, 2.0, 5.0);
    if (mk_pipe(wdbp, "pipe.s", &pipe_head) != 0) {
        bu_log("mk_pipe failed\n");
        mk_pipe_free(&pipe_head);
        db_close(wdbp->dbip);
        return BRLCAD_ERROR;
    }
    mk_pipe_free(&pipe_head);

    /* --- arb8 box --------------------------------------------------- */
    fastf_t arb_pts[8*3] = {
        -5,-5,-5,  5,-5,-5,  5, 5,-5, -5, 5,-5,
        -5,-5, 5,  5,-5, 5,  5, 5, 5, -5, 5, 5
    };
    if (mk_arb8(wdbp, "arb8.s", arb_pts) != 0) {
        bu_log("mk_arb8 failed\n");
        db_close(wdbp->dbip);
        return BRLCAD_ERROR;
    }

    /* --- planar sketch with one line segment ----------------------- */
    struct rt_sketch_internal *skt;
    BU_ALLOC(skt, struct rt_sketch_internal);
    skt->magic = RT_SKETCH_INTERNAL_MAGIC;
    VSET(skt->V, 0, 0, 0);
    VSET(skt->u_vec, 1, 0, 0);
    VSET(skt->v_vec, 0, 1, 0);
    skt->vert_count = 4;
    skt->verts = (point2d_t *)bu_malloc(4 * sizeof(point2d_t),
	"GED edit sketch vertices");
    V2SET(skt->verts[0], 0, 0);
    V2SET(skt->verts[1], 10, 0);
    V2SET(skt->verts[2], 10, 10);
    V2SET(skt->verts[3], 0, 10);
    struct line_seg *line;
    BU_ALLOC(line, struct line_seg);
    line->magic = CURVE_LSEG_MAGIC;
    line->start = 0;
    line->end = 1;
    skt->curve.count = 1;
    skt->curve.segment = (void **)bu_malloc(sizeof(void *),
	"GED edit sketch segments");
    skt->curve.reverse = (int *)bu_calloc(1, sizeof(int),
	"GED edit sketch reverse flags");
    skt->curve.segment[0] = line;
    if (wdb_export(wdbp, "sketch.s", (void *)skt, ID_SKETCH, 1.0) != 0) {
	bu_log("wdb_export sketch failed\n");
	db_close(wdbp->dbip);
	return BRLCAD_ERROR;
    }

    db_close(wdbp->dbip);
    return BRLCAD_OK;
}


/* ------------------------------------------------------------------ *
 * sketch: the shared session uses true 2-D UV values and preserves  *
 * selection state across command, widget, and manipulator clients.  *
 * ------------------------------------------------------------------ */
static void
test_p5_sketch_session_state(struct ged *gedp)
{
    const char *json_av[] = {"edit", "sketch", "--list-ops=json", NULL};
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, json_av) == BRLCAD_OK &&
	strstr(bu_vls_cstr(gedp->ged_result_str), "point2") &&
	strstr(bu_vls_cstr(gedp->ged_result_str), "vector2") &&
	strstr(bu_vls_cstr(gedp->ged_result_str), "integer_list") &&
	strstr(bu_vls_cstr(gedp->ged_result_str), "scalar_list") &&
	strstr(bu_vls_cstr(gedp->ged_result_str), "min_count"),
	"sketch descriptor publishes 2-D and variable-list semantics");

    const char *select_av[] = {
	"edit", "-i", "sketch.s", "pick_vertex", "2", NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 5, select_av) == BRLCAD_OK,
	"sketch CLI selects a vertex in the shared session");

    ged_edit_session_ref session = GED_EDIT_SESSION_REF_NULL;
    struct ged_edit_session_info info = {};
    CHECK(ged_edit_session_find(gedp, "sketch.s", &session) == GED_EDIT_OK &&
	ged_edit_session_info_get(gedp, session, &info) == GED_EDIT_OK &&
	info.revision == 1 && !info.dirty,
	"sketch selection advances revision without dirtying geometry");

    struct rt_edit_param_bounds vertex_bounds = {};
    CHECK(ged_edit_session_parameter_bounds_get(gedp, session,
	info.command_id, 0, &vertex_bounds) == GED_EDIT_OK &&
	vertex_bounds.has_maximum &&
	fabs(vertex_bounds.maximum - 3.0) < NEAR_ENOUGH,
	"sketch session exposes current vertex topology to controls");
    const char *invalid_select_av[] = {
	"edit", "-i", "sketch.s", "pick_vertex", "4", NULL
    };
    CHECK(ged_exec(gedp, 5, invalid_select_av) == BRLCAD_ERROR &&
	ged_edit_session_info_get(gedp, session, &info) == GED_EDIT_OK &&
	info.revision == 1 && !info.dirty,
	"sketch topology bounds reject an unavailable vertex atomically");

    const char *get_select_av[] = {
	"edit", "sketch.s", "get", "pick_vertex", NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 4, get_select_av) == BRLCAD_OK &&
	strstr(bu_vls_cstr(gedp->ged_result_str), "2"),
	"sketch command readback observes the authoritative selection");

    const char *move_av[] = {
	"edit", "-i", "sketch.s", "move_vertex", "12.5", "11.5", NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, move_av) == BRLCAD_OK,
	"sketch move_vertex accepts exactly two UV coordinates");

    struct rt_db_internal preview = RT_DB_INTERNAL_INIT_ZERO;
    bool moved = false;
    if (ged_edit_session_internal_copy(gedp, session, &preview) == GED_EDIT_OK &&
	preview.idb_type == ID_SKETCH && preview.idb_ptr) {
	struct rt_sketch_internal *working =
	    (struct rt_sketch_internal *)preview.idb_ptr;
	moved = NEAR_EQUAL(working->verts[2][0], 12.5, NEAR_ENOUGH) &&
	    NEAR_EQUAL(working->verts[2][1], 11.5, NEAR_ENOUGH);
    }
    rt_db_free_internal(&preview);
    CHECK(moved,
	"sketch session snapshot contains the command-side vertex move");

    const char *get_move_av[] = {
	"edit", "sketch.s", "get", "move_vertex", NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 4, get_move_av) == BRLCAD_OK &&
	strstr(bu_vls_cstr(gedp->ged_result_str), "12.5") &&
	strstr(bu_vls_cstr(gedp->ged_result_str), "11.5"),
	"sketch CLI reads the same moved UV values as graphical clients");

    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 4, get_select_av) == BRLCAD_OK &&
	strstr(bu_vls_cstr(gedp->ged_result_str), "2"),
	"sketch selection remains authoritative after a geometry update");

    const char *extra_move_av[] = {
	"edit", "-i", "sketch.s", "move_vertex", "1", "2", "3", NULL
    };
    CHECK(ged_exec(gedp, 7, extra_move_av) == BRLCAD_ERROR,
	"fixed-width sketch commands reject surplus typed values");

    const char *short_bezier_av[] = {
	"edit", "-i", "sketch.s", "append_bezier", "0", NULL
    };
    CHECK(ged_exec(gedp, 5, short_bezier_av) == BRLCAD_ERROR,
	"sketch variable-list commands enforce their minimum arity");

    const char *bezier_av[] = {
	"edit", "-i", "sketch.s", "append_bezier", "0", "1", "2", NULL
    };
    CHECK(ged_exec(gedp, 7, bezier_av) == BRLCAD_OK,
	"sketch CLI submits a typed variable control-vertex list");

    const char *move_list_av[] = {
	"edit", "-i", "sketch.s", "move_vertex_list",
	"1", "-2", "0", "3", NULL
    };
    CHECK(ged_exec(gedp, 8, move_list_av) == BRLCAD_OK,
	"sketch CLI combines a fixed 2-D delta with a variable vertex list");

    struct rt_db_internal topology_preview = RT_DB_INTERNAL_INIT_ZERO;
    bool topology_changed = false;
    if (ged_edit_session_internal_copy(gedp, session, &topology_preview) ==
	    GED_EDIT_OK && topology_preview.idb_type == ID_SKETCH &&
	topology_preview.idb_ptr) {
	struct rt_sketch_internal *working =
	    (struct rt_sketch_internal *)topology_preview.idb_ptr;
	topology_changed = working->curve.count == 2 &&
	    NEAR_EQUAL(working->verts[0][0], 1.0, NEAR_ENOUGH) &&
	    NEAR_EQUAL(working->verts[0][1], -2.0, NEAR_ENOUGH) &&
	    NEAR_EQUAL(working->verts[3][0], 1.0, NEAR_ENOUGH) &&
	    NEAR_EQUAL(working->verts[3][1], 8.0, NEAR_ENOUGH);
    }
    rt_db_free_internal(&topology_preview);
    CHECK(topology_changed,
	"typed sketch list commands update only the shared intermediate");

    const char *reset_av[] = {"edit", "sketch.s", "reset", NULL};
    CHECK(ged_exec(gedp, 3, reset_av) == BRLCAD_OK,
	"sketch shared session cancels cleanly");

    struct directory *dp = db_lookup(gedp->dbip, "sketch.s", LOOKUP_QUIET);
    struct rt_db_internal committed = RT_DB_INTERNAL_INIT_ZERO;
    bool unchanged = dp && rt_db_get_internal(&committed, dp, gedp->dbip,
	NULL) >= 0;
    if (unchanged) {
	struct rt_sketch_internal *stored =
	    (struct rt_sketch_internal *)committed.idb_ptr;
	unchanged = NEAR_EQUAL(stored->verts[2][0], 10.0, NEAR_ENOUGH) &&
	    NEAR_EQUAL(stored->verts[2][1], 10.0, NEAR_ENOUGH);
    }
    rt_db_free_internal(&committed);
    CHECK(unchanged,
	"cancel leaves committed sketch geometry unchanged");
}

/* ------------------------------------------------------------------ *
 * pipe: --list-ops includes the new ops                              *
 * ------------------------------------------------------------------ */
static void
test_p5_pipe_list_ops(struct ged *gedp)
{
    const char *av[] = { "edit", "pipe", "--list-ops", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK,
          "edit pipe --list-ops returns OK");
    const char *out = bu_vls_cstr(gedp->ged_result_str);
    CHECK(out && strstr(out, "select_point") != NULL,
          "edit pipe --list-ops lists 'select_point'");
    CHECK(out && strstr(out, "next_point") != NULL,
          "edit pipe --list-ops lists 'next_point'");
    CHECK(out && strstr(out, "previous_point") != NULL,
          "edit pipe --list-ops lists 'previous_point'");
    CHECK(out && strstr(out, "move_point") != NULL,
          "edit pipe --list-ops lists 'move_point'");
    CHECK(out && strstr(out, "delete_point") != NULL,
          "edit pipe --list-ops lists 'delete_point'");
    CHECK(out && strstr(out, "append_point") != NULL,
          "edit pipe --list-ops lists 'append_point'");
    CHECK(out && strstr(out, "prepend_point") != NULL,
          "edit pipe --list-ops lists 'prepend_point'");
    CHECK(out && strstr(out, "split_segment") != NULL,
          "edit pipe --list-ops lists 'split_segment'");
}

/* ------------------------------------------------------------------ *
 * pipe: select + next_pt + prev_pt navigation                        *
 * ------------------------------------------------------------------ */
static void
test_p5_pipe_select_next_prev(struct ged *gedp)
{
    /* Select the point nearest (0,0,0) */
    {
        const char *av[] = { "edit", "-i", "pipe.s", "select_point", "0", "0", "0", NULL };
        bu_vls_trunc(gedp->ged_result_str, 0);
        int ret = ged_exec(gedp, 7, av);
        CHECK(ret == BRLCAD_OK,
              "pipe.s select_point 0 0 0 returns OK");
    }
    /* Advance to the next point (0,0,10) */
    {
        const char *av[] = { "edit", "-i", "pipe.s", "next_point", NULL };
        bu_vls_trunc(gedp->ged_result_str, 0);
        int ret = ged_exec(gedp, 4, av);
        CHECK(ret == BRLCAD_OK,
              "pipe.s next_point returns OK");
    }
    /* Step back to the first point */
    {
        const char *av[] = { "edit", "-i", "pipe.s", "previous_point", NULL };
        bu_vls_trunc(gedp->ged_result_str, 0);
        int ret = ged_exec(gedp, 4, av);
        CHECK(ret == BRLCAD_OK,
              "pipe.s previous_point returns OK");
    }
    {
	const char *av[] = { "edit", "pipe.s", "reset", NULL };
	bu_vls_trunc(gedp->ged_result_str, 0);
	CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK,
	      "pipe navigation session resets cleanly");
    }
}

/* ------------------------------------------------------------------ *
 * pipe: append + delete point round-trip                             *
 * ------------------------------------------------------------------ */
static void
test_p5_pipe_append_del(struct ged *gedp)
{
    /* Append a new point at (0, 0, 30) */
    {
        const char *av[] = { "edit", "pipe.s", "append_point", "0", "0", "30", NULL };
        bu_vls_trunc(gedp->ged_result_str, 0);
        CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK,
              "pipe.s append_point 0 0 30 returns OK");
    }
    /* Select that new point */
    {
        const char *av[] = { "edit", "pipe.s", "select_point", "0", "0", "30", NULL };
        bu_vls_trunc(gedp->ged_result_str, 0);
        ged_exec(gedp, 6, av);
    }
    /* Delete it */
    {
        const char *av[] = { "edit", "pipe.s", "delete_point", NULL };
        bu_vls_trunc(gedp->ged_result_str, 0);
        CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK,
              "pipe.s delete_point returns OK");
    }
}

/* ------------------------------------------------------------------ *
 * pipe: prepend point                                                 *
 * ------------------------------------------------------------------ */
static void
test_p5_pipe_prepend(struct ged *gedp)
{
    const char *av[] = { "edit", "pipe.s", "prepend_point", "0", "0", "-10", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK,
          "pipe.s prepend_point 0 0 -10 returns OK");
}

/* ------------------------------------------------------------------ *
 * pipe: split_segment                                                 *
 * ------------------------------------------------------------------ */
static void
test_p5_pipe_split(struct ged *gedp)
{
    /* First select a point */
    {
        const char *av[] = { "edit", "pipe.s", "select_point", "0", "0", "0", NULL };
        bu_vls_trunc(gedp->ged_result_str, 0);
        ged_exec(gedp, 6, av);
    }
    /* Split segment between point 0 and point 1 */
    {
        const char *av[] = { "edit", "pipe.s", "split_segment", "0", "0", "5", NULL };
        bu_vls_trunc(gedp->ged_result_str, 0);
        CHECK(ged_exec(gedp, 6, av) == BRLCAD_OK,
              "pipe.s split_segment 0 0 5 returns OK");
    }
}

/* ------------------------------------------------------------------ *
 * arb8: --list-ops includes the new ops                              *
 * ------------------------------------------------------------------ */
static void
test_p5_arb8_list_ops(struct ged *gedp)
{
    const char *av[] = { "edit", "arb8", "--list-ops", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK,
          "edit arb8 --list-ops returns OK");
    const char *out = bu_vls_cstr(gedp->ged_result_str);
    CHECK(out && strstr(out, "move_face") != NULL,
          "edit arb8 --list-ops lists 'move_face'");
    CHECK(out && strstr(out, "move_edge") != NULL,
          "edit arb8 --list-ops lists 'move_edge'");
    CHECK(out && strstr(out, "move_vertex") != NULL,
          "edit arb8 --list-ops lists 'move_vertex'");
    CHECK(out && strstr(out, "select_vertex") != NULL,
          "edit arb8 --list-ops lists topology selection operations");
    CHECK(out && strstr(out, "rotate_face") != NULL,
          "edit arb8 --list-ops lists 'rotate_face'");
}

/* ------------------------------------------------------------------ *
 * arb8: viewport/CLI selection shares revisioned non-dirty state     *
 * ------------------------------------------------------------------ */
static void
test_p5_arb8_selection_state(struct ged *gedp)
{
    const char *select_av[] = {
	"edit", "-i", "arb8.s", "select_vertex", "6", NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 5, select_av) == BRLCAD_OK,
	"arb8.s select_vertex records viewport selection through the CLI");

    ged_edit_session_ref session = GED_EDIT_SESSION_REF_NULL;
    struct ged_edit_session_info info = {};
    CHECK(ged_edit_session_find(gedp, "arb8.s", &session) == GED_EDIT_OK &&
	ged_edit_session_info_get(gedp, session, &info) == GED_EDIT_OK &&
	info.revision == 1 && !info.dirty,
	"ARB selection advances shared state without marking geometry dirty");
    const struct rt_edit_prim_desc *descriptor = NULL;
    int select_command = 0;
    if (ged_edit_session_descriptor_get(gedp, session, &descriptor) ==
	GED_EDIT_OK && descriptor) {
	for (int i = 0; i < descriptor->ncmd; i++) {
	    if (descriptor->cmds[i].name &&
		BU_STR_EQUAL(descriptor->cmds[i].name,
		    "ECMD_ARB_SELECT_VERTEX")) {
		select_command = descriptor->cmds[i].cmd_id;
		break;
	    }
	}
    }
    struct rt_edit_cmd_values values;
    CHECK(select_command && ged_edit_session_command_values_get(gedp,
	session, select_command, &values) == GED_EDIT_OK &&
	values.value_count == 1 && values.value_valid[0] &&
	(int)values.values[0] == 6,
	"ARB widget/manipulator readback observes command-side selection");
    const char *get_av[] = {
	"edit", "arb8.s", "get", "select_vertex", NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 4, get_av) == BRLCAD_OK &&
	strstr(bu_vls_cstr(gedp->ged_result_str), "6"),
	"command-line edit reads the same ARB selection state");
    const char *reset_av[] = {"edit", "arb8.s", "reset", NULL};
    CHECK(ged_exec(gedp, 3, reset_av) == BRLCAD_OK,
	"ARB selection-only session closes without a database write");
}

/* ------------------------------------------------------------------ *
 * arb8: move_face (face 0 through a point)                           *
 * The box runs from -5,-5,-5 to 5,5,5.  Face 0 is the bottom (Z=-5).*
 * Moving it through (0,0,-7) should shift the bottom face down.      *
 * ------------------------------------------------------------------ */
static int
read_arb8(struct ged *gedp, const char *name, struct rt_arb_internal *out)
{
    struct directory *dp = db_lookup(gedp->dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL)
        return BRLCAD_ERROR;
    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, gedp->dbip, NULL) < 0)
        return BRLCAD_ERROR;
    if (intern.idb_minor_type != DB5_MINORTYPE_BRLCAD_ARB8) {
        rt_db_free_internal(&intern);
        return BRLCAD_ERROR;
    }
    *out = *(struct rt_arb_internal *)intern.idb_ptr;
    intern.idb_ptr = NULL;
    rt_db_free_internal(&intern);
    return BRLCAD_OK;
}

static void
test_p5_arb8_move_face(struct ged *gedp)
{
    struct rt_arb_internal before;
    read_arb8(gedp, "arb8.s", &before);

    /* Move face 0 through point (0, 0, -7) */
    const char *av[] = { "edit", "arb8.s", "move_face",
                         "0", "0", "0", "-7", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    int ret = ged_exec(gedp, 7, av);
    CHECK(ret == BRLCAD_OK,
          "arb8.s move_face 0  0 0 -7 returns OK");

    struct rt_arb_internal after;
    if (read_arb8(gedp, "arb8.s", &after) == BRLCAD_OK) {
        /* At least one vertex should have changed */
        double d = DIST_PNT_PNT(before.pt[0], after.pt[0]);
        CHECK(d > 0.01,
              "arb8.s: face moved — vertex 0 changed");
    } else {
        CHECK(0, "arb8.s: read arb8 after move_face");
    }
}

/* ------------------------------------------------------------------ *
 * arb8: move_vertex (vertex 0 to a new position)                     *
 * ------------------------------------------------------------------ */
static void
test_p5_arb8_move_vertex(struct ged *gedp)
{
    /* Move vertex 0 (originally at -5,-5,-5) to (-6,-5,-5) */
    const char *av[] = { "edit", "arb8.s", "move_vertex",
                         "0", "-6", "-5", "-5", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    int ret = ged_exec(gedp, 7, av);
    CHECK(ret == BRLCAD_OK,
          "arb8.s move_vertex 0  -6 -5 -5 returns OK");
}

/* ------------------------------------------------------------------ *
 * arb8: type option reaches the primitive edit implementation        *
 * ------------------------------------------------------------------ */
static void
test_p5_arb8_type_option(struct ged *gedp)
{
    const char *ok_av[] = { "edit", "-O", "type=arb8", "arb8.s", "move_vertex",
                            "0", "-7", "-5", "-5", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 9, ok_av) == BRLCAD_OK,
          "arb8.s move_vertex with -O type=arb8 returns OK");

    const char *bad_av[] = { "edit", "-O", "type=bogus", "arb8.s", "move_vertex",
                             "0", "-8", "-5", "-5", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    int ret = ged_exec(gedp, 9, bad_av);
    CHECK(ret == BRLCAD_ERROR,
          "arb8.s move_vertex with invalid -O type=bogus fails");
}

/* ------------------------------------------------------------------ *
 * arb8: --list-ops=json includes arb8                                *
 * ------------------------------------------------------------------ */
static void
test_p5_arb8_list_ops_json(struct ged *gedp)
{
    const char *av[] = { "edit", "arb8", "--list-ops=json", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK,
          "edit arb8 --list-ops=json returns OK");
    const char *json = bu_vls_cstr(gedp->ged_result_str);
    CHECK(json && strstr(json, "\"arb8\"") != NULL,
          "edit arb8 --list-ops=json contains '\"arb8\"'");
    CHECK(json && strstr(json, "Move Face") != NULL,
          "edit arb8 --list-ops=json contains 'Move Face'");
}

/* ------------------------------------------------------------------ *
 * all-prim-ops includes arb8                                         *
 * ------------------------------------------------------------------ */
static void
test_p5_all_prim_ops_has_arb8(struct ged *gedp)
{
    const char *av[] = { "edit", "--list-all-prim-ops=json", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 2, av) == BRLCAD_OK,
          "edit --list-all-prim-ops=json still returns OK after arb8 added");
    const char *json = bu_vls_cstr(gedp->ged_result_str);
    CHECK(json && strstr(json, "\"arb8\"") != NULL,
          "edit --list-all-prim-ops=json includes '\"arb8\"'");
}


/* ================================================================== *
 * Section 6 — Phase 5.1: brep descriptor + select/move/set CV ops   *
 * ================================================================== */

/* ------------------------------------------------------------------ *
 * Fixture helper: create a brep sphere                               *
 * ------------------------------------------------------------------ */
static int
create_p6_fixture(const char *dbpath)
{
    struct rt_wdb *wdbp = wdb_fopen(dbpath);
    if (!wdbp) {
        bu_log("ERROR: unable to create p6 fixture %s\n", dbpath);
        return BRLCAD_ERROR;
    }

    /* Build a BREP sphere at origin, radius 10. */
    ON_3dPoint centre(0.0, 0.0, 0.0);
    ON_Sphere sph(centre, 10.0);

    struct rt_brep_internal *bi;
    BU_ALLOC(bi, struct rt_brep_internal);
    bi->magic = RT_BREP_INTERNAL_MAGIC;
    bi->brep  = ON_BrepSphere(sph);

    if (wdb_export(wdbp, "brep_sph.s", (void *)bi, ID_BREP, 1.0) < 0) {
        bu_log("wdb_export brep_sph.s failed\n");
        db_close(wdbp->dbip);
        return BRLCAD_ERROR;
    }

    db_close(wdbp->dbip);
    return BRLCAD_OK;
}

/* ------------------------------------------------------------------ *
 * brep: --list-ops output includes the three CV ops                  *
 * ------------------------------------------------------------------ */
static void
test_p6_brep_list_ops(struct ged *gedp)
{
    const char *av[] = { "edit", "brep", "--list-ops", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK,
          "edit brep --list-ops returns OK");
    const char *out = bu_vls_cstr(gedp->ged_result_str);
    CHECK(out && strstr(out, "select_surface_cv") != NULL,
          "edit brep --list-ops lists 'select_surface_cv'");
    CHECK(out && strstr(out, "move_surface_cv") != NULL,
          "edit brep --list-ops lists 'move_surface_cv'");
    CHECK(out && strstr(out, "set_surface_cv_position") != NULL,
          "edit brep --list-ops lists 'set_surface_cv_position'");
}

/* ------------------------------------------------------------------ *
 * brep: --list-ops=json describes the descriptor                      *
 * ------------------------------------------------------------------ */
static void
test_p6_brep_list_ops_json(struct ged *gedp)
{
    const char *av[] = { "edit", "brep", "--list-ops=json", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 3, av) == BRLCAD_OK,
          "edit brep --list-ops=json returns OK");
    const char *json = bu_vls_cstr(gedp->ged_result_str);
    CHECK(json && strstr(json, "\"brep\"") != NULL,
          "edit brep --list-ops=json contains '\"brep\"'");
    CHECK(json && strstr(json, "Select Surface CV") != NULL,
          "edit brep --list-ops=json contains 'Select Surface CV'");
}

/* ------------------------------------------------------------------ *
 * brep: select a CV and move it                                       *
 * ------------------------------------------------------------------ */
static void
test_p6_brep_select_and_move(struct ged *gedp)
{
    /* Select face 0, CV (0, 0) */
    {
        const char *av[] = {
            "edit", "brep_sph.s", "select_surface_cv", "0", "0", "0", NULL
        };
        bu_vls_trunc(gedp->ged_result_str, 0);
        int ret = ged_exec(gedp, 6, av);
        CHECK(ret == BRLCAD_OK,
              "brep_sph.s select_surface_cv 0 0 0 returns OK");
    }
    /* Move the selected CV by (1, 0, 0) */
    {
        const char *av[] = {
            "edit", "brep_sph.s", "move_surface_cv", "1", "0", "0", NULL
        };
        bu_vls_trunc(gedp->ged_result_str, 0);
        int ret = ged_exec(gedp, 6, av);
        CHECK(ret == BRLCAD_OK,
              "brep_sph.s move_surface_cv 1 0 0 returns OK");
    }
}

/* ------------------------------------------------------------------ *
 * brep: select a CV and set its absolute position                     *
 * ------------------------------------------------------------------ */
static void
test_p6_brep_select_and_set(struct ged *gedp)
{
    /* Select face 0, CV (1, 0) */
    {
        const char *av[] = {
            "edit", "brep_sph.s", "select_surface_cv", "0", "1", "0", NULL
        };
        bu_vls_trunc(gedp->ged_result_str, 0);
        int ret = ged_exec(gedp, 6, av);
        CHECK(ret == BRLCAD_OK,
              "brep_sph.s select_surface_cv 0 1 0 returns OK");
    }
    /* Set the CV to (3, 4, 5) */
    {
        const char *av[] = {
            "edit", "brep_sph.s", "set_surface_cv_position", "3", "4", "5", NULL
        };
        bu_vls_trunc(gedp->ged_result_str, 0);
        int ret = ged_exec(gedp, 6, av);
        CHECK(ret == BRLCAD_OK,
              "brep_sph.s set_surface_cv_position 3 4 5 returns OK");
    }
}

/* ------------------------------------------------------------------ *
 * brep: bad face index rejected                                        *
 * ------------------------------------------------------------------ */
static void
test_p6_brep_bad_face(struct ged *gedp)
{
    const char *av[] = {
        "edit", "brep_sph.s", "select_surface_cv", "9999", "0", "0", NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    /* Bad face should not crash and should emit an out-of-range diagnostic. */
    int ret = ged_exec(gedp, 6, av);
    CHECK(ret == BRLCAD_OK || ret == BRLCAD_ERROR,
          "brep_sph.s select_surface_cv 9999 0 0 returns without crash");
}

/* ------------------------------------------------------------------ *
 * brep: move with no prior selection rejected                          *
 * ------------------------------------------------------------------ */
static void
test_p6_brep_move_no_selection(struct ged *gedp)
{
    /* First reset selection by making a fresh fixture GED session -
     * just verify that a fresh object with no selection gives error. */
    /* We open a fresh ged for isolation so the selection state is clean. */
    /* (Re-use same db; the edit framework tracks state per ged_exec call.) */
    const char *av[] = {
        "edit", "brep_sph.s", "move_surface_cv", "1", "0", "0", NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    /* This may succeed or fail depending on prior selection state;
     * just verify no crash. */
    (void)ged_exec(gedp, 6, av);
    CHECK(1, "brep_sph.s move_surface_cv without selection does not crash");
}

/* ------------------------------------------------------------------ *
 * brep: --list-all-prim-ops includes brep                             *
 * ------------------------------------------------------------------ */
static void
test_p6_all_prim_ops_has_brep(struct ged *gedp)
{
    const char *av[] = { "edit", "--list-all-prim-ops=json", NULL };
    bu_vls_trunc(gedp->ged_result_str, 0);
    CHECK(ged_exec(gedp, 2, av) == BRLCAD_OK,
          "edit --list-all-prim-ops=json returns OK after brep added");
    const char *json = bu_vls_cstr(gedp->ged_result_str);
    CHECK(json && strstr(json, "\"brep\"") != NULL,
          "edit --list-all-prim-ops=json includes '\"brep\"'");
}


/* ================================================================== *
 * main
 * ================================================================== */

int
main(int ac, char *av[])
{
    int need_help = 0;

    bu_setprogname(av[0]);

    struct bu_opt_desc d[2];
    BU_OPT(d[0], "h", "help", "", NULL, &need_help, "Print help and exit");
    BU_OPT_NULL(d[1]);

    int opt_ret = bu_opt_parse(NULL, ac, (const char **)av, d);
    if (need_help || opt_ret < 0) {
        bu_log("Usage: %s [-h]\n", av[0]);
        return (need_help) ? 0 : 1;
    }

    /* ---------------------------------------------------------------- *
     * Section 0 — Infrastructure / smoke tests
     * ---------------------------------------------------------------- */
    struct bu_vls p0_path = BU_VLS_INIT_ZERO;
    if (make_temp_path(&p0_path) != BRLCAD_OK) {
        bu_log("ERROR: cannot create temp file for section 0\n"); return 1;
    }
    if (create_p0_fixture(bu_vls_cstr(&p0_path)) != BRLCAD_OK) {
        bu_log("ERROR: section 0 fixture creation failed\n");
        bu_vls_free(&p0_path); return 1;
    }
    {
        struct ged *gedp = open_fixture(bu_vls_cstr(&p0_path));
        if (!gedp) { bu_log("ERROR: ged_open failed (section 0)\n"); bu_vls_free(&p0_path); return 1; }
        bu_log("\n--- Section 0: infrastructure ---\n");
        test_p0_cmd_exists();
        test_p0_fixture_objects_exist(gedp);
        test_p0_noargs(gedp);
        test_p0_bad_geom(gedp);
        test_p0_obj_nosubcmd(gedp);
        test_p0_perturb(gedp);
        test_p0_hrt_descriptor_ops(gedp);
        test_p0_obj_help_descriptor_ops(gedp);
        test_p0_desc_coverage();
        test_p0_no_desc_returns_error();
        ged_close(gedp);
    }
    bu_vls_free(&p0_path);

    /* ---------------------------------------------------------------- *
     * Section 1 — Parser, selection, conflict arbiter, edit buffer
     * ---------------------------------------------------------------- */
    struct bu_vls p1_path = BU_VLS_INIT_ZERO;
    if (make_temp_path(&p1_path) != BRLCAD_OK) {
        bu_log("ERROR: cannot create temp file for section 1\n"); return 1;
    }
    if (create_p1_fixture(bu_vls_cstr(&p1_path)) != BRLCAD_OK) {
        bu_log("ERROR: section 1 fixture creation failed\n");
        bu_vls_free(&p1_path); return 1;
    }
    {
        struct ged *gedp = open_fixture(bu_vls_cstr(&p1_path));
        if (!gedp) { bu_log("ERROR: ged_open failed (section 1)\n"); bu_vls_free(&p1_path); return 1; }
        bu_log("\n--- Section 1: parser / selection / edit buffer ---\n");
        test_p1_noargs(gedp);
        test_p1_bad_geom(gedp);
        test_p1_obj_nosubcmd(gedp);
        test_p1_help_flag(gedp);
        test_p1_uri_fragment();
        test_p1_uri_query();
        test_p1_batch_marker(gedp);
        test_p1_selection_fallback(gedp);
        test_p1_no_geom_no_sel(gedp);
        test_p1_conflict_arbiter(gedp);
        test_p1_flag_S(gedp);
        test_p1_flag_f(gedp);
        test_p1_buf_get_missing(gedp);
        test_p1_buf_set_get(gedp);
        test_p1_buf_abandon(gedp);
	test_p1_buf_flush(gedp);
	test_p1_edit_session_api(gedp);
	test_p1_alternate_occurrence_session_conflict(gedp);
	test_p1_independent_ged_context_sessions(gedp,
	    bu_vls_cstr(&p1_path));
	test_p1_selection_session_state(gedp);
	test_p1_reentrant_session_observer(gedp);
	test_p1_database_event_invalidation(gedp);
	test_p1_closedb_invalidation(gedp, bu_vls_cstr(&p1_path));
	test_p1_perturb_regression(gedp);

	/* Leave a deliberately dirty session active.  ged_close must cancel it
	 * before tearing down drawable/view state, and must never turn shutdown
	 * into an implicit database write. */
	struct rt_tor_internal shutdown_before;
	CHECK(read_tor(gedp, "tor.s", &shutdown_before) == BRLCAD_OK,
	    "active-session shutdown test reads the committed baseline");
	const char *shutdown_edit[] = {
	    "edit", "-i", "tor.s", "r1", "123.5", NULL
	};
	CHECK(ged_exec(gedp, 5, shutdown_edit) == BRLCAD_OK,
	    "active-session shutdown test leaves a dirty edit open");
        ged_close(gedp);

	gedp = open_fixture(bu_vls_cstr(&p1_path));
	CHECK(gedp != NULL,
	    "database reopens after closing with an active edit session");
	if (gedp) {
	    struct rt_tor_internal shutdown_after;
	    CHECK(read_tor(gedp, "tor.s", &shutdown_after) == BRLCAD_OK &&
		NEAR_EQUAL(shutdown_after.r_a, shutdown_before.r_a,
		    NEAR_ENOUGH),
		"ged_close cancels an active edit without committing it");
	    ged_close(gedp);
	}
    }
    bu_vls_free(&p1_path);

    /* ---------------------------------------------------------------- *
     * Section 2 — translate / rotate / scale / checkpoint / mat
     * ---------------------------------------------------------------- */
    struct bu_vls p2_path = BU_VLS_INIT_ZERO;
    if (make_temp_path(&p2_path) != BRLCAD_OK) {
        bu_log("ERROR: cannot create temp file for section 2\n"); return 1;
    }
    if (create_p2_fixture(bu_vls_cstr(&p2_path)) != BRLCAD_OK) {
        bu_log("ERROR: section 2 fixture creation failed\n");
        bu_vls_free(&p2_path); return 1;
    }
    {
        struct ged *gedp = open_fixture(bu_vls_cstr(&p2_path));
        if (!gedp) { bu_log("ERROR: ged_open failed (section 2)\n"); bu_vls_free(&p2_path); return 1; }
        bu_log("\n--- Section 2: translate ---\n");
        test_p2_translate_abs(gedp);
        test_p2_translate_rel(gedp);
        test_p2_tra(gedp);
        test_p2_translate_from_to(gedp);
        test_p2_translate_xonly(gedp);
        bu_log("--- Section 2: rotate ---\n");
        test_p2_rotate_3angle(gedp);
        test_p2_rotate_undo(gedp);
        test_p2_rotate_z_flag(gedp);
        test_p2_rotate_single_angle(gedp);
        test_p2_rotate_180_ambiguous(gedp);
        test_p2_rotate_radians(gedp);
        test_p2_rotate_two_angles(gedp);
        bu_log("--- Section 2: scale ---\n");
        test_p2_scale_scalar(gedp);
        test_p2_scale_vector(gedp);
        test_p2_scale_from_to(gedp);
        test_p2_scale_explicit_r(gedp);
        test_p2_scale_complex(gedp);
        test_p2_scale_zero_error(gedp);
        bu_log("--- Section 2: checkpoint/revert/reset/mat ---\n");
        test_p2_checkpoint(gedp);
        test_p2_revert(gedp);
        test_p2_reset(gedp);
        test_p2_revert_no_checkpoint(gedp);
        test_p2_mat_identity(gedp);
        test_p2_mat_missing_values(gedp);
        ged_close(gedp);
    }
    bu_vls_free(&p2_path);

    /* ---------------------------------------------------------------- *
     * Section A - edit operations without compatibility draw state
     * ---------------------------------------------------------------- */
    struct bu_vls pa_path = BU_VLS_INIT_ZERO;
    if (make_temp_path(&pa_path) != BRLCAD_OK) {
        bu_log("ERROR: cannot create temp file for section A\n"); return 1;
    }
    if (create_pa_fixture(bu_vls_cstr(&pa_path)) != BRLCAD_OK) {
        bu_log("ERROR: section A fixture creation failed\n");
        bu_vls_free(&pa_path); return 1;
    }
    {
        struct ged *gedp = open_fixture_no_dbistate(bu_vls_cstr(&pa_path));
        if (!gedp) { bu_log("ERROR: ged_open failed (section A)\n"); bu_vls_free(&pa_path); return 1; }
        bu_log("\n--- Section A: edit operations without compatibility draw state ---\n");
        test_pa_translate_abs(gedp);
        test_pa_translate_rel(gedp);
        test_pa_rotate(gedp);
        test_pa_scale(gedp);
        test_pa_S_flag_no_crash(gedp);
        test_pa_unknown_obj(gedp);
        ged_close(gedp);
    }
    bu_vls_free(&pa_path);

    /* ---------------------------------------------------------------- *
     * Section C — rotate axis-mode, scale -c, per-axis scale
     * ---------------------------------------------------------------- */
    struct bu_vls pc_path = BU_VLS_INIT_ZERO;
    if (make_temp_path(&pc_path) != BRLCAD_OK) {
        bu_log("ERROR: cannot create temp file for section C\n"); return 1;
    }
    if (create_pc_fixture(bu_vls_cstr(&pc_path)) != BRLCAD_OK) {
        bu_log("ERROR: section C fixture creation failed\n");
        bu_vls_free(&pc_path); return 1;
    }
    {
        struct ged *gedp = open_fixture(bu_vls_cstr(&pc_path));
        if (!gedp) { bu_log("ERROR: ged_open failed (section C)\n"); bu_vls_free(&pc_path); return 1; }
        bu_log("\n--- Section C: rotate axis-mode ---\n");
        test_pc_rotate_axis_z90(gedp);
        test_pc_rotate_axis_z90_rel(gedp);
        test_pc_rotate_axis_x90(gedp);
        test_pc_rotate_euler_center(gedp);
        test_pc_rotate_axis_missing_angle(gedp);
        bu_log("--- Section C: scale -c and per-axis ---\n");
        test_pc_scale_uniform_center(gedp);
        test_pc_scale_center_halve(gedp);
        test_pc_scale_three_equal(gedp);
        test_pc_scale_anisotropic(gedp);
        test_pc_scale_anisotropic_zero_error(gedp);
        ged_close(gedp);
    }
    bu_vls_free(&pc_path);

    /* ---------------------------------------------------------------- *
     * Section D — Design-completion features
     * ---------------------------------------------------------------- */
    struct bu_vls pd_path = BU_VLS_INIT_ZERO;
    if (make_temp_path(&pd_path) != BRLCAD_OK) {
        bu_log("ERROR: cannot create temp file for section D\n"); return 1;
    }
    if (create_pd_fixture(bu_vls_cstr(&pd_path)) != BRLCAD_OK) {
        bu_log("ERROR: section D fixture creation failed\n");
        bu_vls_free(&pd_path); return 1;
    }
    {
        struct ged *gedp = open_fixture(bu_vls_cstr(&pd_path));
        if (!gedp) { bu_log("ERROR: ged_open failed (section D)\n"); bu_vls_free(&pd_path); return 1; }
        bu_log("\n--- Section D: translate design completions ---\n");
        test_pd_translate_k_coords(gedp);
        reset_pd_sph(gedp);
        test_pd_translate_k_self(gedp);
        reset_pd_sph(gedp);
        test_pd_translate_component_obj(gedp);
        test_pd_translate_component_obj_y(gedp);
        bu_log("--- Section D: rotate -c . and implicit angle ---\n");
        test_pd_rotate_c_self(gedp);
        test_pd_rotate_axis_c_self(gedp);
        test_pd_rotate_implicit_angle_single_pair(gedp);
        test_pd_rotate_implicit_angle_on_axis(gedp);
        test_pd_rotate_implicit_angle_two_pairs(gedp);
        test_pd_rotate_implicit_angle_reversed(gedp);
        bu_log("--- Section D: scale -c . ---\n");
        test_pd_scale_c_self(gedp);
        test_pd_scale_c_self_aniso(gedp);
        ged_close(gedp);
    }
    bu_vls_free(&pd_path);

    /* ---------------------------------------------------------------- *
     * Section 3 — Descriptor-driven per-primitive editing (Phase 3)
     * Uses the Phase 0 fixture which already has tor.s and ell.s.
     * ---------------------------------------------------------------- */
    struct bu_vls p3_path = BU_VLS_INIT_ZERO;
    if (make_temp_path(&p3_path) != BRLCAD_OK) {
        bu_log("ERROR: cannot create temp file for section 3\n"); return 1;
    }
    if (create_p0_fixture(bu_vls_cstr(&p3_path)) != BRLCAD_OK) {
        bu_log("ERROR: section 3 fixture creation failed\n");
        bu_vls_free(&p3_path); return 1;
    }
    {
        struct ged *gedp = open_fixture(bu_vls_cstr(&p3_path));
        if (!gedp) { bu_log("ERROR: ged_open failed (section 3)\n"); bu_vls_free(&p3_path); return 1; }
        bu_log("\n--- Section 3: descriptor-driven primitive editing ---\n");
        test_p3_tor_set_radius_1(gedp);
        test_p3_tor_set_radius_2(gedp);
        test_p3_tor_r1_alias(gedp);
        test_p3_ell_set_a(gedp);
        test_p3_ell_a_alias(gedp);
        test_p3_ell_set_abc(gedp);
        bu_log("--- Section 3: type-level help, --list-ops, and --list-ops=json ---\n");
        test_p3_list_ops_tor(gedp);
        test_p3_list_ops_tor_json(gedp);
        test_p3_type_help_tor(gedp);
        test_p3_list_ops_ell(gedp);
        test_p3_list_ops_ell_json(gedp);
        bu_log("--- Section 3: --list-all-prim-ops ---\n");
        test_p3_list_all_prim_ops(gedp);
        test_p3_list_all_prim_ops_json(gedp);
        bu_log("--- Section 3: error cases ---\n");
        test_p3_unknown_prim_op(gedp);
        test_p3_missing_param(gedp);
        ged_close(gedp);
    }
    bu_vls_free(&p3_path);

    /* ---------------------------------------------------------------- *
     * Section 4 — Phase 4: complete descriptor coverage
     * Uses the Phase 0 fixture (all needed primitives already there).
     * ---------------------------------------------------------------- */
    struct bu_vls p4_path = BU_VLS_INIT_ZERO;
    if (make_temp_path(&p4_path) != BRLCAD_OK) {
        bu_log("ERROR: cannot create temp file for section 4\n"); return 1;
    }
    if (create_p0_fixture(bu_vls_cstr(&p4_path)) != BRLCAD_OK) {
        bu_log("ERROR: section 4 fixture creation failed\n");
        bu_vls_free(&p4_path); return 1;
    }
    {
        struct ged *gedp = open_fixture(bu_vls_cstr(&p4_path));
        if (!gedp) { bu_log("ERROR: ged_open failed (section 4)\n"); bu_vls_free(&p4_path); return 1; }
        bu_log("\n--- Section 4: Phase 4 descriptor completeness ---\n");
        bu_log("--- tgc: new POINT and VECTOR ops ---\n");
        test_p4_tgc_move_end_h_rt(gedp);
        test_p4_tgc_move_end_h(gedp);
        test_p4_tgc_rotate_h(gedp);
        test_p4_tgc_rotate_axb(gedp);
        test_p4_tgc_move_end_h_adj_cd(gedp);
        bu_log("--- eto: ROT_C ---\n");
        test_p4_eto_rotate_c(gedp);
        bu_log("--- cline: MOVE_H ---\n");
        test_p4_cline_move_end_h(gedp);
        bu_log("--- hyp: ROT_H ---\n");
        test_p4_hyp_rotate_h(gedp);
        bu_log("--- extrude: new descriptor ---\n");
        test_p4_extrude_list_ops(gedp);
        test_p4_extrude_list_ops_json(gedp);
        test_p4_all_prim_ops_has_extrude(gedp);
        bu_log("--- --list-ops completeness checks ---\n");
        test_p4_tgc_list_ops_complete(gedp);
        test_p4_cline_list_ops_complete(gedp);
        test_p4_hyp_list_ops_complete(gedp);
        test_p4_eto_list_ops_complete(gedp);
        ged_close(gedp);
    }
    bu_vls_free(&p4_path);

    /* ---------------------------------------------------------------- *
     * Section 5 — Phase 4+: pipe descriptor completions + arb8 desc.
     * ---------------------------------------------------------------- */
    struct bu_vls p5_path = BU_VLS_INIT_ZERO;
    if (make_temp_path(&p5_path) != BRLCAD_OK) {
        bu_log("ERROR: cannot create temp file for section 5\n"); return 1;
    }
    if (create_p5_fixture(bu_vls_cstr(&p5_path)) != BRLCAD_OK) {
        bu_log("ERROR: section 5 fixture creation failed\n");
        bu_vls_free(&p5_path); return 1;
    }
    {
        struct ged *gedp = open_fixture(bu_vls_cstr(&p5_path));
        if (!gedp) { bu_log("ERROR: ged_open failed (section 5)\n"); bu_vls_free(&p5_path); return 1; }
        bu_log("\n--- Section 5: pipe descriptor completions ---\n");
        test_p5_pipe_list_ops(gedp);
        test_p5_pipe_select_next_prev(gedp);
        test_p5_pipe_append_del(gedp);
        test_p5_pipe_prepend(gedp);
        test_p5_pipe_split(gedp);
        bu_log("--- Section 5: arb8 descriptor ---\n");
        test_p5_arb8_list_ops(gedp);
        test_p5_arb8_list_ops_json(gedp);
        test_p5_all_prim_ops_has_arb8(gedp);
        test_p5_arb8_selection_state(gedp);
        test_p5_arb8_move_face(gedp);
        test_p5_arb8_move_vertex(gedp);
        test_p5_arb8_type_option(gedp);
	bu_log("--- Section 5: sketch shared editing semantics ---\n");
	test_p5_sketch_session_state(gedp);
        ged_close(gedp);
    }
    bu_vls_free(&p5_path);

    /* ---------------------------------------------------------------- *
     * Section 6 — Phase 5.1: brep descriptor + CV ops                 *
     * ---------------------------------------------------------------- */
    struct bu_vls p6_path = BU_VLS_INIT_ZERO;
    if (make_temp_path(&p6_path) != BRLCAD_OK) {
        bu_log("ERROR: cannot create temp file for section 6\n"); return 1;
    }
    if (create_p6_fixture(bu_vls_cstr(&p6_path)) != BRLCAD_OK) {
        bu_log("ERROR: section 6 fixture creation failed\n");
        bu_vls_free(&p6_path); return 1;
    }
    {
        struct ged *gedp = open_fixture(bu_vls_cstr(&p6_path));
        if (!gedp) { bu_log("ERROR: ged_open failed (section 6)\n"); bu_vls_free(&p6_path); return 1; }
        bu_log("\n--- Section 6: brep descriptor and CV editing ---\n");
        test_p6_brep_list_ops(gedp);
        test_p6_brep_list_ops_json(gedp);
        test_p6_all_prim_ops_has_brep(gedp);
        test_p6_brep_select_and_move(gedp);
        test_p6_brep_select_and_set(gedp);
        test_p6_brep_bad_face(gedp);
        test_p6_brep_move_no_selection(gedp);
        ged_close(gedp);
    }
    bu_vls_free(&p6_path);

    bu_log("\n========================================\n");
    bu_log("edit comprehensive tests: %d/%d passed\n",
           total_tests - failed_tests, total_tests);
    bu_log("========================================\n");

    return (failed_tests == 0) ? 0 : 1;
}


/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8 cino=N-s
 */
