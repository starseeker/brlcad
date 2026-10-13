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
#include <initializer_list>
#include <vector>

#include "bu/app.h"
#include "bu/env.h"
#include "bu/log.h"
#include "ged.h"
#include "ged/commands.h"
#include "ged/edit.h"
#include "rt/geom.h"
#include "rt/rt_ecmds.h"

static void
check(bool condition, const char *message)
{
    if (!condition)
        bu_exit(EXIT_FAILURE, "edit session: %s\n", message);
}

static int
command(struct ged *gedp, std::initializer_list<const char *> args)
{
    std::vector<const char *> argv(args);
    bu_vls_trunc(gedp->ged_result_str, 0);
    return ged_edit_session_exec(gedp, argv.size(), argv.data());
}

static void
check_center(struct ged *gedp, double expected)
{
    struct rt_db_internal internal;
    RT_DB_INTERNAL_INIT(&internal);
    check(rt_db_get_internal(&internal, db_lookup(gedp->dbip, "sphere", LOOKUP_QUIET),
        gedp->dbip, NULL) == ID_ELL, "could not read saved sphere");
    struct rt_ell_internal *ell = static_cast<struct rt_ell_internal *>(internal.idb_ptr);
    check(NEAR_EQUAL(ell->v[X], expected, SMALL_FASTF), "unexpected saved sphere center");
    rt_db_free_internal(&internal);
}

static double
path_center(struct ged *gedp, const char *name)
{
    struct db_full_path path = DB_FULL_PATH_INIT_ZERO;
    check(db_string_to_path(&path, gedp->dbip, name) == 0, "unable to resolve saved instance");
    mat_t matrix;
    check(db_path_to_mat(gedp->dbip, &path, matrix, 0), "unable to resolve saved instance matrix");
    struct rt_db_internal internal;
    RT_DB_INTERNAL_INIT(&internal);
    check(rt_db_get_internal(&internal, DB_FULL_PATH_CUR_DIR(&path), gedp->dbip, matrix) == ID_ELL,
        "unable to read saved instance");
    double center = static_cast<struct rt_ell_internal *>(internal.idb_ptr)->v[X];
    rt_db_free_internal(&internal);
    db_free_full_path(&path);
    return center;
}

static int
preview(struct ged *gedp, const char *name)
{
    struct db_full_path path = DB_FULL_PATH_INIT_ZERO;
    check(db_string_to_path(&path, gedp->dbip, name) == 0, "unable to resolve preview instance");
    struct bu_list vectors;
    BU_LIST_INIT(&vectors);
    int ret = ged_edit_session_plot(gedp, &path, gedp->ged_gvp, &vectors);
    if (ret == 1) check(!BU_LIST_IS_EMPTY(&vectors), "preview has no vectors");
    bv_vlist_cleanup(&vectors);
    db_free_full_path(&path);
    return ret;
}

int
main(int, char **argv)
{
    bu_setprogname(argv[0]);
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);
    struct ged *gedp = ged_open("inmem", "0", 0);
    struct ged *other = ged_open("inmem", "0", 0);
    check(gedp && other, "unable to create database sessions");
    check(!ged_cmd_same("sed", "accept"), "distinct edit commands were registered as aliases");
    check(command(gedp, {"in", "sphere", "sph", "0", "0", "0", "1"}) == BRLCAD_OK,
        "unable to create sphere");
    check(command(gedp, {"edit", "-i", "sphere", "scale", "2"}) == BRLCAD_OK,
        "unable to create pending scripted edit");
    check(command(gedp, {"sed", "sphere"}) == BRLCAD_ERROR,
        "interactive edit ignored pending scripted geometry");
    check(command(gedp, {"edit", "sphere", "reset"}) == BRLCAD_OK,
        "unable to discard pending scripted edit");
    check(command(gedp, {"sed", "sphere"}) == BRLCAD_OK, "begin failed");
    bu_clbk_t axes_callback = NULL;
    void *axes_data = NULL;
    rt_edit_map_clbk_get(&axes_callback, &axes_data, ged_edit_session_edit(gedp)->m,
        ECMD_EAXES_POS, BU_CLBK_DURING);
    check(axes_callback != NULL, "session has no axes callback");
    check(ged_edit_session_mode(other) == GED_EDIT_VIEW, "editing leaked between sessions");
    check(command(gedp, {"sed", "sphere"}) == BRLCAD_ERROR, "nested edit was accepted");
    check(command(gedp, {"press", "sxy"}) == BRLCAD_OK, "operation selection failed");
    check(command(gedp, {"p", "10", "0", "0"}) == BRLCAD_OK, "working edit failed");
    check_center(gedp, 0);
    check(command(gedp, {"sed_reset"}) == BRLCAD_OK, "reset failed");
    struct rt_edit *edit = ged_edit_session_edit(gedp);
    bu_clbk_t reset_callback = NULL;
    void *reset_data = NULL;
    rt_edit_map_clbk_get(&reset_callback, &reset_data, edit->m, ECMD_EAXES_POS, BU_CLBK_DURING);
    check(reset_callback == axes_callback && reset_data == axes_data, "reset lost session callbacks");
    struct bview *original_view = gedp->ged_gvp;
    gedp->ged_gvp = other->ged_gvp;
    check(command(gedp, {"tra", "-m", "1", "0", "0"}) == BRLCAD_OK,
        "transform after changing views failed");
    check(edit->vp == other->ged_gvp, "working edit retained the previous active view");
    gedp->ged_gvp = original_view;
    check(command(gedp, {"sed_reset"}) == BRLCAD_OK, "view switch reset failed");
    struct rt_ell_internal *ell = static_cast<struct rt_ell_internal *>(edit->es_int.idb_ptr);
    check(ZERO(ell->v[X]), "reset did not restore working geometry");
    check(command(gedp, {"translate", "10", "0", "0"}) == BRLCAD_OK, "translation failed");
    check(command(gedp, {"reject"}) == BRLCAD_OK, "reject failed");
    check_center(gedp, 0);
    check(command(gedp, {"sed", "sphere"}) == BRLCAD_OK, "second begin failed");
    check(command(gedp, {"translate", "10", "0", "0"}) == BRLCAD_OK, "second translation failed");
    gedp->dbip->dbi_read_only = 1;
    check(command(gedp, {"accept"}) == BRLCAD_ERROR, "read-only accept succeeded");
    check(ged_edit_session_mode(gedp) == GED_EDIT_PRIMITIVE, "failed accept discarded session");
    ell = static_cast<struct rt_ell_internal *>(ged_edit_session_edit(gedp)->es_int.idb_ptr);
    check(NEAR_EQUAL(ell->v[X], 10, SMALL_FASTF), "failed accept discarded working geometry");
    gedp->dbip->dbi_read_only = 0;
    check(command(gedp, {"sed_apply"}) == BRLCAD_OK, "apply failed");
    check_center(gedp, 10);
    check(ged_edit_session_mode(gedp) == GED_EDIT_PRIMITIVE, "apply ended the session");
    check(command(gedp, {"translate", "20", "0", "0"}) == BRLCAD_OK, "third translation failed");
    check(command(gedp, {"sed_reset"}) == BRLCAD_OK, "reset after apply failed");
    check(command(gedp, {"accept"}) == BRLCAD_OK, "accept failed");
    check_center(gedp, 10);
    check(ged_edit_session_mode(gedp) == GED_EDIT_VIEW, "accept did not return to viewing");
    check(command(gedp, {"accept"}) == BRLCAD_ERROR, "idle accept succeeded");
    check(command(gedp, {"sed", "sphere"}) == BRLCAD_OK, "rate test begin failed");
    check(command(gedp, {"press", "sxy"}) == BRLCAD_OK, "rate operation failed");
    check(command(gedp, {"knob", "-e", "X", "1"}) == BRLCAD_OK, "rate setup failed");
    edit = ged_edit_session_edit(gedp);
    fastf_t previous = static_cast<struct rt_ell_internal *>(edit->es_int.idb_ptr)->v[X];
    check(ged_edit_session_tick(gedp, 1) == 1, "rate tick did not update geometry");
    ell = static_cast<struct rt_ell_internal *>(edit->es_int.idb_ptr);
    check(ell->v[X] > previous, "rate translation did not move working geometry");
    check(command(gedp, {"knob", "zero"}) == BRLCAD_OK, "rate zero failed");
    check(ged_edit_session_tick(gedp, 1) == 0, "zero left active knob rates");
    check(command(gedp, {"reject"}) == BRLCAD_OK, "rate reject failed");
    check_center(gedp, 10);
    check(command(gedp, {"sed", "sphere"}) == BRLCAD_OK, "script edit begin failed");
    check(command(gedp, {"edit", "-i", "sphere", "scale", "2"}) == BRLCAD_OK,
        "script did not edit active working geometry");
    ell = static_cast<struct rt_ell_internal *>(ged_edit_session_edit(gedp)->es_int.idb_ptr);
    check(NEAR_EQUAL(MAGNITUDE(ell->a), 2, SMALL_FASTF), "script created a separate working copy");
    check(command(gedp, {"edit", "sphere", "scale", "2"}) == BRLCAD_ERROR,
        "script overwrote a primitive with an active edit");
    check(command(gedp, {"reject"}) == BRLCAD_OK, "script edit reject failed");
    check(command(gedp, {"press", "top", "save", "front", "restore"}) == BRLCAD_OK,
        "standard view save/restore failed");
    check(NEAR_EQUAL(gedp->ged_gvp->gv_rotation[0], 1, SMALL_FASTF), "saved view orientation was not restored");
    check(command(gedp, {"in", "member", "sph", "0", "0", "0", "1"}) == BRLCAD_OK, "member creation failed");
    check(command(gedp, {"comb", "inner", "u", "member", "u", "member"}) == BRLCAD_OK, "repeated combination failed");
    check(command(gedp, {"comb", "outer", "u", "inner"}) == BRLCAD_OK, "parent combination failed");
    check(command(gedp, {"arced", "outer/inner", "matrix", "rarc", "xlate", "10", "0", "0"}) == BRLCAD_OK,
        "parent matrix setup failed");
    check(command(gedp, {"oed", "outer/inner", "member@1"}) == BRLCAD_OK, "instance begin failed");
    const char *menu_label = ged_edit_session_menu_label(gedp, 1);
    check(menu_label && command(gedp, {"press", menu_label}) == BRLCAD_OK &&
        ged_edit_session_edit(gedp)->edit_flag == RT_MATRIX_EDIT_ROT,
        "shared menu label did not select its operation");
    check(!ged_edit_session_menu_label(gedp, SIZE_MAX), "menu lookup exceeded its bounds");
    check(ged_edit_session_matches(gedp, "/outer/inner/member@1", GED_EDIT_INSTANCE),
        "selected instance was not recognized");
    check(!ged_edit_session_matches(gedp, "outer/inner/member", GED_EDIT_INSTANCE),
        "another occurrence matched the edit selection");
    check(command(gedp, {"sca", "2", "3", "4"}) == BRLCAD_OK, "nonuniform instance scaling failed");
    edit = ged_edit_session_edit(gedp);
    check(NEAR_EQUAL(edit->model_changes[0], 2, SMALL_FASTF) &&
        NEAR_EQUAL(edit->model_changes[5], 3, SMALL_FASTF) &&
        NEAR_EQUAL(edit->model_changes[10], 4, SMALL_FASTF), "instance scale axes are incorrect");
    check(command(gedp, {"oed_reset"}) == BRLCAD_OK, "scale reset failed");
    check(command(gedp, {"translate", "15", "0", "0"}) == BRLCAD_OK, "instance translation failed");
    check(NEAR_EQUAL(path_center(gedp, "outer/inner/member@1"), 10, SMALL_FASTF), "instance preview changed database");
    check(preview(gedp, "outer/inner/member@1") == 1, "selected occurrence has no preview");
    check(preview(gedp, "outer/inner/member") == 0, "preview affected a different occurrence");
    check(command(gedp, {"oed_apply"}) == BRLCAD_OK, "instance apply failed");
    check(NEAR_EQUAL(path_center(gedp, "outer/inner/member@1"), 15, SMALL_FASTF), "wrong instance matrix committed");
    check(NEAR_EQUAL(path_center(gedp, "outer/inner/member"), 10, SMALL_FASTF), "first occurrence changed");
    check(command(gedp, {"translate", "25", "0", "0"}) == BRLCAD_OK, "second instance translation failed");
    check(command(gedp, {"oed_reset"}) == BRLCAD_OK, "instance reset failed");
    check(command(gedp, {"accept"}) == BRLCAD_OK, "instance accept failed");
    check(NEAR_EQUAL(path_center(gedp, "outer/inner/member@1"), 15, SMALL_FASTF), "instance reset lost apply baseline");

    check(command(gedp, {"draw", "outer"}) == BRLCAD_OK, "combination display failed");
    check(command(gedp, {"oed", "outer/inner"}) == BRLCAD_OK, "combination edit begin failed");
    check(ged_edit_session_matches(gedp, "outer/inner", GED_EDIT_INSTANCE),
        "combination target was confused with its preview leaf");
    check(command(gedp, {"reject"}) == BRLCAD_OK, "combination reject failed");

    check(command(gedp, {"sed", "member"}) == BRLCAD_OK, "definition begin failed");
    check(command(gedp, {"translate", "1", "0", "0"}) == BRLCAD_OK, "definition translation failed");
    check(preview(gedp, "outer/inner/member") == 1 && preview(gedp, "outer/inner/member@1") == 1,
        "definition preview omitted an occurrence");
    check(command(gedp, {"adjust", "member", "V", "2 0 0"}) == BRLCAD_OK, "external edit failed");
    check(command(gedp, {"accept"}) == BRLCAD_ERROR, "stale edit overwrote database");
    check(ged_edit_session_mode(gedp) == GED_EDIT_PRIMITIVE, "stale save discarded working edit");
    check(command(gedp, {"sed_reset"}) == BRLCAD_OK, "stale edit reset failed");
    ell = static_cast<struct rt_ell_internal *>(ged_edit_session_edit(gedp)->es_int.idb_ptr);
    check(NEAR_EQUAL(ell->v[X], 2, SMALL_FASTF), "reset did not import external change");
    check(command(gedp, {"kill", "member"}) == BRLCAD_OK, "external deletion failed");
    check(command(gedp, {"accept"}) == BRLCAD_ERROR, "deleted primitive accept succeeded");
    check(command(gedp, {"reject"}) == BRLCAD_OK, "deleted primitive reject failed");

    check(command(gedp, {"sed", "sphere"}) == BRLCAD_OK, "close test begin failed");
    check(command(gedp, {"closedb"}) == BRLCAD_OK, "database close failed");
    check(ged_edit_session_mode(gedp) == GED_EDIT_VIEW, "database close retained edit state");
    ged_close(other);
    ged_close(gedp);
    return EXIT_SUCCESS;
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
