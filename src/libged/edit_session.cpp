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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>

#include "ged/edit.h"
#include "rt/geom.h"
#include "rt/misc.h"
#include "rt/primitives/arb8.h"
#include "rt/rt_ecmds.h"
#include "./ged_private.h"
#include "./dbi.h"

struct ged_edit_session {
    struct ged *gedp = NULL;
    struct rt_edit *edit = NULL;
    bool owns_edit = false;
    enum ged_edit_mode mode = GED_EDIT_VIEW;
    struct db_full_path path = DB_FULL_PATH_INIT_ZERO;
    size_t path_index = 0;
    struct bn_tol tolerance = BN_TOL_INIT_TOL;
    bool bbox_keypoint = false;
    bool target_view = true;
    bool view_saved = false;
    mat_t saved_rotation = MAT_INIT_IDN;
    mat_t saved_center = MAT_INIT_IDN;
    fastf_t saved_scale = 1;
    std::string selection;
    struct saved_object {
        std::string name;
        std::vector<unsigned char> bytes;
    };
    std::vector<saved_object> baseline;
    struct rt_edit_menu_item *menu = NULL;

    ~ged_edit_session()
    {
        if (owns_edit)
            rt_edit_destroy(edit);
        db_free_full_path(&path);
    }
};

void
ged_edit_session_free(struct ged_edit_session *session)
{
    delete session;
}

static struct ged_edit_session *
edit_session(const struct ged *gedp)
{
    return (gedp && gedp->i && gedp->i->i) ? gedp->i->i->edit_session : NULL;
}

static struct ged_edit_session *
edit_session_create(struct ged *gedp)
{
    struct ged_edit_session *session = edit_session(gedp);
    if (!session) {
        session = new ged_edit_session;
        session->gedp = gedp;
        gedp->i->i->edit_session = session;
    }
    return session;
}

static int
edit_error(struct ged *gedp, const char *message)
{
    if (gedp && gedp->ged_result_str)
        bu_vls_printf(gedp->ged_result_str, "%s\n", message);
    return BRLCAD_ERROR;
}

/* Keep serialized baselines rather than directory addresses: in-memory
 * databases can change an object without changing its address or size. */
static bool
edit_snapshot(struct ged_edit_session *session, bool compare)
{
    std::vector<ged_edit_session::saved_object> objects;
    for (size_t i = 0; i < session->path.fp_len; i++) {
        const std::string name = compare ? session->baseline[i].name :
            session->path.fp_names[i]->d_namep;
        struct directory *dp = db_lookup(session->gedp->dbip, name.c_str(), LOOKUP_QUIET);
        if (!dp)
            return false;
        struct bu_external external = BU_EXTERNAL_INIT_ZERO;
        if (db_get_external(&external, dp, session->gedp->dbip) < 0)
            return false;
        ged_edit_session::saved_object object;
        object.name = name;
        object.bytes.assign(external.ext_buf, external.ext_buf + external.ext_nbytes);
        bu_free_external(&external);
        if (compare && object.bytes != session->baseline[i].bytes)
            return false;
        objects.push_back(std::move(object));
    }
    if (!compare)
        session->baseline = std::move(objects);
    return true;
}

static int
edit_resolve(struct ged_edit_session *session)
{
    struct db_full_path path = DB_FULL_PATH_INIT_ZERO;
    if (db_string_to_path(&path, session->gedp->dbip, session->selection.c_str()) < 0 ||
        path.fp_len != session->path.fp_len) {
        db_free_full_path(&path);
        return edit_error(session->gedp, "The edit selection no longer exists; reject the edit");
    }
    db_free_full_path(&session->path);
    session->path = path;
    return BRLCAD_OK;
}

static int
edit_menu_set(int UNUSED(argc), const char **UNUSED(argv), void *data, void *menu)
{
    static_cast<struct ged_edit_session *>(data)->menu =
        static_cast<struct rt_edit_menu_item *>(menu);
    return BRLCAD_OK;
}

static int
edit_arb_fixed_vertex(int UNUSED(argc), const char **UNUSED(argv), void *data, void *UNUSED(extra))
{
    struct rt_edit *e = static_cast<struct ged_edit_session *>(data)->edit;
    struct rt_arb8_edit *arb = static_cast<struct rt_arb8_edit *>(e->ipe_ptr);
    int type = rt_arb_std_type(&e->es_int, e->tol);
    const int vertices_per_face = 4;
    const short *vertices = &rt_arb_vertices[type - ARB4][arb->edit_menu * vertices_per_face];
    if (type == ARB7)
        for (int i = 0; i < vertices_per_face; i++)
            if (vertices[i] == 5)
                return 5;
    return vertices[0];
}

static int
edit_axes(int UNUSED(argc), const char **UNUSED(argv), void *data, void *extra)
{
    struct rt_edit *e = static_cast<struct ged_edit_session *>(data)->edit;
    int type = e->es_int.idb_type;
    if (EDOBJ[type].ft_e_axes_pos)
        EDOBJ[type].ft_e_axes_pos(e, &e->es_int, e->tol);
    else
        VMOVE(e->curr_e_axes_pos, e->e_keypoint);
    if (extra && *static_cast<int *>(extra))
        VMOVE(e->e_axes_pos, e->curr_e_axes_pos);
    return BRLCAD_OK;
}

static int
edit_keypoint(struct ged_edit_session *session)
{
    struct rt_edit *e = session->edit;
    if (session->bbox_keypoint) {
        struct db_full_path target = DB_FULL_PATH_INIT_ZERO;
        db_dup_full_path(&target, &session->path);
        target.fp_len = session->path_index + 1;
        char *name = db_path_to_string(&target);
        point_t minimum, maximum;
        int ret = rt_obj_bounds(NULL, session->gedp->dbip, 1,
            const_cast<const char **>(&name), 0, minimum, maximum);
        bu_free(name, "edit selection path");
        db_free_full_path(&target);
        if (ret != BRLCAD_OK)
            return edit_error(session->gedp, "Unable to find edit selection bounds");
        VADD2SCALE(e->e_keypoint, minimum, maximum, 0.5);
        e->e_keytag = "bounding-box center";
    }
    VMOVE(e->e_axes_pos, e->e_keypoint);
    VMOVE(e->curr_e_axes_pos, e->e_keypoint);
    return BRLCAD_OK;
}

enum ged_edit_mode
ged_edit_session_mode(const struct ged *gedp)
{
    const struct ged_edit_session *session = edit_session(gedp);
    return session ? session->mode : GED_EDIT_VIEW;
}

struct rt_edit *
ged_edit_session_edit(struct ged *gedp)
{
    struct ged_edit_session *session = edit_session(gedp);
    if (!session || session->mode == GED_EDIT_VIEW)
        return NULL;
    if (session->owns_edit && gedp->ged_gvp) {
        session->edit->vp = gedp->ged_gvp;
        bn_mat_mul(session->edit->model2objview,
            gedp->ged_gvp->gv_model2view, session->edit->model_changes);
    }
    return session->edit;
}

const struct db_full_path *
ged_edit_session_path(const struct ged *gedp)
{
    struct ged_edit_session *session = edit_session(gedp);
    if (!session || session->mode == GED_EDIT_VIEW || !gedp->dbip)
        return NULL;
    return edit_resolve(session) == BRLCAD_OK ? &session->path : NULL;
}

size_t
ged_edit_session_path_index(const struct ged *gedp)
{
    const struct ged_edit_session *session = edit_session(gedp);
    return session && session->mode != GED_EDIT_VIEW ? session->path_index : SIZE_MAX;
}

int
ged_edit_session_matches(const struct ged *gedp, const char *path, enum ged_edit_mode mode)
{
    if (!gedp || !gedp->dbip || !path || mode == GED_EDIT_VIEW ||
        ged_edit_session_mode(gedp) != mode)
        return 0;
    const struct db_full_path *selected = ged_edit_session_path(gedp);
    struct db_full_path requested = DB_FULL_PATH_INIT_ZERO;
    bool matches = db_string_to_path(&requested, gedp->dbip, path) == 0 && selected &&
        requested.fp_len == ged_edit_session_path_index(gedp) + 1 &&
        db_full_path_match_top(&requested, selected);
    db_free_full_path(&requested);
    return matches;
}

int
ged_edit_session_begin(struct ged *gedp, struct rt_edit *borrowed,
    const struct db_full_path *path, enum ged_edit_mode mode, size_t path_index,
    const struct bn_tol *tol, struct bview *view)
{
    if (!gedp || !gedp->dbip || !path || !path->fp_len ||
        path_index >= path->fp_len || (mode != GED_EDIT_PRIMITIVE && mode != GED_EDIT_INSTANCE))
        return edit_error(gedp, "Invalid edit selection");
    if (gedp->dbip->dbi_read_only)
        return edit_error(gedp, "The database is read-only");
    if (!view && !gedp->ged_gvp)
        return edit_error(gedp, "A view is required for interactive editing");
    if (ged_edit_session_mode(gedp) != GED_EDIT_VIEW)
        return edit_error(gedp, "An editing session is already active");
    if (gedp->i && gedp->i->i && !gedp->i->i->edit_buf.empty())
        return edit_error(gedp, "Finish or reset pending scripted edits before starting an interactive edit");
    struct ged_edit_session *session = edit_session_create(gedp);
    if (tol)
        session->tolerance = *tol;
    if (borrowed) {
        if (session->owns_edit && session->edit != borrowed)
            rt_edit_destroy(session->edit);
        session->edit = borrowed;
        session->owns_edit = false;
    } else if (!session->owns_edit) {
        session->edit = rt_edit_create(NULL, NULL, NULL, NULL);
        session->owns_edit = true;
    }
    session->edit->tol = &session->tolerance;
    if (rt_edit_reinit(session->edit, const_cast<struct db_full_path *>(path),
            gedp->dbip, &session->tolerance, view ? view : gedp->ged_gvp) != BRLCAD_OK)
        return edit_error(gedp, "Unable to import selected geometry");
    db_free_full_path(&session->path);
    db_dup_full_path(&session->path, path);
    char *selection = db_path_to_string(path);
    session->selection = selection;
    bu_free(selection, "edit selection");
    session->path_index = path_index;
    if (!edit_snapshot(session, false)) {
        ged_edit_session_end(gedp);
        return edit_error(gedp, "Unable to save edit baseline");
    }
    session->mode = mode;
    session->bbox_keypoint = false;
    session->target_view = false;
    session->edit->mv_context = 1;
    if (session->owns_edit) {
        rt_edit_map_clbk_set(session->edit->m, ECMD_MENU_SET, BU_CLBK_DURING,
            edit_menu_set, session);
        rt_edit_map_clbk_set(session->edit->m, ECMD_EAXES_POS, BU_CLBK_DURING,
            edit_axes, session);
        rt_edit_map_clbk_set(session->edit->m, ECMD_ARB_SETUP_ROTFACE, BU_CLBK_DURING,
            edit_arb_fixed_vertex, session);
    }
    int type = session->edit->es_int.idb_type;
    session->menu = EDOBJ[type].ft_menu_item ?
        EDOBJ[type].ft_menu_item(&session->tolerance) : NULL;
    rt_edit_set_edflag(session->edit, RT_EDIT_IDLE);
    return edit_keypoint(session);
}

void
ged_edit_session_end(struct ged *gedp)
{
    struct ged_edit_session *session = edit_session(gedp);
    if (!session)
        return;
    if (session->edit)
        rt_edit_reset(session->edit);
    if (!session->owns_edit)
        session->edit = NULL;
    session->mode = GED_EDIT_VIEW;
    session->target_view = true;
    session->selection.clear();
    session->baseline.clear();
    session->menu = NULL;
    db_free_full_path(&session->path);
    db_full_path_init(&session->path);
}

int
ged_edit_session_reset(struct ged *gedp)
{
    struct ged_edit_session *session = edit_session(gedp);
    if (!session || session->mode == GED_EDIT_VIEW)
        return edit_error(gedp, "No editing session is active");
    if (edit_resolve(session) != BRLCAD_OK)
        return BRLCAD_ERROR;
    struct rt_edit *e = session->edit;
    struct rt_edit_map *callbacks = rt_edit_map_create();
    rt_edit_map_copy(callbacks, e->m);
    void *host_data = e->u_ptr;
    struct bu_list *vlfree = e->vlfree;
    int context = e->mv_context;
    int update_views = e->update_views;
    int ret = rt_edit_reinit(e, &session->path, gedp->dbip,
        &session->tolerance, e->vp);
    /* Reloading geometry clears rt_edit, including toolkit callbacks. */
    rt_edit_map_copy(e->m, callbacks);
    rt_edit_map_destroy(callbacks);
    e->u_ptr = host_data;
    e->vlfree = vlfree;
    e->mv_context = context;
    e->update_views = update_views;
    if (ret != BRLCAD_OK) {
        ged_edit_session_end(gedp);
        return edit_error(gedp, "Unable to reload selected geometry");
    }
    if (!edit_snapshot(session, false))
        return edit_error(gedp, "Unable to reload edit baseline");
    int type = session->edit->es_int.idb_type;
    session->menu = EDOBJ[type].ft_menu_item ? EDOBJ[type].ft_menu_item(session->edit->tol) : NULL;
    rt_edit_set_edflag(session->edit, RT_EDIT_IDLE);
    return edit_keypoint(session);
}

/* The instance index is the occurrence among leaves with the same name,
 * matching db_full_path::fp_cinst.  Keeping this identity avoids changing
 * the first occurrence when a combination contains repeated members. */
static union tree *
edit_find_leaf(union tree *tree, const char *name, int occurrence, int &index)
{
    if (!tree)
        return NULL;
    if (tree->tr_op == OP_DB_LEAF) {
        if (BU_STR_EQUAL(tree->tr_l.tl_name, name) && index++ == occurrence)
            return tree;
        return NULL;
    }
    union tree *found = edit_find_leaf(tree->tr_b.tb_left, name, occurrence, index);
    if (found || tree->tr_op == OP_NOT || tree->tr_op == OP_GUARD || tree->tr_op == OP_XNOP)
        return found;
    return edit_find_leaf(tree->tr_b.tb_right, name, occurrence, index);
}

static int
edit_write_instance(struct ged_edit_session *session)
{
    struct ged *gedp = session->gedp;
    size_t index = session->path_index;
    struct directory *target = session->path.fp_names[index ? index - 1 : 0];
    struct rt_db_internal internal;
    RT_DB_INTERNAL_INIT(&internal);
    mat_t delta;
    MAT_COPY(delta, session->edit->model_changes);
    if (!index) {
        if (rt_db_get_internal(&internal, target, gedp->dbip, delta) < 0)
            return edit_error(gedp, "Unable to import edited object");
    } else {
        mat_t parent, inverse, intermediate;
        MAT_IDN(parent);
        if (index > 1 && !db_path_to_mat(gedp->dbip, &session->path, parent, index - 1))
            return edit_error(gedp, "Unable to resolve edit parent matrix");
        bn_mat_inv(inverse, parent);
        bn_mat_mul(intermediate, delta, parent);
        bn_mat_mul(delta, inverse, intermediate);
        if (rt_db_get_internal(&internal, target, gedp->dbip, NULL) != ID_COMBINATION) {
            rt_db_free_internal(&internal);
            return edit_error(gedp, "Edit parent is not a combination");
        }
        struct rt_comb_internal *comb = static_cast<struct rt_comb_internal *>(internal.idb_ptr);
        int occurrence = 0;
        union tree *leaf = edit_find_leaf(comb->tree,
            session->path.fp_names[index]->d_namep,
            session->path.fp_cinst[index], occurrence);
        if (!leaf) {
            rt_db_free_internal(&internal);
            return edit_error(gedp, "Selected combination member no longer exists");
        }
        if (!leaf->tr_l.tl_mat) {
            leaf->tr_l.tl_mat = static_cast<matp_t>(bu_malloc(sizeof(mat_t), "edit member matrix"));
            MAT_IDN(leaf->tr_l.tl_mat);
        }
        bn_mat_mul2(delta, leaf->tr_l.tl_mat);
    }
    int ret = rt_db_put_internal(target, gedp->dbip, &internal);
    rt_db_free_internal(&internal);
    return ret < 0 ? edit_error(gedp, "Unable to write edited object") : BRLCAD_OK;
}

int
ged_edit_session_apply(struct ged *gedp, int continue_editing)
{
    struct ged_edit_session *session = edit_session(gedp);
    if (!session || session->mode == GED_EDIT_VIEW)
        return edit_error(gedp, "No editing session is active");
    if (!gedp->dbip || gedp->dbip->dbi_read_only)
        return edit_error(gedp, "The database is unavailable or read-only");
    if (!edit_snapshot(session, true))
        return edit_error(gedp, "The selection changed outside this edit; reset or reject before saving");
    if (edit_resolve(session) != BRLCAD_OK)
        return BRLCAD_ERROR;
    int ret;
    if (session->mode == GED_EDIT_INSTANCE) {
        ret = edit_write_instance(session);
    } else {
        /* Export a copy: rt_db_put_internal consumes its input even on
         * failure, but a failed accept must retain the user's work. */
        struct rt_db_internal copy;
        RT_DB_INTERNAL_INIT(&copy);
        if (rt_generic_xform(&copy, bn_mat_identity, &session->edit->es_int,
                0, gedp->dbip) < 0)
            return edit_error(gedp, "Unable to copy working geometry");
        if (copy.idb_type == ID_BOT) {
            struct rt_bot_internal *bot = static_cast<struct rt_bot_internal *>(copy.idb_ptr);
            if (bot->mode == RT_BOT_SURFACE || bot->mode == RT_BOT_SOLID) {
                bu_free(bot->thickness, "BOT thickness");
                bu_free(bot->face_mode, "BOT face mode");
                bot->thickness = NULL;
                bot->face_mode = NULL;
            } else {
                if (!bot->thickness)
                    bot->thickness = static_cast<fastf_t *>(bu_calloc(bot->num_faces, sizeof(fastf_t), "BOT thickness"));
                if (!bot->face_mode)
                    bot->face_mode = bu_bitv_new(bot->num_faces);
            }
        }
        ret = rt_db_put_internal(DB_FULL_PATH_CUR_DIR(&session->path), gedp->dbip, &copy);
        rt_db_free_internal(&copy);
        if (ret < 0)
            return edit_error(gedp, "Unable to write working geometry");
        ret = BRLCAD_OK;
    }
    if (ret != BRLCAD_OK)
        return ret;
    if (gedp->dbi_state)
        static_cast<DbiState *>(gedp->dbi_state)->update();
    const bool redraw = session->owns_edit;
    if (continue_editing) {
        ret = ged_edit_session_reset(gedp);
    } else {
        ged_edit_session_end(gedp);
    }
    if (redraw) {
        const char *args[] = {"redraw"};
        if (ged_exec(gedp, 1, args) != BRLCAD_OK)
            bu_vls_strcat(gedp->ged_result_str, "Changes were saved, but display regeneration failed\n");
    }
    return ret;
}

int
ged_edit_session_process(struct ged *gedp)
{
    struct rt_edit *e = ged_edit_session_edit(gedp);
    if (!e)
        return edit_error(gedp, "No editing session is active");
    e->base2local = gedp->dbip->dbi_base2local;
    e->local2base = gedp->dbip->dbi_local2base;
    /* Instance keypoints belong to the original selection.  model_changes
     * moves them in the display; primitive parameter processing must not
     * replace an aggregate bounding-box keypoint with a leaf vertex. */
    int fixed = e->e_keyfixed;
    if (ged_edit_session_mode(gedp) == GED_EDIT_INSTANCE)
        e->e_keyfixed = 1;
    int ret = rt_edit_process(e);
    e->e_keyfixed = fixed;
    if (bu_vls_strlen(e->log_str)) {
        bu_vls_vlscat(gedp->ged_result_str, e->log_str);
        bu_vls_trunc(e->log_str, 0);
    }
    return ret;
}

int
ged_edit_session_mouse(struct ged *gedp, const vect_t position)
{
    struct rt_edit *e = ged_edit_session_edit(gedp);
    if (!position || !std::isfinite(position[X]) || !std::isfinite(position[Y]) || !std::isfinite(position[Z]))
        return edit_error(gedp, "Expected finite mouse coordinates");
    if (!e || e->edit_flag <= RT_EDIT_IDLE)
        return edit_error(gedp, "Select an editing operation first");
    if (!EDOBJ[e->es_int.idb_type].ft_edit_xy)
        return edit_error(gedp, "The selected primitive has no mouse editing operation");
    int ret = EDOBJ[e->es_int.idb_type].ft_edit_xy(e, position);
    if (ret != BRLCAD_OK)
        return ret;
    if (ged_edit_session_mode(gedp) == GED_EDIT_INSTANCE) {
        bn_mat_mul(e->model2objview, e->vp->gv_model2view, e->model_changes);
        return BRLCAD_OK;
    }
    return e->edit_flag == RT_PARAMS_EDIT_ROT ? BRLCAD_OK : ged_edit_session_process(gedp);
}


int
ged_edit_session_plot(struct ged *gedp, const struct db_full_path *path,
    struct bview *view, struct bu_list *vlist)
{
    struct ged_edit_session *session = edit_session(gedp);
    if (!session || session->mode == GED_EDIT_VIEW || !path || !path->fp_len)
        return 0;
    if (!vlist || !gedp->dbip)
        return -1;
    if (session->mode == GED_EDIT_PRIMITIVE) {
        if (!BU_STR_EQUAL(DB_FULL_PATH_CUR_DIR(path)->d_namep, session->baseline.back().name.c_str()))
            return 0;
    } else {
        if (path->fp_len <= session->path_index)
            return 0;
        for (size_t i = 0; i <= session->path_index; i++)
            if (!BU_STR_EQUAL(path->fp_names[i]->d_namep, session->baseline[i].name.c_str()) ||
                path->fp_cinst[i] != session->path.fp_cinst[i])
                return 0;
    }
    mat_t matrix;
    if (!db_path_to_mat(gedp->dbip, const_cast<struct db_full_path *>(path), matrix, 0))
        return -1;
    if (session->mode == GED_EDIT_INSTANCE)
        bn_mat_mul2(session->edit->model_changes, matrix);
    struct rt_db_internal internal;
    RT_DB_INTERNAL_INIT(&internal);
    int ret = session->mode == GED_EDIT_PRIMITIVE ?
        rt_generic_xform(&internal, bn_mat_identity, &session->edit->es_int, 0, gedp->dbip) :
        rt_db_get_internal(&internal, DB_FULL_PATH_CUR_DIR(path), gedp->dbip, NULL);
    if (ret < 0)
        return -1;
    const struct bg_tess_tol tessellation = BG_TESS_TOL_INIT_TOL;
    struct bu_list plotted;
    BU_LIST_INIT(&plotted);
    ret = internal.idb_meth->ft_plot ?
        internal.idb_meth->ft_plot(&plotted, &internal, &tessellation, &session->tolerance, view) : -1;
    rt_db_free_internal(&internal);
    if (ret < 0) {
        bv_vlist_cleanup(&plotted);
        return -1;
    }
    /* Matrix edits may distort a primitive beyond what its parameterization
     * can represent.  Transform the plotted vectors, just as MGED does,
     * rather than attempting to import distorted primitive parameters. */
    mat_t inverse, normals;
    bn_mat_inv(inverse, matrix);
    bn_mat_trn(normals, inverse);
    struct bv_vlist *chunk;
    for (BU_LIST_FOR(chunk, bv_vlist, &plotted)) {
        for (size_t i = 0; i < chunk->nused; i++) {
            point_t point;
            VMOVE(point, chunk->pt[i]);
            switch (chunk->cmd[i]) {
                case BV_VLIST_POLY_START: case BV_VLIST_POLY_VERTNORM:
                case BV_VLIST_TRI_START: case BV_VLIST_TRI_VERTNORM:
                    MAT4X3VEC(chunk->pt[i], normals, point);
                    VUNITIZE(chunk->pt[i]);
                    break;
                case BV_VLIST_POINT_SIZE: case BV_VLIST_LINE_WIDTH:
                case BV_VLIST_MODEL_MAT:
                    break;
                default:
                    MAT4X3PNT(chunk->pt[i], matrix, point);
                    break;
            }
        }
    }
    BU_LIST_APPEND_LIST(vlist, &plotted);
    return 1;
}


static bool
edit_number(const char *text, fastf_t &value)
{
    if (!text || !*text)
        return false;
    char *end = NULL;
    value = std::strtod(text, &end);
    return end != text && !*end && std::isfinite(value);
}

static int
edit_start_command(struct ged *gedp, int argc, const char *argv[])
{
    if (argc != 2)
        return edit_error(gedp, "Usage: sed path");
    struct db_full_path path = DB_FULL_PATH_INIT_ZERO;
    if (db_string_to_path(&path, gedp->dbip, argv[1]) < 0 || !path.fp_len) {
        db_free_full_path(&path);
        return edit_error(gedp, "Unable to resolve edit path");
    }
    struct directory *dp = DB_FULL_PATH_CUR_DIR(&path);
    int ret = (dp->d_flags & RT_DIR_COMB) ?
        edit_error(gedp, "Use oed to edit a combination") :
        ged_edit_session_begin(gedp, NULL, &path, GED_EDIT_PRIMITIVE,
            path.fp_len - 1, NULL, gedp->ged_gvp);
    db_free_full_path(&path);
    return ret;
}


static int
edit_start_instance(struct ged *gedp, int argc, const char *argv[])
{
    if (argc != 2 && argc != 3)
        return edit_error(gedp, "Usage: oed path | oed lhs rhs");
    struct db_full_path lhs = DB_FULL_PATH_INIT_ZERO;
    struct db_full_path path = DB_FULL_PATH_INIT_ZERO;
    if (db_string_to_path(&lhs, gedp->dbip, argv[1]) < 0 || (argc == 2 && !lhs.fp_len)) {
        db_free_full_path(&lhs);
        return edit_error(gedp, "Unable to resolve edit path");
    }
    size_t index = argc == 2 ? lhs.fp_len - 1 : lhs.fp_len;
    std::string name = argv[1];
    if (argc == 3) {
        name += "/";
        name += argv[2];
    }
    db_free_full_path(&lhs);
    if (db_string_to_path(&path, gedp->dbip, name.c_str()) < 0 || !path.fp_len) {
        db_free_full_path(&path);
        return edit_error(gedp, "Unable to resolve edit instance");
    }
    if (DB_FULL_PATH_CUR_DIR(&path)->d_flags & RT_DIR_COMB) {
        struct display_list *entry;
        bool found = false;
        for (BU_LIST_FOR(entry, display_list, (struct bu_list *)ged_dl(gedp))) {
            struct bv_scene_obj *object;
            for (BU_LIST_FOR(object, bv_scene_obj, &entry->dl_head_scene_obj)) {
                if (!object->s_u_data)
                    continue;
                struct ged_bv_data *data = static_cast<struct ged_bv_data *>(object->s_u_data);
                if (db_full_path_match_top(&path, &data->s_fullpath)) {
                    db_free_full_path(&path);
                    db_dup_full_path(&path, &data->s_fullpath);
                    found = true;
                    break;
                }
            }
            if (found)
                break;
        }
        if (!found) {
            db_free_full_path(&path);
            return edit_error(gedp, "Draw the selected combination before editing it");
        }
    }
    int ret = ged_edit_session_begin(gedp, NULL, &path, GED_EDIT_INSTANCE,
        index, NULL, gedp->ged_gvp);
    db_free_full_path(&path);
    if (ret == BRLCAD_OK && argc == 2) {
        struct ged_edit_session *session = edit_session(gedp);
        session->bbox_keypoint = true;
        ret = edit_keypoint(session);
    }
    return ret;
}

static int
edit_transform(struct ged *gedp, int argc, const char *argv[])
{
    struct rt_edit *e = ged_edit_session_edit(gedp);
    const bool matrix = ged_edit_session_mode(gedp) == GED_EDIT_INSTANCE;
    const char *name = argv[0];
    if (BU_STR_EQUAL(name, "rot") || BU_STR_EQUAL(name, "arot")) {
        mat_t rotation;
        char coord = e->vp->gv_coord;
        int ret = BU_STR_EQUAL(name, "arot") ? ged_arot_args(gedp, argc, argv, rotation) :
            ged_rot_args(gedp, argc, argv, &coord, rotation);
        if (ret != BRLCAD_OK)
            return ret;
        rt_knob_edit_rot(e, coord, e->vp->gv_rotate_about, matrix, rotation);
        return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(name, "tra")) {
        vect_t translation;
        char coord;
        int ret = ged_tra_args(gedp, argc, argv, &coord, translation);
        if (ret != BRLCAD_OK)
            return ret;
        rt_knob_edit_tran(e, coord, matrix, translation);
        return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(name, "sca") || BU_STR_EQUAL(name, "scale") || BU_STR_EQUAL(name, "oscale")) {
        const bool absolute = BU_STR_EQUAL(name, "oscale");
        if (argc != 2 && (argc != 4 || !matrix || absolute))
            return edit_error(gedp, "Expected one positive scale factor, or three for instance editing");
        vect_t factors;
        for (int i = 0; i < argc - 1; i++)
            if (!edit_number(argv[i + 1], factors[i]) || factors[i] <= SMALL_FASTF)
                return edit_error(gedp, "Expected positive finite scale factors");
        fastf_t scale = factors[0];
        if (BU_STR_EQUAL(name, "oscale") && !matrix)
            return edit_error(gedp, "oscale requires instance editing");
        if (!matrix) {
            int saved_flag = e->edit_flag;
            int saved_mode = e->edit_mode;
            if (e->edit_mode != RT_PARAMS_EDIT_SCALE)
                rt_edit_set_edflag(e, RT_PARAMS_EDIT_SCALE);
            e->es_scale = scale;
            e->acc_sc_sol *= scale;
            int ret = ged_edit_session_process(gedp);
            e->edit_flag = saved_flag;
            e->edit_mode = saved_mode;
            return ret;
        }
        mat_t change, about;
        point_t pivot;
        MAT_IDN(change);
        int axis = e->edit_flag == RT_MATRIX_EDIT_SCALE_X ? X :
            e->edit_flag == RT_MATRIX_EDIT_SCALE_Y ? Y :
            e->edit_flag == RT_MATRIX_EDIT_SCALE_Z ? Z : -1;
        if (argc == 4) {
            for (int i = 0; i < 3; i++) {
                change[i * 5] = factors[i];
                e->acc_sc[i] *= factors[i];
            }
        } else if (axis >= 0) {
            change[axis * 5] = absolute ? scale / e->acc_sc[axis] : scale;
            e->acc_sc[axis] = absolute ? scale : scale * e->acc_sc[axis];
        } else {
            change[15] = absolute ? 1.0 / (scale * e->model_changes[15]) : 1.0 / scale;
        }
        MAT4X3PNT(pivot, e->model_changes, e->e_keypoint);
        bn_mat_xform_about_pnt(about, change, pivot);
        bn_mat_mul2(about, e->model_changes);
        return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(name, "translate") || BU_STR_EQUAL(name, "orot") || BU_STR_EQUAL(name, "rotobj")) {
        bool translation = BU_STR_EQUAL(name, "translate");
        bool incremental = argc == 5 && BU_STR_EQUAL(argv[1], "-i");
        int first = incremental ? 2 : 1;
        if (argc != first + 3 || (translation && incremental) || (!translation && !matrix))
            return edit_error(gedp, "Expected three editing coordinates");
        vect_t values;
        for (int i = 0; i < 3; i++)
            if (!edit_number(argv[first + i], values[i]))
                return edit_error(gedp, "Expected finite numeric editing parameters");
        if (incremental) {
            mat_t rotation;
            bn_mat_angles(rotation, V3ARGS(values));
            rt_knob_edit_rot(e, 'o', 'k', 1, rotation);
            return BRLCAD_OK;
        }
        rt_edit_set_edflag(e, translation ?
            (matrix ? RT_MATRIX_EDIT_TRANS_MODEL_XYZ : RT_PARAMS_EDIT_TRANS) : RT_MATRIX_EDIT_ROT);
        VMOVE(e->e_para, values);
        e->e_inpara = 3;
        return ged_edit_session_process(gedp);
    }
    if (BU_STR_EQUAL(name, "qorot")) {
        fastf_t values[7];
        if (!matrix || argc != 8)
            return edit_error(gedp, "Usage: qorot x y z dx dy dz angle (instance editing)");
        for (int i = 0; i < 7; i++)
            if (!edit_number(argv[i + 1], values[i]))
                return edit_error(gedp, "Expected finite numeric rotation parameters");
        point_t pivot;
        vect_t direction;
        VSCALE(pivot, values, gedp->dbip->dbi_local2base);
        VMOVE(direction, values + 3);
        if (MAGNITUDE(direction) <= SQRT_SMALL_FASTF)
            return edit_error(gedp, "Rotation direction must be nonzero");
        VUNITIZE(direction);
        mat_t rotation;
        bn_mat_arb_rot(rotation, pivot, direction, values[6] * DEG2RAD);
        bn_mat_mul2(rotation, e->model_changes);
        return BRLCAD_OK;
    }
    return edit_error(gedp, "Unknown editing transform");
}

int
ged_edit_session_exec(struct ged *gedp, int argc, const char *argv[])
{
    if (argc && argv && ged_edit_session_mode(gedp) != GED_EDIT_VIEW &&
        !edit_session(gedp)->target_view &&
        (BU_STR_EQUAL(argv[0], "rot") || BU_STR_EQUAL(argv[0], "arot") ||
         BU_STR_EQUAL(argv[0], "tra") || BU_STR_EQUAL(argv[0], "sca") ||
         BU_STR_EQUAL(argv[0], "scale") || BU_STR_EQUAL(argv[0], "keypoint")))
        return ged_edit_session_command(gedp, argc, argv);
    int ret = ged_exec(gedp, argc, argv);
    if (ret == BRLCAD_OK && argc == 3 && BU_STR_EQUAL(argv[0], "cpi") &&
        ged_edit_session_mode(gedp) == GED_EDIT_VIEW) {
        const char *draw[] = {"draw", argv[2]};
        ret = ged_exec(gedp, 2, draw);
        if (ret == BRLCAD_OK) {
            const char *sed[] = {"sed", argv[2]};
            ret = ged_edit_session_command(gedp, 2, sed);
        }
    }
    return ret;
}

static const struct {
    const char *label;
    const char *button;
    int flag;
} instance_operations[] = {
    {"Rotate", "orot", RT_MATRIX_EDIT_ROT},
    {"Scale", "oscale", RT_MATRIX_EDIT_SCALE},
    {"X Move", "ox", RT_MATRIX_EDIT_TRANS_VIEW_X},
    {"Y Move", "oy", RT_MATRIX_EDIT_TRANS_VIEW_Y},
    {"XY Move", "oxy", RT_MATRIX_EDIT_TRANS_VIEW_XY},
    {"Scale X", "oxscale", RT_MATRIX_EDIT_SCALE_X},
    {"Scale Y", "oyscale", RT_MATRIX_EDIT_SCALE_Y},
    {"Scale Z", "ozscale", RT_MATRIX_EDIT_SCALE_Z}
};

const char *
ged_edit_session_menu_label(const struct ged *gedp, size_t row)
{
    const struct ged_edit_session *session = edit_session(gedp);
    if (!session || session->mode == GED_EDIT_VIEW)
        return NULL;
    if (session->mode == GED_EDIT_INSTANCE) {
        if (!row)
            return "Object Edit";
        return row <= sizeof(instance_operations) / sizeof(instance_operations[0]) ?
            instance_operations[row - 1].label : NULL;
    }
    if (session->menu) {
        for (size_t i = 0; session->menu[i].menu_string[0]; i++)
            if (i == row)
                return session->menu[i].menu_string;
    }
    return NULL;
}

static int
edit_press(struct ged *gedp, const char *name)
{
    if (BU_STR_EQUAL(name, "accept"))
        return ged_edit_session_apply(gedp, 0);
    if (BU_STR_EQUAL(name, "reject")) {
        ged_edit_session_end(gedp);
        return BRLCAD_OK;
    }
    struct bview *view = gedp->ged_gvp;
    if (BU_STR_EQUAL(name, "adc")) {
        const char *args[] = {"adc"};
        return ged_exec(gedp, 1, args);
    }
    if (BU_STR_EQUAL(name, "zoomin") || BU_STR_EQUAL(name, "zoomout")) {
        const char *args[] = {"zoom", BU_STR_EQUAL(name, "zoomin") ? "2" : "0.5"};
        return ged_exec(gedp, 2, args);
    }
    if (BU_STR_EQUAL(name, "reset")) {
        const char *args[] = {"autoview"};
        int ret = ged_exec(gedp, 1, args);
        return ret == BRLCAD_OK ? edit_press(gedp, "top") : ret;
    }
    if (BU_STR_EQUAL(name, "save") || BU_STR_EQUAL(name, "restore")) {
        if (!view) return edit_error(gedp, "No view is active");
        struct ged_edit_session *session = edit_session_create(gedp);
        if (BU_STR_EQUAL(name, "save")) {
            MAT_COPY(session->saved_rotation, view->gv_rotation);
            MAT_COPY(session->saved_center, view->gv_center);
            session->saved_scale = view->gv_scale;
            session->view_saved = true;
        } else if (session->view_saved) {
            MAT_COPY(view->gv_rotation, session->saved_rotation);
            MAT_COPY(view->gv_center, session->saved_center);
            bv_view_set_scale(view, session->saved_scale);
            bv_update(view);
        }
        return BRLCAD_OK;
    }
    struct standard_view { const char *name; vect_t angles; };
    static const struct standard_view views[] = {
        {"top", {0, 0, 0}}, {"bottom", {180, 0, 0}},
        {"right", {270, 0, 0}}, {"left", {270, 0, 180}},
        {"front", {270, 0, 270}}, {"rear", {270, 0, 90}}
    };
    for (const auto &standard : views) {
        if (!BU_STR_EQUAL(name, standard.name))
            continue;
        if (!view)
            return edit_error(gedp, "No view is active");
        bn_mat_angles(view->gv_rotation, V3ARGS(standard.angles));
        bv_update(view);
        return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(name, "35,25") || BU_STR_EQUAL(name, "45,45")) {
        if (!view)
            return edit_error(gedp, "No view is active");
        vect_t aet = {35, 25, 0};
        if (BU_STR_EQUAL(name, "45,45")) VSET(aet, 45, 45, 0);
        bv_view_set_aet(view, aet);
        return BRLCAD_OK;
    }
    struct ged_edit_session *session = edit_session(gedp);
    if (!session || session->mode == GED_EDIT_VIEW)
        return edit_error(gedp, "No editing session is active");
    struct rt_edit *e = session->edit;
    if (session->mode == GED_EDIT_INSTANCE) {
        for (const auto &operation : instance_operations) {
            if (BU_STR_EQUAL(name, operation.label) || BU_STR_EQUAL(name, operation.button)) {
                rt_edit_set_edflag(e, operation.flag);
                return BRLCAD_OK;
            }
        }
        return edit_error(gedp, "Unknown instance editing operation");
    }
    if (BU_STR_EQUAL(name, "sedit")) {
        session->menu = EDOBJ[e->es_int.idb_type].ft_menu_item ?
            EDOBJ[e->es_int.idb_type].ft_menu_item(e->tol) : NULL;
        rt_edit_set_edflag(e, RT_EDIT_IDLE);
    } else if (BU_STR_EQUAL(name, "sxy")) {
        rt_edit_set_edflag(e, RT_PARAMS_EDIT_TRANS);
    } else if (BU_STR_EQUAL(name, "srot")) {
        rt_edit_set_edflag(e, RT_PARAMS_EDIT_ROT);
    } else if (BU_STR_EQUAL(name, "sscale")) {
        rt_edit_set_edflag(e, RT_PARAMS_EDIT_SCALE);
    } else if (session->menu) {
        for (int i = 1; session->menu[i].menu_string[0]; i++) {
            struct rt_edit_menu_item item = session->menu[i];
            if (!BU_STR_EQUAL(name, item.menu_string))
                continue;
            if (item.menu_func)
                item.menu_func(e, item.menu_arg, 0, i, session);
            return BRLCAD_OK;
        }
        return edit_error(gedp, "Unknown edit menu operation");
    } else {
        return edit_error(gedp, "Unknown editing operation");
    }
    if (session->owns_edit) {
        int both = 1;
        edit_axes(0, NULL, session, &both);
    }
    return BRLCAD_OK;
}



static int
edit_knob_command(struct ged *gedp, int argc, const char *argv[])
{
    struct bview *view = gedp->ged_gvp;
    if (!view)
        return edit_error(gedp, "No view is active");
    struct rt_edit *e = ged_edit_session_edit(gedp);
    view->gv_local2base = gedp->dbip->dbi_local2base;
    view->gv_base2local = gedp->dbip->dbi_base2local;
    bool incremental = false, force_edit = false, force_view = false, model = false;
    char origin = view->gv_rotate_about;
    int first = 1;
    for (; first < argc && argv[first][0] == '-'; first++) {
        if (BU_STR_EQUAL(argv[first], "-i")) incremental = true;
        else if (BU_STR_EQUAL(argv[first], "-e")) force_edit = true;
        else if (BU_STR_EQUAL(argv[first], "-v")) force_view = true;
        else if (BU_STR_EQUAL(argv[first], "-m")) { model = true; force_view = true; }
        else if (BU_STR_EQUAL(argv[first], "-o") && first + 1 < argc) origin = argv[++first][0];
        else return edit_error(gedp, "Unknown knob option");
    }
    if (force_edit && !e)
        return edit_error(gedp, "No editing session is active");
    if (first == argc) {
        const struct bview_knobs *knobs = e ? &e->k : &view->k;
        bu_vls_printf(gedp->ged_result_str,
            "x = %g\ny = %g\nz = %g\nX = %g\nY = %g\nZ = %g\nS = %g\n"
            "ax = %g\nay = %g\naz = %g\naX = %g\naY = %g\naZ = %g\naS = %g\n",
            V3ARGS(knobs->rot_v), V3ARGS(knobs->tra_v), knobs->sca,
            V3ARGS(knobs->rot_v_abs), V3ARGS(knobs->tra_v_abs), knobs->sca_abs);
        return BRLCAD_OK;
    }
    vect_t view_rotation = VINIT_ZERO, view_translation = VINIT_ZERO;
    vect_t edit_rotation = VINIT_ZERO, edit_translation = VINIT_ZERO;
    int view_rot = 0, view_tra = 0, edit_rot = 0, edit_tra = 0, edit_scale = 0;
    for (int i = first; i < argc; i++) {
        const char *token = argv[i];
        if (BU_STR_EQUAL(token, "zero") || BU_STR_EQUAL(token, "zap")) {
            bv_knobs_reset(&view->k, 0);
            if (e) bv_knobs_reset(&e->k, BV_KNOBS_RATE);
            const char *reset[] = {"adc", "reset"};
            int ret = ged_exec(gedp, 2, reset);
            if (ret != BRLCAD_OK) return ret;
            continue;
        }
        if (BU_STR_EQUAL(token, "calibrate")) {
            VSETALL(view->k.tra_v_abs, 0);
            continue;
        }
        fastf_t value;
        if (++i >= argc || !edit_number(argv[i], value))
            return edit_error(gedp, "Expected knob token/value pairs");
        const char *adc = BU_STR_EQUAL(token, "xadc") ? "x" : BU_STR_EQUAL(token, "yadc") ? "y" :
            BU_STR_EQUAL(token, "ang1") ? "a1" : BU_STR_EQUAL(token, "ang2") ? "a2" :
            BU_STR_EQUAL(token, "distadc") ? "odst" : NULL;
        if (adc) {
            const char *args[] = {"adc", "-i", adc, argv[i]};
            const char *absolute[] = {"adc", adc, argv[i]};
            int ret = incremental ? ged_exec(gedp, 4, args) : ged_exec(gedp, 3, absolute);
            if (ret != BRLCAD_OK) return ret;
            continue;
        }
        char operation = token[0] == 'a' ? token[1] : token[0];
        int edit_mode = operation == 'S' ? RT_PARAMS_EDIT_SCALE :
            operation >= 'X' && operation <= 'Z' ? RT_PARAMS_EDIT_TRANS : RT_PARAMS_EDIT_ROT;
        bool editing = e && !force_view && (force_edit || (!edit_session(gedp)->target_view && e->edit_mode == edit_mode));
        int ret = editing ?
            rt_edit_knob_cmd_process(e, &edit_rotation, &edit_rot, &edit_translation,
                &edit_tra, &edit_scale, view, token, value, origin, incremental, NULL) :
            bv_knobs_cmd_process(&view_rotation, &view_rot, &view_translation,
                &view_tra, view, token, value, origin, model || (!force_view && view->gv_coord == 'm'), incremental);
        if (ret != BRLCAD_OK)
            return edit_error(gedp, "Unknown knob token");
    }
    char coord = force_view ? (model ? 'm' : 'v') : view->gv_coord;
    if (view_tra) {
        char saved = view->gv_coord;
        view->gv_coord = coord;
        bv_knobs_tran(view, view_translation, coord == 'm');
        view->gv_coord = saved;
    }
    if (view_rot)
        bv_knobs_rot(view, view_rotation, origin, coord, e ? e->acc_rot_sol : NULL,
            e ? e->e_keypoint : NULL);
    bool matrix = ged_edit_session_mode(gedp) == GED_EDIT_INSTANCE;
    if (edit_rot) {
        mat_t rotation;
        bn_mat_angles(rotation, V3ARGS(edit_rotation));
        rt_knob_edit_rot(e, view->gv_coord, origin, matrix, rotation);
    }
    if (edit_tra) rt_knob_edit_tran(e, view->gv_coord, matrix, edit_translation);
    if (edit_scale) rt_knob_edit_sca(e, matrix);
    bv_update_rate_flags(view);
    if (e) bv_knobs_update_rate_flags(&e->k);
    return BRLCAD_OK;
}

int
ged_edit_session_tick(struct ged *gedp, fastf_t steps)
{
    if (!gedp || !gedp->dbip || !gedp->ged_gvp || !std::isfinite(steps) || steps <= 0)
        return 0;
    struct bview *view = gedp->ged_gvp;
    struct rt_edit *e = ged_edit_session_edit(gedp);
    const bool matrix = ged_edit_session_mode(gedp) == GED_EDIT_INSTANCE;
    const fastf_t translation_fraction = 0.05;
    const fastf_t edit_scale_fraction = 0.01;
    const fastf_t view_scale_fraction = 0.1;
    int changed = 0;
    char saved_coord = view->gv_coord;
    int saved_flag = e ? e->edit_flag : RT_EDIT_IDLE;
    int saved_mode = e ? e->edit_mode : RT_EDIT_DEFAULT;
    for (int editing = 0; editing <= (e ? 1 : 0); editing++) {
        struct bview_knobs *k = editing ? &e->k : &view->k;
        struct rate { const fastf_t *values; int active; char coord; char origin; bool rotation; };
        const struct rate rates[] = {
            {k->rot_m, k->rot_m_flag, 'm', k->origin_m, true},
            {k->rot_o, editing ? k->rot_o_flag : 0, 'o', k->origin_o, true},
            {k->rot_v, k->rot_v_flag, 'v', k->origin_v, true},
            {k->tra_m, k->tra_m_flag, 'm', 'v', false},
            {k->tra_v, k->tra_v_flag, 'v', 'v', false}
        };
        for (const auto &rate : rates) {
            if (!rate.active) continue;
            view->gv_coord = rate.coord;
            if (editing && !matrix)
                rt_edit_set_edflag(e, rate.rotation ? RT_PARAMS_EDIT_ROT : RT_PARAMS_EDIT_TRANS);
            fastf_t factor = steps * (rate.rotation ? 1 : translation_fraction * view->gv_scale * gedp->dbip->dbi_base2local);
            std::vector<std::string> args = {"knob", "-i", editing ? "-e" : rate.coord == 'm' ? "-m" : "-v", "-o", std::string(1, rate.origin)};
            for (int axis = 0; axis < 3; axis++) {
                args.push_back(std::string("a") + static_cast<char>((rate.rotation ? 'x' : 'X') + axis));
                std::ostringstream value;
                value << std::setprecision(17) << rate.values[axis] * factor;
                args.push_back(value.str());
            }
            std::vector<const char *> argv;
            for (const auto &arg : args) argv.push_back(arg.c_str());
            if (edit_knob_command(gedp, argv.size(), argv.data()) != BRLCAD_OK) {
                changed = -1;
                break;
            }
            changed = 1;
        }
        if (changed < 0) break;
        if (!k->sca_flag) continue;
        if (editing) {
            if (!matrix) rt_edit_set_edflag(e, RT_PARAMS_EDIT_SCALE);
            fastf_t scale = 1 + k->sca * edit_scale_fraction * steps;
            std::ostringstream value;
            value << std::setprecision(17) << scale;
            std::string text = value.str();
            const char *argv[] = {"sca", text.c_str()};
            if (edit_transform(gedp, 2, argv) != BRLCAD_OK) { changed = -1; break; }
        } else {
            fastf_t scale = std::pow(1 - k->sca * view_scale_fraction, steps);
            if (!std::isfinite(scale) || scale <= SMALL_FASTF) { changed = -1; break; }
            bv_view_set_scale(view, view->gv_scale * scale);
            bv_update(view);
        }
        changed = 1;
    }
    view->gv_coord = saved_coord;
    if (e) { e->edit_flag = saved_flag; e->edit_mode = saved_mode; }
    if (changed < 0) {
        bv_knobs_reset(&view->k, BV_KNOBS_RATE);
        if (e) bv_knobs_reset(&e->k, BV_KNOBS_RATE);
        edit_error(gedp, "Knob rates stopped after an editing error");
    }
    return changed;
}

static int
edit_arb_command(struct ged *gedp, int argc, const char *argv[])
{
    struct rt_edit *e = ged_edit_session_edit(gedp);
    if (!e || ged_edit_session_mode(gedp) != GED_EDIT_PRIMITIVE || e->es_int.idb_type != ID_ARB8)
        return edit_error(gedp, "This operation requires an active ARB primitive edit");
    e->local2base = gedp->dbip->dbi_local2base;
    e->base2local = gedp->dbip->dbi_base2local;
    const struct rt_arb_internal original = *static_cast<struct rt_arb_internal *>(e->es_int.idb_ptr);
    int ret = BRLCAD_ERROR;
    if (BU_STR_EQUAL(argv[0], "facedef")) {
        ret = rt_arb_facedef(e, argc, argv);
    } else if (BU_STR_EQUAL(argv[0], "edgedir")) {
        ret = rt_arb_edgedir(e, argc, argv);
    } else if (BU_STR_EQUAL(argv[0], "permute")) {
        if (argc != 2)
            return edit_error(gedp, "Usage: permute vertex_labels");
        ret = arb_permute(static_cast<struct rt_arb_internal *>(e->es_int.idb_ptr), argv[1], e->tol);
    } else {
        fastf_t face;
        if (argc != 3 || !edit_number(argv[1], face) || !EQUAL(std::floor(face), face) || face < 123 || face > 8765)
            return edit_error(gedp, "Expected a face label and extrusion distance or mirror axis");
        struct rt_arb_internal *arb = static_cast<struct rt_arb_internal *>(e->es_int.idb_ptr);
        fastf_t planes[7][4];
        int type = rt_arb_std_type(&e->es_int, e->tol);
        if (rt_arb_calc_planes(e->log_str, arb, type, planes, e->tol))
            return edit_error(gedp, "Unable to calculate ARB planes");
        if (BU_STR_EQUAL(argv[0], "extrude")) {
            fastf_t distance;
            if (!edit_number(argv[2], distance))
                return edit_error(gedp, "Expected a finite extrusion distance");
            ret = arb_extrude(arb, static_cast<int>(face),
                distance * e->e_mat[15] * e->local2base, e->tol, planes);
        } else {
            ret = arb_mirror_face_axis(arb, planes, static_cast<int>(face), argv[2], e->tol);
        }
    }
    if (bu_vls_strlen(e->log_str)) {
        bu_vls_vlscat(gedp->ged_result_str, e->log_str);
        bu_vls_trunc(e->log_str, 0);
    }
    if (ret != BRLCAD_OK) {
        *static_cast<struct rt_arb_internal *>(e->es_int.idb_ptr) = original;
        return edit_error(gedp, "ARB edit failed");
    }
    struct ged_edit_session *session = edit_session(gedp);
    session->menu = EDOBJ[ID_ARB8].ft_menu_item(e->tol);
    if (EDOBJ[ID_ARB8].ft_prim_edit_reset)
        EDOBJ[ID_ARB8].ft_prim_edit_reset(e);
    return BRLCAD_OK;
}

int
ged_edit_session_command(struct ged *gedp, int argc, const char *argv[])
{
    if (!gedp || !gedp->dbip || argc < 1 || !argv || !argv[0])
        return edit_error(gedp, "An open database is required");
    const char *name = argv[0];
    if (BU_STR_EQUAL(name, "rset")) {
        if ((argc != 3 && argc != 4) || !BU_STR_EQUAL(argv[1], "var") || !gedp->ged_gvp)
            return edit_error(gedp, "Usage: rset var coords|transform [value]");
        struct ged_edit_session *session = edit_session_create(gedp);
        bool coords = BU_STR_EQUAL(argv[2], "coords");
        if (!coords && !BU_STR_EQUAL(argv[2], "transform"))
            return edit_error(gedp, "Unknown edit setting");
        if (argc == 4) {
            if (strlen(argv[3]) != 1 || !(coords ? strchr("mvo", argv[3][0]) : strchr("ve", argv[3][0])))
                return edit_error(gedp, "Invalid edit setting value");
            if (coords) gedp->ged_gvp->gv_coord = argv[3][0];
            else session->target_view = argv[3][0] == 'v';
        }
        bu_vls_putc(gedp->ged_result_str, coords ? gedp->ged_gvp->gv_coord : session->target_view ? 'v' : 'e');
        return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(name, "knob"))
        return edit_knob_command(gedp, argc, argv);
    if (BU_STR_EQUAL(name, "oed"))
        return edit_start_instance(gedp, argc, argv);
    if (BU_STR_EQUAL(name, "sed"))
        return edit_start_command(gedp, argc, argv);
    if (BU_STR_EQUAL(name, "accept"))
        return argc == 1 ? ged_edit_session_apply(gedp, 0) : edit_error(gedp, "Usage: accept");
    if (BU_STR_EQUAL(name, "reject")) {
        if (argc != 1)
            return edit_error(gedp, "Usage: reject");
        ged_edit_session_end(gedp);
        return BRLCAD_OK;
    }
    if ((BU_STR_EQUAL(name, "sed_apply") || BU_STR_EQUAL(name, "oed_apply")))
        return argc == 1 ? ged_edit_session_apply(gedp, 1) : edit_error(gedp, "Usage: sed_apply");
    if ((BU_STR_EQUAL(name, "sed_reset") || BU_STR_EQUAL(name, "oed_reset")))
        return argc == 1 ? ged_edit_session_reset(gedp) : edit_error(gedp, "Usage: sed_reset");
    if (BU_STR_EQUAL(name, "press")) {
        if (argc < 2)
            return edit_error(gedp, "Usage: press operation ...");
        for (int i = 1; i < argc; i++) {
            int ret = edit_press(gedp, argv[i]);
            if (ret != BRLCAD_OK)
                return ret;
        }
        return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(name, "status")) {
        if (argc > 2)
            return edit_error(gedp, "Usage: status [state|Viewscale|base2local|local2base|matrix]");
        const char *field = argc == 2 ? argv[1] : "state";
        if (BU_STR_EQUAL(field, "state")) {
            enum ged_edit_mode mode = ged_edit_session_mode(gedp);
            bu_vls_strcat(gedp->ged_result_str, mode == GED_EDIT_PRIMITIVE ? "SOL EDIT" :
                mode == GED_EDIT_INSTANCE ? "OBJ EDIT" : "VIEWING");
            return BRLCAD_OK;
        }
        if (BU_STR_EQUAL(field, "base2local") || BU_STR_EQUAL(field, "local2base")) {
            bu_vls_printf(gedp->ged_result_str, "%.17g", BU_STR_EQUAL(field, "base2local") ?
                gedp->dbip->dbi_base2local : gedp->dbip->dbi_local2base);
            return BRLCAD_OK;
        }
        struct bview *view = gedp->ged_gvp;
        if (!view)
            return edit_error(gedp, "No view is active");
        if (BU_STR_EQUAL(field, "Viewscale")) {
            bu_vls_printf(gedp->ged_result_str, "%.17g", view->gv_scale * gedp->dbip->dbi_base2local);
            return BRLCAD_OK;
        }
        const fastf_t *matrix = NULL;
        mat_t transformed;
        if (BU_STR_EQUAL(field, "toViewcenter")) matrix = view->gv_center;
        if (BU_STR_EQUAL(field, "Viewrot")) matrix = view->gv_rotation;
        if (BU_STR_EQUAL(field, "model2view")) matrix = view->gv_model2view;
        if (BU_STR_EQUAL(field, "view2model")) matrix = view->gv_view2model;
        if (BU_STR_EQUAL(field, "model2objview") || BU_STR_EQUAL(field, "objview2model")) {
            struct rt_edit *active = ged_edit_session_edit(gedp);
            if (active)
                bn_mat_mul(transformed, view->gv_model2view, active->model_changes);
            else
                MAT_COPY(transformed, view->gv_model2view);
            if (BU_STR_EQUAL(field, "objview2model"))
                bn_mat_inv(transformed, transformed);
            matrix = transformed;
        }
        if (!matrix)
            return edit_error(gedp, "Unknown status field");
        bn_encode_mat(gedp->ged_result_str, matrix, 1);
        return BRLCAD_OK;
    }
    struct rt_edit *e = ged_edit_session_edit(gedp);
    if (!e)
        return edit_error(gedp, "No editing session is active");
    if (BU_STR_EQUAL(name, "rot") || BU_STR_EQUAL(name, "arot") ||
        BU_STR_EQUAL(name, "tra") || BU_STR_EQUAL(name, "sca") ||
        BU_STR_EQUAL(name, "scale") || BU_STR_EQUAL(name, "orot") ||
        BU_STR_EQUAL(name, "oscale") || BU_STR_EQUAL(name, "translate") ||
        BU_STR_EQUAL(name, "rotobj") || BU_STR_EQUAL(name, "qorot"))
        return edit_transform(gedp, argc, argv);
    if (BU_STR_EQUAL(name, "keypoint")) {
        if (argc == 1) {
            point_t point;
            VSCALE(point, e->e_keypoint, gedp->dbip->dbi_base2local);
            bu_vls_printf(gedp->ged_result_str, "%s (%g, %g, %g)\n",
                e->e_keytag ? e->e_keytag : "", V3ARGS(point));
            return BRLCAD_OK;
        }
        if (argc == 2 && BU_STR_EQUAL(argv[1], "reset")) {
            e->e_keyfixed = 0;
            e->e_keytag = "";
            rt_get_solid_keypoint(e, &e->e_keypoint, &e->e_keytag, e->e_mat);
            return edit_keypoint(edit_session(gedp));
        }
        point_t point;
        if (argc != 4)
            return edit_error(gedp, "Usage: keypoint [x y z | reset]");
        for (int i = 0; i < 3; i++)
            if (!edit_number(argv[i + 1], point[i]))
                return edit_error(gedp, "Expected finite keypoint coordinates");
        VSCALE(e->e_keypoint, point, gedp->dbip->dbi_local2base);
        e->e_keyfixed = 1;
        e->e_keytag = "user-specified";
        return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(name, "M")) {
        fastf_t button, x, y;
        if (argc != 4 || !edit_number(argv[1], button) ||
            !edit_number(argv[2], x) || !edit_number(argv[3], y))
            return edit_error(gedp, "Usage: M button x y");
        if (ZERO(button))
            return BRLCAD_OK;
        vect_t position = {x * INV_BV, y * INV_BV, 0.0};
        return ged_edit_session_mouse(gedp, position);
    }
    if (BU_STR_EQUAL(name, "extrude") || BU_STR_EQUAL(name, "mirface") ||
        BU_STR_EQUAL(name, "permute") || BU_STR_EQUAL(name, "edgedir") || BU_STR_EQUAL(name, "facedef"))
        return edit_arb_command(gedp, argc, argv);
    if (BU_STR_EQUAL(name, "mmenu_get")) {
        if (argc != 2)
            return edit_error(gedp, "Usage: mmenu_get index");
        size_t row = ged_edit_session_mode(gedp) == GED_EDIT_INSTANCE ? 1 : 0;
        for (const char *label; (label = ged_edit_session_menu_label(gedp, row)); row++)
            bu_vls_printf(gedp->ged_result_str, "{%s} ", label);
        return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(name, "get_sed")) {
        if (argc != 1 && !(argc == 2 && BU_STR_EQUAL(argv[1], "-c")))
            return edit_error(gedp, "Usage: get_sed [-c]");
        bu_vls_printf(gedp->ged_result_str, "%s ", edit_session(gedp)->baseline.back().name.c_str());
        if (argc == 1)
            return e->es_int.idb_meth->ft_get(gedp->ged_result_str, &e->es_int, NULL);
        struct rt_db_internal copy;
        RT_DB_INTERNAL_INIT(&copy);
        if (rt_generic_xform(&copy, e->e_mat, &e->es_int, 0, gedp->dbip) < 0)
            return edit_error(gedp, "Unable to transform working geometry");
        int ret = copy.idb_meth->ft_get(gedp->ged_result_str, &copy, NULL);
        rt_db_free_internal(&copy);
        return ret;
    }
    if (BU_STR_EQUAL(name, "get_sed_menus")) {
        if (argc != 1 || !EDOBJ[e->es_int.idb_type].ft_menu_str)
            return edit_error(gedp, "No primitive edit menus are available");
        return EDOBJ[e->es_int.idb_type].ft_menu_str(gedp->ged_result_str, &e->es_int, e->tol);
    }
    if (BU_STR_EQUAL(name, "p")) {
        if (argc < 2 || argc > RT_EDIT_MAXPARA + 1 || e->edit_flag <= RT_EDIT_IDLE)
            return edit_error(gedp, "Usage: p value [value value]; select an editing operation first");
        fastf_t values[RT_EDIT_MAXPARA] = {0};
        for (int i = 1; i < argc; i++)
            if (!edit_number(argv[i], values[i - 1]))
                return edit_error(gedp, "Expected finite numeric editing parameters");
        e->e_inpara = argc - 1;
        memcpy(e->e_para, values, sizeof(values));
        return ged_edit_session_process(gedp);
    }
    return edit_error(gedp, "Unknown editing command");
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
