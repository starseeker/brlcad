/*                       E D I T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2019-2026 United States Government as represented by
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
/** @file edit.cpp
 *
 * Implementation of .g geometry editing logic not specific to individual
 * primitives.
 *
 * Relationship to ft_get / ft_adjust:
 *   ft_get (OBJ[].ft_get) reads a single named parameter from a solid's
 *   db_internal (e.g. "V", "H") and returns it as a string.  ft_adjust
 *   (OBJ[].ft_adjust) applies a set of key=value string pairs back to the
 *   db_internal.  These are the scripting-level "attr get/set" interface and
 *   are not directly connected to the interactive ECMD editing path implemented
 *   here.  The rt_edit / ECMD approach works at the numeric (fastf_t) level
 *   and drives immediate wireframe feedback; ft_get/ft_adjust are used by the
 *   Tcl "get" and "adjust" commands which operate on a string-serialized form.
 *   Similarly, ft_write_params / ft_read_params are used by the "tedit"
 *   (text-edit) workflow and also operate independently of the ECMD path.
 */

#include "common.h"

#include <cmath>
#include <cctype>
#include <map>
#include <cstring>
#include <set>
#include <string>

extern "C" {
#include "vmath.h"
#include "bv/view.h"
#include "bu/hash.h"
#include "bu/malloc.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "rt/edit.h"
#include "rt/func.h"
#include "rt/functab.h"
}

class RT_Edit_Map_Internal {
    public:
        // Key is ECMD_ type, populated from MGED_Internal map
        std::map<int, std::pair<bu_clbk_t, void *>> cmd_prerun_clbk;
        std::map<int, std::pair<bu_clbk_t, void *>> cmd_during_clbk;
        std::map<int, std::pair<bu_clbk_t, void *>> cmd_postrun_clbk;
        std::map<int, std::pair<bu_clbk_t, void *>> cmd_linger_clbk;

        std::map<bu_clbk_t, int> clbk_recursion_depth_cnt;
        std::map<int, int> cmd_recursion_depth_cnt;
};

struct rt_edit_map {
    RT_Edit_Map_Internal *i;
};

void
rt_edit_knobs_reset(struct rt_edit_knobs *k, int category)
{
    if (!k)
	return;

    if (!category || category == RT_EDIT_KNOBS_RATE) {
	k->rot_m_flag = 0;
	VSETALL(k->rot_m, 0.0);

	k->rot_o_flag = 0;
	VSETALL(k->rot_o, 0.0);

	k->rot_v_flag = 0;
	VSETALL(k->rot_v, 0.0);

	k->tra_m_flag = 0;
	VSETALL(k->tra_m, 0.0);

	k->tra_v_flag = 0;
	VSETALL(k->tra_v, 0.0);

	k->sca_flag = 0;
	k->sca = 0.0;
    }

    if (!category || category == RT_EDIT_KNOBS_ABS) {
	VSETALL(k->rot_m_abs, 0.0);
	VSETALL(k->rot_m_abs_last, 0.0);

	VSETALL(k->rot_o_abs, 0.0);
	VSETALL(k->rot_o_abs_last, 0.0);

	VSETALL(k->rot_v_abs, 0.0);
	VSETALL(k->rot_v_abs_last, 0.0);

	VSETALL(k->tra_m_abs, 0.0);
	VSETALL(k->tra_m_abs_last, 0.0);

	VSETALL(k->tra_v_abs, 0.0);
	VSETALL(k->tra_v_abs_last, 0.0);

	k->sca_abs = 0.0;
    }
}

unsigned long long
rt_edit_knobs_hash(struct rt_edit_knobs *k, struct bu_data_hash_state *state)
{
    if (!k)
	return 0ULL;

    int own_state = 0;
    if (!state) {
	state = bu_data_hash_create();
	if (!state)
	    return 0ULL;
	own_state = 1;
    }

    bu_data_hash_update(state, &k->rot_m, sizeof(k->rot_m));
    bu_data_hash_update(state, &k->rot_m_flag, sizeof(k->rot_m_flag));
    bu_data_hash_update(state, &k->origin_m, sizeof(k->origin_m));
    bu_data_hash_update(state, &k->rot_o, sizeof(k->rot_o));
    bu_data_hash_update(state, &k->rot_o_flag, sizeof(k->rot_o_flag));
    bu_data_hash_update(state, &k->origin_o, sizeof(k->origin_o));
    bu_data_hash_update(state, &k->rot_v, sizeof(k->rot_v));
    bu_data_hash_update(state, &k->rot_v_flag, sizeof(k->rot_v_flag));
    bu_data_hash_update(state, &k->origin_v, sizeof(k->origin_v));

    bu_data_hash_update(state, &k->sca, sizeof(k->sca));
    bu_data_hash_update(state, &k->sca_flag, sizeof(k->sca_flag));

    bu_data_hash_update(state, &k->tra_m, sizeof(k->tra_m));
    bu_data_hash_update(state, &k->tra_m_flag, sizeof(k->tra_m_flag));
    bu_data_hash_update(state, &k->tra_v, sizeof(k->tra_v));
    bu_data_hash_update(state, &k->tra_v_flag, sizeof(k->tra_v_flag));

    bu_data_hash_update(state, &k->rot_m_abs, sizeof(k->rot_m_abs));
    bu_data_hash_update(state, &k->rot_m_abs_last, sizeof(k->rot_m_abs_last));
    bu_data_hash_update(state, &k->rot_o_abs, sizeof(k->rot_o_abs));
    bu_data_hash_update(state, &k->rot_o_abs_last, sizeof(k->rot_o_abs_last));
    bu_data_hash_update(state, &k->rot_v_abs, sizeof(k->rot_v_abs));
    bu_data_hash_update(state, &k->rot_v_abs_last, sizeof(k->rot_v_abs_last));

    bu_data_hash_update(state, &k->sca_abs, sizeof(k->sca_abs));

    bu_data_hash_update(state, &k->tra_m_abs, sizeof(k->tra_m_abs));
    bu_data_hash_update(state, &k->tra_m_abs_last, sizeof(k->tra_m_abs_last));
    bu_data_hash_update(state, &k->tra_v_abs, sizeof(k->tra_v_abs));
    bu_data_hash_update(state, &k->tra_v_abs_last, sizeof(k->tra_v_abs_last));

    if (!own_state)
	return 0ULL;

    unsigned long long hv = bu_data_hash_val(state);
    bu_data_hash_destroy(state);
    return hv;
}

static void
rt_edit_knobs_init(struct rt_edit_knobs *k)
{
    if (!k)
	return;

    rt_edit_knobs_reset(k, RT_EDIT_KNOBS_ALL);
    k->origin_m = '\0';
    k->origin_o = '\0';
    k->origin_v = '\0';
    k->rot_m_udata = NULL;
    k->rot_o_udata = NULL;
    k->rot_v_udata = NULL;
    k->sca_udata = NULL;
    k->tra_m_udata = NULL;
    k->tra_v_udata = NULL;
}

void
rt_edit_view_init(struct rt_edit_view *v)
{
    if (!v)
	return;

    v->gv_scale = 1.0;
    v->gv_base2local = 1.0;
    v->gv_local2base = 1.0;
    v->gv_coord = 'v';
    v->gv_rotate_about = 'v';
    MAT_IDN(v->gv_rotation);
    MAT_IDN(v->gv_center);
    MAT_IDN(v->gv_model2view);
    MAT_IDN(v->gv_view2model);
}

void
rt_edit_view_from_context(struct rt_edit_view *ev, const void *view_ctx)
{
    if (!ev)
	return;

    rt_edit_view_init(ev);
    const struct bv *view =
	bv_context_view_const((const struct bv_context *)view_ctx);
    if (!view)
	return;

    ev->gv_scale = bv_scale_get(view);
    ev->gv_base2local = bv_base2local_get(view);
    ev->gv_local2base = bv_local2base_get(view);
    ev->gv_coord = bv_coord_get(view);
    ev->gv_rotate_about = bv_rotate_about_get(view);
    bv_rotation_get(ev->gv_rotation, view);
    bv_center_mat_get(ev->gv_center, view);
    bv_model2view_get(ev->gv_model2view, view);
    bv_view2model_get(ev->gv_view2model, view);
}

void
rt_edit_set_view(struct rt_edit *s, const struct rt_edit_view *v)
{
    if (!s)
	return;

    if (!v) {
	rt_edit_view_init(&s->view_storage);
	s->vp = NULL;
	return;
    }

    s->view_storage = *v;
    s->vp = &s->view_storage;
}

extern "C" struct rt_edit_map *
rt_edit_map_create(void)
{
    struct rt_edit_map *o = NULL;
    BU_GET(o, struct rt_edit_map);
    o->i = new RT_Edit_Map_Internal;
    return o;
}

extern "C" void
rt_edit_map_destroy(struct rt_edit_map *o)
{
    if (!o)
	return;
    delete o->i;
    BU_PUT(o, struct rt_edit_map);
}

struct rt_edit *
rt_edit_create(struct db_full_path *dfp, struct db_i *dbip, struct bn_tol *tol, const struct rt_edit_view *v)
{
    struct rt_edit *s;
    BU_GET(s, struct rt_edit);
    BU_GET(s->m, struct rt_edit_map);
    s->m->i = new RT_Edit_Map_Internal;

    s->dbip = dbip;
    RT_DB_INTERNAL_INIT(&s->es_int);

    bu_ptbl_init(&s->comb_insts, 8, "comb inst tbl");
    bu_avs_init_empty(&s->options);

    MAT_IDN(s->acc_rot_sol);
    MAT_IDN(s->e_invmat);
    MAT_IDN(s->e_mat);
    MAT_IDN(s->incr_change);
    MAT_IDN(s->model2objview);
    MAT_IDN(s->model_changes);
    VSETALL(s->curr_e_axes_pos , 0);
    VSETALL(s->e_axes_pos , 0);
    VSETALL(s->e_keypoint, 0);
    VSETALL(s->e_mparam, 0);
    memset(s->e_para, 0, sizeof(s->e_para));
    memset(s->e_str,  0, sizeof(s->e_str));

    rt_edit_knobs_init(&s->k);

    s->acc_sc_sol = 1.0;
    s->acc_sc_obj = 1.0;
    VSETALL(s->acc_sc, 1.0);
    s->base2local = 1.0;
    s->e_inpara = 0;
    s->e_nstr = 0;
    s->e_keyfixed = 0;
    s->e_keytag = NULL;
    s->e_mvalid = 0;
    s->edit_flag = 0;
    s->es_scale = 0.0;
    s->ipe_ptr = NULL;
    s->local2base = 1.0;
    s->mv_context = 0;
    s->snap.enabled = 0;
    s->snap.spacing = 1.0;
    BU_EXTERNAL_INIT(&s->es_ckpt);
    s->edit_mode = RT_EDIT_DEFAULT;
    struct bn_tol default_tol = BN_TOL_INIT_TOL;
    s->tol_storage = tol ? *tol : default_tol;
    s->tol = &s->tol_storage;
    s->u_ptr = NULL;
    s->view_update_requested = 0;
    s->vlfree = NULL;
    rt_edit_view_init(&s->view_storage);
    s->vp = NULL;
    if (v)
	rt_edit_set_view(s, v);
    s->dbip = NULL;

    BU_GET(s->log_str, struct bu_vls);
    bu_vls_init(s->log_str);

    if (!dfp || !dbip) {

	// It's preferable to have an existing rt_db_internal to set up the
	// prim specific , but if we don't have it we still need to stub in an
	// empty struct so the logic has *something* to work with.
	if (dfp && DB_FULL_PATH_CUR_DIR(dfp) && EDOBJ[DB_FULL_PATH_CUR_DIR(dfp)->d_minor_type].ft_prim_edit_create)
	    s->ipe_ptr = (*EDOBJ[DB_FULL_PATH_CUR_DIR(dfp)->d_minor_type].ft_prim_edit_create)(NULL);

	return s;
    }

    s->local2base = dbip->dbi_local2base;
    s->base2local = dbip->dbi_base2local;
    s->dbip = dbip;

    if (rt_db_get_internal(&s->es_int, DB_FULL_PATH_CUR_DIR(dfp), dbip, NULL) < 0) {
	rt_edit_destroy(s);
	return NULL;                         /* FAIL */
    }
    RT_CK_DB_INTERNAL(&s->es_int);

    // OK, we have the internal now, set up prim specific struct, if any
    if (dfp && DB_FULL_PATH_CUR_DIR(dfp) && EDOBJ[DB_FULL_PATH_CUR_DIR(dfp)->d_minor_type].ft_prim_edit_create)
	s->ipe_ptr = (*EDOBJ[DB_FULL_PATH_CUR_DIR(dfp)->d_minor_type].ft_prim_edit_create)(s);

    /* Save aggregate path matrix */
    (void)db_path_to_mat(dbip, dfp, s->e_mat, dfp->fp_len-1);

    /* get the inverse matrix */
    bn_mat_inv(s->e_invmat, s->e_mat);

    /* Establish initial keypoint */
    s->e_keytag = "";
    rt_get_solid_keypoint(s, &s->e_keypoint, &s->e_keytag, s->e_mat);

    return s;
}

struct rt_edit *
rt_edit_create_context(struct db_full_path *dfp, struct db_i *dbip,
		       struct bn_tol *tol, const void *view_ctx)
{
    struct rt_edit_view ev;
    rt_edit_view_from_context(&ev, view_ctx);
    return rt_edit_create(dfp, dbip, tol, &ev);
}

void
rt_edit_destroy(struct rt_edit *s)
{
    if (!s)
	return;

    struct rt_db_internal *ip = &s->es_int;

    if (s->ipe_ptr) {
	if (ip->idb_type > 0 && EDOBJ[ip->idb_type].magic == RT_FUNCTAB_MAGIC &&
		EDOBJ[ip->idb_type].ft_prim_edit_destroy)
	    (*EDOBJ[ip->idb_type].ft_prim_edit_destroy)(s->ipe_ptr);
	s->ipe_ptr = NULL;
    }

    bu_ptbl_free(&s->comb_insts);
    bu_avs_free(&s->options);

    rt_db_free_internal(&s->es_int);
    bu_free_external(&s->es_ckpt);
    BU_EXTERNAL_INIT(&s->es_ckpt);

    bu_vls_free(s->log_str);
    BU_PUT(s->log_str, struct bu_vls);

    delete s->m->i;
    BU_PUT(s->m, struct rt_edit_map);
    BU_PUT(s, struct rt_edit);
}


void
rt_edit_reset(struct rt_edit *s)
{
    if (!s)
	return;

    struct rt_db_internal *ip = &s->es_int;

    /* Free primitive-specific private state.  Must happen before
     * rt_db_free_internal() because the destructor needs the correct
     * idb_type to dispatch. */
    if (s->ipe_ptr) {
	if (ip->idb_type > 0 && EDOBJ[ip->idb_type].ft_prim_edit_destroy)
	    (*EDOBJ[ip->idb_type].ft_prim_edit_destroy)(s->ipe_ptr);
	s->ipe_ptr = NULL;
    }

    /* Free stored solid and checkpoint data, clear comb-instance table */
    rt_db_free_internal(&s->es_int);
    bu_free_external(&s->es_ckpt);
    BU_EXTERNAL_INIT(&s->es_ckpt);
    bu_ptbl_reset(&s->comb_insts);
    bu_avs_free(&s->options);
    bu_avs_init_empty(&s->options);

    /* Clear the edit-callback map (will be repopulated by
     * mged_edit_clbk_sync() at the next init_sedit call) */
    rt_edit_map_clear(s->m);

    /* Reset all per-edit-session fields to initial values, matching the
     * state produced by rt_edit_create(NULL, NULL, NULL, NULL). */
    MAT_IDN(s->acc_rot_sol);
    MAT_IDN(s->e_invmat);
    MAT_IDN(s->e_mat);
    MAT_IDN(s->incr_change);
    MAT_IDN(s->model_changes);
    VSETALL(s->curr_e_axes_pos, 0.0);
    VSETALL(s->e_axes_pos, 0.0);
    VSETALL(s->e_keypoint, 0.0);
    VSETALL(s->e_mparam, 0.0);
    memset(s->e_para, 0, sizeof(s->e_para));
    memset(s->e_str, 0, sizeof(s->e_str));

    rt_edit_knobs_init(&s->k);

    s->acc_sc_sol = 1.0;
    s->acc_sc_obj = 1.0;
    VSETALL(s->acc_sc, 1.0);
    s->base2local = 1.0;
    s->e_inpara = 0;
    s->e_nstr = 0;
    s->e_keyfixed = 0;
    s->e_keytag = NULL;
    s->e_mvalid = 0;
    s->edit_flag = 0;
    s->edit_mode = RT_EDIT_DEFAULT;
    s->es_scale = 0.0;
    s->local2base = 1.0;
    s->mv_context = 0;
    s->snap.enabled = 0;
    s->snap.spacing = 1.0;
    s->u_ptr = NULL;
    s->view_update_requested = 0;
    s->vlfree = NULL;
    rt_edit_set_view(s, NULL);

    bu_vls_trunc(s->log_str, 0);
}


int
rt_edit_reinit(struct rt_edit *s, struct db_full_path *dfp, struct db_i *dbip,
               struct bn_tol *tol, const struct rt_edit_view *v)
{
    if (!s)
	return BRLCAD_ERROR;

    /* Reset to a clean idle state first */
    rt_edit_reset(s);

    struct bn_tol default_tol = BN_TOL_INIT_TOL;
    s->tol_storage = tol ? *tol : default_tol;
    s->tol = &s->tol_storage;
    if (v)
	rt_edit_set_view(s, v);

    if (!dfp || !dbip) {
	/* Idle init with no solid — matches rt_edit_create(NULL, NULL, ...) behavior.
	 * ft_prim_edit_create is called with NULL (not s) because there is no solid
	 * internal yet; the callee is expected to create an empty private state. */
	if (dfp && DB_FULL_PATH_CUR_DIR(dfp) && EDOBJ[DB_FULL_PATH_CUR_DIR(dfp)->d_minor_type].ft_prim_edit_create)
	    s->ipe_ptr = (*EDOBJ[DB_FULL_PATH_CUR_DIR(dfp)->d_minor_type].ft_prim_edit_create)(NULL);
	return BRLCAD_OK;
    }

    s->local2base = dbip->dbi_local2base;
    s->base2local = dbip->dbi_base2local;
    s->dbip = dbip;

    if (rt_db_get_internal(&s->es_int, DB_FULL_PATH_CUR_DIR(dfp), dbip, NULL) < 0)
	return BRLCAD_ERROR;

    RT_CK_DB_INTERNAL(&s->es_int);

    /* Set up primitive-specific private edit state */
    if (EDOBJ[DB_FULL_PATH_CUR_DIR(dfp)->d_minor_type].ft_prim_edit_create)
	s->ipe_ptr = (*EDOBJ[DB_FULL_PATH_CUR_DIR(dfp)->d_minor_type].ft_prim_edit_create)(s);

    /* Compute aggregate path matrix and its inverse */
    (void)db_path_to_mat(dbip, dfp, s->e_mat, dfp->fp_len - 1);
    bn_mat_inv(s->e_invmat, s->e_mat);

    /* Establish initial keypoint */
    s->e_keytag = "";
    rt_get_solid_keypoint(s, &s->e_keypoint, &s->e_keytag, s->e_mat);

    return BRLCAD_OK;
}

int
rt_edit_reinit_context(struct rt_edit *s, struct db_full_path *dfp,
		       struct db_i *dbip, struct bn_tol *tol,
		       const void *view_ctx)
{
    struct rt_edit_view ev;
    rt_edit_view_from_context(&ev, view_ctx);
    return rt_edit_reinit(s, dfp, dbip, tol, &ev);
}

void
rt_edit_set_str(struct rt_edit *s, int index, const char *str)
{
    if (!s || index < 0 || index >= RT_EDIT_MAXSTR || !str)
	return;

    bu_strlcpy(s->e_str[index], str, RT_EDIT_MAXSTR_LEN);
    if (index >= s->e_nstr)
	s->e_nstr = index + 1;
}

void
rt_edit_cmd_values_init(struct rt_edit_cmd_values *result)
{
    if (!result)
	return;
    memset(result, 0, sizeof(*result));
}

int
rt_edit_cmd_values_set_value(struct rt_edit_cmd_values *result, int index,
	fastf_t value)
{
    if (!result || index < 0 || index >= RT_EDIT_MAXPARA)
	return BRLCAD_ERROR;
    result->values[index] = value;
    result->value_valid[index] = 1;
    if ((size_t)(index + 1) > result->value_count)
	result->value_count = (size_t)(index + 1);
    return BRLCAD_OK;
}

int
rt_edit_cmd_values_set_string(struct rt_edit_cmd_values *result, int index,
	const char *value)
{
    if (!result || index < 0 || index >= RT_EDIT_MAXSTR || !value)
	return BRLCAD_ERROR;
    bu_strlcpy(result->strings[index], value, RT_EDIT_MAXSTR_LEN);
    result->string_valid[index] = 1;
    if ((size_t)(index + 1) > result->string_count)
	result->string_count = (size_t)(index + 1);
    return BRLCAD_OK;
}

enum rt_edit_value_status
rt_edit_cmd_values_get(struct rt_edit *s, int command_id,
	struct rt_edit_cmd_values *result)
{
    if (!result)
	return RT_EDIT_VALUE_ERROR;
    rt_edit_cmd_values_init(result);
    if (!s || s->es_int.idb_type < 0 ||
	EDOBJ[s->es_int.idb_type].magic != RT_FUNCTAB_MAGIC ||
	!EDOBJ[s->es_int.idb_type].ft_edit_get_values)
	return RT_EDIT_VALUE_UNAVAILABLE;
    const int status = EDOBJ[s->es_int.idb_type].ft_edit_get_values(
	s, command_id, result);
    if (status != RT_EDIT_VALUE_OK)
	rt_edit_cmd_values_init(result);
    if (status < RT_EDIT_VALUE_ERROR || status > RT_EDIT_VALUE_OK)
	return RT_EDIT_VALUE_ERROR;
    return (enum rt_edit_value_status)status;
}

int
rt_edit_map_clbk_set(struct rt_edit_map *em, int ed_cmd, int mode, bu_clbk_t f, void *d)
{
    if (!em)
	return BRLCAD_OK;

    RT_Edit_Map_Internal *i = em->i;

    std::map<int, std::pair<bu_clbk_t, void *>> *mp = NULL;
    switch (mode) {
	case BU_CLBK_DURING:
	    mp = &i->cmd_during_clbk;
	    break;
	case BU_CLBK_POST:
	    mp = &i->cmd_postrun_clbk;
	    break;
	case BU_CLBK_LINGER:
	    mp = &i->cmd_linger_clbk;
	    break;
	default:
	    mp = &i->cmd_prerun_clbk;
	    break;
    }

    if (ed_cmd == ECMD_CLEAR_CLBKS) {
	mp->clear();
	return BRLCAD_OK;
    }

    if (!f) {
	mp->erase(ed_cmd);
	return BRLCAD_OK;
    }

    (*mp)[ed_cmd] = std::make_pair(f, d);

    return BRLCAD_OK;
}

int
rt_edit_map_clbk_get(bu_clbk_t *f, void **d, struct rt_edit_map *em, int ed_cmd, int mode)
{
    // Check for no-op case
    if (!f || !d || !em)
	return BRLCAD_OK;

    RT_Edit_Map_Internal *i = em->i;
    std::map<int, std::pair<bu_clbk_t, void *>> *mp;
    switch (mode) {
	case BU_CLBK_DURING:
	    mp = &i->cmd_during_clbk;
	    break;
	case BU_CLBK_POST:
	    mp = &i->cmd_postrun_clbk;
	    break;
	case BU_CLBK_LINGER:
	    mp = &i->cmd_linger_clbk;
	    break;
	default:
	    mp = &i->cmd_prerun_clbk;
	    break;
    }

    if (mp->find(ed_cmd) == mp->end())
	return BRLCAD_ERROR;

    std::pair<bu_clbk_t, void *> clbk_info = (*mp)[ed_cmd];

    (*f) = clbk_info.first;
    (*d) = clbk_info.second;

    return BRLCAD_OK;
}

int
rt_edit_map_clear(struct rt_edit_map *m)
{
    // Check for no-op case
    if (!m)
	return BRLCAD_OK;

    m->i->cmd_prerun_clbk.clear();
    m->i->cmd_during_clbk.clear();
    m->i->cmd_postrun_clbk.clear();
    m->i->cmd_linger_clbk.clear();
    return BRLCAD_OK;
}

int
rt_edit_map_copy(struct rt_edit_map *om, struct rt_edit_map *im)
{
    // Check for no-op case
    if (!om || !im)
	return BRLCAD_OK;

    const int modes[4] = {BU_CLBK_PRE, BU_CLBK_DURING, BU_CLBK_POST, BU_CLBK_LINGER};
    std::map<int, std::pair<bu_clbk_t, void *>> *ip;
    std::map<int, std::pair<bu_clbk_t, void *>> *op;
    for (int i = 0; i < 4; i++) {
	switch (modes[i]) {
	    case BU_CLBK_DURING:
		ip = &im->i->cmd_during_clbk;
		op = &om->i->cmd_during_clbk;
		break;
	    case BU_CLBK_POST:
		ip = &im->i->cmd_postrun_clbk;
		op = &om->i->cmd_postrun_clbk;
		break;
	    case BU_CLBK_LINGER:
		ip = &im->i->cmd_linger_clbk;
		op = &om->i->cmd_linger_clbk;
		break;
	    default:
		ip = &im->i->cmd_prerun_clbk;
		op = &om->i->cmd_prerun_clbk;
		break;
	}

	std::map<int, std::pair<bu_clbk_t, void *>>::iterator mp_it;
	for (mp_it = ip->begin(); mp_it != ip->end(); ++mp_it) {
	    (*op)[mp_it->first] = mp_it->second;
	}
    }

    return BRLCAD_OK;
}


/*
 * Keypoint in model space is established in "pt".
 * If "str" is set, then that point is used, else default
 * for this solid is selected and set.
 * "str" may be a constant string, in either upper or lower case,
 * or it may be something complex like "(3, 4)" for an ARS or spline
 * to select a particular vertex or control point.
 *
 * The keypoint selection is dispatched through ft_keypoint in the EDOBJ
 * table, which gives each primitive control over which vertex or reference
 * point is highlighted.  A structparse/structprint approach would be an
 * alternative but would require a separate mapping from parameter names to
 * geometric points; the current callback design is simpler.
 */
void
rt_get_solid_keypoint(struct rt_edit *s, point_t *pt, const char **strp, fastf_t *mat)
{
    bu_clbk_t f = NULL;
    void *d = NULL;

    if (!strp)
	return;

    struct rt_db_internal *ip = &s->es_int;
    RT_CK_DB_INTERNAL(ip);

    if (EDOBJ[ip->idb_type].ft_keypoint) {
	bu_vls_trunc(s->log_str, 0);
	*strp = (*EDOBJ[ip->idb_type].ft_keypoint)(pt, *strp, mat, s, s->tol);
	if (bu_vls_strlen(s->log_str)) {
	    rt_edit_map_clbk_get(&f, &d, s->m, ECMD_PRINT_STR, BU_CLBK_DURING);
	    if (f)
		(*f)(0, NULL, d, NULL);
	    bu_vls_trunc(s->log_str, 0);
	}
	return;
    }

    struct bu_vls ltmp = BU_VLS_INIT_ZERO;
    bu_vls_printf(&ltmp, "%s", bu_vls_cstr(s->log_str));
    bu_vls_sprintf(s->log_str, "get_solid_keypoint: unrecognized solid type (setting keypoint to origin)\n");
    rt_edit_map_clbk_get(&f, &d, s->m, ECMD_PRINT_STR, BU_CLBK_DURING);
    if (f)
	(*f)(0, NULL, d, NULL);
    bu_vls_sprintf(s->log_str, "%s", bu_vls_cstr(&ltmp));
    bu_vls_free(&ltmp);

    VSETALL(*pt, 0.0);
    *strp = "(origin)";
}


void
rt_edit_set_edflag(struct rt_edit *s, int edflag)
{
    if (!s)
	return;

    s->edit_flag = edflag;

    // In the case of the four generic (i.e. not primitive data specific) flag
    // settings, we can also set the edit_mode state.  For anything else,
    // it is the responsibility of the primitive specific logic to decode
    // edit_flag (and any other relevant info) into the proper edit_mode.
    // Applications like MGED may use edit_mode to adjust interface
    // behaviors, so it is important to have it properly set, but we can only
    // do so much here.
    switch (edflag) {
	case RT_PARAMS_EDIT_ROT:
	case RT_PARAMS_EDIT_TRANS:
	case RT_PARAMS_EDIT_SCALE:
	case RT_PARAMS_EDIT_PICK:
	case RT_MATRIX_EDIT_ROT:
	case RT_MATRIX_EDIT_TRANS_VIEW_XY:
	case RT_MATRIX_EDIT_TRANS_VIEW_X:
	case RT_MATRIX_EDIT_TRANS_VIEW_Y:
	case RT_MATRIX_EDIT_SCALE:
	case RT_MATRIX_EDIT_SCALE_X:
	case RT_MATRIX_EDIT_SCALE_Y:
	case RT_MATRIX_EDIT_SCALE_Z:
	    s->edit_mode = edflag;
	    break;
	default:
	    s->edit_mode = RT_EDIT_DEFAULT;
	    break;
    }
}

/* Processing of editing knob twists. */
int
rt_edit_knob_cmd_process(
	struct rt_edit *s,
	vect_t *rvec, int *do_rot, vect_t *tvec, int *do_tran, int *do_sca,
	const struct rt_edit_view *v, const char *cmd, fastf_t f,
	char origin, int incr_flag, void *u_data)
{
    if (!s || !rvec || !do_rot || !tvec || !do_tran || !do_sca || !cmd)
	return BRLCAD_ERROR;

    if (v)
	rt_edit_set_view(s, v);

    struct rt_edit_view *ev = s->vp;
    if (!ev)
	return BRLCAD_ERROR;

    char c = (cmd[1] == '\0') ? cmd[0] : cmd[1];

    int ind = -1;
    switch (c) {
	case 'x':
	case 'X':
	    ind = X;
	    break;
	case 'y':
	case 'Y':
	    ind = Y;
	    break;
	case 'z':
	case 'Z':
	    ind = Z;
	    break;
	case 'S':
	    ind = 0;
	    break;
    }

    if (ind < 0)
	return BRLCAD_ERROR;

    if (cmd[1] == '\0') {

	if (cmd[0] == 'x' || cmd[0] == 'y' || cmd[0] == 'z') {

	    fastf_t *rot;
	    char *orig;
	    void **edm;

	    switch (ev->gv_coord) {
		case 'm':
		    rot = &s->k.rot_m[ind];
		    orig = &s->k.origin_m;
		    edm = &s->k.rot_m_udata;
		    break;
		case 'o':
		    rot = &s->k.rot_o[ind];
		    orig = &s->k.origin_o;
		    edm = &s->k.rot_o_udata;
		    break;
		case 'v':
		default:
		    rot = &s->k.rot_v[ind];
		    orig = &s->k.origin_v;
		    edm = &s->k.rot_v_udata;
		    break;
	    }

	    *orig = origin;
	    *edm = u_data;

	    if (incr_flag) {
		*rot += f;
	    } else {
		*rot = f;
	    }

	    return BRLCAD_OK;
	}

	if (cmd[0] == 'X' || cmd[0] == 'Y' || cmd[0] == 'Z') {
	    fastf_t *tra;
	    void **edm;

	    switch (ev->gv_coord) {
		case 'm':
		case 'o':
		    tra = &s->k.tra_m[ind];
		    edm = &s->k.tra_m_udata;
		    break;
		case 'v':
		default:
		    tra = &s->k.tra_v[ind];
		    edm = &s->k.tra_v_udata;
		    break;
	    }

	    *edm = u_data;

	    if (incr_flag) {
		*tra += f;
	    } else {
		*tra = f;
	    }

	    return BRLCAD_OK;
	}

	if (cmd[0] == 'S') {

	    if (incr_flag) {
		s->k.sca += f;
	    } else {
		s->k.sca = f;
	    }

	    return BRLCAD_OK;
	}

    } /* switch cmd[0] */

    if (cmd[0] == 'a' && strlen(cmd) == 2) {

	if (cmd[1] == 'x' || cmd[1] == 'y' || cmd[1] == 'z') {

	    fastf_t *rot_c = NULL;
	    fastf_t *rot_lc = NULL;
	    fastf_t *rvec_c;

	    rvec_c = &(*rvec)[ind];

	    switch (ev->gv_coord) {
		case 'm':
		    rot_c = &s->k.rot_m_abs[ind];
		    rot_lc = &s->k.rot_m_abs_last[ind];
		    break;
		case 'o':
		    rot_c = &s->k.rot_o_abs[ind];
		    rot_lc = &s->k.rot_o_abs_last[ind];
		    break;
		case 'v':
		    rot_c = &s->k.rot_v_abs[ind];
		    rot_lc = &s->k.rot_v_abs_last[ind];
		    break;
	    }
	    if (!rot_c || !rot_lc)
		return BRLCAD_ERROR;

	    if (incr_flag) {
		*rot_c += f;
		*rvec_c = f;
	    } else {
		*rot_c = f;
		*rvec_c = f - *rot_lc;
	    }

	    /* wrap around */
	    fastf_t *arp;
	    fastf_t *larp;

	    switch (ev->gv_coord) {
		case 'm':
		    arp = s->k.rot_m_abs;
		    larp = s->k.rot_m_abs_last;
		    break;
		case 'o':
		    arp = s->k.rot_o_abs;
		    larp = s->k.rot_o_abs_last;
		    break;
		case 'v':
		    arp = s->k.rot_v_abs;
		    larp = s->k.rot_v_abs_last;
		    break;
		default:
		    bu_log("unknown mv_coords\n");
		    arp = larp = NULL;
		    return BRLCAD_ERROR;
	    }

	    if (arp[ind] < -180.0) {
		arp[ind] = arp[ind] + 360.0;
	    } else if (arp[ind] > 180.0) {
		arp[ind] = arp[ind] - 360.0;
	    }

	    larp[ind] = arp[ind];

	    *do_rot = 1;

	    return BRLCAD_OK;
	}

	if (cmd[1] == 'X' || cmd[1] == 'Y' || cmd[1] == 'Z') {
	    fastf_t *eamt = NULL;
	    fastf_t *leamt = NULL;
	    fastf_t *tvec_c;
	    fastf_t sf = f * ev->gv_local2base / ev->gv_scale;

	    tvec_c = &(*tvec)[ind];

	    switch (ev->gv_coord) {
		case 'm':
		case 'o':
		    eamt = &s->k.tra_m_abs[ind];
		    leamt = &s->k.tra_m_abs_last[ind];
		    break;
		case 'v':
		    eamt = &s->k.tra_v_abs[ind];
		    leamt = &s->k.tra_v_abs_last[ind];
		    break;
	    }

	    if (!eamt || !leamt)
		return BRLCAD_ERROR;

	    if (incr_flag) {
		*eamt += sf;
		*tvec_c = f;
	    } else {
		*eamt = sf;
		*tvec_c = f - *leamt * ev->gv_scale * ev->gv_base2local;
	    }
	    *leamt = *eamt;

	    *do_tran = 1;

	    return BRLCAD_OK;
	}

	if (cmd[1] == 'S') {
	    if (incr_flag) {
		s->k.sca_abs += f;
	    } else {
		s->k.sca_abs = f;
	    }

	    *do_sca = 1;

	    return BRLCAD_OK;
	}
    } /* switch (cmd[1]) */

    return BRLCAD_ERROR;
}

int
rt_edit_knob_cmd_process_context(
	struct rt_edit *s,
	vect_t *rvec, int *do_rot, vect_t *tvec, int *do_tran, int *do_sca,
	const void *view_ctx, const char *cmd, fastf_t f,
	char origin, int incr_flag, void *u_data)
{
    if (!view_ctx)
	return BRLCAD_ERROR;

    struct rt_edit_view ev;
    rt_edit_view_from_context(&ev, view_ctx);
    return rt_edit_knob_cmd_process(s, rvec, do_rot, tvec, do_tran, do_sca,
	    &ev, cmd, f, origin, incr_flag, u_data);
}

void
rt_knob_edit_rot(struct rt_edit *s,
	char coords,
	char rotate_about,
	int matrix_edit,
	mat_t newrot)
{
    mat_t temp1, temp2;

    s->view_update_requested = 1;

    switch (coords) {
	case 'm':
	    break;
	case 'o':
	    bn_mat_inv(temp1, s->acc_rot_sol);

	    /* transform into object rotations */
	    bn_mat_mul(temp2, s->acc_rot_sol, newrot);
	    bn_mat_mul(newrot, temp2, temp1);
	    break;
	case 'v':
	    bn_mat_inv(temp1, s->vp->gv_rotation);

	    /* transform into model rotations */
	    bn_mat_mul(temp2, temp1, newrot);
	    bn_mat_mul(newrot, temp2, s->vp->gv_rotation);
	    break;
    }

    if (!matrix_edit) {

	MAT_COPY(s->incr_change, newrot);
	bn_mat_mul2(s->incr_change, s->acc_rot_sol);
	s->e_inpara = 0;

	// Stash state to temporarily put things in the
	// solid edit rotate state (it may already have
	// been there but we're not counting on it)
	char save_rotate_about = s->vp->gv_rotate_about;
	s->vp->gv_rotate_about = rotate_about;

	int save_edflag = s->edit_flag;
	int save_mode = s->edit_mode;

	rt_edit_set_edflag(s, RT_PARAMS_EDIT_ROT);

	rt_edit_process(s);

	// Restore previous state
	s->vp->gv_rotate_about = save_rotate_about;
	s->edit_flag = save_edflag;
	s->edit_mode = save_mode;

    } else {

	point_t point;
	vect_t work;

	bn_mat_mul2(newrot, s->acc_rot_sol);

	/* find point for rotation to take place wrt */
	switch (rotate_about) {
	    case 'v':       /* View Center */
		VSET(work, 0.0, 0.0, 0.0);
		MAT4X3PNT(point, s->vp->gv_view2model, work);
		break;
	    case 'e':       /* Eye */
		VSET(work, 0.0, 0.0, 1.0);
		MAT4X3PNT(point, s->vp->gv_view2model, work);
		break;
	    case 'm':       /* Model Center */
		VSETALL(point, 0.0);
		break;
	    case 'k':
	    default:
		MAT4X3PNT(point, s->model_changes, s->e_keypoint);
	}

	/* Apply newrot to the s->model_changes matrix wrt "point" */
	mat_t t, out;
	bn_mat_xform_about_pnt(t, newrot, point);
	bn_mat_mul(out, t, s->model_changes);
	MAT_COPY(s->model_changes, out);

	/* Update the model2objview matrix */
	bn_mat_mul(s->model2objview, s->vp->gv_model2view, s->model_changes);
    }
}

void
rt_knob_edit_tran(struct rt_edit *s,
        char coords,
        int matrix_edit,
        const vect_t tvec)
{
    point_t p2;
    point_t delta;
    point_t vcenter;
    point_t work;

    /* compute delta */
    switch (coords) {
	case 'm':
	    VSCALE(delta, tvec, s->local2base);
	    break;
	case 'o':
	    VSCALE(p2, tvec, s->local2base);
	    MAT4X3PNT(delta, s->acc_rot_sol, p2);
	    break;
	case 'v':
	default:
	    VSCALE(p2, tvec, s->local2base / s->vp->gv_scale);
	    MAT4X3PNT(work, s->vp->gv_view2model, p2);
	    MAT_DELTAS_GET_NEG(vcenter, s->vp->gv_center);
	    VSUB2(delta, work, vcenter);

	    break;
    }

    if (!matrix_edit) {
	s->e_keyfixed = 0;

	// Get primitive keypoint
	struct rt_db_internal *ip = &s->es_int;
	fastf_t *mat = s->e_mat;
	point_t *pt = &s->e_keypoint;
	RT_CK_DB_INTERNAL(ip);
	if (EDOBJ[ip->idb_type].ft_keypoint) {
	    bu_vls_trunc(s->log_str, 0);
	    s->e_keytag = (*EDOBJ[ip->idb_type].ft_keypoint)(pt, s->e_keytag, mat, s, s->tol);
	    if (bu_vls_strlen(s->log_str)) {
		bu_log("%s", bu_vls_cstr(s->log_str));
		bu_vls_trunc(s->log_str, 0);
	    }
	} else {
	    bu_log("rt_edit_knobs_tra: unrecognized solid type (setting keypoint to origin)\n");
	    VSETALL(*pt, 0.0);
	    s->e_keytag = "(origin)";
	}

	int save_edflag = s->edit_flag;
	int save_mode = s->edit_mode;

	rt_edit_set_edflag(s, RT_PARAMS_EDIT_TRANS);

	VADD2(s->e_para, delta, s->curr_e_axes_pos);
	s->e_inpara = 3;
	rt_edit_process(s);
	s->edit_flag = save_edflag;
	s->edit_mode = save_mode;
    } else {
	mat_t xlatemat;
	MAT_IDN(xlatemat);
	MAT_DELTAS_VEC(xlatemat, delta);
	bn_mat_mul2(xlatemat, s->model_changes);

	/* Update the model2objview matrix */
	bn_mat_mul(s->model2objview, s->vp->gv_model2view, s->model_changes);
    }
}

#define MGED_SMALL_SCALE 1.0e-10
void
rt_knob_edit_sca(struct rt_edit *s, int matrix_edit)
{
   if (!matrix_edit) {

        fastf_t old_acc_sc_sol;

        old_acc_sc_sol = s->acc_sc_sol;

        if (-SMALL_FASTF < s->k.sca_abs && s->k.sca_abs < SMALL_FASTF)
            s->acc_sc_sol = 1.0;
        else if (s->k.sca_abs > 0.0)
            s->acc_sc_sol = 1.0 + s->k.sca_abs * 3.0;
        else {
            if ((s->k.sca_abs - MGED_SMALL_SCALE) < -1.0)
                s->k.sca_abs = -1.0 + MGED_SMALL_SCALE;

            s->acc_sc_sol = 1.0 + s->k.sca_abs;
        }

        s->es_scale = s->acc_sc_sol / old_acc_sc_sol;

	int save_edflag = s->edit_flag;
	int save_mode = s->edit_mode;

	rt_edit_set_edflag(s, RT_PARAMS_EDIT_SCALE);

	rt_edit_process(s);

	s->edit_flag = save_edflag;
	s->edit_mode = save_mode;

   } else {
       fastf_t scale;
       mat_t incr_mat;
       MAT_IDN(incr_mat);

       /* Note: interactive mouse-based matrix scaling (objedit_mouse SARROW)
        * computes scale differently - it uses the raw mousevec[Y] value
        * directly to get an incremental scale (scale = 1 + |mousevec[Y]|).
        * That path goes through edit_mscale_xy() which accepts a mousevec
        * directly.  This function handles the knob-based path where scale is
        * expressed via the k.sca_abs abstraction (0 = no scale, +/- 1 = max
        * scale up/down).  The two paths produce consistent model_changes
        * updates and are both correct for their respective input sources. */

       if (-SMALL_FASTF < s->k.sca_abs && s->k.sca_abs < SMALL_FASTF)
	   scale = 1;
       else if (s->k.sca_abs > 0.0)
	   scale = 1.0 + s->k.sca_abs * 3.0;
       else {
	   if ((s->k.sca_abs - MGED_SMALL_SCALE) < -1.0)
	       s->k.sca_abs = -1.0 + MGED_SMALL_SCALE;

	   scale = 1.0 + s->k.sca_abs;
       }

       /* switch depending on scaling option selected */
       switch (s->edit_flag) {

	   case RT_PARAMS_EDIT_SCALE:
	       /* global scaling */
	       incr_mat[15] = s->acc_sc_obj / scale;
	       s->acc_sc_obj = scale;
	       break;

	   case RT_MATRIX_EDIT_SCALE_X:
	       /* local scaling ... X-axis */
	       incr_mat[0] = scale / s->acc_sc[0];
	       /* accumulate the scale factor */
	       s->acc_sc[0] = scale;
	       break;

	   case RT_MATRIX_EDIT_SCALE_Y:
	       /* local scaling ... Y-axis */
	       incr_mat[5] = scale / s->acc_sc[1];
	       /* accumulate the scale factor */
	       s->acc_sc[1] = scale;
	       break;

	   case RT_MATRIX_EDIT_SCALE_Z:
	       /* local scaling ... Z-axis */
	       incr_mat[10] = scale / s->acc_sc[2];
	       /* accumulate the scale factor */
	       s->acc_sc[2] = scale;
	       break;

	   default:
	       bu_log("Incorrect edit flag for matrix scale:  %d\n", s->edit_flag);
	       return;
       }

       /* Have scaling take place with respect to keypoint, NOT the view
	* center.  model_changes is the matrix that will ultimately be used to
	* alter the geometry on disk. */
       mat_t t, out;
       vect_t pos_model;
       VMOVE(t, s->e_keypoint);
       MAT4X3PNT(pos_model, s->model_changes, t);
       bn_mat_xform_about_pnt(t, incr_mat, pos_model);
       bn_mat_mul(out, t, s->model_changes);
       MAT_COPY(s->model_changes, out);

       /* Update the model2objview matrix */
       bn_mat_mul(s->model2objview, s->vp->gv_model2view, s->model_changes);
   }
}

static int
edit_internal_restore(struct rt_edit *s, const struct bu_external *snapshot,
	int type)
{
    if (!s || !snapshot || !snapshot->ext_buf || type < 0 ||
	EDOBJ[type].magic != RT_FUNCTAB_MAGIC)
	return BRLCAD_ERROR;

    if (s->ipe_ptr && EDOBJ[type].ft_prim_edit_destroy)
	(*EDOBJ[type].ft_prim_edit_destroy)(s->ipe_ptr);
    s->ipe_ptr = NULL;
    rt_db_free_internal(&s->es_int);
    RT_DB_INTERNAL_INIT(&s->es_int);
    s->es_int.idb_minor_type = type;
    mat_t identity;
    MAT_IDN(identity);
    if (rt_obj_import(&s->es_int, snapshot, identity, s->dbip) < 0)
	return BRLCAD_ERROR;
    if (EDOBJ[type].ft_prim_edit_create)
	s->ipe_ptr = (*EDOBJ[type].ft_prim_edit_create)(s);
    s->e_keytag = "";
    rt_get_solid_keypoint(s, &s->e_keypoint, &s->e_keytag, s->e_mat);
    return BRLCAD_OK;
}

/*
 * A great deal of magic takes place here, to accomplish solid editing.
 *
 * Called from mged main loop after parameter entry or mouse events
 * to apply any accumulated edit parameters to the current solid.
 *
 * A lot of processing is deferred to here, so that the "p" command
 * can operate on an equal footing to mouse events.
 */
int
rt_edit_process_result(struct rt_edit *s)
{
    if (!s)
	return BRLCAD_ERROR;

    bu_clbk_t f = NULL;
    void *d = NULL;

    ++s->view_update_requested;

    int had_method = 0;
    const struct rt_db_internal *ip = &s->es_int;
    if (ip->idb_type < 0 || EDOBJ[ip->idb_type].magic != RT_FUNCTAB_MAGIC) {
	bu_vls_printf(s->log_str,
	    "rt_edit_process: invalid primitive type %d.\n", ip->idb_type);
	--s->view_update_requested;
	return BRLCAD_ERROR;
    }
    if (EDOBJ[ip->idb_type].ft_edit) {
	bu_vls_trunc(s->log_str, 0);
	if ((*EDOBJ[ip->idb_type].ft_edit)(s)) {
	    if (bu_vls_strlen(s->log_str)) {
		rt_edit_map_clbk_get(&f, &d, s->m, ECMD_PRINT_STR, BU_CLBK_DURING);
		if (f)
		    (*f)(0, NULL, d, NULL);
	    }
	    --s->view_update_requested;
	    s->e_inpara = 0;
	    s->e_nstr = 0;
	    memset(s->e_str, 0, sizeof(s->e_str));
	    s->e_mvalid = 0;
	    return BRLCAD_ERROR;
	}
	if (bu_vls_strlen(s->log_str)) {
	    rt_edit_map_clbk_get(&f, &d, s->m, ECMD_PRINT_STR, BU_CLBK_DURING);
	    if (f)
		(*f)(0, NULL, d, NULL);
	    bu_vls_trunc(s->log_str, 0);
	}
	had_method = 1;
    }

    switch (s->edit_flag) {

	case RT_EDIT_IDLE:
	    /* do nothing more */
	    --s->view_update_requested;
	    break;
	default:
	    {
		if (had_method)
		    break;
		struct bu_vls tmp_vls = BU_VLS_INIT_ZERO;
		bu_vls_sprintf(&tmp_vls, "%s", bu_vls_cstr(s->log_str));
		bu_vls_sprintf(s->log_str, "rt_edit_process:  unknown edflag = %d.\n", s->edit_flag);
		rt_edit_map_clbk_get(&f, &d, s->m, ECMD_PRINT_STR, BU_CLBK_DURING);
		if (f)
		    (*f)(0, NULL, d, NULL);
		bu_vls_sprintf(s->log_str, "%s", bu_vls_cstr(&tmp_vls));
		bu_vls_free(&tmp_vls);
		rt_edit_map_clbk_get(&f, &d, s->m, ECMD_PRINT_RESULTS, BU_CLBK_DURING);
		if (f)
		    (*f)(0, NULL, d, NULL);
		--s->view_update_requested;
		s->e_inpara = 0;
		s->e_nstr = 0;
		memset(s->e_str, 0, sizeof(s->e_str));
		s->e_mvalid = 0;
		return BRLCAD_ERROR;
	    }
    }

    /* If the keypoint changed location, find about it here */
    if (!s->e_keyfixed)
	rt_get_solid_keypoint(s, &s->e_keypoint, &s->e_keytag, s->e_mat);

    // eaxes_pos callback
    int flag = 0;
    f = NULL; d = NULL;
    rt_edit_map_clbk_get(&f, &d, s->m, ECMD_EAXES_POS, BU_CLBK_DURING);
    if (f)
	(*f)(0, NULL, d, &flag);

    // replot callback
    f = NULL; d = NULL;
    rt_edit_map_clbk_get(&f, &d, s->m, ECMD_REPLOT_EDITING_SOLID, BU_CLBK_DURING);
    if (f)
	(*f)(0, NULL, d, NULL);

    // view update callback
    if (s->view_update_requested) {
	f = NULL; d = NULL;
	rt_edit_map_clbk_get(&f, &d, s->m, ECMD_VIEW_UPDATE, BU_CLBK_DURING);
	if (f)
	    (*f)(0, NULL, d, NULL);
    }

    // Inputs processed, reset
    s->e_inpara = 0;
    s->e_nstr = 0;
    memset(s->e_str, 0, sizeof(s->e_str));
    s->e_mvalid = 0;
    return BRLCAD_OK;
}

void
rt_edit_process(struct rt_edit *s)
{
    (void)rt_edit_process_result(s);
}

void
rt_edit_snap_point(point2d_t pt, const struct rt_edit *s)
{
    if (!s || !s->snap.enabled || s->snap.spacing <= 0.0)
	return;
    fastf_t inv_sp = 1.0 / s->snap.spacing;
    pt[0] = floor(pt[0] * inv_sp + 0.5) * s->snap.spacing;
    pt[1] = floor(pt[1] * inv_sp + 0.5) * s->snap.spacing;
}


int
rt_edit_checkpoint(struct rt_edit *s)
{
    if (!s)
	return BRLCAD_ERROR;

    RT_CK_DB_INTERNAL(&s->es_int);

    /* Release any previous snapshot */
    bu_free_external(&s->es_ckpt);
    BU_EXTERNAL_INIT(&s->es_ckpt);

    if (rt_obj_export(&s->es_ckpt, &s->es_int, 1.0, s->dbip) < 0) {
	bu_vls_printf(s->log_str, "rt_edit_checkpoint: export failed\n");
	return BRLCAD_ERROR;
    }

    return BRLCAD_OK;
}


int
rt_edit_revert(struct rt_edit *s)
{
    if (!s)
	return BRLCAD_ERROR;

    if (!s->es_ckpt.ext_buf) {
	bu_vls_printf(s->log_str, "rt_edit_revert: no checkpoint saved\n");
	return BRLCAD_ERROR;
    }

    int type = s->es_int.idb_type;

    if (edit_internal_restore(s, &s->es_ckpt, type) != BRLCAD_OK) {
	bu_vls_printf(s->log_str, "rt_edit_revert: import failed\n");
	return BRLCAD_ERROR;
    }

    return BRLCAD_OK;
}


/* ======================================================================
 * rt_edit_prim_desc_to_json / rt_edit_type_to_json
 *
 * Serialise the ft_edit_desc() descriptor to JSON using bu_vls so that
 * GUI layers (qged etc.) can consume edit metadata at runtime.
 * ======================================================================
 */

/* Type-code string table */
static const char *
edit_param_type_str(int type)
{
    switch (type) {
	case RT_EDIT_PARAM_SCALAR:  return "scalar";
	case RT_EDIT_PARAM_INTEGER: return "integer";
	case RT_EDIT_PARAM_BOOLEAN: return "boolean";
	case RT_EDIT_PARAM_POINT:   return "point";
	case RT_EDIT_PARAM_VECTOR:  return "vector";
	case RT_EDIT_PARAM_STRING:  return "string";
	case RT_EDIT_PARAM_ENUM:    return "enum";
	case RT_EDIT_PARAM_COLOR:   return "color";
	case RT_EDIT_PARAM_MATRIX:  return "matrix";
	case RT_EDIT_PARAM_POINT2:  return "point2";
	case RT_EDIT_PARAM_VECTOR2: return "vector2";
	case RT_EDIT_PARAM_SCALAR_LIST:  return "scalar_list";
	case RT_EDIT_PARAM_INTEGER_LIST: return "integer_list";
	default:                    return "unknown";
    }
}

static const char *
edit_control_class_str(enum rt_edit_control_class control_class)
{
    switch (control_class) {
	case RT_EDIT_CONTROL_GENERATED:   return "generated";
	case RT_EDIT_CONTROL_ACTION:      return "action";
	case RT_EDIT_CONTROL_CUSTOM:      return "custom";
	case RT_EDIT_CONTROL_UNSUPPORTED: return "unsupported";
	case RT_EDIT_CONTROL_INHERIT:
	default:                          return "unclassified";
    }
}

static const char *
edit_selection_domain_str(enum rt_edit_selection_domain domain)
{
    switch (domain) {
	case RT_EDIT_SELECTION_PRIMITIVE:     return "primitive";
	case RT_EDIT_SELECTION_VERTEX:        return "vertex";
	case RT_EDIT_SELECTION_EDGE:          return "edge";
	case RT_EDIT_SELECTION_FACE:          return "face";
	case RT_EDIT_SELECTION_SEGMENT:       return "segment";
	case RT_EDIT_SELECTION_CONTROL_POINT: return "control_point";
	case RT_EDIT_SELECTION_CUSTOM:        return "custom";
	case RT_EDIT_SELECTION_NONE:
	default:                              return "none";
    }
}

static const char *
edit_coordinate_space_str(enum rt_edit_coordinate_space space)
{
    switch (space) {
	case RT_EDIT_COORDINATE_SCALAR:        return "scalar";
	case RT_EDIT_COORDINATE_OBJECT:        return "object";
	case RT_EDIT_COORDINATE_MODEL:         return "model";
	case RT_EDIT_COORDINATE_VIEW:          return "view";
	case RT_EDIT_COORDINATE_PARAMETRIC_2D: return "parametric_2d";
	case RT_EDIT_COORDINATE_INFER:
	default:                               return "infer";
    }
}

static const char *
edit_manipulator_hint_str(enum rt_edit_manipulator_hint hint)
{
    switch (hint) {
	case RT_EDIT_MANIPULATOR_POINT:         return "point";
	case RT_EDIT_MANIPULATOR_AXIS:          return "axis";
	case RT_EDIT_MANIPULATOR_PLANE:         return "plane";
	case RT_EDIT_MANIPULATOR_RADIUS:        return "radius";
	case RT_EDIT_MANIPULATOR_ROTATION_RING: return "rotation_ring";
	case RT_EDIT_MANIPULATOR_INDEXED_SET:   return "indexed_set";
	case RT_EDIT_MANIPULATOR_CUSTOM:        return "custom";
	case RT_EDIT_MANIPULATOR_NONE:
	default:                                return "none";
    }
}

static const char *
edit_param_semantic_str(enum rt_edit_param_semantic semantic)
{
    switch (semantic) {
	case RT_EDIT_SEMANTIC_INDEX:     return "index";
	case RT_EDIT_SEMANTIC_POSITION:  return "position";
	case RT_EDIT_SEMANTIC_DELTA:     return "delta";
	case RT_EDIT_SEMANTIC_DIRECTION: return "direction";
	case RT_EDIT_SEMANTIC_DISTANCE:  return "distance";
	case RT_EDIT_SEMANTIC_ANGLE:     return "angle";
	case RT_EDIT_SEMANTIC_SCALE:     return "scale";
	case RT_EDIT_SEMANTIC_FRACTION:  return "fraction";
	case RT_EDIT_SEMANTIC_COUNT:     return "count";
	case RT_EDIT_SEMANTIC_PROPERTY:  return "property";
	case RT_EDIT_SEMANTIC_INFER:
	default:                         return "infer";
    }
}

enum rt_edit_control_class
rt_edit_cmd_control_class(const struct rt_edit_prim_desc *desc,
	const struct rt_edit_cmd_desc *command)
{
    if (!desc || !command)
	return RT_EDIT_CONTROL_UNSUPPORTED;
    int command_index = -1;
    for (int ci = 0; ci < desc->ncmd; ci++) {
	if (desc->cmds[ci].cmd_id == command->cmd_id) {
	    command_index = ci;
	    break;
	}
    }
    if (command_index < 0)
	return RT_EDIT_CONTROL_UNSUPPORTED;
    if (desc->command_controls &&
	desc->command_controls[command_index] != RT_EDIT_CONTROL_INHERIT)
	return desc->command_controls[command_index];
    return desc->control_class;
}

struct rt_edit_interaction_desc
rt_edit_cmd_interaction(const struct rt_edit_prim_desc *desc,
	const struct rt_edit_cmd_desc *command)
{
    struct rt_edit_interaction_desc result = {
	RT_EDIT_SELECTION_NONE,
	RT_EDIT_COORDINATE_INFER,
	RT_EDIT_MANIPULATOR_NONE,
	NULL,
	0
    };
    if (!desc || !command || !desc->command_interactions)
	return result;
    for (int ci = 0; ci < desc->ncmd; ci++) {
	if (desc->cmds[ci].cmd_id == command->cmd_id)
	    return desc->command_interactions[ci];
    }
    return result;
}

enum rt_edit_param_semantic
rt_edit_param_semantic_get(const struct rt_edit_prim_desc *desc,
	const struct rt_edit_cmd_desc *command, int parameter_index)
{
    if (!desc || !command || parameter_index < 0 ||
	parameter_index >= command->nparam)
	return RT_EDIT_SEMANTIC_INFER;
    const struct rt_edit_interaction_desc interaction =
	rt_edit_cmd_interaction(desc, command);
    if (!interaction.parameter_semantics ||
	parameter_index >= interaction.parameter_semantic_count)
	return RT_EDIT_SEMANTIC_INFER;
    return interaction.parameter_semantics[parameter_index];
}

static std::string
edit_cmd_slug(const char *label)
{
    std::string slug;
    for (const char *p = label; p && *p; p++) {
	if (isalnum((unsigned char)*p))
	    slug += (char)tolower((unsigned char)*p);
	else if (!slug.empty() && slug.back() != '_')
	    slug += '_';
    }
    while (!slug.empty() && slug.back() == '_')
	slug.pop_back();
    return slug;
}

static std::string
edit_cmd_machine_name(const struct rt_edit_prim_desc *desc,
	const struct rt_edit_cmd_desc *command)
{
    if (!desc || !command)
	return std::string();
    std::string name = edit_cmd_slug(command->name);
    const std::string global_prefix("ecmd_");
    if (name.compare(0, global_prefix.size(), global_prefix) == 0)
	name.erase(0, global_prefix.size());
    const std::string primitive_prefix = edit_cmd_slug(desc->prim_type) + "_";
    if (name.compare(0, primitive_prefix.size(), primitive_prefix) == 0)
	name.erase(0, primitive_prefix.size());
    return name;
}

int
rt_edit_cmd_name(struct bu_vls *out, const struct rt_edit_prim_desc *desc,
		 const struct rt_edit_cmd_desc *command)
{
    if (!out || !desc || !command)
	return BRLCAD_ERROR;
    const std::string name = edit_cmd_machine_name(desc, command);
    if (name.empty())
	return BRLCAD_ERROR;
    bu_vls_strcat(out, name.c_str());
    return BRLCAD_OK;
}

static bool
edit_limit_is_set(fastf_t val)
{
    static_assert(sizeof(fastf_t) == sizeof(double),
	"fastf_t must be double for RT_EDIT_PARAM_NO_LIMIT comparison");
    static const double sentinel = RT_EDIT_PARAM_NO_LIMIT;
    const double value = (double)val;
    return memcmp(&value, &sentinel, sizeof(double)) != 0;
}

int
rt_edit_param_bounds_get(const struct rt_edit *edit, int command_id,
	int parameter_index, struct rt_edit_param_bounds *bounds)
{
    if (bounds)
	memset(bounds, 0, sizeof(*bounds));
    if (!edit || !bounds || parameter_index < 0 ||
	edit->es_int.idb_type < 0 ||
	EDOBJ[edit->es_int.idb_type].magic != RT_FUNCTAB_MAGIC ||
	!EDOBJ[edit->es_int.idb_type].ft_edit_desc)
	return BRLCAD_ERROR;

    const struct rt_edit_prim_desc *desc =
	EDOBJ[edit->es_int.idb_type].ft_edit_desc();
    const struct rt_edit_cmd_desc *command = NULL;
    for (int ci = 0; desc && ci < desc->ncmd; ci++) {
	if (desc->cmds[ci].cmd_id == command_id) {
	    command = &desc->cmds[ci];
	    break;
	}
    }
    if (!command || parameter_index >= command->nparam)
	return BRLCAD_ERROR;

    const struct rt_edit_param_desc *param =
	&command->params[parameter_index];
    if (edit_limit_is_set(param->range_min)) {
	bounds->minimum = param->range_min;
	bounds->has_minimum = 1;
    }
    if (edit_limit_is_set(param->range_max)) {
	bounds->maximum = param->range_max;
	bounds->has_maximum = 1;
    }
    if (desc->parameter_bounds && desc->parameter_bounds(bounds, edit,
	    command_id, parameter_index) != BRLCAD_OK)
	return BRLCAD_ERROR;
    if (bounds->has_minimum && bounds->has_maximum &&
	bounds->minimum > bounds->maximum)
	return BRLCAD_ERROR;
    return BRLCAD_OK;
}

static bool
edit_machine_name_valid(const char *name)
{
    if (!name || !name[0])
	return false;
    for (const unsigned char *p = (const unsigned char *)name; *p; p++) {
	if (!(std::isalnum(*p) || *p == '_'))
	    return false;
    }
    return true;
}

int
rt_edit_prim_desc_validate(struct bu_vls *msg,
			   const struct rt_edit_prim_desc *desc)
{
#define DESC_FAIL(_fmt, ...) do { \
    if (msg) bu_vls_printf(msg, (_fmt), __VA_ARGS__); \
    return BRLCAD_ERROR; \
} while (0)

    if (!desc) {
	if (msg)
	    bu_vls_strcat(msg, "null primitive edit descriptor");
	return BRLCAD_ERROR;
    }
    if (!edit_machine_name_valid(desc->prim_type))
	DESC_FAIL("invalid primitive machine name '%s'",
	    desc->prim_type ? desc->prim_type : "(null)");
    if (!desc->prim_label || !desc->prim_label[0])
	DESC_FAIL("%s: missing primitive label", desc->prim_type);
    if (desc->ncmd < 0 || (desc->ncmd && !desc->cmds))
	DESC_FAIL("%s: invalid command array", desc->prim_type);
    if (desc->nopt < 0 || (desc->nopt && !desc->opts))
	DESC_FAIL("%s: invalid option array", desc->prim_type);
    if (desc->control_class < RT_EDIT_CONTROL_GENERATED ||
	desc->control_class > RT_EDIT_CONTROL_UNSUPPORTED)
	DESC_FAIL("%s: invalid default control class %d", desc->prim_type,
	    (int)desc->control_class);

    std::set<int> command_ids;
    std::set<std::string> command_names;
    std::set<std::string> canonical_names;
    for (int ci = 0; ci < desc->ncmd; ci++) {
	const struct rt_edit_cmd_desc *cmd = &desc->cmds[ci];
	if (!command_ids.insert(cmd->cmd_id).second)
	    DESC_FAIL("%s: duplicate command id %d", desc->prim_type,
		cmd->cmd_id);
	if (!edit_machine_name_valid(cmd->name))
	    DESC_FAIL("%s: command %d has invalid machine name '%s'",
		desc->prim_type, cmd->cmd_id,
		cmd->name ? cmd->name : "(null)");
	if (!command_names.insert(cmd->name).second)
	    DESC_FAIL("%s: duplicate command name '%s'", desc->prim_type,
		cmd->name);
	const std::string canonical_name = edit_cmd_machine_name(desc, cmd);
	if (canonical_name.empty() ||
	    !canonical_names.insert(canonical_name).second)
	    DESC_FAIL("%s: invalid or duplicate canonical command name '%s'",
		desc->prim_type, canonical_name.c_str());
	if (!cmd->label || !cmd->label[0])
	    DESC_FAIL("%s/%s: missing command label", desc->prim_type,
		cmd->name);
	if (!cmd->category || !cmd->category[0])
	    DESC_FAIL("%s/%s: missing command category", desc->prim_type,
		cmd->name);
	if (cmd->interactive != 0 && cmd->interactive != 1)
	    DESC_FAIL("%s/%s: invalid interactive flag %d", desc->prim_type,
		cmd->name, cmd->interactive);
	if (cmd->nparam < 0 || (cmd->nparam && !cmd->params))
	    DESC_FAIL("%s/%s: invalid parameter array", desc->prim_type,
		cmd->name);
	if (desc->command_controls &&
	    (desc->command_controls[ci] < RT_EDIT_CONTROL_INHERIT ||
	     desc->command_controls[ci] > RT_EDIT_CONTROL_UNSUPPORTED))
	    DESC_FAIL("%s/%s: invalid control class %d", desc->prim_type,
		cmd->name, (int)desc->command_controls[ci]);
	if (desc->command_interactions) {
	    const struct rt_edit_interaction_desc *interaction =
		&desc->command_interactions[ci];
	    if (interaction->selection_domain < RT_EDIT_SELECTION_NONE ||
		interaction->selection_domain > RT_EDIT_SELECTION_CUSTOM)
		DESC_FAIL("%s/%s: invalid selection domain %d",
		    desc->prim_type, cmd->name,
		    (int)interaction->selection_domain);
	    if (interaction->coordinate_space < RT_EDIT_COORDINATE_INFER ||
		interaction->coordinate_space >
		RT_EDIT_COORDINATE_PARAMETRIC_2D)
		DESC_FAIL("%s/%s: invalid coordinate space %d",
		    desc->prim_type, cmd->name,
		    (int)interaction->coordinate_space);
	    if (interaction->coordinate_space == RT_EDIT_COORDINATE_INFER)
		DESC_FAIL("%s/%s: interaction coordinate space must be explicit",
		    desc->prim_type, cmd->name);
	    if (interaction->manipulator < RT_EDIT_MANIPULATOR_NONE ||
		interaction->manipulator > RT_EDIT_MANIPULATOR_CUSTOM)
		DESC_FAIL("%s/%s: invalid manipulator hint %d",
		    desc->prim_type, cmd->name,
		    (int)interaction->manipulator);
	    if (interaction->parameter_semantic_count < 0 ||
		(interaction->parameter_semantic_count &&
		 !interaction->parameter_semantics) ||
		interaction->parameter_semantic_count > cmd->nparam)
		DESC_FAIL("%s/%s: invalid parameter semantic array",
		    desc->prim_type, cmd->name);
	    for (int si = 0; si < interaction->parameter_semantic_count; si++) {
		if (interaction->parameter_semantics[si] <
			RT_EDIT_SEMANTIC_INFER ||
		    interaction->parameter_semantics[si] >
			RT_EDIT_SEMANTIC_PROPERTY)
		    DESC_FAIL("%s/%s: invalid parameter semantic %d",
			desc->prim_type, cmd->name,
			(int)interaction->parameter_semantics[si]);
	    }
	}

	std::set<std::string> parameter_names;
	for (int pi = 0; pi < cmd->nparam; pi++) {
	    const struct rt_edit_param_desc *param = &cmd->params[pi];
	    if (!edit_machine_name_valid(param->name))
		DESC_FAIL("%s/%s: parameter %d has invalid name '%s'",
		    desc->prim_type, cmd->name, pi,
		    param->name ? param->name : "(null)");
	    if (!parameter_names.insert(param->name).second)
		DESC_FAIL("%s/%s: duplicate parameter name '%s'",
		    desc->prim_type, cmd->name, param->name);
	    if (!param->label || !param->label[0])
		DESC_FAIL("%s/%s/%s: missing parameter label",
		    desc->prim_type, cmd->name, param->name);

	    int width = 1;
	    const bool list_param =
		param->type == RT_EDIT_PARAM_SCALAR_LIST ||
		param->type == RT_EDIT_PARAM_INTEGER_LIST;
	    switch (param->type) {
		case RT_EDIT_PARAM_POINT2:
		case RT_EDIT_PARAM_VECTOR2:
		    width = 2;
		    break;
		case RT_EDIT_PARAM_POINT:
		case RT_EDIT_PARAM_VECTOR:
		case RT_EDIT_PARAM_COLOR:
		    width = 3;
		    break;
		case RT_EDIT_PARAM_MATRIX:
		    width = 16;
		    break;
		case RT_EDIT_PARAM_SCALAR:
		case RT_EDIT_PARAM_INTEGER:
		case RT_EDIT_PARAM_BOOLEAN:
		case RT_EDIT_PARAM_STRING:
		case RT_EDIT_PARAM_ENUM:
		    break;
		case RT_EDIT_PARAM_SCALAR_LIST:
		case RT_EDIT_PARAM_INTEGER_LIST:
		    if (pi != cmd->nparam - 1)
			DESC_FAIL("%s/%s/%s: a variable list must be the final parameter",
			    desc->prim_type, cmd->name, param->name);
		    if (param->nenum <= 0 ||
			param->index + param->nenum > RT_EDIT_MAXPARA)
			DESC_FAIL("%s/%s/%s: invalid minimum list count %d",
			    desc->prim_type, cmd->name, param->name,
			    param->nenum);
		    width = RT_EDIT_MAXPARA - param->index;
		    break;
		default:
		    DESC_FAIL("%s/%s/%s: invalid parameter type %d",
			desc->prim_type, cmd->name, param->name, param->type);
	    }

	    const int limit = param->type == RT_EDIT_PARAM_STRING ?
		RT_EDIT_MAXSTR : RT_EDIT_MAXPARA;
	    if (param->index < 0 || param->index + width > limit)
		DESC_FAIL("%s/%s/%s: parameter slots [%d,%d) exceed %d",
		    desc->prim_type, cmd->name, param->name, param->index,
		    param->index + width, limit);
	    if (edit_limit_is_set(param->range_min) &&
		!std::isfinite((double)param->range_min))
		DESC_FAIL("%s/%s/%s: nonfinite minimum", desc->prim_type,
		    cmd->name, param->name);
	    if (edit_limit_is_set(param->range_max) &&
		!std::isfinite((double)param->range_max))
		DESC_FAIL("%s/%s/%s: nonfinite maximum", desc->prim_type,
		    cmd->name, param->name);
	    if (!list_param && param->type != RT_EDIT_PARAM_ENUM && param->nenum)
		DESC_FAIL("%s/%s/%s: enum/list count on an unrelated parameter",
		    desc->prim_type, cmd->name, param->name);
	    if (edit_limit_is_set(param->range_min) &&
		edit_limit_is_set(param->range_max) &&
		param->range_min > param->range_max)
		DESC_FAIL("%s/%s/%s: minimum exceeds maximum",
		    desc->prim_type, cmd->name, param->name);
	    if (param->type == RT_EDIT_PARAM_ENUM) {
		if (param->nenum <= 0 || !param->enum_labels || !param->enum_ids)
		    DESC_FAIL("%s/%s/%s: incomplete enum metadata",
			desc->prim_type, cmd->name, param->name);
		std::set<int> enum_ids;
		for (int ei = 0; ei < param->nenum; ei++) {
		    if (!param->enum_labels[ei] || !param->enum_labels[ei][0])
			DESC_FAIL("%s/%s/%s: empty enum label %d",
			    desc->prim_type, cmd->name, param->name, ei);
		    if (!enum_ids.insert(param->enum_ids[ei]).second)
			DESC_FAIL("%s/%s/%s: duplicate enum id %d",
			    desc->prim_type, cmd->name, param->name,
			    param->enum_ids[ei]);
		}
	    }
	    if (param->type == RT_EDIT_PARAM_STRING &&
		(!param->prim_field || !param->prim_field[0]))
		DESC_FAIL("%s/%s/%s: missing string primitive field",
		    desc->prim_type, cmd->name, param->name);
	}
    }

    std::set<std::string> option_names;
    for (int oi = 0; oi < desc->nopt; oi++) {
	const struct rt_edit_opt_desc *opt = &desc->opts[oi];
	if (!edit_machine_name_valid(opt->name) ||
	    !option_names.insert(opt->name).second)
	    DESC_FAIL("%s: invalid or duplicate option name '%s'",
		desc->prim_type, opt->name ? opt->name : "(null)");
	if (!opt->label || !opt->label[0] || !opt->desc || !opt->desc[0])
	    DESC_FAIL("%s/%s: incomplete option metadata", desc->prim_type,
		opt->name);
	if (opt->type < RT_EDIT_PARAM_SCALAR || opt->type > RT_EDIT_PARAM_MATRIX)
	    DESC_FAIL("%s/%s: invalid option type %d", desc->prim_type,
		opt->name, opt->type);
    }

#undef DESC_FAIL
    return BRLCAD_OK;
}

/* Emit a fastf_t value as JSON number or "null" for NO_LIMIT */
static void
emit_limit(struct bu_vls *out, fastf_t val)
{
    /* Use bit-level comparison to avoid -Wfloat-equal.  RT_EDIT_PARAM_NO_LIMIT
     * is defined as (-DBL_MAX), a well-defined sentinel bit-pattern.
     * fastf_t is always double in BRL-CAD (vmath.h), so sizeof(fastf_t) and
     * sizeof(double) are identical; we assert that here for safety. */
    /* fastf_t is always double in BRL-CAD (vmath.h); assert for safety. */
    double v = (double)val;
    if (!edit_limit_is_set(val))
	bu_vls_strcat(out, "null");
    else
	bu_vls_printf(out, "%.17g", v);
}

/* Escape a C string for safe JSON embedding (handles " and \) */
static void
emit_json_str(struct bu_vls *out, const char *s)
{
    bu_vls_putc(out, '"');
    if (s) {
	for (; *s; s++) {
	    if (*s == '"' || *s == '\\')
		bu_vls_putc(out, '\\');
	    bu_vls_putc(out, *s);
	}
    }
    bu_vls_putc(out, '"');
}

int
rt_edit_prim_desc_to_json(struct bu_vls *out,
                          const struct rt_edit_prim_desc *desc)
{
    if (!out || !desc)
	return BRLCAD_ERROR;

    if (rt_edit_prim_desc_validate(NULL, desc) != BRLCAD_OK)
	return BRLCAD_ERROR;

    bu_vls_strcat(out, "{\n");
    bu_vls_strcat(out, "  \"prim_type\": ");
    emit_json_str(out, desc->prim_type);
    bu_vls_strcat(out, ",\n");

    bu_vls_strcat(out, "  \"prim_label\": ");
    emit_json_str(out, desc->prim_label);
    bu_vls_strcat(out, ",\n");

    bu_vls_strcat(out, "  \"commands\": [\n");

    for (int ci = 0; ci < desc->ncmd; ci++) {
	const struct rt_edit_cmd_desc *cmd = &desc->cmds[ci];
	bu_vls_strcat(out, "    {\n");
	bu_vls_printf(out, "      \"cmd_id\": %d,\n", cmd->cmd_id);
	const std::string machine_name = edit_cmd_machine_name(desc, cmd);
	bu_vls_strcat(out, "      \"name\": ");
	emit_json_str(out, machine_name.c_str());
	bu_vls_strcat(out, ",\n");
	bu_vls_strcat(out, "      \"symbol\": ");
	emit_json_str(out, cmd->name);
	bu_vls_strcat(out, ",\n");
	std::string cname = edit_cmd_slug(cmd->label);
	bu_vls_strcat(out, "      \"aliases\": [");
	bool have_alias = false;
	if (cname != machine_name) {
	    emit_json_str(out, cname.c_str());
	    have_alias = true;
	}
	if (cmd->nparam == 1 && cmd->params && cmd->params[0].name) {
	    const char *a = cmd->params[0].name;
	    if (a && !BU_STR_EQUAL(a, machine_name.c_str()) &&
		!BU_STR_EQUAL(a, cname.c_str())) {
		if (have_alias)
		    bu_vls_putc(out, ',');
		emit_json_str(out, a);
		have_alias = true;
	    }
	}
	bu_vls_strcat(out, "],\n");
	bu_vls_strcat(out, "      \"label\": ");
	emit_json_str(out, cmd->label);
	bu_vls_strcat(out, ",\n");
	bu_vls_strcat(out, "      \"category\": ");
	emit_json_str(out, cmd->category);
	bu_vls_strcat(out, ",\n");
	bu_vls_strcat(out, "      \"control\": ");
	emit_json_str(out, edit_control_class_str(
	    rt_edit_cmd_control_class(desc, cmd)));
	bu_vls_strcat(out, ",\n");
	const struct rt_edit_interaction_desc interaction =
	    rt_edit_cmd_interaction(desc, cmd);
	bu_vls_strcat(out, "      \"selection_domain\": ");
	emit_json_str(out, edit_selection_domain_str(
	    interaction.selection_domain));
	bu_vls_strcat(out, ",\n      \"coordinate_space\": ");
	emit_json_str(out, edit_coordinate_space_str(
	    interaction.coordinate_space));
	bu_vls_strcat(out, ",\n      \"manipulator\": ");
	emit_json_str(out, edit_manipulator_hint_str(
	    interaction.manipulator));
	bu_vls_strcat(out, ",\n");
	bu_vls_printf(out, "      \"interactive\": %s,\n",
		      cmd->interactive ? "true" : "false");
	bu_vls_printf(out, "      \"display_order\": %d,\n",
		      cmd->display_order);
	bu_vls_strcat(out, "      \"params\": [\n");

	for (int pi = 0; pi < cmd->nparam; pi++) {
	    const struct rt_edit_param_desc *p = &cmd->params[pi];
	    bu_vls_strcat(out, "        {\n");
	    bu_vls_strcat(out, "          \"name\": ");
	    emit_json_str(out, p->name);
	    bu_vls_strcat(out, ",\n");
	    bu_vls_strcat(out, "          \"label\": ");
	    emit_json_str(out, p->label);
	    bu_vls_strcat(out, ",\n");
	    bu_vls_printf(out, "          \"type\": \"%s\",\n",
			 edit_param_type_str(p->type));
	    bu_vls_strcat(out, "          \"semantic\": ");
	    emit_json_str(out, edit_param_semantic_str(
		rt_edit_param_semantic_get(desc, cmd, pi)));

	    /* index — omitted for STRING */
	    if (p->type != RT_EDIT_PARAM_STRING)
		bu_vls_printf(out, ",\n          \"index\": %d", p->index);

	    /* range — only for numeric types */
	    if (p->type == RT_EDIT_PARAM_SCALAR ||
		p->type == RT_EDIT_PARAM_INTEGER ||
		p->type == RT_EDIT_PARAM_ENUM   ||
		p->type == RT_EDIT_PARAM_COLOR  ||
		p->type == RT_EDIT_PARAM_SCALAR_LIST ||
		p->type == RT_EDIT_PARAM_INTEGER_LIST) {
		bu_vls_strcat(out, ",\n          \"min\": ");
		emit_limit(out, p->range_min);
		bu_vls_strcat(out, ",\n          \"max\": ");
		emit_limit(out, p->range_max);
	    }

	    if (p->type == RT_EDIT_PARAM_SCALAR_LIST ||
		p->type == RT_EDIT_PARAM_INTEGER_LIST) {
		bu_vls_printf(out,
		    ",\n          \"min_count\": %d,\n          \"max_count\": %d",
		    p->nenum, RT_EDIT_MAXPARA - p->index);
	    }

	    /* units */
	    bu_vls_strcat(out, ",\n          \"units\": ");
	    if (p->units)
		emit_json_str(out, p->units);
	    else
		bu_vls_strcat(out, "null");

	    /* enum extras */
	    if (p->type == RT_EDIT_PARAM_ENUM && p->nenum > 0) {
		bu_vls_strcat(out, ",\n          \"enum_labels\": [");
		for (int ei = 0; ei < p->nenum; ei++) {
		    if (ei) bu_vls_putc(out, ',');
		    emit_json_str(out, p->enum_labels[ei]);
		}
		bu_vls_strcat(out, "],\n          \"enum_ids\": [");
		for (int ei = 0; ei < p->nenum; ei++) {
		    if (ei) bu_vls_putc(out, ',');
		    bu_vls_printf(out, "%d", p->enum_ids[ei]);
		}
		bu_vls_putc(out, ']');
	    }

	    /* prim_field for STRING */
	    if (p->type == RT_EDIT_PARAM_STRING) {
		bu_vls_strcat(out, ",\n          \"prim_field\": ");
		emit_json_str(out, p->prim_field);
	    }

	    bu_vls_strcat(out, "\n        }");
	    if (pi < cmd->nparam - 1)
		bu_vls_putc(out, ',');
	    bu_vls_putc(out, '\n');
	}

	bu_vls_strcat(out, "      ]\n");
	bu_vls_strcat(out, "    }");
	if (ci < desc->ncmd - 1)
	    bu_vls_putc(out, ',');
	bu_vls_putc(out, '\n');
    }

    bu_vls_strcat(out, "  ]\n");
    bu_vls_strcat(out, "}\n");

    return BRLCAD_OK;
}


int
rt_edit_type_to_json(struct bu_vls *out, int prim_type_id)
{
    if (!out)
	return BRLCAD_ERROR;

    /* Range-check against the EDOBJ table size */
    extern const struct rt_edit_functab EDOBJ[];
    /* Walk until we find a matching type label or a sentinel entry */
    for (int i = 0; EDOBJ[i].magic == RT_FUNCTAB_MAGIC; i++) {
	if (i == prim_type_id) {
	    if (!EDOBJ[i].ft_edit_desc)
		return BRLCAD_ERROR;
	    return rt_edit_prim_desc_to_json(out, (*EDOBJ[i].ft_edit_desc)());
	}
    }

    return BRLCAD_ERROR;
}

int
rt_edit_set_opt(struct rt_edit *s, const char *key, const char *val)
{
    if (!s || !key || !val) return BRLCAD_ERROR;
    return bu_avs_add(&s->options, key, val);
}

const char *
rt_edit_get_opt(struct rt_edit *s, const char *key)
{
    if (!s || !key) return NULL;
    return bu_avs_get(&s->options, key);
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
