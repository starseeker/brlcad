/*                         D R A W . C P P
 * BRL-CAD
 *
 * Copyright (c) 2008-2026 United States Government as represented by
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
/** @file draw.cpp
 *
 * Drawing routines that generate scene objects.
 *
 */

#include "common.h"

#include <set>
#include <unordered_map>

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include "bsocket.h"

#include "bu/cmd.h"
#include "bu/opt.h"
#include "bu/sort.h"
#include "bu/str.h"

#include "../librt/librt_private.h"
#include "./ged_private.h"

static void
tree_color(struct directory *dp, struct draw_data_t *dd)
{
    struct bu_attribute_value_set c_avs = BU_AVS_INIT_ZERO;

    // Easy answer - if we're overridden, dd color is already set.
    if (dd->vs && dd->vs->color_override)
	return;

    // Not overridden by settings.  Next question - are we under an inherit?
    // If so, dd color is already set.
    if (dd->color_inherit)
	return;

    // Need attributes for the rest of this
    db5_get_attributes(dd->dbip, &c_avs, dp);

    // No inherit.  Do we have a region material table?
    if (db_mater_head(dd->dbip) != MATER_NULL) {
	// If we do, do we have a region id?
	int region_id = -1;
	const char *region_id_val = bu_avs_get(&c_avs, "region_id");
	if (region_id_val) {
	    bu_opt_int(NULL, 1, &region_id_val, (void *)&region_id);
	} else if (dp->d_flags & RT_DIR_REGION) {
	    // If we have a region flag but no region_id, for color table
	    // purposes treat the region_id as 0
	    region_id = 0;
	}
	if (region_id >= 0) {
	    const struct mater *mp;
	    int material_color = 0;
	    for (mp = db_mater_head(dd->dbip); mp != MATER_NULL; mp = mp->mt_forw) {
		if (region_id > mp->mt_high || region_id < mp->mt_low) {
		    continue;
		}
		unsigned char mt[3];
		mt[0] = mp->mt_r;
		mt[1] = mp->mt_g;
		mt[2] = mp->mt_b;
		bu_color_from_rgb_chars(&dd->c, mt);
		material_color = 1;
	    }
	    if (material_color) {
		// Have answer from color table
		bu_avs_free(&c_avs);
		return;
	    }
	}
    }

    // Material table didn't give us the answer - do we have a color or
    // rgb attribute?
    const char *color_val = bu_avs_get(&c_avs, "color");
    if (!color_val) {
	color_val = bu_avs_get(&c_avs, "rgb");
    }
    if (color_val) {
	bu_opt_color(NULL, 1, &color_val, (void *)&dd->c);
    }

    // Check for an inherit flag
    dd->color_inherit = (BU_STR_EQUAL(bu_avs_get(&c_avs, "inherit"), "1")) ? 1 : 0;

    // Done with attributes
    bu_avs_free(&c_avs);
}

/*****************************************************************************

  The primary drawing subtree walk.

  To get an initial idea of scene size and center for adaptive plotting (i.e.
  before we have wireframes drawn) we also have a need to check ft_bbox ahead
  of the vlist generation.  It can be thought of as a "preliminary autoview"
  step.  That mode is also supported by this subtree walk.

******************************************************************************/
static void
draw_walk_tree(struct db_full_path *path, union tree *tp, mat_t *curr_mat,
			 void (*traverse_func) (struct db_full_path *path, mat_t *, void *),
			 void *client_data, void *comb_inst_map)
{
    mat_t om, nm;
    struct directory *dp;
    struct draw_data_t *dd= (struct draw_data_t *)client_data;
    std::unordered_map<std::string, int> *cinst_map = (std::unordered_map<std::string, int> *)comb_inst_map;

    if (!tp)
	return;

    RT_CK_FULL_PATH(path);
    RT_CHECK_DBI(dd->dbip);
    RT_CK_TREE(tp);

    switch (tp->tr_op) {
	case OP_NOT:
	case OP_GUARD:
	case OP_XNOP:
	    draw_walk_tree(path, tp->tr_b.tb_left, curr_mat, traverse_func, client_data, comb_inst_map);
	    break;
	case OP_UNION:
	case OP_INTERSECT:
	case OP_SUBTRACT:
	case OP_XOR:
	    draw_walk_tree(path, tp->tr_b.tb_left, curr_mat, traverse_func, client_data, comb_inst_map);
	    dd->bool_op = tp->tr_op;
	    draw_walk_tree(path, tp->tr_b.tb_right, curr_mat, traverse_func, client_data, comb_inst_map);
	    break;
	case OP_DB_LEAF:
	    if (UNLIKELY(dd->dbip->i->dbi_use_comb_instance_ids && cinst_map))
		(*cinst_map)[std::string(tp->tr_l.tl_name)]++;
	    if ((dp=db_lookup(dd->dbip, tp->tr_l.tl_name, LOOKUP_QUIET)) == RT_DIR_NULL) {
		return;
	    } else {

		/* Update current matrix state to reflect the new branch of
		 * the tree. Either we have a local matrix, or we have an
		 * implicit IDN matrix. */
		MAT_COPY(om, *curr_mat);
		if (tp->tr_l.tl_mat) {
		    MAT_COPY(nm, tp->tr_l.tl_mat);
		} else {
		    MAT_IDN(nm);
		}
		bn_mat_mul(*curr_mat, om, nm);

		// Stash current color settings and see if we're getting new ones
		struct bu_color oc;
		int inherit_old = dd->color_inherit;
		HSET(oc.buc_rgb, dd->c.buc_rgb[0], dd->c.buc_rgb[1], dd->c.buc_rgb[2], dd->c.buc_rgb[3]);
		if (!dd->bound_only) {
		    tree_color(dp, dd);
		}

		// Two things may prevent further processing - a hidden dp, or
		// a cyclic path.  Can check either here or in traverse_func -
		// just do it here since otherwise the logic would have to be
		// duplicated in all traverse functions.
		if (!(dp->d_flags & RT_DIR_HIDDEN)) {
		    db_add_node_to_full_path(path, dp);
		    DB_FULL_PATH_SET_CUR_BOOL(path, tp->tr_op);
		    if (UNLIKELY(dd->dbip->i->dbi_use_comb_instance_ids && cinst_map))
			DB_FULL_PATH_SET_CUR_COMB_INST(path, (*cinst_map)[std::string(tp->tr_l.tl_name)]-1);
		    if (!db_full_path_cyclic(path, NULL, 0)) {
			/* Keep going */
			traverse_func(path, curr_mat, client_data);
		    }
		}

		/* Done with branch - restore path, put back the old matrix state,
		 * and restore previous color settings */
		DB_FULL_PATH_POP(path);
		MAT_COPY(*curr_mat, om);
		if (!dd->bound_only) {
		    dd->color_inherit = inherit_old;
		    HSET(dd->c.buc_rgb, oc.buc_rgb[0], oc.buc_rgb[1], oc.buc_rgb[2], oc.buc_rgb[3]);
		}
		return;
	    }

	default:
	    bu_log("db_functree_subtree: unrecognized operator %d\n", tp->tr_op);
	    bu_bomb("db_functree_subtree: unrecognized operator\n");
    }
}

/**
 * This walker builds a list of db_full_path entries corresponding to
 * the contents of the tree under *path.  It does so while assigning
 * the boolean operation associated with each path entry to the
 * db_full_path structure.  This list is then used for further
 * processing and filtering by the search routines.
 */
extern "C" void
draw_gather_paths(struct db_full_path *path, mat_t *curr_mat, void *client_data)
{
    struct directory *dp;
    struct draw_data_t *dd= (struct draw_data_t *)client_data;
    RT_CK_FULL_PATH(path);
    RT_CK_DBI(dd->dbip);

    dp = DB_FULL_PATH_CUR_DIR(path);
    if (!dp)
	return;

    // If we're skipping subtractions and we have a subtraction op there's no
    // point in going further.
    if (dd->vs && dd->vs->draw_non_subtract_only && dd->bool_op == 4) {
	return;
    }


    if (dp->d_flags & RT_DIR_COMB) {

	struct rt_db_internal in;
	struct rt_comb_internal *comb;

	if (rt_db_get_internal(&in, dp, dd->dbip, NULL) < 0)
	    return;

	comb = (struct rt_comb_internal *)in.idb_ptr;
	if (UNLIKELY(dd->dbip->i->dbi_use_comb_instance_ids)) {
	    std::unordered_map<std::string, int> cinst_map;
	    draw_walk_tree(path, comb->tree, curr_mat, draw_gather_paths, client_data, (void *)&cinst_map);
	} else {
	    draw_walk_tree(path, comb->tree, curr_mat, draw_gather_paths, client_data, NULL);
	}
	rt_db_free_internal(&in);

    } else {

	// If we've got a solid, things get interesting.  There are a lot of
	// potentially relevant options to sort through.  It may be that most
	// will end up getting handled by the object update callbacks, and the
	// job here will just be to set up the key data for later use...

	// Let the object know about its size
	int has_draw_size = 0;
	fastf_t draw_size = 0.0;
	if (dd->s_size && dd->s_size->find(DB_FULL_PATH_CUR_DIR(path)) != dd->s_size->end()) {
	    has_draw_size = 1;
	    draw_size = (*dd->s_size)[DB_FULL_PATH_CUR_DIR(path)];
	}

	unsigned char rgb[3];
	bu_color_to_rgb_chars(&dd->c, rgb);
	(void)ged_draw_source_group_ref_commit_database_leaf_draft(dd->gedp,
		dd->draw_parent_group_ref, dd->view_ctx, dd->dbip, path,
		*curr_mat, dd->tol, dd->ttol, dd->vs, dd->bool_op, rgb,
		has_draw_size, draw_size);

    }
}


extern "C" int
ged_draw_view_context_gobject_create(struct ged *gedp,
				     void *view_ctx,
				     const char *db_path,
				     const char *gobject_name,
				     struct bu_vls *result)
{
    if (!gedp || !gedp->dbip || !view_ctx || !db_path || !db_path[0] ||
	    !gobject_name || !gobject_name[0])
	return 0;

    if (ged_draw_view_context_feature_exists(view_ctx, gobject_name)) {
	if (result)
	    bu_vls_printf(result, "View feature %s already exists\n", gobject_name);
	return 0;
    }

    struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
    if (!wdbp)
	return 0;

    struct db_full_path *fp;
    BU_GET(fp, struct db_full_path);
    db_full_path_init(fp);
    if (db_string_to_path(fp, gedp->dbip, db_path) < 0) {
	db_free_full_path(fp);
	BU_PUT(fp, struct db_full_path);
	if (result)
	    bu_vls_printf(result, "Invalid path: %s\n", db_path);
	return 0;
    }

    mat_t mat;
    MAT_IDN(mat);
    if (!db_path_to_mat(gedp->dbip, fp, mat, fp->fp_len - 1)) {
	db_free_full_path(fp);
	BU_PUT(fp, struct db_full_path);
	if (result)
	    bu_vls_printf(result, "Invalid path matrix: %s\n", db_path);
	return 0;
    }

    struct rt_db_internal *ip;
    BU_GET(ip, struct rt_db_internal);
    RT_DB_INTERNAL_INIT(ip);
    if (rt_db_get_internal(ip, DB_FULL_PATH_CUR_DIR(fp), gedp->dbip, mat) < 0) {
	db_free_full_path(fp);
	BU_PUT(fp, struct db_full_path);
	rt_db_free_internal(ip);
	BU_PUT(ip, struct rt_db_internal);
	return 0;
    }

    ged_draw_group_ref draw_parent_group_ref = GED_DRAW_GROUP_REF_NULL;
    if (!ged_draw_view_context_overlay_internal_create_group_ref(gedp, view_ctx,
	    gobject_name, fp, &ip, &draw_parent_group_ref)) {
	db_free_full_path(fp);
	BU_PUT(fp, struct db_full_path);
	rt_db_free_internal(ip);
	BU_PUT(ip, struct rt_db_internal);
	return 0;
    }

    unsigned char wcolor[3] = {255, 255, 255};
    struct ged_draw_appearance_settings vs = GED_DRAW_APPEARANCE_SETTINGS_INIT;
    std::map<struct directory *, fastf_t> s_size;
    struct draw_data_t dd = {};
    dd.gedp = gedp;
    dd.dbip = gedp->dbip;
    dd.view_ctx = view_ctx;
    dd.tol = &wdbp->wdb_tol;
    dd.ttol = &wdbp->wdb_ttol;
    dd.color_inherit = 0;
    dd.bound_only = 0;
    dd.s_size = &s_size;
    bu_color_from_rgb_chars(&dd.c, wcolor);
    dd.vs = &vs;
    dd.draw_parent_group_ref = draw_parent_group_ref;

    draw_gather_paths(fp, &mat, (void *)&dd);

    db_free_full_path(fp);
    BU_PUT(fp, struct db_full_path);
    return 1;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
