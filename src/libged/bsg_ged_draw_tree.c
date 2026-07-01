/*                B S G _ G E D _ D R A W _ T R E E . C
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
/** @file libged/b_s_g___g_e_d___d_r_a_w___t_r_e_e_._c.c
 *
 * Scene-ref retained draw-tree mutation and group path management.
 */

#include "common.h"

#include <stdlib.h>
#include <string.h>

#include "bu/color.h"
#include "bu/hash.h"
#include "bu/malloc.h"
#include "bu/str.h"
#include "bg/clip.h"

#include "ged.h"
#include "ged/draw.h"
#include "./bsg_ged_draw_private.h"
#include "./ged_private.h"


static int _sg_erase_path(struct ged *gedp, const char *path);
static int _sg_erase_all_paths(struct ged *gedp, const char *path);
static int _sg_erase_path_scoped(struct ged *gedp, const char *path,
				 void *view_ctx, int mode);
static int _sg_erase_all_paths_scoped(struct ged *gedp, const char *path,
				      void *view_ctx, int mode);
static int _sg_erase_nonroot_component_scoped(struct ged *gedp,
					      const char *name,
					      void *view_ctx,
					      int mode);
static int _sg_erase_component_scoped(struct ged *gedp,
				      const char *name,
				      void *view_ctx,
				      int mode);


static int
_sg_erase_path(struct ged *gedp, const char *path)
{
    return ged_draw_scene_root_erase_path(gedp, path);
}


static int
_sg_erase_path_scoped(struct ged *gedp, const char *path,
		      void *view_ctx, int mode)
{
    return ged_draw_source_erase_path_in_active_scope(gedp, path, view_ctx,
	    mode);
}


static void
_sg_erase_all_names(struct ged *gedp, const char *name)
{
    ged_draw_overlay_erase_name(gedp, name);
    (void)ged_draw_scene_root_erase_groups_by_name(gedp, name);
}


static int
_sg_erase_all_paths(struct ged *gedp, const char *path)
{
    return ged_draw_scene_root_erase_path_prefix(gedp, path);
}


static int
_sg_erase_all_paths_scoped(struct ged *gedp, const char *path,
			   void *view_ctx, int mode)
{
    return ged_draw_source_erase_path_prefix_in_active_scope(gedp, path,
	    view_ctx, mode);
}


static int
_sg_erase_component_scoped_impl(struct ged *gedp,
				const char *name,
				void *view_ctx,
				int mode,
				int nonroot_only)
{
    if (!nonroot_only)
	ged_draw_overlay_erase_name(gedp, name);

    return ged_draw_source_erase_component_name_in_active_scope(gedp,
	    name, view_ctx, mode, nonroot_only);
}


static int
_sg_erase_nonroot_component_scoped(struct ged *gedp,
				   const char *name,
				   void *view_ctx,
				   int mode)
{
    return _sg_erase_component_scoped_impl(gedp, name, view_ctx, mode, 1);
}


static int
_sg_erase_component_scoped(struct ged *gedp,
			   const char *name,
			   void *view_ctx,
			   int mode)
{
    return _sg_erase_component_scoped_impl(gedp, name, view_ctx, mode, 0);
}


static int
_sg_bounding_sph(struct ged *gedp, vect_t *min, vect_t *max, int pflag)
{
    return ged_draw_scene_root_subtree_bounds(gedp, min, max, pflag);
}


void
ged_draw_erase_name(struct ged *gedp, const char *name)
{
    if (!gedp || !name)
	return;
    _sg_erase_all_names(gedp, name);
}


int
ged_draw_erase_path_string(struct ged *gedp, const char *path)
{
    if (!gedp || !path)
	return 0;
    uint64_t rev0 = gedp->i->ged_gdp->gd_draw_rev;
    int erased = _sg_erase_path(gedp, ged_draw_dbpath_skip_lead_slash(path));
    return erased ? erased :
	((gedp->i->ged_gdp->gd_draw_rev != rev0) ? 1 : 0);
}


int
ged_draw_erase_path_prefix_string(struct ged *gedp, const char *path)
{
    if (!gedp || !path)
	return 0;
    uint64_t rev0 = gedp->i->ged_gdp->gd_draw_rev;
    int erased = _sg_erase_all_paths(gedp,
	    ged_draw_dbpath_skip_lead_slash(path));
    return erased ? erased :
	((gedp->i->ged_gdp->gd_draw_rev != rev0) ? 1 : 0);
}


int
ged_draw_erase_path_string_scoped(struct ged *gedp,
				  const char *path,
				  void *view_ctx,
				  int mode)
{
    if (!gedp || !path)
	return 0;
    uint64_t rev0 = gedp->i->ged_gdp->gd_draw_rev;
    int erased = _sg_erase_path_scoped(gedp,
	    ged_draw_dbpath_skip_lead_slash(path),
	    view_ctx, mode);
    return erased ? erased :
	((gedp->i->ged_gdp->gd_draw_rev != rev0) ? 1 : 0);
}


int
ged_draw_erase_path_prefix_string_scoped(struct ged *gedp,
					 const char *path,
					 void *view_ctx,
					 int mode)
{
    if (!gedp || !path)
	return 0;
    uint64_t rev0 = gedp->i->ged_gdp->gd_draw_rev;
    int erased = _sg_erase_all_paths_scoped(gedp,
	    ged_draw_dbpath_skip_lead_slash(path),
	    view_ctx, mode);
    return erased ? erased :
	((gedp->i->ged_gdp->gd_draw_rev != rev0) ? 1 : 0);
}


int
ged_draw_erase_nonroot_component_string_scoped(struct ged *gedp,
					       const char *name,
					       void *view_ctx,
					       int mode)
{
    if (!gedp || !name)
	return 0;
    return _sg_erase_nonroot_component_scoped(gedp, name, view_ctx, mode);
}


int
ged_draw_erase_component_string_scoped(struct ged *gedp,
				       const char *name,
				       void *view_ctx,
				       int mode)
{
    if (!gedp || !name)
	return 0;
    return _sg_erase_component_scoped(gedp, name, view_ctx, mode);
}


int
ged_draw_apply_erase_path(struct ged *gedp,
			  const struct db_full_path *dfp)
{
    if (!gedp || !dfp)
	return 0;

    char *s = db_path_to_string(dfp);
    if (!s)
	return 0;

    struct ged_draw_transaction txn =
	ged_draw_transaction_make(GED_DRAW_TXN_ERASE,
				  ged_draw_dbpath_skip_lead_slash(s));
    int ret = ged_draw_apply_transaction(gedp, &txn, NULL);
    bu_free(s, "ged_draw_apply_erase_path: path string");
    return ret;
}


int
ged_draw_apply_erase_path_prefix(struct ged *gedp,
				 const struct db_full_path *dfp)
{
    if (!gedp || !dfp)
	return 0;

    char *s = db_path_to_string(dfp);
    if (!s)
	return 0;

    struct ged_draw_transaction txn =
	ged_draw_transaction_make(GED_DRAW_TXN_ERASE_PREFIX,
				  ged_draw_dbpath_skip_lead_slash(s));
    int ret = ged_draw_apply_transaction(gedp, &txn, NULL);
    bu_free(s, "ged_draw_apply_erase_path_prefix: path string");
    return ret;
}


void
ged_draw_erase_path(struct ged *gedp,
		    const struct db_full_path *dfp)
{
    if (!gedp || !dfp)
	return;

    char *s = db_path_to_string(dfp);
    if (!s)
	return;

    struct ged_draw_transaction txn =
	ged_draw_transaction_make(GED_DRAW_TXN_ERASE,
				  ged_draw_dbpath_skip_lead_slash(s));
    ged_draw_apply_transaction(gedp, &txn, NULL);
    bu_free(s, "ged_draw_erase_path: path string");
}


void
ged_draw_erase_path_prefix(struct ged *gedp,
			   const struct db_full_path *dfp)
{
    if (!gedp || !dfp)
	return;

    char *s = db_path_to_string(dfp);
    if (!s)
	return;

    struct ged_draw_transaction txn =
	ged_draw_transaction_make(GED_DRAW_TXN_ERASE_PREFIX,
				  ged_draw_dbpath_skip_lead_slash(s));
    ged_draw_apply_transaction(gedp, &txn, NULL);
    bu_free(s, "ged_draw_erase_path_prefix: path string");
}


int
ged_draw_bounds(struct ged *gedp, vect_t *min, vect_t *max, int pflag)
{
    if (!gedp || !min || !max)
	return 1;
    return _sg_bounding_sph(gedp, min, max, pflag);
}


void
ged_draw_clear(struct ged *gedp)
{
    if (!gedp)
	return;

    (void)ged_draw_scene_root_clear_all_scope_children(gedp);
    ged_draw_registry_free(gedp);

    gedp->i->ged_gdp->gd_draw_rev = 0;

    struct ged_drawable *gdp = gedp->i->ged_gdp;
    gdp->gd_highlight_token = 0;
    gdp->gd_highlight_scene_rev = 0;
    gdp->gd_highlight_rev++;
}


int
ged_draw_clear_view(struct ged *gedp, void *view_ctx)
{
    if (!gedp)
	return 0;

    return ged_draw_source_clear_db_groups_in_scope(gedp, view_ctx);
}


int
ged_draw_has_groups(struct ged *gedp)
{
    if (!gedp || !gedp->i || !gedp->i->ged_gdp)
	return 0;

    return ged_draw_scene_root_has_groups(gedp);
}


/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
