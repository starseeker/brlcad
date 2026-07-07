/*                        O B J S . C P P
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
/** @file libged/view/objs.c
 *
 * Commands for view features.
 *
 */

#include "common.h"

#include <ctype.h>
#include <cstdlib>
#include <cstring>
#include <set>
#include <string>
#include <vector>

extern "C" {
#include "bu/cmd.h"
#include "bu/color.h"
#include "bu/opt.h"
#include "bu/path.h"
#include "bu/vls.h"
#include "raytrace.h"
#include "ged/draw.h"
}
#include "./ged_view.h"
#include "../ged_private.h"

static struct ged_draw_view_record_query
_view_obj_query(int list_view, int list_db, int local_only, const char *glob)
{
    struct ged_draw_view_record_query query;
    query.flags = 0;
    query.glob = glob;
    query.draw_mode = -1;
    if (list_view)
	query.flags |= GED_DRAW_VIEW_RECORD_QUERY_VIEW_OBJECTS;
    if (list_db)
	query.flags |= GED_DRAW_VIEW_RECORD_QUERY_DB_OBJECTS;
    if (local_only)
	query.flags |= GED_DRAW_VIEW_RECORD_QUERY_LOCAL_ONLY;
    return query;
}

struct view_obj_list_state {
    std::set<std::string> names;
};

static int
_view_obj_list_cb(const struct ged_draw_view_db_object_record *rec, void *ud)
{
    struct view_obj_list_state *ctx = (struct view_obj_list_state *)ud;
    if (!ctx || !rec || !rec->path || !rec->path[0])
	return 1;
    ctx->names.insert(std::string(rec->path));
    return 1;
}

static void
_view_obj_list(struct bu_vls *out, void *view_ctx, int list_view, int list_db, int local_only, const char *glob)
{
    if (!out || !view_ctx || !ged_draw_view_context_scene_attached(view_ctx))
	return;

    struct ged_draw_view_record_query query =
	_view_obj_query(list_view, list_db, local_only, glob);
    struct view_obj_list_state ctx;
    ged_draw_foreach_view_record_query(view_ctx, &query, _view_obj_list_cb, &ctx);

    for (std::set<std::string>::iterator it = ctx.names.begin(); it != ctx.names.end(); ++it)
	bu_vls_printf(out, "%s\n", it->c_str());
}

static void
_view_obj_mode_value_string(struct bu_vls *out, int mode)
{
    if (!out) {
	bu_vls_printf(out, "unknown");
	return;
    }
    switch (mode) {
	case GED_DRAW_MODE_WIRE:
	    bu_vls_printf(out, "wireframe");
	    break;
	case GED_DRAW_MODE_SHADED_BOTS:
	case GED_DRAW_MODE_SHADED:
	    bu_vls_printf(out, "shaded");
	    break;
	case GED_DRAW_MODE_EVAL_WIRE:
	    bu_vls_printf(out, "evaluated_wireframe");
	    break;
	case GED_DRAW_MODE_HIDDEN_LINE:
	    bu_vls_printf(out, "hidden_line");
	    break;
	case GED_DRAW_MODE_EVAL_POINTS:
	    bu_vls_printf(out, "evaluated_points");
	    break;
	default:
	    bu_vls_printf(out, "unknown");
	    break;
    }
}

struct view_obj_record {
    int found;
    std::string type_name;
    int draw_mode;
    int visible;
    unsigned char color[3];
    size_t vlist_structure_count;
};

struct view_obj_find_state {
    const char *name;
    struct view_obj_record *out;
};

static int
_view_obj_find_cb(const struct ged_draw_view_db_object_record *rec, void *ud)
{
    struct view_obj_find_state *ctx = (struct view_obj_find_state *)ud;
    if (!ctx || !ctx->out || !rec || !rec->path || !ctx->name)
	return 1;
    if (!BU_STR_EQUAL(rec->path, ctx->name))
	return 1;

    ctx->out->found = 1;
    ctx->out->type_name = rec->type_name ? rec->type_name : "object";
    ctx->out->draw_mode = rec->draw_mode;
    ctx->out->visible = rec->visible;
    ctx->out->color[0] = rec->color[0];
    ctx->out->color[1] = rec->color[1];
    ctx->out->color[2] = rec->color[2];
    ctx->out->vlist_structure_count = rec->vlist_structure_count;
    return 0;
}

static int
_view_obj_record_find(void *view_ctx,
		      const char *name,
		      int list_view,
		      int list_db,
		      int local_only,
		      struct view_obj_record *out)
{
    if (out) {
	out->found = 0;
	out->type_name.clear();
	out->draw_mode = 0;
	out->visible = 0;
	out->color[0] = 0;
	out->color[1] = 0;
	out->color[2] = 0;
	out->vlist_structure_count = 0;
    }
    if (!view_ctx || !ged_draw_view_context_scene_attached(view_ctx) || !name || !name[0])
	return 0;

    struct ged_draw_view_record_query query =
	_view_obj_query(list_view, list_db, local_only, NULL);
    struct view_obj_find_state ctx;
    ctx.name = name;
    ctx.out = out;
    ged_draw_foreach_view_record_query(view_ctx, &query, _view_obj_find_cb, &ctx);

    return out ? out->found : 0;
}

struct view_obj_shape_ref_find_state {
    struct ged *gedp;
    const char *name;
    ged_draw_shape_ref ref;
};

static int
_view_obj_name_matches_record(const struct ged_draw_shape_record *rec,
			      const char *name)
{
    if (!rec || !name || !name[0])
	return 0;
    if (rec->display_name && BU_STR_EQUAL(rec->display_name, name))
	return 1;
    if (rec->leaf_name && BU_STR_EQUAL(rec->leaf_name, name))
	return 1;
    if (rec->fullpath) {
	char *path = db_path_to_string(rec->fullpath);
	if (path) {
	    int match = BU_STR_EQUAL(path, name) ||
		BU_STR_EQUAL(ged_draw_dbpath_skip_lead_slash(path), name);
	    bu_free(path, "view feature fullpath string");
	    if (match)
		return 1;
	}
    }
    return 0;
}

static int
_view_obj_shape_ref_find_cb(const struct ged_draw_shape_record *rec, void *ud)
{
    struct view_obj_shape_ref_find_state *ctx =
	(struct view_obj_shape_ref_find_state *)ud;
    if (!ctx || !ged_draw_shape_ref_is_null(ctx->ref))
	return 0;
    if (_view_obj_name_matches_record(rec, ctx->name)) {
	ctx->ref = rec->ref;
	return 0;
    }
    return 1;
}

static ged_draw_shape_ref
_view_obj_shape_ref_find(struct ged *gedp,
			 void *view_ctx,
			 const char *name,
			 int list_view,
			 int list_db,
			 int local_only)
{
    ged_draw_shape_ref null_ref = GED_DRAW_SHAPE_REF_NULL;
    struct view_obj_record rec;
    if (!_view_obj_record_find(view_ctx, name, list_view, list_db, local_only, &rec))
	return null_ref;

    struct view_obj_shape_ref_find_state ctx;
    ctx.gedp = gedp;
    ctx.name = name;
    ctx.ref = null_ref;
    ged_draw_foreach_shape_record(gedp, _view_obj_shape_ref_find_cb, &ctx);
    return ctx.ref;
}

int
_objs_cmd_draw(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view feature style set <name> visible [0|1|UP|DOWN]";
    const char *purpose_string = "toggle view feature visibility";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (argc == 0) {
	int visible = 0;
	if (!ged_draw_view_context_managed_feature_visible_get(gedp, gd->cv,
		gd->vobj, gd->shape_ref, &visible, gedp->ged_result_str))
	    return BRLCAD_ERROR;
	bu_vls_printf(gedp->ged_result_str, "%s\n", visible ? "UP" : "DOWN");
	return BRLCAD_OK;
    }

    if (BU_STR_EQUAL(argv[0], "DOWN")) {
	return ged_draw_view_context_managed_feature_visible_set(gedp, gd->cv,
		gd->vobj, gd->shape_ref, 0, gedp->ged_result_str) ?
	    BRLCAD_OK : BRLCAD_ERROR;
    }
    if (BU_STR_EQUAL(argv[0], "UP")) {
	return ged_draw_view_context_managed_feature_visible_set(gedp, gd->cv,
		gd->vobj, gd->shape_ref, 1, gedp->ged_result_str) ?
	    BRLCAD_OK : BRLCAD_ERROR;
    }

    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
    return BRLCAD_ERROR;
}

int
_objs_cmd_delete(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view feature delete <name>";
    const char *purpose_string = "delete view feature";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    return ged_draw_view_context_managed_feature_remove(gedp, gd->cv, gd->vobj,
	    gd->shape_ref, gedp->ged_result_str) ? BRLCAD_OK : BRLCAD_ERROR;
}

int
_objs_cmd_color(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view feature style <get|set> <name> color [r/g/b]";
    const char *purpose_string = "show/set feature color";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    int recurse = 0;

    struct bu_opt_desc d[2];
    BU_OPT(d[0], "r", "recursive",       "",  NULL,  &recurse,  "Report (or set) color of all child objects");
    BU_OPT_NULL(d[1]);

    int ac = bu_opt_parse(NULL, argc, argv, d);

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct ged_draw_view_feature_style style =
	GED_DRAW_VIEW_FEATURE_STYLE_INIT;

    if (ac == 0) {
	if (!ged_draw_view_context_managed_feature_style_get(gedp, gd->cv, gd->vobj,
		gd->shape_ref, &style, gedp->ged_result_str))
	    return BRLCAD_ERROR;
	bu_vls_printf(gedp->ged_result_str, "%d/%d/%d\n",
		style.color[0], style.color[1], style.color[2]);
	return BRLCAD_OK;
    }

    struct bu_color val;
    if (bu_opt_color(NULL, 1, (const char **)&argv[0], (void *)&val) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    bu_color_to_rgb_chars(&val, style.color);
    style.color_valid = 1;
    return ged_draw_view_context_managed_feature_style_apply(gedp, gd->cv, gd->vobj,
	    gd->shape_ref, &style, recurse, gedp->ged_result_str) ?
	BRLCAD_OK : BRLCAD_ERROR;
}

int
_objs_cmd_arrow(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view feature style <get|set> <name> arrow [0|1] [width [#]] [length [#]]";
    const char *purpose_string = "toggle arrow drawing, for features that support it";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct ged_draw_view_feature_style style =
	GED_DRAW_VIEW_FEATURE_STYLE_INIT;
    if (!ged_draw_view_context_managed_feature_style_get(gedp, gd->cv, gd->vobj,
	    gd->shape_ref, &style, gedp->ged_result_str))
	return BRLCAD_ERROR;

    if (style.arrow < 0 && style.arrow_tip_width < 0.0 &&
	    style.arrow_tip_length < 0.0) {
	bu_vls_printf(gedp->ged_result_str,
		"View feature %s does not support arrow settings\n", gd->vobj);
	return BRLCAD_ERROR;
    }

    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", style.arrow > 0 ? 1 : 0);
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(argv[0], "0") || BU_STR_EQUAL(argv[0], "1")) {
	style.arrow = BU_STR_EQUAL(argv[0], "1") ? 1 : 0;
	return ged_draw_view_context_managed_feature_style_apply(gedp, gd->cv,
		gd->vobj, gd->shape_ref, &style, 0,
		gedp->ged_result_str) ? BRLCAD_OK : BRLCAD_ERROR;
    }
    if (BU_STR_EQUAL(argv[0], "width"))  {
	if (argc == 2) {
	    fastf_t width = 0.0;
	    if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[1], (void *)&width) != 1) {
		bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
		return BRLCAD_ERROR;
	    }
	    style.arrow_tip_width = width;
	    return ged_draw_view_context_managed_feature_style_apply(gedp, gd->cv,
		    gd->vobj, gd->shape_ref, &style, 0,
		    gedp->ged_result_str) ? BRLCAD_OK : BRLCAD_ERROR;
	}
	bu_vls_printf(gedp->ged_result_str, "%f\n", style.arrow_tip_width);
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(argv[0], "length"))  {
	if (argc == 2) {
	    fastf_t length = 0.0;
	    if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[1], (void *)&length) != 1) {
		bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
		return BRLCAD_ERROR;
	    }
	    style.arrow_tip_length = length;
	    return ged_draw_view_context_managed_feature_style_apply(gedp, gd->cv,
		    gd->vobj, gd->shape_ref, &style, 0,
		    gedp->ged_result_str) ? BRLCAD_OK : BRLCAD_ERROR;
	}
	bu_vls_printf(gedp->ged_result_str, "%f\n", style.arrow_tip_length);
	return BRLCAD_OK;
    }

    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
    return BRLCAD_ERROR;
}

int
_objs_cmd_lcnt(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view feature info <name> lcnt";
    const char *purpose_string = "print the number of vlist entities";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct view_obj_record rec;
    if (_view_obj_record_find(gd->cv, gd->vobj, 1, 1, gd->local_obj, &rec)) {
	bu_vls_printf(gedp->ged_result_str, "%zu\n", rec.vlist_structure_count);
	return BRLCAD_OK;
    }
    bu_vls_printf(gedp->ged_result_str, "0\n");
    return BRLCAD_OK;
}

int
_objs_cmd_update(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view feature realize <name> [x y]";
    const char *purpose_string = "update feature";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (argc && (argc != 2)) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    int have_xy = 0;
    int x = 0;
    int y = 0;
    if (argc) {
	if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&x) != 1 || x < 0) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	    return BRLCAD_ERROR;
	}
	if (bu_opt_int(NULL, 1, (const char **)&argv[1], (void *)&y) != 1 || y < 0) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[1]);
	    return BRLCAD_ERROR;
	}
	have_xy = 1;
	ged_view_context_mouse_state_set(gd->cv, x, y);
    }

    ged_draw_view_polygon_ref poly_ref =
	ged_draw_view_context_polygon_find_scoped(gd->cv, gd->vobj,
		gd->local_obj);
    if (ged_draw_view_polygon_ref_is_null(poly_ref) && !gd->local_obj) {
	poly_ref = ged_draw_view_context_polygon_find_scoped(gd->cv,
		gd->vobj, 1);
	if (!ged_draw_view_polygon_ref_is_null(poly_ref))
	    gd->local_obj = 1;
    }
    if (!ged_draw_view_polygon_ref_is_null(poly_ref)) {
	if (have_xy) {
	    if (!ged_draw_view_context_polygon_update_screen_pt(poly_ref,
		    gd->cv, x, y, GED_DRAW_VIEW_POLYGON_UPDATE_DEFAULT))
		return BRLCAD_ERROR;
	} else if (!ged_draw_view_context_polygon_update(poly_ref, gd->cv,
		GED_DRAW_VIEW_POLYGON_UPDATE_DEFAULT)) {
	    return BRLCAD_ERROR;
	}
	return BRLCAD_OK;
    }

    return ged_draw_view_context_managed_feature_realize(gedp, gd->cv, gd->vobj,
	    gd->shape_ref, gedp->ged_result_str) ? BRLCAD_OK : BRLCAD_ERROR;
}

static int
_view_object_scene_ready(struct _ged_view_info *gd)
{
    if (!gd || !gd->gedp)
	return 0;
    struct ged *gedp = gd->gedp;
    if (!gd->cv) {
	bu_vls_printf(gedp->ged_result_str, ": no view current in GED");
	return 0;
    }
    if (!ged_draw_view_context_scene_attached(gd->cv)) {
	void *cv = ged_view_active_ctx(gedp);
	ged_view_active_ctx_set(gedp, gd->cv);
	ged_draw_ensure_root_attached(gedp);
	ged_view_active_ctx_set(gedp, cv);
    }
    return 1;
}

static int
_view_object_first_pos(int argc, const char **argv)
{
    for (int i = 0; i < argc; i++) {
	if (!argv[i])
	    return -1;
	if (argv[i][0] == '-')
	    continue;
	return i;
    }
    return -1;
}

static void
_view_object_list_defaults(int annotations_requested, int db_requested,
	int all_requested, int *list_view, int *list_db)
{
    if (all_requested || (!annotations_requested && !db_requested)) {
	*list_view = 1;
	*list_db = 1;
	return;
    }
    *list_view = annotations_requested ? 1 : 0;
    *list_db = db_requested ? 1 : 0;
}

static int
_view_object_info(struct _ged_view_info *gd,
	const char *name,
	const char *field,
	int list_view,
	int list_db)
{
    struct ged *gedp = gd->gedp;
    struct view_obj_record rec;
    if (!_view_obj_record_find(gd->cv, name, list_view, list_db, gd->local_obj, &rec)) {
	bu_vls_printf(gedp->ged_result_str, "No view feature named %s\n", name);
	return BRLCAD_ERROR;
    }

    if (!field) {
	bu_vls_printf(gedp->ged_result_str, "%s %s\n", name, rec.type_name.c_str());
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(field, "mode")) {
	_view_obj_mode_value_string(gedp->ged_result_str, rec.draw_mode);
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(field, "color")) {
	bu_vls_printf(gedp->ged_result_str, "%d/%d/%d\n",
		rec.color[0], rec.color[1], rec.color[2]);
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(field, "visible") || BU_STR_EQUAL(field, "draw")) {
	bu_vls_printf(gedp->ged_result_str, "%s\n", rec.visible ? "UP" : "DOWN");
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(field, "lcnt")) {
	bu_vls_printf(gedp->ged_result_str, "%zu\n", rec.vlist_structure_count);
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(field, "type") || BU_STR_EQUAL(field, "class")) {
	bu_vls_printf(gedp->ged_result_str, "%s\n", rec.type_name.c_str());
	return BRLCAD_OK;
    }

    bu_vls_printf(gedp->ged_result_str, "Unsupported object info field %s", field);
    return BRLCAD_ERROR;
}

static int
_view_object_style_get(struct _ged_view_info *gd, const char *field)
{
    if (!field) {
	bu_vls_printf(gd->gedp->ged_result_str,
		"Usage: view feature style get <name> <color|visible|arrow|arrow_width|arrow_length>");
	return BRLCAD_ERROR;
    }
    if (BU_STR_EQUAL(field, "color")) {
	const char *cargv[2] = {"color", NULL};
	return _objs_cmd_color(gd, 1, cargv);
    }
    if (BU_STR_EQUAL(field, "visible") || BU_STR_EQUAL(field, "draw")) {
	const char *dargv[2] = {"visible", NULL};
	return _objs_cmd_draw(gd, 1, dargv);
    }
    if (BU_STR_EQUAL(field, "arrow")) {
	const char *aargv[2] = {"arrow", NULL};
	return _objs_cmd_arrow(gd, 1, aargv);
    }
    if (BU_STR_EQUAL(field, "arrow_width")) {
	const char *aargv[3] = {"arrow", "width", NULL};
	return _objs_cmd_arrow(gd, 2, aargv);
    }
    if (BU_STR_EQUAL(field, "arrow_length")) {
	const char *aargv[3] = {"arrow", "length", NULL};
	return _objs_cmd_arrow(gd, 2, aargv);
    }
    bu_vls_printf(gd->gedp->ged_result_str, "Unsupported style field %s", field);
    return BRLCAD_ERROR;
}

static int
_view_object_style_set(struct _ged_view_info *gd, int argc, const char **argv)
{
    if (argc < 1) {
	bu_vls_printf(gd->gedp->ged_result_str,
		"Usage: view feature style set <name> <field> <value>");
	return BRLCAD_ERROR;
    }
    if (BU_STR_EQUAL(argv[0], "color")) {
	return _objs_cmd_color(gd, argc, argv);
    }
    if (BU_STR_EQUAL(argv[0], "visible") || BU_STR_EQUAL(argv[0], "draw")) {
	if (argc != 2) {
	    bu_vls_printf(gd->gedp->ged_result_str,
		    "Usage: view feature style set <name> visible <0|1|UP|DOWN>");
	    return BRLCAD_ERROR;
	}
	const char *state = argv[1];
	if (BU_STR_EQUAL(state, "1"))
	    state = "UP";
	else if (BU_STR_EQUAL(state, "0"))
	    state = "DOWN";
	const char *dargv[3] = {"visible", state, NULL};
	return _objs_cmd_draw(gd, 2, dargv);
    }
    if (BU_STR_EQUAL(argv[0], "arrow")) {
	return _objs_cmd_arrow(gd, argc, argv);
    }
    if (BU_STR_EQUAL(argv[0], "arrow_width")) {
	if (argc != 2) {
	    bu_vls_printf(gd->gedp->ged_result_str,
		    "Usage: view feature style set <name> arrow_width <value>");
	    return BRLCAD_ERROR;
	}
	const char *aargv[4] = {"arrow", "width", argv[1], NULL};
	return _objs_cmd_arrow(gd, 3, aargv);
    }
    if (BU_STR_EQUAL(argv[0], "arrow_length")) {
	if (argc != 2) {
	    bu_vls_printf(gd->gedp->ged_result_str,
		    "Usage: view feature style set <name> arrow_length <value>");
	    return BRLCAD_ERROR;
	}
	const char *aargv[4] = {"arrow", "length", argv[1], NULL};
	return _objs_cmd_arrow(gd, 3, aargv);
    }
    bu_vls_printf(gd->gedp->ged_result_str, "Unsupported style field %s", argv[0]);
    return BRLCAD_ERROR;
}

extern "C" int
_view_cmd_feature(void *bs, int argc, const char **argv)
{
    int help = 0;
    int local_only = 0;
    int annotations_requested = 0;
    int db_requested = 0;
    int all_requested = 0;
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;

    const char *usage_string = "view [options] feature [options] <list|info|show|hide|delete|style|realize> ...";
    const char *purpose_string = "manage typed Obol view features";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;
    if (!_view_object_scene_ready(gd))
	return BRLCAD_ERROR;

    argc--; argv++;

    struct bu_opt_desc d[6];
    BU_OPT(d[0], "h", "help",        "",  NULL,  &help,                  "Print help");
    BU_OPT(d[1], "L", "local",       "",  NULL,  &local_only,            "Restrict feature lookup to current/specified view");
    BU_OPT(d[2], "A", "annotations", "",  NULL,  &annotations_requested, "List or query annotation features");
    BU_OPT(d[3], "D", "database",    "",  NULL,  &db_requested,          "List or query database display features");
    BU_OPT(d[4], "a", "all",         "",  NULL,  &all_requested,         "List or query all feature classes");
    BU_OPT_NULL(d[5]);
    gd->gopts = d;

    int first_pos = _view_object_first_pos(argc, argv);
    int acnt = (first_pos >= 0) ? first_pos : argc;
    (void)bu_opt_parse(NULL, acnt, argv, d);
    gd->local_obj = local_only;

    int list_view = 0;
    int list_db = 0;
    _view_object_list_defaults(annotations_requested, db_requested, all_requested, &list_view, &list_db);

    if (help) {
	_ged_cmd_help(gedp, usage_string, d);
	return GED_HELP;
    }

    if (first_pos < 0) {
	_view_obj_list(gedp->ged_result_str, gd->cv, list_view, list_db, gd->local_obj, NULL);
	return BRLCAD_OK;
    }

    const char *cmd = argv[first_pos];
    int sub_argc = argc - first_pos;
    const char **sub_argv = argv + first_pos;

    if (BU_STR_EQUAL(cmd, "list")) {
	if (sub_argc > 2) {
	    bu_vls_printf(gedp->ged_result_str,
		    "Usage: view feature [--all|--annotations|--database] [-L] list [glob_pattern]");
	    return BRLCAD_ERROR;
	}
	const char *glob = (sub_argc > 1) ? sub_argv[1] : NULL;
	_view_obj_list(gedp->ged_result_str, gd->cv, list_view, list_db, gd->local_obj, glob);
	return BRLCAD_OK;
    }

    if (BU_STR_EQUAL(cmd, "info")) {
	if (sub_argc < 2 || sub_argc > 3) {
	    bu_vls_printf(gedp->ged_result_str, "Usage: view feature info <name> [field]");
	    return BRLCAD_ERROR;
	}
	return _view_object_info(gd, sub_argv[1], (sub_argc == 3) ? sub_argv[2] : NULL,
		list_view, list_db);
    }

    if (BU_STR_EQUAL(cmd, "show") || BU_STR_EQUAL(cmd, "hide") ||
	    BU_STR_EQUAL(cmd, "delete") || BU_STR_EQUAL(cmd, "realize")) {
	if (sub_argc < 2) {
	    bu_vls_printf(gedp->ged_result_str, "Usage: view feature %s <name>", cmd);
	    return BRLCAD_ERROR;
	}
	gd->vobj = sub_argv[1];
	gd->shape_ref = _view_obj_shape_ref_find(gedp, gd->cv, gd->vobj, list_view, list_db, gd->local_obj);
	if (BU_STR_EQUAL(cmd, "show")) {
	    const char *dargv[3] = {"visible", "UP", NULL};
	    return _objs_cmd_draw(gd, 2, dargv);
	}
	if (BU_STR_EQUAL(cmd, "hide")) {
	    const char *dargv[3] = {"visible", "DOWN", NULL};
	    return _objs_cmd_draw(gd, 2, dargv);
	}
	if (BU_STR_EQUAL(cmd, "delete")) {
	    if (sub_argc != 2) {
		bu_vls_printf(gedp->ged_result_str, "Usage: view feature delete <name>");
		return BRLCAD_ERROR;
	    }
	    const char *rargv[2] = {"delete", NULL};
	    return _objs_cmd_delete(gd, 1, rargv);
	}
	std::vector<const char *> uargv;
	uargv.push_back("realize");
	for (int i = 2; i < sub_argc; i++)
	    uargv.push_back(sub_argv[i]);
	uargv.push_back(NULL);
	return _objs_cmd_update(gd, (int)uargv.size() - 1, uargv.data());
    }

    if (BU_STR_EQUAL(cmd, "style")) {
	if (sub_argc < 4) {
	    bu_vls_printf(gedp->ged_result_str,
		    "Usage: view feature style <get|set> <name> <field> [value...]");
	    return BRLCAD_ERROR;
	}
	const char *style_cmd = sub_argv[1];
	gd->vobj = sub_argv[2];
	gd->shape_ref = _view_obj_shape_ref_find(gedp, gd->cv, gd->vobj, list_view, list_db, gd->local_obj);
	if (BU_STR_EQUAL(style_cmd, "get")) {
	    if (sub_argc != 4) {
		bu_vls_printf(gedp->ged_result_str,
			"Usage: view feature style get <name> <field>");
		return BRLCAD_ERROR;
	    }
	    return _view_object_style_get(gd, sub_argv[3]);
	}
	if (BU_STR_EQUAL(style_cmd, "set")) {
	    return _view_object_style_set(gd, sub_argc - 3, sub_argv + 3);
	}
	bu_vls_printf(gedp->ged_result_str, "Unsupported style operation %s", style_cmd);
	return BRLCAD_ERROR;
    }

    bu_vls_printf(gedp->ged_result_str,
	    "Unsupported feature subcommand %s (valid: list, info, show, hide, delete, style, realize)",
	    cmd);
    return BRLCAD_ERROR;
}

extern "C" int
_view_cmd_annotation(void *bs, int argc, const char **argv)
{
    int help = 0;
    int local_only = 0;
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;

    const char *usage_string = "view [options] annotation [-L] <line|arrow|label|axes> <verb> <name> [args]";
    const char *purpose_string = "create and edit user annotation objects";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;
    if (!_view_object_scene_ready(gd))
	return BRLCAD_ERROR;

    argc--; argv++;

    struct bu_opt_desc d[3];
    BU_OPT(d[0], "h", "help",  "",  NULL,  &help,       "Print help");
    BU_OPT(d[1], "L", "local", "",  NULL,  &local_only, "Feature is scoped only to current/specified view");
    BU_OPT_NULL(d[2]);
    gd->gopts = d;

    int first_pos = _view_object_first_pos(argc, argv);
    int acnt = (first_pos >= 0) ? first_pos : argc;
    (void)bu_opt_parse(NULL, acnt, argv, d);
    gd->local_obj = local_only;

    if (help) {
	_ged_cmd_help(gedp, usage_string, d);
	return GED_HELP;
    }
    if (first_pos < 0 || argc - first_pos < 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    const char **sub_argv = argv + first_pos;
    int sub_argc = argc - first_pos;
    const char *type = sub_argv[0];
    const char *verb = sub_argv[1];
    gd->vobj = sub_argv[2];
    gd->shape_ref = _view_obj_shape_ref_find(gedp, gd->cv, gd->vobj, 1, 1, gd->local_obj);

    std::vector<const char *> cargv;
    if (BU_STR_EQUAL(type, "arrow")) {
	if (!BU_STR_EQUAL(verb, "create")) {
	    bu_vls_printf(gedp->ged_result_str,
		    "Unsupported arrow annotation operation %s", verb);
	    return BRLCAD_ERROR;
	}
	cargv.push_back("line");
	cargv.push_back("create");
	for (int i = 3; i < sub_argc; i++)
	    cargv.push_back(sub_argv[i]);
	cargv.push_back(NULL);
	int ret = _view_cmd_lines(gd, (int)cargv.size() - 1, cargv.data());
	if (ret != BRLCAD_OK)
	    return ret;
	gd->shape_ref = _view_obj_shape_ref_find(gedp, gd->cv, gd->vobj, 1, 1, gd->local_obj);
	const char *aargv[3] = {"arrow", "1", NULL};
	return _objs_cmd_arrow(gd, 2, aargv);
    }

    cargv.push_back(type);
    cargv.push_back(verb);
    for (int i = 3; i < sub_argc; i++)
	cargv.push_back(sub_argv[i]);
    cargv.push_back(NULL);

    if (BU_STR_EQUAL(type, "line"))
	return _view_cmd_lines(gd, (int)cargv.size() - 1, cargv.data());
    if (BU_STR_EQUAL(type, "label"))
	return _view_cmd_labels(gd, (int)cargv.size() - 1, cargv.data());
    if (BU_STR_EQUAL(type, "axes"))
	return _view_cmd_axes(gd, (int)cargv.size() - 1, cargv.data());

    bu_vls_printf(gedp->ged_result_str,
	    "Unsupported annotation type %s (valid: line, arrow, label, axes)",
	    type);
    return BRLCAD_ERROR;
}

extern "C" int
_view_cmd_polygon(void *bs, int argc, const char **argv)
{
    int help = 0;
    int local_only = 0;
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;

    const char *usage_string = "view [options] polygon [-L] <verb> <name> [args]";
    const char *purpose_string = "create and edit polygon annotation objects";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;
    if (!_view_object_scene_ready(gd))
	return BRLCAD_ERROR;

    argc--; argv++;

    struct bu_opt_desc d[3];
    BU_OPT(d[0], "h", "help",  "",  NULL,  &help,       "Print help");
    BU_OPT(d[1], "L", "local", "",  NULL,  &local_only, "Object is scoped only to current/specified view");
    BU_OPT_NULL(d[2]);
    gd->gopts = d;

    int first_pos = _view_object_first_pos(argc, argv);
    int acnt = (first_pos >= 0) ? first_pos : argc;
    (void)bu_opt_parse(NULL, acnt, argv, d);
    gd->local_obj = local_only;

    if (help) {
	_ged_cmd_help(gedp, usage_string, d);
	return GED_HELP;
    }
    if (first_pos < 0 || argc - first_pos < 2) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    const char **sub_argv = argv + first_pos;
    int sub_argc = argc - first_pos;
    const char *verb = sub_argv[0];
    gd->vobj = sub_argv[1];
    gd->polygon_ref = ged_draw_view_context_polygon_find_scoped(gd->cv,
	    gd->vobj, gd->local_obj);
    gd->shape_ref = _view_obj_shape_ref_find(gedp, gd->cv, gd->vobj, 1, 1, gd->local_obj);

    if (BU_STR_EQUAL(verb, "update")) {
	std::vector<const char *> uargv;
	uargv.push_back("update");
	for (int i = 2; i < sub_argc; i++)
	    uargv.push_back(sub_argv[i]);
	uargv.push_back(NULL);
	return _objs_cmd_update(gd, (int)uargv.size() - 1, uargv.data());
    }

    std::vector<const char *> cargv;
    cargv.push_back("polygon");
    cargv.push_back(verb);
    for (int i = 2; i < sub_argc; i++)
	cargv.push_back(sub_argv[i]);
    cargv.push_back(NULL);
    return _view_cmd_polygons(gd, (int)cargv.size() - 1, cargv.data());
}

extern "C" int
_view_cmd_db_objects(void *bs, int argc, const char **argv)
{
    int help = 0;
    struct bu_vls as_name = BU_VLS_INIT_ZERO;
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;

    const char *usage_string = "view [options] db <add|delete|list> ...";
    const char *purpose_string = "manage Obol database display objects";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;
    if (!_view_object_scene_ready(gd))
	return BRLCAD_ERROR;

    argc--; argv++;

    if (argc <= 0) {
	_view_obj_list(gedp->ged_result_str, gd->cv, 0, 1, 0, NULL);
	return BRLCAD_OK;
    }

    if (BU_STR_EQUAL(argv[0], "list")) {
	const char *glob = (argc > 1) ? argv[1] : NULL;
	if (argc > 2) {
	    bu_vls_printf(gedp->ged_result_str, "Usage: view db list [glob_pattern]");
	    return BRLCAD_ERROR;
	}
	_view_obj_list(gedp->ged_result_str, gd->cv, 0, 1, 0, glob);
	return BRLCAD_OK;
    }

    if (BU_STR_EQUAL(argv[0], "add")) {
	if (argc < 2) {
	    bu_vls_printf(gedp->ged_result_str, "Usage: view db add <dbpath> [--as <name>]");
	    return BRLCAD_ERROR;
	}
	struct bu_opt_desc d[3];
	BU_OPT(d[0], "h", "help",  "",      NULL,        &help,       "Print help");
	BU_OPT(d[1], "",  "as",    "name",  &bu_opt_vls, &as_name,    "View feature name");
	BU_OPT_NULL(d[2]);
	gd->gopts = d;
	const char **add_argv = argv + 1;
	int acnt = bu_opt_parse(NULL, argc - 1, add_argv, d);
	if (help) {
	    _ged_cmd_help(gedp, "view db add <dbpath> [--as <name>]", d);
	    bu_vls_free(&as_name);
	    return GED_HELP;
	}
	if (acnt < 1) {
	    bu_vls_free(&as_name);
	    bu_vls_printf(gedp->ged_result_str, "Usage: view db add <dbpath> [--as <name>]");
	    return BRLCAD_ERROR;
	}
	const char *dbpath = add_argv[0];
	const char *name = bu_vls_strlen(&as_name) ? bu_vls_cstr(&as_name) : dbpath;
	int ret = ged_draw_view_context_gobject_create(gedp, gd->cv, dbpath,
		name, gedp->ged_result_str) ? BRLCAD_OK : BRLCAD_ERROR;
	bu_vls_free(&as_name);
	return ret;
    }

    if (BU_STR_EQUAL(argv[0], "delete")) {
	if (argc != 2) {
	    bu_vls_printf(gedp->ged_result_str, "Usage: view db delete <name>");
	    return BRLCAD_ERROR;
	}
	gd->vobj = argv[1];
	gd->shape_ref = _view_obj_shape_ref_find(gedp, gd->cv, gd->vobj, 0, 1, 0);
	const char *rargv[2] = {"delete", NULL};
	return _objs_cmd_delete(gd, 1, rargv);
    }

    bu_vls_printf(gedp->ged_result_str,
	    "Unsupported db subcommand %s (valid: add, delete, list)", argv[0]);
    return BRLCAD_ERROR;
}


// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
