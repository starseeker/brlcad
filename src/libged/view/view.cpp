/*                     V I E W . C P P
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
/** @file libged/view.c
 *
 * The view command.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "bu/cmd.h"
#include "bu/malloc.h"
#include "bu/vls.h"

#include "BObol/BExportAction.h"
#include "BObol/BSceneController.h"
#include "BObol/BViewStore.h"
#include <Inventor/SoViewport.h>
#include "bv.h"
#include "ged/draw.h"
#include "../ged_bobol_private.hpp"
#include "../ged_private.h"
#include "./ged_view.h"

/* Visit context + callbacks for the "view vZ" autodetect path. */
struct _view_vZ_ctx {
    mat_t model2view;
    int calc_mode;
    double vZ;
    int have_vz;
};

static const char *
_view_name(const struct ged_view_context *view_ctx)
{
    const char *name = bv_name_get(
			   bv_context_view_const((const struct bv_context *)view_ctx));
    return name ? name : "";
}

static int
_view_vZ_consider(struct _view_vZ_ctx *c, fastf_t calc_val)
{
    if (c->calc_mode) {
	if (calc_val > c->vZ) {
	    c->vZ = calc_val;
	    c->have_vz = 1;
	}
    } else {
	if (calc_val < c->vZ) {
	    c->vZ = calc_val;
	    c->have_vz = 1;
	}
    }
    return 1;
}

static void
_view_vZ_consider_model_point(struct _view_vZ_ctx *c, const point_t p)
{
    vect_t vpt;
    MAT4X3PNT(vpt, c->model2view, p);
    _view_vZ_consider(c, vpt[Z]);
}

static void
_view_vZ_consider_feature_record(struct _view_vZ_ctx *ctx,
    const BObolFeatureRecord &record)
{
    if (!ctx)
	return;
    const auto consider = [ctx](const SbVec3f &point) {
	point_t model_point;
	VSET(model_point, point[0], point[1], point[2]);
	_view_vZ_consider_model_point(ctx, model_point);
    };
    for (const SbVec3f &point : record.points)
	consider(point);
    for (const BObolLabel &label : record.labels) {
	consider(label.point);
	if (label.hasLeader)
	    consider(label.target);
    }
    for (const SbVec3f &center : record.axesCenters)
	consider(center);
}

static int
_view_vZ_feature_record_cb(const BObolFeatureRecord &record, void *data)
{
    _view_vZ_consider_feature_record(
	static_cast<struct _view_vZ_ctx *>(data), record);
    return 1;
}

static int
_view_vZ_feature_depth(struct ged *gedp, struct ged_view_context *view_ctx,
    const char *name, int mode, fastf_t *depth)
{
    if (depth)
	*depth = 0.0;
    const struct ged_bobol_feature_binding binding =
	ged_bobol_feature_find(gedp, view_ctx, name);
    BObolFeatureRecord record;
    if (!depth || !binding.controller || !binding.handle.isValid() ||
	!binding.controller->features().record(binding.handle, record))
	return 0;

    struct _view_vZ_ctx ctx;
    if (!bv_model2view_get(ctx.model2view,
	bv_context_view_const(ged_view_context_bv_const(view_ctx))))
	return 0;
    ctx.calc_mode = mode;
    ctx.vZ = mode ? -DBL_MAX : DBL_MAX;
    ctx.have_vz = 0;
    _view_vZ_consider_feature_record(&ctx, record);
    if (!ctx.have_vz)
	return 0;
    *depth = ctx.vZ;
    return 1;
}

static void
_view_vZ_visit_features(struct ged *gedp,
    struct ged_view_context *view_ctx, struct _view_vZ_ctx *ctx)
{
    BObolViewController *local = ged_bobol_view_controller(view_ctx);
    const BObolFeatureOwner owner = ged_bobol_view_feature_owner(view_ctx);
    if (local) {
	local->features().visitRecords(_view_vZ_feature_record_cb, ctx,
	    BOBOL_FEATURE_SCOPE_LOCAL, &owner);
	local->features().visitRecords(_view_vZ_feature_record_cb, ctx,
	    BOBOL_FEATURE_SCOPE_SHARED, nullptr);
    }
    BObolViewController *shared = ged_bobol_shared_view_controller(gedp);
    if (shared && shared != local)
	shared->features().visitRecords(_view_vZ_feature_record_cb, ctx,
	    BOBOL_FEATURE_SCOPE_SHARED, nullptr);
}

static void
_view_vZ_visit_db_exports(struct ged_view_context *view_ctx, struct _view_vZ_ctx *ctx)
{
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (!ctx || !controller || !controller->getViewport() ||
	!controller->getViewport()->getRoot())
	return;

    SoBRLExportAction export_action;
    export_action.setGeometryPolicy(SoBRLExportAction::DISPLAY_LEVEL);
    export_action.apply(controller->getViewport()->getRoot());

    std::vector<SoBRLExportAction::ObjectRecord> records;
    export_action.collectObjectRecords(records,
	SoBRLExportAction::QUERY_VISIBLE_ONLY |
	SoBRLExportAction::QUERY_DATABASE_OBJECTS);
    const auto consider = [ctx](const SbVec3f &point) {
	point_t model_point;
	VSET(model_point, point[0], point[1], point[2]);
	_view_vZ_consider_model_point(ctx, model_point);
    };
    for (const SoBRLExportAction::ObjectRecord &record : records) {
	for (int index : record.lineIndices) {
	    const SoBRLExportAction::LineRecord &line =
		export_action.getLine(index);
	    consider(line.a);
	    consider(line.b);
	}
	for (int index : record.pointIndices)
	    consider(export_action.getPoint(index).point);
	for (int index : record.triangleIndices) {
	    const SoBRLExportAction::TriangleRecord &triangle =
		export_action.getTriangle(index);
	    consider(triangle.a);
	    consider(triangle.b);
	    consider(triangle.c);
	}
	if (record.lineIndices.empty() && record.pointIndices.empty() &&
	    record.triangleIndices.empty() && !record.bounds.isEmpty()) {
	    SbVec3f bmin;
	    SbVec3f bmax;
	    record.bounds.getBounds(bmin, bmax);
	    for (int x = 0; x < 2; x++)
		for (int y = 0; y < 2; y++)
		    for (int z = 0; z < 2; z++)
			consider(SbVec3f(x ? bmax[0] : bmin[0],
			    y ? bmax[1] : bmin[1],
			    z ? bmax[2] : bmin[2]));
	}
    }
}

int
_view_cmd_msgs(void *bs, int argc, const char **argv, const char *us, const char *ps)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    if (argc == 2 && BU_STR_EQUAL(argv[1], HELPFLAG)) {
	bu_vls_printf(gd->gedp->ged_result_str, "%s\n%s\n", us, ps);
	return 1;
    }
    if (argc == 2 && BU_STR_EQUAL(argv[1], PURPOSEFLAG)) {
	bu_vls_printf(gd->gedp->ged_result_str, "%s\n", ps);
	return 1;
    }
    return 0;
}

struct _view_independent_path {
    char *path;
    int mode;
};

typedef int (*view_core_cmd_func)(struct ged *, int, const char **);

static int
_view_call_on_gd_view(struct _ged_view_info *gd, view_core_cmd_func cmd, int argc, const char **argv)
{
    struct ged_view_context *cv = ged_view_active_ctx(gd->gedp);
    ged_view_active_ctx_set(gd->gedp, gd->cv);
    int ret = cmd(gd->gedp, argc, argv);
    ged_view_active_ctx_set(gd->gedp, cv);
    return ret;
}

static void
_view_independent_paths_free(struct _view_independent_path *paths, size_t path_cnt)
{
    if (!paths)
	return;

    for (size_t i = 0; i < path_cnt; i++) {
	if (paths[i].path)
	    bu_free(paths[i].path, "independent path");
    }

    bu_free(paths, "independent paths");
}

static int
_view_independent_paths_add(struct _view_independent_path **paths,
			    size_t *path_cnt,
			    size_t *path_cap,
			    const char *path,
			    int mode)
{
    if (!paths || !path_cnt || !path_cap || !path || !path[0])
	return BRLCAD_ERROR;

    for (size_t i = 0; i < *path_cnt; i++) {
	if ((*paths)[i].mode == mode && BU_STR_EQUAL((*paths)[i].path, path))
	    return BRLCAD_OK;
    }

    if (*path_cnt >= *path_cap) {
	size_t ncap = (*path_cap < 8) ? 8 : (*path_cap * 2);
	struct _view_independent_path *npaths = (struct _view_independent_path *)bu_realloc(
		*paths, ncap * sizeof(struct _view_independent_path), "independent paths");
	if (!npaths)
	    return BRLCAD_ERROR;
	*paths = npaths;
	*path_cap = ncap;
    }

    (*paths)[*path_cnt].path = bu_strdup(path);
    (*paths)[*path_cnt].mode = mode;
    (*path_cnt)++;

    return BRLCAD_OK;
}

static bool
_view_independent_collect_shared(struct ged *gedp,
	struct _view_independent_path **paths, size_t *path_cnt,
	size_t *path_cap)
{
    BObolSceneController *scene = ged_bobol_scene(gedp);
    if (!scene)
	return false;

    for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) ||
	    !ged_bobol_source_in_view(nullptr, summary) || !summary.visible)
	    continue;
	const char *path = summary.parentGroupPath.getLength() &&
	    !BU_STR_EQUAL(summary.parentGroupPath.getString(), "/") ?
	    summary.parentGroupPath.getString() : summary.path.getString();
	while (path && *path == '/')
	    path++;
	const int mode = summary.representationMode >= 0 ?
	    summary.representationMode : GED_DRAW_MODE_WIRE;
	if (_view_independent_paths_add(paths, path_cnt, path_cap, path,
		mode) != BRLCAD_OK)
	    return false;
    }
    return true;
}

int
_view_cmd_aet(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    const char *usage_string = "view [options] aet [vals]";
    const char *purpose_string = "get/set azimuth/elevation/twist of view";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    return _view_call_on_gd_view(gd, ged_aet_core, argc, argv);
}

int
_view_cmd_center(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    const char *usage_string = "view [options] center [vals]";
    const char *purpose_string = "get/set view center";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    return _view_call_on_gd_view(gd, ged_center_core, argc, argv);
}

static int
_view_cmd_cutting(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    const char *usage_string = "view [options] cutting [vals]";
    const char *purpose_string = "control the world-space cutting plane";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;
    return _view_call_on_gd_view(gd, ged_cutting_core, argc, argv);
}

int
_view_cmd_dir(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    const char *usage_string = "view [options] dir [-i] [x y z [twist]]";
    const char *purpose_string = "get/set view direction and optional twist";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    /* A direction query is viewdir; a supplied vector and optional twist are
     * qvrot.  Both use the eye-from-center convention, and -i consistently
     * selects its inverse. */
    return _view_call_on_gd_view(gd, argc < 4 ?
	ged_viewdir_core : ged_qvrot_core, argc, argv);
}

int
_view_cmd_eye(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    const char *usage_string = "view [options] eye [vals]";
    const char *purpose_string = "get/set view eye point";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    return _view_call_on_gd_view(gd, ged_eye_core, argc, argv);
}

int
_view_cmd_faceplate(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    const char *usage_string = "view [options] faceplate [vals]";
    const char *purpose_string = "manage faceplate view elements";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    return _view_call_on_gd_view(gd, ged_faceplate_core, argc, argv);
}

int
_view_cmd_lighting(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    const char *usage_string = "view [options] lighting [vals]";
    const char *purpose_string = "control headlight and in-scene lighting";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    return _view_call_on_gd_view(gd, ged_lighting_core, argc, argv);
}

int
_view_cmd_shading(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    const char *usage_string = "view [options] shading [vals]";
    const char *purpose_string = "control mesh normal presentation";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    return _view_call_on_gd_view(gd, ged_shading_core, argc, argv);
}

/* When a view is "independent", it displays only those objects when have been
 * added to its individual scene storage - the shared objects common to all
 * views will not be drawn.  When shifting a view from shared to independent
 * its local storage is populated with copies of the shared objects to prevent
 * an abrupt change of displayed contents, but once this setup is complete
 * further draw or erase operations in shared views will no longer alter the
 * scene object lists in the independent view.  To modify the independent
 * view's scene, it must be specifically set as the current view in libged.
 * Note also that when a view ceases to be independent, it's local object set
 * is compared to the shared object set and any objects in both are removed
 * from the local set.  However, any object in the independent list that are
 * not present in the shared set will remain, since there is no way for the
 * library to know if the intent is to preserve or remove such objects from the
 * scene.  Removal, as the destructive option, is the responsibility of the
 * application.
 *
 * Note that views may have localized scene objects even when not independent,
 * but they must be defined as view features rather than database objects. */
int
_view_cmd_independent(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    const char *usage_string = "view [options] independent <view> [0|1]";
    const char *purpose_string = "make a view independent (1) or part of the default view set (0)";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    // We know we're the independent command - start processing args
    argc--; argv++;

    struct ged *gedp = gd->gedp;
    if (!argc) {
	bu_vls_printf(gedp->ged_result_str, "no view specified\n");
	return BRLCAD_ERROR;
    }

    struct ged_view_context *view_ctx = ged_view_find_ctx(gedp, argv[0]);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, "view %s not found\n", argv[0]);
	return BRLCAD_ERROR;
    }

    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", ged_view_context_is_independent(view_ctx));
	return BRLCAD_OK;
    }

    if (BU_STR_EQUAL(argv[1], "1")) {
	if (ged_view_context_is_independent(view_ctx))
	    return BRLCAD_OK;

	struct _view_independent_path *paths = NULL;
	size_t path_cnt = 0;
	size_t path_cap = 0;
	if (!_view_independent_collect_shared(gedp, &paths, &path_cnt,
		&path_cap)) {
	    _view_independent_paths_free(paths, path_cnt);
	    bu_vls_printf(gedp->ged_result_str, "failed to snapshot shared draw state\n");
	    return BRLCAD_ERROR;
	}

	struct ged_view_context *cv = ged_view_active_ctx(gedp);
	ged_view_active_ctx_set(gedp, view_ctx);
	ged_draw_ensure_root_attached(gedp);
	ged_view_active_ctx_set(gedp, cv);

	if (!ged_view_context_scene_attached(view_ctx) ||
		ged_view_context_independent_scope_is_null(view_ctx, 1)) {
	    _view_independent_paths_free(paths, path_cnt);
	    bu_vls_printf(gedp->ged_result_str, "failed to create independent draw scope for %s\n",
		    argv[0]);
	    return BRLCAD_ERROR;
	}

	for (size_t i = 0; i < path_cnt; i++) {
	    struct bu_vls mode = BU_VLS_INIT_ZERO;
	    bu_vls_sprintf(&mode, "-m%d", paths[i].mode);
	    const char *av[7] = {"draw", "-R", "-V", NULL, NULL, NULL, NULL};
	    av[3] = _view_name(view_ctx);
	    av[4] = bu_vls_cstr(&mode);
	    av[5] = paths[i].path;
	    ged_exec_draw(gedp, 6, av);
	    bu_vls_free(&mode);
	}
	_view_independent_paths_free(paths, path_cnt);
	return BRLCAD_OK;
    }

    if (BU_STR_EQUAL(argv[1], "0")) {
	if (!ged_view_context_is_independent(view_ctx))
	    return BRLCAD_OK;
	const char *z_av[4] = {"Z", "-V", NULL, "-g"};
	z_av[2] = _view_name(view_ctx);
	ged_exec_Z(gedp, 4, z_av);
	ged_view_context_independent_scope_destroy(view_ctx);
	return BRLCAD_OK;
    }

    bu_vls_printf(gedp->ged_result_str, "invalid value supplied: %s (need 0 or 1)\n", argv[1]);
    return BRLCAD_ERROR;
}

int
_view_cmd_list(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    const char *usage_string = "view [options] ";
    const char *purpose_string = "list available views";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct ged *gedp = gd->gedp;
    struct bu_ptbl *views = ged_view_set_views_ctx(gedp);
    struct ged_view_context *active_view = ged_view_active_ctx(gedp);
    for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	struct ged_view_context *view_ctx =
	    (struct ged_view_context *)BU_PTBL_GET(views, i);
	if (view_ctx != active_view) {
	    bu_vls_printf(gedp->ged_result_str, "  %s\n", _view_name(view_ctx));
	} else {
	    bu_vls_printf(gedp->ged_result_str, "* %s\n", _view_name(view_ctx));
	}
    }

    return BRLCAD_OK;
}

int
_view_cmd_quat(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    const char *usage_string = "view [options] quat [vals]";
    const char *purpose_string = "get/set quaternion of view";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    return _view_call_on_gd_view(gd, ged_quat_core, argc, argv);
}

int
_view_cmd_selections(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    const char *usage_string = "view [options] selections [options] [args]";
    const char *purpose_string = "manipulate view selections";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    argc--; argv++;

    if (!gd->cv) {
	bu_vls_printf(gd->gedp->ged_result_str, "no current view selected\n");
	return BRLCAD_ERROR;
    }

    bu_vls_printf(gd->gedp->ged_result_str, "%zu",
		  ged_view_selection_count(gd->cv));

    return BRLCAD_OK;
}

int
_view_cmd_size(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    const char *usage_string = "view [options] size [vals]";
    const char *purpose_string = "get/set view size";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    return _view_call_on_gd_view(gd, ged_size_core, argc, argv);
}

int
_view_cmd_snap(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    const char *usage_string = "view [options] snap [vals]";
    const char *purpose_string = "snap to view elements";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    return _view_call_on_gd_view(gd, ged_view_snap, argc, argv);
}

int
_view_cmd_ypr(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    const char *usage_string = "view [options] ypr [vals]";
    const char *purpose_string = "get/set yaw/pitch/roll of view";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    return _view_call_on_gd_view(gd, ged_ypr_core, argc, argv);
}


struct vZ_opt {
    int set;
    struct bu_vls vn;
};

static
int vZ_opt_read(struct bu_vls *msg, size_t argc, const char **argv, void *set_var)
{
    struct vZ_opt *vZ = (struct vZ_opt *)set_var;
    vZ->set = 1;
    if (bu_opt_vls(msg, argc, argv, (void *)&vZ->vn) == 1) {
	return 1;
    }
    return 0;
}

int
_view_cmd_vZ(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view [options] vZ [opts] [val|x y z]";
    const char *purpose_string = "get/set/calc view data vZ value";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    int print_help = 0;
    struct vZ_opt calc_near = { 0, BU_VLS_INIT_ZERO };
    struct vZ_opt calc_far = { 0, BU_VLS_INIT_ZERO };
    struct bu_opt_desc d[4];
    BU_OPT(d[0], "h", "help", "",       NULL,  &print_help, "Print help");
    BU_OPT(d[1], "N", "near", "[feature]",  &vZ_opt_read,  &calc_near,  "Find vZ value of closest view feature vertex");
    BU_OPT(d[2], "F", "far",  "[feature]",  &vZ_opt_read,  &calc_far,   "Find vZ value of furthest view feature vertex");
    BU_OPT_NULL(d[3]);

    // We know we're the vZ command - start processing args
    argc--; argv++;

    int ac = bu_opt_parse(NULL, argc, argv, d);
    argc = ac;

    if (print_help || (calc_near.set && calc_far.set)) {
	bu_vls_printf(gedp->ged_result_str, "[WARNING] this command is deprecated - vZ values should be set on data objects\n\nUsage:\n%s", usage_string);
	return GED_HELP;
    }

    int calc_mode = -1;
    struct bu_vls calc_target = BU_VLS_INIT_ZERO;
    if (calc_near.set) {
	calc_mode = 0;
	bu_vls_sprintf(&calc_target, "%s", bu_vls_cstr(&calc_near.vn));
	bu_vls_free(&calc_near.vn);
    }
    if (calc_far.set) {
	calc_mode = 1;
	bu_vls_sprintf(&calc_target, "%s", bu_vls_cstr(&calc_far.vn));
	bu_vls_free(&calc_far.vn);
    }

    if (calc_mode != -1) {
	if (bu_vls_strlen(&calc_target)) {
	    fastf_t vZ = 0.0;
	    if (_view_vZ_feature_depth(gedp, gd->cv,
		    bu_vls_cstr(&calc_target),
			calc_mode, &vZ)) {
		bu_vls_sprintf(gedp->ged_result_str, "%0.15f", vZ);
		return BRLCAD_OK;
	    } else {
		bu_vls_sprintf(gedp->ged_result_str, "View feature %s not found", bu_vls_cstr(&calc_target));
		bu_vls_free(&calc_target);
		return BRLCAD_ERROR;
	    }
	} else {
	    /* Check all drawn view features and database leaves. */
	    struct _view_vZ_ctx ctx;
	    bv_model2view_get(ctx.model2view,
		    bv_context_view_const((const struct bv_context *)gd->cv));
	    ctx.calc_mode = calc_mode;
	    ctx.vZ = (calc_mode) ? -DBL_MAX : DBL_MAX;
	    ctx.have_vz = 0;
	    _view_vZ_visit_features(gedp, gd->cv, &ctx);
	    _view_vZ_visit_db_exports(gd->cv, &ctx);
	    if (ctx.have_vz) {
		bu_vls_sprintf(gedp->ged_result_str, "%0.15f", ctx.vZ);
	    }
	}
	return BRLCAD_OK;
    }

    /* The legacy get/set scratch path that read or wrote `gv_tcl->gv_data_vZ`
     * was removed from libged.  `gv_data_vZ` is a Tcl-mode editing scratch
     * scalar with no meaning to non-Tcl libged callers; callers should set vZ
     * values on data objects directly.  Tcl callers that still need the
     * scratch can use the `data_vZ` command exposed by libtclcad (commands.c),
     * which keeps the gv_tcl-resident scalar consistent with Tcl editing-mode
     * plumbing. */
    bu_vls_printf(gedp->ged_result_str, "[WARNING] this command is deprecated - vZ values should be set on data objects.\n\nUsage:\n%s", usage_string);
    return GED_HELP;
}
int
_view_cmd_width(void *ds, int argc, const char **argv)
{
    const char *usage_string = "view [options] width";
    const char *purpose_string = "report current width in pixels of view.";
    if (_view_cmd_msgs(ds, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    argc--; argv++;

    struct _ged_view_info *gd = (struct _ged_view_info *)ds;
    bu_vls_printf(gd->gedp->ged_result_str, "%d\n",
	    bv_width_get(bv_context_view_const((const struct bv_context *)gd->cv)));
    return BRLCAD_OK;
}

int
_view_cmd_height(void *ds, int argc, const char **argv)
{
    const char *usage_string = "view [options] height";
    const char *purpose_string = "report current height in pixels of view.";
    if (_view_cmd_msgs(ds, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    argc--; argv++;

    struct _ged_view_info *gd = (struct _ged_view_info *)ds;
    bu_vls_printf(gd->gedp->ged_result_str, "%d\n",
	    bv_height_get(bv_context_view_const((const struct bv_context *)gd->cv)));
    return BRLCAD_OK;
}

int
_view_cmd_print(struct ged *gedp, int argc, const char **argv)
{
    // ae
    ged_aet_core(gedp, argc, argv);
    char* ae = bu_vls_strdup(gedp->ged_result_str);

    // dir
    ged_viewdir_core(gedp, argc, argv);
    char* dir = bu_vls_strdup(gedp->ged_result_str);

    // center
    ged_center_core(gedp, argc, argv);
    char* center = bu_vls_strdup(gedp->ged_result_str);

    // eye
    ged_eye_core(gedp, argc, argv);
    char* eye = bu_vls_strdup(gedp->ged_result_str);

    // size
    ged_size_core(gedp, argc, argv);
    char* size = bu_vls_strdup(gedp->ged_result_str);

    // print
    bu_vls_trunc(gedp->ged_result_str, 0);
    bu_vls_printf(gedp->ged_result_str, "    ae: %s\n    dir: <%s>\n    center: (%s)\n    eye: (%s)\n    size: %s\n", ae, dir, center, eye, size);

    // cleanup
    bu_free(ae, "ae str free");
    bu_free(dir, "dir str free");
    bu_free(center, "center str free");
    bu_free(eye, "eye str free");
    bu_free(size, "size str free");

    return BRLCAD_OK;
}

static int
_view_cmd_core(void *bs, int argc, const char **argv, view_core_cmd_func func)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    return _view_call_on_gd_view(gd, func, argc, argv);
}

static int
_view_cmd_auto(void *bs, int argc, const char **argv)
{
    const char *usage_string = "view [options] auto [options] [object ...]";
    const char *purpose_string = "size and center the view to frame geometry";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    return _view_cmd_core(bs, argc, argv, ged_autoview_core);
}

static int
_view_cmd_lookat(void *bs, int argc, const char **argv)
{
    const char *usage_string = "view [options] lookat x y z";
    const char *purpose_string = "point the view at model coordinates";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    return _view_cmd_core(bs, argc, argv, ged_lookat_core);
}

static int
_view_cmd_print_subcmd(void *bs, int argc, const char **argv)
{
    const char *usage_string = "view [options] print";
    const char *purpose_string = "print the current view parameters";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    return _view_cmd_core(bs, argc, argv, _view_cmd_print);
}

static int
_view_cmd_save(void *bs, int argc, const char **argv)
{
    const char *usage_string = "view [options] save [-e command] [-i input] [-l log] [-o output] file [args]";
    const char *purpose_string = "save the current view as a raytrace script";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    return _view_cmd_core(bs, argc, argv, ged_saveview_core);
}

int
_view_cmd_knob(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    const char *usage_string = "view [options] knob [vals]";
    const char *purpose_string = "low level view rotate/translate/scale operations";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    return _view_call_on_gd_view(gd, ged_knob_core, argc, argv);
}

const struct bu_cmdtab _view_cmds[] = {
    { "ae",         _view_cmd_aet},
    { "aet",        _view_cmd_aet},
    { "auto",       _view_cmd_auto},
    { "autoview",   _view_cmd_auto},
    { "center",     _view_cmd_center},
    { "cutting",    _view_cmd_cutting},
    { "dir",        _view_cmd_dir},
    { "eye",        _view_cmd_eye},
    { "annotation", _view_cmd_annotation},
    { "db",         _view_cmd_db_objects},
    { "faceplate",  _view_cmd_faceplate},
    { "height",     _view_cmd_height},
    { "independent",_view_cmd_independent},
    { "knob",       _view_cmd_knob},
    { "lighting",   _view_cmd_lighting},
    { "list",       _view_cmd_list},
    { "lod",        _view_cmd_lod},
    { "lookat",     _view_cmd_lookat},
    { "obj",        _view_cmd_feature},
    { "objs",       _view_cmd_feature},
    { "print",      _view_cmd_print_subcmd},
    { "feature",    _view_cmd_feature},
    { "polygon",    _view_cmd_polygon},
    { "quat",       _view_cmd_quat},
    { "save",       _view_cmd_save},
    { "saveview",   _view_cmd_save},
    { "shading",    _view_cmd_shading},
    { "selections", _view_cmd_selections},
    { "size",       _view_cmd_size},
    { "snap",       _view_cmd_snap},
    { "vZ",         _view_cmd_vZ},
    { "width",      _view_cmd_width},
    { "ypr",        _view_cmd_ypr},
    { (char *)NULL,      NULL}
};

/* Feature stores are intentionally renderer-neutral.  A successful command
 * can therefore change a shared scene node without changing the libbv camera
 * hash that qged traditionally uses to schedule a paint.  Only these command
 * families can mutate retained view content; mutation itself is detected by
 * the shared controller's render-request serial. */
static bool
_view_cmd_can_mutate_obol_content(const char *command)
{
    return command && (BU_STR_EQUAL(command, "annotation") ||
	BU_STR_EQUAL(command, "feature") || BU_STR_EQUAL(command, "obj") ||
	BU_STR_EQUAL(command, "objs") || BU_STR_EQUAL(command, "polygon"));
}

struct _view_shared_presentation_revision {
    uint64_t features = 0;
    uint64_t polygons = 0;

    bool operator==(const _view_shared_presentation_revision &other) const
    {
	return features == other.features && polygons == other.polygons;
    }

    bool operator!=(const _view_shared_presentation_revision &other) const
    {
	return !(*this == other);
    }
};

static _view_shared_presentation_revision
_view_shared_presentation_revision_get(BObolViewController *controller)
{
    _view_shared_presentation_revision revision;
    if (!controller)
	return revision;
    revision.features = controller->features().presentationRevision();
    revision.polygons = controller->polygons().presentationRevision();
    return revision;
}

int
ged_view_core(struct ged *gedp, int argc, const char *argv[])
{
    int help = 0;
    struct _ged_view_info gd;
    memset(&gd, 0, sizeof(gd));
    gd.gedp = gedp;
    gd.cmds = _view_cmds;
    gd.cv = NULL;
    gd.gobj_dbpath = NULL;
    gd.verbosity = 0;

    // Sanity
    if (UNLIKELY(!gedp || !argc || !argv)) {
	return BRLCAD_ERROR;
    }

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    // We know we're the dm command - start processing args
    argc--; argv++;

    // See if we have any high level options set
    struct bu_vls vname = BU_VLS_INIT_ZERO;
    struct bu_opt_desc d[4];
    BU_OPT(d[0], "h", "help",    "",      NULL,               &help,         "Print help");
    BU_OPT(d[1], "v", "verbose", "",      &bu_opt_incr_long,  &gd.verbosity, "Verbose output");
    BU_OPT(d[2], "V", "view",    "name",  &bu_opt_vls,        &vname,        "Specified view (default is GED current)");
    BU_OPT_NULL(d[3]);

    gd.gopts = d;

    // High level options are only defined prior to the subcommand
    int cmd_pos = -1;
    for (int i = 0; i < argc; i++) {
	if (bu_cmd_valid(_view_cmds, argv[i]) == BRLCAD_OK ||
	    _ged_cmd_namespace_has_child("view", argv[i])) {
	    cmd_pos = i;
	    break;
	}
    }

    // Clear out any high level opts prior to subcommand
    int acnt = (cmd_pos >= 0) ? cmd_pos : argc;
    int ac_ret = bu_opt_parse(NULL, acnt, argv, d);
    if (ac_ret) {
	help = 1;
    } else {
	for (int i = 0; i < acnt; i++) {
	    argc--; argv++;
	}
    }

    if (help) {
	if (cmd_pos >= 0) {
	    argc = argc - cmd_pos;
	    argv = &argv[cmd_pos];
	    _ged_subcmd_help(gedp, (struct bu_opt_desc *)d, (const struct bu_cmdtab *)_view_cmds, "view", "[options] subcommand [args]", &gd, argc, argv);
	} else {
	    _ged_subcmd_help(gedp, (struct bu_opt_desc *)d, (const struct bu_cmdtab *)_view_cmds, "view", "[options] subcommand [args]", &gd, 0, NULL);
	}
	_ged_cmd_namespace_help(gedp, "view", _view_cmds);
	bu_vls_free(&vname);
	return BRLCAD_OK;
    }

    // Must have a subcommand
    if (cmd_pos == -1) {
	bu_vls_printf(gedp->ged_result_str, ": no valid subcommand specified\n");
	_ged_subcmd_help(gedp, (struct bu_opt_desc *)d, (const struct bu_cmdtab *)_view_cmds, "view", "[options] subcommand [args]", &gd, 0, NULL);
	_ged_cmd_namespace_help(gedp, "view", _view_cmds);
	bu_vls_free(&vname);
	return BRLCAD_ERROR;
    }

    // Either a view was specified, or we use the current view
    if (bu_vls_strlen(&vname)) {
	gd.cv = ged_view_find_ctx(gedp, bu_vls_cstr(&vname));
	if (!gd.cv) {
	    bu_vls_printf(gedp->ged_result_str, ": invalid view name: %s", bu_vls_cstr(&vname));
	    bu_vls_free(&vname);
	    return BRLCAD_ERROR;
	}
    } else {
	gd.cv = ged_view_active_ctx(gedp);
    }

    /* LoD cache inspection and maintenance are scoped to the open database,
     * not to a camera.  Permit that one command family in headless clients;
     * _view_cmd_lod handles it before consulting gd.cv. */
    const int database_only_lod_cache =
	argc > 1 && BU_STR_EQUAL(argv[0], "lod") &&
	BU_STR_EQUAL(argv[1], "cache");
    if (!gd.cv && !database_only_lod_cache) {
	bu_vls_printf(gedp->ged_result_str, ": no view specified and no view listed as current in GED");
	bu_vls_free(&vname);
	return BRLCAD_ERROR;
    }

    BObolViewController *shared_controller =
	_view_cmd_can_mutate_obol_content(argv[0]) ?
	ged_bobol_shared_view_controller(gedp) : NULL;
    const _view_shared_presentation_revision shared_presentation_revision =
	_view_shared_presentation_revision_get(shared_controller);
    int ret;
    if (bu_cmd_valid(_view_cmds, argv[0]) == BRLCAD_OK &&
	bu_cmd(_view_cmds, argc, argv, 0, (void *)&gd, &ret) == BRLCAD_OK) {
	if (ret == BRLCAD_OK && shared_controller &&
	    _view_shared_presentation_revision_get(shared_controller) !=
		shared_presentation_revision)
	    ged_bobol_shared_view_presentation_request(gedp,
		"ged-shared-view-content");
	bu_vls_free(&vname);
	return ret;
    }

    struct ged_view_context *saved_view = ged_view_active_ctx(gedp);
    ged_view_active_ctx_set(gedp, gd.cv);
    const int plugin_found = _ged_cmd_namespace_exec(gedp, "view",
	argc, argv, &ret);
    ged_view_active_ctx_set(gedp, saved_view);
    if (plugin_found) {
	bu_vls_free(&vname);
	return ret;
    }

    bu_vls_printf(gedp->ged_result_str, "subcommand %s not defined", argv[0]);

    bu_vls_free(&vname);
    return BRLCAD_ERROR;
}

#include "../include/plugin.h"

#define GED_VIEW_COMMANDS(X, XID) \
    X(ae, ged_aet_core, GED_CMD_DEFAULT) \
    X(aet, ged_aet_core, GED_CMD_DEFAULT) \
    X(autoview, ged_autoview_core, GED_CMD_DEFAULT) \
    X(center, ged_center_core, GED_CMD_DEFAULT) \
    X(data_lines, ged_view_data_lines, GED_CMD_DEFAULT) \
    X(eye, ged_eye_core, GED_CMD_DEFAULT) \
    X(eye_pt, ged_eye_core, GED_CMD_DEFAULT) \
    X(lookat, ged_lookat_core, GED_CMD_DEFAULT) \
    X(print, _view_cmd_print, GED_CMD_DEFAULT) \
    X(quat, ged_quat_core, GED_CMD_DEFAULT) \
    X(qvrot, ged_qvrot_core, GED_CMD_DEFAULT) \
    X(saveview, ged_saveview_core, GED_CMD_DEFAULT) \
    X(sdata_lines, ged_view_data_lines, GED_CMD_DEFAULT) \
    X(size, ged_size_core, GED_CMD_DEFAULT) \
    X(view, ged_view_core, GED_CMD_DEFAULT) \
    X(view2, ged_view_core, GED_CMD_DEFAULT) \
    X(view_func, ged_view_core, GED_CMD_DEFAULT) \
    X(viewdir, ged_viewdir_core, GED_CMD_DEFAULT) \
    X(ypr, ged_ypr_core, GED_CMD_DEFAULT) \

GED_DECLARE_COMMAND_SET(GED_VIEW_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_view", 1, GED_VIEW_COMMANDS)

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
