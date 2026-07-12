/*                          M O U S E . C
 * BRL-CAD
 *
 * Copyright (c) 2000-2026 United States Government as represented by
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

#include "bv.h"
#include "bu/path.h"
#include "ged/view.h"
#include "tclcad.h"
#include "brlobol/display_endpoint.h"

/* Private headers */
#include "ged/draw.h"
#include "./tclcad_private.h"
#include "./view/view.h"
#include "./draw_view_move_helpers.h"

static const char *
tclcad_mouse_view_name(const void *view_ctx)
{
    const char *name = bv_context_name_get(
	    (const struct bv_context *)view_ctx);
    return name ? name : "";
}

static struct bv *
tclcad_mouse_bv(void *view_ctx)
{
    return bv_context_view((struct bv_context *)view_ctx);
}

static const struct bv *
tclcad_mouse_bv_const(const void *view_ctx)
{
    return bv_context_view_const((const struct bv_context *)view_ctx);
}

static int
tclcad_mouse_display_width(const void *view_ctx)
{
    brlobol_display_endpoint_t *endpoint =
	ged_view_context_display_endpoint_get(view_ctx);
    struct brlobol_endpoint_property_value property =
	BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    if (endpoint && brlobol_display_endpoint_property_get(endpoint,
	    "endpoint.width", &property) == BRLOBOL_ENDPOINT_PROPERTY_OK &&
	property.type == BRLOBOL_ENDPOINT_PROPERTY_UINT && property.uint_value)
	return (int)property.uint_value;
    return bv_context_width_get((const struct bv_context *)view_ctx);
}

static int
tclcad_mouse_display_height(const void *view_ctx)
{
    brlobol_display_endpoint_t *endpoint =
	ged_view_context_display_endpoint_get(view_ctx);
    struct brlobol_endpoint_property_value property =
	BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    if (endpoint && brlobol_display_endpoint_property_get(endpoint,
	    "endpoint.height", &property) == BRLOBOL_ENDPOINT_PROPERTY_OK &&
	property.type == BRLOBOL_ENDPOINT_PROPERTY_UINT && property.uint_value)
	return (int)property.uint_value;
    return bv_context_height_get((const struct bv_context *)view_ctx);
}

static struct bu_vls *
tclcad_mouse_display_pathname(const void *view_ctx)
{
    return tclcad_view_pathname_vls(view_ctx);
}

static void
tclcad_mouse_sync_dimensions(void *view_ctx)
{
    int width = tclcad_mouse_display_width(view_ctx);
    int height = tclcad_mouse_display_height(view_ctx);
    if (width > 0 && height > 0)
	bv_context_dimensions_set((struct bv_context *)view_ctx, width, height);
}

static void
tclcad_mouse_view_inv_rotation(mat_t inv_rotation, const void *view_ctx)
{
    mat_t view_rotation;

    bv_rotation_get(view_rotation, tclcad_mouse_bv_const(view_ctx));
    bn_mat_inv(inv_rotation, view_rotation);
}

static void
tclcad_mouse_previous_get_set(fastf_t *prev_x, fastf_t *prev_y,
	void *view_ctx, fastf_t x, fastf_t y)
{
    if (prev_x)
	*prev_x = 0.0;
    if (prev_y)
	*prev_y = 0.0;
    (void)bv_previous_mouse_get(prev_x, prev_y,
	    tclcad_mouse_bv_const(view_ctx));
    (void)bv_previous_mouse_set(tclcad_mouse_bv(view_ctx), x, y);
}

static void
tclcad_mouse_delta_settings(struct bv_mouse_delta_settings *settings,
	const void *view_ctx)
{
    struct bv_mouse_delta_settings zero = BV_MOUSE_DELTA_SETTINGS_INIT;

    if (!settings)
	return;

    *settings = zero;
    (void)bv_mouse_delta_settings_get(settings, tclcad_mouse_bv_const(view_ctx));
}

static void
tclcad_mouse_clamp_delta(fastf_t *dx, fastf_t *dy, const void *view_ctx)
{
    struct bv_mouse_delta_settings settings = BV_MOUSE_DELTA_SETTINGS_INIT;

    if (!dx || !dy)
	return;

    tclcad_mouse_delta_settings(&settings, view_ctx);
    if (*dx < settings.min_delta)
	*dx = settings.min_delta;
    else if (settings.max_delta < *dx)
	*dx = settings.max_delta;

    if (*dy < settings.min_delta)
	*dy = settings.min_delta;
    else if (settings.max_delta < *dy)
	*dy = settings.max_delta;
}

static unsigned long long
tclcad_mouse_prepare_snap(void *view_ctx)
{
    struct bv *view = tclcad_mouse_bv(view_ctx);
    int flags;

    if (!view)
	return 0ULL;

    flags = bv_snap_source_flags_get(view);
    (void)bv_snap_source_flags_set(view, flags | BV_SNAP_TCL);
    return bv_snap_kind_mask_get(view);
}

static int
tclcad_mouse_snap_point_2d(void *view_ctx, fastf_t *vx, fastf_t *vy,
	unsigned long long kinds)
{
    return (kinds & BV_SNAP_KIND_GRID) ?
	bv_snap_grid_2d(tclcad_mouse_bv_const(view_ctx), vx, vy) : 0;
}

static fastf_t
tclcad_mouse_rotate_scale(const void *view_ctx)
{
    struct bv_mouse_delta_settings settings = BV_MOUSE_DELTA_SETTINGS_INIT;

    tclcad_mouse_delta_settings(&settings, view_ctx);
    return settings.rotate_scale;
}

static fastf_t
tclcad_mouse_scale_scale(const void *view_ctx)
{
    struct bv_mouse_delta_settings settings = BV_MOUSE_DELTA_SETTINGS_INIT;

    tclcad_mouse_delta_settings(&settings, view_ctx);
    return settings.scale_scale;
}

int
to_get_prev_mouse(struct ged *gedp,
		  int argc,
		  const char *argv[],
		  ged_func_ptr UNUSED(func),
		  const char *usage,
		  int UNUSED(maxargs))
{
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 2) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    (void)bv_previous_mouse_get(&prev_x, &prev_y,
	    tclcad_mouse_bv_const(gdvp));
    bu_vls_printf(gedp->ged_result_str, "%d %d", (int)prev_x, (int)prev_y);
    return BRLCAD_OK;
}


int
to_mouse_append_pnt_common(struct ged *gedp,
			   int argc,
			   const char *argv[],
			   ged_func_ptr func,
			   const char *usage,
			   int UNUSED(maxargs))
{
    int ret;
    const char *av[4];
    point_t view;
    struct bu_vls pt_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 5) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[3], "%lf", &x) != 1 ||
	bu_sscanf(argv[4], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    tclcad_mouse_sync_dimensions(gdvp);
    bv_screen_to_view(&x, &y, tclcad_mouse_bv_const(gdvp), x, y);
    VSET(view, x, y, 0.0);

    tclcad_mouse_sync_dimensions(gdvp);

    ged_view_active_ctx_set(gedp, gdvp);
    {
	unsigned long long snap_kinds = tclcad_mouse_prepare_snap(gdvp);
	if (snap_kinds)
	    tclcad_mouse_snap_point_2d(gdvp, &view[X], &view[Y], snap_kinds);
    }

    bu_vls_printf(&pt_vls, "%lf %lf %lf", view[X], view[Y], view[Z]);

    ged_view_active_ctx_set(gedp, gdvp);
    av[0] = (char *)argv[0];
    av[1] = (char *)argv[2];
    av[2] = bu_vls_addr(&pt_vls);
    av[3] = (char *)0;

    ret = (*func)(gedp, 3, (const char **)av);
    bu_vls_free(&pt_vls);

    if (ret == BRLCAD_OK) {
	av[0] = "draw";
	av[1] = (char *)argv[2];
	av[2] = (char *)0;
	to_edit_redraw(gedp, 2, (const char **)av);
    }

    return BRLCAD_OK;
}


int
to_mouse_brep_selection_append(struct ged *gedp,
			       int argc,
			       const char *argv[],
			       ged_func_ptr UNUSED(func),
			       const char *usage,
			       int maxargs)
{
    const char *cmd_argv[11] = {"brep", NULL, "selection", "append", "active"};
    int ret, cmd_argc = (int)(sizeof(cmd_argv) / sizeof(const char *));
    char *brep_name;
    char *end;
    struct bu_vls bindings = BU_VLS_INIT_ZERO;
    struct bu_vls start[] = {BU_VLS_INIT_ZERO, BU_VLS_INIT_ZERO, BU_VLS_INIT_ZERO};
    struct bu_vls dir[] = {BU_VLS_INIT_ZERO, BU_VLS_INIT_ZERO, BU_VLS_INIT_ZERO};
    point_t screen_pt, view_pt, model_pt;
    vect_t view_dir, model_dir;
    mat_t invRot;
    mat_t view2model, view_rotation;

    if (argc != maxargs) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }


    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* parse args */
    brep_name = bu_path_basename(argv[2], NULL);

    screen_pt[X] = strtol(argv[3], &end, 10);
    if (*end != '\0') {
	bu_vls_printf(gedp->ged_result_str, "ERROR: bad x value %f\n", screen_pt[X]);
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    screen_pt[Y] = strtol(argv[4], &end, 10);
    if (*end != '\0') {
	bu_vls_printf(gedp->ged_result_str, "ERROR: bad y value: %f\n", screen_pt[Y]);
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    /* stash point coordinates for future drag handling */
    bv_previous_mouse_set(tclcad_mouse_bv(gdvp), screen_pt[X], screen_pt[Y]);

    /* convert screen point to model-space start point and direction */
    tclcad_mouse_sync_dimensions(gdvp);
    bv_screen_to_view(&view_pt[X], &view_pt[Y], tclcad_mouse_bv_const(gdvp),
	    screen_pt[X], screen_pt[Y]);
    view_pt[Z] = 1.0;

    bv_view2model_get(view2model, tclcad_mouse_bv_const(gdvp));
    MAT4X3PNT(model_pt, view2model, view_pt);

    VSET(view_dir, 0.0, 0.0, -1.0);
    bv_rotation_get(view_rotation, tclcad_mouse_bv_const(gdvp));
    bn_mat_inv(invRot, view_rotation);
    MAT4X3PNT(model_dir, invRot, view_dir);

    /* brep brep_name selection append selection_name startx starty startz dirx diry dirz */
    bu_vls_printf(&start[X], "%f", model_pt[X]);
    bu_vls_printf(&start[Y], "%f", model_pt[Y]);
    bu_vls_printf(&start[Z], "%f", model_pt[Z]);

    cmd_argv[1] = brep_name;
    cmd_argv[5] = bu_vls_addr(&start[X]);
    cmd_argv[6] = bu_vls_addr(&start[Y]);
    cmd_argv[7] = bu_vls_addr(&start[Z]);

    bu_vls_printf(&dir[X], "%f", model_dir[X]);
    bu_vls_printf(&dir[Y], "%f", model_dir[Y]);
    bu_vls_printf(&dir[Z], "%f", model_dir[Z]);

    cmd_argv[8] = bu_vls_addr(&dir[X]);
    cmd_argv[9] = bu_vls_addr(&dir[Y]);
    cmd_argv[10] = bu_vls_addr(&dir[Z]);

    ged_view_active_ctx_set(gedp, gdvp);
    ret = ged_exec_brep(gedp, cmd_argc, cmd_argv);

    bu_vls_free(&start[X]);
    bu_vls_free(&start[Y]);
    bu_vls_free(&start[Z]);
    bu_vls_free(&dir[X]);
    bu_vls_free(&dir[Y]);
    bu_vls_free(&dir[Z]);

    if (ret != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }

    struct bu_vls *dname = tclcad_mouse_display_pathname(gdvp);
    if (dname && bu_vls_strlen(dname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_brep_selection_translate %s %s %%x %%y; "
		      "%s brep %s plot SCV}",
		      bu_vls_cstr(dname),
		      bu_vls_cstr(&current_top->to_gedp->go_name),
		      tclcad_mouse_view_name(gdvp),
		      brep_name,
		      bu_vls_cstr(&current_top->to_gedp->go_name),
		      brep_name);
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    bu_free((void *)brep_name, "brep_name");

    return BRLCAD_OK;
}


int
to_mouse_brep_selection_translate(struct ged *gedp,
				  int argc,
				  const char *argv[],
				  ged_func_ptr UNUSED(func),
				  const char *usage,
				  int maxargs)
{
    const char *cmd_argv[8] = {"brep", NULL, "selection", "translate", "active"};
    int ret, cmd_argc = (int)(sizeof(cmd_argv) / sizeof(const char *));
    char *brep_name;
    char *end;
    point_t screen_end, view_start, view_end, model_start, model_end;
    vect_t model_delta;
    mat_t view2model;
    struct bu_vls delta[] = {BU_VLS_INIT_ZERO, BU_VLS_INIT_ZERO, BU_VLS_INIT_ZERO};

    if (argc != maxargs) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    brep_name = bu_path_basename(argv[2], NULL);

    screen_end[X] = strtol(argv[3], &end, 10);
    if (*end != '\0') {
	bu_vls_printf(gedp->ged_result_str, "ERROR: bad x value %f\n", screen_end[X]);
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    screen_end[Y] = strtol(argv[4], &end, 10);
    if (*end != '\0') {
	bu_vls_printf(gedp->ged_result_str, "ERROR: bad y value: %f\n", screen_end[Y]);
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    /* convert screen-space delta to model-space delta */
    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    (void)bv_previous_mouse_get(&prev_x, &prev_y,
	    tclcad_mouse_bv_const(gdvp));
    tclcad_mouse_sync_dimensions(gdvp);
    bv_screen_to_view(&view_start[X], &view_start[Y],
	    tclcad_mouse_bv_const(gdvp),
	    prev_x, prev_y);
    view_start[Z] = 1;
    bv_view2model_get(view2model, tclcad_mouse_bv_const(gdvp));
    MAT4X3PNT(model_start, view2model, view_start);

    tclcad_mouse_sync_dimensions(gdvp);
    bv_screen_to_view(&view_end[X], &view_end[Y],
	    tclcad_mouse_bv_const(gdvp),
	    screen_end[X], screen_end[Y]);
    view_end[Z] = 1;
    MAT4X3PNT(model_end, view2model, view_end);

    VSUB2(model_delta, model_end, model_start);

    bu_vls_printf(&delta[X], "%f", model_delta[X]);
    bu_vls_printf(&delta[Y], "%f", model_delta[Y]);
    bu_vls_printf(&delta[Z], "%f", model_delta[Z]);

    cmd_argv[1] = brep_name;
    cmd_argv[5] = bu_vls_addr(&delta[X]);
    cmd_argv[6] = bu_vls_addr(&delta[Y]);
    cmd_argv[7] = bu_vls_addr(&delta[Z]);

    ret = ged_exec_brep(gedp, cmd_argc, cmd_argv);

    bu_free((void *)brep_name, "brep_name");
    bu_vls_free(&delta[X]);
    bu_vls_free(&delta[Y]);
    bu_vls_free(&delta[Z]);

    if (ret != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }

    /* need to tell front-end that we've modified the db */
    tclcad_eval_noresult(current_top->to_interp, "$::ArcherCore::application setSave", 0, NULL);

    bv_previous_mouse_set(tclcad_mouse_bv(gdvp), screen_end[X], screen_end[Y]);

    cmd_argc = 2;
    cmd_argv[0] = "draw";
    cmd_argv[1] = argv[2];
    cmd_argv[2] = NULL;
    ret = to_edit_redraw(gedp, cmd_argc, cmd_argv);

    return ret;
}


int
to_mouse_constrain_rot(struct ged *gedp,
		       int argc,
		       const char *argv[],
		       ged_func_ptr UNUSED(func),
		       const char *usage,
		       int UNUSED(maxargs))
{
    int ret;
    int ac;
    const char *av[4];
    fastf_t dx, dy;
    fastf_t sf;
    struct bu_vls rot_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 5) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }


    if ((argv[2][0] != 'x' && argv[2][0] != 'y' && argv[2][0] != 'z') || argv[2][1] != '\0') {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[3], "%lf", &x) != 1 ||
	bu_sscanf(argv[4], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    tclcad_mouse_previous_get_set(&prev_x, &prev_y, gdvp, x, y);
    dx = x - prev_x;
    dy = prev_y - y;

    tclcad_mouse_clamp_delta(&dx, &dy, gdvp);

    dx *= tclcad_mouse_rotate_scale(gdvp);
    dy *= tclcad_mouse_rotate_scale(gdvp);

    if (fabs(dx) > fabs(dy))
	sf = dx;
    else
	sf = dy;

    switch (argv[2][0]) {
	case 'x':
	    bu_vls_printf(&rot_vls, "%lf 0 0", -sf);
	    break;
	case 'y':
	    bu_vls_printf(&rot_vls, "0 %lf 0", -sf);
	    break;
	case 'z':
	    bu_vls_printf(&rot_vls, "0 0 %lf", -sf);
    }

    ged_view_active_ctx_set(gedp, gdvp);
    ac = 3;
    av[0] = "rot";
    av[1] = "-m";
    av[2] = bu_vls_addr(&rot_vls);
    av[3] = (char *)0;

    ret = ged_exec_rot(gedp, ac, (const char **)av);
    bu_vls_free(&rot_vls);

    if (ret == BRLCAD_OK) {
	struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(gdvp);
	if (0 < bu_vls_strlen(&tvd->gdv_callback)) {
	    tclcad_eval_noresult(current_top->to_interp, bu_vls_addr(&tvd->gdv_callback), 0, NULL);
	}

	to_refresh_view(gdvp);
    }

    return BRLCAD_OK;
}


int
to_mouse_constrain_trans(struct ged *gedp,
			 int argc,
			 const char *argv[],
			 ged_func_ptr UNUSED(func),
			 const char *usage,
			 int UNUSED(maxargs))
{
    int width;
    int ret;
    int ac;
    const char *av[4];
    fastf_t dx, dy;
    fastf_t sf;
    fastf_t inv_width;
    struct bu_vls tran_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 5) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if ((argv[2][0] != 'x' && argv[2][0] != 'y' && argv[2][0] != 'z') || argv[2][1] != '\0') {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[3], "%lf", &x) != 1 ||
	bu_sscanf(argv[4], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    tclcad_mouse_previous_get_set(&prev_x, &prev_y, gdvp, x, y);
    dx = x - prev_x;
    dy = prev_y - y;

    tclcad_mouse_clamp_delta(&dx, &dy, gdvp);

    width = tclcad_mouse_display_width(gdvp);
    inv_width = 1.0 / (fastf_t)width;
    dx *= inv_width * bv_size_get(tclcad_mouse_bv_const(gdvp)) *
	gedp->dbip->dbi_local2base;
    dy *= inv_width * bv_size_get(tclcad_mouse_bv_const(gdvp)) *
	gedp->dbip->dbi_local2base;

    if (fabs(dx) > fabs(dy))
	sf = dx;
    else
	sf = dy;

    switch (argv[2][0]) {
	case 'x':
	    bu_vls_printf(&tran_vls, "%lf 0 0", -sf);
	    break;
	case 'y':
	    bu_vls_printf(&tran_vls, "0 %lf 0", -sf);
	    break;
	case 'z':
	    bu_vls_printf(&tran_vls, "0 0 %lf", -sf);
    }

    ged_view_active_ctx_set(gedp, gdvp);
    ac = 3;
    av[0] = "tra";
    av[1] = "-m";
    av[2] = bu_vls_addr(&tran_vls);
    av[3] = (char *)0;

    ret = ged_exec_tra(gedp, ac, (const char **)av);
    bu_vls_free(&tran_vls);

    if (ret == BRLCAD_OK) {
	struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(gdvp);
	if (0 < bu_vls_strlen(&tvd->gdv_callback)) {
	    tclcad_eval_noresult(current_top->to_interp, bu_vls_addr(&tvd->gdv_callback), 0, NULL);
	}

	to_refresh_view(gdvp);
    }

    return BRLCAD_OK;
}


int
to_mouse_find_arb_edge(struct ged *gedp,
		       int argc,
		       const char *argv[],
		       ged_func_ptr UNUSED(func),
		       const char *usage,
		       int UNUSED(maxargs))
{
    const char *av[6];
    point_t view;
    struct bu_vls pt_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 6) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[3], "%lf", &x) != 1 ||
	bu_sscanf(argv[4], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    tclcad_mouse_sync_dimensions(gdvp);
    bv_screen_to_view(&x, &y, tclcad_mouse_bv_const(gdvp), x, y);
    VSET(view, x, y, 0.0);

    bu_vls_printf(&pt_vls, "%lf %lf %lf", view[X], view[Y], view[Z]);

    ged_view_active_ctx_set(gedp, gdvp);
    av[0] = "find_arb_edge_nearest_pnt";
    av[1] = (char *)argv[2];
    av[2] = bu_vls_addr(&pt_vls);
    av[3] = (char *)argv[5];
    av[4] = (char *)0;

    // TODO - above is not a current GED command - broken
    (void)ged_exec(gedp, 4, (const char **)av);
    bu_vls_free(&pt_vls);

    return BRLCAD_OK;
}


int
to_mouse_find_bot_edge(struct ged *gedp,
		       int argc,
		       const char *argv[],
		       ged_func_ptr UNUSED(func),
		       const char *usage,
		       int UNUSED(maxargs))
{
    const char *av[6];
    point_t view;
    struct bu_vls pt_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 5) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[3], "%lf", &x) != 1 ||
	bu_sscanf(argv[4], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    tclcad_mouse_sync_dimensions(gdvp);
    bv_screen_to_view(&x, &y, tclcad_mouse_bv_const(gdvp), x, y);
    VSET(view, x, y, 0.0);

    bu_vls_printf(&pt_vls, "%lf %lf %lf", view[X], view[Y], view[Z]);

    ged_view_active_ctx_set(gedp, gdvp);
    av[0] = "find_bot_edge_nearest_pnt";
    av[1] = (char *)argv[2];
    av[2] = bu_vls_addr(&pt_vls);
    av[3] = (char *)0;

    // TODO - above is not a current GED command - broken
    (void)ged_exec(gedp, 3, (const char **)av);
    bu_vls_free(&pt_vls);

    return BRLCAD_OK;
}


int
to_mouse_find_bot_pnt(struct ged *gedp,
		      int argc,
		      const char *argv[],
		      ged_func_ptr UNUSED(func),
		      const char *usage,
		      int UNUSED(maxargs))
{
    const char *av[6];
    point_t view;
    struct bu_vls pt_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 5) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[3], "%lf", &x) != 1 ||
	bu_sscanf(argv[4], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    tclcad_mouse_sync_dimensions(gdvp);
    bv_screen_to_view(&x, &y, tclcad_mouse_bv_const(gdvp), x, y);
    VSET(view, x, y, 0.0);

    bu_vls_printf(&pt_vls, "%lf %lf %lf", view[X], view[Y], view[Z]);

    ged_view_active_ctx_set(gedp, gdvp);
    av[0] = "find_bot_pnt_nearest_pnt";
    av[1] = (char *)argv[2];
    av[2] = bu_vls_addr(&pt_vls);
    av[3] = (char *)0;

    // TODO - above is not a current GED command - broken
    (void)ged_exec(gedp, 3, (const char **)av);
    bu_vls_free(&pt_vls);

    return BRLCAD_OK;
}


int
to_mouse_find_metaball_pnt(struct ged *gedp,
			   int argc,
			   const char *argv[],
			   ged_func_ptr UNUSED(func),
			   const char *usage,
			   int UNUSED(maxargs))
{
    const char *av[6];
    point_t model;
    point_t view;
    mat_t view2model;
    struct bu_vls pt_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 5) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[3], "%lf", &x) != 1 ||
	bu_sscanf(argv[4], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    tclcad_mouse_sync_dimensions(gdvp);
    bv_screen_to_view(&x, &y, tclcad_mouse_bv_const(gdvp), x, y);
    VSET(view, x, y, 0.0);
    bv_view2model_get(view2model, tclcad_mouse_bv_const(gdvp));
    MAT4X3PNT(model, view2model, view);

    bu_vls_printf(&pt_vls, "%lf %lf %lf", model[X], model[Y], model[Z]);

    ged_view_active_ctx_set(gedp, gdvp);
    av[0] = "find_metaball_pnt_nearest_pnt";
    av[1] = (char *)argv[2];
    av[2] = bu_vls_addr(&pt_vls);
    av[3] = (char *)0;

    // TODO - above is not a current GED command - broken
    (void)ged_exec(gedp, 3, (const char **)av);
    bu_vls_free(&pt_vls);

    return BRLCAD_OK;
}


int
to_mouse_find_pipe_pnt(struct ged *gedp,
		       int argc,
		       const char *argv[],
		       ged_func_ptr UNUSED(func),
		       const char *usage,
		       int UNUSED(maxargs))
{
    const char *av[6];
    point_t model;
    point_t view;
    mat_t view2model;
    struct bu_vls pt_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 5) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[3], "%lf", &x) != 1 ||
	bu_sscanf(argv[4], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    tclcad_mouse_sync_dimensions(gdvp);
    bv_screen_to_view(&x, &y, tclcad_mouse_bv_const(gdvp), x, y);
    VSET(view, x, y, 0.0);
    bv_view2model_get(view2model, tclcad_mouse_bv_const(gdvp));
    MAT4X3PNT(model, view2model, view);

    bu_vls_printf(&pt_vls, "%lf %lf %lf", model[X], model[Y], model[Z]);

    ged_view_active_ctx_set(gedp, gdvp);
    av[0] = "find_pipe_pnt_nearest_pnt";
    av[1] = (char *)argv[2];
    av[2] = bu_vls_addr(&pt_vls);
    av[3] = (char *)0;

    // TODO - above is not a current GED command - broken
    (void)ged_exec(gedp, 3, (const char **)av);
    bu_vls_free(&pt_vls);

    return BRLCAD_OK;
}


int
to_mouse_joint_select(
    struct ged *gedp,
    int argc,
    const char *argv[],
    ged_func_ptr UNUSED(func),
    const char *usage,
    int maxargs)
{
    const char *cmd_argv[11] = {"joint2", NULL, "selection", "replace", "active"};
    int ret, cmd_argc = (int)(sizeof(cmd_argv) / sizeof(const char *));
    char *joint_name;
    char *end;
    struct bu_vls bindings = BU_VLS_INIT_ZERO;
    struct bu_vls start[] = {BU_VLS_INIT_ZERO, BU_VLS_INIT_ZERO, BU_VLS_INIT_ZERO};
    struct bu_vls dir[] = {BU_VLS_INIT_ZERO, BU_VLS_INIT_ZERO, BU_VLS_INIT_ZERO};
    point_t screen_pt, view_pt, model_pt;
    vect_t view_dir, model_dir;
    mat_t invRot;
    mat_t view2model, view_rotation;

    if (argc != maxargs) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* parse args */
    joint_name = bu_path_basename(argv[2], NULL);

    screen_pt[X] = strtol(argv[3], &end, 10);
    if (*end != '\0') {
	bu_vls_printf(gedp->ged_result_str, "ERROR: bad x value %f\n", screen_pt[X]);
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    screen_pt[Y] = strtol(argv[4], &end, 10);
    if (*end != '\0') {
	bu_vls_printf(gedp->ged_result_str, "ERROR: bad y value: %f\n", screen_pt[Y]);
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    /* stash point coordinates for future drag handling */
    bv_previous_mouse_set(tclcad_mouse_bv(gdvp), screen_pt[X], screen_pt[Y]);

    /* convert screen point to model-space start point and direction */
    tclcad_mouse_sync_dimensions(gdvp);
    bv_screen_to_view(&view_pt[X], &view_pt[Y], tclcad_mouse_bv_const(gdvp),
	    screen_pt[X], screen_pt[Y]);
    view_pt[Z] = 1.0;

    bv_view2model_get(view2model, tclcad_mouse_bv_const(gdvp));
    MAT4X3PNT(model_pt, view2model, view_pt);

    VSET(view_dir, 0.0, 0.0, -1.0);
    bv_rotation_get(view_rotation, tclcad_mouse_bv_const(gdvp));
    bn_mat_inv(invRot, view_rotation);
    MAT4X3PNT(model_dir, invRot, view_dir);

    /* joint2 joint_name selection append selection_name startx starty startz dirx diry dirz */
    bu_vls_printf(&start[X], "%f", model_pt[X]);
    bu_vls_printf(&start[Y], "%f", model_pt[Y]);
    bu_vls_printf(&start[Z], "%f", model_pt[Z]);

    cmd_argv[1] = joint_name;
    cmd_argv[5] = bu_vls_addr(&start[X]);
    cmd_argv[6] = bu_vls_addr(&start[Y]);
    cmd_argv[7] = bu_vls_addr(&start[Z]);

    bu_vls_printf(&dir[X], "%f", model_dir[X]);
    bu_vls_printf(&dir[Y], "%f", model_dir[Y]);
    bu_vls_printf(&dir[Z], "%f", model_dir[Z]);

    cmd_argv[8] = bu_vls_addr(&dir[X]);
    cmd_argv[9] = bu_vls_addr(&dir[Y]);
    cmd_argv[10] = bu_vls_addr(&dir[Z]);

    ged_view_active_ctx_set(gedp, gdvp);
    ret = ged_exec_joint2(gedp, cmd_argc, cmd_argv);

    bu_vls_free(&start[X]);
    bu_vls_free(&start[Y]);
    bu_vls_free(&start[Z]);
    bu_vls_free(&dir[X]);
    bu_vls_free(&dir[Y]);
    bu_vls_free(&dir[Z]);

    if (ret != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }

    struct bu_vls *dname = tclcad_mouse_display_pathname(gdvp);
    if (dname) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_joint_selection_translate %s %s %%x %%y}",
		      bu_vls_cstr(dname),
		      bu_vls_cstr(&current_top->to_gedp->go_name),
		      tclcad_mouse_view_name(gdvp),
		      joint_name);
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    bu_free((void *)joint_name, "joint_name");

    return BRLCAD_OK;
}


int
to_mouse_joint_selection_translate(
    struct ged *gedp,
    int argc,
    const char *argv[],
    ged_func_ptr UNUSED(func),
    const char *usage,
    int maxargs)
{
    const char *cmd_argv[8] = {"joint2", NULL, "selection", "translate", "active"};
    int ret, cmd_argc = (int)(sizeof(cmd_argv) / sizeof(const char *));
    char *joint_name;
    char *end;
    point_t screen_end, view_start, view_end, model_start, model_end;
    vect_t model_delta;
    mat_t view2model;
    struct bu_vls delta[] = {BU_VLS_INIT_ZERO, BU_VLS_INIT_ZERO, BU_VLS_INIT_ZERO};

    if (argc != maxargs) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    joint_name = bu_path_basename(argv[2], NULL);

    screen_end[X] = strtol(argv[3], &end, 10);
    if (*end != '\0') {
	bu_vls_printf(gedp->ged_result_str, "ERROR: bad x value %f\n", screen_end[X]);
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    screen_end[Y] = strtol(argv[4], &end, 10);
    if (*end != '\0') {
	bu_vls_printf(gedp->ged_result_str, "ERROR: bad y value: %f\n", screen_end[Y]);
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    /* convert screen-space delta to model-space delta */
    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    (void)bv_previous_mouse_get(&prev_x, &prev_y,
	    tclcad_mouse_bv_const(gdvp));
    tclcad_mouse_sync_dimensions(gdvp);
    bv_screen_to_view(&view_start[X], &view_start[Y],
	    tclcad_mouse_bv_const(gdvp),
	    prev_x, prev_y);
    view_start[Z] = 1;
    bv_view2model_get(view2model, tclcad_mouse_bv_const(gdvp));
    MAT4X3PNT(model_start, view2model, view_start);

    tclcad_mouse_sync_dimensions(gdvp);
    bv_screen_to_view(&view_end[X], &view_end[Y],
	    tclcad_mouse_bv_const(gdvp),
	    screen_end[X], screen_end[Y]);
    view_end[Z] = 1;
    MAT4X3PNT(model_end, view2model, view_end);

    VSUB2(model_delta, model_end, model_start);

    bu_vls_printf(&delta[X], "%f", model_delta[X]);
    bu_vls_printf(&delta[Y], "%f", model_delta[Y]);
    bu_vls_printf(&delta[Z], "%f", model_delta[Z]);

    cmd_argv[1] = joint_name;
    cmd_argv[5] = bu_vls_addr(&delta[X]);
    cmd_argv[6] = bu_vls_addr(&delta[Y]);
    cmd_argv[7] = bu_vls_addr(&delta[Z]);

    ret = ged_exec_joint2(gedp, cmd_argc, cmd_argv);

    if (ret != BRLCAD_OK) {
	bu_free((void *)joint_name, "joint_name");
	bu_vls_free(&delta[X]);
	bu_vls_free(&delta[Y]);
	bu_vls_free(&delta[Z]);
	return BRLCAD_ERROR;
    }

    /* need to tell front-end that we've modified the db */
    Tcl_Eval(current_top->to_interp, "$::ArcherCore::application setSave");

    bv_previous_mouse_set(tclcad_mouse_bv(gdvp), screen_end[X], screen_end[Y]);

    cmd_argc = 3;
    cmd_argv[0] = "get";
    cmd_argv[1] = joint_name;
    cmd_argv[2] = "RP1";
    cmd_argv[3] = NULL;
    ret = ged_exec_get(gedp, cmd_argc, cmd_argv);

    if (ret == BRLCAD_OK) {
	char *path_name = bu_strdup(bu_vls_cstr(gedp->ged_result_str));
	int draw_mode = 0;
	struct bu_vls path_dmode = BU_VLS_INIT_ZERO;

	/* get current display mode of path */
	cmd_argc = 2;
	cmd_argv[0] = "how";
	cmd_argv[1] = path_name;
	cmd_argv[2] = NULL;
	ret = ged_exec_how(gedp, cmd_argc, cmd_argv);

	if (ret == BRLCAD_OK) {
	    bu_sscanf(bu_vls_cstr(gedp->ged_result_str), "%d", &draw_mode);
	}
	if (draw_mode == 4) {
	    bu_vls_printf(&path_dmode, "-h");
	} else {
	    bu_vls_printf(&path_dmode, "-m%d", draw_mode);
	}

	/* erase path to split it from visible vlists */
	cmd_argc = 2;
	cmd_argv[0] = "erase";
	cmd_argv[1] = path_name;
	cmd_argv[2] = NULL;
	ret = ged_exec_erase(gedp, cmd_argc, cmd_argv);

	if (ret == BRLCAD_OK) {
	    /* redraw path with its previous display mode */
	    cmd_argc = 4;
	    cmd_argv[0] = "draw";
	    cmd_argv[1] = "-R";
	    cmd_argv[2] = bu_vls_cstr(&path_dmode);
	    cmd_argv[3] = path_name;
	    cmd_argv[4] = NULL;
	    ret = ged_exec_draw(gedp, cmd_argc, cmd_argv);

	    to_refresh_all_views(current_top);
	}
	bu_vls_free(&path_dmode);
	bu_free(path_name, "path_name");
    }

    bu_free((void *)joint_name, "joint_name");
    bu_vls_free(&delta[X]);
    bu_vls_free(&delta[Y]);
    bu_vls_free(&delta[Z]);

    return ret;
}


int
to_mouse_move_arb_edge(struct ged *gedp,
		       int argc,
		       const char *argv[],
		       ged_func_ptr UNUSED(func),
		       const char *usage,
		       int UNUSED(maxargs))
{
    int width;
    int ret;
    const char *av[6];
    fastf_t dx, dy;
    fastf_t inv_width;
    point_t model;
    point_t view;
    mat_t inv_rot;
    struct bu_vls pt_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 6) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[4], "%lf", &x) != 1 ||
	bu_sscanf(argv[5], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    tclcad_mouse_previous_get_set(&prev_x, &prev_y, gdvp, x, y);
    dx = x - prev_x;
    dy = prev_y - y;

    tclcad_mouse_clamp_delta(&dx, &dy, gdvp);

    width = tclcad_mouse_display_width(gdvp);
    inv_width = 1.0 / (fastf_t)width;
    /* ged_move_arb_edge expects things to be in local units */
    dx *= inv_width * bv_size_get(tclcad_mouse_bv_const(gdvp)) *
	gedp->dbip->dbi_base2local;
    dy *= inv_width * bv_size_get(tclcad_mouse_bv_const(gdvp)) *
	gedp->dbip->dbi_base2local;
    VSET(view, dx, dy, 0.0);
    tclcad_mouse_view_inv_rotation(inv_rot, gdvp);
    MAT4X3PNT(model, inv_rot, view);

    bu_vls_printf(&pt_vls, "%lf %lf %lf", model[X], model[Y], model[Z]);

    ged_view_active_ctx_set(gedp, gdvp);
    av[0] = "move_arb_edge";
    av[1] = "-r";
    av[2] = (char *)argv[2];
    av[3] = (char *)argv[3];
    av[4] = bu_vls_addr(&pt_vls);
    av[5] = (char *)0;

    ret = ged_exec_move_arb_edge(gedp, 5, (const char **)av);
    bu_vls_free(&pt_vls);

    if (ret == BRLCAD_OK) {
	av[0] = "draw";
	av[1] = (char *)argv[2];
	av[2] = (char *)0;
	to_edit_redraw(gedp, 2, (const char **)av);
    }

    return BRLCAD_OK;
}


int
to_mouse_move_arb_face(struct ged *gedp,
		       int argc,
		       const char *argv[],
		       ged_func_ptr UNUSED(func),
		       const char *usage,
		       int UNUSED(maxargs))
{
    int width;
    int ret;
    const char *av[6];
    fastf_t dx, dy;
    fastf_t inv_width;
    point_t model;
    point_t view;
    mat_t inv_rot;
    struct bu_vls pt_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 6) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[4], "%lf", &x) != 1 ||
	bu_sscanf(argv[5], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    tclcad_mouse_previous_get_set(&prev_x, &prev_y, gdvp, x, y);
    dx = x - prev_x;
    dy = prev_y - y;

    tclcad_mouse_clamp_delta(&dx, &dy, gdvp);

    width = tclcad_mouse_display_width(gdvp);
    inv_width = 1.0 / (fastf_t)width;
    /* ged_move_arb_face expects things to be in local units */
    dx *= inv_width * bv_size_get(tclcad_mouse_bv_const(gdvp)) *
	gedp->dbip->dbi_base2local;
    dy *= inv_width * bv_size_get(tclcad_mouse_bv_const(gdvp)) *
	gedp->dbip->dbi_base2local;
    VSET(view, dx, dy, 0.0);
    tclcad_mouse_view_inv_rotation(inv_rot, gdvp);
    MAT4X3PNT(model, inv_rot, view);

    bu_vls_printf(&pt_vls, "%lf %lf %lf", model[X], model[Y], model[Z]);

    ged_view_active_ctx_set(gedp, gdvp);
    av[0] = "move_arb_face";
    av[1] = "-r";
    av[2] = (char *)argv[2];
    av[3] = (char *)argv[3];
    av[4] = bu_vls_addr(&pt_vls);
    av[5] = (char *)0;

    ret = ged_exec_move_arb_face(gedp, 5, (const char **)av);
    bu_vls_free(&pt_vls);

    if (ret == BRLCAD_OK) {
	av[0] = "draw";
	av[1] = (char *)argv[2];
	av[2] = (char *)0;
	to_edit_redraw(gedp, 2, (const char **)av);
    }

    return BRLCAD_OK;
}


int
to_mouse_move_bot_pnt(struct ged *gedp,
		      int argc,
		      const char *argv[],
		      ged_func_ptr UNUSED(func),
		      const char *usage,
		      int UNUSED(maxargs))
{
    int width;
    int ret;
    int rflag;
    const char *av[6];
    const char *cmd;
    fastf_t dx, dy, dz;
    fastf_t inv_width;
    point_t model;
    point_t view;
    mat_t model2view;
    mat_t v2m_mat;
    struct bu_vls pt_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    cmd = argv[0];

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", cmd, usage);
	return GED_HELP;
    }

    if (argc == 7) {
	if (argv[1][0] != '-' || argv[1][1] != 'r' || argv[1][2] != '\0') {
	    bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", cmd, usage);
	    return BRLCAD_ERROR;
	}

	rflag = 1;
	--argc;
	++argv;
    } else
	rflag = 0;

    if (argc != 6) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", cmd, usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[4], "%lf", &x) != 1 ||
	bu_sscanf(argv[5], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", cmd, usage);
	return BRLCAD_ERROR;
    }

    width = tclcad_mouse_display_width(gdvp);
    inv_width = 1.0 / (fastf_t)width;

    if (rflag) {
	fastf_t prev_x = 0.0;
	fastf_t prev_y = 0.0;
	tclcad_mouse_previous_get_set(&prev_x, &prev_y, gdvp, x, y);
	dx = x - prev_x;
	dy = prev_y - y;
	dz = 0.0;

	tclcad_mouse_clamp_delta(&dx, &dy, gdvp);

	tclcad_mouse_view_inv_rotation(v2m_mat, gdvp);

	dx *= inv_width * bv_size_get(tclcad_mouse_bv_const(gdvp));
	dy *= inv_width * bv_size_get(tclcad_mouse_bv_const(gdvp));
    } else {
	struct rt_db_internal intern;
	struct rt_bot_internal *botip;
	mat_t mat;
	size_t vertex_i;
	const char *last;

	if ((last = strrchr(argv[2], '/')) == NULL)
	    last = (char *)argv[2];
	else
	    ++last;

	if (last[0] == '\0') {
	    bu_vls_printf(gedp->ged_result_str, "%s: illegal input - %s", cmd, argv[2]);
	    return BRLCAD_ERROR;
	}

	if (bu_sscanf(argv[3], "%zu", &vertex_i) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "%s: bad bot vertex index - %s", cmd, argv[3]);
	    return BRLCAD_ERROR;
	}

	struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
	if (wdb_import_from_path2(gedp->ged_result_str, &intern, argv[2], wdbp, mat) & BRLCAD_ERROR) {
	    bu_vls_printf(gedp->ged_result_str, "%s: failed to find %s", cmd, argv[2]);
	    return BRLCAD_ERROR;
	}

	if (intern.idb_major_type != DB5_MAJORTYPE_BRLCAD ||
	    intern.idb_minor_type != DB5_MINORTYPE_BRLCAD_BOT) {
	    bu_vls_printf(gedp->ged_result_str, "Object is not a BOT");
	    rt_db_free_internal(&intern);

	    return BRLCAD_ERROR;
	}

	botip = (struct rt_bot_internal *)intern.idb_ptr;

	if (vertex_i >= botip->num_vertices) {
	    bu_vls_printf(gedp->ged_result_str, "%s: bad bot vertex index - %s", cmd, argv[3]);
	    rt_db_free_internal(&intern);
	    return BRLCAD_ERROR;
	}

	bv_model2view_get(model2view, tclcad_mouse_bv_const(gdvp));
	bv_view2model_get(v2m_mat, tclcad_mouse_bv_const(gdvp));
	MAT4X3PNT(view, model2view, &botip->vertices[vertex_i*3]);

	tclcad_mouse_sync_dimensions(gdvp);
	bv_screen_to_view(&dx, &dy, tclcad_mouse_bv_const(gdvp), x, y);
	dz = view[Z];

	rt_db_free_internal(&intern);
    }

    VSET(view, dx, dy, dz);
    MAT4X3PNT(model, v2m_mat, view);

    /* ged_bot_move_pnt expects things to be in local units */
    VSCALE(model, model, gedp->dbip->dbi_base2local);
    bu_vls_printf(&pt_vls, "%lf %lf %lf", model[X], model[Y], model[Z]);

    ged_view_active_ctx_set(gedp, gdvp);
    av[0] = "bot_move_pnt";
    // TODO - bot_move_pnt is not a current LIBGED cmd - the following are broken
    if (rflag) {
	av[1] = "-r";
	av[2] = (char *)argv[2];
	av[3] = (char *)argv[3];
	av[4] = bu_vls_addr(&pt_vls);
	av[5] = (char *)0;

	ret = ged_exec(gedp, 5, (const char **)av);
    } else {
	av[1] = (char *)argv[2];
	av[2] = (char *)argv[3];
	av[3] = bu_vls_addr(&pt_vls);
	av[4] = (char *)0;

	ret = ged_exec(gedp, 4, (const char **)av);
    }

    bu_vls_free(&pt_vls);

    if (ret == BRLCAD_OK) {
	av[0] = "draw";
	av[1] = (char *)argv[2];
	av[2] = (char *)0;
	to_edit_redraw(gedp, 2, (const char **)av);
    }

    return BRLCAD_OK;
}


int
to_mouse_move_bot_pnts(struct ged *gedp,
		       int argc,
		       const char *argv[],
		       ged_func_ptr UNUSED(func),
		       const char *usage,
		       int UNUSED(maxargs))
{
    int ret, width;
    const char *cmd;
    fastf_t dx, dy, dz;
    fastf_t inv_width;
    point_t model;
    point_t view;
    mat_t v2m_mat;
    struct bu_vls pt_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    cmd = argv[0];

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", cmd, usage);
	return GED_HELP;
    }

    if (argc < 6) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", cmd, usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[2], "%lf", &x) != 1 ||
	bu_sscanf(argv[3], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", cmd, usage);
	return BRLCAD_ERROR;
    }

    width = tclcad_mouse_display_width(gdvp);
    inv_width = 1.0 / (fastf_t)width;

    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    tclcad_mouse_previous_get_set(&prev_x, &prev_y, gdvp, x, y);
    dx = x - prev_x;
    dy = prev_y - y;
    dz = 0.0;

    tclcad_mouse_clamp_delta(&dx, &dy, gdvp);

    tclcad_mouse_view_inv_rotation(v2m_mat, gdvp);

    dx *= inv_width * bv_size_get(tclcad_mouse_bv_const(gdvp));
    dy *= inv_width * bv_size_get(tclcad_mouse_bv_const(gdvp));

    VSET(view, dx, dy, dz);
    MAT4X3PNT(model, v2m_mat, view);

    /* ged_bot_move_pnts expects things to be in local units */
    VSCALE(model, model, gedp->dbip->dbi_base2local);
    bu_vls_printf(&pt_vls, "%lf %lf %lf", model[X], model[Y], model[Z]);

    ged_view_active_ctx_set(gedp, gdvp);

    {
	register int i, j;
	int ac = argc - 2;
	const char **av = (const char **)bu_calloc(ac, sizeof(char *), "to_mouse_move_bot_pnts: av[]");
	av[0] = "bot_move_pnts";
	// TODO - above is not a current GED command - broken

	av[1] = (char *)argv[4];
	av[2] = bu_vls_addr(&pt_vls);
	av[ac-1] = (char *)0;

	for (i=3, j=5; i < ac; ++i, ++j)
	    av[i] = (char *)argv[j];

	ret = ged_exec(gedp, ac, (const char **)av);
	bu_vls_free(&pt_vls);

	if (ret == BRLCAD_OK) {
	    av[0] = "draw";
	    av[1] = (char *)argv[4];
	    av[2] = (char *)0;
	    to_edit_redraw(gedp, 2, (const char **)av);
	}

	bu_free((void *)av, "to_mouse_move_bot_pnts: av[]");
    }

    return BRLCAD_OK;
}


int
to_mouse_move_pnt_common(struct ged *gedp,
			 int argc,
			 const char *argv[],
			 ged_func_ptr func,
			 const char *usage,
			 int UNUSED(maxargs))
{
    int ret, width;
    const char *av[6];
    fastf_t dx, dy;
    fastf_t inv_width;
    point_t model;
    point_t view;
    mat_t inv_rot;
    struct bu_vls pt_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 6) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[4], "%lf", &x) != 1 ||
	bu_sscanf(argv[5], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    tclcad_mouse_previous_get_set(&prev_x, &prev_y, gdvp, x, y);
    dx = x - prev_x;
    dy = prev_y - y;

    tclcad_mouse_clamp_delta(&dx, &dy, gdvp);

    width = tclcad_mouse_display_width(gdvp);
    inv_width = 1.0 / (fastf_t)width;
    /* ged_pipe_move_pnt expects things to be in local units */
    dx *= inv_width * bv_size_get(tclcad_mouse_bv_const(gdvp)) *
	gedp->dbip->dbi_base2local;
    dy *= inv_width * bv_size_get(tclcad_mouse_bv_const(gdvp)) *
	gedp->dbip->dbi_base2local;
    VSET(view, dx, dy, 0.0);
    tclcad_mouse_view_inv_rotation(inv_rot, gdvp);
    MAT4X3PNT(model, inv_rot, view);

    bu_vls_printf(&pt_vls, "%lf %lf %lf", model[X], model[Y], model[Z]);

    ged_view_active_ctx_set(gedp, gdvp);
    av[0] = (char *)argv[0];
    av[1] = "-r";
    av[2] = (char *)argv[2];
    av[3] = (char *)argv[3];
    av[4] = bu_vls_addr(&pt_vls);
    av[5] = (char *)0;

    ret = (*func)(gedp, 5, (const char **)av);
    bu_vls_free(&pt_vls);

    if (ret == BRLCAD_OK) {
	av[0] = "draw";
	av[1] = (char *)argv[2];
	av[2] = (char *)0;
	to_edit_redraw(gedp, 2, (const char **)av);
    }

    return BRLCAD_OK;
}


int
to_mouse_orotate(struct ged *gedp,
		 int argc,
		 const char *argv[],
		 ged_func_ptr UNUSED(func),
		 const char *usage,
		 int UNUSED(maxargs))
{
    fastf_t dx, dy;
    point_t model;
    point_t view;
    mat_t inv_rot;
    struct bu_vls rot_x_vls = BU_VLS_INIT_ZERO;
    struct bu_vls rot_y_vls = BU_VLS_INIT_ZERO;
    struct bu_vls rot_z_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 5) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[3], "%lf", &x) != 1 ||
	bu_sscanf(argv[4], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    tclcad_mouse_previous_get_set(&prev_x, &prev_y, gdvp, x, y);
    dx = y - prev_y;
    dy = x - prev_x;

    tclcad_mouse_clamp_delta(&dx, &dy, gdvp);

    dx *= tclcad_mouse_rotate_scale(gdvp);
    dy *= tclcad_mouse_rotate_scale(gdvp);

    VSET(view, dx, dy, 0.0);
    tclcad_mouse_view_inv_rotation(inv_rot, gdvp);
    MAT4X3PNT(model, inv_rot, view);

    bu_vls_printf(&rot_x_vls, "%lf", model[X]);
    bu_vls_printf(&rot_y_vls, "%lf", model[Y]);
    bu_vls_printf(&rot_z_vls, "%lf", model[Z]);

    ged_view_active_ctx_set(gedp, gdvp);

    struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(gdvp);
    if (0 < bu_vls_strlen(&tvd->gdv_edit_motion_delta_callback)) {
	const char *command = bu_vls_addr(&tvd->gdv_edit_motion_delta_callback);
	const char *args[4];
	args[0] = "orotate";
	args[1] = bu_vls_addr(&rot_x_vls);
	args[2] = bu_vls_addr(&rot_y_vls);
	args[3] = bu_vls_addr(&rot_z_vls);
	tclcad_eval(current_top->to_interp, command, sizeof(args) / sizeof(args[0]), args);
    } else {
	const char *av[6];

	av[0] = "orotate";
	av[1] = (char *)argv[2];
	av[2] = bu_vls_addr(&rot_x_vls);
	av[3] = bu_vls_addr(&rot_y_vls);
	av[4] = bu_vls_addr(&rot_z_vls);
	av[5] = (char *)0;

	if (ged_exec_orotate(gedp, 5, (const char **)av) == BRLCAD_OK) {
	    av[0] = "draw";
	    av[1] = (char *)argv[2];
	    av[2] = (char *)0;
	    to_edit_redraw(gedp, 2, (const char **)av);
	}
    }

    bu_vls_free(&rot_x_vls);
    bu_vls_free(&rot_y_vls);
    bu_vls_free(&rot_z_vls);

    return BRLCAD_OK;
}


int
to_mouse_oscale(struct ged *gedp,
		int argc,
		const char *argv[],
		ged_func_ptr UNUSED(func),
		const char *usage,
		int UNUSED(maxargs))
{
    int width;
    fastf_t dx, dy;
    fastf_t sf;
    fastf_t inv_width;
    struct bu_vls sf_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 5) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[3], "%lf", &x) != 1 ||
	bu_sscanf(argv[4], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    tclcad_mouse_previous_get_set(&prev_x, &prev_y, gdvp, x, y);
    dx = x - prev_x;
    dy = prev_y - y;

    tclcad_mouse_clamp_delta(&dx, &dy, gdvp);

    width = tclcad_mouse_display_width(gdvp);
    inv_width = 1.0 / (fastf_t)width;
    dx *= inv_width * tclcad_mouse_scale_scale(gdvp);
    dy *= inv_width * tclcad_mouse_scale_scale(gdvp);

    if (fabs(dx) < fabs(dy))
	sf = 1.0 + dy;
    else
	sf = 1.0 + dx;

    bu_vls_printf(&sf_vls, "%lf", sf);

    ged_view_active_ctx_set(gedp, gdvp);

    struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(gdvp);
    if (0 < bu_vls_strlen(&tvd->gdv_edit_motion_delta_callback)) {
	struct bu_vls tcl_cmd;

	bu_vls_init(&tcl_cmd);
	bu_vls_printf(&tcl_cmd, "%s oscale %s", bu_vls_addr(&tvd->gdv_edit_motion_delta_callback), bu_vls_addr(&sf_vls));
	Tcl_Eval(current_top->to_interp, bu_vls_addr(&tcl_cmd));
	bu_vls_free(&tcl_cmd);
    } else {
	const char *av[6];

	av[0] = "oscale";
	av[1] = (char *)argv[2];
	av[2] = bu_vls_addr(&sf_vls);
	av[3] = (char *)0;

	if (ged_exec_oscale(gedp, 3, (const char **)av) == BRLCAD_OK) {
	    av[0] = "draw";
	    av[1] = (char *)argv[2];
	    av[2] = (char *)0;
	    to_edit_redraw(gedp, 2, (const char **)av);
	}
    }

    bu_vls_free(&sf_vls);

    return BRLCAD_OK;
}


int
to_mouse_otranslate(struct ged *gedp,
		    int argc,
		    const char *argv[],
		    ged_func_ptr UNUSED(func),
		    const char *usage,
		    int UNUSED(maxargs))
{
    int width;
    fastf_t dx, dy;
    fastf_t inv_width;
    point_t model = VINIT_ZERO;
    point_t view = VINIT_ZERO;
    mat_t inv_rot;
    struct bu_vls tran_x_vls = BU_VLS_INIT_ZERO;
    struct bu_vls tran_y_vls = BU_VLS_INIT_ZERO;
    struct bu_vls tran_z_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 5) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[3], "%lf", &x) != 1 ||
	bu_sscanf(argv[4], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    tclcad_mouse_previous_get_set(&prev_x, &prev_y, gdvp, x, y);
    dx = x - prev_x;
    dy = prev_y - y;

    tclcad_mouse_clamp_delta(&dx, &dy, gdvp);

    width = tclcad_mouse_display_width(gdvp);
    inv_width = 1.0 / (fastf_t)width;
    /* ged_otranslate expects things to be in local units */
    dx *= inv_width * bv_size_get(tclcad_mouse_bv_const(gdvp)) *
	gedp->dbip->dbi_base2local;
    dy *= inv_width * bv_size_get(tclcad_mouse_bv_const(gdvp)) *
	gedp->dbip->dbi_base2local;

    VSET(view, dx, dy, 0.0);
    bu_vls_printf(&tran_x_vls, "%lf", model[X]);
    bu_vls_printf(&tran_y_vls, "%lf", model[Y]);
    bu_vls_printf(&tran_z_vls, "%lf", model[Z]);

    ged_view_active_ctx_set(gedp, gdvp);

    struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(gdvp);
    struct tclcad_ged_data *tgd = (struct tclcad_ged_data *)current_top->to_gedp->u_data;
    if (0 < bu_vls_strlen(&tvd->gdv_edit_motion_delta_callback)) {
	const char *path_string = argv[2];
	vect_t dvec;
	struct dm_path_edit_params *params = (struct dm_path_edit_params *)bu_hash_get(tgd->go_dmv.edited_paths,
										 (uint8_t *)path_string,
										 sizeof(char) * strlen(path_string) + 1);

	if (!params) {
	    BU_GET(params, struct dm_path_edit_params);
	    params->edit_mode = tclcad_view_polygon_mode_from_view_ctx(gdvp);
	    params->dx = params->dy = 0.0;
	    (void)bu_hash_set(tgd->go_dmv.edited_paths,
			      (uint8_t *)path_string,
			      sizeof(char) * strlen(path_string) + 1, (void *)params);
	}

	params->dx += dx;
	params->dy += dy;
	VSET(view, params->dx, params->dy, 0.0);
	tclcad_mouse_view_inv_rotation(inv_rot, gdvp);
	MAT4X3PNT(model, inv_rot, view);

	MAT_IDN(params->edit_mat);
	MAT4X3PNT(model, inv_rot, view);
	VSCALE(dvec, model, gedp->dbip->dbi_local2base);
	MAT_DELTAS_VEC(params->edit_mat, dvec);

	to_refresh_view(gdvp);
    } else {
	const char *av[6];

	av[0] = "otranslate";
	av[1] = (char *)argv[2];
	av[2] = bu_vls_addr(&tran_x_vls);
	av[3] = bu_vls_addr(&tran_y_vls);
	av[4] = bu_vls_addr(&tran_z_vls);
	av[5] = (char *)0;

	if (ged_exec_otranslate(gedp, 5, (const char **)av) == BRLCAD_OK) {
	    av[0] = "draw";
	    av[1] = (char *)argv[2];
	    av[2] = (char *)0;
	    to_edit_redraw(gedp, 2, (const char **)av);
	}
    }

    bu_vls_free(&tran_x_vls);
    bu_vls_free(&tran_y_vls);
    bu_vls_free(&tran_z_vls);

    return BRLCAD_OK;
}


int
go_mouse_poly_circ(Tcl_Interp *interp,
		   struct ged *gedp,
		   void *draw_view_ctx,
		   int argc,
		   const char *argv[],
		   const char *usage)
{
    void *view_ctx = draw_view_ctx;
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }
    to_refresh_suppress_all_begin(current_top);

    int ret = to_mouse_poly_circ_func(interp, gedp, view_ctx, argc, argv, usage);
    to_refresh_suppress_all_end(current_top);
    return ret;
}


int
to_mouse_poly_circ(struct ged *gedp,
		   int argc,
		   const char *argv[],
		   ged_func_ptr UNUSED(func),
		   const char *usage,
		   int UNUSED(maxargs))
{
    int ret;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *view_ctx = ged_view_find_ctx(gedp, argv[1]);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* shift the command name to argv[1] before calling to_mouse_poly_circ_func */
    argv[1] = argv[0];
    ret = to_mouse_poly_circ_func(current_top->to_interp, gedp, view_ctx, argc-1, argv+1, usage);
#if 0
    if (ret == BRLCAD_ERROR)
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
#endif

    to_refresh_view(view_ctx);

    return ret;
}


int
to_mouse_poly_circ_func(Tcl_Interp *interp,
			struct ged *gedp,
			void *view_ctx,
			int UNUSED(argc),
			const char *argv[],
			const char *usage)
{
    int ac;
    const char *av[5];
    int x, y;
    fastf_t fx, fy;
    point_t v_pt, m_pt;
    mat_t view2model;
    struct bu_vls plist = BU_VLS_INIT_ZERO;
    struct bu_vls i_vls = BU_VLS_INIT_ZERO;
    tclcad_polygon_state *gdpsp =
	tclcad_view_polygon_state_from_view_ctx(view_ctx, argv[0][0] == 's');
    if (!gdpsp)
	return BRLCAD_ERROR;

    if (bu_sscanf(argv[1], "%d", &x) != 1 ||
	bu_sscanf(argv[2], "%d", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_mouse_bv(view_ctx), x, y);

    tclcad_mouse_sync_dimensions(view_ctx);
    bv_screen_to_view(&fx, &fy, tclcad_mouse_bv_const(view_ctx), x, y);
    bv_view2model_get(view2model, tclcad_mouse_bv_const(view_ctx));

    {
	unsigned long long snap_kinds = tclcad_mouse_prepare_snap(view_ctx);
	if (snap_kinds)
	    tclcad_mouse_snap_point_2d(view_ctx, &fx, &fy, snap_kinds);
    }

    bu_vls_printf(&plist, "{0 ");

    {
	vect_t vdiff;
	fastf_t r, arc;
	fastf_t curr_fx, curr_fy;
	register int nsegs, n;

	VSET(v_pt, fx, fy, gdpsp->gdps_data_vZ);
	VSUB2(vdiff, v_pt, gdpsp->gdps_prev_point);
	r = MAGNITUDE(vdiff);

	/* use a variable number of segments based on the size of the
	 * circle being created so small circles have few segments and
	 * large ones are nice and smooth.  select a chord length that
	 * results in segments approximately 4 pixels in length.
	 *
	 * circumference / 4 = PI * diameter / 4
	 *
	 */
	nsegs = M_PI_2 * r * bv_scale_get(tclcad_mouse_bv_const(view_ctx));

	if (nsegs < 32)
	    nsegs = 32;

	arc = 360.0 / nsegs;
	for (n = 0; n < nsegs; ++n) {
	    fastf_t ang = n * arc;

	    curr_fx = cos(ang*DEG2RAD) * r + gdpsp->gdps_prev_point[X];
	    curr_fy = sin(ang*DEG2RAD) * r + gdpsp->gdps_prev_point[Y];
	    VSET(v_pt, curr_fx, curr_fy, gdpsp->gdps_data_vZ);
	    MAT4X3PNT(m_pt, view2model, v_pt);
	    bu_vls_printf(&plist, " {%lf %lf %lf}", V3ARGS(m_pt));
	}
    }

    bu_vls_printf(&plist, " }");
    bu_vls_printf(&i_vls, "%zu", gdpsp->gdps_curr_polygon_i);

    ged_view_active_ctx_set(gedp, view_ctx);
    ac = 4;
    av[0] = "data_polygons";
    av[1] = "replace_poly";
    av[2] = bu_vls_addr(&i_vls);
    av[3] = bu_vls_addr(&plist);
    av[4] = (char *)0;

    (void)to_data_polygons_func(interp, gedp, view_ctx, ac, (const char **)av);
    bu_vls_free(&plist);
    bu_vls_free(&i_vls);

    return BRLCAD_OK;
}


int
go_mouse_poly_cont(Tcl_Interp *interp,
		   struct ged *gedp,
		   void *draw_view_ctx,
		   int argc,
		   const char *argv[],
		   const char *usage)
{
    void *view_ctx = draw_view_ctx;
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }
    to_refresh_suppress_all_begin(current_top);

    int ret = to_mouse_poly_cont_func(interp, gedp, view_ctx, argc, argv, usage);
    to_refresh_suppress_all_end(current_top);
    return ret;
}


int
to_mouse_poly_cont(struct ged *gedp,
		   int argc,
		   const char *argv[],
		   ged_func_ptr UNUSED(func),
		   const char *usage,
		   int UNUSED(maxargs))
{
    int ret;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *view_ctx = ged_view_find_ctx(gedp, argv[1]);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* shift the command name to argv[1] before calling to_mouse_poly_cont_func */
    argv[1] = argv[0];
    ret = to_mouse_poly_cont_func(current_top->to_interp, gedp, view_ctx, argc-1, argv+1, usage);
#if 0
    if (ret == BRLCAD_ERROR)
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
#endif

    to_refresh_view(view_ctx);

    return ret;
}


int
to_mouse_poly_cont_func(Tcl_Interp *interp,
			struct ged *gedp,
			void *view_ctx,
			int UNUSED(argc),
			const char *argv[],
			const char *usage)
{
    int ac;
    const char *av[7];
    int x, y;
    fastf_t fx, fy;
    point_t v_pt, m_pt;
    mat_t view2model;
    tclcad_polygon_state *gdpsp =
	tclcad_view_polygon_state_from_view_ctx(view_ctx, argv[0][0] == 's');
    if (!gdpsp)
	return BRLCAD_ERROR;

    if (bu_sscanf(argv[1], "%d", &x) != 1 ||
	bu_sscanf(argv[2], "%d", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_mouse_bv(view_ctx), x, y);

    tclcad_mouse_sync_dimensions(view_ctx);
    bv_screen_to_view(&fx, &fy, tclcad_mouse_bv_const(view_ctx), x, y);
    VSET(v_pt, fx, fy, gdpsp->gdps_data_vZ);

    bv_view2model_get(view2model, tclcad_mouse_bv_const(view_ctx));
    MAT4X3PNT(m_pt, view2model, v_pt);
    ged_view_active_ctx_set(gedp, view_ctx);

    {
	struct bu_vls i_vls = BU_VLS_INIT_ZERO;
	struct bu_vls k_vls = BU_VLS_INIT_ZERO;
	struct bu_vls plist = BU_VLS_INIT_ZERO;

	bu_vls_printf(&i_vls, "%zu", gdpsp->gdps_curr_polygon_i);
	bu_vls_printf(&k_vls, "%zu", gdpsp->gdps_curr_point_i);
	bu_vls_printf(&plist, "%lf %lf %lf", V3ARGS(m_pt));

	ac = 6;
	av[0] = "data_polygons";
	av[1] = "replace_point";
	av[2] = bu_vls_addr(&i_vls);
	av[3] = "0";
	av[4] = bu_vls_addr(&k_vls);
	av[5] = bu_vls_addr(&plist);
	av[6] = (char *)0;

	(void)to_data_polygons_func(interp, gedp, view_ctx, ac, (const char **)av);
	bu_vls_free(&i_vls);
	bu_vls_free(&k_vls);
	bu_vls_free(&plist);
    }

    return BRLCAD_OK;
}


int
go_mouse_poly_ell(Tcl_Interp *interp,
		  struct ged *gedp,
		  void *draw_view_ctx,
		  int argc,
		  const char *argv[],
		  const char *usage)
{
    void *view_ctx = draw_view_ctx;
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }
    to_refresh_suppress_all_begin(current_top);

    int ret = to_mouse_poly_ell_func(interp, gedp, view_ctx, argc, argv, usage);
    to_refresh_suppress_all_end(current_top);
    return ret;
}


int
to_mouse_poly_ell(struct ged *gedp,
		  int argc,
		  const char *argv[],
		  ged_func_ptr UNUSED(func),
		  const char *usage,
		  int UNUSED(maxargs))
{
    int ret;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *view_ctx = ged_view_find_ctx(gedp, argv[1]);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* shift the command name to argv[1] before calling to_mouse_poly_ell_func */
    argv[1] = argv[0];
    ret = to_mouse_poly_ell_func(current_top->to_interp, gedp, view_ctx, argc-1, argv+1, usage);
#if 0
    if (ret == BRLCAD_ERROR)
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
#endif

    to_refresh_view(view_ctx);

    return ret;
}


int
to_mouse_poly_ell_func(Tcl_Interp *interp,
		       struct ged *gedp,
		       void *view_ctx,
		       int UNUSED(argc),
		       const char *argv[],
		       const char *usage)
{
    int ac;
    const char *av[5];
    int x, y;
    fastf_t fx, fy;
    point_t m_pt;
    mat_t view2model;
    struct bu_vls plist = BU_VLS_INIT_ZERO;
    struct bu_vls i_vls = BU_VLS_INIT_ZERO;
    tclcad_polygon_state *gdpsp =
	tclcad_view_polygon_state_from_view_ctx(view_ctx, argv[0][0] == 's');
    if (!gdpsp)
	return BRLCAD_ERROR;

    if (bu_sscanf(argv[1], "%d", &x) != 1 ||
	bu_sscanf(argv[2], "%d", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_mouse_bv(view_ctx), x, y);


    tclcad_mouse_sync_dimensions(view_ctx);
    bv_screen_to_view(&fx, &fy, tclcad_mouse_bv_const(view_ctx), x, y);
    bv_view2model_get(view2model, tclcad_mouse_bv_const(view_ctx));

    {
	unsigned long long snap_kinds = tclcad_mouse_prepare_snap(view_ctx);
	if (snap_kinds)
	    tclcad_mouse_snap_point_2d(view_ctx, &fx, &fy, snap_kinds);
    }

    bu_vls_printf(&plist, "{0 ");

    {
	fastf_t a, b, arc;
	point_t ellout;
	point_t A, B;
	register int nsegs, n;

	a = fx - gdpsp->gdps_prev_point[X];
	b = fy - gdpsp->gdps_prev_point[Y];

	/*
	 * For angle alpha, compute surface point as
	 *
	 * V + cos(alpha) * A + sin(alpha) * B
	 *
	 * note that sin(alpha) is cos(90-alpha).
	 */

	VSET(A, a, 0, gdpsp->gdps_data_vZ);
	VSET(B, 0, b, gdpsp->gdps_data_vZ);

	/* use a variable number of segments based on the size of the
	 * circle being created so small circles have few segments and
	 * large ones are nice and smooth.  select a chord length that
	 * results in segments approximately 4 pixels in length.
	 *
	 * circumference / 4 = PI * diameter / 4
	 *
	 */
	nsegs = M_PI_2 * FMAX(a, b) *
	    bv_scale_get(tclcad_mouse_bv_const(view_ctx));

	if (nsegs < 32)
	    nsegs = 32;

	arc = 360.0 / nsegs;
	for (n = 0; n < nsegs; ++n) {
	    fastf_t cosa = cos(n * arc * DEG2RAD);
	    fastf_t sina = sin(n * arc * DEG2RAD);

	    VJOIN2(ellout, gdpsp->gdps_prev_point, cosa, A, sina, B);
	    MAT4X3PNT(m_pt, view2model, ellout);
	    bu_vls_printf(&plist, " {%lf %lf %lf}", V3ARGS(m_pt));
	}
    }

    bu_vls_printf(&plist, " }");
    bu_vls_printf(&i_vls, "%zu", gdpsp->gdps_curr_polygon_i);

    ged_view_active_ctx_set(gedp, view_ctx);
    ac = 4;
    av[0] = "data_polygons";
    av[1] = "replace_poly";
    av[2] = bu_vls_addr(&i_vls);
    av[3] = bu_vls_addr(&plist);
    av[4] = (char *)0;

    (void)to_data_polygons_func(interp, gedp, view_ctx, ac, (const char **)av);
    bu_vls_free(&plist);
    bu_vls_free(&i_vls);

    return BRLCAD_OK;
}


int
go_mouse_poly_rect(Tcl_Interp *interp,
		   struct ged *gedp,
		   void *draw_view_ctx,
		   int argc,
		   const char *argv[],
		   const char *usage)
{
    void *view_ctx = draw_view_ctx;
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }
    to_refresh_suppress_all_begin(current_top);

    int ret = to_mouse_poly_rect_func(interp, gedp, view_ctx, argc, argv, usage);
    to_refresh_suppress_all_end(current_top);
    return ret;
}


int
to_mouse_poly_rect(struct ged *gedp,
		   int argc,
		   const char *argv[],
		   ged_func_ptr UNUSED(func),
		   const char *usage,
		   int UNUSED(maxargs))
{
    int ret;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *view_ctx = ged_view_find_ctx(gedp, argv[1]);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* shift the command name to argv[1] before calling to_mouse_poly_rect_func */
    argv[1] = argv[0];
    ret = to_mouse_poly_rect_func(current_top->to_interp, gedp, view_ctx, argc-1, argv+1, usage);
#if 0
    if (ret == BRLCAD_ERROR)
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
#endif

    to_refresh_view(view_ctx);

    return ret;
}


int
to_mouse_poly_rect_func(Tcl_Interp *interp,
			struct ged *gedp,
			void *view_ctx,
			int UNUSED(argc),
			const char *argv[],
			const char *usage)
{
    int ac;
    const char *av[5];
    int x, y;
    fastf_t fx, fy;
    point_t v_pt, m_pt;
    mat_t view2model;
    struct bu_vls plist = BU_VLS_INIT_ZERO;
    struct bu_vls i_vls = BU_VLS_INIT_ZERO;
    tclcad_polygon_state *gdpsp =
	tclcad_view_polygon_state_from_view_ctx(view_ctx, argv[0][0] == 's');
    if (!gdpsp)
	return BRLCAD_ERROR;

    if (bu_sscanf(argv[1], "%d", &x) != 1 ||
	bu_sscanf(argv[2], "%d", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_mouse_bv(view_ctx), x, y);

    tclcad_mouse_sync_dimensions(view_ctx);
    bv_screen_to_view(&fx, &fy, tclcad_mouse_bv_const(view_ctx), x, y);
    bv_view2model_get(view2model, tclcad_mouse_bv_const(view_ctx));

    {
	unsigned long long snap_kinds = tclcad_mouse_prepare_snap(view_ctx);
	if (snap_kinds)
	    tclcad_mouse_snap_point_2d(view_ctx, &fx, &fy, snap_kinds);
    }


    if (tclcad_view_polygon_mode_from_view_ctx(view_ctx) == TCLCAD_POLY_SQUARE_MODE) {
	fastf_t dx, dy;

	dx = fx - gdpsp->gdps_prev_point[X];
	dy = fy - gdpsp->gdps_prev_point[Y];

	if (fabs(dx) > fabs(dy)) {
	    if (dy < 0.0)
		fy = gdpsp->gdps_prev_point[Y] - fabs(dx);
	    else
		fy = gdpsp->gdps_prev_point[Y] + fabs(dx);
	} else {
	    if (dx < 0.0)
		fx = gdpsp->gdps_prev_point[X] - fabs(dy);
	    else
		fx = gdpsp->gdps_prev_point[X] + fabs(dy);
	}
    }

    MAT4X3PNT(m_pt, view2model, gdpsp->gdps_prev_point);
    bu_vls_printf(&plist, "{0 {%lf %lf %lf} ",  V3ARGS(m_pt));

    VSET(v_pt, gdpsp->gdps_prev_point[X], fy, gdpsp->gdps_data_vZ);
    MAT4X3PNT(m_pt, view2model, v_pt);
    bu_vls_printf(&plist, "{%lf %lf %lf} ",  V3ARGS(m_pt));

    VSET(v_pt, fx, fy, gdpsp->gdps_data_vZ);
    MAT4X3PNT(m_pt, view2model, v_pt);
    bu_vls_printf(&plist, "{%lf %lf %lf} ",  V3ARGS(m_pt));
    VSET(v_pt, fx, gdpsp->gdps_prev_point[Y], gdpsp->gdps_data_vZ);
    MAT4X3PNT(m_pt, view2model, v_pt);
    bu_vls_printf(&plist, "{%lf %lf %lf} }",  V3ARGS(m_pt));

    bu_vls_printf(&i_vls, "%zu", gdpsp->gdps_curr_polygon_i);

    ged_view_active_ctx_set(gedp, view_ctx);
    ac = 4;
    av[0] = "data_polygons";
    av[1] = "replace_poly";
    av[2] = bu_vls_addr(&i_vls);
    av[3] = bu_vls_addr(&plist);
    av[4] = (char *)0;

    (void)to_data_polygons_func(interp, gedp, view_ctx, ac, (const char **)av);
    bu_vls_free(&plist);
    bu_vls_free(&i_vls);

    return BRLCAD_OK;
}


int
to_mouse_ray(struct ged *gedp,
	     int argc,
	     const char *argv[],
	     ged_func_ptr UNUSED(func),
	     const char *usage,
	     int UNUSED(maxargs))
{
    bu_vls_trunc(gedp->ged_result_str, 0);
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }
    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *view_ctx = ged_view_find_ctx(gedp, argv[1]);
    double screen_x = 0.0;
    double screen_y = 0.0;
    if (!view_ctx || bu_sscanf(argv[2], "%lf", &screen_x) != 1 ||
	bu_sscanf(argv[3], "%lf", &screen_y) != 1)
	return BRLCAD_ERROR;

    tclcad_mouse_sync_dimensions(view_ctx);
    bv_screen_to_view(&screen_x, &screen_y,
	tclcad_mouse_bv_const(view_ctx), screen_x, screen_y);
    struct bv_grid_state grid = BV_GRID_STATE_INIT;
    if (bv_grid_state_get(&grid, tclcad_mouse_bv_const(view_ctx)) &&
	grid.snap)
	(void)bv_snap_grid_2d(tclcad_mouse_bv_const(view_ctx),
	    &screen_x, &screen_y);
    mat_t view2model;
    bv_view2model_get(view2model, tclcad_mouse_bv_const(view_ctx));
    point_t view_point;
    point_t ray_start;
    point_t ray_target;
    VSET(view_point, screen_x, screen_y, BV_VIEW_MAX);
    MAT4X3PNT(ray_start, view2model, view_point);
    VSET(view_point, screen_x, screen_y, 0.0);
    MAT4X3PNT(ray_target, view2model, view_point);
    bu_vls_printf(gedp->ged_result_str,
	"{%0.17g %0.17g %0.17g} {%0.17g %0.17g %0.17g}",
	V3ARGS(ray_start), V3ARGS(ray_target));
    return BRLCAD_OK;
}

static int
tclcad_pick_path_matches(const char *path, const char *target)
{
    if (!target || !target[0])
	return 1;
    while (path && path[0] == '/')
	path++;
    while (target[0] == '/')
	target++;
    if (!path)
	return 0;
    if (BU_STR_EQUAL(path, target))
	return 1;
    const size_t path_len = strlen(path);
    const size_t target_len = strlen(target);
    return path_len > target_len &&
	path[path_len - target_len - 1] == '/' &&
	BU_STR_EQUAL(path + path_len - target_len, target);
}

int
to_mouse_pick_detail(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    bu_vls_trunc(gedp->ged_result_str, 0);
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }
    if (argc != 4 && argc != 5) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *view_ctx = ged_view_find_ctx(gedp, argv[1]);
    int x = 0;
    int y = 0;
    if (!view_ctx || bu_sscanf(argv[2], "%d", &x) != 1 ||
	bu_sscanf(argv[3], "%d", &y) != 1)
	return BRLCAD_ERROR;

    const char *target = argc == 5 ? argv[4] : NULL;
    struct ged_draw_pick_result *result =
	ged_draw_view_context_pick_point(view_ctx, x, y, target ? 0 : 1);
    const size_t count = ged_draw_pick_result_count(result);
    for (size_t i = 0; i < count; i++) {
	struct bu_vls path = BU_VLS_INIT_ZERO;
	struct ged_draw_pick_detail detail = GED_DRAW_PICK_DETAIL_INIT;
	if (!ged_draw_pick_result_path(result, i, &path) ||
	    !tclcad_pick_path_matches(bu_vls_cstr(&path), target) ||
	    !ged_draw_pick_result_detail(result, i, &detail)) {
	    bu_vls_free(&path);
	    continue;
	}
	bu_vls_printf(gedp->ged_result_str,
	    "path {%s} hit_distance %.17g primitive_kind %d "
	    "primitive_index %d source_id %u material_id %d "
	    "face_vertices {%d %d %d} nearest_face_vertex %d",
	    bu_vls_cstr(&path), ged_draw_pick_result_hit_dist(result, i),
	    detail.primitive_kind, detail.primitive_index, detail.source_id,
	    detail.material_id, detail.face_vertex_index[0],
	    detail.face_vertex_index[1], detail.face_vertex_index[2],
	    detail.nearest_face_vertex_index);
	if (detail.model_point_valid)
	    bu_vls_printf(gedp->ged_result_str, " model_point {%0.17g %0.17g %0.17g}",
		V3ARGS(detail.model_point));
	bu_vls_free(&path);
	break;
    }
    ged_draw_pick_result_free(result);
    return BRLCAD_OK;
}


int
to_mouse_rect(struct ged *gedp,
	      int argc,
	      const char *argv[],
	      ged_func_ptr UNUSED(func),
	      const char *usage,
	      int UNUSED(maxargs))
{
    int ret;
    int ac;
    const char *av[5];
    int x, y;
    int dx, dy;
    struct bu_vls dx_vls = BU_VLS_INIT_ZERO;
    struct bu_vls dy_vls = BU_VLS_INIT_ZERO;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[2], "%d", &x) != 1 ||
	bu_sscanf(argv[3], "%d", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    (void)bv_previous_mouse_get(&prev_x, &prev_y,
	    tclcad_mouse_bv_const(gdvp));
    dx = x - prev_x;
    dy = tclcad_mouse_display_height(gdvp) - y - prev_y;

    bu_vls_printf(&dx_vls, "%d", dx);
    bu_vls_printf(&dy_vls, "%d", dy);
    ged_view_active_ctx_set(gedp, gdvp);
    ac = 4;
    av[0] = "rect";
    av[1] = "dim";
    av[2] = bu_vls_addr(&dx_vls);
    av[3] = bu_vls_addr(&dy_vls);
    av[4] = (char *)0;

    ret = ged_exec_rect(gedp, ac, (const char **)av);
    bu_vls_free(&dx_vls);
    bu_vls_free(&dy_vls);

    if (ret == BRLCAD_OK)
	to_refresh_view(gdvp);

    return BRLCAD_OK;
}


int
to_mouse_rot(struct ged *gedp,
	     int argc,
	     const char *argv[],
	     ged_func_ptr UNUSED(func),
	     const char *usage,
	     int UNUSED(maxargs))
{
    int ret;
    int ac;
    const char *av[4];
    fastf_t dx, dy;
    struct bu_vls rot_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[2], "%lf", &x) != 1 ||
	bu_sscanf(argv[3], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    tclcad_mouse_previous_get_set(&prev_x, &prev_y, gdvp, x, y);
    dx = prev_y - y;
    dy = prev_x - x;

    tclcad_mouse_clamp_delta(&dx, &dy, gdvp);

    dx *= tclcad_mouse_rotate_scale(gdvp);
    dy *= tclcad_mouse_rotate_scale(gdvp);

    bu_vls_printf(&rot_vls, "%lf %lf 0", dx, dy);

    ged_view_active_ctx_set(gedp, gdvp);
    ac = 3;
    av[0] = "rot";
    av[1] = "-v";
    av[2] = bu_vls_addr(&rot_vls);
    av[3] = (char *)0;

    ret = ged_exec_rot(gedp, ac, (const char **)av);
    bu_vls_free(&rot_vls);

    if (ret == BRLCAD_OK) {
	struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(gdvp);
	if (0 < bu_vls_strlen(&tvd->gdv_callback)) {
	    Tcl_Eval(current_top->to_interp, bu_vls_addr(&tvd->gdv_callback));
	}

	to_refresh_view(gdvp);
    }

    return BRLCAD_OK;
}


int
to_mouse_rotate_arb_face(struct ged *gedp,
			 int argc,
			 const char *argv[],
			 ged_func_ptr UNUSED(func),
			 const char *usage,
			 int UNUSED(maxargs))
{
    int ret;
    const char *av[6];
    fastf_t dx, dy;
    point_t model;
    point_t view;
    mat_t inv_rot;
    struct bu_vls pt_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 7) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[5], "%lf", &x) != 1 ||
	bu_sscanf(argv[6], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    tclcad_mouse_previous_get_set(&prev_x, &prev_y, gdvp, x, y);
    dx = y - prev_y;
    dy = x - prev_x;

    tclcad_mouse_clamp_delta(&dx, &dy, gdvp);

    dx *= tclcad_mouse_rotate_scale(gdvp);
    dy *= tclcad_mouse_rotate_scale(gdvp);

    VSET(view, dx, dy, 0.0);
    tclcad_mouse_view_inv_rotation(inv_rot, gdvp);
    MAT4X3PNT(model, inv_rot, view);

    bu_vls_printf(&pt_vls, "%lf %lf %lf", model[X], model[Y], model[Z]);

    ged_view_active_ctx_set(gedp, gdvp);
    av[0] = "rotate_arb_face";
    av[1] = (char *)argv[2];
    av[2] = (char *)argv[3];
    av[3] = (char *)argv[4];
    av[4] = bu_vls_addr(&pt_vls);
    av[5] = (char *)0;

    ret = ged_exec_rotate_arb_face(gedp, 5, (const char **)av);
    bu_vls_free(&pt_vls);

    if (ret == BRLCAD_OK) {
	av[0] = "draw";
	av[1] = (char *)argv[2];
	av[2] = (char *)0;
	to_edit_redraw(gedp, 2, (const char **)av);
    }

    return BRLCAD_OK;
}


#define TO_COMMON_MOUSE_SCALE(_gdvp, _zoom_vls, _argc, _argv, _usage) { \
	int _width; \
	fastf_t _dx, _dy; \
	fastf_t _inv_width; \
	fastf_t _prev_x, _prev_y; \
	fastf_t _sf; \
 \
	/* must be double for scanf */ \
	double _x, _y; \
 \
	/* initialize result */ \
	bu_vls_trunc(gedp->ged_result_str, 0); \
 \
	/* must be wanting help */ \
	if ((_argc) == 1) { \
	    bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", (_argv)[0], (_usage)); \
	    return GED_HELP; \
	} \
 \
	if ((_argc) != 4) { \
	    bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", (_argv)[0], (_usage)); \
	    return BRLCAD_ERROR; \
	} \
 \
        (_gdvp) = ged_view_find_ctx(gedp, (_argv)[1]); \
        if (!(_gdvp)) { \
	    bu_vls_printf(gedp->ged_result_str, "View not found - %s", (_argv)[1]); \
	    return BRLCAD_ERROR; \
	} \
 \
	if (bu_sscanf((_argv)[2], "%lf", &_x) != 1 || \
	    bu_sscanf((_argv)[3], "%lf", &_y) != 1) { \
	    bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", (_argv)[0], (_usage)); \
	    return BRLCAD_ERROR; \
	} \
 \
	(void)bv_previous_mouse_get(&_prev_x, &_prev_y, \
		tclcad_mouse_bv_const((_gdvp))); \
	_dx = _x - _prev_x; \
	_dy = _prev_y - _y; \
 \
	bv_previous_mouse_set(tclcad_mouse_bv((_gdvp)), _x, _y); \
 \
	tclcad_mouse_clamp_delta(&_dx, &_dy, (_gdvp)); \
 \
	_width = tclcad_mouse_display_width((_gdvp)); \
	_inv_width = 1.0 / (fastf_t)_width; \
	_dx *= _inv_width * tclcad_mouse_scale_scale((_gdvp)); \
	_dy *= _inv_width * tclcad_mouse_scale_scale((_gdvp)); \
 \
	if (fabs(_dx) > fabs(_dy)) \
	    _sf = 1.0 + _dx; \
	else \
	    _sf = 1.0 + _dy; \
 \
	bu_vls_printf(&(_zoom_vls), "%lf", _sf);	\
    }

/*
 * Usage: data_scale vname dtype sf
 */
static int
to_data_scale(struct ged *gedp,
	      int argc,
	      const char *argv[],
	      ged_func_ptr UNUSED(func),
	      const char *usage,
	      int UNUSED(maxargs))
{
    register int i;
    fastf_t sf;
    mat_t model2view, view2model;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    usage = "vname dtype sf";
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[2], "%lf", &sf) != 1 || sf < 0) {
	bu_vls_printf(gedp->ged_result_str, "Invalid scale factor - %s", argv[2]);
	return BRLCAD_ERROR;
    }
    bv_model2view_get(model2view, tclcad_mouse_bv_const(gdvp));
    bv_view2model_get(view2model, tclcad_mouse_bv_const(gdvp));

    /* scale data arrows through typed data-arrow facades */
    {
	const char *feature_name = "_tcl_data_arrows";
	point_t vcenter = {0, 0, 0};
	point_t *_pts = NULL;
	int _npts = _tclcad_draw_view_data_arrows_points_copy(gdvp, feature_name, &_pts);

	/* Arrows are stored as MOVE/DRAW pairs; need even count. */
	if (_npts >= 2 && (_npts % 2) == 0) {
	    /* Scale the length of each arrow (even-indexed endpoints = shaft starts) */
	    for (i = 0; i < _npts; i += 2) {
		vect_t diff;
		point_t vpoint;
		MAT4X3PNT(vpoint, model2view, _pts[i]);
		vcenter[Z] = vpoint[Z];
		VSUB2(diff, vpoint, vcenter);
		VSCALE(diff, diff, sf);
		VADD2(vpoint, vcenter, diff);
		MAT4X3PNT(_pts[i], view2model, vpoint);
	    }

	    int _color[3]; int _lw, _tl, _tw, _vis;
	    _tclcad_draw_view_data_arrows_style_read(gdvp, feature_name, _color, &_lw, &_tl, &_tw, &_vis);
	    _tclcad_draw_view_data_arrows_replace(gdvp, feature_name, _pts, _npts,
			_color, _lw, _tl, _tw, _vis);
	}
	if (_pts)
	    bu_free(_pts, "TclCAD draw-view points");
    }

    /* scale data labels through typed data-label facades */
    {
	const char *label_name = "_tcl_data_labels";
	size_t _child_cnt = _tclcad_draw_view_data_labels_count(gdvp, label_name);
	if (_child_cnt > 0) {
	    point_t vcenter = {0, 0, 0};
	    point_t vpoint;

	    for (size_t _k = 0; _k < _child_cnt; _k++) {
		vect_t diff;
		point_t label_pt;

		if (!_tclcad_draw_view_data_label_copy(gdvp, label_name, _k, NULL, label_pt, NULL))
		    continue;
		MAT4X3PNT(vpoint, model2view, label_pt);
		vcenter[Z] = vpoint[Z];
		VSUB2(diff, vpoint, vcenter);
		VSCALE(diff, diff, sf);
		VADD2(vpoint, vcenter, diff);
		MAT4X3PNT(label_pt, view2model, vpoint);
		(void)_tclcad_draw_view_data_label_point_set(gdvp, label_name, _k, label_pt);
	    }
	}
    }


    to_refresh_view(gdvp);
    return BRLCAD_OK;
}

int
to_mouse_data_scale(struct ged *gedp,
		    int argc,
		    const char *argv[],
		    ged_func_ptr UNUSED(func),
		    const char *usage,
		    int UNUSED(maxargs))
{
    int ret;
    const char *av[4];
    struct bu_vls scale_vls = BU_VLS_INIT_ZERO;
    void *gdvp;

    TO_COMMON_MOUSE_SCALE(gdvp, scale_vls, argc, argv, usage);
    ged_view_active_ctx_set(gedp, gdvp);

    av[0] = "to_data_scale";
    av[1] = (char *)argv[1];
    av[2] = bu_vls_addr(&scale_vls);
    av[3] = (char *)0;

    ret = to_data_scale(gedp, 3, (const char **)av, (ged_func_ptr)NULL, NULL, 4);

    bu_vls_free(&scale_vls);

    return ret;
}


int
to_mouse_scale(struct ged *gedp,
	       int argc,
	       const char *argv[],
	       ged_func_ptr UNUSED(func),
	       const char *usage,
	       int UNUSED(maxargs))
{
    int ret;
    const char *av[3];
    struct bu_vls zoom_vls = BU_VLS_INIT_ZERO;
    void *gdvp;

    TO_COMMON_MOUSE_SCALE(gdvp, zoom_vls, argc, argv, usage);
    ged_view_active_ctx_set(gedp, gdvp);

    av[0] = "zoom";
    av[1] = bu_vls_addr(&zoom_vls);
    av[2] = (char *)0;
    ret = ged_exec_zoom(gedp, 2, (const char **)av);
    bu_vls_free(&zoom_vls);

    if (ret == BRLCAD_OK) {
	struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(gdvp);
	if (0 < bu_vls_strlen(&tvd->gdv_callback)) {
	    Tcl_Eval(current_top->to_interp, bu_vls_addr(&tvd->gdv_callback));
	}

	to_refresh_view(gdvp);
    }

    return BRLCAD_OK;
}


int
to_mouse_protate(struct ged *gedp,
		 int argc,
		 const char *argv[],
		 ged_func_ptr UNUSED(func),
		 const char *usage,
		 int UNUSED(maxargs))
{
    int ret;
    const char *av[6];
    fastf_t dx, dy;
    point_t model;
    point_t view;
    mat_t inv_rot;
    struct bu_vls mrot_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 6) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[4], "%lf", &x) != 1 ||
	bu_sscanf(argv[5], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    tclcad_mouse_previous_get_set(&prev_x, &prev_y, gdvp, x, y);
    dx = y - prev_y;
    dy = x - prev_x;

    tclcad_mouse_clamp_delta(&dx, &dy, gdvp);

    dx *= tclcad_mouse_rotate_scale(gdvp);
    dy *= tclcad_mouse_rotate_scale(gdvp);

    VSET(view, dx, dy, 0.0);
    tclcad_mouse_view_inv_rotation(inv_rot, gdvp);
    MAT4X3PNT(model, inv_rot, view);

    bu_vls_printf(&mrot_vls, "%lf %lf %lf", V3ARGS(model));

    ged_view_active_ctx_set(gedp, gdvp);
    av[0] = "protate";
    av[1] = (char *)argv[2];
    av[2] = (char *)argv[3];
    av[3] = bu_vls_addr(&mrot_vls);
    av[4] = (char *)0;

    ret = ged_exec_protate(gedp, 4, (const char **)av);
    bu_vls_free(&mrot_vls);

    if (ret == BRLCAD_OK) {
	av[0] = "draw";
	av[1] = (char *)argv[2];
	av[2] = (char *)0;
	to_edit_redraw(gedp, 2, (const char **)av);
    }

    return BRLCAD_OK;
}


int
to_mouse_pscale(struct ged *gedp,
		int argc,
		const char *argv[],
		ged_func_ptr UNUSED(func),
		const char *usage,
		int UNUSED(maxargs))
{
    int ret, width;
    const char *av[6];
    fastf_t dx, dy;
    fastf_t sf;
    fastf_t inv_width;
    struct bu_vls sf_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 6) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[4], "%lf", &x) != 1 ||
	bu_sscanf(argv[5], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    tclcad_mouse_previous_get_set(&prev_x, &prev_y, gdvp, x, y);
    dx = x - prev_x;
    dy = prev_y - y;

    tclcad_mouse_clamp_delta(&dx, &dy, gdvp);

    width = tclcad_mouse_display_width(gdvp);
    inv_width = 1.0 / (fastf_t)width;
    dx *= inv_width * tclcad_mouse_scale_scale(gdvp);
    dy *= inv_width * tclcad_mouse_scale_scale(gdvp);

    if (fabs(dx) < fabs(dy))
	sf = 1.0 + dy;
    else
	sf = 1.0 + dx;

    bu_vls_printf(&sf_vls, "%lf", sf);

    ged_view_active_ctx_set(gedp, gdvp);
    av[0] = "pscale";
    av[1] = "-r";
    av[2] = (char *)argv[2];
    av[3] = (char *)argv[3];
    av[4] = bu_vls_addr(&sf_vls);
    av[5] = (char *)0;

    ret = ged_exec_pscale(gedp, 5, (const char **)av);
    bu_vls_free(&sf_vls);

    if (ret == BRLCAD_OK) {
	av[0] = "draw";
	av[1] = (char *)argv[2];
	av[2] = (char *)0;
	to_edit_redraw(gedp, 2, (const char **)av);
    }

    return BRLCAD_OK;
}


int
to_mouse_ptranslate(struct ged *gedp,
		    int argc,
		    const char *argv[],
		    ged_func_ptr UNUSED(func),
		    const char *usage,
		    int UNUSED(maxargs))
{
    int ret, width;
    const char *av[6];
    fastf_t dx, dy;
    point_t model;
    point_t view;
    fastf_t inv_width;
    mat_t inv_rot;
    struct bu_vls tvec_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 6) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[4], "%lf", &x) != 1 ||
	bu_sscanf(argv[5], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    tclcad_mouse_previous_get_set(&prev_x, &prev_y, gdvp, x, y);
    dx = x - prev_x;
    dy = prev_y - y;

    tclcad_mouse_clamp_delta(&dx, &dy, gdvp);

    width = tclcad_mouse_display_width(gdvp);
    inv_width = 1.0 / (fastf_t)width;
    /* ged_ptranslate expects things to be in local units */
    dx *= inv_width * bv_size_get(tclcad_mouse_bv_const(gdvp)) *
	gedp->dbip->dbi_base2local;
    dy *= inv_width * bv_size_get(tclcad_mouse_bv_const(gdvp)) *
	gedp->dbip->dbi_base2local;
    VSET(view, dx, dy, 0.0);
    tclcad_mouse_view_inv_rotation(inv_rot, gdvp);
    MAT4X3PNT(model, inv_rot, view);

    bu_vls_printf(&tvec_vls, "%lf %lf %lf", V3ARGS(model));

    ged_view_active_ctx_set(gedp, gdvp);
    av[0] = "ptranslate";
    av[1] = "-r";
    av[2] = (char *)argv[2];
    av[3] = (char *)argv[3];
    av[4] = bu_vls_addr(&tvec_vls);
    av[5] = (char *)0;

    ret = ged_exec_ptranslate(gedp, 5, (const char **)av);
    bu_vls_free(&tvec_vls);

    if (ret == BRLCAD_OK) {
	av[0] = "draw";
	av[1] = (char *)argv[2];
	av[2] = (char *)0;
	to_edit_redraw(gedp, 2, (const char **)av);
    }

    return BRLCAD_OK;
}


int
to_mouse_trans(struct ged *gedp,
	       int argc,
	       const char *argv[],
	       ged_func_ptr UNUSED(func),
	       const char *usage,
	       int UNUSED(maxargs))
{
    int ret, width;
    int ac;
    const char *av[4];
    fastf_t dx, dy;
    fastf_t inv_width;
    struct bu_vls trans_vls = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[2], "%lf", &x) != 1 ||
	bu_sscanf(argv[3], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    tclcad_mouse_previous_get_set(&prev_x, &prev_y, gdvp, x, y);
    dx = prev_x - x;
    dy = y - prev_y;

    tclcad_mouse_clamp_delta(&dx, &dy, gdvp);

    width = tclcad_mouse_display_width(gdvp);
    inv_width = 1.0 / (fastf_t)width;
    dx *= inv_width * bv_size_get(tclcad_mouse_bv_const(gdvp)) *
	gedp->dbip->dbi_local2base;
    dy *= inv_width * bv_size_get(tclcad_mouse_bv_const(gdvp)) *
	gedp->dbip->dbi_local2base;

    bu_vls_printf(&trans_vls, "%lf %lf 0", dx, dy);

    ged_view_active_ctx_set(gedp, gdvp);
    ac = 3;
    av[0] = "tra";
    av[1] = "-v";
    av[2] = bu_vls_addr(&trans_vls);
    av[3] = (char *)0;

    ret = ged_exec_tra(gedp, ac, (const char **)av);
    bu_vls_free(&trans_vls);

    if (ret == BRLCAD_OK) {
	struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(gdvp);
	if (0 < bu_vls_strlen(&tvd->gdv_callback)) {
	    Tcl_Eval(current_top->to_interp, bu_vls_addr(&tvd->gdv_callback));
	}

	to_refresh_view(gdvp);
    }

    return BRLCAD_OK;
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
