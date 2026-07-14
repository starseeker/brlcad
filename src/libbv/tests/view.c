/*                   L I B B V / V I E W . C
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

#include "common.h"

#include "bu/app.h"

#include <stdio.h>

#include "bv.h"
#include "bu/str.h"

static int
near_fastf(fastf_t a, fastf_t b)
{
    return NEAR_EQUAL(a, b, 1.0e-6);
}

static int
near_point(const point_t a, fastf_t x, fastf_t y, fastf_t z)
{
    return near_fastf(a[X], x) && near_fastf(a[Y], y) && near_fastf(a[Z], z);
}

static int
fail(const char *msg)
{
    fprintf(stderr, "%s\n", msg);
    return 1;
}

struct callback_state {
    int count;
    uint64_t flags;
    struct bv_context *ctx;
};

static void
context_callback(struct bv_context *ctx, uint64_t changed_flags,
	void *client_data)
{
    struct callback_state *state = (struct callback_state *)client_data;

    if (!state)
	return;

    state->count++;
    state->flags = changed_flags;
    state->ctx = ctx;
}

int
main(int UNUSED(argc), const char **UNUSED(argv))
{
    bu_setprogname("bv_test");
    struct bv *v = bv_create();
    struct bv *v2 = bv_create();
    struct bv_set *set = bv_set_create();
    struct bv_context_set *ctx_set = bv_context_set_create();
    struct bv_context ctx;
    struct bv_context *owned_ctx = NULL;
    struct callback_state cb_state = {0, 0, NULL};
    point_t min;
    point_t max;
    point_t center;
    point_t model;
    point_t keypoint = VINIT_ZERO;
    struct bv_grid_state grid;
    struct bv_axes_state axes;
    struct bv_other_state overlay = BV_OTHER_STATE_INIT;
    struct bv_params_state params;
    struct bv_snap_state snap;
    struct bv_mouse_delta_settings mouse_settings = BV_MOUSE_DELTA_SETTINGS_INIT;
    struct bv_view_info view_info = BV_VIEW_INFO_INIT;
    struct bv_lod_policy lod_policy;
    vect_t rvec = VINIT_ZERO;
    vect_t tvec = VINIT_ZERO;
    fastf_t sx = 1.0;
    fastf_t sy = 1.0;
    void *sentinel = (void *)v;
    int do_rot = 0;
    int do_tran = 0;

    if (!bv_is_valid(v))
	return fail("new view is not valid");
    if (!near_fastf(v->scale, BV_DEFAULT_SCALE) || !near_fastf(v->size, 2.0 * BV_DEFAULT_SCALE))
	return fail("new view scale/size defaults are incorrect");
    if (!bv_refresh_enabled_get(v) || bv_refresh_dirty_get(v))
	return fail("new view refresh defaults are incorrect");
    if (!bv_grid_state_get(&grid, v) ||
	    !near_fastf(grid.res_h, 1.0) || !near_fastf(grid.res_v, 1.0) ||
	    grid.res_major_h != 5 || grid.res_major_v != 5 ||
	    grid.color[0] != 255 || grid.color[1] != 255 || grid.color[2] != 255)
	return fail("new view grid defaults are incorrect");
    if (!bv_view_axes_state_get(&axes, v) ||
	    !near_point(axes.axes_pos, 0.80, -0.80, 0.0) ||
	    !near_fastf(axes.axes_size, 0.2) ||
	    !axes.pos_only || !axes.label_flag || !axes.triple_color)
	return fail("new view axes defaults are incorrect");
    if (!bv_params_state_get(&params, v) ||
	    !params.draw_size || !params.draw_center || !params.draw_az ||
	    !params.draw_el || !params.draw_tw || params.font_size != 20)
	return fail("new view params defaults are incorrect");
    if (!bv_snap_state_get(&snap, v) || snap.lines ||
	    snap.source_flags || !near_fastf(snap.tolerance_factor, 10.0) ||
	    snap.kind_mask != BV_SNAP_KIND_DEFAULT_MASK)
	return fail("new view snap defaults are incorrect");
    view_info.width = 200;
    view_info.height = 100;
    view_info.size = 20.0;
    view_info.lod.curve_scale = 3.0;
    view_info.lod.point_scale = 2.0;
    view_info.lod.bot_threshold = 42;
    if (!near_fastf(bv_view_lod_curve_scale(&view_info), 3.0) ||
	    bv_view_lod_bot_threshold(&view_info) != 42 ||
	    !near_fastf(bv_view_avg_sample_spacing(&view_info), 0.1))
	return fail("view info helpers failed");
    bv_view_info_sanitize(NULL);
    view_info.width = 0;
    view_info.height = -1;
    view_info.size = 0.0;
    view_info.lod.scale = -1.0;
    view_info.lod.curve_scale = 0.0;
    view_info.lod.point_scale = 0.0;
    bv_view_info_sanitize(&view_info);
    if (view_info.width != 1 || view_info.height != 1 ||
	    !near_fastf(view_info.size, 1.0) ||
	    !near_fastf(view_info.lod.scale, 1.0) ||
	    !near_fastf(view_info.lod.curve_scale, 1.0) ||
	    !near_fastf(view_info.lod.point_scale, 1.0))
	return fail("view info sanitize failed");
    bv_lod_policy_init(NULL);
    bv_lod_policy_init(&lod_policy);
    if (lod_policy.policy != BV_LOD_AUTO ||
	    !near_fastf(lod_policy.scale, 1.0) ||
	    !near_fastf(lod_policy.curve_scale, 1.0) ||
	    !near_fastf(lod_policy.point_scale, 1.0))
	return fail("lod policy defaults failed");

    if (!bv_dimensions_set(v, 800, 400))
	return fail("failed to set view dimensions");
    if (!bv_unit_conversion_set(v, 2.0, 0.5) ||
	    !near_fastf(bv_local2base_get(v), 2.0) ||
	    !near_fastf(bv_base2local_get(v), 0.5))
	return fail("view unit conversion accessors failed");
    if (!bv_perspective_set(v, 35.0) ||
	    !near_fastf(bv_perspective_get(v), 35.0) ||
	    !bv_coord_set(v, 'm') || bv_coord_get(v) != 'm' ||
	    !bv_rotate_about_set(v, 'e') || bv_rotate_about_get(v) != 'e')
	return fail("view scalar setters failed");
    if (bv_width_get(v) != 800 || bv_height_get(v) != 400 ||
	    !near_fastf(bv_scale_get(v), BV_DEFAULT_SCALE) ||
	    !near_fastf(bv_size_get(v), BV_DEFAULT_SCALE * 2.0) ||
	    !near_fastf(bv_perspective_get(v), 35.0) ||
	    !near_fastf(bv_local2base_get(v), 2.0) ||
	    !near_fastf(bv_base2local_get(v), 0.5))
	return fail("view scalar accessors failed");
    if (!bv_mouse_delta_settings_get(&mouse_settings, v) ||
	!near_fastf(mouse_settings.min_delta, -20.0) ||
	!near_fastf(mouse_settings.max_delta, 20.0) ||
	!near_fastf(mouse_settings.rotate_scale, 0.4) ||
	!near_fastf(mouse_settings.scale_scale, 2.0))
	return fail("mouse delta defaults are incorrect");
    mouse_settings.min_delta = -5.0;
    mouse_settings.max_delta = 5.0;
    mouse_settings.rotate_scale = 0.5;
    mouse_settings.scale_scale = 4.0;
    if (!bv_mouse_delta_settings_set(v, &mouse_settings))
	return fail("mouse delta settings set failed");
    fastf_t mouse_dx = 99.0;
    fastf_t mouse_dy = -99.0;
    if (!bv_mouse_delta_clamp(&mouse_dx, &mouse_dy, v) ||
	!near_fastf(mouse_dx, 5.0) || !near_fastf(mouse_dy, -5.0))
	return fail("mouse delta clamp did not use view settings");
    const fastf_t scale_before_mouse_adjust = bv_scale_get(v);
    if (!bv_mouse_delta_adjust(v, 99, 0, keypoint, BV_ADJUST_SCALE) ||
	!near_fastf(bv_scale_get(v), scale_before_mouse_adjust / 1.05))
	return fail("mouse delta scale adjustment did not use view settings");
    mouse_settings.min_delta = 1.0;
    if (bv_mouse_delta_settings_set(v, &mouse_settings))
	return fail("invalid mouse delta settings were accepted");
    if (!bv_user_data_set(v, sentinel) || bv_user_data_get(v) != sentinel)
	return fail("failed to set view user data");
    bv_context_init(&ctx, v);
    if (!bv_context_is_valid(&ctx) || bv_context_view(&ctx) != v ||
	    bv_context_view_const(&ctx) != v)
	return fail("view context did not wrap the view");
    owned_ctx = bv_context_create();
    if (!bv_context_is_valid(owned_ctx) || !bv_context_view(owned_ctx) ||
	    bv_context_view(owned_ctx) == v)
	return fail("owned view context create failed");
    if (!bv_context_name_set(owned_ctx, "owned") ||
	    !BU_STR_EQUAL(bv_context_name_get(owned_ctx), "owned"))
	return fail("view context name helpers failed");
    if (!bv_context_user_data_set(owned_ctx, sentinel) ||
	    bv_context_user_data_get(owned_ctx) != sentinel)
	return fail("view context user-data helpers failed");
    if (!bv_context_dimensions_set(owned_ctx, 320, 240) ||
	    bv_context_width_get(owned_ctx) != 320 ||
	    bv_context_height_get(owned_ctx) != 240)
	return fail("view context dimension helpers failed");
    if (!bv_context_callback_add(&ctx, context_callback, &cb_state) ||
	    !bv_context_notify(&ctx, BV_CONTEXT_CHANGED_REFRESH) ||
	    cb_state.count != 1 || cb_state.flags != BV_CONTEXT_CHANGED_REFRESH ||
	    cb_state.ctx != &ctx)
	return fail("view context notify failed");
    if (!bv_context_update(&ctx, 0) || cb_state.count != 2 ||
	    cb_state.flags != BV_CONTEXT_CHANGED_VIEW)
	return fail("view context update callback failed");
    if (!bv_context_callback_remove(&ctx, context_callback, &cb_state) ||
	    !bv_context_notify(&ctx, BV_CONTEXT_CHANGED_ALL) ||
	    cb_state.count != 2)
	return fail("view context callback remove failed");
    if (!bv_context_refresh_request(owned_ctx, 0x5) ||
	    !bv_refresh_dirty_get(bv_context_view(owned_ctx)) ||
	    !bv_context_refresh_complete(owned_ctx))
	return fail("view context refresh helpers failed");
    if (!bv_screen_to_view(&sx, &sy, v, 400.0, 200.0))
	return fail("failed to convert screen center to view");
    if (!near_fastf(sx, 0.0) || !near_fastf(sy, 0.0))
	return fail("screen center did not map to view origin");

    VSET(min, -1.0, -2.0, -3.0);
    VSET(max, 3.0, 2.0, 1.0);
    if (!bv_autoview_bounds(v, BV_AUTOVIEW_SCALE_DEFAULT, min, max))
	return fail("autoview bounds failed");
    if (!bv_center_get(center, v) || !near_point(center, 1.0, 0.0, -1.0))
	return fail("autoview center is incorrect");
    if (!near_fastf(v->scale, 2.0) || !near_fastf(v->size, 4.0))
	return fail("autoview scale/size are incorrect");
    if (!bv_screen_to_model(model, v, 400.0, 200.0) || !near_point(model, 1.0, 0.0, -1.0))
	return fail("screen center did not map to model center");

    if (bv_knobs_cmd_process(&rvec, &do_rot, &tvec, &do_tran, v, "aX", 10.0, 'v', 0, 0) != BRLCAD_OK)
	return fail("absolute X knob command failed");
    if (!do_tran || do_rot || !near_fastf(tvec[X], 10.0))
	return fail("absolute X knob command did not report the expected translation");

    if (!bv_refresh_request(v, 0x3) || !bv_refresh_dirty_get(v))
	return fail("refresh request did not mark view dirty");
    if (bv_refresh_consume(v) != 0x3 || bv_refresh_dirty_get(v))
	return fail("refresh consume did not return and clear dirty flags");
    if (!bv_refresh_suppress_begin(v) || !bv_refresh_request(v, 0x4) ||
	    bv_refresh_dirty_get(v) || !bv_refresh_suppressed_get(v))
	return fail("refresh suppression did not hold dirty flags");
    if (!bv_refresh_suppress_end(v) || bv_refresh_suppressed_get(v))
	return fail("refresh suppression did not clear");
    if (!bv_refresh_enabled_set(v, 0) || bv_refresh_enabled_get(v) ||
	    !bv_refresh_request(v, 0x8) || bv_refresh_dirty_get(v))
	return fail("disabled refresh did not hold dirty flags");
    if (!bv_refresh_enabled_set(v, 1) || !bv_refresh_drawn_count_set(v, 7) ||
	    bv_refresh_drawn_count_get(v) != 7)
	return fail("refresh drawn count state failed");
    if (!bv_frametime_set(v, 1234) || bv_frametime_get(v) != 1234 ||
	    !bv_zclip_set(v, 1) || !bv_zclip_get(v) ||
	    !bv_framebuffer_mode_set(v, 2) || bv_framebuffer_mode_get(v) != 2 ||
	    !bv_cleared_set(v, 1) || !bv_cleared_get(v))
	return fail("view bookkeeping state failed");
    grid.snap = 1;
    grid.res_h = 2.0;
    grid.res_major_h = 7;
    if (!bv_grid_state_set(v, &grid) || !bv_grid_state_get(&grid, v) ||
	    !grid.snap || !near_fastf(grid.res_h, 2.0) ||
	    grid.res_major_h != 7)
	return fail("grid state set/get failed");
    overlay.gos_draw = 1;
    overlay.gos_line_color[0] = 12;
    overlay.gos_text_color[1] = 34;
    overlay.gos_font_size = 18;
    if (!bv_scale_overlay_state_set(v, &overlay) ||
	    !bv_scale_overlay_state_get(&overlay, v) ||
	    !overlay.gos_draw || overlay.gos_line_color[0] != 12 ||
	    overlay.gos_text_color[1] != 34 || overlay.gos_font_size != 18)
	return fail("overlay state set/get failed");
    snap.lines = 4;
    snap.source_flags = BV_SNAP_DB | BV_SNAP_TCL;
    snap.kind_mask = BV_SNAP_KIND_ENDPOINT;
    snap.tolerance_factor = -5.0;
    if (!bv_snap_state_set(v, &snap) || !bv_snap_state_get(&snap, v) ||
	    snap.lines != 1 || snap.source_flags != (BV_SNAP_DB | BV_SNAP_TCL) ||
	    snap.kind_mask != BV_SNAP_KIND_ENDPOINT ||
	    !near_fastf(snap.tolerance_factor, 1.0))
	return fail("snap state set/get failed");
    if (!bv_snap_tolerance_factor_set(v, 3.5) ||
	    !near_fastf(bv_snap_tolerance_factor_get(v), 3.5) ||
	    !bv_snap_lines_set(v, 0) || bv_snap_lines_get(v))
	return fail("snap scalar accessors failed");

    if (!bv_name_set(v, "primary") || !bv_set_add(set, v))
	return fail("failed to add named view to set");
    if (bv_set_find(set, "primary") != v)
	return fail("failed to find named view in set");
    if (!bv_set_remove(set, v) || bv_set_find(set, "primary") != NULL)
	return fail("failed to remove view from set");

    struct bu_ptbl *views = NULL;

    if (!bv_context_set_add(ctx_set, &ctx) ||
	    bv_context_set_find(ctx_set, "primary") != &ctx ||
	    !(views = bv_context_set_views(ctx_set)) ||
	    BU_PTBL_LEN(views) != 1)
	return fail("failed to add named context to context set");
    if (!bv_name_set(bv_context_view(owned_ctx), "owned") ||
	    !bv_context_set_add(ctx_set, owned_ctx) ||
	    bv_context_set_find(ctx_set, "owned") != owned_ctx ||
	    !(views = bv_context_set_views(ctx_set)) ||
	    BU_PTBL_LEN(views) != 2)
	return fail("failed to add owned context to context set");
    if (!bv_context_set_remove(ctx_set, &ctx) ||
	    bv_context_set_find(ctx_set, "primary") != NULL ||
	    !(views = bv_context_set_views(ctx_set)) ||
	    BU_PTBL_LEN(views) != 1)
	return fail("failed to remove context from context set");
    if (!bv_context_set_attach(ctx_set, &ctx) ||
	    ctx.view_set != ctx_set ||
	    !(views = bv_context_set_views(ctx_set)) ||
	    BU_PTBL_LEN(views) != 1 ||
	    !bv_context_set_attach(NULL, &ctx) ||
	    ctx.view_set != NULL)
	return fail("failed to attach context to context set without insertion");

    if (!bv_copy(v2, v) || !bv_center_get(center, v2) || !near_point(center, 1.0, 0.0, -1.0))
	return fail("view copy did not preserve model center");
    mat_t center_mat;
    if (!bv_center_mat_get(center_mat, v2) ||
	    !near_fastf(center_mat[MDX], -1.0) ||
	    !near_fastf(center_mat[MDY], 0.0) ||
	    !near_fastf(center_mat[MDZ], 1.0) ||
	    !bv_center_mat_set(v2, center_mat))
	return fail("view center matrix accessors failed");
    if (bv_user_data_get(v2) != sentinel ||
	    bv_refresh_drawn_count_get(v2) != 7 ||
	    bv_frametime_get(v2) != 1234 || !bv_zclip_get(v2) ||
	    bv_framebuffer_mode_get(v2) != 2 ||
	    !bv_cleared_get(v2))
	return fail("view copy did not preserve bookkeeping state");
    if (!bv_grid_state_get(&grid, v2) || !near_fastf(grid.res_h, 2.0) ||
	    grid.res_major_h != 7 ||
	    !bv_snap_state_get(&snap, v2) ||
	    snap.kind_mask != BV_SNAP_KIND_ENDPOINT ||
	    !near_fastf(snap.tolerance_factor, 3.5) ||
	    !bv_scale_overlay_state_get(&overlay, v2) ||
	    overlay.gos_font_size != 18)
	return fail("view copy did not preserve passive view state");

    bv_set_destroy(set);
    bv_context_set_destroy(ctx_set);
    bv_context_destroy(owned_ctx);
    bv_context_free(&ctx);
    bv_destroy(v2);
    bv_destroy(v);

    return 0;
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
