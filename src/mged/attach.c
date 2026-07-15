/*                        A T T A C H . C
 * BRL-CAD
 *
 * Copyright (c) 1985-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file mged/attach.c
 *
 */

#include "common.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef HAVE_SYS_TIME_H
#  include <sys/time.h>		/* for struct timeval */
#endif

#include "tcl.h"
#ifdef HAVE_TK
#  include "tk.h"
#endif

#include "bnetwork.h"

/* Make sure this comes after bio.h (for Windows) */
#ifdef HAVE_GL_GL_H
#  include <GL/gl.h>
#endif

#include "vmath.h"
#include "bu/env.h"
#include "bu/ptbl.h"
#include "ged.h"
#include "ged/draw_obol.h"
#include "ged/view.h"
#include "brlobol/display_endpoint.h"
#include "tclcad.h"

#include "./mged.h"
#include "./sedit.h"
#include "./mged_display.h"

// FIXME: Globals
/* Geometry display instances used by MGED */
struct bu_ptbl active_display_set = BU_PTBL_INIT_ZERO;  /* set of active graphics views */
struct mged_display *mged_initial_display = NULL;


extern struct _color_scheme default_color_scheme;
static unsigned long mged_tkobol_count = 0;

/* attach accepts only the built-in Tk Obol host.  Its only pre-Tk option is
 * the target X display, so do not instantiate a generic DM just to parse it. */
static const char *
mged_tkobol_display_arg(int argc, const char *argv[])
{
    const char *display = NULL;

    for (int i = 1; i < argc - 1; i++) {
	if (!BU_STR_EQUAL(argv[i], "-d"))
	    continue;
	if (++i >= argc - 1)
	    return NULL;
	display = argv[i];
    }

    return display;
}

static void
mged_obol_input_refresh(struct mged_state *s)
{
    if (!s || s->mged_curr_display == MGED_DISPLAY_NULL)
	return;
    mged_refresh_request_view(s, view_state, GED_VIEW_REFRESH_VIEW);
    mged_display_repaint_request(s->mged_curr_display, MGED_REPAINT_INTERACTION);
}

static void
mged_obol_input_update_gui(struct mged_state *s, const char *name, int value)
{
    const char *pathname = s ? mged_display_pathname(s->mged_curr_display) : NULL;
    if (!s || !s->interp || !pathname || !name)
	return;

    Tcl_DString command;
    Tcl_DStringInit(&command);
    Tcl_DStringAppendElement(&command, "update_gui");
    Tcl_DStringAppendElement(&command, pathname);
    Tcl_DStringAppendElement(&command, name);
    char value_string[32] = {0};
    snprintf(value_string, sizeof(value_string), "%d", value);
    Tcl_DStringAppendElement(&command, value_string);
    Tcl_Obj *saved_result = Tcl_GetObjResult(s->interp);
    Tcl_IncrRefCount(saved_result);
    (void)Tcl_EvalEx(s->interp, Tcl_DStringValue(&command), -1,
	TCL_EVAL_GLOBAL);
    Tcl_SetObjResult(s->interp, saved_result);
    Tcl_DecrRefCount(saved_result);
    Tcl_DStringFree(&command);
}

static int
mged_obol_input_forwarding(const struct mged_state *s)
{
    const char *pathname = s ? mged_display_pathname(s->mged_curr_display) : NULL;
    if (!s || !s->interp || !pathname)
	return 0;
    const char *forwarding = Tcl_GetVar2(s->interp, "forwarding_key",
	pathname, TCL_GLOBAL_ONLY);
    return forwarding && atoi(forwarding) != 0;
}

static int
mged_obol_input_pan_begin(struct mged_state *s,
	const BRLObolInputEvent *event)
{
    if (!s || !event || event->type != BRLOBOL_INPUT_POINTER_PRESS ||
	event->button != 0 || s->mged_curr_display == MGED_DISPLAY_NULL ||
	s->global_editing_state != ST_VIEW ||
	mged_variables->mv_transform != 'v' || mged_variables->mv_rateknobs ||
	!s->dbip)
	return 0;

    struct bv_grid_state grid = BV_GRID_STATE_INIT;
    if (!mged_display_grid_state_get(s->mged_curr_display, &grid) || !grid.snap)
	return 0;

    /* This is the same baseline the historical "dm am t" path establishes
     * before applying snapped view translations. */
    snap_view_center_to_grid(s);
    mouse_dx = 0;
    mouse_dy = 0;
    s->mged_curr_display->display_input_snap_pan_active = 1;
    (void)bv_context_update(view_state->vs_gvp, BV_CONTEXT_CHANGED_VIEW);
    mged_obol_input_refresh(s);
    return 1;
}

static int
mged_obol_input_camera_adjust(struct mged_state *s,
	BRLObolInputAction action, const BRLObolInputEvent *event)
{
    const int wheel_zoom = event && event->type == BRLOBOL_INPUT_WHEEL &&
	action == BRLOBOL_ACTION_VIEW_ZOOM;
    if (!s || !event ||
	(!wheel_zoom && event->type != BRLOBOL_INPUT_POINTER_MOTION) ||
	s->mged_curr_display == MGED_DISPLAY_NULL ||
	s->global_editing_state != ST_VIEW ||
	mged_variables->mv_transform != 'v' ||
	mged_variables->mv_rateknobs ||
	(!wheel_zoom && !(event->buttons & (1u << 0))))
	return 0;

    struct bv_grid_state grid = BV_GRID_STATE_INIT;
	if (!mged_display_grid_state_get(s->mged_curr_display, &grid))
	return 0;

    if (grid.snap) {
	if (action != BRLOBOL_ACTION_VIEW_PAN ||
	    !s->mged_curr_display->display_input_snap_pan_active)
	    return 0;
	const int width = bv_width_get(mged_view_state_view(view_state));
	const int height = bv_height_get(mged_view_state_view(view_state));
	if (width <= 0 || height <= 0)
	    return 0;
	s->mged_curr_display->display_mouse_dx += event->dx;
	s->mged_curr_display->display_mouse_dy += event->dy;
	const fastf_t aspect = (fastf_t)width / (fastf_t)height;
	snap_view_to_grid(s,
	    s->mged_curr_display->display_mouse_dx / (fastf_t)width * 2.0,
	    -s->mged_curr_display->display_mouse_dy / (fastf_t)height / aspect * 2.0);
	(void)bv_context_update(view_state->vs_gvp, BV_CONTEXT_CHANGED_VIEW);
	mged_obol_input_refresh(s);
	return 1;
    }

    unsigned long long flags = BV_ADJUST_IDLE;
    if (action == BRLOBOL_ACTION_VIEW_ROTATE)
	flags = BV_ADJUST_ROT;
    else if (action == BRLOBOL_ACTION_VIEW_PAN)
	flags = BV_ADJUST_TRANS;
    else if (action == BRLOBOL_ACTION_VIEW_ZOOM)
	flags = BV_ADJUST_SCALE;
    else
	return 0;

    struct bv *view = mged_view_state_view(view_state);
    if (!view)
	return 0;

    if (wheel_zoom) {
	point_t origin = VINIT_ZERO;
	/* Keep the toolkit-neutral wheel convention identical to QgCanvasInput. */
	if (bv_adjust(view, 100 + event->wheelDelta, 100,
		origin, 0, BV_ADJUST_SCALE)) {
	    (void)bv_context_update(view_state->vs_gvp, BV_CONTEXT_CHANGED_VIEW);
	    mged_obol_input_refresh(s);
	}
	return 1;
    }

    point_t center;
    mat_t center_mat;
    bv_center_mat_get(center_mat, view);
    MAT_DELTAS_GET_NEG(center, center_mat);

    if (bv_mouse_delta_adjust(view, event->dx, event->dy, center, flags)) {
	(void)bv_context_update(view_state->vs_gvp, BV_CONTEXT_CHANGED_VIEW);
	mged_obol_input_refresh(s);
    }

    /* A zero first motion still belongs to this semantic drag, not the
     * legacy X-motion path that would otherwise apply a second adjustment. */
    return 1;
}

static int
mged_obol_input_action_apply(struct mged_state *s,
	BRLObolInputAction action, const BRLObolInputEvent *event)
{
    if (!s || !event || s->mged_curr_display == MGED_DISPLAY_NULL)
	return 0;
    if (mged_obol_input_forwarding(s))
	return 0;

    if (action == BRLOBOL_ACTION_VIEW_ROTATE ||
	action == BRLOBOL_ACTION_VIEW_PAN ||
	action == BRLOBOL_ACTION_VIEW_ZOOM)
	return mged_obol_input_camera_adjust(s, action, event);

    switch (action) {
	case BRLOBOL_ACTION_TOGGLE_ADC: {
	    if (!brlobol_display_endpoint_input_faceplate_toggle_apply(
		ged_view_context_display_endpoint_get(view_state->vs_gvp),
		view_state->vs_gvp, action, NULL))
		return 0;
	    mged_obol_input_refresh(s);
	    return 1;
	}
	case BRLOBOL_ACTION_TOGGLE_MODEL_AXES: {
	    int visible = 0;
	    if (!brlobol_display_endpoint_input_faceplate_toggle_apply(
		ged_view_context_display_endpoint_get(view_state->vs_gvp),
		view_state->vs_gvp, action, &visible))
		return 0;
	    axes_state->ax_model_draw = visible;
	    mged_obol_input_update_gui(s, "model_draw", visible);
	    mged_obol_input_refresh(s);
	    return 1;
	}
	case BRLOBOL_ACTION_TOGGLE_VIEW_AXES: {
	    int visible = 0;
	    if (!brlobol_display_endpoint_input_faceplate_toggle_apply(
		ged_view_context_display_endpoint_get(view_state->vs_gvp),
		view_state->vs_gvp, action, &visible))
		return 0;
	    axes_state->ax_view_draw = visible;
	    mged_obol_input_update_gui(s, "view_draw", visible);
	    mged_obol_input_refresh(s);
	    return 1;
	}
	case BRLOBOL_ACTION_VIEW_2:
	case BRLOBOL_ACTION_VIEW_3:
	case BRLOBOL_ACTION_VIEW_4:
	case BRLOBOL_ACTION_VIEW_5:
	case BRLOBOL_ACTION_VIEW_6:
	case BRLOBOL_ACTION_VIEW_7:
	case BRLOBOL_ACTION_VIEW_FRONT:
	case BRLOBOL_ACTION_VIEW_TOP:
	case BRLOBOL_ACTION_VIEW_BOTTOM:
	case BRLOBOL_ACTION_VIEW_LEFT:
	case BRLOBOL_ACTION_VIEW_REAR:
	case BRLOBOL_ACTION_VIEW_RIGHT:
	    if (!brlobol_input_view_orientation_apply(view_state->vs_gvp, action))
		return 0;
	    mged_obol_input_refresh(s);
	    return 1;
	case BRLOBOL_ACTION_VIEW_PAN_BEGIN:
	    return mged_obol_input_pan_begin(s, event);
	case BRLOBOL_ACTION_VIEW_PRIMARY_RELEASE:
	case BRLOBOL_ACTION_VIEW_PAN_END:
	    s->mged_curr_display->display_input_snap_pan_active = 0;
	    return 0;
	default:
	    return 0;
    }
}

static int
mged_obol_input_action(void *data, BRLObolInputAction action,
	const BRLObolInputEvent *event)
{
    struct mged_display *mdmp = (struct mged_display *)data;
    struct mged_state *s = mdmp ? mdmp->display_input_state : NULL;
    if (!s || !mdmp || !event)
	return 0;

    struct mged_display *saved = s->mged_curr_display;
    if (saved != mdmp)
	mged_current_display_set(s, mdmp);
    const int ret = mged_obol_input_action_apply(s, action, event);
    if (ret > 0 && event->type == BRLOBOL_INPUT_POINTER_MOTION) {
	mdmp->display_input_motion_pending = 1;
	mdmp->display_input_motion_timestamp = (unsigned long)event->timestamp;
	mdmp->display_input_motion_x = event->x;
	mdmp->display_input_motion_y = event->y;
    }
    if (saved != mdmp)
	mged_current_display_set(s, saved);
    return ret;
}

static const BRLObolInputProfile *
mged_obol_active_mode_profile(void)
{
    static const BRLObolInputBinding bindings[] = {
	{BRLOBOL_INPUT_POINTER_MOTION, BRLOBOL_INPUT_ANY, BRLOBOL_INPUT_ANY,
	 0, 0, 10, BRLOBOL_ACTION_APP_MGED_ACTIVE_MODE_MOTION}
    };
    static const BRLObolInputProfile profile = {
	"mged-active-mode", bindings, sizeof(bindings) / sizeof(bindings[0])
    };
    return &profile;
}

static int
mged_obol_input_semantic_action(void *data, BRLObolInputAction action,
	const BRLObolInputEvent *event)
{
    struct mged_display *mdmp = (struct mged_display *)data;
    struct mged_state *s = mdmp ? mdmp->display_input_state : NULL;
    if (!s || !mdmp || !event ||
	action != BRLOBOL_ACTION_APP_MGED_ACTIVE_MODE_MOTION ||
	event->type != BRLOBOL_INPUT_POINTER_MOTION)
	return BRLOBOL_INPUT_RESULT_UNHANDLED;

    struct mged_display *saved = s->mged_curr_display;
    if (saved != mdmp)
	mged_current_display_set(s, mdmp);
    if (am_mode == AMM_IDLE) {
	if (saved != mdmp)
	    mged_current_display_set(s, saved);
	return BRLOBOL_INPUT_RESULT_UNHANDLED;
    }
    const int handled = mged_display_motion(s, event->x, event->y) > 0;
    if (handled) {
	mdmp->display_input_motion_pending = 1;
	mdmp->display_input_motion_timestamp = (unsigned long)event->timestamp;
	mdmp->display_input_motion_x = event->x;
	mdmp->display_input_motion_y = event->y;
    }
    if (saved != mdmp)
	mged_current_display_set(s, saved);
    return handled ? BRLOBOL_INPUT_RESULT_HANDLED :
	BRLOBOL_INPUT_RESULT_UNHANDLED;
}

void
mged_obol_input_semantic_mode_sync(struct mged_state *s)
{
    if (!s || s->mged_curr_display == MGED_DISPLAY_NULL ||
	!s->mged_curr_display->display_view_state)
	return;
    struct mged_display *mdmp = s->mged_curr_display;
    brlobol_display_endpoint_t *endpoint =
	ged_view_context_display_endpoint_get(mdmp->display_view_state->vs_gvp);
    if (!endpoint || !(brlobol_display_endpoint_host_capabilities(endpoint) &
	BRLOBOL_HOST_CAP_INPUT))
	return;

    if (am_mode == AMM_IDLE) {
	(void)brlobol_display_endpoint_input_semantic_action_handler_clear_if(
	    endpoint, mged_obol_input_semantic_action, mdmp);
	(void)brlobol_display_endpoint_input_semantic_profile_clear_if(endpoint,
	    mdmp);
	return;
    }

    if (!brlobol_display_endpoint_input_semantic_profile_set(endpoint,
	mged_obol_active_mode_profile(), mdmp))
	return;
    if (!brlobol_display_endpoint_input_semantic_action_handler_set(endpoint,
	mged_obol_input_semantic_action, mdmp)) {
	(void)brlobol_display_endpoint_input_semantic_profile_clear_if(endpoint,
	    mdmp);
    }
}

int
mged_obol_input_motion_consumed(struct mged_display *mdmp,
	unsigned long timestamp, int x, int y)
{
    if (!mdmp || !mdmp->display_input_motion_pending ||
	mdmp->display_input_motion_timestamp != timestamp ||
	mdmp->display_input_motion_x != x || mdmp->display_input_motion_y != y)
	return 0;
    mdmp->display_input_motion_pending = 0;
    return 1;
}

/* If we changed the active dm, need to update GEDP as well.. */
void mged_current_display_set(struct mged_state *s, struct mged_display *nc)
{
    // Normally we can assume mged_state is present, since it is allocated early
    // in the application startup, but mged_current_display_set is called from doEvent which
    // gets triggered even during shutdown.  So we need to sanity check in this
    // instance.
    if (!s)
	return;

    // Make sure the magic is non-zero.  We don't want to bomb if we're
    // shutting down and some pending function callback still has the
    // non-nulled MGED_STATE value, so we just do the check.
    if (UNLIKELY(( ((uintptr_t)(s) == 0) /* non-zero pointer */
		    || ((uintptr_t)(s) & (sizeof((uintptr_t)(s))-1)) /* aligned ptr */
		    || (*((const uint32_t *)(s)) != (uint32_t)(MGED_STATE_MAGIC)) /* matches value */
		 ))) {
	return;
    }

    s->mged_curr_display = nc;
    if (s->gedp) {
	if (nc != MGED_DISPLAY_NULL && nc->display_view_state) {
	    ged_view_active_ctx_set(s->gedp, nc->display_view_state->vs_gvp);
	} else {
	    ged_view_active_ctx_set(s->gedp, NULL);
	}
    }
}

void
mged_display_adc_state_set(struct mged_display *dm, const struct bv_adc_state *adc)
{
    struct bv_adc_state bv_adc;

    if (!dm || !dm->display_view_state || !dm->display_view_state->vs_gvp)
	return;
    memcpy(&bv_adc, adc, sizeof(bv_adc));
    bv_adc_state_set(mged_view_state_view(dm->display_view_state), &bv_adc);
}

int
mged_display_adc_visibility_set(struct mged_display *dm, int enabled)
{
    if (!dm || !dm->display_view_state || !dm->display_view_state->vs_gvp)
	return 0;

    void *view_ctx = dm->display_view_state->vs_gvp;
    if (ged_view_context_display_endpoint_get(view_ctx)) {
	struct brlobol_endpoint_property_value value =
	    BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	value.type = BRLOBOL_ENDPOINT_PROPERTY_BOOL;
	value.bool_value = enabled ? 1 : 0;
	return ged_view_context_display_property_set(view_ctx,
	    "view.faceplate.adc.visible", &value) ==
	    BRLOBOL_ENDPOINT_PROPERTY_OK;
    }

    struct bv_adc_state adc = {0};
    if (!mged_display_adc_state_get(dm, &adc))
	return 0;
    adc.draw = enabled ? 1 : 0;
    mged_display_adc_state_set(dm, &adc);
    return 1;
}

int
mged_display_init(
	struct mged_state *s,
	struct mged_display *o_dm,
	const char *host_type,
	int argc,
	const char *argv[])
{
    struct bu_vls vls = BU_VLS_INIT_ZERO;
    const char *failure = "unknown direct endpoint error";

    if (!BU_STR_EQUAL(host_type, "tkobol"))
	return TCL_ERROR;

    mged_display_var_init(s, o_dm);

    /* register application provided routines */
    cmd_hook = mged_display_command;

    void *ctx = view_state->vs_gvp;

    int width = bv_width_get(mged_view_state_view(view_state));
    int height = bv_height_get(mged_view_state_view(view_state));
    int toplevel = 1;
    int software = 0;
    struct bu_vls pathname = BU_VLS_INIT_ZERO;
    for (int i = 1; i < argc - 1; i++) {
	if (BU_STR_EQUAL(argv[i], "sw")) {
	    software = 1;
	    continue;
	}
	if (BU_STR_EQUAL(argv[i], "hw")) {
	    software = 0;
	    continue;
	}
	if (i + 1 >= argc - 1)
	    continue;
	if (BU_STR_EQUAL(argv[i], "-W"))
	    width = atoi(argv[++i]);
	else if (BU_STR_EQUAL(argv[i], "-N"))
	    height = atoi(argv[++i]);
	else if (BU_STR_EQUAL(argv[i], "-S") ||
		BU_STR_EQUAL(argv[i], "-s")) {
	    width = height = atoi(argv[++i]);
	} else if (BU_STR_EQUAL(argv[i], "-t"))
	    toplevel = atoi(argv[++i]) ? 1 : 0;
	else if (BU_STR_EQUAL(argv[i], "-n")) {
	    const char *name = argv[++i];
	    if (name[0] == '.')
		bu_vls_strcpy(&pathname, name);
	    else
		bu_vls_printf(&pathname, ".%s", name);
	} else if (BU_STR_EQUAL(argv[i], "-d") ||
		BU_STR_EQUAL(argv[i], "-i")) {
	    i++;
	}
    }
    if (width <= 0)
	width = 512;
    if (height <= 0)
	height = 512;
    if (!bu_vls_strlen(&pathname))
	bu_vls_printf(&pathname, ".tkobol%lu", ++mged_tkobol_count);
    bu_vls_strcpy(&s->mged_curr_display->display_pathname, bu_vls_cstr(&pathname));
    (void)bv_dimensions_set(mged_view_state_view(view_state), width, height);


    if (!tclcad_obol_host_factories_register()) {
	failure = "Tk Obol host factory registration failed";
	goto attach_fail;
    }
    brlobol_display_endpoint_t *endpoint =
	brlobol_display_endpoint_create(NULL, 0);

    if (!endpoint) {
	failure = "Obol display endpoint creation failed";
	goto attach_fail;
    }
    if (!brlobol_display_endpoint_render_engine_set(endpoint,
	    software ? BRLOBOL_RENDER_ENGINE_SW : BRLOBOL_RENDER_ENGINE_HW)) {
	failure = "Obol render-engine selection failed";
	brlobol_display_endpoint_destroy(endpoint);
	goto attach_fail;
    }
    if (!ged_view_context_display_endpoint_set(ctx, endpoint, 1)) {
	failure = "GED view endpoint attachment failed";
	brlobol_display_endpoint_destroy(endpoint);
	goto attach_fail;
    }

    struct brlobol_host_desc desc = {0};
    desc.struct_size = sizeof(desc);
    desc.mode = toplevel ? BRLOBOL_HOST_MODE_TOPLEVEL :
	BRLOBOL_HOST_MODE_EMBEDDED;
    desc.width = (unsigned int)width;
    desc.height = (unsigned int)height;
    desc.device_pixel_ratio = 1.0;
    desc.visible = 1;
    desc.required_capabilities = software ?
	BRLOBOL_HOST_CAP_PIXEL_PRESENT : BRLOBOL_HOST_CAP_SYSTEM_GL;
    desc.title = bu_vls_cstr(&pathname);
    desc.native_id_hint = bu_vls_cstr(&pathname);
    desc.application_context = s->interp;
    if (!brlobol_display_endpoint_host_open(endpoint,
	    software ? "tk-photo" : "tk-gl", &desc)) {
	failure = software ? "TkPhoto host open failed" :
	    "Tk OpenGL host open failed";
	(void)ged_view_context_display_endpoint_set(ctx, NULL, 0);
	goto attach_fail;
    }

#ifdef HAVE_TK
    struct bu_vls widget_path = BU_VLS_INIT_ZERO;
    if (toplevel)
	bu_vls_printf(&widget_path, "%s.__obol", bu_vls_cstr(&pathname));
    else
	bu_vls_strcpy(&widget_path, bu_vls_cstr(&pathname));
    Tk_Window host_window = Tk_NameToWindow(s->interp,
	bu_vls_cstr(&widget_path), Tk_MainWindow(s->interp));
    if (!host_window) {
	failure = "Tk host window lookup failed";
	bu_vls_free(&widget_path);
	(void)ged_view_context_display_endpoint_set(ctx, NULL, 0);
	goto attach_fail;
    }
    Tk_MakeWindowExist(host_window);
    s->mged_curr_display->display_native_id = Tk_WindowId(host_window);
    bu_vls_free(&widget_path);
#endif
    s->mged_curr_display->display_hosted = 1;
    s->mged_curr_display->display_input_state = s;
    s->mged_curr_display->display_input_motion_pending = 0;
    s->mged_curr_display->display_input_snap_pan_active = 0;
    (void)brlobol_display_endpoint_input_profile_set(endpoint,
	brlobol_input_default_view_profile());
    (void)brlobol_display_endpoint_input_action_handler_set(endpoint,
	mged_obol_input_action, s->mged_curr_display);
	/* A mode may have been selected before this view was attached. */
	mged_obol_input_semantic_mode_sync(s);

#ifdef HAVE_TK
    if (tkwin != NULL) {
	Tk_DeleteGenericHandler(doEvent, (ClientData)s);
	Tk_CreateGenericHandler(doEvent, (ClientData)s);
    }
#endif

    const char *logical_path = mged_display_pathname(s->mged_curr_display);
    if (logical_path) {
	(void)Tcl_SetVar2(s->interp, "mged_obol_semantic_input", logical_path,
	    "1", TCL_GLOBAL_ONLY);
	bu_vls_printf(&vls, "mged_bind_endpoint %s", logical_path);
	Tcl_Eval(s->interp, bu_vls_cstr(&vls));
    }
    bu_vls_free(&vls);
    bu_vls_free(&pathname);

    return TCL_OK;

attach_fail:
    Tcl_AppendResult(s->interp, "attach(tkobol): ", failure, "\n",
	    (char *)NULL);
    bu_vls_free(&pathname);
    return TCL_ERROR;
}

int
mged_obol_framebuffer_ensure(struct mged_state *s)
{
    if (!s || !s->gedp || !s->mged_curr_display || !view_state ||
	!view_state->vs_gvp)
	return 0;
    return ged_draw_obol_framebuffer_backend_ensure_for_view(s->gedp,
	view_state->vs_gvp) == BRLCAD_OK;
}


void
mged_slider_init_vls(struct mged_display *p)
{
    bu_vls_init(&p->display_fps_name);
    bu_vls_init(&p->display_aet_name);
    bu_vls_init(&p->display_ang_name);
    bu_vls_init(&p->display_center_name);
    bu_vls_init(&p->display_size_name);
    bu_vls_init(&p->display_adc_name);
}


void
mged_slider_free_vls(struct mged_display *p)
{
    if (BU_VLS_IS_INITIALIZED(&p->display_fps_name)) {
	bu_vls_free(&p->display_fps_name);
	bu_vls_free(&p->display_aet_name);
	bu_vls_free(&p->display_ang_name);
	bu_vls_free(&p->display_center_name);
	bu_vls_free(&p->display_size_name);
	bu_vls_free(&p->display_adc_name);
    }
    if (BU_VLS_IS_INITIALIZED(&p->display_pathname))
	bu_vls_free(&p->display_pathname);
}


void
mged_obol_display_detach(struct mged_state *s, struct mged_display *mdmp)
{
    if (!s || !s->gedp || !mdmp || !mdmp->display_view_state)
	return;

    const char *logical_path = mged_display_pathname(mdmp);
    brlobol_display_endpoint_t *endpoint =
	ged_view_context_display_endpoint_get(mdmp->display_view_state->vs_gvp);
    if (endpoint) {
	(void)brlobol_display_endpoint_input_action_handler_clear_if(endpoint,
	    mged_obol_input_action, mdmp);
	    (void)brlobol_display_endpoint_input_semantic_action_handler_clear_if(
		endpoint, mged_obol_input_semantic_action, mdmp);
	    (void)brlobol_display_endpoint_input_semantic_profile_clear_if(endpoint,
		mdmp);
	}
    if (s->interp && logical_path)
	(void)Tcl_UnsetVar2(s->interp, "mged_obol_semantic_input", logical_path,
	    TCL_GLOBAL_ONLY);
    mdmp->display_input_state = NULL;
    mdmp->display_input_motion_pending = 0;
    mdmp->display_input_snap_pan_active = 0;
    (void)ged_view_context_display_endpoint_set(mdmp->display_view_state->vs_gvp,
	    NULL, 0);
}

static int
release(struct mged_state *s, char *name, int need_close)
{
    struct mged_display *save_display = MGED_DISPLAY_NULL;

    if (name != NULL) {
	struct mged_display *p = MGED_DISPLAY_NULL;

	if (BU_STR_EQUAL("nu", name))
	    return TCL_OK;  /* Ignore */

	for (size_t i = 0; i < BU_PTBL_LEN(&active_display_set); i++) {
	    struct mged_display *m_dmp = (struct mged_display *)BU_PTBL_GET(&active_display_set, i);
	    if (!m_dmp)
		continue;

	    const char *display_path = mged_display_pathname(m_dmp);
	    if (!display_path || !BU_STR_EQUAL(name, display_path))
		continue;

	    /* found it */
	    p = m_dmp;
	    if (p != s->mged_curr_display) {
		save_display = s->mged_curr_display;
		mged_current_display_set(s, p);
	    }
	    break;
	}

	if (p == MGED_DISPLAY_NULL) {
	    Tcl_AppendResult(s->interp, "release: ", name, " not found\n", (char *)NULL);
	    return TCL_ERROR;
	}
    } else {
	const char *display_path = mged_display_pathname(s->mged_curr_display);
	if (display_path && BU_STR_EQUAL("nu", display_path))
	    return TCL_OK;  /* Ignore */
    }

    if (mged_variables->mv_listen) {
	/* Drop all clients before detaching the endpoint-owned backend. */
	mged_variables->mv_listen = 0;
	fbserv_set_port(NULL, NULL, NULL, NULL, s);
    }

    /*
     * This saves the state of the resources to the "nu" display
     * manager, which is beneficial only if closing the last display
	 * view. So when another graphics view is opened, it looks
     * like the last one the user had open.
     */
    usurp_all_resources(mged_initial_display, s->mged_curr_display);

    /* If this display is being referenced by a command window, then
     * remove the reference.
     */
    if (s->mged_curr_display->display_tie != NULL)
	s->mged_curr_display->display_tie->cl_tie = (struct mged_display *)NULL;

    if (need_close && s->gedp)
	ged_draw_obol_framebuffer_release(s->gedp);

    if (need_close) {
	mged_obol_display_detach(s, s->mged_curr_display);
    }

    bu_ptbl_rm(&active_display_set, (long *)s->mged_curr_display);
    mged_slider_free_vls(s->mged_curr_display);
    bu_free((void *)s->mged_curr_display, "release: s->mged_curr_display");

    if (save_display != MGED_DISPLAY_NULL)
	mged_current_display_set(s, save_display);
    else {
	if (BU_PTBL_LEN(&active_display_set) > 0) {
	    mged_current_display_set(s, (struct mged_display *)BU_PTBL_GET(&active_display_set, 0));
	} else {
	    mged_current_display_set(s, MGED_DISPLAY_NULL);
	}
    }
    return TCL_OK;
}


int
f_release(ClientData clientData, Tcl_Interp *interpreter, int argc, const char *argv[])
{
    struct cmdtab *ctp = (struct cmdtab *)clientData;
    MGED_CK_CMD(ctp);
    struct mged_state *s = ctp->s;

    struct bu_vls vls = BU_VLS_INIT_ZERO;

    if (argc < 1 || 2 < argc) {
	bu_vls_printf(&vls, "help release");
	Tcl_Eval(interpreter, bu_vls_addr(&vls));
	bu_vls_free(&vls);

	return TCL_ERROR;
    }

    if (argc == 2) {
	int status;

	if (*argv[1] != '.')
	    bu_vls_printf(&vls, ".%s", argv[1]);
	else
	    bu_vls_strcpy(&vls, argv[1]);

	status = release(s, bu_vls_addr(&vls), 1);

	bu_vls_free(&vls);
	return status;
    } else
	return release(s, (char *)NULL, 1);
}


static void
print_valid_dm(Tcl_Interp *interpreter)
{
    Tcl_AppendResult(interpreter,
	    "\tThe following display host types are valid: nu tkobol\n",
	    (char *)NULL);
}


int
f_attach(ClientData clientData, Tcl_Interp *interpreter, int argc, const char *argv[])
{
    struct cmdtab *ctp = (struct cmdtab *)clientData;
    MGED_CK_CMD(ctp);
    struct mged_state *s = ctp->s;

    if (argc < 2) {
	struct bu_vls vls = BU_VLS_INIT_ZERO;

	bu_vls_printf(&vls, "help attach");
	Tcl_Eval(interpreter, bu_vls_addr(&vls));
	bu_vls_free(&vls);
	print_valid_dm(interpreter);

	return TCL_ERROR;
    }

    if (BU_STR_EQUAL(argv[argc-1], "nu")) {
	/* nothing to do */
	return TCL_OK;
    }

    if (!BU_STR_EQUAL(argv[argc-1], "tkobol")) {
	Tcl_AppendResult(interpreter, "attach(", argv[argc - 1], "): BAD\n", (char *)NULL);
	print_valid_dm(interpreter);
	return TCL_ERROR;
    }

    return mged_attach(s, argv[argc - 1], argc, argv);
}


int
gui_setup(struct mged_state *s, const char *dstr)
{
    if (!s)
	return TCL_ERROR;

#ifdef HAVE_TK
    Tk_GenericProc *handler = doEvent;
#endif
    /* initialize only once */
    if (tkwin != NULL)
	return TCL_OK;

    Tcl_ResetResult(s->interp);

    /* set DISPLAY to dstr */
    if (dstr != (char *)NULL) {
	Tcl_SetVar(s->interp, "env(DISPLAY)", dstr, TCL_GLOBAL_ONLY);
    }

#ifdef HAVE_TK
    /* This runs the tk.tcl script */
    if (Tk_Init(s->interp) == TCL_ERROR) {
	const char *result = Tcl_GetStringResult(s->interp);
	/* hack to avoid a stupid Tk error */
	if (bu_strncmp(result, "this isn't a Tk applicationcouldn't", 35) == 0) {
	    result = (result + 27);
	    Tcl_ResetResult(s->interp);
	    Tcl_AppendResult(s->interp, result, (char *)NULL);
	}
	return TCL_ERROR;
    }

    /* Initialize [incr Tk] */
    if (Tcl_Eval(s->interp, "package require Itk") != TCL_OK) {
      return TCL_ERROR;
    }

    /* Import [incr Tk] commands into the global namespace */
    if (Tcl_Import(s->interp, Tcl_GetGlobalNamespace(s->interp),
		   "::itk::*", /* allowOverwrite */ 1) != TCL_OK) {
	return TCL_ERROR;
    }
#endif

    /* Initialize the Iwidgets package */
    if (Tcl_Eval(s->interp, "package require Iwidgets") != TCL_OK) {
	return TCL_ERROR;
    }

    /* Import iwidgets into the global namespace */
    if (Tcl_Import(s->interp, Tcl_GetGlobalNamespace(s->interp),
		   "::iwidgets::*", /* allowOverwrite */ 1) != TCL_OK) {
	return TCL_ERROR;
    }

    /* Initialize libtclcad */
#ifdef HAVE_TK
    (void)tclcad_init(s->interp, 1, NULL);
#else
    (void)tclcad_init(s->interp, 0, NULL);
#endif

#ifdef HAVE_TK
    if ((tkwin = Tk_MainWindow(s->interp)) == NULL) {
	return TCL_ERROR;
    }

    /* create the event handler */
    Tk_CreateGenericHandler(handler, (ClientData)s);

    Tcl_Eval(s->interp, "wm withdraw .");
    Tcl_Eval(s->interp, "tk appname mged");
#endif

    return TCL_OK;
}


int
mged_attach(struct mged_state *s, const char *wp_name, int argc, const char *argv[])
{
    struct mged_display *o_dm;

    if (!wp_name) {
	return TCL_ERROR;
    }

    o_dm = s->mged_curr_display;
    BU_ALLOC(s->mged_curr_display, struct mged_display);

    /* Only need to do this once */
    if (tkwin == NULL && BU_STR_EQUAL(wp_name, "tkobol")) {
	if (gui_setup(s, mged_tkobol_display_arg(argc, argv)) == TCL_ERROR) {
	    bu_free((void *)s->mged_curr_display, "f_attach: current display");
	    mged_current_display_set(s, o_dm);
	    return TCL_ERROR;
	}
    }

    bu_ptbl_ins(&active_display_set, (long *)s->mged_curr_display);

    if (!wp_name) {
	return TCL_ERROR;
    }

    if (mged_display_init(s, o_dm, wp_name, argc, argv) == TCL_ERROR) {
	goto Bad;
    }

    /* initialize the background color */
    {
	/* need dummy values for func signature--they are unused in the func */
	const struct bu_structparse *sdp = 0;
	const char name[] = "name";
	void *base = 0;
	const char value[] = "value";
	cs_set_bg(sdp, name, base, value, s);
    }

    mged_link_vars(s->mged_curr_display);

    Tcl_ResetResult(s->interp);
    Tcl_AppendResult(s->interp, "ATTACHING tkobol (Tk Obol display endpoint)\n",
	    (char *)NULL);

    (void)mged_obol_framebuffer_ensure(s);

    ged_view_active_ctx_set(s->gedp, s->mged_curr_display->display_view_state->vs_gvp);

    return TCL_OK;

 Bad:
    Tcl_AppendResult(s->interp, "attach(", argv[argc - 1], "): BAD\n", (char *)NULL);

    release(s, (char *)NULL, 1);

    return TCL_ERROR;
}


#define MAX_ATTACH_RETRIES 100

void
get_attached(struct mged_state *s)
{
    int inflimit = MAX_ATTACH_RETRIES;
    int ret;
    struct bu_vls wanted_type = BU_VLS_INIT_ZERO;
    struct bu_vls prompt = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&prompt, "attach (nu tkobol)[nu]? ");

    while (inflimit > 0) {
	bu_log("%s", bu_vls_cstr(&prompt));

	ret = bu_vls_gets(&wanted_type, stdin);
	if (ret < 0) {
	    /* handle EOF */
	    bu_log("\n");
	    bu_vls_free(&wanted_type);
	    bu_vls_free(&prompt);
	    return;
	}

	if (bu_vls_strlen(&wanted_type) == 0 || BU_STR_EQUAL(bu_vls_addr(&wanted_type), "nu")) {
	    /* Nothing more to do. */
	    bu_vls_free(&wanted_type);
	    bu_vls_free(&prompt);
	    return;
	}

	/* trim whitespace before comparisons (but not before checking empty) */
	bu_vls_trimspace(&wanted_type);

	if (BU_STR_EQUAL(bu_vls_cstr(&wanted_type), "tkobol")) {
	    break;
	}

	/* Not a valid choice, loop. */
	inflimit--;
    }

    bu_vls_free(&prompt);

    if (inflimit <= 0) {
	bu_log("\nInfinite loop protection, attach aborted!\n");
	bu_vls_free(&wanted_type);
	return;
    }

    bu_log("Starting an %s graphics endpoint\n", bu_vls_cstr(&wanted_type));

    int argc = 1;
    const char *argv[3];
    argv[0] = "";
    argv[1] = "";
    argv[2] = (char *)NULL;
    (void)mged_attach(s, bu_vls_cstr(&wanted_type), argc, argv);
    bu_vls_free(&wanted_type);
}


/*
 * Run graphics endpoint-specific command(s).
 */
int
f_dm(ClientData clientData, Tcl_Interp *interpreter, int argc, const char *argv[])
{
    struct cmdtab *ctp = (struct cmdtab *)clientData;
    MGED_CK_CMD(ctp);
    struct mged_state *s = ctp->s;
    struct bu_vls vls = BU_VLS_INIT_ZERO;

    if (argc < 2) {
	bu_vls_printf(&vls, "help dm");
	Tcl_Eval(interpreter, bu_vls_addr(&vls));
	bu_vls_free(&vls);
	return TCL_ERROR;
    }

    if (BU_STR_EQUAL(argv[1], "valid")) {
	if (argc < 3) {
	    bu_vls_printf(&vls, "help dm");
	    Tcl_Eval(interpreter, bu_vls_addr(&vls));
	    bu_vls_free(&vls);
	    return TCL_ERROR;
	}
	if (BU_STR_EQUAL(argv[argc-1], "nu") ||
	    BU_STR_EQUAL(argv[argc-1], "tkobol")) {
	    Tcl_AppendResult(interpreter, argv[argc-1], (char *)NULL);
	}
	return TCL_OK;
    }

    if (!cmd_hook) {
	Tcl_AppendResult(interpreter,
		"The current display host does not support local commands.\n",
		(char *)NULL);
	return TCL_ERROR;
    }


    return cmd_hook(argc-1, argv+1, (void *)s);
}

void
mged_display_var_init(struct mged_state *s, struct mged_display *target_display)
{
    s->mged_curr_display->display_input_state = NULL;
    s->mged_curr_display->display_input_motion_pending = 0;
    s->mged_curr_display->display_input_snap_pan_active = 0;
    bu_vls_init(&s->mged_curr_display->display_pathname);
    BU_ALLOC(menu_state, struct _menu_state);
    *menu_state = *target_display->display_menu_state;		/* struct copy */
    menu_state->ms_rc = 1;

    BU_ALLOC(rubber_band, struct _rubber_band);
    *rubber_band = *target_display->display_rubber_band;		/* struct copy */
    rubber_band->rb_rc = 1;

    BU_ALLOC(mged_variables, struct _mged_variables);
    *mged_variables = *target_display->display_variables;	/* struct copy */
    mged_variables->mv_rc = 1;
    mged_variables->mv_listen = 0;
    mged_variables->mv_port = 0;
    mged_variables->mv_fb = 0;

    BU_ALLOC(color_scheme, struct _color_scheme);

    /* Initialize from the initial graphics view. */
    if (mged_initial_display && mged_initial_display->display_color_scheme) {
	*color_scheme = *mged_initial_display->display_color_scheme;
    }

    color_scheme->cs_rc = 1;
    s->mged_curr_display->display_color_scheme_dirty = 1;

    BU_ALLOC(axes_state, struct _axes_state);
    *axes_state = *target_display->display_axes_state;		/* struct copy */
    axes_state->ax_rc = 1;
    s->mged_curr_display->display_axes_state_dirty = 1;
    s->mged_curr_display->display_adc_style_dirty = 1;

    BU_ALLOC(view_state, struct _view_state);
    *view_state = *target_display->display_view_state;			/* struct copy */
    void *target_view_ctx = target_display->display_view_state->vs_gvp;
    void *view_set_ctx = ged_view_set_ctx(s->gedp);
    void *view_ctx = ged_view_context_create_copy_with_set(target_view_ctx, view_set_ctx);
    view_state->vs_gvp = view_ctx;

    ged_view_set_context_add(view_set_ctx, view_ctx);
    ged_view_context_owned_add(s->gedp, view_ctx);
    ged_view_context_update_callback_set(view_ctx,
	    mged_view_callback, (void *)view_state);
    ged_draw_view_lod_policy lod_policy = BV_LOD_POLICY_INIT;
    if (ged_draw_view_context_lod_policy_get(&lod_policy, view_ctx)) {
	lod_policy.csg_enabled = 0;
	lod_policy.zoom_refresh = 0;
	lod_policy.point_scale = 1.0;
	lod_policy.curve_scale = 1.0;
	ged_draw_view_context_lod_policy_apply(view_ctx, &lod_policy);
    }
    view_state->vs_rc = 1;
    view_ring_init(s->mged_curr_display->display_view_state, (struct _view_state *)NULL);

    mged_display_repaint_request(target_display, MGED_REPAINT_NATIVE_EVENT);
    mapped = 1;
    s->mged_curr_display->display_owner = 1;
    am_mode = AMM_IDLE;
    adc_auto = 1;
    grid_auto_size = 1;
}


void
mged_link_vars(struct mged_display *p)
{
    mged_slider_init_vls(p);
    const char *pn = mged_display_pathname(p);
    if (pn) {
	bu_vls_printf(&p->display_fps_name, "%s(%s,fps)", MGED_DISPLAY_VAR, pn);
	bu_vls_printf(&p->display_aet_name, "%s(%s,aet)", MGED_DISPLAY_VAR, pn);
	bu_vls_printf(&p->display_ang_name, "%s(%s,ang)", MGED_DISPLAY_VAR, pn);
	bu_vls_printf(&p->display_center_name, "%s(%s,center)", MGED_DISPLAY_VAR, pn);
	bu_vls_printf(&p->display_size_name, "%s(%s,size)", MGED_DISPLAY_VAR, pn);
	bu_vls_printf(&p->display_adc_name, "%s(%s,adc)", MGED_DISPLAY_VAR, pn);
    }
}


int
f_get_display_list(ClientData UNUSED(clientData), Tcl_Interp *interpreter, int argc, const char *argv[])
{
    if (argc != 1 || !argv) {
	struct bu_vls vls = BU_VLS_INIT_ZERO;

	bu_vls_printf(&vls, "helpdevel get_display_list");
	Tcl_Eval(interpreter, bu_vls_addr(&vls));
	bu_vls_free(&vls);

	return TCL_ERROR;
    }

    for (size_t i = 0; i < BU_PTBL_LEN(&active_display_set); i++) {
	struct mged_display *dlp = (struct mged_display *)BU_PTBL_GET(&active_display_set, i);
	const char *pn = mged_display_pathname(dlp);
	if (pn)
	    Tcl_AppendElement(interpreter, pn);
    }
    return TCL_OK;
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
