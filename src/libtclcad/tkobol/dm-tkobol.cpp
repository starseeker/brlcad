/*                      D M - T K O B O L . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file libtclcad/tkobol/dm-tkobol.cpp
 *
 * Tk-native host for an Obol view.  Togl owns the Tk window and OpenGL
 * context; this host owns the retained controller directly.  The struct dm
 * surface is a compatibility attachment for classic Tcl applications.
 */

#include "common.h"

#include <cstring>

#include "tcl.h"
#include "tk.h"

#ifdef TKOBOL_X11
#  include <X11/Xlib.h>
#  include <GL/gl.h>
#  include <GL/glx.h>
#endif

extern "C" {
#include "bu/log.h"
#include "bu/malloc.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "dm.h"
#include "../../libdm/include/private.h"
#include "../../libdm/null/dm-Null.h"
}

#include "brlobol/init.h"
#include "brlobol/scene_group.h"
#include "brlobol/view_controller.h"
#include <Inventor/SoRenderManager.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include "vendor/togl/togl.h"

enum tkobol_callback_kind {
    TKOBOL_DISPLAY_CALLBACK = 1,
    TKOBOL_RESHAPE_CALLBACK = 2
};

struct tkobol_vars;

struct tkobol_callback {
    struct tkobol_vars *vars;
    enum tkobol_callback_kind kind;
};

struct tkobol_vars {
    struct dm *dmp;
    Tcl_Interp *interp;
    Tk_Window tkwin;
    Tk_Window container;
    BRLObolViewController *controller;
    struct bu_vls widget_path;
    struct bu_vls widget_command;
    struct bu_vls container_command;
    struct bu_vls display_command;
    struct bu_vls reshape_command;
    struct tkobol_callback display_callback;
    struct tkobol_callback reshape_callback;
    int widget_command_private;
    int closing;
    fastf_t viewscale;
    int (*original_close)(struct dm *);
};

static int tkobol_close(struct dm *dmp);
static int tkobol_configure(struct dm *dmp, int force);
static int tkobol_draw_end(struct dm *dmp);
static int tkobol_get_display_image(struct dm *dmp, unsigned char **image,
	int flip, int alpha);
static int tkobol_make_current(struct dm *dmp);
static int tkobol_reshape(struct dm *dmp, int width, int height);
static int tkobol_swap_buffers(struct dm *dmp);

static struct tkobol_vars *
tkobol_pvars(struct dm *dmp)
{
    return (dmp && dmp->i) ?
	static_cast<struct tkobol_vars *>(dmp->i->dm_udata) : NULL;
}

static int
tkobol_widget_command(struct tkobol_vars *tv, const char *command)
{
    if (!tv || !tv->interp || !tv->dmp || !command || tv->closing)
	return TCL_ERROR;

    Tcl_Obj *objv[2];
    objv[0] = Tcl_NewStringObj(tv->widget_command_private ?
	bu_vls_cstr(&tv->widget_command) : bu_vls_cstr(&tv->widget_path), -1);
    objv[1] = Tcl_NewStringObj(command, -1);
    Tcl_IncrRefCount(objv[0]);
    Tcl_IncrRefCount(objv[1]);

    Tcl_InterpState state = Tcl_SaveInterpState(tv->interp, TCL_OK);
    int ret = Tcl_EvalObjv(tv->interp, 2, objv, TCL_EVAL_GLOBAL);
    if (ret != TCL_OK)
	bu_log("tkobol: %s failed for %s: %s\n", command,
		bu_vls_cstr(&tv->dmp->i->dm_pathName),
		Tcl_GetStringResult(tv->interp));
    Tcl_RestoreInterpState(tv->interp, state);

    Tcl_DecrRefCount(objv[1]);
    Tcl_DecrRefCount(objv[0]);
    return ret;
}

static int
tkobol_rename_command(Tcl_Interp *interp, const char *from, const char *to)
{
    Tcl_Obj *objv[3];
    const char *argv[] = {"rename", from, to};
    for (int i = 0; i < 3; i++) {
	objv[i] = Tcl_NewStringObj(argv[i], -1);
	Tcl_IncrRefCount(objv[i]);
    }
    const int ret = Tcl_EvalObjv(interp, 3, objv, TCL_EVAL_GLOBAL);
    for (int i = 2; i >= 0; i--)
	Tcl_DecrRefCount(objv[i]);
    return ret;
}

static int
tkobol_sync_view(struct tkobol_vars *tv, int request_render)
{
    if (!tv || !tv->dmp || !tv->controller)
	return BRLCAD_ERROR;

    int width = tv->tkwin ? Tk_Width(tv->tkwin) : 0;
    int height = tv->tkwin ? Tk_Height(tv->tkwin) : 0;
    if (width < 2)
	width = tv->dmp->i->dm_width > 1 ? tv->dmp->i->dm_width : 512;
    if (height < 2)
	height = tv->dmp->i->dm_height > 1 ? tv->dmp->i->dm_height : 512;

    const int resized = width != tv->dmp->i->dm_width ||
	height != tv->dmp->i->dm_height;
    tv->dmp->i->dm_width = width;
    tv->dmp->i->dm_height = height;
    tv->dmp->i->dm_aspect = static_cast<fastf_t>(width) /
	static_cast<fastf_t>(height);
    if (tv->dmp->i->dm_ctx)
	dm_view_context_dimensions_set(tv->dmp->i->dm_ctx, width, height);

    tv->controller->setViewportSize(static_cast<unsigned int>(width),
	static_cast<unsigned int>(height));
    if (tv->dmp->i->dm_ctx &&
	!tv->controller->syncCameraFromViewContext(tv->dmp->i->dm_ctx))
	return BRLCAD_ERROR;

    if (request_render || resized) {
	tv->controller->requestRender(resized ? "Tk Obol resize" :
		"Tk Obol refresh");
	dm_set_native_repaint_pending(tv->dmp, 1);
	if (tv->dmp->i->dm_ctx)
	    dm_view_context_refresh_request(tv->dmp->i->dm_ctx,
		    DM_VIEW_REFRESH_VIEW | DM_VIEW_REFRESH_FORCE);
    }
    return BRLCAD_OK;
}

static int
tkobol_render_frame(struct tkobol_vars *tv, GLenum render_buffer,
	GLboolean double_buffered)
{
    if (!tv || !tv->controller)
	return BRLCAD_ERROR;

    tv->controller->getRenderManager()->setDoubleBuffer(double_buffered ?
	TRUE : FALSE);
    glDrawBuffer(render_buffer);
    (void)tv->controller->realizePending();
    (void)tv->controller->advanceProgressiveWork();
    glClearColor(static_cast<GLfloat>(tv->dmp->i->dm_bg1[0]) / 255.0f,
	static_cast<GLfloat>(tv->dmp->i->dm_bg1[1]) / 255.0f,
	static_cast<GLfloat>(tv->dmp->i->dm_bg1[2]) / 255.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (tv->controller->getCamera() && tv->controller->getRenderRoot()) {
	SoGLRenderAction *action =
	    tv->controller->getRenderManager()->getGLRenderAction();
	action->setViewportRegion(tv->controller->getViewportRegion());
	action->apply(tv->controller->getRenderRoot());
    }
    tv->controller->clearRenderRequest();
    return BRLCAD_OK;
}

static int
tkobol_render_current(struct tkobol_vars *tv, int make_current)
{
    if (!tv || !tv->controller)
	return BRLCAD_ERROR;
    if (make_current && tkobol_widget_command(tv, "makecurrent") != TCL_OK)
	return BRLCAD_ERROR;
    if (tkobol_sync_view(tv, 0) != BRLCAD_OK)
	return BRLCAD_ERROR;

    GLboolean double_buffered = GL_FALSE;
    glGetBooleanv(GL_DOUBLEBUFFER, &double_buffered);
    if (tkobol_render_frame(tv, double_buffered ? GL_BACK : GL_FRONT,
		double_buffered) != BRLCAD_OK)
	return BRLCAD_ERROR;

    if (double_buffered &&
	    tkobol_widget_command(tv, "swapbuffers") != TCL_OK)
	return BRLCAD_ERROR;

    const int pending = tv->controller->isRenderRequested() ||
	tv->controller->hasProgressiveWorkPending();
    dm_set_native_repaint_pending(tv->dmp, pending);
    if (pending && tv->dmp->i->dm_ctx)
	dm_view_context_refresh_request(tv->dmp->i->dm_ctx,
		DM_VIEW_REFRESH_VIEW);
    return BRLCAD_OK;
}

static int
tkobol_callback_command(ClientData client_data, Tcl_Interp *interp,
	int objc, Tcl_Obj *const *objv)
{
    struct tkobol_callback *callback =
	static_cast<struct tkobol_callback *>(client_data);
    if (!callback || !callback->vars || callback->vars->closing)
	return TCL_OK;
    if (objc != 2) {
	Tcl_WrongNumArgs(interp, 1, objv, "widget");
	return TCL_ERROR;
    }

    if (callback->kind == TKOBOL_RESHAPE_CALLBACK) {
	if (tkobol_sync_view(callback->vars, 1) != BRLCAD_OK)
	    return TCL_ERROR;
	return tkobol_widget_command(callback->vars, "postredisplay");
    }

    return tkobol_render_current(callback->vars, 0) == BRLCAD_OK ?
	TCL_OK : TCL_ERROR;
}

static int
tkobol_create_widget(struct tkobol_vars *tv)
{
    if (!tv || !tv->dmp || !tv->interp)
	return TCL_ERROR;

    static unsigned long callback_count = 0;
    const unsigned long id = ++callback_count;
    bu_vls_sprintf(&tv->display_command,
	"::brlcad::tkobol::display_%lu", id);
    bu_vls_sprintf(&tv->reshape_command,
	"::brlcad::tkobol::reshape_%lu", id);
    bu_vls_sprintf(&tv->widget_command,
	"::brlcad::tkobol::widget_%lu", id);
    bu_vls_sprintf(&tv->container_command,
	"::brlcad::tkobol::container_%lu", id);
    tv->display_callback.vars = tv;
    tv->display_callback.kind = TKOBOL_DISPLAY_CALLBACK;
    tv->reshape_callback.vars = tv;
    tv->reshape_callback.kind = TKOBOL_RESHAPE_CALLBACK;

    Tcl_CreateObjCommand(tv->interp, bu_vls_cstr(&tv->display_command),
	tkobol_callback_command, &tv->display_callback, NULL);
    Tcl_CreateObjCommand(tv->interp, bu_vls_cstr(&tv->reshape_command),
	tkobol_callback_command, &tv->reshape_callback, NULL);

    bu_vls_sprintf(&tv->widget_path, "%s",
	bu_vls_cstr(&tv->dmp->i->dm_pathName));
    if (tv->dmp->i->dm_top) {
	const char *top_argv[] = {
	    "toplevel", bu_vls_cstr(&tv->dmp->i->dm_pathName)
	};
	Tcl_Obj *top_objv[2];
	for (int i = 0; i < 2; i++) {
	    top_objv[i] = Tcl_NewStringObj(top_argv[i], -1);
	    Tcl_IncrRefCount(top_objv[i]);
	}
	const int top_ret = Tcl_EvalObjv(tv->interp, 2, top_objv,
		TCL_EVAL_GLOBAL);
	for (int i = 1; i >= 0; i--)
	    Tcl_DecrRefCount(top_objv[i]);
	if (top_ret != TCL_OK) {
	    bu_log("tkobol: unable to create toplevel %s: %s\n",
		    bu_vls_cstr(&tv->dmp->i->dm_pathName),
		    Tcl_GetStringResult(tv->interp));
	    return TCL_ERROR;
	}
	tv->container = Tk_NameToWindow(tv->interp,
		bu_vls_cstr(&tv->dmp->i->dm_pathName),
		Tk_MainWindow(tv->interp));
	if (!tv->container)
	    return TCL_ERROR;
	if (tkobol_rename_command(tv->interp,
		bu_vls_cstr(&tv->dmp->i->dm_pathName),
		bu_vls_cstr(&tv->container_command)) != TCL_OK)
	    return TCL_ERROR;
	bu_vls_printf(&tv->widget_path, ".__obol");
    }

    char width[32];
    char height[32];
    snprintf(width, sizeof(width), "%d", tv->dmp->i->dm_width);
    snprintf(height, sizeof(height), "%d", tv->dmp->i->dm_height);
    const char *argv[] = {
	"::brlcad::tkobol::host",
	bu_vls_cstr(&tv->widget_path),
	"-width", width,
	"-height", height,
	"-rgba", "true",
	"-double", "true",
	"-depth", "true",
	"-depthsize", "24",
	"-major", "2",
	"-minor", "1",
	"-coreprofile", "false",
	"-swapinterval", "1",
	"-displayproc", bu_vls_cstr(&tv->display_command),
	"-reshapeproc", bu_vls_cstr(&tv->reshape_command)
    };
    constexpr int argc = static_cast<int>(sizeof(argv) / sizeof(argv[0]));
    Tcl_Obj *objv[argc];
    for (int i = 0; i < argc; i++) {
	objv[i] = Tcl_NewStringObj(argv[i], -1);
	Tcl_IncrRefCount(objv[i]);
    }
    const int ret = Tcl_EvalObjv(tv->interp, argc, objv, TCL_EVAL_GLOBAL);
    for (int i = argc - 1; i >= 0; i--)
	Tcl_DecrRefCount(objv[i]);
    if (ret != TCL_OK) {
	bu_log("tkobol: unable to create %s: %s\n",
		bu_vls_cstr(&tv->dmp->i->dm_pathName),
		Tcl_GetStringResult(tv->interp));
	return TCL_ERROR;
    }
    if (tkobol_rename_command(tv->interp, bu_vls_cstr(&tv->widget_path),
	    bu_vls_cstr(&tv->widget_command)) != TCL_OK)
	return TCL_ERROR;
    tv->widget_command_private = 1;

    tv->tkwin = Tk_NameToWindow(tv->interp,
	bu_vls_cstr(&tv->widget_path), Tk_MainWindow(tv->interp));
    if (!tv->tkwin)
	return TCL_ERROR;
    if (!tv->container)
	tv->container = tv->tkwin;

    if (tv->dmp->i->dm_top) {
	const char *pack_argv[] = {
	    "pack", bu_vls_cstr(&tv->widget_path),
	    "-expand", "true", "-fill", "both"
	};
	Tcl_Obj *pack_objv[6];
	for (int i = 0; i < 6; i++) {
	    pack_objv[i] = Tcl_NewStringObj(pack_argv[i], -1);
	    Tcl_IncrRefCount(pack_objv[i]);
	}
	const int pack_ret = Tcl_EvalObjv(tv->interp, 6, pack_objv,
		TCL_EVAL_GLOBAL);
	for (int i = 5; i >= 0; i--)
	    Tcl_DecrRefCount(pack_objv[i]);
	if (pack_ret != TCL_OK)
	    return TCL_ERROR;
    }

    Tk_MapWindow(tv->tkwin);

    /* Windows ignores an early deiconify before the native window exists. */
    if (tv->dmp->i->dm_top) {
	Tk_MakeWindowExist(tv->container);
	Tk_MapWindow(tv->container);
	const char *wm_argv[] = {
	    "wm", "deiconify", bu_vls_cstr(&tv->dmp->i->dm_pathName)
	};
	Tcl_Obj *wm_objv[3];
	for (int i = 0; i < 3; i++) {
	    wm_objv[i] = Tcl_NewStringObj(wm_argv[i], -1);
	    Tcl_IncrRefCount(wm_objv[i]);
	}
	Tcl_InterpState state = Tcl_SaveInterpState(tv->interp, TCL_OK);
	const int wm_ret = Tcl_EvalObjv(tv->interp, 3, wm_objv,
		TCL_EVAL_GLOBAL);
	Tcl_RestoreInterpState(tv->interp, state);
	for (int i = 2; i >= 0; i--)
	    Tcl_DecrRefCount(wm_objv[i]);
	if (wm_ret != TCL_OK)
	    return TCL_ERROR;
    }

    tv->dmp->i->dm_id = Tk_WindowId(tv->tkwin);
    bu_vls_sprintf(&tv->dmp->i->dm_tkName, "%s", Tk_Name(tv->tkwin));
    return TCL_OK;
}

static int
tkobol_close(struct dm *dmp)
{
    struct tkobol_vars *tv = tkobol_pvars(dmp);
    int (*original_close)(struct dm *) = NULL;
    if (tv) {
	tv->closing = 1;
	original_close = tv->original_close;
	if (tv->container)
	    Tk_DestroyWindow(tv->container);
	if (tv->interp) {
	    Tcl_DeleteCommand(tv->interp,
		    bu_vls_cstr(&tv->display_command));
	    Tcl_DeleteCommand(tv->interp,
		    bu_vls_cstr(&tv->reshape_command));
	}
	bu_vls_free(&tv->display_command);
	bu_vls_free(&tv->reshape_command);
	bu_vls_free(&tv->widget_path);
	bu_vls_free(&tv->widget_command);
	bu_vls_free(&tv->container_command);
	delete tv->controller;
	dmp->i->dm_udata = NULL;
	bu_free(tv, "tkobol private vars");
    }
    return original_close ? original_close(dmp) : BRLCAD_OK;
}

static int
tkobol_configure(struct dm *dmp, int force)
{
    struct tkobol_vars *tv = tkobol_pvars(dmp);
    if (!tv)
	return BRLCAD_ERROR;
    int ret = tkobol_sync_view(tv, force ? 1 : 0);
    if (ret == BRLCAD_OK)
	(void)tkobol_widget_command(tv, "postredisplay");
    return ret;
}

static int
tkobol_draw_end(struct dm *dmp)
{
    return tkobol_render_current(tkobol_pvars(dmp), 1);
}

static int
tkobol_get_display_image(struct dm *dmp, unsigned char **image, int flip,
	int alpha)
{
    struct tkobol_vars *tv = tkobol_pvars(dmp);
    if (!tv || !tv->controller || !image)
	return BRLCAD_ERROR;
    *image = NULL;

    if (tkobol_widget_command(tv, "makecurrent") != TCL_OK ||
	    tkobol_sync_view(tv, 0) != BRLCAD_OK)
	return BRLCAD_ERROR;

    GLboolean double_buffered = GL_FALSE;
    glGetBooleanv(GL_DOUBLEBUFFER, &double_buffered);
    const GLenum render_buffer = double_buffered ? GL_BACK : GL_FRONT;
    if (tkobol_render_frame(tv, render_buffer, double_buffered) != BRLCAD_OK)
	return BRLCAD_ERROR;

    const int width = dmp->i->dm_width;
    const int height = dmp->i->dm_height;
    const int components = alpha ? 4 : 3;
    if (width <= 0 || height <= 0)
	return BRLCAD_ERROR;

    GLint old_alignment = 4;
    GLint old_read_buffer = GL_BACK;
    glGetIntegerv(GL_PACK_ALIGNMENT, &old_alignment);
    glGetIntegerv(GL_READ_BUFFER, &old_read_buffer);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadBuffer(render_buffer);

    const size_t stride = static_cast<size_t>(width) *
	static_cast<size_t>(components);
    unsigned char *pixels = static_cast<unsigned char *>(bu_malloc(
	static_cast<size_t>(height) * stride, "Tk Obol display image"));
    glReadPixels(0, 0, width, height, alpha ? GL_RGBA : GL_RGB,
	GL_UNSIGNED_BYTE, pixels);
    glFinish();
    glReadBuffer(old_read_buffer);
    glPixelStorei(GL_PACK_ALIGNMENT, old_alignment);

    if (flip) {
	unsigned char *row = static_cast<unsigned char *>(
	    bu_malloc(stride, "Tk Obol display image row"));
	for (int y = 0; y < height / 2; y++) {
	    unsigned char *lower = pixels + static_cast<size_t>(y) * stride;
	    unsigned char *upper = pixels +
		static_cast<size_t>(height - y - 1) * stride;
	    std::memcpy(row, lower, stride);
	    std::memcpy(lower, upper, stride);
	    std::memcpy(upper, row, stride);
	}
	bu_free(row, "Tk Obol display image row");
    }

    *image = pixels;
    return BRLCAD_OK;
}

static int
tkobol_make_current(struct dm *dmp)
{
    return tkobol_widget_command(tkobol_pvars(dmp), "makecurrent") == TCL_OK ?
	BRLCAD_OK : BRLCAD_ERROR;
}

static int
tkobol_reshape(struct dm *dmp, int width, int height)
{
    struct tkobol_vars *tv = tkobol_pvars(dmp);
    if (!tv || width <= 0 || height <= 0)
	return BRLCAD_ERROR;
    dmp->i->dm_width = width;
    dmp->i->dm_height = height;
    if (tv->tkwin)
	Tk_GeometryRequest(tv->tkwin, width, height);
    return tkobol_sync_view(tv, 1);
}

static int
tkobol_swap_buffers(struct dm *dmp)
{
    return tkobol_widget_command(tkobol_pvars(dmp), "swapbuffers") == TCL_OK ?
	BRLCAD_OK : BRLCAD_ERROR;
}

static int
tkobol_doevent(struct dm *dmp, void *UNUSED(client_data),
	void *UNUSED(event_ptr))
{
    struct tkobol_vars *tv = tkobol_pvars(dmp);
    if (!tv)
	return TCL_ERROR;
    (void)tkobol_sync_view(tv, 1);
    return tkobol_widget_command(tv, "postredisplay");
}

static int
tkobol_viable(const char *display_name)
{
#ifdef TKOBOL_X11
    Display *display = XOpenDisplay(display_name && display_name[0] ?
	display_name : NULL);
    if (!display)
	return 0;
    int error_base = 0;
    int event_base = 0;
    const int viable = glXQueryExtension(display, &error_base, &event_base) ?
	1 : 0;
    XCloseDisplay(display);
    return viable;
#else
    (void)display_name;
    return 1;
#endif
}

static struct dm *
tkobol_open(void *ctx, void *vinterp, int argc, const char **argv)
{
    static int count = 0;
    Tcl_Interp *interp = static_cast<Tcl_Interp *>(vinterp);
    if (!interp || !Tk_MainWindow(interp))
	return DM_NULL;
    if (BrlcadTkObolHost_Init(interp) != TCL_OK)
	return DM_NULL;

    struct dm *dmp = NULL;
    struct dm_impl *dmpi = NULL;
    BU_GET(dmp, struct dm);
    BU_GET(dmpi, struct dm_impl);
    dmp->magic = DM_MAGIC;
    dmp->start_time = 0;
    *dmpi = *dm_null.i;
    dmp->i = dmpi;
    dmp->i->dm_ctx = ctx;
    dmp->i->dm_interp = vinterp;
    dmp->i->dm_lineWidth = 1;
    dmp->i->dm_depthMask = 1;
    dmp->i->dm_top = 1;
    dmp->i->dm_bytes_per_pixel = 3;
    dmp->i->dm_bits_per_channel = 8;
    dmp->i->dm_fontsize = 20;
    dmp->i->dm_bg1[0] = dmp->i->dm_bg1[1] = dmp->i->dm_bg1[2] = 0;
    dmp->i->dm_fg[0] = dmp->i->dm_fg[1] = dmp->i->dm_fg[2] = 255;
    bu_vls_init(&dmp->i->dm_pathName);
    bu_vls_init(&dmp->i->dm_tkName);
    bu_vls_init(&dmp->i->dm_dName);
    bu_vls_init(&dmp->i->dm_log);

    struct bu_vls init_proc = BU_VLS_INIT_ZERO;
    if (argc > 0)
	dm_processOptions(dmp, &init_proc, argc - 1, argv + 1);
    bu_vls_free(&init_proc);
    const int view_width = ctx ? dm_view_context_width_get(ctx) : 0;
    const int view_height = ctx ? dm_view_context_height_get(ctx) : 0;
    if (dmp->i->dm_width <= 0)
	dmp->i->dm_width = view_width > 0 ? view_width : 512;
    if (dmp->i->dm_height <= 0)
	dmp->i->dm_height = view_height > 0 ? view_height : 512;
    dmp->i->dm_aspect = static_cast<fastf_t>(dmp->i->dm_width) /
	static_cast<fastf_t>(dmp->i->dm_height);
    if (bu_vls_strlen(&dmp->i->dm_pathName) == 0)
	bu_vls_printf(&dmp->i->dm_pathName, ".dm_tkobol%d", count);
    if (bu_vls_strlen(&dmp->i->dm_dName) == 0)
	bu_vls_strcpy(&dmp->i->dm_dName, "tkobol");
    ++count;

    struct tkobol_vars *tv = NULL;
    BU_ALLOC(tv, struct tkobol_vars);
    tv->dmp = dmp;
    tv->interp = interp;
    brlobol_init(brlobol_headless_context_manager());
    tv->controller = new BRLObolViewController(new SoBRLSceneGroup);
    tv->viewscale = 1.0;
    tv->original_close = null_close;
    bu_vls_init(&tv->widget_path);
    bu_vls_init(&tv->widget_command);
    bu_vls_init(&tv->container_command);
    bu_vls_init(&tv->display_command);
    bu_vls_init(&tv->reshape_command);
    dmp->i->dm_udata = tv;
    dmp->i->p_vars = tv->controller;
    dmp->i->m_vars = tv;
    dmp->i->dm_vp = &tv->viewscale;

    if (!tv->controller || tkobol_create_widget(tv) != TCL_OK) {
	tkobol_close(dmp);
	return DM_NULL;
    }

    dmp->i->dm_close = tkobol_close;
    dmp->i->dm_viable = tkobol_viable;
    dmp->i->dm_drawEnd = tkobol_draw_end;
    dmp->i->dm_getDisplayImage = tkobol_get_display_image;
    dmp->i->dm_configureWin = tkobol_configure;
    dmp->i->dm_reshape = tkobol_reshape;
    dmp->i->dm_makeCurrent = tkobol_make_current;
    dmp->i->dm_SwapBuffers = tkobol_swap_buffers;
    dmp->i->dm_doevent = tkobol_doevent;
    dmp->i->dm_graphical = 1;
    dmp->i->graphics_system = "Tk/OpenGL";
    dmp->i->dm_name = "tkobol";
    dmp->i->dm_lname = "Tk Obol graphics";
    (void)tkobol_sync_view(tv, 1);
    return dmp;
}

static struct dm_impl dm_tkobol_impl;
static struct dm dm_tkobol = {DM_MAGIC, &dm_tkobol_impl, 0};

extern "C" void *
tclcad_tkobol_controller(struct dm *dmp)
{
    if (!dmp || !BU_STR_EQUIV(dm_get_type(dmp), "tkobol"))
	return NULL;
    struct tkobol_vars *tv = tkobol_pvars(dmp);
    return tv ? static_cast<void *>(tv->controller) : NULL;
}

extern "C" const struct dm *
tclcad_tkobol_dm(void)
{
    dm_tkobol_impl.dm_open = tkobol_open;
    dm_tkobol_impl.dm_close = tkobol_close;
    dm_tkobol_impl.dm_viable = tkobol_viable;
    dm_tkobol_impl.dm_graphical = 1;
    dm_tkobol_impl.dm_getDisplayImage = tkobol_get_display_image;
    dm_tkobol_impl.graphics_system = "Tk/OpenGL";
    dm_tkobol_impl.dm_name = "tkobol";
    dm_tkobol_impl.dm_lname = "Tk Obol graphics";
    return &dm_tkobol;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
