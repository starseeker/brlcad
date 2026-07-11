/*                       D M - O B O L . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file libdm/dm-obol.cpp
 *
 * Headless Obol/Coin display manager.
 *
 * This backend is intentionally retained-scene first: libged mirrors draw
 * transactions into the backend-owned BRLObolViewController, and libdm renders
 * that Obol scene on demand for screengrab/classic refresh.  The historical
 * immediate-mode libdm drawing callbacks are no-ops here because translating
 * those legacy streams back into Obol would reintroduce the BSG-era adapter
 * path this backend is meant to retire.
 */

#include "common.h"

#include <stddef.h>
#include <stdio.h>
#include <cstring>

extern "C" {
#include "vmath.h"
#include "bu.h"
#include "dm.h"
#include "dm/util.h"
#include "rt/view.h"

#include "../include/private.h"
#include "../null/dm-Null.h"
}

#include "brlobol/init.h"
#include "brlobol/scene_group.h"
#include "brlobol/view_controller.h"

#include <Inventor/SbColor.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/SoDB.h>
#include <Inventor/SoViewport.h>
#include <Inventor/nodes/SoSeparator.h>

#include "./dm-obol.h"

struct obol_vars {
    struct dm *this_dm;
    void *view_ctx;
    SoDB::ContextManager *manager;
    BRLObolViewController *controller;
    int lighting;
    int zbuffer;
    int zclip;
    int transparency;
    int depthmask;
    double bound;
    int bound_flag;
    fastf_t viewscale;
};

static struct dm *obol_open(void *ctx, void *interp, int argc, const char **argv);
static int obol_close(struct dm *dmp);
static int obol_viable(const char *dpy_string);
static int obol_drawBegin(struct dm *dmp);
static int obol_drawEnd(struct dm *dmp);
static int obol_noop(struct dm *dmp);
static int obol_getDisplayImage(struct dm *dmp, unsigned char **image, int flip, int alpha);
static int obol_configureWin(struct dm *dmp, int force);
static int obol_reshape(struct dm *dmp, int width, int height);
static int obol_setFGColor(struct dm *dmp, unsigned char r, unsigned char g, unsigned char b, int strict, fastf_t transparency);
static int obol_setBGColor(struct dm *dmp, unsigned char r1, unsigned char g1, unsigned char b1, unsigned char r2, unsigned char g2, unsigned char b2);
static int obol_setLineAttr(struct dm *dmp, int width, int style);
static int obol_setWinBounds(struct dm *dmp, fastf_t *w);
static int obol_setLight(struct dm *dmp, int light_on);
static int obol_getLight(struct dm *dmp);
static int obol_setTransparency(struct dm *dmp, int transparency_on);
static int obol_getTransparency(struct dm *dmp);
static int obol_setDepthMask(struct dm *dmp, int depthMask_on);
static int obol_setZBuffer(struct dm *dmp, int zbuffer_on);
static int obol_getZBuffer(struct dm *dmp);
static int obol_setZClip(struct dm *dmp, int zclip);
static int obol_getZClip(struct dm *dmp);
static int obol_setBound(struct dm *dmp, double bound);
static double obol_getBound(struct dm *dmp);
static int obol_setBoundFlag(struct dm *dmp, int bound);
static int obol_getBoundFlag(struct dm *dmp);
static int obol_debug(struct dm *dmp, int lvl);
static int obol_logfile(struct dm *dmp, const char *filename);
static int obol_internal_var(struct bu_vls *result, struct dm *dmp, const char *key);
static int obol_makeCurrent(struct dm *dmp);
static int obol_openFb(struct dm *dmp);
static int obol_write_image(struct bu_vls *msgs, FILE *fp, struct dm *dmp);
static int obol_event_cmp(struct dm *dmp, dm_event_t type, int event);

#define OBOL_MV_O(_m) offsetof(struct obol_vars, _m)

static void
obol_lighting_hook(const struct bu_structparse *sdp, const char *name,
	void *base, const char *value, void *data)
{
    struct obol_vars *pv = (struct obol_vars *)base;
    if (pv && pv->this_dm)
	(void)dm_set_light(pv->this_dm, pv->lighting);
    dm_generic_hook(sdp, name, base, value, data);
}

static void
obol_zbuffer_hook(const struct bu_structparse *sdp, const char *name,
	void *base, const char *value, void *data)
{
    struct obol_vars *pv = (struct obol_vars *)base;
    if (pv && pv->this_dm)
	(void)dm_set_zbuffer(pv->this_dm, pv->zbuffer);
    dm_generic_hook(sdp, name, base, value, data);
}

static void
obol_zclip_hook(const struct bu_structparse *sdp, const char *name,
	void *base, const char *value, void *data)
{
    struct obol_vars *pv = (struct obol_vars *)base;
    if (pv && pv->this_dm)
	dm_set_zclip(pv->this_dm, pv->zclip);
    dm_generic_hook(sdp, name, base, value, data);
}

static void
obol_transparency_hook(const struct bu_structparse *sdp, const char *name,
	void *base, const char *value, void *data)
{
    struct obol_vars *pv = (struct obol_vars *)base;
    if (pv && pv->this_dm)
	(void)dm_set_transparency(pv->this_dm, pv->transparency);
    dm_generic_hook(sdp, name, base, value, data);
}

static void
obol_bound_hook(const struct bu_structparse *sdp, const char *name,
	void *base, const char *value, void *data)
{
    struct obol_vars *pv = (struct obol_vars *)base;
    if (pv && pv->this_dm)
	dm_set_bound(pv->this_dm, pv->bound);
    dm_generic_hook(sdp, name, base, value, data);
}

static void
obol_bound_flag_hook(const struct bu_structparse *sdp, const char *name,
	void *base, const char *value, void *data)
{
    struct obol_vars *pv = (struct obol_vars *)base;
    if (pv && pv->this_dm)
	dm_set_bound_flag(pv->this_dm, pv->bound_flag);
    dm_generic_hook(sdp, name, base, value, data);
}

static struct bu_structparse obol_vparse[] = {
    {"%d", 1, "zclip",        OBOL_MV_O(zclip),       obol_zclip_hook, NULL, NULL},
    {"%d", 1, "zbuffer",      OBOL_MV_O(zbuffer),     obol_zbuffer_hook, NULL, NULL},
    {"%d", 1, "lighting",     OBOL_MV_O(lighting),    obol_lighting_hook, NULL, NULL},
    {"%d", 1, "transparency", OBOL_MV_O(transparency), obol_transparency_hook, NULL, NULL},
    {"%g", 1, "bound",        OBOL_MV_O(bound),       obol_bound_hook, NULL, NULL},
    {"%d", 1, "useBound",     OBOL_MV_O(bound_flag),  obol_bound_flag_hook, NULL, NULL},
    {"",   0, (char *)0,      0,                      BU_STRUCTPARSE_FUNC_NULL, NULL, NULL}
};

static SoDB::ContextManager *
obol_context_manager(void)
{
    static SoDB::ContextManager *manager = SoDB::createOSMesaContextManager();
    return manager;
}

static struct obol_vars *
obol_pvars(struct dm *dmp)
{
    return (dmp && dmp->i) ?
	(struct obol_vars *)dmp->i->dm_vars.priv_vars : NULL;
}

static int
obol_dimensions(struct dm *dmp, int *width, int *height)
{
    if (!dmp)
	return BRLCAD_ERROR;

    struct obol_vars *pv = obol_pvars(dmp);
    int vw = pv && pv->view_ctx ? dm_view_context_width_get(pv->view_ctx) : 0;
    int vh = pv && pv->view_ctx ? dm_view_context_height_get(pv->view_ctx) : 0;
    int w = dmp->i->dm_width > 0 ? dmp->i->dm_width : (vw > 0 ? vw : 512);
    int h = dmp->i->dm_height > 0 ? dmp->i->dm_height : (vh > 0 ? vh : 512);

    if (w <= 0 || h <= 0)
	return BRLCAD_ERROR;

    dmp->i->dm_width = w;
    dmp->i->dm_height = h;
    dmp->i->dm_aspect = (fastf_t)w / (fastf_t)h;
    if (pv && pv->view_ctx)
	dm_view_context_dimensions_set(pv->view_ctx, w, h);

    if (width)
	*width = w;
    if (height)
	*height = h;
    return BRLCAD_OK;
}

static int
obol_sync_view(struct dm *dmp)
{
    struct obol_vars *pv = obol_pvars(dmp);
    if (!pv || !pv->controller)
	return BRLCAD_ERROR;

    int width = 0;
    int height = 0;
    if (obol_dimensions(dmp, &width, &height) != BRLCAD_OK)
	return BRLCAD_ERROR;

    pv->controller->setViewportSize((unsigned int)width,
	(unsigned int)height);
    if (pv->view_ctx &&
	    !pv->controller->syncCameraFromViewContext(pv->view_ctx))
	return BRLCAD_ERROR;

    return BRLCAD_OK;
}

static int
obol_render_to_image(struct dm *dmp, unsigned char **image, int flip,
	int alpha)
{
    struct obol_vars *pv = obol_pvars(dmp);
    if (!pv || !pv->manager || !pv->controller ||
	    !pv->controller->getViewport() || !image)
	return BRLCAD_ERROR;
    *image = NULL;

    if (obol_sync_view(dmp) != BRLCAD_OK)
	return BRLCAD_ERROR;

    SbColor background(
	(float)dmp->i->dm_bg1[0] / 255.0f,
	(float)dmp->i->dm_bg1[1] / 255.0f,
	(float)dmp->i->dm_bg1[2] / 255.0f);
    BRLObolProgressiveStatus progressive_status;
    int ret = pv->controller->renderToImage(image, flip, alpha, &background,
	    pv->manager, &progressive_status);
    if (progressive_status.hasMore)
	dm_set_native_repaint_pending(dmp, 1);

    return ret;
}

static int
obol_render(struct dm *dmp)
{
    unsigned char *image = NULL;
    int ret = obol_render_to_image(dmp, &image, 0, 0);
    if (image)
	bu_free(image, "dm-obol display image");
    return ret;
}

static int
obol_viable(const char *UNUSED(dpy_string))
{
    return obol_context_manager() ? 1 : 0;
}

static struct dm *
obol_open(void *ctx, void *interp, int argc, const char **argv)
{
    static int count = 0;
    struct dm *dmp = NULL;
    struct dm_impl *dmpi = NULL;
    struct obol_vars *privars = NULL;
    struct bu_vls init_proc_vls = BU_VLS_INIT_ZERO;

    SoDB::ContextManager *manager = obol_context_manager();
    if (!manager)
	return DM_NULL;
    brlobol_init(manager);

    BU_GET(dmp, struct dm);
    dmp->magic = DM_MAGIC;
    dmp->start_time = 0;

    BU_GET(dmpi, struct dm_impl);
    *dmpi = *dm_obol.i;
    dmp->i = dmpi;

    BU_ALLOC(dmp->i->dm_vars.priv_vars, struct obol_vars);
    privars = (struct obol_vars *)dmp->i->dm_vars.priv_vars;
    privars->this_dm = dmp;
    privars->view_ctx = ctx;
    privars->manager = manager;
    privars->controller = new BRLObolViewController(new SoBRLSceneGroup);
    privars->lighting = 1;
    privars->zbuffer = 1;
    privars->zclip = 0;
    privars->transparency = 0;
    privars->depthmask = 1;
    privars->bound = 1.0;
    privars->bound_flag = 1;
    privars->viewscale = 1.0;

    dmp->i->dm_vars.pub_vars = NULL;
    dmp->i->m_vars = privars;
    dmp->i->p_vars = privars->controller;
    dmp->i->dm_ctx = ctx;
    dmp->i->dm_interp = interp;
    dmp->i->dm_vp = &privars->viewscale;
    dmp->i->dm_lineWidth = 1;
    dmp->i->dm_lineStyle = 0;
    dmp->i->dm_depthMask = 1;
    dmp->i->dm_bytes_per_pixel = 3;
    dmp->i->dm_bits_per_channel = 8;
    dmp->i->dm_fontsize = 20;

    bu_vls_init(&dmp->i->dm_pathName);
    bu_vls_init(&dmp->i->dm_tkName);
    bu_vls_init(&dmp->i->dm_dName);
    bu_vls_init(&dmp->i->dm_log);

    dm_processOptions(dmp, &init_proc_vls, --argc, ++argv);
    if (bu_vls_strlen(&dmp->i->dm_pathName) == 0)
	bu_vls_printf(&dmp->i->dm_pathName, ".dm_obol%d", count);
    if (bu_vls_strlen(&dmp->i->dm_dName) == 0)
	bu_vls_strcpy(&dmp->i->dm_dName, "obol");
    ++count;
    bu_vls_free(&init_proc_vls);

    if (obol_dimensions(dmp, NULL, NULL) != BRLCAD_OK ||
	    obol_sync_view(dmp) != BRLCAD_OK) {
	obol_close(dmp);
	return DM_NULL;
    }
    (void)obol_setBGColor(dmp, 0, 0, 0, 0, 0, 0);
    (void)obol_setFGColor(dmp, 255, 255, 255, 1, 1.0);

    return dmp;
}

static int
obol_close(struct dm *dmp)
{
    if (!dmp)
	return BRLCAD_OK;

    struct obol_vars *pv = obol_pvars(dmp);
    if (pv) {
	delete pv->controller;
	bu_free(pv, "dm-obol private vars");
    }

    bu_vls_free(&dmp->i->dm_pathName);
    bu_vls_free(&dmp->i->dm_tkName);
    bu_vls_free(&dmp->i->dm_dName);
    bu_vls_free(&dmp->i->dm_log);
    BU_PUT(dmp->i, struct dm_impl);
    BU_PUT(dmp, struct dm);
    return BRLCAD_OK;
}

static int
obol_drawBegin(struct dm *UNUSED(dmp))
{
    return BRLCAD_OK;
}

static int
obol_drawEnd(struct dm *dmp)
{
    return obol_render(dmp);
}

static int
obol_noop(struct dm *UNUSED(dmp))
{
    return BRLCAD_OK;
}

static int
obol_getDisplayImage(struct dm *dmp, unsigned char **image, int flip, int alpha)
{
    if (!image)
	return BRLCAD_ERROR;
    return obol_render_to_image(dmp, image, flip, alpha);
}

static int
obol_configureWin(struct dm *dmp, int UNUSED(force))
{
    return obol_sync_view(dmp);
}

static int
obol_reshape(struct dm *dmp, int width, int height)
{
    if (!dmp || width <= 0 || height <= 0)
	return BRLCAD_ERROR;
    dmp->i->dm_width = width;
    dmp->i->dm_height = height;
    return obol_sync_view(dmp);
}

static int
obol_setFGColor(struct dm *dmp, unsigned char r, unsigned char g,
	unsigned char b, int UNUSED(strict), fastf_t UNUSED(transparency))
{
    if (!dmp)
	return BRLCAD_ERROR;
    dmp->i->dm_fg[0] = r;
    dmp->i->dm_fg[1] = g;
    dmp->i->dm_fg[2] = b;
    return BRLCAD_OK;
}

static int
obol_setBGColor(struct dm *dmp, unsigned char r1, unsigned char g1,
	unsigned char b1, unsigned char r2, unsigned char g2,
	unsigned char b2)
{
    if (!dmp)
	return BRLCAD_ERROR;
    dmp->i->dm_bg1[0] = r1;
    dmp->i->dm_bg1[1] = g1;
    dmp->i->dm_bg1[2] = b1;
    dmp->i->dm_bg2[0] = r2;
    dmp->i->dm_bg2[1] = g2;
    dmp->i->dm_bg2[2] = b2;
    if (dmp->i->dm_ctx) {
	struct bv_background_state background = BV_BACKGROUND_STATE_INIT;
	VSET(background.bottom, r1, g1, b1);
	VSET(background.top, r2, g2, b2);
	(void)bv_background_state_set(
	    bv_context_view(static_cast<struct bv_context *>(dmp->i->dm_ctx)),
	    &background);
    }
    struct obol_vars *pv = obol_pvars(dmp);
    if (pv && pv->controller)
	pv->controller->setBackgroundColors(
	    SbColor(r1 / 255.0f, g1 / 255.0f, b1 / 255.0f),
	    SbColor(r2 / 255.0f, g2 / 255.0f, b2 / 255.0f));
    return BRLCAD_OK;
}

static int
obol_setLineAttr(struct dm *dmp, int width, int style)
{
    if (!dmp)
	return BRLCAD_ERROR;
    dmp->i->dm_lineWidth = width;
    dmp->i->dm_lineStyle = style;
    return BRLCAD_OK;
}

static int
obol_setWinBounds(struct dm *dmp, fastf_t *w)
{
    if (!dmp || !w)
	return BRLCAD_ERROR;
    dmp->i->dm_clipmin[0] = w[0];
    dmp->i->dm_clipmin[1] = w[2];
    dmp->i->dm_clipmin[2] = w[4];
    dmp->i->dm_clipmax[0] = w[1];
    dmp->i->dm_clipmax[1] = w[3];
    dmp->i->dm_clipmax[2] = w[5];
    return BRLCAD_OK;
}

static int
obol_setLight(struct dm *dmp, int light_on)
{
    struct obol_vars *pv = obol_pvars(dmp);
    if (!pv)
	return BRLCAD_ERROR;
    pv->lighting = light_on ? 1 : 0;
    return BRLCAD_OK;
}

static int
obol_getLight(struct dm *dmp)
{
    struct obol_vars *pv = obol_pvars(dmp);
    return pv ? pv->lighting : 0;
}

static int
obol_setTransparency(struct dm *dmp, int transparency_on)
{
    struct obol_vars *pv = obol_pvars(dmp);
    if (!pv)
	return BRLCAD_ERROR;
    pv->transparency = transparency_on ? 1 : 0;
    return BRLCAD_OK;
}

static int
obol_getTransparency(struct dm *dmp)
{
    struct obol_vars *pv = obol_pvars(dmp);
    return pv ? pv->transparency : 0;
}

static int
obol_setDepthMask(struct dm *dmp, int depthMask_on)
{
    struct obol_vars *pv = obol_pvars(dmp);
    if (!pv || !dmp)
	return BRLCAD_ERROR;
    pv->depthmask = depthMask_on ? 1 : 0;
    dmp->i->dm_depthMask = pv->depthmask;
    return BRLCAD_OK;
}

static int
obol_setZBuffer(struct dm *dmp, int zbuffer_on)
{
    struct obol_vars *pv = obol_pvars(dmp);
    if (!pv)
	return BRLCAD_ERROR;
    pv->zbuffer = zbuffer_on ? 1 : 0;
    return BRLCAD_OK;
}

static int
obol_getZBuffer(struct dm *dmp)
{
    struct obol_vars *pv = obol_pvars(dmp);
    return pv ? pv->zbuffer : 0;
}

static int
obol_setZClip(struct dm *dmp, int zclip)
{
    struct obol_vars *pv = obol_pvars(dmp);
    if (!pv)
	return BRLCAD_ERROR;
    pv->zclip = zclip ? 1 : 0;
    return BRLCAD_OK;
}

static int
obol_getZClip(struct dm *dmp)
{
    struct obol_vars *pv = obol_pvars(dmp);
    return pv ? pv->zclip : 0;
}

static int
obol_setBound(struct dm *dmp, double bound)
{
    struct obol_vars *pv = obol_pvars(dmp);
    if (!pv)
	return BRLCAD_ERROR;
    pv->bound = bound;
    return BRLCAD_OK;
}

static double
obol_getBound(struct dm *dmp)
{
    struct obol_vars *pv = obol_pvars(dmp);
    return pv ? pv->bound : 0.0;
}

static int
obol_setBoundFlag(struct dm *dmp, int bound)
{
    struct obol_vars *pv = obol_pvars(dmp);
    if (!pv)
	return BRLCAD_ERROR;
    pv->bound_flag = bound ? 1 : 0;
    return BRLCAD_OK;
}

static int
obol_getBoundFlag(struct dm *dmp)
{
    struct obol_vars *pv = obol_pvars(dmp);
    return pv ? pv->bound_flag : 0;
}

static int
obol_debug(struct dm *dmp, int lvl)
{
    if (!dmp)
	return BRLCAD_ERROR;
    dmp->i->dm_debugLevel = lvl;
    return BRLCAD_OK;
}

static int
obol_logfile(struct dm *dmp, const char *filename)
{
    if (!dmp || !filename)
	return BRLCAD_ERROR;
    bu_vls_sprintf(&dmp->i->dm_log, "%s", filename);
    return BRLCAD_OK;
}

static int
obol_internal_var(struct bu_vls *result, struct dm *UNUSED(dmp),
	const char *UNUSED(key))
{
    if (result)
	bu_vls_printf(result, "dm-obol has no backend-specific variables");
    return BRLCAD_OK;
}

static int
obol_makeCurrent(struct dm *UNUSED(dmp))
{
    return BRLCAD_OK;
}

static int
obol_openFb(struct dm *UNUSED(dmp))
{
    return BRLCAD_ERROR;
}

static int
obol_write_image(struct bu_vls *UNUSED(msgs), FILE *UNUSED(fp),
	struct dm *UNUSED(dmp))
{
    return BRLCAD_ERROR;
}

static int
obol_event_cmp(struct dm *UNUSED(dmp), dm_event_t UNUSED(type),
	int UNUSED(event))
{
    return -1;
}

struct dm_impl dm_obol_impl = {
    obol_open,
    obol_close,
    obol_viable,
    obol_drawBegin,
    obol_drawEnd,
    obol_noop,
    obol_noop,
    null_loadMatrix,
    null_loadPMatrix,
    null_popPMatrix,
    null_drawString2D,
    NULL,
    null_String2DBBox,
    null_drawLine2D,
    null_drawLine3D,
    null_drawLines3D,
    null_drawPoint2D,
    null_drawPoint3D,
    null_drawPoints3D,
    null_drawVList,
    null_drawVListHiddenLine,
    null_draw,
    obol_setFGColor,
    obol_setBGColor,
    obol_setLineAttr,
    obol_configureWin,
    obol_setWinBounds,
    obol_setLight,
    obol_getLight,
    obol_setTransparency,
    obol_getTransparency,
    obol_setDepthMask,
    obol_setZBuffer,
    obol_getZBuffer,
    obol_setZClip,
    obol_getZClip,
    obol_setBound,
    obol_getBound,
    obol_setBoundFlag,
    obol_getBoundFlag,
    obol_debug,
    obol_logfile,
    obol_getDisplayImage,
    obol_reshape,
    obol_makeCurrent,
    null_SwapBuffers,
    null_doevent,
    obol_openFb,
    NULL,
    NULL,
    NULL,
    obol_internal_var,
    obol_write_image,
    NULL,
    NULL,
    obol_event_cmp,
    NULL,
    0,
    1,				/* graphical, but headless */
    "Obol",
    1,				/* retained scene cache */
    0,				/* no stereo */
    "obol",
    "Obol/Coin headless graphics",
    1,				/* top */
    0,				/* width */
    0,				/* height */
    0,				/* dirty */
    3,				/* bytes per pixel */
    8,				/* bits per channel */
    1,				/* line width */
    0,				/* line style */
    1.0,			/* aspect ratio */
    0,
    {0, 0},
    NULL,
    NULL,
    BU_VLS_INIT_ZERO,		/* path name */
    BU_VLS_INIT_ZERO,		/* tk name */
    BU_VLS_INIT_ZERO,		/* display name */
    BU_VLS_INIT_ZERO,		/* logfile */
    {0, 0, 0},			/* bg1 color */
    {0, 0, 0},			/* bg2 color */
    {255, 255, 255},		/* fg color */
    {255, 0, 0},		/* geometry default color */
    {RT_VIEW_MIN, RT_VIEW_MIN, RT_VIEW_MIN},
    {RT_VIEW_MAX, RT_VIEW_MAX, RT_VIEW_MAX},
    0,				/* debugging */
    0,				/* perspective */
    1,				/* depth buffer is writable */
    1,				/* clear after drawing */
    0,				/* font override */
    obol_vparse,
    FB_NULL,
    0,				/* interpreter */
    NULL,			/* drawing context */
    NULL			/* app data */
};

struct dm dm_obol = { DM_MAGIC, &dm_obol_impl, 0 };

#ifdef DM_PLUGIN
static const struct dm_plugin pinfo = { DM_API, &dm_obol };
extern "C" {
COMPILER_DLLEXPORT const struct dm_plugin *dm_plugin_info(void)
{
    return &pinfo;
}
}
#endif

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
