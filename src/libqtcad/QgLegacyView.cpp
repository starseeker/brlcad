/*                    Q G L E G A C Y V I E W . C P P
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
/** @file QgLegacyView.cpp
 *
 * Neutral qtcad helpers for the staged legacy view handle.
 */

#include "common.h"

#include "bu/malloc.h"
#include "dm.h"
#include "ged.h"
#include "ged/draw.h"
#include "QgLegacyViewContext.h"
#include "QgLegacyViewDm.h"
#include "qtcad/QgLegacyView.h"

static struct dm *
qg_legacy_dm_to_dm(qg_legacy_dm *dmp)
{
    return reinterpret_cast<struct dm *>(dmp);
}

static qg_legacy_dm *
qg_legacy_dm_from_dm(struct dm *dmp)
{
    return reinterpret_cast<qg_legacy_dm *>(dmp);
}

static struct fb *
qg_legacy_fb_to_fb(qg_legacy_fb *ifp)
{
    return reinterpret_cast<struct fb *>(ifp);
}

static qg_legacy_fb *
qg_legacy_fb_from_fb(struct fb *ifp)
{
    return reinterpret_cast<qg_legacy_fb *>(ifp);
}

qg_legacy_view *
qg_legacy_view_local_create(const char *name)
{
    void *view_ctx = rt_view_context_create();
    if (name)
	rt_view_context_name_set(view_ctx, name);
    return qg_legacy_view_from_context(view_ctx);
}

void
qg_legacy_view_local_free(qg_legacy_view *view)
{
    /* Canvas-local views may outlive GED scene cleanup in tests and app
     * shutdown; keep the historical qtcad release behavior until the lower
     * retained view source is retired. */
    rt_view_context_release_storage(qg_legacy_view_to_context(view));
}

void
qg_legacy_view_local_destroy(qg_legacy_view *view)
{
    void *view_ctx = qg_legacy_view_to_context(view);
    if (!view_ctx)
	return;

    rt_view_context_free(view_ctx);
}

int
qg_legacy_view_dimensions_set(qg_legacy_view *view, int width, int height)
{
    return rt_view_context_dimensions_set(qg_legacy_view_to_context(view),
	    width, height);
}

int
qg_legacy_view_unit_conversion_set(qg_legacy_view *view,
	double local2base,
	double base2local)
{
    return rt_view_context_unit_conversion_set(qg_legacy_view_to_context(view),
	    local2base, base2local);
}

int
qg_legacy_view_name_set(qg_legacy_view *view, const char *name)
{
    return rt_view_context_name_set(qg_legacy_view_to_context(view), name);
}

qg_legacy_view *
qg_legacy_view_ged_active_get(struct ged *gedp)
{
    return qg_legacy_view_from_context(ged_view_active_ctx(gedp));
}

void
qg_legacy_view_ged_active_set(struct ged *gedp, qg_legacy_view *view)
{
    ged_view_active_ctx_set(gedp, qg_legacy_view_to_context(view));
}

struct db_i *
qg_legacy_view_ged_database(struct ged *gedp)
{
    return gedp ? gedp->dbip : nullptr;
}

int
qg_legacy_view_ged_view_set_add(struct ged *gedp, qg_legacy_view *view)
{
    if (!gedp || !view)
	return 0;

    return rt_view_set_context_add(ged_view_set_ctx(gedp),
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_ged_view_set_remove(struct ged *gedp, qg_legacy_view *view)
{
    if (!gedp || !view)
	return 0;

    return rt_view_set_context_remove(ged_view_set_ctx(gedp),
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_ged_view_set_attach(struct ged *gedp, qg_legacy_view *view)
{
    if (!gedp || !view)
	return 0;

    return rt_view_context_view_set_attach(qg_legacy_view_to_context(view),
	    ged_view_set_ctx(gedp));
}

int
qg_legacy_view_display_manager_set(qg_legacy_view *view, qg_legacy_dm *dmp)
{
    return rt_view_context_display_manager_set(qg_legacy_view_to_context(view),
	    qg_legacy_dm_to_dm(dmp));
}

int
qg_legacy_view_dm_bind(qg_legacy_view *view, qg_legacy_dm *dmp)
{
    if (!view || !dmp)
	return 0;

    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    dm_set_vp(dm, rt_view_context_scale_storage_get(
		qg_legacy_view_to_context(view)));
    return qg_legacy_view_display_manager_set(view, dmp);
}

int
qg_legacy_view_dm_sync_dimensions(qg_legacy_view *view, qg_legacy_dm *dmp)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    if (!view || !dm)
	return 0;

    return qg_legacy_view_dimensions_set(view, dm_get_width(dm),
	    dm_get_height(dm));
}

unsigned long long
qg_legacy_view_dm_hash(qg_legacy_dm *dmp)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    return dm ? dm_hash(dm) : 0ULL;
}

int
qg_legacy_view_dm_native_repaint_pending_get(qg_legacy_dm *dmp)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    return dm ? dm_get_native_repaint_pending(dm) : 0;
}

void
qg_legacy_view_dm_native_repaint_pending_set(qg_legacy_dm *dmp, int pending)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    if (dm)
	dm_set_native_repaint_pending(dm, pending);
}

int
qg_legacy_view_dm_configure_window(qg_legacy_dm *dmp, int force)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    return dm ? dm_configure_win(dm, force) : 0;
}

int
qg_legacy_view_dm_width_get(qg_legacy_dm *dmp)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    return dm ? dm_get_width(dm) : 0;
}

int
qg_legacy_view_dm_height_get(qg_legacy_dm *dmp)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    return dm ? dm_get_height(dm) : 0;
}

void
qg_legacy_view_dm_dimensions_set(qg_legacy_dm *dmp, int width, int height)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    if (!dm)
	return;

    dm_set_width(dm, width);
    dm_set_height(dm, height);
}

int
qg_legacy_view_dm_framebuffer_setup_existing(qg_legacy_fb *ifp,
	qg_legacy_dm *dmp)
{
    struct fb *fb = qg_legacy_fb_to_fb(ifp);
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    if (!fb || !dm)
	return 0;

    struct fb_platform_specific *fbps = fb_get_platform_specific(FB_QTGL_MAGIC);
    if (!fbps)
	return 0;

    fbps->data = (void *)dm;
    fb_setup_existing(fb, dm_get_width(dm), dm_get_height(dm), fbps);
    fb_put_platform_specific(fbps);
    return 1;
}

qg_legacy_dm *
qg_legacy_view_dm_open_qtgl(void *context)
{
    const char *acmd = "attach";
    return qg_legacy_dm_from_dm(dm_open(context, nullptr, "qtgl", 1, &acmd));
}

qg_legacy_dm *
qg_legacy_view_dm_open_swrast(qg_legacy_view *view, void *context)
{
    const char *acmd = "attach";
    struct dm *dmp = dm_open(qg_legacy_view_to_context(view), nullptr,
	    "swrast", 1, &acmd);
    if (dmp)
	dm_set_udata(dmp, context);
    return qg_legacy_dm_from_dm(dmp);
}

int
qg_legacy_view_dm_close(qg_legacy_dm *dmp)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    return dm ? dm_close(dm) : 0;
}

int
qg_legacy_view_dm_setup_qtgl(qg_legacy_dm *dmp, fastf_t zmin, fastf_t zmax)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    if (!dm)
	return 0;

    dm_set_pathname(dm, "QTDM");
    fastf_t windowbounds[6] = { -1, 1, -1, 1, zmin, zmax };
    dm_set_win_bounds(dm, windowbounds);
    return 1;
}

int
qg_legacy_view_dm_setup_swrast(qg_legacy_dm *dmp, fastf_t zmin, fastf_t zmax)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    if (!dm)
	return 0;

    dm_configure_win(dm, 0);
    dm_set_pathname(dm, "SWDM");
    dm_set_zbuffer(dm, 1);
    fastf_t windowbounds[6] = { -1, 1, -1, 1, zmin, zmax };
    dm_set_win_bounds(dm, windowbounds);
    return 1;
}

int
qg_legacy_view_dm_background_get(qg_legacy_dm *dmp, unsigned char bg1[3],
	unsigned char bg2[3])
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    if (!dm || !bg1 || !bg2)
	return 0;

    unsigned char *dm_bg1 = nullptr;
    unsigned char *dm_bg2 = nullptr;
    dm_get_bg(&dm_bg1, &dm_bg2, dm);
    if (!dm_bg1 || !dm_bg2)
	return 0;

    bg1[0] = dm_bg1[0];
    bg1[1] = dm_bg1[1];
    bg1[2] = dm_bg1[2];
    bg2[0] = dm_bg2[0];
    bg2[1] = dm_bg2[1];
    bg2[2] = dm_bg2[2];
    return 1;
}

int
qg_legacy_view_dm_background_set(qg_legacy_dm *dmp, const unsigned char bg1[3],
	const unsigned char bg2[3])
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    if (!dm || !bg1 || !bg2)
	return 0;

    dm_set_bg(dm, bg1[0], bg1[1], bg1[2], bg2[0], bg2[1], bg2[2]);
    return 1;
}

int
qg_legacy_view_dm_background_restore(qg_legacy_dm *dmp)
{
    unsigned char bg1[3] = {0, 0, 0};
    unsigned char bg2[3] = {0, 0, 0};
    if (!qg_legacy_view_dm_background_get(dmp, bg1, bg2))
	return 0;

    return qg_legacy_view_dm_background_set(dmp, bg1, bg2);
}

int
qg_legacy_view_dm_draw_begin(qg_legacy_dm *dmp)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    return dm ? dm_draw_begin(dm) : 0;
}

int
qg_legacy_view_dm_draw_end(qg_legacy_dm *dmp)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    return dm ? dm_draw_end(dm) : 0;
}

int
qg_legacy_view_dm_display_image_get(qg_legacy_dm *dmp, unsigned char **image,
	int copy, int release)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    return (dm && image) ? dm_get_display_image(dm, image, copy, release) : 1;
}

void
qg_legacy_view_dm_display_image_release(unsigned char *image)
{
    if (image)
	bu_free(image, "copy of backend image");
}

int
qg_legacy_view_dm_load_current_model2view(qg_legacy_dm *dmp,
	const qg_legacy_view *view, int which_eye_or_mode)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    if (!dm || !view)
	return 0;

    mat_t model2view;
    if (!rt_view_context_model2view_get(model2view,
	    qg_legacy_view_to_context(view)))
	return 0;

    dm_loadmatrix(dm, model2view, which_eye_or_mode);
    return 1;
}

qg_legacy_fb *
qg_legacy_view_framebuffer_raw_create(const char *type)
{
    struct fb *ifp = fb_raw(type);
    if (ifp)
	fb_set_standalone(ifp, 0);
    return qg_legacy_fb_from_fb(ifp);
}

qg_legacy_fb *
qg_legacy_view_framebuffer_handle_from_raw(void *ifp)
{
    return reinterpret_cast<qg_legacy_fb *>(ifp);
}

QgGL *
qg_legacy_view_framebuffer_qtgl_canvas_get(qg_legacy_fb *ifp)
{
    struct fb *fb = qg_legacy_fb_to_fb(ifp);
    struct dm *dm = fb ? fb_get_dm(fb) : nullptr;
    return dm ? reinterpret_cast<QgGL *>(dm_get_ctx(dm)) : nullptr;
}

QgSW *
qg_legacy_view_framebuffer_swrast_canvas_get(qg_legacy_fb *ifp)
{
    struct fb *fb = qg_legacy_fb_to_fb(ifp);
    struct dm *dm = fb ? fb_get_dm(fb) : nullptr;
    return dm ? reinterpret_cast<QgSW *>(dm_get_udata(dm)) : nullptr;
}

int
qg_legacy_view_framebuffer_release(qg_legacy_fb *ifp, int initialized)
{
    struct fb *fb = qg_legacy_fb_to_fb(ifp);
    if (!fb || fb_get_standalone(fb))
	return 0;

    if (initialized)
	return fb_close_existing(fb);

    fb_put(fb);
    return 1;
}

int
qg_legacy_view_framebuffer_standalone_get(qg_legacy_fb *ifp)
{
    struct fb *fb = qg_legacy_fb_to_fb(ifp);
    return fb ? fb_get_standalone(fb) : 0;
}

int
qg_legacy_view_framebuffer_configure(qg_legacy_fb *ifp, int width, int height)
{
    struct fb *fb = qg_legacy_fb_to_fb(ifp);
    return fb ? fb_configure_window(fb, width, height) : 0;
}

int
qg_legacy_view_dm_draw(qg_legacy_view *view)
{
    void *view_ctx = qg_legacy_view_to_context(view);
    if (!view_ctx)
	return 0;

    dm_draw_objs(view_ctx);
    return 1;
}

const char *
qg_legacy_view_dm_init_messages(void)
{
    return dm_init_msgs();
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
