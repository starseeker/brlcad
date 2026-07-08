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

#include "bv.h"
#include "bu/malloc.h"
#include "dm.h"
#include "ged.h"
#include "ged/draw.h"
#include "ged/view.h"
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

qg_legacy_view *
qg_legacy_view_local_create(const char *name)
{
    void *view_ctx = ged_view_context_create();
    if (name)
	bv_context_name_set(reinterpret_cast<struct bv_context *>(view_ctx),
		name);
    return qg_legacy_view_from_context(view_ctx);
}

void
qg_legacy_view_local_free(qg_legacy_view *view)
{
    qg_legacy_view_local_destroy(view);
}

void
qg_legacy_view_local_destroy(qg_legacy_view *view)
{
    void *view_ctx = qg_legacy_view_to_context(view);
    if (!view_ctx)
	return;

    ged_view_context_free(view_ctx);
}

int
qg_legacy_view_dimensions_set(qg_legacy_view *view, int width, int height)
{
    return bv_context_dimensions_set(
	reinterpret_cast<struct bv_context *>(qg_legacy_view_to_context(view)),
	width, height);
}

int
qg_legacy_view_unit_conversion_set(qg_legacy_view *view,
	double local2base,
	double base2local)
{
    return bv_unit_conversion_set(bv_context_view(
	reinterpret_cast<struct bv_context *>(qg_legacy_view_to_context(view))),
	local2base, base2local);
}

int
qg_legacy_view_name_set(qg_legacy_view *view, const char *name)
{
    return bv_context_name_set(
	reinterpret_cast<struct bv_context *>(qg_legacy_view_to_context(view)),
	name);
}

qg_legacy_view *
qg_legacy_view_ged_active_get(struct ged *gedp)
{
    return qg_legacy_view_from_context(ged_view_active_ctx(gedp));
}

void
qg_legacy_view_ged_active_set(struct ged *gedp, qg_legacy_view *view)
{
    void *view_ctx = qg_legacy_view_to_context(view);
    ged_view_active_ctx_set(gedp, view_ctx);
    (void)ged_view_context_host_attach(gedp, view_ctx);
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

    void *view_ctx = qg_legacy_view_to_context(view);
    if (!ged_view_set_context_add(ged_view_set_ctx(gedp), view_ctx))
	return 0;
    return ged_view_context_host_attach(gedp, view_ctx);
}

int
qg_legacy_view_ged_view_set_remove(struct ged *gedp, qg_legacy_view *view)
{
    if (!gedp || !view)
	return 0;

    return ged_view_set_context_remove(ged_view_set_ctx(gedp),
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_ged_view_set_attach(struct ged *gedp, qg_legacy_view *view)
{
    if (!gedp || !view)
	return 0;

    void *view_ctx = qg_legacy_view_to_context(view);
    if (!ged_view_context_view_set_attach(view_ctx, ged_view_set_ctx(gedp)))
	return 0;
    return ged_view_context_host_attach(gedp, view_ctx);
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
