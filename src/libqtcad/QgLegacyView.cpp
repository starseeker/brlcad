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
#include "ged.h"
#include "ged/draw.h"
#include "ged/view.h"
#include "QgLegacyViewContext.h"
#include "qtcad/QgLegacyView.h"

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

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
