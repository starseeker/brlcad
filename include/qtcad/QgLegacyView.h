/*                    Q G L E G A C Y V I E W . H
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
/** @file QgLegacyView.h
 *
 * Opaque transitional qtcad view handle.
 */

#ifndef QGLEGACYVIEW_H
#define QGLEGACYVIEW_H

#include "qtcad/defines.h"

struct db_i;
struct ged;
struct qg_legacy_view;

QTCAD_EXPORT extern qg_legacy_view *qg_legacy_view_local_create(
	const char *name);

QTCAD_EXPORT extern void qg_legacy_view_local_free(qg_legacy_view *view);

QTCAD_EXPORT extern void qg_legacy_view_local_destroy(qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_dimensions_set(qg_legacy_view *view,
	int width,
	int height);

QTCAD_EXPORT extern int qg_legacy_view_unit_conversion_set(qg_legacy_view *view,
	double local2base,
	double base2local);

QTCAD_EXPORT extern int qg_legacy_view_name_set(qg_legacy_view *view,
	const char *name);

QTCAD_EXPORT extern qg_legacy_view *qg_legacy_view_ged_active_get(
	struct ged *gedp);

QTCAD_EXPORT extern void qg_legacy_view_ged_active_set(struct ged *gedp,
	qg_legacy_view *view);

QTCAD_EXPORT extern struct db_i *qg_legacy_view_ged_database(
	struct ged *gedp);

QTCAD_EXPORT extern int qg_legacy_view_ged_view_set_add(struct ged *gedp,
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_ged_view_set_remove(struct ged *gedp,
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_ged_view_set_attach(struct ged *gedp,
	qg_legacy_view *view);

#endif /* QGLEGACYVIEW_H */

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
