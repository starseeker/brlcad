/*              Q G L E G A C Y V I E W S K E T C H . H
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
/** @file QgLegacyViewSketch.h
 *
 * Private qtcad sketch helpers for callers that need qsketch's temporary
 * retained edit/wireframe path without directly depending on the backend.
 */

#ifndef QGLEGACYVIEWSKETCH_H
#define QGLEGACYVIEWSKETCH_H

#include <stddef.h>

#include "qtcad/QgLegacyView.h"

struct bn_tol;
struct db_full_path;
struct db_i;
struct rt_edit;
struct qg_legacy_view_sketch_lines;

QTCAD_EXPORT extern struct rt_edit *qg_legacy_view_sketch_edit_create(
	qg_legacy_view *view,
	struct db_full_path *dfp,
	struct db_i *dbip,
	struct bn_tol *tol);

QTCAD_EXPORT extern struct qg_legacy_view_sketch_lines *
qg_legacy_view_sketch_lines_create(qg_legacy_view *view,
	const char *name,
	unsigned char r,
	unsigned char g,
	unsigned char b);

QTCAD_EXPORT extern int qg_legacy_view_sketch_lines_is_null(
	const struct qg_legacy_view_sketch_lines *lines);

QTCAD_EXPORT extern int qg_legacy_view_sketch_lines_set_points(
	struct qg_legacy_view_sketch_lines *lines,
	const point_t *points,
	const int *commands,
	size_t point_count);

QTCAD_EXPORT extern void qg_legacy_view_sketch_lines_destroy(
	struct qg_legacy_view_sketch_lines *lines);

#endif /* QGLEGACYVIEWSKETCH_H */

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
