/*                  Q G V I E W F I L T E R . C P P
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
/** @file QgViewFilter.cpp
 *
 * Shared Qt event-filter helpers for QgView-based interaction filters.
 */

#include "common.h"

#include "qtcad/QgObolSnap.h"
#include "qtcad/QgLegacyView.h"
#include "qtcad/QgView.h"
#include "qtcad/QgViewFilter.h"

static uint32_t
qg_obol_snap_kinds_from_legacy_mask(unsigned kinds)
{
	uint32_t obol_kinds = 0;

	if (kinds & QG_LEGACY_VIEW_SNAP_KIND_ENDPOINT)
		obol_kinds |= QgObolSnapRecord::ENDPOINT |
			QgObolSnapRecord::MIDPOINT |
			QgObolSnapRecord::LINE_NEAREST;
	if (kinds & QG_LEGACY_VIEW_SNAP_KIND_MIDPOINT)
		obol_kinds |= QgObolSnapRecord::MIDPOINT;
	if (kinds & QG_LEGACY_VIEW_SNAP_KIND_OVERLAY_HANDLE)
		obol_kinds |= QgObolSnapRecord::CENTER;

	return obol_kinds;
}

static int
qg_obol_db_snap_enabled(const qg_legacy_view *v)
{
	return qg_legacy_view_db_snap_enabled(v);
}

static float
qg_obol_snap_tolerance(const qg_legacy_view *v)
{
	if (!v)
		return 0.0f;

	int width = qg_legacy_view_width_get(v);
	int height = qg_legacy_view_height_get(v);
	if (width <= 0 || height <= 0)
		return 0.0f;

	double lavg = ((double)width + (double)height) * 0.5;
	double lratio = 1.0 / lavg;
	double tol = qg_legacy_view_size_get(v) * lratio *
		qg_legacy_view_snap_tolerance_factor_get(v);

	return tol > 0.0 ? (float)tol : 0.0f;
}

static void
qg_obol_refine_db_snap(QgView *display, qg_legacy_view *v)
{
	if (!display || !qg_obol_db_snap_enabled(v))
		return;

	uint32_t obol_kinds = qg_obol_snap_kinds_from_legacy_mask(
		qg_legacy_view_snap_kind_mask_get(v));
	float tolerance = qg_obol_snap_tolerance(v);
	if (!obol_kinds || tolerance <= 0.0f)
		return;

	QgObolSnapRecord record;
	point_t query_point = VINIT_ZERO;
	if (!qg_legacy_view_current_point_get(v, query_point))
		return;

	SbVec3f query((float)query_point[X],
		(float)query_point[Y],
		(float)query_point[Z]);
	if (!qg_obol_snap_point(display, query, tolerance, obol_kinds, record))
		return;

	point_t snapped_point = VINIT_ZERO;
	VSET(snapped_point, (fastf_t)record.point[0],
		(fastf_t)record.point[1], (fastf_t)record.point[2]);
	(void)qg_legacy_view_current_point_set(v, snapped_point);
}

class QgViewFilter::QgViewFilterPrivate {
public:
	QgView *display = nullptr;
};

QgViewFilter::QgViewFilter(QObject *parent)
	: QObject(parent), m(new QgViewFilterPrivate)
{
}

QgViewFilter::~QgViewFilter()
{
	delete m;
}

void
QgViewFilter::set_view_widget(QgView *display)
{
	m->display = display;
}

QgView *
QgViewFilter::view_widget() const
{
	return m->display;
}

QMouseEvent *
QgViewFilter::view_sync(QEvent *e)
{
	qg_legacy_view *lv = m->display ? m->display->view() : nullptr;
	if (!lv)
		return nullptr;

	/* If this is not a supported mouse event, there is nothing to do. */
	QMouseEvent *m_e = nullptr;
	if (e->type() == QEvent::MouseButtonPress || e->type() == QEvent::MouseButtonRelease || e->type() == QEvent::MouseButtonDblClick || e->type() == QEvent::MouseMove)
		m_e = (QMouseEvent *)e;
	if (!m_e)
		return nullptr;

	/* Capture current mouse coordinates from the Qt event. */
	int e_x = 0, e_y = 0;
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
	e_x = m_e->x();
	e_y = m_e->y();
#else
	e_x = (int)m_e->position().x();
	e_y = (int)m_e->position().y();
#endif

	/* Keep retained legacy view state synchronized with the event stream. */
	qg_legacy_view_mouse_state_set(lv, e_x, e_y);
	qg_obol_refine_db_snap(m->display, lv);

	/* Modifier keys are typically view-nav gestures, not edit operations. */
	if (m_e->modifiers() != Qt::NoModifier)
		return nullptr;

	return m_e;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
