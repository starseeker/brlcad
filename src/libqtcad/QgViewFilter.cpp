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

extern "C" {
#include "bsg.h"
#include "bsg/view_state.h"
}

#include "qtcad/QgObolSnap.h"
#include "qtcad/QgViewFilter.h"

static uint32_t
qg_obol_snap_kinds_from_bsg(bsg_snap_kind_mask kinds)
{
	uint32_t obol_kinds = 0;

	if (kinds & BSG_SNAP_KIND_ENDPOINT)
		obol_kinds |= QgObolSnapRecord::ENDPOINT |
			QgObolSnapRecord::MIDPOINT |
			QgObolSnapRecord::LINE_NEAREST;
	if (kinds & BSG_SNAP_KIND_MIDPOINT)
		obol_kinds |= QgObolSnapRecord::MIDPOINT;
	if (kinds & BSG_SNAP_KIND_OVERLAY_HANDLE)
		obol_kinds |= QgObolSnapRecord::CENTER;

	return obol_kinds;
}

static int
qg_obol_db_snap_enabled(const struct bsg_view *v)
{
	if (!v || !bsg_view_snap_lines(v))
		return 0;

	int flags = bsg_view_snap_source_flags(v);
	if (flags == BSG_SNAP_TCL)
		return 0;
	return !flags || (flags & BSG_SNAP_DB);
}

static float
qg_obol_snap_tolerance(const struct bsg_view *v)
{
	if (!v || v->gv_width <= 0 || v->gv_height <= 0)
		return 0.0f;

	double lavg = ((double)v->gv_width + (double)v->gv_height) * 0.5;
	double lratio = 1.0 / lavg;
	double tol = v->gv_size * lratio * bsg_view_snap_tolerance_factor(v);

	return tol > 0.0 ? (float)tol : 0.0f;
}

static void
qg_obol_refine_db_snap(QgView *display, struct bsg_view *v)
{
	if (!display || !qg_obol_db_snap_enabled(v))
		return;

	uint32_t obol_kinds = qg_obol_snap_kinds_from_bsg(
		bsg_view_snap_kind_mask(v));
	float tolerance = qg_obol_snap_tolerance(v);
	if (!obol_kinds || tolerance <= 0.0f)
		return;

	QgObolSnapRecord record;
	SbVec3f query((float)v->gv_point[X],
		(float)v->gv_point[Y],
		(float)v->gv_point[Z]);
	if (!qg_obol_snap_point(display, query, tolerance, obol_kinds, record))
		return;

	v->gv_point[X] = (fastf_t)record.point[0];
	v->gv_point[Y] = (fastf_t)record.point[1];
	v->gv_point[Z] = (fastf_t)record.point[2];
}

class QgView;

class QgViewFilter::QgViewFilterPrivate {
public:
	struct bsg_view *v = nullptr;
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
QgViewFilter::set_view(struct bsg_view *nv)
{
	m->v = nv;
}

struct bsg_view *
QgViewFilter::view() const
{
	return m->v;
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
	if (!m->v)
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

	/* Keep bsg_view mouse state synchronized with the event stream. */
	m->v->gv_prevMouseX = m->v->gv_mouse_x;
	m->v->gv_prevMouseY = m->v->gv_mouse_y;
	m->v->gv_mouse_x = e_x;
	m->v->gv_mouse_y = e_y;
	bsg_screen_pt(&m->v->gv_point, (fastf_t)e_x, (fastf_t)e_y, m->v);
	qg_obol_refine_db_snap(m->display, m->v);

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
