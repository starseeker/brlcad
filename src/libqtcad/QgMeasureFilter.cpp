/*                 Q G M E A S U R E F I L T E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2021-2026 United States Government as represented by
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
/** @file QgMeasureFilter.cpp
 *
 * Measurement tool for Qt views.
 *
 */

#include "common.h"

#include "bv.h"

#include <string>
#include "qtcad/QgMeasureFilter.h"
#include "qtcad/QgObolMeasure.h"
#include "qtcad/QgSignalFlags.h"
#include "qtcad/QgView.h"

static struct ged_view_context *
qg_measure_filter_view(const QgMeasureFilter *filter)
{
	QgView *display = filter ? filter->view_widget() : nullptr;
	return display ? ged_view_context_from_bv(display->viewContext()) : nullptr;
}

void
QgMeasureFilter::update_current_mouse(QMouseEvent *m_e)
{
	current_mouse_valid = false;
	if (!m_e)
		return;

#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
	current_mouse_x = m_e->x();
	current_mouse_y = m_e->y();
#else
	current_mouse_x = (int)m_e->position().x();
	current_mouse_y = (int)m_e->position().y();
#endif
	current_mouse_valid = true;
}

bool
QgMeasureFilter::current_mouse_xy(int *sx, int *sy) const
{
	if (!sx || !sy || !current_mouse_valid)
		return false;

	*sx = current_mouse_x;
	*sy = current_mouse_y;
	return true;
}

void
QgMeasureFilter::clear_measure_overlay()
{
	obol_overlay_point_count = 0;
	(void)qg_obol_measure_clear_overlay(view_widget(), oname.c_str());
}

void
QgMeasureFilter::update_measure_overlay(int point_count)
{
	if (point_count <= 0)
		return;

	if (point_count > 3)
		point_count = 3;
	obol_overlay_point_count = point_count;

	SbVec3f points[3] = {
		SbVec3f((float)p1[X], (float)p1[Y], (float)p1[Z]),
		SbVec3f((float)p2[X], (float)p2[Y], (float)p2[Z]),
		SbVec3f((float)p3[X], (float)p3[Y], (float)p3[Z])
	};

	(void)qg_obol_measure_update_overlay(view_widget(), oname.c_str(),
		points, point_count, &measure_color);
}

double
QgMeasureFilter::length1()
{
	return DIST_PNT_PNT(p1, p2);
}

double
QgMeasureFilter::length2()
{
	if (mode < 3)
		return 0.0;

	return DIST_PNT_PNT(p2, p3);
}

double
QgMeasureFilter::angle(bool radians)
{
	if (mode < 3)
		return 0.0;

	vect_t v1, v2;
	VSUB2(v1, p1, p2);
	VSUB2(v2, p3, p2);
	VUNITIZE(v1);
	VUNITIZE(v2);
	double a = acos(VDOT(v1, v2));
	if (radians)
		return a*180/M_PI;
	return a;
}

void
QgMeasureFilter::update_color(struct bu_color *c)
{
	if (!c)
		return;

	BU_COLOR_CPY(&measure_color, c);
	if (obol_overlay_point_count > 0)
		update_measure_overlay(obol_overlay_point_count);
}

bool
QgMeasureFilter::eventFilter(QObject *, QEvent *e)
{
	QMouseEvent *m_e = view_sync(e);
	if (!m_e)
		return false;
	update_current_mouse(m_e);

	struct ged_view_context *view_ctx = qg_measure_filter_view(this);
	if (!view_ctx)
		return false;

	if (e->type() == QEvent::MouseButtonPress) {
		if (m_e->button() == Qt::RightButton) {
			clear_measure_overlay();
			mode = 0;
			VSETALL(p1, 0.0);
			VSETALL(p2, 0.0);
			VSETALL(p3, 0.0);
			emit view_updated(QG_VIEW_REFRESH);
			return true;
		}
		if (mode == 4) {
			clear_measure_overlay();
			mode = 0;
			emit view_updated(QG_VIEW_REFRESH);
			return true;
		}
		if (!mode) {
			if (!get_point())
				return true;

			VSETALL(p1, 0.0);
			VSETALL(p2, 0.0);
			VSETALL(p3, 0.0);

			mode = 1;
			VMOVE(p1, mpnt);
			VMOVE(p2, mpnt);
			update_measure_overlay(1);
			emit view_updated(QG_VIEW_REFRESH);
			return true;
		}
		if (mode == 1) {
			if (!get_point())
				return true;

			VMOVE(p2, mpnt);
			return true;
		}
		if (mode == 2) {
			if (!get_point())
				return true;
			mode = 3;
			VMOVE(p3, mpnt);
			update_measure_overlay(3);
			emit view_updated(QG_VIEW_REFRESH);
		}
		return true;
	}

	if (e->type() == QEvent::MouseMove) {
		if (!mode)
			return false;
		if (mode == 1) {
			if (!get_point())
				return true;

			VMOVE(p2, mpnt);
			update_measure_overlay(2);
			emit view_updated(QG_VIEW_REFRESH);
		}
		if (mode == 3) {
			if (!get_point())
				return true;

			VMOVE(p3, mpnt);
			update_measure_overlay(3);
			emit view_updated(QG_VIEW_REFRESH);
		}
		return true;
	}

	if (e->type() == QEvent::MouseButtonRelease) {
		if (m_e->button() == Qt::RightButton) {
			mode = 0;
			clear_measure_overlay();
			emit view_updated(QG_VIEW_REFRESH);
			return true;
		}
		if (!mode)
			return false;
		if (mode == 1 && DIST_PNT_PNT(p1, p2) < SMALL_FASTF) {
			return true;
		}
		if (mode == 1) {
			if (!get_point())
				return true;

			if (length_only) {
				// Angle measurement disabled, starting over
				mode = 0;
				emit view_updated(QG_VIEW_REFRESH);
				return true;
			}

			mode = 2;
			VMOVE(p2, mpnt);
			update_measure_overlay(2);
			emit view_updated(QG_VIEW_REFRESH);
			return true;
		}
		if (mode == 3) {
			if (!get_point())
				return true;
			mode = 4;
			VMOVE(p3, mpnt);
			update_measure_overlay(3);
			emit view_updated(QG_VIEW_REFRESH);
			return true;
		}

		return true;
	}

	// Shouldn't get here
	return false;
}

bool
QMeasure2DFilter::get_point()
{
	struct ged_view_context *view_ctx = qg_measure_filter_view(this);
	if (!view_ctx)
		return false;
	int sx = 0, sy = 0;
	if (!current_mouse_xy(&sx, &sy))
		return false;
	const struct bv *view = bv_context_view_const(
	    ged_view_context_bv_const(view_ctx));
	fastf_t vx, vy;
	if (!bv_screen_to_view(&vx, &vy, view, sx, sy))
		return false;
	point_t vpnt;
	mat_t view2model;
	VSET(vpnt, vx, vy, 0);
	if (!bv_view2model_get(view2model, view))
		return false;
	MAT4X3PNT(mpnt, view2model, vpnt);
	return true;
}

bool
QMeasure2DFilter::eventFilter(QObject *o, QEvent *e)
{
	return QgMeasureFilter::eventFilter(o, e);
}

QMeasure3DFilter::QMeasure3DFilter()
{
}

QMeasure3DFilter::~QMeasure3DFilter()
{
}

bool
QMeasure3DFilter::get_point()
{
	struct ged_view_context *view_ctx = qg_measure_filter_view(this);
	if (!view_ctx)
		return false;

	SbVec3f obolPoint;
	int sx = 0, sy = 0;
	if (!current_mouse_xy(&sx, &sy))
		return false;
	if (qg_obol_measure_pick_point(view_widget(), sx, sy, obolPoint, NULL)) {
		VSET(mpnt, obolPoint[0], obolPoint[1], obolPoint[2]);
		return true;
	}

	return false;
}

bool
QMeasure3DFilter::eventFilter(QObject *o, QEvent *e)
{
	return QgMeasureFilter::eventFilter(o, e);
}




// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
