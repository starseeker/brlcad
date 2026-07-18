/*                   Q G P O L Y F I L T E R . C P P
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
/** @file QgPolyFilter.cpp
 *
 * Polygon interaction logic for Qt views.
 *
 */

#include "common.h"

#include "bv.h"

extern "C" {
#include "bg/polygon.h"
#include "raytrace.h" // For finalize polygon sketch export functionality (TODO - need to move...)
}

#include "qtcad/QgPolyFilter.h"
#include "qtcad/QgSignalFlags.h"
#include "qtcad/QgView.h"

static void *
qg_poly_filter_view(const QgPolyFilter *filter)
{
	QgView *display = filter ? filter->view_widget() : nullptr;
	return display ? display->viewContext() : nullptr;
}

static void
qg_poly_mouse_xy(QMouseEvent *m_e, int *sx, int *sy)
{
	if (!m_e || !sx || !sy)
		return;

#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
	*sx = m_e->x();
	*sy = m_e->y();
#else
	*sx = (int)m_e->position().x();
	*sy = (int)m_e->position().y();
#endif
}

static void
qg_poly_screen_point(void *v, QMouseEvent *m_e, point_t point)
{
	VSETALL(point, 0.0);
	if (!v || !m_e)
		return;

	int sx = 0;
	int sy = 0;
	qg_poly_mouse_xy(m_e, &sx, &sy);
	(void)bv_screen_to_model(point, bv_context_view_const(static_cast<const struct bv_context *>(v)),
		(fastf_t)sx, (fastf_t)sy);
}

static int
qg_poly_update_from_event(qg_polygon_ref ref, void *v,
	QMouseEvent *m_e, qg_polygon_update_mode utype)
{
	if (!m_e)
		return 0;

	int sx = 0;
	int sy = 0;
	qg_poly_mouse_xy(m_e, &sx, &sy);
	return ged_polygon_update_screen_pt(ref,
		ged_view_context_from_bv(static_cast<struct bv_context *>(v)),
		sx, sy, utype);
}

bool
QgPolyFilter::close_polygon()
{
	// Close the general polygon - if that's what we're creating,
	// at this point it will still be open.
	qg_polygon_record rec;
	if (ged_polygon_record_get(polygon, &rec) && rec.first_contour_open) {
		if (!ged_polygon_close(polygon)) {
			polygon = GED_POLYGON_REF_NULL;
			return false;
		}
	}

	return true;
}

bool
QgPolyCreateFilter::eventFilter(QObject *, QEvent *e)
{
	QMouseEvent *m_e = view_sync(e);
	if (!m_e)
		return false;

	void *v = qg_poly_filter_view(this);
	if (!v)
		return false;

	// Handle Left Click
	if (m_e->type() == QEvent::MouseButtonPress && m_e->buttons().testFlag(Qt::LeftButton)) {

		if (ged_polygon_ref_is_null(polygon)) {

			point_t current_point = VINIT_ZERO;
			qg_poly_screen_point(v, m_e, current_point);

			polygon = ged_polygon_create(
				ged_view_context_from_bv(
				    static_cast<struct bv_context *>(v)),
				"_tmp_view_polygon", 0,
				ptype, current_point);
			ged_polygon_set_visual(polygon, &edge_color, &fill_color, fill_slope_x, fill_slope_y, fill_density, vZ, fill_poly ? 1 : 0);

			emit view_updated(QG_VIEW_REFRESH);
			return true;
		}

		// If we don't have a polygon at this point, we're done - subsequent logic assumes it
		if (ged_polygon_ref_is_null(polygon))
			return true;

		// If we are in the process of creating a general polygon, after the initial creation
		// left clicks will append new points
		qg_polygon_record rec;
		if (ged_polygon_record_get(polygon, &rec) && rec.type == QG_POLYGON_GENERAL) {
			qg_poly_update_from_event(polygon, v, m_e, QG_POLYGON_UPDATE_PT_APPEND);
			emit view_updated(QG_VIEW_REFRESH);
			return true;
		}

		// When we're dealing with polygons stray left clicks shouldn't zoom - just
		// consume them if we're not using them above.
		return true;
	}

	if (m_e->type() == QEvent::MouseButtonPress && m_e->buttons().testFlag(Qt::RightButton)) {
		// No-op if we're not in the process of creating a polygon
		if (ged_polygon_ref_is_null(polygon))
			return true;

		// Non-general polygon creation doesn't use right click.
		qg_polygon_record rec;
		if (!ged_polygon_record_get(polygon, &rec) || rec.type != QG_POLYGON_GENERAL)
			return true;

		// When creating a general polygon, right click indicates we're done.
		finalize(true);

		return true;
	}

	if (m_e->type() == QEvent::MouseButtonPress) {
		// We don't want other stray mouse clicks to do something surprising
		return true;
	}

	// During initial add/creation of non-general polygons, mouse movement
	// adjusts the shape
	if (m_e->type() == QEvent::MouseMove) {
		// No-op if no current polygon is defined
		if (ged_polygon_ref_is_null(polygon))
			return true;

		// General polygon creation doesn't use mouse movement.
		qg_polygon_record rec;
		if (ged_polygon_record_get(polygon, &rec) && rec.type == QG_POLYGON_GENERAL)
			return true;

		// For every other polygon type, call the libbv update routine
		// with the view's x,y coordinates
		if (m_e->buttons().testFlag(Qt::LeftButton) && m_e->modifiers() == Qt::NoModifier) {
			qg_poly_update_from_event(polygon, v, m_e, QG_POLYGON_UPDATE_DEFAULT);
			emit view_updated(QG_VIEW_REFRESH);
			return true;
		}
	}

	// For the constrained polygon shapes, we're done creating once we release
	// the mouse button (i.e. a "click, hold and move" creation paradigm)
	if (m_e->type() == QEvent::MouseButtonRelease) {

		// No-op if no current polygon is defined
		if (ged_polygon_ref_is_null(polygon))
			return true;

		// General polygons are finalized by a right-click close, since
		// appending multiple points requires multiple mouse click-and-release
		// operations
		qg_polygon_record rec;
		if (ged_polygon_record_get(polygon, &rec) && rec.type == QG_POLYGON_GENERAL)
			return true;

		// For all non-general polygons, mouse release is the signal
		// to finish up.
		finalize(true);
		polygon = GED_POLYGON_REF_NULL;

		return true;
	}

	return false;
}

void
QgPolyCreateFilter::finalize(bool)
{
	int icnt = 0;

	if (ged_polygon_ref_is_null(polygon))
		return;

	if (!close_polygon())
		return;

	if (op == bg_None || bool_objs.empty()) {
		// No interactions, so we're keeping it - assign a proper name
		ged_polygon_set_name(polygon, vname.c_str());
	}
	else {

		for (auto target : bool_objs) {
			icnt += ged_polygon_csg(target, polygon, op);
		}

		// When doing boolean operations, the convention is if there were one
		// or more interactions with other polygons, the original polygon is
		// not retained
		if (icnt || op == bg_Difference || op == bg_Intersection) {
			ged_polygon_remove(polygon);
			polygon = GED_POLYGON_REF_NULL;
		}
		else {
			// No interactions, so we're keeping it - assign a proper name
			ged_polygon_set_name(polygon, vname.c_str());
		}
	}

	emit view_updated(QG_VIEW_REFRESH);
	emit finalized((icnt > 0) ? true : false);
}

bool
QgPolyUpdateFilter::eventFilter(QObject *, QEvent *e)
{
	QMouseEvent *m_e = view_sync(e);
	if (!m_e)
		return false;

	// The update filter needs an active polygon to operate on
	if (ged_polygon_ref_is_null(polygon))
		return false;

	// We don't want other stray mouse clicks to do something surprising
	if (m_e->type() == QEvent::MouseButtonPress || m_e->type() == QEvent::MouseButtonRelease) {
		return true;
	}

	if (m_e->type() == QEvent::MouseMove) {

		// General polygon creation doesn't use mouse movement.
		qg_polygon_record rec;
		if (ged_polygon_record_get(polygon, &rec) && rec.type == QG_POLYGON_GENERAL)
			return true;

		// For every other polygon type, call the libbv update routine
		// with the view's x,y coordinates
		if (m_e->buttons().testFlag(Qt::LeftButton) && m_e->modifiers() == Qt::NoModifier) {
			qg_poly_update_from_event(polygon, qg_poly_filter_view(this), m_e, QG_POLYGON_UPDATE_DEFAULT);
			emit view_updated(QG_VIEW_REFRESH);
			return true;
		}

	}

	return false;
}

bool
QgPolySelectFilter::eventFilter(QObject *, QEvent *e)
{
	QMouseEvent *m_e = view_sync(e);
	if (!m_e)
		return false;

	void *v = qg_poly_filter_view(this);
	if (!v)
		return false;

	// Handle Left Click
	if (m_e->type() == QEvent::MouseButtonPress && m_e->buttons().testFlag(Qt::LeftButton)) {
		/* Use typed polygon selection records rather than the legacy
		 * feature-table and ptbl compatibility path. */
		point_t current_point = VINIT_ZERO;
		qg_poly_screen_point(v, m_e, current_point);
		polygon = ged_polygon_select(
			ged_view_context_from_bv(
			    static_cast<struct bv_context *>(v)), current_point);
		if (ged_polygon_ref_is_null(polygon))
			return true;
		qg_polygon_record rec;
		if (ged_polygon_record_get(polygon, &rec)) {
			ptype = rec.type;
			close_general_poly = rec.first_contour_open;
		}

		return true;
	}

	// We also don't want other stray mouse clicks to do something surprising
	if (m_e->type() == QEvent::MouseButtonPress || m_e->type() == QEvent::MouseButtonRelease)
		return true;

	return false;
}

bool
QgPolyPointFilter::eventFilter(QObject *, QEvent *e)
{
	QMouseEvent *m_e = view_sync(e);
	if (!m_e)
		return false;

	// The point filter needs an active general polygon to operate on
	if (ged_polygon_ref_is_null(polygon) || ptype != QG_POLYGON_GENERAL)
		return false;

	qg_polygon_record rec;
	if (!ged_polygon_record_get(polygon, &rec))
		return false;

	// If we have a Left release, clear point selection
	if (m_e->type() == QEvent::MouseButtonRelease) {
		qg_poly_update_from_event(polygon, qg_poly_filter_view(this), m_e, QG_POLYGON_UPDATE_PT_SELECT_CLEAR);
		emit view_updated(QG_VIEW_REFRESH);
		return true;
	}

	// Left press selects a point
	if (m_e->type() == QEvent::MouseButtonPress && m_e->buttons().testFlag(Qt::LeftButton)) {
		qg_poly_update_from_event(polygon, qg_poly_filter_view(this), m_e, QG_POLYGON_UPDATE_PT_SELECT);
		emit view_updated(QG_VIEW_REFRESH);
		return true;
	}

	// We also don't want other stray mouse clicks to do something surprising
	if (m_e->type() == QEvent::MouseButtonPress || m_e->type() == QEvent::MouseButtonRelease)
		return true;

	// Handle Mouse Move - move selected point with a left button click-and-hold
	if (m_e->type() == QEvent::MouseMove) {
		if (rec.curr_point_i < 0) {
			// No selected point
			return true;
		}
		if (m_e->buttons().testFlag(Qt::LeftButton) && m_e->modifiers() == Qt::NoModifier) {
			qg_poly_update_from_event(polygon, qg_poly_filter_view(this), m_e, QG_POLYGON_UPDATE_PT_MOVE);
			emit view_updated(QG_VIEW_REFRESH);
			return true;
		}

		return true;
	}

	return false;
}

bool
QgPolyMoveFilter::eventFilter(QObject *, QEvent *e)
{
	QMouseEvent *m_e = view_sync(e);
	if (!m_e)
		return false;

	void *v = qg_poly_filter_view(this);
	if (!v)
		return false;

	// The move filter needs an active polygon to operate on
	if (ged_polygon_ref_is_null(polygon) && move_objs.empty())
		return false;

	point_t current_point = VINIT_ZERO;
	qg_poly_screen_point(v, m_e, current_point);

	// We don't want other stray mouse clicks to do something surprising
	if (m_e->type() == QEvent::MouseButtonPress || m_e->type() == QEvent::MouseButtonRelease) {
		VMOVE(m_prev_point, current_point);
		m_prev_point_valid = true;
		return true;
	}

	// If we're clicking-and-holding, it's time to move
	if (m_e->type() == QEvent::MouseMove) {
		if (!m_prev_point_valid) {
			VMOVE(m_prev_point, current_point);
			m_prev_point_valid = true;
		}
		if (m_e->buttons().testFlag(Qt::LeftButton) && m_e->modifiers() == Qt::NoModifier) {
			if (!move_objs.empty()) {
				for (auto mpoly : move_objs)
					ged_polygon_move(mpoly, &current_point, &m_prev_point);
			}
			else {
				ged_polygon_move(polygon, &current_point, &m_prev_point);
			}
			emit view_updated(QG_VIEW_REFRESH);
		}
		VMOVE(m_prev_point, current_point);
		return true;
	}

	return false;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
