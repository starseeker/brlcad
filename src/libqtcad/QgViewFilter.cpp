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

#include "bv.h"
#include "BObol/BSnapAction.h"
#include "ged/view_types.h"
#include "qtcad/QgObolSnap.h"
#include "qtcad/QgView.h"
#include "qtcad/QgViewFilter.h"

#include <float.h>

static uint32_t
qg_obol_snap_kinds_from_bv_mask(unsigned long long kinds)
{
    uint32_t obol_kinds = 0;

    if (kinds & BV_SNAP_KIND_ENDPOINT)
	obol_kinds |= QgObolSnapRecord::ENDPOINT |
		      QgObolSnapRecord::MIDPOINT |
		      QgObolSnapRecord::LINE_NEAREST;
    if (kinds & BV_SNAP_KIND_MIDPOINT)
	obol_kinds |= QgObolSnapRecord::MIDPOINT;
    if (kinds & BV_SNAP_KIND_OVERLAY_HANDLE)
	obol_kinds |= QgObolSnapRecord::CENTER;

    return obol_kinds;
}

static uint32_t
qg_obol_source_filter_from_bv(int source_flags)
{
    if (!source_flags)
	return SoBRLSnapAction::ALL_SOURCES;

    uint32_t sources = 0;
    if (source_flags & BV_SNAP_DB)
	sources |= SoBRLSnapAction::DATABASE_SOURCES;
    if (source_flags & BV_SNAP_VIEW)
	sources |= SoBRLSnapAction::VIEW_SOURCES;
    if (source_flags & BV_SNAP_LOCAL)
	sources |= SoBRLSnapAction::LOCAL_SOURCES;
    if (source_flags & BV_SNAP_SHARED)
	sources |= SoBRLSnapAction::SHARED_SOURCES;
    return sources;
}

static const struct bv *
qg_obol_snap_bv_const(const struct ged_view_context *view_ctx)
{
    return view_ctx ? bv_context_view_const(
	ged_view_context_bv_const(view_ctx)) : nullptr;
}

static struct bv *
qg_obol_snap_bv(struct ged_view_context *view_ctx)
{
    return view_ctx ? bv_context_view(
	ged_view_context_bv(view_ctx)) : nullptr;
}

static int
qg_obol_snap_enabled(const struct ged_view_context *view_ctx)
{
    const struct bv *view = qg_obol_snap_bv_const(view_ctx);
    if (!view || !bv_snap_lines_get(view))
	return 0;

    int flags = bv_snap_source_flags_get(view);
    if (flags == BV_SNAP_TCL)
	return 0;

    return !flags || (flags & (BV_SNAP_DB | BV_SNAP_VIEW |
			       BV_SNAP_LOCAL | BV_SNAP_SHARED));
}

static int
qg_grid_snap_enabled(const struct bv *view)
{
    struct bv_grid_state grid = BV_GRID_STATE_INIT;
    if (!view || !bv_grid_state_get(&grid, view) || !grid.snap ||
	!(bv_snap_kind_mask_get(view) & BV_SNAP_KIND_GRID))
	return 0;

    const int flags = bv_snap_source_flags_get(view);
    return flags != BV_SNAP_TCL &&
	(!flags || (flags & (BV_SNAP_VIEW | BV_SNAP_LOCAL |
			    BV_SNAP_SHARED)));
}

static int
qg_grid_snap_model_point(point_t snapped, const struct bv *view,
	const point_t query)
{
    mat_t model2view;
    mat_t view2model;
    point_t view_point = VINIT_ZERO;
    if (!snapped || !view || !query ||
	!bv_model2view_get(model2view, view) ||
	!bv_view2model_get(view2model, view))
	return 0;

    MAT4X3PNT(view_point, model2view, query);
    if (!bv_snap_grid_2d(view, &view_point[X], &view_point[Y]))
	return 0;
    MAT4X3PNT(snapped, view2model, view_point);
    return 1;
}

static float
qg_obol_snap_tolerance(const struct ged_view_context *view_ctx)
{
	if (!view_ctx)
	return 0.0f;

    const struct bv *view = qg_obol_snap_bv_const(view_ctx);
    int width = bv_width_get(view);
    int height = bv_height_get(view);
    if (width <= 0 || height <= 0)
	return 0.0f;

    double lavg = ((double)width + (double)height) * 0.5;
    double lratio = 1.0 / lavg;
    double tol = bv_size_get(view) * lratio *
		 bv_snap_tolerance_factor_get(view);

    return tol > 0.0 ? (float)tol : 0.0f;
}

static void
qg_obol_refine_snap(QgView *display, struct ged_view_context *view_ctx)
{
    if (!display || !view_ctx)
	return;

    struct bv *view = qg_obol_snap_bv(view_ctx);
    const float tolerance = qg_obol_snap_tolerance(view_ctx);
    if (tolerance <= 0.0f)
	return;

    point_t query_point = VINIT_ZERO;
    if (!bv_current_point_get(query_point, view))
	return;

    point_t best_point = VINIT_ZERO;
    VMOVE(best_point, query_point);
    fastf_t best_distance = FLT_MAX;
    int found = 0;

    if (qg_grid_snap_enabled(view)) {
	point_t grid_point = VINIT_ZERO;
	if (qg_grid_snap_model_point(grid_point, view, query_point)) {
	    const fastf_t distance = DIST_PNT_PNT(query_point, grid_point);
	    if (distance <= (fastf_t)tolerance) {
		VMOVE(best_point, grid_point);
		best_distance = distance;
		found = 1;
	    }
	}
    }

    if (qg_obol_snap_enabled(view_ctx)) {
	const int source_flags = bv_snap_source_flags_get(view);
	const uint32_t source_filter =
	    qg_obol_source_filter_from_bv(source_flags);
	const uint32_t obol_kinds = qg_obol_snap_kinds_from_bv_mask(
	    bv_snap_kind_mask_get(view));
	QgObolSnapRecord record;
	if ((source_filter || !source_flags) && obol_kinds) {
	    SbVec3f query((float)query_point[X], (float)query_point[Y],
		    (float)query_point[Z]);
	    if (qg_obol_snap_point_filtered(display, query, tolerance,
		    obol_kinds, source_filter, record) &&
		    (!found || record.distance <= best_distance + 1.0e-6f)) {
		VSET(best_point, (fastf_t)record.point[0],
		     (fastf_t)record.point[1], (fastf_t)record.point[2]);
		found = 1;
	    }
	}
    }

    if (found)
	(void)bv_current_point_set(view, best_point);
}

class QgViewFilter::QgViewFilterPrivate
{
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
	struct ged_view_context *view_ctx = m->display ?
	    ged_view_context_from_bv(m->display->viewContext()) : nullptr;
	if (!view_ctx)
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

    /* Keep neutral view state synchronized with the event stream. */
	bv_mouse_state_set(qg_obol_snap_bv(view_ctx), e_x, e_y);
	qg_obol_refine_snap(m->display, view_ctx);

    /* Modifier keys are typically view-nav gestures, not edit operations. */
    if (m_e->modifiers() != Qt::NoModifier)
	return nullptr;

    return m_e;
}

bool
QgViewFilter::view_sync(const BObolInputEvent *event)
{
    struct ged_view_context *view_ctx = m->display ?
	ged_view_context_from_bv(m->display->viewContext()) : nullptr;
    if (!view_ctx || !event)
	return false;
    if (event->type != BOBOL_INPUT_POINTER_PRESS &&
	event->type != BOBOL_INPUT_POINTER_RELEASE &&
	event->type != BOBOL_INPUT_POINTER_MOTION)
	return false;

    bv_mouse_state_set(qg_obol_snap_bv(view_ctx), event->x, event->y);
    qg_obol_refine_snap(m->display, view_ctx);

    return event->modifiers == BOBOL_INPUT_MOD_NONE;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
