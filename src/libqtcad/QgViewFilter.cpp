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
#include "qtcad/QgObolSnap.h"
#include "qtcad/QgLegacyView.h"
#include "qtcad/QgView.h"
#include "qtcad/QgViewFilter.h"
#include "QgLegacyViewContext.h"

static uint32_t
qg_obol_snap_kinds_from_bv_mask(unsigned long long kinds, int source_flags)
{
    uint32_t obol_kinds = 0;
    const int all_sources = !source_flags;
    const int allow_db = all_sources || (source_flags & BV_SNAP_DB);
    const int allow_view = all_sources ||
			   (source_flags & (BV_SNAP_VIEW | BV_SNAP_LOCAL));

    if (allow_view && (kinds & BV_SNAP_KIND_GRID))
	obol_kinds |= QgObolSnapRecord::GRID;
    if (allow_db && (kinds & BV_SNAP_KIND_ENDPOINT))
	obol_kinds |= QgObolSnapRecord::ENDPOINT |
		      QgObolSnapRecord::MIDPOINT |
		      QgObolSnapRecord::LINE_NEAREST;
    if (allow_db && (kinds & BV_SNAP_KIND_MIDPOINT))
	obol_kinds |= QgObolSnapRecord::MIDPOINT;
    if (allow_view && (kinds & BV_SNAP_KIND_OVERLAY_HANDLE))
	obol_kinds |= QgObolSnapRecord::CENTER;

    return obol_kinds;
}

static const struct bv *
qg_obol_snap_bv_const(const qg_legacy_view *v)
{
    const void *view_ctx = qg_legacy_view_to_context(v);
    return view_ctx ? bv_context_view_const(
	reinterpret_cast<const struct bv_context *>(view_ctx)) : nullptr;
}

static struct bv *
qg_obol_snap_bv(qg_legacy_view *v)
{
    void *view_ctx = qg_legacy_view_to_context(v);
    return view_ctx ? bv_context_view(
	reinterpret_cast<struct bv_context *>(view_ctx)) : nullptr;
}

static int
qg_obol_snap_enabled(const qg_legacy_view *v)
{
    const struct bv *view = qg_obol_snap_bv_const(v);
    if (!view || !bv_snap_lines_get(view))
	return 0;

    int flags = bv_snap_source_flags_get(view);
    if (flags == BV_SNAP_TCL)
	return 0;

    return !flags || (flags & (BV_SNAP_DB | BV_SNAP_VIEW |
			       BV_SNAP_LOCAL));
}

static float
qg_obol_snap_tolerance(const qg_legacy_view *v)
{
    if (!v)
	return 0.0f;

    const struct bv *view = qg_obol_snap_bv_const(v);
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
qg_obol_refine_db_snap(QgView *display, qg_legacy_view *v)
{
    if (!display || !qg_obol_snap_enabled(v))
	return;

    struct bv *view = qg_obol_snap_bv(v);
    const int source_flags = bv_snap_source_flags_get(view);
    uint32_t obol_kinds = qg_obol_snap_kinds_from_bv_mask(
			      bv_snap_kind_mask_get(view), source_flags);
    float tolerance = qg_obol_snap_tolerance(v);
    if (!obol_kinds || tolerance <= 0.0f)
	return;

    QgObolSnapRecord record;
    point_t query_point = VINIT_ZERO;
    if (!bv_current_point_get(query_point, view))
	return;

    SbVec3f query((float)query_point[X],
		  (float)query_point[Y],
		  (float)query_point[Z]);
    if (!qg_obol_snap_point(display, query, tolerance, obol_kinds, record))
	return;

    point_t snapped_point = VINIT_ZERO;
    VSET(snapped_point, (fastf_t)record.point[0],
	 (fastf_t)record.point[1], (fastf_t)record.point[2]);
    (void)bv_current_point_set(view, snapped_point);
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

    /* Keep neutral view state synchronized with the event stream. */
    bv_mouse_state_set(qg_obol_snap_bv(lv), e_x, e_y);
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
