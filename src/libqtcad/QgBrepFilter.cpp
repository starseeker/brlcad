/*                  Q G B R E P F I L T E R . C P P
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this file; see the file named COPYING for more information.
 */
/** @file QgBrepFilter.cpp
 *
 * Qt B-rep NURBS CV editing filters.
 */

#include "common.h"

#include <algorithm>
#include <cfloat>

#include "bn/mat.h"
#include "bv.h"
#include "brep/edit.h"
#include "raytrace.h"
#include "rt/geom.h"
#include "rt/rt_ecmds.h"

#include "qtcad/QgBrepFilter.h"
#include "qtcad/QgSignalFlags.h"
#include "qtcad/QgView.h"

static void
qg_brep_mouse_xy(QMouseEvent *event, int *sx, int *sy)
{
    if (!event || !sx || !sy)
	return;
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    *sx = event->x();
    *sy = event->y();
#else
    *sx = (int)event->position().x();
    *sy = (int)event->position().y();
#endif
}

static const struct bv *
qg_brep_view(const QgBrepFilter *filter)
{
    QgView *display = filter ? filter->view_widget() : nullptr;
    return display
	? bv_context_view_const(
	    static_cast<const struct bv_context *>(display->viewContext()))
	: nullptr;
}

static bool
qg_brep_selected_point(const struct rt_edit *es, int face, int cv_u,
	int cv_v, ON_3dPoint *point)
{
    if (!es || !point)
	return false;
    const struct rt_brep_internal *bip =
	(const struct rt_brep_internal *)es->es_int.idb_ptr;
    if (!bip || !bip->brep || face < 0 || face >= bip->brep->m_F.Count())
	return false;
    const ON_NurbsSurface *surface =
	dynamic_cast<const ON_NurbsSurface *>(
	    bip->brep->m_F[face].SurfaceOf());
    return surface && surface->GetCV(cv_u, cv_v, *point);
}

QgBrepFilter::QgBrepFilter(QObject *parent)
    : QgViewFilter(parent)
{
}

bool
QgBrepFilter::selected_cv(int *face, int *cv_u, int *cv_v) const
{
    if (!es || !face || !cv_u || !cv_v)
	return false;

    fastf_t values[3] = {0.0, 0.0, 0.0};
    int count = EDOBJ[es->es_int.idb_type].ft_edit_get_params(
	    es, ECMD_BREP_SRF_SELECT, values);
    if (count != 3 || values[0] < 0.0 || values[1] < 0.0 || values[2] < 0.0)
	return false;

    *face = (int)values[0];
    *cv_u = (int)values[1];
    *cv_v = (int)values[2];
    return true;
}

bool
QgBrepFilter::selected_cv_is_topology_safe() const
{
    int face = -1;
    int cv_u = -1;
    int cv_v = -1;
    if (!selected_cv(&face, &cv_u, &cv_v))
	return false;

    const struct rt_brep_internal *bip =
	(const struct rt_brep_internal *)es->es_int.idb_ptr;
    return bip && bip->brep
	&& brep_face_cv_is_topology_safe(bip->brep, face, cv_u, cv_v);
}

bool
QgBrepFilter::selected_cv_can_translate() const
{
    int face = -1;
    int cv_u = -1;
    int cv_v = -1;
    if (!selected_cv(&face, &cv_u, &cv_v))
	return false;

    const struct rt_brep_internal *bip =
	(const struct rt_brep_internal *)es->es_int.idb_ptr;
    return bip && bip->brep
	&& brep_face_cv_can_translate(bip->brep, face, cv_u, cv_v);
}

bool
QgBrepFilter::pick_cv(int sx, int sy, fastf_t max_px)
{
    const struct bv *view = qg_brep_view(this);
    if (!es || !view || max_px <= 0.0)
	return false;

    const struct rt_brep_internal *bip =
	(const struct rt_brep_internal *)es->es_int.idb_ptr;
    if (!bip || !bip->brep)
	return false;

    fastf_t cursor_x = 0.0;
    fastf_t cursor_y = 0.0;
    if (!bv_screen_to_view(&cursor_x, &cursor_y, view,
		(fastf_t)sx, (fastf_t)sy))
	return false;

    const int width = bv_width_get(view);
    const fastf_t px_to_view = width > 0 ? 2.0 / (fastf_t)width : 0.005;
    const fastf_t max_dist_sq =
	(max_px * px_to_view) * (max_px * px_to_view);

    mat_t model2view;
    mat_t local2view;
    if (!bv_model2view_get(model2view, view))
	return false;
    bn_mat_mul(local2view, model2view, es->model_changes);

    int best_face = -1;
    int best_u = -1;
    int best_v = -1;
    fastf_t best_dist_sq = max_dist_sq;

    for (int fi = 0; fi < bip->brep->m_F.Count(); ++fi) {
	const ON_BrepFace &face = bip->brep->m_F[fi];
	const ON_NurbsSurface *surface =
	    dynamic_cast<const ON_NurbsSurface *>(face.SurfaceOf());
	if (!surface)
	    continue;

	for (int i = 0; i < surface->CVCount(0); ++i) {
	    for (int j = 0; j < surface->CVCount(1); ++j) {
		ON_3dPoint cv;
		if (!surface->GetCV(i, j, cv))
		    continue;
		point_t local = VINIT_ZERO;
		point_t projected = VINIT_ZERO;
		VSET(local, cv.x, cv.y, cv.z);
		MAT4X3PNT(projected, local2view, local);
		const fastf_t dx = projected[X] - cursor_x;
		const fastf_t dy = projected[Y] - cursor_y;
		const fastf_t dist_sq = dx * dx + dy * dy;
		if (dist_sq < best_dist_sq) {
		    best_dist_sq = dist_sq;
		    best_face = fi;
		    best_u = i;
		    best_v = j;
		}
	    }
	}
    }

    if (best_face < 0)
	return false;

    EDOBJ[es->es_int.idb_type].ft_set_edit_mode(
	    es, ECMD_BREP_SRF_SELECT);
    es->e_inpara = 3;
    es->e_para[0] = (fastf_t)best_face;
    es->e_para[1] = (fastf_t)best_u;
    es->e_para[2] = (fastf_t)best_v;
    rt_edit_process(es);

    struct brep_face_cv_constraint status = {};
    const bool have_status = brep_face_cv_constraint_status(
	    bip->brep, best_face, best_u, best_v, &status);
    const bool safe = have_status && status.topology_safe;
    const bool can_translate = have_status && status.can_translate;
    emit brep_selection_changed(best_face, best_u, best_v, safe,
	    can_translate, have_status ? status.edit_backend
				      : BREP_CV_EDIT_BACKEND_NONE,
	    have_status ? status.constraint_edge_count : 0,
	    have_status ? status.constraint_face_count : 0);
    emit view_updated(QG_VIEW_REFRESH);
    return true;
}

bool
QgBrepFilter::screen_delta_to_local(int from_x, int from_y,
	int to_x, int to_y, vect_t delta) const
{
    const struct bv *view = qg_brep_view(this);
    if (!es || !view)
	return false;

    fastf_t from_vx = 0.0;
    fastf_t from_vy = 0.0;
    fastf_t to_vx = 0.0;
    fastf_t to_vy = 0.0;
    if (!bv_screen_to_view(&from_vx, &from_vy, view,
		(fastf_t)from_x, (fastf_t)from_y)
	    || !bv_screen_to_view(&to_vx, &to_vy, view,
		(fastf_t)to_x, (fastf_t)to_y))
	return false;

    mat_t view2model;
    mat_t model2local;
    if (!bv_view2model_get(view2model, view)
	    || !bn_mat_inverse(model2local, es->model_changes))
	return false;

    vect_t view_delta = VINIT_ZERO;
    vect_t model_delta = VINIT_ZERO;
    VSET(view_delta, to_vx - from_vx, to_vy - from_vy, 0.0);
    MAT4X3VEC(model_delta, view2model, view_delta);
    MAT4X3VEC(delta, model2local, model_delta);
    return true;
}

bool
QgBrepPickCVFilter::eventFilter(QObject *, QEvent *event)
{
    QMouseEvent *mouse = view_sync(event);
    if (!mouse)
	return false;

    if (mouse->type() == QEvent::MouseButtonPress
	    && mouse->buttons().testFlag(Qt::LeftButton)) {
	int sx = 0;
	int sy = 0;
	qg_brep_mouse_xy(mouse, &sx, &sy);
	(void)pick_cv(sx, sy, pick_radius_px);
	return true;
    }

    if (mouse->type() == QEvent::MouseButtonPress
	    || mouse->type() == QEvent::MouseButtonRelease)
	return true;
    return false;
}

bool
QgBrepMoveCVFilter::eventFilter(QObject *, QEvent *event)
{
    QMouseEvent *mouse = view_sync(event);
    if (!mouse || !es)
	return false;

    int sx = 0;
    int sy = 0;
    qg_brep_mouse_xy(mouse, &sx, &sy);

    if (mouse->type() == QEvent::MouseButtonPress
	    && mouse->buttons().testFlag(Qt::LeftButton)) {
	int face = -1;
	int cv_u = -1;
	int cv_v = -1;
	if (pick_on_press)
	    (void)pick_cv(sx, sy, pick_radius_px);
	if (!selected_cv(&face, &cv_u, &cv_v))
	    return true;
	if (!selected_cv_can_translate()) {
	    emit edit_rejected(
		    QStringLiteral("This control vertex does not have a supported "
			"topology-preserving boundary constraint."));
	    return true;
	}
	m_prev_x = sx;
	m_prev_y = sy;
	m_dragging = true;
	return true;
    }

    if (mouse->type() == QEvent::MouseMove && m_dragging
	    && mouse->buttons().testFlag(Qt::LeftButton)) {
	vect_t delta = VINIT_ZERO;
	if (!screen_delta_to_local(m_prev_x, m_prev_y, sx, sy, delta))
	    return true;
	int face = -1;
	int cv_u = -1;
	int cv_v = -1;
	ON_3dPoint before;
	const bool have_before = selected_cv(&face, &cv_u, &cv_v)
	    && qg_brep_selected_point(es, face, cv_u, cv_v, &before);

	EDOBJ[es->es_int.idb_type].ft_set_edit_mode(
		es, ECMD_BREP_SRF_CV_MOVE);
	es->e_inpara = 3;
	es->e_para[0] = delta[X] * es->base2local;
	es->e_para[1] = delta[Y] * es->base2local;
	es->e_para[2] = delta[Z] * es->base2local;
	rt_edit_process(es);
	m_prev_x = sx;
	m_prev_y = sy;
	ON_3dPoint after;
	const double requested_length = MAGNITUDE(delta);
	if (have_before && requested_length > SMALL_FASTF
		&& qg_brep_selected_point(es, face, cv_u, cv_v, &after)
		&& before.DistanceTo(after) <= 1.0e-12
		    * std::max(1.0, requested_length)) {
	    emit edit_rejected(QStringLiteral(
		    "The requested drag failed its C0 constraint or B-rep "
		    "validation checks."));
	    emit view_updated(QG_VIEW_REFRESH);
	    return true;
	}
	emit brep_changed();
	emit view_updated(QG_VIEW_REFRESH);
	return true;
    }

    if (mouse->type() == QEvent::MouseButtonRelease
	    && mouse->button() == Qt::LeftButton) {
	m_dragging = false;
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
