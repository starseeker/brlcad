/*                 Q G S E L E C T F I L T E R . C P P
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
/** @file QgSelectFilter.cpp
 *
 * Graphical selection tool for Qt views.
 *
 */

#include "common.h"

extern "C" {
#include "bu/malloc.h"
#include "bu/vls.h"
#include "raytrace.h"
}

#include <algorithm>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "qtcad/QgObolPick.h"
#include "qtcad/QgSelectFilter.h"
#include "qtcad/QgSignalFlags.h"
#include "qtcad/QgView.h"
#include "rt/view.h"
#include "QgLegacyViewContext.h"

static qg_legacy_view *
qg_select_filter_legacy_view(const QgSelectFilter *filter)
{
    QgView *display = filter ? filter->view_widget() : nullptr;
    return display ? display->view() : nullptr;
}

static void
qg_select_mouse_xy(QMouseEvent *m_e, int *sx, int *sy)
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

class QgSelectFilter::QgSelectFilterPrivate {
public:
    void *selected_result = nullptr;
    std::vector<std::string> selected_path_strings;
};

static std::string
_qg_normalize_path(const char *path)
{
    if (!path)
	return std::string();
    while (*path == '/')
	path++;
    return std::string(path);
}

static void
_qg_append_unique_path(std::vector<std::string> &paths, const char *path)
{
    std::string normalized = _qg_normalize_path(path);
    if (normalized.empty())
	return;
    if (std::find(paths.begin(), paths.end(), normalized) == paths.end())
	paths.push_back(normalized);
}

static std::string
_qg_pick_result_path(const void *result, size_t index)
{
    struct bu_vls path = BU_VLS_INIT_ZERO;
    std::string ret;
    if (rt_view_pick_result_context_path(result, index, &path))
	ret = bu_vls_cstr(&path);
    bu_vls_free(&path);
    return ret;
}

static void
_qg_append_unique_path_cb(const char *path, void *data)
{
    std::vector<std::string> *paths =
	static_cast<std::vector<std::string> *>(data);
    if (paths)
	_qg_append_unique_path(*paths, path);
}

static int
_qg_apply_obol_pick_selection(QgView *display,
	const std::vector<QgObolPickRecord> &picks,
	std::vector<std::string> &database_paths)
{
    std::vector<QgObolPickRecord> feature_picks;
    for (const QgObolPickRecord &pick : picks) {
	if (pick.featurePickResolved) {
	    feature_picks.push_back(pick);
	    continue;
	}
	_qg_append_unique_path(database_paths, pick.path.c_str());
    }

    return qg_obol_pick_apply_feature_states(display, feature_picks, true,
	    true);
}

QgSelectFilter::QgSelectFilter() :
    m(new QgSelectFilterPrivate)
{
}

QgSelectFilter::~QgSelectFilter()
{
    clear_selected_result();
    delete m;
    m = nullptr;
}

const std::vector<std::string> &
QgSelectFilter::selected_paths() const
{
    static const std::vector<std::string> empty_paths;
    return m ? m->selected_path_strings : empty_paths;
}

void
QgSelectFilter::clear_selected_result()
{
    if (!m)
	return;

    m->selected_path_strings.clear();
    if (m->selected_result) {
	rt_view_pick_result_context_free(m->selected_result);
	m->selected_result = nullptr;
    }
}

void
QgSelectFilter::set_selected_result(void *res)
{
    clear_selected_result();
    if (!m)
	return;

    m->selected_result = res;
    qg_legacy_view *v = qg_select_filter_legacy_view(this);
    if (v)
	rt_view_context_selection_set_pick_result_context(
		qg_legacy_view_to_context(v), m->selected_result,
		_qg_append_unique_path_cb, &m->selected_path_strings);
}

void
QgSelectFilter::set_selected_paths(const std::vector<std::string> &paths)
{
    clear_selected_result();
    if (!m)
	return;

    for (const std::string &path : paths)
	_qg_append_unique_path(m->selected_path_strings, path.c_str());

    qg_legacy_view *v = qg_select_filter_legacy_view(this);
    if (v)
	rt_view_context_selection_clear(qg_legacy_view_to_context(v));
}

bool
QgSelectPntFilter::eventFilter(QObject *, QEvent *e)
{
    QMouseEvent *m_e = view_sync(e);
    if (!m_e)
	return false;

    qg_legacy_view *v = qg_select_filter_legacy_view(this);

    if (e->type() != QEvent::MouseButtonRelease)
	return true;
    if (m_e->button() != Qt::LeftButton)
	return true;
    if (!v)
	return true;

    int sx = 0, sy = 0;
    qg_select_mouse_xy(m_e, &sx, &sy);

    std::vector<QgObolPickRecord> obolPicks;
    int submittedSourceRequests = 0;
    if (qg_obol_pick_point(view_widget(), sx, sy,
	    6.0f, !first_only, obolPicks,
	    &submittedSourceRequests) > 0) {
	std::vector<std::string> paths;
	int feature_count = _qg_apply_obol_pick_selection(view_widget(),
		obolPicks, paths);
	if (!paths.empty() || feature_count > 0) {
	    set_selected_paths(paths);
	    return true;
	}
    }
    if (submittedSourceRequests > 0) {
	std::vector<std::string> paths;
	set_selected_paths(paths);
	return true;
    }

    void *res = first_only ?
	rt_view_context_pick_nearest(qg_legacy_view_to_context(v), sx, sy) :
	rt_view_context_pick_point(qg_legacy_view_to_context(v), sx, sy, 0);
    set_selected_result(res);

    return true;
}

bool
QgSelectBoxFilter::eventFilter(QObject *, QEvent *e)
{
    QMouseEvent *m_e = view_sync(e);
    if (!m_e)
	return false;

    qg_legacy_view *v = qg_select_filter_legacy_view(this);
    if (!v)
	return false;
    void *view_ctx = qg_legacy_view_to_context(v);

    if (e->type() == QEvent::MouseButtonDblClick)
	return true;
    if (m_e->button() != Qt::LeftButton && e->type() != QEvent::MouseMove)
	return true;

    int sx = 0, sy = 0;
    qg_select_mouse_xy(m_e, &sx, &sy);

    if (e->type() == QEvent::MouseButtonPress) {
	px = sx;
	py = sy;
	int view_width = rt_view_context_width_get(view_ctx);
	int view_height = rt_view_context_height_get(view_ctx);
	struct rt_view_interactive_rect_state rect;
	if (!rt_view_context_interactive_rect_state_get(&rect, view_ctx))
	    return true;
	rect.line_width = 1;
	rect.dim[0] = 0;
	rect.dim[1] = 0;
	rect.x = px;
	rect.y = view_height - py;
	rect.pos[0] = rect.x;
	rect.pos[1] = rect.y;
	rect.cdim[0] = view_width;
	rect.cdim[1] = view_height;
	rect.aspect = (fastf_t)rect.cdim[X] / rect.cdim[Y];
	rt_view_context_interactive_rect_state_set(view_ctx, &rect);
	emit view_updated(QG_VIEW_DRAWN);
	return true;
    }

    if (e->type() == QEvent::MouseMove) {
	struct rt_view_interactive_rect_state rect;
	if (!rt_view_context_interactive_rect_state_get(&rect, view_ctx))
	    return true;
	int view_height = rt_view_context_height_get(view_ctx);
	rect.draw = 1;
	rect.dim[0] = sx - px;
	rect.dim[1] = (view_height - sy) - rect.pos[1];
	rect.x = (rect.pos[X] / (fastf_t)rect.cdim[X] - 0.5) * 2.0;
	rect.y = ((0.5 - (rect.cdim[Y] - rect.pos[Y]) / (fastf_t)rect.cdim[Y]) / rect.aspect * 2.0);
	rect.width = rect.dim[X] * 2.0 / (fastf_t)rect.cdim[X];
	rect.height = rect.dim[Y] * 2.0 / (fastf_t)rect.cdim[X];
	rt_view_context_interactive_rect_state_set(view_ctx, &rect);
	emit view_updated(QG_VIEW_DRAWN);
	return true;
    }

    if (e->type() == QEvent::MouseButtonRelease) {
	int ipx = (int)px;
	int ipy = (int)py;
	std::vector<QgObolPickRecord> obolPicks;
	int submittedSourceRequests = 0;
	if (qg_obol_pick_rect(view_widget(), ipx, ipy,
		sx, sy, 6.0f, first_only,
		obolPicks, &submittedSourceRequests) > 0) {
	    std::vector<std::string> paths;
	    int feature_count = _qg_apply_obol_pick_selection(view_widget(),
		    obolPicks, paths);
	    if (!paths.empty() || feature_count > 0) {
		set_selected_paths(paths);

		struct rt_view_interactive_rect_state rect;
		if (!rt_view_context_interactive_rect_state_get(&rect, view_ctx))
		    return true;
		rect.draw = 0;
		rect.line_width = 0;
		rect.pos[0] = 0;
		rect.pos[1] = 0;
		rect.dim[0] = 0;
		rect.dim[1] = 0;
		rt_view_context_interactive_rect_state_set(view_ctx, &rect);
		emit view_updated(QG_VIEW_DRAWN);
		return true;
	    }
	}

	if (submittedSourceRequests > 0) {
	    std::vector<std::string> paths;
	    set_selected_paths(paths);
	    struct rt_view_interactive_rect_state rect;
	    if (!rt_view_context_interactive_rect_state_get(&rect, view_ctx))
		return true;
	    rect.draw = 0;
	    rect.line_width = 0;
	    rect.pos[0] = 0;
	    rect.pos[1] = 0;
	    rect.dim[0] = 0;
	    rect.dim[1] = 0;
	    rt_view_context_interactive_rect_state_set(view_ctx, &rect);
	    emit view_updated(QG_VIEW_DRAWN);
	    return true;
	}

	void *res = rt_view_context_pick_rect(view_ctx, ipx, ipy, sx, sy);
	if (first_only && res && rt_view_pick_result_context_count(res) > 1) {
	    void *nearest = rt_view_pick_result_context_filter_first(res);
	    rt_view_pick_result_context_free(res);
	    res = nearest;
	}
	set_selected_result(res);

	struct rt_view_interactive_rect_state rect;
	if (!rt_view_context_interactive_rect_state_get(&rect, view_ctx))
	    return true;
	rect.draw = 0;
	rect.line_width = 0;
	rect.pos[0] = 0;
	rect.pos[1] = 0;
	rect.dim[0] = 0;
	rect.dim[1] = 0;
	rt_view_context_interactive_rect_state_set(view_ctx, &rect);
	emit view_updated(QG_VIEW_DRAWN);
	return true;
    }

    return false;
}

struct select_rec_state {
    std::unordered_map<std::string, fastf_t> hits;
    int rec_all;
    double cdist;
    std::string closest;
};

static void
_select_record_hit(struct select_rec_state *rc, const char *name, fastf_t hit_dist)
{
    if (!rc || !name || !name[0])
	return;

    std::string key = _qg_normalize_path(name);
    if (key.empty())
	return;

    std::unordered_map<std::string, fastf_t>::iterator h_it = rc->hits.find(key);
    if (h_it == rc->hits.end() || hit_dist < h_it->second)
	rc->hits[key] = hit_dist;

    if (hit_dist < rc->cdist) {
	rc->closest = key;
	rc->cdist = hit_dist;
    }
}

static int
_obj_record(struct application *ap, struct partition *p_hp, struct seg *UNUSED(segs))
{
    struct select_rec_state *rc = (struct select_rec_state *)ap->a_uptr;
    for (struct partition *pp = p_hp->pt_forw; pp != p_hp; pp = pp->pt_forw) {
	if (rc->rec_all) {
	    _select_record_hit(rc, pp->pt_regionp->reg_name, pp->pt_inhit->hit_dist);
	} else {
	    _select_record_hit(rc, pp->pt_regionp->reg_name, pp->pt_inhit->hit_dist);
	}
    }
    return 1;
}

static int
_ovlp_record(struct application *ap, struct partition *pp, struct region *reg1, struct region *reg2, struct partition *UNUSED(ihp))
{
    struct select_rec_state *rc = (struct select_rec_state *)ap->a_uptr;
    if (rc->rec_all) {
	_select_record_hit(rc, reg1->reg_name, pp->pt_inhit->hit_dist);
	_select_record_hit(rc, reg2->reg_name, pp->pt_inhit->hit_dist);
    } else {
	_select_record_hit(rc, reg1->reg_name, pp->pt_inhit->hit_dist);
    }
    return 1;
}

static void *
_qg_pick_result_from_ray_hits(const void *candidates,
			      const struct select_rec_state *rc,
			      int first_only)
{
    void *res = rt_view_pick_result_context_create();
    if (!res || !candidates || !rc)
	return res;

    std::unordered_set<std::string> seen_paths;
    for (size_t i = 0; i < rt_view_pick_result_context_count(candidates); i++) {
	std::string key = _qg_normalize_path(_qg_pick_result_path(candidates, i).c_str());
	if (key.empty())
	    continue;

	if (first_only) {
	    if (key != rc->closest)
		continue;
	} else {
	    if (rc->hits.find(key) == rc->hits.end())
		continue;
	}

	if (!seen_paths.insert(key).second)
	    continue;

	fastf_t hit_dist = rt_view_pick_result_context_hit_dist(candidates, i);
	std::unordered_map<std::string, fastf_t>::const_iterator h_it = rc->hits.find(key);
	if (h_it != rc->hits.end())
	    hit_dist = h_it->second;

	rt_view_pick_result_context_append_copy(res, candidates, i, hit_dist);

	if (first_only)
	    break;
    }

    return res;
}

static bool
_qg_select_ray_from_view(qg_legacy_view *v, int sx, int sy, point_t origin, vect_t direction)
{
    void *view_ctx = qg_legacy_view_to_context(v);
    if (!view_ctx || !origin || !direction)
	return false;

    fastf_t vx = 0.0;
    fastf_t vy = 0.0;
    if (!rt_view_context_screen_to_view(&vx, &vy, view_ctx, sx, sy))
	return false;

    point_t vpnt;
    point_t mpnt;
    VSET(vpnt, vx, vy, 0.0);

    mat_t view2model;
    if (!rt_view_context_view2model_get(view2model, view_ctx))
	return false;
    MAT4X3PNT(mpnt, view2model, vpnt);

    mat_t view_rotation;
    if (!rt_view_context_rotation_get(view_rotation, view_ctx))
	return false;
    VMOVEN(direction, view_rotation + 8, 3);
    VUNITIZE(direction);
    VSCALE(direction, direction, rt_view_context_radius_get(view_ctx));
    VADD2(origin, mpnt, direction);
    VUNITIZE(direction);
    VSCALE(direction, direction, -1.0);
    return true;
}

bool
QgSelectRayFilter::eventFilter(QObject *, QEvent *e)
{
    QMouseEvent *m_e = view_sync(e);
    if (!m_e)
	return false;

    qg_legacy_view *v = qg_select_filter_legacy_view(this);
    if (!v)
	return false;
    if (e->type() != QEvent::MouseButtonRelease)
	return true;
    if (m_e->button() != Qt::LeftButton)
	return true;

    int sx = 0, sy = 0;
    qg_select_mouse_xy(m_e, &sx, &sy);

    point_t rayOrigin;
    vect_t rayDirection;
    if (_qg_select_ray_from_view(v, sx, sy, rayOrigin, rayDirection)) {
	std::vector<QgObolPickRecord> obolRayPicks;
	int submittedSourceRequests = 0;
	if (qg_obol_pick_ray(view_widget(),
		SbVec3f(rayOrigin[0], rayOrigin[1], rayOrigin[2]),
		SbVec3f(rayDirection[0], rayDirection[1], rayDirection[2]),
		!first_only, obolRayPicks,
		&submittedSourceRequests) > 0) {
	    std::vector<std::string> paths;
	    int feature_count = _qg_apply_obol_pick_selection(view_widget(),
		    obolRayPicks, paths);
	    if (!paths.empty() || feature_count > 0) {
		set_selected_paths(paths);
		return true;
	    }
	}
	if (submittedSourceRequests > 0) {
	    std::vector<std::string> paths;
	    set_selected_paths(paths);
	    return true;
	}
    }

    std::vector<QgObolPickRecord> obolPicks;
    int submittedSourceRequests = 0;
    if (qg_obol_pick_point(view_widget(), sx, sy,
	    6.0f, !first_only, obolPicks,
	    &submittedSourceRequests) > 0) {
	std::vector<std::string> paths;
	int feature_count = _qg_apply_obol_pick_selection(view_widget(),
		obolPicks, paths);
	if (!paths.empty() || feature_count > 0) {
	    set_selected_paths(paths);
	    return true;
	}
    }
    if (submittedSourceRequests > 0) {
	std::vector<std::string> paths;
	set_selected_paths(paths);
	return true;
    }

    if (!dbip) {
	std::vector<std::string> paths;
	set_selected_paths(paths);
	return true;
    }

    void *candidates = rt_view_context_pick_point(qg_legacy_view_to_context(v),
	    sx, sy, 0);
    size_t candidate_count = candidates ?
	rt_view_pick_result_context_count(candidates) : 0;
    if (!candidates || !candidate_count) {
	set_selected_result(candidates);
	return true;
    }

    struct application *ap;
    BU_GET(ap, struct application);
    RT_APPLICATION_INIT(ap);
    ap->a_onehit = 0;
    ap->a_hit = _obj_record;
    ap->a_miss = nullptr;
    ap->a_overlap = _ovlp_record;
    ap->a_logoverlap = nullptr;

    struct rt_i *rtip = rt_i_create(dbip);
    struct resource *resp = nullptr;
    BU_GET(resp, struct resource);
    rt_init_resource(resp, 0, rtip);
    ap->a_resource = resp;
    ap->a_rt_i = rtip;
    std::vector<std::string> candidate_paths;
    candidate_paths.reserve(candidate_count);
    const char **objs = (const char **)bu_calloc(candidate_count + 1, sizeof(char *), "objs");
    for (size_t i = 0; i < candidate_count; i++) {
	candidate_paths.push_back(_qg_pick_result_path(candidates, i));
	objs[i] = candidate_paths.back().empty() ? nullptr : candidate_paths.back().c_str();
    }
    if (rt_gettrees_and_attrs(rtip, nullptr, (int)candidate_count, objs, 1)) {
	bu_free(objs, "objs");
	rt_i_destroy(rtip);
	BU_PUT(resp, struct resource);
	BU_PUT(ap, struct application);
	rt_view_pick_result_context_free(candidates);
	return false;
    }
    size_t ncpus = bu_avail_cpus();
    rt_prep_parallel(rtip, (int)ncpus);

    struct select_rec_state rc;
    rc.cdist = INFINITY;
    if (!first_only) {
	rc.rec_all = 1;
    } else {
	rc.rec_all = 0;
    }
    ap->a_uptr = (void *)&rc;

    if (_qg_select_ray_from_view(v, sx, sy, rayOrigin, rayDirection)) {
	VMOVE(ap->a_ray.r_pt, rayOrigin);
	VMOVE(ap->a_ray.r_dir, rayDirection);
	(void)rt_shootray(ap);
    }
    bu_free(objs, "objs");
    rt_i_destroy(rtip);
    BU_PUT(resp, struct resource);
    BU_PUT(ap, struct application);

    void *res =
	_qg_pick_result_from_ray_hits(candidates, &rc, first_only);
    rt_view_pick_result_context_free(candidates);
    set_selected_result(res);

    return true;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
