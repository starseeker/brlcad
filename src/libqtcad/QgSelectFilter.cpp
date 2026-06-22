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
#include "rt/view_legacy_bsg.h"
}

#include <algorithm>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "qtcad/QgLegacyView.h"
#include "qtcad/QgLegacyViewBsg.h"
#include "qtcad/QgObolPick.h"
#include "qtcad/QgSelectFilter.h"
#include "qtcad/QgSignalFlags.h"
#include "qtcad/QgView.h"

static qg_legacy_view *
qg_select_filter_legacy_view(const QgSelectFilter *filter)
{
    QgView *display = filter ? filter->view_widget() : nullptr;
    return display ? display->view() : nullptr;
}

static struct bsg_view *
qg_select_filter_view(const QgSelectFilter *filter)
{
    return qg_legacy_view_to_bsg(qg_select_filter_legacy_view(filter));
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
    struct rt_view_pick_result_bsg *selected_result = nullptr;
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
_qg_pick_result_path(const struct rt_view_pick_result_bsg *result, size_t index)
{
    struct bu_vls path = BU_VLS_INIT_ZERO;
    std::string ret;
    if (rt_view_pick_result_path_bsg(result, index, &path))
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
	rt_view_pick_result_free_bsg(m->selected_result);
	m->selected_result = nullptr;
    }
}

void
QgSelectFilter::set_selected_result(struct rt_view_pick_result_bsg *res)
{
    clear_selected_result();
    if (!m)
	return;

    m->selected_result = res;
    struct bsg_view *v = qg_select_filter_view(this);
    if (v)
	rt_view_selection_set_pick_result_ref_bsg(v, m->selected_result,
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

    struct bsg_view *v = qg_select_filter_view(this);
    if (v)
	rt_view_selection_clear_bsg(v);
}

bool
QgSelectPntFilter::eventFilter(QObject *, QEvent *e)
{
    QMouseEvent *m_e = view_sync(e);
    if (!m_e)
	return false;

    struct bsg_view *v = qg_select_filter_view(this);

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
	for (const QgObolPickRecord &pick : obolPicks)
	    _qg_append_unique_path(paths, pick.path.c_str());
	set_selected_paths(paths);
	return true;
    }
    if (submittedSourceRequests > 0) {
	std::vector<std::string> paths;
	set_selected_paths(paths);
	return true;
    }

    struct rt_view_pick_result_bsg *res = first_only ?
	rt_view_pick_nearest_bsg(v, sx, sy) :
	rt_view_pick_point_bsg(v, sx, sy, 0);
    set_selected_result(res);

    return true;
}

bool
QgSelectBoxFilter::eventFilter(QObject *, QEvent *e)
{
    QMouseEvent *m_e = view_sync(e);
    if (!m_e)
	return false;

    struct bsg_view *v = qg_select_filter_view(this);
    if (!v)
	return false;

    if (e->type() == QEvent::MouseButtonDblClick)
	return true;
    if (m_e->button() != Qt::LeftButton && e->type() != QEvent::MouseMove)
	return true;

    int sx = 0, sy = 0;
    qg_select_mouse_xy(m_e, &sx, &sy);

    if (e->type() == QEvent::MouseButtonPress) {
	px = sx;
	py = sy;
	qg_legacy_view *lv = qg_select_filter_legacy_view(this);
	int view_width = qg_legacy_view_width_get(lv);
	int view_height = qg_legacy_view_height_get(lv);
	struct rt_view_interactive_rect_state rect;
	if (!rt_view_interactive_rect_state_from_bsg(&rect, v))
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
	rt_view_interactive_rect_state_set_bsg(v, &rect);
	emit view_updated(QG_VIEW_DRAWN);
	return true;
    }

    if (e->type() == QEvent::MouseMove) {
	struct rt_view_interactive_rect_state rect;
	if (!rt_view_interactive_rect_state_from_bsg(&rect, v))
	    return true;
	qg_legacy_view *lv = qg_select_filter_legacy_view(this);
	int view_height = qg_legacy_view_height_get(lv);
	rect.draw = 1;
	rect.dim[0] = sx - px;
	rect.dim[1] = (view_height - sy) - rect.pos[1];
	rect.x = (rect.pos[X] / (fastf_t)rect.cdim[X] - 0.5) * 2.0;
	rect.y = ((0.5 - (rect.cdim[Y] - rect.pos[Y]) / (fastf_t)rect.cdim[Y]) / rect.aspect * 2.0);
	rect.width = rect.dim[X] * 2.0 / (fastf_t)rect.cdim[X];
	rect.height = rect.dim[Y] * 2.0 / (fastf_t)rect.cdim[X];
	rt_view_interactive_rect_state_set_bsg(v, &rect);
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
	    for (const QgObolPickRecord &pick : obolPicks)
		_qg_append_unique_path(paths, pick.path.c_str());
	    if (!paths.empty()) {
		set_selected_paths(paths);

		struct rt_view_interactive_rect_state rect;
		if (!rt_view_interactive_rect_state_from_bsg(&rect, v))
		    return true;
		rect.draw = 0;
		rect.line_width = 0;
		rect.pos[0] = 0;
		rect.pos[1] = 0;
		rect.dim[0] = 0;
		rect.dim[1] = 0;
		rt_view_interactive_rect_state_set_bsg(v, &rect);
		emit view_updated(QG_VIEW_DRAWN);
		return true;
	    }
	}

	if (submittedSourceRequests > 0) {
	    std::vector<std::string> paths;
	    set_selected_paths(paths);
	    struct rt_view_interactive_rect_state rect;
	    if (!rt_view_interactive_rect_state_from_bsg(&rect, v))
		return true;
	    rect.draw = 0;
	    rect.line_width = 0;
	    rect.pos[0] = 0;
	    rect.pos[1] = 0;
	    rect.dim[0] = 0;
	    rect.dim[1] = 0;
	    rt_view_interactive_rect_state_set_bsg(v, &rect);
	    emit view_updated(QG_VIEW_DRAWN);
	    return true;
	}

	struct rt_view_pick_result_bsg *res =
	    rt_view_pick_rect_bsg(v, ipx, ipy, sx, sy);
	if (first_only && res && rt_view_pick_result_count_bsg(res) > 1) {
	    struct rt_view_pick_result_bsg *nearest =
		rt_view_pick_result_filter_first_bsg(res);
	    rt_view_pick_result_free_bsg(res);
	    res = nearest;
	}
	set_selected_result(res);

	struct rt_view_interactive_rect_state rect;
	if (!rt_view_interactive_rect_state_from_bsg(&rect, v))
	    return true;
	rect.draw = 0;
	rect.line_width = 0;
	rect.pos[0] = 0;
	rect.pos[1] = 0;
	rect.dim[0] = 0;
	rect.dim[1] = 0;
	rt_view_interactive_rect_state_set_bsg(v, &rect);
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

static struct rt_view_pick_result_bsg *
_qg_pick_result_from_ray_hits(const struct rt_view_pick_result_bsg *candidates,
			      const struct select_rec_state *rc,
			      int first_only)
{
    struct rt_view_pick_result_bsg *res = rt_view_pick_result_create_bsg();
    if (!res || !candidates || !rc)
	return res;

    std::unordered_set<std::string> seen_paths;
    for (size_t i = 0; i < rt_view_pick_result_count_bsg(candidates); i++) {
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

	fastf_t hit_dist = rt_view_pick_result_hit_dist_bsg(candidates, i);
	std::unordered_map<std::string, fastf_t>::const_iterator h_it = rc->hits.find(key);
	if (h_it != rc->hits.end())
	    hit_dist = h_it->second;

	rt_view_pick_result_append_copy_bsg(res, candidates, i, hit_dist);

	if (first_only)
	    break;
    }

    return res;
}

static bool
_qg_select_ray_from_view(struct bsg_view *v, int sx, int sy, point_t origin, vect_t direction)
{
    if (!v)
	return false;

    fastf_t vx = -FLT_MAX;
    fastf_t vy = -FLT_MAX;
    if (!rt_view_screen_to_view_from_bsg(&vx, &vy, v, sx, sy))
	return false;
    point_t vpnt, mpnt;
    VSET(vpnt, vx, vy, 0);
    mat_t view2model;
    rt_view_view2model_from_bsg(view2model, v);
    MAT4X3PNT(mpnt, view2model, vpnt);
    mat_t view_rotation;
    rt_view_rotation_from_bsg(view_rotation, v);
    VMOVEN(direction, view_rotation + 8, 3);
    VUNITIZE(direction);
    VSCALE(direction, direction, rt_view_radius_from_bsg(v));
    VADD2(origin, mpnt, direction);
    VUNITIZE(direction);
    VSCALE(direction, direction, -1);
    return true;
}

bool
QgSelectRayFilter::eventFilter(QObject *, QEvent *e)
{
    QMouseEvent *m_e = view_sync(e);
    if (!m_e)
	return false;

    struct bsg_view *v = qg_select_filter_view(this);
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
	    for (const QgObolPickRecord &pick : obolRayPicks)
		_qg_append_unique_path(paths, pick.path.c_str());
	    if (!paths.empty()) {
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
	for (const QgObolPickRecord &pick : obolPicks)
	    _qg_append_unique_path(paths, pick.path.c_str());
	if (!paths.empty()) {
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

    struct rt_view_pick_result_bsg *candidates =
	rt_view_pick_point_bsg(v, sx, sy, 0);
    size_t candidate_count = rt_view_pick_result_count_bsg(candidates);
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
	rt_view_pick_result_free_bsg(candidates);
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

    struct rt_view_pick_result_bsg *res =
	_qg_pick_result_from_ray_hits(candidates, &rc, first_only);
    rt_view_pick_result_free_bsg(candidates);
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
