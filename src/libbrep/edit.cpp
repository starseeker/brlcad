/*                        E D I T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2023-2026 United States Government as represented by
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
/** @file edit.cpp
 *
 * Implementation of edit support for brep.
 *
 */

#include "brep/edit.h"
#include "bu/log.h"
#include "bu/malloc.h"
#include <Eigen/Core>
#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <functional>
#include <set>
#include <utility>
#include <vector>

namespace {

static bool
valid_brep(const ON_Brep *brep)
{
    if (!brep) {
	bu_log("brep is NULL\n");
	return false;
    }
    return true;
}

static bool
curve3d_referenced(const ON_Brep *brep, int curve_id)
{
    return brep && brep->EdgeCurveUseCount(curve_id, 1) > 0;
}

static bool
surface_referenced(const ON_Brep *brep, int surface_id)
{
    return brep && brep->SurfaceUseCount(surface_id, 1) > 0;
}

static bool
replace_curve3d(ON_Brep *brep, int curve_id, ON_Curve *replacement)
{
    if (!brep || !replacement || curve_id < 0 || curve_id >= brep->m_C3.Count())
	return false;

    ON_Curve *old_curve = brep->m_C3[curve_id];
    brep->m_C3[curve_id] = replacement;
    for (int ei = 0; ei < brep->m_E.Count(); ++ei) {
	ON_BrepEdge &edge = brep->m_E[ei];
	if (edge.m_c3i == curve_id) {
	    ON_Interval proxy_domain = edge.ProxyCurveDomain();
	    ON_Interval edge_domain = edge.Domain();
	    const bool proxy_reversed = edge.ProxyCurveIsReversed();
	    edge.SetProxyCurve(replacement, proxy_domain);
	    if (proxy_reversed != edge.ProxyCurveIsReversed())
		edge.ON_CurveProxy::Reverse();
	    edge.SetDomain(edge_domain);
	}
    }
    delete old_curve;
    brep->DestroyRuntimeCache(true);
    return true;
}

static bool
replace_surface(ON_Brep *brep, int surface_id, ON_Surface *replacement)
{
    if (!brep || !replacement || surface_id < 0 || surface_id >= brep->m_S.Count())
	return false;

    ON_Surface *old_surface = brep->m_S[surface_id];
    brep->m_S[surface_id] = replacement;
    for (int fi = 0; fi < brep->m_F.Count(); ++fi) {
	ON_BrepFace &face = brep->m_F[fi];
	if (face.m_si == surface_id) {
	    const bool proxy_transposed = face.ProxySurfaceIsTransposed();
	    face.DestroyRuntimeCache(true);
	    face.SetProxySurface(replacement);
	    if (proxy_transposed != face.ProxySurfaceIsTransposed())
		face.ON_SurfaceProxy::Transpose();
	}
    }
    delete old_surface;
    brep->DestroyRuntimeCache(true);
    return true;
}

static bool
edit_nurbs_curve(ON_Brep *brep, int curve_id, bool allow_referenced,
	const std::function<bool(ON_NurbsCurve *)> &operation)
{
    if (!valid_brep(brep) || curve_id < 0 || curve_id >= brep->m_C3.Count()) {
	bu_log("curve_id is out of range\n");
	return false;
    }
    if (!allow_referenced && curve3d_referenced(brep, curve_id)) {
	bu_log("curve %d is referenced by a brep edge\n", curve_id);
	return false;
    }

    const ON_NurbsCurve *curve = dynamic_cast<const ON_NurbsCurve *>(brep->m_C3[curve_id]);
    if (!curve) {
	bu_log("curve %d is not a NURBS curve\n", curve_id);
	return false;
    }

    ON_Curve *curve_copy = curve->DuplicateCurve();
    ON_NurbsCurve *nurbs_copy = dynamic_cast<ON_NurbsCurve *>(curve_copy);
    if (!nurbs_copy || !operation(nurbs_copy)) {
	delete curve_copy;
	return false;
    }
    return replace_curve3d(brep, curve_id, curve_copy);
}

static bool
edit_nurbs_surface(ON_Brep *brep, int surface_id, bool allow_referenced,
	const std::function<bool(ON_NurbsSurface *)> &operation)
{
    if (!valid_brep(brep) || surface_id < 0 || surface_id >= brep->m_S.Count()) {
	bu_log("surface_id is out of range\n");
	return false;
    }
    if (!allow_referenced && surface_referenced(brep, surface_id)) {
	bu_log("surface %d is referenced by a brep face\n", surface_id);
	return false;
    }

    const ON_NurbsSurface *surface = dynamic_cast<const ON_NurbsSurface *>(brep->m_S[surface_id]);
    if (!surface) {
	bu_log("surface %d is not a NURBS surface\n", surface_id);
	return false;
    }

    ON_Surface *surface_copy = surface->DuplicateSurface();
    ON_NurbsSurface *nurbs_copy = dynamic_cast<ON_NurbsSurface *>(surface_copy);
    if (!nurbs_copy || !operation(nurbs_copy)) {
	delete surface_copy;
	return false;
    }
    return replace_surface(brep, surface_id, surface_copy);
}

static bool
edit_face_nurbs_surface(ON_Brep *brep, int face_id,
	const std::function<bool(ON_NurbsSurface *)> &operation)
{
    if (!valid_brep(brep) || face_id < 0 || face_id >= brep->m_F.Count()) {
	bu_log("face is out of range\n");
	return false;
    }

    ON_BrepFace &face = brep->m_F[face_id];
    const int surface_id = face.m_si;
    if (surface_id < 0 || surface_id >= brep->m_S.Count())
	return false;

    const ON_NurbsSurface *surface = dynamic_cast<const ON_NurbsSurface *>(brep->m_S[surface_id]);
    if (!surface) {
	bu_log("face %d does not reference a NURBS surface\n", face_id);
	return false;
    }

    ON_Surface *surface_copy = surface->DuplicateSurface();
    ON_NurbsSurface *nurbs_copy = dynamic_cast<ON_NurbsSurface *>(surface_copy);
    if (!nurbs_copy || !operation(nurbs_copy)) {
	delete surface_copy;
	return false;
    }

    if (brep->SurfaceUseCount(surface_id, 2) > 1) {
	const int new_surface_id = brep->AddSurface(surface_copy);
	if (new_surface_id < 0) {
	    delete surface_copy;
	    return false;
	}
	face.m_si = new_surface_id;
	face.DestroyRuntimeCache(true);
	face.SetProxySurface(surface_copy);
	brep->DestroyRuntimeCache(true);
	return true;
    }

    return replace_surface(brep, surface_id, surface_copy);
}

static bool
nurbs_surface_cv_support(const ON_NurbsSurface *surface, int dir, int cv_id,
	ON_Interval *support)
{
    if (!surface || !support || dir < 0 || dir > 1
	    || cv_id < 0 || cv_id >= surface->CVCount(dir))
	return false;

    const int order = surface->Order(dir);
    const int knot_count = surface->KnotCount(dir);
    if (order < 2 || knot_count < 2)
	return false;

    /*
     * openNURBS stores the standard knot vector without its first and last
     * repeated entries.  Convert the full-vector support indices
     * [cv_id, cv_id + order] to openNURBS knot indices.
     */
    const int full_left = cv_id;
    const int full_right = cv_id + order;
    const double left = (full_left <= 0)
	? surface->Knot(dir, 0)
	: surface->Knot(dir, full_left - 1);
    const double right = (full_right >= knot_count + 1)
	? surface->Knot(dir, knot_count - 1)
	: surface->Knot(dir, full_right - 1);

    support->Set(left, right);
    return support->IsIncreasing();
}

static bool
support_intersects_interval(double interval_min, double interval_max,
	const ON_Interval &support, const ON_Interval &domain,
	bool clamped_start_cv, bool clamped_end_cv)
{
    /*
     * A basis function is zero at the open ends of its support.  The
     * first/last basis of a clamped vector is the exception at the
     * corresponding surface-domain end.
     */
    if (interval_max > support.Min() + ON_ZERO_TOLERANCE
	    && interval_min < support.Max() - ON_ZERO_TOLERANCE)
	return true;
    if (clamped_start_cv
	    && interval_min <= domain.Min() + ON_ZERO_TOLERANCE
	    && interval_max >= domain.Min() - ON_ZERO_TOLERANCE)
	return true;
    if (clamped_end_cv
	    && interval_min <= domain.Max() + ON_ZERO_TOLERANCE
	    && interval_max >= domain.Max() - ON_ZERO_TOLERANCE)
	return true;
    return false;
}

static bool
face_cv_influenced_trims(const ON_Brep *brep, int face_id, int cv_id_u,
	int cv_id_v, std::vector<int> *trim_ids)
{
    if (!brep || !trim_ids || face_id < 0 || face_id >= brep->m_F.Count())
	return false;
    trim_ids->clear();

    const ON_BrepFace &face = brep->m_F[face_id];
    if (face.ProxySurfaceIsTransposed())
	return false;
    const ON_NurbsSurface *surface =
	dynamic_cast<const ON_NurbsSurface *>(face.SurfaceOf());
    if (!surface)
	return false;

    ON_Interval u_support;
    ON_Interval v_support;
    if (!nurbs_surface_cv_support(surface, 0, cv_id_u, &u_support)
	    || !nurbs_surface_cv_support(surface, 1, cv_id_v, &v_support))
	return false;

    const ON_Interval u_domain = surface->Domain(0);
    const ON_Interval v_domain = surface->Domain(1);
    const bool u_clamped_start_cv =
	cv_id_u == 0 && surface->IsClamped(0, 0);
    const bool u_clamped_end_cv =
	cv_id_u == surface->CVCount(0) - 1 && surface->IsClamped(0, 1);
    const bool v_clamped_start_cv =
	cv_id_v == 0 && surface->IsClamped(1, 0);
    const bool v_clamped_end_cv =
	cv_id_v == surface->CVCount(1) - 1 && surface->IsClamped(1, 1);

    for (int fli = 0; fli < face.m_li.Count(); ++fli) {
	const int loop_id = face.m_li[fli];
	if (loop_id < 0 || loop_id >= brep->m_L.Count())
	    return false;
	const ON_BrepLoop &loop = brep->m_L[loop_id];
	for (int lti = 0; lti < loop.m_ti.Count(); ++lti) {
	    const int trim_id = loop.m_ti[lti];
	    if (trim_id < 0 || trim_id >= brep->m_T.Count())
		return false;
	    const ON_BrepTrim &trim = brep->m_T[trim_id];
	    if (trim.m_c2i < 0 || trim.m_c2i >= brep->m_C2.Count()
		    || !brep->m_C2[trim.m_c2i])
		return false;

	    ON_BoundingBox trim_bbox;
	    if (!brep->m_C2[trim.m_c2i]->GetBoundingBox(trim_bbox, false)
		    || !trim_bbox.IsValid())
		return false;

	    const bool intersects_u = support_intersects_interval(
		    trim_bbox.m_min.x, trim_bbox.m_max.x, u_support, u_domain,
		    u_clamped_start_cv, u_clamped_end_cv);
	    const bool intersects_v = support_intersects_interval(
		    trim_bbox.m_min.y, trim_bbox.m_max.y, v_support, v_domain,
		    v_clamped_start_cv, v_clamped_end_cv);
	    if (intersects_u && intersects_v)
		trim_ids->push_back(trim_id);
	}
    }

    return true;
}

static bool
iso_flag(ON_Surface::ISO iso)
{
    return iso == ON_Surface::x_iso || iso == ON_Surface::y_iso
	|| iso == ON_Surface::W_iso || iso == ON_Surface::S_iso
	|| iso == ON_Surface::E_iso || iso == ON_Surface::N_iso;
}

static bool
natural_iso_flag(ON_Surface::ISO iso)
{
    return iso == ON_Surface::W_iso || iso == ON_Surface::S_iso
	|| iso == ON_Surface::E_iso || iso == ON_Surface::N_iso;
}

static int
trim_face_index(const ON_Brep *brep, int trim_id)
{
    if (!brep || trim_id < 0 || trim_id >= brep->m_T.Count())
	return -1;
    const int loop_id = brep->m_T[trim_id].m_li;
    if (loop_id < 0 || loop_id >= brep->m_L.Count())
	return -1;
    return brep->m_L[loop_id].m_fi;
}

struct iso_trim_trace {
    int trim_id = -1;
    int face_id = -1;
    int edge_id = -1;
    int varying_dir = -1;
    double constant_parameter = ON_UNSET_VALUE;
    double varying_min = ON_UNSET_VALUE;
    double varying_max = ON_UNSET_VALUE;
    bool reverse_to_edge = false;
};

static bool
get_iso_trim_trace(const ON_Brep *brep, int trim_id,
	struct iso_trim_trace *trace)
{
    if (!brep || !trace || trim_id < 0 || trim_id >= brep->m_T.Count())
	return false;

    const ON_BrepTrim &trim = brep->m_T[trim_id];
    if (!iso_flag(trim.m_iso) || trim.m_ei < 0
	    || trim.m_ei >= brep->m_E.Count())
	return false;
    const int face_id = trim_face_index(brep, trim_id);
    if (face_id < 0 || face_id >= brep->m_F.Count())
	return false;
    const ON_BrepFace &face = brep->m_F[face_id];
    if (face.ProxySurfaceIsTransposed())
	return false;
    const ON_NurbsSurface *surface =
	dynamic_cast<const ON_NurbsSurface *>(face.SurfaceOf());
    if (!surface)
	return false;

    const ON_3dPoint start = trim.PointAtStart();
    const ON_3dPoint end = trim.PointAtEnd();
    const bool constant_u = trim.m_iso == ON_Surface::x_iso
	|| trim.m_iso == ON_Surface::W_iso
	|| trim.m_iso == ON_Surface::E_iso;
    const double constant_start = constant_u ? start.x : start.y;
    const double constant_end = constant_u ? end.x : end.y;
    const double varying_start = constant_u ? start.y : start.x;
    const double varying_end = constant_u ? end.y : end.x;
    const ON_Interval constant_domain = surface->Domain(constant_u ? 0 : 1);
    const ON_Interval varying_domain = surface->Domain(constant_u ? 1 : 0);
    const double parameter_tol = 1.0e-10
	* std::max(1.0, std::max(constant_domain.Length(),
		varying_domain.Length()));
    if (std::fabs(constant_start - constant_end) > parameter_tol
	    || constant_start < constant_domain.Min() - parameter_tol
	    || constant_start > constant_domain.Max() + parameter_tol
	    || varying_start < varying_domain.Min() - parameter_tol
	    || varying_start > varying_domain.Max() + parameter_tol
	    || varying_end < varying_domain.Min() - parameter_tol
	    || varying_end > varying_domain.Max() + parameter_tol
	    || std::fabs(varying_end - varying_start) <= parameter_tol)
	return false;

    /*
     * Isoparametric flags describe the locus, but an unusual pcurve can
     * still backtrack along that locus.  The exact trace backend requires a
     * monotone parameter correspondence.
     */
    const int sample_count = 16;
    const double sign = varying_end > varying_start ? 1.0 : -1.0;
    double previous = varying_start;
    for (int sample = 1; sample <= sample_count; ++sample) {
	const ON_3dPoint p = trim.PointAt(
		trim.Domain().ParameterAt((double)sample / sample_count));
	const double constant = constant_u ? p.x : p.y;
	const double varying = constant_u ? p.y : p.x;
	if (std::fabs(constant - constant_start) > parameter_tol
		|| sign * (varying - previous) < -parameter_tol)
	    return false;
	previous = varying;
    }

    trace->trim_id = trim_id;
    trace->face_id = face_id;
    trace->edge_id = trim.m_ei;
    trace->varying_dir = constant_u ? 1 : 0;
    trace->constant_parameter = 0.5 * (constant_start + constant_end);
    trace->varying_min = std::min(varying_start, varying_end);
    trace->varying_max = std::max(varying_start, varying_end);
    trace->reverse_to_edge = (varying_end < varying_start) != trim.m_bRev3d;
    return true;
}

static ON_NurbsCurve *
extract_iso_trace(const ON_NurbsSurface *surface,
	const struct iso_trim_trace &trace, const ON_Interval &edge_domain)
{
    if (!surface || !edge_domain.IsIncreasing())
	return NULL;

    ON_Curve *curve = surface->IsoCurve(trace.varying_dir,
	    trace.constant_parameter);
    if (!curve)
	return NULL;
    ON_NurbsCurve *nurbs = dynamic_cast<ON_NurbsCurve *>(curve);
    if (!nurbs) {
	ON_NurbsCurve *converted = ON_NurbsCurve::New();
	if (curve->GetNurbForm(*converted, 0.0) <= 0) {
	    delete converted;
	    delete curve;
	    return NULL;
	}
	delete curve;
	nurbs = converted;
    }

    const ON_Interval keep(trace.varying_min, trace.varying_max);
    const ON_Interval curve_domain = nurbs->Domain();
    const double domain_tol = 1.0e-12
	* std::max(1.0, curve_domain.Length());
    if (keep.Min() > curve_domain.Min() + domain_tol
	    || keep.Max() < curve_domain.Max() - domain_tol) {
	if (!nurbs->Trim(keep)) {
	    delete nurbs;
	    return NULL;
	}
    }
    if (trace.reverse_to_edge && !nurbs->Reverse()) {
	delete nurbs;
	return NULL;
    }
    if (!nurbs->SetDomain(edge_domain.Min(), edge_domain.Max())) {
	delete nurbs;
	return NULL;
    }
    return nurbs;
}

static bool
same_curve_space(const ON_NurbsCurve *a, const ON_NurbsCurve *b,
	double tolerance = 1.0e-10)
{
    if (!a || !b || a->Dimension() != b->Dimension()
	    || a->IsRational() != b->IsRational()
	    || a->Order() != b->Order() || a->CVCount() != b->CVCount()
	    || a->KnotCount() != b->KnotCount())
	return false;

    const double scale = std::max(1.0,
	    std::max(a->Domain().Length(), b->Domain().Length()));
    for (int knot = 0; knot < a->KnotCount(); ++knot) {
	if (std::fabs(a->Knot(knot) - b->Knot(knot))
		> tolerance * scale)
	    return false;
    }
    for (int cv = 0; cv < a->CVCount(); ++cv) {
	if (std::fabs(a->Weight(cv) - b->Weight(cv)) > tolerance
		* std::max(1.0, std::max(std::fabs(a->Weight(cv)),
			std::fabs(b->Weight(cv)))))
	    return false;
    }
    return true;
}

static bool
translate_surface_cv(ON_NurbsSurface *surface, int cv_u, int cv_v,
	const ON_3dVector &delta)
{
    if (!surface)
	return false;
    ON_4dPoint cv;
    if (!surface->GetCV(cv_u, cv_v, ON::euclidean_rational, &cv.x))
	return false;
    cv.x += delta.x;
    cv.y += delta.y;
    cv.z += delta.z;
    return surface->SetCV(cv_u, cv_v, ON::euclidean_rational, &cv.x);
}

static bool
trace_cv_influence(const ON_NurbsSurface *surface, int cv_u, int cv_v,
	const struct iso_trim_trace &trace, const ON_Interval &edge_domain,
	const ON_NurbsCurve *baseline, std::vector<double> *influence,
	bool *moves_endpoint)
{
    if (!surface || !baseline || !influence || !moves_endpoint)
	return false;

    ON_Surface *surface_copy = surface->DuplicateSurface();
    ON_NurbsSurface *nurbs_copy =
	dynamic_cast<ON_NurbsSurface *>(surface_copy);
    if (!nurbs_copy
	    || !translate_surface_cv(nurbs_copy, cv_u, cv_v,
		ON_3dVector(1.0, 0.0, 0.0))) {
	delete surface_copy;
	return false;
    }
    ON_NurbsCurve *changed =
	extract_iso_trace(nurbs_copy, trace, edge_domain);
    delete surface_copy;
    if (!changed || !same_curve_space(baseline, changed)) {
	delete changed;
	return false;
    }

    influence->assign((size_t)baseline->CVCount(), 0.0);
    for (int cv = 0; cv < baseline->CVCount(); ++cv) {
	ON_4dPoint before;
	ON_4dPoint after;
	if (!baseline->GetCV(cv, ON::euclidean_rational, &before.x)
		|| !changed->GetCV(cv, ON::euclidean_rational, &after.x)) {
	    delete changed;
	    return false;
	}
	(*influence)[(size_t)cv] = after.x - before.x;
	if (std::fabs(after.y - before.y) > 1.0e-10
		|| std::fabs(after.z - before.z) > 1.0e-10) {
	    delete changed;
	    return false;
	}
    }

    const ON_3dPoint before_start = baseline->PointAtStart();
    const ON_3dPoint before_end = baseline->PointAtEnd();
    const ON_3dPoint after_start = changed->PointAtStart();
    const ON_3dPoint after_end = changed->PointAtEnd();
    *moves_endpoint = before_start.DistanceTo(after_start) > 1.0e-9
	|| before_end.DistanceTo(after_end) > 1.0e-9;
    delete changed;
    return true;
}

struct coupled_iso_plan {
    int source_face = -1;
    int source_trim = -1;
    int mate_face = -1;
    int mate_trim = -1;
    int edge_id = -1;
    int source_cv_u = -1;
    int source_cv_v = -1;
    struct iso_trim_trace source_trace;
    struct iso_trim_trace mate_trace;
    std::vector<std::pair<int, int>> mate_cvs;
    std::vector<double> coefficients;
};

static bool
prepare_coupled_iso_plan(const ON_Brep *brep, int face_id, int cv_u,
	int cv_v, struct coupled_iso_plan *plan)
{
    if (!brep || !plan || !brep->IsValid(NULL)
	    || face_id < 0 || face_id >= brep->m_F.Count())
	return false;

    std::vector<int> influenced;
    if (!face_cv_influenced_trims(brep, face_id, cv_u, cv_v, &influenced)
	    || influenced.size() != 1)
	return false;
    const int source_trim_id = influenced[0];
    const ON_BrepTrim &source_trim = brep->m_T[source_trim_id];
    if (source_trim.m_type != ON_BrepTrim::mated
	    || source_trim.m_ei < 0
	    || source_trim.m_ei >= brep->m_E.Count())
	return false;
    const ON_BrepEdge &edge = brep->m_E[source_trim.m_ei];
    if (edge.m_ti.Count() != 2 || edge.m_c3i < 0
	    || edge.m_c3i >= brep->m_C3.Count())
	return false;

    const int mate_trim_id = edge.m_ti[0] == source_trim_id
	? edge.m_ti[1] : edge.m_ti[0];
    if (mate_trim_id < 0 || mate_trim_id >= brep->m_T.Count())
	return false;
    const ON_BrepTrim &mate_trim = brep->m_T[mate_trim_id];
    if (mate_trim.m_type != ON_BrepTrim::mated)
	return false;
    const int mate_face_id = trim_face_index(brep, mate_trim_id);
    if (mate_face_id < 0 || mate_face_id == face_id)
	return false;

    struct iso_trim_trace source_trace;
    struct iso_trim_trace mate_trace;
    if (!get_iso_trim_trace(brep, source_trim_id, &source_trace)
	    || !get_iso_trim_trace(brep, mate_trim_id, &mate_trace))
	return false;

    const ON_BrepFace &source_face = brep->m_F[face_id];
    const ON_BrepFace &mate_face = brep->m_F[mate_face_id];
    if (source_face.m_si < 0 || source_face.m_si >= brep->m_S.Count()
	    || mate_face.m_si < 0 || mate_face.m_si >= brep->m_S.Count()
	    || source_face.m_si == mate_face.m_si
	    || brep->SurfaceUseCount(source_face.m_si, 2) > 1
	    || brep->SurfaceUseCount(mate_face.m_si, 2) > 1)
	return false;
    const ON_NurbsSurface *source_surface =
	dynamic_cast<const ON_NurbsSurface *>(source_face.SurfaceOf());
    const ON_NurbsSurface *mate_surface =
	dynamic_cast<const ON_NurbsSurface *>(mate_face.SurfaceOf());
    if (!source_surface || !mate_surface)
	return false;

    ON_NurbsCurve *source_curve =
	extract_iso_trace(source_surface, source_trace, edge.Domain());
    ON_NurbsCurve *mate_curve =
	extract_iso_trace(mate_surface, mate_trace, edge.Domain());
    if (!source_curve || !mate_curve
	    || !same_curve_space(source_curve, mate_curve)) {
	delete source_curve;
	delete mate_curve;
	return false;
    }

    std::vector<double> source_influence;
    bool source_moves_endpoint = false;
    if (!trace_cv_influence(source_surface, cv_u, cv_v, source_trace,
	    edge.Domain(), source_curve, &source_influence,
	    &source_moves_endpoint)
	    || source_moves_endpoint) {
	delete source_curve;
	delete mate_curve;
	return false;
    }

    std::vector<std::pair<int, int>> candidates;
    std::vector<std::vector<double>> columns;
    for (int i = 0; i < mate_surface->CVCount(0); ++i) {
	for (int j = 0; j < mate_surface->CVCount(1); ++j) {
	    std::vector<int> candidate_trims;
	    if (!face_cv_influenced_trims(brep, mate_face_id, i, j,
		    &candidate_trims)
		    || candidate_trims.size() != 1
		    || candidate_trims[0] != mate_trim_id)
		continue;
	    std::vector<double> influence;
	    bool moves_endpoint = false;
	    if (!trace_cv_influence(mate_surface, i, j, mate_trace,
		    edge.Domain(), mate_curve, &influence, &moves_endpoint))
		continue;
	    double norm = 0.0;
	    for (double value : influence)
		norm += value * value;
	    if (norm <= 1.0e-24)
		continue;
	    candidates.push_back(std::make_pair(i, j));
	    columns.push_back(influence);
	}
    }
    if (candidates.empty()) {
	delete source_curve;
	delete mate_curve;
	return false;
    }

    Eigen::MatrixXd matrix(source_curve->CVCount(), (int)candidates.size());
    Eigen::VectorXd rhs(source_curve->CVCount());
    for (int row = 0; row < source_curve->CVCount(); ++row) {
	rhs(row) = source_influence[(size_t)row];
	for (int col = 0; col < (int)candidates.size(); ++col)
	    matrix(row, col) = columns[(size_t)col][(size_t)row];
    }
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> decomposition(matrix);
    const Eigen::VectorXd solution = decomposition.solve(rhs);
    const double residual = (matrix * solution - rhs).norm();
    const double residual_limit = 1.0e-10 * std::max(1.0, rhs.norm());
    if (!solution.allFinite() || residual > residual_limit) {
	delete source_curve;
	delete mate_curve;
	return false;
    }
    for (int i = 0; i < solution.size(); ++i) {
	if (std::fabs(solution(i)) > 1.0e6) {
	    delete source_curve;
	    delete mate_curve;
	    return false;
	}
    }

    plan->source_face = face_id;
    plan->source_trim = source_trim_id;
    plan->mate_face = mate_face_id;
    plan->mate_trim = mate_trim_id;
    plan->edge_id = source_trim.m_ei;
    plan->source_cv_u = cv_u;
    plan->source_cv_v = cv_v;
    plan->source_trace = source_trace;
    plan->mate_trace = mate_trace;
    plan->mate_cvs = candidates;
    plan->coefficients.resize(candidates.size());
    for (int i = 0; i < solution.size(); ++i)
	plan->coefficients[(size_t)i] = solution(i);

    delete source_curve;
    delete mate_curve;
    return true;
}

static bool
replace_edge_curve(ON_Brep *brep, int edge_id, ON_NurbsCurve *curve)
{
    if (!brep || !curve || edge_id < 0 || edge_id >= brep->m_E.Count())
	return false;
    ON_BrepEdge &edge = brep->m_E[edge_id];
    const int old_curve_id = edge.m_c3i;
    if (old_curve_id < 0 || old_curve_id >= brep->m_C3.Count())
	return false;

    if (brep->EdgeCurveUseCount(old_curve_id, 2) == 1) {
	ON_Curve *old_curve = brep->m_C3[old_curve_id];
	brep->m_C3[old_curve_id] = curve;
	if (!edge.ChangeEdgeCurve(old_curve_id)) {
	    brep->m_C3[old_curve_id] = old_curve;
	    delete curve;
	    return false;
	}
	delete old_curve;
	return true;
    }

    const int new_curve_id = brep->AddEdgeCurve(curve);
    return new_curve_id >= 0 && edge.ChangeEdgeCurve(new_curve_id);
}

}

void *brep_create()
{
    ON_Brep *brep = new ON_Brep();
    return (void *)brep;
}

bool brep_is_valid(const ON_Brep *brep)
{
    return brep && brep->IsValid(NULL);
}

int brep_vertex_create(ON_Brep *brep, ON_3dPoint point)
{
    if (!valid_brep(brep))
	return -1;
    ON_BrepVertex& v = brep->NewVertex(point);
    v.m_tolerance = 0.0;
    return brep->m_V.Count() - 1;
}

bool brep_vertex_remove(ON_Brep *brep, int v_id)
{
    if (!valid_brep(brep) || v_id < 0 || v_id >= brep->m_V.Count()) {
	bu_log("v_id is out of range\n");
	return false;
    }

    for (int ei = 0; ei < brep->m_E.Count(); ++ei) {
	const ON_BrepEdge &edge = brep->m_E[ei];
	if (edge.m_vi[0] == v_id || edge.m_vi[1] == v_id) {
	    bu_log("vertex %d is referenced by brep edge %d\n", v_id, ei);
	    return false;
	}
    }
    for (int ti = 0; ti < brep->m_T.Count(); ++ti) {
	const ON_BrepTrim &trim = brep->m_T[ti];
	if (trim.m_vi[0] == v_id || trim.m_vi[1] == v_id) {
	    bu_log("vertex %d is referenced by brep trim %d\n", v_id, ti);
	    return false;
	}
    }

    brep->m_V.Remove(v_id);
    for (int vi = v_id; vi < brep->m_V.Count(); ++vi)
	brep->m_V[vi].m_vertex_index = vi;
    for (int ei = 0; ei < brep->m_E.Count(); ++ei) {
	ON_BrepEdge &edge = brep->m_E[ei];
	if (edge.m_vi[0] > v_id) --edge.m_vi[0];
	if (edge.m_vi[1] > v_id) --edge.m_vi[1];
    }
    for (int ti = 0; ti < brep->m_T.Count(); ++ti) {
	ON_BrepTrim &trim = brep->m_T[ti];
	if (trim.m_vi[0] > v_id) --trim.m_vi[0];
	if (trim.m_vi[1] > v_id) --trim.m_vi[1];
    }
    brep->DestroyRuntimeCache(true);
    return true;
}

int brep_curve2d_make_line(ON_Brep *brep, const ON_2dPoint &from, const ON_2dPoint &to)
{
    if (!valid_brep(brep))
	return -1;
    ON_Curve* c2 = new ON_LineCurve(from, to);
    c2->SetDomain(0.0, 1.0);
    brep->m_C2.Append(c2);
    return brep->m_C2.Count() - 1;
}

bool brep_curve2d_remove(ON_Brep *brep, int curve_id)
{
    if (!valid_brep(brep) || curve_id < 0 || curve_id >= brep->m_C2.Count()) {
	bu_log("curve_id is out of range\n");
	return false;
    }

    if (brep->TrimCurveUseCount(curve_id, 1) > 0) {
	bu_log("curve2d %d is referenced by a brep trim\n", curve_id);
	return false;
    }

    delete brep->m_C2[curve_id];
    brep->m_C2.Remove(curve_id);
    for (int ti = 0; ti < brep->m_T.Count(); ++ti) {
	if (brep->m_T[ti].m_c2i > curve_id)
	    --brep->m_T[ti].m_c2i;
    }
    brep->DestroyRuntimeCache(true);
    return true;
}

int brep_curve_make(ON_Brep *brep, const ON_3dPoint &position)
{
    if (!valid_brep(brep))
	return -1;
    ON_NurbsCurve *curve = ON_NurbsCurve::New(3, true, 3, 4);
    curve->SetCV(0, ON_3dPoint(-0.1, -1.5, 0));
    curve->SetCV(1, ON_3dPoint(0.1, -0.5, 0));
    curve->SetCV(2, ON_3dPoint(0.1, 0.5, 0));
    curve->SetCV(3, ON_3dPoint(-0.1, 1.5, 0));
    curve->MakeClampedUniformKnotVector();
    curve->Translate(position);
    return brep->AddEdgeCurve(curve);
}

// TODO: add more options about knot vector
int brep_curve_in(ON_Brep *brep, bool is_rational, int order, int cv_count, std::vector<ON_4dPoint> cv)
{
    int dim = 3;
    if (cv.size() != (size_t)cv_count) {
	bu_log("cv_count is not equal to cv.size()\n");
	return -1;
    }
    ON_NurbsCurve *curve = ON_NurbsCurve::New(dim, is_rational, order, cv_count);

    for (int i = 0; i < cv_count; i++) {
	curve->SetCV(i, cv[i]);
    }

    // make uniform knot vector
    curve->MakeClampedUniformKnotVector();
    return brep->AddEdgeCurve(curve);
}

ON_NurbsCurve *brep_get_nurbs_curve(ON_Brep *brep, int curve_id)
{
    if (!valid_brep(brep) || curve_id < 0 || curve_id >= brep->m_C3.Count()) {
	bu_log("curve_id is out of range\n");
	return NULL;
    }
    ON_NurbsCurve *curve = dynamic_cast<ON_NurbsCurve *>(brep->m_C3[curve_id]);

    if (!curve) {
	bu_log("curve %d is not a NURBS curve\n", curve_id);
	return NULL;
    }
    return curve;
}

ON_NurbsSurface *brep_get_nurbs_surface(ON_Brep *brep, int surface_id)
{
    if (!valid_brep(brep) || surface_id < 0 || surface_id >= brep->m_S.Count()) {
	bu_log("surface_id is out of range\n");
	return NULL;
    }
    ON_NurbsSurface *surface = dynamic_cast<ON_NurbsSurface *>(brep->m_S[surface_id]);

    if (!surface) {
	bu_log("surface %d is not a NURBS surface\n", surface_id);
	return NULL;
    }
    return surface;
}

/**
 * calculate tangent vectors for each point
 */
std::vector<ON_3dVector> calculateTangentVectors(const std::vector<ON_3dPoint> &points);

/**
 * calculate control points and knot vector for interpolating a curve
 */
void calcuBsplineCVsKnots(std::vector<ON_3dPoint> &cvs, std::vector<double> &knots, const std::vector<ON_3dPoint> &endPoints, const std::vector<ON_3dVector> &tangentVectors);

int brep_curve_interpCrv(ON_Brep *brep, std::vector<ON_3dPoint> points)
{
    std::vector<ON_3dVector> tangentVectors = calculateTangentVectors(points);
    std::vector<ON_3dPoint> controlPs;
    std::vector<double> knotVector;

    calcuBsplineCVsKnots(controlPs, knotVector, points, tangentVectors);
    ON_NurbsCurve *curve = ON_NurbsCurve::New(3, false, 4, controlPs.size());
    for (size_t i = 0; i < controlPs.size(); i++) {
	curve->SetCV(i, controlPs[i]);
    }
    for(size_t i = 0; i < knotVector.size(); i++) {
	curve->SetKnot(i, knotVector[i]);
    }
    return brep->AddEdgeCurve(curve);
}

int brep_curve_copy(ON_Brep *brep, int curve_id)
{
    if (curve_id < 0 || curve_id >= brep->m_C3.Count()) {
	bu_log("curve_id is out of range\n");
	return -1;
    }
    ON_Curve *curve = brep->m_C3[curve_id];
    if (!curve) {
	return -1;
    }
    ON_Curve *curve_copy = curve->DuplicateCurve();
    return brep->AddEdgeCurve(curve_copy);
}

bool brep_curve_remove(ON_Brep *brep, int curve_id)
{
    if (!valid_brep(brep) || curve_id < 0 || curve_id >= brep->m_C3.Count()) {
	bu_log("curve_id is out of range\n");
	return false;
    }

    if (curve3d_referenced(brep, curve_id)) {
	bu_log("curve %d is referenced by a brep edge\n", curve_id);
	return false;
    }

    delete brep->m_C3[curve_id];
    brep->m_C3.Remove(curve_id);
    for (int ei = 0; ei < brep->m_E.Count(); ++ei) {
	if (brep->m_E[ei].m_c3i > curve_id)
	    --brep->m_E[ei].m_c3i;
    }
    brep->DestroyRuntimeCache(true);
    return true;
}

bool brep_curve_move(ON_Brep *brep, int curve_id, const ON_3dVector &point)
{
    if (!valid_brep(brep) || curve_id < 0 || curve_id >= brep->m_C3.Count()) {
	bu_log("curve_id is out of range\n");
	return false;
    }
    if (curve3d_referenced(brep, curve_id)) {
	bu_log("curve %d is referenced by a brep edge\n", curve_id);
	return false;
    }
    ON_Curve *curve = brep->m_C3[curve_id];
    if (!curve) {
	return false;
    }
    ON_Curve *curve_copy = curve->DuplicateCurve();
    if (!curve_copy || !curve_copy->Translate(point)) {
	delete curve_copy;
	return false;
    }
    return replace_curve3d(brep, curve_id, curve_copy);
}

bool brep_curve_set_cv(ON_Brep *brep, int curve_id, int cv_id, const ON_4dPoint &point)
{
    return edit_nurbs_curve(brep, curve_id, false,
	    [cv_id, &point](ON_NurbsCurve *curve) {
		return curve->SetCV(cv_id, point);
	    });
}

bool brep_curve_set_weight(ON_Brep *brep, int curve_id, int cv_id, double weight)
{
    return edit_nurbs_curve(brep, curve_id, false,
	    [cv_id, weight](ON_NurbsCurve *curve) {
		return ON_IsValid(weight) && weight > 0.0 && curve->SetWeight(cv_id, weight);
	    });
}

bool brep_curve_reverse(ON_Brep *brep, int curve_id)
{
    return edit_nurbs_curve(brep, curve_id, false,
	    [](ON_NurbsCurve *curve) {
		return curve->Reverse();
	    });
}

bool brep_curve_insert_knot(ON_Brep *brep, int curve_id, double knot, int multiplicity)
{
    return edit_nurbs_curve(brep, curve_id, true,
	    [knot, multiplicity](ON_NurbsCurve *curve) {
		return curve->InsertKnot(knot, multiplicity);
	    });
}

bool brep_curve_clamp(ON_Brep *brep, int curve_id, int end)
{
    if (end < 0 || end > 2)
	return false;
    return edit_nurbs_curve(brep, curve_id, true,
	    [end](ON_NurbsCurve *curve) {
		return curve->ClampEnd(end);
	    });
}

bool brep_curve_elevate_degree(ON_Brep *brep, int curve_id, int desired_degree)
{
    return edit_nurbs_curve(brep, curve_id, true,
	    [desired_degree](ON_NurbsCurve *curve) {
		return desired_degree >= curve->Degree()
		    && curve->IncreaseDegree(desired_degree);
	    });
}

bool brep_curve_trim(ON_Brep *brep, int curve_id, double t0, double t1)
{
    return edit_nurbs_curve(brep, curve_id, false,
	    [t0, t1](ON_NurbsCurve *curve) {
		ON_Interval interval(t0, t1);
		return interval.IsIncreasing() && curve->Trim(interval);
	    });
}

bool brep_curve_split(ON_Brep *brep, int curve_id, double t)
{
    ON_NurbsCurve *curve = brep_get_nurbs_curve(brep, curve_id);
    if (!curve || curve3d_referenced(brep, curve_id)) {
	return false;
    }
    ON_Curve *curve1 = NULL;
    ON_Curve *curve2 = NULL;
    bool flag = curve->Split(t, curve1, curve2);
    if (flag) {
	brep_curve_remove(brep, curve_id);
	brep->AddEdgeCurve(curve1);
	brep->AddEdgeCurve(curve2);
	bu_log("old curve removed, id: %d, new curve id: %d, %d\n", curve_id, brep->m_C3.Count() - 2, brep->m_C3.Count() - 1);
    }
    return flag;
}

int brep_curve_join(ON_Brep *brep, int curve_id_1, int curve_id_2)
{
    ON_NurbsCurve *curve1 = brep_get_nurbs_curve(brep, curve_id_1);
    ON_NurbsCurve *curve2 = brep_get_nurbs_curve(brep, curve_id_2);
    if (!curve1 || !curve2 || curve_id_1 == curve_id_2
	    || curve3d_referenced(brep, curve_id_1)
	    || curve3d_referenced(brep, curve_id_2)) {
	return -1;
    }

    ON_NurbsCurve curve1_copy(*curve1);
    ON_NurbsCurve curve2_copy(*curve2);

    /// force move ends of the two curves to the same location
    if (!ON_ForceMatchCurveEnds(curve1_copy, 1, curve2_copy, 0)) {
	bu_log("ON_ForceMatchCurveEnds failed\n");
	return -1;
    }

    ON_SimpleArray<const ON_Curve *> in_curves;
    in_curves.Append(&curve1_copy);
    in_curves.Append(&curve2_copy);

    ON_SimpleArray<ON_Curve *> out_curves;
    /// join the two curves
    int joined_num = ON_JoinCurves(in_curves, out_curves, 0.0f, 0.0f, false);
    if (joined_num != 1) {
	bu_log("ON_JoinCurves failed\n");
	return -1;
    }

    /// remove the two curves.
    /// Remark: the index of m_C3 is massed up after removing the curves!
    int high_id = curve_id_1 > curve_id_2 ? curve_id_1 : curve_id_2;
    int low_id = curve_id_1 < curve_id_2 ? curve_id_1 : curve_id_2;
    if (!brep_curve_remove(brep, high_id) || !brep_curve_remove(brep, low_id)) {
	for (int i = 0; i < out_curves.Count(); ++i)
	    delete out_curves[i];
	return -1;
    }

    return brep->AddEdgeCurve(out_curves[0]);
}

int brep_surface_make(ON_Brep *brep, const ON_3dPoint &position)
{
    ON_NurbsSurface *surface = ON_NurbsSurface::New(3, true, 3, 3, 4, 4);
    surface->SetCV(0, 0, ON_4dPoint(-1.5, -1.5, 0, 1));
    surface->SetCV(0, 1, ON_4dPoint(-0.5, -1.5, 0, 1));
    surface->SetCV(0, 2, ON_4dPoint(0.5, -1.5, 0, 1));
    surface->SetCV(0, 3, ON_4dPoint(1.5, -1.5, 0, 1));
    surface->SetCV(1, 0, ON_4dPoint(-1.5, -0.5, 0, 1));
    surface->SetCV(1, 1, ON_4dPoint(-0.5, -0.5, 0.5, 1));
    surface->SetCV(1, 2, ON_4dPoint(0.5, -0.5, 0.5, 1));
    surface->SetCV(1, 3, ON_4dPoint(1.5, -0.5, 0, 1));
    surface->SetCV(2, 0, ON_4dPoint(-1.5, 0.5, 0, 1));
    surface->SetCV(2, 1, ON_4dPoint(-0.5, 0.5, 0.5, 1));
    surface->SetCV(2, 2, ON_4dPoint(0.5, 0.5, 0.5, 1));
    surface->SetCV(2, 3, ON_4dPoint(1.5, 0.5, 0, 1));
    surface->SetCV(3, 0, ON_4dPoint(-1.5, 1.5, 0, 1));
    surface->SetCV(3, 1, ON_4dPoint(-0.5, 1.5, 0, 1));
    surface->SetCV(3, 2, ON_4dPoint(0.5, 1.5, 0, 1));
    surface->SetCV(3, 3, ON_4dPoint(1.5, 1.5, 0, 1));
    surface->MakeClampedUniformKnotVector(0, 1);
    surface->MakeClampedUniformKnotVector(1, 1);
    surface->Translate(position);
    return brep->AddSurface(surface);
}

int brep_surface_extract_vertex(ON_Brep *brep, int surface_id, double u, double v)
{
    ON_NurbsSurface *surface = brep_get_nurbs_surface(brep, surface_id);
    if (!surface) {
	return -1;
    }
    ON_3dPoint point;
    bool res = surface->Evaluate(u, v, 0, 3, point);
    if(!res) {
	return -1;
    }
    ON_BrepVertex& vertex = brep->NewVertex(point);
    vertex.m_tolerance = 0.0;
    return brep->m_V.Count() - 1;
}

int brep_surface_extract_curve(ON_Brep *brep, int surface_id, int dir, double t)
{
    ON_NurbsSurface *surface = brep_get_nurbs_surface(brep, surface_id);
    if (!surface) {
	return -1;
    }
    ON_Curve *curve = surface->IsoCurve(dir, t);
    if(!curve) {
	return -1;
    }
    return brep->AddEdgeCurve(curve);
}

int brep_surface_copy(ON_Brep *brep, int surface_id)
{
    if (surface_id < 0 || surface_id >= brep->m_S.Count()) {
	bu_log("surface_id is out of range\n");
	return -1;
    }
    ON_Surface *surface = brep->m_S[surface_id];
    if (!surface) {
	return -1;
    }
    ON_Surface *surface_copy = surface->DuplicateSurface();
    return brep->AddSurface(surface_copy);
}

/**
 * calculate parameter values of each point
 */
int SurfMeshParams(int n, int m, std::vector<ON_3dPoint> points, std::vector<double> &uk, std::vector<double> &ul);

/**
 * global cubic curve interpolation with C2 continuity
 * input: n, knots, Q
 * output: P
 * reference: The NURBS Book (2nd Edition), chapter 9.2.3
 * TODO: while n <= 3, the special case should be considered
 */
void globalCurveInterp(int n, std::vector<double> &knots, const std::vector<ON_3dPoint> &Q, std::vector<ON_3dPoint> &P);

int brep_surface_interpCrv(ON_Brep *brep, int cv_count_x, int cv_count_y, std::vector<ON_3dPoint> points)
{
    cv_count_x = cv_count_x < 2 ? 2 : cv_count_x;
    cv_count_y = cv_count_y < 2 ? 2 : cv_count_y;
    if (points.size() != (size_t)(cv_count_x * cv_count_y)) {
	return -1;
    }
    int n = cv_count_x - 1;
    int m = cv_count_y - 1;

    /// calculate parameter values of each point
    std::vector<double> uk, ul;
    SurfMeshParams(n, m, points, uk, ul);

    /// calculate knots of the cubic B-spline surface
    std::vector<double> knots_u, knots_v;
    for (size_t i = 0; i < 4; i++) {
	knots_u.push_back(0);
	knots_v.push_back(0);
    }
    for (size_t i = 1; i < uk.size() - 1; i++) {
	knots_u.push_back(uk[i]);
    }
    for (size_t i = 1; i < ul.size() - 1; i++) {
	knots_v.push_back(ul[i]);
    }
    for (int i = 0; i < 4; i++) {
	knots_u.push_back(1);
	knots_v.push_back(1);
    }

    /// curve interpolation in v direction
    // temporary control points
    std::vector<std::vector<ON_3dPoint>> R(n + 1, std::vector<ON_3dPoint>(m + 3));

    for (int l = 0; l <= n; l++) {
	std::vector<ON_3dPoint> Q(m + 1);
	for (int k = 0; k <= m; k++) {
	    Q[k] = points[l * (m + 1) + k];
	}
	globalCurveInterp(m, knots_v, Q, R[l]);
    }

    ON_NurbsSurface *surface = ON_NurbsSurface::New(3, false, 4, 4, n + 3, m + 3);
    for (int i = 0; i < m + 3; i++) {
	std::vector<ON_3dPoint> Q(n + 1);
	for (int j = 0; j < n + 1; j++) {
	    Q[j] = R[j][i];
	}
	std::vector<ON_3dPoint> P(n + 3);
	globalCurveInterp(n, knots_u, Q, P);
	for (int j = 0; j < n + 3; j++) {
	    surface->SetCV(j, i, P[j]);
	}
    }
    /// the knot vector of openNURBS is different from the NURBS book
    for (size_t i = 1; i < knots_u.size() - 1; i++) {
	surface->SetKnot(0, i - 1, knots_u[i]);
    }
    for (size_t i = 1; i < knots_v.size() - 1; i++) {
	surface->SetKnot(1, i - 1, knots_v[i]);
    }
    return brep->AddSurface(surface);
}

bool brep_surface_move(ON_Brep *brep, int surface_id, const ON_3dVector &point)
{
    if (!valid_brep(brep) || surface_id < 0 || surface_id >= brep->m_S.Count()) {
	bu_log("surface_id is out of range\n");
	return false;
    }
    if (surface_referenced(brep, surface_id)) {
	bu_log("surface %d is referenced by a brep face\n", surface_id);
	return false;
    }
    ON_Surface *surface = brep->m_S[surface_id];
    if (!surface) {
	return false;
    }
    ON_Surface *surface_copy = surface->DuplicateSurface();
    if (!surface_copy || !surface_copy->Translate(point)) {
	delete surface_copy;
	return false;
    }
    return replace_surface(brep, surface_id, surface_copy);
}

bool brep_surface_set_cv(ON_Brep *brep, int surface_id, int cv_id_u, int cv_id_v, const ON_4dPoint &point)
{
    return edit_nurbs_surface(brep, surface_id, false,
	    [cv_id_u, cv_id_v, &point](ON_NurbsSurface *surface) {
		return surface->SetCV(cv_id_u, cv_id_v, point);
	    });
}

bool brep_surface_set_weight(ON_Brep *brep, int surface_id, int cv_id_u, int cv_id_v, double weight)
{
    return edit_nurbs_surface(brep, surface_id, false,
	    [cv_id_u, cv_id_v, weight](ON_NurbsSurface *surface) {
		return ON_IsValid(weight) && weight > 0.0
		    && surface->SetWeight(cv_id_u, cv_id_v, weight);
	    });
}

bool brep_surface_reverse(ON_Brep *brep, int surface_id, int dir)
{
    if (dir < 0 || dir > 1)
	return false;
    return edit_nurbs_surface(brep, surface_id, false,
	    [dir](ON_NurbsSurface *surface) {
		return surface->Reverse(dir);
	    });
}

bool brep_surface_transpose(ON_Brep *brep, int surface_id)
{
    return edit_nurbs_surface(brep, surface_id, false,
	    [](ON_NurbsSurface *surface) {
		return surface->Transpose();
	    });
}

bool brep_surface_insert_knot(ON_Brep *brep, int surface_id, int dir, double knot, int multiplicity)
{
    if (dir < 0 || dir > 1)
	return false;
    return edit_nurbs_surface(brep, surface_id, true,
	    [dir, knot, multiplicity](ON_NurbsSurface *surface) {
		return surface->InsertKnot(dir, knot, multiplicity);
	    });
}

bool brep_surface_clamp(ON_Brep *brep, int surface_id, int dir, int end)
{
    if (dir < 0 || dir > 1 || end < 0 || end > 2)
	return false;
    return edit_nurbs_surface(brep, surface_id, true,
	    [dir, end](ON_NurbsSurface *surface) {
		return surface->ClampEnd(dir, end);
	    });
}

bool brep_surface_elevate_degree(ON_Brep *brep, int surface_id, int dir, int desired_degree)
{
    if (dir < 0 || dir > 1)
	return false;
    return edit_nurbs_surface(brep, surface_id, true,
	    [dir, desired_degree](ON_NurbsSurface *surface) {
		return desired_degree >= surface->Degree(dir)
		    && surface->IncreaseDegree(dir, desired_degree);
	    });
}

bool brep_surface_trim(ON_Brep *brep, int surface_id, int dir, double t0, double t1)
{
    if (dir < 0 || dir > 1)
	return false;
    return edit_nurbs_surface(brep, surface_id, false,
	    [dir, t0, t1](ON_NurbsSurface *surface) {
		ON_Interval interval(t0, t1);
		return interval.IsIncreasing() && surface->Trim(dir, interval);
	    });
}

bool brep_surface_split(ON_Brep *brep, int surface_id, int dir, double t)
{
    ON_NurbsSurface *surface = brep_get_nurbs_surface(brep, surface_id);
    if (!surface || dir < 0 || dir > 1 || surface_referenced(brep, surface_id)) {
	return false;
    }
    ON_Surface *surf1=NULL, *surf2=NULL;
    bool flag = surface->Split(dir, t, surf1, surf2);
    if (flag) {
	brep_surface_remove(brep, surface_id);
	brep->AddSurface(surf1);
	brep->AddSurface(surf2);
    bu_log("brep_surface_split: split surface %d into %d and %d\n", surface_id, brep->m_S.Count()-2, brep->m_S.Count()-1);
    }
    return flag;
}

int brep_surface_create_ruled(ON_Brep *brep, int curve_id0, int curve_id1)
{
    ON_NurbsCurve *curve0 = brep_get_nurbs_curve(brep, curve_id0);
    ON_NurbsCurve *curve1 = brep_get_nurbs_curve(brep, curve_id1);
    if (!curve0 || !curve1) {
	return -1;
    }
    ON_NurbsSurface *surface = ON_NurbsSurface::New();
    if (!surface->CreateRuledSurface(*curve0, *curve1)) {
	delete surface;
	return -1;
    }
    return brep->AddSurface(surface);
}

int brep_surface_tensor_product(ON_Brep *brep, int curve_id0, int curve_id1)
{
    ON_NurbsCurve *curve0 = brep_get_nurbs_curve(brep, curve_id0);
    ON_NurbsCurve *curve1 = brep_get_nurbs_curve(brep, curve_id1);
    if (!curve0 || !curve1) {
	return -1;
    }
    ON_SumSurface *surface = ON_SumSurface::New();
    if(!surface->Create(*curve0, *curve1)) {
	delete surface;
	return -1;
    }
    ON_NurbsSurface *nurbs_surface = surface->NurbsSurface();
    delete surface;
    if(nurbs_surface == NULL) {
	return -1;
    }
    return brep->AddSurface(nurbs_surface);
}

int brep_surface_revolution(ON_Brep *brep, int curve_id0, ON_3dPoint line_start, ON_3dPoint line_end, double angle)
{
    ON_NurbsCurve *curve0 = brep_get_nurbs_curve(brep, curve_id0);
    ON_Line line(line_start, line_end);
    if (!curve0 || !line.IsValid() || !ON_IsValid(angle)
	    || std::fabs(angle) <= ON_ZERO_TOLERANCE) {
	return -1;
    }

    if (angle < 0.0) {
	angle = -angle;
	line.Reverse();
    }
    if (angle > 2 * ON_PI)
	angle = 2 * ON_PI;

    ON_RevSurface *rev_surf = ON_RevSurface::New();
    rev_surf->m_curve = curve0->DuplicateCurve();
    if (!rev_surf->m_curve) {
	delete rev_surf;
	return -1;
    }
    rev_surf->m_axis = line;
    rev_surf->m_angle = ON_Interval(0, angle);

    ON_NurbsSurface *nurbs_surface = ON_NurbsSurface::New();
    const int nurbs_status = rev_surf->GetNurbForm(*nurbs_surface, 0.0);
    delete rev_surf;
    if (nurbs_status <= 0 || !nurbs_surface->IsValid(NULL)) {
	delete nurbs_surface;
	return -1;
    }
    return brep->AddSurface(nurbs_surface);
}

bool brep_surface_remove(ON_Brep *brep, int surface_id)
{
    if (!valid_brep(brep) || surface_id < 0 || surface_id >= brep->m_S.Count()) {
	bu_log("surface_id is out of range\n");
	return false;
    }

    if (surface_referenced(brep, surface_id)) {
	bu_log("surface %d is referenced by a brep face\n", surface_id);
	return false;
    }

    delete brep->m_S[surface_id];
    brep->m_S.Remove(surface_id);
    for (int fi = 0; fi < brep->m_F.Count(); ++fi) {
	if (brep->m_F[fi].m_si > surface_id)
	    --brep->m_F[fi].m_si;
    }
    brep->DestroyRuntimeCache(true);
    return true;
}

int brep_edge_create(ON_Brep *brep, int from, int to, int curve)
{
    ON_BrepVertex& v0 = brep->m_V[from];
    ON_BrepVertex& v1 = brep->m_V[to];
    ON_BrepEdge& edge = brep->NewEdge(v0, v1, curve);
    edge.m_tolerance = 0.0;
    return brep->m_E.Count() - 1;
}

int brep_face_create(ON_Brep *brep, int surface)
{
    if(surface < 0 || surface >= brep->m_S.Count()) {
	bu_log("surface is out of range\n");
	return -1;
    }
    brep->NewFace(surface);
    return brep->m_F.Count() - 1;
}

bool brep_face_reverse(ON_Brep *brep, int face)
{
    if(face < 0 || face >= brep->m_F.Count()) {
	bu_log("face is out of range\n");
	return false;
    }
    ON_BrepFace& f = brep->m_F[face];
    f.m_bRev = !f.m_bRev;
    return true;
}

bool brep_face_cv_is_topology_safe(const ON_Brep *brep, int face, int cv_id_u, int cv_id_v)
{
    std::vector<int> trim_ids;
    return face_cv_influenced_trims(brep, face, cv_id_u, cv_id_v, &trim_ids)
	&& trim_ids.empty();
}

int brep_edge_constraint_type(const ON_Brep *brep, int edge_id)
{
    if (!brep || edge_id < 0 || edge_id >= brep->m_E.Count())
	return BREP_EDGE_CONSTRAINT_INVALID;
    const ON_BrepEdge &edge = brep->m_E[edge_id];
    if (edge.m_ti.Count() == 1) {
	const int trim_id = edge.m_ti[0];
	if (trim_id < 0 || trim_id >= brep->m_T.Count())
	    return BREP_EDGE_CONSTRAINT_INVALID;
	const ON_BrepTrim &trim = brep->m_T[trim_id];
	return trim.m_type == ON_BrepTrim::seam
	    || trim.m_type == ON_BrepTrim::singular
	    ? BREP_EDGE_CONSTRAINT_SPECIAL : BREP_EDGE_CONSTRAINT_OPEN;
    }
    if (edge.m_ti.Count() != 2)
	return BREP_EDGE_CONSTRAINT_SPECIAL;

    const int trim0_id = edge.m_ti[0];
    const int trim1_id = edge.m_ti[1];
    if (trim0_id < 0 || trim0_id >= brep->m_T.Count()
	    || trim1_id < 0 || trim1_id >= brep->m_T.Count())
	return BREP_EDGE_CONSTRAINT_INVALID;
    const ON_BrepTrim &trim0 = brep->m_T[trim0_id];
    const ON_BrepTrim &trim1 = brep->m_T[trim1_id];
    if (trim0.m_type != ON_BrepTrim::mated
	    || trim1.m_type != ON_BrepTrim::mated
	    || trim_face_index(brep, trim0_id) == trim_face_index(brep, trim1_id))
	return BREP_EDGE_CONSTRAINT_SPECIAL;

    const bool iso0 = iso_flag(trim0.m_iso);
    const bool iso1 = iso_flag(trim1.m_iso);
    if (natural_iso_flag(trim0.m_iso) && natural_iso_flag(trim1.m_iso))
	return BREP_EDGE_CONSTRAINT_NATURAL;
    if (iso0 && iso1)
	return BREP_EDGE_CONSTRAINT_ISOPARAMETRIC;
    if (iso0 || iso1)
	return BREP_EDGE_CONSTRAINT_ONE_ISOPARAMETRIC;
    return BREP_EDGE_CONSTRAINT_GENERAL;
}

const char *
brep_edge_constraint_type_name(int classification)
{
    switch (classification) {
	case BREP_EDGE_CONSTRAINT_OPEN:
	    return "open";
	case BREP_EDGE_CONSTRAINT_NATURAL:
	    return "natural";
	case BREP_EDGE_CONSTRAINT_ISOPARAMETRIC:
	    return "isoparametric";
	case BREP_EDGE_CONSTRAINT_ONE_ISOPARAMETRIC:
	    return "one-isoparametric";
	case BREP_EDGE_CONSTRAINT_GENERAL:
	    return "general";
	case BREP_EDGE_CONSTRAINT_SPECIAL:
	    return "special";
	default:
	    return "invalid";
    }
}

const char *
brep_cv_constraint_type_name(int classification)
{
    switch (classification) {
	case BREP_CV_CONSTRAINT_INTERIOR:
	    return "interior";
	case BREP_CV_CONSTRAINT_NATURAL:
	    return "natural";
	case BREP_CV_CONSTRAINT_ISOPARAMETRIC:
	    return "isoparametric";
	case BREP_CV_CONSTRAINT_GENERAL:
	    return "general";
	case BREP_CV_CONSTRAINT_SPECIAL:
	    return "special";
	default:
	    return "invalid";
    }
}

bool
brep_face_cv_constraint_status(const ON_Brep *brep, int face_id,
	int cv_id_u, int cv_id_v, struct brep_face_cv_constraint *status)
{
    if (!status)
	return false;
    status->classification = BREP_CV_CONSTRAINT_INVALID;
    status->trim_count = 0;
    status->edge_count = 0;
    status->other_face_count = 0;
    status->natural_trim_count = 0;
    status->isoparametric_trim_count = 0;
    status->general_trim_count = 0;
    status->topology_safe = false;
    status->can_translate = false;

    std::vector<int> trim_ids;
    if (!face_cv_influenced_trims(brep, face_id, cv_id_u, cv_id_v,
	    &trim_ids))
	return false;
    status->trim_count = (int)trim_ids.size();
    status->topology_safe = trim_ids.empty();
    if (trim_ids.empty()) {
	status->classification = BREP_CV_CONSTRAINT_INTERIOR;
	status->can_translate = true;
	return true;
    }

    std::set<int> edge_ids;
    std::set<int> other_face_ids;
    bool special = false;
    for (int trim_id : trim_ids) {
	const ON_BrepTrim &trim = brep->m_T[trim_id];
	if (natural_iso_flag(trim.m_iso))
	    ++status->natural_trim_count;
	if (iso_flag(trim.m_iso))
	    ++status->isoparametric_trim_count;
	else
	    ++status->general_trim_count;
	if (trim.m_ei >= 0 && trim.m_ei < brep->m_E.Count()) {
	    edge_ids.insert(trim.m_ei);
	    const ON_BrepEdge &edge = brep->m_E[trim.m_ei];
	    for (int eti = 0; eti < edge.m_ti.Count(); ++eti) {
		const int other_trim = edge.m_ti[eti];
		const int other_face = trim_face_index(brep, other_trim);
		if (other_face >= 0 && other_face != face_id)
		    other_face_ids.insert(other_face);
	    }
	    const int edge_class = brep_edge_constraint_type(brep, trim.m_ei);
	    if (edge_class == BREP_EDGE_CONSTRAINT_SPECIAL
		    || edge_class == BREP_EDGE_CONSTRAINT_INVALID)
		special = true;
	} else {
	    special = true;
	}
    }
    status->edge_count = (int)edge_ids.size();
    status->other_face_count = (int)other_face_ids.size();
    if (special)
	status->classification = BREP_CV_CONSTRAINT_SPECIAL;
    else if (status->general_trim_count)
	status->classification = BREP_CV_CONSTRAINT_GENERAL;
    else if (status->natural_trim_count == status->trim_count)
	status->classification = BREP_CV_CONSTRAINT_NATURAL;
    else
	status->classification = BREP_CV_CONSTRAINT_ISOPARAMETRIC;

    struct coupled_iso_plan plan;
    status->can_translate =
	prepare_coupled_iso_plan(brep, face_id, cv_id_u, cv_id_v, &plan);
    return true;
}

bool
brep_face_cv_can_translate(const ON_Brep *brep, int face_id,
	int cv_id_u, int cv_id_v)
{
    struct brep_face_cv_constraint status;
    return brep_face_cv_constraint_status(brep, face_id, cv_id_u, cv_id_v,
	    &status) && status.can_translate;
}

bool brep_face_set_cv(ON_Brep *brep, int face, int cv_id_u, int cv_id_v, const ON_4dPoint &point)
{
    if (!brep_face_cv_is_topology_safe(brep, face, cv_id_u, cv_id_v)) {
	bu_log("face %d CV (%d,%d) influences a trimming edge\n",
		face, cv_id_u, cv_id_v);
	return false;
    }
    return edit_face_nurbs_surface(brep, face,
	    [cv_id_u, cv_id_v, &point](ON_NurbsSurface *surface) {
		return surface->SetCV(cv_id_u, cv_id_v,
			ON::euclidean_rational, &point.x);
	    });
}

bool brep_face_translate_cv(ON_Brep *brep, int face, int cv_id_u, int cv_id_v, const ON_3dVector &delta)
{
    if (!brep_face_cv_is_topology_safe(brep, face, cv_id_u, cv_id_v)) {
	bu_log("face %d CV (%d,%d) influences a trimming edge\n",
		face, cv_id_u, cv_id_v);
	return false;
    }
    return edit_face_nurbs_surface(brep, face,
	    [cv_id_u, cv_id_v, &delta](ON_NurbsSurface *surface) {
		ON_4dPoint cv;
		if (!surface->GetCV(cv_id_u, cv_id_v, ON::euclidean_rational, &cv.x))
		    return false;
		cv.x += delta.x;
		cv.y += delta.y;
		cv.z += delta.z;
		return surface->SetCV(cv_id_u, cv_id_v, ON::euclidean_rational, &cv.x);
	    });
}

bool
brep_face_translate_cv_constrained(ON_Brep *brep, int face_id,
	int cv_id_u, int cv_id_v, const ON_3dVector &delta)
{
    if (!brep || !delta.IsValid())
	return false;
    if (brep_face_cv_is_topology_safe(brep, face_id, cv_id_u, cv_id_v))
	return brep_face_translate_cv(brep, face_id, cv_id_u, cv_id_v, delta);

    struct coupled_iso_plan plan;
    if (!prepare_coupled_iso_plan(brep, face_id, cv_id_u, cv_id_v, &plan)) {
	bu_log("face %d CV (%d,%d) does not satisfy the exact C0 "
		"isoparametric coupling requirements\n",
		face_id, cv_id_u, cv_id_v);
	return false;
    }

    ON_Brep trial(*brep);
    ON_BrepFace &source_face = trial.m_F[plan.source_face];
    ON_BrepFace &mate_face = trial.m_F[plan.mate_face];
    ON_NurbsSurface *source_surface =
	dynamic_cast<ON_NurbsSurface *>(trial.m_S[source_face.m_si]);
    ON_NurbsSurface *mate_surface =
	dynamic_cast<ON_NurbsSurface *>(trial.m_S[mate_face.m_si]);
    if (!source_surface || !mate_surface
	    || !translate_surface_cv(source_surface, plan.source_cv_u,
		plan.source_cv_v, delta))
	return false;

    for (size_t candidate = 0; candidate < plan.mate_cvs.size(); ++candidate) {
	const ON_3dVector mate_delta =
	    plan.coefficients[candidate] * delta;
	if (!translate_surface_cv(mate_surface,
		plan.mate_cvs[candidate].first,
		plan.mate_cvs[candidate].second, mate_delta))
	    return false;
    }
    source_surface->DestroyRuntimeCache(true);
    mate_surface->DestroyRuntimeCache(true);
    source_face.DestroyRuntimeCache(true);
    mate_face.DestroyRuntimeCache(true);

    const ON_BrepEdge &old_edge = brep->m_E[plan.edge_id];
    ON_NurbsCurve *old_source_curve = extract_iso_trace(
	    dynamic_cast<const ON_NurbsSurface *>(
		brep->m_F[plan.source_face].SurfaceOf()),
	    plan.source_trace, old_edge.Domain());
    ON_NurbsCurve *old_mate_curve = extract_iso_trace(
	    dynamic_cast<const ON_NurbsSurface *>(
		brep->m_F[plan.mate_face].SurfaceOf()),
	    plan.mate_trace, old_edge.Domain());
    ON_NurbsCurve *new_source_curve = extract_iso_trace(source_surface,
	    plan.source_trace, old_edge.Domain());
    ON_NurbsCurve *new_mate_curve = extract_iso_trace(mate_surface,
	    plan.mate_trace, old_edge.Domain());
    if (!old_source_curve || !old_mate_curve || !new_source_curve
	    || !new_mate_curve
	    || !same_curve_space(old_source_curve, old_mate_curve)
	    || !same_curve_space(old_source_curve, new_source_curve)
	    || !same_curve_space(old_mate_curve, new_mate_curve)) {
	delete old_source_curve;
	delete old_mate_curve;
	delete new_source_curve;
	delete new_mate_curve;
	return false;
    }

    double max_delta_error = 0.0;
    for (int cv = 0; cv < new_source_curve->CVCount(); ++cv) {
	ON_4dPoint old_source_cv;
	ON_4dPoint old_mate_cv;
	ON_4dPoint new_source_cv;
	ON_4dPoint new_mate_cv;
	if (!old_source_curve->GetCV(cv, ON::euclidean_rational,
		&old_source_cv.x)
		|| !old_mate_curve->GetCV(cv, ON::euclidean_rational,
		    &old_mate_cv.x)
		|| !new_source_curve->GetCV(cv, ON::euclidean_rational,
		    &new_source_cv.x)
		|| !new_mate_curve->GetCV(cv, ON::euclidean_rational,
		    &new_mate_cv.x)) {
	    delete old_source_curve;
	    delete old_mate_curve;
	    delete new_source_curve;
	    delete new_mate_curve;
	    return false;
	}
	const ON_3dVector source_change(
		new_source_cv.x - old_source_cv.x,
		new_source_cv.y - old_source_cv.y,
		new_source_cv.z - old_source_cv.z);
	const ON_3dVector mate_change(
		new_mate_cv.x - old_mate_cv.x,
		new_mate_cv.y - old_mate_cv.y,
		new_mate_cv.z - old_mate_cv.z);
	max_delta_error = std::max(max_delta_error,
		(source_change - mate_change).Length());
    }
    const double edit_scale = std::max(1.0, delta.Length());
    if (max_delta_error > 1.0e-8 * edit_scale) {
	delete old_source_curve;
	delete old_mate_curve;
	delete new_source_curve;
	delete new_mate_curve;
	return false;
    }

    const ON_3dPoint old_start = old_source_curve->PointAtStart();
    const ON_3dPoint old_end = old_source_curve->PointAtEnd();
    const ON_3dPoint new_start = new_source_curve->PointAtStart();
    const ON_3dPoint new_end = new_source_curve->PointAtEnd();
    delete old_source_curve;
    delete old_mate_curve;
    delete new_mate_curve;
    if (old_start.DistanceTo(new_start) > 1.0e-8 * edit_scale
	    || old_end.DistanceTo(new_end) > 1.0e-8 * edit_scale) {
	delete new_source_curve;
	return false;
    }

    if (!replace_edge_curve(&trial, plan.edge_id, new_source_curve))
	return false;
    trial.DestroyRuntimeCache(true);
    if (!trial.IsValid(NULL))
	return false;

    *brep = trial;
    return true;
}

bool
brep_face_set_cv_constrained(ON_Brep *brep, int face_id, int cv_id_u,
	int cv_id_v, const ON_4dPoint &point)
{
    if (!brep || face_id < 0 || face_id >= brep->m_F.Count()
	    || !point.IsValid())
	return false;
    const ON_NurbsSurface *surface =
	dynamic_cast<const ON_NurbsSurface *>(brep->m_F[face_id].SurfaceOf());
    if (!surface)
	return false;
    ON_4dPoint current;
    if (!surface->GetCV(cv_id_u, cv_id_v, ON::euclidean_rational,
	    &current.x))
	return false;
    if (brep_face_cv_is_topology_safe(brep, face_id, cv_id_u, cv_id_v))
	return brep_face_set_cv(brep, face_id, cv_id_u, cv_id_v, point);
    const double weight_scale =
	std::max(1.0, std::max(std::fabs(current.w), std::fabs(point.w)));
    if (std::fabs(current.w - point.w) > 1.0e-12 * weight_scale)
	return false;
    return brep_face_translate_cv_constrained(brep, face_id, cv_id_u,
	    cv_id_v, ON_3dVector(point.x - current.x,
		point.y - current.y, point.z - current.z));
}

bool brep_face_set_weight(ON_Brep *brep, int face, int cv_id_u, int cv_id_v, double weight)
{
    if (!brep_face_cv_is_topology_safe(brep, face, cv_id_u, cv_id_v)) {
	bu_log("face %d CV (%d,%d) influences a trimming edge\n",
		face, cv_id_u, cv_id_v);
	return false;
    }
    return edit_face_nurbs_surface(brep, face,
	    [cv_id_u, cv_id_v, weight](ON_NurbsSurface *surface) {
		return ON_IsValid(weight) && weight > 0.0
		    && surface->SetWeight(cv_id_u, cv_id_v, weight);
	    });
}

bool brep_face_reverse_parameter(ON_Brep *brep, int face, int dir)
{
    if (!valid_brep(brep) || face < 0 || face >= brep->m_F.Count()
	    || dir < 0 || dir > 1)
	return false;

    ON_Brep trial(*brep);
    if (!trial.m_F[face].Reverse(dir))
	return false;
    if (brep->IsValid(NULL) && !trial.IsValid(NULL))
	return false;
    *brep = trial;
    return true;
}

bool brep_face_transpose(ON_Brep *brep, int face)
{
    if (!valid_brep(brep) || face < 0 || face >= brep->m_F.Count())
	return false;

    ON_Brep trial(*brep);
    if (!trial.m_F[face].Transpose())
	return false;
    if (brep->IsValid(NULL) && !trial.IsValid(NULL))
	return false;
    *brep = trial;
    return true;
}

ON_Curve* getEdgeCurve(const ON_Surface& s, int side)
{
    ON_2dPoint from, to;
    double u0, u1, v0, v1;
    s.GetDomain(0, &u0, &u1);
    s.GetDomain(1, &v0, &v1);

    switch (side) {
	case 0:
	    from.x = u0; from.y = v0;
	    to.x   = u1; to.y   = v0;
	    break;
	case 1:
	    from.x = u1; from.y = v0;
	    to.x   = u1; to.y   = v1;
	    break;
	case 2:
	    from.x = u1; from.y = v1;
	    to.x   = u0; to.y   = v1;
	    break;
	case 3:
	    from.x = u0; from.y = v1;
	    to.x   = u0; to.y   = v0;
	    break;
	default:
	    return NULL;
    }
    ON_Curve* c2d = new ON_LineCurve(from, to);
    c2d->SetDomain(0.0, 1.0);
    return c2d;
}

int brep_loop_create(ON_Brep *brep, int face_id)
{
    if(face_id < 0 || face_id >= brep->m_F.Count()) {
	bu_log("face_id is out of range\n");
	return -1;
    }
    ON_BrepFace& face = brep->m_F[face_id];
    ON_BrepLoop& loop = brep->NewLoop(ON_BrepLoop::outer, face);
    return loop.m_loop_index;
}

int brep_trim_create(ON_Brep *brep, int loop_id, int edge_id, int orientation, int para_curve_id)
{
    if(loop_id < 0 || loop_id >= brep->m_L.Count()) {
	bu_log("loop_id is out of range\n");
	return -1;
    }
    if(edge_id < 0 || edge_id >= brep->m_E.Count()) {
	bu_log("edge_id is out of range\n");
	return -1;
    }
    if(para_curve_id < 0 || para_curve_id >= brep->m_C2.Count()) {
	bu_log("para_curve_id is out of range\n");
	return -1;
    }
    ON_BrepLoop& loop = brep->m_L[loop_id];
    ON_BrepTrim& trim = brep->NewTrim(brep->m_E[edge_id], orientation, loop, para_curve_id);
    trim.m_type = ON_BrepTrim::mated;
    const ON_Curve* c2 = brep->m_C2[trim.m_c2i];
    const ON_Surface* s = loop.SurfaceOf();
    ON_Interval PD = trim.ProxyCurveDomain();
    trim.m_iso = s->IsIsoparametric(*c2, &PD);
    trim.m_tolerance[0] = 0.0;
    trim.m_tolerance[1] = 0.0;
    return trim.m_trim_index;
}

std::vector<ON_3dVector> calculateTangentVectors(const std::vector<ON_3dPoint> &points)
{
    int n = points.size() - 1;

    /// calculate qk
    ON_3dVectorArray qk(n + 3);
    ON_3dVector q_;
    for (size_t i = 1; i < points.size(); ++i) {
	qk[(int)i] = ON_3dVector(points[i] - points[i - 1]);
    }
    qk[0] = 2 * qk[1] - qk[2];
    q_ = 2 * qk[0] - qk[1];
    qk[n + 1] = 2 * qk[n] - qk[n - 1];
    qk[n + 2] = 2 * qk[n + 1] - qk[n];

    /// tk is a middle variable for calculating ak. tk=|qk-1 x qk|
    double *tk = (double *)bu_calloc(n + 3, sizeof(double), "tk");
    tk[0] = ON_CrossProduct(q_, qk[0]).Length();
    for (int i = 1; i < n + 3; ++i) {
	tk[i] = ON_CrossProduct(qk[i - 1], qk[i]).Length();
    }

    /// calculate ak
    double *ak = (double *)bu_calloc(n + 1, sizeof(double), "ak");
    for (int i = 0; i < n + 1; ++i) {
	ak[i] = tk[i] / (tk[i] + tk[i + 2]);
    }
    bu_free(tk, "tk");

    /// calculate vk
    ON_3dVectorArray vk(n + 1);
    for (int i = 0; i < n + 1; ++i) {
	vk[i] = (1 - ak[i]) * qk[i] + ak[i] * qk[i];
    }
    bu_free(ak, "ak");

    std::vector<ON_3dVector> tangentVectors;

    for (int i = 0; i < n + 1; ++i) {
	tangentVectors.push_back(vk[i].UnitVector());
    }
    return tangentVectors;
}

double getPosRoot(const double a, const double b,const double c)
{
    double delta = b * b - 4 * a * c;
    if (delta < 0) {
	return -1;
    }
    else if (NEAR_ZERO(delta, 1e-6)) {
	return -b / (2 * a);
    }
    else {
	double x1 = (-b + sqrt(delta)) / (2 * a);
	if (x1 > 0) {
	    return x1;
	}
	else {
	    return -1;
	}
    }

}

void calcuBsplineCVsKnots(std::vector<ON_3dPoint> &cvs, std::vector<double> &knots, const std::vector<ON_3dPoint> &endPoints, const std::vector<ON_3dVector> &tangentVectors)
{
    std::vector<double> uk;
    uk.push_back(0);
    for (size_t i = 0; i < endPoints.size() - 1; i++) {
	cvs.push_back(endPoints[i]);
	double a = 16 - (tangentVectors[i] + tangentVectors[i + 1]).LengthSquared();
	double b = 12 * (tangentVectors[i] + tangentVectors[i + 1])*(endPoints[i + 1] - endPoints[i]);
	double c = -36 * (endPoints[i + 1] - endPoints[i]).LengthSquared();
	double d = getPosRoot(a, b, c);
	ON_3dPoint p1 = endPoints[i] + tangentVectors[i] * d / 3.0;
	ON_3dPoint p2 = endPoints[i + 1] - tangentVectors[i + 1] * d / 3.0;
	cvs.push_back(p1);
	cvs.push_back(p2);
	uk.push_back(uk[i] + 3 * (p1.DistanceTo(p2)));
    }
    cvs.push_back(endPoints[endPoints.size() - 1]);


    knots.clear();
    for (size_t i = 0; i < 3; i++)
	knots.push_back(0);
    for (size_t i = 1; i < uk.size() - 1; i++) {
	knots.push_back(uk[i] / uk[uk.size() - 1]);
	knots.push_back(uk[i] / uk[uk.size() - 1]);
	knots.push_back(uk[i] / uk[uk.size() - 1]);
    }
    for (size_t i = 0; i < 3; i++)
	knots.push_back(1);
}

// calculate parameters of the surface
// input: n, m, points
// output: uk, ul
int SurfMeshParams(int n, int m, std::vector<ON_3dPoint> points, std::vector<double> &uk, std::vector<double> &ul)
{
    std::vector<double> cds;
    n > m ? cds.resize(n + 1) : cds.resize(m + 1);
    int num = m + 1; // number of nondegenerate rows
    uk.resize(n + 1);
    uk[0] = 0.0f;
    uk[n] = 1.0f;
    for (int i = 1; i < n; i++)
	uk[i] = 0;
    for (int i = 0; i <= m; i++) {
	double sum = 0;
	for (int j = 1; j <= n; j++) {
	    cds[j] = points[i * (n + 1) + j].DistanceTo(points[i * (n + 1) + j - 1]);
	    sum += cds[j];
	}
	if (sum <= 0)
	    num--;
	else
	{
	    double d = 0.0f;
	    for (int j = 1; j < n; j++)
	    {
		d += cds[j];
		uk[j] += d / sum;
	    }
	}
    }
    if (num == 0)
	return -1;
    for (int i = 1; i < n; i++)
	uk[i] /= num;

    num = n + 1;
    ul.resize(m + 1);
    ul[0] = 0.0f;
    ul[m] = 1.0f;
    for (int i = 1; i < m; i++)
	ul[i] = 0;
    for (int i = 0; i <= n; i++) {
	double sum = 0;
	for (int j = 1; j <= m; j++) {
	    cds[j] = points[j * (n + 1) + i].DistanceTo(points[(j - 1) * (n + 1) + i]);
	    sum += cds[j];
	}
	if (sum <= 0)
	    num--;
	else
	{
	    double d = 0.0f;
	    for (int j = 1; j < m; j++)
	    {
		d += cds[j];
		ul[j] += d / sum;
	    }
	}
    }
    if (num == 0)
	return -1;
    for (int i = 1; i < m; i++)
	ul[i] /= num;
    return 0;
}

void bsplineBasisFuns(int i, double u, int p, std::vector<double> U,  std::vector<double> &N)
{
    double *left = (double *)bu_calloc(p + 1, sizeof(double), "left");
    double *right = (double *)bu_calloc(p + 1, sizeof(double), "right");
    double saved, temp;
    int j, r;
    N[0] = 1.0;

    for (j = 1; j <= p; j++) {
	left[j] = u - U[i + 1 - j];
	right[j] = U[i + j] - u;
	saved = 0.0;

	for (r = 0; r < j; r++) {
	    temp = N[r] / (right[r + 1] + left[j - r]);
	    N[r] = saved + right[r + 1] * temp;
	    saved = left[j - r] * temp;
	}
	N[j] = saved;
    }
    bu_free(left, "left");
    bu_free(right, "right");
}

/**
 * solve tridiagonal equation for cubic spline interpolation
 */
int solveTridiagonalint(int n, std::vector<ON_3dPoint> Q, std::vector<double> U, std::vector<ON_3dPoint>& P)
{
    std::vector<double> abc(4);
    std::vector<ON_3dPoint> R(n + 1);
    std::vector<double> dd(n + 1);
    double den;
    int i;

    for (i = 3; i < n; i++) {
	R[i] = Q[i - 1];
    }

    bsplineBasisFuns(4, U[4], 3, U, abc);
    den = abc[1];

    /* P[2] */
    P[2] = (Q[1] - abc[0] * P[1]) / den;

    for (i = 3; i < n; i++) {
	dd[i] = abc[2] / den;

	bsplineBasisFuns(i + 2, U[i + 2], 3, U, abc);
	den = abc[1] - abc[0] * dd[i];
	P[i] = (R[i] - abc[0] * P[i - 1]) / den;
    }

    dd[n] = abc[2] / den;

    bsplineBasisFuns(n + 2, U[n + 2], 3, U, abc);
    den = abc[1] - abc[0] * dd[n];

    P[n] = (Q[n - 1] - abc[2] * P[n + 1] - abc[0] * P[n - 1]) / den;

    for (i = (n - 1); i >= 2; i--) {
	P[i] = P[i] - dd[i + 1] * P[i + 1];
    }

    if (n == 2) {
	P[2] /= 4.0f / 3.0f;
    }
    return 1;
}

void globalCurveInterp(int n, std::vector<double> &knots, const std::vector<ON_3dPoint> &Q, std::vector<ON_3dPoint> &P)
{
    /// initialize control points of P[0], P[1], P[n+1], P[n+2]
    P.resize(n + 3);
    P[0] = Q[0];
    P[1] = Q[0] + (Q[1] - Q[0]) / 3.0f * knots[4];
    P[n + 2] = Q[n];
    P[n + 1] = Q[n] - (Q[n] - Q[n - 1]) / 3.0f * (1.0f - knots[n + 2]);

    solveTridiagonalint(n, Q, knots, P);
}
// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
