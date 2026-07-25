/*                         E D I T . C P P
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
/** @file edit.cpp
 *
 * Tests for reference-safe and face-aware libbrep NURBS editing.
 */

#include "common.h"

#include <cmath>

#include "brep/edit.h"
#include "bu/app.h"
#include "bu/log.h"

static bool
near_point(const ON_3dPoint &a, const ON_3dPoint &b, double tol = 1.0e-9)
{
    return a.DistanceTo(b) <= tol;
}

static ON_4dPoint
surface_cv(const ON_Brep &brep, int face_id, int i, int j)
{
    const ON_Surface *surface = brep.m_F[face_id].SurfaceOf();
    const ON_NurbsSurface *nurbs = dynamic_cast<const ON_NurbsSurface *>(surface);
    if (!nurbs)
	bu_exit(1, "face %d does not reference a NURBS surface\n", face_id);

    ON_4dPoint cv;
    if (!nurbs->GetCV(i, j, ON::euclidean_rational, &cv.x))
	bu_exit(1, "failed to get face %d CV (%d, %d)\n", face_id, i, j);
    return cv;
}

static double
surface_greville(const ON_NurbsSurface &surface, int dir, int cv)
{
    const int degree = surface.Degree(dir);
    const int knot_count = surface.KnotCount(dir);
    double value = 0.0;
    for (int k = 1; k <= degree; ++k) {
	const int full_knot = cv + k;
	const double knot = full_knot <= 0
	    ? surface.Knot(dir, 0)
	    : (full_knot >= knot_count + 1
		? surface.Knot(dir, knot_count - 1)
		: surface.Knot(dir, full_knot - 1));
	value += knot;
    }
    return value / degree;
}

static ON_NurbsSurface *
planar_edit_surface(int cv_count, double x_offset, bool rational)
{
    ON_NurbsSurface *surface =
	ON_NurbsSurface::New(3, rational, 4, 4, cv_count, cv_count);
    surface->MakeClampedUniformKnotVector(0, 1.0);
    surface->MakeClampedUniformKnotVector(1, 1.0);
    for (int i = 0; i < cv_count; ++i) {
	for (int j = 0; j < cv_count; ++j) {
	    ON_4dPoint cv(
		    x_offset + surface_greville(*surface, 0, i),
		    surface_greville(*surface, 1, j), 0.0,
		    rational ? 1.0 + 0.1 * (j % 3) : 1.0);
	    surface->SetCV(i, j, ON::euclidean_rational, &cv.x);
	}
    }
    return surface;
}

static int
append_line_edge(ON_Brep &brep, int from, int to)
{
    ON_LineCurve *curve = new ON_LineCurve(
	    brep.m_V[from].point, brep.m_V[to].point);
    curve->SetDomain(0.0, 1.0);
    const int curve_id = brep.AddEdgeCurve(curve);
    ON_BrepEdge &edge = brep.NewEdge(
	    brep.m_V[from], brep.m_V[to], curve_id, NULL, 1.0e-8);
    return edge.m_edge_index;
}

static int
append_rect_face(ON_Brep &brep, int surface_id,
	const ON_Interval &u, const ON_Interval &v,
	const int edges[4], const bool reversed[4], int mated_side)
{
    ON_BrepFace &face = brep.NewFace(surface_id);
    ON_BrepLoop &loop = brep.NewLoop(ON_BrepLoop::outer, face);
    const ON_2dPoint starts[4] = {
	ON_2dPoint(u.Min(), v.Min()),
	ON_2dPoint(u.Max(), v.Min()),
	ON_2dPoint(u.Max(), v.Max()),
	ON_2dPoint(u.Min(), v.Max())
    };
    const ON_2dPoint ends[4] = {
	ON_2dPoint(u.Max(), v.Min()),
	ON_2dPoint(u.Max(), v.Max()),
	ON_2dPoint(u.Min(), v.Max()),
	ON_2dPoint(u.Min(), v.Min())
    };
    const ON_Surface *surface = brep.m_S[surface_id];

    for (int side = 0; side < 4; ++side) {
	ON_LineCurve *curve = new ON_LineCurve(starts[side], ends[side]);
	curve->SetDomain(0.0, 1.0);
	const int curve_id = brep.AddTrimCurve(curve);
	ON_BrepTrim &trim = brep.NewTrim(
		brep.m_E[edges[side]], reversed[side], loop, curve_id);
	trim.m_iso = surface->IsIsoparametric(*curve);
	trim.m_type = side == mated_side
	    ? ON_BrepTrim::mated : ON_BrepTrim::boundary;
	trim.m_tolerance[0] = 0.0;
	trim.m_tolerance[1] = 0.0;
    }
    return face.m_face_index;
}

static ON_Brep
two_face_iso_fixture(bool interior_iso, bool rational = false)
{
    ON_Brep brep;
    const int cv_count = interior_iso ? 9 : 7;
    ON_NurbsSurface *surface0 =
	planar_edit_surface(cv_count, 0.0, rational);
    ON_NurbsSurface *surface1 = interior_iso
	? new ON_NurbsSurface(*surface0)
	: planar_edit_surface(cv_count,
	    surface0->Domain(0).Length(), rational);
    const int surface0_id = brep.AddSurface(surface0);
    const int surface1_id = brep.AddSurface(surface1);

    const ON_Interval full_u = surface0->Domain(0);
    const ON_Interval full_v = surface0->Domain(1);
    const double split = interior_iso
	? full_u.Mid() : full_u.Max();
    const ON_Interval face0_u(full_u.Min(), split);
    const ON_Interval face1_u(interior_iso
	? split : surface1->Domain(0).Min(),
	interior_iso ? full_u.Max() : surface1->Domain(0).Max());

    const ON_3dPoint p0 = surface0->PointAt(face0_u.Min(), full_v.Min());
    const ON_3dPoint p1 = surface0->PointAt(face0_u.Max(), full_v.Min());
    const ON_3dPoint p2 = surface0->PointAt(face0_u.Max(), full_v.Max());
    const ON_3dPoint p3 = surface0->PointAt(face0_u.Min(), full_v.Max());
    const ON_3dPoint p4 = surface1->PointAt(face1_u.Max(), full_v.Min());
    const ON_3dPoint p5 = surface1->PointAt(face1_u.Max(), full_v.Max());
    const int vertices[6] = {
	brep.NewVertex(p0, 1.0e-8).m_vertex_index,
	brep.NewVertex(p1, 1.0e-8).m_vertex_index,
	brep.NewVertex(p2, 1.0e-8).m_vertex_index,
	brep.NewVertex(p3, 1.0e-8).m_vertex_index,
	brep.NewVertex(p4, 1.0e-8).m_vertex_index,
	brep.NewVertex(p5, 1.0e-8).m_vertex_index
    };

    const int edges[7] = {
	append_line_edge(brep, vertices[0], vertices[1]),
	append_line_edge(brep, vertices[1], vertices[2]),
	append_line_edge(brep, vertices[2], vertices[3]),
	append_line_edge(brep, vertices[3], vertices[0]),
	append_line_edge(brep, vertices[1], vertices[4]),
	append_line_edge(brep, vertices[4], vertices[5]),
	append_line_edge(brep, vertices[5], vertices[2])
    };
    const int face0_edges[4] = {edges[0], edges[1], edges[2], edges[3]};
    const bool face0_reversed[4] = {false, false, false, false};
    const int face1_edges[4] = {edges[4], edges[5], edges[6], edges[1]};
    const bool face1_reversed[4] = {false, false, false, true};
    append_rect_face(brep, surface0_id, face0_u, full_v,
	    face0_edges, face0_reversed, 1);
    append_rect_face(brep, surface1_id, face1_u, full_v,
	    face1_edges, face1_reversed, 3);
    return brep;
}

static void
test_vertex_reference_safety()
{
    ON_Brep brep;
    int c3 = brep_curve_make(&brep, ON_3dPoint::Origin);
    int v0 = brep_vertex_create(&brep, ON_3dPoint(0.0, 0.0, 0.0));
    int unused = brep_vertex_create(&brep, ON_3dPoint(5.0, 5.0, 5.0));
    int v1 = brep_vertex_create(&brep, ON_3dPoint(1.0, 0.0, 0.0));
    int edge = brep_edge_create(&brep, v0, v1, c3);

    if (brep_vertex_remove(&brep, v0))
	bu_exit(1, "removed a referenced vertex\n");
    if (!brep_vertex_remove(&brep, unused))
	bu_exit(1, "failed to remove an unreferenced vertex\n");
    if (brep.m_E[edge].m_vi[0] != 0 || brep.m_E[edge].m_vi[1] != 1)
	bu_exit(1, "vertex removal did not remap edge vertex indices\n");
}

static void
test_curve_reference_safety()
{
    ON_Brep brep;
    int c0 = brep_curve_make(&brep, ON_3dPoint(0.0, 0.0, 0.0));
    int unused = brep_curve_make(&brep, ON_3dPoint(1.0, 0.0, 0.0));
    int c2 = brep_curve_make(&brep, ON_3dPoint(2.0, 0.0, 0.0));
    int v0 = brep_vertex_create(&brep, ON_3dPoint(0.0, 0.0, 0.0));
    int v1 = brep_vertex_create(&brep, ON_3dPoint(1.0, 0.0, 0.0));
    int edge = brep_edge_create(&brep, v0, v1, c2);
    const ON_Curve *edge_curve = brep.m_E[edge].ProxyCurve();

    if (!brep_curve_remove(&brep, unused))
	bu_exit(1, "failed to remove an unreferenced 3D curve\n");
    if (brep.m_E[edge].m_c3i != 1 || brep.m_E[edge].ProxyCurve() != edge_curve)
	bu_exit(1, "3D curve removal did not preserve the edge reference\n");
    if (brep_curve_remove(&brep, 1))
	bu_exit(1, "removed a referenced 3D curve\n");

    brep.m_E[edge].ON_CurveProxy::Reverse();
    const bool edge_reversed = brep.m_E[edge].ProxyCurveIsReversed();
    const ON_Interval edge_domain = brep.m_E[edge].Domain();
    const ON_3dPoint edge_midpoint =
	brep.m_E[edge].PointAt(edge_domain.Mid());
    const ON_NurbsCurve *referenced_curve =
	dynamic_cast<const ON_NurbsCurve *>(brep.m_C3[1]);
    if (!referenced_curve
	    || !brep_curve_insert_knot(&brep, 1,
		referenced_curve->Domain().Mid(), 1))
	bu_exit(1, "shape-preserving referenced curve edit failed\n");
    if (brep.m_E[edge].ProxyCurveIsReversed() != edge_reversed
	    || brep.m_E[edge].Domain() != edge_domain
	    || !near_point(brep.m_E[edge].PointAt(edge_domain.Mid()),
		edge_midpoint, 1.0e-8))
	bu_exit(1, "referenced curve edit changed edge proxy semantics\n");

    const ON_NurbsCurve *curve =
	dynamic_cast<const ON_NurbsCurve *>(brep.m_C3[c0]);
    if (!curve)
	bu_exit(1, "test curve is not NURBS\n");

    ON_Interval domain = curve->Domain();
    ON_3dPoint samples_before[5];
    for (int i = 0; i < 5; ++i)
	samples_before[i] = curve->PointAt(domain.ParameterAt(0.25 * i));

    if (!brep_curve_insert_knot(&brep, c0, domain.Mid(), 1))
	bu_exit(1, "curve knot insertion failed\n");
    curve = dynamic_cast<const ON_NurbsCurve *>(brep.m_C3[c0]);
    if (!brep_curve_elevate_degree(&brep, c0, curve->Degree() + 1))
	bu_exit(1, "curve degree elevation failed\n");
    curve = dynamic_cast<const ON_NurbsCurve *>(brep.m_C3[c0]);
    for (int i = 0; i < 5; ++i) {
	ON_3dPoint sample = curve->PointAt(domain.ParameterAt(0.25 * i));
	if (!near_point(samples_before[i], sample, 1.0e-8))
	    bu_exit(1, "shape-preserving curve edit changed the curve locus\n");
    }
}

static void
test_surface_reference_safety()
{
    ON_Brep brep;
    int s0 = brep_surface_make(&brep, ON_3dPoint(0.0, 0.0, 0.0));
    int unused = brep_surface_make(&brep, ON_3dPoint(0.0, 0.0, 1.0));
    int s2 = brep_surface_make(&brep, ON_3dPoint(0.0, 0.0, 2.0));
    int face = brep_face_create(&brep, s2);
    const ON_Surface *face_surface = brep.m_F[face].SurfaceOf();

    if (!brep_surface_remove(&brep, unused))
	bu_exit(1, "failed to remove an unreferenced surface\n");
    if (brep.m_F[face].m_si != 1 || brep.m_F[face].SurfaceOf() != face_surface)
	bu_exit(1, "surface removal did not preserve the face reference\n");
    if (brep_surface_remove(&brep, 1))
	bu_exit(1, "removed a referenced surface\n");

    brep.m_F[face].ON_SurfaceProxy::Transpose();
    const bool face_transposed = brep.m_F[face].ProxySurfaceIsTransposed();
    const ON_3dPoint face_midpoint = brep.m_F[face].PointAt(
	    brep.m_F[face].Domain(0).Mid(), brep.m_F[face].Domain(1).Mid());
    const ON_NurbsSurface *referenced_surface =
	dynamic_cast<const ON_NurbsSurface *>(brep.m_S[1]);
    if (!referenced_surface
	    || !brep_surface_insert_knot(&brep, 1, 0,
		referenced_surface->Domain(0).Mid(), 1))
	bu_exit(1, "shape-preserving referenced surface edit failed\n");
    if (brep.m_F[face].ProxySurfaceIsTransposed() != face_transposed
	    || !near_point(brep.m_F[face].PointAt(
		    brep.m_F[face].Domain(0).Mid(),
		    brep.m_F[face].Domain(1).Mid()),
		face_midpoint, 1.0e-8))
	bu_exit(1, "referenced surface edit changed face proxy semantics\n");

    const ON_NurbsSurface *surface =
	dynamic_cast<const ON_NurbsSurface *>(brep.m_S[s0]);
    if (!surface)
	bu_exit(1, "test surface is not NURBS\n");

    ON_Interval udom = surface->Domain(0);
    ON_Interval vdom = surface->Domain(1);
    ON_3dPoint samples_before[3][3];
    for (int i = 0; i < 3; ++i) {
	for (int j = 0; j < 3; ++j) {
	    samples_before[i][j] = surface->PointAt(
		    udom.ParameterAt(0.5 * i), vdom.ParameterAt(0.5 * j));
	}
    }

    if (!brep_surface_insert_knot(&brep, s0, 0, udom.Mid(), 1))
	bu_exit(1, "surface knot insertion failed\n");
    surface = dynamic_cast<const ON_NurbsSurface *>(brep.m_S[s0]);
    if (!brep_surface_elevate_degree(&brep, s0, 1, surface->Degree(1) + 1))
	bu_exit(1, "surface degree elevation failed\n");
    surface = dynamic_cast<const ON_NurbsSurface *>(brep.m_S[s0]);
    for (int i = 0; i < 3; ++i) {
	for (int j = 0; j < 3; ++j) {
	    ON_3dPoint sample = surface->PointAt(
		    udom.ParameterAt(0.5 * i), vdom.ParameterAt(0.5 * j));
	    if (!near_point(samples_before[i][j], sample, 1.0e-8))
		bu_exit(1, "shape-preserving surface edit changed the surface locus\n");
	}
    }
}

static void
test_shared_surface_isolation()
{
    ON_Brep brep;
    int surface = brep_surface_make(&brep, ON_3dPoint::Origin);
    int face0 = brep_face_create(&brep, surface);
    int face1 = brep_face_create(&brep, surface);
    ON_4dPoint before0 = surface_cv(brep, face0, 1, 1);
    ON_4dPoint before1 = surface_cv(brep, face1, 1, 1);

    if (!brep_face_translate_cv(&brep, face0, 1, 1, ON_3dVector(1.0, 2.0, 3.0)))
	bu_exit(1, "face-aware CV translation failed\n");
    if (brep.m_F[face0].m_si == brep.m_F[face1].m_si)
	bu_exit(1, "face-aware CV translation did not isolate a shared surface\n");

    ON_4dPoint after0 = surface_cv(brep, face0, 1, 1);
    ON_4dPoint after1 = surface_cv(brep, face1, 1, 1);
    if (std::fabs(after0.x - before0.x - 1.0) > 1.0e-9
	    || std::fabs(after0.y - before0.y - 2.0) > 1.0e-9
	    || std::fabs(after0.z - before0.z - 3.0) > 1.0e-9
	    || std::fabs(after0.w - before0.w) > 1.0e-9)
	bu_exit(1, "face-aware CV translation produced the wrong CV\n");
    if (std::fabs(after1.x - before1.x) > 1.0e-9
	    || std::fabs(after1.y - before1.y) > 1.0e-9
	    || std::fabs(after1.z - before1.z) > 1.0e-9
	    || std::fabs(after1.w - before1.w) > 1.0e-9)
	bu_exit(1, "face-aware CV translation changed another face\n");

    if (!brep_face_set_weight(&brep, face0, 1, 1, 2.0))
	bu_exit(1, "face-aware weight edit failed\n");
    after0 = surface_cv(brep, face0, 1, 1);
    after1 = surface_cv(brep, face1, 1, 1);
    if (std::fabs(after0.w - 2.0) > 1.0e-9
	    || std::fabs(after1.w - before1.w) > 1.0e-9)
	bu_exit(1, "face-aware weight edit did not remain isolated\n");
}

static void
test_trim_support_lock()
{
    ON_Brep brep;
    ON_NurbsSurface surface(3, false, 4, 4, 7, 7);
    for (int i = 0; i < 7; ++i) {
	for (int j = 0; j < 7; ++j)
	    surface.SetCV(i, j, ON_3dPoint((double)i, (double)j, 0.0));
    }
    surface.MakeClampedUniformKnotVector(0, 1.0);
    surface.MakeClampedUniformKnotVector(1, 1.0);
    if (!brep.NewFace(surface) || !brep.IsValid(NULL))
	bu_exit(1, "failed to create a valid trimmed NURBS face\n");

    if (brep_face_cv_is_topology_safe(&brep, 0, 0, 0))
	bu_exit(1, "corner CV was incorrectly classified topology-safe\n");
    if (!brep_face_cv_is_topology_safe(&brep, 0, 3, 3))
	bu_exit(1, "interior CV was incorrectly classified trim-influencing\n");

    ON_4dPoint before = surface_cv(brep, 0, 0, 0);
    if (brep_face_translate_cv(&brep, 0, 0, 0,
		ON_3dVector(0.0, 0.0, 1.0)))
	bu_exit(1, "trim-influencing corner CV edit was accepted\n");
    ON_4dPoint after = surface_cv(brep, 0, 0, 0);
    if (std::fabs(after.x - before.x) > 1.0e-12
	    || std::fabs(after.y - before.y) > 1.0e-12
	    || std::fabs(after.z - before.z) > 1.0e-12
	    || std::fabs(after.w - before.w) > 1.0e-12)
	bu_exit(1, "rejected trim-influencing edit changed the B-rep\n");
}

static void
test_coupled_isoparametric_edit(bool interior_iso, bool rational = false)
{
    ON_Brep brep = two_face_iso_fixture(interior_iso, rational);
    if (!brep.IsValid(NULL))
	bu_exit(1, "failed to create a valid two-face isoparametric fixture\n");

    const int shared_edge = 1;
    const int expected_edge_class = interior_iso
	? BREP_EDGE_CONSTRAINT_ISOPARAMETRIC
	: BREP_EDGE_CONSTRAINT_NATURAL;
    if (brep_edge_constraint_type(&brep, shared_edge)
	    != expected_edge_class)
	bu_exit(1, "shared edge received the wrong constraint class\n");

    const ON_NurbsSurface *source =
	dynamic_cast<const ON_NurbsSurface *>(brep.m_F[0].SurfaceOf());
    const int cv_u = interior_iso ? 4 : source->CVCount(0) - 1;
    const int cv_v = source->CVCount(1) / 2;
    struct brep_face_cv_constraint status;
    if (!brep_face_cv_constraint_status(&brep, 0, cv_u, cv_v, &status)
	    || status.topology_safe || !status.can_translate
	    || status.trim_count != 1)
	bu_exit(1, "eligible isoparametric CV was not classified editable\n");

    const int vertex_count = brep.m_V.Count();
    const int edge_count = brep.m_E.Count();
    const int trim_count = brep.m_T.Count();
    const int face_count = brep.m_F.Count();
    const ON_3dPoint edge_start = brep.m_E[shared_edge].PointAtStart();
    const ON_3dPoint edge_end = brep.m_E[shared_edge].PointAtEnd();
    if (!brep_face_translate_cv_constrained(&brep, 0, cv_u, cv_v,
	    ON_3dVector(0.0, 0.0, 0.5)))
	bu_exit(1, "eligible isoparametric C0 edit was rejected\n");
    if (!brep.IsValid(NULL))
	bu_exit(1, "coupled isoparametric edit invalidated the B-rep\n");
    if (brep.m_V.Count() != vertex_count || brep.m_E.Count() != edge_count
	    || brep.m_T.Count() != trim_count || brep.m_F.Count() != face_count)
	bu_exit(1, "coupled edit changed B-rep topology\n");
    if (!near_point(edge_start, brep.m_E[shared_edge].PointAtStart(), 1.0e-8)
	    || !near_point(edge_end, brep.m_E[shared_edge].PointAtEnd(), 1.0e-8))
	bu_exit(1, "non-corner coupled edit moved an edge endpoint\n");

    const ON_BrepTrim &source_trim = brep.m_T[
	    brep.m_E[shared_edge].m_ti[0]];
    const ON_BrepTrim &mate_trim = brep.m_T[
	    brep.m_E[shared_edge].m_ti[1]];
    for (int sample = 0; sample <= 20; ++sample) {
	const double f = (double)sample / 20.0;
	const ON_3dPoint source_uv = source_trim.PointAt(
		source_trim.Domain().ParameterAt(f));
	const ON_3dPoint mate_uv = mate_trim.PointAt(
		mate_trim.Domain().ParameterAt(1.0 - f));
	const ON_3dPoint source_point =
	    brep.m_F[0].PointAt(source_uv.x, source_uv.y);
	const ON_3dPoint mate_point =
	    brep.m_F[1].PointAt(mate_uv.x, mate_uv.y);
	if (!near_point(source_point, mate_point, 1.0e-7))
	    bu_exit(1, "coupled faces do not share the edited trace\n");
    }

    ON_Brep rejected(brep);
    const ON_NurbsSurface *rejected_surface =
	dynamic_cast<const ON_NurbsSurface *>(rejected.m_F[0].SurfaceOf());
    const int corner_u = interior_iso ? cv_u : rejected_surface->CVCount(0) - 1;
    if (brep_face_translate_cv_constrained(&rejected, 0, corner_u, 0,
	    ON_3dVector(0.0, 0.0, 0.25)))
	bu_exit(1, "multi-trim/endpoint coupled edit was accepted\n");
    if (!near_point(rejected.m_E[shared_edge].PointAtStart(),
	    brep.m_E[shared_edge].PointAtStart(), 1.0e-12))
	bu_exit(1, "rejected coupled edit changed the B-rep\n");
}

static void
test_edge_constraint_classification()
{
    ON_Brep brep = two_face_iso_fixture(false);
    const int shared_edge = 1;
    if (brep_edge_constraint_type(&brep, 0)
	    != BREP_EDGE_CONSTRAINT_OPEN
	    || brep_edge_constraint_type(&brep, shared_edge)
	    != BREP_EDGE_CONSTRAINT_NATURAL)
	bu_exit(1, "open or natural edge classification failed\n");

    ON_BrepEdge &edge = brep.m_E[shared_edge];
    brep.m_T[edge.m_ti[0]].m_iso = ON_Surface::x_iso;
    brep.m_T[edge.m_ti[1]].m_iso = ON_Surface::not_iso;
    if (brep_edge_constraint_type(&brep, shared_edge)
	    != BREP_EDGE_CONSTRAINT_ONE_ISOPARAMETRIC)
	bu_exit(1, "one-isoparametric edge classification failed\n");

    brep.m_T[edge.m_ti[0]].m_iso = ON_Surface::not_iso;
    if (brep_edge_constraint_type(&brep, shared_edge)
	    != BREP_EDGE_CONSTRAINT_GENERAL)
	bu_exit(1, "general edge classification failed\n");

    brep.m_T[edge.m_ti[0]].m_type = ON_BrepTrim::boundary;
    if (brep_edge_constraint_type(&brep, shared_edge)
	    != BREP_EDGE_CONSTRAINT_SPECIAL)
	bu_exit(1, "special edge classification failed\n");
}

static void
test_revolution_preserves_source_curve()
{
    ON_Brep brep;
    const int curve_id = brep_curve_make(&brep, ON_3dPoint::Origin);
    const ON_NurbsCurve *curve =
	dynamic_cast<const ON_NurbsCurve *>(brep.m_C3[curve_id]);
    if (!curve)
	bu_exit(1, "revolution source curve is not NURBS\n");
    const ON_Interval domain = curve->Domain();
    const ON_3dPoint midpoint = curve->PointAt(domain.Mid());

    const int surface_id = brep_surface_revolution(&brep, curve_id,
	    ON_3dPoint(0.0, 0.0, -1.0),
	    ON_3dPoint(0.0, 0.0, 1.0), ON_PI);
    if (surface_id < 0 || surface_id >= brep.m_S.Count()
	    || !brep.m_S[surface_id])
	bu_exit(1, "surface revolution failed\n");

    curve = dynamic_cast<const ON_NurbsCurve *>(brep.m_C3[curve_id]);
    if (!curve || !near_point(curve->PointAt(domain.Mid()), midpoint))
	bu_exit(1, "surface revolution consumed or changed its source curve\n");
}

static void
test_face_parameter_edits_preserve_validity()
{
    ON_Brep brep;
    ON_NurbsSurface surface(3, false, 4, 4, 7, 7);
    for (int i = 0; i < 7; ++i) {
	for (int j = 0; j < 7; ++j)
	    surface.SetCV(i, j,
		    ON_3dPoint((double)i, (double)j, 0.1 * i * j));
    }
    surface.MakeClampedUniformKnotVector(0, 1.0);
    surface.MakeClampedUniformKnotVector(1, 1.0);
    if (!brep.NewFace(surface) || !brep.IsValid(NULL))
	bu_exit(1, "failed to create face parameter-edit fixture\n");

    const ON_Interval udom = brep.m_F[0].Domain(0);
    const ON_Interval vdom = brep.m_F[0].Domain(1);
    const double u = udom.ParameterAt(0.25);
    const double v = vdom.ParameterAt(0.75);
    const ON_3dPoint before = brep.m_F[0].PointAt(u, v);

    if (!brep_face_reverse_parameter(&brep, 0, 0)
	    || !brep.IsValid(NULL))
	bu_exit(1, "face parameter reversal failed or invalidated the B-rep\n");
    if (!near_point(before, brep.m_F[0].PointAt(-u, v), 1.0e-8))
	bu_exit(1, "face parameter reversal changed its locus mapping\n");

    if (!brep_face_transpose(&brep, 0) || !brep.IsValid(NULL))
	bu_exit(1, "face transposition failed or invalidated the B-rep\n");
    if (!near_point(before, brep.m_F[0].PointAt(v, -u), 1.0e-8))
	bu_exit(1, "face transposition changed its locus mapping\n");
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    if (argc != 1)
	bu_exit(1, "Usage: %s\n", argv[0]);

    test_vertex_reference_safety();
    test_curve_reference_safety();
    test_surface_reference_safety();
    test_shared_surface_isolation();
    test_trim_support_lock();
    test_coupled_isoparametric_edit(false);
    test_coupled_isoparametric_edit(true);
    test_coupled_isoparametric_edit(false, true);
    test_edge_constraint_classification();
    test_revolution_preserves_source_curve();
    test_face_parameter_edits_preserve_validity();
    bu_log("libbrep edit tests passed\n");
    return 0;
}
