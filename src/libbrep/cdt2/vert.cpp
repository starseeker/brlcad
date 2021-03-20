/*                        C D T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2007-2021 United States Government as represented by
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
/** @addtogroup libbrep */
/** @{ */
/** @file cdt.cpp
 *
 * Constrained Delaunay Triangulation of NURBS B-Rep objects.
 *
 */

#include "common.h"
#include "bu/malloc.h"
#include "bu/vls.h"
#include "brep/pullback.h"
#include "./cdt.h"

static ON_3dVector
vert_trim_vnorm(ON_BrepVertex& v, ON_BrepTrim *trim)
{
    ON_3dPoint t1, t2;
    ON_3dVector v1 = ON_3dVector::UnsetVector;
    ON_3dVector v2 = ON_3dVector::UnsetVector;
    ON_3dVector trim_norm = ON_3dVector::UnsetVector;

    ON_Interval trange = trim->Domain();
    ON_3dPoint t_2d1 = trim->PointAt(trange[0]);
    ON_3dPoint t_2d2 = trim->PointAt(trange[1]);

    ON_Plane fplane;
    const ON_Surface *s = trim->SurfaceOf();
    double ptol = s->BoundingBox().Diagonal().Length()*0.001;
    ptol = (ptol < BREP_PLANAR_TOL) ? ptol : BREP_PLANAR_TOL;
    if (s->IsPlanar(&fplane, ptol)) {
        trim_norm = fplane.Normal();
        if (trim->Face()->m_bRev) {
            trim_norm = trim_norm * -1;
        }
    } else {
        int ev1 = 0;
        int ev2 = 0;
        if (surface_EvNormal(s, t_2d1.x, t_2d1.y, t1, v1)) {
            if (trim->Face()->m_bRev) {
                v1 = v1 * -1;
            }
            ev1 = 1;
        }
        if (surface_EvNormal(s, t_2d2.x, t_2d2.y, t2, v2)) {
            if (trim->Face()->m_bRev) {
                v2 = v2 * -1;
            }
            ev2 = 1;
        }
        // If we got both of them, go with the closest one
        if (ev1 && ev2) {
            trim_norm = (v.Point().DistanceTo(t1) < v.Point().DistanceTo(t2)) ? v1 : v2;
        }

        if (ev1 && !ev2) {
            trim_norm = v1;
        }

        if (!ev1 && ev2) {
            trim_norm = v2;
        }
    }

    return trim_norm;
}

ON_3dVector
singular_vert_norm(ON_Brep *brep, int index)
{
    ON_BrepVertex &v = brep->m_V[index];
    ON_3dVector vnrml = ON_3dVector::UnsetVector;
    bool have_calculated = false;
    for (int eind = 0; eind != v.EdgeCount(); eind++) {
	ON_3dVector trim1_norm = ON_3dVector::UnsetVector;
	ON_3dVector trim2_norm = ON_3dVector::UnsetVector;
	ON_BrepEdge& edge = brep->m_E[v.m_ei[eind]];
	if (edge.TrimCount() != 2) {
	    // Don't know what to do with this yet... skip.
	    continue;
	}
	ON_BrepTrim *trim1 = edge.Trim(0);
	ON_BrepTrim *trim2 = edge.Trim(1);

	if (trim1->m_type != ON_BrepTrim::singular) {
	    trim1_norm = vert_trim_vnorm(v, trim1);
	}
	if (trim2->m_type != ON_BrepTrim::singular) {
	    trim2_norm = vert_trim_vnorm(v, trim2);
	}

	// If one of the normals is unset and the other comes from a plane, use it
	if (trim1_norm == ON_3dVector::UnsetVector && trim2_norm != ON_3dVector::UnsetVector) {
	    const ON_Surface *s2 = trim2->SurfaceOf();
	    if (!s2->IsPlanar(NULL, ON_ZERO_TOLERANCE)) {
		continue;
	    }
	    trim1_norm = trim2_norm;
	}
	if (trim1_norm != ON_3dVector::UnsetVector && trim2_norm == ON_3dVector::UnsetVector) {
	    const ON_Surface *s1 = trim1->SurfaceOf();
	    if (!s1->IsPlanar(NULL, ON_ZERO_TOLERANCE)) {
		continue;
	    }
	    trim2_norm = trim1_norm;
	}

	// If we have disagreeing normals and one of them is from a planar surface, go
	// with that one
	if (NEAR_EQUAL(ON_DotProduct(trim1_norm, trim2_norm), -1, VUNITIZE_TOL)) {
	    const ON_Surface *s1 = trim1->SurfaceOf();
	    const ON_Surface *s2 = trim2->SurfaceOf();
	    if (!s1->IsPlanar(NULL, ON_ZERO_TOLERANCE) && !s2->IsPlanar(NULL, ON_ZERO_TOLERANCE)) {
		// Normals severely disagree, no planar surface to fall back on - can't use this
		continue;
	    }
	    if (s1->IsPlanar(NULL, ON_ZERO_TOLERANCE) && s2->IsPlanar(NULL, ON_ZERO_TOLERANCE)) {
		// Two disagreeing planes - can't use this
		continue;
	    }
	    if (s1->IsPlanar(NULL, ON_ZERO_TOLERANCE)) {
		trim2_norm = trim1_norm;
	    }
	    if (s2->IsPlanar(NULL, ON_ZERO_TOLERANCE)) {
		trim1_norm = trim2_norm;
	    }
	}

	// Add the normals to the vnrml total
	vnrml += trim1_norm;
	vnrml += trim2_norm;
	have_calculated = 1;
    }

    if (!have_calculated) {
	return ON_3dVector::UnsetVector;
    }

    // Average all the successfully calculated normals into a new unit normal
    vnrml.Unitize();

    return vnrml;
}

ON_3dVector
mesh_point_t::norm(int f_ind)
{
    // No point type info, no normal - return default (unset) normal
    if (type == B_PT || !cdt)
       	return ON_3dVector::UnsetVector;

    // Surface point normals will be calculated unambiguously when B_SURF is
    // assigned - we can just return the normal for them.
    if (type == B_SURF)
       	return n_surf;

    // For other types of points, we need more context
    if (f_ind == -1)
       	return ON_3dVector::UnsetVector;

    // Check for pre-calculated stashed values - if we don't have
    // any, we need to calculate them.
    f_ind = (singular) ? -1 : f_ind;
    std::unordered_map<int, ON_3dVector>::iterator n_it;
    n_it = fnorms.find(f_ind);
    if (n_it != fnorms.end()){
	ON_3dVector vn = n_it->second;
	return vn;
    }

    // Don't have it, need to calculate it.
    ON_Brep *brep = cdt->i->s.brep;

    if (type == B_VERT) {
	// First question - is it a singularity?
	// TODO - data type wise, singularities are only singular per face -
	// should we attempt to treat them that way for normal calculations?
	// How reliable would that be?
	if (singular) {
	    ON_3dVector vn = singular_vert_norm(brep, vert_index);
	    fnorms[f_ind] = vn;
	    return vn;
	}

	// If not singular, the face index supplied as context
	// tells us which normal to calculate
	const ON_Surface *s = brep->m_F[f_ind].SurfaceOf();
	double ptol = s->BoundingBox().Diagonal().Length()*0.001;
	ON_Plane fplane;
	ptol = (ptol < BREP_PLANAR_TOL) ? ptol : BREP_PLANAR_TOL;

	if (s->IsPlanar(&fplane, ptol)) {
	    // Planar case is easy - just use the plane normal
	    ON_3dVector vn = fplane.Normal();
	    fnorms[f_ind] = vn;
	    return vn;

	} else {
	    // Need to evaluate at the 2D surface point associated with this
	    // vertex.  To find that, check the trims of the associated edges
	    // to find a trim with the same face index as ind, and get the
	    // 2D point from the trim
	    ON_BrepTrim *trim = NULL;
	    for (int i = 0; i < brep->m_V[vert_index].m_ei.Count(); i++) {
		ON_BrepEdge &edge = brep->m_E[brep->m_V[vert_index].m_ei[i]];
		if (edge.Trim(0)->FaceIndexOf() == f_ind) {
		    trim = edge.Trim(0);
		}
		if (edge.Trim(1)->FaceIndexOf() == f_ind) {
		    trim = edge.Trim(1);
		}
		if (trim)
		    break;
	    }
	    // If we couldn't find a trim, there was an error
	    if (!trim) {
		std::cerr << "No trim associated with vertex and face " << f_ind << "\n";
		return ON_3dVector::UnsetVector;
	    }

	    // Found a trim associated with the user specified face - now,
	    // determine which end is the vertex end and get the 2D surface
	    // point.
	    int tind = (trim->Vertex(0)->m_vertex_index == vert_index) ? 0 : 1;
	    ON_3dPoint trimpt = trim->PointAt(trim->Domain().ParameterAt(tind));
	    ON_3dPoint tmp1;
	    ON_3dVector vn;
	    surface_EvNormal(trim->SurfaceOf(), trimpt.x, trimpt.y, tmp1, vn);
	    fnorms[f_ind] = vn;
	    return vn;
	}
    }

    if (type == B_EDGE) {
	// If we're an edge point, we have some additional context from the edge
	// association
	int tind = -1;
	if (trim_index[0] != -1 && brep->m_T[trim_index[0]].FaceIndexOf() == f_ind) {
	    tind = 0;
	}
	if (trim_index[1] != -1 && brep->m_T[trim_index[1]].FaceIndexOf() == f_ind) {
	    tind = 1;
	}
	if (tind == -1) {
	    std::cerr << "No edge trim for this point associated with face " << f_ind << "\n";
	    return ON_3dVector::UnsetVector;
	}

	const ON_Surface *s = brep->m_T[tind].SurfaceOf();
	double ptol = s->BoundingBox().Diagonal().Length()*0.001;
	ON_Plane fplane;
	ptol = (ptol < BREP_PLANAR_TOL) ? ptol : BREP_PLANAR_TOL;
	if (s->IsPlanar(&fplane, ptol)) {
	    // Planar case is easy - just use the plane normal
	    ON_3dVector vn = fplane.Normal();
	    fnorms[f_ind] = vn;
	    return vn;
	}

	// Non-planar - do the calculation
	ON_3dPoint trimpt = brep->m_T[tind].PointAt(trim_t[tind]);
	ON_3dPoint tmp1;
	ON_3dVector vn;
	surface_EvNormal(brep->m_T[tind].SurfaceOf(), trimpt.x, trimpt.y, tmp1, vn);
	fnorms[f_ind] = vn;
	return vn;
    }

    return ON_3dVector::UnsetVector;
}

void
mesh_point_t::update()
{
    // Before updating, remove old box (if any) from RTree
    double p1[3];
    p1[0] = bb.Min().x - 2*ON_ZERO_TOLERANCE;
    p1[1] = bb.Min().y - 2*ON_ZERO_TOLERANCE;
    p1[2] = bb.Min().z - 2*ON_ZERO_TOLERANCE;
    double p2[3];
    p2[0] = bb.Max().x + 2*ON_ZERO_TOLERANCE;
    p2[1] = bb.Max().y + 2*ON_ZERO_TOLERANCE;
    p2[2] = bb.Max().z + 2*ON_ZERO_TOLERANCE;
    cdt->i->s.b_pnts_tree.Remove(p1, p2, vect_ind);

    // Start with a minimal box around the point - anything else is up to
    // the edge lengths
    ON_3dPoint pztol(ON_ZERO_TOLERANCE, ON_ZERO_TOLERANCE, ON_ZERO_TOLERANCE);
    bb = ON_BoundingBox(p,p);
    bb.m_max = bb.m_max + pztol;
    bb.m_min = bb.m_min - pztol;

    // With no edges assigned this is a no-op, but otherwise bump the bbox
    // dimensions using half the smallest connected edge length.
    double slen = DBL_MAX;
    std::set<mesh_uedge_t *>::iterator ue_it;
    for (ue_it = uedges.begin(); ue_it != uedges.end(); ue_it++) {
	if ((*ue_it)->type == B_SINGULAR) continue;
	slen = (slen > (*ue_it)->len) ? (*ue_it)->len : slen;
    }
    slen = (slen < DBL_MAX) ? slen : 0;

    ON_3dPoint nbbp;
    ON_3dPoint bdelta(slen*0.5, slen*0.5, slen*0.5);
    nbbp = p + bdelta;
    bb.Set(nbbp, true);
    nbbp = p - bdelta;
    bb.Set(nbbp, true);

    // Reinsert vert into RTree
    p1[0] = bb.Min().x;
    p1[1] = bb.Min().y;
    p1[2] = bb.Min().z;
    p2[0] = bb.Max().x;
    p2[1] = bb.Max().y;
    p2[2] = bb.Max().z;
    cdt->i->s.b_pnts_tree.Insert(p1, p2, vect_ind);
}

/** @} */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8

