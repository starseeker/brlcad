/*                         B R E P _ F A S T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by the
 * U.S. Army Research Laboratory.
 *
 * Distributed under the terms of the GNU Lesser General Public License
 * (LGPL), version 2.1.
 */

#include "common.h"

#include "bu/app.h"
#include "bu/malloc.h"
#include <cstring>
#include <vector>

#include "brep/cdt.h"
#include "cdt/test_api.h"
#include "raytrace.h"
#include "rt/geom.h"
#include "rt/primitives/brep.h"

struct fast_output {
    std::vector<int> faces;
    std::vector<fastf_t> normals;
    std::vector<fastf_t> points;
    struct brep_cdt_fast_report report = {};
};

static int
run_fast(fast_output *output, const ON_Brep *brep, size_t workers,
	size_t max_points, size_t max_working_bytes = 0,
	size_t max_triangles = 0, bool adaptive_quality = true, int face_index = -1)
{
    int *faces = NULL;
    int face_count = 0;
    vect_t *normals = NULL;
    point_t *points = NULL;
    int point_count = 0;
    struct bg_tess_tol ttol = BG_TESS_TOL_INIT_TOL;
    struct bn_tol tol = BN_TOL_INIT_TOL;
    struct brep_cdt_fast_options options;
    brep_cdt_fast_options_default(&options);
    options.max_workers = workers;
    options.max_points = max_points;
    if (max_working_bytes)
	options.max_working_bytes = max_working_bytes;
    if (max_triangles)
	options.max_triangles = max_triangles;
    options.adaptive_quality = adaptive_quality ? 1 : 0;

    int ret = brep_cdt_fast_ex(&faces, &face_count, &normals, &points,
	&point_count, brep, face_index, &ttol, &tol, &options, &output->report);
    if (faces)
	output->faces.assign(faces, faces + (size_t)face_count * 3);
    if (normals) {
	const fastf_t *values = (const fastf_t *)normals;
	output->normals.assign(values, values + (size_t)face_count * 9);
    }
    if (points) {
	const fastf_t *values = (const fastf_t *)points;
	output->points.assign(values, values + (size_t)point_count * 3);
    }

    bu_free(faces, "fast test faces");
    bu_free(normals, "fast test normals");
    bu_free(points, "fast test points");
    return ret;
}

static bool
triangle_geometric_normal(const fast_output &output, size_t triangle,
	ON_3dVector &normal)
{
    ON_3dPoint points[3];
    for (int corner = 0; corner < 3; ++corner) {
	const int vertex = output.faces[triangle * 3 + corner];
	if (vertex < 0 || (size_t)vertex >= output.points.size() / 3)
	    return false;
	points[corner] = ON_3dPoint(&output.points[(size_t)vertex * 3]);
    }
    normal = ON_CrossProduct(points[1] - points[0], points[2] - points[0]);
    return normal.IsValid() && normal.Unitize();
}

static bool
normals_match_winding(const fast_output &output)
{
    if (output.faces.empty() || output.faces.size() % 3 ||
	    output.normals.size() != output.faces.size() * 3)
	return false;
    for (size_t triangle = 0; triangle < output.faces.size() / 3; ++triangle) {
	ON_3dVector geometric;
	if (!triangle_geometric_normal(output, triangle, geometric))
	    return false;
	for (int corner = 0; corner < 3; ++corner) {
	    ON_3dVector normal(&output.normals[triangle * 9 + corner * 3]);
	    if (!normal.IsValid() || !normal.Unitize() ||
		    !(geometric * normal > 0.0))
		return false;
	}
    }
    return true;
}

static bool
planar_winding_matches_source(const ON_Brep &brep)
{
    /* Check source orientation independently of the emitted normals, so
     * the display fallback cannot hide a reversed triangulator winding. */
    for (int face_index = 0; face_index < brep.m_F.Count(); ++face_index) {
	const ON_BrepFace &face = brep.m_F[face_index];
	const ON_Surface *surface = face.SurfaceOf();
	if (!surface || !surface->IsPlanar())
	    continue;
	ON_3dPoint point;
	ON_3dVector expected;
	if (!surface->EvNormal(surface->Domain(0).Mid(),
		surface->Domain(1).Mid(), point, expected))
	    return false;
	if (face.m_bRev)
	    expected = -expected;
	fast_output output;
	if (run_fast(&output, &brep, 1, 0, 0, 0, true, face_index) !=
		BREP_CDT_FAST_OK || output.faces.empty() ||
		!output.report.peak_working_bytes)
	    return false;
	for (size_t triangle = 0; triangle < output.faces.size() / 3; ++triangle) {
	    ON_3dVector geometric;
	    if (!triangle_geometric_normal(output, triangle, geometric) ||
		    !(geometric * expected > 0.0))
		return false;
	}
    }
    return true;
}

int
main(int argc, const char **argv)
{
    bu_setprogname(argv[0]);
    if (argc != 3)
	return 2;
    if (cdt_test_fast_display_normals())
	return 1;

    ON_Brep empty_brep;
    fast_output empty;
    fast_output nonexistent_face;
    if (run_fast(&empty, &empty_brep, 1, 0) != BREP_CDT_FAST_OK ||
	    !empty.faces.empty() || !empty.points.empty() || !empty.normals.empty() ||
	    empty.report.requested_faces || empty.report.failed_faces ||
	    empty.report.completed_faces ||
	    run_fast(&nonexistent_face, &empty_brep, 1, 0, 0, 0, true, 0) !=
	    BREP_CDT_FAST_ERROR)
	return 1;

    struct db_i *dbip = db_open(argv[1], DB_OPEN_READONLY);
    if (dbip == DBI_NULL || db_dirbuild(dbip) < 0)
	return 2;
    struct directory *dp = db_lookup(dbip, argv[2], LOOKUP_QUIET);
    if (dp == RT_DIR_NULL) {
	db_close(dbip);
	return 2;
    }

    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, dbip, NULL) < 0 ||
	    intern.idb_minor_type != DB5_MINORTYPE_BRLCAD_BREP) {
	db_close(dbip);
	return 2;
    }
    struct rt_brep_internal *bi = (struct rt_brep_internal *)intern.idb_ptr;
    RT_BREP_CK_MAGIC(bi);

    std::vector<ON_Interval> domains;
    for (int i = 0; i < bi->brep->m_F.Count(); i++) {
	const ON_Surface *surface = bi->brep->m_F[i].SurfaceOf();
	domains.push_back(surface->Domain(0));
	domains.push_back(surface->Domain(1));
    }
    std::vector<void *> trim_users;
    for (int i = 0; i < bi->brep->m_T.Count(); i++)
	trim_users.push_back(bi->brep->m_T[i].m_trim_user.p);

    fast_output serial;
    fast_output parallel;
    const size_t working_budget = (size_t)64 * 1024 * 1024;
    if (run_fast(&serial, bi->brep, 1, 16 * 1024 * 1024) !=
	    BREP_CDT_FAST_OK ||
	    run_fast(&parallel, bi->brep, 4, 16 * 1024 * 1024,
		working_budget) !=
	    BREP_CDT_FAST_OK) {
	rt_db_free_internal(&intern);
	db_close(dbip);
	return 1;
    }

    bool same = serial.faces == parallel.faces &&
	serial.normals == parallel.normals && serial.points == parallel.points;
    ON_Brep reversed_brep(*bi->brep);
    reversed_brep.Flip();
    fast_output reversed;
    const bool oriented_normals = normals_match_winding(serial) &&
	run_fast(&reversed, &reversed_brep, 1, 16 * 1024 * 1024) ==
	    BREP_CDT_FAST_OK && normals_match_winding(reversed) &&
	planar_winding_matches_source(*bi->brep) &&
	planar_winding_matches_source(reversed_brep);
    bool complete = serial.report.failed_faces == 0 &&
	parallel.report.failed_faces == 0 &&
	serial.report.completed_faces == bi->brep->m_F.Count() &&
	parallel.report.completed_faces == bi->brep->m_F.Count();
    bool working_bounded = parallel.report.peak_working_bytes > 0 &&
	parallel.report.peak_working_bytes <= working_budget;
    bool adaptive_reported = serial.report.triangle_budget > 0 &&
	serial.report.triangle_budget <= (size_t)256 * 1024 &&
	serial.report.refinement_passes == parallel.report.refinement_passes &&
	serial.report.approximated_faces ==
	parallel.report.approximated_faces &&
	serial.report.boundary_envelope_incomplete_faces == 0 &&
	parallel.report.boundary_envelope_incomplete_faces == 0;

    fast_output limited;
    int limit_ret = run_fast(&limited, bi->brep, 4, 1);
    bool limited_cleanly = limit_ret == BREP_CDT_FAST_LIMIT &&
	limited.report.hit_point_limit && limited.faces.empty() &&
	limited.normals.empty() && limited.points.empty();

    fast_output triangle_targeted;
    int triangle_target_ret = run_fast(&triangle_targeted, bi->brep, 4,
	16 * 1024 * 1024, 0, 1);
    bool authoritative_boundaries_retained =
	triangle_target_ret == BREP_CDT_FAST_OK &&
	triangle_targeted.report.triangle_budget == 1 &&
	triangle_targeted.report.triangle_budget_limited_faces > 0 &&
	triangle_targeted.report.boundary_envelope_incomplete_faces == 0 &&
	!triangle_targeted.faces.empty();

    bool unchanged = true;
    size_t domain_index = 0;
    for (int i = 0; i < bi->brep->m_F.Count(); i++) {
	const ON_Surface *surface = bi->brep->m_F[i].SurfaceOf();
	unchanged = unchanged &&
	    surface->Domain(0) == domains[domain_index++] &&
	    surface->Domain(1) == domains[domain_index++];
    }
    for (int i = 0; i < bi->brep->m_T.Count(); i++)
	unchanged = unchanged &&
	    bi->brep->m_T[i].m_trim_user.p == trim_users[(size_t)i];

    rt_db_free_internal(&intern);
    db_close(dbip);
    return (same && oriented_normals && complete && working_bounded &&
	adaptive_reported &&
	unchanged && limited_cleanly && authoritative_boundaries_retained) ? 0 : 1;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
