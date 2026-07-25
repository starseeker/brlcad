/*                       B R E P . C P P
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
/** @file brep.cpp
 *
 * Unit tests for BREP primitive editing via edbrep.cpp.
 *
 * Reference BREP: a 7x7 clamped planar NURBS surface with its standard
 * outer B-rep loop.  The larger control net has CVs whose basis support
 * does not intersect the four boundary trims, permitting topology-safe
 * interactive edit tests.
 *
 * Tests verify:
 *   - ECMD_BREP_SRF_SELECT stores face/i/j correctly
 *   - ECMD_BREP_SRF_CV_MOVE translates the CV by the given delta
 *   - ECMD_BREP_SRF_CV_SET places the CV at the given absolute position
 *   - ECMD_BREP_SRF_CV_WEIGHT changes an interior CV weight
 *   - CV edits whose support reaches an outer trim are rejected
 *   - rt_edit_brep_get_params returns sensible values
 *   - Invalid inputs are rejected gracefully (wrong e_inpara, bad face index)
 */

#include "common.h"

#include <cmath>
#include <cstring>

#include "vmath.h"
#include "bu/app.h"
#include "bu/log.h"
#include "bu/malloc.h"
#include "bu/vls.h"
#include "raytrace.h"
#include "edit_test_view.h"
#include "rt/geom.h"
#include "rt/rt_ecmds.h"
#include "brep/util.h"

/* Mirror the private struct from edbrep.cpp so the test can read state. */
struct rt_brep_edit_local {
    int face_index;
    int srf_cv_i;
    int srf_cv_j;
};


/* ------------------------------------------------------------------ *
 * Fixture helpers
 * ------------------------------------------------------------------ */

static struct directory *
make_brep_surface(struct rt_wdb *wdbp)
{
    const char *objname = "brep_surface";

    struct rt_brep_internal *bi;
    BU_ALLOC(bi, struct rt_brep_internal);
    bi->magic = RT_BREP_INTERNAL_MAGIC;
    bi->brep = new ON_Brep;

    ON_NurbsSurface surface(3, false, 4, 4, 7, 7);
    for (int i = 0; i < 7; ++i) {
	for (int j = 0; j < 7; ++j)
	    surface.SetCV(i, j, ON_3dPoint((double)i, (double)j, 0.0));
    }
    surface.MakeClampedUniformKnotVector(0, 1.0);
    surface.MakeClampedUniformKnotVector(1, 1.0);
    if (!bi->brep->NewFace(surface))
	bu_exit(1, "ERROR: unable to create trimmed NURBS face\n");

    wdb_export(wdbp, objname, (void *)bi, ID_BREP, 1.0);

    struct directory *dp = db_lookup(wdbp->dbip, objname, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL)
	bu_exit(1, "ERROR: unable to create brep sphere object\n");

    return dp;
}

static struct rt_edit *
open_edit(struct directory *dp, struct db_i *dbip)
{
    struct bn_tol tol = BN_TOL_INIT_TOL;
    struct db_full_path fp;
    db_full_path_init(&fp);
    db_add_node_to_full_path(&fp, dp);

    struct rt_edit_view v;
    rt_edit_test_view_init_identity_size(&v, 1000.0);

    struct rt_edit *s = rt_edit_create(&fp, dbip, &tol, &v);
    s->mv_context = 1;
    return s;
}


/* ------------------------------------------------------------------ *
 * Helpers to read CV position from the edit state
 * ------------------------------------------------------------------ */

static void
get_cv_pos(struct rt_edit *s, int face_index, int cv_i, int cv_j,
	   double *x, double *y, double *z)
{
    struct rt_brep_internal *bip =
	(struct rt_brep_internal *)s->es_int.idb_ptr;
    RT_BREP_CK_MAGIC(bip);
    ON_Brep *brep = bip->brep;
    ON_BrepFace *face = brep->Face(face_index);
    const ON_Surface *surf = face->SurfaceOf();
    const ON_NurbsSurface *ns = dynamic_cast<const ON_NurbsSurface *>(surf);
    if (!ns)
	bu_exit(1, "ERROR: face %d has no NURBS surface\n", face_index);
    double *cv = ns->CV(cv_i, cv_j);
    *x = cv[0];
    *y = cv[1];
    *z = cv[2];
}

static double
get_cv_weight(struct rt_edit *s, int face_index, int cv_i, int cv_j)
{
    struct rt_brep_internal *bip =
	(struct rt_brep_internal *)s->es_int.idb_ptr;
    RT_BREP_CK_MAGIC(bip);
    const ON_NurbsSurface *ns =
	dynamic_cast<const ON_NurbsSurface *>(
	    bip->brep->m_F[face_index].SurfaceOf());
    if (!ns)
	bu_exit(1, "ERROR: face %d has no NURBS surface\n", face_index);
    ON_4dPoint cv;
    if (!ns->GetCV(cv_i, cv_j, ON::euclidean_rational, &cv.x))
	bu_exit(1, "ERROR: unable to read face %d CV (%d,%d)\n",
		face_index, cv_i, cv_j);
    return cv.w;
}


/* ------------------------------------------------------------------ *
 * Tests
 * ------------------------------------------------------------------ */

/* 1. SELECT stores indices correctly */
static void
test_brep_select(struct rt_edit *s)
{
    struct rt_brep_edit_local *b = (struct rt_brep_edit_local *)s->ipe_ptr;

    /* Select face 0, CV (0, 0) */
    EDOBJ[ID_BREP].ft_set_edit_mode(s, ECMD_BREP_SRF_SELECT);
    s->e_inpara = 3;
    s->e_para[0] = 0.0; /* face_index */
    s->e_para[1] = 0.0; /* cv_i */
    s->e_para[2] = 0.0; /* cv_j */

    rt_edit_process(s);

    if (b->face_index != 0 || b->srf_cv_i != 0 || b->srf_cv_j != 0)
	bu_exit(1,
		"ERROR: ECMD_BREP_SRF_SELECT: expected (0,0,0) got (%d,%d,%d)\n",
		b->face_index, b->srf_cv_i, b->srf_cv_j);

    bu_log("ECMD_BREP_SRF_SELECT (0,0,0) PASS\n");
}

/* 2. SELECT rejects invalid face index */
static void
test_brep_select_invalid_face(struct rt_edit *s)
{
    struct rt_brep_edit_local *b = (struct rt_brep_edit_local *)s->ipe_ptr;
    /* Reset to known state */
    b->face_index = -1; b->srf_cv_i = -1; b->srf_cv_j = -1;

    EDOBJ[ID_BREP].ft_set_edit_mode(s, ECMD_BREP_SRF_SELECT);
    s->e_inpara = 3;
    s->e_para[0] = 9999.0; /* out-of-range face */
    s->e_para[1] = 0.0;
    s->e_para[2] = 0.0;

    bu_vls_trunc(s->log_str, 0);
    rt_edit_process(s);

    if (b->face_index != -1)
	bu_exit(1,
		"ERROR: ECMD_BREP_SRF_SELECT invalid face did not reject: "
		"face_index=%d\n", b->face_index);

    bu_log("ECMD_BREP_SRF_SELECT invalid face PASS (rejected as expected)\n");
}

/* 3. SELECT rejects when e_inpara != 3 */
static void
test_brep_select_wrong_inpara(struct rt_edit *s)
{
    struct rt_brep_edit_local *b = (struct rt_brep_edit_local *)s->ipe_ptr;
    b->face_index = -1; b->srf_cv_i = -1; b->srf_cv_j = -1;

    EDOBJ[ID_BREP].ft_set_edit_mode(s, ECMD_BREP_SRF_SELECT);
    s->e_inpara = 1; /* wrong */
    s->e_para[0] = 0.0;

    rt_edit_process(s);

    if (b->face_index != -1)
	bu_exit(1,
		"ERROR: ECMD_BREP_SRF_SELECT wrong e_inpara did not reject\n");

    bu_log("ECMD_BREP_SRF_SELECT wrong e_inpara PASS\n");
}

/* 4. MOVE translates the selected CV */
static void
test_brep_cv_move(struct rt_edit *s)
{
    struct rt_brep_edit_local *b = (struct rt_brep_edit_local *)s->ipe_ptr;

    /* First select face 0 CV (3, 3), whose support misses all trims. */
    EDOBJ[ID_BREP].ft_set_edit_mode(s, ECMD_BREP_SRF_SELECT);
    s->e_inpara = 3;
    s->e_para[0] = 0.0; /* face 0 */
    s->e_para[1] = 3.0; /* cv_i = 3 */
    s->e_para[2] = 3.0; /* cv_j = 3 */
    rt_edit_process(s);

    if (b->face_index != 0 || b->srf_cv_i != 3 || b->srf_cv_j != 3)
	bu_exit(1, "ERROR: ECMD_BREP_SRF_SELECT (0,3,3) failed\n");

    double x0, y0, z0;
    get_cv_pos(s, 0, 3, 3, &x0, &y0, &z0);

    /* Translate by (1, 2, 3) */
    EDOBJ[ID_BREP].ft_set_edit_mode(s, ECMD_BREP_SRF_CV_MOVE);
    s->e_inpara = 3;
    s->e_para[0] = 1.0;
    s->e_para[1] = 2.0;
    s->e_para[2] = 3.0;
    rt_edit_process(s);

    double x1, y1, z1;
    get_cv_pos(s, 0, 3, 3, &x1, &y1, &z1);

    /* local2base == 1.0 so delta in model units */
    double ex = x0 + 1.0;
    double ey = y0 + 2.0;
    double ez = z0 + 3.0;

    if (fabs(x1 - ex) > 1e-6 || fabs(y1 - ey) > 1e-6 || fabs(z1 - ez) > 1e-6)
	bu_exit(1,
		"ERROR: ECMD_BREP_SRF_CV_MOVE: expected (%.6f,%.6f,%.6f) "
		"got (%.6f,%.6f,%.6f)\n",
		ex, ey, ez, x1, y1, z1);

    bu_log("ECMD_BREP_SRF_CV_MOVE PASS: delta (1,2,3) applied correctly\n");
}

/* 5. SET places the CV at an absolute position */
static void
test_brep_cv_set(struct rt_edit *s)
{
    struct rt_brep_edit_local *b = (struct rt_brep_edit_local *)s->ipe_ptr;

    /* Select face 0 CV (3, 3), whose support misses all trims. */
    EDOBJ[ID_BREP].ft_set_edit_mode(s, ECMD_BREP_SRF_SELECT);
    s->e_inpara = 3;
    s->e_para[0] = 0.0;
    s->e_para[1] = 3.0;
    s->e_para[2] = 3.0;
    rt_edit_process(s);

    if (b->face_index != 0 || b->srf_cv_i != 3 || b->srf_cv_j != 3)
	bu_exit(1, "ERROR: ECMD_BREP_SRF_SELECT (0,3,3) failed\n");

    /* Place at (5, 5, 5) */
    EDOBJ[ID_BREP].ft_set_edit_mode(s, ECMD_BREP_SRF_CV_SET);
    s->e_inpara = 3;
    s->e_para[0] = 5.0;
    s->e_para[1] = 5.0;
    s->e_para[2] = 5.0;
    rt_edit_process(s);

    double xc, yc, zc;
    get_cv_pos(s, 0, 3, 3, &xc, &yc, &zc);

    if (fabs(xc - 5.0) > 1e-6 || fabs(yc - 5.0) > 1e-6 || fabs(zc - 5.0) > 1e-6)
	bu_exit(1,
		"ERROR: ECMD_BREP_SRF_CV_SET: expected (5,5,5) "
		"got (%.6f,%.6f,%.6f)\n", xc, yc, zc);

    bu_log("ECMD_BREP_SRF_CV_SET PASS: CV placed at (5,5,5)\n");
}

/* 6. MOVE rejected when no CV is selected */
static void
test_brep_cv_move_no_selection(struct rt_edit *s)
{
    struct rt_brep_edit_local *b = (struct rt_brep_edit_local *)s->ipe_ptr;
    /* Force clear selection */
    b->face_index = -1; b->srf_cv_i = -1; b->srf_cv_j = -1;

    EDOBJ[ID_BREP].ft_set_edit_mode(s, ECMD_BREP_SRF_CV_MOVE);
    s->e_inpara = 3;
    s->e_para[0] = s->e_para[1] = s->e_para[2] = 1.0;

    bu_vls_trunc(s->log_str, 0);
    rt_edit_process(s);

    /* Verify the operation was a no-op (no crash, CV state unchanged).
     * log_str is cleared by rt_edit_process so we cannot inspect it here. */
    if (b->face_index != -1 || b->srf_cv_i != -1 || b->srf_cv_j != -1)
	bu_exit(1,
		"ERROR: ECMD_BREP_SRF_CV_MOVE no-selection: state changed "
		"unexpectedly: face=%d i=%d j=%d\n",
		b->face_index, b->srf_cv_i, b->srf_cv_j);

    bu_log("ECMD_BREP_SRF_CV_MOVE no-selection PASS\n");
}

/* 7. get_params returns the current selection */
static void
test_brep_get_params_select(struct rt_edit *s)
{
    /* Select face 0 CV (3, 3) */
    EDOBJ[ID_BREP].ft_set_edit_mode(s, ECMD_BREP_SRF_SELECT);
    s->e_inpara = 3;
    s->e_para[0] = 0.0;
    s->e_para[1] = 3.0;
    s->e_para[2] = 3.0;
    rt_edit_process(s);

    fastf_t vals[3] = {0, 0, 0};
    int n = EDOBJ[ID_BREP].ft_edit_get_params(s, ECMD_BREP_SRF_SELECT, vals);

    if (n != 3 || (int)vals[0] != 0 || (int)vals[1] != 3 || (int)vals[2] != 3)
	bu_exit(1,
		"ERROR: get_params(SELECT) returned n=%d vals=(%.0f,%.0f,%.0f)\n",
		n, vals[0], vals[1], vals[2]);

    bu_log("get_params(ECMD_BREP_SRF_SELECT) PASS: (0,3,3)\n");
}

/* 8. A trim-influencing CV is selectable but cannot be moved. */
static void
test_brep_boundary_cv_locked(struct rt_edit *s)
{
    EDOBJ[ID_BREP].ft_set_edit_mode(s, ECMD_BREP_SRF_SELECT);
    s->e_inpara = 3;
    s->e_para[0] = 0.0;
    s->e_para[1] = 0.0;
    s->e_para[2] = 0.0;
    rt_edit_process(s);

    double x0, y0, z0;
    get_cv_pos(s, 0, 0, 0, &x0, &y0, &z0);
    EDOBJ[ID_BREP].ft_set_edit_mode(s, ECMD_BREP_SRF_CV_MOVE);
    s->e_inpara = 3;
    s->e_para[0] = 1.0;
    s->e_para[1] = 2.0;
    s->e_para[2] = 3.0;
    rt_edit_process(s);

    double x1, y1, z1;
    get_cv_pos(s, 0, 0, 0, &x1, &y1, &z1);
    if (fabs(x1 - x0) > 1e-12 || fabs(y1 - y0) > 1e-12
	    || fabs(z1 - z0) > 1e-12)
	bu_exit(1, "ERROR: trim-influencing boundary CV was moved\n");
    bu_log("ECMD_BREP_SRF_CV_MOVE boundary lock PASS\n");
}

/* 9. An interior CV weight can be changed without moving another CV. */
static void
test_brep_cv_weight(struct rt_edit *s)
{
    EDOBJ[ID_BREP].ft_set_edit_mode(s, ECMD_BREP_SRF_SELECT);
    s->e_inpara = 3;
    s->e_para[0] = 0.0;
    s->e_para[1] = 3.0;
    s->e_para[2] = 3.0;
    rt_edit_process(s);

    EDOBJ[ID_BREP].ft_set_edit_mode(s, ECMD_BREP_SRF_CV_WEIGHT);
    s->e_inpara = 1;
    s->e_para[0] = 2.5;
    rt_edit_process(s);
    if (fabs(get_cv_weight(s, 0, 3, 3) - 2.5) > 1e-12)
	bu_exit(1, "ERROR: ECMD_BREP_SRF_CV_WEIGHT did not set weight\n");

    fastf_t val = 0.0;
    int n = EDOBJ[ID_BREP].ft_edit_get_params(
	    s, ECMD_BREP_SRF_CV_WEIGHT, &val);
    if (n != 1 || fabs(val - 2.5) > 1e-12)
	bu_exit(1, "ERROR: get_params(WEIGHT) returned n=%d val=%.17g\n",
		n, val);
    bu_log("ECMD_BREP_SRF_CV_WEIGHT PASS\n");
}

/* 10. Descriptor is well-formed */
static void
test_brep_edit_desc(struct directory *dp)
{
    const struct rt_edit_prim_desc *desc =
	EDOBJ[dp->d_minor_type].ft_edit_desc();

    if (!desc)
	bu_exit(1, "ERROR: rt_edit_brep_edit_desc() returned NULL\n");

    if (!desc->prim_type || !BU_STR_EQUAL(desc->prim_type, "brep"))
	bu_exit(1, "ERROR: prim_type is '%s', expected 'brep'\n",
		desc->prim_type ? desc->prim_type : "(null)");

    if (desc->ncmd != 4)
	bu_exit(1, "ERROR: ncmd=%d, expected 4\n", desc->ncmd);

    bu_log("rt_edit_brep_edit_desc PASS: prim_type='%s' ncmd=%d\n",
	   desc->prim_type, desc->ncmd);
}


/* ------------------------------------------------------------------ *
 * main
 * ------------------------------------------------------------------ */

int
main(int argc, char *argv[])
{
    bu_setprogname(argv[0]);
    if (argc != 1)
	return BRLCAD_ERROR;

    struct db_i *dbip = db_open_inmem();
    if (dbip == DBI_NULL)
	bu_exit(1, "ERROR: db_open_inmem failed\n");

    struct rt_wdb *wdbp = wdb_dbopen(dbip, RT_WDB_TYPE_DB_INMEM);
    struct directory *dp = make_brep_surface(wdbp);

    struct rt_edit *s = open_edit(dp, dbip);

    test_brep_edit_desc(dp);
    test_brep_select(s);
    test_brep_select_invalid_face(s);
    test_brep_select_wrong_inpara(s);
    test_brep_cv_move(s);
    test_brep_cv_set(s);
    test_brep_cv_move_no_selection(s);
    test_brep_get_params_select(s);
    test_brep_boundary_cv_locked(s);
    test_brep_cv_weight(s);

    rt_edit_destroy(s);
    db_close(dbip);

    bu_log("ALL brep edit tests PASSED\n");
    return 0;
}


/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
