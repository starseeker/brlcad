/*              I N D E X E D _ F A C E _ S E T . C
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/* Validate direct primitive indexed-face providers. */

#include "common.h"

#include "bu/app.h"

#include <math.h>
#include <string.h>

#include "bu/malloc.h"
#include "raytrace.h"
#include "rt/geom.h"

static void
face_set_free(struct rt_primitive_indexed_face_set *set)
{
    rt_primitive_indexed_face_set_free(set);
}

static int
face_set_closed_volume(const struct rt_primitive_indexed_face_set *set,
	double *volume)
{
    size_t *edges;
    int face[1024];
    int face_count = 0;
    double signed_volume = 0.0;

    if (!set || !set->points || !set->indices || !set->point_count ||
	set->point_count > 1024 || !volume)
	return 0;
    if (set->point_count > SIZE_MAX / set->point_count)
	return 0;
    edges = (size_t *)bu_calloc(set->point_count * set->point_count,
	sizeof(size_t), "indexed-face test edges");

    for (size_t i = 0; i <= set->index_count; i++) {
	const int index = i < set->index_count ? set->indices[i] : -1;
	if (index >= 0) {
	    if ((size_t)index >= set->point_count || face_count >= 1024) {
		bu_free(edges, "indexed-face test edges");
		return 0;
	    }
	    face[face_count++] = index;
	    continue;
	}
	if (!face_count)
	    continue;
	if (face_count < 3) {
	    bu_free(edges, "indexed-face test edges");
	    return 0;
	}
	for (int j = 0; j < face_count; j++) {
	    size_t a = (size_t)face[j];
	    size_t b = (size_t)face[(j + 1) % face_count];
	    if (a > b) {
		const size_t tmp = a;
		a = b;
		b = tmp;
	    }
	    edges[a * set->point_count + b]++;
	}
	for (int j = 1; j + 1 < face_count; j++) {
	    vect_t cross;
	    VCROSS(cross, set->points[face[j]], set->points[face[j + 1]]);
	    signed_volume += VDOT(set->points[face[0]], cross) / 6.0;
	}
	face_count = 0;
    }

    for (size_t a = 0; a < set->point_count; a++) {
	for (size_t b = a + 1; b < set->point_count; b++) {
	    const size_t count = edges[a * set->point_count + b];
	    if (count && count != 2) {
		bu_free(edges, "indexed-face test edges");
		return 0;
	    }
	}
    }
    bu_free(edges, "indexed-face test edges");
    *volume = signed_volume;
    return signed_volume > 0.0;
}

static int
check_arb(void)
{
    struct rt_arb_internal arb;
    struct rt_db_internal intern;
    struct rt_primitive_indexed_face_set set = {0};
    struct bg_tess_tol ttol = BG_TESS_TOL_INIT_ZERO;
    struct bn_tol tol = BN_TOL_INIT_TOL;
    struct bv_view_info view = BV_VIEW_INFO_INIT;
    double volume = 0.0;

    memset(&arb, 0, sizeof(arb));
    arb.magic = RT_ARB_INTERNAL_MAGIC;
    VSET(arb.pt[0], 0.0, 0.0, 0.0);
    VSET(arb.pt[1], 1.0, 0.0, 0.0);
    VSET(arb.pt[2], 1.0, 1.0, 0.0);
    VSET(arb.pt[3], 0.0, 1.0, 0.0);
    VSET(arb.pt[4], 0.0, 0.0, 1.0);
    VSET(arb.pt[5], 1.0, 0.0, 1.0);
    VSET(arb.pt[6], 1.0, 1.0, 1.0);
    VSET(arb.pt[7], 0.0, 1.0, 1.0);
    RT_DB_INTERNAL_INIT(&intern);
    intern.idb_type = ID_ARB8;
    intern.idb_meth = &OBJ[ID_ARB8];
    intern.idb_ptr = &arb;
    if (!intern.idb_meth->ft_indexed_face_set ||
	intern.idb_meth->ft_indexed_face_set(&set, &intern, &ttol, &tol,
	    &view) != BRLCAD_OK ||
	!face_set_closed_volume(&set, &volume) || fabs(volume - 1.0) > 1.0e-9) {
	face_set_free(&set);
	return 0;
    }
    face_set_free(&set);

    for (int i = 4; i < 8; i++)
	VSET(arb.pt[i], 0.5, 0.5, 1.0);
    volume = 0.0;
    if (intern.idb_meth->ft_indexed_face_set(&set, &intern, &ttol, &tol,
	    &view) != BRLCAD_OK ||
	!face_set_closed_volume(&set, &volume) ||
	fabs(volume - 1.0 / 3.0) > 1.0e-9) {
	face_set_free(&set);
	return 0;
    }
    face_set_free(&set);
    return 1;
}

static int
check_tgc(int apex)
{
    struct rt_tgc_internal tgc;
    struct rt_db_internal intern;
    struct rt_primitive_indexed_face_set set = {0};
    struct bg_tess_tol ttol = BG_TESS_TOL_INIT_ZERO;
    struct bn_tol tol = BN_TOL_INIT_TOL;
    struct bv_view_info view = BV_VIEW_INFO_INIT;
    double volume = 0.0;
    const double expected = apex ? 8.0 * M_PI / 3.0 : 8.0 * M_PI;

    memset(&tgc, 0, sizeof(tgc));
    tgc.magic = RT_TGC_INTERNAL_MAGIC;
    VSET(tgc.v, 0.0, 0.0, 0.0);
    VSET(tgc.h, 0.0, 0.0, 4.0);
    VSET(tgc.a, 2.0, 0.0, 0.0);
    VSET(tgc.b, 0.0, 1.0, 0.0);
    if (!apex) {
	VMOVE(tgc.c, tgc.a);
	VMOVE(tgc.d, tgc.b);
    }
    ttol.rel = 0.01;
    RT_DB_INTERNAL_INIT(&intern);
    intern.idb_type = ID_TGC;
    intern.idb_meth = &OBJ[ID_TGC];
    intern.idb_ptr = &tgc;
    if (!intern.idb_meth->ft_indexed_face_set ||
	intern.idb_meth->ft_indexed_face_set(&set, &intern, &ttol, &tol,
	    &view) != BRLCAD_OK ||
	!face_set_closed_volume(&set, &volume) ||
	fabs(volume - expected) / expected > 0.02) {
	face_set_free(&set);
	return 0;
    }
    face_set_free(&set);
    return 1;
}

static int
check_general_tgc(void)
{
    struct rt_tgc_internal tgc;
    struct rt_db_internal intern;
    struct rt_primitive_indexed_face_set set = {0};
    struct bg_tess_tol ttol = BG_TESS_TOL_INIT_ZERO;
    struct bn_tol tol = BN_TOL_INIT_TOL;
    struct bv_view_info view = BV_VIEW_INFO_INIT;
    double volume = 0.0;
    const double expected = 16.0 * M_PI / 3.0;

    memset(&tgc, 0, sizeof(tgc));
    tgc.magic = RT_TGC_INTERNAL_MAGIC;
    VSET(tgc.h, 0.0, 0.0, 4.0);
    VSET(tgc.a, 2.0, 0.0, 0.0);
    VSET(tgc.b, 0.0, 1.0, 0.0);
    VSET(tgc.c, 1.0, 0.0, 0.0);
    VSET(tgc.d, 0.0, 0.75, 0.0);
    ttol.rel = 0.01;
    RT_DB_INTERNAL_INIT(&intern);
    intern.idb_type = ID_TGC;
    intern.idb_meth = &OBJ[ID_TGC];
    intern.idb_ptr = &tgc;
    const int provider_status = intern.idb_meth->ft_indexed_face_set ?
	intern.idb_meth->ft_indexed_face_set(&set, &intern, &ttol, &tol,
	    &view) : BRLCAD_ERROR;
    const int closed = provider_status == BRLCAD_OK ?
	face_set_closed_volume(&set, &volume) : 0;
    if (provider_status != BRLCAD_OK || !closed ||
	fabs(volume - expected) / expected > 0.02) {
	bu_log("general TGC provider status=%d closed=%d volume=%.17g expected=%.17g points=%zu indices=%zu\n",
	    provider_status, closed, volume, expected, set.point_count,
	    set.index_count);
	face_set_free(&set);
	return 0;
    }
    face_set_free(&set);

    VSET(tgc.c, 1.0, 0.25, 0.0);
    if (intern.idb_meth->ft_indexed_face_set(&set, &intern, &ttol, &tol,
	    &view) == BRLCAD_OK || set.points || set.indices) {
	face_set_free(&set);
	return 0;
    }
    VSET(tgc.c, 1.0, 0.0, 0.0);
    VSETALL(tgc.d, 0.0);
    if (intern.idb_meth->ft_indexed_face_set(&set, &intern, &ttol, &tol,
	    &view) == BRLCAD_OK || set.points || set.indices) {
	face_set_free(&set);
	return 0;
    }
    return 1;
}

static int
check_ell(void)
{
    struct rt_ell_internal ell;
    struct rt_db_internal intern;
    struct rt_primitive_indexed_face_set set = {0};
    struct bg_tess_tol ttol = BG_TESS_TOL_INIT_ZERO;
    struct bn_tol tol = BN_TOL_INIT_TOL;
    struct bv_view_info view = BV_VIEW_INFO_INIT;
    double volume = 0.0;
    const double expected = 8.0 * M_PI;

    memset(&ell, 0, sizeof(ell));
    ell.magic = RT_ELL_INTERNAL_MAGIC;
    VSET(ell.v, 1.0, -2.0, 3.0);
    VSET(ell.a, 2.0, 0.0, 0.0);
    VSET(ell.b, 0.0, 1.0, 0.0);
    VSET(ell.c, 0.0, 0.0, 3.0);
    ttol.rel = 0.01;
    RT_DB_INTERNAL_INIT(&intern);
    intern.idb_type = ID_ELL;
    intern.idb_meth = &OBJ[ID_ELL];
    intern.idb_ptr = &ell;
    const int provider_status = intern.idb_meth->ft_indexed_face_set ?
	intern.idb_meth->ft_indexed_face_set(&set, &intern, &ttol, &tol,
	    &view) : BRLCAD_ERROR;
    const int closed = provider_status == BRLCAD_OK ?
	face_set_closed_volume(&set, &volume) : 0;
    if (provider_status != BRLCAD_OK || !closed ||
	fabs(volume - expected) / expected > 0.03) {
	bu_log("ELL provider status=%d closed=%d volume=%.17g expected=%.17g points=%zu indices=%zu\n",
	    provider_status, closed, volume, expected, set.point_count,
	    set.index_count);
	face_set_free(&set);
	return 0;
    }
    face_set_free(&set);
    return 1;
}

static int
check_tor(void)
{
    struct rt_tor_internal tor;
    struct rt_db_internal intern;
    struct rt_primitive_indexed_face_set set = {0};
    struct bg_tess_tol ttol = BG_TESS_TOL_INIT_ZERO;
    struct bn_tol tol = BN_TOL_INIT_TOL;
    struct bv_view_info view = BV_VIEW_INFO_INIT;
    double volume = 0.0;
    double expected;

    memset(&tor, 0, sizeof(tor));
    tor.magic = RT_TOR_INTERNAL_MAGIC;
    VSET(tor.v, -1.0, 2.0, 3.0);
    VSET(tor.h, 0.0, 0.0, 1.0);
    VSET(tor.a, 3.0, 0.0, 0.0);
    VSET(tor.b, 0.0, 3.0, 0.0);
    tor.r_a = 3.0;
    tor.r_b = 3.0;
    tor.r_h = 0.75;
    expected = 2.0 * M_PI * M_PI * tor.r_a * tor.r_h * tor.r_h;
    ttol.rel = 0.01;
    RT_DB_INTERNAL_INIT(&intern);
    intern.idb_type = ID_TOR;
    intern.idb_meth = &OBJ[ID_TOR];
    intern.idb_ptr = &tor;
    const int provider_status = intern.idb_meth->ft_indexed_face_set ?
	intern.idb_meth->ft_indexed_face_set(&set, &intern, &ttol, &tol,
	    &view) : BRLCAD_ERROR;
    const int closed = provider_status == BRLCAD_OK ?
	face_set_closed_volume(&set, &volume) : 0;
    if (provider_status != BRLCAD_OK || !closed ||
	fabs(volume - expected) / expected > 0.03) {
	bu_log("TOR provider status=%d closed=%d volume=%.17g expected=%.17g points=%zu indices=%zu\n",
	    provider_status, closed, volume, expected, set.point_count,
	    set.index_count);
	face_set_free(&set);
	return 0;
    }
    face_set_free(&set);

    tor.r_h = tor.r_a * 1.25;
    if (intern.idb_meth->ft_indexed_face_set(&set, &intern, &ttol, &tol,
	    &view) == BRLCAD_OK || set.points || set.indices) {
	face_set_free(&set);
	return 0;
    }
    return 1;
}

static int
check_database_object(const char *database, const char *object)
{
    struct db_i *dbip = DBI_NULL;
    struct directory *dp = RT_DIR_NULL;
    struct rt_db_internal intern;
    struct rt_primitive_indexed_face_set set = {0};
    struct bg_tess_tol ttol = BG_TESS_TOL_INIT_ZERO;
    struct bn_tol tol = BN_TOL_INIT_TOL;
    struct bv_view_info view = BV_VIEW_INFO_INIT;
    point_t mesh_min = VINIT_ZERO;
    point_t mesh_max = VINIT_ZERO;
    point_t source_min = VINIT_ZERO;
    point_t source_max = VINIT_ZERO;
    unsigned char *referenced = NULL;
    size_t referenced_count = 0;
    double volume = 0.0;
    int closed_check = -1;
    int result = 1;

    if (!database || !object)
	return 0;

    dbip = db_open(database, DB_OPEN_READONLY);
    if (dbip == DBI_NULL || db_dirbuild(dbip) < 0) {
	bu_log("Unable to open database %s\n", database);
	goto cleanup;
    }
    dp = db_lookup(dbip, object, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL) {
	bu_log("Unable to find object %s\n", object);
	goto cleanup;
    }

    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, dbip, NULL) < 0) {
	bu_log("Unable to import object %s\n", object);
	goto cleanup;
    }
    if (!intern.idb_meth || !intern.idb_meth->ft_indexed_face_set) {
	bu_log("Object %s has no indexed-face-set provider\n", object);
	rt_db_free_internal(&intern);
	goto cleanup;
    }

    ttol.rel = 0.01;
    if (intern.idb_meth->ft_indexed_face_set(&set, &intern, &ttol, &tol,
	    &view) != BRLCAD_OK || !set.points || !set.indices ||
	    !set.point_count || !set.index_count) {
	bu_log("Object %s did not produce an indexed face set\n", object);
	rt_db_free_internal(&intern);
	goto cleanup;
    }

    VMOVE(mesh_min, set.points[0]);
    VMOVE(mesh_max, set.points[0]);
    for (size_t i = 0; i < set.point_count; i++) {
	if (!isfinite(set.points[i][X]) || !isfinite(set.points[i][Y]) ||
	    !isfinite(set.points[i][Z])) {
	    bu_log("Object %s produced a non-finite point\n", object);
	    rt_db_free_internal(&intern);
	    goto cleanup;
	}
	VMINMAX(mesh_min, mesh_max, set.points[i]);
    }
    for (size_t i = 0; i < set.index_count; i++) {
	if (set.indices[i] < -1 ||
	    (set.indices[i] >= 0 && (size_t)set.indices[i] >= set.point_count)) {
	    bu_log("Object %s produced an invalid index\n", object);
	    rt_db_free_internal(&intern);
	    goto cleanup;
	}
    }

    referenced = (unsigned char *)bu_calloc(set.point_count,
	sizeof(unsigned char), "indexed-face referenced points");
    for (size_t i = 0; i < set.index_count; i++) {
	if (set.indices[i] >= 0)
	    referenced[set.indices[i]] = 1;
    }
    int have_referenced_bounds = 0;
    for (size_t i = 0; i < set.point_count; i++) {
	if (!referenced[i])
	    continue;
	if (!have_referenced_bounds) {
	    VMOVE(mesh_min, set.points[i]);
	    VMOVE(mesh_max, set.points[i]);
	    have_referenced_bounds = 1;
	} else {
	    VMINMAX(mesh_min, mesh_max, set.points[i]);
	}
	referenced_count++;
    }

    /* The simple closure helper is intentionally quadratic and bounded to
     * small unit-test meshes.  Report an unevaluated check distinctly from a
     * mesh which was actually tested and found open. */
    if (set.point_count <= 1024)
	closed_check = face_set_closed_volume(&set, &volume);
    const int source_bounds = rt_bound_internal(dbip, dp, source_min,
	source_max) == 0;
    bu_log("%s: points=%zu referenced=%zu indices=%zu closed_check=%d "
	"volume=%.17g "
	"mesh_min=(%.17g %.17g %.17g) mesh_max=(%.17g %.17g %.17g)",
	object, set.point_count, referenced_count, set.index_count, closed_check,
	volume,
	V3ARGS(mesh_min), V3ARGS(mesh_max));
    if (source_bounds)
	bu_log(" source_min=(%.17g %.17g %.17g) "
	    "source_max=(%.17g %.17g %.17g)",
	    V3ARGS(source_min), V3ARGS(source_max));
    bu_log("\n");

    rt_db_free_internal(&intern);
    result = 0;

cleanup:
    if (referenced)
	bu_free(referenced, "indexed-face referenced points");
    face_set_free(&set);
    if (dbip != DBI_NULL)
	db_close(dbip);
    return result;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    if (argc == 3)
	return check_database_object(argv[1], argv[2]);
    if (argc != 1) {
	bu_log("Usage: %s [file.g object]\n", argv[0]);
	return 1;
    }
    if (!check_arb()) {
	bu_log("ARB indexed-face validation failed\n");
	return 1;
    }
    if (!check_tgc(0) || !check_tgc(1)) {
	bu_log("TGC indexed-face validation failed\n");
	return 1;
    }
    if (!check_general_tgc()) {
	bu_log("general TGC indexed-face validation failed\n");
	return 1;
    }
    if (!check_ell()) {
	bu_log("ELL indexed-face validation failed\n");
	return 1;
    }
    if (!check_tor()) {
	bu_log("TOR indexed-face validation failed\n");
	return 1;
    }
    return 0;
}
