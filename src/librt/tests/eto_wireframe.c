/*                   E T O _ W I R E F R A M E . C
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/* Exercise numerically delicate ETO renderer-facing line providers. */

#include "common.h"

#include "bu/app.h"
#include "raytrace.h"
#include "rt/geom.h"

#include <math.h>
#include <string.h>

static int
line_set_is_finite(const struct rt_primitive_lod_realization *realization)
{
    if (!realization || !realization->has_line_set ||
	!realization->line_points || !realization->line_count)
	return 0;
    for (size_t i = 0; i < realization->line_count; i++)
	if (!isfinite(realization->line_points[i][X]) ||
	    !isfinite(realization->line_points[i][Y]) ||
	    !isfinite(realization->line_points[i][Z]))
	    return 0;
    return 1;
}

static void
report_line_set(const char *label, int result,
	const struct rt_primitive_lod_realization *realization)
{
    bu_log("%s result=%d has=%d count=%zu points=%p\n", label, result,
	realization->has_line_set, realization->line_count,
	(void *)realization->line_points);
    for (size_t i = 0; i < realization->line_count; i++) {
	if (!isfinite(realization->line_points[i][X]) ||
	    !isfinite(realization->line_points[i][Y]) ||
	    !isfinite(realization->line_points[i][Z])) {
	    bu_log("%s first non-finite point %zu: %g %g %g\n", label, i,
		realization->line_points[i][X], realization->line_points[i][Y],
		realization->line_points[i][Z]);
	    break;
	}
    }
}

static int
check_parallel_axis_eto(void)
{
    struct rt_eto_internal eto;
    struct rt_db_internal intern;
    struct rt_primitive_lod_realization realization;
    struct bg_tess_tol ttol = BG_TESS_TOL_INIT_ZERO;
    struct bn_tol tol = BN_TOL_INIT_TOL;
    struct bv_view_info view = BV_VIEW_INFO_INIT;
    int failed = 0;

    memset(&eto, 0, sizeof(eto));
    eto.eto_magic = RT_ETO_INTERNAL_MAGIC;
    VSET(eto.eto_V, 0.0, 0.0, 25.0);
    VSET(eto.eto_N, 0.0, 0.0, 0.98633320096251198);
    VSET(eto.eto_C, 0.0, 0.0, 4.93166600481256);
    eto.eto_r = 156.0;
    eto.eto_rd = 2.958999602887538;

    RT_DB_INTERNAL_INIT(&intern);
    intern.idb_type = ID_ETO;
    intern.idb_meth = &OBJ[ID_ETO];
    intern.idb_ptr = &eto;
    ttol.rel = 0.01;

    memset(&realization, 0, sizeof(realization));
    int result = intern.idb_meth->ft_wireframe_line_set ?
	intern.idb_meth->ft_wireframe_line_set(&realization, &intern, &ttol,
	    &tol) : BRLCAD_ERROR;
    if (result != BRLCAD_OK || !line_set_is_finite(&realization)) {
	bu_log("ETO standard wireframe line set was invalid\n");
	report_line_set("standard", result, &realization);
	failed = 1;
    }
    rt_primitive_lod_realization_free(&realization);

    memset(&realization, 0, sizeof(realization));
    result = intern.idb_meth->ft_lod_realize ?
	intern.idb_meth->ft_lod_realize(&realization, &intern, &tol, &view,
	    2.0 * (eto.eto_r + MAGNITUDE(eto.eto_C))) : BRLCAD_ERROR;
    if (result != BRLCAD_OK || !line_set_is_finite(&realization)) {
	bu_log("ETO LoD wireframe line set was invalid\n");
	report_line_set("LoD", result, &realization);
	failed = 1;
    }
    rt_primitive_lod_realization_free(&realization);
    return failed;
}

int
main(int UNUSED(argc), const char **argv)
{
    bu_setprogname(argv[0]);
    if (check_parallel_axis_eto()) {
	bu_log("ETO parallel-axis wireframe realization failed\n");
	return 1;
    }
    return 0;
}
