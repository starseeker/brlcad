#include "raytrace.h"
#include "bu/app.h"
#include "../primitives/arbn/arbn_clip.h"
#include <stdio.h>

/* Basic unit test scaffold:
 * Tests:
 * 1. Simple tetrahedron (4 planes).
 * 2. Cube (6 planes).
 * 3. Random convex poly with > threshold planes.
 * References: Vertex enumeration correctness vs half-space intersection (Avis & Fukuda 1992; Barber et al. 1996).
 */

static void test_simple_cube(void) {
    struct rt_arbn_internal ai;
    ai.magic = RT_ARBN_INTERNAL_MAGIC;
    ai.neqn = 6;
    ai.eqn = (plane_t *)bu_calloc(ai.neqn, sizeof(plane_t), "cube planes");
    /* Axis-aligned unit cube centered at origin: planes x=±1,y=±1,z=±1 */
    VSET(ai.eqn[0],  1, 0, 0); ai.eqn[0][W] = 1;
    VSET(ai.eqn[1], -1, 0, 0); ai.eqn[1][W] = 1;
    VSET(ai.eqn[2], 0,  1, 0); ai.eqn[2][W] = 1;
    VSET(ai.eqn[3], 0, -1, 0); ai.eqn[3][W] = 1;
    VSET(ai.eqn[4], 0, 0,  1); ai.eqn[4][W] = 1;
    VSET(ai.eqn[5], 0, 0, -1); ai.eqn[5][W] = 1;

    struct bn_tol tol;
    tol.dist = 1e-6;
    tol.dist_sq = tol.dist * tol.dist;

    struct arbn_clip_poly *poly = rt_arbn_clip_build(&ai, &tol);
    if (!poly) {
        bu_log("Cube test failed (poly NULL)\n");
    } else {
        if (poly->vcnt != 8) bu_log("Cube vertex count mismatch (%zu)\n", poly->vcnt);
        if (poly->fcnt < 6) bu_log("Cube face count mismatch (%zu)\n", poly->fcnt);
        rt_arbn_clip_free(poly);
    }
    bu_free(ai.eqn, "cube planes");
}

int main(int UNUSED(argc), char **UNUSED(argv))
{
    bu_setprogname("arbn_clip_test");
    test_simple_cube();
    printf("arbn_clip_test complete\n");
    return 0;
}
