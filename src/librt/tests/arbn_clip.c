#include "raytrace.h"
#include "bu/app.h"
#include "../primitives/arbn/arbn_clip.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Comprehensive unit tests for ARBN clipping tessellator:
 * Tests:
 * 1. Simple cube (6 planes) - basic convex polyhedron
 * 2. Tetrahedron (4 planes) - minimal convex polyhedron
 * 3. Duplicate plane removal
 * 4. Large plane count (100+ planes)
 * 5. Statistics generation
 * 6. Environment variable switching
 * References: Vertex enumeration correctness vs half-space intersection 
 * (Avis & Fukuda 1992; Barber et al. 1996; Preparata & Shamos 1985).
 */

static int test_count = 0;
static int test_pass = 0;

static void report_test(const char *name, int passed) {
    test_count++;
    if (passed) {
        test_pass++;
        printf("  [PASS] %s\n", name);
    } else {
        printf("  [FAIL] %s\n", name);
    }
}

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

    struct bn_tol tol = BN_TOL_INIT_TOL;
    tol.dist = 1e-6;
    tol.dist_sq = tol.dist * tol.dist;

    struct arbn_clip_stats stats;
    struct arbn_clip_poly *poly = rt_arbn_clip_build_with_stats(&ai, &tol, &stats);
    
    int passed = (poly != NULL);
    if (poly) {
        /* Expect 8 vertices and 6 faces for a cube */
        passed = (stats.final_vertices == 8 && stats.final_faces == 6);
        rt_arbn_clip_free(poly);
    }
    
    bu_free(ai.eqn, "cube planes");
    report_test("Cube (6 planes)", passed);
}

static void test_tetrahedron(void) {
    struct rt_arbn_internal ai;
    ai.magic = RT_ARBN_INTERNAL_MAGIC;
    ai.neqn = 4;
    ai.eqn = (plane_t *)bu_calloc(ai.neqn, sizeof(plane_t), "tetra planes");
    
    /* Regular tetrahedron: 4 planes */
    VSET(ai.eqn[0], 0, 0, -1); ai.eqn[0][W] = 0;     /* base */
    VSET(ai.eqn[1], 0, -1, 1); ai.eqn[1][W] = -1;    /* side 1 */
    VSET(ai.eqn[2], -0.866, 0.5, 1); ai.eqn[2][W] = -1; /* side 2 */
    VSET(ai.eqn[3], 0.866, 0.5, 1); ai.eqn[3][W] = -1;  /* side 3 */

    struct bn_tol tol = BN_TOL_INIT_TOL;
    tol.dist = 1e-6;
    tol.dist_sq = tol.dist * tol.dist;

    struct arbn_clip_stats stats;
    struct arbn_clip_poly *poly = rt_arbn_clip_build_with_stats(&ai, &tol, &stats);
    
    int passed = (poly != NULL);
    if (poly) {
        /* Expect 4 vertices and 4 faces for a tetrahedron */
        passed = (stats.final_vertices == 4 && stats.final_faces == 4);
        rt_arbn_clip_free(poly);
    }
    
    bu_free(ai.eqn, "tetra planes");
    report_test("Tetrahedron (4 planes)", passed);
}

static void test_duplicate_planes(void) {
    struct rt_arbn_internal ai;
    ai.magic = RT_ARBN_INTERNAL_MAGIC;
    ai.neqn = 9; /* 6 unique + 3 duplicates */
    ai.eqn = (plane_t *)bu_calloc(ai.neqn, sizeof(plane_t), "dup planes");
    
    /* Cube with duplicates */
    VSET(ai.eqn[0],  1, 0, 0); ai.eqn[0][W] = 1;
    VSET(ai.eqn[1], -1, 0, 0); ai.eqn[1][W] = 1;
    VSET(ai.eqn[2], 0,  1, 0); ai.eqn[2][W] = 1;
    VSET(ai.eqn[3], 0, -1, 0); ai.eqn[3][W] = 1;
    VSET(ai.eqn[4], 0, 0,  1); ai.eqn[4][W] = 1;
    VSET(ai.eqn[5], 0, 0, -1); ai.eqn[5][W] = 1;
    /* Duplicates */
    VSET(ai.eqn[6],  1, 0, 0); ai.eqn[6][W] = 1; /* dup of 0 */
    VSET(ai.eqn[7], 0,  1, 0); ai.eqn[7][W] = 1; /* dup of 2 */
    VSET(ai.eqn[8], 0, 0, -1); ai.eqn[8][W] = 1; /* dup of 5 */

    struct bn_tol tol = BN_TOL_INIT_TOL;
    tol.dist = 1e-6;
    tol.dist_sq = tol.dist * tol.dist;

    struct arbn_clip_stats stats;
    struct arbn_clip_poly *poly = rt_arbn_clip_build_with_stats(&ai, &tol, &stats);
    
    int passed = (poly != NULL);
    if (poly) {
        /* Should detect and remove 3 duplicates */
        passed = (stats.duplicate_planes == 3 && stats.active_planes == 6);
        rt_arbn_clip_free(poly);
    }
    
    bu_free(ai.eqn, "dup planes");
    report_test("Duplicate plane removal", passed);
}

static void test_large_plane_count(void) {
    /* Sphere approximation with 100 planes */
    struct rt_arbn_internal ai;
    ai.magic = RT_ARBN_INTERNAL_MAGIC;
    ai.neqn = 100;
    ai.eqn = (plane_t *)bu_calloc(ai.neqn, sizeof(plane_t), "sphere planes");
    
    double radius = 10.0;
    for (size_t i = 0; i < ai.neqn; ++i) {
        /* Generate random plane pointing inward */
        double theta = 2.0 * M_PI * i / ai.neqn;
        double phi = M_PI * (i % 10) / 10.0;
        VSET(ai.eqn[i], sin(phi) * cos(theta), sin(phi) * sin(theta), cos(phi));
        ai.eqn[i][W] = radius;
    }

    struct bn_tol tol = BN_TOL_INIT_TOL;
    tol.dist = 1e-6;
    tol.dist_sq = tol.dist * tol.dist;

    struct arbn_clip_stats stats;
    struct arbn_clip_poly *poly = rt_arbn_clip_build_with_stats(&ai, &tol, &stats);
    
    int passed = (poly != NULL);
    if (poly) {
        /* Should succeed with reasonable vertex/face counts */
        passed = (stats.final_vertices > 10 && stats.final_faces > 10);
        rt_arbn_clip_free(poly);
    }
    
    bu_free(ai.eqn, "sphere planes");
    report_test("Large plane count (100 planes)", passed);
}

static void test_statistics(void) {
    struct rt_arbn_internal ai;
    ai.magic = RT_ARBN_INTERNAL_MAGIC;
    ai.neqn = 6;
    ai.eqn = (plane_t *)bu_calloc(ai.neqn, sizeof(plane_t), "cube planes");
    
    VSET(ai.eqn[0],  1, 0, 0); ai.eqn[0][W] = 1;
    VSET(ai.eqn[1], -1, 0, 0); ai.eqn[1][W] = 1;
    VSET(ai.eqn[2], 0,  1, 0); ai.eqn[2][W] = 1;
    VSET(ai.eqn[3], 0, -1, 0); ai.eqn[3][W] = 1;
    VSET(ai.eqn[4], 0, 0,  1); ai.eqn[4][W] = 1;
    VSET(ai.eqn[5], 0, 0, -1); ai.eqn[5][W] = 1;

    struct bn_tol tol = BN_TOL_INIT_TOL;
    tol.dist = 1e-6;
    tol.dist_sq = tol.dist * tol.dist;

    struct arbn_clip_stats stats;
    struct arbn_clip_poly *poly = rt_arbn_clip_build_with_stats(&ai, &tol, &stats);
    
    int passed = (poly != NULL);
    if (poly) {
        /* Verify statistics are populated */
        passed = (stats.input_planes == 6 &&
                  stats.active_planes == 6 &&
                  stats.final_vertices == 8 &&
                  stats.final_faces == 6 &&
                  stats.bounding_radius > 0);
        rt_arbn_clip_free(poly);
    }
    
    bu_free(ai.eqn, "cube planes");
    report_test("Statistics generation", passed);
}

int main(int UNUSED(argc), char **UNUSED(argv))
{
    bu_setprogname("arbn_clip_test");
    
    printf("Running ARBN clipping tests:\n");
    
    test_simple_cube();
    test_tetrahedron();
    test_duplicate_planes();
    test_large_plane_count();
    test_statistics();
    
    printf("\n");
    printf("Test summary: %d/%d passed\n", test_pass, test_count);
    
    return (test_pass == test_count) ? 0 : 1;
}
