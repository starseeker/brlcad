/*                   A R B N _ E D G E _ C A S E S . C
 * BRL-CAD
 *
 * Copyright (c) 2025 United States Government as represented by
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

#include "common.h"

#include <stdio.h>

#include "bu/app.h"
#include "bu/log.h"
#include "vmath.h"
#include "raytrace.h"


/**
 * Test that rt_arbn_bbox handles edge cases with invalid plane counts
 * without segfaulting. This test specifically checks for the unsigned
 * underflow issue when neqn < 3.
 */
int
main(int ac, char *av[])
{
    struct rt_db_internal intern;
    struct rt_arbn_internal *arbn;
    point_t min, max;
    struct bn_tol tol;
    int ret;
    int test_passed = 1;

    bu_setprogname(av[0]);

    if (ac > 1) {
	bu_exit(1, "Usage: %s\n", av[0]);
    }

    BN_TOL_INIT(&tol);

    /* Test 1: ARBN with 0 planes (should not crash) */
    bu_log("Test 1: ARBN with 0 planes\n");
    RT_DB_INTERNAL_INIT(&intern);
    intern.idb_major_type = DB5_MAJORTYPE_BRLCAD;
    intern.idb_type = ID_ARBN;
    intern.idb_meth = &OBJ[ID_ARBN];
    
    BU_ALLOC(arbn, struct rt_arbn_internal);
    arbn->magic = RT_ARBN_INTERNAL_MAGIC;
    arbn->neqn = 0;
    arbn->eqn = NULL;
    intern.idb_ptr = (void *)arbn;

    ret = rt_arbn_bbox(&intern, &min, &max, &tol);
    if (ret < 0) {
	bu_log("  PASS: rt_arbn_bbox correctly rejected 0 planes\n");
    } else {
	bu_log("  FAIL: rt_arbn_bbox should have rejected 0 planes\n");
	test_passed = 0;
    }
    rt_db_free_internal(&intern);

    /* Test 2: ARBN with 1 plane (should not crash) */
    bu_log("Test 2: ARBN with 1 plane\n");
    RT_DB_INTERNAL_INIT(&intern);
    intern.idb_major_type = DB5_MAJORTYPE_BRLCAD;
    intern.idb_type = ID_ARBN;
    intern.idb_meth = &OBJ[ID_ARBN];
    
    BU_ALLOC(arbn, struct rt_arbn_internal);
    arbn->magic = RT_ARBN_INTERNAL_MAGIC;
    arbn->neqn = 1;
    arbn->eqn = (plane_t *)bu_calloc(1, sizeof(plane_t), "arbn eqn");
    HSET(arbn->eqn[0], 1, 0, 0, 100);
    intern.idb_ptr = (void *)arbn;

    ret = rt_arbn_bbox(&intern, &min, &max, &tol);
    if (ret < 0) {
	bu_log("  PASS: rt_arbn_bbox correctly rejected 1 plane\n");
    } else {
	bu_log("  FAIL: rt_arbn_bbox should have rejected 1 plane\n");
	test_passed = 0;
    }
    rt_db_free_internal(&intern);

    /* Test 3: ARBN with 2 planes (should not crash) */
    bu_log("Test 3: ARBN with 2 planes\n");
    RT_DB_INTERNAL_INIT(&intern);
    intern.idb_major_type = DB5_MAJORTYPE_BRLCAD;
    intern.idb_type = ID_ARBN;
    intern.idb_meth = &OBJ[ID_ARBN];
    
    BU_ALLOC(arbn, struct rt_arbn_internal);
    arbn->magic = RT_ARBN_INTERNAL_MAGIC;
    arbn->neqn = 2;
    arbn->eqn = (plane_t *)bu_calloc(2, sizeof(plane_t), "arbn eqn");
    HSET(arbn->eqn[0], 1, 0, 0, 100);
    HSET(arbn->eqn[1], -1, 0, 0, 100);
    intern.idb_ptr = (void *)arbn;

    ret = rt_arbn_bbox(&intern, &min, &max, &tol);
    if (ret < 0) {
	bu_log("  PASS: rt_arbn_bbox correctly rejected 2 planes\n");
    } else {
	bu_log("  FAIL: rt_arbn_bbox should have rejected 2 planes\n");
	test_passed = 0;
    }
    rt_db_free_internal(&intern);

    /* Test 4: Valid ARBN with 8 planes forming a cube (should work) */
    bu_log("Test 4: Valid ARBN with 8 planes (cube)\n");
    RT_DB_INTERNAL_INIT(&intern);
    intern.idb_major_type = DB5_MAJORTYPE_BRLCAD;
    intern.idb_type = ID_ARBN;
    intern.idb_meth = &OBJ[ID_ARBN];
    
    BU_ALLOC(arbn, struct rt_arbn_internal);
    arbn->magic = RT_ARBN_INTERNAL_MAGIC;
    arbn->neqn = 8;
    arbn->eqn = (plane_t *)bu_calloc(8, sizeof(plane_t), "arbn eqn");
    /* Define a cube centered at origin with faces at ±100 */
    HSET(arbn->eqn[0],  1,  0,  0, 100);  /* x >= -100 */
    HSET(arbn->eqn[1], -1,  0,  0, 100);  /* x <= 100 */
    HSET(arbn->eqn[2],  0,  1,  0, 100);  /* y >= -100 */
    HSET(arbn->eqn[3],  0, -1,  0, 100);  /* y <= 100 */
    HSET(arbn->eqn[4],  0,  0,  1, 100);  /* z >= -100 */
    HSET(arbn->eqn[5],  0,  0, -1, 100);  /* z <= 100 */
    /* Add two diagonal planes to make it properly bounded */
    HSET(arbn->eqn[6],  0.57735,  0.57735,  0.57735, 200);
    HSET(arbn->eqn[7], -0.57735, -0.57735, -0.57735, 200);
    intern.idb_ptr = (void *)arbn;

    ret = rt_arbn_bbox(&intern, &min, &max, &tol);
    if (ret == 0) {
	bu_log("  PASS: rt_arbn_bbox accepted valid 8-plane ARBN\n");
	bu_log("  Bounding box: min=(%.2f, %.2f, %.2f) max=(%.2f, %.2f, %.2f)\n",
	       V3ARGS(min), V3ARGS(max));
    } else {
	bu_log("  FAIL: rt_arbn_bbox should have accepted valid 8-plane ARBN\n");
	test_passed = 0;
    }
    rt_db_free_internal(&intern);

    if (test_passed) {
	bu_log("\nAll tests PASSED\n");
	return 0;
    } else {
	bu_log("\nSome tests FAILED\n");
	return 1;
    }
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
