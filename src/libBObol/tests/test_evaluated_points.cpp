/*              T E S T _ E V A L U A T E D _ P O I N T S . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */

#include "common.h"

#include "BObol/BEvaluatedPoints.h"

#include "bu/app.h"
#include "bu/file.h"
#include "raytrace.h"
#include "rt/wdb.h"
#include "wdb.h"

#include <cstdio>
#include <cstring>

static int
make_path_transform_db(const char *dbpath)
{
    rt_wdb *wdbp = wdb_fopen_v(dbpath, 5);
    if (!wdbp)
	return 1;

    point_t bmin = {-1.0, -1.0, -1.0};
    point_t bmax = { 1.0,  1.0,  1.0};
    if (mk_rpp(wdbp, "leaf.s", bmin, bmax) != 0) {
	wdb_close(wdbp);
	return 1;
    }

    struct wmember child_members;
    BU_LIST_INIT(&child_members.l);
    if (!mk_addmember("leaf.s", &child_members.l, NULL, WMOP_UNION) ||
	mk_comb(wdbp, "child.c", &child_members.l, 0, NULL, NULL,
		NULL, 0, 0, 0, 0, 0, 0, 0) != 0) {
	wdb_close(wdbp);
	return 1;
    }

    mat_t xlate;
    MAT_IDN(xlate);
    MAT_DELTAS(xlate, 50.0, 0.0, 0.0);
    struct wmember parent_members;
    BU_LIST_INIT(&parent_members.l);
    if (!mk_addmember("child.c", &parent_members.l, xlate, WMOP_UNION) ||
	mk_comb(wdbp, "parent.c", &parent_members.l, 0, NULL, NULL,
		NULL, 0, 0, 0, 0, 0, 0, 0) != 0) {
	wdb_close(wdbp);
	return 1;
    }

    wdb_close(wdbp);
    return 0;
}

static int
test_parent_path_matrix(void)
{
    char dbpath[MAXPATHLEN] = {0};
    db_i *dbip = NULL;
    struct rt_primitive_indexed_face_set face_set;
    const point_t expected_min = {49.0, -1.0, -1.0};
    const point_t expected_max = {51.0, 1.0, 1.0};
    memset(&face_set, 0, sizeof(face_set));

    FILE *fp = bu_temp_file(dbpath, MAXPATHLEN);
    if (!fp) {
	std::printf("FAIL: evaluated-points temp file\n");
	return 1;
    }
    std::fclose(fp);

    int ret = 0;
    if (make_path_transform_db(dbpath)) {
	std::printf("FAIL: evaluated-points fixture database\n");
	ret = 1;
	goto cleanup;
    }

    dbip = db_open(dbpath, DB_OPEN_READONLY);
    if (!dbip || db_dirbuild(dbip) < 0) {
	std::printf("FAIL: evaluated-points db open\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_evaluated_points_evaluate_path_face_set(dbip,
	    "parent.c/child.c", &face_set) != BRLCAD_OK ||
	!face_set.points || !face_set.point_count) {
	std::printf("FAIL: evaluated-points provider output\n");
	ret = 1;
	goto cleanup_db;
    }

    for (size_t i = 0; i < face_set.point_count; i++) {
	if (face_set.points[i][X] < 45.0) {
	    std::printf("FAIL: evaluated-points full-path matrix not applied\n");
	    ret = 1;
	    break;
	}
    }

    if (!ret && (!face_set.source_bounds_valid ||
	!VNEAR_EQUAL(face_set.source_bounds_min, expected_min, 0.000001) ||
	!VNEAR_EQUAL(face_set.source_bounds_max, expected_max, 0.000001))) {
	std::printf("FAIL: evaluated-points authoritative bounds not transformed\n");
	ret = 1;
    }

cleanup_db:
    bobol_evaluated_points_face_set_free(&face_set);
    if (dbip)
	db_close(dbip);
cleanup:
    if (dbpath[0])
	bu_file_delete(dbpath);
    return ret;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    (void)argc;
    return test_parent_path_matrix();
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
