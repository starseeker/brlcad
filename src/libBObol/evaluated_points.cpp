/*             E V A L U A T E D _ P O I N T S . C P P
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
/** @addtogroup libBObol */
/** @{ */
/** @file libBObol/evaluated_points.cpp
 *
 * libBObol evaluated-points provider.
 *
 */
/** @} */

#include "common.h"

#include "BObol/BEvaluatedPoints.h"

#include "bn/mat.h"
#include "bu/malloc.h"
#include "raytrace.h"
#include "rt/db_fullpath.h"
#include "rt/func.h"

#include <string.h>


void
bobol_evaluated_points_face_set_free(
    struct rt_primitive_indexed_face_set *face_set)
{
    if (!face_set)
	return;

    if (face_set->points)
	bu_free(face_set->points, "evaluated-points face-set points");
    if (face_set->normals)
	bu_free(face_set->normals, "evaluated-points face-set normals");
    if (face_set->indices)
	bu_free(face_set->indices, "evaluated-points face-set indices");
    memset(face_set, 0, sizeof(*face_set));
}


static int
bobol_evaluated_points_transform_face_set(
    struct rt_primitive_indexed_face_set *face_set,
    const mat_t matrix)
{
    if (!face_set || !matrix)
	return BRLCAD_ERROR;
    if (bn_mat_is_identity(matrix))
	return BRLCAD_OK;

    for (size_t i = 0; i < face_set->point_count; i++) {
	point_t transformed;
	MAT4X3PNT(transformed, matrix, face_set->points[i]);
	VMOVE(face_set->points[i], transformed);
    }

    if (!face_set->normals || !face_set->normal_count)
	return BRLCAD_OK;

    mat_t inverse;
    mat_t normal_matrix;
    if (!bn_mat_inverse(inverse, matrix))
	return BRLCAD_ERROR;
    bn_mat_trn(normal_matrix, inverse);

    for (size_t i = 0; i < face_set->normal_count; i++) {
	vect_t transformed;
	MAT4X3VEC(transformed, normal_matrix, face_set->normals[i]);
	if (MAGNITUDE(transformed) > SMALL_FASTF)
	    VUNITIZE(transformed);
	VMOVE(face_set->normals[i], transformed);
    }

    return BRLCAD_OK;
}


int
bobol_evaluated_points_evaluate_path_face_set(
    struct db_i *dbip,
    const char *path,
    struct rt_primitive_indexed_face_set *face_set)
{
    struct db_full_path full_path;
    mat_t path_matrix;
    struct rt_db_internal intern;
    int ret = BRLCAD_ERROR;

    if (face_set)
	memset(face_set, 0, sizeof(*face_set));
    if (!dbip || !path || !path[0] || !face_set)
	return BRLCAD_ERROR;

    db_full_path_init(&full_path);
    if (db_string_to_path(&full_path, dbip, path) != 0)
	return BRLCAD_ERROR;

    MAT_IDN(path_matrix);
    if (!db_path_to_mat(dbip, &full_path, path_matrix,
	    (int)full_path.fp_len - 1)) {
	db_free_full_path(&full_path);
	return BRLCAD_ERROR;
    }

    struct directory *dp = DB_FULL_PATH_CUR_DIR(&full_path);
    if (!dp) {
	db_free_full_path(&full_path);
	return BRLCAD_ERROR;
    }

    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, dbip, NULL) < 0) {
	db_free_full_path(&full_path);
	return BRLCAD_ERROR;
    }

    ret = rt_obj_sampled_face_set(face_set, &intern);
    rt_db_free_internal(&intern);
    db_free_full_path(&full_path);

    if (ret == BRLCAD_OK) {
	ret = bobol_evaluated_points_transform_face_set(face_set,
		path_matrix);
	if (ret != BRLCAD_OK)
	    bobol_evaluated_points_face_set_free(face_set);
    }

    return ret;
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
