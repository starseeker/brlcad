/*               E V A L U A T E D _ P O I N T S . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/evaluated_points.h
 *
 * libbrlobol evaluated-points provider.
 *
 * This API resolves a database full path and returns a source-local sampled
 * indexed-face representation suitable for Obol database-source realization.
 * Callers own successful output arrays and must release them with
 * brlobol_evaluated_points_face_set_free.
 */

#ifndef BRLOBOL_EVALUATED_POINTS_H
#define BRLOBOL_EVALUATED_POINTS_H

#include "common.h"
#include "brlobol/defines.h"

struct db_i;
struct rt_primitive_indexed_face_set;

__BEGIN_DECLS

BRLOBOL_EXPORT extern int
brlobol_evaluated_points_evaluate_path_face_set(
	struct db_i *dbip,
	const char *path,
	struct rt_primitive_indexed_face_set *face_set);

BRLOBOL_EXPORT extern void
brlobol_evaluated_points_face_set_free(
	struct rt_primitive_indexed_face_set *face_set);

__END_DECLS

#endif /* BRLOBOL_EVALUATED_POINTS_H */
