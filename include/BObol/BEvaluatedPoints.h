/*            B E V A L U A T E D P O I N T S . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BEvaluatedPoints.h
 *
 * libBObol evaluated-points provider.
 *
 * This API resolves a database full path and returns a source-local sampled
 * indexed-face representation suitable for Obol database-source realization.
 * Callers own successful output arrays and must release them with
 * bobol_evaluated_points_face_set_free.
 */

#ifndef BOBOL_BEVALUATEDPOINTS_H
#define BOBOL_BEVALUATEDPOINTS_H

#include "common.h"
#include "BObol/BDefines.h"

struct db_i;
struct rt_primitive_indexed_face_set;

__BEGIN_DECLS

BOBOL_EXPORT extern int
bobol_evaluated_points_evaluate_path_face_set(
	struct db_i *dbip,
	const char *path,
	struct rt_primitive_indexed_face_set *face_set);

BOBOL_EXPORT extern void
bobol_evaluated_points_face_set_free(
	struct rt_primitive_indexed_face_set *face_set);

__END_DECLS

#endif /* BOBOL_BEVALUATEDPOINTS_H */
