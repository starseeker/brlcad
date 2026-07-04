/*                 E V A L U A T E D _ W I R E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/evaluated_wire.h
 *
 * libbrlobol evaluated-wire provider boundary.
 *
 * This API returns source-local neutral BG line-command arrays.  Callers own
 * successful output arrays and must release them with
 * brlobol_evaluated_wire_line_set_free.
 */

#ifndef BRLOBOL_EVALUATED_WIRE_H
#define BRLOBOL_EVALUATED_WIRE_H

#include "common.h"
#include "brlobol/defines.h"
#include "vmath.h"

#include <stddef.h>

struct bg_tess_tol;
struct bn_tol;
struct db_i;

__BEGIN_DECLS

BRLOBOL_EXPORT extern int
brlobol_evaluated_wire_evaluate_path_line_set(
	struct db_i *dbip,
	const char *path,
	const struct bn_tol *tol,
	const struct bg_tess_tol *ttol,
	point_t **points_out,
	int **commands_out,
	size_t *count_out);

BRLOBOL_EXPORT extern void
brlobol_evaluated_wire_line_set_free(
	point_t *points,
	int *commands);

__END_DECLS

#endif /* BRLOBOL_EVALUATED_WIRE_H */
