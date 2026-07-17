/*              B E V A L U A T E D W I R E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BEvaluatedWire.h
 *
 * libBObol compatibility adapter for librt evaluated-wire output.
 *
 * New evaluated-wire correctness work belongs in rt_eval_wireframe().  This
 * API keeps existing Obol/GED publication code on source-local neutral BG
 * line-command arrays.  Callers own successful output arrays and must release
 * them with bobol_evaluated_wire_line_set_free.
 */

#ifndef BOBOL_BEVALUATEDWIRE_H
#define BOBOL_BEVALUATEDWIRE_H

#include "common.h"
#include "BObol/BDefines.h"
#include "vmath.h"

#include <stddef.h>

struct bg_tess_tol;
struct bn_tol;
struct db_i;

__BEGIN_DECLS

BOBOL_EXPORT extern int
bobol_evaluated_wire_evaluate_path_line_set(
	struct db_i *dbip,
	const char *path,
	const struct bn_tol *tol,
	const struct bg_tess_tol *ttol,
	point_t **points_out,
	int **commands_out,
	size_t *count_out);

BOBOL_EXPORT extern void
bobol_evaluated_wire_line_set_free(
	point_t *points,
	int *commands);

__END_DECLS

#endif /* BOBOL_BEVALUATEDWIRE_H */
