/*                 E V A L _ W I R E F R A M E . H
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
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file rt/eval_wireframe.h
 *
 * Evaluated CSG wireframe generation.
 */

#ifndef RT_EVAL_WIREFRAME_H
#define RT_EVAL_WIREFRAME_H

#include "common.h"

#include <stdint.h>

#include "rt/defines.h"
#include "rt/vlist.h"

struct bg_tess_tol;
struct bn_tol;
struct db_i;

__BEGIN_DECLS

/**
 * Use the current default evaluated-wireframe method.  The current default is
 * the symbolic boolweave-backed evaluator.  Callers needing the historical
 * candidate evaluator for comparison or fallback should request
 * RT_EVAL_WIREFRAME_F_BIGE explicitly.
 */
#define RT_EVAL_WIREFRAME_F_DEFAULT 0x00000000u

/**
 * Force the historical Big-E candidate evaluation path.
 */
#define RT_EVAL_WIREFRAME_F_BIGE 0x00000001u

/**
 * Force the symbolic boolweave-backed candidate evaluation path.  This is
 * also the default when no method flag is supplied.
 */
#define RT_EVAL_WIREFRAME_F_BOOLWEAVE 0x00000002u

/**
 * Emit evaluated-wireframe profiling summaries to bu_log.  This is intended
 * for implementation tuning and should not be enabled by normal interactive
 * callers.
 */
#define RT_EVAL_WIREFRAME_F_PROFILE 0x00010000u

#define RT_EVAL_WIREFRAME_F_METHOD_MASK \
    (RT_EVAL_WIREFRAME_F_BIGE | RT_EVAL_WIREFRAME_F_BOOLWEAVE)

struct rt_eval_wireframe_opts {
    uint32_t flags;
    int ncpu;
};

#define RT_EVAL_WIREFRAME_OPTS_INIT {RT_EVAL_WIREFRAME_F_DEFAULT, 1}

/**
 * Generate an evaluated wireframe for a database path.
 *
 * The caller owns both list heads and must initialize them with BU_LIST_INIT.
 * Output RT_VLIST_LINE_MOVE/RT_VLIST_LINE_DRAW chunks are appended to vhead
 * and allocated from vlfree.  Release the returned geometry with
 * RT_FREE_VLIST(vlfree, vhead), then clean vlfree using the caller's normal
 * vlist-free-list policy.
 */
RT_EXPORT extern int
rt_eval_wireframe(struct bu_list *vhead,
		  struct bu_list *vlfree,
		  struct db_i *dbip,
		  const char *path,
		  const struct bn_tol *tol,
		  const struct bg_tess_tol *ttol,
		  const struct rt_eval_wireframe_opts *opts);

__END_DECLS

#endif /* RT_EVAL_WIREFRAME_H */
