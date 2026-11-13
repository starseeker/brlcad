/*              D I S P L A Y _ B O U N D S . C P P
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

#include "raytrace.h"
#include "rt/display_bounds.h"

#include <cmath>
#include <unordered_set>

namespace {

struct display_bound {
    point_t minimum = {INFINITY, INFINITY, INFINITY};
    point_t maximum = {-INFINITY, -INFINITY, -INFINITY};
    bool empty = true;
};

struct display_bound_context {
    struct db_i *dbip = NULL;
    struct bu_vls *messages = NULL;
    struct bg_tess_tol tessellationTolerance = BG_TESS_TOL_INIT_TOL;
    struct bn_tol tolerance = BN_TOL_INIT_TOL;
    std::unordered_set<struct directory *> activeCombinations;
};

static bool
display_bound_finite(const display_bound &bounds)
{
    if (bounds.empty)
	return true;
    for (size_t axis = 0; axis < 3; ++axis) {
	if (!std::isfinite(bounds.minimum[axis]) ||
	    !std::isfinite(bounds.maximum[axis]) ||
	    bounds.minimum[axis] > bounds.maximum[axis])
	    return false;
    }
    return true;
}

static void
display_bound_union(display_bound &result, const display_bound &right)
{
    if (right.empty)
	return;
    if (result.empty) {
	VMOVE(result.minimum, right.minimum);
	VMOVE(result.maximum, right.maximum);
	result.empty = false;
	return;
    }
    VMIN(result.minimum, right.minimum);
    VMAX(result.maximum, right.maximum);
}

static void
display_bound_intersect(display_bound &result, const display_bound &right)
{
    if (result.empty || right.empty) {
	result = display_bound();
	return;
    }
    VMAX(result.minimum, right.minimum);
    VMIN(result.maximum, right.maximum);
    for (size_t axis = 0; axis < 3; ++axis) {
	if (result.minimum[axis] > result.maximum[axis]) {
	    result = display_bound();
	    return;
	}
    }
}

static bool display_bound_object(display_bound_context &context,
	struct directory *dp, const mat_t matrix, display_bound &bounds);

static bool
display_bound_tree(display_bound_context &context, const union tree *tree,
	const mat_t matrix, display_bound &bounds)
{
    if (!tree)
	return false;

    switch (tree->tr_op) {
	case OP_NOP:
	    bounds = display_bound();
	    return true;
	case OP_DB_LEAF: {
	    struct directory *dp = db_lookup(context.dbip,
		tree->tr_l.tl_name, LOOKUP_QUIET);
	    if (dp == RT_DIR_NULL)
		return false;
	    mat_t leafMatrix;
	    if (tree->tr_l.tl_mat)
		bn_mat_mul(leafMatrix, matrix, tree->tr_l.tl_mat);
	    else
		MAT_COPY(leafMatrix, matrix);
	    return display_bound_object(context, dp, leafMatrix, bounds);
	}
	case OP_UNION:
	case OP_XOR:
	case OP_INTERSECT: {
	    display_bound left;
	    display_bound right;
	    if (!display_bound_tree(context, tree->tr_b.tb_left, matrix,
		    left) ||
		!display_bound_tree(context, tree->tr_b.tb_right, matrix,
		    right))
		return false;
	    bounds = left;
	    if (tree->tr_op == OP_INTERSECT)
		display_bound_intersect(bounds, right);
	    else
		display_bound_union(bounds, right);
	    return display_bound_finite(bounds);
	}
	case OP_SUBTRACT:
	    /* A - B cannot extend beyond A.  Avoid importing B solely to prove
	     * the conservative display extent, matching rt_bound_tree. */
	    return display_bound_tree(context, tree->tr_b.tb_left, matrix,
		bounds);
	case OP_GUARD:
	case OP_XNOP:
	    return display_bound_tree(context, tree->tr_b.tb_left, matrix,
		bounds);
	default:
	    if (context.messages)
		bu_vls_printf(context.messages,
		    "unsupported Boolean operation %d while bounding\n",
		    tree->tr_op);
	    return false;
    }
}

static bool
display_bound_object(display_bound_context &context, struct directory *dp,
	const mat_t matrix, display_bound &bounds)
{
    if (!dp)
	return false;
    if (!(dp->d_flags & RT_DIR_COMB)) {
	mat_t primitiveMatrix;
	MAT_COPY(primitiveMatrix, matrix);
	point_t minimum;
	point_t maximum;
	if (rt_bound_instance(&minimum, &maximum, dp, context.dbip,
		&context.tessellationTolerance, &context.tolerance,
		&primitiveMatrix) < 0)
	    return false;
	VMOVE(bounds.minimum, minimum);
	VMOVE(bounds.maximum, maximum);
	bounds.empty = false;
	return display_bound_finite(bounds);
    }

    if (!context.activeCombinations.insert(dp).second) {
	if (context.messages)
	    bu_vls_printf(context.messages,
		"cyclic combination reference at %s\n",
		dp->d_namep ? dp->d_namep : "?");
	return false;
    }

    struct rt_db_internal internal;
    RT_DB_INTERNAL_INIT(&internal);
    const int imported = rt_db_get_internal(&internal, dp, context.dbip,
	NULL);
    bool success = false;
    if (imported >= 0 && internal.idb_type == ID_COMBINATION &&
	internal.idb_ptr) {
	const struct rt_comb_internal *combination =
	    static_cast<const struct rt_comb_internal *>(internal.idb_ptr);
	success = combination->tree && display_bound_tree(context,
	    combination->tree, matrix, bounds);
    }
    if (imported >= 0)
	rt_db_free_internal(&internal);
    context.activeCombinations.erase(dp);
    return success;
}

} /* namespace */

int
rt_display_bounds(struct bu_vls *messages, struct db_i *dbip, int argc,
	const char *argv[], point_t bounds_min, point_t bounds_max)
{
    if (!dbip || argc <= 0 || !argv || !bounds_min || !bounds_max)
	return BRLCAD_ERROR;

    display_bound_context context;
    context.dbip = dbip;
    context.messages = messages;
    display_bound aggregate;
    for (int index = 0; index < argc; ++index) {
	if (!argv[index] || !argv[index][0])
	    return BRLCAD_ERROR;
	struct db_full_path path;
	db_full_path_init(&path);
	if (db_string_to_path(&path, dbip, argv[index]) || !path.fp_len) {
	    if (messages)
		bu_vls_printf(messages, "invalid database path %s\n",
		    argv[index]);
	    db_free_full_path(&path);
	    return BRLCAD_ERROR;
	}
	mat_t pathMatrix;
	const bool matrixValid = db_path_to_mat(dbip, &path, pathMatrix,
	    static_cast<int>(path.fp_len - 1)) != 0;
	struct directory *target = DB_FULL_PATH_CUR_DIR(&path);
	display_bound pathBounds;
	const bool bounded = matrixValid && target &&
	    display_bound_object(context, target, pathMatrix, pathBounds);
	db_free_full_path(&path);
	if (!bounded)
	    return BRLCAD_ERROR;
	display_bound_union(aggregate, pathBounds);
    }

    if (aggregate.empty || !display_bound_finite(aggregate))
	return BRLCAD_ERROR;
    VMOVE(bounds_min, aggregate.minimum);
    VMOVE(bounds_max, aggregate.maximum);
    return BRLCAD_OK;
}
