/*                 D I S P L A Y _ B O U N D S . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file rt/display_bounds.h
 *
 * Read-only, Boolean-aware bounds intended for interactive display setup.
 * Unlike rt_obj_bounds(), this path does not construct regions, soltabs, or
 * raytracing acceleration structures.  It evaluates the authored combination
 * expression directly from primitive bounding boxes.
 */

#ifndef RT_DISPLAY_BOUNDS_H
#define RT_DISPLAY_BOUNDS_H

#include "common.h"

#include "vmath.h"
#include "rt/defines.h"

__BEGIN_DECLS

struct bu_vls;
struct db_i;

/**
 * Compute the conservative semantic bounds of one or more database paths.
 *
 * Union and xor combine both operands, intersection clips their boxes, and
 * subtraction retains the left operand's box.  The result is therefore safe
 * for view framing but is not claimed to be a tight bound of evaluated CSG.
 * Primitive imports and combination traversal are local to this call; the
 * database must remain read-only until it returns.
 *
 * Returns BRLCAD_OK on success and BRLCAD_ERROR when no finite bound can be
 * established.  Optional diagnostics are appended to @p messages.
 */
RT_EXPORT extern int
rt_display_bounds(struct bu_vls *messages,
	struct db_i *dbip,
	int argc,
	const char *argv[],
	point_t bounds_min,
	point_t bounds_max);

__END_DECLS

#endif /* RT_DISPLAY_BOUNDS_H */
