/*              S K E T C H _ L E G A C Y _ B S G . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file rt/primitives/sketch_legacy_bsg.h
 *
 * Transitional adapters for the BSG-backed view polygon workflow.
 */

#ifndef RT_PRIMITIVES_SKETCH_LEGACY_BSG_H
#define RT_PRIMITIVES_SKETCH_LEGACY_BSG_H

#include "common.h"
#include "rt/defines.h"
#include "rt/primitives/sketch.h"
#include "rt/view.h"

__BEGIN_DECLS

RT_EXPORT extern rt_view_polygon_ref
db_sketch_to_view_polygon_ref(const char *sname, struct db_i *dbip, struct directory *dp, void *view_ctx);

RT_EXPORT extern rt_view_polygon_ref
db_sketch_to_view_polygon_scoped_ref(const char *sname, struct db_i *dbip, struct directory *dp, void *view_ctx, int local);

RT_EXPORT extern struct directory *
db_view_polygon_ref_to_sketch(struct db_i *dbip, const char *sname, rt_view_polygon_ref ref);

__END_DECLS

#endif /* RT_PRIMITIVES_SKETCH_LEGACY_BSG_H */
