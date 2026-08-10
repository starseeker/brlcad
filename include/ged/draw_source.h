/* D R A W _ S O U R C E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged/draw_source.h */

#ifndef GED_DRAW_SOURCE_H
#define GED_DRAW_SOURCE_H

#include "ged/draw_types.h"

__BEGIN_DECLS

GED_EXPORT extern int
ged_draw_source_mark_changed(struct ged *gedp,
			      const char *path,
			      ged_draw_stale_reason reason);

GED_EXPORT extern int
ged_draw_source_apply_update(struct ged *gedp,
			       const char *path,
			       int removed,
			       int redraw);

GED_EXPORT extern const char *
ged_draw_source_stale_reason_name(ged_draw_stale_reason reason);

/**
 * Return non-zero when @p gedp has an initialized draw scene.
 */

GED_EXPORT extern int
ged_draw_source_snapshot(
    struct ged *gedp,
    ged_draw_shape_ref ref,
    struct ged_draw_shape_source_snapshot *out);


__END_DECLS

#endif /* GED_DRAW_SOURCE_H */
