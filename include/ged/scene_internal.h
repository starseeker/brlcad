/*                 S C E N E _ I N T E R N A L . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file ged/scene_internal.h
 *
 * Uninstalled migration surface for in-tree clients that still navigate
 * realized records or animate a transient vlist directly.  This is not a
 * semantic scene API.  New code must use edit previews, semantic identities,
 * or typed scene operations.
 */

#ifndef GED_SCENE_INTERNAL_H
#define GED_SCENE_INTERNAL_H

#include "common.h"

#include "ged/defines.h"
#include "ged/draw_types.h"
#include "ged/scene_record_api_internal.h"
#include "vmath.h"

__BEGIN_DECLS

GED_EXPORT extern int
ged_scene_internal_shape_translate_geometry(struct ged *gedp,
	ged_draw_shape_ref ref, const vect_t translation);

GED_EXPORT extern int
ged_scene_internal_shape_set_center(struct ged *gedp,
	ged_draw_shape_ref ref, const point_t center);

GED_EXPORT extern int
ged_scene_internal_shape_geometry_clear(struct ged *gedp,
	ged_draw_shape_ref ref);

GED_EXPORT extern int
ged_scene_internal_shape_bounds_update(struct ged *gedp,
	ged_draw_shape_ref ref, int *bad_command);

__END_DECLS

#endif /* GED_SCENE_INTERNAL_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
