/*                      D R A W _ T Y P E S . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged_draw_types_private.h
 *
 * Private draw-command migration and realization value types.  Public clients
 * use ged_scene_draw_mode, ged_scene_style, and ged_scene_realization_policy.
 */

#ifndef GED_DRAW_TYPES_PRIVATE_H
#define GED_DRAW_TYPES_PRIVATE_H

#include "vmath.h"
#include "ged/defines.h"

__BEGIN_DECLS

typedef enum ged_draw_stale_reason {
    GED_DRAW_STALE_NONE = 0,
    GED_DRAW_STALE_SOURCE_CHANGED,
    GED_DRAW_STALE_VIEW_INPUT_CHANGED,
    GED_DRAW_STALE_SETTINGS_CHANGED,
    GED_DRAW_STALE_FORCED,
    GED_DRAW_STALE_UPDATE_FAILED
} ged_draw_stale_reason;

typedef enum ged_draw_mode {
    GED_DRAW_MODE_WIRE         = 0,
    GED_DRAW_MODE_SHADED_BOTS  = 1,
    GED_DRAW_MODE_SHADED       = 2,
    GED_DRAW_MODE_EVAL_WIRE    = 3,
    GED_DRAW_MODE_HIDDEN_LINE  = 4,
    GED_DRAW_MODE_EVAL_POINTS  = 5
} ged_draw_mode;

struct ged_draw_appearance_settings {
    int draw_mode;
    int mixed_modes;
    /* Legacy draw-option opacity: 0 is clear and 1 is opaque. */
    fastf_t transparency;
    int color_override;
    unsigned char color[3];
    int s_line_width;
    fastf_t s_arrow_tip_length;
    fastf_t s_arrow_tip_width;
    int draw_solid_lines_only;
    int draw_non_subtract_only;
    int strict_fallback;
    int defer_leaf_expansion;
};

#define GED_DRAW_APPEARANCE_SETTINGS_INIT { \
    GED_DRAW_MODE_WIRE, 0, 1.0, 0, {255, 0, 0}, 1, 0.0, 0.0, 0, 0, 0, 0 }

__END_DECLS

#endif /* GED_DRAW_TYPES_PRIVATE_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
