/*                         H U D . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */

#ifndef MGED_HUD_H
#define MGED_HUD_H

#include "common.h"
#include "ged/draw.h"

__BEGIN_DECLS

struct mged_hud_line_layer {
    int color[3];
    int line_width;
    int line_style;
    point_t *points;
    int *commands;
    size_t point_count;
    size_t point_capacity;
    char name[32];
};

struct mged_hud_builder {
    struct ged_view_context *view_ctx;
    char *prefix;
    struct ged_annotation_label *labels;
    char **label_text;
    size_t label_count;
    size_t label_capacity;
    struct mged_hud_line_layer *layers;
    size_t layer_count;
    size_t layer_capacity;
    int color[3];
    fastf_t font_size;
    int line_width;
    int line_style;
};

void mged_hud_builder_init(struct mged_hud_builder *builder,
	struct ged_view_context *view_ctx, const char *prefix);
void mged_hud_builder_free(struct mged_hud_builder *builder);
void mged_hud_color_set(struct mged_hud_builder *builder,
	const int color[3]);
void mged_hud_font_size_set(struct mged_hud_builder *builder,
	fastf_t font_size);
void mged_hud_line_style_set(struct mged_hud_builder *builder,
	int line_width, int line_style);
int mged_hud_line_add(struct mged_hud_builder *builder,
	fastf_t x1, fastf_t y1, fastf_t x2, fastf_t y2);
int mged_hud_label_add(struct mged_hud_builder *builder,
	const char *text, fastf_t x, fastf_t y, fastf_t font_size,
	int anchor);
int mged_hud_builder_publish(struct mged_hud_builder *builder);

__END_DECLS

#endif /* MGED_HUD_H */
