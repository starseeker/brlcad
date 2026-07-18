/*                           H U D . C
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

#include <limits.h>
#include <string.h>

#include "bu/malloc.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "vmath.h"

#include "./hud.h"

static int
mged_hud_capacity_grow(void **storage, size_t *capacity, size_t count,
	size_t item_size, const char *label)
{
    if (count < *capacity)
	return 1;
    if (item_size == 0 || *capacity > SIZE_MAX / 2 ||
	(*capacity ? *capacity * 2 : 8) > SIZE_MAX / item_size)
	return 0;
    size_t next = *capacity ? *capacity * 2 : 8;
    *storage = bu_realloc(*storage, next * item_size, label);
    *capacity = next;
    return 1;
}

void
mged_hud_builder_init(struct mged_hud_builder *builder,
	struct ged_view_context *view_ctx, const char *prefix)
{
    if (!builder)
	return;
    memset(builder, 0, sizeof(*builder));
    builder->view_ctx = view_ctx;
    builder->prefix = bu_strdup(prefix && prefix[0] ? prefix : "_faceplate/mged");
    VSET(builder->color, 255, 255, 255);
    builder->font_size = 20.0;
    builder->line_width = 1;
}

void
mged_hud_builder_free(struct mged_hud_builder *builder)
{
    if (!builder)
	return;
    for (size_t i = 0; i < builder->label_count; i++)
	bu_free(builder->label_text[i], "MGED HUD label text");
    for (size_t i = 0; i < builder->layer_count; i++) {
	bu_free(builder->layers[i].points, "MGED HUD line points");
	bu_free(builder->layers[i].commands, "MGED HUD line commands");
    }
    bu_free(builder->labels, "MGED HUD labels");
    bu_free(builder->label_text, "MGED HUD label strings");
    bu_free(builder->layers, "MGED HUD line layers");
    bu_free(builder->prefix, "MGED HUD prefix");
    memset(builder, 0, sizeof(*builder));
}

void
mged_hud_color_set(struct mged_hud_builder *builder, const int color[3])
{
    if (!builder || !color)
	return;
    VMOVE(builder->color, color);
}

void
mged_hud_font_size_set(struct mged_hud_builder *builder, fastf_t font_size)
{
    if (!builder || font_size <= 0.0)
	return;
    builder->font_size = font_size;
}

void
mged_hud_line_style_set(struct mged_hud_builder *builder,
	int line_width, int line_style)
{
    if (!builder)
	return;
    builder->line_width = line_width > 0 ? line_width : 1;
    builder->line_style = line_style ? 1 : 0;
}

static struct mged_hud_line_layer *
mged_hud_line_layer_get(struct mged_hud_builder *builder)
{
    for (size_t i = 0; i < builder->layer_count; i++) {
	struct mged_hud_line_layer *layer = &builder->layers[i];
	if (VEQUAL(layer->color, builder->color) &&
	    layer->line_width == builder->line_width &&
	    layer->line_style == builder->line_style)
	    return layer;
    }
    if (!mged_hud_capacity_grow((void **)&builder->layers,
	    &builder->layer_capacity, builder->layer_count,
	    sizeof(struct mged_hud_line_layer), "MGED HUD line layers"))
	return NULL;
    struct mged_hud_line_layer *layer =
	&builder->layers[builder->layer_count++];
    memset(layer, 0, sizeof(*layer));
    VMOVE(layer->color, builder->color);
    layer->line_width = builder->line_width;
    layer->line_style = builder->line_style;
    snprintf(layer->name, sizeof(layer->name), "style-%zu",
	    builder->layer_count - 1);
    return layer;
}

int
mged_hud_line_add(struct mged_hud_builder *builder,
	fastf_t x1, fastf_t y1, fastf_t x2, fastf_t y2)
{
    if (!builder || !builder->view_ctx)
	return 0;
    struct mged_hud_line_layer *layer = mged_hud_line_layer_get(builder);
    if (!layer)
	return 0;
    if (!mged_hud_capacity_grow((void **)&layer->points,
	    &layer->point_capacity, layer->point_count + 1,
	    sizeof(point_t), "MGED HUD line points"))
	return 0;
    size_t command_capacity = layer->point_capacity;
    layer->commands = bu_realloc(layer->commands,
	command_capacity * sizeof(int), "MGED HUD line commands");
    VSET(layer->points[layer->point_count], x1, y1, 0.0);
    layer->commands[layer->point_count++] = GED_DRAW_VIEW_LINE_MOVE;
    VSET(layer->points[layer->point_count], x2, y2, 0.0);
    layer->commands[layer->point_count++] = GED_DRAW_VIEW_LINE_DRAW;
    return 1;
}

int
mged_hud_label_add(struct mged_hud_builder *builder,
	const char *text, fastf_t x, fastf_t y, fastf_t font_size, int anchor)
{
    if (!builder || !builder->view_ctx || !text || !text[0])
	return 0;
    if (!mged_hud_capacity_grow((void **)&builder->labels,
	    &builder->label_capacity, builder->label_count,
	    sizeof(struct ged_annotation_label), "MGED HUD labels"))
	return 0;
    builder->label_text = bu_realloc(builder->label_text,
	builder->label_capacity * sizeof(char *), "MGED HUD label strings");
    size_t index = builder->label_count++;
    builder->label_text[index] = bu_strdup(text);
    struct ged_annotation_label label = GED_ANNOTATION_LABEL_INIT;
    label.text = builder->label_text[index];
    VSET(label.point, x, y, 0.0);
    label.color_valid = 1;
    VMOVE(label.color, builder->color);
    label.font_size = font_size > 0.0 ? font_size : builder->font_size;
    label.anchor = anchor;
    builder->labels[index] = label;
    return 1;
}

int
mged_hud_builder_publish(struct mged_hud_builder *builder)
{
    if (!builder || !builder->view_ctx || !builder->prefix)
	return 0;

    struct bu_vls name = BU_VLS_INIT_ZERO;
    bu_vls_printf(&name, "%s/labels", builder->prefix);
    int labels_ok = ged_annotation_hud_labels_replace(
	builder->view_ctx, bu_vls_cstr(&name), builder->labels,
	builder->label_count, NULL);

    struct ged_annotation_line_layer *layers = NULL;
    if (builder->layer_count)
	layers = (struct ged_annotation_line_layer *)bu_calloc(
	    builder->layer_count, sizeof(*layers), "MGED HUD publish layers");
    for (size_t i = 0; i < builder->layer_count; i++) {
	layers[i].name = builder->layers[i].name;
	layers[i].points = (const point_t *)builder->layers[i].points;
	layers[i].commands = builder->layers[i].commands;
	layers[i].point_count = builder->layers[i].point_count;
	layers[i].style = (struct ged_view_feature_style)
	    GED_VIEW_FEATURE_STYLE_INIT;
	layers[i].style.visible = 1;
	layers[i].style.selectable = 0;
	layers[i].style.color_valid = 1;
	VMOVE(layers[i].style.color, builder->layers[i].color);
	layers[i].style.line_width = builder->layers[i].line_width;
	layers[i].style.line_style = builder->layers[i].line_style;
    }
    bu_vls_sprintf(&name, "%s/lines", builder->prefix);
    int lines_ok = ged_annotation_hud_line_layers_replace(
	builder->view_ctx, bu_vls_cstr(&name), layers,
	builder->layer_count, NULL);
    bu_free(layers, "MGED HUD publish layers");
    bu_vls_free(&name);
    return labels_ok && lines_ok;
}
