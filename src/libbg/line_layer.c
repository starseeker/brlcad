/*                    L I N E _ L A Y E R . C
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
/** @file libbg/line_layer.c
 *
 * Typed multi-color line layer builder.
 */

#include "common.h"

#include <string.h>

#include "bu/malloc.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "bg/line_layer.h"
#include "bg/plot3.h"


struct bg_line_layer {
    long rgb;
    char *name;
    point_t *points;
    int *commands;
    size_t point_count;
    size_t point_capacity;
};

struct bg_line_layer_builder {
    struct bg_line_layer *layers;
    size_t layer_count;
    size_t layer_capacity;
};


static long
_line_layer_rgb(int r, int g, int b)
{
    return ((r & 0xFF) << 16) | ((g & 0xFF) << 8) | (b & 0xFF);
}


static int
_line_layer_geometry_cmd_valid(int command)
{
    switch (command) {
	case BG_GEOMETRY_LINE_MOVE:
	case BG_GEOMETRY_LINE_DRAW:
	case BG_GEOMETRY_POINT_DRAW:
	    return 1;
	default:
	    break;
    }
    return 0;
}


static void
_line_layer_free(struct bg_line_layer *layer)
{
    if (!layer)
	return;

    if (layer->name)
	bu_free(layer->name, "bg line layer name");
    if (layer->points)
	bu_free(layer->points, "bg line layer points");
    if (layer->commands)
	bu_free(layer->commands, "bg line layer commands");
    memset(layer, 0, sizeof(*layer));
}


static int
_line_layer_reserve_points(struct bg_line_layer *layer, size_t needed)
{
    if (!layer)
	return 0;
    if (needed <= layer->point_capacity)
	return 1;

    size_t capacity = layer->point_capacity ? layer->point_capacity : 64;
    while (capacity < needed)
	capacity *= 2;

    layer->points = (point_t *)bu_realloc(layer->points, capacity * sizeof(point_t),
	    "bg line layer points");
    layer->commands = (int *)bu_realloc(layer->commands, capacity * sizeof(int),
	    "bg line layer commands");
    layer->point_capacity = capacity;
    return 1;
}


static struct bg_line_layer *
_line_layer_find_or_create(struct bg_line_layer_builder *builder, int r, int g, int b)
{
    if (!builder)
	return NULL;

    long rgb = _line_layer_rgb(r, g, b);
    for (size_t i = 0; i < builder->layer_count; i++) {
	if (builder->layers[i].rgb == rgb)
	    return &builder->layers[i];
    }

    if (builder->layer_count == builder->layer_capacity) {
	size_t capacity = builder->layer_capacity ? builder->layer_capacity * 2 : 8;
	builder->layers = (struct bg_line_layer *)bu_realloc(builder->layers,
		capacity * sizeof(struct bg_line_layer),
		"bg line layer builder layers");
	for (size_t i = builder->layer_capacity; i < capacity; i++)
	    memset(&builder->layers[i], 0, sizeof(builder->layers[i]));
	builder->layer_capacity = capacity;
    }

    struct bg_line_layer *layer = &builder->layers[builder->layer_count++];
    memset(layer, 0, sizeof(*layer));
    layer->rgb = rgb;

    struct bu_vls lname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&lname, "rgb_%03ld_%03ld_%03ld",
	    (rgb >> 16) & 0xFF, (rgb >> 8) & 0xFF, rgb & 0xFF);
    layer->name = bu_strdup(bu_vls_cstr(&lname));
    bu_vls_free(&lname);

    return layer;
}


struct bg_line_layer_builder *
bg_line_layer_builder_create(void)
{
    struct bg_line_layer_builder *builder = NULL;
    BU_ALLOC(builder, struct bg_line_layer_builder);
    return builder;
}


void
bg_line_layer_builder_clear(struct bg_line_layer_builder *builder)
{
    if (!builder)
	return;

    for (size_t i = 0; i < builder->layer_count; i++)
	_line_layer_free(&builder->layers[i]);
    builder->layer_count = 0;
}


void
bg_line_layer_builder_free(struct bg_line_layer_builder *builder)
{
    if (!builder)
	return;

    bg_line_layer_builder_clear(builder);
    if (builder->layers)
	bu_free(builder->layers, "bg line layer builder layers");
    bu_free(builder, "bg line layer builder");
}


size_t
bg_line_layer_builder_layer_count(const struct bg_line_layer_builder *builder)
{
    return builder ? builder->layer_count : 0;
}


size_t
bg_line_layer_builder_point_count(const struct bg_line_layer_builder *builder)
{
    size_t count = 0;
    if (!builder)
	return 0;

    for (size_t i = 0; i < builder->layer_count; i++)
	count += builder->layers[i].point_count;
    return count;
}


const struct bg_line_layer *
bg_line_layer_builder_layer_at(const struct bg_line_layer_builder *builder,
			       size_t idx)
{
    if (!builder || idx >= builder->layer_count)
	return NULL;
    return &builder->layers[idx];
}


struct bg_line_layer *
bg_line_layer_builder_find(struct bg_line_layer_builder *builder,
			   int r,
			   int g,
			   int b)
{
    return _line_layer_find_or_create(builder, r, g, b);
}


int
bg_line_layer_add(struct bg_line_layer *layer,
		  const point_t point,
		  int command)
{
    if (!layer || !point || !_line_layer_geometry_cmd_valid(command))
	return 0;

    if (!_line_layer_reserve_points(layer, layer->point_count + 1))
	return 0;

    VMOVE(layer->points[layer->point_count], point);
    layer->commands[layer->point_count] = command;
    layer->point_count++;
    return 1;
}


int
bg_line_layer_builder_add(struct bg_line_layer_builder *builder,
			  int r,
			  int g,
			  int b,
			  const point_t point,
			  int command)
{
    struct bg_line_layer *layer = bg_line_layer_builder_find(builder, r, g, b);
    return bg_line_layer_add(layer, point, command);
}


const char *
bg_line_layer_name(const struct bg_line_layer *layer)
{
    return layer ? layer->name : NULL;
}


const point_t *
bg_line_layer_points(const struct bg_line_layer *layer)
{
    return layer ? (const point_t *)layer->points : NULL;
}


const int *
bg_line_layer_commands(const struct bg_line_layer *layer)
{
    return layer ? (const int *)layer->commands : NULL;
}


int
bg_line_layer_color(const struct bg_line_layer *layer,
		    unsigned char *r,
		    unsigned char *g,
		    unsigned char *b)
{
    if (!layer)
	return 0;
    if (r)
	*r = (unsigned char)((layer->rgb >> 16) & 0xFF);
    if (g)
	*g = (unsigned char)((layer->rgb >> 8) & 0xFF);
    if (b)
	*b = (unsigned char)(layer->rgb & 0xFF);
    return 1;
}


size_t
bg_line_layer_point_count(const struct bg_line_layer *layer)
{
    return layer ? layer->point_count : 0;
}


int
bg_line_layer_point_at(const struct bg_line_layer *layer,
		       size_t idx,
		       point_t point)
{
    if (!layer || !point || idx >= layer->point_count)
	return 0;
    VMOVE(point, layer->points[idx]);
    return 1;
}


int
bg_line_layer_command_at(const struct bg_line_layer *layer,
			 size_t idx,
			 int *command)
{
    if (!layer || !command || idx >= layer->point_count)
	return 0;
    *command = layer->commands[idx];
    return 1;
}


int
bg_line_layer_write_plot3(FILE *plotfp, const struct bg_line_layer *layer)
{
    if (!plotfp || !layer)
	return 0;

    unsigned char r = 0;
    unsigned char g = 0;
    unsigned char b = 0;
    if (!bg_line_layer_color(layer, &r, &g, &b))
	return 0;

    pl_color(plotfp, r, g, b);
    for (size_t i = 0; i < layer->point_count; i++) {
	switch (layer->commands[i]) {
	    case BG_GEOMETRY_LINE_MOVE:
		pdv_3move(plotfp, layer->points[i]);
		break;
	    case BG_GEOMETRY_LINE_DRAW:
		pdv_3cont(plotfp, layer->points[i]);
		break;
	    case BG_GEOMETRY_POINT_DRAW:
		pdv_3point(plotfp, layer->points[i]);
		break;
	    default:
		return 0;
	}
    }

    return ferror(plotfp) ? 0 : 1;
}


int
bg_line_layer_builder_write_plot3(FILE *plotfp,
				  const struct bg_line_layer_builder *builder)
{
    if (!plotfp || !builder)
	return 0;

    for (size_t i = 0; i < builder->layer_count; i++) {
	if (!bg_line_layer_write_plot3(plotfp, &builder->layers[i]))
	    return 0;
    }

    return ferror(plotfp) ? 0 : 1;
}


/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
