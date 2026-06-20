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
/** @file libbsg/line_layer.c
 *
 * BSG compatibility wrappers for the neutral bg_line_layer builder.
 */

#include "common.h"

#include "bu/malloc.h"
#include "bg/line_layer.h"
#include "bsg/feature.h"
#include "bsg/geometry.h"


static int
_bsg_line_layer_command_to_bg(int command, int *bg_command)
{
    switch (command) {
	case BSG_GEOMETRY_LINE_MOVE:
	    *bg_command = BG_GEOMETRY_LINE_MOVE;
	    return 1;
	case BSG_GEOMETRY_LINE_DRAW:
	    *bg_command = BG_GEOMETRY_LINE_DRAW;
	    return 1;
	case BSG_GEOMETRY_POINT_DRAW:
	    *bg_command = BG_GEOMETRY_POINT_DRAW;
	    return 1;
	default:
	    break;
    }
    return 0;
}


static int
_bsg_line_layer_command_from_bg(int bg_command, int *command)
{
    switch (bg_command) {
	case BG_GEOMETRY_LINE_MOVE:
	    *command = BSG_GEOMETRY_LINE_MOVE;
	    return 1;
	case BG_GEOMETRY_LINE_DRAW:
	    *command = BSG_GEOMETRY_LINE_DRAW;
	    return 1;
	case BG_GEOMETRY_POINT_DRAW:
	    *command = BSG_GEOMETRY_POINT_DRAW;
	    return 1;
	default:
	    break;
    }
    return 0;
}


struct bsg_line_layer_builder *
bsg_line_layer_builder_create(void)
{
    return (struct bsg_line_layer_builder *)bg_line_layer_builder_create();
}


void
bsg_line_layer_builder_clear(struct bsg_line_layer_builder *builder)
{
    bg_line_layer_builder_clear((struct bg_line_layer_builder *)builder);
}


void
bsg_line_layer_builder_free(struct bsg_line_layer_builder *builder)
{
    bg_line_layer_builder_free((struct bg_line_layer_builder *)builder);
}


size_t
bsg_line_layer_builder_layer_count(const struct bsg_line_layer_builder *builder)
{
    return bg_line_layer_builder_layer_count((const struct bg_line_layer_builder *)builder);
}


size_t
bsg_line_layer_builder_point_count(const struct bsg_line_layer_builder *builder)
{
    return bg_line_layer_builder_point_count((const struct bg_line_layer_builder *)builder);
}


const struct bsg_line_layer *
bsg_line_layer_builder_layer_at(const struct bsg_line_layer_builder *builder,
				size_t idx)
{
    return (const struct bsg_line_layer *)bg_line_layer_builder_layer_at(
	    (const struct bg_line_layer_builder *)builder, idx);
}


struct bsg_line_layer *
bsg_line_layer_builder_find(struct bsg_line_layer_builder *builder,
			    int r,
			    int g,
			    int b)
{
    return (struct bsg_line_layer *)bg_line_layer_builder_find(
	    (struct bg_line_layer_builder *)builder, r, g, b);
}


int
bsg_line_layer_add(struct bsg_line_layer *layer,
		   const point_t point,
		   int command)
{
    int bg_command = -1;
    if (!_bsg_line_layer_command_to_bg(command, &bg_command))
	return 0;

    return bg_line_layer_add((struct bg_line_layer *)layer, point, bg_command);
}


int
bsg_line_layer_builder_add(struct bsg_line_layer_builder *builder,
			   int r,
			   int g,
			   int b,
			   const point_t point,
			   int command)
{
    int bg_command = -1;
    if (!_bsg_line_layer_command_to_bg(command, &bg_command))
	return 0;

    return bg_line_layer_builder_add((struct bg_line_layer_builder *)builder,
	    r, g, b, point, bg_command);
}


int
bsg_line_layer_color(const struct bsg_line_layer *layer,
		     unsigned char *r,
		     unsigned char *g,
		     unsigned char *b)
{
    return bg_line_layer_color((const struct bg_line_layer *)layer, r, g, b);
}


size_t
bsg_line_layer_point_count(const struct bsg_line_layer *layer)
{
    return bg_line_layer_point_count((const struct bg_line_layer *)layer);
}


int
bsg_line_layer_point_at(const struct bsg_line_layer *layer,
			size_t idx,
			point_t point)
{
    return bg_line_layer_point_at((const struct bg_line_layer *)layer, idx, point);
}


int
bsg_line_layer_command_at(const struct bsg_line_layer *layer,
			  size_t idx,
			  int *command)
{
    int bg_command = -1;
    if (!bg_line_layer_command_at((const struct bg_line_layer *)layer, idx, &bg_command))
	return 0;

    return _bsg_line_layer_command_from_bg(bg_command, command);
}


bsg_feature_ref
bsg_feature_replace_line_layer_builder(struct bsg_view *v,
				       const char *name,
				       int local,
				       const struct bsg_line_layer_builder *builder,
				       const struct bsg_feature_style *style)
{
    bsg_feature_ref null_ref = BSG_FEATURE_REF_NULL_INIT;
    if (!v || !name)
	return null_ref;

    const struct bg_line_layer_builder *bg_builder = (const struct bg_line_layer_builder *)builder;
    size_t layer_count = bg_line_layer_builder_layer_count(bg_builder);
    if (!builder || !layer_count)
	return bsg_feature_replace_line_layers(v, name, local, NULL, 0, style);

    struct bsg_feature_line_layer *layers = (struct bsg_feature_line_layer *)bu_calloc(
	    layer_count, sizeof(struct bsg_feature_line_layer), "bsg feature line layers");

    for (size_t i = 0; i < layer_count; i++) {
	const struct bg_line_layer *src = bg_line_layer_builder_layer_at(bg_builder, i);
	struct bsg_feature_line_layer init = BSG_FEATURE_LINE_LAYER_INIT;
	unsigned char r = 0;
	unsigned char g = 0;
	unsigned char b = 0;
	layers[i] = init;
	layers[i].name = NULL;
	layers[i].points = bg_line_layer_points(src);
	layers[i].commands = bg_line_layer_commands(src);
	layers[i].point_count = bg_line_layer_point_count(src);
	bg_line_layer_color(src, &r, &g, &b);
	layers[i].style.color_valid = 1;
	layers[i].style.color[0] = r;
	layers[i].style.color[1] = g;
	layers[i].style.color[2] = b;
    }

    bsg_feature_ref ref = bsg_feature_replace_line_layers(v, name, local,
	    (const struct bsg_feature_line_layer *)layers, layer_count, style);
    bu_free(layers, "bsg feature line layers");
    return ref;
}
