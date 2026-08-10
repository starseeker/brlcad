/*                    L I N E _ L A Y E R . H
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
/** @file bg/line_layer.h
 *
 * Low-level color-keyed point/line command layer storage.
 */

#ifndef BG_LINE_LAYER_H
#define BG_LINE_LAYER_H

#include "common.h"

#include <stddef.h>
#include <stdio.h>

#include "vmath.h"
#include "bg/defines.h"

__BEGIN_DECLS

enum bg_geometry_line_command {
    BG_GEOMETRY_LINE_MOVE = 0,
    BG_GEOMETRY_LINE_DRAW = 1,
    BG_GEOMETRY_POINT_DRAW = 12
};

struct bg_line_layer_builder;
struct bg_line_layer;

BG_EXPORT extern struct bg_line_layer_builder *
bg_line_layer_builder_create(void);

BG_EXPORT extern void
bg_line_layer_builder_clear(struct bg_line_layer_builder *builder);

BG_EXPORT extern void
bg_line_layer_builder_free(struct bg_line_layer_builder *builder);

BG_EXPORT extern size_t
bg_line_layer_builder_layer_count(const struct bg_line_layer_builder *builder);

BG_EXPORT extern size_t
bg_line_layer_builder_point_count(const struct bg_line_layer_builder *builder);

BG_EXPORT extern const struct bg_line_layer *
bg_line_layer_builder_layer_at(const struct bg_line_layer_builder *builder,
			       size_t idx);

BG_EXPORT extern int
bg_line_layer_builder_add(struct bg_line_layer_builder *builder,
			  int r,
			  int g,
			  int b,
			  const point_t point,
			  int command);

BG_EXPORT extern struct bg_line_layer *
bg_line_layer_builder_find(struct bg_line_layer_builder *builder,
			   int r,
			   int g,
			   int b);

BG_EXPORT extern int
bg_line_layer_add(struct bg_line_layer *layer,
		  const point_t point,
		  int command);

BG_EXPORT extern const char *
bg_line_layer_name(const struct bg_line_layer *layer);

BG_EXPORT extern const point_t *
bg_line_layer_points(const struct bg_line_layer *layer);

BG_EXPORT extern const int *
bg_line_layer_commands(const struct bg_line_layer *layer);

BG_EXPORT extern int
bg_line_layer_color(const struct bg_line_layer *layer,
		    unsigned char *r,
		    unsigned char *g,
		    unsigned char *b);

BG_EXPORT extern size_t
bg_line_layer_point_count(const struct bg_line_layer *layer);

BG_EXPORT extern int
bg_line_layer_point_at(const struct bg_line_layer *layer,
		       size_t idx,
		       point_t point);

BG_EXPORT extern int
bg_line_layer_command_at(const struct bg_line_layer *layer,
			 size_t idx,
			 int *command);

BG_EXPORT extern int
bg_line_layer_write_plot3(FILE *plotfp,
			  const struct bg_line_layer *layer);

BG_EXPORT extern int
bg_line_layer_builder_write_plot3(FILE *plotfp,
				  const struct bg_line_layer_builder *builder);

__END_DECLS

#endif /* BG_LINE_LAYER_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
