/*                        A N N O T . H
 * BRL-CAD
 *
 * Copyright (c) 2017-2026 United States Government as represented by
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
/** @addtogroup rt_annotation */
/** @{ */
/** @file rt/primitives/annot.h */

#ifndef RT_PRIMITIVES_ANNOT_H
#define RT_PRIMITIVES_ANNOT_H

#include "common.h"
#include "vmath.h"
#include "bu/list.h"
#include "bu/vls.h"
#include "bn/tol.h"
#include "rt/defines.h"

/** Nominal display density for stored screen-plane annotation coordinates.
 * Explicit authoring DPI is already reflected in stored text dimensions. */
#define RT_ANNOT_SCREEN_DPI 96.0
#define RT_ANNOT_DISPLAY_PIXELS_PER_MM (RT_ANNOT_SCREEN_DPI / 25.4)

__BEGIN_DECLS

/** Primitive ranges emitted for one annotation segment.  Line and triangle
 * indices refer to their independent output streams, after fill backgrounds
 * have been ordered ahead of strokes.  The style pointer remains owned by
 * the annotation and is valid only for the duration of the call. */
struct rt_annot_plot_range {
    size_t segment;
    size_t first_line;
    size_t line_count;
    size_t first_triangle;
    size_t triangle_count;
    const struct rt_annot_seg_style *style;
};

typedef void (*rt_annot_plot_range_callback)(
    const struct rt_annot_plot_range *range, void *data);

RT_EXPORT extern struct rt_annot_internal *rt_copy_annot(const struct rt_annot_internal *annot_ip);

/** Validate annotation topology, model-space placement, and optional segment
 * presentation data.  Returns zero when valid. */
RT_EXPORT extern int rt_annot_validate(const struct rt_annot_internal *annot_ip,
	struct bu_vls *messages);

/** Plot an annotation and report the source style for each emitted primitive
 * range.  This leaves the generic vlist ABI unchanged while retained
 * consumers preserve segment provenance. */
RT_EXPORT extern int rt_annot_plot_with_styles(struct bu_list *vhead,
	struct rt_db_internal *ip, const struct bg_tess_tol *ttol,
	rt_annot_plot_range_callback callback, void *data);

__END_DECLS

/** @} */

#endif /* RT_PRIMITIVES_ANNOT_H */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
