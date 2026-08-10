/*                          V I E W . H
 * BRL-CAD
 *
 * Copyright (c) 1993-2026 United States Government as represented by
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
/** @file rt/view.h
 *
 */

#ifndef RT_VIEW_H
#define RT_VIEW_H

#include "common.h"
#include "vmath.h"
#include "bg/polygon.h"
#include "bu/color.h"
#include "bu/list.h"
#include "bu/hash.h"
#include "bu/ptbl.h"
#include "bv/view.h"
#include "bn/tol.h"
#include "rt/defines.h"

__BEGIN_DECLS

struct db_i;
struct directory;
struct bu_vls;

/* Historical normalized view-coordinate range. */
#define RT_VIEW_MAX 2047.0
#define RT_VIEW_MIN -2048.0
#define RT_VIEW_RANGE 4095.0
#define RT_INV_VIEW 0.00048828125
#define RT_INV_4096 0.000244140625
#define RT_VIEW_MIN_SIZE 0.0001
#define RT_VIEW_MIN_SCALE 0.00005
#define RT_VIEW_AUTOVIEW_SCALE_DEFAULT -1
#define RT_VIEW_ADJUST_IDLE   0x000ULL
#define RT_VIEW_ADJUST_ROT    0x001ULL
#define RT_VIEW_ADJUST_TRANS  0x002ULL
#define RT_VIEW_ADJUST_SCALE  0x004ULL
#define RT_VIEW_ADJUST_CENTER 0x008ULL
#define RT_VIEW_ADJUST_CON_X  0x010ULL
#define RT_VIEW_ADJUST_CON_Y  0x020ULL
#define RT_VIEW_ADJUST_CON_Z  0x040ULL
#define RT_VIEW_ADJUST_CON_GRID 0x080ULL
#define RT_VIEW_ADJUST_CON_LINES 0x100ULL
#define RT_VIEW_REFRESH_VIEW        0x00000001u
#define RT_VIEW_REFRESH_DRAW        0x00000002u
#define RT_VIEW_REFRESH_EDIT        0x00000004u
#define RT_VIEW_REFRESH_FRAMEBUFFER 0x00000008u
#define RT_VIEW_REFRESH_OVERLAY     0x00000010u
#define RT_VIEW_REFRESH_FORCE       0x80000000u
#define RT_VIEW_REFRESH_ALL         0xffffffffu
#define RT_VIEW_CLEAR_DB    0x01
#define RT_VIEW_CLEAR_VIEW  0x02
#define RT_VIEW_CLEAR_LOCAL 0x04

/**
 * Specifies a subset of a primitive's geometry as the target for an
 * operation.
 *
 * TODO: This structure is tentative and subject to change or removal
 *       without notice.
 */
struct rt_selection {
    void *obj; /**< @brief primitive-specific selection object */
};

/**
 * TODO: This structure is tentative and subject to change or removal
 *       without notice.
 */
struct rt_selection_set {
    struct bu_ptbl selections; /**< @brief holds struct rt_selection */

    /** selection-object-specific routine that will free all memory
     *  associated with any of the stored selections
     */
    void (*free_selection)(struct rt_selection *);
};

/**
 * Stores selections associated with an object. There is an entry in
 * the selections table for each kind of selection (e.g. "active",
 * "option"). The table entries are sets to allow more than one
 * selection of the same type (e.g. multiple "option" selections).
 *
 * TODO: This structure is tentative and subject to change or removal
 *       without notice.
 */
struct rt_object_selections {
    /** selection type -> struct rt_selection_set */
    struct bu_hash_tbl *sets;
};

/**
 * Analogous to a database query. Specifies how to filter and sort the
 * selectable components of a primitive in order to find the most
 * relevant selections for a particular application.
 *
 * TODO: This structure is tentative and subject to change or removal
 *       without notice.
 */
struct rt_selection_query {
    point_t start;     /**< @brief start point of query ray */
    vect_t dir;        /**< @brief direction of query ray */

#define RT_SORT_UNSORTED         0
#define RT_SORT_CLOSEST_TO_START 1
    int sorting;
};

/**
 * Parameters of a translation applied to a selection.
 *
 * TODO: This structure is tentative and subject to change or removal
 *       without notice.
 */
struct rt_selection_translation {
    fastf_t dx;
    fastf_t dy;
    fastf_t dz;
};

/**
 * Describes an operation that can be applied to a selection.
 *
 * TODO: This structure is tentative and subject to change or removal
 *       without notice.
 */
struct rt_selection_operation {
#define RT_SELECTION_NOP         0
#define RT_SELECTION_TRANSLATION 1
    int type;
    union {
	struct rt_selection_translation tran;
    } parameters;
};

__END_DECLS

#endif /* RT_VIEW_H */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
