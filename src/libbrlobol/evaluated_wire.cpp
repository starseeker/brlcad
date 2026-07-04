/*              E V A L U A T E D _ W I R E . C P P
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
/** @addtogroup libbrlobol */
/** @{ */
/** @file libbrlobol/evaluated_wire.cpp
 *
 * libbrlobol compatibility adapter for librt evaluated-wireframe output.
 *
 */
/** @} */

#include "common.h"

#include <stddef.h>

#include "bg/line_layer.h"
#include "bg/vlist.h"
#include "brlobol/evaluated_wire.h"
#include "bu/list.h"
#include "bu/malloc.h"
#include "rt/eval_wireframe.h"
#include "rt/vlist.h"


static int
brlobol_vlist_line_count(const struct bu_list *vhead, size_t *count_out)
{
    rt_vlist *vp = NULL;

    if (!vhead || !count_out)
	return BRLCAD_ERROR;

    *count_out = 0;
    BU_LIST_EACH(vhead, vp, rt_vlist) {
	for (size_t i = 0; i < vp->nused; i++) {
	    if (vp->cmd[i] == RT_VLIST_LINE_MOVE ||
		    vp->cmd[i] == RT_VLIST_LINE_DRAW)
		(*count_out)++;
	}
    }

    return BRLCAD_OK;
}


static int
brlobol_vlist_export_line_set(const struct bu_list *vhead,
			      point_t **points_out,
			      int **commands_out,
			      size_t *count_out)
{
    rt_vlist *vp = NULL;
    size_t count = 0;
    size_t used = 0;

    if (!vhead || !points_out || !commands_out || !count_out)
	return BRLCAD_ERROR;

    *points_out = NULL;
    *commands_out = NULL;
    *count_out = 0;

    if (brlobol_vlist_line_count(vhead, &count) != BRLCAD_OK)
	return BRLCAD_ERROR;
    if (!count)
	return BRLCAD_OK;

    point_t *points = (point_t *)bu_calloc(count, sizeof(point_t),
	    "evaluated-wire typed line points");
    int *commands = (int *)bu_calloc(count, sizeof(int),
	    "evaluated-wire typed line commands");

    BU_LIST_EACH(vhead, vp, rt_vlist) {
	for (size_t i = 0; i < vp->nused; i++) {
	    switch (vp->cmd[i]) {
		case RT_VLIST_LINE_MOVE:
		    VMOVE(points[used], vp->pt[i]);
		    commands[used] = BG_GEOMETRY_LINE_MOVE;
		    used++;
		    break;
		case RT_VLIST_LINE_DRAW:
		    VMOVE(points[used], vp->pt[i]);
		    commands[used] = BG_GEOMETRY_LINE_DRAW;
		    used++;
		    break;
		default:
		    break;
	    }
	}
    }

    *points_out = points;
    *commands_out = commands;
    *count_out = used;
    return BRLCAD_OK;
}


int
brlobol_evaluated_wire_evaluate_path_line_set(
	struct db_i *dbip,
	const char *path,
	const struct bn_tol *tol,
	const struct bg_tess_tol *ttol,
	point_t **points_out,
	int **commands_out,
	size_t *count_out)
{
    int ret = BRLCAD_OK;
    struct bu_list vhead;
    struct bu_list vlfree;

    if (points_out)
	*points_out = NULL;
    if (commands_out)
	*commands_out = NULL;
    if (count_out)
	*count_out = 0;
    if (!dbip || !path || !path[0] || !tol || !ttol ||
	    !points_out || !commands_out || !count_out)
	return BRLCAD_ERROR;

    BU_LIST_INIT(&vhead);
    BU_LIST_INIT(&vlfree);

    struct rt_eval_wireframe_opts opts = RT_EVAL_WIREFRAME_OPTS_INIT;
    ret = rt_eval_wireframe(&vhead, &vlfree, dbip, path, tol, ttol, &opts);
    if (ret == BRLCAD_OK)
	ret = brlobol_vlist_export_line_set(&vhead, points_out, commands_out,
		count_out);

    RT_FREE_VLIST(&vlfree, &vhead);
    bg_vlist_cleanup(&vlfree);
    return ret;
}


void
brlobol_evaluated_wire_line_set_free(point_t *points, int *commands)
{
    if (points)
	bu_free(points, "evaluated-wire typed line points");
    if (commands)
	bu_free(commands, "evaluated-wire typed line commands");
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
