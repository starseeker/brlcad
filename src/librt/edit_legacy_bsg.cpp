/*                E D I T _ L E G A C Y _ B S G . C P P
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
/** @file edit_legacy_bsg.cpp
 *
 * Transitional BSG view adapters for the librt edit API.
 */

#include "common.h"

extern "C" {
#include "vmath.h"
#include "bsg/view_state.h"
#include "rt/edit_legacy_bsg.h"
}

void
rt_edit_view_from_bsg(struct rt_edit_view *ev, const struct bsg_view *v)
{
    if (!ev)
	return;

    if (!v) {
	rt_edit_view_init(ev);
	return;
    }

    ev->gv_scale = v->gv_scale;
    ev->gv_base2local = v->gv_base2local;
    ev->gv_local2base = v->gv_local2base;
    ev->gv_coord = v->gv_coord;
    ev->gv_rotate_about = v->gv_rotate_about;
    MAT_COPY(ev->gv_rotation, v->gv_rotation);
    MAT_COPY(ev->gv_center, v->gv_center);
    MAT_COPY(ev->gv_model2view, v->gv_model2view);
    MAT_COPY(ev->gv_view2model, v->gv_view2model);
}

struct rt_edit *
rt_edit_create_bsg(struct db_full_path *dfp, struct db_i *dbip,
		   struct bn_tol *tol, const struct bsg_view *v)
{
    struct rt_edit_view ev;

    if (!v)
	return rt_edit_create(dfp, dbip, tol, NULL);

    rt_edit_view_from_bsg(&ev, v);
    return rt_edit_create(dfp, dbip, tol, &ev);
}

int
rt_edit_reinit_bsg(struct rt_edit *s, struct db_full_path *dfp,
		   struct db_i *dbip, struct bn_tol *tol,
		   const struct bsg_view *v)
{
    struct rt_edit_view ev;

    if (!v)
	return rt_edit_reinit(s, dfp, dbip, tol, NULL);

    rt_edit_view_from_bsg(&ev, v);
    return rt_edit_reinit(s, dfp, dbip, tol, &ev);
}

int
rt_edit_knob_cmd_process_bsg(
	struct rt_edit *s,
	vect_t *rvec, int *do_rot, vect_t *tvec, int *do_tran, int *do_sca,
	const struct bsg_view *v, const char *cmd, fastf_t f,
	char origin, int incr_flag, void *u_data)
{
    struct rt_edit_view ev;

    if (!v)
	return rt_edit_knob_cmd_process(s, rvec, do_rot, tvec, do_tran,
	    do_sca, NULL, cmd, f, origin, incr_flag, u_data);

    rt_edit_view_from_bsg(&ev, v);
    return rt_edit_knob_cmd_process(s, rvec, do_rot, tvec, do_tran,
	do_sca, &ev, cmd, f, origin, incr_flag, u_data);
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
