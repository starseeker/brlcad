/*          B S G _ G E D _ D R A W _ N A V I G A T I O N . C
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
/** @file bsg_ged_draw_navigation.c
 *
 * Draw-tree shape and group navigation helpers.
 */

#include "common.h"

#include "bu/ptbl.h"
#include "ged/draw.h"
#include "./ged_private.h"
#include "./bsg_ged_draw_private.h"


static int
_any_shape_cb(bsg_scene_ref UNUSED(ref), void *ud)
{
    *(int *)ud = 1;
    return 0;
}


int
ged_draw_has_shapes(struct ged *gedp)
{
    if (!gedp)
	return 0;
    int found = 0;
    ged_draw_scene_root_foreach_shape(gedp, 0, _any_shape_cb, &found);
    return found;
}


struct _first_shape_data {
    bsg_scene_ref result;
};


static int
_first_shape_cb(bsg_scene_ref ref, void *ud)
{
    struct _first_shape_data *d = (struct _first_shape_data *)ud;
    d->result = ref;
    return 0;
}


bsg_scene_ref
ged_draw_first_shape_scene_ref(struct ged *gedp)
{
    if (!gedp)
	return ged_draw_scene_ref_null();

    struct _first_shape_data d;
    d.result = ged_draw_scene_ref_null();
    ged_draw_scene_root_foreach_shape(gedp, 1, _first_shape_cb, &d);
    return d.result;
}


static int
_snap_shape_cb(bsg_scene_ref ref, void *ud)
{
    bu_ptbl_ins((struct bu_ptbl *)ud,
	    (long *)ged_draw_scene_ref_context(ref));
    return 1;
}


static void
_sg_build_shape_snapshot(struct ged *gedp, struct bu_ptbl *out)
{
    ged_draw_scene_root_foreach_shape(gedp, 1, _snap_shape_cb, (void *)out);
}


int
ged_draw_shape_count(struct ged *gedp)
{
    if (!gedp)
	return 0;
    struct bu_ptbl snap = BU_PTBL_INIT_ZERO;
    _sg_build_shape_snapshot(gedp, &snap);
    int n = (int)BU_PTBL_LEN(&snap);
    bu_ptbl_free(&snap);
    return n;
}


static bsg_scene_ref
_shape_scene_ref_at(struct ged *gedp, int idx)
{
    if (!gedp)
	return ged_draw_scene_ref_null();
    struct bu_ptbl snap = BU_PTBL_INIT_ZERO;
    _sg_build_shape_snapshot(gedp, &snap);
    bsg_scene_ref shape_ref = ged_draw_scene_ref_null();
    int n = (int)BU_PTBL_LEN(&snap);
    if (n > 0) {
	idx = ((idx % n) + n) % n;
	shape_ref = ged_draw_scene_ref_from_context(
		(void *)BU_PTBL_GET(&snap, idx));
    }
    bu_ptbl_free(&snap);
    return shape_ref;
}


ged_draw_shape_ref
ged_draw_first_shape_ref(struct ged *gedp)
{
    return ged_draw_shape_ref_from_scene_ref(gedp,
	    ged_draw_first_shape_scene_ref(gedp));
}


ged_draw_shape_ref
ged_draw_shape_ref_at(struct ged *gedp, int idx)
{
    return ged_draw_shape_ref_from_scene_ref(gedp,
	    _shape_scene_ref_at(gedp, idx));
}


int
ged_draw_shape_ref_index(struct ged *gedp, ged_draw_shape_ref ref)
{
    if (!gedp || ged_draw_shape_ref_is_null(ref))
	return -1;
    bsg_scene_ref target = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(target))
	return -1;

    struct bu_ptbl snap = BU_PTBL_INIT_ZERO;
    _sg_build_shape_snapshot(gedp, &snap);
    int found = -1;
    for (int i = 0; i < (int)BU_PTBL_LEN(&snap); i++) {
	if (ged_draw_scene_ref_equal(ged_draw_scene_ref_from_context(
		    (void *)BU_PTBL_GET(&snap, i)), target)) {
	    found = i;
	    break;
	}
    }
    bu_ptbl_free(&snap);
    return found;
}


ged_draw_shape_ref
ged_draw_advance_shape_ref(struct ged *gedp, ged_draw_shape_ref ref, int delta)
{
    if (!gedp)
	return GED_DRAW_SHAPE_REF_NULL;

    bsg_scene_ref target = ged_draw_registry_shape_scene_ref(gedp, ref);
    struct bu_ptbl snap = BU_PTBL_INIT_ZERO;
    _sg_build_shape_snapshot(gedp, &snap);

    ged_draw_shape_ref result = GED_DRAW_SHAPE_REF_NULL;
    int n = (int)BU_PTBL_LEN(&snap);
    if (n > 0) {
	int idx = 0;
	if (!ged_draw_scene_ref_is_null(target)) {
	    for (int i = 0; i < n; i++) {
		if (ged_draw_scene_ref_equal(
			ged_draw_scene_ref_from_context(
			    (void *)BU_PTBL_GET(&snap, i)),
			target)) {
		    idx = i;
		    break;
		}
	    }
	}
	int new_idx = (((idx + delta) % n) + n) % n;
	result = ged_draw_shape_ref_from_scene_ref(gedp,
		ged_draw_scene_ref_from_context(
		    (void *)BU_PTBL_GET(&snap, new_idx)));
    }

    bu_ptbl_free(&snap);
    return result;
}


ged_draw_group_ref
ged_draw_group_ref_of_shape(struct ged *gedp, ged_draw_shape_ref ref)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    struct ged_draw_shape_record_summary shape_summary;
    if (!ged_draw_scene_ref_shape_record_summary(shape_ref, &shape_summary))
	return GED_DRAW_GROUP_REF_NULL;
    return ged_draw_group_ref_from_scene_ref(gedp,
	    shape_summary.owning_group_ref);
}


/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
