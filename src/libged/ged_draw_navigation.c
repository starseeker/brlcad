/*              G E D _ D R A W _ N A V I G A T I O N . C
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
/** @file ged_draw_navigation.c
 *
 * Draw-tree shape and group navigation helpers.
 */

#include "common.h"

#include <limits.h>

#include "bu/malloc.h"
#include "ged/draw.h"
#include "./ged_private.h"
#include "./ged_draw_private.h"


static int
_any_shape_cb(ged_draw_shape_ref UNUSED(ref), void *ud)
{
    *(int *)ud = 1;
    return 0;
}


int
ged_draw_has_shapes(struct ged *gedp)
{
    if (!gedp)
	return 0;

    /* A synchronized Obol scene already owns an authoritative top-level
     * database-source index.  Walking every compact occurrence to answer a
     * yes/no question is both semantically unnecessary (all occurrences of
     * one retained source resolve to the same shape ref) and disastrous for
     * large scenes: the 150k GUI stress spent measurable wall time allocating
     * transient path/context records just to rediscover that one draw root
     * exists.  Keep the traversal as the representation-independent fallback
     * for mixed/legacy scene backends. */
    if (ged_draw_obol_scene_controller_full_synced(gedp)) {
	size_t source_count = 0;
	if (ged_draw_obol_database_source_count(gedp, 0, &source_count))
	    return source_count != 0;
    }
    int found = 0;
    ged_draw_source_root_foreach_shape_ref(gedp, 0, _any_shape_cb, &found);
    return found;
}


struct _first_shape_data {
    ged_draw_shape_ref result;
};


static int
_first_shape_cb(ged_draw_shape_ref ref, void *ud)
{
    struct _first_shape_data *d = (struct _first_shape_data *)ud;
    d->result = ref;
    return 0;
}


ged_draw_shape_ref
ged_draw_first_shape_ref(struct ged *gedp)
{
    if (!gedp)
	return GED_DRAW_SHAPE_REF_NULL;

    struct _first_shape_data d;
    d.result = GED_DRAW_SHAPE_REF_NULL;
    ged_draw_source_root_foreach_shape_ref(gedp, 1, _first_shape_cb, &d);
    return d.result;
}


struct _shape_ref_snapshot {
    ged_draw_shape_ref *refs;
    size_t count;
    size_t capacity;
};


static void
_shape_ref_snapshot_append(struct _shape_ref_snapshot *snap,
			   ged_draw_shape_ref ref)
{
    if (!snap || ged_draw_shape_ref_is_null(ref))
	return;
    if (snap->count + 1 > snap->capacity) {
	size_t new_capacity = snap->capacity ? snap->capacity * 2 : 16;
	snap->refs = (ged_draw_shape_ref *)bu_realloc(snap->refs,
		new_capacity * sizeof(ged_draw_shape_ref),
		"draw navigation shape-ref snapshot");
	snap->capacity = new_capacity;
    }
    snap->refs[snap->count++] = ref;
}


static int
_snap_shape_cb(ged_draw_shape_ref ref, void *ud)
{
    _shape_ref_snapshot_append((struct _shape_ref_snapshot *)ud, ref);
    return 1;
}


static void
_sg_build_shape_snapshot(struct ged *gedp, struct _shape_ref_snapshot *out)
{
    if (!out)
	return;
    out->refs = NULL;
    out->count = 0;
    out->capacity = 0;
    ged_draw_source_root_foreach_shape_ref(gedp, 1, _snap_shape_cb,
	    (void *)out);
}


static void
_shape_ref_snapshot_free(struct _shape_ref_snapshot *snap)
{
    if (!snap)
	return;
    if (snap->refs)
	bu_free(snap->refs, "draw navigation shape-ref snapshot");
    snap->refs = NULL;
    snap->count = 0;
    snap->capacity = 0;
}


int
ged_draw_shape_count(struct ged *gedp)
{
    if (!gedp)
	return 0;

    /* Compact source records are leaf-level semantic detail beneath one
     * top-level shape ref.  The Obol source index therefore supplies the
     * exact result without materializing and deduplicating every leaf path.
     * Saturate the historical int return type defensively. */
    if (ged_draw_obol_scene_controller_full_synced(gedp)) {
	size_t source_count = 0;
	if (ged_draw_obol_database_source_count(gedp, 0, &source_count))
	    return source_count > (size_t)INT_MAX ? INT_MAX :
		(int)source_count;
    }
    struct _shape_ref_snapshot snap;
    _sg_build_shape_snapshot(gedp, &snap);
    int n = (int)snap.count;
    _shape_ref_snapshot_free(&snap);
    return n;
}


ged_draw_shape_ref
ged_draw_shape_ref_at(struct ged *gedp, int idx)
{
    if (!gedp)
	return GED_DRAW_SHAPE_REF_NULL;
    struct _shape_ref_snapshot snap;
    _sg_build_shape_snapshot(gedp, &snap);
    ged_draw_shape_ref ref = GED_DRAW_SHAPE_REF_NULL;
    int n = (int)snap.count;
    if (n > 0) {
	idx = ((idx % n) + n) % n;
	ref = snap.refs[idx];
    }
    _shape_ref_snapshot_free(&snap);
    return ref;
}


int
ged_draw_shape_ref_index(struct ged *gedp, ged_draw_shape_ref ref)
{
    if (!gedp || ged_draw_shape_ref_is_null(ref))
	return -1;
    if (ref.scene_revision && ref.scene_revision != ged_draw_scene_revision(gedp))
	return -1;

    struct _shape_ref_snapshot snap;
    _sg_build_shape_snapshot(gedp, &snap);
    int found = -1;
    for (int i = 0; i < (int)snap.count; i++) {
	if (snap.refs[i].token == ref.token) {
	    found = i;
	    break;
	}
    }
    _shape_ref_snapshot_free(&snap);
    return found;
}


ged_draw_shape_ref
ged_draw_advance_shape_ref(struct ged *gedp, ged_draw_shape_ref ref, int delta)
{
    if (!gedp)
	return GED_DRAW_SHAPE_REF_NULL;

    struct _shape_ref_snapshot snap;
    _sg_build_shape_snapshot(gedp, &snap);

    ged_draw_shape_ref result = GED_DRAW_SHAPE_REF_NULL;
    int n = (int)snap.count;
    if (n > 0) {
	int idx = 0;
	if (!ged_draw_shape_ref_is_null(ref) &&
		(!ref.scene_revision ||
		 ref.scene_revision == ged_draw_scene_revision(gedp))) {
	    for (int i = 0; i < n; i++) {
		if (snap.refs[i].token == ref.token) {
		    idx = i;
		    break;
		}
	    }
	}
	int new_idx = (((idx + delta) % n) + n) % n;
	result = snap.refs[new_idx];
    }

    _shape_ref_snapshot_free(&snap);
    return result;
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
