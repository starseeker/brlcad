/*             B S G _ G E D _ D R A W _ B I N D I N G . C
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file libged/bsg_ged_draw_binding.c
 *
 * GED draw-source database binding and realization-state mirrors.
 */

#include "common.h"

#include "bu/malloc.h"
#include "bsg/appearance.h"
#include "bsg/database_source.h"
#include "bsg/scene_builder.h"
#include "ged/draw.h"
#include "rt/view.h"
#include "./ged_private.h"
#include "./bsg_ged_draw_private.h"


bsg_scene_ref
ged_draw_scene_ref_source_owner(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return bsg_scene_ref_null();
    if (bsg_database_source_ref_is_container(
		bsg_database_source_ref_from_scene(ref)))
	return ref;

    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);
    if (shape_data && !ged_draw_scene_ref_is_null(shape_data->source_ref))
	return shape_data->source_ref;

    bsg_scene_ref parent_ref = ged_draw_scene_ref_parent(ref);
    while (!ged_draw_scene_ref_is_null(parent_ref)) {
	if (bsg_database_source_ref_is_container(
		    bsg_database_source_ref_from_scene(parent_ref)))
	    return parent_ref;
	parent_ref = ged_draw_scene_ref_parent(parent_ref);
    }

    return ref;
}


bsg_scene_ref
ged_draw_view_context_database_source_create(void *view_ctx, const char *name)
{
    struct bsg_view *v = (struct bsg_view *)view_ctx;
    bsg_database_source_ref source =
	bsg_database_source_ref_create(v, name ? name : "_db_source");
    return bsg_database_source_ref_as_scene(source);
}


int
ged_draw_scene_ref_is_database_source(bsg_scene_ref ref)
{
    bsg_database_source_ref source =
	bsg_database_source_ref_from_scene(ref);
    if (bsg_database_source_ref_is_null(source))
	return 0;
    return bsg_scene_is_database_source(ref) ||
	bsg_database_source_ref_is_container(source);
}


static bsg_database_source_stale_reason
ged_draw_source_stale_reason_to_bsg(ged_draw_stale_reason reason)
{
    switch (reason) {
	case GED_DRAW_STALE_SOURCE_CHANGED:
	    return BSG_DATABASE_SOURCE_STALE_SOURCE_CHANGED;
	case GED_DRAW_STALE_VIEW_INPUT_CHANGED:
	    return BSG_DATABASE_SOURCE_STALE_VIEW_INPUT_CHANGED;
	case GED_DRAW_STALE_SETTINGS_CHANGED:
	    return BSG_DATABASE_SOURCE_STALE_SETTINGS_CHANGED;
	case GED_DRAW_STALE_FORCED:
	    return BSG_DATABASE_SOURCE_STALE_FORCED;
	case GED_DRAW_STALE_UPDATE_FAILED:
	    return BSG_DATABASE_SOURCE_STALE_UPDATE_FAILED;
	case GED_DRAW_STALE_NONE:
	default:
	    return BSG_DATABASE_SOURCE_STALE_NONE;
    }
}


static ged_draw_stale_reason
ged_draw_source_stale_reason_from_bsg(bsg_database_source_stale_reason reason)
{
    switch (reason) {
	case BSG_DATABASE_SOURCE_STALE_SOURCE_CHANGED:
	    return GED_DRAW_STALE_SOURCE_CHANGED;
	case BSG_DATABASE_SOURCE_STALE_VIEW_INPUT_CHANGED:
	    return GED_DRAW_STALE_VIEW_INPUT_CHANGED;
	case BSG_DATABASE_SOURCE_STALE_SETTINGS_CHANGED:
	    return GED_DRAW_STALE_SETTINGS_CHANGED;
	case BSG_DATABASE_SOURCE_STALE_FORCED:
	    return GED_DRAW_STALE_FORCED;
	case BSG_DATABASE_SOURCE_STALE_UPDATE_FAILED:
	    return GED_DRAW_STALE_UPDATE_FAILED;
	case BSG_DATABASE_SOURCE_STALE_NONE:
	default:
	    return GED_DRAW_STALE_NONE;
    }
}


static bsg_database_source_material_policy
ged_draw_source_material_policy_to_bsg(
	ged_draw_database_source_material_policy policy)
{
    switch (policy) {
	case GED_DRAW_DATABASE_SOURCE_MATERIAL_DATABASE:
	    return BSG_DATABASE_SOURCE_MATERIAL_DATABASE;
	case GED_DRAW_DATABASE_SOURCE_MATERIAL_INHERIT:
	default:
	    return BSG_DATABASE_SOURCE_MATERIAL_INHERIT;
    }
}


static ged_draw_database_source_material_policy
ged_draw_source_material_policy_from_bsg(
	bsg_database_source_material_policy policy)
{
    switch (policy) {
	case BSG_DATABASE_SOURCE_MATERIAL_DATABASE:
	    return GED_DRAW_DATABASE_SOURCE_MATERIAL_DATABASE;
	case BSG_DATABASE_SOURCE_MATERIAL_INHERIT:
	default:
	    return GED_DRAW_DATABASE_SOURCE_MATERIAL_INHERIT;
    }
}


static bsg_database_source_realization_status
ged_draw_source_realization_status_to_bsg(
	ged_draw_database_source_realization_status status)
{
    switch (status) {
	case GED_DRAW_DATABASE_SOURCE_REALIZATION_CURRENT:
	    return BSG_DATABASE_SOURCE_REALIZATION_CURRENT;
	case GED_DRAW_DATABASE_SOURCE_REALIZATION_STALE:
	default:
	    return BSG_DATABASE_SOURCE_REALIZATION_STALE;
    }
}


static ged_draw_database_source_realization_status
ged_draw_source_realization_status_from_bsg(
	bsg_database_source_realization_status status)
{
    switch (status) {
	case BSG_DATABASE_SOURCE_REALIZATION_CURRENT:
	    return GED_DRAW_DATABASE_SOURCE_REALIZATION_CURRENT;
	case BSG_DATABASE_SOURCE_REALIZATION_STALE:
	default:
	    return GED_DRAW_DATABASE_SOURCE_REALIZATION_STALE;
    }
}


static int
ged_draw_source_realization_roles_to_bsg(int roles)
{
    int out = BSG_DATABASE_SOURCE_REALIZATION_ROLE_NONE;
    if (roles & GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_CSG)
	out |= BSG_DATABASE_SOURCE_REALIZATION_ROLE_CSG;
    if (roles & GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_MESH)
	out |= BSG_DATABASE_SOURCE_REALIZATION_ROLE_MESH;
    return out;
}


static int
ged_draw_source_realization_roles_from_bsg(int roles)
{
    int out = GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_NONE;
    if (roles & BSG_DATABASE_SOURCE_REALIZATION_ROLE_CSG)
	out |= GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_CSG;
    if (roles & BSG_DATABASE_SOURCE_REALIZATION_ROLE_MESH)
	out |= GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_MESH;
    return out;
}


int
ged_draw_scene_ref_database_source_record_get(
	bsg_scene_ref ref,
	struct ged_draw_database_source_record *record)
{
    if (!record)
	return 0;

    memset(record, 0, sizeof(*record));
    bsg_database_source_ref source =
	bsg_database_source_ref_from_scene(ref);
    if (bsg_database_source_ref_is_null(source))
	return 0;

    struct bsg_database_source_record bsg_record;
    memset(&bsg_record, 0, sizeof(bsg_record));
    if (!bsg_database_source_record_get(source, &bsg_record))
	return 0;

    record->database_path = bsg_record.database_path;
    record->draw_mode = (int)bsg_record.draw_mode;
    record->material_policy =
	ged_draw_source_material_policy_from_bsg(
		bsg_record.material_policy);
    record->source_revision = bsg_record.source_revision;
    record->inputs_revision = bsg_record.inputs_revision;
    record->realized_source_revision =
	bsg_record.realized_source_revision;
    record->realized_inputs_revision =
	bsg_record.realized_inputs_revision;
    record->stale_reason =
	ged_draw_source_stale_reason_from_bsg(bsg_record.stale_reason);
    record->realization_identity = bsg_record.realization_identity;
    record->realization_status =
	ged_draw_source_realization_status_from_bsg(
		bsg_record.realization_status);
    record->realization_role_flags =
	ged_draw_source_realization_roles_from_bsg(
		bsg_record.realization_role_flags);
    record->realization_view_dependent =
	bsg_record.realization_view_dependent;
    record->realization_view_scale =
	(fastf_t)bsg_record.realization_view_scale;
    record->realization_bot_threshold =
	bsg_record.realization_bot_threshold;
    record->realization_curve_scale =
	(fastf_t)bsg_record.realization_curve_scale;
    record->realization_point_scale =
	(fastf_t)bsg_record.realization_point_scale;
    return 1;
}


int
ged_draw_scene_ref_database_source_record_apply(
	bsg_scene_ref ref,
	const struct ged_draw_database_source_record *record)
{
    if (!record)
	return 0;

    bsg_database_source_ref source =
	bsg_database_source_ref_from_scene(ref);
    if (bsg_database_source_ref_is_null(source))
	return 0;

    struct bsg_database_source_record bsg_record;
    memset(&bsg_record, 0, sizeof(bsg_record));
    bsg_record.database_path = record->database_path;
    bsg_record.draw_mode = (bsg_draw_mode)record->draw_mode;
    bsg_record.material_policy =
	ged_draw_source_material_policy_to_bsg(record->material_policy);
    bsg_record.source_revision = record->source_revision;
    bsg_record.inputs_revision = record->inputs_revision;
    bsg_record.realized_source_revision =
	record->realized_source_revision;
    bsg_record.realized_inputs_revision =
	record->realized_inputs_revision;
    bsg_record.stale_reason =
	ged_draw_source_stale_reason_to_bsg(record->stale_reason);
    bsg_record.realization_identity = record->realization_identity;
    bsg_record.realization_status =
	ged_draw_source_realization_status_to_bsg(
		record->realization_status);
    bsg_record.realization_role_flags =
	ged_draw_source_realization_roles_to_bsg(
		record->realization_role_flags);
    bsg_record.realization_view_dependent =
	record->realization_view_dependent;
    bsg_record.realization_view_scale =
	(double)record->realization_view_scale;
    bsg_record.realization_bot_threshold =
	record->realization_bot_threshold;
    bsg_record.realization_curve_scale =
	(double)record->realization_curve_scale;
    bsg_record.realization_point_scale =
	(double)record->realization_point_scale;

    return bsg_database_source_record_apply(source, &bsg_record);
}


struct ged_draw_source_revision_record {
    uint64_t source_revision;
    uint64_t inputs_revision;
    uint64_t realized_source_revision;
    uint64_t realized_inputs_revision;
    ged_draw_stale_reason stale_reason;
};


static void
ged_draw_scene_ref_database_source_sync(bsg_scene_ref ref,
					const struct ged_draw_source_state *source_data,
					const ged_draw_shape_state *shape_data);


static void
ged_draw_source_revision_record_from_state(
	struct ged_draw_source_revision_record *record,
	const struct ged_draw_source_state *source_data,
	const ged_draw_shape_state *shape_data)
{
    if (!record)
	return;

    memset(record, 0, sizeof(*record));
    if (shape_data) {
	record->source_revision = shape_data->source_revision;
	record->inputs_revision = shape_data->inputs_revision;
	record->realized_source_revision =
	    shape_data->realized_source_revision;
	record->realized_inputs_revision =
	    shape_data->realized_inputs_revision;
	record->stale_reason = shape_data->stale_reason;
	return;
    }

    if (source_data) {
	record->source_revision = source_data->source_revision;
	record->inputs_revision = source_data->inputs_revision;
	record->realized_source_revision =
	    source_data->realized_source_revision;
	record->realized_inputs_revision =
	    source_data->realized_inputs_revision;
	record->stale_reason = source_data->stale_reason;
    }
}


static void
ged_draw_scene_ref_database_source_sync(bsg_scene_ref ref,
					const struct ged_draw_source_state *source_data,
					const ged_draw_shape_state *shape_data)
{
    if (ged_draw_scene_ref_is_null(ref))
	return;

    const struct db_full_path *fp = NULL;
    char *path = NULL;
    bsg_scene_ref source_ref = ged_draw_scene_ref_source_owner(ref);
    struct ged_draw_database_source_record record = {0};

    (void)ged_draw_scene_ref_database_source_record_get(source_ref, &record);

    if (shape_data)
	fp = &shape_data->s_fullpath;
    else if (source_data)
	fp = source_data->fp;

    if (fp)
	path = db_path_to_string(fp);

    record.database_path = path ? path : "";
    struct ged_draw_scene_draw_state_summary draw_summary;
    if (ged_draw_scene_ref_draw_state_summary(ref, &draw_summary))
	record.draw_mode = draw_summary.draw_mode;
    record.material_policy = fp ?
	GED_DRAW_DATABASE_SOURCE_MATERIAL_DATABASE :
	GED_DRAW_DATABASE_SOURCE_MATERIAL_INHERIT;

    if (shape_data || source_data) {
	struct ged_draw_source_revision_record revisions;
	ged_draw_source_revision_record_from_state(&revisions, source_data,
		shape_data);
	record.source_revision = revisions.source_revision;
	record.inputs_revision = revisions.inputs_revision;
	record.realized_source_revision = revisions.realized_source_revision;
	record.realized_inputs_revision = revisions.realized_inputs_revision;
	record.stale_reason = revisions.stale_reason;
	if (shape_data)
	    record.realization_identity = shape_data->path_hash;
    }

    (void)ged_draw_scene_ref_database_source_record_apply(source_ref,
	    &record);
    if (!ged_draw_scene_ref_equal(source_ref, ref))
	(void)ged_draw_scene_ref_database_source_record_apply(ref, &record);

    if (path)
	bu_free(path, "db_path_to_string");
}


void
ged_draw_source_data_free(void *data)
{
    struct ged_draw_source_state *d = (struct ged_draw_source_state *)data;
    if (!d)
	return;
    if (d->fp) {
	db_free_full_path(d->fp);
	bu_free(d->fp, "ged draw source full path");
	d->fp = NULL;
    }
    if (d->rt_mesh_lod) {
	rt_mesh_lod_destroy(d->rt_mesh_lod);
	d->rt_mesh_lod = NULL;
    }
    BU_PUT(d, struct ged_draw_source_state);
}


static void
ged_draw_scene_ref_source_data_set(bsg_scene_ref ref,
				   struct ged_draw_source_state *data)
{
    if (ged_draw_scene_ref_is_null(ref) || !data)
	return;
    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);
    if (!shape_data) {
	ged_draw_source_data_free(data);
	return;
    }
    if (shape_data->source_data && shape_data->source_data != data)
	ged_draw_source_data_free(shape_data->source_data);
    shape_data->source_data = data;
    ged_draw_scene_ref_database_source_sync(ref, data, shape_data);
}


int
ged_draw_scene_ref_attach_source_state_record(
	bsg_scene_ref ref,
	const struct ged_draw_source_state_record *record)
{
    if (ged_draw_scene_ref_is_null(ref) || !record || !record->dbip ||
	    !record->fullpath || record->fullpath->fp_len <= 0)
	return 0;

    struct ged_draw_source_state *data;
    BU_GET(data, struct ged_draw_source_state);
    memset(data, 0, sizeof(*data));

    BU_GET(data->fp, struct db_full_path);
    db_full_path_init(data->fp);
    db_dup_full_path(data->fp, record->fullpath);

    data->dbip = record->dbip;
    data->tol = record->tol;
    data->ttol = record->ttol;
    if (!data->tol || !data->ttol) {
	struct rt_wdb *wdbp = wdb_dbopen(data->dbip, RT_WDB_TYPE_DB_DEFAULT);
	if (wdbp) {
	    if (!data->tol)
		data->tol = &wdbp->wdb_tol;
	    if (!data->ttol)
		data->ttol = &wdbp->wdb_ttol;
	}
    }
    data->stale_reason = record->stale_reason;

    ged_draw_scene_ref_source_data_set(ref, data);
    return ged_draw_shape_state_get_scene_ref(ref) ? 1 : 0;
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
