/*                           L O D . C
 * BRL-CAD
 *
 * Copyright (c) 2013-2026 United States Government as represented by
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

#include "common.h"

#include "vmath.h"
#include "bu/app.h"
#include "bu/time.h"
#include "raytrace.h"
#include "rt/view.h"

int
main(int argc, char *argv[])
{
    int64_t start, elapsed;
    fastf_t seconds;
    struct db_i *dbip;
    struct directory *dp;

    bu_setprogname(argv[0]);

    if (argc != 3) {
	bu_exit(1, "Usage: %s file.g [object]", argv[0]);
    }

    start = bu_gettime();

    dbip = db_open(argv[1], DB_OPEN_READWRITE);
    if (dbip == DBI_NULL) {
	bu_exit(1, "ERROR: Unable to read from %s\n", argv[1]);
    }

    if (db_dirbuild(dbip) < 0) {
	bu_exit(1, "ERROR: Unable to read from %s\n", argv[1]);
    }

    db_update_nref(dbip);

    dp = db_lookup(dbip, argv[2], LOOKUP_QUIET);
    if (dp == RT_DIR_NULL) {
	bu_exit(1, "ERROR: Unable to look up object %s\n", argv[2]);
    }

    // Unpack bot
    struct rt_db_internal intern;
    mat_t s_mat;
    RT_DB_INTERNAL_INIT(&intern);
    MAT_IDN(s_mat);
    if (rt_db_get_internal(&intern, dp, dbip, s_mat) < 0)
	bu_exit(1, "ERROR: %s internal get failed\n", argv[2]);
    if (intern.idb_minor_type != DB5_MINORTYPE_BRLCAD_BOT)
	bu_exit(1, "ERROR: %s is not a BoT primitive\n", argv[2]);

    struct rt_bot_internal *bot = (struct rt_bot_internal *)intern.idb_ptr;
    RT_BOT_CK_MAGIC(bot);

    if (!bot->num_faces)
	bu_exit(1, "ERROR: %s - no faces found\n", argv[2]);

    struct rt_mesh_lod_cache_status cache_status = RT_MESH_LOD_CACHE_STATUS_INIT;
    if (db_mesh_lod_status(dbip, argv[2], &cache_status) != BRLCAD_OK ||
	    !cache_status.directory_found || !cache_status.is_bot ||
	    cache_status.has_cache_key || cache_status.has_cached_payload ||
	    cache_status.stale_cache_entry)
	bu_exit(1, "ERROR: %s - unexpected initial lod cache status\n", argv[2]);

    if (db_mesh_lod_refresh(dbip, argv[2], &cache_status) != BRLCAD_OK)
	bu_exit(1, "ERROR: %s - lod cache update failed\n", argv[2]);
    if (!cache_status.generated_cache_entry ||
	    !cache_status.has_cache_key ||
	    !cache_status.has_cached_payload ||
	    cache_status.stale_cache_entry)
	bu_exit(1, "ERROR: %s - lod cache refresh status failed\n", argv[2]);

    struct rt_mesh_lod *mlod = db_mesh_lod_get(dbip, argv[2]);
    if (!mlod)
	bu_exit(1, "ERROR: %s - lod cache get failed\n", argv[2]);

    elapsed = bu_gettime() - start;
    seconds = elapsed / 1000000.0;
    bu_log("Initialization time: %f seconds\n", seconds);

    start = bu_gettime();

    int coarse_fcnt = -1;
    int prev_fcnt = -1;
    int saw_face_increase = 0;
    int saw_snapped_points = 0;

    for (int i = 0; i < 16; i++) {
	struct rt_mesh_lod_data data;
	struct rt_mesh_lod_info info;
	int have_data;
	size_t face_count;
	int level = rt_mesh_lod_load_level(mlod, i, 0);
	if (level != i)
	    bu_exit(1, "ERROR: requested LoD level %d, got %d\n", i, level);
	if (rt_mesh_lod_current_level(mlod) != i)
	    bu_exit(1, "ERROR: current LoD level query did not report %d\n", i);

	have_data = rt_mesh_lod_data_get(mlod, &data);
	if (rt_mesh_lod_has_active_data(mlod) != have_data)
	    bu_exit(1, "ERROR: LoD level %d active data status mismatch\n", i);
	if (rt_mesh_lod_info_get(mlod, &info) != have_data)
	    bu_exit(1, "ERROR: LoD level %d info/data availability mismatch\n", i);
	if (have_data && (info.active_level != i ||
		info.face_count != data.face_count ||
		info.point_count != data.point_count ||
		info.point_orig_count != data.point_orig_count))
	    bu_exit(1, "ERROR: LoD level %d info/data mismatch\n", i);
	face_count = have_data ? data.face_count : 0;
	if (face_count > bot->num_faces)
	    bu_exit(1, "ERROR: LoD level %d reports invalid face count %zu (full mesh has %zu)\n", i, face_count, bot->num_faces);

	if (i == 0)
	    coarse_fcnt = (int)face_count;
	if (prev_fcnt >= 0 && face_count > (size_t)prev_fcnt)
	    saw_face_increase = 1;
	if (prev_fcnt >= 0 && face_count < (size_t)prev_fcnt)
	    bu_exit(1, "ERROR: LoD level %d has fewer faces than the previous finer step baseline (%zu < %d)\n", i, face_count, prev_fcnt);
	prev_fcnt = (int)face_count;

	if (have_data && !saw_snapped_points && data.points && data.points_orig) {
	    size_t pmax = (data.point_count < data.point_orig_count) ? data.point_count : data.point_orig_count;
	    for (size_t p = 0; p < pmax; p++) {
		if (!VNEAR_EQUAL(data.points[p], data.points_orig[p], SQRT_SMALL_FASTF)) {
		    saw_snapped_points = 1;
		    break;
		}
	    }
	}
    }

    {
	int active_level;
	struct rt_mesh_lod_data shrink_data;
	struct rt_mesh_lod_info shrink_info = RT_MESH_LOD_INFO_INIT;
	if (rt_mesh_lod_load_level(mlod, 0, 0) != 0 ||
		!rt_mesh_lod_data_get(mlod, &shrink_data))
	    bu_exit(1, "ERROR: unable to load coarse LoD before memshrink\n");
	active_level = rt_mesh_lod_current_level(mlod);
	if (active_level < 0)
	    bu_exit(1, "ERROR: current LoD level query failed before memshrink\n");
	if (!rt_mesh_lod_has_active_data(mlod))
	    bu_exit(1, "ERROR: active data status failed before memshrink\n");
	rt_mesh_lod_memshrink(mlod);
	if (rt_mesh_lod_current_level(mlod) != active_level)
	    bu_exit(1, "ERROR: memshrink did not preserve current LoD level\n");
	if (rt_mesh_lod_has_active_data(mlod))
	    bu_exit(1, "ERROR: memshrink did not clear active data status\n");
	if (rt_mesh_lod_data_get(mlod, &shrink_data))
	    bu_exit(1, "ERROR: memshrink did not release active LoD data\n");
	if (rt_mesh_lod_info_get(mlod, &shrink_info) ||
		shrink_info.active_level != active_level)
	    bu_exit(1, "ERROR: memshrink status did not report unloaded active level\n");
	if (rt_mesh_lod_load_level(mlod, active_level, 0) != active_level ||
		!rt_mesh_lod_data_get(mlod, &shrink_data))
	    bu_exit(1, "ERROR: LoD reload after memshrink failed\n");
	if (!rt_mesh_lod_has_active_data(mlod))
	    bu_exit(1, "ERROR: active data status failed after reload\n");
    }

    if (bot->num_faces > 100) {
	if ((size_t)coarse_fcnt >= bot->num_faces)
	    bu_exit(1, "ERROR: coarse LoD did not reduce face count (%d >= %zu)\n", coarse_fcnt, bot->num_faces);
	if (!saw_face_increase)
	    bu_exit(1, "ERROR: finer LoD levels did not add geometric information\n");
	if (!saw_snapped_points)
	    bu_exit(1, "ERROR: POP LoD levels did not snap points to quantized positions\n");
    }

    elapsed = bu_gettime() - start;
    seconds = elapsed / 1000000.0;
    bu_log("lod level setting: %f sec\n", seconds);

    rt_mesh_lod_destroy(mlod);
    mlod = NULL;
    if (db_mesh_lod_invalidate(dbip, argv[2], &cache_status) != BRLCAD_OK ||
	    !cache_status.cleared_cache_entry ||
	    !cache_status.cleared_cache_key ||
	    cache_status.has_cache_key ||
	    cache_status.has_cached_payload)
	bu_exit(1, "ERROR: %s - lod cache invalidation status failed\n", argv[2]);
    if (db_mesh_lod_get(dbip, argv[2]))
	bu_exit(1, "ERROR: %s - lod cache get succeeded after invalidation\n", argv[2]);
    if (db_mesh_lod_update(dbip, argv[2]) != BRLCAD_OK ||
	    !(mlod = db_mesh_lod_get(dbip, argv[2])))
	bu_exit(1, "ERROR: %s - lod cache compatibility update failed\n", argv[2]);
    rt_mesh_lod_destroy(mlod);
    rt_db_free_internal(&intern);
    db_close(dbip);

    return 0;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
