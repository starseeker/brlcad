/*                  G E D _ D R A W _ B R E P . C P P
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this file; see the file named COPYING for more information.
 */
/** @file ged_draw_brep.cpp
 *
 * Private BREP mesh-LoD detail setup for libged drawing paths.
 */

#include "common.h"

#include <string.h>

#include "brep/cdt.h"
#include "bn/tol.h"
#include "bu/hash.h"
#include "bu/malloc.h"
#include "raytrace.h"
#include "rt/geom.h"
#include "rt/primitives/brep.h"
#include "rt/view.h"
#include "../librt/librt_private.h"
#include "./ged_private.h"
#include "./ged_draw_private.h"


struct ged_brep_lod_detail_clbk_data {
    struct db_i *dbip;
    struct directory *dp;
    struct rt_db_internal *intern;
    struct bg_tess_tol ttol_storage;
    struct bn_tol tol_storage;
    const struct bg_tess_tol *ttol;
    const struct bn_tol *tol;
    int *faces;
    int face_cnt;
    vect_t *normals;
    point_t *pnts;
    int pnt_cnt;
};


static int
_ged_draw_mesh_lod_bounds_from_points(point_t bmin,
				      point_t bmax,
				      const point_t *points,
				      int point_count)
{
    if (!points || point_count <= 0)
	return 0;

    VMOVE(bmin, points[0]);
    VMOVE(bmax, points[0]);
    for (int i = 1; i < point_count; i++) {
	VMIN(bmin, points[i]);
	VMAX(bmax, points[i]);
    }

    return 1;
}


static int
_ged_draw_brep_lod_bounds_prepare(point_t bmin,
				  point_t bmax,
				  struct db_i *dbip,
				  struct directory *dp,
				  const struct bg_tess_tol *ttol,
				  const struct bn_tol *tol)
{
    if (!dbip || !dp)
	return 0;

    struct rt_db_internal dbintern;
    RT_DB_INTERNAL_INIT(&dbintern);
    if (rt_db_get_internal(&dbintern, dp, dbip, NULL) < 0)
	return 0;
    if (dbintern.idb_minor_type != DB5_MINORTYPE_BRLCAD_BREP) {
	rt_db_free_internal(&dbintern);
	return 0;
    }

    struct rt_brep_internal *bi = (struct rt_brep_internal *)dbintern.idb_ptr;
    RT_BREP_CK_MAGIC(bi);

    int *faces = NULL;
    int face_cnt = 0;
    vect_t *normals = NULL;
    point_t *points = NULL;
    int point_count = 0;
    int ret = brep_cdt_fast(&faces, &face_cnt, &normals, &points, &point_count,
	    bi->brep, -1, ttol, tol);
    int ok = 0;
    if (ret == BRLCAD_OK)
	ok = _ged_draw_mesh_lod_bounds_from_points(bmin, bmax,
		(const point_t *)points, point_count);

    bu_free(faces, "faces");
    bu_free(normals, "normals");
    bu_free(points, "pnts");
    rt_db_free_internal(&dbintern);

    return ok;
}


extern "C" int
ged_draw_brep_mesh_lod_cache_prepare(struct BObolMeshLod **lod,
				     point_t bmin,
				     point_t bmax,
				     int *bounds_valid,
				     struct db_i *dbip,
				     struct directory *dp,
				     const struct bg_tess_tol *ttol,
				     const struct bn_tol *tol)
{
    if (lod)
	*lod = NULL;
    if (bounds_valid)
	*bounds_valid = 0;
    if (!lod || !bounds_valid || !dbip || !dp)
	return 0;

    struct BObolMeshLod *mesh_lod = NULL;
    struct BObolMeshLodCacheStatus status = BOBOL_MESH_LOD_CACHE_STATUS_INIT;

    if (bobol_mesh_lod_cache_status(dbip, dp->d_namep, &status) != BRLCAD_OK)
	return 0;
    if (status.has_cache_key && status.has_cached_payload &&
	    !status.stale_cache_entry)
	mesh_lod = bobol_mesh_lod_get(dbip, dp->d_namep);

    if (!mesh_lod) {
	struct bu_external ext = BU_EXTERNAL_INIT_ZERO;
	if (db_get_external(&ext, dp, dbip))
	    return 0;
	unsigned long long key = bu_data_hash((void *)ext.ext_buf,
		ext.ext_nbytes);
	bu_free_external(&ext);
	if (!key)
	    return 0;

	struct rt_db_internal dbintern;
	RT_DB_INTERNAL_INIT(&dbintern);
	struct rt_db_internal *ip = &dbintern;
	int ret = rt_db_get_internal(ip, dp, dbip, NULL);
	if (ret < 0)
	    return 0;
	if (ip->idb_minor_type != DB5_MINORTYPE_BRLCAD_BREP) {
	    rt_db_free_internal(&dbintern);
	    return 0;
	}
	struct rt_brep_internal *bi = (struct rt_brep_internal *)ip->idb_ptr;
	RT_BREP_CK_MAGIC(bi);

	int *faces = NULL;
	int face_cnt = 0;
	vect_t *normals = NULL;
	point_t *pnts = NULL;
	int pnt_cnt = 0;

	ret = brep_cdt_fast(&faces, &face_cnt, &normals, &pnts, &pnt_cnt,
		bi->brep, -1, ttol, tol);
	if (ret != BRLCAD_OK) {
	    bu_free(faces, "faces");
	    bu_free(normals, "normals");
	    bu_free(pnts, "pnts");
	    rt_db_free_internal(&dbintern);
	    return 0;
	}

	if (_ged_draw_mesh_lod_bounds_from_points(bmin, bmax,
		(const point_t *)pnts, pnt_cnt))
	    *bounds_valid = 1;

	ret = bobol_mesh_lod_cache_store_mesh(dbip, dp->d_namep, (const point_t *)pnts,
		(size_t)pnt_cnt, normals, faces, (size_t)face_cnt, key, 1.0,
		&status);

	rt_db_free_internal(&dbintern);
	bu_free(faces, "faces");
	bu_free(normals, "normals");
	bu_free(pnts, "pnts");
	if (ret != BRLCAD_OK)
	    return 0;

	mesh_lod = bobol_mesh_lod_get(dbip, dp->d_namep);
    }
    if (!mesh_lod)
	return 0;

    if (!*bounds_valid) {
	if (_ged_draw_brep_lod_bounds_prepare(bmin, bmax, dbip, dp, ttol, tol)) {
	    *bounds_valid = 1;
	} else {
	    struct BObolMeshLodInfo info = BOBOL_MESH_LOD_INFO_INIT;
	    if (bobol_mesh_lod_info_get(mesh_lod, &info)) {
		VMOVE(bmin, info.bmin);
		VMOVE(bmax, info.bmax);
		*bounds_valid = 1;
	    }
	}
    }

    *lod = mesh_lod;
    return 1;
}


static int
_ged_draw_brep_mesh_info_clbk(struct BObolMeshLodDetail *detail, void *cb_data)
{
    if (!detail || !cb_data)
	return -1;

    struct ged_brep_lod_detail_clbk_data *cd = (struct ged_brep_lod_detail_clbk_data *)cb_data;
    BU_GET(cd->intern, struct rt_db_internal);
    RT_DB_INTERNAL_INIT(cd->intern);
    struct rt_db_internal *ip = cd->intern;
    int ret = rt_db_get_internal(ip, cd->dp, cd->dbip, NULL);
    if (ret < 0) {
	BU_PUT(cd->intern, struct rt_db_internal);
	cd->intern = NULL;
	return -1;
    }
    if (ip->idb_minor_type != DB5_MINORTYPE_BRLCAD_BREP || !ip->idb_ptr) {
	rt_db_free_internal(cd->intern);
	BU_PUT(cd->intern, struct rt_db_internal);
	cd->intern = NULL;
	return -1;
    }

    struct rt_brep_internal *bi = (struct rt_brep_internal *)ip->idb_ptr;
    RT_BREP_CK_MAGIC(bi);

    ret = brep_cdt_fast(&cd->faces, &cd->face_cnt, &cd->normals,
	    &cd->pnts, &cd->pnt_cnt, bi->brep, -1, cd->ttol, cd->tol);
    rt_db_free_internal(cd->intern);
    BU_PUT(cd->intern, struct rt_db_internal);
    cd->intern = NULL;
    if (ret != BRLCAD_OK) {
	bu_free(cd->faces, "faces");
	bu_free(cd->normals, "normals");
	bu_free(cd->pnts, "pnts");
	cd->faces = NULL;
	cd->face_cnt = 0;
	cd->normals = NULL;
	cd->pnts = NULL;
	cd->pnt_cnt = 0;
	return -1;
    }

    detail->faces = cd->faces;
    detail->face_count = (cd->face_cnt > 0) ? (size_t)cd->face_cnt : 0;
    detail->points = (const point_t *)cd->pnts;
    detail->point_count = (cd->pnt_cnt > 0) ? (size_t)cd->pnt_cnt : 0;
    detail->points_orig = (const point_t *)cd->pnts;
    detail->point_orig_count = (cd->pnt_cnt > 0) ? (size_t)cd->pnt_cnt : 0;
    detail->normals = cd->normals;
    detail->normal_count = (cd->normals && detail->face_count) ?
	detail->face_count * 3 : 0;

    return 0;
}


static int
_ged_draw_brep_mesh_info_clear_clbk(void *cb_data)
{
    struct ged_brep_lod_detail_clbk_data *cd = (struct ged_brep_lod_detail_clbk_data *)cb_data;
    if (!cd)
	return 0;

    if (cd->intern) {
	rt_db_free_internal(cd->intern);
	BU_PUT(cd->intern, struct rt_db_internal);
	cd->intern = NULL;
    }
    bu_free(cd->faces, "faces");
    bu_free(cd->normals, "normals");
    bu_free(cd->pnts, "pnts");
    cd->faces = NULL;
    cd->face_cnt = 0;
    cd->normals = NULL;
    cd->pnts = NULL;
    cd->pnt_cnt = 0;

    return 0;
}


static int
_ged_draw_brep_mesh_info_free_clbk(void *cb_data)
{
    _ged_draw_brep_mesh_info_clear_clbk(cb_data);
    struct ged_brep_lod_detail_clbk_data *cd = (struct ged_brep_lod_detail_clbk_data *)cb_data;
    if (!cd)
	return 0;

    BU_PUT(cd, struct ged_brep_lod_detail_clbk_data);
    return 0;
}


extern "C" int
ged_draw_brep_mesh_lod_detail_setup(struct BObolMeshLod *lod,
				    struct db_i *dbip,
				    struct directory *dp,
				    const struct bg_tess_tol *ttol,
				    const struct bn_tol *tol)
{
    if (!lod || !dbip || !dp)
	return 0;

    struct ged_brep_lod_detail_clbk_data *cbd;
    BU_GET(cbd, struct ged_brep_lod_detail_clbk_data);
    memset(cbd, 0, sizeof(*cbd));
    cbd->dbip = dbip;
    cbd->dp = dp;
    if (ttol) {
	cbd->ttol_storage = *ttol;
	cbd->ttol = &cbd->ttol_storage;
    }
    if (tol) {
	cbd->tol_storage = *tol;
	cbd->tol = &cbd->tol_storage;
    }
    if (!bobol_mesh_lod_detail_callbacks_set(lod,
	    &_ged_draw_brep_mesh_info_clbk,
	    &_ged_draw_brep_mesh_info_clear_clbk,
	    &_ged_draw_brep_mesh_info_free_clbk,
	    (void *)cbd)) {
	_ged_draw_brep_mesh_info_free_clbk((void *)cbd);
	return 0;
    }

    return 1;
}
