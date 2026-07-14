/*              E V A L U A T E D _ W I R E . C P P
 * BRL-CAD
 *
 * Copyright (c) 1997-2026 United States Government as represented by
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
/** @addtogroup librt */
/** @{ */
/** @file librt/eval_wireframe.cpp
 *
 * This module implements "evaluated" wireframe drawing: it takes seed nmg
 * primitive wireframes and evaluates the CSG booleans to produce more minimal
 * wireframes - a.k.a drawing mode 3, a.k.a the "bigE" mode originally exposed
 * by the "E" command.
 *
 */
/** @} */

#include "common.h"

#include "bu/str.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <algorithm>
#include <atomic>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "bu/parallel.h"
#include "bu/debug.h"
#include "bu/getopt.h"
#include "bu/malloc.h"
#include "bu/time.h"
#include "vmath.h"
#include "nmg.h"
#include "rt/geom.h"
#include "rt/boolweave.h"
#include "raytrace.h"
#include "rt/func.h"
#include "rt/misc.h"

#include "rt/eval_wireframe.h"
#include "rt/vlist.h"
#include "bg/vlist.h"
#include "../libbg/RTree.h"

static void
bigE_get_fu_plane(plane_t n, const struct faceuse *fu)
{
    const struct face_g_plane *fg = NULL;

    NMG_CK_FACEUSE(fu);
    NMG_CK_FACE(fu->f_p);
    fg = fu->f_p->g.plane_p;
    NMG_CK_FACE_G_PLANE(fg);
    if ((fu->orientation != OT_SAME) != (fu->f_p->flip != 0)) {
	HREVERSE(n, fg->N);
    } else {
	HMOVE(n, fg->N);
    }
}

struct bigE_line_record {
    point_t point;
    int command;
};

struct bigE_line_set {
    std::vector<bigE_line_record> records;
};

struct bigE_csg_eval_node {
    int op;
    int left;
    int right;
    int local_index;

    bigE_csg_eval_node() :
	op(OP_NOP),
	left(-1),
	right(-1),
	local_index(-1)
    {}
};

struct bigE_face_edge_record {
    point_t start;
    vect_t dir;
    vect_t unit_dir;

    bigE_face_edge_record()
    {
	VSETALL(start, 0.0);
	VSETALL(dir, 0.0);
	VSETALL(unit_dir, 0.0);
    }
};

struct bigE_face_record {
    struct faceuse *fu;
    struct face *face;
    point_t min_pt;
    point_t max_pt;
    plane_t plane;
    std::vector<bigE_face_edge_record> edges;
    bool edges_cached;
    size_t order;

    bigE_face_record() :
	fu(NULL),
	face(NULL),
	edges_cached(false),
	order(0)
    {
	VSETALL(min_pt, 0.0);
	VSETALL(max_pt, 0.0);
	HSETALL(plane, 0.0);
    }
};

typedef RTree<int, double, 3> bigE_face_rtree;

static const size_t BIGE_FACE_RTREE_MIN_FACES = 32;

struct bigE_leaf_face_index {
    std::vector<bigE_face_record> faces;
    std::unique_ptr<bigE_face_rtree> rtree;

    bigE_leaf_face_index() :
	faces(),
	rtree()
    {}
};

struct bigE_profile_region_summary {
    size_t index;
    std::string name;
    int64_t usec;
    int64_t edge_tab_usec;
    int64_t edge_candidate_usec;
    int64_t intersection_candidate_usec;
    int64_t intersection_face_index_usec;
    int64_t intersection_pair_collect_usec;
    int64_t intersection_pair_eval_usec;
    int64_t candidate_eval_usec;
    int64_t boolweave_shoot_usec;
    int64_t boolweave_weave_usec;
    int64_t boolweave_partition_usec;
    size_t leaves;
    size_t active_soltabs;
    size_t candidates;
    size_t edge_candidates;
    size_t intersection_candidates;
    size_t intersection_face_pairs;
    size_t degenerate_candidates;
    size_t bige_segments;
    size_t boolweave_segments;
    size_t boolweave_bbox_tests;
    size_t boolweave_bbox_hits;
    size_t boolweave_source_intervals;
    size_t boolweave_source_only;
    size_t boolweave_shot_calls;
    size_t boolweave_shot_hits;
    size_t boolweave_shot_misses;
    size_t boolweave_partitions;

    bigE_profile_region_summary() :
	index(0),
	name(),
	usec(0),
	edge_tab_usec(0),
	edge_candidate_usec(0),
	intersection_candidate_usec(0),
	intersection_face_index_usec(0),
	intersection_pair_collect_usec(0),
	intersection_pair_eval_usec(0),
	candidate_eval_usec(0),
	boolweave_shoot_usec(0),
	boolweave_weave_usec(0),
	boolweave_partition_usec(0),
	leaves(0),
	active_soltabs(0),
	candidates(0),
	edge_candidates(0),
	intersection_candidates(0),
	intersection_face_pairs(0),
	degenerate_candidates(0),
	bige_segments(0),
	boolweave_segments(0),
	boolweave_bbox_tests(0),
	boolweave_bbox_hits(0),
	boolweave_source_intervals(0),
	boolweave_source_only(0),
	boolweave_shot_calls(0),
	boolweave_shot_hits(0),
	boolweave_shot_misses(0),
	boolweave_partitions(0)
    {}
};

static void bigE_vlfree_clear(struct bu_list *vlfree);
static void bigE_context_cleanup(struct bigE_data *dgcdp);

struct bigE_data {
    struct db_i *dbip;
    struct db_full_path *fp;
    const struct bn_tol *tol;
    const struct bg_tess_tol *ttol;

    struct application *ap;
    struct resource resource;
    struct bu_ptbl leaf_list;
    struct region *active_region;
    struct rt_i *rtip;
    struct bu_list vlfree_storage;
    struct bu_list *vlfree;
    struct bigE_line_set lines;
    std::vector<struct soltab *> active_region_soltabs;
    std::vector<int> active_soltab_local_index;
    std::vector<bigE_csg_eval_node> active_csg_nodes;
    std::vector<std::vector<int>> active_source_indices_by_leaf;
    std::vector<std::vector<int>> active_candidate_indices_by_leaf;
    std::unordered_map<uint64_t, std::vector<int>> active_candidate_indices_by_pair;
    std::unordered_map<std::string, unsigned char> active_csg_eval_cache;
    int active_csg_root;
    time_t start_time;
    time_t etime;
    long nvectors;
    int do_polysolids;
    int num_halfs;
    int resource_initialized;
    int leaf_list_initialized;
    int bige_stats_enabled;
    int boolweave_eval_enabled;
    int boolweave_eval_filter_enabled;
    int boolweave_diag_enabled;
    int profile_enabled;
    int profile_region_report_enabled;
    int profile_progress_enabled;
    int profile_reported;
    uint32_t eval_flags;
    size_t boolweave_diag_limit;
    size_t profile_top_limit;
    size_t boolweave_diag_logged;
    size_t boolweave_diag_candidates;
    size_t boolweave_diag_matches;
    size_t boolweave_diag_count_mismatches;
    size_t boolweave_diag_length_mismatches;
    size_t boolweave_diag_skipped;
    size_t candidate_index;
    fastf_t boolweave_diag_max_abs_len_delta;
    fastf_t boolweave_diag_total_abs_len_delta;
    int64_t boolweave_diag_bige_usec;
    int64_t boolweave_diag_bw_usec;
    int64_t boolweave_diag_bw_shot_usec;
    int64_t boolweave_diag_bw_bool_usec;
    size_t bige_stats_candidates;
    size_t bige_stats_segments;
    int64_t bige_stats_usec;
    size_t boolweave_eval_candidates;
    size_t boolweave_eval_segments;
    int64_t boolweave_eval_usec;
    int64_t profile_total_usec;
    int64_t profile_edge_tab_usec;
    int64_t profile_edge_candidate_usec;
    int64_t profile_intersection_candidate_usec;
    int64_t profile_intersection_face_index_usec;
    int64_t profile_intersection_pair_collect_usec;
    int64_t profile_intersection_pair_eval_usec;
    int64_t profile_candidate_eval_usec;
    int64_t profile_boolweave_shoot_usec;
    int64_t profile_boolweave_weave_usec;
    int64_t profile_boolweave_partition_usec;
    size_t profile_regions;
    size_t profile_total_leaves;
    size_t profile_max_leaves;
    size_t profile_total_active_soltabs;
    size_t profile_max_active_soltabs;
    size_t profile_total_candidates;
    size_t profile_edge_candidates;
    size_t profile_intersection_candidates;
    size_t profile_intersection_face_pairs;
    size_t profile_degenerate_candidates;
    size_t profile_bige_segments;
    size_t profile_boolweave_segments;
    size_t profile_boolweave_bbox_tests;
    size_t profile_boolweave_bbox_hits;
    size_t profile_boolweave_source_intervals;
    size_t profile_boolweave_source_only;
    size_t profile_boolweave_shot_calls;
    size_t profile_boolweave_shot_hits;
    size_t profile_boolweave_shot_misses;
    size_t profile_boolweave_partitions;
    int64_t profile_current_start_usec;
    bigE_profile_region_summary profile_current;
    std::vector<bigE_profile_region_summary> profile_top_regions;

    bigE_data();
    ~bigE_data();
    bigE_data(const bigE_data &) = delete;
    bigE_data &operator=(const bigE_data &) = delete;
};


bigE_data::bigE_data() :
    dbip(NULL),
    fp(NULL),
    tol(NULL),
    ttol(NULL),
    ap(NULL),
    resource(),
    leaf_list(),
    active_region(NULL),
    rtip(NULL),
    vlfree_storage(),
    vlfree(&vlfree_storage),
    lines(),
    active_region_soltabs(),
    active_soltab_local_index(),
    active_csg_nodes(),
    active_source_indices_by_leaf(),
    active_candidate_indices_by_leaf(),
    active_candidate_indices_by_pair(),
    active_csg_eval_cache(),
    active_csg_root(-1),
    start_time(0),
    etime(0),
    nvectors(0),
    do_polysolids(0),
    num_halfs(0),
    resource_initialized(0),
    leaf_list_initialized(0),
    bige_stats_enabled(0),
    boolweave_eval_enabled(0),
    boolweave_eval_filter_enabled(0),
    boolweave_diag_enabled(0),
    profile_enabled(0),
    profile_region_report_enabled(0),
    profile_progress_enabled(0),
    profile_reported(0),
    eval_flags(RT_EVAL_WIREFRAME_F_DEFAULT),
    boolweave_diag_limit(0),
    profile_top_limit(10),
    boolweave_diag_logged(0),
    boolweave_diag_candidates(0),
    boolweave_diag_matches(0),
    boolweave_diag_count_mismatches(0),
    boolweave_diag_length_mismatches(0),
    boolweave_diag_skipped(0),
    candidate_index(0),
    boolweave_diag_max_abs_len_delta(0.0),
    boolweave_diag_total_abs_len_delta(0.0),
    boolweave_diag_bige_usec(0),
    boolweave_diag_bw_usec(0),
    boolweave_diag_bw_shot_usec(0),
    boolweave_diag_bw_bool_usec(0),
    bige_stats_candidates(0),
    bige_stats_segments(0),
    bige_stats_usec(0),
    boolweave_eval_candidates(0),
    boolweave_eval_segments(0),
    boolweave_eval_usec(0),
    profile_total_usec(0),
    profile_edge_tab_usec(0),
    profile_edge_candidate_usec(0),
    profile_intersection_candidate_usec(0),
    profile_intersection_face_index_usec(0),
    profile_intersection_pair_collect_usec(0),
    profile_intersection_pair_eval_usec(0),
    profile_candidate_eval_usec(0),
    profile_boolweave_shoot_usec(0),
    profile_boolweave_weave_usec(0),
    profile_boolweave_partition_usec(0),
    profile_regions(0),
    profile_total_leaves(0),
    profile_max_leaves(0),
    profile_total_active_soltabs(0),
    profile_max_active_soltabs(0),
    profile_total_candidates(0),
    profile_edge_candidates(0),
    profile_intersection_candidates(0),
    profile_intersection_face_pairs(0),
    profile_degenerate_candidates(0),
    profile_bige_segments(0),
    profile_boolweave_segments(0),
    profile_boolweave_bbox_tests(0),
    profile_boolweave_bbox_hits(0),
    profile_boolweave_source_intervals(0),
    profile_boolweave_source_only(0),
    profile_boolweave_shot_calls(0),
    profile_boolweave_shot_hits(0),
    profile_boolweave_shot_misses(0),
    profile_boolweave_partitions(0),
    profile_current_start_usec(0),
    profile_current(),
    profile_top_regions()
{
    memset(&resource, 0, sizeof(resource));
    memset(&leaf_list, 0, sizeof(leaf_list));
    BU_LIST_INIT(&vlfree_storage);
}


bigE_data::~bigE_data()
{
    bigE_context_cleanup(this);
}


struct bigE_region_result {
    size_t index;
    int ret;
    std::string debug_name;
    struct bigE_line_set lines;

    bigE_region_result() : index(0), ret(BRLCAD_OK), debug_name(), lines() {}
};

struct bigE_segment_summary {
    size_t count;
    fastf_t length;

    bigE_segment_summary() : count(0), length(0.0) {}
};

/* #define debug 1 */

static union E_tree *build_etree(union tree *tp, struct bigE_data *dgcdp);
static void Eplot(union E_tree *eptr, struct bigE_data *dgcdp);
static void free_etree(union E_tree *eptr, struct bigE_data *dgcdp);
static void fix_halfs(struct bigE_data *dgcdp);
static void bigE_collect_region_soltabs(union tree *tp,
					std::vector<struct soltab *> *soltabs);
static void bigE_build_active_soltab_index(struct bigE_data *dgcdp);
static int bigE_compile_active_csg_tree(struct bigE_data *dgcdp,
					union tree *tp);
static void bigE_build_active_source_index(struct bigE_data *dgcdp);

/* segment types (stored in the "seg_stp" field of the (struct seg) */
#define ON_SURF	(struct soltab *)0x1
#define IN_SOL  (struct soltab *)0x2
#define ON_INT	(struct soltab *)0x3

#define NOT_SEG_OVERLAP(_a, _b)	((_a->seg_out.hit_dist <= _b->seg_in.hit_dist) || (_b->seg_out.hit_dist <= _a->seg_in.hit_dist))

/* stolen from g_half.c */
struct half_specific {
    plane_t half_eqn;              /* Plane equation, outward normal */
    vect_t half_Xbase;             /* "X" basis direction */
    vect_t half_Ybase;             /* "Y" basis direction */
};
#define HALF_NULL ((struct half_specific *)0)

/* structures for building a tree corresponding to the region to be
 * drawn uses the same "op" values as "union tree"
 */
union E_tree {
    uint32_t magic;

    struct E_node {
	/* the operator nodes */
	uint32_t magic;
	int op;
	union E_tree *left;
	union E_tree *right;
    } n;

    struct E_leaf {
	/* the leaf nodes */
	uint32_t magic;
	int op;
	struct model *m;		 /* NMG version of this leaf solid */
	struct bu_list seghead;		 /* head of list of segments for this leaf solid */
	struct bu_ptbl edge_list;	 /* list of edges from above NMG */
	struct soltab *stp;		 /* the usual soltab pointer */
	unsigned char do_not_free_model; /* A flag indicating that the
					  * NMG model pointer is a
					  * reference to the NMG model
					  * in the soltab structure.
					  */
	unsigned char edge_list_initialized;
    } l;
};


#define E_TREE_MAGIC 0x45545245
#define CK_ETREE(_p) BU_CKMAG(_p, E_TREE_MAGIC, "struct E_tree")


static void
bigE_vlfree_clear(struct bu_list *vlfree)
{
    bg_vlist *vp;

    if (!vlfree || !BU_LIST_IS_INITIALIZED(vlfree))
	return;

    while (BU_LIST_WHILE (vp, bg_vlist, vlfree)) {
	BU_LIST_DEQUEUE(&vp->l);
	bu_free(vp, "bigE vlist free node");
    }
}


static struct resource *
bigE_resource(struct bigE_data *dgcdp)
{
    return (dgcdp && dgcdp->ap) ? dgcdp->ap->a_resource : RESOURCE_NULL;
}


static struct seg *
bigE_seg_get(struct bigE_data *dgcdp)
{
    struct seg *seg = NULL;
    RT_GET_SEG(seg, bigE_resource(dgcdp));
    return seg;
}


static void
bigE_seg_free(struct bigE_data *dgcdp, struct seg *seg)
{
    if (seg)
	RT_FREE_SEG(seg, bigE_resource(dgcdp));
}


static void
bigE_seg_list_free(struct bu_list *seghead, struct bigE_data *dgcdp)
{
    struct seg *seg;

    if (!seghead)
	return;

    while (BU_LIST_WHILE (seg, seg, seghead)) {
	BU_LIST_DEQUEUE(&seg->l);
	bigE_seg_free(dgcdp, seg);
    }
}


static struct seg *
bigE_ray_seghead_create(void)
{
    struct seg *seghead;

    BU_ALLOC(seghead, struct seg);
    BU_LIST_INIT(&seghead->l);
    return seghead;
}


static void
bigE_ray_seghead_free(struct bigE_data *dgcdp, struct seg *seghead)
{
    if (!seghead)
	return;

    bigE_seg_list_free(&seghead->l, dgcdp);
    bu_free(seghead, "ray_data seghead");
}


static size_t
bigE_leaf_count(const struct bigE_data *dgcdp)
{
    return dgcdp ? BU_PTBL_LEN(&dgcdp->leaf_list) : 0;
}


static union E_tree *
bigE_leaf_at(struct bigE_data *dgcdp, size_t index)
{
    if (!dgcdp || index >= bigE_leaf_count(dgcdp))
	return NULL;
    return (union E_tree *)BU_PTBL_GET(&dgcdp->leaf_list, index);
}


static union E_tree *
bigE_tree_create(void)
{
    union E_tree *eptr = NULL;

    BU_ALLOC(eptr, union E_tree);
    memset(eptr, 0, sizeof(union E_tree));
    eptr->magic = E_TREE_MAGIC;
    return eptr;
}


static struct soltab *
bigE_soltab_create(const struct directory *dp, matp_t mat)
{
    struct soltab *stp = NULL;

    BU_ALLOC(stp, struct soltab);
    memset(stp, 0, sizeof(struct soltab));
    stp->l.magic = RT_SOLTAB_MAGIC;
    stp->l2.magic = RT_SOLTAB2_MAGIC;
    stp->st_dp = dp;
    stp->st_matp = mat;
    return stp;
}


static int
bigE_soltabs_match_instance(const struct soltab *a,
			    const struct soltab *b,
			    const struct bn_tol *tol)
{
    if (!a || !b || a->st_dp != b->st_dp)
	return 0;

    if (!a->st_matp && !b->st_matp)
	return 1;
    if (!a->st_matp)
	return bn_mat_is_equal(b->st_matp, bn_mat_identity, tol);
    if (!b->st_matp)
	return bn_mat_is_equal(a->st_matp, bn_mat_identity, tol);

    return bn_mat_is_equal(a->st_matp, b->st_matp, tol);
}


static int
bigE_should_skip_leaf(struct bigE_data *dgcdp,
		      size_t shoot_leaf,
		      int skip_leaf1,
		      int skip_leaf2)
{
    union E_tree *shoot = bigE_leaf_at(dgcdp, shoot_leaf);
    int skip_ids[2] = {skip_leaf1, skip_leaf2};

    if (!shoot)
	return 1;

    for (size_t i = 0; i < 2; i++) {
	int skip_leaf = skip_ids[i];
	if (skip_leaf < 0)
	    continue;
	if ((size_t)skip_leaf == shoot_leaf)
	    return 1;

	union E_tree *source = bigE_leaf_at(dgcdp, (size_t)skip_leaf);
	if (source && bigE_soltabs_match_instance(source->l.stp, shoot->l.stp,
		dgcdp->tol))
	    return 1;
    }

    return 0;
}


static int
bigE_soltab_is_candidate_source(struct bigE_data *dgcdp,
				const struct soltab *stp,
				int skip_leaf1,
				int skip_leaf2)
{
    int skip_ids[2] = {skip_leaf1, skip_leaf2};

    if (!dgcdp || !stp)
	return 0;

    for (size_t i = 0; i < 2; i++) {
	int skip_leaf = skip_ids[i];
	if (skip_leaf < 0)
	    continue;
	union E_tree *source = bigE_leaf_at(dgcdp, (size_t)skip_leaf);
	if (source && bigE_soltabs_match_instance(source->l.stp, stp,
		dgcdp->tol))
	    return 1;
    }

    return 0;
}


static void
bigE_line_append(struct bigE_data *dgcdp, const point_t point, int command)
{
    if (!dgcdp || !point)
	return;

    struct bigE_line_record rec;
    VMOVE(rec.point, point);
    rec.command = command;
    dgcdp->lines.records.push_back(rec);
}


static void
bigE_line_free(struct bigE_data *dgcdp)
{
    if (!dgcdp)
	return;
    dgcdp->lines.records.clear();
    dgcdp->lines.records.shrink_to_fit();
}


static int
bigE_line_set_export_vlist(struct bigE_line_set *lines,
			   struct bu_list *vhead,
			   struct bu_list *vlfree)
{
    if (!lines || !vhead || !vlfree)
	return BRLCAD_ERROR;

    const size_t count = lines->records.size();
    if (!count)
	return BRLCAD_OK;

    for (size_t i = 0; i < count; i++) {
	RT_ADD_VLIST(vlfree, vhead, lines->records[i].point,
		lines->records[i].command);
    }

    lines->records.clear();
    return BRLCAD_OK;
}


static void
bigE_boolweave_diag_reset(struct bigE_data *dgcdp)
{
    if (!dgcdp)
	return;

    dgcdp->boolweave_diag_logged = 0;
    dgcdp->boolweave_diag_candidates = 0;
    dgcdp->boolweave_diag_matches = 0;
    dgcdp->boolweave_diag_count_mismatches = 0;
    dgcdp->boolweave_diag_length_mismatches = 0;
    dgcdp->boolweave_diag_skipped = 0;
    dgcdp->boolweave_diag_max_abs_len_delta = 0.0;
    dgcdp->boolweave_diag_total_abs_len_delta = 0.0;
    dgcdp->boolweave_diag_bige_usec = 0;
    dgcdp->boolweave_diag_bw_usec = 0;
    dgcdp->boolweave_diag_bw_shot_usec = 0;
    dgcdp->boolweave_diag_bw_bool_usec = 0;
    dgcdp->bige_stats_candidates = 0;
    dgcdp->bige_stats_segments = 0;
    dgcdp->bige_stats_usec = 0;
    dgcdp->boolweave_eval_candidates = 0;
    dgcdp->boolweave_eval_segments = 0;
    dgcdp->boolweave_eval_usec = 0;
}


static void
bigE_boolweave_diag_report(struct bigE_data *dgcdp)
{
    if (!dgcdp || !dgcdp->boolweave_diag_enabled ||
	    !dgcdp->boolweave_diag_candidates)
	return;

    const double bige_sec = ((double)dgcdp->boolweave_diag_bige_usec) / 1000000.0;
    const double bw_sec = ((double)dgcdp->boolweave_diag_bw_usec) / 1000000.0;
    const double bw_shot_sec =
	((double)dgcdp->boolweave_diag_bw_shot_usec) / 1000000.0;
    const double bw_bool_sec =
	((double)dgcdp->boolweave_diag_bw_bool_usec) / 1000000.0;
    const double ratio = (dgcdp->boolweave_diag_bige_usec > 0) ?
	((double)dgcdp->boolweave_diag_bw_usec /
	 (double)dgcdp->boolweave_diag_bige_usec) : 0.0;

    bu_log("evaluated-wire boolweave summary: candidates=%zu matched=%zu "
	    "count_mismatch=%zu length_mismatch=%zu skipped=%zu "
	    "detail_logged=%zu max_abs_len_delta=%.12g total_abs_len_delta=%.12g "
	    "bige_time=%.6fs boolweave_time=%.6fs ratio=%.6g "
	    "boolweave_shot=%.6fs boolweave_bool=%.6fs\n",
	    dgcdp->boolweave_diag_candidates,
	    dgcdp->boolweave_diag_matches,
	    dgcdp->boolweave_diag_count_mismatches,
	    dgcdp->boolweave_diag_length_mismatches,
	    dgcdp->boolweave_diag_skipped,
	    dgcdp->boolweave_diag_logged,
	    dgcdp->boolweave_diag_max_abs_len_delta,
	    dgcdp->boolweave_diag_total_abs_len_delta,
	    bige_sec,
	    bw_sec,
	    ratio,
	    bw_shot_sec,
	    bw_bool_sec);

    bigE_boolweave_diag_reset(dgcdp);
}


static void
bigE_boolweave_eval_report(struct bigE_data *dgcdp)
{
    if (!dgcdp || !dgcdp->boolweave_eval_enabled ||
	    !dgcdp->boolweave_eval_candidates)
	return;

    const double eval_sec = ((double)dgcdp->boolweave_eval_usec) / 1000000.0;
    bu_log("evaluated-wire boolweave eval summary: candidates=%zu "
	    "segments=%zu eval_time=%.6fs\n",
	    dgcdp->boolweave_eval_candidates,
	    dgcdp->boolweave_eval_segments,
	    eval_sec);
}


static void
bigE_stats_report(struct bigE_data *dgcdp)
{
    if (!dgcdp || !dgcdp->bige_stats_enabled || !dgcdp->bige_stats_candidates)
	return;

    const double bige_sec = ((double)dgcdp->bige_stats_usec) / 1000000.0;
    bu_log("evaluated-wire Big-E summary: candidates=%zu segments=%zu "
	    "eval_time=%.6fs\n",
	    dgcdp->bige_stats_candidates,
	    dgcdp->bige_stats_segments,
	    bige_sec);
}


static int
bigE_env_enabled(const char *name)
{
    const char *value = getenv(name);
    return (value && value[0] && bu_strcmp(value, "0")) ? 1 : 0;
}


static size_t
bigE_env_size(const char *name, size_t default_value)
{
    const char *value = getenv(name);
    if (!value || !value[0])
	return default_value;
    return (size_t)strtoull(value, NULL, 10);
}


static const char *
bigE_method_name(uint32_t eval_flags)
{
    const uint32_t method = eval_flags & RT_EVAL_WIREFRAME_F_METHOD_MASK;
    if (method == RT_EVAL_WIREFRAME_F_BOOLWEAVE)
	return "boolweave";
    if (method == RT_EVAL_WIREFRAME_F_BIGE)
	return "bige";
    return "boolweave";
}


static void
bigE_profile_reset(struct bigE_data *dgcdp)
{
    if (!dgcdp)
	return;

    dgcdp->profile_enabled = 0;
    dgcdp->profile_region_report_enabled = 0;
    dgcdp->profile_progress_enabled = 0;
    dgcdp->profile_reported = 0;
    dgcdp->profile_top_limit = 10;
    dgcdp->profile_total_usec = 0;
    dgcdp->profile_edge_tab_usec = 0;
    dgcdp->profile_edge_candidate_usec = 0;
    dgcdp->profile_intersection_candidate_usec = 0;
    dgcdp->profile_intersection_face_index_usec = 0;
    dgcdp->profile_intersection_pair_collect_usec = 0;
    dgcdp->profile_intersection_pair_eval_usec = 0;
    dgcdp->profile_candidate_eval_usec = 0;
    dgcdp->profile_boolweave_shoot_usec = 0;
    dgcdp->profile_boolweave_weave_usec = 0;
    dgcdp->profile_boolweave_partition_usec = 0;
    dgcdp->profile_regions = 0;
    dgcdp->profile_total_leaves = 0;
    dgcdp->profile_max_leaves = 0;
    dgcdp->profile_total_active_soltabs = 0;
    dgcdp->profile_max_active_soltabs = 0;
    dgcdp->profile_total_candidates = 0;
    dgcdp->profile_edge_candidates = 0;
    dgcdp->profile_intersection_candidates = 0;
    dgcdp->profile_intersection_face_pairs = 0;
    dgcdp->profile_degenerate_candidates = 0;
    dgcdp->profile_bige_segments = 0;
    dgcdp->profile_boolweave_segments = 0;
    dgcdp->profile_boolweave_bbox_tests = 0;
    dgcdp->profile_boolweave_bbox_hits = 0;
    dgcdp->profile_boolweave_source_intervals = 0;
    dgcdp->profile_boolweave_source_only = 0;
    dgcdp->profile_boolweave_shot_calls = 0;
    dgcdp->profile_boolweave_shot_hits = 0;
    dgcdp->profile_boolweave_shot_misses = 0;
    dgcdp->profile_boolweave_partitions = 0;
    dgcdp->profile_current_start_usec = 0;
    dgcdp->profile_current = bigE_profile_region_summary();
    dgcdp->profile_top_regions.clear();
}


static void
bigE_profile_region_start(struct bigE_data *dgcdp,
			  size_t index,
			  const struct region *rp)
{
    if (!dgcdp || !dgcdp->profile_enabled)
	return;

    dgcdp->profile_current = bigE_profile_region_summary();
    dgcdp->profile_current.index = index;
    if (rp && rp->reg_name)
	dgcdp->profile_current.name = rp->reg_name;
    dgcdp->profile_current_start_usec = bu_gettime();
}


static void
bigE_profile_top_insert(struct bigE_data *dgcdp,
			const bigE_profile_region_summary &summary)
{
    if (!dgcdp || !dgcdp->profile_enabled || dgcdp->profile_top_limit == 0)
	return;

    dgcdp->profile_top_regions.push_back(summary);
    std::sort(dgcdp->profile_top_regions.begin(),
	    dgcdp->profile_top_regions.end(),
	    [](const bigE_profile_region_summary &a,
	       const bigE_profile_region_summary &b) {
		return a.usec > b.usec;
	    });
    if (dgcdp->profile_top_regions.size() > dgcdp->profile_top_limit)
	dgcdp->profile_top_regions.resize(dgcdp->profile_top_limit);
}


static void
bigE_profile_region_finish(struct bigE_data *dgcdp)
{
    if (!dgcdp || !dgcdp->profile_enabled || !dgcdp->profile_current_start_usec)
	return;

    bigE_profile_region_summary summary = dgcdp->profile_current;
    summary.usec = bu_gettime() - dgcdp->profile_current_start_usec;

    dgcdp->profile_regions++;
    dgcdp->profile_total_usec += summary.usec;
    dgcdp->profile_edge_tab_usec += summary.edge_tab_usec;
    dgcdp->profile_edge_candidate_usec += summary.edge_candidate_usec;
    dgcdp->profile_intersection_candidate_usec +=
	summary.intersection_candidate_usec;
    dgcdp->profile_intersection_face_index_usec +=
	summary.intersection_face_index_usec;
    dgcdp->profile_intersection_pair_collect_usec +=
	summary.intersection_pair_collect_usec;
    dgcdp->profile_intersection_pair_eval_usec +=
	summary.intersection_pair_eval_usec;
    dgcdp->profile_candidate_eval_usec += summary.candidate_eval_usec;
    dgcdp->profile_boolweave_shoot_usec += summary.boolweave_shoot_usec;
    dgcdp->profile_boolweave_weave_usec += summary.boolweave_weave_usec;
    dgcdp->profile_boolweave_partition_usec += summary.boolweave_partition_usec;
    dgcdp->profile_total_leaves += summary.leaves;
    if (summary.leaves > dgcdp->profile_max_leaves)
	dgcdp->profile_max_leaves = summary.leaves;
    dgcdp->profile_total_active_soltabs += summary.active_soltabs;
    if (summary.active_soltabs > dgcdp->profile_max_active_soltabs)
	dgcdp->profile_max_active_soltabs = summary.active_soltabs;
    dgcdp->profile_total_candidates += summary.candidates;
    dgcdp->profile_edge_candidates += summary.edge_candidates;
    dgcdp->profile_intersection_candidates += summary.intersection_candidates;
    dgcdp->profile_intersection_face_pairs +=
	summary.intersection_face_pairs;
    dgcdp->profile_degenerate_candidates += summary.degenerate_candidates;
    dgcdp->profile_bige_segments += summary.bige_segments;
    dgcdp->profile_boolweave_segments += summary.boolweave_segments;
    dgcdp->profile_boolweave_bbox_tests += summary.boolweave_bbox_tests;
    dgcdp->profile_boolweave_bbox_hits += summary.boolweave_bbox_hits;
    dgcdp->profile_boolweave_source_intervals +=
	summary.boolweave_source_intervals;
    dgcdp->profile_boolweave_source_only += summary.boolweave_source_only;
    dgcdp->profile_boolweave_shot_calls += summary.boolweave_shot_calls;
    dgcdp->profile_boolweave_shot_hits += summary.boolweave_shot_hits;
    dgcdp->profile_boolweave_shot_misses += summary.boolweave_shot_misses;
    dgcdp->profile_boolweave_partitions += summary.boolweave_partitions;

    if (dgcdp->profile_region_report_enabled) {
	bu_log("evaluated-wire profile region: index=%zu name=\"%s\" "
		"time=%.6fs leaves=%zu active_soltabs=%zu candidates=%zu "
		"edge=%zu intersection=%zu face_pairs=%zu "
		"degenerate=%zu "
		"segments=(bige:%zu boolweave:%zu) "
		"intersection_time=(index:%.6fs collect:%.6fs eval:%.6fs) "
		"boolweave=(bbox:%zu/%zu source:%zu source_only:%zu "
		"shot:%zu hit:%zu miss:%zu "
		"partitions:%zu shoot_time:%.6fs weave_time:%.6fs "
		"partition_time:%.6fs eval_time:%.6fs)\n",
		summary.index,
		summary.name.empty() ? "(null)" : summary.name.c_str(),
		((double)summary.usec) / 1000000.0,
		summary.leaves,
		summary.active_soltabs,
		summary.candidates,
		summary.edge_candidates,
		summary.intersection_candidates,
		summary.intersection_face_pairs,
		summary.degenerate_candidates,
		summary.bige_segments,
		summary.boolweave_segments,
		((double)summary.intersection_face_index_usec) / 1000000.0,
		((double)summary.intersection_pair_collect_usec) / 1000000.0,
		((double)summary.intersection_pair_eval_usec) / 1000000.0,
		summary.boolweave_bbox_hits,
		summary.boolweave_bbox_tests,
		summary.boolweave_source_intervals,
		summary.boolweave_source_only,
		summary.boolweave_shot_calls,
		summary.boolweave_shot_hits,
		summary.boolweave_shot_misses,
		summary.boolweave_partitions,
		((double)summary.boolweave_shoot_usec) / 1000000.0,
		((double)summary.boolweave_weave_usec) / 1000000.0,
		((double)summary.boolweave_partition_usec) / 1000000.0,
		((double)summary.candidate_eval_usec) / 1000000.0);
    }

    bigE_profile_top_insert(dgcdp, summary);
    dgcdp->profile_current = bigE_profile_region_summary();
    dgcdp->profile_current_start_usec = 0;
}


static void
bigE_profile_report(struct bigE_data *dgcdp)
{
    if (!dgcdp || !dgcdp->profile_enabled || dgcdp->profile_reported ||
	    !dgcdp->profile_regions)
	return;

    dgcdp->profile_reported = 1;
    const double total_sec = ((double)dgcdp->profile_total_usec) / 1000000.0;
    const double edge_tab_sec =
	((double)dgcdp->profile_edge_tab_usec) / 1000000.0;
    const double edge_candidate_sec =
	((double)dgcdp->profile_edge_candidate_usec) / 1000000.0;
    const double intersection_candidate_sec =
	((double)dgcdp->profile_intersection_candidate_usec) / 1000000.0;
    const double intersection_face_index_sec =
	((double)dgcdp->profile_intersection_face_index_usec) / 1000000.0;
    const double intersection_pair_collect_sec =
	((double)dgcdp->profile_intersection_pair_collect_usec) / 1000000.0;
    const double intersection_pair_eval_sec =
	((double)dgcdp->profile_intersection_pair_eval_usec) / 1000000.0;
    const double candidate_eval_sec =
	((double)dgcdp->profile_candidate_eval_usec) / 1000000.0;
    const double boolweave_shoot_sec =
	((double)dgcdp->profile_boolweave_shoot_usec) / 1000000.0;
    const double boolweave_weave_sec =
	((double)dgcdp->profile_boolweave_weave_usec) / 1000000.0;
    const double boolweave_partition_sec =
	((double)dgcdp->profile_boolweave_partition_usec) / 1000000.0;
    const double avg_leaves =
	((double)dgcdp->profile_total_leaves) / (double)dgcdp->profile_regions;
    const double avg_active =
	((double)dgcdp->profile_total_active_soltabs) /
	(double)dgcdp->profile_regions;
    const double bbox_keep =
	dgcdp->profile_boolweave_bbox_tests ?
	((double)dgcdp->profile_boolweave_bbox_hits /
	 (double)dgcdp->profile_boolweave_bbox_tests) : 0.0;

    bu_log("evaluated-wire profile summary: method=%s regions=%zu "
	    "time=%.6fs candidates=%zu edge=%zu intersection=%zu "
	    "face_pairs=%zu degenerate=%zu "
	    "segments=(bige:%zu boolweave:%zu) "
	    "leaves=(avg:%.3f max:%zu) active_soltabs=(avg:%.3f max:%zu)\n",
	    bigE_method_name(dgcdp->eval_flags),
	    dgcdp->profile_regions,
	    total_sec,
	    dgcdp->profile_total_candidates,
	    dgcdp->profile_edge_candidates,
	    dgcdp->profile_intersection_candidates,
	    dgcdp->profile_intersection_face_pairs,
	    dgcdp->profile_degenerate_candidates,
	    dgcdp->profile_bige_segments,
	    dgcdp->profile_boolweave_segments,
	    avg_leaves,
	    dgcdp->profile_max_leaves,
	    avg_active,
	    dgcdp->profile_max_active_soltabs);

    bu_log("evaluated-wire profile stages: edge_tab=%.6fs "
	    "edge_candidates=%.6fs intersection_candidates=%.6fs "
	    "intersection_face_index=%.6fs intersection_pair_collect=%.6fs "
	    "intersection_pair_eval=%.6fs candidate_eval=%.6fs\n",
	    edge_tab_sec,
	    edge_candidate_sec,
	    intersection_candidate_sec,
	    intersection_face_index_sec,
	    intersection_pair_collect_sec,
	    intersection_pair_eval_sec,
	    candidate_eval_sec);

    if (dgcdp->profile_boolweave_bbox_tests ||
	    dgcdp->profile_boolweave_shot_calls ||
	    dgcdp->profile_boolweave_partitions) {
	bu_log("evaluated-wire profile boolweave: bbox_hits=%zu "
		"bbox_tests=%zu bbox_keep=%.6f source_intervals=%zu "
		"source_only=%zu shot_calls=%zu shot_hits=%zu "
		"shot_misses=%zu partitions=%zu shoot_time=%.6fs "
		"weave_time=%.6fs partition_time=%.6fs\n",
		dgcdp->profile_boolweave_bbox_hits,
		dgcdp->profile_boolweave_bbox_tests,
		bbox_keep,
		dgcdp->profile_boolweave_source_intervals,
		dgcdp->profile_boolweave_source_only,
		dgcdp->profile_boolweave_shot_calls,
		dgcdp->profile_boolweave_shot_hits,
		dgcdp->profile_boolweave_shot_misses,
		dgcdp->profile_boolweave_partitions,
		boolweave_shoot_sec,
		boolweave_weave_sec,
		boolweave_partition_sec);
    }

    for (size_t i = 0; i < dgcdp->profile_top_regions.size(); i++) {
	const bigE_profile_region_summary &summary =
	    dgcdp->profile_top_regions[i];
	bu_log("evaluated-wire profile top-region: rank=%zu index=%zu "
		"name=\"%s\" time=%.6fs candidates=%zu edge=%zu "
		"intersection=%zu face_pairs=%zu "
		"segments=(bige:%zu boolweave:%zu) "
		"leaves=%zu active_soltabs=%zu "
		"intersection_time=(index:%.6fs collect:%.6fs eval:%.6fs) "
		"candidate_eval=%.6fs\n",
		i + 1,
		summary.index,
		summary.name.empty() ? "(null)" : summary.name.c_str(),
		((double)summary.usec) / 1000000.0,
		summary.candidates,
		summary.edge_candidates,
		summary.intersection_candidates,
		summary.intersection_face_pairs,
		summary.bige_segments,
		summary.boolweave_segments,
		summary.leaves,
		summary.active_soltabs,
		((double)summary.intersection_face_index_usec) / 1000000.0,
		((double)summary.intersection_pair_collect_usec) / 1000000.0,
		((double)summary.intersection_pair_eval_usec) / 1000000.0,
		((double)summary.candidate_eval_usec) / 1000000.0);
    }
}


static void
bigE_context_cleanup(struct bigE_data *dgcdp)
{
    if (!dgcdp)
	return;

    bigE_profile_report(dgcdp);
    bigE_stats_report(dgcdp);
    bigE_boolweave_eval_report(dgcdp);
    bigE_boolweave_diag_report(dgcdp);
    bigE_line_free(dgcdp);
    dgcdp->active_region_soltabs.clear();
    dgcdp->active_soltab_local_index.clear();
    dgcdp->active_csg_nodes.clear();
    dgcdp->active_source_indices_by_leaf.clear();
    dgcdp->active_candidate_indices_by_leaf.clear();
    dgcdp->active_candidate_indices_by_pair.clear();
    dgcdp->active_csg_eval_cache.clear();
    dgcdp->active_csg_root = -1;

    if (dgcdp->leaf_list_initialized) {
	bu_ptbl_free(&dgcdp->leaf_list);
	dgcdp->leaf_list_initialized = 0;
    }
    if (dgcdp->resource_initialized) {
	rt_clean_resource(dgcdp->rtip, &dgcdp->resource);
	dgcdp->resource_initialized = 0;
    }
    if (dgcdp->rtip) {
	rt_i_destroy(dgcdp->rtip);
	dgcdp->rtip = NULL;
    }
    if (dgcdp->ap) {
	bu_free(dgcdp->ap, "struct application");
	dgcdp->ap = NULL;
    }
    bigE_vlfree_clear(&dgcdp->vlfree_storage);
}


static int
bigE_context_prepare(struct bigE_data *dgcdp,
		     struct db_i *dbip,
		     const char *path,
		     const struct bn_tol *tol,
		     const struct bg_tess_tol *ttol,
		     uint32_t eval_flags)
{
    if (!dgcdp || !dbip || !path || !path[0] || !tol || !ttol)
	return BRLCAD_ERROR;

    bigE_context_cleanup(dgcdp);
    bigE_profile_reset(dgcdp);

    dgcdp->dbip = dbip;
    dgcdp->do_polysolids = 0;
    dgcdp->fp = NULL;
    dgcdp->tol = tol;
    dgcdp->ttol = ttol;
    dgcdp->vlfree = &dgcdp->vlfree_storage;
    BU_LIST_INIT(&dgcdp->vlfree_storage);
    dgcdp->active_region = NULL;
    bigE_boolweave_diag_reset(dgcdp);
    dgcdp->candidate_index = 0;
    dgcdp->eval_flags = eval_flags;

    const char *bw_diag = getenv("BRLOBOL_EVAL_WIRE_BOOLWEAVE_DIAG");
    dgcdp->boolweave_diag_enabled = (bw_diag && bw_diag[0] && bu_strcmp(bw_diag, "0"));
    const char *bige_stats = getenv("BRLOBOL_EVAL_WIRE_BIGE_STATS");
    dgcdp->bige_stats_enabled = (bige_stats && bige_stats[0] && bu_strcmp(bige_stats, "0"));
    const uint32_t method = eval_flags & RT_EVAL_WIREFRAME_F_METHOD_MASK;
    dgcdp->boolweave_eval_enabled =
	(method != RT_EVAL_WIREFRAME_F_BIGE);
    dgcdp->boolweave_eval_filter_enabled =
	bigE_env_enabled("BRLOBOL_EVAL_WIRE_BOOLWEAVE_FILTER");
    dgcdp->profile_enabled =
	((eval_flags & RT_EVAL_WIREFRAME_F_PROFILE) ||
	 bigE_env_enabled("RT_EVAL_WIREFRAME_PROFILE") ||
	 bigE_env_enabled("BRLOBOL_EVAL_WIRE_PROFILE"));
    dgcdp->profile_region_report_enabled =
	bigE_env_enabled("RT_EVAL_WIREFRAME_PROFILE_REGIONS") ||
	bigE_env_enabled("BRLOBOL_EVAL_WIRE_PROFILE_REGIONS");
    dgcdp->profile_progress_enabled =
	dgcdp->profile_enabled &&
	bigE_env_enabled("RT_EVAL_WIREFRAME_PROFILE_PROGRESS");
    dgcdp->profile_top_limit =
	bigE_env_size("RT_EVAL_WIREFRAME_PROFILE_TOP", dgcdp->profile_top_limit);
    dgcdp->boolweave_diag_limit = 64;
    if (dgcdp->boolweave_diag_enabled) {
	const char *bw_limit = getenv("BRLOBOL_EVAL_WIRE_BOOLWEAVE_DIAG_LIMIT");
	if (bw_limit && bw_limit[0])
	    dgcdp->boolweave_diag_limit = (size_t)strtoull(bw_limit, NULL, 10);
    }

    BU_ALLOC(dgcdp->ap, struct application);
    RT_APPLICATION_INIT(dgcdp->ap);

    dgcdp->rtip = rt_i_create(dgcdp->dbip);
    if (!dgcdp->rtip)
	return BRLCAD_ERROR;

    dgcdp->rtip->rti_tol = *dgcdp->tol;
    dgcdp->rtip->useair = 1;
    dgcdp->ap->a_rt_i = dgcdp->rtip;
    rt_init_resource(&dgcdp->resource, 0, dgcdp->rtip);
    dgcdp->resource_initialized = 1;
    dgcdp->ap->a_resource = &dgcdp->resource;
    dgcdp->nvectors = 0;

    bu_ptbl_init(&dgcdp->leaf_list, 8, "leaf_list");
    dgcdp->leaf_list_initialized = 1;

    const char *paths[1] = {path};
    if (rt_gettrees(dgcdp->rtip, 1, paths, 1))
	return BRLCAD_ERROR;

    return BRLCAD_OK;
}


static struct region *
bigE_region_at(struct rt_i *rtip, size_t index)
{
    struct region *rp;
    size_t i = 0;

    if (!rtip)
	return NULL;

    for (BU_LIST_FOR (rp, region, &(rtip->HeadRegion))) {
	if (i == index)
	    return rp;
	i++;
    }

    return NULL;
}


static const char *
bigE_region_debug_name(const struct region *rp)
{
    if (!rp || !rp->reg_name)
	return NULL;
    return rp->reg_name;
}


static int
bigE_region_count_discover(struct db_i *dbip,
			   const char *path,
			   const struct bn_tol *tol,
			   const struct bg_tess_tol *ttol,
			   uint32_t eval_flags,
			   size_t *region_count)
{
    if (!region_count)
	return BRLCAD_ERROR;

    *region_count = 0;

    struct bigE_data discovery;
    if (bigE_context_prepare(&discovery, dbip, path, tol, ttol,
	    eval_flags) != BRLCAD_OK)
	return BRLCAD_ERROR;

    struct region *rp;
    for (BU_LIST_FOR (rp, region, &(discovery.rtip->HeadRegion))) {
	(void)rp;
	(*region_count)++;
    }

    return BRLCAD_OK;
}


static int
bigE_region_evaluate(struct bigE_data *dgcdp,
		     struct region *rp,
		     size_t region_index,
		     struct bigE_line_set *lines_out)
{
    if (!dgcdp || !rp || !rp->reg_treetop || !lines_out)
	return BRLCAD_ERROR;

    bigE_profile_region_start(dgcdp, region_index, rp);
    dgcdp->lines.records.clear();
    dgcdp->active_region_soltabs.clear();
    dgcdp->num_halfs = 0;
    dgcdp->active_region = rp;
    bigE_collect_region_soltabs(rp->reg_treetop,
	    &dgcdp->active_region_soltabs);
    bigE_build_active_soltab_index(dgcdp);
    dgcdp->active_csg_nodes.clear();
    dgcdp->active_csg_eval_cache.clear();
    dgcdp->active_csg_root =
	bigE_compile_active_csg_tree(dgcdp, rp->reg_treetop);
    if (dgcdp->profile_enabled)
	dgcdp->profile_current.active_soltabs =
	    dgcdp->active_region_soltabs.size();

    union E_tree *eptr = build_etree(rp->reg_treetop, dgcdp);
    if (!eptr) {
	dgcdp->active_region_soltabs.clear();
	dgcdp->active_soltab_local_index.clear();
	dgcdp->active_csg_nodes.clear();
	dgcdp->active_source_indices_by_leaf.clear();
	dgcdp->active_candidate_indices_by_leaf.clear();
	dgcdp->active_candidate_indices_by_pair.clear();
	dgcdp->active_csg_eval_cache.clear();
	dgcdp->active_csg_root = -1;
	dgcdp->active_region = NULL;
	bigE_profile_region_finish(dgcdp);
	return BRLCAD_ERROR;
    }
    if (dgcdp->profile_enabled)
	dgcdp->profile_current.leaves = bigE_leaf_count(dgcdp);

    if (dgcdp->num_halfs)
	fix_halfs(dgcdp);
    bigE_build_active_source_index(dgcdp);

    Eplot(eptr, dgcdp);
    free_etree(eptr, dgcdp);
    bu_ptbl_reset(&dgcdp->leaf_list);
    dgcdp->active_region_soltabs.clear();
    dgcdp->active_soltab_local_index.clear();
    dgcdp->active_csg_nodes.clear();
    dgcdp->active_source_indices_by_leaf.clear();
    dgcdp->active_candidate_indices_by_leaf.clear();
    dgcdp->active_candidate_indices_by_pair.clear();
    dgcdp->active_csg_eval_cache.clear();
    dgcdp->active_csg_root = -1;
    dgcdp->active_region = NULL;

    lines_out->records = dgcdp->lines.records;
    dgcdp->lines.records.clear();
    bigE_profile_region_finish(dgcdp);
    return BRLCAD_OK;
}


static union E_tree *
add_solid(const struct directory *dp,
	  matp_t mat,
	  struct bigE_data *dgcdp)
{
    union E_tree *eptr = NULL;
    struct nmgregion *r = NULL;
    struct rt_db_internal intern;
    int id;
    int solid_is_plate_mode_bot=0;

    eptr = bigE_tree_create();

    id = rt_db_get_internal(&intern, dp, dgcdp->dbip, mat);
    if (id < 0) {
	eptr->l.m = (struct model *)NULL;
	return eptr;
    }
    if (id == ID_COMBINATION) {
	/* do explicit expansion of referenced combinations */

	struct rt_comb_internal *comb;

	bu_free((char *)eptr, "eptr");

	comb = (struct rt_comb_internal *)intern.idb_ptr;
	RT_CK_COMB(comb);

	eptr = build_etree(comb->tree, dgcdp);
	rt_db_free_internal(&intern);
	return eptr;
    }
    if (id == ID_HALF) {
	eptr->l.m = NULL;
	dgcdp->num_halfs++;
    } else if (id == ID_NMG) {
	/* steal the nmg model */
	eptr->l.m = (struct model *)intern.idb_ptr;
	eptr->l.do_not_free_model = 1;
    } else {
	/* create the NMG version of this solid */
	eptr->l.m = nmg_mm();

	if (!OBJ[id].ft_tessellate ||
	    OBJ[id].ft_tessellate(&r, eptr->l.m, &intern,
					 dgcdp->ttol,
					 dgcdp->tol) < 0)
	{
	    nmg_km(eptr->l.m);
	    eptr->l.m = NULL;
	}
    }

    /* get the soltab stuff */
    eptr->l.stp = bigE_soltab_create(dp, mat);

    struct rt_bot_internal *bot = NULL;

    if (dgcdp->do_polysolids) {
	struct shell *s=(struct shell *)NULL;

	/* create and prep a BoT version of this solid */
	if (eptr->l.m) {
	    r = BU_LIST_FIRST(nmgregion, &eptr->l.m->r_hd);
	    s = BU_LIST_FIRST(shell, &r->s_hd);
	}

	if (solid_is_plate_mode_bot
		|| !eptr->l.m
		|| (bot = nmg_bot(s, dgcdp->vlfree, dgcdp->tol)) == (struct rt_bot_internal *)NULL)
	{
	    eptr->l.stp->st_id = id;
	    eptr->l.stp->st_meth = &OBJ[id];
	    rt_obj_prep(eptr->l.stp, &intern, dgcdp->rtip);
	} else {
	    struct rt_db_internal intern2;
	    RT_DB_INTERNAL_INIT(&intern2);
	    intern2.idb_major_type = DB5_MAJORTYPE_BRLCAD;
	    intern2.idb_type = ID_BOT;
	    intern2.idb_meth = &OBJ[ID_BOT];
	    intern2.idb_ptr = (void *)bot;
	    eptr->l.stp->st_id = ID_BOT;
	    eptr->l.stp->st_meth = &OBJ[ID_BOT];
	    rt_obj_prep(eptr->l.stp, &intern2, dgcdp->rtip);
	    rt_db_free_internal(&intern2);
	}

    } else {
	/* prep this solid */

	eptr->l.stp->st_id = id;
	eptr->l.stp->st_meth = &OBJ[id];
	rt_obj_prep(eptr->l.stp, &intern, dgcdp->rtip);
    }

    if (id != ID_NMG)
	rt_db_free_internal(&intern);

    /* add this leaf to the leaf list */
    bu_ptbl_ins(&dgcdp->leaf_list, (long *)eptr);

    return eptr;
}


/* build an E_tree corresponding to the region tree (tp) */
static union E_tree *
build_etree(union tree *tp, struct bigE_data *dgcdp)
{
    union E_tree *eptr = NULL;
    struct soltab *stp;
    struct directory *dp;

    RT_CK_TREE(tp);

    switch (tp->tr_op) {
	case OP_UNION:
	case OP_SUBTRACT:
	case OP_INTERSECT:
	    eptr = bigE_tree_create();
	    eptr->n.op = tp->tr_op;
	    eptr->n.left = build_etree(tp->tr_b.tb_left, dgcdp);
	    eptr->n.right = build_etree(tp->tr_b.tb_right, dgcdp);
	    break;
	case OP_SOLID:
	    stp = tp->tr_a.tu_stp;
	    eptr = add_solid(stp->st_dp, stp->st_matp, dgcdp);
	    eptr->l.op = tp->tr_op;
	    BU_LIST_INIT(&eptr->l.seghead);
	    break;
	case OP_DB_LEAF:
	    dp = db_lookup(dgcdp->dbip, tp->tr_l.tl_name, LOOKUP_NOISY);
	    if (dp == RT_DIR_NULL) {
	      break;
	    }
	    eptr = add_solid(dp, tp->tr_l.tl_mat, dgcdp);
	    eptr->l.op = tp->tr_op;
	    BU_LIST_INIT(&eptr->l.seghead);
	    break;
	case OP_NOP:
	    /* add a NULL solid */
	    eptr = bigE_tree_create();
	    eptr->l.m = (struct model *)NULL;
	    break;
	default:
	    bu_bomb("build_etree() Unknown tr_op\n");
    }
    return eptr;
}


/* a handy routine (for debugging) that prints a segment list */
void
show_seg(struct bu_list *seg, int str)
{
    struct seg *ptr;

    if (!seg)
	bu_log("%d - NULL seg list\n", str);
    else {
	if (BU_LIST_IS_EMPTY(seg))
	    bu_log("%d - empty\n", str);
	else {
	    bu_log("%d:\n", str);
	    for (BU_LIST_FOR (ptr, seg, seg)) {
		if (ptr->seg_stp == ON_SURF)
		    bu_log("\t %g to %g (ON_SURF)\n", ptr->seg_in.hit_dist, ptr->seg_out.hit_dist);
		else if (ptr->seg_stp == ON_INT)
		    bu_log("\t %g to %g (ON_INT)\n", ptr->seg_in.hit_dist, ptr->seg_out.hit_dist);
		else if (ptr->seg_stp == IN_SOL)
		    bu_log("\t %g to %g (IN)\n", ptr->seg_in.hit_dist, ptr->seg_out.hit_dist);
		else
		    bu_log("\t %g to %g (?)\n", ptr->seg_in.hit_dist, ptr->seg_out.hit_dist);
		bu_log("\t\t(%g %g %g) <-> (%g %g %g)\n", V3ARGS(ptr->seg_in.hit_point),
		       V3ARGS(ptr->seg_out.hit_point));
	    }
	}
    }
}


/* given a segment list, eliminate any overlaps in the segments */
static void
eliminate_overlaps(struct bu_list *seghead,
		   struct bigE_data *dgcdp)
{
    struct seg *a, *b, *nextb;

    a = BU_LIST_FIRST(seg, seghead);
    while (BU_LIST_NOT_HEAD(&a->l, seghead)) {
	b = BU_LIST_PNEXT(seg, &a->l);
	if (BU_LIST_IS_HEAD(&b->l, seghead))
	    break;

	while (BU_LIST_NOT_HEAD(&b->l, seghead)) {
	    nextb = BU_LIST_PNEXT(seg, &b->l);
	    if (NOT_SEG_OVERLAP(a, b))
		break;

	    if (b->seg_in.hit_dist < a->seg_out.hit_dist) {
		if (b->seg_out.hit_dist > a->seg_out.hit_dist)
		    a->seg_out.hit_dist = b->seg_out.hit_dist;

		BU_LIST_DEQUEUE(&b->l);
		bigE_seg_free(dgcdp, b);
		b = nextb;
		continue;
	    }

	    b = nextb;
	}

	a = BU_LIST_PNEXT(seg, &a->l);
    }
}


/* perform the intersection of two segments the result is assigned the
 * provided type
 */
static void
do_intersect(struct seg *A,
	     struct seg *B,
	     struct bu_list *seghead,
	     struct soltab *type,
	     struct bigE_data *dgcdp)
{
    struct seg *tmp=(struct seg *)NULL;

    if (NOT_SEG_OVERLAP(A, B))
	return;

    tmp = bigE_seg_get(dgcdp);
    if (A->seg_in.hit_dist <= B->seg_in.hit_dist) {
	if (B->seg_out.hit_dist <= A->seg_out.hit_dist) {
	    *tmp = *B;
	    tmp->seg_stp = type;
	} else {
	    tmp->seg_in.hit_dist = B->seg_in.hit_dist;
	    tmp->seg_out.hit_dist = A->seg_out.hit_dist;
	    tmp->seg_stp = type;
	}
    } else {
	if (B->seg_out.hit_dist >= A->seg_out.hit_dist) {
	    *tmp = *A;
	    tmp->seg_stp = type;
	} else {
	    tmp->seg_in.hit_dist = A->seg_in.hit_dist;
	    tmp->seg_out.hit_dist = B->seg_out.hit_dist;
	    tmp->seg_stp = type;
	}
    }
    if (tmp) {
	BU_LIST_INSERT(seghead, &tmp->l);
    }

    return;
}


/* perform the subtraction of one segment from another the result is
 * assigned the type from segment A
 */
static void
do_subtract(struct seg *A,
	    struct seg *B,
	    struct bu_list *seghead,
	    struct bigE_data *dgcdp)
{
    struct seg *tmp=(struct seg *)NULL;

    if (NOT_SEG_OVERLAP(A, B)) {
	tmp = bigE_seg_get(dgcdp);
	*tmp = *A;
	BU_LIST_INSERT(seghead, &tmp->l);
	return;
    }

    if (B->seg_in.hit_dist<= A->seg_in.hit_dist) {
	if (B->seg_out.hit_dist < A->seg_out.hit_dist) {
	    tmp = bigE_seg_get(dgcdp);
	    *tmp = *A;
	    tmp->seg_in.hit_dist = B->seg_out.hit_dist;
	    BU_LIST_INSERT(seghead, &tmp->l);
	    return;
	} else {
	    return;
	}
    } else {
	if (B->seg_out.hit_dist >= A->seg_out.hit_dist) {
	    tmp = bigE_seg_get(dgcdp);
	    *tmp = *A;
	    tmp->seg_out.hit_dist = B->seg_in.hit_dist;
	    BU_LIST_INSERT(seghead, &tmp->l);
	    return;
	} else {
	    tmp = bigE_seg_get(dgcdp);
	    tmp->seg_in.hit_dist = A->seg_in.hit_dist;
	    tmp->seg_out.hit_dist = B->seg_in.hit_dist;
	    tmp->seg_stp = A->seg_stp;
	    BU_LIST_INSERT(seghead, &tmp->l);
	    tmp = bigE_seg_get(dgcdp);
	    tmp->seg_in.hit_dist = B->seg_out.hit_dist;
	    tmp->seg_out.hit_dist = A->seg_out.hit_dist;
	    tmp->seg_stp = A->seg_stp;
	    BU_LIST_INSERT(seghead, &tmp->l);
	    return;
	}
    }
}


static void
promote_ints(struct bu_list *head,
	     struct bigE_data *dgcdp)
{
    struct seg *a, *b, *tmp;

    a = BU_LIST_FIRST(seg, head);
    while (BU_LIST_NOT_HEAD(&a->l, head)) {
	b = BU_LIST_PNEXT(seg, &a->l);
	while (BU_LIST_NOT_HEAD(&b->l, head)) {
	    if (a->seg_stp == ON_INT && b->seg_stp == ON_SURF) {
		if (NOT_SEG_OVERLAP(a, b)) {
		    b = BU_LIST_PNEXT(seg, &b->l);
		    continue;
		}

		if (ZERO(a->seg_in.hit_dist - b->seg_in.hit_dist)
		    && ZERO(a->seg_out.hit_dist - b->seg_out.hit_dist))
		{
		    a->seg_stp = ON_SURF;
		    tmp = b;
		    b = BU_LIST_PNEXT(seg, &b->l);
		    BU_LIST_DEQUEUE(&tmp->l);
		    bigE_seg_free(dgcdp, tmp);
		    continue;
		}

		if (ZERO(a->seg_out.hit_dist - b->seg_out.hit_dist))
		    a->seg_out.hit_dist = b->seg_in.hit_dist;
		else if (a->seg_out.hit_dist < b->seg_out.hit_dist) {
		    if (b->seg_in.hit_dist > a->seg_in.hit_dist)
			a->seg_out.hit_dist = b->seg_in.hit_dist;
		    else {
			tmp = a;
			a  = BU_LIST_PLAST(seg, &a->l);
			BU_LIST_DEQUEUE(&tmp->l);
			bigE_seg_free(dgcdp, tmp);
			break;
		    }
		} else if (ZERO(a->seg_in.hit_dist - b->seg_in.hit_dist)) {
		    fastf_t tmp_dist;

		    tmp_dist = a->seg_out.hit_dist;
		    a->seg_out.hit_dist = b->seg_out.hit_dist;
		    b->seg_in.hit_dist = a->seg_out.hit_dist;
		    b->seg_out.hit_dist = tmp_dist;
		    a->seg_stp = ON_SURF;
		    b->seg_stp = ON_INT;
		} else {
		    tmp = bigE_seg_get(dgcdp);
		    *tmp = *a;
		    tmp->seg_in.hit_dist = b->seg_out.hit_dist;
		    a->seg_out.hit_dist = b->seg_in.hit_dist;
		    BU_LIST_APPEND(&b->l, &tmp->l);
		}
	    } else if (b->seg_stp == ON_INT && a->seg_stp == ON_SURF) {
		if (NOT_SEG_OVERLAP(b, a)) {
		    b = BU_LIST_PNEXT(seg, &b->l);
		    continue;
		}

		if (ZERO(b->seg_in.hit_dist - a->seg_in.hit_dist)
		    && ZERO(b->seg_out.hit_dist - a->seg_out.hit_dist))
		{
		    b->seg_stp = ON_SURF;
		    tmp = a;
		    a = BU_LIST_PLAST(seg, &a->l);
		    BU_LIST_DEQUEUE(&tmp->l);
		    bigE_seg_free(dgcdp, tmp);
		    break;
		}

		if (ZERO(b->seg_out.hit_dist - a->seg_out.hit_dist)) {
		    tmp = b;
		    b = BU_LIST_PNEXT(seg, &b->l);
		    BU_LIST_DEQUEUE(&tmp->l);
		    bigE_seg_free(dgcdp, tmp);
		} else if (b->seg_out.hit_dist < a->seg_out.hit_dist) {
		    if (a->seg_in.hit_dist > b->seg_in.hit_dist)
			b->seg_out.hit_dist = a->seg_in.hit_dist;
		    else {
			tmp = b;
			b = BU_LIST_PNEXT(seg, &b->l);
			BU_LIST_DEQUEUE(&tmp->l);
			bigE_seg_free(dgcdp, tmp);
			continue;
		    }
		} else if (ZERO(b->seg_in.hit_dist - a->seg_in.hit_dist)) {
		    b->seg_in.hit_dist = a->seg_out.hit_dist;
		} else {
		    tmp = bigE_seg_get(dgcdp);
		    *tmp = *b;
		    tmp->seg_in.hit_dist = a->seg_out.hit_dist;
		    b->seg_out.hit_dist = a->seg_in.hit_dist;
		    BU_LIST_APPEND(&a->l, &tmp->l);
		}
	    }

	    if ((a->seg_stp != ON_INT) || (b->seg_stp != ON_INT)) {
		b = BU_LIST_PNEXT(seg, &b->l);
		continue;
	    }

	    if (NOT_SEG_OVERLAP(a, b)) {
		b = BU_LIST_PNEXT(seg, &b->l);
		continue;
	    }

	    if (ZERO(a->seg_in.hit_dist - b->seg_in.hit_dist)
		&& ZERO(a->seg_out.hit_dist - b->seg_out.hit_dist))
	    {
		a->seg_stp = ON_SURF;
		BU_LIST_DEQUEUE(&b->l);
		bigE_seg_free(dgcdp, b);
		break;
	    }

	    if (ZERO(a->seg_out.hit_dist - b->seg_out.hit_dist)) {
		b->seg_stp = ON_SURF;
		a->seg_out.hit_dist = b->seg_in.hit_dist;
	    } else if (a->seg_out.hit_dist < b->seg_out.hit_dist) {
		if (b->seg_in.hit_dist > a->seg_in.hit_dist) {
		    tmp = bigE_seg_get(dgcdp);
		    tmp->seg_stp = ON_SURF;
		    tmp->seg_in.hit_dist = b->seg_in.hit_dist;
		    tmp->seg_out.hit_dist = a->seg_out.hit_dist;
		    b->seg_in.hit_dist = a->seg_out.hit_dist;
		    a->seg_out.hit_dist = tmp->seg_in.hit_dist;
		    BU_LIST_INSERT(&b->l, &tmp->l);
		} else {
		    b->seg_in.hit_dist = a->seg_out.hit_dist;
		    a->seg_stp = ON_SURF;
		}
	    } else {
		if (ZERO(a->seg_in.hit_dist - b->seg_in.hit_dist)) {
		    fastf_t tmp_dist;

		    tmp_dist = a->seg_out.hit_dist;
		    a->seg_out.hit_dist = b->seg_out.hit_dist;
		    a->seg_stp = ON_SURF;
		    b->seg_in.hit_dist = a->seg_out.hit_dist;
		    b->seg_out.hit_dist = tmp_dist;
		} else {
		    tmp = bigE_seg_get(dgcdp);
		    *tmp = *a;
		    tmp->seg_in.hit_dist = b->seg_out.hit_dist;
		    a->seg_out.hit_dist = b->seg_in.hit_dist;
		    b->seg_stp = ON_SURF;
		    BU_LIST_APPEND(&b->l, &tmp->l);
		}
	    }
	    b = BU_LIST_PNEXT(seg, &b->l);
	}
	a = BU_LIST_PNEXT(seg, &a->l);
    }

#ifdef debug
    bu_log("Results of promote_ints()\n");
    show_seg(head, "SEGS");
#endif
}


/* Evaluate an operation on the operands (segment lists) */
static struct bu_list *
eval_op(struct bu_list *A,
	int op,
	struct bu_list *B,
	struct bigE_data *dgcdp)
{
    struct seg *sega, *segb, *tmp, *next;
    struct bu_list ret, ons, ins;
    int inserted;

    BU_LIST_INIT(&ret);

    switch (op) {
	case OP_SUBTRACT:
	    if (BU_LIST_IS_EMPTY(A)) {
		bigE_seg_list_free(B, dgcdp);
		bu_free((char *)B, "bu_list");
		return A;
	    } else if (BU_LIST_IS_EMPTY(B)) {
		bu_free((char *)B, "bu_list");
		return A;
	    }

	    /* A - B:
	     * keep segments:
	     * ON_A - IN_B
	     * ON_A - ON_B
	     * ON_B + IN_A
	     * IN_A - IN_B
	     */
	    for (BU_LIST_FOR (sega, seg, A)) {
		for (BU_LIST_FOR (segb, seg, B)) {
		    if (sega->seg_stp == ON_INT && segb->seg_stp == ON_INT) {
			do_intersect(sega, segb, &ret, ON_SURF, dgcdp);
		    } else if (sega->seg_stp == ON_SURF || sega->seg_stp == ON_INT) {
			do_subtract(sega, segb, &ret, dgcdp);
		    } else if (segb->seg_stp == ON_SURF ||  segb->seg_stp == ON_INT) {
			do_intersect(segb, sega, &ret, segb->seg_stp, dgcdp);
		    }
		}
	    }
	    bigE_seg_list_free(B, dgcdp);
	    bu_free((char *)B, "bu_list");
	    bigE_seg_list_free(A, dgcdp);
	    BU_LIST_INSERT_LIST(A, &ret);
	    return A;
	case OP_INTERSECT:
	    if (BU_LIST_IS_EMPTY(A) || BU_LIST_IS_EMPTY(B)) {
		bigE_seg_list_free(A, dgcdp);
		bigE_seg_list_free(B, dgcdp);
		bu_free((char *)B, "bu_list");
		return A;
	    }
	    /* A + B
	     *
	     * This is merely the intersection of segments from A with
	     * those from B.  The two different calls to
	     * "do_intersect" get the types (IN, ON) right
	     */
	    for (BU_LIST_FOR (sega, seg, A)) {
		for (BU_LIST_FOR (segb, seg, B)) {
		    if (sega->seg_stp == ON_INT && segb->seg_stp == ON_INT)
			do_intersect(sega, segb, &ret, ON_SURF, dgcdp);
		    else if (sega->seg_stp == ON_SURF || sega->seg_stp == ON_INT)
			do_intersect(sega, segb, &ret, sega->seg_stp, dgcdp);
		    else
			do_intersect(segb, sega, &ret, segb->seg_stp, dgcdp);
		}
	    }
	    bigE_seg_list_free(B, dgcdp);
	    bu_free((char *)B, "bu_list");
	    bigE_seg_list_free(A, dgcdp);
	    BU_LIST_INSERT_LIST(A, &ret)
	    return A;
	case OP_UNION:
	    if (BU_LIST_IS_EMPTY(A)) {
		bu_free((char *)A, "bu_list");
		return B;
	    }
	    if (BU_LIST_IS_EMPTY(B)) {
		bu_free((char *)B, "bu_list");
		return A;
	    }
	    /* A u B:
	     * keep segments:
	     * ON_A - IN_B (ON)
	     * IN_B + ON_A (IN)
	     * ON_B - IN_A (ON)
	     * IN_A + ON_B (IN)
	     * all remaining unique ON or IN segments
	     */

	    /* create two new lists, one with all the ON segments, the
	     * other with all the IN segments
	     */
	    BU_LIST_INIT(&ons);
	    BU_LIST_INIT(&ins);

	    /* Put the A operand segments on the lists */
	    while (BU_LIST_WHILE (sega, seg, A)) {
		BU_LIST_DEQUEUE(&sega->l);

		if (sega->seg_stp == ON_SURF || sega->seg_stp == ON_INT) {
		    BU_LIST_INSERT(&ons, &sega->l);
		} else {
		    BU_LIST_INSERT(&ins, &sega->l);
		}
	    }

	    /* insert the B operand segments in the lists (maintaining
	     * order from smaller starting hit distance to larger
	     */
	    while (BU_LIST_WHILE (segb, seg, B)) {
		BU_LIST_DEQUEUE(&segb->l);

		if (segb->seg_stp == IN_SOL) {
		    inserted = 0;
		    for (BU_LIST_FOR (tmp, seg, &ins)) {
			if (tmp->seg_in.hit_dist >= segb->seg_in.hit_dist) {
			    inserted = 1;
			    BU_LIST_INSERT(&tmp->l, &segb->l);
			    break;
			}
		    }
		    if (!inserted)
			BU_LIST_INSERT(&ins, &segb->l);
		} else {
		    inserted = 0;
		    for (BU_LIST_FOR (tmp, seg, &ons)) {
			if (tmp->seg_in.hit_dist >= segb->seg_in.hit_dist) {
			    inserted = 1;
			    BU_LIST_INSERT(&tmp->l, &segb->l);
			    break;
			}
		    }
		    if (!inserted)
			BU_LIST_INSERT(&ons, &segb->l);
		}
	    }

	    /* promote intersecting ON_INT's to ON_SURF */
	    promote_ints(&ons, dgcdp);

	    /* make sure the segments are unique */
	    eliminate_overlaps(&ins, dgcdp);
	    eliminate_overlaps(&ons, dgcdp);

	    /* subtract INS from ONS */
	    sega = BU_LIST_FIRST(seg, &ons);
	    while (BU_LIST_NOT_HEAD(&sega->l, &ons)) {
		next = BU_LIST_PNEXT(seg, &sega->l);
		for (BU_LIST_FOR (segb, seg, &ins)) {
		    if (NOT_SEG_OVERLAP(sega, segb)) {
			continue;
		    }

		    if (segb->seg_in.hit_dist <= sega->seg_in.hit_dist &&
			segb->seg_out.hit_dist >= sega->seg_out.hit_dist) {
			/* eliminate sega */
			BU_LIST_DEQUEUE(&sega->l);
			bigE_seg_free(dgcdp, sega);
			break;
		    }

		    if (segb->seg_in.hit_dist > sega->seg_in.hit_dist &&
			segb->seg_out.hit_dist < sega->seg_out.hit_dist)
		    {
			/* split sega */
			tmp = bigE_seg_get(dgcdp);
			*tmp = *sega;
			tmp->seg_in.hit_dist = segb->seg_out.hit_dist;
			sega->seg_out.hit_dist = segb->seg_in.hit_dist;
			BU_LIST_APPEND(&sega->l, &tmp->l);
			next = tmp;
		    } else {
			/* subtract edges */
			if (segb->seg_in.hit_dist > sega->seg_in.hit_dist)
			    sega->seg_out.hit_dist = segb->seg_in.hit_dist;
			if (segb->seg_out.hit_dist < sega->seg_out.hit_dist)
			    sega->seg_in.hit_dist = segb->seg_out.hit_dist;
		    }
		}
		sega = next;
	    }

	    /* put the resulting ONS list on the result list */
	    BU_LIST_INSERT_LIST(&ret, &ons);

	    /* add INS to the return list (maintain order) */
	    while (BU_LIST_WHILE (sega, seg, &ins)) {
		BU_LIST_DEQUEUE(&sega->l);

		inserted = 0;
		for (BU_LIST_FOR (segb, seg, &ret)) {
		    if (sega->seg_in.hit_dist < segb->seg_in.hit_dist) {
			BU_LIST_INSERT(&segb->l, &sega->l);
			inserted = 1;
			break;
		    }
		}

		if (!inserted)
		    BU_LIST_INSERT(&ret, &sega->l)
			}

	    bigE_seg_list_free(B, dgcdp);
	    bu_free((char *)B, "bu_list");
	    bigE_seg_list_free(A, dgcdp);
	    BU_LIST_INSERT_LIST(A, &ret)
	    return A;
    }

    /* should never get here */
    bigE_seg_list_free(A, dgcdp);
    bigE_seg_list_free(B, dgcdp);
    bu_free((char *)B, "bu_list");
    return A;
}


/* evaluate an E-tree */
static struct bu_list *
eval_etree(union E_tree *eptr,
	   struct bigE_data *dgcdp)

{
    struct bu_list *A, *B;

    CK_ETREE(eptr);

    switch (eptr->l.op) {
	case OP_DB_LEAF:
	case OP_SOLID:
	    BU_ALLOC(A, struct bu_list);
	    BU_LIST_INIT(A);
	    BU_LIST_INSERT_LIST(A, &eptr->l.seghead);
	    return A;
	case OP_SUBTRACT:
	case OP_INTERSECT:
	case OP_UNION:
	    A = eval_etree(eptr->n.left, dgcdp);
	    B = eval_etree(eptr->n.right, dgcdp);
	    return eval_op(A, eptr->n.op, B, dgcdp);
    }

    /* should never get here */
    return (struct bu_list *)NULL;	/* for the compilers */
}


static struct soltab *
classify_seg(struct seg *segp, struct soltab *shoot, struct xray *rp, struct bigE_data *dgcdp)
{
    fastf_t mid_dist;
    struct xray new_rp;
    struct ray_data rd;
    struct soltab *ret = IN_SOL;

    memset(&rd, 0, sizeof(struct ray_data));

    rd.seghead = bigE_ray_seghead_create();

    mid_dist = (segp->seg_in.hit_dist + segp->seg_out.hit_dist) / 2.0;
    VJOIN1(new_rp.r_pt, rp->r_pt, mid_dist, rp->r_dir);
    bn_vec_ortho(new_rp.r_dir, rp->r_dir);
    /* Compute the inverse of the direction cosines */
    VINVDIR(rd.rd_invdir, new_rp.r_dir);

    /* set up "ray_data" structure for nmg raytrace */
    rd.rp = &new_rp;
    rd.tol = dgcdp->tol;
    rd.ap = dgcdp->ap;
    rd.magic = NMG_RAY_DATA_MAGIC;
    rd.classifying_ray = 0;
    rd.hitmiss = (struct hitmiss **)NULL;
    rd.stp = shoot;

    if (OBJ[shoot->st_id].ft_shot && OBJ[shoot->st_id].ft_shot(shoot, &new_rp, dgcdp->ap, rd.seghead)) {
	struct seg *seg;

	while (BU_LIST_WHILE (seg, seg, &rd.seghead->l)) {
	    BU_LIST_DEQUEUE(&seg->l);
	    if (ret != ON_SURF) {
		if (NEAR_ZERO(seg->seg_in.hit_dist, rd.tol->dist)) {
		    ret = ON_SURF;
		}
		if (NEAR_ZERO(seg->seg_out.hit_dist, rd.tol->dist)) {
		    ret = ON_SURF;
		}
	    }
	    bigE_seg_free(dgcdp, seg);
	}
    }

    if (ret != ON_SURF) {
	vect_t new_dir;

	VCROSS(new_dir, new_rp.r_dir, rp->r_dir);
	VMOVE(new_rp.r_dir, new_dir);
	/* Compute the inverse of the direction cosines */
	VINVDIR(rd.rd_invdir, new_rp.r_dir);

	if (OBJ[shoot->st_id].ft_shot && OBJ[shoot->st_id].ft_shot(shoot, &new_rp, dgcdp->ap, rd.seghead)) {
	    struct seg *seg;

	    while (BU_LIST_WHILE (seg, seg, &rd.seghead->l)) {
		BU_LIST_DEQUEUE(&seg->l);
		if (ret != ON_SURF) {
		    if (NEAR_ZERO(seg->seg_in.hit_dist, rd.tol->dist)) {
			ret = ON_SURF;
		    }
		    if (NEAR_ZERO(seg->seg_out.hit_dist, rd.tol->dist)) {
			ret = ON_SURF;
		    }
		}
		bigE_seg_free(dgcdp, seg);
	    }
	}
    }
    bigE_ray_seghead_free(dgcdp, rd.seghead);
    return ret;
}


static const char *
bigE_seg_type_name(struct soltab *type)
{
    if (type == ON_SURF)
	return "ON_SURF";
    if (type == ON_INT)
	return "ON_INT";
    if (type == IN_SOL)
	return "IN_SOL";
    return "UNKNOWN";
}


static void
bigE_collect_region_soltabs(union tree *tp, std::vector<struct soltab *> *soltabs)
{
    if (!tp || !soltabs)
	return;

    switch (tp->tr_op) {
	case OP_UNION:
	case OP_SUBTRACT:
	case OP_INTERSECT:
	    bigE_collect_region_soltabs(tp->tr_b.tb_left, soltabs);
	    bigE_collect_region_soltabs(tp->tr_b.tb_right, soltabs);
	    break;
	case OP_SOLID: {
	    struct soltab *stp = tp->tr_a.tu_stp;
	    if (!stp)
		break;
	    for (size_t i = 0; i < soltabs->size(); i++) {
		if ((*soltabs)[i] == stp)
		    return;
	    }
	    soltabs->push_back(stp);
	    break;
	}
	default:
	    break;
    }
}


static void
bigE_build_active_soltab_index(struct bigE_data *dgcdp)
{
    if (!dgcdp)
	return;

    dgcdp->active_soltab_local_index.clear();
    size_t index_size = dgcdp->rtip ? (size_t)dgcdp->rtip->stats.nsolids : 0;
    for (size_t i = 0; i < dgcdp->active_region_soltabs.size(); i++) {
	struct soltab *stp = dgcdp->active_region_soltabs[i];
	if (stp && stp->st_bit >= 0 && (size_t)stp->st_bit + 1 > index_size)
	    index_size = (size_t)stp->st_bit + 1;
    }

    dgcdp->active_soltab_local_index.assign(index_size, -1);
    for (size_t i = 0; i < dgcdp->active_region_soltabs.size(); i++) {
	struct soltab *stp = dgcdp->active_region_soltabs[i];
	if (!stp || stp->st_bit < 0)
	    continue;
	if ((size_t)stp->st_bit >= dgcdp->active_soltab_local_index.size())
	    continue;
	dgcdp->active_soltab_local_index[(size_t)stp->st_bit] = (int)i;
    }
}


static int
bigE_active_soltab_local_index(const struct bigE_data *dgcdp,
			       const struct soltab *stp)
{
    if (!dgcdp || !stp || stp->st_bit < 0 ||
	    (size_t)stp->st_bit >= dgcdp->active_soltab_local_index.size())
	return -1;
    return dgcdp->active_soltab_local_index[(size_t)stp->st_bit];
}


static int
bigE_compile_active_csg_tree(struct bigE_data *dgcdp, union tree *tp)
{
    if (!dgcdp || !tp)
	return -1;

    bigE_csg_eval_node node;
    node.op = tp->tr_op;

    switch (tp->tr_op) {
	case OP_UNION:
	case OP_SUBTRACT:
	case OP_INTERSECT:
	    node.left = bigE_compile_active_csg_tree(dgcdp, tp->tr_b.tb_left);
	    node.right = bigE_compile_active_csg_tree(dgcdp, tp->tr_b.tb_right);
	    break;
	case OP_SOLID:
	    node.local_index =
		bigE_active_soltab_local_index(dgcdp, tp->tr_a.tu_stp);
	    break;
	case OP_NOP:
	    break;
	default:
	    node.op = OP_NOP;
	    break;
    }

    dgcdp->active_csg_nodes.push_back(node);
    return (int)dgcdp->active_csg_nodes.size() - 1;
}


static void
bigE_face_bounds(const bigE_face_record &face,
		 double min_pt[3],
		 double max_pt[3],
		 fastf_t grow)
{
    if (!face.face || !min_pt || !max_pt)
	return;

    for (int i = 0; i < 3; i++) {
	min_pt[i] = face.min_pt[i] - grow;
	max_pt[i] = face.max_pt[i] + grow;
    }
}


struct bigE_face_pair_search_context {
    std::vector<std::pair<int, int>> *pairs;
    int face_index;
    bool query_is_first;

    bigE_face_pair_search_context() :
	pairs(NULL),
	face_index(-1),
	query_is_first(true)
    {}
};


static bool
bigE_face_pair_search_callback(const int &match, void *context)
{
    bigE_face_pair_search_context *ctx =
	static_cast<bigE_face_pair_search_context *>(context);

    if (!ctx || !ctx->pairs)
	return false;

    if (ctx->query_is_first)
	ctx->pairs->push_back(std::make_pair(ctx->face_index, match));
    else
	ctx->pairs->push_back(std::make_pair(match, ctx->face_index));

    return true;
}


static void
bigE_build_leaf_face_index(union E_tree *leaf,
			   bigE_leaf_face_index *index)
{
    if (!leaf || !leaf->l.m || !index)
	return;

    index->faces.clear();
    index->rtree.reset();

    struct nmgregion *r = BU_LIST_FIRST(nmgregion, &leaf->l.m->r_hd);
    if (!r)
	return;
    struct shell *s = BU_LIST_FIRST(shell, &r->s_hd);
    if (!s)
	return;

    size_t order = 0;
    struct faceuse *fu;
    for (BU_LIST_FOR (fu, faceuse, &s->fu_hd)) {
	if (fu->orientation != OT_SAME)
	    continue;

	bigE_face_record rec;
	rec.fu = fu;
	rec.face = fu->f_p;
	rec.order = order++;
	if (!rec.face)
	    continue;
	VMOVE(rec.min_pt, rec.face->min_pt);
	VMOVE(rec.max_pt, rec.face->max_pt);
	bigE_get_fu_plane(rec.plane, fu);
	index->faces.push_back(rec);
    }

    if (index->faces.size() < BIGE_FACE_RTREE_MIN_FACES)
	return;

    index->rtree.reset(new bigE_face_rtree);
    for (size_t i = 0; i < index->faces.size(); i++) {
	double min_pt[3] = {0.0, 0.0, 0.0};
	double max_pt[3] = {0.0, 0.0, 0.0};
	bigE_face_bounds(index->faces[i], min_pt, max_pt, 0.0);
	index->rtree->Insert(min_pt, max_pt, (int)i);
    }
}


static void
bigE_cache_face_edges(bigE_face_record *face)
{
    if (!face || face->edges_cached)
	return;

    face->edges_cached = true;
    if (!face->fu)
	return;

    struct loopuse *lu;
    for (BU_LIST_FOR (lu, loopuse, &face->fu->lu_hd)) {
	if (BU_LIST_FIRST_MAGIC(&lu->down_hd) != NMG_EDGEUSE_MAGIC)
	    continue;
	struct edgeuse *eu;
	for (BU_LIST_FOR (eu, edgeuse, &lu->down_hd)) {
	    struct vertex_g *v1 = eu->vu_p->v_p->vg_p;
	    struct vertex_g *v2 = eu->eumate_p->vu_p->v_p->vg_p;
	    if (!v1 || !v2)
		continue;
	    bigE_face_edge_record edge;
	    VMOVE(edge.start, v1->coord);
	    VSUB2(edge.dir, v2->coord, v1->coord);
	    VMOVE(edge.unit_dir, edge.dir);
	    VUNITIZE(edge.unit_dir);
	    face->edges.push_back(edge);
	}
    }
}


static int
bigE_isect_cached_edge_plane(fastf_t *dist,
			     const bigE_face_edge_record &edge,
			     const fastf_t *plane,
			     const struct bn_tol *tol)
{
    const fastf_t norm_dist = plane[3] - VDOT(plane, edge.start);
    const fastf_t slant_factor = VDOT(plane, edge.dir);
    const fastf_t dot = VDOT(plane, edge.unit_dir);

    if (slant_factor < -SMALL_FASTF && dot < -tol->perp) {
	*dist = norm_dist / slant_factor;
	return 1;
    }
    if (slant_factor > SMALL_FASTF && dot > tol->perp) {
	*dist = norm_dist / slant_factor;
	return 2;
    }

    *dist = 0.0;
    if (norm_dist < -tol->dist)
	return -2;
    if (norm_dist > tol->dist)
	return -1;
    return 0;
}


static void
bigE_collect_face_pairs(const bigE_leaf_face_index *a,
			const bigE_leaf_face_index *b,
			const struct bn_tol *tol,
			std::vector<std::pair<int, int>> *pairs)
{
    if (!pairs)
	return;

    pairs->clear();
    if (!a || !b || !tol || a->faces.empty() || b->faces.empty())
	return;

    bool tree_ordered = false;
    if (a->rtree && b->rtree) {
	a->rtree->Overlaps(*b->rtree, pairs, tol->dist);
	tree_ordered = true;
    } else if (b->rtree) {
	for (size_t i = 0; i < a->faces.size(); i++) {
	    double min_pt[3] = {0.0, 0.0, 0.0};
	    double max_pt[3] = {0.0, 0.0, 0.0};
	    bigE_face_bounds(a->faces[i], min_pt, max_pt, tol->dist);
	    bigE_face_pair_search_context ctx;
	    ctx.pairs = pairs;
	    ctx.face_index = (int)i;
	    ctx.query_is_first = true;
	    b->rtree->Search(min_pt, max_pt, bigE_face_pair_search_callback, &ctx);
	}
	tree_ordered = true;
    } else if (a->rtree) {
	for (size_t i = 0; i < b->faces.size(); i++) {
	    double min_pt[3] = {0.0, 0.0, 0.0};
	    double max_pt[3] = {0.0, 0.0, 0.0};
	    bigE_face_bounds(b->faces[i], min_pt, max_pt, tol->dist);
	    bigE_face_pair_search_context ctx;
	    ctx.pairs = pairs;
	    ctx.face_index = (int)i;
	    ctx.query_is_first = false;
	    a->rtree->Search(min_pt, max_pt, bigE_face_pair_search_callback, &ctx);
	}
	tree_ordered = true;
    } else {
	for (size_t i = 0; i < a->faces.size(); i++) {
	    const bigE_face_record &face1 = a->faces[i];
	    if (!face1.face)
		continue;
	    for (size_t j = 0; j < b->faces.size(); j++) {
		const bigE_face_record &face2 = b->faces[j];
		if (face2.face && V3RPP_OVERLAP_TOL(face1.min_pt,
			face1.max_pt, face2.min_pt,
			face2.max_pt, tol->dist))
		    pairs->push_back(std::make_pair((int)i, (int)j));
	    }
	}
    }

    if (tree_ordered) {
	std::sort(pairs->begin(), pairs->end(),
		[](const std::pair<int, int> &pa,
		    const std::pair<int, int> &pb) {
		    if (pa.first != pb.first)
			return pa.first < pb.first;
		    return pa.second < pb.second;
		});
    }
}


static void
bigE_build_active_source_index(struct bigE_data *dgcdp)
{
    if (!dgcdp)
	return;

    dgcdp->active_source_indices_by_leaf.clear();
    dgcdp->active_candidate_indices_by_leaf.clear();
    dgcdp->active_candidate_indices_by_pair.clear();
    dgcdp->active_source_indices_by_leaf.resize(bigE_leaf_count(dgcdp));
    dgcdp->active_candidate_indices_by_leaf.resize(bigE_leaf_count(dgcdp));
    for (size_t leaf_no = 0; leaf_no < bigE_leaf_count(dgcdp); leaf_no++) {
	union E_tree *leaf = bigE_leaf_at(dgcdp, leaf_no);
	if (!leaf || !leaf->l.stp)
	    continue;

	std::vector<int> &sources =
	    dgcdp->active_source_indices_by_leaf[leaf_no];
	std::vector<int> &candidates =
	    dgcdp->active_candidate_indices_by_leaf[leaf_no];
	for (size_t i = 0; i < dgcdp->active_region_soltabs.size(); i++) {
	    struct soltab *stp = dgcdp->active_region_soltabs[i];
	    if (bigE_soltabs_match_instance(leaf->l.stp, stp, dgcdp->tol))
		sources.push_back((int)i);
	    if (stp && V3RPP_OVERLAP_TOL(stp->st_min, stp->st_max,
		    leaf->l.stp->st_min, leaf->l.stp->st_max,
		    dgcdp->tol->dist))
		candidates.push_back((int)i);
	}
    }
}


static uint64_t
bigE_leaf_pair_key(int leaf1, int leaf2)
{
    uint32_t a = (uint32_t)((leaf1 < leaf2) ? leaf1 : leaf2);
    uint32_t b = (uint32_t)((leaf1 < leaf2) ? leaf2 : leaf1);
    return (((uint64_t)a) << 32) | (uint64_t)b;
}


static const std::vector<int> *
bigE_candidate_indices_for_sources(struct bigE_data *dgcdp,
				   int skip_leaf1,
				   int skip_leaf2)
{
    if (!dgcdp || skip_leaf1 < 0)
	return NULL;

    if (skip_leaf2 < 0) {
	if ((size_t)skip_leaf1 >= dgcdp->active_candidate_indices_by_leaf.size())
	    return NULL;
	return &dgcdp->active_candidate_indices_by_leaf[(size_t)skip_leaf1];
    }

    uint64_t key = bigE_leaf_pair_key(skip_leaf1, skip_leaf2);
    auto it = dgcdp->active_candidate_indices_by_pair.find(key);
    if (it != dgcdp->active_candidate_indices_by_pair.end())
	return &it->second;

    std::vector<int> candidates;
    union E_tree *leaf1 = bigE_leaf_at(dgcdp, (size_t)skip_leaf1);
    union E_tree *leaf2 = bigE_leaf_at(dgcdp, (size_t)skip_leaf2);
    if (leaf1 && leaf2 && leaf1->l.stp && leaf2->l.stp) {
	point_t min_pt;
	point_t max_pt;
	VSET(min_pt,
		std::max(leaf1->l.stp->st_min[X], leaf2->l.stp->st_min[X]),
		std::max(leaf1->l.stp->st_min[Y], leaf2->l.stp->st_min[Y]),
		std::max(leaf1->l.stp->st_min[Z], leaf2->l.stp->st_min[Z]));
	VSET(max_pt,
		std::min(leaf1->l.stp->st_max[X], leaf2->l.stp->st_max[X]),
		std::min(leaf1->l.stp->st_max[Y], leaf2->l.stp->st_max[Y]),
		std::min(leaf1->l.stp->st_max[Z], leaf2->l.stp->st_max[Z]));
	for (size_t i = 0; i < dgcdp->active_region_soltabs.size(); i++) {
	    struct soltab *stp = dgcdp->active_region_soltabs[i];
	    if (stp && V3RPP_OVERLAP_TOL(stp->st_min, stp->st_max,
		    min_pt, max_pt, dgcdp->tol->dist))
		candidates.push_back((int)i);
	}
    }

    auto inserted = dgcdp->active_candidate_indices_by_pair.emplace(key,
	    std::move(candidates));
    return &inserted.first->second;
}


struct bigE_saved_regions {
    struct soltab *stp;
    std::vector<long *> regions;

    bigE_saved_regions() : stp(NULL), regions() {}
};


static void
bigE_restrict_soltabs_to_region(const std::vector<struct soltab *> &soltabs,
				struct region *region,
				std::vector<struct bigE_saved_regions> *saved)
{
    if (!region || !saved)
	return;

    saved->clear();
    saved->reserve(soltabs.size());
    for (size_t i = 0; i < soltabs.size(); i++) {
	struct soltab *stp = soltabs[i];
	if (!stp)
	    continue;

	struct bigE_saved_regions entry;
	entry.stp = stp;
	entry.regions.reserve(BU_PTBL_LEN(&stp->st_regions));
	for (size_t j = 0; j < BU_PTBL_LEN(&stp->st_regions); j++)
	    entry.regions.push_back(BU_PTBL_GET(&stp->st_regions, j));
	saved->push_back(entry);

	bu_ptbl_reset(&stp->st_regions);
	bu_ptbl_ins_unique(&stp->st_regions, (long *)region);
    }
}


static void
bigE_restore_soltab_regions(std::vector<struct bigE_saved_regions> *saved)
{
    if (!saved)
	return;

    for (size_t i = 0; i < saved->size(); i++) {
	struct bigE_saved_regions &entry = (*saved)[i];
	if (!entry.stp)
	    continue;
	bu_ptbl_reset(&entry.stp->st_regions);
	for (size_t j = 0; j < entry.regions.size(); j++)
	    bu_ptbl_ins_unique(&entry.stp->st_regions, entry.regions[j]);
    }
    saved->clear();
}


static struct bigE_segment_summary
bigE_final_surface_summary(struct bu_list *final_segs)
{
    struct bigE_segment_summary summary;
    struct seg *seg;

    if (!final_segs)
	return summary;

    for (BU_LIST_FOR (seg, seg, final_segs)) {
	if (seg->seg_stp != ON_SURF)
	    continue;
	summary.count++;
	if (seg->seg_out.hit_dist > seg->seg_in.hit_dist)
	    summary.length += seg->seg_out.hit_dist - seg->seg_in.hit_dist;
    }

    return summary;
}


static struct bigE_segment_summary
bigE_partition_summary(struct partition *part_head, struct region *region)
{
    struct bigE_segment_summary summary;
    struct partition *pp;

    if (!part_head)
	return summary;

    for (BU_LIST_FOR (pp, partition, (struct bu_list *)part_head)) {
	if (region && pp->pt_regionp != region)
	    continue;
	if (!pp->pt_inhit || !pp->pt_outhit)
	    continue;
	summary.count++;
	if (pp->pt_outhit->hit_dist > pp->pt_inhit->hit_dist)
	    summary.length += pp->pt_outhit->hit_dist - pp->pt_inhit->hit_dist;
    }

    return summary;
}


static fastf_t
bigE_boolweave_length_tol(const struct bigE_data *dgcdp,
			  const struct bigE_segment_summary *bigE_summary,
			  const struct bigE_segment_summary *bw_summary)
{
    fastf_t tol = SMALL_FASTF;
    size_t nsegs = 1;

    if (dgcdp && dgcdp->tol && dgcdp->tol->dist > tol)
	tol = dgcdp->tol->dist;
    if (bigE_summary && bigE_summary->count > nsegs)
	nsegs = bigE_summary->count;
    if (bw_summary && bw_summary->count > nsegs)
	nsegs = bw_summary->count;

    return tol * (fastf_t)nsegs;
}


static void
bigE_partition_list_free(struct partition *head, struct resource *res)
{
    if (!head || !res)
	return;

    struct partition *pp = head->pt_forw;
    while (pp && pp != head) {
	struct partition *zap = pp;
	pp = pp->pt_forw;
	DEQUEUE_PT(zap);
	if (zap->pt_overlap_reg) {
	    bu_free((void *)zap->pt_overlap_reg, "pt_overlap_reg");
	    zap->pt_overlap_reg = NULL;
	}
	BU_LIST_APPEND(&res->re_parthead, (struct bu_list *)zap);
	res->re_partfree++;
    }
    head->pt_forw = head->pt_back = head;
}


struct bigE_point_inside_state {
    struct region *region;
    fastf_t tol;
    fastf_t target_dist;
    int inside;

    bigE_point_inside_state() :
	region(NULL),
	tol(SMALL_FASTF),
	target_dist(0.0),
	inside(0)
    {}
};


static int
bigE_point_inside_hit(struct application *ap,
		      struct partition *PartHeadp,
		      struct seg *UNUSED(segs))
{
    struct bigE_point_inside_state *state =
	(struct bigE_point_inside_state *)ap->a_uptr;
    struct partition *pp;

    if (!state || !PartHeadp)
	return 0;

    for (BU_LIST_FOR (pp, partition, (struct bu_list *)PartHeadp)) {
	if (state->region && pp->pt_regionp != state->region)
	    continue;
	if (!pp->pt_inhit || !pp->pt_outhit)
	    continue;
	if (pp->pt_inhit->hit_dist <= state->target_dist + state->tol &&
		pp->pt_outhit->hit_dist > state->target_dist + state->tol) {
	    state->inside = 1;
	    break;
	}
    }

    return 0;
}


static int
bigE_point_inside_miss(struct application *UNUSED(ap))
{
    return 0;
}


static int
bigE_region_point_inside(struct bigE_data *dgcdp, const point_t pt)
{
    if (!dgcdp || !dgcdp->rtip || !dgcdp->active_region || !pt)
	return 0;

    vect_t model_diag;
    fastf_t tol = (dgcdp->tol && dgcdp->tol->dist > SMALL_FASTF) ?
	dgcdp->tol->dist : SMALL_FASTF;
    VSUB2(model_diag, dgcdp->rtip->mdl_max, dgcdp->rtip->mdl_min);
    fastf_t target_dist = 2.0 * MAGNITUDE(model_diag) + 100.0 * tol;
    if (target_dist <= tol)
	target_dist = 1.0;

    static const vect_t ray_dirs[] = {
	{ 0.533,  0.667, 0.521},
	{-0.371,  0.812, 0.451},
	{ 0.707, -0.203, 0.678}
    };

    for (size_t i = 0; i < sizeof(ray_dirs) / sizeof(ray_dirs[0]); i++) {
	vect_t ray_dir;
	struct bigE_point_inside_state state;
	state.region = dgcdp->active_region;
	state.tol = tol;
	state.target_dist = target_dist;
	VMOVE(ray_dir, ray_dirs[i]);
	VUNITIZE(ray_dir);

	struct application ap;
	RT_APPLICATION_INIT(&ap);
	ap.a_rt_i = dgcdp->rtip;
	ap.a_resource = &dgcdp->resource;
	ap.a_hit = bigE_point_inside_hit;
	ap.a_miss = bigE_point_inside_miss;
	ap.a_logoverlap = rt_silent_logoverlap;
	ap.a_onehit = 0;
	ap.a_uptr = (void *)&state;
	VJOIN1(ap.a_ray.r_pt, pt, -state.target_dist, ray_dir);
	VMOVE(ap.a_ray.r_dir, ray_dir);
	ap.a_ray.r_min = 0.0;
	ap.a_ray.r_max = state.target_dist + 100.0 * state.tol;
	(void)rt_shootray(&ap);
	if (state.inside)
	    return 1;
    }

    return 0;
}


static fastf_t
bigE_boundary_epsilon(struct bigE_data *dgcdp, fastf_t edge_len)
{
    fastf_t tol = (dgcdp && dgcdp->tol && dgcdp->tol->dist > SMALL_FASTF) ?
	dgcdp->tol->dist : SMALL_FASTF;
    fastf_t eps = tol * 10.0;
    fastf_t scaled = edge_len * 1.0e-6;

    if (scaled > eps)
	eps = scaled;
    if (edge_len > tol && eps > edge_len * 0.02)
	eps = edge_len * 0.02;
    if (eps <= SMALL_FASTF)
	eps = SMALL_FASTF;
    return eps;
}


static int
bigE_candidate_interval_is_boundary(struct bigE_data *dgcdp,
				    const point_t start_pt,
				    const vect_t dir,
				    fastf_t edge_len,
				    fastf_t dist)
{
    point_t mid;
    vect_t u, v;
    int seen_inside = 0;
    int seen_outside = 0;
    const fastf_t eps = bigE_boundary_epsilon(dgcdp, edge_len);
    static const fastf_t offsets[8][2] = {
	{ 1.0,  0.0},
	{-1.0,  0.0},
	{ 0.0,  1.0},
	{ 0.0, -1.0},
	{ 1.0,  1.0},
	{ 1.0, -1.0},
	{-1.0,  1.0},
	{-1.0, -1.0}
    };

    if (!dgcdp || !start_pt || !dir)
	return 0;

    VJOIN1(mid, start_pt, dist, dir);
    bn_vec_ortho(u, dir);
    VUNITIZE(u);
    VCROSS(v, dir, u);
    VUNITIZE(v);

    for (size_t i = 0; i < sizeof(offsets) / sizeof(offsets[0]); i++) {
	vect_t off;
	point_t pt;
	VBLEND2(off, offsets[i][0], u, offsets[i][1], v);
	VUNITIZE(off);
	VJOIN1(pt, mid, eps, off);
	if (bigE_region_point_inside(dgcdp, pt))
	    seen_inside = 1;
	else
	    seen_outside = 1;
	if (seen_inside && seen_outside)
	    return 1;
    }

    return 0;
}


static void
bigE_append_candidate_interval(struct bigE_data *dgcdp,
			       const point_t start_pt,
			       const vect_t dir,
			       fastf_t in,
			       fastf_t out)
{
    point_t pt;

    if (!dgcdp || !start_pt || !dir || out <= in)
	return;

    VJOIN1(pt, start_pt, in, dir);
    bigE_line_append(dgcdp, pt, RT_VLIST_LINE_MOVE);
    VJOIN1(pt, start_pt, out, dir);
    bigE_line_append(dgcdp, pt, RT_VLIST_LINE_DRAW);
    dgcdp->nvectors++;
    dgcdp->boolweave_eval_segments++;
    if (dgcdp->profile_enabled)
	dgcdp->profile_current.boolweave_segments++;
}


enum bigE_csg_state {
    BIGE_CSG_OUT = 0,
    BIGE_CSG_IN = 1,
    BIGE_CSG_ON = 2
};


static void
bigE_partition_soltab_bits(const struct bigE_data *dgcdp,
			   const struct partition *pp,
			   std::vector<unsigned char> *bits)
{
    struct seg **segpp;

    if (!dgcdp || !pp || !bits)
	return;

    std::fill(bits->begin(), bits->end(), 0);
    for (BU_PTBL_FOR (segpp, (struct seg **), &pp->pt_seglist)) {
	struct soltab *stp = (*segpp) ? (*segpp)->seg_stp : NULL;
	int local_index = bigE_active_soltab_local_index(dgcdp, stp);
	if (local_index < 0 || (size_t)local_index >= bits->size())
	    continue;
	(*bits)[(size_t)local_index] = 1;
    }
}


static enum bigE_csg_state
bigE_boolweave_tree_state(const struct bigE_data *dgcdp,
			  int node_index,
			  const std::vector<unsigned char> &partition_bits,
			  const std::vector<unsigned char> &source_bits)
{
    if (!dgcdp || node_index < 0 ||
	    (size_t)node_index >= dgcdp->active_csg_nodes.size())
	return BIGE_CSG_OUT;

    const bigE_csg_eval_node &node =
	dgcdp->active_csg_nodes[(size_t)node_index];

    switch (node.op) {
	case OP_UNION: {
	    enum bigE_csg_state a = bigE_boolweave_tree_state(dgcdp,
		    node.left, partition_bits, source_bits);
	    if (a == BIGE_CSG_IN)
		return BIGE_CSG_IN;
	    enum bigE_csg_state b = bigE_boolweave_tree_state(dgcdp,
		    node.right, partition_bits, source_bits);
	    if (b == BIGE_CSG_IN)
		return BIGE_CSG_IN;
	    if (a == BIGE_CSG_ON || b == BIGE_CSG_ON)
		return BIGE_CSG_ON;
	    return BIGE_CSG_OUT;
	}
	case OP_INTERSECT: {
	    enum bigE_csg_state a = bigE_boolweave_tree_state(dgcdp,
		    node.left, partition_bits, source_bits);
	    if (a == BIGE_CSG_OUT)
		return BIGE_CSG_OUT;
	    enum bigE_csg_state b = bigE_boolweave_tree_state(dgcdp,
		    node.right, partition_bits, source_bits);
	    if (b == BIGE_CSG_OUT)
		return BIGE_CSG_OUT;
	    if (a == BIGE_CSG_ON || b == BIGE_CSG_ON)
		return BIGE_CSG_ON;
	    return BIGE_CSG_IN;
	}
	case OP_SUBTRACT: {
	    enum bigE_csg_state a = bigE_boolweave_tree_state(dgcdp,
		    node.left, partition_bits, source_bits);
	    if (a == BIGE_CSG_OUT)
		return BIGE_CSG_OUT;
	    enum bigE_csg_state b = bigE_boolweave_tree_state(dgcdp,
		    node.right, partition_bits, source_bits);
	    if (a == BIGE_CSG_ON)
		return (b == BIGE_CSG_IN) ? BIGE_CSG_OUT : BIGE_CSG_ON;
	    if (b == BIGE_CSG_OUT)
		return BIGE_CSG_IN;
	    if (b == BIGE_CSG_ON)
		return BIGE_CSG_ON;
	    return BIGE_CSG_OUT;
	}
	case OP_SOLID: {
	    if (node.local_index < 0 ||
		    (size_t)node.local_index >= source_bits.size() ||
		    (size_t)node.local_index >= partition_bits.size())
		return BIGE_CSG_OUT;
	    if (source_bits[(size_t)node.local_index])
		return BIGE_CSG_ON;
	    return partition_bits[(size_t)node.local_index] ?
		BIGE_CSG_IN : BIGE_CSG_OUT;
	}
	default:
	    return BIGE_CSG_OUT;
    }
}


static enum bigE_csg_state
bigE_boolweave_tree_state_cached(struct bigE_data *dgcdp,
				 const std::vector<unsigned char> &partition_bits,
				 const std::vector<unsigned char> &source_bits)
{
    if (!dgcdp)
	return BIGE_CSG_OUT;

    if (dgcdp->active_csg_nodes.size() < 64)
	return bigE_boolweave_tree_state(dgcdp, dgcdp->active_csg_root,
		partition_bits, source_bits);

    std::string key;
    key.resize(source_bits.size() + partition_bits.size());
    for (size_t i = 0; i < source_bits.size(); i++)
	key[i] = (char)source_bits[i];
    for (size_t i = 0; i < partition_bits.size(); i++)
	key[source_bits.size() + i] = (char)partition_bits[i];

    auto it = dgcdp->active_csg_eval_cache.find(key);
    if (it != dgcdp->active_csg_eval_cache.end())
	return (enum bigE_csg_state)it->second;

    enum bigE_csg_state state = bigE_boolweave_tree_state(dgcdp,
	    dgcdp->active_csg_root, partition_bits, source_bits);
    if (dgcdp->active_csg_eval_cache.size() < 262144)
	dgcdp->active_csg_eval_cache.emplace(std::move(key),
		(unsigned char)state);
    return state;
}


static void
bigE_boolweave_candidate_plot(struct bigE_data *dgcdp,
			      const point_t start_pt,
			      const vect_t dir,
			      fastf_t edge_len,
			      int skip_leaf1,
			      int skip_leaf2)
{
    if (!dgcdp || !dgcdp->active_region || !dgcdp->active_region->reg_treetop ||
	    !start_pt || !dir || edge_len <= SMALL_FASTF)
	return;

    const int64_t eval_start = bu_gettime();
    dgcdp->boolweave_eval_candidates++;

    const std::vector<struct soltab *> &soltabs =
	dgcdp->active_region_soltabs;
    if (soltabs.empty()) {
	if (dgcdp->profile_enabled)
	    dgcdp->profile_current.candidate_eval_usec +=
		bu_gettime() - eval_start;
	return;
    }

    std::vector<unsigned char> source_bits(soltabs.size(), 0);

    struct application ap;
    RT_APPLICATION_INIT(&ap);
    ap.a_rt_i = dgcdp->rtip;
    ap.a_resource = &dgcdp->resource;
    VMOVE(ap.a_ray.r_pt, start_pt);
    VMOVE(ap.a_ray.r_dir, dir);
    ap.a_ray.r_min = 0.0;
    ap.a_ray.r_max = edge_len;
    ap.a_ray_length = edge_len;
    VINVDIR(ap.a_inv_dir, ap.a_ray.r_dir);

    struct seg waiting_segs;
    struct seg finished_segs;
    struct partition initial_part;
    size_t non_source_segments = 0;
    int64_t profile_stage_start = 0;
    BU_LIST_INIT(&waiting_segs.l);
    BU_LIST_INIT(&finished_segs.l);
    initial_part.pt_forw = initial_part.pt_back = &initial_part;
    initial_part.pt_magic = PT_HD_MAGIC;

    if (dgcdp->profile_enabled)
	profile_stage_start = bu_gettime();
    auto inject_source = [&](int local_index, struct soltab *stp) {
	if (local_index < 0 || (size_t)local_index >= soltabs.size())
	    return;
	if (source_bits[(size_t)local_index])
	    return;
	if (!stp)
	    return;
	struct seg *seg = bigE_seg_get(dgcdp);
	if (dgcdp->profile_enabled)
	    dgcdp->profile_current.boolweave_source_intervals++;
	source_bits[(size_t)local_index] = 1;
	seg->l.magic = RT_SEG_MAGIC;
	seg->seg_in.hit_dist = 0.0;
	seg->seg_out.hit_dist = edge_len;
	seg->seg_in.hit_rayp = &ap.a_ray;
	seg->seg_out.hit_rayp = &ap.a_ray;
	seg->seg_stp = stp;
	BU_LIST_INSERT(&waiting_segs.l, &seg->l);
    };

    auto shoot_soltab = [&](struct soltab *stp) {
	if (!stp || !stp->st_meth || !stp->st_meth->ft_shot)
	    return;

	struct seg new_segs;
	BU_LIST_INIT(&new_segs.l);
	if (dgcdp->profile_enabled)
	    dgcdp->profile_current.boolweave_shot_calls++;
	if (stp->st_meth->ft_shot(stp, &ap.a_ray, &ap, &new_segs) <= 0) {
	    if (dgcdp->profile_enabled)
		dgcdp->profile_current.boolweave_shot_misses++;
	    return;
	}
	if (dgcdp->profile_enabled)
	    dgcdp->profile_current.boolweave_shot_hits++;

	struct seg *seg;
	while (BU_LIST_WHILE (seg, seg, &new_segs.l)) {
	    BU_LIST_DEQUEUE(&seg->l);
	    if (seg->seg_in.hit_dist >= edge_len || seg->seg_out.hit_dist <= 0.0) {
		bigE_seg_free(dgcdp, seg);
		continue;
	    }
	    if (seg->seg_in.hit_dist < 0.0)
		seg->seg_in.hit_dist = 0.0;
	    if (seg->seg_out.hit_dist > edge_len)
		seg->seg_out.hit_dist = edge_len;
	    seg->seg_in.hit_rayp = &ap.a_ray;
	    seg->seg_out.hit_rayp = &ap.a_ray;
	    BU_LIST_INSERT(&waiting_segs.l, &seg->l);
	    non_source_segments++;
	}
    };

    const std::vector<int> *source_candidate_indices =
	bigE_candidate_indices_for_sources(dgcdp, skip_leaf1, skip_leaf2);
    const int use_source_candidates = (soltabs.size() >= 32 &&
	    source_candidate_indices);
    if (use_source_candidates) {
	int skip_ids[2] = {skip_leaf1, skip_leaf2};
	for (size_t si = 0; si < 2; si++) {
	    int skip_leaf = skip_ids[si];
	    if (skip_leaf < 0 ||
		    (size_t)skip_leaf >= dgcdp->active_source_indices_by_leaf.size())
		continue;
	    const std::vector<int> &source_indices =
		dgcdp->active_source_indices_by_leaf[(size_t)skip_leaf];
	    for (size_t i = 0; i < source_indices.size(); i++) {
		int local_index = source_indices[i];
		if (local_index < 0 || (size_t)local_index >= soltabs.size())
		    continue;
		inject_source(local_index, soltabs[(size_t)local_index]);
	    }
	}

	for (size_t i = 0; i < source_candidate_indices->size(); i++) {
	    int local_index = (*source_candidate_indices)[i];
	    if (local_index < 0 || (size_t)local_index >= soltabs.size())
		continue;
	    if (source_bits[(size_t)local_index])
		continue;
	    struct soltab *stp = soltabs[(size_t)local_index];
	    if (!stp)
		continue;
	    if (dgcdp->profile_enabled)
		dgcdp->profile_current.boolweave_bbox_tests++;
	    if (!rt_in_rpp(&ap.a_ray, ap.a_inv_dir, stp->st_min, stp->st_max))
		continue;
	    if (dgcdp->profile_enabled)
		dgcdp->profile_current.boolweave_bbox_hits++;
	    shoot_soltab(soltabs[(size_t)local_index]);
	}
    } else {
	for (size_t i = 0; i < soltabs.size(); i++) {
	    struct soltab *stp = soltabs[i];
	    if (!stp || !stp->st_meth || !stp->st_meth->ft_shot)
		continue;

	    if (bigE_soltab_is_candidate_source(dgcdp, stp, skip_leaf1,
		    skip_leaf2)) {
		inject_source(bigE_active_soltab_local_index(dgcdp, stp), stp);
		continue;
	    }

	    if (dgcdp->profile_enabled)
		dgcdp->profile_current.boolweave_bbox_tests++;
	    if (!rt_in_rpp(&ap.a_ray, ap.a_inv_dir, stp->st_min, stp->st_max))
		continue;
	    if (dgcdp->profile_enabled)
		dgcdp->profile_current.boolweave_bbox_hits++;
	    shoot_soltab(stp);
	}
    }
    if (dgcdp->profile_enabled)
	dgcdp->profile_current.boolweave_shoot_usec +=
	    bu_gettime() - profile_stage_start;

    if (non_source_segments == 0) {
	std::vector<unsigned char> partition_bits(soltabs.size(), 0);
	if (dgcdp->profile_enabled)
	    dgcdp->profile_current.boolweave_source_only++;
	if (bigE_boolweave_tree_state_cached(dgcdp, partition_bits,
		source_bits) == BIGE_CSG_ON)
	    bigE_append_candidate_interval(dgcdp, start_pt, dir, 0.0, edge_len);
	bigE_seg_list_free(&waiting_segs.l, dgcdp);
	int64_t eval_elapsed = bu_gettime() - eval_start;
	dgcdp->boolweave_eval_usec += eval_elapsed;
	if (dgcdp->profile_enabled)
	    dgcdp->profile_current.candidate_eval_usec += eval_elapsed;
	return;
    }

    if (dgcdp->profile_enabled)
	profile_stage_start = bu_gettime();
    if (BU_LIST_NON_EMPTY(&waiting_segs.l))
	rt_boolweave(&finished_segs, &waiting_segs, &initial_part, &ap);
    if (dgcdp->profile_enabled)
	dgcdp->profile_current.boolweave_weave_usec +=
	    bu_gettime() - profile_stage_start;

    const fastf_t tol = (dgcdp->tol && dgcdp->tol->dist > SMALL_FASTF) ?
	dgcdp->tol->dist : SMALL_FASTF;
    std::vector<unsigned char> partition_bits(soltabs.size(), 0);
    struct partition *pp;
    if (dgcdp->profile_enabled)
	profile_stage_start = bu_gettime();
    for (BU_LIST_FOR (pp, partition, (struct bu_list *)&initial_part)) {
	if (dgcdp->profile_enabled)
	    dgcdp->profile_current.boolweave_partitions++;
	if (!pp->pt_inhit || !pp->pt_outhit)
	    continue;
	fastf_t in = pp->pt_inhit->hit_dist;
	fastf_t out = pp->pt_outhit->hit_dist;
	if (in < 0.0)
	    in = 0.0;
	if (out > edge_len)
	    out = edge_len;
	if (out - in <= tol)
	    continue;
	bigE_partition_soltab_bits(dgcdp, pp, &partition_bits);
	if (bigE_boolweave_tree_state_cached(dgcdp, partition_bits,
		source_bits) != BIGE_CSG_ON)
	    continue;
	if (dgcdp->boolweave_eval_filter_enabled) {
	    fastf_t mid = (in + out) * 0.5;
	    if (!bigE_candidate_interval_is_boundary(dgcdp, start_pt, dir,
		    edge_len, mid))
		continue;
	}
	bigE_append_candidate_interval(dgcdp, start_pt, dir, in, out);
    }
    if (dgcdp->profile_enabled)
	dgcdp->profile_current.boolweave_partition_usec +=
	    bu_gettime() - profile_stage_start;

    struct resource *res = &dgcdp->resource;
    bigE_partition_list_free(&initial_part, res);
    bigE_seg_list_free(&finished_segs.l, dgcdp);
    bigE_seg_list_free(&waiting_segs.l, dgcdp);

    int64_t eval_elapsed = bu_gettime() - eval_start;
    dgcdp->boolweave_eval_usec += eval_elapsed;
    if (dgcdp->profile_enabled)
	dgcdp->profile_current.candidate_eval_usec += eval_elapsed;
}


static void
bigE_boolweave_candidate_diag(struct bigE_data *dgcdp,
			      size_t candidate_id,
			      const point_t start_pt,
			      const vect_t dir,
			      fastf_t edge_len,
			      struct soltab *candidate_type,
			      int skip_leaf1,
			      int skip_leaf2,
			      const struct bigE_segment_summary *bigE_summary)
{
    if (!dgcdp || !dgcdp->boolweave_diag_enabled ||
	    !dgcdp->active_region || !dgcdp->active_region->reg_treetop ||
	    !bigE_summary)
	return;

    const int64_t diag_start = bu_gettime();
    dgcdp->boolweave_diag_candidates++;

    const std::vector<struct soltab *> &soltabs =
	dgcdp->active_region_soltabs;
    if (soltabs.empty()) {
	dgcdp->boolweave_diag_skipped++;
	return;
    }

    struct application ap;
    RT_APPLICATION_INIT(&ap);
    ap.a_rt_i = dgcdp->rtip;
    ap.a_resource = &dgcdp->resource;
    VMOVE(ap.a_ray.r_pt, start_pt);
    VMOVE(ap.a_ray.r_dir, dir);
    ap.a_ray.r_min = 0.0;
    ap.a_ray.r_max = edge_len;
    ap.a_ray_length = edge_len;
    VINVDIR(ap.a_inv_dir, ap.a_ray.r_dir);

    struct seg waiting_segs;
    struct seg finished_segs;
    struct partition initial_part;
    struct partition final_part;
    struct bu_ptbl regionbits;
    struct bu_bitv *solidbits = NULL;
    int shot_hits = 0;
    int shot_misses = 0;
    int source_injected = 0;

    BU_LIST_INIT(&waiting_segs.l);
    BU_LIST_INIT(&finished_segs.l);
    initial_part.pt_forw = initial_part.pt_back = &initial_part;
    initial_part.pt_magic = PT_HD_MAGIC;
    final_part.pt_forw = final_part.pt_back = &final_part;
    final_part.pt_magic = PT_HD_MAGIC;
    bu_ptbl_init(&regionbits, 7, "bigE boolweave diagnostic regionbits");
    solidbits = rt_get_solidbitv(dgcdp->rtip->stats.nsolids, &dgcdp->resource);

    std::vector<struct bigE_saved_regions> saved_regions;
    bigE_restrict_soltabs_to_region(soltabs, dgcdp->active_region, &saved_regions);

    const int64_t shot_start = bu_gettime();
    for (size_t i = 0; i < soltabs.size(); i++) {
	struct soltab *stp = soltabs[i];
	if (!stp || !stp->st_meth || !stp->st_meth->ft_shot)
	    continue;
	BU_BITSET(solidbits, stp->st_bit);
	if (bigE_soltab_is_candidate_source(dgcdp, stp, skip_leaf1, skip_leaf2)) {
	    /* Candidate generation already proved this source surface.  Inject
	     * that interval directly instead of asking ft_shot to graze it. */
	    struct seg *seg = bigE_seg_get(dgcdp);
	    seg->l.magic = RT_SEG_MAGIC;
	    seg->seg_in.hit_dist = 0.0;
	    seg->seg_out.hit_dist = edge_len;
	    seg->seg_in.hit_rayp = &ap.a_ray;
	    seg->seg_out.hit_rayp = &ap.a_ray;
	    seg->seg_stp = stp;
	    BU_LIST_INSERT(&waiting_segs.l, &seg->l);
	    source_injected++;
	    continue;
	}
	if (!rt_in_rpp(&ap.a_ray, ap.a_inv_dir, stp->st_min, stp->st_max)) {
	    shot_misses++;
	    continue;
	}

	struct seg new_segs;
	BU_LIST_INIT(&new_segs.l);
	if (stp->st_meth->ft_shot(stp, &ap.a_ray, &ap, &new_segs) <= 0) {
	    shot_misses++;
	    continue;
	}
	shot_hits++;

	struct seg *seg;
	while (BU_LIST_WHILE (seg, seg, &new_segs.l)) {
	    BU_LIST_DEQUEUE(&seg->l);
	    if (seg->seg_in.hit_dist >= edge_len || seg->seg_out.hit_dist <= 0.0) {
		bigE_seg_free(dgcdp, seg);
		continue;
	    }
	    if (seg->seg_in.hit_dist < 0.0)
		seg->seg_in.hit_dist = 0.0;
	    if (seg->seg_out.hit_dist > edge_len)
		seg->seg_out.hit_dist = edge_len;
	    seg->seg_in.hit_rayp = &ap.a_ray;
	    seg->seg_out.hit_rayp = &ap.a_ray;
	    BU_LIST_INSERT(&waiting_segs.l, &seg->l);
	}
    }
    dgcdp->boolweave_diag_bw_shot_usec += bu_gettime() - shot_start;

    const int64_t bool_start = bu_gettime();
    if (BU_LIST_NON_EMPTY(&waiting_segs.l))
	rt_boolweave(&finished_segs, &waiting_segs, &initial_part, &ap);

    if (BU_LIST_NON_EMPTY(&finished_segs.l))
	(void)rt_boolfinal(&initial_part, &final_part, 0.0, edge_len,
		&regionbits, &ap, solidbits);
    dgcdp->boolweave_diag_bw_bool_usec += bu_gettime() - bool_start;

    struct bigE_segment_summary bw_summary =
	bigE_partition_summary(&final_part, dgcdp->active_region);

    const fastf_t len_delta = bw_summary.length - bigE_summary->length;
    const fastf_t abs_len_delta = fabs(len_delta);
    const fastf_t len_tol =
	bigE_boolweave_length_tol(dgcdp, bigE_summary, &bw_summary);
    const int count_match = (bigE_summary->count == bw_summary.count);
    const int length_match = (abs_len_delta <= len_tol);

    if (count_match && length_match)
	dgcdp->boolweave_diag_matches++;
    if (!count_match)
	dgcdp->boolweave_diag_count_mismatches++;
    if (!length_match)
	dgcdp->boolweave_diag_length_mismatches++;
    dgcdp->boolweave_diag_total_abs_len_delta += abs_len_delta;
    if (abs_len_delta > dgcdp->boolweave_diag_max_abs_len_delta)
	dgcdp->boolweave_diag_max_abs_len_delta = abs_len_delta;

    if (dgcdp->boolweave_diag_logged < dgcdp->boolweave_diag_limit) {
	bu_log("evaluated-wire boolweave diag: status=%s region=\"%s\" "
		"candidate=%zu type=%s len=%.12g skip=(%d,%d) leaves=%zu "
		"source=%d shot=(%d hit,%d miss) bigE=(%zu seg, %.12g len) "
		"boolweave=(%zu part, %.12g len) len_delta=%.12g len_tol=%.12g\n",
		(count_match && length_match) ? "match" : "mismatch",
		dgcdp->active_region->reg_name ? dgcdp->active_region->reg_name : "(null)",
		candidate_id, bigE_seg_type_name(candidate_type), edge_len,
		skip_leaf1, skip_leaf2, soltabs.size(), source_injected,
		shot_hits, shot_misses,
		bigE_summary->count, bigE_summary->length,
		bw_summary.count, bw_summary.length, len_delta, len_tol);
	dgcdp->boolweave_diag_logged++;
    }

    struct resource *res = &dgcdp->resource;
    bigE_restore_soltab_regions(&saved_regions);
    bigE_partition_list_free(&initial_part, res);
    bigE_partition_list_free(&final_part, res);
    bigE_seg_list_free(&finished_segs.l, dgcdp);
    bigE_seg_list_free(&waiting_segs.l, dgcdp);
    if (solidbits)
	BU_LIST_APPEND(&res->re_solid_bitv, &solidbits->l);
    bu_ptbl_free(&regionbits);

    dgcdp->boolweave_diag_bw_usec += bu_gettime() - diag_start;
}


/* Shoot rays (corresponding to possible edges in the result) at the
 * solids, put the results in the E-tree leaves as type IN_SOL.  Call
 * eval_etree() and plot the results
 */
static void
shoot_and_plot(point_t start_pt,
	       vect_t dir,
	       fastf_t edge_len,
	       int skip_leaf1,
	       int skip_leaf2,
	       union E_tree *eptr,
	       struct soltab *type,
	       struct bigE_data *dgcdp)
{
    struct xray rp;
    struct ray_data rd;
    size_t shoot_leaf;
    struct bu_list *final_segs;
    size_t candidate_id;
    int64_t bige_start = 0;

    CK_ETREE(eptr);
    if (dgcdp->profile_enabled) {
	if (type == ON_INT)
	    dgcdp->profile_current.intersection_candidates++;
	else
	    dgcdp->profile_current.edge_candidates++;
    }
    if (edge_len <= dgcdp->tol->dist) {
	if (dgcdp->profile_enabled)
	    dgcdp->profile_current.degenerate_candidates++;
	return;
    }
    if (dgcdp->profile_enabled)
	dgcdp->profile_current.candidates++;

    candidate_id = dgcdp->candidate_index++;
    if (dgcdp->boolweave_eval_enabled) {
	bigE_boolweave_candidate_plot(dgcdp, start_pt, dir, edge_len,
		skip_leaf1, skip_leaf2);
	return;
    }
    if (dgcdp->boolweave_diag_enabled || dgcdp->bige_stats_enabled ||
	    dgcdp->profile_enabled)
	bige_start = bu_gettime();

    memset(&rd, 0, sizeof(struct ray_data));

    rd.seghead = bigE_ray_seghead_create();

    VMOVE(rp.r_pt, start_pt);
    VMOVE(rp.r_dir, dir);
    /* Compute the inverse of the direction cosines */
    VINVDIR(rd.rd_invdir, rp.r_dir);

    /* set up "ray_data" structure for nmg raytrace */
    rd.rp = &rp;
    rd.tol = dgcdp->tol;
    rd.ap = dgcdp->ap;
    rd.magic = NMG_RAY_DATA_MAGIC;
    rd.classifying_ray = 0;
    rd.hitmiss = (struct hitmiss **)NULL;

    /* shoot this ray at every leaf solid except the one this edge
     * came from (or the two that this intersection line came from
     */
    for (shoot_leaf=0; shoot_leaf < bigE_leaf_count(dgcdp); shoot_leaf++) {
	union E_tree *shoot;

	shoot = bigE_leaf_at(dgcdp, shoot_leaf);
	if (!shoot)
	    continue;

	if (BU_LIST_NON_EMPTY(&shoot->l.seghead)) {
	    bigE_seg_list_free(&shoot->l.seghead, dgcdp);
	}
	BU_LIST_INIT(&shoot->l.seghead);

	if (bigE_should_skip_leaf(dgcdp, shoot_leaf, skip_leaf1, skip_leaf2)) {
	    struct seg *seg;

	    /* put entire edge in seg list and mark it as ON the
	     * surface.
	     */
	    seg = bigE_seg_get(dgcdp);
	    seg->l.magic = RT_SEG_MAGIC;
	    seg->seg_in.hit_dist = 0.0;
	    seg->seg_out.hit_dist = edge_len;
	    seg->seg_stp = type;
	    BU_LIST_INSERT(&shoot->l.seghead, &seg->l);
	    continue;
	}

	/* initialize the lists of things that have been hit/missed */
	rd.rd_m = shoot->l.m;
	BU_LIST_INIT(&rd.rd_hit);
	BU_LIST_INIT(&rd.rd_miss);

	rd.stp = shoot->l.stp;

	/* actually shoot the ray, assign segments to the leaf, and
	 * mark them as IN_SOL.
	 */
	if (rt_in_rpp(&rp, rd.rd_invdir, shoot->l.stp->st_min, shoot->l.stp->st_max)) {
	    if (OBJ[shoot->l.stp->st_id].ft_shot && OBJ[shoot->l.stp->st_id].ft_shot(shoot->l.stp, &rp, dgcdp->ap, rd.seghead)) {
		struct seg *seg;

		/* put the segments in the lead solid structure */
		while (BU_LIST_WHILE (seg, seg, &rd.seghead->l)) {
		    BU_LIST_DEQUEUE(&seg->l);
		    /* clip segments to the edge being considered */
		    if (seg->seg_in.hit_dist >= edge_len || seg->seg_out.hit_dist <= 0) {
			bigE_seg_free(dgcdp, seg);
		    } else {
			if (seg->seg_in.hit_dist < 0.0)
			    seg->seg_in.hit_dist = 0.0;
			if (seg->seg_out.hit_dist > edge_len)
			    seg->seg_out.hit_dist = edge_len;
			seg->seg_stp = classify_seg(seg, shoot->l.stp, &rp, dgcdp);
			BU_LIST_INSERT(&shoot->l.seghead, &seg->l);
		    }
		}
	    }
	}
    }

    /* Evaluate the Boolean tree to get the "final" segments which are
     * to be plotted.
     */
    final_segs = eval_etree(eptr, dgcdp);
    struct bigE_segment_summary final_summary =
	bigE_final_surface_summary(final_segs);
    if (final_segs) {
	struct seg *seg;

	/* add the segments to the typed line accumulator */
	for (BU_LIST_FOR (seg, seg, final_segs)) {
	    point_t pt;

	    /* only plot the resulting segments that are ON the
	     * SURFace.
	     */
	    if (seg->seg_stp != ON_SURF)
		continue;

	    dgcdp->nvectors++;
	    VJOIN1(pt, rp.r_pt, seg->seg_in.hit_dist, rp.r_dir);
	    bigE_line_append(dgcdp, pt, RT_VLIST_LINE_MOVE);
	    VJOIN1(pt, rp.r_pt, seg->seg_out.hit_dist, rp.r_dir);
	    bigE_line_append(dgcdp, pt, RT_VLIST_LINE_DRAW);
	}

    }

    if (dgcdp->boolweave_diag_enabled || dgcdp->bige_stats_enabled ||
	    dgcdp->profile_enabled) {
	int64_t bige_elapsed = bu_gettime() - bige_start;
	if (dgcdp->boolweave_diag_enabled)
	    dgcdp->boolweave_diag_bige_usec += bige_elapsed;
	if (dgcdp->bige_stats_enabled) {
	    dgcdp->bige_stats_candidates++;
	    dgcdp->bige_stats_segments += final_summary.count;
	    dgcdp->bige_stats_usec += bige_elapsed;
	}
	if (dgcdp->profile_enabled) {
	    dgcdp->profile_current.bige_segments += final_summary.count;
	    dgcdp->profile_current.candidate_eval_usec += bige_elapsed;
	}
    }

    bigE_boolweave_candidate_diag(dgcdp, candidate_id, start_pt, dir,
	    edge_len, type, skip_leaf1, skip_leaf2, &final_summary);

    if (final_segs)
	bigE_seg_list_free(final_segs, dgcdp);
    bu_free((char *)final_segs, "bu_list");
    bigE_ray_seghead_free(dgcdp, rd.seghead);
}


#define HITS_BLOCK 20

static void
Eplot(union E_tree *eptr,
      struct bigE_data *dgcdp)
{
    point_t start_pt;
    size_t leaf_no;
    union E_tree *leaf_ptr;
    int hit_count1=0, hit_count2=0;
    point_t *hits1=NULL, *hits2=NULL;
    int hits_avail1=0, hits_avail2=0;
    int i;
    struct bu_list *result;
    const struct bn_tol *tol = dgcdp->tol;
    struct bu_list *vlfree = dgcdp->vlfree;
    int64_t stage_start = 0;
    std::vector<bigE_leaf_face_index> face_indices;
    std::vector<std::pair<int, int>> face_pairs;

    CK_ETREE(eptr);

    /* create an edge list for each leaf solid */
    if (dgcdp->profile_enabled)
	stage_start = bu_gettime();
    for (leaf_no=0; leaf_no < bigE_leaf_count(dgcdp); leaf_no++) {
	leaf_ptr = bigE_leaf_at(dgcdp, leaf_no);
	if (!leaf_ptr)
	    continue;
	CK_ETREE(leaf_ptr);
	if (leaf_ptr->l.op != OP_DB_LEAF && leaf_ptr->l.op != OP_SOLID) {
	    continue;
	}

	if (leaf_ptr->l.m) {
	    nmg_edge_tabulate(&leaf_ptr->l.edge_list, &leaf_ptr->l.m->magic, vlfree);
	    leaf_ptr->l.edge_list_initialized = 1;
	} else {
	    bu_ptbl_init(&leaf_ptr->l.edge_list, 1, "edge_list");
	    leaf_ptr->l.edge_list_initialized = 1;
	}
    }
    if (dgcdp->profile_enabled)
	dgcdp->profile_current.edge_tab_usec += bu_gettime() - stage_start;

    /* now plot appropriate parts of each solid */

    /* loop through every leaf solid */
    if (dgcdp->profile_enabled)
	stage_start = bu_gettime();
    for (leaf_no=0; leaf_no < bigE_leaf_count(dgcdp); leaf_no++) {
	size_t edge_no;

	leaf_ptr = bigE_leaf_at(dgcdp, leaf_no);
	if (!leaf_ptr)
	    continue;

	if (!leaf_ptr->l.m)
	    continue;

	/* do each edge of the current leaf solid */
	for (edge_no=0; edge_no < BU_PTBL_LEN(&leaf_ptr->l.edge_list); edge_no++) {
	    struct edge *e;
	    struct vertex_g *vg;
	    struct vertex_g *vg2;
	    vect_t dir;
	    fastf_t edge_len;
	    fastf_t inv_len;

	    e = (struct edge *)BU_PTBL_GET(&leaf_ptr->l.edge_list, edge_no);
	    NMG_CK_EDGE(e);
	    vg = e->eu_p->vu_p->v_p->vg_p;
	    NMG_CK_VERTEX_G(vg);

	    vg2 = e->eu_p->eumate_p->vu_p->v_p->vg_p;
	    NMG_CK_VERTEX_G(vg2);

	    /* set up a ray from vg towards vg2 */
	    VSUB2(dir, vg2->coord, vg->coord);
	    edge_len = MAGNITUDE(dir);
	    if (edge_len < tol->dist)
		continue;
	    inv_len = 1.0/edge_len;
	    VSCALE(dir, dir, inv_len);
	    shoot_and_plot(vg->coord, dir, edge_len, (int)leaf_no, -1, eptr, ON_SURF, dgcdp);

	}
    }
    if (dgcdp->profile_enabled)
	dgcdp->profile_current.edge_candidate_usec += bu_gettime() - stage_start;

    hits1 = (point_t *)bu_calloc(HITS_BLOCK, sizeof(point_t), "hits");
    hits_avail1 = HITS_BLOCK;
    hits2 = (point_t *)bu_calloc(HITS_BLOCK, sizeof(point_t), "hits");
    hits_avail2 = HITS_BLOCK;

    /* Now draw solid intersection lines */
    if (dgcdp->profile_enabled)
	stage_start = bu_gettime();
    int64_t face_index_start = dgcdp->profile_enabled ? bu_gettime() : 0;
    face_indices.resize(bigE_leaf_count(dgcdp));
    for (leaf_no=0; leaf_no < bigE_leaf_count(dgcdp); leaf_no++) {
	leaf_ptr = bigE_leaf_at(dgcdp, leaf_no);
	if (leaf_ptr && leaf_ptr->l.m)
	    bigE_build_leaf_face_index(leaf_ptr, &face_indices[leaf_no]);
    }
    if (dgcdp->profile_enabled)
	dgcdp->profile_current.intersection_face_index_usec +=
	    bu_gettime() - face_index_start;
    for (leaf_no=0; leaf_no < bigE_leaf_count(dgcdp); leaf_no++) {
	size_t leaf2;

	leaf_ptr = bigE_leaf_at(dgcdp, leaf_no);
	if (!leaf_ptr)
	    continue;
	if (!leaf_ptr->l.m)
	    continue;
	if (face_indices[leaf_no].faces.empty())
	    continue;

	for (leaf2=leaf_no+1; leaf2 < bigE_leaf_count(dgcdp); leaf2++) {
	    union E_tree *leaf2_ptr;
	    struct faceuse *fu1, *fu2;
	    struct face *f1, *f2;
	    plane_t pl1, pl2;
	    struct bu_list *A, *B;

	    leaf2_ptr = bigE_leaf_at(dgcdp, leaf2);
	    if (!leaf2_ptr)
		continue;
	    if (!leaf2_ptr->l.m)
		continue;
	    if (face_indices[leaf2].faces.empty())
		continue;
	    if (!leaf_ptr->l.stp || !leaf2_ptr->l.stp)
		continue;
	    if (!V3RPP_OVERLAP_TOL(leaf2_ptr->l.stp->st_min,
		    leaf2_ptr->l.stp->st_max,
		    leaf_ptr->l.stp->st_min,
		    leaf_ptr->l.stp->st_max,
		    tol->dist))
		continue;

	    /* find intersection lines between these two NMG's */
	    int64_t pair_collect_start =
		dgcdp->profile_enabled ? bu_gettime() : 0;
	    bigE_collect_face_pairs(&face_indices[leaf_no],
		    &face_indices[leaf2], tol, &face_pairs);
	    int64_t pair_collect_usec = 0;
	    if (dgcdp->profile_enabled) {
		pair_collect_usec = bu_gettime() - pair_collect_start;
		dgcdp->profile_current.intersection_pair_collect_usec +=
		    pair_collect_usec;
		dgcdp->profile_current.intersection_face_pairs +=
		    face_pairs.size();
	    }
	    int64_t pair_eval_start = dgcdp->profile_enabled ? bu_gettime() : 0;
	    size_t pair_eval_start_candidates =
		dgcdp->profile_current.intersection_candidates;
	    size_t pair_eval_start_bige_segments =
		dgcdp->profile_current.bige_segments;
	    size_t pair_eval_start_boolweave_segments =
		dgcdp->profile_current.boolweave_segments;
	    if (dgcdp->profile_progress_enabled && !face_pairs.empty()) {
		bu_log("evaluated-wire profile leaf-pair-start: region=\"%s\" "
			"leaf=%zu leaf2=%zu faces=(%zu,%zu) pairs=%zu "
			"collect=%.6fs\n",
			dgcdp->profile_current.name.empty() ? "(null)" :
			dgcdp->profile_current.name.c_str(),
			leaf_no,
			leaf2,
			face_indices[leaf_no].faces.size(),
			face_indices[leaf2].faces.size(),
			face_pairs.size(),
			((double)pair_collect_usec) / 1000000.0);
	    }
	    for (size_t face_pair_no = 0;
		    face_pair_no < face_pairs.size();
		    face_pair_no++) {
		fastf_t dist;
		vect_t dir;
		vect_t vdiff;
		fastf_t *dists1, *dists2;
		fastf_t min_dist, max_dist;
		int min_hit, max_hit;
		int done;
		struct seg *aseg;
		int face1_index = face_pairs[face_pair_no].first;
		int face2_index = face_pairs[face_pair_no].second;

		if (face1_index < 0 ||
			(size_t)face1_index >= face_indices[leaf_no].faces.size())
		    continue;
		if (face2_index < 0 ||
			(size_t)face2_index >= face_indices[leaf2].faces.size())
		    continue;

		bigE_face_record &face1 =
		    face_indices[leaf_no].faces[(size_t)face1_index];
		fu1 = face1.fu;
		f1 = face1.face;
		if (!fu1 || !f1)
		    continue;
		HMOVE(pl1, face1.plane);

		bigE_face_record &face2 =
		    face_indices[leaf2].faces[(size_t)face2_index];
		fu2 = face2.fu;
		f2 = face2.face;
		if (!fu2 || !f2)
		    continue;

		if (!V3RPP_OVERLAP_TOL(face2.min_pt, face2.max_pt,
			face1.min_pt, face1.max_pt, tol->dist))
		    continue;

		HMOVE(pl2, face2.plane);

		    if (bg_coplanar(pl1, pl2, tol)) {
			continue;
		    }
		    bigE_cache_face_edges(&face1);
		    bigE_cache_face_edges(&face2);

		    hit_count1=0;
		    hit_count2=0;
		    for (const bigE_face_edge_record &edge : face1.edges) {
			if (bigE_isect_cached_edge_plane(&dist, edge, pl2, tol) < 1)
			    continue;
			if (dist < -tol->dist || dist > 1.0 + tol->dist)
			    continue;

			if (hit_count1 >= hits_avail1) {
			    hits_avail1 += HITS_BLOCK;
			    hits1 = (point_t *)bu_realloc(hits1,
				    hits_avail1 * sizeof(point_t), "hits1");
			}
			VJOIN1(hits1[hit_count1], edge.start, dist, edge.dir);
			hit_count1++;
		    }
		    for (const bigE_face_edge_record &edge : face2.edges) {
			if (bigE_isect_cached_edge_plane(&dist, edge, pl1, tol) < 1)
			    continue;
			if (dist < -tol->dist || dist > 1.0 + tol->dist)
			    continue;

			if (hit_count2 >= hits_avail2) {
			    hits_avail2 += HITS_BLOCK;
			    hits2 = (point_t *)bu_realloc(hits2,
				    hits_avail2 * sizeof(point_t), "hits2");
			}
			VJOIN1(hits2[hit_count2], edge.start, dist, edge.dir);
			hit_count2++;
		    }

		    if (hit_count1 < 2 || hit_count2 < 2) {
			/* nothing to plot */
			continue;
		    }

		    /* sort the hits on face 1 */
		    dists1 = (fastf_t *)bu_calloc(hit_count1,
						  sizeof(fastf_t), "dists1");
		    dists2 = (fastf_t *)bu_calloc(hit_count2,
						  sizeof(fastf_t), "dists2");
		    VMOVE(start_pt, hits1[0]);
		    dists1[0] = 0.0;
		    min_dist = 0.0;
		    min_hit = 0;
		    VSUB2(dir, hits1[1], hits1[0]);
		    dists1[1] = MAGNITUDE(dir);
		    VUNITIZE(dir);
		    max_dist = dists1[1];
		    max_hit = 1;
		    for (i = 2; i < hit_count1; i++) {
			VSUB2(vdiff, hits1[i], start_pt);
			dists1[i] = MAGNITUDE(vdiff);
			if (VDOT(dir, vdiff) < 0.0)
			    dists1[i] = -dists1[i];
			if (dists1[i] > max_dist) {
			    max_dist = dists1[i];
			    max_hit = i;
			}
			if (dists1[i] < min_dist) {
			    min_dist = dists1[i];
			    min_hit = i;
			}
		    }

		    /* recalculate dir */
		    VSUB2(dir, hits1[max_hit], hits1[min_hit]);
		    VUNITIZE(dir);

		    done = 0;
		    while (!done) {
			done = 1;
			for (i = 1; i < hit_count1; i++) {
			    if (dists1[i-1] > dists1[i]) {
				fastf_t tmp;
				point_t tmp_pt;

				done = 0;
				tmp = dists1[i];
				VMOVE(tmp_pt, hits1[i]);
				dists1[i] = dists1[i-1];
				VMOVE(hits1[i], hits1[i-1]);
				dists1[i-1] = tmp;
				VMOVE(hits1[i-1], tmp_pt);
			    }
			}
		    }

		    /* sort the hits on face 2 */
		    min_dist = MAX_FASTF;
		    max_dist = -min_dist;
		    for (i = 0; i < hit_count2; i++) {
			VSUB2(vdiff, hits2[i], start_pt);
			dists2[i] = MAGNITUDE(vdiff);
			if (VDOT(dir, vdiff) < 0.0)
			    dists2[i] = -dists2[i];
			if (dists2[i] > max_dist) {
			    max_dist = dists2[i];
			}
			if (dists2[i] < min_dist) {
			    min_dist = dists2[i];
			}
		    }

		    done = 0;
		    while (!done) {
			done = 1;
			for (i = 1; i < hit_count2; i++) {
			    if (dists2[i-1] > dists2[i]) {
				fastf_t tmp;
				point_t tmp_pt;

				done = 0;
				tmp = dists2[i];
				VMOVE(tmp_pt, hits2[i]);
				dists2[i] = dists2[i-1];
				VMOVE(hits2[i], hits2[i-1]);
				dists2[i-1] = tmp;
				VMOVE(hits2[i-1], tmp_pt);
			    }
			}
		    }

		    /* build a segment list for each solid */
		    BU_ALLOC(A, struct bu_list);
		    BU_ALLOC(B, struct bu_list);
		    BU_LIST_INIT(A);
		    BU_LIST_INIT(B);

		    for (i = 1; i < hit_count1; i += 2) {
			fastf_t diff;

			diff = dists1[i] - dists1[i-1];
			if (NEAR_ZERO(diff, tol->dist)) {
			    continue;
			}
			aseg = bigE_seg_get(dgcdp);
			aseg->l.magic = RT_SEG_MAGIC;
			aseg->seg_stp = ON_INT;
			VMOVE(aseg->seg_in.hit_point, hits1[i-1]);
			aseg->seg_in.hit_dist = dists1[i-1];
			VMOVE(aseg->seg_out.hit_point, hits1[i]);
			aseg->seg_out.hit_dist = dists1[i];

			BU_LIST_APPEND(A, &aseg->l);
		    }

		    for (i = 1; i < hit_count2; i += 2) {
			fastf_t diff;

			diff = dists2[i] - dists2[i-1];
			if (NEAR_ZERO(diff, tol->dist)) {
			    continue;
			}
			aseg = bigE_seg_get(dgcdp);
			aseg->l.magic = RT_SEG_MAGIC;
			aseg->seg_stp = ON_INT;
			VMOVE(aseg->seg_in.hit_point, hits2[i-1]);
			aseg->seg_in.hit_dist = dists2[i-1];
			VMOVE(aseg->seg_out.hit_point, hits2[i]);
			aseg->seg_out.hit_dist = dists2[i];

			BU_LIST_APPEND(B, &aseg->l);
		    }

		    result = eval_op(A, OP_INTERSECT, B, dgcdp);

		    for (BU_LIST_FOR (aseg, seg, result)) {
			point_t ray_start;
			fastf_t seg_len = aseg->seg_out.hit_dist -
			    aseg->seg_in.hit_dist;

			if (seg_len <= tol->dist)
			    continue;
			VJOIN1(ray_start, start_pt, aseg->seg_in.hit_dist, dir);
			shoot_and_plot(ray_start, dir, seg_len,
				       (int)leaf_no, (int)leaf2, eptr, ON_INT, dgcdp);
		    }
		    bigE_seg_list_free(result, dgcdp);
		    bu_free((char *)result, "bu_list");

		    bu_free((char *)dists1, "dists1");
		    bu_free((char *)dists2, "dists2");
		}
	    if (dgcdp->profile_enabled) {
		int64_t pair_eval_usec = bu_gettime() - pair_eval_start;
		dgcdp->profile_current.intersection_pair_eval_usec +=
		    pair_eval_usec;
		if (dgcdp->profile_progress_enabled && !face_pairs.empty()) {
		    bu_log("evaluated-wire profile leaf-pair-end: region=\"%s\" "
			    "leaf=%zu leaf2=%zu eval=%.6fs "
			    "intersection_candidates=%zu "
			    "segments=(bige:%zu boolweave:%zu)\n",
			    dgcdp->profile_current.name.empty() ? "(null)" :
			    dgcdp->profile_current.name.c_str(),
			    leaf_no,
			    leaf2,
			    ((double)pair_eval_usec) / 1000000.0,
			    dgcdp->profile_current.intersection_candidates -
			    pair_eval_start_candidates,
			    dgcdp->profile_current.bige_segments -
			    pair_eval_start_bige_segments,
			    dgcdp->profile_current.boolweave_segments -
			    pair_eval_start_boolweave_segments);
		}
	    }
	    }
	}
    if (dgcdp->profile_enabled)
	dgcdp->profile_current.intersection_candidate_usec +=
	    bu_gettime() - stage_start;

    bu_free((char *)hits1, "hits1");
    bu_free((char *)hits2, "hits2");
}


static void
free_etree(union E_tree *eptr,
	   struct bigE_data *dgcdp)
{
    CK_ETREE(eptr);

    switch (eptr->l.op) {
	case OP_UNION:
	case OP_SUBTRACT:
	case OP_INTERSECT:
	    free_etree(eptr->n.left, dgcdp);
	    free_etree(eptr->n.right, dgcdp);
	    bu_free((char *)eptr, "node pointer");
	    break;
	case OP_DB_LEAF:
	case OP_SOLID:
	    if (eptr->l.m && !eptr->l.do_not_free_model) {
		nmg_km(eptr->l.m);
		eptr->l.m = (struct model *)NULL;
	    }
	    if (BU_LIST_NON_EMPTY(&eptr->l.seghead)) {
		bigE_seg_list_free(&eptr->l.seghead, dgcdp);
	    }
	    if (eptr->l.edge_list_initialized) {
		bu_ptbl_free(&eptr->l.edge_list);
		eptr->l.edge_list_initialized = 0;
	    }
	    if (eptr->l.stp) {
		if (eptr->l.stp->st_specific && OBJ[eptr->l.stp->st_id].ft_free)
		    OBJ[eptr->l.stp->st_id].ft_free(eptr->l.stp);
		bu_free((char *)eptr->l.stp, "struct soltab");
	    }

	    bu_free((char *)eptr, "leaf pointer");
	    break;
    }
}


/* convert all "half" solids to polysolids */
static void
fix_halfs(struct bigE_data *dgcdp)
{
    point_t max, min;
    size_t i;
    size_t count=0;
    const struct bn_tol *tol = dgcdp->tol;

    VSETALL(max, -INFINITY);
    VSETALL(min, INFINITY);

    for (i = 0; i < BU_PTBL_LEN(&dgcdp->leaf_list); i++) {
	union E_tree *tp;

	tp = (union E_tree *)BU_PTBL_GET(&dgcdp->leaf_list, i);

	if (tp->l.stp->st_id == ID_HALF)
	    continue;

	VMINMAX(min, max, tp->l.stp->st_min);
	VMINMAX(min, max, tp->l.stp->st_max);
	count++;
    }

    if (!count) {
	return;
    }

    for (i = 0; i < BU_PTBL_LEN(&dgcdp->leaf_list); i++) {
	union E_tree *tp;
	struct vertex *v[8];
	struct vertex **vp[4];
	struct nmgregion *r;
	struct shell *s;
	struct rt_pg_internal *pg;
	struct faceuse *fu;
	plane_t haf_pl;
	struct half_specific *hp;
	int j;
	struct bu_list *vlfree = dgcdp->vlfree;

	tp = (union E_tree *)BU_PTBL_GET(&dgcdp->leaf_list, i);

	if (tp->l.stp->st_id != ID_HALF)
	    continue;

	hp = (struct half_specific *)tp->l.stp->st_specific;

	HMOVE(haf_pl, hp->half_eqn);

	if (DIST_PNT_PLANE(max, haf_pl) >= -tol->dist &&
	    DIST_PNT_PLANE(min, haf_pl) >= -tol->dist)
	    continue;

	/* make an NMG the size of our model bounding box */
	tp->l.m = nmg_mm();
	r = nmg_mrsv(tp->l.m);
	s = BU_LIST_FIRST(shell, &r->s_hd);

	for (j = 0; j < 8; j++)
	    v[j] = (struct vertex *)NULL;

	vp[0] = &v[0];
	vp[1] = &v[1];
	vp[2] = &v[2];
	vp[3] = &v[3];
	fu = nmg_cmface(s, vp, 4);
	nmg_vertex_g(v[0], max[X], min[Y], min[Z]);
	nmg_vertex_g(v[1], max[X], max[Y], min[Z]);
	nmg_vertex_g(v[2], max[X], max[Y], max[Z]);
	nmg_vertex_g(v[3], max[X], min[Y], max[Z]);
	nmg_calc_face_g(fu, vlfree);

	vp[0] = &v[4];
	vp[1] = &v[5];
	vp[2] = &v[6];
	vp[3] = &v[7];
	fu = nmg_cmface(s, vp, 4);
	nmg_vertex_g(v[4], min[X], min[Y], min[Z]);
	nmg_vertex_g(v[5], min[X], min[Y], max[Z]);
	nmg_vertex_g(v[6], min[X], max[Y], max[Z]);
	nmg_vertex_g(v[7], min[X], max[Y], min[Z]);
	nmg_calc_face_g(fu, vlfree);

	vp[0] = &v[0];
	vp[1] = &v[3];
	vp[2] = &v[5];
	vp[3] = &v[4];
	fu = nmg_cmface(s, vp, 4);
	nmg_calc_face_g(fu, vlfree);

	vp[0] = &v[1];
	vp[1] = &v[7];
	vp[2] = &v[6];
	vp[3] = &v[2];
	fu = nmg_cmface(s, vp, 4);
	nmg_calc_face_g(fu, vlfree);

	vp[0] = &v[3];
	vp[1] = &v[2];
	vp[2] = &v[6];
	vp[3] = &v[5];
	fu = nmg_cmface(s, vp, 4);
	nmg_calc_face_g(fu, vlfree);

	vp[0] = &v[1];
	vp[1] = &v[0];
	vp[2] = &v[4];
	vp[3] = &v[7];
	fu = nmg_cmface(s, vp, 4);
	nmg_calc_face_g(fu, vlfree);

	nmg_region_a(r, tol);

	for (BU_LIST_FOR (fu, faceuse, &s->fu_hd)) {
	    struct edgeuse *eu, *new_eu;
	    struct loopuse *lu, *new_lu;
	    plane_t pl;
	    struct vertexuse *vcut[2];
	    point_t pt[2];
	    struct edgeuse *eu_split[2];

	    if (fu->orientation != OT_SAME)
		continue;

	    bigE_get_fu_plane(pl, fu);

	    if (bg_coplanar(pl, haf_pl, tol) > 0)
		continue;

	    lu = BU_LIST_FIRST(loopuse, &fu->lu_hd);

	    count = 0;
	    for (BU_LIST_FOR (eu, edgeuse, &lu->down_hd)) {
		vect_t dir;
		struct vertex_g *v1g, *v2g;
		fastf_t dist;

		v1g = eu->vu_p->v_p->vg_p;
		v2g = eu->eumate_p->vu_p->v_p->vg_p;

		VSUB2(dir, v2g->coord, v1g->coord);

		if (bg_isect_line3_plane(&dist, v1g->coord, dir, haf_pl, tol) < 1)
		    continue;

		if (dist < 0.0 || dist >=1.0)
		    continue;

		VJOIN1(pt[count], v1g->coord, dist, dir);
		eu_split[count] = eu;

		count++;
		if (count == 2)
		    break;
	    }

	    if (count != 2)
		continue;

	    new_eu = nmg_eusplit((struct vertex *)NULL, eu_split[0], 1);
	    vcut[0] = new_eu->vu_p;
	    nmg_vertex_gv(vcut[0]->v_p, pt[0]);

	    new_eu = nmg_eusplit((struct vertex *)NULL, eu_split[1], 1);
	    vcut[1] = new_eu->vu_p;
	    nmg_vertex_gv(vcut[1]->v_p, pt[1]);

	    new_lu = nmg_cut_loop(vcut[0], vcut[1], vlfree);
	    nmg_lu_reorient(lu);
	    nmg_lu_reorient(new_lu);

	    for (BU_LIST_FOR (eu, edgeuse, &lu->down_hd)) {
		if (eu->vu_p->v_p == vcut[0]->v_p || eu->vu_p->v_p == vcut[1]->v_p)
		    continue;

		if (DIST_PNT_PLANE(eu->vu_p->v_p->vg_p->coord, haf_pl) > tol->dist) {
		    nmg_klu(lu);
		    break;
		} else {
		    nmg_klu(new_lu);
		    break;
		}
	    }
	}

	/* kill any faces outside the half */
	fu = BU_LIST_FIRST(faceuse, &s->fu_hd);
	if (fu->orientation != OT_SAME)
	    fu = BU_LIST_PNEXT(faceuse, &fu->l);
	while (BU_LIST_NOT_HEAD(&fu->l, &s->fu_hd)) {
	    struct faceuse *next_fu;
	    struct loopuse *lu;
	    int killfu=0;

	    next_fu = BU_LIST_PNEXT(faceuse, &fu->l);
	    if (fu->fumate_p == next_fu)
		next_fu = BU_LIST_PNEXT(faceuse, &next_fu->l);

	    if (fu->orientation != OT_SAME) {
		fu = next_fu;
		continue;
	    }

	    lu = BU_LIST_FIRST(loopuse, &fu->lu_hd);
	    while (BU_LIST_NOT_HEAD(&lu->l, &fu->lu_hd)) {
		struct loopuse *next_lu;
		struct edgeuse *eu;
		int killit;

		next_lu = BU_LIST_PNEXT(loopuse, &lu->l);

		killit = 0;
		for (BU_LIST_FOR (eu, edgeuse, &lu->down_hd)) {
		    struct vertex_g *vg;

		    vg = eu->vu_p->v_p->vg_p;

		    if (DIST_PNT_PLANE(vg->coord, haf_pl) > tol->dist) {
			killit = 1;
			break;
		    }
		}

		if (killit) {
		    if (nmg_klu(lu)) {
			killfu = 1;
			break;
		    }
		}
		lu = next_lu;
	    }

	    if (killfu)
		nmg_kfu(fu);

	    fu = next_fu;
	}

	nmg_rebound(tp->l.m, tol);
	nmg_model_fuse(tp->l.m, vlfree, tol);
	nmg_close_shell(s, vlfree, tol);
	nmg_rebound(tp->l.m, tol);

	BU_ALLOC(pg, struct rt_pg_internal);

	if (!nmg_to_poly(tp->l.m, pg, vlfree, tol)) {
	    bu_free((char *)pg, "rt_pg_internal");
	} else {
	    struct rt_db_internal intern2;

	    RT_DB_INTERNAL_INIT(&intern2);
	    intern2.idb_major_type = DB5_MAJORTYPE_BRLCAD;
	    intern2.idb_type = ID_POLY;
	    intern2.idb_meth = &OBJ[ID_POLY];
	    intern2.idb_ptr = (void *)pg;
	    if (OBJ[tp->l.stp->st_id].ft_free)
		OBJ[tp->l.stp->st_id].ft_free(tp->l.stp);
	    tp->l.stp->st_specific = NULL;
	    tp->l.stp->st_id = ID_POLY;
	    VSETALL(tp->l.stp->st_max, -INFINITY);
	    VSETALL(tp->l.stp->st_min,  INFINITY);
	    rt_obj_prep(tp->l.stp, &intern2, dgcdp->rtip);

	    rt_db_free_internal(&intern2);
	}
    }
}


static const char *
bigE_skip_leading_slash(const char *path)
{
    if (!path)
	return NULL;
    while (*path == '/')
	path++;
    return path;
}


struct bigE_region_parallel_state {
    size_t region_count;
    std::vector<bigE_region_result> *results;
    std::vector<bigE_data *> *workers;
    std::atomic<size_t> next_index;
    std::atomic<int> failed;

    bigE_region_parallel_state() :
	region_count(0),
	results(NULL),
	workers(NULL),
	next_index(0),
	failed(0)
    {}
};


static void
bigE_region_parallel_worker(int cpu, void *data)
{
    struct bigE_region_parallel_state *state =
	(struct bigE_region_parallel_state *)data;
    if (!state || !state->results ||
	    !state->workers || cpu < 0)
	return;

    size_t worker_index = (size_t)cpu;
    if (worker_index >= state->workers->size()) {
	state->failed.store(1);
	return;
    }

    struct bigE_data *worker = (*state->workers)[worker_index];
    if (!worker) {
	state->failed.store(1);
	return;
    }

    while (!state->failed.load()) {
	size_t index = state->next_index.fetch_add(1);
	if (index >= state->region_count)
	    break;

	struct bigE_region_result result;
	result.index = index;
	struct region *rp = bigE_region_at(worker->rtip, index);
	const char *debug_name = bigE_region_debug_name(rp);
	if (debug_name)
	    result.debug_name = debug_name;
	if (!rp) {
	    result.ret = BRLCAD_ERROR;
	    (*state->results)[index] = result;
	    state->failed.store(1);
	    break;
	}

	result.ret = bigE_region_evaluate(worker, rp, index, &result.lines);
	(*state->results)[index] = result;
	if (result.ret != BRLCAD_OK) {
	    state->failed.store(1);
	    break;
	}
    }
}


static int
bigE_regions_evaluate(struct db_i *dbip,
		      const char *path,
		      const struct bn_tol *tol,
		      const struct bg_tess_tol *ttol,
		      uint32_t eval_flags,
		      int ncpu,
		      size_t region_count,
		      struct bigE_line_set *lines_out)
{
    if (!dbip || !path || !tol || !ttol || !lines_out)
	return BRLCAD_ERROR;

    lines_out->records.clear();
    if (region_count == 0)
	return BRLCAD_OK;

    /*
     * Keep this seam serial until the parallel worker path is proven
     * image-equivalent.  Multiple prepared worker contexts are equivalent when
     * used in serial region order, but dynamic out-of-order region execution
     * perturbs the evaluated-wire baseline.  The next parallel layer should
     * preserve region-order semantics and move only isolated ray/classification
     * jobs behind per-thread resources and deterministic merge.
     */
    (void)ncpu;
    size_t nworkers = 1;
    if (nworkers < 1)
	nworkers = 1;

    std::vector<bigE_region_result> results(region_count);
    std::vector<bigE_data *> workers(nworkers, NULL);
    int ret = BRLCAD_OK;
    size_t merged_record_count = 0;
    const int profile_pipeline =
	(eval_flags & RT_EVAL_WIREFRAME_F_PROFILE) ||
	bigE_env_enabled("RT_EVAL_WIREFRAME_PROFILE");
    const int64_t profile_start = profile_pipeline ? bu_gettime() : 0;
    int64_t profile_prepare = profile_start;
    int64_t profile_work = profile_start;
    int64_t profile_merge = profile_start;
    int64_t profile_worker_cleanup = profile_start;

    for (size_t i = 0; i < nworkers; i++) {
	workers[i] = new bigE_data;
	if (bigE_context_prepare(workers[i], dbip, path, tol, ttol,
		eval_flags) != BRLCAD_OK) {
	    ret = BRLCAD_ERROR;
	    goto cleanup;
	}
    }
    if (profile_pipeline)
	profile_prepare = bu_gettime();

    if (nworkers == 1) {
	for (size_t i = 0; i < region_count; i++) {
	    results[i].index = i;
	    struct region *rp = bigE_region_at(workers[0]->rtip, i);
	    const char *debug_name = bigE_region_debug_name(rp);
	    if (debug_name)
		results[i].debug_name = debug_name;
	    if (!rp ||
		bigE_region_evaluate(workers[0], rp, i, &results[i].lines) != BRLCAD_OK) {
		ret = BRLCAD_ERROR;
		goto cleanup;
	    }
	}
    } else {
	struct bigE_region_parallel_state state;
	state.region_count = region_count;
	state.results = &results;
	state.workers = &workers;
	bu_parallel(bigE_region_parallel_worker, nworkers, &state);
	if (state.failed.load()) {
	    ret = BRLCAD_ERROR;
	    goto cleanup;
	}
    }
    if (profile_pipeline)
	profile_work = bu_gettime();

    merged_record_count = lines_out->records.size();
    for (size_t i = 0; i < results.size(); i++) {
	if (results[i].ret != BRLCAD_OK || results[i].index != i) {
	    ret = BRLCAD_ERROR;
	    goto cleanup;
	}
	const size_t result_count = results[i].lines.records.size();
	if (result_count >
	    std::numeric_limits<size_t>::max() - merged_record_count) {
	    ret = BRLCAD_ERROR;
	    goto cleanup;
	}
	merged_record_count += result_count;
    }
    lines_out->records.reserve(merged_record_count);
    for (size_t i = 0; i < results.size(); i++) {
	std::vector<bigE_line_record> &records = results[i].lines.records;
	lines_out->records.insert(lines_out->records.end(),
	    std::make_move_iterator(records.begin()),
	    std::make_move_iterator(records.end()));
    }
    if (profile_pipeline)
	profile_merge = bu_gettime();

cleanup:
    for (size_t i = 0; i < workers.size(); i++) {
	delete workers[i];
	workers[i] = NULL;
    }
    if (profile_pipeline)
	profile_worker_cleanup = bu_gettime();
    results.clear();
    if (profile_pipeline && ret == BRLCAD_OK) {
	const int64_t profile_result_cleanup = bu_gettime();
	bu_log("evaluated-wire evaluate stages: prepare=%.6fs work=%.6fs "
	    "merge=%.6fs worker_cleanup=%.6fs result_cleanup=%.6fs "
	    "records=%zu\n",
	    (profile_prepare - profile_start) / 1000000.0,
	    (profile_work - profile_prepare) / 1000000.0,
	    (profile_merge - profile_work) / 1000000.0,
	    (profile_worker_cleanup - profile_merge) / 1000000.0,
	    (profile_result_cleanup - profile_worker_cleanup) / 1000000.0,
	    lines_out->records.size());
    }
    if (ret != BRLCAD_OK)
	lines_out->records.clear();
    return ret;
}


int
rt_eval_wireframe(struct bu_list *vhead,
		  struct bu_list *vlfree,
		  struct db_i *dbip,
		  const char *path,
		  const struct bn_tol *tol,
		  const struct bg_tess_tol *ttol,
		  const struct rt_eval_wireframe_opts *opts)
{
    int ret = BRLCAD_OK;
    int ncpu = 1;
    size_t region_count = 0;
    size_t output_record_count = 0;
    uint32_t eval_flags = RT_EVAL_WIREFRAME_F_DEFAULT;
    struct bigE_line_set lines;
    int64_t profile_start = 0;
    int64_t profile_discover = 0;
    int64_t profile_evaluate = 0;
    int profile_pipeline = 0;

    if (!vhead || !vlfree || !dbip || !path || !tol || !ttol)
	return BRLCAD_ERROR;
    path = bigE_skip_leading_slash(path);
    if (!path[0])
	return BRLCAD_ERROR;
    if (opts) {
	eval_flags = opts->flags;
	if (opts->ncpu > 0)
	    ncpu = opts->ncpu;
    }
    if ((eval_flags & RT_EVAL_WIREFRAME_F_BIGE) &&
	    (eval_flags & RT_EVAL_WIREFRAME_F_BOOLWEAVE))
	return BRLCAD_ERROR;
    if (ncpu < 1)
	ncpu = 1;
    profile_pipeline = (eval_flags & RT_EVAL_WIREFRAME_F_PROFILE) ||
	bigE_env_enabled("RT_EVAL_WIREFRAME_PROFILE");
    if (profile_pipeline)
	profile_start = bu_gettime();

    ret = bigE_region_count_discover(dbip, path, tol, ttol, eval_flags,
	    &region_count);
    if (ret != BRLCAD_OK)
	return ret;
    if (profile_pipeline)
	profile_discover = bu_gettime();

    ret = bigE_regions_evaluate(dbip, path, tol, ttol, eval_flags, ncpu,
	    region_count, &lines);
    if (ret != BRLCAD_OK)
	return ret;
    if (profile_pipeline)
	profile_evaluate = bu_gettime();

    output_record_count = lines.records.size();
    ret = bigE_line_set_export_vlist(&lines, vhead, vlfree);
    if (profile_pipeline) {
	const int64_t profile_end = bu_gettime();
	bu_log("evaluated-wire pipeline stages: discover=%.6fs "
	    "evaluate=%.6fs export=%.6fs total=%.6fs regions=%zu records=%zu\n",
	    (profile_discover - profile_start) / 1000000.0,
	    (profile_evaluate - profile_discover) / 1000000.0,
	    (profile_end - profile_evaluate) / 1000000.0,
	    (profile_end - profile_start) / 1000000.0,
	    region_count, output_record_count);
    }
    return ret;
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
