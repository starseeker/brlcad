/*               G E D _ D R A W _ S C E N E . C P P
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
/** @file ged_draw_scene.cpp
 *
 * Exercises the ged_draw_* scene API declared in include/ged/draw.h.
 * The API is backed by the GED BSG draw tree and must remain usable without
 * an attached display manager.
 *
 * Usage: ged_test_ged_draw_scene <directory-containing-moss.g>
 */

#include "common.h"

#include <fstream>
#include <cctype>
#include <cstdio>
#include <cstring>

#include <bu.h>
#include "dm.h"
#include <ged.h>
#include "ged/draw.h"
#include "ged/db_index.h"
#include "ged/selection_state.h"
#include "bg/plot3.h"
#include "bg/line_layer.h"
#include "nmg.h"
#include "rt/db_attr.h"
#include "wdb.h"
#include "../../bsg_ged_draw_private.h"
#include "../../ged_private.h"

typedef void test_scene_context;

static void *
test_active_view_ctx(struct ged *gedp)
{
    return gedp ? ged_view_active_ctx(gedp) : NULL;
}


#define ASSERT(cond) do { \
    nchecks++; \
    if (!(cond)) { \
	bu_log("FAIL [%s:%d] %s\n", __FILE__, __LINE__, #cond); \
	nfails++; \
    } \
} while (0)

static int nchecks = 0;

static int
test_uplot_stream_process_args(struct ged_uplot_stream *stream, int command,
			       const char *args)
{
    FILE *fp = tmpfile();
    if (!fp)
	return BRLCAD_ERROR;
    if (args)
	std::fputs(args, fp);
    std::rewind(fp);
    int ret = _ged_uplot_stream_process(stream, fp, command);
    std::fclose(fp);
    return ret;
}

static int nfails  = 0;

#define TEST_LINE_VIEW_MAX_SAMPLES 8

struct draw_observer_test_state {
    struct ged *gedp;
    ged_draw_observer_token self_token;
    int calls;
    ged_draw_transaction_kind last_kind;
    int last_status;
    uint64_t last_before;
    uint64_t last_after;
    struct bu_vls names;
};

struct event_observer_test_state {
    struct ged *gedp;
    ged_event_observer_token self_token;
    int calls;
    ged_event_kind last_kind;
    size_t last_event_count;
    size_t last_coalesced_count;
    int last_status;
    int last_redraw;
    uint64_t last_before;
    uint64_t last_after;
    struct bu_vls names;
};

struct event_observer_followup_state {
    ged_event_observer_token self_token;
    int calls;
    ged_event_kind kinds[4];
};

static void
draw_observer_test_cb(struct ged *gedp,
		      const struct ged_draw_transaction *txn,
		      const struct ged_draw_transaction_result *result,
		      void *client_data)
{
    (void)gedp;

    struct draw_observer_test_state *state =
	(struct draw_observer_test_state *)client_data;
    if (!state)
	return;

    state->calls++;
    state->last_kind = txn ? txn->kind : GED_DRAW_TXN_NONE;
    state->last_status = result ? result->status : -999;
    state->last_before = result ? result->scene_revision_before : 0;
    state->last_after = result ? result->scene_revision_after : 0;
    if (BU_VLS_IS_INITIALIZED(&state->names)) {
	bu_vls_trunc(&state->names, 0);
	if (result && BU_VLS_IS_INITIALIZED(&result->names))
	    bu_vls_strcat(&state->names, bu_vls_cstr(&result->names));
    }
}

static void
event_observer_test_cb(struct ged *gedp,
		       const struct ged_event *events,
		       size_t event_count,
		       const struct ged_event_txn_result *result,
		       void *client_data)
{
    (void)gedp;

    struct event_observer_test_state *state =
	(struct event_observer_test_state *)client_data;
    if (!state)
	return;

    state->calls++;
    state->last_kind = (events && event_count) ? events[0].kind : GED_EVENT_NONE;
    state->last_event_count = result ? result->event_count : event_count;
    state->last_coalesced_count = result ? result->coalesced_event_count : 0;
    state->last_status = result ? result->status : -999;
    state->last_redraw = (events && event_count) ? events[0].redraw : 0;
    state->last_before = result ? result->draw_scene_revision_before : 0;
    state->last_after = result ? result->draw_scene_revision_after : 0;
    if (BU_VLS_IS_INITIALIZED(&state->names)) {
	bu_vls_trunc(&state->names, 0);
	if (result && BU_VLS_IS_INITIALIZED(&result->affected_names))
	    bu_vls_strcat(&state->names, bu_vls_cstr(&result->affected_names));
    }
}

static void
event_observer_followup_cb(struct ged *gedp,
			   const struct ged_event *events,
			   size_t event_count,
			   const struct ged_event_txn_result *result,
			   void *client_data)
{
    (void)result;

    struct event_observer_followup_state *state =
	(struct event_observer_followup_state *)client_data;
    if (!state)
	return;

    if (state->calls < 4)
	state->kinds[state->calls] =
	    (events && event_count) ? events[0].kind : GED_EVENT_NONE;
    state->calls++;

    if (gedp && state->calls == 1)
	ged_event_notify_attribute_changed(gedp, "all.g", 0, NULL);
}

static void
event_observer_self_remove_cb(struct ged *gedp,
			      const struct ged_event *events,
			      size_t event_count,
			      const struct ged_event_txn_result *result,
			      void *client_data)
{
    (void)events;
    (void)event_count;
    (void)result;

    struct event_observer_test_state *state =
	(struct event_observer_test_state *)client_data;
    if (!state)
	return;

    state->calls++;
    if (gedp && state->self_token)
	ged_event_observer_remove(gedp, state->self_token);
}

static void
draw_observer_self_remove_cb(struct ged *gedp,
			     const struct ged_draw_transaction *txn,
			     const struct ged_draw_transaction_result *result,
			     void *client_data)
{
    (void)txn;
    (void)result;

    struct draw_observer_test_state *state =
	(struct draw_observer_test_state *)client_data;
    if (!state)
	return;

    state->calls++;
    if (gedp && state->self_token)
	ged_draw_observer_remove(gedp, state->self_token);
}

static int
_scene_group_count_cb(const struct ged_draw_group_record * /*rec*/, void *ud)
{
    int *n = (int *)ud;
    (*n)++;
    return 1;
}

static int
scene_group_count(struct ged *gedp)
{
    int n = 0;
    ged_draw_foreach_group_record(gedp, _scene_group_count_cb, &n);
    return n;
}


static int
vls_has_line(const struct bu_vls *vls, const char *line)
{
    if (!vls || !line)
	return 0;

    size_t line_len = strlen(line);
    const char *start = bu_vls_cstr(vls);
    while (start && *start) {
	const char *end = strchr(start, '\n');
	size_t len = end ? (size_t)(end - start) : strlen(start);
	if (len == line_len && bu_strncmp(start, line, line_len) == 0)
	    return 1;
	if (!end)
	    break;
	start = end + 1;
    }

    return 0;
}


static int
vls_has_word(const struct bu_vls *vls, const char *word)
{
    if (!vls || !word)
	return 0;

    size_t word_len = strlen(word);
    const char *start = bu_vls_cstr(vls);
    while (start && *start) {
	while (*start && std::isspace((unsigned char)*start))
	    start++;
	if (!*start)
	    break;
	const char *end = start;
	while (*end && !std::isspace((unsigned char)*end))
	    end++;
	size_t len = (size_t)(end - start);
	if (len == word_len && bu_strncmp(start, word, word_len) == 0)
	    return 1;
	start = end;
    }

    return 0;
}

static int
test_db_path_names_equal(const char *lhs, const char *rhs)
{
    if (!lhs || !rhs)
	return 0;
    if (lhs[0] == '/')
	lhs++;
    if (rhs[0] == '/')
	rhs++;
    return BU_STR_EQUAL(lhs, rhs);
}


struct test_shape_record_path_ctx {
    const char *path;
    struct ged_draw_shape_record *out;
    int found;
};


static int
test_shape_record_by_path_cb(const struct ged_draw_shape_record *rec, void *ud)
{
    struct test_shape_record_path_ctx *ctx =
	(struct test_shape_record_path_ctx *)ud;
    if (!ctx || !ctx->path || !rec)
	return 1;

    int matched = 0;
    if (rec->display_name)
	matched = test_db_path_names_equal(rec->display_name, ctx->path);
    if (!matched && rec->fullpath) {
	char *fp_path = db_path_to_string(rec->fullpath);
	if (fp_path) {
	    matched = test_db_path_names_equal(fp_path, ctx->path);
	    bu_free(fp_path, "shape record db path");
	}
    }

    if (!matched)
	return 1;

    if (ctx->out)
	*ctx->out = *rec;
    ctx->found = 1;
    return 0;
}


static int
test_shape_record_by_path(struct ged *gedp,
			  const char *path,
			  struct ged_draw_shape_record *out)
{
    struct test_shape_record_path_ctx ctx;
    ctx.path = path;
    ctx.out = out;
    ctx.found = 0;
    ged_draw_foreach_shape_record(gedp, test_shape_record_by_path_cb, &ctx);
    return ctx.found;
}


static test_scene_context *
test_scene_root(struct ged *gedp)
{
    return gedp ? (test_scene_context *)ged_draw_view_context_scene_root(
	    test_active_view_ctx(gedp)) : NULL;
}


static size_t
test_scene_child_count(test_scene_context *node)
{
    struct ged_draw_scene_tree_summary tree_summary;
    return ged_draw_scene_context_tree_summary(node, &tree_summary) &&
	tree_summary.valid ? tree_summary.child_count : 0;
}


static test_scene_context *
test_scene_child_at(test_scene_context *node, size_t idx)
{
    return (test_scene_context *)ged_draw_scene_context_child_at(node, idx);
}


static test_scene_context *
test_scene_parent(test_scene_context *node)
{
    return (test_scene_context *)ged_draw_scene_context_parent(node);
}


static const char *
test_scene_name(test_scene_context *node)
{
    return ged_draw_scene_context_name(node);
}


static int
test_scene_is_group(test_scene_context *node)
{
    struct ged_draw_scene_tree_summary tree_summary;
    return ged_draw_scene_context_tree_summary(node, &tree_summary) &&
	tree_summary.valid && tree_summary.is_group;
}


static int
test_scene_is_shape(test_scene_context *node)
{
    struct ged_draw_scene_tree_summary tree_summary;
    return ged_draw_scene_context_tree_summary(node, &tree_summary) &&
	tree_summary.valid && tree_summary.is_shape;
}


static test_scene_context *
test_find_shape_by_dbpath(test_scene_context *node, const char *path)
{
    if (!node || !path)
	return NULL;

    if (test_scene_is_shape(node)) {
	const struct db_full_path *fp = ged_draw_scene_context_fullpath(node);
	if (fp && fp->fp_len > 0) {
	    char *fp_path = db_path_to_string(fp);
	    int matched = 0;
	    if (fp_path) {
		const char *lhs = fp_path;
		const char *rhs = path;
		if (lhs[0] == '/')
		    lhs++;
		if (rhs[0] == '/')
		    rhs++;
		matched = BU_STR_EQUAL(lhs, rhs);
		bu_free(fp_path, "test shape dbpath string");
	    }
	    if (matched)
		return node;
	}
    }

    for (size_t i = 0; i < test_scene_child_count(node); i++) {
	test_scene_context *found =
	    test_find_shape_by_dbpath(test_scene_child_at(node, i), path);
	if (found)
	    return found;
    }
    return NULL;
}


struct segment_count_state {
    size_t count;
};


static int
count_segment_cb(const point_t UNUSED(a), const point_t UNUSED(b), void *ud)
{
    struct segment_count_state *state = (struct segment_count_state *)ud;

    if (state)
	state->count++;
    return 1;
}


struct point_count_state {
    size_t count;
};


static int
count_point_cb(const point_t UNUSED(pt), void *ud)
{
    struct point_count_state *state = (struct point_count_state *)ud;

    if (state)
	state->count++;
    return 1;
}


static int
source_name_matches_drawn_prefix(const char *name, const char *prefix)
{
    size_t nlen = 0;

    if (!name || !prefix)
	return 0;
    if (name[0] == '/')
	name++;
    if (prefix[0] == '/')
	prefix++;
    nlen = strlen(prefix);
    return BU_STR_EQUAL(name, prefix) ||
	(bu_strncmp(name, prefix, nlen) == 0 && name[nlen] == '/');
}


static test_scene_context *
ged_draw_first_shape(struct ged *gedp)
{
    return (test_scene_context *)ged_draw_first_shape_context(gedp);
}


static int
test_shape_record_at(struct ged *gedp, int idx, struct ged_draw_shape_record *rec)
{
    if (!rec)
	return 0;

    memset(rec, 0, sizeof(*rec));
    ged_draw_shape_ref ref = ged_draw_shape_ref_at(gedp, idx);
    if (ged_draw_shape_ref_is_null(ref))
	return 0;

    return ged_draw_shape_record_get(gedp, ref, rec);
}


static int
test_first_shape_group_records(struct ged *gedp,
			       struct ged_draw_shape_record *shape_rec,
			       struct ged_draw_group_record *group_rec)
{
    if (!shape_rec || !group_rec)
	return 0;

    memset(shape_rec, 0, sizeof(*shape_rec));
    memset(group_rec, 0, sizeof(*group_rec));
    ged_draw_shape_ref shape_ref = ged_draw_first_shape_ref(gedp);
    if (ged_draw_shape_ref_is_null(shape_ref) ||
	    !ged_draw_shape_record_get(gedp, shape_ref, shape_rec))
	return 0;

    ged_draw_group_ref group_ref = ged_draw_group_ref_of_shape(gedp, shape_ref);
    if (ged_draw_group_ref_is_null(group_ref))
	return 0;

    return ged_draw_group_record_get(gedp, group_ref, group_rec);
}


struct test_selected_view_record_ctx {
    const char *path;
    int saw_selected;
    int saw_path_selected;
};


static int
test_selected_view_record_cb(const struct ged_draw_view_db_object_record *rec,
			     void *ud)
{
    struct test_selected_view_record_ctx *ctx =
	(struct test_selected_view_record_ctx *)ud;

    if (!ctx || !rec || !rec->is_database_source || !rec->selected)
	return 1;

    ctx->saw_selected = 1;
    if (test_db_path_names_equal(rec->path, ctx->path)) {
	ctx->saw_path_selected = 1;
	return 0;
    }

    return 1;
}


static int
test_selected_visible_view_record(struct ged *gedp, const char *path)
{
    struct test_selected_view_record_ctx ctx = {path, 0, 0};

    ged_draw_foreach_visible_view_db_object_record(
	    test_active_view_ctx(gedp),
	    test_selected_view_record_cb,
	    &ctx);

    return ctx.saw_path_selected || ctx.saw_selected;
}


struct test_qray_view_record_ctx {
    const char *prefix;
    size_t items;
    size_t segments;
    int saw_red;
    int saw_blue;
};


static int
test_qray_view_record_cb(const struct ged_draw_view_db_object_record *rec,
			 void *ud)
{
    struct test_qray_view_record_ctx *ctx =
	(struct test_qray_view_record_ctx *)ud;
    const char *path = rec ? rec->path : NULL;

    if (!ctx || !ctx->prefix || !rec || !path || !rec->type_name ||
	    bu_strncmp(path, ctx->prefix, strlen(ctx->prefix)) != 0 ||
	    !BU_STR_EQUAL(rec->type_name, "line") ||
	    rec->vlist_point_count == 0)
	return 1;

    ASSERT(rec->is_view_source);
    ASSERT(!rec->is_database_source);
    ASSERT(rec->non_database_source);

    struct segment_count_state segs = {0};
    ASSERT(ged_draw_view_db_object_record_foreach_segment(rec,
	    count_segment_cb, &segs) == 1);
    ctx->items++;
    ctx->segments += segs.count;
    if (rec->color[0] == 255 && rec->color[1] == 0 && rec->color[2] == 0)
	ctx->saw_red = 1;
    if (rec->color[0] == 0 && rec->color[1] == 0 && rec->color[2] == 255)
	ctx->saw_blue = 1;

    return 1;
}


static int
test_qray_view_records(struct ged *gedp, const char *prefix)
{
    struct test_qray_view_record_ctx ctx = {prefix, 0, 0, 0, 0};
    struct ged_draw_view_record_query query;

    memset(&query, 0, sizeof(query));
    query.flags = GED_DRAW_VIEW_RECORD_QUERY_VIEW_OBJECTS;
    query.draw_mode = -1;
    ged_draw_foreach_view_record_query(
	    test_active_view_ctx(gedp),
	    &query,
	    test_qray_view_record_cb,
	    &ctx);

    return ctx.items == 2 && ctx.segments == 2 && ctx.saw_red && ctx.saw_blue;
}


struct test_view_record_geometry_count_ctx {
    const char *prefix;
    const char *geometry_name;
    size_t count;
};


static int
test_view_record_geometry_count_cb(
	const struct ged_draw_view_db_object_record *rec,
	void *ud)
{
    struct test_view_record_geometry_count_ctx *ctx =
	(struct test_view_record_geometry_count_ctx *)ud;

    if (!ctx || !rec || !ctx->prefix || !ctx->geometry_name ||
	    !rec->path || !rec->geometry_name)
	return 1;

    if (source_name_matches_drawn_prefix(rec->path, ctx->prefix) &&
	    BU_STR_EQUAL(rec->geometry_name, ctx->geometry_name))
	ctx->count++;

    return 1;
}


static size_t
test_visible_db_view_record_geometry_count(struct ged *gedp,
					   const char *prefix,
					   const char *geometry_name)
{
    struct test_view_record_geometry_count_ctx ctx = {prefix, geometry_name, 0};

    ged_draw_foreach_visible_view_db_object_record(
	    test_active_view_ctx(gedp),
	    test_view_record_geometry_count_cb,
	    &ctx);

    return ctx.count;
}


struct test_annotation_view_record_ctx {
    const char *path;
    const char *text;
    uint64_t *cache_identity;
    uint64_t *source_identity;
    int found;
};


static int
test_annotation_view_record_cb(
	const struct ged_draw_view_db_object_record *rec,
	void *ud)
{
    struct test_annotation_view_record_ctx *ctx =
	(struct test_annotation_view_record_ctx *)ud;

    if (!ctx || !rec || !ctx->path || !ctx->text || !rec->path ||
	    !rec->geometry_name ||
	    !source_name_matches_drawn_prefix(rec->path, ctx->path) ||
	    !BU_STR_EQUAL(rec->geometry_name, "annotation"))
	return 1;

    ASSERT(rec->is_database_source);
    ASSERT(!rec->non_database_source);

    struct ged_draw_view_annotation_summary summary;
    ASSERT(ged_draw_view_db_object_record_annotation_summary(rec, 1, &summary));
    ASSERT(summary.valid);
    ASSERT(summary.display_space);
    ASSERT(summary.point_count == 2);
    ASSERT(summary.segment_count == 2);
    ASSERT(summary.has_point);
    if (summary.has_point) {
	point_t annot_local_expected;
	VSET(annot_local_expected, 0.2, 0.3, 0.0);
	ASSERT(VNEAR_EQUAL(summary.point, annot_local_expected, SMALL_FASTF));
    }
    ASSERT(summary.line_segment_valid);
    ASSERT(summary.line_start == 0);
    ASSERT(summary.line_end == 1);
    ASSERT(summary.text_segment_valid);
    ASSERT(summary.text_ref_point == 1);
    ASSERT(summary.text != NULL && BU_STR_EQUAL(summary.text, ctx->text));
    ASSERT(summary.cache_identity != 0);
    if (ctx->cache_identity)
	*ctx->cache_identity = summary.cache_identity;
    if (ctx->source_identity)
	*ctx->source_identity = summary.source_identity;
    ctx->found++;
    return 0;
}


static int
test_visible_annotation_view_record(struct ged *gedp,
				    const char *path,
				    const char *text,
				    uint64_t *cache_identity,
				    uint64_t *source_identity)
{
    struct test_annotation_view_record_ctx ctx =
	{path, text, cache_identity, source_identity, 0};

    ged_draw_foreach_visible_view_db_object_record(
	    test_active_view_ctx(gedp),
	    test_annotation_view_record_cb,
	    &ctx);

    return ctx.found == 1;
}


struct test_surface_view_record {
    struct ged_draw_view_surface_summary summary;
    fastf_t transparency;
    int draw_mode;
    int visible;
    int highlighted;
    int line_style;
    unsigned char color[3];
    size_t segment_count;
    size_t point_count;
    size_t index_sample_count;
    int index_sample[4];
};


struct test_line_view_record {
    struct ged_draw_view_line_summary summary;
    fastf_t transparency;
    int draw_mode;
    int visible;
    int highlighted;
    int line_style;
    unsigned char color[3];
    size_t segment_count;
    size_t point_count;
    size_t sample_count;
    point_t sample_points[TEST_LINE_VIEW_MAX_SAMPLES];
    int sample_commands[TEST_LINE_VIEW_MAX_SAMPLES];
    uint64_t cache_identity;
    uint64_t source_identity;
};


struct test_line_view_record_ctx {
    const char *path;
    int draw_mode;
    size_t sample_count;
    struct test_line_view_record *out;
    int found;
};


static int
test_line_view_record_cb(
	const struct ged_draw_view_db_object_record *rec,
	void *ud)
{
    struct test_line_view_record_ctx *ctx =
	(struct test_line_view_record_ctx *)ud;

    if (!ctx || !ctx->out || !rec || !ctx->path || !rec->path ||
	    !rec->geometry_name ||
	    !source_name_matches_drawn_prefix(rec->path, ctx->path) ||
	    !BU_STR_EQUAL(rec->geometry_name, "line-set") ||
	    (ctx->draw_mode >= 0 && rec->draw_mode != ctx->draw_mode))
	return 1;

    ASSERT(rec->is_database_source);
    ASSERT(rec->is_database_intent);
    ASSERT(!rec->non_database_source);
    ASSERT(!rec->is_view_source);
    ASSERT(!rec->is_local_source);

    memset(ctx->out, 0, sizeof(*ctx->out));
    ASSERT(ged_draw_view_db_object_record_line_summary(
	    rec, &ctx->out->summary) == 1);
    ASSERT(ctx->out->summary.valid);
    ctx->out->transparency = rec->transparency;
    ctx->out->draw_mode = rec->draw_mode;
    ctx->out->visible = rec->visible;
    ctx->out->highlighted = rec->highlighted;
    ctx->out->line_style = rec->line_style;
    ctx->out->color[0] = rec->color[0];
    ctx->out->color[1] = rec->color[1];
    ctx->out->color[2] = rec->color[2];
    ctx->out->cache_identity = rec->cache_identity;
    ctx->out->source_identity = rec->source_identity;
    ASSERT(ctx->out->summary.cache_identity == rec->cache_identity);
    ASSERT(ctx->out->summary.source_identity == rec->source_identity);

    struct segment_count_state segments = {0};
    ASSERT(ged_draw_view_db_object_record_foreach_segment(
	    rec, count_segment_cb, &segments) >= 0);
    ctx->out->segment_count = segments.count;

    struct point_count_state points = {0};
    ASSERT(ged_draw_view_db_object_record_foreach_point(
	    rec, count_point_cb, &points) >= 0);
    ctx->out->point_count = points.count;
    ASSERT(ctx->out->summary.point_count == points.count);

    ASSERT(ctx->sample_count <= ctx->out->summary.point_count);
    ctx->out->sample_count = ctx->sample_count;
    if (ctx->out->sample_count > TEST_LINE_VIEW_MAX_SAMPLES)
	ctx->out->sample_count = TEST_LINE_VIEW_MAX_SAMPLES;
    if (ctx->out->sample_count > ctx->out->summary.point_count)
	ctx->out->sample_count = ctx->out->summary.point_count;
    for (size_t i = 0; i < ctx->out->sample_count; i++) {
	ASSERT(ged_draw_view_db_object_record_line_point_at(
		rec, i, ctx->out->sample_points[i]) == 1);
	ASSERT(ged_draw_view_db_object_record_line_command_at(
		rec, i, &ctx->out->sample_commands[i]) == 1);
    }

    ctx->found++;
    return 0;
}


static int
test_visible_line_view_record(struct ged *gedp,
			      const char *path,
			      int draw_mode,
			      size_t sample_count,
			      struct test_line_view_record *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    struct test_line_view_record_ctx ctx = {path, draw_mode, sample_count, out, 0};

    ged_draw_foreach_visible_view_db_object_record(
	    test_active_view_ctx(gedp),
	    test_line_view_record_cb,
	    &ctx);

    return ctx.found == 1;
}


struct test_surface_view_record_ctx {
    const char *path;
    int draw_mode;
    size_t index_sample_count;
    struct test_surface_view_record *out;
    int found;
};


static int
test_surface_view_record_cb(
	const struct ged_draw_view_db_object_record *rec,
	void *ud)
{
    struct test_surface_view_record_ctx *ctx =
	(struct test_surface_view_record_ctx *)ud;

    if (!ctx || !ctx->out || !rec || !ctx->path || !rec->path ||
	    !rec->geometry_name ||
	    !source_name_matches_drawn_prefix(rec->path, ctx->path) ||
	    !BU_STR_EQUAL(rec->geometry_name, "indexed-face-set") ||
	    (ctx->draw_mode >= 0 && rec->draw_mode != ctx->draw_mode))
	return 1;

    ASSERT(rec->is_database_source);
    ASSERT(rec->is_database_intent);
    ASSERT(!rec->non_database_source);
    ASSERT(!rec->is_view_source);
    ASSERT(!rec->is_local_source);

    memset(ctx->out, 0, sizeof(*ctx->out));
    ASSERT(ged_draw_view_db_object_record_surface_summary(
	    rec, &ctx->out->summary));
    ASSERT(ctx->out->summary.valid);
    ctx->out->transparency = rec->transparency;
    ctx->out->draw_mode = rec->draw_mode;
    ctx->out->visible = rec->visible;
    ctx->out->highlighted = rec->highlighted;
    ctx->out->line_style = rec->line_style;
    ctx->out->color[0] = rec->color[0];
    ctx->out->color[1] = rec->color[1];
    ctx->out->color[2] = rec->color[2];

    struct segment_count_state segments = {0};
    ASSERT(ged_draw_view_db_object_record_foreach_segment(
	    rec, count_segment_cb, &segments) >= 0);
    ctx->out->segment_count = segments.count;

    struct point_count_state points = {0};
    ASSERT(ged_draw_view_db_object_record_foreach_point(
	    rec, count_point_cb, &points) >= 0);
    ctx->out->point_count = points.count;

    for (size_t i = 0; i < ctx->index_sample_count &&
	    i < sizeof(ctx->out->index_sample) / sizeof(ctx->out->index_sample[0]);
	    i++) {
	ASSERT(ged_draw_view_db_object_record_surface_index_at(
		rec, i, &ctx->out->index_sample[i]));
	ctx->out->index_sample_count++;
    }

    ctx->found++;
    return 0;
}


static int
test_visible_surface_view_record(struct ged *gedp,
				 const char *path,
				 int draw_mode,
				 size_t index_sample_count,
				 struct test_surface_view_record *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    struct test_surface_view_record_ctx ctx =
	{path, draw_mode, index_sample_count, out, 0};

    ged_draw_foreach_visible_view_db_object_record(
	    test_active_view_ctx(gedp),
	    test_surface_view_record_cb,
	    &ctx);

    return ctx.found == 1;
}


static test_scene_context *
ged_draw_shape_at(struct ged *gedp, int idx)
{
    return (test_scene_context *)ged_draw_shape_ref_context(gedp,
	    ged_draw_shape_ref_at(gedp, idx));
}


static int
ged_draw_shape_index(struct ged *gedp, test_scene_context *shape)
{
    return ged_draw_shape_ref_index(gedp,
	    ged_draw_shape_ref_from_context(gedp, shape));
}


static test_scene_context *
ged_draw_advance_shape(struct ged *gedp, test_scene_context *shape, int delta)
{
    ged_draw_shape_ref ref = ged_draw_shape_ref_from_context(gedp, shape);
    return (test_scene_context *)ged_draw_shape_ref_context(gedp,
	    ged_draw_advance_shape_ref(gedp, ref, delta));
}


static test_scene_context *
ged_draw_next_shape(struct ged *gedp, test_scene_context *shape)
{
    return ged_draw_advance_shape(gedp, shape, 1);
}


static test_scene_context *
ged_draw_group_of_shape(struct ged *gedp, test_scene_context *shape)
{
    ged_draw_shape_ref ref = ged_draw_shape_ref_from_context(gedp, shape);
    return (test_scene_context *)ged_draw_group_ref_context(gedp,
	    ged_draw_group_ref_of_shape(gedp, ref));
}


static int
ged_draw_group_dbpath(struct ged *gedp, test_scene_context *group, struct db_full_path *out)
{
    return ged_draw_group_context_dbpath(gedp, group, out);
}


static int
ged_draw_group_is_overlay(test_scene_context *group)
{
    return ged_draw_group_context_is_overlay(group);
}


static const struct db_full_path *
ged_draw_shape_fullpath(test_scene_context *shape)
{
    return ged_draw_scene_context_fullpath(shape);
}


static ged_draw_shape_ref
ged_draw_shape_ref_from_node(struct ged *gedp, test_scene_context *shape)
{
    return ged_draw_shape_ref_from_context(gedp, shape);
}


static test_scene_context *
ged_draw_shape_node_from_cache_ref(struct ged *gedp, ged_draw_shape_ref ref)
{
    return (test_scene_context *)ged_draw_shape_ref_cache_context(gedp, ref);
}


static ged_draw_shape_ref
publish_nmg_test_shape(struct ged *gedp, const struct nmgregion *r, int style)
{
    ged_draw_shape_draft *draft = NULL;
    ged_draw_group_ref group_ref = GED_DRAW_GROUP_REF_NULL;
    ged_draw_shape_ref shape_ref = GED_DRAW_SHAPE_REF_NULL;
    struct db_full_path dfp;

    db_full_path_init(&dfp);
    ASSERT(db_string_to_path(&dfp, gedp->dbip, "all.g") == 0);
    group_ref = ged_draw_group_ref_lookup_or_create(gedp, &dfp);
    db_free_full_path(&dfp);
    ASSERT(!ged_draw_group_ref_is_null(group_ref));
    if (ged_draw_group_ref_is_null(group_ref))
	return shape_ref;

    draft = ged_draw_shape_draft_create_context(gedp,
	    test_active_view_ctx(gedp), 0);
    ASSERT(draft != NULL);
    if (!draft)
	return shape_ref;

    {
	int published = ged_draw_shape_draft_publish_nmg_region(draft, r, style);
	ASSERT(published == 1);
	if (published)
	    shape_ref = ged_draw_shape_draft_commit_to_group(draft, group_ref);
	else
	    ged_draw_shape_draft_destroy(draft);
    }

    ASSERT(!ged_draw_shape_ref_is_null(shape_ref));
    return shape_ref;
}


static ged_draw_shape_ref
publish_nmg_test_shape_with_late_state(struct ged *gedp,
				       const struct nmgregion *r,
				       int style,
				       const unsigned char material_rgb[3],
				       int line_style,
				       int line_width,
				       fastf_t transparency,
				       int draw_mode)
{
    ged_draw_shape_draft *draft = NULL;
    ged_draw_group_ref group_ref = GED_DRAW_GROUP_REF_NULL;
    ged_draw_shape_ref shape_ref = GED_DRAW_SHAPE_REF_NULL;
    struct db_full_path dfp;

    db_full_path_init(&dfp);
    ASSERT(db_string_to_path(&dfp, gedp->dbip, "all.g") == 0);
    group_ref = ged_draw_group_ref_lookup_or_create(gedp, &dfp);
    ASSERT(!ged_draw_group_ref_is_null(group_ref));
    if (ged_draw_group_ref_is_null(group_ref)) {
	db_free_full_path(&dfp);
	return shape_ref;
    }

    draft = ged_draw_shape_draft_create_context(gedp,
	    test_active_view_ctx(gedp), 1);
    ASSERT(draft != NULL);
    if (!draft) {
	db_free_full_path(&dfp);
	return shape_ref;
    }

    {
	int published = ged_draw_shape_draft_publish_nmg_region(draft, r, style);
	ASSERT(published == 1);
	if (!published) {
	    ged_draw_shape_draft_destroy(draft);
	    db_free_full_path(&dfp);
	    return shape_ref;
	}
    }

    struct ged_draw_appearance_settings appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    appearance.transparency = transparency;
    appearance.draw_mode = draw_mode;
    appearance.s_line_width = line_width;
    ASSERT(ged_draw_shape_draft_apply_late_display_state(draft, &dfp,
	    line_style, &appearance, material_rgb, 1) == 1);

    shape_ref = ged_draw_shape_draft_commit_to_group(draft, group_ref);
    db_free_full_path(&dfp);

    ASSERT(!ged_draw_shape_ref_is_null(shape_ref));
    return shape_ref;
}


static int
nmg_test_vertex_uv(struct vertex *v,
		   struct vertex **verts,
		   const point_t *uvs,
		   size_t vert_count,
		   point_t uv)
{
    for (size_t i = 0; i < vert_count; i++) {
	if (v == verts[i]) {
	    VMOVE(uv, uvs[i]);
	    return 1;
	}
    }
    return 0;
}


static int
nmg_test_assign_cnurb_uv(struct vertexuse *vu,
			 struct vertex **verts,
			 const point_t *uvs,
			 size_t vert_count)
{
    point_t uv = VINIT_ZERO;

    if (!vu || !vu->v_p)
	return 0;
    if (vu->a.magic_p)
	return 1;
    if (!nmg_test_vertex_uv(vu->v_p, verts, uvs, vert_count, uv))
	return 0;
    nmg_vertexuse_a_cnurb(vu, uv);
    return 1;
}


static int
nmg_test_make_loop_edges_plinear_cnurb(struct loopuse *lu,
				       struct vertex **verts,
				       const point_t *uvs,
				       size_t vert_count)
{
    struct edgeuse *eu = NULL;

    if (!lu || BU_LIST_FIRST_MAGIC(&lu->down_hd) != NMG_EDGEUSE_MAGIC)
	return 0;

    for (BU_LIST_FOR(eu, edgeuse, &lu->down_hd)) {
	if (!nmg_test_assign_cnurb_uv(eu->vu_p, verts, uvs, vert_count) ||
		!nmg_test_assign_cnurb_uv(eu->eumate_p->vu_p, verts, uvs,
		    vert_count))
	    return 0;
    }

    for (BU_LIST_FOR(eu, edgeuse, &lu->down_hd)) {
	if (!eu->g.magic_p)
	    nmg_edge_g_cnurb_plinear(eu);
    }

    return 1;
}


static int
test_group_shape_cb_depth(test_scene_context *node, test_scene_context **out,
			  int depth)
{
    if (!node || depth > 128)
	return 1;

    size_t child_count = test_scene_child_count(node);
    if (child_count == 0) {
	*out = node;
	return 0;
    }
    for (size_t i = 0; i < child_count; i++) {
	test_scene_context *child = test_scene_child_at(node, i);
	if (child == node)
	    continue;
	if (!test_group_shape_cb_depth(child, out, depth + 1))
	    return 0;
    }
    return 1;
}


static int
test_group_shape_cb(test_scene_context *node, test_scene_context **out)
{
    return test_group_shape_cb_depth(node, out, 0);
}


static int
test_group_last_shape_cb_depth(test_scene_context *node,
			       test_scene_context **out,
			       int depth)
{
    if (!node || depth > 128)
	return 1;

    size_t child_count = test_scene_child_count(node);
    if (child_count == 0) {
	*out = node;
	return 1;
    }
    for (size_t i = 0; i < child_count; i++) {
	test_scene_context *child = test_scene_child_at(node, i);
	if (child == node)
	    continue;
	(void)test_group_last_shape_cb_depth(child, out, depth + 1);
    }
    return 1;
}


static int
test_group_last_shape_cb(test_scene_context *node, test_scene_context **out)
{
    return test_group_last_shape_cb_depth(node, out, 0);
}


static test_scene_context *
ged_draw_group_first_shape(test_scene_context *group)
{
    test_scene_context *shape = NULL;
    (void)test_group_shape_cb(group, &shape);
    return shape;
}


static test_scene_context *
ged_draw_group_last_shape(test_scene_context *group)
{
    test_scene_context *shape = NULL;
    (void)test_group_last_shape_cb(group, &shape);
    return shape;
}


static int
ged_draw_group_has_shapes(test_scene_context *group)
{
    return ged_draw_group_first_shape(group) ? 1 : 0;
}


static void
assert_shaded_surface_payload_for_source_prefix(struct ged *gedp,
						const char *source_prefix,
						int draw_mode)
{
    struct test_surface_view_record surface;

    ASSERT(test_visible_surface_view_record(gedp, source_prefix, draw_mode,
	    0, &surface));
    ASSERT(surface.summary.valid);
    ASSERT(surface.draw_mode == draw_mode);
    ASSERT(surface.summary.point_count > 0);
    ASSERT(surface.summary.normal_count > 0);
    ASSERT(surface.summary.index_count > 0);
    ASSERT(surface.summary.face_count > 0);
    ASSERT(surface.summary.normals_per_index);
    ASSERT(surface.summary.material_valid == 1);
    ASSERT(surface.summary.material_draw_mode == draw_mode);
    ASSERT(NEAR_EQUAL(surface.summary.material_transparency,
	    surface.transparency, SMALL_FASTF));
    ASSERT(surface.segment_count == 3 * surface.summary.face_count);
    ASSERT(surface.segment_count > 0);
    ASSERT(surface.point_count == surface.summary.point_count);
    ASSERT(surface.summary.cache_identity != 0);
    ASSERT(surface.summary.source_identity != 0);

    struct ged_draw_view_rendered_object_summary rendered;
    ASSERT(ged_draw_view_rendered_object_summary(test_active_view_ctx(gedp),
	    surface.summary.cache_identity, &rendered));
    ASSERT(rendered.found);
    ASSERT(rendered.is_indexed_face_set);
    ASSERT(rendered.surface_point_count == surface.summary.point_count);
    ASSERT(rendered.surface_normal_count == surface.summary.normal_count);
    ASSERT(rendered.surface_index_count == surface.summary.index_count);
    ASSERT(rendered.surface_face_count == surface.summary.face_count);
    ASSERT(rendered.surface_normals_per_index);
    ASSERT(rendered.surface_material_valid == 1);
    ASSERT(rendered.source_identity == surface.summary.source_identity);
    ASSERT(rendered.material_draw_mode == draw_mode);
    ASSERT(rendered.selection_visible == surface.visible);
}


enum draw_scene_shard {
    DRAW_SCENE_SHARD_ALL = 0,
    DRAW_SCENE_SHARD_TYPED,
    DRAW_SCENE_SHARD_CORE,
    DRAW_SCENE_SHARD_INDEX,
    DRAW_SCENE_SHARD_TREE_VIEW,
    DRAW_SCENE_SHARD_TREE,
    DRAW_SCENE_SHARD_HIGHLIGHT_MATERIAL,
    DRAW_SCENE_SHARD_VIEW_SCENE,
    DRAW_SCENE_SHARD_BOUNDS_FRAME,
    DRAW_SCENE_SHARD_HIGHLIGHT_REVISION,
    DRAW_SCENE_SHARD_RETAINED,
    DRAW_SCENE_SHARD_PATH_API,
    DRAW_SCENE_SHARD_RENAME,
    DRAW_SCENE_SHARD_REF_REMOVE,
    DRAW_SCENE_SHARD_SOURCE_UPDATE,
    DRAW_SCENE_SHARD_SELECTION_FANOUT,
    DRAW_SCENE_SHARD_COMB_UPDATE,
    DRAW_SCENE_SHARD_SCALE,
    DRAW_SCENE_SHARD_BACKEND,
    DRAW_SCENE_SHARD_UNKNOWN
};


static const char *
draw_scene_shard_name(enum draw_scene_shard shard)
{
    switch (shard) {
	case DRAW_SCENE_SHARD_ALL:
	    return "all";
	case DRAW_SCENE_SHARD_TYPED:
	    return "typed";
	case DRAW_SCENE_SHARD_CORE:
	    return "core";
	case DRAW_SCENE_SHARD_INDEX:
	    return "index";
	case DRAW_SCENE_SHARD_TREE_VIEW:
	    return "tree_view";
	case DRAW_SCENE_SHARD_TREE:
	    return "tree";
	case DRAW_SCENE_SHARD_HIGHLIGHT_MATERIAL:
	    return "highlight_material";
	case DRAW_SCENE_SHARD_VIEW_SCENE:
	    return "view_scene";
	case DRAW_SCENE_SHARD_BOUNDS_FRAME:
	    return "bounds_frame";
	case DRAW_SCENE_SHARD_HIGHLIGHT_REVISION:
	    return "highlight_revision";
	case DRAW_SCENE_SHARD_RETAINED:
	    return "retained";
	case DRAW_SCENE_SHARD_PATH_API:
	    return "path_api";
	case DRAW_SCENE_SHARD_RENAME:
	    return "rename";
	case DRAW_SCENE_SHARD_REF_REMOVE:
	    return "ref_remove";
	case DRAW_SCENE_SHARD_SOURCE_UPDATE:
	    return "source_update";
	case DRAW_SCENE_SHARD_SELECTION_FANOUT:
	    return "selection_fanout";
	case DRAW_SCENE_SHARD_COMB_UPDATE:
	    return "comb_update";
	case DRAW_SCENE_SHARD_SCALE:
	    return "scale";
	case DRAW_SCENE_SHARD_BACKEND:
	    return "backend";
	default:
	    return "unknown";
    }
}


static enum draw_scene_shard
draw_scene_shard_parse(const char *name)
{
    if (!name || !*name || BU_STR_EQUAL(name, "all"))
	return DRAW_SCENE_SHARD_ALL;
    if (BU_STR_EQUAL(name, "typed"))
	return DRAW_SCENE_SHARD_TYPED;
    if (BU_STR_EQUAL(name, "core"))
	return DRAW_SCENE_SHARD_CORE;
    if (BU_STR_EQUAL(name, "index"))
	return DRAW_SCENE_SHARD_INDEX;
    if (BU_STR_EQUAL(name, "tree_view"))
	return DRAW_SCENE_SHARD_TREE_VIEW;
    if (BU_STR_EQUAL(name, "tree"))
	return DRAW_SCENE_SHARD_TREE;
    if (BU_STR_EQUAL(name, "highlight_material"))
	return DRAW_SCENE_SHARD_HIGHLIGHT_MATERIAL;
    if (BU_STR_EQUAL(name, "view_scene"))
	return DRAW_SCENE_SHARD_VIEW_SCENE;
    if (BU_STR_EQUAL(name, "bounds_frame"))
	return DRAW_SCENE_SHARD_BOUNDS_FRAME;
    if (BU_STR_EQUAL(name, "highlight_revision"))
	return DRAW_SCENE_SHARD_HIGHLIGHT_REVISION;
    if (BU_STR_EQUAL(name, "retained"))
	return DRAW_SCENE_SHARD_RETAINED;
    if (BU_STR_EQUAL(name, "path_api"))
	return DRAW_SCENE_SHARD_PATH_API;
    if (BU_STR_EQUAL(name, "rename"))
	return DRAW_SCENE_SHARD_RENAME;
    if (BU_STR_EQUAL(name, "ref_remove"))
	return DRAW_SCENE_SHARD_REF_REMOVE;
    if (BU_STR_EQUAL(name, "source_update"))
	return DRAW_SCENE_SHARD_SOURCE_UPDATE;
    if (BU_STR_EQUAL(name, "selection_fanout"))
	return DRAW_SCENE_SHARD_SELECTION_FANOUT;
    if (BU_STR_EQUAL(name, "comb_update"))
	return DRAW_SCENE_SHARD_COMB_UPDATE;
    if (BU_STR_EQUAL(name, "scale"))
	return DRAW_SCENE_SHARD_SCALE;
    if (BU_STR_EQUAL(name, "backend"))
	return DRAW_SCENE_SHARD_BACKEND;

    return DRAW_SCENE_SHARD_UNKNOWN;
}


static int
draw_scene_shard_is_tree_family(enum draw_scene_shard shard)
{
    return shard == DRAW_SCENE_SHARD_ALL ||
	shard == DRAW_SCENE_SHARD_TREE_VIEW ||
	shard == DRAW_SCENE_SHARD_TREE ||
	shard == DRAW_SCENE_SHARD_HIGHLIGHT_MATERIAL ||
	shard == DRAW_SCENE_SHARD_VIEW_SCENE ||
	shard == DRAW_SCENE_SHARD_BOUNDS_FRAME ||
	shard == DRAW_SCENE_SHARD_HIGHLIGHT_REVISION;
}


static int
draw_scene_shard_is_retained_family(enum draw_scene_shard shard)
{
    return shard == DRAW_SCENE_SHARD_ALL ||
	shard == DRAW_SCENE_SHARD_RETAINED ||
	shard == DRAW_SCENE_SHARD_PATH_API ||
	shard == DRAW_SCENE_SHARD_RENAME ||
	shard == DRAW_SCENE_SHARD_REF_REMOVE ||
	shard == DRAW_SCENE_SHARD_SOURCE_UPDATE ||
	shard == DRAW_SCENE_SHARD_SELECTION_FANOUT ||
	shard == DRAW_SCENE_SHARD_COMB_UPDATE ||
	shard == DRAW_SCENE_SHARD_SCALE;
}


static int
draw_scene_shard_is_tree_child(enum draw_scene_shard shard)
{
    return shard == DRAW_SCENE_SHARD_TREE ||
	shard == DRAW_SCENE_SHARD_HIGHLIGHT_MATERIAL ||
	shard == DRAW_SCENE_SHARD_VIEW_SCENE ||
	shard == DRAW_SCENE_SHARD_BOUNDS_FRAME ||
	shard == DRAW_SCENE_SHARD_HIGHLIGHT_REVISION;
}


static int
draw_scene_shard_is_retained_child(enum draw_scene_shard shard)
{
    return shard == DRAW_SCENE_SHARD_PATH_API ||
	shard == DRAW_SCENE_SHARD_RENAME ||
	shard == DRAW_SCENE_SHARD_REF_REMOVE ||
	shard == DRAW_SCENE_SHARD_SOURCE_UPDATE ||
	shard == DRAW_SCENE_SHARD_SELECTION_FANOUT ||
	shard == DRAW_SCENE_SHARD_COMB_UPDATE ||
	shard == DRAW_SCENE_SHARD_SCALE;
}


static int
draw_scene_shard_runs(enum draw_scene_shard selected,
	enum draw_scene_shard candidate)
{
    if (selected == DRAW_SCENE_SHARD_ALL || selected == candidate)
	return 1;
    if (selected == DRAW_SCENE_SHARD_TREE_VIEW &&
	    draw_scene_shard_is_tree_child(candidate))
	return 1;
    if (selected == DRAW_SCENE_SHARD_RETAINED &&
	    draw_scene_shard_is_retained_child(candidate))
	return 1;
    return 0;
}


static void
draw_scene_usage(const char *prog)
{
    bu_log("Usage: %s <directory-containing-moss.g> "
	    "[all|typed|core|index|tree_view|tree|highlight_material|"
	    "view_scene|bounds_frame|highlight_revision|retained|path_api|"
	    "rename|ref_remove|source_update|selection_fanout|comb_update|"
	    "scale|backend]\n", prog);
}


int
main(int ac, char *av[])
{
    bu_setprogname(av[0]);

    if (ac < 2 || ac > 3) {
	draw_scene_usage(av[0]);
	return 1;
    }
    if (!bu_file_directory(av[1])) {
	bu_log("ERROR: [%s] is not a directory.\n", av[1]);
	return 2;
    }

    enum draw_scene_shard shard = draw_scene_shard_parse(ac == 3 ? av[2] : NULL);
    if (shard == DRAW_SCENE_SHARD_UNKNOWN) {
	draw_scene_usage(av[0]);
	return 1;
    }

    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    /* Local cache dir (mirrors the convention of basic.cpp). */
    char lcache[MAXPATHLEN] = {0};
    char cache_leaf[MAXPATHLEN] = {0};
    std::snprintf(cache_leaf, sizeof(cache_leaf), "ged_draw_scene_cache_%s",
	    draw_scene_shard_name(shard));
    bu_dir(lcache, MAXPATHLEN, BU_DIR_CURR, cache_leaf, NULL);
    bu_mkdir(lcache);
    bu_setenv("BU_DIR_CACHE", lcache, 1);

    /* Copy moss.g into a working file. */
    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s/moss.g", av[1]);
    char tmpg_name[MAXPATHLEN] = {0};
    std::snprintf(tmpg_name, sizeof(tmpg_name), "ged_draw_scene_%s.g",
	    draw_scene_shard_name(shard));
    std::ifstream orig(bu_vls_cstr(&fname), std::ios::binary);
    std::ofstream tmpg(tmpg_name, std::ios::binary);
    tmpg << orig.rdbuf();
    orig.close();
    tmpg.close();
    bu_vls_free(&fname);

    struct ged *gedp = ged_open("db", tmpg_name, 1);
    if (!gedp) {
	bu_log("ged_open failed\n");
	return 1;
    }
    bu_log("=== ged_draw_* API regression (%s) ===\n",
	    draw_scene_shard_name(shard));

    if (draw_scene_shard_runs(shard, DRAW_SCENE_SHARD_TYPED)) {
    /* ---------------------------------------------------------------- *
     * 1. NULL-arg safety: every helper must tolerate NULL gedp/path.   *
     * ---------------------------------------------------------------- */
    bu_log("[1] NULL-arg safety...\n");
    ged_draw_erase_name(NULL, "x");              /* no crash */
    ged_draw_erase_name(gedp, NULL);             /* no crash */
    ged_draw_set_highlight_state(NULL, 0);                       /* no crash */
    ged_draw_refresh_material_colors(NULL);                  /* no crash */
    ASSERT(ged_draw_scene_hash(NULL) == 0);

    /* db_full_path-keyed entry points must also be NULL-safe. */
    ASSERT(ged_draw_group_ref_is_null(ged_draw_group_ref_lookup_or_create(NULL, NULL)));
    ASSERT(ged_draw_group_ref_is_null(ged_draw_group_ref_lookup_or_create(gedp, NULL)));
    ged_draw_erase_path(NULL, NULL);              /* no crash */
    ged_draw_erase_path(gedp, NULL);              /* no crash */
    ged_draw_erase_path_prefix(NULL, NULL);            /* no crash */
    ged_draw_erase_path_prefix(gedp, NULL);            /* no crash */
    ASSERT(ged_draw_group_ref_set_dbpath(NULL, GED_DRAW_GROUP_REF_NULL, NULL) == 0);
    {
	struct db_full_path tmp;
	db_full_path_init(&tmp);
	ASSERT(ged_draw_group_dbpath(NULL, NULL, &tmp) != 0);
	ASSERT(ged_draw_group_dbpath(gedp, NULL, &tmp) != 0);
	db_free_full_path(&tmp);
    }
    {
	vect_t mn, mx;
	int empty = ged_draw_bounds(NULL, &mn, &mx, 0);
	ASSERT(empty == 1);
	empty = ged_draw_bounds(gedp, NULL, &mx, 0);
	ASSERT(empty == 1);
    }

    /* ---------------------------------------------------------------- *
     * 2. Empty state: nothing drawn yet.                                *
     * ---------------------------------------------------------------- */
    bu_log("[2] Empty draw set...\n");
    ASSERT(scene_group_count(gedp) == 0);
    ASSERT(ged_draw_scene_hash(gedp) == 0);
    {
	vect_t mn, mx;
	int empty = ged_draw_bounds(gedp, &mn, &mx, 0);
	ASSERT(empty == 1);
    }

    /* ---------------------------------------------------------------- *
     * 2a. Direct NMG vector output should publish typed line geometry.  *
     * ---------------------------------------------------------------- */
    bu_log("[2a] NMG typed line-set publisher...\n");
    {
	struct model *m = nmg_mmr();
	struct nmgregion *r = BU_LIST_FIRST(nmgregion, &m->r_hd);
	struct shell *s = nmg_msv(r);
	point_t point = {1.0, 2.0, 3.0};
	ged_draw_shape_ref shape_ref = GED_DRAW_SHAPE_REF_NULL;

	ASSERT(s != NULL);
	nmg_vertex_gv(s->vu_p->v_p, point);

	shape_ref = publish_nmg_test_shape(gedp, r, GED_DRAW_NMG_STYLE_VECTOR);

	if (!ged_draw_shape_ref_is_null(shape_ref)) {
	    struct ged_draw_view_line_summary line_summary;
	    point_t got = VINIT_ZERO;
	    int command = -1;

	    ASSERT(ged_draw_shape_ref_line_summary(gedp, shape_ref,
		    &line_summary));
	    ASSERT(line_summary.valid);
	    ASSERT(line_summary.point_count == 2);
	    ASSERT(ged_draw_shape_ref_line_point_at(gedp, shape_ref, 0, got) &&
		    VNEAR_EQUAL(got, point, SMALL_FASTF));
	    ASSERT(ged_draw_shape_ref_line_command_at(gedp, shape_ref, 0,
		    &command) && command == GED_DRAW_VIEW_LINE_MOVE);
	    ASSERT(ged_draw_shape_ref_line_point_at(gedp, shape_ref, 1, got) &&
		    VNEAR_EQUAL(got, point, SMALL_FASTF));
	    ASSERT(ged_draw_shape_ref_line_command_at(gedp, shape_ref, 1,
		    &command) && command == GED_DRAW_VIEW_LINE_DRAW);
	}

	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
	nmg_km(m);
    }

    /* ---------------------------------------------------------------- *
     * 2b. NMG polygon NO_SURFACES should publish typed outline lines.  *
     * ---------------------------------------------------------------- */
    bu_log("[2b] NMG polygon no-surfaces typed line-set publisher...\n");
    {
	struct bn_tol tol = BN_TOL_INIT_TOL;
	struct model *m = nmg_mm();
	struct nmgregion *r = nmg_mrsv(m);
	struct shell *s = BU_LIST_FIRST(shell, &r->s_hd);
	struct vertex *verts[3] = {NULL, NULL, NULL};
	struct vertex **fv[3] = {&verts[0], &verts[1], &verts[2]};
	struct faceuse *fu = nmg_cmface(s, fv, 3);
	point_t p0 = {0.0, 0.0, 0.0};
	point_t p1 = {1.0, 0.0, 0.0};
	point_t p2 = {0.0, 1.0, 0.0};
	ged_draw_shape_ref shape_ref = GED_DRAW_SHAPE_REF_NULL;

	ASSERT(fu != NULL);
	nmg_vertex_gv(verts[0], p0);
	nmg_vertex_gv(verts[1], p1);
	nmg_vertex_gv(verts[2], p2);
	ASSERT(nmg_fu_planeeqn(fu, &tol) == 0);

	shape_ref = publish_nmg_test_shape(gedp, r,
		GED_DRAW_NMG_STYLE_POLYGON | GED_DRAW_NMG_STYLE_NO_SURFACES);

	if (!ged_draw_shape_ref_is_null(shape_ref)) {
	    const point_t expected_points[4] = {
		{0.0, 0.0, 0.0},
		{1.0, 0.0, 0.0},
		{0.0, 1.0, 0.0},
		{0.0, 0.0, 0.0}
	    };
	    const int expected_commands[4] = {
		GED_DRAW_VIEW_LINE_MOVE,
		GED_DRAW_VIEW_LINE_DRAW,
		GED_DRAW_VIEW_LINE_DRAW,
		GED_DRAW_VIEW_LINE_DRAW
	    };
	    struct ged_draw_view_line_summary line_summary;

	    ASSERT(ged_draw_shape_ref_line_summary(gedp, shape_ref,
		    &line_summary));
	    ASSERT(line_summary.valid);
	    ASSERT(line_summary.point_count == 4);
	    for (size_t i = 0; i < 4; i++) {
		point_t got = VINIT_ZERO;
		int command = -1;
		ASSERT(ged_draw_shape_ref_line_point_at(gedp, shape_ref, i,
			got) && VNEAR_EQUAL(got, expected_points[i],
			    SMALL_FASTF));
		ASSERT(ged_draw_shape_ref_line_command_at(gedp, shape_ref, i,
			&command) && command == expected_commands[i]);
	    }
	}

	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);

	shape_ref = publish_nmg_test_shape(gedp, r,
		GED_DRAW_NMG_STYLE_POLYGON |
		GED_DRAW_NMG_STYLE_NO_SURFACES |
		GED_DRAW_NMG_STYLE_VISUALIZE_NORMALS);

	if (!ged_draw_shape_ref_is_null(shape_ref)) {
	    point_t centroid = {1.0 / 3.0, 1.0 / 3.0, 0.0};
	    struct ged_draw_view_line_summary line_summary;
	    point_t got = VINIT_ZERO;
	    point_t tip = VINIT_ZERO;
	    int command = -1;

	    ASSERT(ged_draw_shape_ref_line_summary(gedp, shape_ref,
		    &line_summary));
	    ASSERT(line_summary.valid);
	    ASSERT(line_summary.point_count == 6);
	    ASSERT(ged_draw_shape_ref_line_command_at(gedp, shape_ref, 4,
		    &command) && command == GED_DRAW_VIEW_LINE_MOVE);
	    ASSERT(ged_draw_shape_ref_line_point_at(gedp, shape_ref, 4, got) &&
		    VNEAR_EQUAL(got, centroid, SMALL_FASTF));
	    ASSERT(ged_draw_shape_ref_line_command_at(gedp, shape_ref, 5,
		    &command) && command == GED_DRAW_VIEW_LINE_DRAW);
	    ASSERT(ged_draw_shape_ref_line_point_at(gedp, shape_ref, 5, tip));
	    ASSERT(NEAR_EQUAL(tip[X], centroid[X], SMALL_FASTF));
	    ASSERT(NEAR_EQUAL(tip[Y], centroid[Y], SMALL_FASTF));
	    ASSERT(!NEAR_ZERO(tip[Z] - centroid[Z], SMALL_FASTF));
	}

	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
	nmg_km(m);
    }

    /* ---------------------------------------------------------------- *
     * 2c. NMG surface normals publish as auxiliary typed line geometry. *
     * ---------------------------------------------------------------- */
    bu_log("[2c] NMG surface normal auxiliary typed geometry...\n");
    {
	struct bn_tol tol = BN_TOL_INIT_TOL;
	struct model *m = nmg_mm();
	struct nmgregion *r = nmg_mrsv(m);
	struct shell *s = BU_LIST_FIRST(shell, &r->s_hd);
	struct vertex *verts[3] = {NULL, NULL, NULL};
	struct vertex **fv[3] = {&verts[0], &verts[1], &verts[2]};
	struct faceuse *fu = nmg_cmface(s, fv, 3);
	point_t p0 = {0.0, 0.0, 0.0};
	point_t p1 = {1.0, 0.0, 0.0};
	point_t p2 = {0.0, 1.0, 0.0};
	const unsigned char material_rgb[3] = {12, 34, 56};
	ged_draw_shape_ref shape_ref = GED_DRAW_SHAPE_REF_NULL;

	ASSERT(fu != NULL);
	nmg_vertex_gv(verts[0], p0);
	nmg_vertex_gv(verts[1], p1);
	nmg_vertex_gv(verts[2], p2);
	ASSERT(nmg_fu_planeeqn(fu, &tol) == 0);

	shape_ref = publish_nmg_test_shape_with_late_state(gedp, r,
		GED_DRAW_NMG_STYLE_POLYGON |
		GED_DRAW_NMG_STYLE_VISUALIZE_NORMALS,
		material_rgb, 1, 5, 0.25, GED_DRAW_MODE_SHADED);

	if (!ged_draw_shape_ref_is_null(shape_ref)) {
	    test_scene_context *shape = ged_draw_shape_node_from_cache_ref(gedp, shape_ref);
	    struct ged_draw_shape_geometry_summary face_geometry;
	    struct ged_draw_scene_display_summary face_display;
	    test_scene_context *source = test_scene_parent(shape);
	    struct ged_draw_scene_display_summary source_display;
	    test_scene_context *normal_shape = NULL;

	    ASSERT(shape != NULL);
	    ASSERT(ged_draw_shape_context_geometry_summary(shape, &face_geometry));
	    ASSERT(face_geometry.valid);
	    ASSERT(BU_STR_EQUAL(face_geometry.geometry_name, "indexed-face-set"));
	    ASSERT(ged_draw_scene_context_display_summary(shape, &face_display));
	    ASSERT(face_display.valid);
	    ASSERT(face_display.has_draw_intent);
	    ASSERT(face_display.intent_path &&
		    BU_STR_EQUAL(face_display.intent_path, "all.g"));
	    ASSERT(face_display.line_style == 1);
	    ASSERT(face_display.line_width == 5);
	    ASSERT(NEAR_EQUAL(face_display.transparency, 0.25, SMALL_FASTF));
	    ASSERT(face_display.draw_mode == GED_DRAW_MODE_SHADED);
	    ASSERT(face_geometry.point_count == 3);
	    ASSERT(face_geometry.index_count == 4);
	    ASSERT(source != NULL);
	    ASSERT(ged_draw_scene_context_display_summary(source,
		    &source_display));
	    ASSERT(source_display.valid);
	    ASSERT(source_display.is_database_source);
	    ASSERT(test_scene_child_count(source) == 2);

	    for (size_t i = 0; i < test_scene_child_count(source); i++) {
		test_scene_context *child = test_scene_child_at(source, i);
		if (child != shape &&
			BU_STR_EQUAL(test_scene_name(child), "nmg_surface_normals")) {
		    normal_shape = child;
		    break;
		}
	    }

	    ASSERT(normal_shape != NULL);
	    if (normal_shape) {
		struct ged_draw_view_line_summary line_summary;
		struct ged_draw_scene_display_summary normal_display;
		point_t centroid = {1.0 / 3.0, 1.0 / 3.0, 0.0};
		point_t got = VINIT_ZERO;
		point_t tip = VINIT_ZERO;
		int command = -1;

		ASSERT(ged_draw_scene_context_display_summary(normal_shape,
			&normal_display));
		ASSERT(normal_display.valid);
		ASSERT(normal_display.is_database_source);
		ASSERT(normal_display.has_draw_intent);
		ASSERT(normal_display.intent_path &&
			BU_STR_EQUAL(normal_display.intent_path, "all.g"));
		ASSERT(normal_display.intent_draw_mode ==
			GED_DRAW_MODE_SHADED);
		ASSERT(normal_display.visible == face_display.visible);
		ASSERT(normal_display.line_style == face_display.line_style);
		ASSERT(normal_display.line_width == face_display.line_width);
		ASSERT(NEAR_EQUAL(normal_display.transparency,
			face_display.transparency, SMALL_FASTF));
		ASSERT(normal_display.draw_mode == face_display.draw_mode);
		ASSERT(normal_display.highlighted == face_display.highlighted);
		ASSERT(normal_display.material_valid);
		ASSERT(normal_display.material_color[0] == material_rgb[0] &&
			normal_display.material_color[1] == material_rgb[1] &&
			normal_display.material_color[2] == material_rgb[2]);
		ASSERT(ged_draw_shape_context_line_summary(normal_shape,
			&line_summary));
		ASSERT(line_summary.valid);
		ASSERT(line_summary.point_count == 2);
		ASSERT(ged_draw_shape_context_line_command_at(normal_shape,
			0, &command) && command == GED_DRAW_VIEW_LINE_MOVE);
		ASSERT(ged_draw_shape_context_line_point_at(normal_shape, 0,
			got) &&
			VNEAR_EQUAL(got, centroid, SMALL_FASTF));
		ASSERT(ged_draw_shape_context_line_command_at(normal_shape,
			1, &command) && command == GED_DRAW_VIEW_LINE_DRAW);
		ASSERT(ged_draw_shape_context_line_point_at(normal_shape, 1,
			tip));
		ASSERT(NEAR_EQUAL(tip[X], centroid[X], SMALL_FASTF));
		ASSERT(NEAR_EQUAL(tip[Y], centroid[Y], SMALL_FASTF));
		ASSERT(!NEAR_ZERO(tip[Z] - centroid[Z], SMALL_FASTF));
	    }

	    {
		struct test_surface_view_record face_record;
		struct test_line_view_record normal_record;

		ASSERT(test_visible_surface_view_record(gedp, "all.g",
			GED_DRAW_MODE_SHADED, 0, &face_record));
		ASSERT(test_visible_line_view_record(gedp, "all.g",
			GED_DRAW_MODE_SHADED, 0, &normal_record));

		ASSERT(face_record.summary.valid);
		ASSERT(face_record.summary.point_count == 3);
		ASSERT(face_record.summary.normal_count == 3);
		ASSERT(face_record.summary.index_count == 4);
		ASSERT(face_record.summary.face_count == 1);
		ASSERT(face_record.summary.normals_per_index);
		ASSERT(face_record.summary.material_valid == 1);
		ASSERT(face_record.summary.material_color[0] == material_rgb[0] &&
			face_record.summary.material_color[1] == material_rgb[1] &&
			face_record.summary.material_color[2] == material_rgb[2]);
		ASSERT(NEAR_EQUAL(face_record.summary.material_transparency,
			0.25, SMALL_FASTF));
		ASSERT(face_record.summary.material_draw_mode ==
			GED_DRAW_MODE_SHADED);
		ASSERT(face_record.summary.material_highlighted == 1);
		ASSERT(face_record.color[0] == material_rgb[0] &&
			face_record.color[1] == material_rgb[1] &&
			face_record.color[2] == material_rgb[2]);
		ASSERT(face_record.line_style == 1);
		ASSERT(NEAR_EQUAL(face_record.transparency, 0.25, SMALL_FASTF));
		ASSERT(face_record.draw_mode == GED_DRAW_MODE_SHADED);
		ASSERT(face_record.highlighted == 1);
		ASSERT(face_record.segment_count == 3);
		ASSERT(face_record.point_count == 3);
		ASSERT(face_record.summary.cache_identity != 0);
		ASSERT(face_record.summary.source_identity != 0);

		ASSERT(normal_record.color[0] == material_rgb[0] &&
			normal_record.color[1] == material_rgb[1] &&
			normal_record.color[2] == material_rgb[2]);
		ASSERT(normal_record.line_style == 1);
		ASSERT(NEAR_EQUAL(normal_record.transparency, 0.25,
			SMALL_FASTF));
		ASSERT(normal_record.draw_mode == GED_DRAW_MODE_SHADED);
		ASSERT(normal_record.highlighted == 1);
		ASSERT(normal_record.segment_count == 1);
		ASSERT(normal_record.point_count == 2);
		ASSERT(normal_record.cache_identity != 0);
		ASSERT(normal_record.source_identity != 0);

		struct ged_draw_view_rendered_object_summary rendered;
		ASSERT(ged_draw_view_rendered_object_summary(test_active_view_ctx(gedp),
			face_record.summary.cache_identity, &rendered));
		ASSERT(rendered.found);
		ASSERT(rendered.is_indexed_face_set);
		ASSERT(rendered.surface_point_count ==
			face_record.summary.point_count);
		ASSERT(rendered.surface_normal_count ==
			face_record.summary.normal_count);
		ASSERT(rendered.surface_index_count ==
			face_record.summary.index_count);
		ASSERT(rendered.surface_face_count ==
			face_record.summary.face_count);
		ASSERT(rendered.surface_normals_per_index);
		ASSERT(rendered.surface_material_valid == 1);
		ASSERT(rendered.source_identity ==
			face_record.summary.source_identity);
		ASSERT(rendered.material_draw_mode ==
			GED_DRAW_MODE_SHADED);
		ASSERT(NEAR_EQUAL(rendered.material_transparency,
			0.25, SMALL_FASTF));
		ASSERT(rendered.selection_highlighted == 1);
	    }

	    ASSERT(ged_draw_shape_count(gedp) == 1);
	}

	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
	nmg_km(m);
    }

    /* ---------------------------------------------------------------- *
     * 2d. NMG mixed face/wire regions split into typed geometry.        *
     * ---------------------------------------------------------------- */
    bu_log("[2d] NMG mixed face/wire typed auxiliary geometry...\n");
    {
	struct bn_tol tol = BN_TOL_INIT_TOL;
	struct model *m = nmg_mm();
	struct nmgregion *r = nmg_mrsv(m);
	struct shell *s = BU_LIST_FIRST(shell, &r->s_hd);
	struct vertex *verts[3] = {NULL, NULL, NULL};
	struct vertex **fv[3] = {&verts[0], &verts[1], &verts[2]};
	struct faceuse *fu = nmg_cmface(s, fv, 3);
	struct edgeuse *wire_eu = NULL;
	point_t p0 = {0.0, 0.0, 0.0};
	point_t p1 = {1.0, 0.0, 0.0};
	point_t p2 = {0.0, 1.0, 0.0};
	point_t w0 = {2.0, 0.0, 0.0};
	point_t w1 = {2.0, 1.0, 0.0};
	const unsigned char material_rgb[3] = {64, 96, 128};
	ged_draw_shape_ref shape_ref = GED_DRAW_SHAPE_REF_NULL;

	ASSERT(fu != NULL);
	nmg_vertex_gv(verts[0], p0);
	nmg_vertex_gv(verts[1], p1);
	nmg_vertex_gv(verts[2], p2);
	ASSERT(nmg_fu_planeeqn(fu, &tol) == 0);
	wire_eu = nmg_me(NULL, NULL, s);
	ASSERT(wire_eu != NULL && wire_eu->eumate_p != NULL);
	if (wire_eu && wire_eu->eumate_p) {
	    nmg_vertex_gv(wire_eu->vu_p->v_p, w0);
	    nmg_vertex_gv(wire_eu->eumate_p->vu_p->v_p, w1);
	}

	shape_ref = publish_nmg_test_shape_with_late_state(gedp, r,
		GED_DRAW_NMG_STYLE_POLYGON,
		material_rgb, 1, 4, 0.5, GED_DRAW_MODE_SHADED);

	if (!ged_draw_shape_ref_is_null(shape_ref)) {
	    test_scene_context *shape = ged_draw_shape_node_from_cache_ref(gedp, shape_ref);
	    struct ged_draw_shape_geometry_summary face_geometry;
	    test_scene_context *source = test_scene_parent(shape);
	    struct ged_draw_scene_display_summary source_display;
	    test_scene_context *wire_shape = NULL;

	    ASSERT(shape != NULL);
	    ASSERT(ged_draw_shape_context_geometry_summary(shape, &face_geometry));
	    ASSERT(face_geometry.valid);
	    ASSERT(BU_STR_EQUAL(face_geometry.geometry_name, "indexed-face-set"));
	    ASSERT(face_geometry.point_count == 3);
	    ASSERT(face_geometry.index_count == 4);
	    ASSERT(source != NULL);
	    ASSERT(ged_draw_scene_context_display_summary(source,
		    &source_display));
	    ASSERT(source_display.valid);
	    ASSERT(source_display.is_database_source);
	    ASSERT(test_scene_child_count(source) == 2);

	    for (size_t i = 0; i < test_scene_child_count(source); i++) {
		test_scene_context *child = test_scene_child_at(source, i);
		if (child != shape &&
			BU_STR_EQUAL(test_scene_name(child), "nmg_mixed_wire")) {
		    wire_shape = child;
		    break;
		}
	    }

	    ASSERT(wire_shape != NULL);
	    if (wire_shape) {
		struct ged_draw_view_line_summary line_summary;
		struct ged_draw_scene_display_summary wire_display;
		point_t got = VINIT_ZERO;
		int command = -1;

		ASSERT(ged_draw_scene_context_display_summary(wire_shape,
			&wire_display));
		ASSERT(wire_display.valid);
		ASSERT(wire_display.is_database_source);
		ASSERT(wire_display.line_style == 1);
		ASSERT(wire_display.line_width == 4);
		ASSERT(NEAR_EQUAL(wire_display.transparency, 0.5,
			SMALL_FASTF));
		ASSERT(wire_display.draw_mode == GED_DRAW_MODE_SHADED);
		ASSERT(ged_draw_shape_context_line_summary(wire_shape,
			&line_summary));
		ASSERT(line_summary.valid);
		ASSERT(line_summary.point_count == 4);
		ASSERT(ged_draw_shape_context_line_command_at(wire_shape,
			0, &command) && command == GED_DRAW_VIEW_LINE_MOVE);
		ASSERT(ged_draw_shape_context_line_point_at(wire_shape, 0,
			got) &&
			VNEAR_EQUAL(got, w0, SMALL_FASTF));
		ASSERT(ged_draw_shape_context_line_command_at(wire_shape,
			1, &command) && command == GED_DRAW_VIEW_LINE_DRAW);
		ASSERT(ged_draw_shape_context_line_point_at(wire_shape, 1,
			got) &&
			VNEAR_EQUAL(got, w1, SMALL_FASTF));
		ASSERT(ged_draw_shape_context_line_command_at(wire_shape,
			2, &command) && command == GED_DRAW_VIEW_LINE_MOVE);
		ASSERT(ged_draw_shape_context_line_point_at(wire_shape, 2,
			got) &&
			VNEAR_EQUAL(got, w1, SMALL_FASTF));
		ASSERT(ged_draw_shape_context_line_command_at(wire_shape,
			3, &command) && command == GED_DRAW_VIEW_LINE_DRAW);
		ASSERT(ged_draw_shape_context_line_point_at(wire_shape, 3,
			got) &&
			VNEAR_EQUAL(got, w0, SMALL_FASTF));
	    }

	    ASSERT(ged_draw_shape_count(gedp) == 1);
	}

	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
	nmg_km(m);
    }

    /* ---------------------------------------------------------------- *
     * 2e. NMG snurb/cnurb outlines publish as typed line geometry.     *
     * ---------------------------------------------------------------- */
    bu_log("[2e] NMG snurb/cnurb typed outline geometry...\n");
    {
	struct model *m = nmg_mm();
	struct nmgregion *r = nmg_mrsv(m);
	struct shell *s = BU_LIST_FIRST(shell, &r->s_hd);
	struct vertex *verts[4] = {NULL, NULL, NULL, NULL};
	struct vertex **fv[4] = {&verts[0], &verts[1], &verts[2], &verts[3]};
	struct faceuse *fu = nmg_cmface(s, fv, 4);
	struct loopuse *lu = NULL;
	point_t p0 = {0.0, 0.0, 0.0};
	point_t p1 = {1.0, 0.0, 0.0};
	point_t p2 = {1.0, 1.0, 0.5};
	point_t p3 = {0.0, 1.0, 0.0};
	point_t uvs[4] = {
	    {0.0, 0.0, 1.0},
	    {1.0, 0.0, 1.0},
	    {1.0, 1.0, 1.0},
	    {0.0, 1.0, 1.0}
	};
	fastf_t *ukv = (fastf_t *)bu_malloc(4 * sizeof(fastf_t),
		"test snurb u knots");
	fastf_t *vkv = (fastf_t *)bu_malloc(4 * sizeof(fastf_t),
		"test snurb v knots");
	fastf_t *mesh = (fastf_t *)bu_malloc(12 * sizeof(fastf_t),
		"test snurb mesh");
	ged_draw_shape_ref shape_ref = GED_DRAW_SHAPE_REF_NULL;

	ASSERT(fu != NULL);
	nmg_vertex_gv(verts[0], p0);
	nmg_vertex_gv(verts[1], p1);
	nmg_vertex_gv(verts[2], p2);
	nmg_vertex_gv(verts[3], p3);
	ukv[0] = 0.0; ukv[1] = 0.0; ukv[2] = 1.0; ukv[3] = 1.0;
	vkv[0] = 0.0; vkv[1] = 0.0; vkv[2] = 1.0; vkv[3] = 1.0;
	VMOVE(&mesh[0], p0);
	VMOVE(&mesh[3], p1);
	VMOVE(&mesh[6], p3);
	VMOVE(&mesh[9], p2);
	nmg_face_g_snurb(fu, 2, 2, 4, 4, ukv, vkv, 2, 2,
		RT_NURB_MAKE_PT_TYPE(3, RT_NURB_PT_XYZ, RT_NURB_PT_NONRAT),
		mesh);

	lu = BU_LIST_FIRST(loopuse, &fu->lu_hd);
	ASSERT(lu != NULL);
	ASSERT(nmg_test_make_loop_edges_plinear_cnurb(lu, verts, uvs, 4) == 1);

	shape_ref = publish_nmg_test_shape(gedp, r,
		GED_DRAW_NMG_STYLE_POLYGON | GED_DRAW_NMG_STYLE_NO_SURFACES);

	if (!ged_draw_shape_ref_is_null(shape_ref)) {
	    const size_t expected_count = 4 * (10 + 2) - 3;
	    struct ged_draw_view_line_summary line_summary;
	    point_t got = VINIT_ZERO;
	    int command = -1;
	    int saw_surface_sample = 0;

	    ASSERT(ged_draw_shape_ref_line_summary(gedp, shape_ref,
		    &line_summary));
	    ASSERT(line_summary.valid);
	    ASSERT(line_summary.point_count == expected_count);
	    ASSERT(ged_draw_shape_ref_line_command_at(gedp, shape_ref, 0,
		    &command) && command == GED_DRAW_VIEW_LINE_MOVE);
	    ASSERT(ged_draw_shape_ref_line_point_at(gedp, shape_ref, 0, got) &&
		    VNEAR_EQUAL(got, p0, SMALL_FASTF));
	    ASSERT(ged_draw_shape_ref_line_point_at(gedp, shape_ref,
		    expected_count - 1, got) && VNEAR_EQUAL(got, p0,
			SMALL_FASTF));
	    for (size_t i = 1; i < expected_count; i++) {
		ASSERT(ged_draw_shape_ref_line_command_at(gedp, shape_ref, i,
			&command) && command == GED_DRAW_VIEW_LINE_DRAW);
		if (ged_draw_shape_ref_line_point_at(gedp, shape_ref, i, got) &&
			got[Z] > 0.01 && got[Z] < 0.5)
		    saw_surface_sample = 1;
	    }
	    ASSERT(saw_surface_sample == 1);
	    ASSERT(ged_draw_shape_count(gedp) == 1);
	}

	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
	nmg_km(m);
    }

    /* ---------------------------------------------------------------- *
     * 2f. SUBMODEL draw publishes typed realized child geometry.        *
     * ---------------------------------------------------------------- */
    bu_log("[2f] SUBMODEL typed child-shape publisher...\n");
    {
	struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
	const char *submodel_name = "_ged_test_submodel.s";
	const char *draw_av[3] = {"draw", submodel_name, NULL};

	ASSERT(wdbp != RT_WDB_NULL);
	ASSERT(db_lookup(gedp->dbip, "all.g", LOOKUP_QUIET) != RT_DIR_NULL);
	if (wdbp)
	    ASSERT(mk_submodel(wdbp, submodel_name, NULL, "all.g", 0) == 0);
	ASSERT(db_lookup(gedp->dbip, submodel_name, LOOKUP_QUIET) != RT_DIR_NULL);

	ASSERT(ged_exec(gedp, 2, draw_av) == BRLCAD_OK);

	test_scene_context *submodel_shape =
	    test_find_shape_by_dbpath(test_scene_root(gedp), submodel_name);
	ASSERT(submodel_shape != NULL);
	if (submodel_shape) {
	    struct ged_draw_scene_display_summary submodel_display;
	    test_scene_context *submodel_source = test_scene_parent(submodel_shape);
	    struct ged_draw_scene_display_summary submodel_source_display;
	    test_scene_context *leaf_source = NULL;
	    test_scene_context *leaf_shape = NULL;
	    size_t submodel_source_child_count = 0;

	    ASSERT(ged_draw_scene_context_display_summary(submodel_shape,
		    &submodel_display));
	    ASSERT(submodel_display.valid);
	    ASSERT(submodel_display.is_database_source);
	    ASSERT(submodel_source != NULL);
	    ASSERT(ged_draw_scene_context_display_summary(submodel_source,
		    &submodel_source_display));
	    ASSERT(submodel_source_display.valid);
	    ASSERT(submodel_source_display.is_database_source);
	    submodel_source_child_count = test_scene_child_count(submodel_source);
	    ASSERT(submodel_source_child_count > 1);

	    for (size_t i = 0; i < submodel_source_child_count; i++) {
		test_scene_context *child = test_scene_child_at(submodel_source, i);
		struct ged_draw_scene_display_summary candidate_source_display;
		struct ged_draw_shape_geometry_summary direct_geometry_summary;
		int have_direct_geometry_summary =
		    ged_draw_shape_context_geometry_summary(child,
			    &direct_geometry_summary);
		int direct_typed_geometry =
		    have_direct_geometry_summary &&
		    (BU_STR_EQUAL(direct_geometry_summary.geometry_name, "line-set") ||
		     BU_STR_EQUAL(direct_geometry_summary.geometry_name, "point-set") ||
		     BU_STR_EQUAL(direct_geometry_summary.geometry_name, "indexed-line-set") ||
		     BU_STR_EQUAL(direct_geometry_summary.geometry_name, "indexed-face-set"));
		if (ged_draw_shape_context_has_state(child) &&
			direct_typed_geometry &&
			direct_geometry_summary.point_count > 0) {
		    leaf_source = submodel_source;
		    leaf_shape = child;
		    break;
		}

		if (child == submodel_shape ||
			!ged_draw_scene_context_display_summary(child,
				&candidate_source_display) ||
			!candidate_source_display.is_database_source)
		    continue;

		if (!candidate_source_display.has_draw_intent ||
			!candidate_source_display.intent_path ||
			!source_name_matches_drawn_prefix(
			    candidate_source_display.intent_path, "all.g"))
		    continue;

		for (size_t j = 0; j < test_scene_child_count(child); j++) {
		    test_scene_context *candidate_shape = test_scene_child_at(child, j);
		    struct ged_draw_shape_geometry_summary geometry_summary;
		    int have_geometry_summary =
			ged_draw_shape_context_geometry_summary(candidate_shape,
				&geometry_summary);
		    int typed_geometry =
			have_geometry_summary &&
			(BU_STR_EQUAL(geometry_summary.geometry_name, "line-set") ||
			 BU_STR_EQUAL(geometry_summary.geometry_name, "point-set") ||
			 BU_STR_EQUAL(geometry_summary.geometry_name, "indexed-line-set") ||
			 BU_STR_EQUAL(geometry_summary.geometry_name, "indexed-face-set"));

		    if (ged_draw_shape_context_has_state(candidate_shape) &&
			    typed_geometry &&
			    geometry_summary.point_count > 0) {
			leaf_source = child;
			leaf_shape = candidate_shape;
			break;
		    }
		}
		if (leaf_shape)
		    break;
	    }

	    ASSERT(leaf_source != NULL);
	    ASSERT(leaf_shape != NULL);
	    if (leaf_source && leaf_shape) {
		struct ged_draw_scene_display_summary leaf_source_display;
		struct ged_draw_scene_display_summary leaf_shape_display;
		vect_t mn, mx;

		ASSERT(test_scene_child_count(leaf_source) >= 1);
		ASSERT(ged_draw_scene_context_display_summary(leaf_source,
			&leaf_source_display));
		ASSERT(leaf_source_display.valid);
		ASSERT(leaf_source_display.has_draw_intent);
		ASSERT(leaf_source_display.intent_path &&
			(source_name_matches_drawn_prefix(
			    leaf_source_display.intent_path, "all.g") ||
			 source_name_matches_drawn_prefix(
			    leaf_source_display.intent_path, submodel_name)));

		ASSERT(ged_draw_shape_context_has_state(leaf_shape));
		ASSERT(ged_draw_shape_context_source(leaf_shape) != NULL);
		ASSERT(ged_draw_scene_context_display_summary(leaf_shape,
			&leaf_shape_display));
		ASSERT(leaf_shape_display.valid);
		ASSERT(leaf_shape_display.is_database_source);
		ASSERT(leaf_shape_display.has_draw_intent);
		ASSERT(leaf_shape_display.intent_path &&
			(source_name_matches_drawn_prefix(
			    leaf_shape_display.intent_path, "all.g") ||
			 source_name_matches_drawn_prefix(
			    leaf_shape_display.intent_path, submodel_name)));
		ASSERT(ged_draw_bounds(gedp, &mn, &mx, 0) == 0);
	    }

	    struct ged_draw_transaction redraw_op =
		ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, submodel_name);
	    ASSERT(ged_draw_apply_transaction(gedp, &redraw_op, NULL) > 0);
	    submodel_shape =
		test_find_shape_by_dbpath(test_scene_root(gedp), submodel_name);
	    ASSERT(submodel_shape != NULL);
	    if (submodel_shape) {
		submodel_source = test_scene_parent(submodel_shape);
		ASSERT(submodel_source != NULL);
		ASSERT(test_scene_child_count(submodel_source) ==
			submodel_source_child_count);
	    }
	}

	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);

	{
	    struct bu_vls external_file = BU_VLS_INIT_ZERO;
	    const char *external_submodel_name = "_ged_test_ext_submodel.s";
	    const char *external_draw_av[3] = {"draw", external_submodel_name, NULL};
	    bu_vls_sprintf(&external_file, "%s/moss.g", av[1]);

	    if (wdbp)
		ASSERT(mk_submodel(wdbp, external_submodel_name,
			    bu_vls_cstr(&external_file), "all.g", 0) == 0);
	    ASSERT(db_lookup(gedp->dbip, external_submodel_name,
		    LOOKUP_QUIET) != RT_DIR_NULL);
	    ASSERT(ged_exec(gedp, 2, external_draw_av) == BRLCAD_OK);

	    test_scene_context *external_shape =
		test_find_shape_by_dbpath(test_scene_root(gedp),
			external_submodel_name);
	    ASSERT(external_shape != NULL);
	    if (external_shape) {
		test_scene_context *external_source = test_scene_parent(external_shape);
		struct ged_draw_scene_display_summary external_source_display;
		size_t nested_sources = 0;

		ASSERT(external_source != NULL);
		ASSERT(ged_draw_scene_context_display_summary(external_source,
			&external_source_display));
		ASSERT(external_source_display.valid);
		ASSERT(external_source_display.is_database_source);
		ASSERT(test_scene_child_count(external_source) > 1);
		for (size_t i = 0; i < test_scene_child_count(external_source);
			i++) {
		    test_scene_context *child = test_scene_child_at(external_source, i);
		    struct ged_draw_scene_display_summary child_display;
		    if (child != external_shape &&
			    ged_draw_scene_context_display_summary(child,
				&child_display) &&
			    child_display.is_database_source &&
			    test_scene_child_count(child) > 0)
			nested_sources++;
		}
		ASSERT(nested_sources > 0);
	    }

	    ged_draw_clear(gedp);
	    ASSERT(scene_group_count(gedp) == 0);
	    bu_vls_free(&external_file);
	}
    }

    /* ---------------------------------------------------------------- *
     * 2g. ANNOT draw publishes a typed annotation record.               *
     * ---------------------------------------------------------------- */
    bu_log("[2g] ANNOT typed annotation publisher...\n");
    {
	struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
	const char *annot_name = "_ged_test_annot.s";
	const char *draw_av[3] = {"draw", annot_name, NULL};
	struct rt_annot_internal ann;
	memset(&ann, 0, sizeof(ann));
	ann.magic = RT_ANNOT_INTERNAL_MAGIC;
	VSET(ann.V, 0.25, 0.5, 0.0);
	ann.vert_count = 2;
	ann.verts = (point2d_t *)bu_calloc(2, sizeof(point2d_t),
		"test annot verts");
	V2SET(ann.verts[0], 0.0, 0.0);
	V2SET(ann.verts[1], 0.2, 0.3);
	ann.ant.count = 2;
	ann.ant.reverse = (int *)bu_calloc(2, sizeof(int),
		"test annot reverse");
	ann.ant.segments = (void **)bu_calloc(2, sizeof(void *),
		"test annot segments");
	struct line_seg *lsg;
	BU_ALLOC(lsg, struct line_seg);
	lsg->magic = CURVE_LSEG_MAGIC;
	lsg->start = 0;
	lsg->end = 1;
	ann.ant.segments[0] = (void *)lsg;
	struct txt_seg *tsg;
	BU_ALLOC(tsg, struct txt_seg);
	tsg->magic = ANN_TSEG_MAGIC;
	tsg->ref_pt = 1;
	tsg->rel_pos = RT_TXT_POS_BL;
	bu_vls_init(&tsg->label);
	bu_vls_strcpy(&tsg->label, "typed note");
	tsg->txt_size = 12.0;
	tsg->txt_rot_angle = 15.0;
	ann.ant.segments[1] = (void *)tsg;

	ASSERT(wdbp != RT_WDB_NULL);
	if (wdbp)
	    ASSERT(mk_annot(wdbp, annot_name, &ann) == 0);
	bu_vls_free(&tsg->label);
	BU_PUT(tsg, struct txt_seg);
	BU_PUT(lsg, struct line_seg);
	bu_free(ann.ant.segments, "test annot segments");
	bu_free(ann.ant.reverse, "test annot reverse");
	bu_free(ann.verts, "test annot verts");
	ASSERT(db_lookup(gedp->dbip, annot_name, LOOKUP_QUIET) != RT_DIR_NULL);

	ASSERT(ged_exec(gedp, 2, draw_av) == BRLCAD_OK);
	test_scene_context *annot_shape =
	    test_find_shape_by_dbpath(test_scene_root(gedp), annot_name);
	struct ged_draw_shape_geometry_summary annot_geometry;
	ASSERT(annot_shape != NULL);
	if (annot_shape) {
	    ASSERT(ged_draw_shape_context_geometry_summary(annot_shape,
		    &annot_geometry));
	    ASSERT(annot_geometry.valid);
	    ASSERT(BU_STR_EQUAL(annot_geometry.geometry_name, "annotation"));

	    uint64_t annot_cache_identity = 0;
	    uint64_t annot_source_identity = 0;
	    ASSERT(test_visible_annotation_view_record(gedp, annot_name,
		    "typed note", &annot_cache_identity,
		    &annot_source_identity));
	    ASSERT(annot_cache_identity != 0);
	    ASSERT(annot_source_identity != 0);

	    struct ged_draw_view_rendered_object_summary rendered;
	    ASSERT(ged_draw_view_rendered_object_summary(test_active_view_ctx(gedp),
		    annot_cache_identity, &rendered));
	    ASSERT(rendered.found);
	    ASSERT(rendered.is_annotation);
	    ASSERT(rendered.annotation_segment_count == 2);
	    ASSERT(rendered.source_identity == annot_source_identity);

	    struct ged_draw_transaction redraw_op =
		ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, annot_name);
	    ASSERT(ged_draw_apply_transaction(gedp, &redraw_op, NULL) > 0);
	    annot_shape =
		test_find_shape_by_dbpath(test_scene_root(gedp), annot_name);
	    ASSERT(annot_shape != NULL);
	    if (annot_shape) {
		ASSERT(ged_draw_shape_context_geometry_summary(annot_shape,
			&annot_geometry));
		ASSERT(annot_geometry.valid);
		ASSERT(BU_STR_EQUAL(annot_geometry.geometry_name,
			"annotation"));
	    }

	    annot_cache_identity = 0;
	    annot_source_identity = 0;
	    ASSERT(test_visible_annotation_view_record(gedp, annot_name,
		    "typed note", &annot_cache_identity,
		    &annot_source_identity));
	    ASSERT(annot_cache_identity != 0);
	    ASSERT(annot_source_identity != 0);
	}

	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
	ASSERT(test_find_shape_by_dbpath(test_scene_root(gedp), annot_name) == NULL);
	{
	    const char *fallback_av[4] = {"draw", "-m2", annot_name, NULL};
	    ASSERT(ged_exec(gedp, 3, fallback_av) == BRLCAD_OK);
	    annot_shape =
		test_find_shape_by_dbpath(test_scene_root(gedp), annot_name);
	    ASSERT(annot_shape != NULL);
	    if (annot_shape) {
		ASSERT(ged_draw_shape_context_geometry_summary(annot_shape,
			&annot_geometry));
		ASSERT(annot_geometry.valid);
		ASSERT(BU_STR_EQUAL(annot_geometry.geometry_name,
			"annotation"));
	    }
	}

	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
	{
	    const char *strict_av[5] = {"draw", "-m2", "--strict", annot_name, NULL};
	    ASSERT(ged_exec(gedp, 4, strict_av) == BRLCAD_OK);
	    annot_shape =
		test_find_shape_by_dbpath(test_scene_root(gedp), annot_name);
	    ASSERT(annot_shape != NULL);
	    if (annot_shape) {
		ASSERT(ged_draw_shape_context_geometry_summary(annot_shape,
			&annot_geometry));
		ASSERT(annot_geometry.valid);
		ASSERT(BU_STR_EQUAL(annot_geometry.geometry_name,
			"annotation"));
	    }
	}

	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
    }

    /* ---------------------------------------------------------------- *
     * 2h. Shaded BoT draw publishes a typed surface payload.            *
     * ---------------------------------------------------------------- */
    bu_log("[2h] shaded BoT command path typed surface payload...\n");
    {
	const char *view_av[5] = {"view", "lod", "mesh", "0", NULL};
	const char *facetize_av[5] = {"facetize", "-r", "all.g", "all.bot", NULL};
	const char *draw_av[4] = {"draw", "-m1", "all.bot", NULL};

	ASSERT(ged_exec(gedp, 4, view_av) == BRLCAD_OK);
	ASSERT(ged_exec(gedp, 4, facetize_av) == BRLCAD_OK);
	ASSERT(db_lookup(gedp->dbip, "all.bot", LOOKUP_QUIET) != RT_DIR_NULL);
	ASSERT(ged_exec(gedp, 3, draw_av) == BRLCAD_OK);
	ASSERT(scene_group_count(gedp) > 0);
	assert_shaded_surface_payload_for_source_prefix(gedp, "all.bot",
		GED_DRAW_MODE_SHADED_BOTS);
	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
    }

    /* ---------------------------------------------------------------- *
     * 2h1. Shaded BoT direct face-set honors CW winding.                *
     * ---------------------------------------------------------------- */
    bu_log("[2h1] shaded BoT CW winding typed surface payload...\n");
    {
	const char *cw_bot = "_ged_draw_cw_winding.bot";
	fastf_t vertices[9] = {
	    0.0, 0.0, 0.0,
	    1.0, 0.0, 0.0,
	    0.0, 1.0, 0.0
	};
	int faces[3] = {0, 1, 2};
	const char *draw_av[4] = {"draw", "-m1", cw_bot, NULL};
	struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);

	ASSERT(wdbp != RT_WDB_NULL);
	if (wdbp)
	    ASSERT(mk_bot(wdbp, cw_bot, RT_BOT_SURFACE, RT_BOT_CW, 0,
			3, 1, vertices, faces, NULL, NULL) == 0);
	ASSERT(db_lookup(gedp->dbip, cw_bot, LOOKUP_QUIET) != RT_DIR_NULL);
	ASSERT(ged_exec(gedp, 3, draw_av) == BRLCAD_OK);
	{
	    struct test_surface_view_record surface;
	    ASSERT(test_visible_surface_view_record(gedp, cw_bot,
		    GED_DRAW_MODE_SHADED_BOTS, 4, &surface));
	    ASSERT(surface.summary.valid);
	    ASSERT(surface.summary.index_count == 4);
	    ASSERT(surface.index_sample_count == 4);
	    if (surface.index_sample_count == 4) {
		ASSERT(surface.index_sample[0] == 0);
		ASSERT(surface.index_sample[1] == 2);
		ASSERT(surface.index_sample[2] == 1);
		ASSERT(surface.index_sample[3] == -1);
	    }
	}
	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
    }

    /* ---------------------------------------------------------------- *
     * 2h2. Direct face-set failure falls back unless strict is set.     *
     * ---------------------------------------------------------------- */
    bu_log("[2h2] shaded direct face-set failure fallback...\n");
    {
	const char *fallback_av[4] = {"draw", "-m1", "all.bot", NULL};
	const char *strict_av[5] = {"draw", "-m1", "--strict", "all.bot", NULL};

	ged_draw_test_force_primitive_face_set_failure(1);
	ASSERT(ged_exec(gedp, 3, fallback_av) == BRLCAD_OK);
	ASSERT(scene_group_count(gedp) > 0);
	assert_shaded_surface_payload_for_source_prefix(gedp, "all.bot",
		GED_DRAW_MODE_SHADED_BOTS);
	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);

	ASSERT(ged_exec(gedp, 4, strict_av) == BRLCAD_OK);
	ASSERT(test_visible_db_view_record_geometry_count(gedp, "all.bot",
		"indexed-face-set") == 0);
	ged_draw_test_force_primitive_face_set_failure(0);
	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
    }

    /* ---------------------------------------------------------------- *
     * 2i. Shaded BREP/CDT draw publishes a typed surface payload.       *
     * ---------------------------------------------------------------- */
    bu_log("[2i] shaded BREP/CDT command path typed surface payload...\n");
    {
	const char *tol_av[4] = {"tol", "rel", "0.01", NULL};
	const char *brep_av[5] = {"brep", "all.g", "brep", "all.brep", NULL};
	const char *draw_av[4] = {"draw", "-m1", "all.brep", NULL};

	ASSERT(ged_exec(gedp, 3, tol_av) == BRLCAD_OK);
	ASSERT(ged_exec(gedp, 4, brep_av) == BRLCAD_OK);
	ASSERT(db_lookup(gedp->dbip, "all.brep", LOOKUP_QUIET) != RT_DIR_NULL);
	ASSERT(ged_exec(gedp, 3, draw_av) == BRLCAD_OK);
	ASSERT(scene_group_count(gedp) > 0);
	assert_shaded_surface_payload_for_source_prefix(gedp, "all.brep",
		GED_DRAW_MODE_SHADED_BOTS);
	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
    }
    }

    if (draw_scene_shard_runs(shard, DRAW_SCENE_SHARD_CORE)) {

    /* ---------------------------------------------------------------- *
     * 3. Draw via ged_exec, then probe API.                             *
     * ---------------------------------------------------------------- */
    bu_log("[3] Draw and probe...\n");
    {
	const char *s_av[3] = {"draw", "all.g", NULL};
	ged_exec(gedp, 2, s_av);
    }
    int after_draw = scene_group_count(gedp);
    ASSERT(after_draw > 0);
    ASSERT(ged_draw_path_state(NULL, NULL, "all.g", -1) == 0);
    ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), NULL, -1) == 0);
    ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "definitely_no_such_obj", -1) == 0);
    ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 1);
    ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "/all.g/", -1) == 1);
    ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g/platform.r", -1) == 1);
    ASSERT(ged_draw_has_paths(gedp, test_active_view_ctx(gedp), -1) == 1);
    {
	struct bu_vls group_paths = BU_VLS_INIT_ZERO;
	struct bu_vls leaf_paths = BU_VLS_INIT_ZERO;
	ASSERT(ged_draw_list_paths(gedp, test_active_view_ctx(gedp), -1, 0, &group_paths) > 0);
	ASSERT(vls_has_line(&group_paths, "all.g"));
	ASSERT(ged_draw_list_paths(gedp, test_active_view_ctx(gedp), -1, 1, &leaf_paths) > 0);
	ASSERT(vls_has_line(&leaf_paths, "all.g/platform.r/platform.s"));
	bu_vls_free(&group_paths);
	bu_vls_free(&leaf_paths);
    }
    {
	struct db_full_path dfp;
	db_full_path_init(&dfp);
	ASSERT(db_string_to_path(&dfp, gedp->dbip, "all.g/box.r") == 0);
	ged_draw_erase_path(gedp, &dfp);
	db_free_full_path(&dfp);
    }
    ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 2);
    ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g/box.r", -1) == 0);
    ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g/platform.r", -1) == 1);
    {
	struct bu_vls group_paths = BU_VLS_INIT_ZERO;
	struct bu_vls leaf_paths = BU_VLS_INIT_ZERO;
	ASSERT(ged_draw_list_paths(gedp, test_active_view_ctx(gedp), -1, 0, &group_paths) > 0);
	ASSERT(!vls_has_line(&group_paths, "all.g"));
	ASSERT(vls_has_line(&group_paths, "all.g/platform.r"));
	ASSERT(ged_draw_list_paths(gedp, test_active_view_ctx(gedp), -1, 1, &leaf_paths) > 0);
	ASSERT(!vls_has_line(&leaf_paths, "all.g/box.r/box.s"));
	ASSERT(vls_has_line(&leaf_paths, "all.g/platform.r/platform.s"));
	bu_vls_free(&group_paths);
	bu_vls_free(&leaf_paths);
    }
    {
	const char *who_av[2] = {"who", NULL};
	const char *who_real_av[3] = {"who", "real", NULL};
	const char *who_invalid_av[3] = {"who", "not-a-who-mode", NULL};

	ASSERT(ged_exec_who(gedp, 1, who_av) == BRLCAD_OK);
	ASSERT(!vls_has_line(gedp->ged_result_str, "all.g"));
	ASSERT(vls_has_line(gedp->ged_result_str, "all.g/platform.r"));

	ASSERT(ged_exec_who(gedp, 2, who_real_av) == BRLCAD_OK);
	ASSERT(!vls_has_word(gedp->ged_result_str, "all.g"));
	ASSERT(vls_has_word(gedp->ged_result_str, "all.g/platform.r"));

	ASSERT(ged_exec_who(gedp, 2, who_invalid_av) == BRLCAD_ERROR);
    }
    ged_draw_clear(gedp);
    {
	const char *child_av[3] = {"draw", "all.g/box.r", NULL};
	ASSERT(ged_exec(gedp, 2, child_av) == BRLCAD_OK);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g/box.r", -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 2);

	struct ged_draw_transaction txn =
	    ged_draw_transaction_make(GED_DRAW_TXN_ERASE, "all.g/box.r");
	txn.view = test_active_view_ctx(gedp);
	ASSERT(ged_draw_apply_transaction(gedp, &txn, NULL) > 0);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g/box.r", -1) == 0);
    }
    ged_draw_clear(gedp);
    {
	const char *s_av[3] = {"draw", "all.g", NULL};
	ged_exec(gedp, 2, s_av);
    }

    /* lookup_or_add_dbpath on an already-drawn path must return non-NULL
     * and must not insert a duplicate. */
    int before_lookup = scene_group_count(gedp);
    ged_draw_group_ref h = GED_DRAW_GROUP_REF_NULL;
    {
	struct db_full_path dfp;
	db_full_path_init(&dfp);
	if (db_string_to_path(&dfp, gedp->dbip, "all.g") == 0)
	    h = ged_draw_group_ref_lookup_or_create(gedp, &dfp);
	db_free_full_path(&dfp);
    }
    ASSERT(!ged_draw_group_ref_is_null(h));
    ASSERT(scene_group_count(gedp) == before_lookup);

    /* lookup_or_add_dbpath on a non-existent leaf must return NULL.  We
     * test the path-string variant for the legacy noisy-log fallback
     * via a pragma-guarded call below in section [15]. */
    ged_draw_group_ref h_missing = GED_DRAW_GROUP_REF_NULL;
    {
	struct db_full_path dfp;
	db_full_path_init(&dfp);
	if (db_string_to_path(&dfp, gedp->dbip, "definitely_no_such_obj") == 0)
	    h_missing = ged_draw_group_ref_lookup_or_create(gedp, &dfp);
	db_free_full_path(&dfp);
    }
    ASSERT(ged_draw_group_ref_is_null(h_missing));

    /* bounds must report non-empty after a draw. */
    {
	vect_t mn, mx;
	int empty = ged_draw_bounds(gedp, &mn, &mx, 0);
	ASSERT(empty == 0);
	/* Sanity: min <= max on each axis. */
	ASSERT(mn[X] <= mx[X]);
	ASSERT(mn[Y] <= mx[Y]);
	ASSERT(mn[Z] <= mx[Z]);
    }

    /* scene_hash must be non-zero with content drawn. */
    unsigned long long h1 = ged_draw_scene_hash(gedp);
    ASSERT(h1 != 0);

    /* highlight flag must propagate to every drawn scene obj. */
    ged_draw_set_highlight_state(gedp, 1);
    {
	int all_up = 1;
	auto cb = +[](const struct ged_draw_shape_record *rec, void *ud) -> int {
	    int *ok = (int *)ud;
	    if (!rec || !rec->highlighted) *ok = 0;
	    return 1;
	};
	ged_draw_foreach_shape_record(gedp, cb, &all_up);
	ASSERT(all_up);
    }
    ged_draw_set_highlight_state(gedp, 0);
    {
	int all_down = 1;
	auto cb = +[](const struct ged_draw_shape_record *rec, void *ud) -> int {
	    int *ok = (int *)ud;
	    if (rec && rec->highlighted) *ok = 0;
	    return 1;
	};
	ged_draw_foreach_shape_record(gedp, cb, &all_down);
	ASSERT(all_down);
    }

    /* material color refresh must run cleanly. */
    ged_draw_refresh_material_colors(gedp);

    /* ---------------------------------------------------------------- *
     * 4. insert an overlay shape.                           *
     * ---------------------------------------------------------------- */
    bu_log("[4] overlay insert...\n");
    {
	point_t p1 = {0, 0, 0}, p2 = {100, 100, 100};
	point_t points[2];
	int commands[2] = {BG_GEOMETRY_LINE_MOVE, BG_GEOMETRY_LINE_DRAW};
	VMOVE(points[0], p1);
	VMOVE(points[1], p2);
	struct ged_draw_overlay_geometry geometry = {};
	geometry.kind = GED_DRAW_OVERLAY_GEOMETRY_LINE_SET;
	geometry.points = (const point_t *)points;
	geometry.point_count = 2;
	geometry.commands = commands;
	geometry.command_count = 2;

	/* Save current count so we can verify a new entry was added.
	 * The overlay insert creates the _overlays group as a new root child. */
	int before_overlay = scene_group_count(gedp);
	ged_draw_shape_ref overlay_ref = GED_DRAW_SHAPE_REF_NULL;
	int rc = ged_draw_overlay_geometry_insert(gedp,
				     "_ged_test_overlay", &geometry,
				     0xFF8800,  /* orange */
				     1.0, 0, 0, &overlay_ref);
	ASSERT(rc == 0);
	ASSERT(!ged_draw_shape_ref_is_null(overlay_ref));
	ASSERT(scene_group_count(gedp) > before_overlay);

	/* The _overlays group must be present as a root child and be identified as overlay storage. */
	{
	    test_scene_context *root = test_scene_root(gedp);
	    test_scene_context *overlays_grp = NULL;
	    for (size_t i = 0; i < test_scene_child_count(root); i++) {
		test_scene_context *g = test_scene_child_at(root, i);
		if (BU_STR_EQUAL("_overlays", test_scene_name(g))) {
		    overlays_grp = g;
		    break;
		}
	    }
	    ASSERT(overlays_grp != NULL);
	    ASSERT(ged_draw_group_is_overlay(overlays_grp));

	    /* The overlay shape must be retained under overlay storage. */
	    if (overlays_grp && test_scene_child_count(overlays_grp) > 0) {
		test_scene_context *sp = test_scene_child_at(overlays_grp, 0);
		ASSERT(sp != NULL);
		/* No database entry should exist for this name. */
		ASSERT(db_lookup(gedp->dbip, "_ged_test_overlay", LOOKUP_QUIET)
		       == RT_DIR_NULL);
	    }
	}

	/* Erase the overlay shape by name. */
	ged_draw_erase_name(gedp, "_ged_test_overlay");
	/* _overlays group should be gone (empty → freed). */
	ASSERT(scene_group_count(gedp) == before_overlay);
    }

    /* ---------------------------------------------------------------- *
     * 4b. NIRT/qray uplot import publishes typed view line layers.      *
     * ---------------------------------------------------------------- */
    bu_log("[4b] NIRT/qray uplot feature import...\n");
    {
	const char *qray_name = "_ged_test_qray_uplot";
	struct ged_uplot_stream *stream =
	    _ged_uplot_stream_create(1.0, PL_OUTPUT_MODE_TEXT);
	ASSERT(stream != NULL);
	if (stream) {
	    ASSERT(test_uplot_stream_process_args(stream, 'C',
		    "255 0 0") == BRLCAD_OK);
	    ASSERT(test_uplot_stream_process_args(stream, 'L',
		    "0 0 0 100 0 0") == BRLCAD_OK);
	    ASSERT(test_uplot_stream_process_args(stream, 'C',
		    "0 0 255") == BRLCAD_OK);
	    ASSERT(test_uplot_stream_process_args(stream, 'L',
		    "0 0 0 0 100 0") == BRLCAD_OK);
	    ASSERT(_ged_uplot_stream_publish_feature(gedp, stream,
		    qray_name) == BRLCAD_OK);

	    ASSERT(db_lookup(gedp->dbip, qray_name, LOOKUP_QUIET) ==
		    RT_DIR_NULL);

	    struct ged_draw_view_feature_summary feature_summary =
		GED_DRAW_VIEW_FEATURE_SUMMARY_INIT;
	    ASSERT(ged_draw_view_context_feature_summary(test_active_view_ctx(gedp), qray_name,
		    &feature_summary));
	    ASSERT(feature_summary.exists);
	    ASSERT(feature_summary.is_overlay);
	    ASSERT(feature_summary.child_count == 2);
	    ASSERT(feature_summary.geometry_command_count == 4);

	    ASSERT(test_qray_view_records(gedp, qray_name));

	    ASSERT(ged_draw_view_context_feature_remove(test_active_view_ctx(gedp),
		    qray_name) == 1);
	    ASSERT(ged_draw_view_context_feature_summary(test_active_view_ctx(gedp), qray_name,
		    &feature_summary));
	    ASSERT(!feature_summary.exists);
	    ASSERT(!test_qray_view_records(gedp, qray_name));
	    _ged_uplot_stream_free(stream);
	}
    }

    /* ---------------------------------------------------------------- *
     * 5. erase_by_path / erase_all_paths semantics.                     *
     * ---------------------------------------------------------------- */
    bu_log("[5] erase_*...\n");
    {
	/* Make sure we still have all.g drawn. */
	const char *s_av[3] = {"draw", "all.g", NULL};
	ged_exec(gedp, 2, s_av);

	int before = scene_group_count(gedp);
	ASSERT(before > 0);

	/* erase_by_dbpath on the exact drawn name must remove that entry. */
	{
	    struct db_full_path dfp;
	    db_full_path_init(&dfp);
	    if (db_string_to_path(&dfp, gedp->dbip, "all.g") == 0)
		ged_draw_erase_path(gedp, &dfp);
	    db_free_full_path(&dfp);
	}
	ASSERT(scene_group_count(gedp) < before);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 0);

	/* Re-draw and try erase_all_dbpaths. */
	ged_exec(gedp, 2, s_av);
	int before2 = scene_group_count(gedp);
	ASSERT(before2 > 0);
	{
	    struct db_full_path dfp;
	    db_full_path_init(&dfp);
	    if (db_string_to_path(&dfp, gedp->dbip, "all.g") == 0)
		ged_draw_erase_path_prefix(gedp, &dfp);
	    db_free_full_path(&dfp);
	}
	/* Note: erase_all_dbpaths matches subset paths, so should clear
	 * everything that has all.g as a prefix component. */
	ASSERT(scene_group_count(gedp) <= before2);
    }
    {
	const char *draw_av[3] = {"draw", "all.g", NULL};
	const char *erase_av[3] = {"erase", "all.g", NULL};
	const char *erase_recursive_av[3] = {"erase", "-r", "all.g"};
	const char *erase_attr_av[4] = {"erase", "-oA", "ged_draw_attr_test", "yes"};
	const char *zap_av[2] = {"zap", NULL};
	const char *autoview_av[2] = {"autoview", NULL};
	const char *lod_av[2] = {"lod", NULL};
	const char *draw_hidden_av[3] = {"draw", "-h", "all.g"};
	const char *draw_attr_av[4] = {"draw", "-A", "ged_draw_attr_test", "yes"};
	ASSERT(ged_exec(gedp, 2, draw_av) == BRLCAD_OK);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 1);
	ASSERT(ged_draw_has_paths(gedp, test_active_view_ctx(gedp), -1) == 1);
	ASSERT(ged_exec_autoview(gedp, 1, autoview_av) == BRLCAD_OK);
	ASSERT(ged_exec(gedp, 1, lod_av) == BRLCAD_OK);
	ASSERT(ged_exec_erase(gedp, 2, erase_av) == BRLCAD_OK);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 0);
	ASSERT(ged_draw_has_paths(gedp, test_active_view_ctx(gedp), -1) == 0);
	ASSERT(ged_exec(gedp, 2, draw_av) == BRLCAD_OK);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 1);
	ASSERT(ged_exec_erase(gedp, 3, erase_recursive_av) == BRLCAD_OK);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 0);
	ASSERT(ged_exec(gedp, 2, draw_av) == BRLCAD_OK);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 1);
	ASSERT(ged_exec_zap(gedp, 1, zap_av) == BRLCAD_OK);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 0);
	ASSERT(ged_draw_has_paths(gedp, test_active_view_ctx(gedp), -1) == 0);
	ASSERT(ged_exec(gedp, 3, draw_hidden_av) == BRLCAD_OK);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", GED_DRAW_MODE_HIDDEN_LINE) == 1);
	ASSERT(ged_exec_zap(gedp, 1, zap_av) == BRLCAD_OK);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 0);
	{
	    struct directory *attr_dp = db_lookup(gedp->dbip, "platform.r", LOOKUP_QUIET);
	    ASSERT(attr_dp != RT_DIR_NULL);
	    if (attr_dp != RT_DIR_NULL) {
		struct bu_attribute_value_set avs = BU_AVS_INIT_ZERO;
		ASSERT(db5_get_attributes(gedp->dbip, &avs, attr_dp) == 0);
		bu_avs_add(&avs, "ged_draw_attr_test", "yes");
		ASSERT(db5_update_attributes(attr_dp, &avs, gedp->dbip) == 0);
		ASSERT(ged_exec(gedp, 4, draw_attr_av) == BRLCAD_OK);
		ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "platform.r", -1) == 1);
		ASSERT(ged_exec_erase(gedp, 4, erase_attr_av) == BRLCAD_OK);
		ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "platform.r", -1) == 0);
	    }
	}
    }

    /* ---------------------------------------------------------------- *
     * 5b. command-level operation intents.                             *
     * ---------------------------------------------------------------- */
    bu_log("[5b] draw operation intent dispatch...\n");
    {
	const char *s_av[3] = {"draw", "all.g", NULL};
	ged_exec(gedp, 2, s_av);
	ASSERT(scene_group_count(gedp) > 0);

	struct ged_draw_transaction op =
	    ged_draw_transaction_make(GED_DRAW_TXN_ERASE, "all.g");
	struct ged_draw_transaction_result result;
	ged_draw_transaction_result_init(&result);
	ASSERT(ged_draw_apply_transaction(gedp, &op, &result) == 1);
	ASSERT(result.status == 1);
	ASSERT(result.affected_groups == 1);
	ASSERT(result.affected_shapes == 1);
	ASSERT(result.scene_revision_after > result.scene_revision_before);
	ASSERT(strstr(bu_vls_cstr(&result.names), "all.g") != NULL);
	ged_draw_transaction_result_free(&result);
	ASSERT(scene_group_count(gedp) == 0);

	ged_exec(gedp, 2, s_av);
	ASSERT(scene_group_count(gedp) > 0);
	op = ged_draw_transaction_make(GED_DRAW_TXN_ERASE_PREFIX, "all.g");
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 1);
	ASSERT(scene_group_count(gedp) == 0);

	ged_exec(gedp, 2, s_av);
	ASSERT(scene_group_count(gedp) > 0);
	op = ged_draw_transaction_make(GED_DRAW_TXN_CLEAR, NULL);
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 1);
	ASSERT(scene_group_count(gedp) == 0);

	ged_exec(gedp, 2, s_av);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 1);
	op = ged_draw_transaction_make(GED_DRAW_TXN_CLEAR_SCOPE, NULL);
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 0);
	ASSERT(scene_group_count(gedp) == 0);

	op = ged_draw_transaction_make(GED_DRAW_TXN_NONE, NULL);
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 0);
	ASSERT(ged_draw_apply_transaction(NULL, &op, NULL) == 0);
	ASSERT(ged_draw_apply_transaction(gedp, NULL, NULL) == 0);

	const char *draw_paths[3] = {"all.g", "box.r", NULL};
	struct ged_draw_appearance_settings draw_vs =
	    GED_DRAW_APPEARANCE_SETTINGS_INIT;
	draw_vs.draw_mode = GED_DRAW_MODE_HIDDEN_LINE;
	draw_vs.s_line_width = 5;
	draw_vs.color_override = 1;
	draw_vs.color[0] = 10;
	draw_vs.color[1] = 20;
	draw_vs.color[2] = 30;
	op = ged_draw_transaction_make(GED_DRAW_TXN_DRAW, NULL);
	op.view = test_active_view_ctx(gedp);
	op.paths = draw_paths;
	op.path_count = 2;
	op.appearance = &draw_vs;
	ged_draw_transaction_result_init(&result);
	ASSERT(ged_draw_apply_transaction(gedp, &op, &result) == 2);
	ASSERT(result.status == 2);
	ASSERT(result.affected_groups == 2);
	ASSERT(result.affected_shapes > 0);
	ASSERT(result.scene_revision_after > result.scene_revision_before);
	ASSERT(strstr(bu_vls_cstr(&result.names), "all.g") != NULL);
	ASSERT(strstr(bu_vls_cstr(&result.names), "box.r") != NULL);
	ged_draw_transaction_result_free(&result);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", GED_DRAW_MODE_HIDDEN_LINE) == 1);
	{
	    ged_draw_shape_ref first_ref = ged_draw_first_shape_ref(gedp);
	    struct ged_draw_shape_record rec;
	    ASSERT(!ged_draw_shape_ref_is_null(first_ref));
	    ASSERT(ged_draw_shape_record_get(gedp, first_ref, &rec) == 1);
	    ASSERT(rec.draw_mode == GED_DRAW_MODE_HIDDEN_LINE);
	    ASSERT(rec.line_width == 5);
	}
	op = ged_draw_transaction_make(GED_DRAW_TXN_CLEAR, NULL);
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 1);

	const char *child_draw_paths[3] = {"all.g/box.r", "all.g/platform.r", NULL};
	op = ged_draw_transaction_make(GED_DRAW_TXN_DRAW, NULL);
	op.view = test_active_view_ctx(gedp);
	op.paths = child_draw_paths;
	op.path_count = 2;
	ged_draw_transaction_result_init(&result);
	ASSERT(ged_draw_apply_transaction(gedp, &op, &result) == 2);
	ASSERT(result.status == 2);
	ASSERT(result.affected_groups == 2);
	ASSERT(result.affected_shapes > 0);
	ASSERT(result.scene_revision_after > result.scene_revision_before);
	ASSERT(strstr(bu_vls_cstr(&result.names), "all.g/box.r") != NULL);
	ASSERT(strstr(bu_vls_cstr(&result.names), "all.g/platform.r") != NULL);
	ged_draw_transaction_result_free(&result);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g/box.r", -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g/platform.r", -1) == 1);
	op = ged_draw_transaction_make(GED_DRAW_TXN_CLEAR, NULL);
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 1);

	ged_draw_clear(gedp);
	struct draw_observer_test_state obs;
	memset(&obs, 0, sizeof(obs));
	obs.gedp = gedp;
	BU_VLS_INIT(&obs.names);
	ged_draw_observer_token token =
	    ged_draw_observer_add(gedp, draw_observer_test_cb, &obs);
	ASSERT(token != 0);

	op = ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "all.g");
	op.view = test_active_view_ctx(gedp);
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 1);
	ASSERT(obs.calls == 1);
	ASSERT(obs.last_kind == GED_DRAW_TXN_DRAW);
	ASSERT(obs.last_status == 1);
	ASSERT(obs.last_after > obs.last_before);
	ASSERT(strstr(bu_vls_cstr(&obs.names), "all.g") != NULL);

	int calls_before_noop = obs.calls;
	op = ged_draw_transaction_make(GED_DRAW_TXN_ERASE, "definitely_no_such_object");
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 0);
	ASSERT(obs.calls == calls_before_noop);

	op = ged_draw_transaction_make(GED_DRAW_TXN_ERASE, "all.g");
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 1);
	ASSERT(obs.calls == calls_before_noop + 1);
	ASSERT(obs.last_kind == GED_DRAW_TXN_ERASE);
	ASSERT(obs.last_status == 1);

	struct draw_observer_test_state self_obs;
	memset(&self_obs, 0, sizeof(self_obs));
	self_obs.gedp = gedp;
	BU_VLS_INIT(&self_obs.names);
	self_obs.self_token =
	    ged_draw_observer_add(gedp, draw_observer_self_remove_cb, &self_obs);
	ASSERT(self_obs.self_token != 0);

	op = ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "all.g");
	op.view = test_active_view_ctx(gedp);
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 1);
	ASSERT(self_obs.calls == 1);
	op = ged_draw_transaction_make(GED_DRAW_TXN_ERASE, "all.g");
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 1);
	ASSERT(self_obs.calls == 1);

	int calls_before_remove = obs.calls;
	ASSERT(ged_draw_observer_remove(gedp, token) == 1);
	op = ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "all.g");
	op.view = test_active_view_ctx(gedp);
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 1);
	ASSERT(obs.calls == calls_before_remove);
	ASSERT(ged_draw_observer_remove(gedp, token) == 0);
	bu_vls_free(&self_obs.names);
	bu_vls_free(&obs.names);
	ged_draw_clear(gedp);

	ASSERT(ged_event_txn_available(gedp) == 1);
	ged_exec(gedp, 2, s_av);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 1);

	struct event_observer_test_state internal_obs;
	struct event_observer_test_state post_obs;
	memset(&internal_obs, 0, sizeof(internal_obs));
	memset(&post_obs, 0, sizeof(post_obs));
	BU_VLS_INIT(&internal_obs.names);
	BU_VLS_INIT(&post_obs.names);

	ged_event_observer_token internal_token =
	    ged_event_observer_add(gedp, GED_EVENT_OBSERVER_INTERNAL,
		    event_observer_test_cb, &internal_obs);
	ged_event_observer_token post_token =
	    ged_event_observer_add(gedp, GED_EVENT_OBSERVER_POST_RECONCILE,
		    event_observer_test_cb, &post_obs);
	ASSERT(internal_token != 0);
	ASSERT(post_token != 0);

	struct ged_event_txn_result eresult;
	ged_event_txn_result_init(&eresult);
	ASSERT(ged_event_batch_begin(gedp) == 1);
	ASSERT(ged_event_notify_object_modified(gedp, "all.g", 1, NULL) == 0);
	ASSERT(ged_event_notify_object_modified(gedp, "all.g", 1, NULL) == 0);
	ASSERT(ged_event_batch_end(gedp, &eresult) > 0);
	ASSERT(eresult.event_count == 2);
	ASSERT(eresult.coalesced_event_count == 1);
	ASSERT(eresult.draw_scene_revision_after > eresult.draw_scene_revision_before);
	ASSERT(strstr(bu_vls_cstr(&eresult.affected_names), "all.g") != NULL);
	ASSERT(internal_obs.calls == 1);
	ASSERT(post_obs.calls == 1);
	ASSERT(internal_obs.last_kind == GED_EVENT_OBJECT_MODIFIED);
	ASSERT(post_obs.last_kind == GED_EVENT_OBJECT_MODIFIED);
	ASSERT(post_obs.last_coalesced_count == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 1);
	ged_event_txn_result_free(&eresult);

	ged_event_txn_result_init(&eresult);
	int post_calls_before = post_obs.calls;
	ASSERT(ged_event_batch_begin(gedp) == 1);
	ASSERT(ged_event_notify_object_modified(gedp, "all.g", 0, NULL) == 0);
	ASSERT(ged_event_notify_object_modified(gedp, "all.g", 1, NULL) == 0);
	ASSERT(ged_event_batch_end(gedp, &eresult) > 0);
	ASSERT(eresult.event_count == 2);
	ASSERT(eresult.coalesced_event_count == 1);
	ASSERT(post_obs.calls == post_calls_before + 1);
	ASSERT(post_obs.last_kind == GED_EVENT_OBJECT_MODIFIED);
	ASSERT(post_obs.last_redraw == 1);
	ged_event_txn_result_free(&eresult);

	ged_event_txn_result_init(&eresult);
	post_calls_before = post_obs.calls;
	ASSERT(ged_event_notify_object_added(gedp, "all.g", &eresult) >= 0);
	ASSERT(post_obs.calls == post_calls_before + 1);
	ASSERT(post_obs.last_kind == GED_EVENT_OBJECT_ADDED);
	ASSERT(eresult.coalesced_event_count == 1);
	ASSERT(strstr(bu_vls_cstr(&eresult.affected_names), "all.g") != NULL);
	ged_event_txn_result_free(&eresult);

	ged_event_txn_result_init(&eresult);
	post_calls_before = post_obs.calls;
	ASSERT(ged_event_notify_comb_tree_changed(gedp, "all.g", 1, &eresult) >= 0);
	ASSERT(post_obs.calls == post_calls_before + 1);
	ASSERT(post_obs.last_kind == GED_EVENT_COMB_TREE_CHANGED);
	ASSERT(eresult.coalesced_event_count == 1);
	ASSERT(strstr(bu_vls_cstr(&eresult.affected_names), "all.g") != NULL);
	ged_event_txn_result_free(&eresult);

	ged_event_txn_result_init(&eresult);
	post_calls_before = post_obs.calls;
	ASSERT(ged_event_notify_attribute_changed(gedp, "all.g", 1, &eresult) >= 0);
	ASSERT(post_obs.calls == post_calls_before + 1);
	ASSERT(post_obs.last_kind == GED_EVENT_ATTRIBUTE_CHANGED);
	ASSERT(eresult.coalesced_event_count == 1);
	ASSERT(strstr(bu_vls_cstr(&eresult.affected_names), "all.g") != NULL);
	ged_event_txn_result_free(&eresult);

	ged_event_txn_result_init(&eresult);
	post_calls_before = post_obs.calls;
	ASSERT(ged_event_notify_batch_rebuild(gedp, &eresult) >= 0);
	ASSERT(post_obs.calls == post_calls_before + 1);
	ASSERT(post_obs.last_kind == GED_EVENT_BATCH_REBUILD);
	ASSERT(eresult.coalesced_event_count == 1);
	ged_event_txn_result_free(&eresult);

	struct event_observer_followup_state follow_obs;
	memset(&follow_obs, 0, sizeof(follow_obs));
	follow_obs.self_token =
	    ged_event_observer_add(gedp, GED_EVENT_OBSERVER_POST_RECONCILE,
		    event_observer_followup_cb, &follow_obs);
	ASSERT(follow_obs.self_token != 0);
	post_calls_before = post_obs.calls;
	ASSERT(ged_event_notify_object_added(gedp, "all.g", NULL) >= 0);
	ASSERT(follow_obs.calls == 2);
	ASSERT(follow_obs.kinds[0] == GED_EVENT_OBJECT_ADDED);
	ASSERT(follow_obs.kinds[1] == GED_EVENT_ATTRIBUTE_CHANGED);
	ASSERT(post_obs.calls == post_calls_before + 2);
	ASSERT(post_obs.last_kind == GED_EVENT_ATTRIBUTE_CHANGED);
	ASSERT(ged_event_observer_remove(gedp, follow_obs.self_token) == 1);

	struct event_observer_test_state event_self_obs;
	memset(&event_self_obs, 0, sizeof(event_self_obs));
	BU_VLS_INIT(&event_self_obs.names);
	event_self_obs.self_token =
	    ged_event_observer_add(gedp, GED_EVENT_OBSERVER_POST_RECONCILE,
		    event_observer_self_remove_cb, &event_self_obs);
	ASSERT(event_self_obs.self_token != 0);
	ASSERT(ged_event_notify_object_modified(gedp, "all.g", 0, NULL) >= 0);
	ASSERT(event_self_obs.calls == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 1);
	ASSERT(ged_event_notify_object_modified(gedp, "all.g", 0, NULL) >= 0);
	ASSERT(event_self_obs.calls == 1);

	ASSERT(ged_event_observer_remove(gedp, internal_token) == 1);
	ASSERT(ged_event_observer_remove(gedp, post_token) == 1);
	ASSERT(ged_event_observer_remove(gedp, event_self_obs.self_token) == 0);
	bu_vls_free(&event_self_obs.names);
	bu_vls_free(&post_obs.names);
	bu_vls_free(&internal_obs.names);
	ged_draw_clear(gedp);
    }

    /* ---------------------------------------------------------------- *
     * 5c. property-update operation intents.                            *
     * ---------------------------------------------------------------- */
    bu_log("[5c] draw property operation intents...\n");
    {
	const char *s_av[3] = {"draw", "all.g", NULL};
	ged_exec(gedp, 2, s_av);

	struct ged_draw_shape_record first_rec;
	ASSERT(test_shape_record_at(gedp, 0, &first_rec));
	ASSERT(ZERO(first_rec.transparency - 1.0));

	struct ged_draw_transaction op =
	    ged_draw_transaction_make_value(GED_DRAW_TXN_TRANSPARENCY,
					  "all.g", 0.35);
	struct ged_draw_transaction_result result;
	ged_draw_transaction_result_init(&result);
	ASSERT(ged_draw_apply_transaction(gedp, &op, &result) > 0);
	ASSERT(result.affected_shapes > 0);
	ASSERT(strstr(bu_vls_cstr(&result.names), "all.g") != NULL);
	ged_draw_transaction_result_free(&result);
	ASSERT(test_shape_record_at(gedp, 0, &first_rec));
	ASSERT(ZERO(first_rec.transparency - 0.35));
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 0);

	const char *tr_av[4] = {"set_transparency", "all.g", "0.65", NULL};
	ged_exec(gedp, 3, tr_av);
	ASSERT(test_shape_record_at(gedp, 0, &first_rec));
	ASSERT(ZERO(first_rec.transparency - 0.65));

	op = ged_draw_transaction_make_value(GED_DRAW_TXN_TRANSPARENCY,
					   NULL, 0.2);
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 0);

	op = ged_draw_transaction_make_value(GED_DRAW_TXN_TRANSPARENCY,
					   "no_such_object", 0.2);
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 0);
    }

    /* ---------------------------------------------------------------- *
     * 5d. default render-style and material operation intents.          *
     * ---------------------------------------------------------------- */
    bu_log("[5d] draw default/material operation intents...\n");
    {
	ASSERT(ged_draw_default_mode(gedp) == GED_DRAW_MODE_WIRE);

	struct ged_draw_transaction op =
	    ged_draw_transaction_make_value(GED_DRAW_TXN_DEFAULT_DRAW_MODE,
					  NULL, (double)GED_DRAW_MODE_SHADED);
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 1);
	ASSERT(ged_draw_default_mode(gedp) == GED_DRAW_MODE_SHADED);
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 0);

	const char *sm_av[3] = {"shaded_mode", "1", NULL};
	ged_exec(gedp, 2, sm_av);
	ASSERT(ged_draw_default_mode(gedp) == GED_DRAW_MODE_SHADED_BOTS);

	uint64_t rev0 = ged_draw_material_revision(gedp);
	op = ged_draw_transaction_make(GED_DRAW_TXN_MATERIAL_CHANGED, NULL);
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 1);
	uint64_t rev1 = ged_draw_material_revision(gedp);
	ASSERT(rev1 > rev0);

	op = ged_draw_transaction_make(GED_DRAW_TXN_REFRESH_MATERIAL_COLORS, NULL);
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 1);
	ASSERT(ged_draw_material_revision(gedp) == rev1);

	op = ged_draw_transaction_make_value(GED_DRAW_TXN_DEFAULT_DRAW_MODE,
					   NULL, (double)GED_DRAW_MODE_WIRE);
	ged_draw_apply_transaction(gedp, &op, NULL);
	ASSERT(ged_draw_default_mode(gedp) == GED_DRAW_MODE_WIRE);
    }

    /* ---------------------------------------------------------------- *
     * 5e. graph visibility operation intents.                           *
     * ---------------------------------------------------------------- */
    bu_log("[5e] draw visibility operation intents...\n");
    {
	struct ged_draw_shape_record first_rec;
	struct ged_draw_group_record group_rec;
	ASSERT(test_first_shape_group_records(gedp, &first_rec, &group_rec));
	ASSERT(group_rec.visible);
	ASSERT(first_rec.visible);

	struct ged_draw_transaction op =
	    ged_draw_transaction_make_value(GED_DRAW_TXN_VISIBILITY,
					  "all.g", 0.0);
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) > 0);
	ASSERT(test_first_shape_group_records(gedp, &first_rec, &group_rec));
	ASSERT(!group_rec.visible);
	ASSERT(!first_rec.visible);

	op = ged_draw_transaction_make_value(GED_DRAW_TXN_VISIBILITY,
					   "all.g", 1.0);
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) > 0);
	ASSERT(test_first_shape_group_records(gedp, &first_rec, &group_rec));
	ASSERT(group_rec.visible);
	ASSERT(first_rec.visible);

	const char *leaf = first_rec.leaf_name;
	ASSERT(leaf != NULL);
	if (leaf) {
	    op = ged_draw_transaction_make_value(GED_DRAW_TXN_VISIBILITY,
					       leaf, 0.0);
	    ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) > 0);
	    ASSERT(test_shape_record_at(gedp, 0, &first_rec));
	    ASSERT(!first_rec.visible);

	    op = ged_draw_transaction_make_value(GED_DRAW_TXN_VISIBILITY,
					       leaf, 1.0);
	    ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) > 0);
	    ASSERT(test_shape_record_at(gedp, 0, &first_rec));
	    ASSERT(first_rec.visible);
	}

	const char *hide_av[3] = {"hide", "all.g", NULL};
	ged_exec(gedp, 2, hide_av);
	ASSERT(test_first_shape_group_records(gedp, &first_rec, &group_rec));
	ASSERT(!group_rec.visible);
	const char *unhide_av[3] = {"unhide", "all.g", NULL};
	ged_exec(gedp, 2, unhide_av);
	ASSERT(test_first_shape_group_records(gedp, &first_rec, &group_rec));
	ASSERT(group_rec.visible);
	ASSERT(first_rec.visible);
    }

    /* ---------------------------------------------------------------- *
     * 5f. redraw operation intent.                                      *
     * ---------------------------------------------------------------- */
    bu_log("[5f] draw redraw operation intent...\n");
    {
	struct ged_draw_transaction op =
	    ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, "all.g");
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) > 0);

	const char *redraw_av[3] = {"redraw", "all.g", NULL};
	ASSERT(ged_exec(gedp, 2, redraw_av) == BRLCAD_OK);

	op = ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, "no_such_object");
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) < 0);
    }

    /* ---------------------------------------------------------------- *
     * 5g. highlight operation intent.                                   *
     * ---------------------------------------------------------------- */
    bu_log("[5g] draw highlight operation intent...\n");
    {
	struct ged_draw_shape_record first_rec;
	ASSERT(test_shape_record_at(gedp, 0, &first_rec));
	const char *leaf = first_rec.leaf_name;
	ASSERT(leaf != NULL);
	if (leaf) {
	    struct ged_draw_transaction op =
		ged_draw_transaction_make_value(GED_DRAW_TXN_HIGHLIGHT,
					      leaf, 1.0);
	    ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) > 0);
	    ASSERT(test_shape_record_at(gedp, 0, &first_rec));
	    ASSERT(first_rec.highlighted);

	    const char *illum_av[4] = {"illum", "-n", leaf, NULL};
	    ASSERT(ged_exec(gedp, 3, illum_av) == BRLCAD_OK);
	    ASSERT(test_shape_record_at(gedp, 0, &first_rec));
	    ASSERT(!first_rec.highlighted);
	}

	struct ged_draw_transaction op =
	    ged_draw_transaction_make_value(GED_DRAW_TXN_HIGHLIGHT,
					  "no_such_object", 1.0);
	ASSERT(ged_draw_apply_transaction(gedp, &op, NULL) == 0);
    }

    /* ---------------------------------------------------------------- *
     * 5h. GED path-hash selection bridges to render-item selection.     *
     * ---------------------------------------------------------------- */
    bu_log("[5h] GED selection render bridge...\n");
    {
	test_scene_context *first = ged_draw_first_shape(gedp);
	ASSERT(first != NULL);
	const struct db_full_path *fp = ged_draw_shape_fullpath(first);
	ASSERT(fp != NULL);
	char *path = fp ? db_path_to_string(fp) : NULL;
	ASSERT(path != NULL);
	if (path) {
	    const char *clear_av[3] = {"select", "clear", NULL};
	    ASSERT(ged_exec_select(gedp, 2, clear_av) == BRLCAD_OK);
	    ged_draw_shape_ref first_ref = ged_draw_shape_ref_from_node(gedp, first);
	    ASSERT(!ged_draw_shape_ref_is_null(first_ref));
	    ASSERT(ged_draw_shape_set_highlighted(gedp, first_ref, 0) == 1);
	    const char *select_av[4] = {"select", "add", path, NULL};
	    ASSERT(ged_exec_select(gedp, 3, select_av) == BRLCAD_OK);
	    first_ref = ged_draw_shape_ref_from_node(gedp, first);
	    ASSERT(!ged_draw_shape_ref_is_null(first_ref));
	    struct ged_draw_shape_record first_selected_record;
	    ASSERT(ged_draw_shape_record_get(gedp, first_ref, &first_selected_record) == 1);
	    ASSERT(first_selected_record.selected);
	    ASSERT(test_selected_visible_view_record(gedp, path));

	    ASSERT(ged_exec_select(gedp, 2, clear_av) == BRLCAD_OK);
	    bu_free(path, "selection bridge path");
	}
    }
    }

    if (draw_scene_shard_runs(shard, DRAW_SCENE_SHARD_INDEX)) {

    /* ---------------------------------------------------------------- *
     * 6. zap then verify hash returns to zero.                          *
     * ---------------------------------------------------------------- */
    bu_log("[6] zap and rehash...\n");
    {
	const char *s_av[2] = {"zap", NULL};
	ged_exec(gedp, 1, s_av);
	ASSERT(ged_draw_scene_hash(gedp) == 0);
	ASSERT(scene_group_count(gedp) == 0);
    }

    /* ---------------------------------------------------------------- *
     * 7. shape_count / shape_at / shape_index — snapshotted DFS index. *
     * ---------------------------------------------------------------- */
    bu_log("[7] shape_count/at/index...\n");
    {
	/* Draw moss scene. */
	const char *s_av[3] = {"draw", "all.g", NULL};
	ged_exec(gedp, 2, s_av);

	int count = ged_draw_shape_count(gedp);
	ASSERT(count > 0);

	/* shape_at(0) must be non-NULL and equal to first_shape. */
	test_scene_context *first = ged_draw_first_shape(gedp);
	test_scene_context *at0 = ged_draw_shape_at(gedp, 0);
	ASSERT(at0 != NULL);
	ASSERT(at0 == first);

	ged_draw_shape_ref first_ref = ged_draw_first_shape_ref(gedp);
	ASSERT(!ged_draw_shape_ref_is_null(first_ref));
	ASSERT(ged_draw_shape_ref_index(gedp, first_ref) == 0);
	ged_draw_shape_ref at0_ref = ged_draw_shape_ref_at(gedp, 0);
	ASSERT(!ged_draw_shape_ref_is_null(at0_ref));
	ASSERT(ged_draw_shape_ref_index(gedp, at0_ref) == 0);

	struct ged_draw_shape_record first_rec;
	ASSERT(ged_draw_shape_record_get(gedp, first_ref, &first_rec) == 1);
	ASSERT(first_rec.fullpath != NULL);
	ASSERT(first_rec.leaf_name != NULL);
	ASSERT(first_rec.path_hash != 0);
	ASSERT(first_rec.visible == 1);
	ASSERT(first_rec.stale == 0);
	ASSERT(first_rec.stale_reason != NULL);
	auto db_source_self_or_ancestor = [](test_scene_context *n) -> test_scene_context * {
	    for (test_scene_context *p = n;
		    p != NULL;
		    p = test_scene_parent(p)) {
		struct ged_draw_database_source_summary source_summary;
		if (ged_draw_scene_context_source_summary(p, &source_summary) &&
			source_summary.valid &&
			source_summary.is_database_source)
		    return p;
	    }
	    return NULL;
	};
	test_scene_context *first_source_node = db_source_self_or_ancestor(first);
	ASSERT(first_source_node != NULL);
	struct ged_draw_database_source_summary first_source_summary;
	ASSERT(ged_draw_scene_context_source_summary(first,
		&first_source_summary));
	ASSERT(first_source_summary.valid);
	ASSERT(first_source_summary.is_database_source);
	ASSERT(first_source_summary.has_state);
	struct ged_draw_database_source_summary first_source_container_summary;
	ASSERT(ged_draw_scene_context_source_summary(first_source_node,
		&first_source_container_summary));
	ASSERT(first_source_container_summary.valid);
	ASSERT(first_source_container_summary.is_database_source);
	ASSERT(first_source_container_summary.has_state);
	char *first_path = db_path_to_string(first_rec.fullpath);
	ASSERT(first_source_summary.database_path != NULL);
	ASSERT(BU_STR_EQUAL(first_source_summary.database_path, first_path));
	ASSERT(BU_STR_EQUAL(first_source_container_summary.database_path,
		first_path));
	ASSERT(first_source_summary.source_revision == first_rec.source_revision);
	ASSERT(first_source_container_summary.source_revision ==
		first_rec.source_revision);
	ASSERT(first_source_summary.realized_source_revision ==
		first_rec.realized_source_revision);
	ASSERT(first_source_container_summary.realized_source_revision ==
		first_rec.realized_source_revision);
	ASSERT(!first_source_summary.stale);
	ASSERT(!first_source_container_summary.stale);
	ASSERT(BU_STR_EQUAL(first_source_summary.stale_reason, "current"));
	ASSERT(BU_STR_EQUAL(first_source_container_summary.stale_reason,
		"current"));
	bu_free(first_path, "db_path_to_string");
	ASSERT(ged_draw_shape_ref_set_evaluated_region(gedp, first_ref, 1));
	struct ged_draw_shape_record eval_rec;
	ASSERT(ged_draw_shape_record_get(gedp, first_ref, &eval_rec) == 1);
	ASSERT(eval_rec.evaluated_region == 1);
	ASSERT(ged_draw_shape_ref_set_evaluated_region(gedp, first_ref, 0));

	ged_draw_group_ref first_group_ref = ged_draw_group_ref_of_shape(gedp, first_ref);
	ASSERT(!ged_draw_group_ref_is_null(first_group_ref));
	struct ged_draw_group_record first_group_rec;
	ASSERT(ged_draw_group_record_get(gedp, first_group_ref, &first_group_rec) == 1);
	ASSERT(first_group_rec.path != NULL);
	ASSERT(first_group_rec.is_overlay == 0);

	struct {
	    int seen;
	    int visible;
	    unsigned long long hash;
	} record_iter = {0, 1, 0};
	ged_draw_foreach_shape_record(gedp,
	    [](const struct ged_draw_shape_record *rec, void *ud) -> int {
		auto *ctx = (decltype(record_iter) *)ud;
		if (!rec || ged_draw_shape_ref_is_null(rec->ref) || !rec->fullpath)
		    ctx->visible = 0;
		if (!rec || !rec->visible)
		    ctx->visible = 0;
		ctx->hash ^= rec ? rec->path_hash : 0;
		ctx->seen++;
		return 1;
	    }, &record_iter);
	ASSERT(record_iter.seen == count);
	ASSERT(record_iter.visible == 1);
	ASSERT(record_iter.hash != 0);

	ASSERT(ged_draw_mark_database_change(gedp, first_rec.leaf_name,
		GED_DRAW_STALE_SOURCE_CHANGED) > 0);
	first_ref = ged_draw_first_shape_ref(gedp);
	ASSERT(!ged_draw_shape_ref_is_null(first_ref));
	struct ged_draw_shape_record stale_rec;
	ASSERT(ged_draw_shape_record_get(gedp, first_ref, &stale_rec) == 1);
	ASSERT(stale_rec.stale == 1);
	ASSERT(stale_rec.source_revision > stale_rec.realized_source_revision);
	ASSERT(BU_STR_EQUAL(stale_rec.stale_reason, "source-changed"));
	struct ged_draw_database_source_summary stale_source_summary;
	ASSERT(ged_draw_scene_context_source_summary(first,
		&stale_source_summary));
	ASSERT(stale_source_summary.valid);
	ASSERT(stale_source_summary.source_revision == stale_rec.source_revision);
	ASSERT(stale_source_summary.realized_source_revision ==
		stale_rec.realized_source_revision);
	ASSERT(stale_source_summary.stale);
	ASSERT(BU_STR_EQUAL(stale_source_summary.stale_reason,
		"source-changed"));
	struct ged_draw_database_source_summary stale_source_container_summary;
	ASSERT(ged_draw_scene_context_source_summary(first_source_node,
		&stale_source_container_summary));
	ASSERT(stale_source_container_summary.valid);
	ASSERT(stale_source_container_summary.source_revision ==
		stale_rec.source_revision);
	ASSERT(stale_source_container_summary.realized_source_revision ==
		stale_rec.realized_source_revision);
	ASSERT(stale_source_container_summary.stale);
	ASSERT(BU_STR_EQUAL(stale_source_container_summary.stale_reason,
		"source-changed"));
	struct ged_draw_transaction redraw_op =
	    ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, NULL);
	ASSERT(ged_draw_apply_transaction(gedp, &redraw_op, NULL) > 0);
	first_ref = ged_draw_first_shape_ref(gedp);
	ASSERT(!ged_draw_shape_ref_is_null(first_ref));
	first = ged_draw_shape_node_from_cache_ref(gedp, first_ref);
	ASSERT(first != NULL);
	first_source_node = db_source_self_or_ancestor(first);
	ASSERT(first_source_node != NULL);
	struct ged_draw_shape_record current_rec;
	ASSERT(ged_draw_shape_record_get(gedp, first_ref, &current_rec) == 1);
	ASSERT(current_rec.stale == 0);
	ASSERT(current_rec.source_revision == current_rec.realized_source_revision);
	ASSERT(BU_STR_EQUAL(current_rec.stale_reason, "current"));
	struct ged_draw_database_source_summary current_source_summary;
	ASSERT(ged_draw_scene_context_source_summary(first,
		&current_source_summary));
	ASSERT(current_source_summary.valid);
	ASSERT(current_source_summary.source_revision == current_rec.source_revision);
	ASSERT(current_source_summary.realized_source_revision ==
		current_rec.realized_source_revision);
	ASSERT(!current_source_summary.stale);
	ASSERT(BU_STR_EQUAL(current_source_summary.stale_reason, "current"));
	struct ged_draw_database_source_summary current_source_container_summary;
	ASSERT(ged_draw_scene_context_source_summary(first_source_node,
		&current_source_container_summary));
	ASSERT(current_source_container_summary.valid);
	ASSERT(current_source_container_summary.source_revision ==
		current_rec.source_revision);
	ASSERT(current_source_container_summary.realized_source_revision ==
		current_rec.realized_source_revision);
	ASSERT(!current_source_container_summary.stale);
	ASSERT(BU_STR_EQUAL(current_source_container_summary.stale_reason,
		"current"));

	/* shape_index must round-trip with shape_at. */
	int idx_first = ged_draw_shape_index(gedp, first);
	ASSERT(idx_first == 0);
	ASSERT(ged_draw_shape_ref_index(gedp, first_ref) == 0);

	/* last shape: shape_at(-1) should wrap to count-1. */
	test_scene_context *last = ged_draw_shape_at(gedp, -1);
	ASSERT(last != NULL);
	int idx_last = ged_draw_shape_index(gedp, last);
	ASSERT(idx_last == count - 1);
	ged_draw_shape_ref last_ref = ged_draw_shape_ref_at(gedp, -1);
	ASSERT(ged_draw_shape_ref_index(gedp, last_ref) == count - 1);

	/* advance_shape wraps correctly: last+1 == first. */
	test_scene_context *wrap_fwd = ged_draw_advance_shape(gedp, last, 1);
	ASSERT(wrap_fwd == first);
	ged_draw_shape_ref wrap_fwd_ref = ged_draw_advance_shape_ref(gedp, last_ref, 1);
	ASSERT(ged_draw_shape_ref_index(gedp, wrap_fwd_ref) == 0);

	/* advance_shape backward: first-1 == last. */
	test_scene_context *wrap_bwd = ged_draw_advance_shape(gedp, first, -1);
	ASSERT(wrap_bwd == last);
	ged_draw_shape_ref wrap_bwd_ref = ged_draw_advance_shape_ref(gedp, first_ref, -1);
	ASSERT(ged_draw_shape_ref_index(gedp, wrap_bwd_ref) == count - 1);

	/* Non-drawn pointer returns -1 from shape_index. */
	ASSERT(ged_draw_shape_index(gedp, NULL) == -1);
	ASSERT(ged_draw_shape_ref_index(gedp, GED_DRAW_SHAPE_REF_NULL) == -1);

	/* Overlay shapes should NOT appear in the snapshot. */
	{
	    point_t p1 = {0, 0, 0};
	    point_t points[1];
	    int commands[1] = {BG_GEOMETRY_LINE_MOVE};
	    VMOVE(points[0], p1);
	    struct ged_draw_overlay_geometry geometry = {};
	    geometry.kind = GED_DRAW_OVERLAY_GEOMETRY_LINE_SET;
	    geometry.points = (const point_t *)points;
	    geometry.point_count = 1;
	    geometry.commands = commands;
	    geometry.command_count = 1;
	    ged_draw_shape_ref overlay_ref = GED_DRAW_SHAPE_REF_NULL;
	    ged_draw_overlay_geometry_insert(gedp, "_snap_test_overlay",
			       &geometry, 0xFF0000, 1.0, 0, 0, &overlay_ref);
	    ASSERT(!ged_draw_shape_ref_is_null(overlay_ref));

	    /* count must not have changed */
	    ASSERT(ged_draw_shape_count(gedp) == count);
	    /* Clean up */
	    ged_draw_erase_name(gedp, "_snap_test_overlay");
	}
    }

    /* ---------------------------------------------------------------- *
     * 8. draw_rev / scene_hash revision counter (Step 4 / B7).          *
     * ---------------------------------------------------------------- */
    bu_log("[8] draw_rev revision counter...\n");
    {
	/* Zap → rev must be 0, hash must be 0. */
	{
	    const char *s_av[2] = {"zap", NULL};
	    ged_exec(gedp, 1, s_av);
	}
	ASSERT(ged_draw_scene_revision(gedp) == 0);
	ASSERT(ged_draw_scene_hash(gedp) == 0);

	/* Draw something — rev must be non-zero. */
	{
	    const char *s_av[3] = {"draw", "all.g", NULL};
	    ged_exec(gedp, 2, s_av);
	}
	uint64_t rev_after_draw = ged_draw_scene_revision(gedp);
	ASSERT(rev_after_draw != 0);
	ASSERT(ged_draw_scene_hash(gedp) == (unsigned long long)rev_after_draw);

	/* Erase something — rev must have increased again. */
	{
	    struct db_full_path dfp;
	    db_full_path_init(&dfp);
	    if (db_string_to_path(&dfp, gedp->dbip, "all.g") == 0)
		ged_draw_erase_path(gedp, &dfp);
	    db_free_full_path(&dfp);
	}
	uint64_t rev_after_erase = ged_draw_scene_revision(gedp);
	ASSERT(rev_after_erase > rev_after_draw);

	/* Zap → rev reset to 0. */
	{
	    const char *s_av[2] = {"zap", NULL};
	    ged_exec(gedp, 1, s_av);
	}
	ASSERT(ged_draw_scene_revision(gedp) == 0);

	/* Invent an overlay → rev bumped. */
	{
	    const char *s_av[3] = {"draw", "all.g", NULL};
	    ged_exec(gedp, 2, s_av);
	    uint64_t rev_pre = ged_draw_scene_revision(gedp);
	    test_scene_context *first_before_overlay = ged_draw_first_shape(gedp);
	    ASSERT(first_before_overlay != NULL);
	    ged_draw_shape_ref cache_ref = ged_draw_shape_ref_from_node(gedp, first_before_overlay);

	    point_t p = {1, 1, 1};
	    point_t points[1];
	    int commands[1] = {BG_GEOMETRY_LINE_MOVE};
	    VMOVE(points[0], p);
	    struct ged_draw_overlay_geometry geometry = {};
	    geometry.kind = GED_DRAW_OVERLAY_GEOMETRY_LINE_SET;
	    geometry.points = (const point_t *)points;
	    geometry.point_count = 1;
	    geometry.commands = commands;
	    geometry.command_count = 1;
	    ged_draw_shape_ref overlay_ref = GED_DRAW_SHAPE_REF_NULL;
	    ged_draw_overlay_geometry_insert(gedp, "_rev_test_ov",
			       &geometry, 0x00FF00, 1.0, 0, 0, &overlay_ref);
	    ASSERT(!ged_draw_shape_ref_is_null(overlay_ref));
	    ASSERT(ged_draw_scene_revision(gedp) > rev_pre);
	    struct ged_draw_shape_record cached_rec;
	    ASSERT(ged_draw_shape_record_get(gedp, cache_ref, &cached_rec) == 1);
	    ASSERT(cached_rec.fullpath != NULL);

	    ged_draw_erase_name(gedp, "_rev_test_ov");
	}
    }
    }

    if (draw_scene_shard_is_tree_family(shard)) {

    if (draw_scene_shard_runs(shard, DRAW_SCENE_SHARD_TREE)) {
    /* ---------------------------------------------------------------- *
     * 9. Nested group tree structure (Step 5 — A2+B1+B2).             *
     * ---------------------------------------------------------------- */
    bu_log("[9] nested group tree structure...\n");
    {
	/* Start with a clean slate */
	{
	    const char *s_av[2] = {"zap", NULL};
	    ged_exec(gedp, 1, s_av);
	}

	/* Draw "all.g" */
	{
	    const char *s_av[3] = {"draw", "all.g", NULL};
	    ged_exec(gedp, 2, s_av);
	}

	test_scene_context *root = test_scene_root(gedp);
	ASSERT(root != NULL);

	/* Root should have exactly one non-_overlays child after drawing "all.g" */
	int real_groups = 0;
	test_scene_context *all_g_group = NULL;
	for (size_t i = 0; i < test_scene_child_count(root); i++) {
	    test_scene_context *g = test_scene_child_at(root, i);
	    if (!BU_STR_EQUAL("_overlays", test_scene_name(g))) {
		real_groups++;
		if (!all_g_group)
		    all_g_group = g;
	    }
	}
	ASSERT(real_groups == 1);
	ASSERT(all_g_group != NULL);

	/* Root child must be named "all.g" (single component, not "all.g/hull.r") */
	ASSERT(BU_STR_EQUAL("all.g", test_scene_name(all_g_group)));

	/* Root child must contain sub-groups (not just flat shapes) for any
	 * multi-level hierarchy in moss.g */
	int has_subgroup = 0;
	for (size_t i = 0; i < test_scene_child_count(all_g_group); i++) {
	    test_scene_context *c = test_scene_child_at(all_g_group, i);
	    if (test_scene_is_group(c)) {
		has_subgroup = 1;
		break;
	    }
	}
	ASSERT(has_subgroup);

	/* group_first_shape and group_last_shape must return SHAPE nodes
	 * (not GROUP nodes) even when children include sub-groups */
	test_scene_context *fs = ged_draw_group_first_shape(all_g_group);
	ASSERT(fs != NULL);
	ASSERT(test_scene_is_shape(fs));

	test_scene_context *ls = ged_draw_group_last_shape(all_g_group);
	ASSERT(ls != NULL);
	ASSERT(test_scene_is_shape(ls));

	/* group_is_nonempty must return 1 when shapes exist in sub-tree */
	ASSERT(ged_draw_group_has_shapes(all_g_group) == 1);

	/* group_of_shape must return the root child, not the immediate parent.
	 * Obol-backed context adapters may be re-created, so compare semantic
	 * group identity rather than raw wrapper pointer identity. */
	test_scene_context *fs_group = ged_draw_group_of_shape(gedp, fs);
	ASSERT(fs_group != NULL);
	ASSERT(fs_group && BU_STR_EQUAL(test_scene_name(fs_group),
		test_scene_name(all_g_group)));

	/* group_of_shape on the last shape also returns the root child */
	test_scene_context *ls_group = ged_draw_group_of_shape(gedp, ls);
	ASSERT(ls_group != NULL);
	ASSERT(ls_group && BU_STR_EQUAL(test_scene_name(ls_group),
		test_scene_name(all_g_group)));

	/* Erase "all.g" should clean up cleanly */
	{
	    struct db_full_path dfp;
	    db_full_path_init(&dfp);
	    if (db_string_to_path(&dfp, gedp->dbip, "all.g") == 0)
		ged_draw_erase_path(gedp, &dfp);
	    db_free_full_path(&dfp);
	}
	ASSERT(ged_draw_shape_count(gedp) == 0);
	ASSERT(scene_group_count(gedp) == 0);
    }
    }

    if (draw_scene_shard_runs(shard, DRAW_SCENE_SHARD_HIGHLIGHT_MATERIAL)) {
    /* ---------------------------------------------------------------- *
     * [10] B5: highlight tracking / highlight flag O(1) / material revision (B4 counter).   *
     * ---------------------------------------------------------------- */
    {
	bu_log("[10] highlight/material revision...\n");

	/* Draw all.g again to have some shapes in the tree. */
	{
	    const char *s_av[3] = {"draw", "all.g", NULL};
	    ged_exec(gedp, 2, s_av);
	}
	int ns = ged_draw_shape_count(gedp);
	ASSERT(ns > 0);

	/* highlighted-shape ref lookup returns NULL initially. */
	ASSERT(ged_draw_shape_ref_is_null(ged_draw_highlighted_shape_ref(gedp)));

	/* Highlight the first shape. */
	test_scene_context *s0 = ged_draw_shape_at(gedp, 0);
	ASSERT(s0 != NULL);
	ged_draw_shape_ref s0_ref = ged_draw_shape_ref_at(gedp, 0);
	ASSERT(ged_draw_shape_ref_index(gedp, s0_ref) == 0);
	struct ged_draw_shape_record s0_rec;
	ASSERT(ged_draw_shape_record_get(gedp, s0_ref, &s0_rec) == 1);
	ASSERT(s0_rec.fullpath != NULL);
	ged_draw_set_highlighted_shape_ref(gedp, s0_ref);

	/* highlighted-shape lookup returns s0's ref and s0 is highlighted. */
	ASSERT(ged_draw_highlighted_shape_ref(gedp).token == s0_ref.token);
	ASSERT(test_shape_record_at(gedp, 0, &s0_rec));
	ASSERT(s0_rec.highlighted);

	/* highlight flag(DOWN) should run in O(1) — s0 is the tracked shape. */
	ged_draw_set_highlight_state(gedp, 0);
	ASSERT(test_shape_record_at(gedp, 0, &s0_rec));
	ASSERT(!s0_rec.highlighted);
	ASSERT(ged_draw_shape_ref_is_null(ged_draw_highlighted_shape_ref(gedp)));

	/* Path-prefix highlighting can mark multiple records and deliberately
	 * leaves the single highlighted-ref cache invalid. */
	int prefix_matches = ged_draw_set_highlighted_path_prefix(gedp,
		s0_rec.fullpath, 0, 1);
	ASSERT(prefix_matches > 0);
	ASSERT(ged_draw_shape_ref_is_null(ged_draw_highlighted_shape_ref(gedp)));
	struct ged_draw_shape_record s0_prefix_rec;
	ASSERT(test_shape_record_at(gedp, 0, &s0_prefix_rec));
	ASSERT(s0_prefix_rec.highlighted == 1);
	ged_draw_set_highlight_state(gedp, 0);
	ASSERT(test_shape_record_at(gedp, 0, &s0_rec));
	ASSERT(!s0_rec.highlighted);

	/* highlight(s0) then highlight(s1) clears s0 and highlights s1. */
	if (ns >= 2) {
	    struct ged_draw_shape_record s1_rec;
	    ged_draw_shape_ref s1_ref = ged_draw_shape_ref_at(gedp, 1);
	    ASSERT(!ged_draw_shape_ref_is_null(s1_ref));
	    s0_ref = ged_draw_shape_ref_at(gedp, 0);
	    ASSERT(!ged_draw_shape_ref_is_null(s0_ref));
	    ged_draw_set_highlighted_shape_ref(gedp, s0_ref);
	    ASSERT(test_shape_record_at(gedp, 0, &s0_rec));
	    ASSERT(s0_rec.highlighted);
	    ged_draw_set_highlighted_shape_ref(gedp, s1_ref);
	    ASSERT(test_shape_record_at(gedp, 0, &s0_rec));
	    ASSERT(!s0_rec.highlighted);
	    ASSERT(test_shape_record_at(gedp, 1, &s1_rec));
	    ASSERT(s1_rec.highlighted);
	    ASSERT(ged_draw_highlighted_shape_ref(gedp).token == s1_ref.token);
	    /* Clean up */
	    ged_draw_set_highlight_state(gedp, 0);
	    ASSERT(test_shape_record_at(gedp, 1, &s1_rec));
	    ASSERT(!s1_rec.highlighted);
	}

	/* highlight(NULL) invalidates tracking — subsequent highlight flag(DOWN)
	 * falls back to O(N) sweep (both paths yield correct result). */
	ASSERT(ged_draw_shape_set_highlighted(gedp, s0_ref, 1) == 1);
	ged_draw_highlighted_shape_ref_invalidate(gedp);
	ASSERT(ged_draw_shape_ref_is_null(ged_draw_highlighted_shape_ref(gedp)));
	ged_draw_set_highlight_state(gedp, 0);  /* O(N) fallback */
	/* After O(N) sweep, s0 must not be highlighted. */
	ASSERT(test_shape_record_at(gedp, 0, &s0_rec));
	ASSERT(!s0_rec.highlighted);

	/* B4 activated: material color refresh does NOT bump material revision by itself.
	 * The counter is event-driven: only ged_draw_bump_material_revision() moves it.
	 *
	 * Freshly drawn shapes have material revision 0.
	 * The GED material revision is initialized to 1 so new shapes are always stale.
	 * The first material color refresh call colors them and stamps revision 1.
	 * The counter itself stays at 1 (no self-bump). */
	uint64_t rev0 = ged_draw_material_revision(gedp);
	ged_draw_refresh_material_colors(gedp);
	/* Counter must be unchanged — no material-change event occurred. */
	ASSERT(ged_draw_material_revision(gedp) == rev0);

	/* Verify the first shape was stamped with rev0. */
	ASSERT(s0 != NULL);
	struct ged_draw_shape_material_summary s0_material;
	ASSERT(ged_draw_shape_ref_material_summary(gedp, s0_ref, &s0_material));
	ASSERT(s0_material.valid);
	ASSERT(s0_material.material_revision == rev0);

	/* Simulate a material-change event: bump the counter. */
	ged_draw_bump_material_revision(gedp);
	uint64_t rev1 = ged_draw_material_revision(gedp);
	ASSERT(rev1 == rev0 + 1);

	/* After a bump, material color refresh recolors stale shapes and stamps
	 * them with the new rev, but the counter itself stays put. */
	ged_draw_refresh_material_colors(gedp);
	ASSERT(ged_draw_material_revision(gedp) == rev1);  /* unchanged */
	ASSERT(ged_draw_shape_ref_material_summary(gedp, s0_ref, &s0_material));
	ASSERT(s0_material.material_revision == rev1);      /* stamped at rev1 */

	/* A second call without a bump must leave the material revision current.
	 * Obol owner-state synchronization may re-apply source material color. */
	unsigned char forced_material_rgb[3] = {123, 45, 67};
	ASSERT(ged_draw_shape_ref_set_material_color(gedp, s0_ref,
		forced_material_rgb));
	s0_ref = ged_draw_shape_ref_at(gedp, 0);
	ASSERT(!ged_draw_shape_ref_is_null(s0_ref));
	ged_draw_refresh_material_colors(gedp);  /* skip: material revision is current */
	ASSERT(ged_draw_shape_ref_material_summary(gedp, s0_ref, &s0_material));
	ASSERT(s0_material.material_revision == rev1);  /* stamp unchanged */

	/* highlight pointer is cleared by zap. */
	ged_draw_shape_ref zap_highlight_ref = ged_draw_shape_ref_at(gedp, 0);
	ASSERT(!ged_draw_shape_ref_is_null(zap_highlight_ref));
	ged_draw_set_highlighted_shape_ref(gedp, zap_highlight_ref);
	ASSERT(ged_draw_highlighted_shape_ref(gedp).token == zap_highlight_ref.token);
	ged_draw_clear(gedp);
	ASSERT(ged_draw_shape_ref_is_null(ged_draw_highlighted_shape_ref(gedp)));
	ASSERT(ged_draw_shape_record_get(gedp, s0_ref, &s0_rec) == 0);
    }
    }

    if (draw_scene_shard_runs(shard, DRAW_SCENE_SHARD_VIEW_SCENE)) {
    /* ---------------------------------------------------------------- *
     * [11] A3: view scene registration.                                *
     *      The public scene ref carries the retained scene handle; sync *
     *      remains a no-op because draw/erase mutations maintain it.    *
     * ---------------------------------------------------------------- */
    {
	bu_log("[11] A3: view scene registration...\n");

	/* Draw one object to populate the tree. */
	{
	    const char *dav[3] = {"draw", "all.g", NULL};
	    ged_exec(gedp, 2, dav);
	}

	void *view_ctx = test_active_view_ctx(gedp);
	ASSERT(view_ctx != NULL);

	/* The view scene ref must be set now (registered by draw-scene setup). */
	ASSERT(ged_draw_view_context_scene_attached(view_ctx));

	/* The public view carries only an opaque handle. */
	ASSERT(ged_draw_view_context_scene_attached(view_ctx));

	/* Internal BSG traversal smoke test. */
	test_scene_context *draw_root = test_scene_root(gedp);
	ASSERT(draw_root != NULL);
	ASSERT(ged_draw_view_context_scene_root(view_ctx) == draw_root);

	/* The draw root must have at least one child group (from the draw) */
	test_scene_context *dr = (test_scene_context *)draw_root;
	struct ged_draw_scene_tree_summary draw_root_tree;
	ASSERT(ged_draw_scene_context_tree_summary(draw_root, &draw_root_tree));
	ASSERT(draw_root_tree.valid);
	ASSERT(draw_root_tree.child_count > 0);
	ASSERT(test_scene_child_count(draw_root) == draw_root_tree.child_count);

	/* The draw root should have depth 0 (no parent). */
	ASSERT(draw_root_tree.draw_tree_depth == 0);
	ASSERT(!draw_root_tree.has_parent);

	/* A child's depth should be 1. */
	test_scene_context *first_child = test_scene_child_at(dr, 0);
	ASSERT(first_child != NULL);
	test_scene_context *first_child_parent = test_scene_parent(first_child);
	ASSERT(first_child_parent == draw_root);
	struct ged_draw_scene_tree_summary first_child_tree;
	ASSERT(ged_draw_scene_context_tree_summary(first_child, &first_child_tree));
	ASSERT(first_child_tree.valid);
	if (first_child_tree.is_group || first_child_tree.is_shape) {
	    ASSERT(first_child_tree.draw_tree_depth == 1);
	    ASSERT(first_child_tree.has_parent);
	}

	/* The retained scene is maintained by draw/erase mutation paths. */
	test_scene_context *retained_root = test_scene_root(gedp);
	ASSERT(retained_root == dr);
	ASSERT(test_scene_child_count(retained_root) ==
	       test_scene_child_count(dr));

	/* Obol-backed context wrappers are allowed to be re-created; compare
	 * semantic child identity rather than adapter pointer identity. */
	for (size_t i = 0; i < test_scene_child_count(dr); i++) {
	    test_scene_context *retained_child =
		test_scene_child_at(retained_root, i);
	    test_scene_context *draw_child = test_scene_child_at(dr, i);
	    ASSERT(retained_child != NULL);
	    ASSERT(draw_child != NULL);
	    ASSERT(test_scene_name(retained_child) != NULL);
	    ASSERT(test_scene_name(draw_child) != NULL);
	    if (retained_child && draw_child &&
		    test_scene_name(retained_child) &&
		    test_scene_name(draw_child))
		ASSERT(BU_STR_EQUAL(test_scene_name(retained_child),
			test_scene_name(draw_child)));
	}

	/* After zap the draw root has no children. */
	ged_draw_clear(gedp);
	ASSERT(ged_draw_scene_context_tree_summary(retained_root, &draw_root_tree));
	ASSERT(draw_root_tree.child_count == 0);

	/* The view scene ref itself remains valid after zap (scene anchor not freed). */
	ASSERT(ged_draw_view_context_scene_attached(view_ctx));
    }
    }

    if (draw_scene_shard_runs(shard, DRAW_SCENE_SHARD_BOUNDS_FRAME)) {
    /* ---------------------------------------------------------------- *
     * [12] Phase 9.1 (B3): cached aggregate bbox.                       *
     *      Verify that the public draw-bounds query and scene-ref subtree *
     *      walk agree, and that erase invalidates the semantic result.    *
     * ---------------------------------------------------------------- */
    {
	bu_log("[12] Phase 9.1: cached aggregate bbox...\n");

	/* Start clean. */
	{
	    const char *s_av[2] = {"zap", NULL};
	    ged_exec(gedp, 1, s_av);
	}

	/* Empty draw tree: bounds report empty. */
	{
	    vect_t emin, emax;
	    int empty = ged_draw_bounds(gedp, &emin, &emax, 0);
	    ASSERT(empty == 1);
	}

	/* Draw two paths so the tree has multiple groups. */
	{
	    const char *dav[3] = {"draw", "all.g", NULL};
	    ged_exec(gedp, 2, dav);
	}

	test_scene_context *root = test_scene_root(gedp);
	ASSERT(root != NULL);

	/* First query: cache is cold, must compute. */
	vect_t min1, max1;
	int empty1 = ged_draw_bounds(gedp, &min1, &max1, 0);
	ASSERT(empty1 == 0);

	/* Second query: cache hit, must return the identical bbox. */
	vect_t min2, max2;
	int empty2 = ged_draw_bounds(gedp, &min2, &max2, 0);
	ASSERT(empty2 == 0);
	ASSERT(VNEAR_EQUAL(min1, min2, SMALL_FASTF));
	ASSERT(VNEAR_EQUAL(max1, max2, SMALL_FASTF));

	/* Cross-check vs an explicit walk when the scene-context adapter can
	 * produce subtree bounds; the public draw-bounds query is authoritative
	 * for Obol-primary draw state. */
	vect_t min3, max3;
	int empty3 = ged_draw_scene_context_subtree_bounds(root, &min3, &max3, 1);
	if (empty3 == 0) {
	    ASSERT(VNEAR_EQUAL(min1, min3, SMALL_FASTF));
	    ASSERT(VNEAR_EQUAL(max1, max3, SMALL_FASTF));
	}

	/* Erase invalidates the semantic bounds result. */
	{
	    const char *eav[3] = {"erase", "all.g", NULL};
	    ged_exec(gedp, 2, eav);
	}

	/* Re-query reports empty again (no shapes left). */
	{
	    vect_t emin, emax;
	    int e = ged_draw_bounds(gedp, &emin, &emax, 0);
	    ASSERT(e == 1);
	}
    }


    /* ---------------------------------------------------------------- *
     * [13] Phase 9.2 (B5): frame revision / draw-record stability.       *
     *      Verify frame-revision bumps and zap/redraw do not accumulate   *
     *      stale semantic draw records.                                  *
     * ---------------------------------------------------------------- */
    {
	bu_log("[13] Phase 9.2: frame revision / draw records...\n");

	/* Start clean. */
	{
	    const char *s_av[2] = {"zap", NULL};
	    ged_exec(gedp, 1, s_av);
	}

	/* Initial view frame revision is 0; nothing has been drawn yet. */
	void *view_ctx = test_active_view_ctx(gedp);
	ASSERT(view_ctx != NULL);
	ASSERT(ged_draw_view_context_frame_revision(view_ctx) == 0);

	/* Draw something so we have shapes to stamp. */
	{
	    const char *dav[3] = {"draw", "all.g", NULL};
	    ged_exec(gedp, 2, dav);
	}

	int counted = ged_draw_shape_count(gedp);
	ASSERT(counted > 0);

	/* Frame revision changes must not mutate the semantic draw record set.
	 * The off-screen swrast test exercises the render-time drawn-revision
	 * stamping path; this command-layer test stays on public GED records. */
	ged_draw_view_context_bump_frame_revision(view_ctx);
	ASSERT(ged_draw_shape_count(gedp) == counted);
	ged_draw_view_context_bump_frame_revision(view_ctx);
	ASSERT(ged_draw_shape_count(gedp) == counted);

	/* Zap/redraw should rebuild the same semantic record count, not
	 * accumulate stale records from the previous tree. */
	{
	    const char *zav[2] = {"zap", NULL};
	    ged_exec(gedp, 1, zav);
	}
	ASSERT(ged_draw_shape_count(gedp) == 0);
	{
	    const char *dav[3] = {"draw", "all.g", NULL};
	    ged_exec(gedp, 2, dav);
	}
	ASSERT(ged_draw_shape_count(gedp) == counted);

	/* Erase to leave a clean state. */
	{
	    const char *eav[3] = {"erase", "all.g", NULL};
	    ged_exec(gedp, 2, eav);
	}
    }
    }


    if (draw_scene_shard_runs(shard, DRAW_SCENE_SHARD_HIGHLIGHT_REVISION)) {
    /* ---------------------------------------------------------------- *
     * [14] Highlight draw-ref tracker + highlight revision.            *
     *      Verify the highlight-state revision counter bumps on        *
     *      draw-ref transitions, and that visibility mutations are not  *
     *      authoritative highlight state.                              *
     * ---------------------------------------------------------------- */
    {
	bu_log("[14] highlight draw-ref tracker + highlight revision...\n");

	{
	    const char *s_av[2] = {"zap", NULL};
	    ged_exec(gedp, 1, s_av);
	}
	{
	    const char *dav[3] = {"draw", "all.g", NULL};
	    ged_exec(gedp, 2, dav);
	}

	/* Locate a shape under the draw root. */
	test_scene_context *root = test_scene_root(gedp);
	ASSERT(root != NULL);
	test_scene_context *target = ged_draw_first_shape(gedp);
	ASSERT(target != NULL);
	ged_draw_shape_ref target_ref = ged_draw_first_shape_ref(gedp);
	ASSERT(!ged_draw_shape_ref_is_null(target_ref));

	/* Snapshot the current highlight rev. */
	uint64_t r0 = ged_draw_highlight_revision(gedp);

	/* Transition NULL -> target: rev bumps. */
	ged_draw_set_highlighted_shape_ref(gedp, target_ref);
	uint64_t r1 = ged_draw_highlight_revision(gedp);
	ASSERT(r1 > r0);
	ASSERT(ged_draw_highlighted_shape_ref(gedp).token == target_ref.token);
	struct ged_draw_shape_record target_rec;
	ASSERT(test_shape_record_at(gedp, 0, &target_rec));
	ASSERT(target_rec.highlighted);

	/* Visibility changes on the highlighted shape do not change semantic
	 * highlight identity. */
	ASSERT(ged_draw_shape_ref_set_visible(gedp, target_ref, 0));
	uint64_t r2 = ged_draw_highlight_revision(gedp);
	ASSERT(r2 == r1);

	/* Visibility changes on a non-highlighted node do NOT bump. */
	test_scene_context *other = ged_draw_next_shape(gedp, target);
	if (other && other != target) {
	    ged_draw_shape_ref other_ref = ged_draw_shape_ref_from_node(gedp, other);
	    ASSERT(!ged_draw_shape_ref_is_null(other_ref));
	    ASSERT(ged_draw_shape_ref_set_visible(gedp, other_ref, 0));
	    uint64_t r3 = ged_draw_highlight_revision(gedp);
	    ASSERT(r3 == r2);
	}

	/* Transition target -> NULL: rev bumps. */
	ged_draw_set_highlighted_shape_ref(gedp, GED_DRAW_SHAPE_REF_NULL);
	uint64_t r4 = ged_draw_highlight_revision(gedp);
	ASSERT(r4 > r2);
	ASSERT(ged_draw_shape_ref_is_null(ged_draw_highlighted_shape_ref(gedp)));

	/* After teardown, changing the previously-highlighted shape visibility still
	 * does not bump highlight revision. */
	target_ref = ged_draw_shape_ref_from_node(gedp, target);
	ASSERT(!ged_draw_shape_ref_is_null(target_ref));
	ASSERT(ged_draw_shape_ref_set_visible(gedp, target_ref, 1));
	uint64_t r5 = ged_draw_highlight_revision(gedp);
	ASSERT(r5 == r4);

	{
	    const char *eav[3] = {"erase", "all.g", NULL};
	    ged_exec(gedp, 2, eav);
	}
    }
    }
    }


    if (draw_scene_shard_is_retained_family(shard)) {
    if (draw_scene_shard_runs(shard, DRAW_SCENE_SHARD_PATH_API)) {
    /* ---------------------------------------------------------------- *
     * [15] Phase 10: db_full_path-keyed entry points.                   *
     *      Verify lookup_or_add_dbpath / erase_by_dbpath /              *
     *      erase_all_dbpaths / group_dbpath / group_set_dbpath behave   *
     *      identically to their path-string counterparts.               *
     * ---------------------------------------------------------------- */
    {
	bu_log("[15] Phase 10: db_full_path-keyed entry points...\n");

	{
	    const char *s_av[2] = {"zap", NULL};
	    ged_exec(gedp, 1, s_av);
	}

	struct db_full_path dfp;
	db_full_path_init(&dfp);
	ASSERT(db_string_to_path(&dfp, gedp->dbip, "all.g") == 0);

	/* lookup_or_add via dbpath. */
	ged_draw_group_ref g = ged_draw_group_ref_lookup_or_create(gedp, &dfp);
	ASSERT(!ged_draw_group_ref_is_null(g));
	ASSERT(scene_group_count(gedp) == 1);
	struct ged_draw_group_record grec;
	memset(&grec, 0, sizeof(grec));
	ASSERT(ged_draw_group_record_get(gedp, g, &grec) == 1);
	ASSERT(BU_STR_EQUAL(grec.path, "all.g"));
	ASSERT(grec.fullpath && grec.fullpath->fp_len == dfp.fp_len &&
		db_full_path_match_top(grec.fullpath, &dfp) &&
		db_full_path_match_top(&dfp, grec.fullpath));
	ASSERT(grec.draw_mode == GED_DRAW_MODE_WIRE);

	/* Draw-intent metadata is the canonical path/mode source. */
	ASSERT(ged_draw_group_ref_set_mode(gedp, g, GED_DRAW_MODE_HIDDEN_LINE) == 1);
	memset(&grec, 0, sizeof(grec));
	ASSERT(ged_draw_group_record_get(gedp, g, &grec) == 1);
	ASSERT(grec.draw_mode == GED_DRAW_MODE_HIDDEN_LINE);
	ASSERT(BU_STR_EQUAL(grec.path, "all.g"));

	/* group_set_dbpath keeps the draw-intent path synchronized. */
	ASSERT(ged_draw_group_ref_set_dbpath(gedp, g, &dfp) == 1);
	memset(&grec, 0, sizeof(grec));
	ASSERT(ged_draw_group_record_get(gedp, g, &grec) == 1);
	ASSERT(BU_STR_EQUAL(grec.path, "all.g"));
	ASSERT(grec.fullpath && grec.fullpath->fp_len == dfp.fp_len &&
		db_full_path_match_top(grec.fullpath, &dfp) &&
		db_full_path_match_top(&dfp, grec.fullpath));

	/* erase_by_dbpath removes it. */
	ged_draw_erase_path(gedp, &dfp);
	ASSERT(scene_group_count(gedp) == 0);

	/* erase_all_dbpaths is a no-op on an empty set (no crash). */
	ged_draw_erase_path_prefix(gedp, &dfp);
	ASSERT(scene_group_count(gedp) == 0);

	db_free_full_path(&dfp);
    }
    }


    if (draw_scene_shard_runs(shard, DRAW_SCENE_SHARD_RENAME)) {
    /* ---------------------------------------------------------------- *
     * [15b] Rename transaction updates retained draw intent paths.      *
     * ---------------------------------------------------------------- */
    {
	bu_log("[15b] Retained draw rename transaction...\n");

	const char *zap_av[2] = {"zap", NULL};
	const char *draw_av[3] = {"draw", "all.g", NULL};
	const char *move_av[4] = {"move", "all.g", "all_renamed.g", NULL};
	const char *move_back_av[4] = {"move", "all_renamed.g", "all.g", NULL};

	ASSERT(ged_exec(gedp, 1, zap_av) == BRLCAD_OK);
	ASSERT(ged_exec(gedp, 2, draw_av) == BRLCAD_OK);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 1);

	ASSERT(ged_exec(gedp, 3, move_av) == BRLCAD_OK);
	ASSERT(db_lookup(gedp->dbip, "all.g", LOOKUP_QUIET) == RT_DIR_NULL);
	ASSERT(db_lookup(gedp->dbip, "all_renamed.g", LOOKUP_QUIET) != RT_DIR_NULL);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 0);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all_renamed.g", -1) == 1);
	{
	    struct bu_vls paths = BU_VLS_INIT_ZERO;
	    ASSERT(ged_draw_list_paths(gedp, test_active_view_ctx(gedp), -1, 0, &paths) > 0);
	    ASSERT(!vls_has_line(&paths, "all.g"));
	    ASSERT(vls_has_line(&paths, "all_renamed.g"));
	    bu_vls_free(&paths);
	}

	ASSERT(ged_exec(gedp, 3, move_back_av) == BRLCAD_OK);
	ASSERT(db_lookup(gedp->dbip, "all.g", LOOKUP_QUIET) != RT_DIR_NULL);
	ASSERT(db_lookup(gedp->dbip, "all_renamed.g", LOOKUP_QUIET) == RT_DIR_NULL);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g", -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all_renamed.g", -1) == 0);
	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);

	const char *draw_child_av[3] = {"draw", "all.g/platform.r", NULL};
	const char *mvall_av[4] = {"move_all", "platform.r", "platform_renamed.r", NULL};
	const char *mvall_back_av[4] = {"move_all", "platform_renamed.r", "platform.r", NULL};
	ASSERT(ged_exec(gedp, 2, draw_child_av) == BRLCAD_OK);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g/platform.r", -1) == 1);

	ASSERT(ged_exec(gedp, 3, mvall_av) == BRLCAD_OK);
	ASSERT(db_lookup(gedp->dbip, "platform.r", LOOKUP_QUIET) == RT_DIR_NULL);
	ASSERT(db_lookup(gedp->dbip, "platform_renamed.r", LOOKUP_QUIET) != RT_DIR_NULL);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g/platform.r", -1) == 0);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g/platform_renamed.r", -1) == 1);
	{
	    struct bu_vls paths = BU_VLS_INIT_ZERO;
	    ASSERT(ged_draw_list_paths(gedp, test_active_view_ctx(gedp), -1, 0, &paths) > 0);
	    ASSERT(!vls_has_line(&paths, "all.g/platform.r"));
	    ASSERT(vls_has_line(&paths, "all.g/platform_renamed.r"));
	    bu_vls_free(&paths);
	}

	ASSERT(ged_exec(gedp, 3, mvall_back_av) == BRLCAD_OK);
	ASSERT(db_lookup(gedp->dbip, "platform.r", LOOKUP_QUIET) != RT_DIR_NULL);
	ASSERT(db_lookup(gedp->dbip, "platform_renamed.r", LOOKUP_QUIET) == RT_DIR_NULL);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g/platform.r", -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g/platform_renamed.r", -1) == 0);
	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);

	const char *prefix_av[4] = {"prefix", "pre_", "platform.r", NULL};
	const char *prefix_restore_av[4] = {"move_all", "pre_platform.r", "platform.r", NULL};
	ASSERT(ged_exec(gedp, 2, draw_child_av) == BRLCAD_OK);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g/platform.r", -1) == 1);

	ASSERT(ged_exec(gedp, 3, prefix_av) == BRLCAD_OK);
	ASSERT(db_lookup(gedp->dbip, "platform.r", LOOKUP_QUIET) == RT_DIR_NULL);
	ASSERT(db_lookup(gedp->dbip, "pre_platform.r", LOOKUP_QUIET) != RT_DIR_NULL);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g/platform.r", -1) == 0);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g/pre_platform.r", -1) == 1);
	{
	    struct bu_vls paths = BU_VLS_INIT_ZERO;
	    ASSERT(ged_draw_list_paths(gedp, test_active_view_ctx(gedp), -1, 0, &paths) > 0);
	    ASSERT(!vls_has_line(&paths, "all.g/platform.r"));
	    ASSERT(vls_has_line(&paths, "all.g/pre_platform.r"));
	    bu_vls_free(&paths);
	}

	ASSERT(ged_exec(gedp, 3, prefix_restore_av) == BRLCAD_OK);
	ASSERT(db_lookup(gedp->dbip, "platform.r", LOOKUP_QUIET) != RT_DIR_NULL);
	ASSERT(db_lookup(gedp->dbip, "pre_platform.r", LOOKUP_QUIET) == RT_DIR_NULL);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g/platform.r", -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), "all.g/pre_platform.r", -1) == 0);
	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
    }
    }


    if (draw_scene_shard_runs(shard, DRAW_SCENE_SHARD_REF_REMOVE)) {
    /* ---------------------------------------------------------------- *
     * [15c] Reference-removal transaction preserves direct draws.       *
     * ---------------------------------------------------------------- */
    {
	bu_log("[15c] Retained draw reference-removal transaction...\n");

	const char *kr_parent = "_ged_killrefs_parent.c";
	const char *kr_child = "_ged_killrefs_child.s";
	const char *kr_sibling = "_ged_killrefs_sibling.s";
	const char *kr_child_path = "_ged_killrefs_parent.c/_ged_killrefs_child.s";
	const char *kr_sibling_path = "_ged_killrefs_parent.c/_ged_killrefs_sibling.s";
	struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
	point_t child_center = {100.0, 0.0, 0.0};
	point_t sibling_center = {104.0, 0.0, 0.0};

	ASSERT(wdbp != RT_WDB_NULL);
	if (wdbp) {
	    struct wmember wm;
	    BU_LIST_INIT(&wm.l);
	    ASSERT(mk_sph(wdbp, kr_child, child_center, 1.0) == 0);
	    ASSERT(mk_sph(wdbp, kr_sibling, sibling_center, 1.0) == 0);
	    ASSERT(mk_addmember(kr_child, &wm.l, NULL, WMOP_UNION) != NULL);
	    ASSERT(mk_addmember(kr_sibling, &wm.l, NULL, WMOP_UNION) != NULL);
	    ASSERT(mk_comb(wdbp, kr_parent, &wm.l, 0, NULL, NULL, NULL,
		    0, 0, 0, 0, 0, 0, 0) == 0);
	}

	const char *zap_av[2] = {"zap", NULL};
	const char *draw_parent_av[3] = {"draw", kr_parent, NULL};
	const char *draw_child_av[3] = {"draw", kr_child, NULL};
	const char *killrefs_av[3] = {"killrefs", kr_child, NULL};

	ASSERT(ged_exec(gedp, 1, zap_av) == BRLCAD_OK);
	ASSERT(ged_exec(gedp, 2, draw_parent_av) == BRLCAD_OK);
	ASSERT(ged_exec(gedp, 2, draw_child_av) == BRLCAD_OK);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), kr_parent, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), kr_child_path, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), kr_sibling_path, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), kr_child, -1) == 1);

	ASSERT(ged_exec(gedp, 2, killrefs_av) == BRLCAD_OK);
	ASSERT(db_lookup(gedp->dbip, kr_child, LOOKUP_QUIET) != RT_DIR_NULL);
	ASSERT(db_lookup(gedp->dbip, kr_parent, LOOKUP_QUIET) != RT_DIR_NULL);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), kr_child_path, -1) == 0);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), kr_child, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), kr_sibling_path, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), kr_parent, -1) == 1);

	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
    }


    /* ---------------------------------------------------------------- *
     * [15c2] Comb-child removal reconciles only the affected branch.    *
     * ---------------------------------------------------------------- */
    {
	bu_log("[15c2] Retained comb-child remove reconciliation...\n");

	const char *cr_parent = "_ged_remove_parent.c";
	const char *cr_child = "_ged_remove_child.s";
	const char *cr_sibling = "_ged_remove_sibling.s";
	const char *cr_child_path = "_ged_remove_parent.c/_ged_remove_child.s";
	const char *cr_sibling_path =
	    "_ged_remove_parent.c/_ged_remove_sibling.s";
	struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
	point_t child_center = {112.0, 0.0, 0.0};
	point_t sibling_center = {116.0, 0.0, 0.0};

	ASSERT(wdbp != RT_WDB_NULL);
	if (wdbp) {
	    struct wmember wm;
	    BU_LIST_INIT(&wm.l);
	    ASSERT(mk_sph(wdbp, cr_child, child_center, 1.0) == 0);
	    ASSERT(mk_sph(wdbp, cr_sibling, sibling_center, 1.0) == 0);
	    ASSERT(mk_addmember(cr_child, &wm.l, NULL, WMOP_UNION) != NULL);
	    ASSERT(mk_addmember(cr_sibling, &wm.l, NULL, WMOP_UNION) != NULL);
	    ASSERT(mk_comb(wdbp, cr_parent, &wm.l, 0, NULL, NULL, NULL,
		    0, 0, 0, 0, 0, 0, 0) == 0);
	}

	const char *zap_av[2] = {"zap", NULL};
	const char *draw_parent_av[3] = {"draw", cr_parent, NULL};
	const char *draw_child_av[3] = {"draw", cr_child, NULL};
	const char *remove_child_av[4] = {"remove", cr_parent, cr_child, NULL};

	ASSERT(ged_exec(gedp, 1, zap_av) == BRLCAD_OK);
	ASSERT(ged_exec(gedp, 2, draw_parent_av) == BRLCAD_OK);
	ASSERT(ged_exec(gedp, 2, draw_child_av) == BRLCAD_OK);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), cr_parent, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), cr_child_path, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), cr_sibling_path, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), cr_child, -1) == 1);

	ASSERT(ged_exec(gedp, 3, remove_child_av) == BRLCAD_OK);
	ASSERT(db_lookup(gedp->dbip, cr_child, LOOKUP_QUIET) != RT_DIR_NULL);
	ASSERT(db_lookup(gedp->dbip, cr_parent, LOOKUP_QUIET) != RT_DIR_NULL);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), cr_child_path, -1) == 0);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), cr_child, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), cr_sibling_path, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), cr_parent, -1) == 1);
	ASSERT(test_shape_record_by_path(gedp, cr_child_path, NULL) == 0);

	struct ged_draw_shape_record sibling_rec;
	memset(&sibling_rec, 0, sizeof(sibling_rec));
	ASSERT(test_shape_record_by_path(gedp, cr_sibling_path, &sibling_rec) == 1);
	ASSERT(sibling_rec.stale == 0);

	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
    }


    /* ---------------------------------------------------------------- *
     * [15d] Removed source reconciliation preserves valid siblings.     *
     * ---------------------------------------------------------------- */
    {
	bu_log("[15d] Retained draw source removal reconciliation...\n");

	const char *rm_parent = "_ged_kill_parent.c";
	const char *rm_child = "_ged_kill_child.s";
	const char *rm_sibling = "_ged_kill_sibling.s";
	const char *rm_child_path = "_ged_kill_parent.c/_ged_kill_child.s";
	const char *rm_sibling_path = "_ged_kill_parent.c/_ged_kill_sibling.s";
	struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
	point_t child_center = {120.0, 0.0, 0.0};
	point_t sibling_center = {124.0, 0.0, 0.0};

	ASSERT(wdbp != RT_WDB_NULL);
	if (wdbp) {
	    struct wmember wm;
	    BU_LIST_INIT(&wm.l);
	    ASSERT(mk_sph(wdbp, rm_child, child_center, 1.0) == 0);
	    ASSERT(mk_sph(wdbp, rm_sibling, sibling_center, 1.0) == 0);
	    ASSERT(mk_addmember(rm_child, &wm.l, NULL, WMOP_UNION) != NULL);
	    ASSERT(mk_addmember(rm_sibling, &wm.l, NULL, WMOP_UNION) != NULL);
	    ASSERT(mk_comb(wdbp, rm_parent, &wm.l, 0, NULL, NULL, NULL,
		    0, 0, 0, 0, 0, 0, 0) == 0);
	}

	const char *zap_av[2] = {"zap", NULL};
	const char *draw_parent_av[3] = {"draw", rm_parent, NULL};
	const char *draw_child_av[3] = {"draw", rm_child, NULL};
	const char *kill_child_av[3] = {"kill", rm_child, NULL};

	ASSERT(ged_exec(gedp, 1, zap_av) == BRLCAD_OK);
	ASSERT(ged_exec(gedp, 2, draw_parent_av) == BRLCAD_OK);
	ASSERT(ged_exec(gedp, 2, draw_child_av) == BRLCAD_OK);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), rm_child_path, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), rm_sibling_path, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), rm_child, -1) == 1);

	ASSERT(ged_exec(gedp, 2, kill_child_av) == BRLCAD_OK);
	ASSERT(db_lookup(gedp->dbip, rm_child, LOOKUP_QUIET) == RT_DIR_NULL);
	ASSERT(db_lookup(gedp->dbip, rm_parent, LOOKUP_QUIET) != RT_DIR_NULL);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), rm_child_path, -1) == 0);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), rm_child, -1) == 0);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), rm_sibling_path, -1) == 1);

	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
    }
    }


    if (draw_scene_shard_runs(shard, DRAW_SCENE_SHARD_SOURCE_UPDATE)) {
    /* ---------------------------------------------------------------- *
     * [15e] Source update redraws retained component instances.         *
     * ---------------------------------------------------------------- */
    {
	bu_log("[15e] Retained source update redraw fan-out...\n");

	const char *up_parent = "_ged_update_parent.c";
	const char *up_child = "_ged_update_child.s";
	const char *up_sibling = "_ged_update_sibling.s";
	const char *up_child_path = "_ged_update_parent.c/_ged_update_child.s";
	const char *up_sibling_path = "_ged_update_parent.c/_ged_update_sibling.s";
	struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
	point_t child_center = {140.0, 0.0, 0.0};
	point_t sibling_center = {144.0, 0.0, 0.0};

	ASSERT(wdbp != RT_WDB_NULL);
	if (wdbp) {
	    struct wmember wm;
	    BU_LIST_INIT(&wm.l);
	    ASSERT(mk_sph(wdbp, up_child, child_center, 1.0) == 0);
	    ASSERT(mk_sph(wdbp, up_sibling, sibling_center, 1.0) == 0);
	    ASSERT(mk_addmember(up_child, &wm.l, NULL, WMOP_UNION) != NULL);
	    ASSERT(mk_addmember(up_sibling, &wm.l, NULL, WMOP_UNION) != NULL);
	    ASSERT(mk_comb(wdbp, up_parent, &wm.l, 0, NULL, NULL, NULL,
		    0, 0, 0, 0, 0, 0, 0) == 0);
	}

	const char *zap_av[2] = {"zap", NULL};
	const char *draw_parent_av[3] = {"draw", up_parent, NULL};
	ASSERT(ged_exec(gedp, 1, zap_av) == BRLCAD_OK);
	ASSERT(ged_exec(gedp, 2, draw_parent_av) == BRLCAD_OK);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), up_child_path, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), up_sibling_path, -1) == 1);

	struct ged_draw_shape_record child_rec;
	memset(&child_rec, 0, sizeof(child_rec));
	ASSERT(test_shape_record_by_path(gedp, up_child_path, &child_rec) == 1);
	ASSERT(child_rec.stale == 0);

	struct ged_draw_transaction stale_txn =
	    ged_draw_transaction_make(GED_DRAW_TXN_SOURCE_UPDATED, up_child);
	stale_txn.redraw = 0;
	struct ged_draw_transaction_result result;
	ged_draw_transaction_result_init(&result);
	ASSERT(ged_draw_apply_transaction(gedp, &stale_txn, &result) > 0);
	ASSERT(result.stale_count > 0);
	ged_draw_transaction_result_free(&result);

	memset(&child_rec, 0, sizeof(child_rec));
	ASSERT(test_shape_record_by_path(gedp, up_child_path, &child_rec) == 1);
	ASSERT(child_rec.stale == 1);
	ASSERT(BU_STR_EQUAL(child_rec.stale_reason, "source-changed"));

	struct ged_draw_shape_record sibling_rec;
	memset(&sibling_rec, 0, sizeof(sibling_rec));
	ASSERT(test_shape_record_by_path(gedp, up_sibling_path, &sibling_rec) == 1);
	ASSERT(sibling_rec.stale == 0);

	struct ged_draw_transaction redraw_txn =
	    ged_draw_transaction_make(GED_DRAW_TXN_SOURCE_UPDATED, up_child);
	redraw_txn.redraw = 1;
	ged_draw_transaction_result_init(&result);
	ASSERT(ged_draw_apply_transaction(gedp, &redraw_txn, &result) > 0);
	ASSERT(result.redrawn_count > 0);
	ged_draw_transaction_result_free(&result);

	memset(&child_rec, 0, sizeof(child_rec));
	ASSERT(test_shape_record_by_path(gedp, up_child_path, &child_rec) == 1);
	ASSERT(child_rec.stale == 0);
	ASSERT(child_rec.source_revision == child_rec.realized_source_revision);
	ASSERT(BU_STR_EQUAL(child_rec.stale_reason, "current"));
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), up_sibling_path, -1) == 1);
	memset(&sibling_rec, 0, sizeof(sibling_rec));
	ASSERT(test_shape_record_by_path(gedp, up_sibling_path, &sibling_rec) == 1);
	ASSERT(sibling_rec.stale == 0);

	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
    }


    /* ---------------------------------------------------------------- *
     * [15f] Primitive edits redraw all retained component instances.    *
     * ---------------------------------------------------------------- */
    {
	bu_log("[15f] Native edit redraw fan-out...\n");

	const char *ei_child = "_ged_edit_instance_child.s";
	const char *ei_parent_a = "_ged_edit_instance_parent_a.c";
	const char *ei_parent_b = "_ged_edit_instance_parent_b.c";
	const char *ei_child_path_a =
	    "_ged_edit_instance_parent_a.c/_ged_edit_instance_child.s";
	const char *ei_child_path_b =
	    "_ged_edit_instance_parent_b.c/_ged_edit_instance_child.s";
	struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
	point_t child_center = {180.0, 0.0, 0.0};

	ASSERT(wdbp != RT_WDB_NULL);
	if (wdbp) {
	    struct wmember wm_a;
	    struct wmember wm_b;
	    BU_LIST_INIT(&wm_a.l);
	    BU_LIST_INIT(&wm_b.l);
	    ASSERT(mk_sph(wdbp, ei_child, child_center, 1.0) == 0);
	    ASSERT(mk_addmember(ei_child, &wm_a.l, NULL, WMOP_UNION) != NULL);
	    ASSERT(mk_addmember(ei_child, &wm_b.l, NULL, WMOP_UNION) != NULL);
	    ASSERT(mk_comb(wdbp, ei_parent_a, &wm_a.l, 0, NULL, NULL, NULL,
		    0, 0, 0, 0, 0, 0, 0) == 0);
	    ASSERT(mk_comb(wdbp, ei_parent_b, &wm_b.l, 0, NULL, NULL, NULL,
		    0, 0, 0, 0, 0, 0, 0) == 0);
	}

	const char *zap_av[2] = {"zap", NULL};
	const char *draw_parent_a_av[3] = {"draw", ei_parent_a, NULL};
	const char *draw_parent_b_av[3] = {"draw", ei_parent_b, NULL};
	const char *edit_child_av[8] = {
	    "edit", ei_child, "translate", "-r", "5", "0", "0", NULL
	};

	ASSERT(ged_exec(gedp, 1, zap_av) == BRLCAD_OK);
	ASSERT(ged_exec(gedp, 2, draw_parent_a_av) == BRLCAD_OK);
	ASSERT(ged_exec(gedp, 2, draw_parent_b_av) == BRLCAD_OK);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), ei_child_path_a, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), ei_child_path_b, -1) == 1);

	struct ged_draw_shape_record before_a;
	struct ged_draw_shape_record before_b;
	memset(&before_a, 0, sizeof(before_a));
	memset(&before_b, 0, sizeof(before_b));
	ASSERT(test_shape_record_by_path(gedp, ei_child_path_a, &before_a) == 1);
	ASSERT(test_shape_record_by_path(gedp, ei_child_path_b, &before_b) == 1);
	ASSERT(before_a.stale == 0);
	ASSERT(before_b.stale == 0);
	ASSERT(before_a.source_revision == before_a.realized_source_revision);
	ASSERT(before_b.source_revision == before_b.realized_source_revision);

	ged_draw_index_stats_reset(gedp);
	struct ged_event_txn_result indexed_update_result;
	ged_event_txn_result_init(&indexed_update_result);
	ASSERT(ged_event_notify_object_modified(gedp, ei_child, 1,
		&indexed_update_result) > 0);
	ASSERT(vls_has_word(&indexed_update_result.affected_names,
		ei_child_path_a));
	ASSERT(vls_has_word(&indexed_update_result.affected_names,
		ei_child_path_b));
	ged_event_txn_result_free(&indexed_update_result);
	struct ged_draw_index_stats ei_index_stats;
	memset(&ei_index_stats, 0, sizeof(ei_index_stats));
	ged_draw_index_stats_get(gedp, &ei_index_stats);
	ASSERT(ei_index_stats.fallback_shape_scans == 0);
	ASSERT(ei_index_stats.fallback_group_scans == 0);

	struct ged_draw_shape_record event_a;
	struct ged_draw_shape_record event_b;
	memset(&event_a, 0, sizeof(event_a));
	memset(&event_b, 0, sizeof(event_b));
	ASSERT(test_shape_record_by_path(gedp, ei_child_path_a, &event_a) == 1);
	ASSERT(test_shape_record_by_path(gedp, ei_child_path_b, &event_b) == 1);
	ASSERT(event_a.source_revision > before_a.source_revision);
	ASSERT(event_b.source_revision > before_b.source_revision);
	ASSERT(event_a.source_revision == event_a.realized_source_revision);
	ASSERT(event_b.source_revision == event_b.realized_source_revision);
	ASSERT(event_a.stale == 0);
	ASSERT(event_b.stale == 0);
	before_a = event_a;
	before_b = event_b;

	bu_vls_trunc(gedp->ged_result_str, 0);
	ASSERT(ged_exec(gedp, 7, edit_child_av) == BRLCAD_OK);

	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), ei_child_path_a, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), ei_child_path_b, -1) == 1);
	struct ged_draw_shape_record after_a;
	struct ged_draw_shape_record after_b;
	memset(&after_a, 0, sizeof(after_a));
	memset(&after_b, 0, sizeof(after_b));
	ASSERT(test_shape_record_by_path(gedp, ei_child_path_a, &after_a) == 1);
	ASSERT(test_shape_record_by_path(gedp, ei_child_path_b, &after_b) == 1);
	ASSERT(after_a.source_revision > before_a.source_revision);
	ASSERT(after_b.source_revision > before_b.source_revision);
	ASSERT(after_a.source_revision == after_a.realized_source_revision);
	ASSERT(after_b.source_revision == after_b.realized_source_revision);
	ASSERT(after_a.stale == 0);
	ASSERT(after_b.stale == 0);
	ASSERT(BU_STR_EQUAL(after_a.stale_reason, "current"));
	ASSERT(BU_STR_EQUAL(after_b.stale_reason, "current"));

	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
    }
    }


    if (draw_scene_shard_runs(shard, DRAW_SCENE_SHARD_SELECTION_FANOUT)) {
    /* ---------------------------------------------------------------- *
     * [15g] Selection fan-out tracks all active source instances.       *
     * ---------------------------------------------------------------- */
    {
	bu_log("[15g] Native selection fan-out...\n");

	const char *sf_child = "_ged_select_instance_child.s";
	const char *sf_parent_a = "_ged_select_instance_parent_a.c";
	const char *sf_parent_b = "_ged_select_instance_parent_b.c";
	const char *sf_child_path_a =
	    "_ged_select_instance_parent_a.c/_ged_select_instance_child.s";
	const char *sf_child_path_b =
	    "_ged_select_instance_parent_b.c/_ged_select_instance_child.s";
	struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
	point_t child_center = {188.0, 0.0, 0.0};

	ASSERT(wdbp != RT_WDB_NULL);
	if (wdbp) {
	    struct wmember wm_a;
	    struct wmember wm_b;
	    BU_LIST_INIT(&wm_a.l);
	    BU_LIST_INIT(&wm_b.l);
	    ASSERT(mk_sph(wdbp, sf_child, child_center, 1.0) == 0);
	    ASSERT(mk_addmember(sf_child, &wm_a.l, NULL, WMOP_UNION) != NULL);
	    ASSERT(mk_addmember(sf_child, &wm_b.l, NULL, WMOP_UNION) != NULL);
	    ASSERT(mk_comb(wdbp, sf_parent_a, &wm_a.l, 0, NULL, NULL, NULL,
		    0, 0, 0, 0, 0, 0, 0) == 0);
	    ASSERT(mk_comb(wdbp, sf_parent_b, &wm_b.l, 0, NULL, NULL, NULL,
		    0, 0, 0, 0, 0, 0, 0) == 0);
	}

	const char *zap_av[2] = {"zap", NULL};
	const char *draw_parent_a_av[3] = {"draw", sf_parent_a, NULL};
	const char *draw_parent_b_av[3] = {"draw", sf_parent_b, NULL};
	const char *select_child_av[4] = {"select", "add", sf_child, NULL};
	const char *clear_select_av[3] = {"select", "clear", NULL};

	ASSERT(ged_exec(gedp, 1, zap_av) == BRLCAD_OK);
	ASSERT(ged_exec(gedp, 2, draw_parent_a_av) == BRLCAD_OK);
	ASSERT(ged_exec(gedp, 2, draw_parent_b_av) == BRLCAD_OK);
	ASSERT(ged_exec_select(gedp, 2, clear_select_av) == BRLCAD_OK);

	ged_db_index_id sf_path_ids[2] = {0, 0};
	ASSERT(ged_db_index_path_resolve(gedp, sf_child_path_a,
		sf_path_ids, 2) == 2);
	ASSERT(ged_db_index_valid_id(gedp, sf_path_ids[0]) == 1);
	ASSERT(ged_db_index_valid_id(gedp, sf_path_ids[1]) == 1);
	ged_db_index_id sf_path_hash =
	    ged_db_index_path_hash(gedp, sf_path_ids, 2, 0);
	ged_db_index_id sf_parent_path_hash =
	    ged_db_index_path_hash(gedp, sf_path_ids, 2, 1);
	ASSERT(sf_path_hash != 0);
	ASSERT(sf_parent_path_hash != 0);

	struct bu_vls sf_printed_path = BU_VLS_INIT_ZERO;
	ASSERT(ged_db_index_path_print(gedp, &sf_printed_path, sf_path_ids,
		2, 0, 0) == 1);
	ASSERT(BU_STR_EQUAL(bu_vls_cstr(&sf_printed_path), sf_child_path_a));
	bu_vls_free(&sf_printed_path);

	ASSERT(ged_selection_select_path_ids(gedp, NULL, sf_path_ids, 2, 1) == 1);
	ASSERT(ged_selection_count(gedp, NULL) == 1);
	ASSERT(ged_selection_is_path_selected(gedp, NULL, sf_path_hash) == 1);
	ASSERT(ged_selection_is_path_active(gedp, NULL, sf_path_hash) == 1);
	ASSERT(ged_selection_is_path_active_parent(gedp, NULL,
		sf_parent_path_hash) == 1);
	ASSERT(ged_selection_is_object_immediate_parent(gedp, NULL,
		sf_path_ids[0]) == 1);
	ASSERT(ged_selection_is_object_parent(gedp, NULL, sf_path_ids[0]) == 1);
	ASSERT(ged_selection_deselect_path_ids(gedp, NULL, sf_path_ids, 2, 1) == 1);
	ASSERT(ged_selection_count(gedp, NULL) == 0);
	ASSERT(ged_selection_is_path_selected(gedp, NULL, sf_path_hash) == 0);

	ASSERT(ged_exec_select(gedp, 3, select_child_av) == BRLCAD_OK);
	struct ged_draw_shape_record sel_a;
	struct ged_draw_shape_record sel_b;
	memset(&sel_a, 0, sizeof(sel_a));
	memset(&sel_b, 0, sizeof(sel_b));
	ASSERT(test_shape_record_by_path(gedp, sf_child_path_a, &sel_a) == 1);
	ASSERT(test_shape_record_by_path(gedp, sf_child_path_b, &sel_b) == 1);
	ASSERT(sel_a.selected == 1);
	ASSERT(sel_b.selected == 1);

	ASSERT(ged_exec_select(gedp, 2, clear_select_av) == BRLCAD_OK);
	memset(&sel_a, 0, sizeof(sel_a));
	memset(&sel_b, 0, sizeof(sel_b));
	ASSERT(test_shape_record_by_path(gedp, sf_child_path_a, &sel_a) == 1);
	ASSERT(test_shape_record_by_path(gedp, sf_child_path_b, &sel_b) == 1);
	ASSERT(sel_a.selected == 0);
	ASSERT(sel_b.selected == 0);

	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
    }
    }


    if (draw_scene_shard_runs(shard, DRAW_SCENE_SHARD_COMB_UPDATE)) {
    /* ---------------------------------------------------------------- *
     * [15h] Comb source updates re-expand retained draw intents.        *
     * ---------------------------------------------------------------- */
    {
	bu_log("[15h] Retained comb update re-expansion...\n");

	const char *cu_parent = "_ged_comb_update_parent.c";
	const char *cu_child_a = "_ged_comb_update_a.s";
	const char *cu_child_b = "_ged_comb_update_b.s";
	const char *cu_child_a_path =
	    "_ged_comb_update_parent.c/_ged_comb_update_a.s";
	const char *cu_child_b_path =
	    "_ged_comb_update_parent.c/_ged_comb_update_b.s";
	struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
	point_t child_a_center = {160.0, 0.0, 0.0};
	point_t child_b_center = {164.0, 0.0, 0.0};

	ASSERT(wdbp != RT_WDB_NULL);
	if (wdbp) {
	    struct wmember wm;
	    BU_LIST_INIT(&wm.l);
	    ASSERT(mk_sph(wdbp, cu_child_a, child_a_center, 1.0) == 0);
	    ASSERT(mk_sph(wdbp, cu_child_b, child_b_center, 1.0) == 0);
	    ASSERT(mk_addmember(cu_child_a, &wm.l, NULL, WMOP_UNION) != NULL);
	    ASSERT(mk_comb(wdbp, cu_parent, &wm.l, 0, NULL, NULL, NULL,
		    0, 0, 0, 0, 0, 0, 0) == 0);
	}

	const char *zap_av[2] = {"zap", NULL};
	const char *group_child_av[4] = {"g", cu_parent, cu_child_b, NULL};

	ASSERT(ged_exec(gedp, 1, zap_av) == BRLCAD_OK);
	struct ged_draw_appearance_settings cu_settings =
	    GED_DRAW_APPEARANCE_SETTINGS_INIT;
	cu_settings.s_line_width = 7;
	struct ged_draw_transaction draw_txn =
	    ged_draw_transaction_make(GED_DRAW_TXN_DRAW, cu_parent);
	draw_txn.view = test_active_view_ctx(gedp);
	draw_txn.appearance = &cu_settings;
	ASSERT(ged_draw_apply_transaction(gedp, &draw_txn, NULL) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), cu_child_a_path, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), cu_child_b_path, -1) == 0);
	struct ged_draw_shape_record child_a_rec;
	memset(&child_a_rec, 0, sizeof(child_a_rec));
	ASSERT(test_shape_record_by_path(gedp, cu_child_a_path, &child_a_rec) == 1);
	ASSERT(child_a_rec.line_width == 7);

	ged_draw_index_stats_reset(gedp);
	ASSERT(ged_exec(gedp, 3, group_child_av) == BRLCAD_OK);
	struct ged_draw_index_stats cu_index_stats;
	memset(&cu_index_stats, 0, sizeof(cu_index_stats));
	ged_draw_index_stats_get(gedp, &cu_index_stats);
	ASSERT(cu_index_stats.fallback_shape_scans == 0);
	ASSERT(cu_index_stats.fallback_group_scans == 0);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), cu_parent, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), cu_child_a_path, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp), cu_child_b_path, -1) == 1);
	struct ged_draw_shape_record child_b_rec;
	memset(&child_b_rec, 0, sizeof(child_b_rec));
	ASSERT(test_shape_record_by_path(gedp, cu_child_b_path, &child_b_rec) == 1);
	ASSERT(child_b_rec.stale == 0);
	ASSERT(child_b_rec.line_width == 7);

	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
    }
    }


    if (draw_scene_shard_runs(shard, DRAW_SCENE_SHARD_SCALE)) {
    /* ---------------------------------------------------------------- *
     * [15i] Source update fan-out is indexed at retained-scene scale.   *
     * ---------------------------------------------------------------- */
    {
	bu_log("[15i] Source update reverse-index scale gate...\n");

	const int scale_count = 12;
	const char *scale_target = "_ged_index_scale_target.s";
	const char *scale_parent_a = "_ged_index_scale_parent_a.c";
	const char *scale_parent_b = "_ged_index_scale_parent_b.c";
	const char *scale_target_path_a =
	    "_ged_index_scale_parent_a.c/_ged_index_scale_target.s";
	const char *scale_target_path_b =
	    "_ged_index_scale_parent_b.c/_ged_index_scale_target.s";
	char unrelated0_path[128] = {0};
	std::snprintf(unrelated0_path, sizeof(unrelated0_path),
	    "_ged_index_scale_unrelated_000.c/_ged_index_scale_unrelated_000.s");

	struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
	ASSERT(wdbp != RT_WDB_NULL);
	if (wdbp) {
	    point_t target_center = {220.0, 0.0, 0.0};
	    struct wmember wm_a;
	    struct wmember wm_b;
	    BU_LIST_INIT(&wm_a.l);
	    BU_LIST_INIT(&wm_b.l);
	    ASSERT(mk_sph(wdbp, scale_target, target_center, 1.0) == 0);
	    ASSERT(mk_addmember(scale_target, &wm_a.l, NULL, WMOP_UNION) != NULL);
	    ASSERT(mk_addmember(scale_target, &wm_b.l, NULL, WMOP_UNION) != NULL);
	    ASSERT(mk_comb(wdbp, scale_parent_a, &wm_a.l, 0, NULL, NULL, NULL,
		    0, 0, 0, 0, 0, 0, 0) == 0);
	    ASSERT(mk_comb(wdbp, scale_parent_b, &wm_b.l, 0, NULL, NULL, NULL,
		    0, 0, 0, 0, 0, 0, 0) == 0);

	    for (int i = 0; i < scale_count; i++) {
		char child[64] = {0};
		char parent[64] = {0};
		std::snprintf(child, sizeof(child),
		    "_ged_index_scale_unrelated_%03d.s", i);
		std::snprintf(parent, sizeof(parent),
		    "_ged_index_scale_unrelated_%03d.c", i);
		point_t center = {240.0 + (fastf_t)i * 3.0, 0.0, 0.0};
		struct wmember wm;
		BU_LIST_INIT(&wm.l);
		ASSERT(mk_sph(wdbp, child, center, 1.0) == 0);
		ASSERT(mk_addmember(child, &wm.l, NULL, WMOP_UNION) != NULL);
		ASSERT(mk_comb(wdbp, parent, &wm.l, 0, NULL, NULL, NULL,
			0, 0, 0, 0, 0, 0, 0) == 0);
	    }
	}

	const char *zap_av[2] = {"zap", NULL};
	ASSERT(ged_exec(gedp, 1, zap_av) == BRLCAD_OK);

	struct ged_draw_transaction draw_a =
	    ged_draw_transaction_make(GED_DRAW_TXN_DRAW, scale_parent_a);
	struct ged_draw_transaction draw_b =
	    ged_draw_transaction_make(GED_DRAW_TXN_DRAW, scale_parent_b);
	ASSERT(ged_draw_apply_transaction(gedp, &draw_a, NULL) == 1);
	ASSERT(ged_draw_apply_transaction(gedp, &draw_b, NULL) == 1);

	for (int i = 0; i < scale_count; i++) {
	    char parent[64] = {0};
	    std::snprintf(parent, sizeof(parent),
		"_ged_index_scale_unrelated_%03d.c", i);
	    struct ged_draw_transaction draw_unrelated =
		ged_draw_transaction_make(GED_DRAW_TXN_DRAW, parent);
	    ASSERT(ged_draw_apply_transaction(gedp, &draw_unrelated, NULL) == 1);
	}

	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp),
		scale_target_path_a, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp),
		scale_target_path_b, -1) == 1);
	ASSERT(ged_draw_path_state(gedp, test_active_view_ctx(gedp),
		unrelated0_path, -1) == 1);

	ASSERT(ged_selection_clear(gedp, NULL) == 1);
	ASSERT(ged_selection_select_path(gedp, NULL, scale_target_path_a, 1) == 1);
	ged_draw_index_stats_reset(gedp);
	ASSERT(ged_selection_draw_sync(gedp, NULL) == 1);
	struct ged_draw_index_stats selection_scale_stats;
	memset(&selection_scale_stats, 0, sizeof(selection_scale_stats));
	ged_draw_index_stats_get(gedp, &selection_scale_stats);
	ASSERT(selection_scale_stats.fallback_shape_scans == 0);
	ASSERT(selection_scale_stats.fallback_group_scans == 0);

	struct ged_draw_shape_record selected_target;
	struct ged_draw_shape_record unselected_target;
	memset(&selected_target, 0, sizeof(selected_target));
	memset(&unselected_target, 0, sizeof(unselected_target));
	ASSERT(test_shape_record_by_path(gedp, scale_target_path_a,
		&selected_target) == 1);
	ASSERT(test_shape_record_by_path(gedp, scale_target_path_b,
		&unselected_target) == 1);
	ASSERT(selected_target.selected == 1);
	ASSERT(unselected_target.selected == 0);
	ASSERT(ged_selection_clear(gedp, NULL) == 1);
	ASSERT(ged_selection_draw_sync(gedp, NULL) == 1);

	struct ged_draw_shape_record unrelated_before;
	memset(&unrelated_before, 0, sizeof(unrelated_before));
	ASSERT(test_shape_record_by_path(gedp, unrelated0_path,
		&unrelated_before) == 1);

	ged_draw_index_stats_reset(gedp);
	struct ged_event_txn_result scale_result;
	ged_event_txn_result_init(&scale_result);
	ASSERT(ged_event_notify_object_modified(gedp, scale_target, 1,
		&scale_result) > 0);
	ASSERT(vls_has_word(&scale_result.affected_names,
		scale_target_path_a));
	ASSERT(vls_has_word(&scale_result.affected_names,
		scale_target_path_b));
	ged_event_txn_result_free(&scale_result);

	struct ged_draw_index_stats scale_stats;
	memset(&scale_stats, 0, sizeof(scale_stats));
	ged_draw_index_stats_get(gedp, &scale_stats);
	ASSERT(scale_stats.fallback_shape_scans == 0);
	ASSERT(scale_stats.fallback_group_scans == 0);

	struct ged_draw_shape_record unrelated_after;
	memset(&unrelated_after, 0, sizeof(unrelated_after));
	ASSERT(test_shape_record_by_path(gedp, unrelated0_path,
		&unrelated_after) == 1);
	ASSERT(unrelated_after.source_revision ==
		unrelated_before.source_revision);
	ASSERT(unrelated_after.realized_source_revision ==
		unrelated_before.realized_source_revision);

	ged_draw_clear(gedp);
	ASSERT(scene_group_count(gedp) == 0);
    }
    }
    }


    if (draw_scene_shard_runs(shard, DRAW_SCENE_SHARD_BACKEND)) {
    /* ---------------------------------------------------------------- *
     * [16] Backend-owned resource cache.                                *
     *      Resources are keyed by render-item cache/source identity and  *
     *      owned by the display manager, not scene node.                   *
     * ---------------------------------------------------------------- */
    {
	bu_log("[16] Backend-owned resource cache...\n");

	struct phase16_state {
	    int free_calls;
	    uint64_t last_cache_identity;
	    uint64_t last_source_identity;
	} st;
	memset(&st, 0, sizeof(st));

	struct phase16_helpers {
	    static void free_resource(struct dm *UNUSED(dmp), struct dm_backend_resource *r) {
		struct phase16_state *p = (struct phase16_state *)r->handle;
		p->free_calls++;
		p->last_cache_identity = r->cache_identity;
		p->last_source_identity = r->source_identity;
		r->handle = NULL;
	    }
	};

	const char *av0 = "attach";
	struct dm *dmp = dm_open(NULL, NULL, "nu", 1, &av0);
	ASSERT(dmp != NULL);

	const uint64_t backend_identity = 0x9001;
	ASSERT(dm_backend_begin_frame(dmp) == 0);
	struct dm_backend_resource *r = dm_backend_resource_get(dmp,
		0x1001, 0x501, backend_identity, 1,
		phase16_helpers::free_resource);
	ASSERT(r != NULL);
	r->handle = &st;
	dm_backend_resource_touch(dmp, r);

	dm_backend_resource_invalidate(dmp, 0x501);
	ASSERT(r->stale == 1);
	dm_backend_release_source(dmp, 0x501);
	ASSERT(st.free_calls == 1);
	ASSERT(st.last_cache_identity == 0x1001);
	ASSERT(st.last_source_identity == 0x501);

	r = dm_backend_resource_get(dmp, 0x1002, 0x502, backend_identity, 1,
		phase16_helpers::free_resource);
	ASSERT(r != NULL);
	r->handle = &st;
	dm_backend_resource_touch(dmp, r);
	ASSERT(dm_backend_begin_frame(dmp) == 0);
	dm_backend_resource_reclaim_unseen(dmp);
	ASSERT(st.free_calls == 2);
	ASSERT(st.last_cache_identity == 0x1002);
	ASSERT(st.last_source_identity == 0x502);

	dm_backend_resource_invalidate(NULL, 0x501);
	dm_backend_release_source(NULL, 0x501);
	dm_close(dmp);
    }
    }


    /* Final zap to leave clean state. */
    {
	const char *s_av[2] = {"zap", NULL};
	ged_exec(gedp, 1, s_av);
	ASSERT(ged_draw_scene_hash(gedp) == 0);
	ASSERT(scene_group_count(gedp) == 0);
    }

    /* ---------------------------------------------------------------- *
     * Summary.                                                          *
     * ---------------------------------------------------------------- */
    ged_close(gedp);

    bu_log("=== ged_draw_* API regression: %d/%d checks passed ===\n",
	   nchecks - nfails, nchecks);
    return (nfails == 0) ? 0 : 1;
}


/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
