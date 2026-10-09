/*                         A E T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2018-2026 United States Government as represented by
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
/** @file aet.cpp
 *
 * Tests shared-view progressive autoview and the stability of view state and
 * retained rendering across resize events.
 */

#include "common.h"

#include <stdio.h>
#include <fstream>

#include <bu.h>
#include <ged.h>
#include <ged/draw.h>
#include <icv.h>
#include <rt/view.h>

#include "view_test_util.h"

extern "C" int draw_test_obol_screengrab_view_if_enabled(struct ged *gedp,
	struct ged_view_context *view_ctx, int id, const char *filename);

static const size_t AET_VIEW_COUNT = 4;
static const int AET_INITIAL_DIMENSION = 512;
static const int AET_LARGE_DIMENSION = 600;
static const int AET_CYCLE_FIRST_DIMENSION = 513;
static const int AET_CYCLE_LAST_DIMENSION = 599;
static const fastf_t AET_SENTINEL_CENTER_BASE = 1000.0;
static const fastf_t AET_SENTINEL_CENTER_STEP = 100.0;
static const fastf_t AET_SENTINEL_SIZE_BASE = 1000.0;
static const fastf_t AET_SENTINEL_SIZE_STEP = 125.0;
/* Resize updates should reproduce retained camera state, not merely remain
 * visually close.  This tolerance allows only floating-point roundoff from a
 * deterministic matrix refresh. */
static const fastf_t AET_STATE_TOLERANCE = 1.0e-12;

enum aet_capture_stage {
    AET_CAPTURE_INITIAL_512 = 0,
    AET_CAPTURE_INITIAL_600,
    AET_CAPTURE_SHRINK_512,
    AET_CAPTURE_POST_CYCLE_600,
    AET_CAPTURE_FINAL_512,
    AET_CAPTURE_STAGE_COUNT
};

static const char *aet_capture_stage_name[AET_CAPTURE_STAGE_COUNT] = {
    "initial_512",
    "initial_600",
    "shrink_512",
    "post_cycle_600",
    "final_512"
};

struct aet_view_state {
    point_t center;
    point_t eye;
    vect_t aet;
    mat_t model2view;
    fastf_t size;
    fastf_t scale;
    fastf_t perspective;
};

static int
aet_view_state_get(struct aet_view_state *state,
	const struct ged_view_context *view_ctx)
{
    const struct bv *view = DRAW_TEST_BV_CONST(view_ctx);
    if (!state || !view ||
	!bv_center_get(state->center, view) ||
	!bv_eye_pos_get(state->eye, view) ||
	!bv_aet_get(state->aet, view) ||
	!bv_model2view_get(state->model2view, view))
	return 0;

    state->size = bv_size_get(view);
    state->scale = bv_scale_get(view);
    state->perspective = bv_perspective_get(view);
    return 1;
}

static int
aet_view_state_matches(const struct aet_view_state *expected,
	const struct ged_view_context *view_ctx, size_t view_index,
	const char *phase)
{
    struct aet_view_state actual;
    if (!expected || !aet_view_state_get(&actual, view_ctx)) {
	bu_log("FAIL: V%zu could not read view state after %s\n",
	    view_index, phase);
	return 0;
    }

    int matches = 1;
    if (!VNEAR_EQUAL(expected->center, actual.center,
	    AET_STATE_TOLERANCE)) {
	bu_log("FAIL: V%zu center changed after %s: "
	       "(%0.17g %0.17g %0.17g) -> (%0.17g %0.17g %0.17g)\n",
	       view_index, phase, V3ARGS(expected->center),
	       V3ARGS(actual.center));
	matches = 0;
    }
    if (!VNEAR_EQUAL(expected->eye, actual.eye, AET_STATE_TOLERANCE)) {
	bu_log("FAIL: V%zu eye changed after %s: "
	       "(%0.17g %0.17g %0.17g) -> (%0.17g %0.17g %0.17g)\n",
	       view_index, phase, V3ARGS(expected->eye), V3ARGS(actual.eye));
	matches = 0;
    }
    if (!VNEAR_EQUAL(expected->aet, actual.aet, AET_STATE_TOLERANCE)) {
	bu_log("FAIL: V%zu AET changed after %s: "
	       "(%0.17g %0.17g %0.17g) -> (%0.17g %0.17g %0.17g)\n",
	       view_index, phase, V3ARGS(expected->aet), V3ARGS(actual.aet));
	matches = 0;
    }
    if (!NEAR_EQUAL(expected->size, actual.size, AET_STATE_TOLERANCE) ||
	!NEAR_EQUAL(expected->scale, actual.scale, AET_STATE_TOLERANCE) ||
	!NEAR_EQUAL(expected->perspective, actual.perspective,
	    AET_STATE_TOLERANCE)) {
	bu_log("FAIL: V%zu scale state changed after %s: "
	       "size=%0.17g/%0.17g scale=%0.17g/%0.17g "
	       "perspective=%0.17g/%0.17g\n",
	       view_index, phase, expected->size, actual.size,
	       expected->scale, actual.scale, expected->perspective,
	       actual.perspective);
	matches = 0;
    }
    for (size_t i = 0; i < 16; i++) {
	if (NEAR_EQUAL(expected->model2view[i], actual.model2view[i],
		AET_STATE_TOLERANCE))
	    continue;
	bu_log("FAIL: V%zu model2view[%zu] changed after %s: "
	       "%0.17g -> %0.17g\n", view_index, i, phase,
	       expected->model2view[i], actual.model2view[i]);
	matches = 0;
	break;
    }
    return matches;
}

static int
aet_verify_shared_autoview(struct ged_view_context *views[AET_VIEW_COUNT],
	const struct aet_view_state sentinels[AET_VIEW_COUNT],
	struct aet_view_state fitted[AET_VIEW_COUNT])
{
    int failures = 0;
    for (size_t i = 0; i < AET_VIEW_COUNT; i++) {
	if (!aet_view_state_get(&fitted[i], views[i])) {
	    bu_log("FAIL: V%zu could not read the progressive autoview result\n",
		i);
	    fitted[i] = sentinels[i];
	    failures++;
	    continue;
	}
	if (VNEAR_EQUAL(fitted[i].center, sentinels[i].center,
		AET_STATE_TOLERANCE)) {
	    bu_log("FAIL: V%zu retained its distinct pre-draw center\n", i);
	    failures++;
	}
	if (NEAR_EQUAL(fitted[i].size, sentinels[i].size,
		AET_STATE_TOLERANCE)) {
	    bu_log("FAIL: V%zu retained its distinct pre-draw size\n", i);
	    failures++;
	}
	if (!VNEAR_EQUAL(fitted[i].aet, sentinels[i].aet,
		AET_STATE_TOLERANCE)) {
	    bu_log("FAIL: V%zu progressive autoview changed its AET\n", i);
	    failures++;
	}
    }

    for (size_t i = 1; i < AET_VIEW_COUNT; i++) {
	if (!VNEAR_EQUAL(fitted[0].center, fitted[i].center,
		AET_STATE_TOLERANCE) ||
	    !NEAR_EQUAL(fitted[0].size, fitted[i].size,
		AET_STATE_TOLERANCE)) {
	    bu_log("FAIL: shared progressive autoview diverged: "
		   "V0 center=(%0.17g %0.17g %0.17g) size=%0.17g, "
		   "V%zu center=(%0.17g %0.17g %0.17g) size=%0.17g\n",
		   V3ARGS(fitted[0].center), fitted[0].size, i,
		   V3ARGS(fitted[i].center), fitted[i].size);
	    failures++;
	}
    }
    if (!failures)
	bu_log("PASS: progressive automatic autoview reached all shared views\n");
    return failures;
}

static int
aet_resize_views(struct ged *gedp,
	struct ged_view_context *views[AET_VIEW_COUNT], int dimension,
	const struct aet_view_state fitted[AET_VIEW_COUNT], const char *phase)
{
    int failures = 0;
    for (size_t i = 0; i < AET_VIEW_COUNT; i++) {
	if (draw_test_obol_view_init(gedp, views[i], dimension, dimension) !=
		BRLCAD_OK) {
	    bu_log("FAIL: V%zu endpoint resize to %dx%d failed during %s\n",
		i, dimension, dimension, phase);
	    failures++;
	    continue;
	}
	(void)ged_view_context_update(views[i]);
	const struct bv *view = DRAW_TEST_BV_CONST(views[i]);
	if (bv_width_get(view) != dimension ||
	    bv_height_get(view) != dimension) {
	    bu_log("FAIL: V%zu dimensions after %s are %dx%d, expected %dx%d\n",
		i, phase, bv_width_get(view), bv_height_get(view), dimension,
		dimension);
	    failures++;
	}
	if (!aet_view_state_matches(&fitted[i], views[i], i, phase))
	    failures++;
    }
    return failures;
}

static int
aet_validate_capture(const char *filename, int dimension, int require_nonblack)
{
    icv_image_t *image = icv_read(filename, BU_MIME_IMAGE_PNG, 0, 0);
    if (!image) {
	bu_log("FAIL: could not read capture %s\n", filename);
	return 0;
    }

    int valid = 1;
    if (image->width != static_cast<size_t>(dimension) ||
	image->height != static_cast<size_t>(dimension)) {
	bu_log("FAIL: capture %s is %zux%zu, expected %dx%d\n", filename,
	    image->width, image->height, dimension, dimension);
	valid = 0;
    }
    if (require_nonblack) {
	size_t nonblack = 0;
	const size_t pixel_count = image->width * image->height;
	for (size_t i = 0; i < pixel_count; i++) {
	    if (image->data[3 * i] > 0.001 ||
		image->data[3 * i + 1] > 0.001 ||
		image->data[3 * i + 2] > 0.001)
		nonblack++;
	}
	if (!nonblack) {
	    bu_log("FAIL: initial capture %s contains no geometry\n", filename);
	    valid = 0;
	}
    }
    icv_destroy(image);
    return valid;
}

static int
aet_capture_views(struct ged *gedp,
	struct ged_view_context *views[AET_VIEW_COUNT],
	const struct aet_view_state fitted[AET_VIEW_COUNT],
	enum aet_capture_stage stage, int dimension,
	char captures[AET_CAPTURE_STAGE_COUNT][AET_VIEW_COUNT][MAXPATHLEN])
{
    int failures = 0;
    for (size_t i = 0; i < AET_VIEW_COUNT; i++) {
	if (!draw_test_obol_progressive_drain(gedp, views[i], 2000, 1)) {
	    bu_log("FAIL: V%zu did not settle before %s capture\n", i,
		aet_capture_stage_name[stage]);
	    failures++;
	}
	if (!aet_view_state_matches(&fitted[i], views[i], i,
		aet_capture_stage_name[stage]))
	    failures++;
	const int capture_id = static_cast<int>(stage * AET_VIEW_COUNT + i + 1);
	if (draw_test_obol_screengrab_view_if_enabled(gedp, views[i],
		capture_id, captures[stage][i]) != 1) {
	    bu_log("FAIL: could not capture V%zu at %s\n", i,
		aet_capture_stage_name[stage]);
	    failures++;
	    continue;
	}
	const int initial = stage == AET_CAPTURE_INITIAL_512 ||
	    stage == AET_CAPTURE_INITIAL_600;
	if (!aet_validate_capture(captures[stage][i], dimension, initial))
	    failures++;
    }
    return failures;
}

static int
aet_images_match_exactly(const char *first, const char *second)
{
    icv_image_t *first_image = icv_read(first, BU_MIME_IMAGE_PNG, 0, 0);
    icv_image_t *second_image = icv_read(second, BU_MIME_IMAGE_PNG, 0, 0);
    if (!first_image || !second_image) {
	bu_log("FAIL: could not read same-run captures %s and %s\n", first,
	    second);
	if (first_image)
	    icv_destroy(first_image);
	if (second_image)
	    icv_destroy(second_image);
	return 0;
    }

    int exact = first_image->width == second_image->width &&
	first_image->height == second_image->height &&
	first_image->color_space == second_image->color_space;
    int matching = 0;
    int off_by_one = 0;
    int off_by_many = 0;
    if (exact && icv_diff(&matching, &off_by_one, &off_by_many,
	    first_image, second_image))
	exact = 0;
    if (!exact) {
	bu_log("FAIL: same-run captures differ: %s vs %s "
	       "(matching=%d off-by-one=%d off-by-many=%d)\n",
	       first, second, matching, off_by_one, off_by_many);
    }

    icv_destroy(first_image);
    icv_destroy(second_image);
    return exact;
}

static void
aet_capture_names_init(
	char captures[AET_CAPTURE_STAGE_COUNT][AET_VIEW_COUNT][MAXPATHLEN])
{
    for (size_t stage = 0; stage < AET_CAPTURE_STAGE_COUNT; stage++) {
	for (size_t view = 0; view < AET_VIEW_COUNT; view++) {
	    snprintf(captures[stage][view], MAXPATHLEN, "aet_v%zu_%s.png",
		view, aet_capture_stage_name[stage]);
	    /* A failed capture must not accidentally compare a prior run's file. */
	    bu_file_delete(captures[stage][view]);
	}
    }
}

static void
aet_capture_cleanup(
	char captures[AET_CAPTURE_STAGE_COUNT][AET_VIEW_COUNT][MAXPATHLEN])
{
    for (size_t stage = 0; stage < AET_CAPTURE_STAGE_COUNT; stage++) {
	for (size_t view = 0; view < AET_VIEW_COUNT; view++)
	    bu_file_delete(captures[stage][view]);
    }
}

int
main(int ac, char *av[])
{
    int need_help = 0;
    int soft_fail = 0;
    int keep_images = 0;

    bu_setprogname(av[0]);

    struct bu_opt_desc options[4];
    BU_OPT(options[0], "h", "help", "", NULL, &need_help,
	"Print help and exit");
    BU_OPT(options[1], "c", "continue", "", NULL, &soft_fail,
	"Continue testing if a failure is encountered.");
    BU_OPT(options[2], "k", "keep", "", NULL, &keep_images,
	"Keep images generated by the run.");
    BU_OPT_NULL(options[3]);

    const int unused_count = bu_opt_parse(NULL, ac, (const char **)av,
	options);
    if (unused_count != 2 || need_help)
	bu_exit(EXIT_FAILURE, "%s [-h] [-c] [-k] <data directory>", av[0]);
    const char *data_directory = av[1];
    if (!bu_file_directory(data_directory)) {
	bu_log("ERROR: [%s] is not a data directory\n", data_directory);
	return 2;
    }

    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    char cache_root[MAXPATHLEN] = {0};
    char runtime_cache[MAXPATHLEN] = {0};
    bu_dir(cache_root, MAXPATHLEN, BU_DIR_CURR, "ged_aet_test_cache", NULL);
    bu_mkdir(cache_root);
    bu_dir(runtime_cache, MAXPATHLEN, BU_DIR_CURR, "ged_aet_test_cache",
	"cache", NULL);
    bu_mkdir(runtime_cache);
    bu_setenv("BU_DIR_CACHE", runtime_cache, 1);

    struct bu_vls source_name = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&source_name, "%s/moss.g", data_directory);
    std::ifstream source(bu_vls_cstr(&source_name), std::ios::binary);
    std::ofstream destination("moss_aet_tmp.g", std::ios::binary);
    if (!source || !destination) {
	bu_log("ERROR: could not prepare the AET test database from %s\n",
	    bu_vls_cstr(&source_name));
	bu_vls_free(&source_name);
	return 2;
    }
    destination << source.rdbuf();
    source.close();
    destination.close();
    bu_vls_free(&source_name);

    struct ged *gedp = ged_open("db", "moss_aet_tmp.g", 1);
    if (!gedp) {
	bu_log("ERROR: could not open the AET test database\n");
	bu_file_delete("moss_aet_tmp.g");
	return 2;
    }

    char captures[AET_CAPTURE_STAGE_COUNT][AET_VIEW_COUNT][MAXPATHLEN];
    aet_capture_names_init(captures);

    struct ged_view_set *view_set_ctx = ged_view_set_ctx(gedp);
    ged_view_set_context_remove(view_set_ctx, NULL);

    struct ged_view_context *views[AET_VIEW_COUNT] = {NULL};
    int failures = 0;
    for (size_t i = 0; i < AET_VIEW_COUNT; i++) {
	char view_name[16] = {0};
	snprintf(view_name, sizeof(view_name), "V%zu", i);
	views[i] = ged_view_context_create_with_set(view_set_ctx);
	if (!views[i]) {
	    bu_log("FAIL: could not create %s\n", view_name);
	    failures++;
	    break;
	}
	if (!i)
	    ged_view_active_ctx_set(gedp, views[i]);
	bv_name_set(DRAW_TEST_BV(views[i]), view_name);
	ged_view_set_context_add(view_set_ctx, views[i]);
	ged_view_context_owned_add(gedp, views[i]);
	if (draw_test_obol_view_init(gedp, views[i], AET_INITIAL_DIMENSION,
		AET_INITIAL_DIMENSION) != BRLCAD_OK) {
	    bu_log("FAIL: could not initialize the Obol endpoint for %s\n",
		view_name);
	    failures++;
	    break;
	}
    }

    struct aet_view_state sentinels[AET_VIEW_COUNT];
    struct aet_view_state fitted[AET_VIEW_COUNT];
    const vect_t requested_aet[AET_VIEW_COUNT] = {
	{0.0, 0.0, 90.0},
	{90.0, 90.0, 180.0},
	{-90.0, 270.0, -90.0},
	{270.0, -180.0, 90.0}
    };

    do {
	if (failures)
	    break;
	for (size_t i = 0; i < AET_VIEW_COUNT; i++) {
	    point_t sentinel_center;
	    const fastf_t center_offset = AET_SENTINEL_CENTER_BASE +
		static_cast<fastf_t>(i) * AET_SENTINEL_CENTER_STEP;
	    VSET(sentinel_center, center_offset, -2.0 * center_offset,
		3.0 * center_offset);
	    const fastf_t sentinel_size = AET_SENTINEL_SIZE_BASE +
		static_cast<fastf_t>(i) * AET_SENTINEL_SIZE_STEP;
	    if (!bv_aet_set(DRAW_TEST_BV(views[i]), requested_aet[i]) ||
		!bv_center_set(DRAW_TEST_BV(views[i]), sentinel_center) ||
		!bv_size_set(DRAW_TEST_BV(views[i]), sentinel_size) ||
		!ged_view_context_update(views[i]) ||
		!aet_view_state_get(&sentinels[i], views[i])) {
		bu_log("FAIL: could not establish the V%zu pre-draw sentinel\n",
		    i);
		failures++;
	    }
	}
	/* No later assertion is meaningful without complete sentinel snapshots. */
	if (failures)
	    break;

	bu_log("Testing shared progressive automatic autoview...\n");
	const char *draw_argv[] = {
	    "draw", "--defer-leaf-expansion", "-m0", "all.g"
	};
	if (ged_exec_draw(gedp, 4, draw_argv) != BRLCAD_OK) {
	    bu_log("FAIL: deferred draw command failed\n");
	    failures++;
	}
	if (!draw_test_obol_progressive_drain(gedp, views[0], 2000, 1)) {
	    bu_log("FAIL: progressive realization owner V0 did not settle\n");
	    failures++;
	}
	failures += aet_verify_shared_autoview(views, sentinels, fitted);
	if (failures && !soft_fail)
	    break;

	failures += aet_capture_views(gedp, views, fitted,
	    AET_CAPTURE_INITIAL_512, AET_INITIAL_DIMENSION, captures);
	if (failures && !soft_fail)
	    break;

	failures += aet_resize_views(gedp, views, AET_LARGE_DIMENSION,
	    fitted, "initial resize to 600");
	failures += aet_capture_views(gedp, views, fitted,
	    AET_CAPTURE_INITIAL_600, AET_LARGE_DIMENSION, captures);
	if (failures && !soft_fail)
	    break;

	failures += aet_resize_views(gedp, views, AET_INITIAL_DIMENSION,
	    fitted, "shrink to 512");
	failures += aet_capture_views(gedp, views, fitted,
	    AET_CAPTURE_SHRINK_512, AET_INITIAL_DIMENSION, captures);
	if (failures && !soft_fail)
	    break;

	bu_log("Cycling all views through dimensions 513 to 599...\n");
	for (int dimension = AET_CYCLE_FIRST_DIMENSION;
		dimension <= AET_CYCLE_LAST_DIMENSION; dimension++) {
	    char phase[64] = {0};
	    snprintf(phase, sizeof(phase), "resize cycle at %d", dimension);
	    failures += aet_resize_views(gedp, views, dimension, fitted,
		phase);
	    if (failures && !soft_fail)
		break;
	}
	if (failures && !soft_fail)
	    break;

	failures += aet_resize_views(gedp, views, AET_LARGE_DIMENSION,
	    fitted, "post-cycle resize to 600");
	failures += aet_capture_views(gedp, views, fitted,
	    AET_CAPTURE_POST_CYCLE_600, AET_LARGE_DIMENSION, captures);
	if (failures && !soft_fail)
	    break;

	failures += aet_resize_views(gedp, views, AET_INITIAL_DIMENSION,
	    fitted, "final resize to 512");
	failures += aet_capture_views(gedp, views, fitted,
	    AET_CAPTURE_FINAL_512, AET_INITIAL_DIMENSION, captures);
	if (failures && !soft_fail)
	    break;

	for (size_t i = 0; i < AET_VIEW_COUNT; i++) {
	    if (!aet_images_match_exactly(captures[AET_CAPTURE_INITIAL_512][i],
		    captures[AET_CAPTURE_SHRINK_512][i]))
		failures++;
	    if (!aet_images_match_exactly(captures[AET_CAPTURE_INITIAL_600][i],
		    captures[AET_CAPTURE_POST_CYCLE_600][i]))
		failures++;
	    if (!aet_images_match_exactly(captures[AET_CAPTURE_INITIAL_512][i],
		    captures[AET_CAPTURE_FINAL_512][i]))
		failures++;
	}
    } while (0);

    ged_close(gedp);
    bu_file_delete("moss_aet_tmp.g");

    if (!failures && !keep_images)
	aet_capture_cleanup(captures);
    else if (failures)
	bu_log("AET test failed; preserving same-run captures for inspection\n");

    if (!failures)
	bu_log("PASS: shared-view AET and resize state remained stable\n");
    return failures ? BRLCAD_ERROR : BRLCAD_OK;
}


// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
