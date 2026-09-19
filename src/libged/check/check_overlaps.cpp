/*            C H E C K _ O V E R L A P S . C P P
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

#include "common.h"
#include "ged.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "BObol/BViewStore.h"
#include "bg/line_layer.h"
#include "ged/draw.h"

#include "../ged_bobol_private.hpp"
#include "../ged_private.h"
#include "./check_private.h"

struct ged_check_plot_segment {
    char *region_a;
    char *region_b;
    double depth_mm;
    point_t entry_mm;
    point_t exit_mm;
};

struct ged_check_plot {
    struct bg_line_layer_builder *builder;
    struct ged_check_plot_segment *segments;
    size_t segment_count;
    size_t segment_capacity;
};

static int
ged_check_plot_append(struct ged_check_plot *plot,
    const char *region_a,
    const char *region_b,
    double depth_mm,
    const point_t a,
    const point_t b)
{
    if (!plot)
	return 0;

    if (!plot->builder)
	return 0;

    if (!bg_line_layer_builder_add(plot->builder, 255, 255, 0, a, BG_GEOMETRY_LINE_MOVE))
	return 0;
    if (!bg_line_layer_builder_add(plot->builder, 255, 255, 0, b, BG_GEOMETRY_LINE_DRAW))
	return 0;

    if (plot->segment_count >= plot->segment_capacity) {
	size_t capacity = plot->segment_capacity ?
	    plot->segment_capacity * 2 : 64;
	plot->segments = (struct ged_check_plot_segment *)bu_realloc(
	    plot->segments, capacity * sizeof(struct ged_check_plot_segment),
	    "check overlap plot segments");
	memset(plot->segments + plot->segment_capacity, 0,
	    (capacity - plot->segment_capacity) *
	    sizeof(struct ged_check_plot_segment));
	plot->segment_capacity = capacity;
    }

    struct ged_check_plot_segment *segment =
	&plot->segments[plot->segment_count++];
    segment->region_a = bu_strdup(region_a ? region_a : "");
    segment->region_b = bu_strdup(region_b ? region_b : "");
    segment->depth_mm = depth_mm;
    VMOVE(segment->entry_mm, a);
    VMOVE(segment->exit_mm, b);
    return 1;
}

static void
ged_check_plot_free(struct ged_check_plot *plot)
{
    if (!plot)
	return;

    if (plot->builder)
	bg_line_layer_builder_free(plot->builder);
    for (size_t i = 0; i < plot->segment_count; i++) {
	bu_free(plot->segments[i].region_a, "check overlap region a");
	bu_free(plot->segments[i].region_b, "check overlap region b");
    }
    if (plot->segments)
	bu_free(plot->segments, "check overlap plot segments");
    memset(plot, 0, sizeof(*plot));
}


struct overlap_list {
    struct bu_list l;
    char *reg1; 	/* overlapping region 1 */
    char *reg2; 	/* overlapping region 2 */
    size_t count;	/* number of time reported */
    double maxdepth;	/* maximum overlap depth */
};


struct overlaps_context {
    struct overlap_list *overlapList;
    size_t noverlaps;		/* Number of overlaps seen */
    size_t overlap_count;	/* Number of overlap pairs seen */
    size_t unique_overlap_count; /* Number of unique overlap pairs seen */
    int overlaps_overlay_flag;
    int rpt_overlap_flag;
    FILE *plot_overlaps;
    int overlap_color[3];
    struct ged_check_plot *overlaps_overlay_plot;
    int sem_stats;
    int sem_lists;
    struct ged *gedp;
};

static BObolFeatureOwner
check_feature_owner(void)
{
    BObolFeatureOwner owner;
    owner.ownerId = "check";
    owner.ownerRole = "command-result";
    return owner;
}

static BObolOverlayInfo
check_overlay_info(struct ged_view_context *view_ctx, const char *name)
{
    BObolOverlayInfo overlay;
    overlay.isOverlay = TRUE;
    overlay.ownerToken = view_ctx;
    overlay.role = BObolOverlayRole::Model;
    overlay.overlayClass = BObolOverlayClass::CommandResult;
    overlay.lifecycle = BObolOverlayLifecycle::PerCommand;
    overlay.order = BObolOverlayOrder::PostTransparent;
    overlay.sourcePath = name ? name : "";
    return overlay;
}

static BObolFeatureMetadata
check_metadata(const char *key, const char *value)
{
    BObolFeatureMetadata item;
    item.key = key ? key : "";
    item.value = value ? value : "";
    return item;
}

static void
check_publish_overlap_label(struct ged *gedp, struct ged_view_context *view_ctx,
	const struct overlaps_context *context)
{
    if (!gedp || !context)
	return;

    struct bu_vls text = BU_VLS_INIT_ZERO;
    bu_vls_printf(&text, "check overlaps: %zu total, %zu unique",
	    context->noverlaps, context->unique_overlap_count);

    struct ged_diagnostic_hud_label label = GED_DIAGNOSTIC_HUD_LABEL_INIT;
    label.label_id = "check::overlaps::summary";
    label.text = bu_vls_cstr(&text);
    label.position[0] = 8.0;
    label.position[1] = 42.0;
    label.color[0] = 255;
    label.color[1] = 255;
    label.color[2] = 0;
    label.font_size = 12.0;
    label.source_id = (uint32_t)context->noverlaps;

    BObolViewController *controller =
	ged_bobol_shared_view_controller(gedp);
    if (view_ctx && controller) {
	BObolFeatureOwner owner = check_feature_owner();
	controller->features().markCommandOwnerGeneration(owner);

	BObolLabel item;
	item.text = label.text;
	item.point = SbVec3f(static_cast<float>(label.position[0]),
	    static_cast<float>(label.position[1]), 0.0f);
	item.hasColor = TRUE;
	item.color = SbColor(1.0f, 1.0f, 0.0f);
	item.fontSize = static_cast<float>(label.font_size);
	item.sourceId = label.source_id;
	std::vector<BObolLabel> labels(1, item);
	BObolFeatureStyle style;
	style.hasVisible = TRUE;
	style.visible = TRUE;
	style.hasSelectable = TRUE;
	style.selectable = TRUE;
	BObolFeatureHandle handle = controller->features().publishHudLabels(
	    label.label_id, BObolFeatureScope::Shared, labels, &style, &owner);

	char total[64] = {0};
	char unique[64] = {0};
	snprintf(total, sizeof(total), "%zu", context->noverlaps);
	snprintf(unique, sizeof(unique), "%zu",
	    context->unique_overlap_count);
	std::vector<BObolFeatureMetadata> metadata;
	metadata.push_back(check_metadata("result.feature", label.label_id));
	metadata.push_back(check_metadata("result.owner", "check"));
	metadata.push_back(check_metadata("result.kind", "overlap-summary"));
	metadata.push_back(check_metadata("result.overlap_count", total));
	metadata.push_back(check_metadata("result.unique_overlap_count", unique));
	metadata.push_back(check_metadata("result.schema",
	    "brlcad.check.overlap-summary.v1"));
	metadata.push_back(check_metadata("result.severity",
	    context->noverlaps ? "error" : "info"));
	if (handle.isValid() &&
	    controller->features().replaceMetadata(handle, metadata)) {
	    (void)controller->features().setOverlayInfo(handle,
		check_overlay_info(view_ctx, label.label_id));
	    bu_vls_free(&text);
	    return;
	}
    }

    if (ged_diagnostic_hud_label_handler_available(gedp))
	(void)ged_diagnostic_hud_label_publish(gedp, &label);

    bu_vls_free(&text);
}

static int
check_publish_overlap_lines(struct ged *gedp, struct ged_view_context *view_ctx,
	const struct ged_check_plot *plot)
{
    const struct bg_line_layer_builder *builder = plot ? plot->builder : NULL;
    if (!gedp || !builder)
	return 0;

    BObolViewController *controller =
	ged_bobol_shared_view_controller(gedp);
    if (view_ctx && controller) {
	BObolFeatureOwner owner = check_feature_owner();
	controller->features().markCommandOwnerGeneration(owner);
	(void)controller->features().removePrefix("check::overlaps",
	    BOBOL_FEATURE_SCOPE_SHARED, &owner);
	if (!bg_line_layer_builder_point_count(builder))
	    return 1;

	BObolFeatureStyle style;
	style.hasVisible = TRUE;
	style.visible = TRUE;
	style.hasSelectable = TRUE;
	style.selectable = TRUE;
	BObolFeatureHandle handle =
	    controller->features().publishLineLayerBuilder(
		"check::overlaps", BObolFeatureScope::Shared,
		builder, &style, &owner);
	if (handle.isValid()) {
	    char layer_count[64] = {0};
	    char point_count[64] = {0};
	    char primitive_count[64] = {0};
	    snprintf(layer_count, sizeof(layer_count), "%zu",
		    bg_line_layer_builder_layer_count(builder));
	    snprintf(point_count, sizeof(point_count), "%zu",
		    bg_line_layer_builder_point_count(builder));
	    snprintf(primitive_count, sizeof(primitive_count), "%zu",
		    plot->segment_count);
	    std::vector<BObolFeatureMetadata> metadata;
	    metadata.push_back(check_metadata("result.feature", "check::overlaps"));
	    metadata.push_back(check_metadata("result.format", "line-layer-builder"));
	    metadata.push_back(check_metadata("result.layer_count", layer_count));
	    metadata.push_back(check_metadata("result.point_count", point_count));
	    metadata.push_back(check_metadata("result.owner", "check"));
	    metadata.push_back(check_metadata("result.kind", "overlap"));
	    metadata.push_back(check_metadata("result.schema", "brlcad.check.overlap.v1"));
	    metadata.push_back(check_metadata("result.severity", "error"));
	    metadata.push_back(check_metadata("result.primitive_count", primitive_count));
	    metadata.push_back(check_metadata("result.units", "mm"));
	    bool metadata_published =
		controller->features().replaceMetadata(handle, metadata);
	    for (size_t i = 0; metadata_published &&
		    i < plot->segment_count; i++) {
		const struct ged_check_plot_segment *segment =
		    &plot->segments[i];
		char primitive[64] = {0};
		char depth[64] = {0};
		char entry[192] = {0};
		char exit[192] = {0};
		snprintf(primitive, sizeof(primitive), "%zu", i);
		snprintf(depth, sizeof(depth), "%.17g", segment->depth_mm);
		snprintf(entry, sizeof(entry), "%.17g %.17g %.17g",
		    V3ARGS(segment->entry_mm));
		snprintf(exit, sizeof(exit), "%.17g %.17g %.17g",
		    V3ARGS(segment->exit_mm));
		std::vector<BObolFeatureMetadata> primitive_metadata;
		primitive_metadata.push_back(check_metadata("result.schema", "brlcad.check.overlap.v1"));
		primitive_metadata.push_back(check_metadata("result.primitive", primitive));
		primitive_metadata.push_back(check_metadata("result.primitive.kind", "overlap"));
		primitive_metadata.push_back(check_metadata("result.severity", "error"));
		primitive_metadata.push_back(check_metadata("overlap.region_a", segment->region_a));
		primitive_metadata.push_back(check_metadata("overlap.region_b", segment->region_b));
		primitive_metadata.push_back(check_metadata("overlap.depth_mm", depth));
		primitive_metadata.push_back(check_metadata("hit.entry_mm", entry));
		primitive_metadata.push_back(check_metadata("hit.exit_mm", exit));
		primitive_metadata.push_back(check_metadata("result.units", "mm"));
		metadata_published = controller->features().replacePrimitiveMetadata(
		    handle, static_cast<int32_t>(i), primitive_metadata);
	    }
	    if (metadata_published) {
		(void)controller->features().setOverlayInfo(handle,
		    check_overlay_info(view_ctx, "check::overlaps"));
		return 1;
	    }
	}
    }

    int handled = ged_diagnostic_line_layer_publish(gedp,
	    "check::overlaps", builder);
    return handled;
}

static void
check_log_overlaps(struct ged *gedp, const char *reg1, const char *reg2, double depth, vect_t ihit, vect_t ohit, void *context)
{
    struct overlaps_context *callbackdata = (struct overlaps_context*)context;
    struct overlap_list *olist= (struct overlap_list *)callbackdata->overlapList;
    struct overlap_list *new_op;
    struct overlap_list *op;

    bu_semaphore_acquire(callbackdata->sem_stats);
    callbackdata->noverlaps++;
    bu_semaphore_release(callbackdata->sem_stats);

    if (!callbackdata->rpt_overlap_flag) {
	bu_vls_printf(gedp->ged_result_str,
		      "OVERLAP %zu: %s\nOVERLAP %zu: %s\nOVERLAP %zu: depth %gmm\nOVERLAP %zu: in_hit_point (%g, %g, %g) mm\nOVERLAP %zu: out_hit_point (%g, %g, %g) mm\n------------------------------------------------------------\n",
		      callbackdata->noverlaps, reg1,
		      callbackdata->noverlaps, reg2,
		      callbackdata->noverlaps, depth,
		      callbackdata->noverlaps, V3ARGS(ihit),
		      callbackdata->noverlaps, V3ARGS(ohit));
	/* If we report overlaps, don't print if already noted once.
	 * Build up a linked list of known overlapping regions and compare
	 * against it.
	 */
    } else {
	BU_GET(new_op, struct overlap_list);
	bu_semaphore_acquire(callbackdata->sem_stats);
	for (BU_LIST_FOR(op, overlap_list, &(olist->l))) {
	    if ((BU_STR_EQUAL(reg1, op->reg1)) && (BU_STR_EQUAL(reg2, op->reg2))) {
		op->count++;
		if (depth > op->maxdepth)
		    op->maxdepth = depth;
		bu_semaphore_release(callbackdata->sem_stats);
		bu_free((char *) new_op, "overlap list");
		return;
	    }
	}

	for (BU_LIST_FOR(op, overlap_list, &(olist->l))) {
	    /* if this pair was seen in reverse, decrease the unique counter */
	    if ((BU_STR_EQUAL(reg1, op->reg2)) && (BU_STR_EQUAL(reg2, op->reg1))) {
		callbackdata->unique_overlap_count--;
		break;
	    }
	}

	/* we have a new overlapping region pair */
	callbackdata->overlap_count++;
	callbackdata->unique_overlap_count++;
	new_op->reg1 = (char *)bu_malloc(strlen(reg1)+1, "reg1");
	new_op->reg2 = (char *)bu_malloc(strlen(reg2)+1, "reg2");
	bu_strlcpy(new_op->reg1, reg1, strlen(reg1)+1);
	bu_strlcpy(new_op->reg2, reg2, strlen(reg2)+1);
	new_op->maxdepth = depth;
	new_op->count = 1;
	BU_LIST_INSERT(&(olist->l), &(new_op->l));
	bu_semaphore_release(callbackdata->sem_stats);
    }
}


static void
printOverlaps(struct ged *gedp, void *context, struct check_parameters *options)
{
    struct overlaps_context *overlapData = (struct overlaps_context*)context;
    struct overlap_list *olist= (struct overlap_list *)overlapData->overlapList;
    struct overlap_list *op;
    struct overlap_list *backop;
    struct overlap_list *nextop;
    size_t object_counter = 0;

    if (overlapData->rpt_overlap_flag) {
	struct bu_vls str = BU_VLS_INIT_ZERO;
	/* using counters instead of the actual variables to be able to
	 * summarize after checking for matching pairs
	 */
	size_t overlap_counter = overlapData->overlap_count;
	size_t unique_overlap_counter = overlapData->unique_overlap_count;

	/* iterate over the overlap pairs and output one OVERLAP section
	 * per unordered pair.  a summary is output at the end.
	 */
	bu_vls_printf(&str, "OVERLAP PAIRS\n------------------------------------------\n");
	for (BU_LIST_FOR(op, overlap_list, &(olist->l))) {

	    for (BU_LIST_FOR_BACKWARDS(backop, overlap_list, &(op->l))) {
		if ((BU_STR_EQUAL(op->reg2, backop->reg1)) && (BU_STR_EQUAL(op->reg1, backop->reg2))) break;
		if (BU_LIST_IS_HEAD(&(backop->l), &(olist->l))) break;
	    }
	    if (backop && BU_LIST_NOT_HEAD(&(backop->l), &(olist->l))) continue;

	    bu_vls_printf(&str, "%s and %s overlap\n", op->reg1, op->reg2);

	    nextop = (struct overlap_list *)NULL;
	    /* if there are still matching pairs to search for */
	    if (overlap_counter > unique_overlap_counter) {
		/* iterate until end of pairs or we find a
		 * reverse matching pair (done inside loop
		 * explicitly)*/
		BU_LIST_TAIL(&(olist->l), op, nextop, struct overlap_list) {
		    if ((BU_STR_EQUAL(op->reg1, nextop->reg2)) && (BU_STR_EQUAL(op->reg2, nextop->reg1)))
			break;
		}
		/* when we leave the loop, nextop is either
		 * null (hit tail of list) or the matching
		 * reverse pair */
	    }

	    bu_vls_printf(&str, "\t<%s, %s>: %zu overlap%c detected, maximum depth is %g %s\n",
			  op->reg1, op->reg2, op->count, op->count>1 ? 's' : (char) 0, op->maxdepth/options->units[LINE]->val, options->units[LINE]->name);
	    if (nextop && BU_LIST_NOT_HEAD(nextop, &(olist->l))) {
		    bu_vls_printf(&str,
				  "\t<%s, %s>: %zu overlap%c detected, maximum depth is %g %s\n",
				  nextop->reg1, nextop->reg2, nextop->count,
				  nextop->count > 1 ? 's' : (char)0, nextop->maxdepth/options->units[LINE]->val, options->units[LINE]->name);
		    /* counter the decrement below to account for
		 * the matched reverse pair
		 */
		    unique_overlap_counter++;
	    }

	    /* decrement so we may stop scanning for unique overlaps asap */
	    unique_overlap_counter--;
	    overlap_counter--;
	}

	if (overlapData->noverlaps) {
	    bu_vls_printf(&str, "==========================================\n");
	    bu_vls_printf(&str, "SUMMARY\n");

	    bu_vls_printf(&str, "\t%zu overlap%s detected\n",
			  overlapData->noverlaps, (overlapData->noverlaps == 1) ? "" : "s");
	    bu_vls_printf(&str, "\t%zu unique overlapping pair%s (%zu ordered pair%s)\n",
			  overlapData->unique_overlap_count,
			  (overlapData->unique_overlap_count == 1) ? "" : "s",
			  overlapData->overlap_count, (overlapData->overlap_count == 1) ? "" : "s");

	    if (olist) {
		bu_vls_printf(&str, "\tOverlapping objects: ");

		for (BU_LIST_FOR(op, overlap_list, &(olist->l))) {
		    /* iterate over the list and see if we already printed this one */
		    for (BU_LIST_FOR_BACKWARDS(backop, overlap_list, &(op->l))) {
			if ((BU_STR_EQUAL(op->reg1, backop->reg1)) || (BU_STR_EQUAL(op->reg1, backop->reg2)))
			    break;
			if (BU_LIST_IS_HEAD(&(backop->l), &(olist->l))) break;
		    }
		    /* if we got to the head of the list (backop points to the match) */
		    if (BU_LIST_IS_HEAD(&(backop->l), &(olist->l))) {
			bu_vls_printf(&str, "%s  ", op->reg1);
			object_counter++;
		    }

		   /* iterate over the list again up to where we are to see if the second
		    * region was already printed */
		    for (BU_LIST_FOR_BACKWARDS(backop, overlap_list, &(op->l))) {
			if ((BU_STR_EQUAL(op->reg2, backop->reg1)) || (BU_STR_EQUAL(op->reg2, backop->reg2)))
			    break;
			if (BU_LIST_IS_HEAD(&(backop->l), &(olist->l))) break;
		    }
		    /* if we got to the head of the list (backop points to the match) */
		    if (BU_LIST_IS_HEAD(&(backop->l), &(olist->l))) {
			bu_vls_printf(&str, "%s  ", op->reg2);
			object_counter++;
		    }
		}
		bu_vls_printf(&str, "\n\t%zu unique overlapping object%s detected\n",
			      object_counter, (object_counter == 1) ? "" : "s");
	    }
	} else {
		bu_vls_printf(&str, "%zu overlap%s detected\n\n",
			      overlapData->noverlaps, (overlapData->noverlaps==1)?"":"s");
	}
	bu_vls_vlscat(gedp->ged_result_str, &str);
	bu_vls_free(&str);
    }
}


static void
overlap(const struct xray *ray,
	const struct partition *pp,
	const struct region *reg1,
	const struct region *reg2,
	double depth,
	void* callback_data)
{
    struct overlaps_context *context = (struct overlaps_context*) callback_data;
    struct ged *gedp = context->gedp;
    struct hit *ihitp = pp->pt_inhit;
    struct hit *ohitp = pp->pt_outhit;
    vect_t ihit;
    vect_t ohit;

    VJOIN1(ihit, ray->r_pt, ihitp->hit_dist, ray->r_dir);
    VJOIN1(ohit, ray->r_pt, ohitp->hit_dist, ray->r_dir);

    if (context->overlaps_overlay_flag) {
	bu_semaphore_acquire(context->sem_stats);
	(void)ged_check_plot_append(context->overlaps_overlay_plot,
	    reg1->reg_name, reg2->reg_name, depth, ihit, ohit);
	bu_semaphore_release(context->sem_stats);
    }

    bu_semaphore_acquire(context->sem_lists);
    check_log_overlaps(gedp, reg1->reg_name, reg2->reg_name, depth, ihit, ohit, context);
    bu_semaphore_release(context->sem_lists);

    if (context->plot_overlaps) {
	pl_color(context->plot_overlaps, V3ARGS(context->overlap_color));
	pdv_3line(context->plot_overlaps, ihit, ohit);
    }
}


int check_overlaps(struct ged *gedp, struct current_state *state,
		   struct db_i *dbip,
		   char **tobjtab,
		   int tnobjs,
		   struct check_parameters *options)
{
    struct ged_check_plot check_plot;
    struct overlaps_context callbackdata;
    callbackdata.gedp = gedp;
    struct overlap_list overlapList;
    struct overlap_list *op;

    FILE *plot_overlaps = NULL;
    const char *name = "overlaps.plot3";
    int overlap_color[3] = { 255, 255, 0 };	/* yellow */
    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    int overlay_enabled = options->overlaps_overlay_flag &&
	(ged_bobol_shared_view_controller(gedp) ||
	 ged_diagnostic_line_layer_handler_available(gedp));

    /* init overlaps list */
    BU_LIST_INIT(&(overlapList.l));
    overlapList.count = 0;
    overlapList.reg1 = NULL;
    overlapList.reg2 = NULL;
    overlapList.maxdepth = 0;

    if (options->plot_files) {
	plot_overlaps = fopen(name, "wb");
	if (plot_overlaps == (FILE *)NULL) {
	    bu_vls_printf(gedp->ged_result_str, "cannot open plot file %s\n", name);
	    return BRLCAD_ERROR;
	}
    }

    if (options->overlaps_overlay_flag) {
	memset(&check_plot, 0, sizeof(check_plot));
	if (overlay_enabled)
	    check_plot.builder = bg_line_layer_builder_create();
	else
	    bu_vls_printf(gedp->ged_result_str, "overlap overlay requested, but no view is available\n");
    }

    callbackdata.noverlaps = 0;
    callbackdata.overlap_count = 0;
    callbackdata.unique_overlap_count = 0;
    callbackdata.rpt_overlap_flag = options->rpt_overlap_flag;
    VMOVE(callbackdata.overlap_color,overlap_color);
    callbackdata.plot_overlaps = plot_overlaps;
    callbackdata.overlapList = &overlapList;
    callbackdata.overlaps_overlay_flag = overlay_enabled;
    callbackdata.overlaps_overlay_plot = &check_plot;
    callbackdata.sem_stats = bu_semaphore_register("check_stats");
    callbackdata.sem_lists = bu_semaphore_register("check_lists");

    /* register callback */
    analyze_register_overlaps_callback(state, overlap, &callbackdata);

    if (perform_raytracing(state, dbip, tobjtab, tnobjs, ANALYSIS_OVERLAPS)) {
	while (BU_LIST_WHILE(op, overlap_list, &(overlapList.l))){
	    bu_free(op->reg1, "reg1 name");
	    bu_free(op->reg2, "reg1 name");

	    BU_LIST_DEQUEUE(&(op->l));
	    BU_PUT(op, struct overlap_list);
	}
	if (overlay_enabled)
	    ged_check_plot_free(&check_plot);
	return BRLCAD_ERROR;
    }

    print_verbose_debug(gedp, options);

    printOverlaps(gedp, &callbackdata, options);

    if (overlay_enabled) {
	(void)check_publish_overlap_lines(gedp, view_ctx, &check_plot);
	ged_check_plot_free(&check_plot);
    }
    if (options->overlaps_overlay_flag)
	check_publish_overlap_label(gedp, view_ctx, &callbackdata);

    if (plot_overlaps) {
	fclose(plot_overlaps);
	bu_vls_printf(gedp->ged_result_str, "\nplot file saved as %s",name);
    }

    while (BU_LIST_WHILE(op, overlap_list, &(overlapList.l))) {
	bu_free(op->reg1, "reg1 name");
	bu_free(op->reg2, "reg1 name");

	BU_LIST_DEQUEUE(&(op->l));
	BU_PUT(op, struct overlap_list);
    }

    return BRLCAD_OK;
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
