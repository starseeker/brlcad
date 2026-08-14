/*          Q G P R O G R E S S I V E D I A G N O S T I C S . C P P
 * BRL-CAD
 *
 * Copyright (c) 2014-2026 United States Government as represented by
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
/** @file QgProgressiveDiagnostics.cpp
 *
 * White-box progressive CAD diagnostics used by qged's GUI test driver.
 *
 */

#include "common.h"

#include <algorithm>
#include <climits>
#include <limits>
#include <unordered_set>
#include <vector>

#include <QJsonArray>
#include <QJsonObject>
#include <QRegularExpression>

#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoOrthographicCamera.h>

#include <Obol/cad/SoCADAssembly.h>

#include "BObol/BLodService.h"
#include "BObol/BDatabaseSource.h"
#include "BObol/BSourceRealization.h"
#include "BObol/BViewController.h"
#include "BObol/BViewLod.h"
#include "bu/vls.h"
#include "ged/scene.h"
#include "ged/selection.h"
#include "ged/view_feature.h"
#include "qtcad/QgCanvasBase.h"

#include "QgEdApp.h"
#include "QgGuiTestDriver.h"

static void
qged_test_collect_database_source_roots(
    SoNode *node, std::vector<SoBRLDatabaseSource *> &sources)
{
    if (!node)
	return;
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	SoBRLDatabaseSource *source =
	    static_cast<SoBRLDatabaseSource *>(node);
	if (std::find(sources.begin(), sources.end(), source) == sources.end())
	    sources.push_back(source);
	return;
    }
    if (!node->isOfType(SoGroup::getClassTypeId()))
	return;
    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++)
	qged_test_collect_database_source_roots(group->getChild(i), sources);
}

static QgView *
qged_test_event_view(QgEdApp &app, const QgTestEvent &event)
{
    QObject *target = app.w ?
	QgEventRecorder::resolveObject(app.w, event.target) : nullptr;
    for (QObject *object = target; object; object = object->parent()) {
	if (QgView *view = qobject_cast<QgView *>(object))
	    return view;
    }
    return app.w ? app.w->CurrentDisplay() : nullptr;
}

static BObolViewController *
qged_test_event_controller(QgEdApp &app, const QgTestEvent &event)
{
    QgView *view = qged_test_event_view(app, event);
    return view ? view->obolViewController() : nullptr;
}

QJsonObject
qged_collect_progressive_sample(QgEdApp &app, int eventIndex,
    const QgTestEvent &event, qint64 elapsedMilliseconds,
    qint64 eventMicroseconds, bool collectDeepLodDiagnostics)
{
    QJsonObject sample;
    sample.insert(QStringLiteral("event_index"), eventIndex);
    sample.insert(QStringLiteral("action"), event.action);
    sample.insert(QStringLiteral("target"), event.target);
    if (event.action == QLatin1String("qged_command"))
	sample.insert(QStringLiteral("command"),
	    event.arguments.value(QStringLiteral("command")));
    if (event.action == QLatin1String("qged_command_batch"))
	sample.insert(QStringLiteral("commands"),
	    event.arguments.value(QStringLiteral("commands")));
    if (event.action == QLatin1String("checkpoint"))
	sample.insert(QStringLiteral("checkpoint"),
	    event.arguments.value(QStringLiteral("name")));
    if (event.action == QLatin1String("resize")) {
	sample.insert(QStringLiteral("requested_width"),
	    event.arguments.value(QStringLiteral("width")));
	sample.insert(QStringLiteral("requested_height"),
	    event.arguments.value(QStringLiteral("height")));
    }
    if (event.action == QLatin1String("window_state"))
	sample.insert(QStringLiteral("requested_window_state"),
	    event.arguments.value(QStringLiteral("state")));
    sample.insert(QStringLiteral("elapsed_ms"), elapsedMilliseconds);
    sample.insert(QStringLiteral("event_duration_us"), eventMicroseconds);
    if (app.w) {
	sample.insert(QStringLiteral("window_width"), app.w->width());
	sample.insert(QStringLiteral("window_height"), app.w->height());
	sample.insert(QStringLiteral("window_minimized"),
	    app.w->isMinimized());
	sample.insert(QStringLiteral("window_maximized"),
	    app.w->isMaximized());
	sample.insert(QStringLiteral("window_fullscreen"),
	    app.w->isFullScreen());
	sample.insert(QStringLiteral("window_visible"), app.w->isVisible());
    }
    if (QgView *view = qged_test_event_view(app, event)) {
	sample.insert(QStringLiteral("view_width"), view->width());
	sample.insert(QStringLiteral("view_height"), view->height());
	struct ged_view_context *sampleViewContext =
	    ged_view_context_from_bv(view->viewContext());
	if (sampleViewContext) {
	    sample.insert(QStringLiteral("lod_progress_track_present"),
		ged_view_feature_exists(sampleViewContext,
		    "_faceplate/lod_progress_track") != 0);
	    sample.insert(QStringLiteral("lod_progress_fill_present"),
		ged_view_feature_exists(sampleViewContext,
		    "_faceplate/lod_progress_fill") != 0);
	    sample.insert(QStringLiteral("lod_progress_label_present"),
		ged_view_feature_exists(sampleViewContext,
		    "_faceplate/lod_progress_label") != 0);
	}
	if (QgCanvasBase *canvas = view->canvasBase()) {
	    if (QWidget *widget = canvas->canvasWidget()) {
		sample.insert(QStringLiteral("canvas_width"), widget->width());
		sample.insert(QStringLiteral("canvas_height"), widget->height());
		sample.insert(QStringLiteral("canvas_device_pixel_ratio"),
		    static_cast<double>(widget->devicePixelRatioF()));
	    }
	}
    }
    /* Full draw-frontier and selection listings cross the GED/Obol hierarchy
     * and are structural diagnostics, not frame telemetry.  Sampling them
     * after every 8 ms motion interval made the validation harness dominate
     * its own perf captures. */
    const bool collectStructuralDiagnostics =
	collectDeepLodDiagnostics ||
	event.action == QLatin1String("qged_command") ||
	event.action == QLatin1String("qged_command_batch") ||
	event.action == QLatin1String("set_current") ||
	event.action == QLatin1String("clear_selection");

    QgModel *model = app.mdl;
    sample.insert(QStringLiteral("tree_model_available"), model != nullptr);
    if (model) {
	sample.insert(QStringLiteral("tree_loaded_items"),
	    static_cast<qint64>(model->allItems().size()));
	sample.insert(QStringLiteral("tree_top_items"),
	    static_cast<qint64>(model->topItems().size()));
	sample.insert(QStringLiteral("tree_top_rows"), model->rowCount());
	sample.insert(QStringLiteral("tree_flattened"),
	    model->flattenHierarchy());

	const QgModel::NotificationStats notification =
	    model->notificationStats();
	sample.insert(QStringLiteral("tree_notify_path_queries"),
	    static_cast<qint64>(notification.path_queries));
	sample.insert(QStringLiteral("tree_notify_path_candidates"),
	    static_cast<qint64>(notification.path_candidates));
	sample.insert(QStringLiteral("tree_notify_fallback_scans"),
	    static_cast<qint64>(notification.path_fallback_scans));
	sample.insert(QStringLiteral("tree_notify_items"),
	    static_cast<qint64>(notification.items_notified));
	sample.insert(QStringLiteral("tree_notify_full_items"),
	    static_cast<qint64>(notification.full_items_notified));
	sample.insert(QStringLiteral("tree_notify_subtree_items"),
	    static_cast<qint64>(notification.subtree_items_notified));
	sample.insert(QStringLiteral("tree_notify_path_us"),
	    static_cast<qint64>(notification.path_notify_us));

	const QgModel::DrawTimingStats draw = model->drawTimingStats();
	sample.insert(QStringLiteral("tree_draw_calls"),
	    static_cast<qint64>(draw.draw_calls));
	sample.insert(QStringLiteral("tree_draw_successes"),
	    static_cast<qint64>(draw.successful_draw_calls));
	sample.insert(QStringLiteral("tree_draw_blank_slate_us"),
	    static_cast<qint64>(draw.blank_slate_check_us));
	sample.insert(QStringLiteral("tree_draw_transaction_us"),
	    static_cast<qint64>(draw.transaction_us));
	sample.insert(QStringLiteral("tree_draw_observer_us"),
	    static_cast<qint64>(draw.observer_callback_us));
	sample.insert(QStringLiteral("tree_draw_notify_path_us"),
	    static_cast<qint64>(draw.notify_path_us));
	sample.insert(QStringLiteral("tree_draw_notify_all_us"),
	    static_cast<qint64>(draw.notify_all_us));
	sample.insert(QStringLiteral("tree_draw_view_signal_us"),
	    static_cast<qint64>(draw.view_signal_us));

	const QgModel::FetchMoreStats fetch = model->fetchMoreStats();
	sample.insert(QStringLiteral("tree_fetch_calls"),
	    static_cast<qint64>(fetch.calls));
	sample.insert(QStringLiteral("tree_fetch_populated_calls"),
	    static_cast<qint64>(fetch.populated_calls));
	sample.insert(QStringLiteral("tree_fetch_inserted_rows"),
	    static_cast<qint64>(fetch.inserted_rows));
	sample.insert(QStringLiteral("tree_fetch_inserted_ranges"),
	    static_cast<qint64>(fetch.inserted_row_ranges));
	sample.insert(QStringLiteral("tree_fetch_max_children"),
	    static_cast<qint64>(fetch.max_child_count));
	sample.insert(QStringLiteral("tree_fetch_elapsed_us"),
	    static_cast<qint64>(fetch.elapsed_us));

	const QgModel::HierarchyDeltaStats hierarchy =
	    model->hierarchyDeltaStats();
	sample.insert(QStringLiteral("tree_delta_attempts"),
	    static_cast<qint64>(hierarchy.attempts));
	sample.insert(QStringLiteral("tree_delta_applied"),
	    static_cast<qint64>(hierarchy.applied));
	sample.insert(QStringLiteral("tree_delta_reset_fallbacks"),
	    static_cast<qint64>(hierarchy.reset_fallbacks));
	sample.insert(QStringLiteral("tree_delta_changed_parents"),
	    static_cast<qint64>(hierarchy.changed_parents));
	sample.insert(QStringLiteral("tree_delta_inserted_rows"),
	    static_cast<qint64>(hierarchy.inserted_rows));
	sample.insert(QStringLiteral("tree_delta_removed_rows"),
	    static_cast<qint64>(hierarchy.removed_rows));
	sample.insert(QStringLiteral("tree_delta_moved_rows"),
	    static_cast<qint64>(hierarchy.moved_rows));
	sample.insert(QStringLiteral("tree_delta_elapsed_us"),
	    static_cast<qint64>(hierarchy.elapsed_us));

	struct ged *gedp = model->ged();
	sample.insert(QStringLiteral("selection_state_available"),
	    ged_selection_state_available(gedp) != 0);
	if (gedp && ged_selection_state_available(gedp)) {
	    const size_t selectionCount = ged_selection_count(gedp, nullptr);
	    sample.insert(QStringLiteral("selection_count"),
		static_cast<qint64>(selectionCount));
	    sample.insert(QStringLiteral("selection_hash"),
		QString::number(ged_selection_state_hash(gedp, nullptr)));
	    if (collectStructuralDiagnostics) {
		struct bu_vls selectionListing = BU_VLS_INIT_ZERO;
		(void)ged_selection_list_paths(gedp, nullptr, &selectionListing);
		const QStringList selectionPaths =
		    QString::fromLocal8Bit(bu_vls_cstr(&selectionListing))
		    .split(QRegularExpression(QStringLiteral("[\\r\\n]+")),
			Qt::SkipEmptyParts);
		QJsonArray selected;
		const qsizetype selectedLimit =
		    std::min<qsizetype>(selectionPaths.size(), 64);
		for (qsizetype i = 0; i < selectedLimit; ++i)
		    selected.append(selectionPaths[i]);
		sample.insert(QStringLiteral("selection_paths"), selected);
		sample.insert(QStringLiteral("selection_paths_truncated"),
		    selectionPaths.size() > selectedLimit);
		bu_vls_free(&selectionListing);
	    }
	}

	if (gedp && ged_scene_available(gedp)) {
	    sample.insert(QStringLiteral("draw_scene_revision"),
		QString::number(ged_scene_revision(gedp)));
	    if (collectStructuralDiagnostics) {
		sample.insert(QStringLiteral("draw_occurrence_count"),
		    static_cast<qint64>(ged_scene_occurrence_count(gedp)));
		struct bu_vls drawListing = BU_VLS_INIT_ZERO;
		const size_t drawCount = ged_scene_paths_append(gedp,
		    ged_view_active_ctx(gedp), GED_SCENE_DRAW_DEFAULT,
		    GED_SCENE_PATHS_DRAW_INTENTS, &drawListing);
		const QStringList drawPaths =
		    QString::fromLocal8Bit(bu_vls_cstr(&drawListing))
		    .split(QRegularExpression(QStringLiteral("[\\r\\n]+")),
			Qt::SkipEmptyParts);
		QJsonArray frontier;
		const qsizetype frontierLimit =
		    std::min<qsizetype>(drawPaths.size(), 64);
		for (qsizetype i = 0; i < frontierLimit; ++i)
		    frontier.append(drawPaths[i]);
		sample.insert(QStringLiteral("draw_frontier_count"),
		    static_cast<qint64>(drawCount));
		sample.insert(QStringLiteral("draw_frontier_paths"), frontier);
		sample.insert(QStringLiteral("draw_frontier_paths_truncated"),
		    drawPaths.size() > frontierLimit);
		bu_vls_free(&drawListing);
	    }
	}

	if (QgView *sampleView = qged_test_event_view(app, event)) {
	    ged_view_lod_policy policy;
	    if (ged_view_lod_policy_get(&policy,
		    ged_view_context_from_bv(sampleView->viewContext()))) {
		sample.insert(QStringLiteral("view_lod_policy"), policy.policy);
		sample.insert(QStringLiteral("view_lod_mesh_enabled"),
		    policy.mesh_enabled != 0);
		sample.insert(QStringLiteral("view_lod_csg_enabled"),
		    policy.csg_enabled != 0);
		sample.insert(QStringLiteral("view_lod_zoom_refresh"),
		    policy.zoom_refresh != 0);
		sample.insert(QStringLiteral("view_lod_scale"),
		    static_cast<double>(policy.scale));
	    }
	}
    }

    BObolViewController *controller =
	qged_test_event_controller(app, event);
    if (!controller) {
	sample.insert(QStringLiteral("controller_available"), false);
	return sample;
    }

    sample.insert(QStringLiteral("controller_available"), true);
    sample.insert(QStringLiteral("render_request_serial"),
	static_cast<qint64>(controller->getRenderRequestSerial()));
    sample.insert(QStringLiteral("render_completion_serial"),
	static_cast<qint64>(controller->getRenderCompletionSerial()));
    sample.insert(QStringLiteral("presented_frame_serial"),
	static_cast<qint64>(controller->getPresentedFrameSerial()));
    sample.insert(QStringLiteral("lod_settle_after_render_serial"),
	static_cast<qint64>(controller->getLodSettleAfterRenderSerial()));
    sample.insert(QStringLiteral("lod_refinement_resume_after_render_serial"),
	static_cast<qint64>(
	    controller->getLodRefinementResumeAfterRenderSerial()));
    struct bv_view_info viewInfo = BV_VIEW_INFO_INIT;
    if (controller->getViewInfo(&viewInfo)) {
	sample.insert(QStringLiteral("model_view_size"),
	    static_cast<double>(viewInfo.size));
	sample.insert(QStringLiteral("model_view_width"), viewInfo.width);
	sample.insert(QStringLiteral("model_view_height"), viewInfo.height);
    }
    const BObolViewController::LightingProfile lightingProfile =
	controller->getLightingProfile();
    sample.insert(QStringLiteral("lighting_profile"),
	lightingProfile == BObolViewController::LIGHTING_MGED ?
	    QStringLiteral("mged") : QStringLiteral("studio"));
    sample.insert(QStringLiteral("lighting_ambient_intensity"),
	static_cast<double>(controller->getLightingAmbientIntensity()));
    std::vector<BObolSceneLightRealization> cameraLights;
    controller->getCameraLights(cameraLights);
    sample.insert(QStringLiteral("lighting_camera_light_count"),
	static_cast<qint64>(cameraLights.size()));
    BObolSourceRealizationCoordinator &sourceCoordinator =
	BObolSourceRealizationCoordinator::global();
    sample.insert(QStringLiteral("source_realization_workers"),
	static_cast<qint64>(
	    sourceCoordinator.workerCountForDiagnostics()));
    sample.insert(QStringLiteral("source_realization_queued_items"),
	static_cast<qint64>(
	    sourceCoordinator.queuedItemCountForDiagnostics()));
    sample.insert(QStringLiteral("source_realization_active_items"),
	static_cast<qint64>(
	    sourceCoordinator.activeItemCountForDiagnostics()));
    sample.insert(QStringLiteral("source_realization_active_bytes"),
	static_cast<qint64>(
	    sourceCoordinator.activeWorkingSetBytesForDiagnostics()));
    sample.insert(QStringLiteral("source_realization_limit_bytes"),
	static_cast<qint64>(
	    sourceCoordinator.workingSetLimitBytesForDiagnostics()));
    if (SoCamera *camera = controller->getCamera()) {
	const SbVec3f position = camera->position.getValue();
	const SbRotation orientation = camera->orientation.getValue();
	SbVec3f orientationAxis;
	float orientationAngle = 0.0f;
	orientation.getValue(orientationAxis, orientationAngle);
	sample.insert(QStringLiteral("camera_position_x"),
	    static_cast<double>(position[0]));
	sample.insert(QStringLiteral("camera_position_y"),
	    static_cast<double>(position[1]));
	sample.insert(QStringLiteral("camera_position_z"),
	    static_cast<double>(position[2]));
	sample.insert(QStringLiteral("camera_focal_distance"),
	    static_cast<double>(camera->focalDistance.getValue()));
	sample.insert(QStringLiteral("camera_near_distance"),
	    static_cast<double>(camera->nearDistance.getValue()));
	sample.insert(QStringLiteral("camera_far_distance"),
	    static_cast<double>(camera->farDistance.getValue()));
	sample.insert(QStringLiteral("camera_orientation_axis_x"),
	    static_cast<double>(orientationAxis[0]));
	sample.insert(QStringLiteral("camera_orientation_axis_y"),
	    static_cast<double>(orientationAxis[1]));
	sample.insert(QStringLiteral("camera_orientation_axis_z"),
	    static_cast<double>(orientationAxis[2]));
	sample.insert(QStringLiteral("camera_orientation_angle"),
	    static_cast<double>(orientationAngle));
	const SbVec2s viewportSize =
	    controller->getViewportRegion().getViewportSizePixels();
	sample.insert(QStringLiteral("viewport_width"), viewportSize[0]);
	sample.insert(QStringLiteral("viewport_height"), viewportSize[1]);
	float aspect = viewportSize[1] > 0 ?
	    static_cast<float>(viewportSize[0]) /
		static_cast<float>(viewportSize[1]) : 1.0f;
	if (camera->isOfType(SoOrthographicCamera::getClassTypeId()))
	    sample.insert(QStringLiteral("camera_orthographic_height"),
		static_cast<double>(
		    static_cast<SoOrthographicCamera *>(camera)->
			height.getValue()));
	if (collectDeepLodDiagnostics) {
	    const SbMatrix viewProjection =
		camera->getViewVolume(aspect).getMatrix();
	    QJsonArray matrix;
	    const float *values = viewProjection[0];
	    for (size_t valueIndex = 0; valueIndex < 16; ++valueIndex)
		matrix.append(static_cast<double>(values[valueIndex]));
	    sample.insert(QStringLiteral("camera_view_projection"), matrix);
	}
    }
    const uint64_t lastRender = controller->getLastRenderTimeNanoseconds();
    const uint64_t smoothRender =
	controller->getSmoothedRenderTimeNanoseconds();
    const uint64_t smoothPresentation =
	controller->getSmoothedPresentationIntervalNanoseconds();
    const uint64_t smoothInteractivePresentation =
	controller->getSmoothedInteractivePresentationIntervalNanoseconds();
    const uint64_t displayedPresentation =
	controller->getDisplayedPresentationIntervalNanoseconds();
    sample.insert(QStringLiteral("last_render_ms"),
	static_cast<double>(lastRender) / 1000000.0);
    sample.insert(QStringLiteral("smoothed_render_ms"),
	static_cast<double>(smoothRender) / 1000000.0);
    sample.insert(QStringLiteral("last_background_render_ms"),
	static_cast<double>(
	    controller->getLastBackgroundRenderTimeNanoseconds()) /
	    1000000.0);
    sample.insert(QStringLiteral("last_scene_render_ms"),
	static_cast<double>(
	    controller->getLastSceneRenderTimeNanoseconds()) /
	    1000000.0);
    sample.insert(QStringLiteral("smoothed_present_ms"),
	static_cast<double>(smoothPresentation) / 1000000.0);
    sample.insert(QStringLiteral("smoothed_present_fps"),
	smoothPresentation ?
	1000000000.0 / static_cast<double>(smoothPresentation) : 0.0);
    sample.insert(QStringLiteral("smoothed_interactive_present_ms"),
	static_cast<double>(smoothInteractivePresentation) / 1000000.0);
    sample.insert(QStringLiteral("smoothed_interactive_present_fps"),
	smoothInteractivePresentation ?
	1000000000.0 /
	    static_cast<double>(smoothInteractivePresentation) : 0.0);
    sample.insert(QStringLiteral("displayed_present_ms"),
	static_cast<double>(displayedPresentation) / 1000000.0);
    sample.insert(QStringLiteral("displayed_present_fps"),
	displayedPresentation ?
	1000000000.0 / static_cast<double>(displayedPresentation) : 0.0);
    sample.insert(QStringLiteral("lod_interactive_target_fps"),
	static_cast<double>(controller->getLodInteractiveTargetFps()));
    sample.insert(QStringLiteral("lod_stable_target_fps"),
	static_cast<double>(controller->getLodStableTargetFps()));
    sample.insert(QStringLiteral("presentation_deadline_interactive_ms"),
	static_cast<double>(
	    controller->getInteractivePresentationFrameDeadline()) /
	    1000000.0);
    sample.insert(QStringLiteral("presentation_deadline_stable_ms"),
	static_cast<double>(controller->getStablePresentationFrameDeadline()) /
	    1000000.0);
    sample.insert(QStringLiteral("presentation_deadline_current_ms"),
	static_cast<double>(controller->getCurrentPresentationFrameDeadline()) /
	    1000000.0);
    sample.insert(QStringLiteral("presentation_interrupted_frames"),
	static_cast<qint64>(
	    controller->getInterruptedPresentationFrameCount()));
    sample.insert(QStringLiteral("presentation_last_interrupted_ms"),
	static_cast<double>(
	    controller->getLastInterruptedPresentationTimeNanoseconds()) /
	    1000000.0);
    sample.insert(QStringLiteral("active_lod_scene_faces"),
	static_cast<qint64>(controller->getActiveLodFaceCount()));
    sample.insert(QStringLiteral("active_lod_scene_render_cost"),
	static_cast<qint64>(controller->getActiveLodRenderCost()));
    sample.insert(QStringLiteral("lod_scene_render_cost_budget"),
	static_cast<qint64>(controller->getCurrentLodRenderCostBudget()));
    sample.insert(QStringLiteral("lod_calibrated_render_cost_per_second"),
	controller->getCalibratedLodRenderCostPerSecond());
    sample.insert(QStringLiteral(
	"lod_interactive_calibrated_render_cost_per_second"),
	controller->getInteractiveCalibratedLodRenderCostPerSecond());
    sample.insert(QStringLiteral(
	"lod_stable_calibrated_render_cost_per_second"),
	controller->getStableCalibratedLodRenderCostPerSecond());
    sample.insert(QStringLiteral("last_progressive_advance_ms"),
	static_cast<double>(
	    controller->getLastProgressiveAdvanceTimeNanoseconds()) /
	    1000000.0);
    sample.insert(QStringLiteral("last_lod_result_processing_ms"),
	static_cast<double>(
	    controller->getLastLodResultProcessingTimeNanoseconds()) /
	    1000000.0);
    sample.insert(QStringLiteral("last_progressive_provider_ms"),
	static_cast<double>(
	    controller->getLastProgressiveProviderTimeNanoseconds()) /
	    1000000.0);
    sample.insert(QStringLiteral("last_lod_submission_ms"),
	static_cast<double>(
	    controller->getLastLodSubmissionTimeNanoseconds()) /
	    1000000.0);
    sample.insert(QStringLiteral("last_presentation_sync_ms"),
	static_cast<double>(
	    controller->getLastPresentationSyncTimeNanoseconds()) /
	    1000000.0);
    sample.insert(QStringLiteral("render_requested"),
	controller->isRenderRequested() ? true : false);
    sample.insert(QStringLiteral("render_reason"),
	QString::fromLocal8Bit(controller->getRenderReason().getString()));
    sample.insert(QStringLiteral("progressive_pending"),
	controller->hasProgressiveWorkPending() ? true : false);
    sample.insert(QStringLiteral("lod_results_pending"),
	controller->hasPendingLodResults() ? true : false);
    sample.insert(QStringLiteral("lod_submissions_pending"),
	controller->hasPendingLodSubmissions() ? true : false);
    sample.insert(QStringLiteral("lod_refinement_frame_pending"),
	controller->hasPendingLodRefinementFrame() ? true : false);
    sample.insert(QStringLiteral("lod_view_revision"),
	static_cast<qint64>(controller->getLodViewRevision()));
    sample.insert(QStringLiteral("lod_policy_revision"),
	static_cast<qint64>(controller->getLodPolicyRevision()));
    sample.insert(QStringLiteral("lod_interactive"),
	controller->isLodInteractionActive() ? true : false);
    sample.insert(QStringLiteral("lod_gesture_active"),
	controller->isLodGestureActive() ? true : false);
    sample.insert(QStringLiteral("lod_scale_changing_interaction"),
	controller->isLodScaleChangingInteraction() ? true : false);
    sample.insert(QStringLiteral("lod_target_pixel_error"),
	static_cast<double>(controller->getLodTargetPixelError()));
    sample.insert(QStringLiteral("lod_interactive_progressive_ceiling"),
	controller->getLodInteractiveProgressiveCeiling());
    sample.insert(QStringLiteral("lod_visited_meshes"),
	static_cast<int>(controller->getLastLodVisitedMeshCount()));
    sample.insert(QStringLiteral("lod_submitted_tasks"),
	static_cast<int>(controller->getLastLodSubmittedTaskCount()));
    sample.insert(QStringLiteral("lod_updated_cuts"),
	static_cast<int>(controller->getLastLodUpdatedCutCount()));
    sample.insert(QStringLiteral("lod_skipped_meshes"),
	static_cast<int>(controller->getLastLodSkippedMeshCount()));
    sample.insert(QStringLiteral("lod_results"),
	static_cast<qint64>(controller->getLastLodResultCount()));
    sample.insert(QStringLiteral("lod_matched_results"),
	static_cast<int>(controller->getLastLodMatchedResultCount()));
    sample.insert(QStringLiteral("lod_applied_results"),
	static_cast<int>(controller->getLastLodAppliedResultCount()));
    sample.insert(QStringLiteral("lod_rejected_results"),
	static_cast<int>(controller->getLastLodRejectedResultCount()));
    sample.insert(QStringLiteral("lod_unmatched_results"),
	static_cast<int>(controller->getLastLodUnmatchedResultCount()));
    sample.insert(QStringLiteral("active_lod_mesh_payloads"),
	static_cast<qint64>(controller->getActiveLodMeshPayloadCount()));
    sample.insert(QStringLiteral("active_lod_proxy_payloads"),
	static_cast<qint64>(controller->getActiveLodProxyPayloadCount(
	    BOBOL_LOD_PROXY_NONE)));
    sample.insert(QStringLiteral("active_lod_aabb_payloads"),
	static_cast<qint64>(controller->getActiveLodProxyPayloadCount(
	    BOBOL_LOD_PROXY_AABB)));
    sample.insert(QStringLiteral("active_lod_obb_payloads"),
	static_cast<qint64>(controller->getActiveLodProxyPayloadCount(
	    BOBOL_LOD_PROXY_OBB)));
    sample.insert(QStringLiteral("active_lod_sphere_payloads"),
	static_cast<qint64>(controller->getActiveLodProxyPayloadCount(
	    BOBOL_LOD_PROXY_SPHERE)));
    sample.insert(QStringLiteral("active_lod_cad_payloads"),
	static_cast<qint64>(controller->getActiveLodCadPayloadCount()));
    BObolLodConvergenceStatus convergence;
    controller->getLodConvergenceStatus(convergence);
    sample.insert(QStringLiteral("lod_convergence_phase"),
	convergence.phase);
    sample.insert(QStringLiteral("lod_coordinator_phase"),
	convergence.coordinatorPhase);
    sample.insert(QStringLiteral("lod_coordinator_event"),
	convergence.coordinatorEvent);
    sample.insert(QStringLiteral("lod_coordinator_transition_serial"),
	static_cast<qint64>(convergence.coordinatorTransitionSerial));
    sample.insert(QStringLiteral("lod_coordinator_progress_sequence"),
	static_cast<qint64>(convergence.coordinatorProgressSequence));
    sample.insert(QStringLiteral("lod_coordinator_dispatch_serial"),
	static_cast<qint64>(convergence.coordinatorDispatchSerial));
    sample.insert(QStringLiteral("lod_coordinator_stagnant_dispatches"),
	static_cast<qint64>(convergence.coordinatorStagnantDispatchCount));
    sample.insert(QStringLiteral("lod_coordinator_invariant_violations"),
	static_cast<qint64>(convergence.coordinatorInvariantViolationCount));
    sample.insert(QStringLiteral("lod_coordinator_invariant_mask"),
	static_cast<qint64>(convergence.coordinatorInvariantMask));
    sample.insert(QStringLiteral("lod_coordinator_invariant_history_mask"),
	static_cast<qint64>(convergence.coordinatorInvariantHistoryMask));
    sample.insert(QStringLiteral("lod_convergence_fraction"),
	static_cast<double>(convergence.fraction));
    sample.insert(QStringLiteral("lod_convergence_view_ready"),
	convergence.viewReady ? true : false);
    sample.insert(QStringLiteral("lod_convergence_background_pending"),
	convergence.backgroundPending ? true : false);
    sample.insert(QStringLiteral("lod_convergence_performance_limited"),
	convergence.performanceLimited ? true : false);
    sample.insert(QStringLiteral("lod_convergence_expected_leaves"),
	static_cast<qint64>(convergence.expectedLeafCount));
    sample.insert(QStringLiteral("lod_convergence_available_leaves"),
	static_cast<qint64>(convergence.availableLeafCount));
    sample.insert(QStringLiteral("lod_convergence_visible_targets"),
	static_cast<qint64>(convergence.visibleTargetCount));
    sample.insert(QStringLiteral("lod_convergence_active_payloads"),
	static_cast<qint64>(convergence.activePayloadCount));
    sample.insert(QStringLiteral("lod_convergence_satisfied_payloads"),
	static_cast<qint64>(convergence.satisfiedPayloadCount));
    sample.insert(QStringLiteral("lod_convergence_presented_subpixel_occurrences"),
	static_cast<qint64>(std::min<size_t>(
	    convergence.presentedSubpixelOccurrenceCount,
	    static_cast<size_t>(std::numeric_limits<qint64>::max()))));
    sample.insert(QStringLiteral("lod_convergence_presented_structural_boxes"),
	static_cast<qint64>(std::min<size_t>(
	    convergence.presentedStructuralBoxCount,
	    static_cast<size_t>(std::numeric_limits<qint64>::max()))));
    sample.insert(QStringLiteral("lod_convergence_terminal_occurrence_failures"),
	static_cast<qint64>(std::min<size_t>(
	    convergence.terminalOccurrenceFailureCount,
	    static_cast<size_t>(std::numeric_limits<qint64>::max()))));
    sample.insert(QStringLiteral("lod_convergence_resident_mesh_bytes"),
	static_cast<qint64>(convergence.residentMeshBytes));
    sample.insert(QStringLiteral("lod_convergence_stable_resident_mesh_bytes"),
	static_cast<qint64>(convergence.stableResidentMeshBytes));
    sample.insert(QStringLiteral("lod_convergence_reserved_mesh_growth_bytes"),
	static_cast<qint64>(convergence.reservedResidentMeshGrowthBytes));
    sample.insert(QStringLiteral("lod_convergence_resident_mesh_limit_bytes"),
	static_cast<qint64>(convergence.residentMeshLimitBytes));
    sample.insert(QStringLiteral("lod_convergence_memory_limited_payloads"),
	static_cast<qint64>(convergence.memoryLimitedPayloadCount));
    sample.insert(QStringLiteral("lod_convergence_memory_limited"),
	convergence.memoryLimited ? true : false);
    sample.insert(QStringLiteral("lod_convergence_refinement_frame_pending"),
	convergence.refinementFramePending ? true : false);
    sample.insert(QStringLiteral("lod_active_generation"),
	static_cast<qint64>(convergence.activeGeneration));
    sample.insert(QStringLiteral("lod_submission_source_index"),
	static_cast<qint64>(convergence.submissionSourceIndex));
    sample.insert(QStringLiteral("lod_submission_entry_offset"),
	static_cast<qint64>(convergence.submissionEntryOffset));
    sample.insert(QStringLiteral("lod_convergence_budget_calibration_pending"),
	convergence.budgetCalibrationPending ? true : false);
    sample.insert(QStringLiteral(
	"lod_convergence_stable_handoff_pending"),
	convergence.stablePresentationHandoffPending ? true : false);
    sample.insert(QStringLiteral(
	"lod_convergence_point_proxy_calibration_pending"),
	convergence.pointProxyCalibrationPending ? true : false);
    sample.insert(QStringLiteral(
	"lod_convergence_resident_growth_reallocation_pending"),
	convergence.residentGrowthReallocationPending ? true : false);
    sample.insert(QStringLiteral("lod_convergence_publication_frame_pending"),
	convergence.publicationFramePending ? true : false);
    sample.insert(QStringLiteral("lod_convergence_active_working_set_bytes"),
	static_cast<qint64>(convergence.activeWorkingSetBytes));
    sample.insert(QStringLiteral("lod_convergence_peak_working_set_bytes"),
	static_cast<qint64>(convergence.peakWorkingSetBytes));
    sample.insert(QStringLiteral("lod_convergence_compactions"),
	static_cast<qint64>(convergence.residentCompactionCount));
    sample.insert(QStringLiteral("lod_gpu_tracked_buffer_bytes"),
	static_cast<qint64>(convergence.gpuTrackedBufferBytes));
    sample.insert(QStringLiteral("lod_gpu_ordinary_part_buffer_bytes"),
	static_cast<qint64>(convergence.gpuOrdinaryPartBufferBytes));
    sample.insert(QStringLiteral("lod_gpu_progressive_cut_buffer_bytes"),
	static_cast<qint64>(convergence.gpuProgressiveCutBufferBytes));
    sample.insert(QStringLiteral("lod_gpu_progressive_active_cut_bytes"),
	static_cast<qint64>(convergence.gpuProgressiveActiveCutBytes));
    sample.insert(QStringLiteral("lod_gpu_batch_buffer_bytes"),
	static_cast<qint64>(convergence.gpuBatchBufferBytes));
    sample.insert(QStringLiteral("lod_gpu_atlas_allocated_bytes"),
	static_cast<qint64>(convergence.gpuTriangleAtlasAllocatedBytes));
    sample.insert(QStringLiteral("lod_gpu_atlas_live_bytes"),
	static_cast<qint64>(convergence.gpuTriangleAtlasLiveBytes));
    sample.insert(QStringLiteral("lod_gpu_atlas_configured_capacity_bytes"),
	static_cast<qint64>(
	    convergence.gpuTriangleAtlasConfiguredCapacityBytes));
    sample.insert(QStringLiteral("lod_gpu_atlas_parts"),
	static_cast<qint64>(convergence.gpuTriangleAtlasPartCount));
    sample.insert(QStringLiteral("lod_gpu_atlas_pages"),
	static_cast<qint64>(convergence.gpuTriangleAtlasPageCount));
    sample.insert(QStringLiteral("lod_gpu_ordinary_full_upload_bytes"),
	static_cast<qint64>(convergence.gpuOrdinaryPartFullUploadBytes));
    sample.insert(QStringLiteral("lod_gpu_ordinary_suffix_upload_bytes"),
	static_cast<qint64>(convergence.gpuOrdinaryPartSuffixUploadBytes));
    sample.insert(QStringLiteral("lod_gpu_ordinary_copy_bytes"),
	static_cast<qint64>(convergence.gpuOrdinaryPartGpuCopyBytes));
    sample.insert(QStringLiteral("lod_gpu_ordinary_lineage_reuses"),
	static_cast<qint64>(convergence.gpuOrdinaryPartLineageReuseCount));
    sample.insert(QStringLiteral("lod_gpu_atlas_full_upload_bytes"),
	static_cast<qint64>(convergence.gpuTriangleAtlasFullUploadBytes));
    sample.insert(QStringLiteral("lod_gpu_atlas_suffix_upload_bytes"),
	static_cast<qint64>(convergence.gpuTriangleAtlasSuffixUploadBytes));
    sample.insert(QStringLiteral("lod_gpu_atlas_lineage_reuses"),
	static_cast<qint64>(convergence.gpuTriangleAtlasLineageReuseCount));
    sample.insert(QStringLiteral("lod_gpu_pressure_proxies"),
	static_cast<qint64>(convergence.gpuPressureProxyCount));
    sample.insert(QStringLiteral("lod_gpu_progressive_evictions"),
	static_cast<qint64>(convergence.gpuProgressiveEvictionCount));
    sample.insert(QStringLiteral("lod_gpu_atlas_reclamations"),
	static_cast<qint64>(convergence.gpuTriangleAtlasReclamationCount));
    sample.insert(QStringLiteral("lod_gpu_resource_sample_serial"),
	static_cast<qint64>(convergence.gpuResourceSampleSerial));
    sample.insert(QStringLiteral("lod_gpu_memory_pressure"),
	convergence.gpuMemoryPressure ? true : false);
    BObolViewLodState *viewLodState = controller->getViewLodState();
    qint64 compactEntries = 0;
    qint64 compactSelectedEntries = 0;
    qint64 cadSelectedInstances = 0;
    qint64 compactLodEntries = 0;
    qint64 compactFallbackEntries = 0;
    qint64 compactLodEntriesWithPayload = 0;
    qint64 compactFallbackEntriesWithPayload = 0;
    qint64 cadPayloadsWithoutEntry = 0;
    qint64 supersededFallbackPresentations = 0;
    qint64 activeStructuralFallbackPresentations = 0;
    qint64 activeProgressiveCadPayloads = 0;
    qint64 activeProgressiveCadCullPayloads = 0;
    double activeProgressiveCadFaces = 0.0;
    double activeProgressiveCadCullFaces = 0.0;
    double activeProgressiveCadPoints = 0.0;
    qint64 projectedDemandCadPayloads = 0;
    qint64 prominentCadPayloads = 0;
    qint64 cadQualityFloorViolations = 0;
    qint64 prominentCadQualityFloorViolations = 0;
    double maximumCadProjectedErrorPixels = 0.0;
    double maximumCadNormalizedError = 0.0;
    double maximumCadVisualFootprint = 0.0;
    double cadVisualImportanceDebt = 0.0;
    struct CadVisualOutlier {
	std::string sourcePath;
	std::string sourceName;
	std::string occurrenceKey;
	int activeCut;
	int residentCut;
	int requestedCut;
	double diameter;
	double footprint;
	double targetError;
	double projectedError;
	double normalizedError;
	uint64_t faces;
	uint64_t points;
    };
    std::vector<CadVisualOutlier> cadVisualOutliers;
    qint64 activeCadSubpixelProxyPoints = 0;
    qint64 visibleStructuralFallbackBoxes = 0;
    qint64 cadOccurrenceTerminalFailures = 0;
    qint64 presentedCadFaces = 0;
    qint64 presentedCadLines = 0;
    qint64 presentedCadPositions = 0;
    qint64 presentedCadNormals = 0;
    qint64 presentedCadOccurrences = 0;
    bool presentedCadWorkExact = true;
    qint64 measuredCadGpuFaces = 0;
    uint64_t measuredCadGpuNanoseconds = 0;
    uint64_t measuredCadGpuSampleSerial = 0;
    double measuredCadGpuPointProxyPixelThreshold = 1.0;
    qint64 cadFramePlanBuildCountMax = 0;
    qint64 cadFramePlanInstanceRecordCountMax = 0;
    bool cadPreparedReplayAll = true;
    qint64 cadPresentationCount = 0;
    double cadPointProxyPixelThresholdMax = 1.0;
    int cadRenderTierMin = INT_MAX;
    int cadRenderTierMax = -1;
    int cadIndirectStatusMin = INT_MAX;
    int cadIndirectStatusMax = -1;
    std::unordered_set<const SoCADAssembly *> sampledCadPresentations;
    int activeProgressiveCadCutMin = -1;
    int activeProgressiveCadCutMax = -1;
    int requestedProgressiveCadCutMin = -1;
    int requestedProgressiveCadCutMax = -1;
    uint64_t activeProgressiveCadOccurrenceHash = 0;
    QJsonArray supersededFallbackPaths;
    QJsonArray databaseSourceBounds;
    QJsonArray compactPlanningBounds;
    QJsonArray compactMissingPayloadBounds;
    /*
     * Payload aggregation is intentionally checkpoint-only.  At distinct
     * asset scale findCadPayloadsUnordered plus the occurrence digest walks
     * every resident payload; doing that after each synthetic mouse move
     * added roughly 35 ms of GUI-thread work and made the observer manufacture
     * the low presentation rate it was meant to measure.  Controller timing
     * and aggregate scene-budget counters above remain available on every
     * event.  The terminal event also requests deep diagnostics.
     */
    const bool collectPresentationDiagnostics =
	collectDeepLodDiagnostics ||
	event.action == QLatin1String("checkpoint");
    if (collectPresentationDiagnostics) {
	if (viewLodState)
	    cadOccurrenceTerminalFailures = static_cast<qint64>(
		std::min<size_t>(
		    viewLodState->cadOccurrenceTerminalFailureCount(),
		    static_cast<size_t>(
			std::numeric_limits<qint64>::max())));
	std::vector<SoBRLDatabaseSource *> renderSources;
	qged_test_collect_database_source_roots(
	    controller->getRenderSceneRoot(), renderSources);
	if (renderSources.empty())
	    qged_test_collect_database_source_roots(controller->getSceneRoot(),
		renderSources);
	for (SoBRLDatabaseSource *source : renderSources) {
	if (!source || !source->hasCompactInstanceIndex())
	    continue;
	compactSelectedEntries +=
	    source->getCompactSelectedInstanceCount();
	if (viewLodState) {
	    SoCADAssembly *presentation =
		viewLodState->findCadPresentation(source);
	    if (presentation &&
		sampledCadPresentations.insert(presentation).second) {
		cadPresentationCount++;
		cadPreparedReplayAll =
		    cadPreparedReplayAll &&
		    presentation->lastRenderUsedPreparedReplay();
		cadFramePlanBuildCountMax = std::max(
		    cadFramePlanBuildCountMax,
		    static_cast<qint64>(std::min<uint64_t>(
			presentation->framePlanBuildCount(),
			static_cast<uint64_t>(
			    std::numeric_limits<qint64>::max()))));
		cadFramePlanInstanceRecordCountMax = std::max(
		    cadFramePlanInstanceRecordCountMax,
		    static_cast<qint64>(std::min<size_t>(
			presentation->framePlanInstanceRecordCount(),
			static_cast<size_t>(
			    std::numeric_limits<qint64>::max()))));
		cadSelectedInstances += static_cast<qint64>(
		    presentation->selectedInstanceCount());
		activeCadSubpixelProxyPoints += static_cast<qint64>(
		    presentation->lastSubpixelProxyCount());
		cadPointProxyPixelThresholdMax = std::max(
		    cadPointProxyPixelThresholdMax,
		    static_cast<double>(
			presentation->pointProxyPixelThreshold.getValue()));
		visibleStructuralFallbackBoxes += static_cast<qint64>(
		    presentation->lastUncollapsedStructuralProxyCount());
		const Obol::CadRenderedWork renderedWork =
		    presentation->lastRenderedWork();
		presentedCadWorkExact =
		    presentedCadWorkExact && renderedWork.exact;
		const auto addRenderedWork = [](qint64 current, uint64_t value) {
		    return value > static_cast<uint64_t>(
			std::numeric_limits<qint64>::max() - current) ?
			std::numeric_limits<qint64>::max() :
			current + static_cast<qint64>(value);
		};
		presentedCadFaces = addRenderedWork(
		    presentedCadFaces, renderedWork.triangleCount);
		presentedCadLines = addRenderedWork(
		    presentedCadLines, renderedWork.lineCount);
		presentedCadPositions = addRenderedWork(
		    presentedCadPositions, renderedWork.positionCount);
		presentedCadNormals = addRenderedWork(
		    presentedCadNormals, renderedWork.normalCount);
		presentedCadOccurrences = addRenderedWork(
		    presentedCadOccurrences, renderedWork.occurrenceCount);
		const uint64_t gpuSampleSerial =
		    presentation->gpuTimerSampleSerial();
		const uint64_t gpuNanoseconds =
		    presentation->lastGpuRenderNanoseconds();
		const uint64_t gpuTriangles =
		    presentation->lastGpuRenderedTriangleCount();
		const double gpuPointThreshold =
		    presentation->lastGpuPointProxyPixelThreshold();
		if (gpuSampleSerial && gpuNanoseconds) {
		    measuredCadGpuSampleSerial =
			std::max(measuredCadGpuSampleSerial, gpuSampleSerial);
		    measuredCadGpuNanoseconds =
			gpuNanoseconds >
			    UINT64_MAX - measuredCadGpuNanoseconds ?
			UINT64_MAX :
			measuredCadGpuNanoseconds + gpuNanoseconds;
		    measuredCadGpuFaces =
			gpuTriangles > static_cast<uint64_t>(
			    std::numeric_limits<qint64>::max() -
				measuredCadGpuFaces) ?
			std::numeric_limits<qint64>::max() :
			measuredCadGpuFaces +
			    static_cast<qint64>(gpuTriangles);
		    measuredCadGpuPointProxyPixelThreshold = std::max(
			measuredCadGpuPointProxyPixelThreshold,
			gpuPointThreshold);
		}
		const int tier = presentation->lastRenderTier();
		if (tier >= 0) {
		    cadRenderTierMin = std::min(cadRenderTierMin, tier);
		    cadRenderTierMax = std::max(cadRenderTierMax, tier);
		}
		const int indirectStatus =
		    presentation->lastIndirectStatus();
		if (indirectStatus >= 0) {
		    cadIndirectStatusMin = std::min(
			cadIndirectStatusMin, indirectStatus);
		    cadIndirectStatusMax = std::max(
			cadIndirectStatusMax, indirectStatus);
		}
	    }
	}
	if (collectDeepLodDiagnostics && databaseSourceBounds.size() < 16) {
	    QJsonObject boundsSample;
	    boundsSample.insert(QStringLiteral("path"),
		QString::fromLocal8Bit(source->path.getValue().getString()));
	    SbBox3f bounds;
	    const bool valid = source->getSourceBounds(bounds) &&
		!bounds.isEmpty();
	    boundsSample.insert(QStringLiteral("valid"), valid);
	    boundsSample.insert(QStringLiteral("exact"),
		source->hasExactSourceBounds() ? true : false);
	    boundsSample.insert(QStringLiteral("realization_status"),
		source->realizationStatus.getValue());
	    boundsSample.insert(QStringLiteral("stale"),
		source->stale.getValue() ? true : false);
	    boundsSample.insert(QStringLiteral("stale_reason"),
		static_cast<int>(source->staleReason.getValue()));
	    boundsSample.insert(QStringLiteral("diagnostic"),
		QString::fromLocal8Bit(
		    source->realizationDiagnostic.getValue().getString()));
	    if (valid) {
		const SbVec3f minimum = bounds.getMin();
		const SbVec3f maximum = bounds.getMax();
		boundsSample.insert(QStringLiteral("min_x"), minimum[0]);
		boundsSample.insert(QStringLiteral("min_y"), minimum[1]);
		boundsSample.insert(QStringLiteral("min_z"), minimum[2]);
		boundsSample.insert(QStringLiteral("max_x"), maximum[0]);
		boundsSample.insert(QStringLiteral("max_y"), maximum[1]);
		boundsSample.insert(QStringLiteral("max_z"), maximum[2]);
	    }
	    databaseSourceBounds.append(boundsSample);
	}
	std::vector<const BObolViewLodState::CadPayload *> payloads;
	if (viewLodState)
	    viewLodState->findCadPayloadsUnordered(source, payloads);
	for (const BObolViewLodState::CadPayload *payload : payloads) {
	    if (!payload || !payload->progressiveMesh)
		continue;
	    activeProgressiveCadPayloads++;
	    /*
	     * Keep this order-independent so the unordered payload store does
	     * not perturb reports.  Hash each key independently and combine the
	     * digests commutatively.  The former copy/sort/hash implementation
	     * allocated and sorted 50,000 strings after every motion event,
	     * making the test observer a material GUI-thread perf sample.
	     */
	    uint64_t occurrenceDigest = UINT64_C(1469598103934665603);
	    for (const unsigned char *c =
		    reinterpret_cast<const unsigned char *>(
			payload->sourceInstanceKey.getString());
		 *c; ++c) {
		occurrenceDigest ^= *c;
		occurrenceDigest *= UINT64_C(1099511628211);
	    }
	    occurrenceDigest ^= occurrenceDigest >> 33;
	    occurrenceDigest *= UINT64_C(0xff51afd7ed558ccd);
	    occurrenceDigest ^= occurrenceDigest >> 33;
	    occurrenceDigest *= UINT64_C(0xc4ceb9fe1a85ec53);
	    occurrenceDigest ^= occurrenceDigest >> 33;
	    activeProgressiveCadOccurrenceHash ^= occurrenceDigest;
	    if (payload->shadedCullBackfaces) {
		activeProgressiveCadCullPayloads++;
		activeProgressiveCadCullFaces +=
		    static_cast<double>(payload->counts.faceCount);
	    }
	    activeProgressiveCadFaces +=
		static_cast<double>(payload->counts.faceCount);
	    activeProgressiveCadPoints +=
		static_cast<double>(payload->counts.pointCount);
	    if (payload->activeCut >= 0 &&
		std::isfinite(payload->projectedPixelDiameter) &&
		payload->projectedPixelDiameter > 0.0f) {
		projectedDemandCadPayloads++;
		const double diameter = payload->projectedPixelDiameter;
		const double area = std::max(0.0,
		    static_cast<double>(payload->projectedPixelArea));
		const double perimeter = std::max(0.0,
		    static_cast<double>(payload->projectedPixelPerimeter));
		const double footprint = std::max(std::sqrt(area),
		    std::max(perimeter * 0.25, diameter * 0.25));
		/* A cut is an opaque producer-authored schedule position, not an
		 * isotropic quantization-bit count.  The canonical PoP schedule may
		 * refine only one object axis at a cut and may retain an unchanged
		 * topology population solely to improve coordinate precision.  Ask
		 * the retained hierarchy for its actual object-error projection;
		 * diameter / 2^cut was an obsolete fixed-level approximation that
		 * falsely reported converged anisotropic meshes as several pixels
		 * coarse and consequently corrupted the HUD and matrix oracle. */
		const double projectedError =
		    payload->progressiveMesh->projectedErrorAtCut(
			payload->activeCut, diameter);
		const double target = std::max(1.0,
		    static_cast<double>(payload->targetPixelError));
		const double normalizedError = projectedError / target;
		const bool prominent = footprint >= 16.0;
		const bool floorViolation = normalizedError > 3.0;
		prominentCadPayloads += prominent ? 1 : 0;
		cadQualityFloorViolations += floorViolation ? 1 : 0;
		prominentCadQualityFloorViolations +=
		    prominent && floorViolation ? 1 : 0;
		maximumCadProjectedErrorPixels = std::max(
		    maximumCadProjectedErrorPixels, projectedError);
		maximumCadNormalizedError = std::max(
		    maximumCadNormalizedError, normalizedError);
		maximumCadVisualFootprint = std::max(
		    maximumCadVisualFootprint, footprint);
		cadVisualImportanceDebt += footprint *
		    std::max(0.0, normalizedError - 1.0);
		/* Keep a fixed-size, deterministic witness set for perceptual
		 * regressions.  Scene aggregates revealed that Hubble wire drawing
		 * was unfair, but not which occurrence retained a minimum cut.  A
		 * bounded top-N avoids turning 50k/150k checkpoint telemetry into an
		 * unbounded JSON or string-allocation workload. */
		const char *occurrenceKey =
		    payload->sourceInstanceKey.getString();
		const auto ranksBefore = [](double errorA, double footprintA,
			const char *keyA, const CadVisualOutlier &b) {
		    if (errorA > b.normalizedError)
			return true;
		    if (errorA < b.normalizedError)
			return false;
		    if (footprintA > b.footprint)
			return true;
		    if (footprintA < b.footprint)
			return false;
		    return std::string(keyA ? keyA : "") < b.occurrenceKey;
		};
		const bool entersOutlierSet = cadVisualOutliers.size() < 16 ||
		    ranksBefore(normalizedError, footprint, occurrenceKey,
			cadVisualOutliers.back());
		if (prominent && normalizedError > 1.0 && entersOutlierSet) {
		    CadVisualOutlier outlier;
		    outlier.sourcePath = payload->sourcePath.getString();
		    outlier.sourceName = payload->sourceName.getString();
		    outlier.occurrenceKey = occurrenceKey ? occurrenceKey : "";
		    outlier.activeCut = payload->activeCut;
		    outlier.residentCut = payload->residentCut;
		    outlier.requestedCut = payload->requestedCut;
		    outlier.diameter = diameter;
		    outlier.footprint = footprint;
		    outlier.targetError = target;
		    outlier.projectedError = projectedError;
		    outlier.normalizedError = normalizedError;
		    outlier.faces = payload->counts.faceCount;
		    outlier.points = payload->counts.pointCount;
		    cadVisualOutliers.push_back(std::move(outlier));
		    std::sort(cadVisualOutliers.begin(), cadVisualOutliers.end(),
			[](const CadVisualOutlier &a,
			    const CadVisualOutlier &b) {
			    if (a.normalizedError > b.normalizedError)
				return true;
			    if (a.normalizedError < b.normalizedError)
				return false;
			    if (a.footprint > b.footprint)
				return true;
			    if (a.footprint < b.footprint)
				return false;
			    return a.occurrenceKey < b.occurrenceKey;
			});
		    if (cadVisualOutliers.size() > 16)
			cadVisualOutliers.resize(16);
		}
	    }
	    if (activeProgressiveCadCutMin < 0 ||
		payload->activeCut < activeProgressiveCadCutMin)
		activeProgressiveCadCutMin = payload->activeCut;
	    activeProgressiveCadCutMax = std::max(
		activeProgressiveCadCutMax, payload->activeCut);
	    if (requestedProgressiveCadCutMin < 0 ||
		payload->requestedCut < requestedProgressiveCadCutMin)
		requestedProgressiveCadCutMin = payload->requestedCut;
	    requestedProgressiveCadCutMax = std::max(
		requestedProgressiveCadCutMax, payload->requestedCut);
	}
	if (collectDeepLodDiagnostics) {
	    std::unordered_set<std::string> payloadKeys;
	    payloadKeys.reserve(payloads.size());
	    for (const BObolViewLodState::CadPayload *payload : payloads)
		if (payload && payload->sourceInstanceKey.getLength() > 0)
		    payloadKeys.insert(payload->sourceInstanceKey.getString());
	    std::unordered_set<std::string> entryKeys;
	    const int entryCount = source->getCompactInstanceCount();
	    entryKeys.reserve(static_cast<size_t>(entryCount));
	    for (int entryIndex = 0; entryIndex < entryCount; entryIndex++) {
		if (compactPlanningBounds.size() < 16) {
		    BObolCompactLodPlanningSummary planning;
		    if (source->getCompactLodPlanningSummary(entryIndex,
			    planning) && planning.valid) {
			QJsonObject planningSample;
			planningSample.insert(QStringLiteral("entry"),
			    entryIndex);
			planningSample.insert(QStringLiteral("lod_backed"),
			    planning.lodBacked ? true : false);
			planningSample.insert(
			    QStringLiteral("source_request_valid"),
			    planning.sourceMeshRequestValid ? true : false);
			const SbVec3f minimum = planning.localBounds.getMin();
			const SbVec3f maximum = planning.localBounds.getMax();
			planningSample.insert(QStringLiteral("min_x"),
			    minimum[0]);
			planningSample.insert(QStringLiteral("min_y"),
			    minimum[1]);
			planningSample.insert(QStringLiteral("min_z"),
			    minimum[2]);
			planningSample.insert(QStringLiteral("max_x"),
			    maximum[0]);
			planningSample.insert(QStringLiteral("max_y"),
			    maximum[1]);
			planningSample.insert(QStringLiteral("max_z"),
			    maximum[2]);
			SbVec3f origin;
			planning.localToSource.multVecMatrix(
			    SbVec3f(0.0f, 0.0f, 0.0f), origin);
			planningSample.insert(QStringLiteral("origin_x"),
			    origin[0]);
			planningSample.insert(QStringLiteral("origin_y"),
			    origin[1]);
			planningSample.insert(QStringLiteral("origin_z"),
			    origin[2]);
			compactPlanningBounds.append(planningSample);
		    }
		}
		BObolCompactInstanceHandle handle;
		BObolCompactInstanceSummary summary;
		if (!source->getCompactInstanceHandle(entryIndex, handle) ||
		    !source->getCompactInstanceSummary(handle, summary))
		    continue;
		compactEntries++;
		const std::string key = summary.sourceInstanceKey.getString();
		if (!key.empty())
		    entryKeys.insert(key);
		const bool hasPayload = !key.empty() &&
		    payloadKeys.find(key) != payloadKeys.end();
		if (!hasPayload && compactMissingPayloadBounds.size() < 16) {
		    BObolCompactLodPlanningSummary planning;
		    if (source->getCompactLodPlanningSummary(entryIndex,
			    planning) && planning.valid) {
			SbBox3f worldBounds;
			worldBounds.makeEmpty();
			const SbVec3f minimum = planning.localBounds.getMin();
			const SbVec3f maximum = planning.localBounds.getMax();
			for (int corner = 0; corner < 8; ++corner) {
			    const SbVec3f local(
				(corner & 1) ? maximum[0] : minimum[0],
				(corner & 2) ? maximum[1] : minimum[1],
				(corner & 4) ? maximum[2] : minimum[2]);
			    SbVec3f world;
			    planning.localToSource.multVecMatrix(local, world);
			    worldBounds.extendBy(world);
			}
			QJsonObject missingSample;
			missingSample.insert(QStringLiteral("entry"),
			    entryIndex);
			missingSample.insert(QStringLiteral("key"),
			    QString::fromStdString(key));
			missingSample.insert(QStringLiteral("eligible"),
			    planning.lodBacked &&
				planning.sourceMeshRequestValid);
			if (!worldBounds.isEmpty()) {
			    const SbVec3f worldMinimum =
				worldBounds.getMin();
			    const SbVec3f worldMaximum =
				worldBounds.getMax();
			    missingSample.insert(QStringLiteral("min_x"),
				worldMinimum[0]);
			    missingSample.insert(QStringLiteral("min_y"),
				worldMinimum[1]);
			    missingSample.insert(QStringLiteral("min_z"),
				worldMinimum[2]);
			    missingSample.insert(QStringLiteral("max_x"),
				worldMaximum[0]);
			    missingSample.insert(QStringLiteral("max_y"),
				worldMaximum[1]);
			    missingSample.insert(QStringLiteral("max_z"),
				worldMaximum[2]);
			}
			compactMissingPayloadBounds.append(missingSample);
		    }
		}
		if (summary.lodBacked && summary.sourceMeshRequestValid) {
		    compactLodEntries++;
		    if (hasPayload)
			compactLodEntriesWithPayload++;
		}
		if (BU_STR_EQUAL(summary.geometryKind.getString(), "aabb") ||
		    BU_STR_EQUAL(summary.geometryKind.getString(), "obb")) {
		    compactFallbackEntries++;
		    if (hasPayload)
			compactFallbackEntriesWithPayload++;
		}
	    }
	    for (const std::string &key : payloadKeys)
		if (entryKeys.find(key) == entryKeys.end())
		    cadPayloadsWithoutEntry++;
	    std::vector<SbString> supersededPaths;
	    supersededFallbackPresentations +=
		source->getCompactViewLodSupersededFallbackCount(viewLodState,
		&supersededPaths);
	    activeStructuralFallbackPresentations +=
		source->getCompactViewLodActiveFallbackCount(viewLodState);
	    for (const SbString &path : supersededPaths)
		if (supersededFallbackPaths.size() < 16)
		    supersededFallbackPaths.append(
			QString::fromLocal8Bit(path.getString()));
	}
    }
    }
    sample.insert(QStringLiteral("deep_lod_diagnostics"),
	collectDeepLodDiagnostics);
    sample.insert(QStringLiteral("compact_entries"), compactEntries);
    sample.insert(QStringLiteral("compact_selected_entries"),
	compactSelectedEntries);
    sample.insert(QStringLiteral("cad_selected_instances"),
	cadSelectedInstances);
    sample.insert(QStringLiteral("compact_lod_entries"), compactLodEntries);
    sample.insert(QStringLiteral("compact_fallback_entries"),
	compactFallbackEntries);
    sample.insert(QStringLiteral("compact_lod_entries_with_payload"),
	compactLodEntriesWithPayload);
    sample.insert(QStringLiteral("compact_fallback_entries_with_payload"),
	compactFallbackEntriesWithPayload);
    sample.insert(QStringLiteral("cad_payloads_without_entry"),
	cadPayloadsWithoutEntry);
    sample.insert(QStringLiteral("superseded_fallback_presentations"),
	supersededFallbackPresentations);
    sample.insert(QStringLiteral("active_structural_fallback_presentations"),
	activeStructuralFallbackPresentations);
    sample.insert(QStringLiteral("superseded_fallback_paths"),
	supersededFallbackPaths);
    if (collectDeepLodDiagnostics) {
	sample.insert(QStringLiteral("database_source_bounds"),
	    databaseSourceBounds);
	sample.insert(QStringLiteral("compact_planning_bounds"),
	    compactPlanningBounds);
	sample.insert(QStringLiteral("compact_missing_payload_bounds"),
	    compactMissingPayloadBounds);
    }
    sample.insert(QStringLiteral("active_progressive_cad_payloads"),
	activeProgressiveCadPayloads);
    sample.insert(QStringLiteral("lod_projected_demand_cad_payloads"),
	projectedDemandCadPayloads);
    sample.insert(QStringLiteral("lod_prominent_cad_payloads"),
	prominentCadPayloads);
    sample.insert(QStringLiteral("lod_cad_quality_floor_violations"),
	cadQualityFloorViolations);
    sample.insert(
	QStringLiteral("lod_prominent_cad_quality_floor_violations"),
	prominentCadQualityFloorViolations);
    sample.insert(QStringLiteral("lod_max_cad_projected_error_pixels"),
	maximumCadProjectedErrorPixels);
    sample.insert(QStringLiteral("lod_max_cad_normalized_error"),
	maximumCadNormalizedError);
    sample.insert(QStringLiteral("lod_max_cad_visual_footprint_pixels"),
	maximumCadVisualFootprint);
    sample.insert(QStringLiteral("lod_cad_visual_importance_debt"),
	cadVisualImportanceDebt);
    QJsonArray cadVisualOutlierSamples;
    for (const CadVisualOutlier &outlier : cadVisualOutliers) {
	QJsonObject item;
	item.insert(QStringLiteral("source_path"),
	    QString::fromStdString(outlier.sourcePath));
	item.insert(QStringLiteral("source_name"),
	    QString::fromStdString(outlier.sourceName));
	item.insert(QStringLiteral("occurrence_key"),
	    QString::fromStdString(outlier.occurrenceKey));
	item.insert(QStringLiteral("active_cut"), outlier.activeCut);
	item.insert(QStringLiteral("resident_cut"), outlier.residentCut);
	item.insert(QStringLiteral("requested_cut"), outlier.requestedCut);
	item.insert(QStringLiteral("projected_diameter_pixels"),
	    outlier.diameter);
	item.insert(QStringLiteral("visual_footprint_pixels"),
	    outlier.footprint);
	item.insert(QStringLiteral("target_error_pixels"),
	    outlier.targetError);
	item.insert(QStringLiteral("projected_error_pixels"),
	    outlier.projectedError);
	item.insert(QStringLiteral("normalized_error"),
	    outlier.normalizedError);
	item.insert(QStringLiteral("faces"),
	    static_cast<qint64>(std::min<uint64_t>(outlier.faces,
		static_cast<uint64_t>(std::numeric_limits<qint64>::max()))));
	item.insert(QStringLiteral("points"),
	    static_cast<qint64>(std::min<uint64_t>(outlier.points,
		static_cast<uint64_t>(std::numeric_limits<qint64>::max()))));
	cadVisualOutlierSamples.append(item);
    }
    sample.insert(QStringLiteral("lod_cad_visual_importance_outliers"),
	cadVisualOutlierSamples);
    sample.insert(QStringLiteral("active_progressive_cad_occurrence_hash"),
	QString::number(
	    static_cast<qulonglong>(activeProgressiveCadOccurrenceHash), 16));
    sample.insert(QStringLiteral("active_progressive_cad_cull_payloads"),
	activeProgressiveCadCullPayloads);
    sample.insert(QStringLiteral("active_progressive_cad_faces"),
	activeProgressiveCadFaces);
    sample.insert(QStringLiteral("active_progressive_cad_cull_faces"),
	activeProgressiveCadCullFaces);
    sample.insert(QStringLiteral("active_progressive_cad_points"),
	activeProgressiveCadPoints);
    sample.insert(QStringLiteral("active_cad_subpixel_proxy_points"),
	activeCadSubpixelProxyPoints);
    sample.insert(QStringLiteral("cad_point_proxy_pixel_threshold_max"),
	cadPointProxyPixelThresholdMax);
    sample.insert(QStringLiteral("visible_structural_fallback_boxes"),
	visibleStructuralFallbackBoxes);
    sample.insert(QStringLiteral("cad_occurrence_terminal_failures"),
	cadOccurrenceTerminalFailures);
    sample.insert(QStringLiteral("presented_cad_faces"),
	presentedCadFaces);
    sample.insert(QStringLiteral("presented_cad_lines"),
	presentedCadLines);
    sample.insert(QStringLiteral("presented_cad_positions"),
	presentedCadPositions);
    sample.insert(QStringLiteral("presented_cad_normals"),
	presentedCadNormals);
    sample.insert(QStringLiteral("presented_cad_occurrences"),
	presentedCadOccurrences);
    sample.insert(QStringLiteral("presented_cad_work_exact"),
	cadPresentationCount > 0 && presentedCadWorkExact);
    sample.insert(QStringLiteral("measured_cad_gpu_faces"),
	measuredCadGpuFaces);
    sample.insert(QStringLiteral("measured_cad_gpu_ms"),
	static_cast<double>(measuredCadGpuNanoseconds) / 1000000.0);
    sample.insert(QStringLiteral("measured_cad_gpu_sample_serial"),
	QString::number(measuredCadGpuSampleSerial));
    sample.insert(
	QStringLiteral("measured_cad_gpu_point_proxy_pixel_threshold"),
	measuredCadGpuPointProxyPixelThreshold);
    sample.insert(QStringLiteral("cad_frame_plan_build_count_max"),
	cadFramePlanBuildCountMax);
    sample.insert(QStringLiteral("cad_frame_plan_instance_record_count_max"),
	cadFramePlanInstanceRecordCountMax);
    sample.insert(QStringLiteral("cad_prepared_replay_all"),
	cadPresentationCount > 0 && cadPreparedReplayAll);
    sample.insert(QStringLiteral("cad_render_tier_min"),
	cadRenderTierMin == INT_MAX ? -1 : cadRenderTierMin);
    sample.insert(QStringLiteral("cad_render_tier_max"), cadRenderTierMax);
    sample.insert(QStringLiteral("cad_indirect_status_min"),
	cadIndirectStatusMin == INT_MAX ? -1 : cadIndirectStatusMin);
    sample.insert(QStringLiteral("cad_indirect_status_max"),
	cadIndirectStatusMax);
    sample.insert(QStringLiteral("active_progressive_cad_cut_min"),
	activeProgressiveCadCutMin);
    sample.insert(QStringLiteral("active_progressive_cad_cut_max"),
	activeProgressiveCadCutMax);
    sample.insert(QStringLiteral("requested_progressive_cad_cut_min"),
	requestedProgressiveCadCutMin);
    sample.insert(QStringLiteral("requested_progressive_cad_cut_max"),
	requestedProgressiveCadCutMax);
    BObolLodService *service = controller->getLodService();
    if (service) {
	sample.insert(QStringLiteral("lod_service_in_flight_tasks"),
	    static_cast<qint64>(service->inFlightCount()));
	sample.insert(QStringLiteral("lod_service_result_reservations"),
	    static_cast<qint64>(
		service->resultReservationCountForDiagnostics()));
	sample.insert(QStringLiteral("lod_service_available_result_capacity"),
	    static_cast<qint64>(service->availableResultTaskCapacity()));
	sample.insert(QStringLiteral("lod_service_pending_tasks"),
	    static_cast<qint64>(service->pendingTaskCountForDiagnostics()));
	sample.insert(QStringLiteral("lod_service_active_requests"),
	    static_cast<qint64>(service->activeRequestCountForDiagnostics()));
	sample.insert(QStringLiteral("lod_service_queued_results"),
	    static_cast<qint64>(service->queuedResultCountForDiagnostics()));
	sample.insert(QStringLiteral("lod_service_queued_cache_writes"),
	    static_cast<qint64>(
		service->queuedCacheWriteCountForDiagnostics()));
	sample.insert(QStringLiteral("lod_service_completed_tasks"),
	    static_cast<qint64>(service->completedTaskCountForDiagnostics()));
	sample.insert(QStringLiteral("lod_service_resident_assets"),
	    static_cast<qint64>(
		service->residentMeshAssetCountForDiagnostics()));
	sample.insert(QStringLiteral("lod_service_resident_bytes"),
	    static_cast<qint64>(
		service->residentMeshBytesForDiagnostics()));
	sample.insert(QStringLiteral("lod_service_stable_resident_bytes"),
	    static_cast<qint64>(
		service->stableResidentMeshBytesForDiagnostics()));
	sample.insert(QStringLiteral("lod_service_reserved_growth_bytes"),
	    static_cast<qint64>(
		service->reservedResidentMeshGrowthBytesForDiagnostics()));
	sample.insert(QStringLiteral("lod_service_resident_admission_revision"),
	    static_cast<qint64>(
		service->residentMeshAdmissionRevision()));
	sample.insert(QStringLiteral("lod_service_resident_limit_bytes"),
	    static_cast<qint64>(service->getResidentMeshLimit()));
	sample.insert(QStringLiteral("lod_service_working_set_limit_bytes"),
	    static_cast<qint64>(service->getWorkingSetLimit()));
	sample.insert(QStringLiteral("lod_service_active_working_set_bytes"),
	    static_cast<qint64>(
		service->activeWorkingSetBytesForDiagnostics()));
	sample.insert(QStringLiteral("lod_service_executing_tasks"),
	    static_cast<qint64>(
		service->executingTaskCountForDiagnostics()));
	sample.insert(QStringLiteral("lod_service_peak_working_set_bytes"),
	    static_cast<qint64>(
		service->peakWorkingSetBytesForDiagnostics()));
	sample.insert(QStringLiteral("lod_service_peak_executing_tasks"),
	    static_cast<qint64>(
		service->peakExecutingTaskCountForDiagnostics()));
	sample.insert(QStringLiteral("lod_service_cache_loads"),
	    static_cast<qint64>(
		service->residentMeshCacheLoadCountForDiagnostics()));
	sample.insert(QStringLiteral("lod_service_resident_hits"),
	    static_cast<qint64>(
		service->residentMeshHitCountForDiagnostics()));
	sample.insert(QStringLiteral("lod_service_compactions"),
	    static_cast<qint64>(
		service->residentMeshCompactionCountForDiagnostics()));
	sample.insert(QStringLiteral("lod_service_evictions"),
	    static_cast<qint64>(
		service->residentMeshEvictionCountForDiagnostics()));
    }
    sample.insert(QStringLiteral("lod_global_working_set_limit_bytes"),
	static_cast<qint64>(bobol_lod_working_set_global_limit()));
    sample.insert(QStringLiteral("lod_global_active_working_set_bytes"),
	static_cast<qint64>(bobol_lod_working_set_global_active_bytes()));
    sample.insert(QStringLiteral("lod_global_active_working_set_tasks"),
	static_cast<qint64>(bobol_lod_working_set_global_active_tasks()));
    sample.insert(QStringLiteral("lod_global_peak_working_set_bytes"),
	static_cast<qint64>(bobol_lod_working_set_global_peak_bytes()));
    sample.insert(QStringLiteral("lod_global_peak_working_set_tasks"),
	static_cast<qint64>(bobol_lod_working_set_global_peak_tasks()));
    sample.insert(QStringLiteral("visited_sources"),
	static_cast<int>(controller->getLastVisitedSourceCount()));
    sample.insert(QStringLiteral("realized_sources"),
	static_cast<int>(controller->getLastRealizedSourceCount()));
    sample.insert(QStringLiteral("failed_sources"),
	static_cast<int>(controller->getLastFailedSourceCount()));
    sample.insert(QStringLiteral("last_mesh_budget_final_bytes"),
	static_cast<qint64>(
	    controller->getLastMeshBudgetFinalResidentBytes()));
    sample.insert(QStringLiteral("diagnostics"),
	QString::fromLocal8Bit(controller->getLastDiagnostics().getString()));
    sample.insert(QStringLiteral("lod_diagnostics"),
	QString::fromLocal8Bit(controller->getLastLodDiagnostics().getString()));
    return sample;
}


/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
