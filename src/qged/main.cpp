/*                        M A I N . C X X
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
/** @file main.cxx
 *
 * Command line parsing and main application launching for qbrlcad
 *
 */

#include "common.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <unordered_set>
#include <vector>

#include <QApplication>
#include <QCoreApplication>
#include <QElapsedTimer>
#include <QImage>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QOpenGLWidget>
#include <QPixmap>
#include <QRegularExpression>
#include <QSaveFile>
#include <QSurfaceFormat>
#include <QTimer>

#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoGroup.h>

#include "BObol/BLodService.h"
#include "BObol/BDatabaseSource.h"
#include "BObol/BViewController.h"
#include "BObol/BViewLod.h"
#include "bu/app.h"
#include "bu/log.h"
#include "bu/opt.h"
#include "bu/vls.h"
#include "brlcad_version.h"
#include "ged/scene.h"
#include "ged/selection_state.h"

#include "QgEdApp.h"
#include "fbserv.h"
#include "qtcad/QgTestEvent.h"

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

static QJsonObject
qged_test_controller_sample(QgEdApp &app, int eventIndex,
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
    if (event.action == QLatin1String("checkpoint"))
	sample.insert(QStringLiteral("checkpoint"),
	    event.arguments.value(QStringLiteral("name")));
    sample.insert(QStringLiteral("elapsed_ms"), elapsedMilliseconds);
    sample.insert(QStringLiteral("event_duration_us"), eventMicroseconds);
    /* Full draw-frontier and selection listings cross the GED/Obol hierarchy
     * and are structural diagnostics, not frame telemetry.  Sampling them
     * after every 8 ms motion interval made the validation harness dominate
     * its own perf captures. */
    const bool collectStructuralDiagnostics =
	collectDeepLodDiagnostics ||
	event.action == QLatin1String("qged_command") ||
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

	if (gedp && ged_draw_scene_available(gedp)) {
	    sample.insert(QStringLiteral("draw_scene_revision"),
		QString::number(ged_draw_scene_revision(gedp)));
	    if (collectStructuralDiagnostics) {
		sample.insert(QStringLiteral("draw_shape_count"),
		    ged_draw_shape_count(gedp));
		struct bu_vls drawListing = BU_VLS_INIT_ZERO;
		const size_t drawCount = ged_draw_list_paths(gedp,
		    ged_view_active_ctx(gedp), -1, 0, &drawListing);
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
    }

    QgView *view = app.w ? app.w->CurrentDisplay() : nullptr;
    BObolViewController *controller =
	view ? view->obolViewController() : nullptr;
    if (!controller) {
	sample.insert(QStringLiteral("controller_available"), false);
	return sample;
    }

    sample.insert(QStringLiteral("controller_available"), true);
    if (SoCamera *camera = controller->getCamera()) {
	const SbVec3f position = camera->position.getValue();
	sample.insert(QStringLiteral("camera_position_x"),
	    static_cast<double>(position[0]));
	sample.insert(QStringLiteral("camera_position_y"),
	    static_cast<double>(position[1]));
	sample.insert(QStringLiteral("camera_position_z"),
	    static_cast<double>(position[2]));
	sample.insert(QStringLiteral("camera_focal_distance"),
	    static_cast<double>(camera->focalDistance.getValue()));
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
    sample.insert(QStringLiteral("active_lod_scene_faces"),
	static_cast<qint64>(controller->getActiveLodFaceCount()));
    sample.insert(QStringLiteral("lod_scene_face_budget"),
	static_cast<qint64>(controller->getCurrentLodFaceBudget()));
    sample.insert(QStringLiteral("lod_calibrated_faces_per_second"),
	controller->getCalibratedLodFacesPerSecond());
    sample.insert(QStringLiteral(
	"lod_interactive_calibrated_faces_per_second"),
	controller->getInteractiveCalibratedLodFacesPerSecond());
    sample.insert(QStringLiteral("lod_stable_calibrated_faces_per_second"),
	controller->getStableCalibratedLodFacesPerSecond());
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
    BObolViewLodState *viewLodState = controller->getViewLodState();
    qint64 compactEntries = 0;
    qint64 compactLodEntries = 0;
    qint64 compactFallbackEntries = 0;
    qint64 compactLodEntriesWithPayload = 0;
    qint64 compactFallbackEntriesWithPayload = 0;
    qint64 cadPayloadsWithoutEntry = 0;
    qint64 supersededFallbackPresentations = 0;
    qint64 activeProgressiveCadPayloads = 0;
    qint64 activeProgressiveCadCullPayloads = 0;
    double activeProgressiveCadFaces = 0.0;
    double activeProgressiveCadCullFaces = 0.0;
    double activeProgressiveCadPoints = 0.0;
    int activeProgressiveCadLevelMin = -1;
    int activeProgressiveCadLevelMax = -1;
    int requestedProgressiveCadLevelMin = -1;
    int requestedProgressiveCadLevelMax = -1;
    std::vector<std::string> activeProgressiveCadOccurrenceKeys;
    QJsonArray supersededFallbackPaths;
    std::vector<SoBRLDatabaseSource *> renderSources;
    qged_test_collect_database_source_roots(
	controller->getRenderSceneRoot(), renderSources);
    if (renderSources.empty())
	qged_test_collect_database_source_roots(controller->getSceneRoot(),
	    renderSources);
    for (SoBRLDatabaseSource *source : renderSources) {
	if (!source || !source->hasCompactInstanceIndex())
	    continue;
	std::vector<const BObolViewLodState::CadPayload *> payloads;
	if (viewLodState)
	    viewLodState->findCadPayloadsUnordered(source, payloads);
	for (const BObolViewLodState::CadPayload *payload : payloads) {
	    if (!payload || !payload->progressiveMesh)
		continue;
	    activeProgressiveCadPayloads++;
	    activeProgressiveCadOccurrenceKeys.push_back(
		payload->sourceInstanceKey.getString());
	    if (payload->shadedCullBackfaces) {
		activeProgressiveCadCullPayloads++;
		activeProgressiveCadCullFaces +=
		    static_cast<double>(payload->counts.faceCount);
	    }
	    activeProgressiveCadFaces +=
		static_cast<double>(payload->counts.faceCount);
	    activeProgressiveCadPoints +=
		static_cast<double>(payload->counts.pointCount);
	    if (activeProgressiveCadLevelMin < 0 ||
		payload->activeLevel < activeProgressiveCadLevelMin)
		activeProgressiveCadLevelMin = payload->activeLevel;
	    activeProgressiveCadLevelMax = std::max(
		activeProgressiveCadLevelMax, payload->activeLevel);
	    if (requestedProgressiveCadLevelMin < 0 ||
		payload->requestedLevel < requestedProgressiveCadLevelMin)
		requestedProgressiveCadLevelMin = payload->requestedLevel;
	    requestedProgressiveCadLevelMax = std::max(
		requestedProgressiveCadLevelMax, payload->requestedLevel);
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
	    for (const SbString &path : supersededPaths)
		if (supersededFallbackPaths.size() < 16)
		    supersededFallbackPaths.append(
			QString::fromLocal8Bit(path.getString()));
	}
    }
    std::sort(activeProgressiveCadOccurrenceKeys.begin(),
	activeProgressiveCadOccurrenceKeys.end());
    uint64_t activeProgressiveCadOccurrenceHash =
	UINT64_C(1469598103934665603);
    for (const std::string &key : activeProgressiveCadOccurrenceKeys) {
	for (unsigned char c : key) {
	    activeProgressiveCadOccurrenceHash ^= c;
	    activeProgressiveCadOccurrenceHash *= UINT64_C(1099511628211);
	}
	activeProgressiveCadOccurrenceHash ^= 0xff;
	activeProgressiveCadOccurrenceHash *= UINT64_C(1099511628211);
    }
    sample.insert(QStringLiteral("deep_lod_diagnostics"),
	collectDeepLodDiagnostics);
    sample.insert(QStringLiteral("compact_entries"), compactEntries);
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
    sample.insert(QStringLiteral("superseded_fallback_paths"),
	supersededFallbackPaths);
    sample.insert(QStringLiteral("active_progressive_cad_payloads"),
	activeProgressiveCadPayloads);
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
    sample.insert(QStringLiteral("active_progressive_cad_level_min"),
	activeProgressiveCadLevelMin);
    sample.insert(QStringLiteral("active_progressive_cad_level_max"),
	activeProgressiveCadLevelMax);
    sample.insert(QStringLiteral("requested_progressive_cad_level_min"),
	requestedProgressiveCadLevelMin);
    sample.insert(QStringLiteral("requested_progressive_cad_level_max"),
	requestedProgressiveCadLevelMax);
    BObolLodService *service = controller->getLodService();
    if (service) {
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

static bool
qged_write_test_report(const QString &fileName, const QJsonObject &report,
    QString *error)
{
    if (fileName.isEmpty())
	return true;
    QSaveFile file(fileName);
    if (!file.open(QIODevice::WriteOnly)) {
	if (error)
	    *error = file.errorString();
	return false;
    }
    const QByteArray data = QJsonDocument(report).toJson(
	QJsonDocument::Indented);
    if (file.write(data) != data.size() || !file.commit()) {
	if (error)
	    *error = file.errorString();
	return false;
    }
    return true;
}

static bool
qged_test_wait_progressive_idle(QgEdApp &app, int timeoutMilliseconds,
    int quietMilliseconds, QString *error)
{
    QgView *view = app.w ? app.w->CurrentDisplay() : nullptr;
    BObolViewController *controller =
	view ? view->obolViewController() : nullptr;
    if (!controller) {
	if (error)
	    *error = QStringLiteral(
		"wait_progressive_idle requires an Obol view controller");
	return false;
    }

    timeoutMilliseconds = std::max(0, timeoutMilliseconds);
    quietMilliseconds = std::max(0, quietMilliseconds);
    QElapsedTimer elapsed;
    QElapsedTimer quiet;
    elapsed.start();

    while (elapsed.elapsed() <= timeoutMilliseconds) {
	BObolLodService *service = controller->getLodService();
	const bool serviceIdle =
	    !service ||
	    (service->pendingTaskCountForDiagnostics() == 0 &&
	     service->delayedTaskCountForDiagnostics() == 0 &&
	     service->inFlightCount() == 0 &&
	     service->activeRequestCountForDiagnostics() == 0 &&
	     service->queuedResultCountForDiagnostics() == 0 &&
	     service->queuedCacheWriteCountForDiagnostics() == 0);
	const bool idle =
	    !controller->hasProgressiveWorkPending() &&
	    !controller->hasPendingLodResults() &&
	    !controller->hasPendingLodSubmissions() &&
	    !controller->hasPendingLodRefinementFrame() &&
	    serviceIdle;

	if (idle) {
	    if (!quiet.isValid())
		quiet.start();
	    if (quiet.elapsed() >= quietMilliseconds)
		return true;
	} else {
	    quiet.invalidate();
	}

	/* Let queued paints, worker wakeups, and the 16 ms progressive frame
	 * pump run, while retaining a hard deadline for a defective pipeline. */
	QEventLoop loop;
	const int remaining = timeoutMilliseconds -
	    static_cast<int>(elapsed.elapsed());
	QTimer::singleShot(std::max(1, std::min(16, remaining)),
	    &loop, &QEventLoop::quit);
	loop.exec(QEventLoop::AllEvents);
    }

    if (error) {
	BObolLodService *service = controller->getLodService();
	*error = QStringLiteral(
	    "progressive pipeline did not become idle within %1 ms "
	    "(progressive=%2 results=%3 submissions=%4 refinement_frame=%5 "
	    "pending=%6 delayed=%7 in_flight=%8 active=%9 queued=%10 "
	    "cache_writes=%11)")
	    .arg(timeoutMilliseconds)
	    .arg(controller->hasProgressiveWorkPending() ? 1 : 0)
	    .arg(controller->hasPendingLodResults() ? 1 : 0)
	    .arg(controller->hasPendingLodSubmissions() ? 1 : 0)
	    .arg(controller->hasPendingLodRefinementFrame() ? 1 : 0)
	    .arg(service ?
		static_cast<qulonglong>(
		    service->pendingTaskCountForDiagnostics()) : 0)
	    .arg(service ?
		static_cast<qulonglong>(
		    service->delayedTaskCountForDiagnostics()) : 0)
	    .arg(service ?
		static_cast<qulonglong>(service->inFlightCount()) : 0)
	    .arg(service ?
		static_cast<qulonglong>(
		    service->activeRequestCountForDiagnostics()) : 0)
	    .arg(service ?
		static_cast<qulonglong>(
		    service->queuedResultCountForDiagnostics()) : 0)
	    .arg(service ?
		static_cast<qulonglong>(
		    service->queuedCacheWriteCountForDiagnostics()) : 0);
    }
    return false;
}

static void
qged_usage(const char *cmd, struct bu_opt_desc *d) {
    struct bu_vls str = BU_VLS_INIT_ZERO;
    char *option_help = bu_opt_describe(d, NULL);
    bu_vls_sprintf(&str, "Usage: %s [options] [file.g]\n", cmd);
    if (option_help) {
	bu_vls_printf(&str, "Options:\n%s\n", option_help);
    }
    bu_free(option_help, "help str");
    bu_log("%s", bu_vls_cstr(&str));
    bu_vls_free(&str);
}

#ifdef HAVE_WINDOWS_H
int APIENTRY
WinMain(HINSTANCE hInstance,
	HINSTANCE hPrevInstance,
	LPSTR lpszCmdLine,
	int nCmdShow)
{
    int argc = __argc;
    char **argv = __argv;
#else
int
main(int argc, char **argv)
{
#endif

    int console_mode = 0;
    int swrast_mode = 0;
    int quad_mode = 0;
    int print_help = 0;
    struct bu_vls msg = BU_VLS_INIT_ZERO;
    const char *startup_commands = NULL;
    const char *test_script = NULL;
    const char *test_report = NULL;
    const char *exec_name = argv[0];

    // All BRL-CAD programs need to set this in order for relative path lookups
    // to work reliably.
    bu_setprogname(argv[0]);

    /* Done with command name argv[0] */
    argc-=(argc>0); argv+=(argc>0);

    /* Handle top level application options */
    struct bu_opt_desc d[9];
    BU_OPT(d[0],  "h", "help",   "", NULL, &print_help,    "Print help and exit");
    BU_OPT(d[1],  "?", "",       "", NULL, &print_help,    "");
    BU_OPT(d[2],  "c", "no-gui", "", NULL, &console_mode,  "Run without GUI");
    BU_OPT(d[3],  "s", "swrast", "", NULL, &swrast_mode,   "Use offscreen rendering for 3D view");
    BU_OPT(d[4],  "4", "quad",   "", NULL, &quad_mode,     "Launch using quad view");
    BU_OPT(d[5],  "e", "exec",   "commands", &bu_opt_str, &startup_commands,
	    "Run semicolon-separated GED commands after GUI initialization");
    BU_OPT(d[6],  "", "test-script", "events.json", &bu_opt_str, &test_script,
	    "Replay a qtcad GUI event stream and exit with test status");
    BU_OPT(d[7],  "", "test-report", "report.json", &bu_opt_str, &test_report,
	    "Write timing and LoD samples while replaying a GUI test");
    BU_OPT_NULL(d[8]);

    // High level options are only defined prior to the file argument (if there
    // is one).  See if we need to limit our processing
    int acmax = 0;
    for (int i = 0; i < argc; i++) {
	if (argv[i][0] == '-') {
	    acmax++;
	    if ((BU_STR_EQUAL(argv[i], "-e") ||
		    BU_STR_EQUAL(argv[i], "--exec") ||
		    BU_STR_EQUAL(argv[i], "--test-script") ||
		    BU_STR_EQUAL(argv[i], "--test-report")) && i + 1 < argc) {
		acmax++;
		i++;
	    }
	} else {
	    break;
	}
    }

    // Process high level args
    int opt_ac = bu_opt_parse(&msg, acmax, (const char **)argv, d);
    if (opt_ac < 0) {
	bu_log("%s\n", bu_vls_cstr(&msg));
	bu_vls_free(&msg);
	return BRLCAD_ERROR;
    }

    // Shift everything not processed by bu_opt_parse down in argv to the last
    // option left behind by bu_opt_parse (or the beginning of the array if
    // opt_ac == 0
    int opt_rem = opt_ac;
    for (int i = acmax; i < argc; i++) {
	argv[opt_ac + (i - acmax)] = argv[i];
	opt_rem++;
    }

    // Set argc to the full count of whatever is left
    argc = opt_rem;

    if (print_help) {
	qged_usage(exec_name, d);
	bu_vls_free(&msg);
	return BRLCAD_OK;
    }

    if (argc > 1 && !console_mode) {
	bu_log("For qged GUI mode need either zero or one .g files specified\n");
	return BRLCAD_ERROR;
    }

    if (console_mode) {
	bu_log("Unimplemented\n");
	return BRLCAD_ERROR;
    }

    if (startup_commands)
	bu_log("qged: queued startup commands: %s\n", startup_commands);

    // Qt6 requires QSurfaceFormat::setDefaultFormat() to be called BEFORE
    // QApplication is constructed.  Without specifying QSurfaceFormat::OpenGL
    // here, Qt's RHI layer may fall back to OpenGL ES (QRhiGles2) even on a
    // desktop system, causing context-creation failures at runtime.
#ifdef BRLCAD_OPENGL
    if (!swrast_mode) {
	QSurfaceFormat fmt;
	fmt.setRenderableType(QSurfaceFormat::OpenGL);
	fmt.setDepthBufferSize(24);
	fmt.setStencilBufferSize(8);
	/* QOpenGLWidget renders into a Qt-owned FBO which the window-system
	 * compositor presents later.  V-syncing that intermediate context adds a
	 * second, unrelated presentation gate.  On GLX this reduced an otherwise
	 * 45 FPS retained view to about 7 FPS and made an 8 ms input interval wait
	 * nearly a second for the compositor; it does not prevent tearing in the
	 * final composed window.  Let Qt/the compositor own pacing and keep the
	 * off-screen widget context unpaced on every platform.  PartialUpdate in
	 * QgGL preserves the last complete FBO, so this no longer exposes the
	 * invalid-buffer flash that originally motivated the Unix interval of 1.
	 * The override remains available for driver experiments. */
	int swapInterval = 0;
	const char *swapOverride = std::getenv("QGED_SWAP_INTERVAL");
	if (swapOverride && swapOverride[0]) {
	    char *end = NULL;
	    const long value = std::strtol(swapOverride, &end, 10);
	    if (end && end[0] == '\0' && value >= -1 && value <= 1)
		swapInterval = static_cast<int>(value);
	    else
		bu_log("qged: ignoring invalid QGED_SWAP_INTERVAL=%s "
		       "(expected -1, 0, or 1)\n", swapOverride);
	}
	fmt.setSwapInterval(swapInterval);
	QSurfaceFormat::setDefaultFormat(fmt);
    }
#endif

    const QString startup = startup_commands ?
	QString::fromLocal8Bit(startup_commands) : QString();

    // QApplication requires argv[0] to be the executable name.  The qged
    // option parser above has removed that entry and left only the optional
    // database argument, so give Qt a valid argument vector and pass the
    // database path separately.
    const char *db_file = argc ? argv[0] : NULL;
    int qt_argc = 1;
    char *qt_argv[] = {const_cast<char *>(exec_name), NULL};
    QgEdApp app(qt_argc, qt_argv, db_file, swrast_mode, quad_mode);
    if (!startup.isEmpty()) {
	const QStringList commands = startup.split(';', Qt::SkipEmptyParts);
	for (const QString &command : commands) {
	    const QString trimmed = command.trimmed();
	    bu_log("qged: executing startup command: %s\n",
		    trimmed.toLocal8Bit().constData());
	    app.run_qcmd(trimmed);
	}
    }
    if (test_script) {
	const QString script = QString::fromLocal8Bit(test_script);
	const QString reportFile = test_report ?
	    QString::fromLocal8Bit(test_report) : QString();
	const bool software = swrast_mode != 0;
	QTimer::singleShot(0, &app, [&app, script, reportFile, software]() {
	    QVector<QgTestEvent> events;
	    QString error;
	    QgEventPlayer player(app.w);
	    QElapsedTimer elapsed;
	    QJsonArray samples;
	    QJsonObject report;
	    report.insert(QStringLiteral("schema"),
		QStringLiteral("brlcad.qged.gui-report"));
	    report.insert(QStringLiteral("version"), 1);
	    report.insert(QStringLiteral("event_script"), script);
	    report.insert(QStringLiteral("backend"),
		software ? QStringLiteral("osmesa") :
		QStringLiteral("system_gl"));
	    report.insert(QStringLiteral("cache_directory"),
		QString::fromLocal8Bit(std::getenv("BU_DIR_CACHE") ?
		    std::getenv("BU_DIR_CACHE") : ""));
	    report.insert(QStringLiteral("swap_interval"),
		QString::fromLocal8Bit(std::getenv("QGED_SWAP_INTERVAL") ?
		    std::getenv("QGED_SWAP_INTERVAL") : "default"));
	    player.setCheckpointHandler([](QWidget *widget,
		    const QString &name, QString *checkpointError) {
		if (!widget || name.isEmpty()) {
		    if (checkpointError)
			*checkpointError =
			    QStringLiteral("checkpoint needs a widget and output path");
		    return false;
		}
		bool saved = false;
		if (QOpenGLWidget *glWidget =
			qobject_cast<QOpenGLWidget *>(widget)) {
		    const QImage frame = glWidget->grabFramebuffer();
		    saved = !frame.isNull() && frame.save(name);
		} else {
		    const QPixmap frame = widget->grab();
		    saved = !frame.isNull() && frame.save(name);
		}
		if (!saved) {
		    if (checkpointError)
			*checkpointError =
			    QStringLiteral("unable to save checkpoint image: %1")
			    .arg(name);
		    return false;
		}
		return true;
	    });
	    elapsed.start();
	    bool success = QgEventRecorder::load(script, events, &error);
	    for (int i = 0; success && i < events.size(); ++i) {
		QElapsedTimer eventElapsed;
		eventElapsed.start();
		if (events[i].action == QLatin1String("qged_command")) {
		    const QString command = events[i].arguments.value(
			QStringLiteral("command")).toString();
		    if (command.isEmpty()) {
			error = QStringLiteral(
			    "qged_command requires a command argument");
			success = false;
		    } else {
			app.run_qcmd(command);
		    }
		} else if (events[i].action ==
		    QLatin1String("wait_progressive_idle")) {
		    success = qged_test_wait_progressive_idle(app,
			events[i].arguments.value(
			    QStringLiteral("timeout_ms")).toInt(30000),
			events[i].arguments.value(
			    QStringLiteral("quiet_ms")).toInt(100),
			&error);
		} else {
		    success = player.play(events[i], &error);
		}
		/* Do not drain the event queue between scripted inputs.  Progressive
		 * rendering deliberately queues another frame while work remains,
		 * so even a nominally bounded processEvents call can consume a
		 * backlog of complete frames and turn one input into a one-second
		 * pause.  Commands and synthetic inputs dispatch synchronously;
		 * explicit wait events provide the real event-loop intervals. */
		const qint64 eventMicroseconds =
		    eventElapsed.nsecsElapsed() / 1000;
		const bool deepLodDiagnostics =
		    i + 1 == events.size() ||
		    std::getenv("QGED_TEST_DEEP_LOD_REPORT") != NULL;
		QElapsedTimer sampleElapsed;
		sampleElapsed.start();
		QJsonObject sample = qged_test_controller_sample(app, i,
		    events[i], elapsed.elapsed(), eventMicroseconds,
		    deepLodDiagnostics);
		sample.insert(QStringLiteral("sample_duration_us"),
		    sampleElapsed.nsecsElapsed() / 1000);
		samples.append(sample);
	    }
	    report.insert(QStringLiteral("success"), success);
	    report.insert(QStringLiteral("elapsed_ms"), elapsed.elapsed());
	    report.insert(QStringLiteral("samples"), samples);
	    if (!success)
		report.insert(QStringLiteral("error"), error);
	    QString reportError;
	    if (!qged_write_test_report(reportFile, report, &reportError)) {
		bu_log("qged: unable to write GUI test report: %s\n",
		    reportError.toLocal8Bit().constData());
		success = false;
	    }
	    if (!success) {
		bu_log("qged: GUI test replay failed: %s\n",
		    error.toLocal8Bit().constData());
		app.exit(BRLCAD_ERROR);
		return;
	    }
	    app.exit(BRLCAD_OK);
	});
    }
    bu_vls_free(&msg);

    // Setup complete - time to enter the interactive event loop
    return app.exec();
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
