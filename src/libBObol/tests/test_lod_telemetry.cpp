/*                    T E S T _ L O D _ T E L E M E T R Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol.h"
#include "BObol/BInit.h"
#include "BObol/BLodService.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "raytrace.h"
#include "wdb.h"

#include <Obol/cad/CadGeometry.h>
#include <Obol/cad/CadGeometryValidation.h>

#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoSeparator.h>

#include "../../libbu/json.hpp"

#include <cstdio>
#include <chrono>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

static const char *telemetry_file_environment =
    "BOBOL_LOD_TELEMETRY_FILE";
static const char *telemetry_sanitize_environment =
    "BOBOL_LOD_TELEMETRY_SANITIZE";
static const char *private_combination = "classified-root-private.c";
static const char *private_solid = "classified-private.bot";
static const char *private_second_solid = "classified-second-private.bot";
static const char *private_instance = "classified-private-instance";
static const char *private_result_diagnostic =
    "classified-private-result-diagnostic";

static int
make_database(const char *path)
{
    struct rt_wdb *wdbp = wdb_fopen_v(path, 5);
    if (!wdbp)
	return 0;

    fastf_t vertices[12] = {
	0.0, 0.0, 0.0,
	1.0, 0.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 0.0, 1.0
    };
    int faces[12] = {
	0, 1, 2,
	0, 3, 1,
	1, 3, 2,
	2, 3, 0
    };
    int result = mk_bot(wdbp, private_solid, RT_BOT_SOLID,
	RT_BOT_UNORIENTED, 0, 4, 4, vertices, faces, NULL, NULL);
    struct wmember members;
    BU_LIST_INIT(&members.l);
    if (result == 0 &&
	(!mk_addmember(private_solid, &members.l, NULL, WMOP_UNION) ||
	 mk_lcomb(wdbp, private_combination, &members, 0, NULL, NULL,
	     NULL, 0) != 0))
	result = -1;
    wdb_close(wdbp);
    return result == 0;
}

static std::shared_ptr<const Obol::PartGeometry>
make_geometry(void)
{
    Obol::PartGeometryBuilder geometry;
    Obol::TriMesh mesh;
    mesh.positions.push_back(SbVec3f(0.0f, 0.0f, 0.0f));
    mesh.positions.push_back(SbVec3f(1.0f, 0.0f, 0.0f));
    mesh.positions.push_back(SbVec3f(0.0f, 1.0f, 0.0f));
    mesh.indices.push_back(0);
    mesh.indices.push_back(1);
    mesh.indices.push_back(2);
    mesh.bounds = SbBox3f(SbVec3f(0.0f, 0.0f, 0.0f),
	SbVec3f(1.0f, 1.0f, 0.0f));
    geometry.shaded = std::move(mesh);
    const Obol::CadGeometryAdmission admitted =
	Obol::cadAdmitPartGeometry(std::move(geometry));
    return admitted ? admitted.geometry.shared() :
	std::shared_ptr<const Obol::PartGeometry>();
}

static BObolLodResult
cancelled_result(const BObolLodRequest &request, void *)
{
    BObolLodResult result;
    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.providerStatus = BOBOL_LOD_PROVIDER_CANCELLED;
    result.terminal = TRUE;
    result.diagnostic = private_result_diagnostic;
    return result;
}

static std::string
read_file(const char *path)
{
    std::ifstream input(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(input),
	std::istreambuf_iterator<char>());
}

static std::string
path_basename(const char *path)
{
    const std::string value = path ? path : "";
    const size_t separator = value.find_last_of("/\\");
    return separator == std::string::npos ? value :
	value.substr(separator + 1);
}

static int
validate_sanitized_log(const std::string &log,
    const std::string &databaseBasename)
{
    if (log.find(private_combination) != std::string::npos ||
	log.find(private_solid) != std::string::npos ||
	log.find(private_second_solid) != std::string::npos ||
	log.find(private_instance) != std::string::npos ||
	log.find(private_result_diagnostic) != std::string::npos ||
	(!databaseBasename.empty() &&
	 log.find(databaseBasename) != std::string::npos)) {
	std::fprintf(stderr,
	    "FAIL: sanitized telemetry disclosed database identity\n");
	return 1;
    }
    if (log.find("\"render_reason\":") != std::string::npos) {
	std::fprintf(stderr,
	    "FAIL: sanitized telemetry disclosed an arbitrary diagnostic string\n");
	return 1;
    }

    bool schema = false;
    bool transition = false;
    bool source = false;
    bool object = false;
    bool hiddenObject = false;
    bool initialInventoryReset = false;
    bool incrementalInventoryDelta = false;
    bool secondObject = false;
    bool publication = false;
    bool nonFiniteCameraValid = false;
    bool signalsConsistent = true;
    std::istringstream lines(log);
    std::string line;
    try {
	while (std::getline(lines, line)) {
	    if (line.empty())
		continue;
	    const nlohmann::json record = nlohmann::json::parse(line);
	    const std::string kind = record.at("record").get<std::string>();
	    const nlohmann::json &data = record.at("data");
	    if (kind == "schema") {
		schema = data.at("sanitized").get<bool>() &&
		    data.at("version").get<int>() == 1;
	    } else if (kind == "transition") {
		const nlohmann::json &camera = data.at("camera");
		if (camera.at("type").get<std::string>() == "orthographic") {
		    const nlohmann::json &position = camera.at("position");
		    nonFiniteCameraValid = nonFiniteCameraValid ||
			(position.size() == 3 && position.at(0).is_null() &&
			 position.at(1).is_null() && position.at(2).is_null());
		}
		const nlohmann::json &after = data.at("after");
		const nlohmann::json &signals = after.at("signals");
		const bool noProgressRoute = signals.at(
		    "handoff_without_progress_route").get<bool>();
		if (noProgressRoute &&
		    (after.at("host_work").at("pump_pending").get<bool>() ||
		     after.at("host_work").at("render_pending").get<bool>() ||
		     after.at("host_work").at("frame_claimed").get<bool>() ||
		     !signals.at("producer_idle").get<bool>() ||
		     !signals.at("renderer_prepared").get<bool>() ||
		     after.at("pending").at(
			 "refinement_cooldown_remaining_us").get<uint64_t>() > 0))
		    signalsConsistent = false;
		transition = data.at("transition_serial").get<uint64_t>() > 0 &&
		    after.at("control").contains("facts") &&
		    after.at("control").contains("owner_name") &&
		    after.at("progress").at("display").contains("class_name") &&
		    after.at("population").contains(
			"presented_structural_boxes") &&
		    after.at("allocation").contains("render_cost_budget") &&
		    after.at("allocation").contains(
			"interactive_calibrated_render_cost_per_second") &&
		    after.at("allocation").contains(
			"interactive_entry_render_cost_budget") &&
		    after.at("allocation").contains(
			"interactive_entry_progressive_cut_ceiling") &&
		    after.at("allocation").contains("retained_render_cost") &&
		    after.at("allocation").contains(
			"allocation_managed_render_cost") &&
		    after.at("allocation").contains(
			"allocation_unmanaged_render_cost") &&
		    after.at("allocation").contains(
			"maximum_nonaggregated_progressive_cut") &&
		after.at("allocation").contains("certificate_current") &&
		after.at("allocation").contains("cuts_applied") &&
		after.at("allocation").contains("active_mismatch_count") &&
		after.at("allocation").contains("presentation_realized") &&
		    after.at("work").contains("producer_stage") &&
		    after.at("pending").contains(
			"refinement_cooldown_remaining_us") &&
		    after.at("submission").contains("delta_entry_count") &&
		    after.at("submission").contains(
			"structural_frontier_count") &&
		    after.at("submission").contains(
			"structural_terminal_proxy") &&
		    after.at("submission").contains(
			"pass_missing_mesh_budget_blocked") &&
		    after.at("submission").contains("last_submitted_tasks") &&
		    after.at("submission").contains("last_updated_cuts") &&
		    signals.contains(
			"handoff_waiting_for_cooldown") &&
		    signals.contains(
			"handoff_without_render_or_producer") &&
		    signals.contains(
			"handoff_without_progress_route") &&
		    data.contains("camera");
	    } else if (kind == "inventory_source") {
		source = data.at("database").get<std::string>() == "db1.g" &&
		    data.at("path").get<std::string>() == "c1.c" &&
		    data.at("object_kind").get<std::string>() == "combination" &&
		    data.contains("compact_population_complete") &&
		    !data.at("compact_population_complete").get<bool>();
		initialInventoryReset = initialInventoryReset ||
		    data.at("inventory_reset").get<bool>();
		incrementalInventoryDelta = incrementalInventoryDelta ||
		    (!data.at("inventory_reset").get<bool>() &&
		     data.at("occurrence_count").get<uint64_t>() == 2);
	    } else if (kind == "inventory_object") {
		if (data.at("operation").get<std::string>() == "upsert") {
		    const bool expectedObject =
			data.at("path").get<std::string>() == "c1.c/s1.s" &&
			data.at("object").get<std::string>() == "s1.s" &&
			data.at("object_type").get<std::string>() == "bot" &&
			data.at("faces").get<uint64_t>() == 4 &&
			data.at("vertices").get<uint64_t>() == 4;
		    object = object || expectedObject;
		    hiddenObject = hiddenObject || (expectedObject &&
			!data.at("visible").get<bool>());
		    secondObject = secondObject ||
			(data.at("path").get<std::string>() == "c1.c/s2.s" &&
			 data.at("object").get<std::string>() == "s2.s");
		}
	    } else if (kind == "publication_results") {
		const nlohmann::json &provider = data.at("provider_status");
		const nlohmann::json &disposition =
		    data.at("publication_disposition");
		const nlohmann::json &replay = data.at("replay");
		const nlohmann::json &cursor =
		    data.at("submission_cursor");
		const nlohmann::json &samples = data.at("result_samples");
		const bool sampleValid = samples.size() == 1 &&
		    samples.at(0).at("source_entry_index").get<uint64_t>() == 0 &&
		    samples.at(0).at("source_routing_id").get<uint64_t>() > 0 &&
		    samples.at(0).at("request").at(
			"submission_reason_name").get<std::string>() ==
			    "spatial_presentation_repair" &&
		    samples.at(0).at("request").at(
			"draw_mode").get<int>() == BOBOL_LOD_DRAW_WIRE &&
		    samples.at(0).at("request").at(
			"required_chunk_count").get<uint64_t>() == 2 &&
		    samples.at(0).at("request").at(
			"required_chunk_hash").get<uint64_t>() != 0 &&
		    !samples.at(0).at("before").at("resident").get<bool>() &&
		    samples.at(0).at("outcome").at(
			"retry_current_demand").get<bool>() &&
		    !samples.at(0).at("outcome").at(
			"semantic_state_changed").get<bool>() &&
		    !samples.at(0).at("after").at("resident").get<bool>();
		publication = data.at("processed").get<uint64_t>() == 1 &&
		    data.at("matched").get<uint64_t>() == 1 &&
		    data.at("rejected").get<uint64_t>() == 1 &&
		    provider.at("cancelled").get<uint64_t>() == 1 &&
		    disposition.at("retry_current_demand").get<uint64_t>() == 1 &&
		    data.at("retry_source_entry_count").get<uint64_t>() == 1 &&
		    data.at("retry_source_entry_indices").size() == 1 &&
		    data.at("retry_source_entry_indices").at(0).get<uint64_t>() == 0 &&
		    data.at("result_sample_count").get<uint64_t>() == 1 &&
		    !data.at("result_samples_truncated").get<bool>() &&
		    sampleValid &&
		    replay.at("source").get<bool>() &&
		    !replay.at("requested").get<bool>() &&
		    !cursor.at("rescan_pending_after").get<bool>();
	    }
	}
    } catch (const std::exception &error) {
	std::fprintf(stderr, "FAIL: invalid telemetry JSONL: %s\n", error.what());
	return 1;
    }
    if (!nonFiniteCameraValid) {
	std::fprintf(stderr,
	    "FAIL: telemetry did not encode non-finite camera values as JSON null\n");
	return 1;
    }
    if (!schema || !transition || !source || !object || !hiddenObject ||
	!initialInventoryReset || !incrementalInventoryDelta || !secondObject ||
	!publication || !signalsConsistent) {
	std::fprintf(stderr,
	    "FAIL: telemetry omitted schema, control, or sanitized inventory/visibility facts\n");
	return 1;
    }
    return 0;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    if (argc != 1) {
	std::fprintf(stderr, "FAIL: unexpected arguments\n");
	return 1;
    }
    bobol_init(NULL);

    char databasePath[MAXPATHLEN] = {0};
    char sanitizedPath[MAXPATHLEN] = {0};
    char unsanitizedPath[MAXPATHLEN] = {0};
    FILE *temporary = bu_temp_file(databasePath, sizeof(databasePath));
    if (!temporary)
	return 1;
    std::fclose(temporary);
    temporary = bu_temp_file(sanitizedPath, sizeof(sanitizedPath));
    if (!temporary) {
	(void)bu_file_delete(databasePath);
	return 1;
    }
    std::fclose(temporary);
    temporary = bu_temp_file(unsanitizedPath, sizeof(unsanitizedPath));
    if (!temporary) {
	(void)bu_file_delete(databasePath);
	(void)bu_file_delete(sanitizedPath);
	return 1;
    }
    std::fclose(temporary);

    int result = 0;
    if (!make_database(databasePath)) {
	std::fprintf(stderr, "FAIL: could not create telemetry database\n");
	result = 1;
	goto cleanup_files;
    }
    struct db_i *database = db_open(databasePath, DB_OPEN_READONLY);
    if (!database || db_dirbuild(database) < 0) {
	std::fprintf(stderr, "FAIL: could not open telemetry database\n");
	if (database)
	    db_close(database);
	result = 1;
	goto cleanup_files;
    }

    {
	SoSeparator *root = new SoSeparator;
	root->ref();
	SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
	source->setDatabase(database);
	source->path = private_combination;
	source->instanceKey = private_instance;
	source->representationMode =
	    SoBRLDatabaseSource::REPRESENTATION_SHADED;
	source->drawMode = SoBRLDatabaseSource::SHADED;
	BObolCompactOccurrence occurrence;
	occurrence.geometry = make_geometry();
	occurrence.summary.valid = TRUE;
	occurrence.summary.shapeKind = BObolRealizedShapeSummary::SHAPE_MESH;
	const std::string occurrencePath =
	    std::string(private_combination) + "/" + private_solid;
	occurrence.summary.path = occurrencePath.c_str();
	occurrence.summary.sourceName = private_solid;
	occurrence.summary.sourceType = "bot";
	occurrence.summary.geometryKind = "aabb";
	occurrence.summary.visible = TRUE;
	occurrence.lodBacked = TRUE;
	occurrence.sourceMeshRequestValid = TRUE;
	occurrence.sourceMeshRequest.path = occurrence.summary.path;
	occurrence.sourceMeshRequest.sourceName = private_solid;
	occurrence.sourceMeshRequest.sourceType = "bot";
	occurrence.sourceMeshRequest.meshAssetPath = private_solid;
	occurrence.sourceMeshRequest.meshAssetName = private_solid;
	occurrence.sourceMeshRequest.faceCount = 4;
	occurrence.sourceMeshRequest.pointCount = 4;
	occurrence.sourceMeshRequest.bounds = SbBox3f(
	    SbVec3f(0.0f, 0.0f, 0.0f), SbVec3f(1.0f, 1.0f, 1.0f));
	occurrence.sourceMeshRequest.meshAssetBounds =
	    occurrence.sourceMeshRequest.bounds;
	if (!occurrence.geometry || source->setCompactOccurrenceRegistry(
		std::vector<BObolCompactOccurrence>(1, occurrence)) != 1) {
	    std::fprintf(stderr, "FAIL: could not publish telemetry fixture\n");
	    root->unref();
	    db_close(database);
	    result = 1;
	    goto cleanup_files;
	}
	root->addChild(source);

	bu_setenv(telemetry_file_environment, sanitizedPath, 1);
	bu_setenv(telemetry_sanitize_environment, "1", 1);
	{
	    BObolViewController controller(root, NULL);
	    SoOrthographicCamera *camera = new SoOrthographicCamera;
	    controller.setCamera(camera);
	    camera->position = SbVec3f(
		std::numeric_limits<float>::quiet_NaN(),
		std::numeric_limits<float>::infinity(),
		-std::numeric_limits<float>::infinity());
	    controller.clearRenderRequest();
	    camera->position = SbVec3f(0.0f, 0.0f, 10.0f);
	    controller.setLodControlTransitionTracing(TRUE, 16);
	    controller.setViewportSize(320, 240);
	    if (source->setCompactInstanceVisibilityOverrideForPathMatch(
		    occurrencePath.c_str(), BOBOL_COMPACT_PATH_EXACT, FALSE) != 1) {
		std::fprintf(stderr,
		    "FAIL: could not publish telemetry visibility fixture\n");
		result = 1;
	    }
	    BObolCompactOccurrence secondOccurrence = occurrence;
	    const std::string secondPath =
		std::string(private_combination) + "/" + private_second_solid;
	    secondOccurrence.summary.path = secondPath.c_str();
	    secondOccurrence.summary.sourceName = private_second_solid;
	    secondOccurrence.sourceMeshRequest.path = secondPath.c_str();
	    secondOccurrence.sourceMeshRequest.sourceName =
		private_second_solid;
	    secondOccurrence.sourceMeshRequest.meshAssetPath =
		private_second_solid;
	    secondOccurrence.sourceMeshRequest.meshAssetName =
		private_second_solid;
	    if (source->mergeCompactOccurrences({secondOccurrence}, TRUE) != 1) {
		std::fprintf(stderr,
		    "FAIL: could not publish incremental telemetry fixture\n");
		result = 1;
	    }
	    controller.setViewportSize(321, 240);

	    BObolCompactInstanceHandle handle;
	    BObolCompactInstanceSummary summary;
	    BObolLodService service;
	    if (!source->getCompactInstanceHandle(0, handle) ||
		!source->getCompactInstanceSummary(handle, summary) ||
		summary.sourceInstanceKey.getLength() == 0 ||
		!service.start(1, FALSE)) {
		std::fprintf(stderr,
		    "FAIL: could not prepare telemetry publication fixture\n");
		result = 1;
	    } else {
		controller.setLodAutoSubmit(FALSE);
		controller.setLodService(&service);
		BObolLodTask task;
		task.generation = controller.beginLodGeneration();
		task.request.databaseId = "private-database-id";
		task.request.databaseRevision = 1;
		task.request.sourceRevision = source->sourceRevision.getValue();
		task.request.sourceContentHash = 1234;
		task.request.objectPath = occurrencePath.c_str();
		task.request.objectName = private_solid;
		task.request.occurrenceKey = summary.sourceInstanceKey;
		task.request.sourceRoutingId =
		    source->getCompactSourceRoutingId();
		task.request.sourcePopulationEpoch =
		    source->getCompactPopulationEpoch();
		task.request.sourceEntryIndex = 0;
		task.request.viewRevision = controller.getLodViewRevision();
		task.request.policyRevision = controller.getLodPolicyRevision();
		task.request.submissionReason =
		    BOBOL_LOD_SUBMISSION_SPATIAL_PRESENTATION_REPAIR;
		task.request.drawMode = BOBOL_LOD_DRAW_WIRE;
		task.request.requiredChunks = {7, 11};
		task.request.providerId = "private-provider";
		task.request.providerVersion = "private-version";
		task.realize = cancelled_result;
		if (service.submit(task) == 0) {
		    result = 1;
		} else {
		    for (int attempt = 0;
			 attempt < 100 &&
			 service.queuedResultCountForDiagnostics() == 0;
			 ++attempt)
			std::this_thread::sleep_for(
			    std::chrono::milliseconds(5));
		    if (service.queuedResultCountForDiagnostics() == 0 ||
			controller.applyLodResults(
			    &service, 1, 0, task.generation) != 0 ||
			controller.getLastLodMatchedResultCount() != 1 ||
			controller.getLastLodRejectedResultCount() != 1) {
			std::fprintf(stderr,
			    "FAIL: telemetry publication fixture did not drain a retry\n");
			result = 1;
		    }
		}
		controller.setLodService(NULL);
		service.stop();
	    }
	    std::vector<BObolLodControlTransitionRecord> records;
	    controller.drainLodControlTransitions(records);
	    if (records.empty() || records.front().serial != 1 ||
		records.front().event != BOBOL_LOD_CONTROL_TRANSITION_INITIAL) {
		std::fprintf(stderr,
		    "FAIL: telemetry disturbed the bounded transition journal\n");
		result = 1;
	    }
	}
	const std::string sanitized = read_file(sanitizedPath);
	if (!result)
	    result = validate_sanitized_log(sanitized,
		path_basename(databasePath));

	bu_setenv(telemetry_file_environment, unsanitizedPath, 1);
	bu_setenv(telemetry_sanitize_environment, "0", 1);
	{
	    BObolViewController controller(root, NULL);
	    controller.setViewportSize(640, 480);
	}
	const std::string unsanitized = read_file(unsanitizedPath);
	if (!result &&
	    (unsanitized.find(private_combination) == std::string::npos ||
	     unsanitized.find(private_solid) == std::string::npos ||
	     unsanitized.find(private_instance) == std::string::npos ||
	     unsanitized.find(path_basename(databasePath)) == std::string::npos ||
	     unsanitized.find("\"render_reason\":") == std::string::npos)) {
	    std::fprintf(stderr,
		"FAIL: unsanitized telemetry omitted diagnostic identity\n");
	    result = 1;
	}

	bu_setenv(telemetry_file_environment, "", 1);
	bu_setenv(telemetry_sanitize_environment, "", 1);
	root->unref();
    }
    db_close(database);

cleanup_files:
    (void)bu_file_delete(databasePath);
    (void)bu_file_delete(sanitizedPath);
    (void)bu_file_delete(unsanitizedPath);
    return result;
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8 cino=N-s
 */
