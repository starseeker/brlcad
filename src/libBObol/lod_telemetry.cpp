/*                    L O D _ T E L E M E T R Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/log.h"
#include "bu/str.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BViewController.h"
#include "lod_control_private.h"
#include "lod_telemetry_private.h"
#include "raytrace.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <locale>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <Inventor/SbRotation.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>

namespace {

static constexpr int bobol_lod_telemetry_schema_version = 1;
static const char *bobol_lod_telemetry_file_environment =
    "BOBOL_LOD_TELEMETRY_FILE";
static const char *bobol_lod_telemetry_sanitize_environment =
    "BOBOL_LOD_TELEMETRY_SANITIZE";

static void
telemetry_append_json_string(std::string &output, const char *value)
{
    static const char hexadecimal[] = "0123456789abcdef";
    output.push_back('"');
    const unsigned char *current = reinterpret_cast<const unsigned char *>(
	value ? value : "");
    while (*current) {
	const unsigned char c = *current++;
	switch (c) {
	    case '"': output += "\\\""; break;
	    case '\\': output += "\\\\"; break;
	    case '\b': output += "\\b"; break;
	    case '\f': output += "\\f"; break;
	    case '\n': output += "\\n"; break;
	    case '\r': output += "\\r"; break;
	    case '\t': output += "\\t"; break;
	    default:
		if (c < 0x20) {
		    output += "\\u00";
		    output.push_back(hexadecimal[(c >> 4) & 0xf]);
		    output.push_back(hexadecimal[c & 0xf]);
		} else {
		    output.push_back(static_cast<char>(c));
		}
		break;
	}
    }
    output.push_back('"');
}

static void
telemetry_append_json_real(std::string &output, double value, int precision)
{
    if (!std::isfinite(value)) {
	output += "null";
	return;
    }
    std::ostringstream formatted;
    formatted.imbue(std::locale::classic());
    formatted << std::setprecision(precision) << value;
    output += formatted.str();
}

class TelemetryJsonObject {
public:
    TelemetryJsonObject()
    {
	this->value.push_back('{');
    }

    void text(const char *name, const char *fieldValue)
    {
	this->key(name);
	telemetry_append_json_string(this->value, fieldValue);
    }

    void text(const char *name, const std::string &fieldValue)
    {
	this->text(name, fieldValue.c_str());
    }

    template <typename Integer>
    void integer(const char *name, Integer fieldValue)
    {
	this->key(name);
	this->value += std::to_string(fieldValue);
    }

    void boolean(const char *name, bool fieldValue)
    {
	this->key(name);
	this->value += fieldValue ? "true" : "false";
    }

    void real(const char *name, double fieldValue)
    {
	this->key(name);
	telemetry_append_json_real(this->value, fieldValue,
	    std::numeric_limits<double>::max_digits10);
    }

    void raw(const char *name, const std::string &fieldValue)
    {
	this->key(name);
	this->value += fieldValue;
    }

    std::string take()
    {
	this->value.push_back('}');
	return std::move(this->value);
    }

private:
    void key(const char *name)
    {
	if (!this->first)
	    this->value.push_back(',');
	this->first = false;
	telemetry_append_json_string(this->value, name);
	this->value.push_back(':');
    }

    std::string value;
    bool first = true;
};

static std::string
telemetry_size_array(const size_t *values, size_t count)
{
    std::string result("[");
    for (size_t i = 0; i < count; ++i) {
	if (i)
	    result.push_back(',');
	result += std::to_string(values[i]);
    }
    result.push_back(']');
    return result;
}

static std::string
telemetry_u32_array(const std::vector<uint32_t> &values)
{
    std::string result("[");
    for (size_t i = 0; i < values.size(); ++i) {
	if (i)
	    result.push_back(',');
	result += std::to_string(values[i]);
    }
    result.push_back(']');
    return result;
}

static std::string
telemetry_provider_status_counts(
    const BObolLodPublicationTelemetryRecord &publication)
{
    TelemetryJsonObject result;
    result.integer("unknown", publication.providerStatusCounts[
	BOBOL_LOD_PROVIDER_UNKNOWN]);
    result.integer("ready", publication.providerStatusCounts[
	BOBOL_LOD_PROVIDER_READY]);
    result.integer("cache_miss", publication.providerStatusCounts[
	BOBOL_LOD_PROVIDER_CACHE_MISS]);
    result.integer("stale", publication.providerStatusCounts[
	BOBOL_LOD_PROVIDER_STALE]);
    result.integer("running", publication.providerStatusCounts[
	BOBOL_LOD_PROVIDER_RUNNING]);
    result.integer("terminal", publication.providerStatusCounts[
	BOBOL_LOD_PROVIDER_TERMINAL]);
    result.integer("fallback", publication.providerStatusCounts[
	BOBOL_LOD_PROVIDER_FALLBACK]);
    result.integer("error", publication.providerStatusCounts[
	BOBOL_LOD_PROVIDER_ERROR]);
    result.integer("cancelled", publication.providerStatusCounts[
	BOBOL_LOD_PROVIDER_CANCELLED]);
    result.integer("superseded", publication.providerStatusCounts[
	BOBOL_LOD_PROVIDER_SUPERSEDED]);
    result.integer("out_of_range", publication.unknownProviderStatusCount);
    return result.take();
}

static const char *
telemetry_submission_reason_name(int reason)
{
    switch (reason) {
	case BOBOL_LOD_SUBMISSION_INITIAL: return "initial";
	case BOBOL_LOD_SUBMISSION_ASSET_REPLACEMENT:
	    return "asset_replacement";
	case BOBOL_LOD_SUBMISSION_NORMAL_PRESENTATION_REPAIR:
	    return "normal_presentation_repair";
	case BOBOL_LOD_SUBMISSION_SPATIAL_PRESENTATION_REPAIR:
	    return "spatial_presentation_repair";
	case BOBOL_LOD_SUBMISSION_PRESENTATION_REFINEMENT:
	    return "presentation_refinement";
	case BOBOL_LOD_SUBMISSION_RESIDENT_PREFETCH:
	    return "resident_prefetch";
	case BOBOL_LOD_SUBMISSION_PRESENTATION_AND_RESIDENT_PREFETCH:
	    return "presentation_and_resident_prefetch";
	case BOBOL_LOD_SUBMISSION_FORCED: return "forced";
	case BOBOL_LOD_SUBMISSION_RESET: return "reset";
	case BOBOL_LOD_SUBMISSION_TERMINAL_PROMOTION:
	    return "terminal_promotion";
	case BOBOL_LOD_SUBMISSION_UNSPECIFIED:
	default: return "unspecified";
    }
}

static std::string
telemetry_result_sample(
    const BObolLodPublicationTelemetryRecord::ResultSample &sample)
{
    TelemetryJsonObject result;
    result.integer("source_routing_id", sample.sourceRoutingId);
    result.integer("source_entry_index", sample.sourceEntryIndex);

    TelemetryJsonObject request;
    request.integer("submission_reason", sample.submissionReason);
    request.text("submission_reason_name",
	telemetry_submission_reason_name(sample.submissionReason));
    request.integer("draw_mode", sample.drawMode);
    request.integer("normal_style", sample.normalStyle);
    request.integer("requested_cut", sample.requestCut);
    request.integer("required_chunk_count", sample.requestChunkCount);
    request.integer("required_chunk_hash", sample.requestChunkHash);
    result.raw("request", request.take());

    TelemetryJsonObject incoming;
    incoming.integer("resolved_cut", sample.resolvedCut);
    incoming.integer("active_cut", sample.incomingActiveCut);
    incoming.integer("resident_cut", sample.incomingResidentCut);
    incoming.integer("presentation_admission_cut",
	sample.incomingPresentationAdmissionCut);
    incoming.integer("presentation_layer_count",
	sample.incomingLayerCount);
    incoming.integer("mesh_revision", sample.incomingMeshRevision);
    incoming.integer("prepared_revision",
	sample.incomingPreparedRevision);
    incoming.boolean("prepared_geometry",
	sample.incomingPreparedGeometry);
    incoming.integer("faces", sample.incomingFaceCount);
    incoming.integer("points", sample.incomingPointCount);
    incoming.boolean("terminal", sample.incomingTerminal);
    incoming.boolean("memory_limited", sample.incomingMemoryLimited);
    result.raw("incoming", incoming.take());

    TelemetryJsonObject before;
    before.boolean("resident", sample.residentBefore);
    before.integer("active_cut", sample.beforeActiveCut);
    before.integer("resident_cut", sample.beforeResidentCut);
    before.integer("requested_cut", sample.beforeRequestedCut);
    before.integer("allocated_cut", sample.beforeAllocatedCut);
    before.integer("normal_style", sample.beforeNormalStyle);
    before.integer("required_chunk_count",
	sample.beforeRequiredChunkCount);
    before.integer("required_chunk_hash",
	sample.beforeRequiredChunkHash);
    before.integer("presented_chunk_count",
	sample.beforePresentedChunkCount);
    before.integer("presented_chunk_hash",
	sample.beforePresentedChunkHash);
    before.integer("presentation_layer_count", sample.beforeLayerCount);
    before.integer("mesh_revision", sample.beforeMeshRevision);
    before.integer("prepared_revision", sample.beforePreparedRevision);
    before.boolean("prepared_geometry", sample.beforePreparedGeometry);
    before.integer("faces", sample.beforeFaceCount);
    before.integer("points", sample.beforePointCount);
    before.boolean("normal_presentation_matches",
	sample.beforeNormalPresentationMatches);
    before.boolean("requested_cut_drawable",
	sample.beforeRequestedCutDrawable);
    before.boolean("allocated_cut_drawable",
	sample.beforeAllocatedCutDrawable);
    result.raw("before", before.take());

    TelemetryJsonObject outcome;
    outcome.boolean("accepted", sample.accepted);
    outcome.boolean("unchanged", sample.unchanged);
    outcome.boolean("retry_current_demand", sample.retryCurrentDemand);
    outcome.boolean("semantic_state_changed", sample.semanticStateChanged);
    outcome.boolean("same_progressive_asset",
	sample.sameProgressiveAsset);
    result.raw("outcome", outcome.take());

    TelemetryJsonObject after;
    after.boolean("resident", sample.residentAfter);
    after.integer("active_cut", sample.afterActiveCut);
    after.integer("resident_cut", sample.afterResidentCut);
    after.integer("requested_cut", sample.afterRequestedCut);
    after.integer("allocated_cut", sample.afterAllocatedCut);
    after.integer("normal_style", sample.afterNormalStyle);
    after.integer("required_chunk_count", sample.afterRequiredChunkCount);
    after.integer("required_chunk_hash", sample.afterRequiredChunkHash);
    after.integer("presented_chunk_count",
	sample.afterPresentedChunkCount);
    after.integer("presented_chunk_hash",
	sample.afterPresentedChunkHash);
    after.integer("presentation_layer_count", sample.afterLayerCount);
    after.integer("mesh_revision", sample.afterMeshRevision);
    after.integer("prepared_revision", sample.afterPreparedRevision);
    after.boolean("prepared_geometry", sample.afterPreparedGeometry);
    after.integer("faces", sample.afterFaceCount);
    after.integer("points", sample.afterPointCount);
    after.boolean("normal_presentation_matches",
	sample.afterNormalPresentationMatches);
    after.boolean("requested_cut_drawable",
	sample.afterRequestedCutDrawable);
    after.boolean("allocated_cut_drawable",
	sample.afterAllocatedCutDrawable);
    result.raw("after", after.take());
    return result.take();
}

static std::string
telemetry_result_samples(
    const std::vector<BObolLodPublicationTelemetryRecord::ResultSample> &samples)
{
    std::string result("[");
    for (size_t i = 0; i < samples.size(); ++i) {
	if (i)
	    result.push_back(',');
	result += telemetry_result_sample(samples[i]);
    }
    result.push_back(']');
    return result;
}

struct TelemetryBitName {
    uint32_t bit;
    const char *name;
};

template <size_t Count>
static std::string
telemetry_named_bits(uint32_t mask,
    const TelemetryBitName (&names)[Count])
{
    std::string result("[");
    bool first = true;
    for (const TelemetryBitName &entry : names) {
	if (!(mask & entry.bit))
	    continue;
	if (!first)
	    result.push_back(',');
	first = false;
	telemetry_append_json_string(result, entry.name);
    }
    result.push_back(']');
    return result;
}

static const char *
telemetry_control_owner_name(int owner)
{
    switch (owner) {
	case BOBOL_LOD_CONTROL_OWNER_INTERACTION: return "interaction";
	case BOBOL_LOD_CONTROL_OWNER_INVENTORY: return "inventory";
	case BOBOL_LOD_CONTROL_OWNER_AVAILABILITY: return "availability";
	case BOBOL_LOD_CONTROL_OWNER_PUBLICATION: return "publication";
	case BOBOL_LOD_CONTROL_OWNER_PLANNING: return "planning";
	case BOBOL_LOD_CONTROL_OWNER_PRESENTATION: return "presentation";
	case BOBOL_LOD_CONTROL_OWNER_HANDOFF: return "handoff";
	case BOBOL_LOD_CONTROL_OWNER_COMPACTION: return "compaction";
	case BOBOL_LOD_CONTROL_OWNER_CACHE_WRITE: return "cache_write";
	case BOBOL_LOD_CONTROL_OWNER_NONE:
	default: return "none";
    }
}

static const char *
telemetry_phase_name(int phase)
{
    switch (phase) {
	case BOBOL_LOD_CONVERGENCE_IDLE: return "idle";
	case BOBOL_LOD_CONVERGENCE_DISCOVERING: return "discovering";
	case BOBOL_LOD_CONVERGENCE_PREPARING: return "preparing";
	case BOBOL_LOD_CONVERGENCE_INTERACTIVE: return "interactive";
	case BOBOL_LOD_CONVERGENCE_REFINING: return "refining";
	case BOBOL_LOD_CONVERGENCE_CALIBRATING: return "calibrating";
	case BOBOL_LOD_CONVERGENCE_BACKGROUND: return "background";
	case BOBOL_LOD_CONVERGENCE_ERROR: return "error";
	default: return "unknown";
    }
}

static const char *
telemetry_outcome_name(int outcome)
{
    switch (outcome) {
	case BOBOL_LOD_PRESENTATION_ACTIVE: return "active";
	case BOBOL_LOD_PRESENTATION_READY: return "ready";
	case BOBOL_LOD_PRESENTATION_CONSTRAINED: return "constrained";
	case BOBOL_LOD_PRESENTATION_ERROR: return "error";
	default: return "unknown";
    }
}

static const char *
telemetry_progress_display_class_name(BObolLodProgressDisplayClass displayClass)
{
    switch (displayClass) {
	case BOBOL_LOD_PROGRESS_DISPLAY_IDLE: return "idle";
	case BOBOL_LOD_PROGRESS_DISPLAY_DISCOVERING: return "discovering";
	case BOBOL_LOD_PROGRESS_DISPLAY_PREPARING: return "preparing";
	case BOBOL_LOD_PROGRESS_DISPLAY_INTERACTIVE: return "interactive";
	case BOBOL_LOD_PROGRESS_DISPLAY_SETTLING: return "settling";
	case BOBOL_LOD_PROGRESS_DISPLAY_BACKGROUND: return "background";
	case BOBOL_LOD_PROGRESS_DISPLAY_ERROR: return "error";
	case BOBOL_LOD_PROGRESS_DISPLAY_TERMINAL_ERROR:
	    return "terminal_error";
	default: return "unknown";
    }
}

static const char *
telemetry_producer_stage_name(int stage)
{
    switch (stage) {
	case BOBOL_LOD_PRODUCER_STAGE_CACHE_LOOKUP: return "cache_lookup";
	case BOBOL_LOD_PRODUCER_STAGE_SOURCE_PREPARATION:
	    return "source_preparation";
	case BOBOL_LOD_PRODUCER_STAGE_SOURCE_HASHING: return "source_hashing";
	case BOBOL_LOD_PRODUCER_STAGE_BOUNDS_ANALYSIS:
	    return "bounds_analysis";
	case BOBOL_LOD_PRODUCER_STAGE_FACE_CLASSIFICATION:
	    return "face_classification";
	case BOBOL_LOD_PRODUCER_STAGE_PREFIX_MATERIALIZATION:
	    return "prefix_materialization";
	case BOBOL_LOD_PRODUCER_STAGE_SPATIAL_CONSTRUCTION:
	    return "spatial_construction";
	case BOBOL_LOD_PRODUCER_STAGE_CACHE_PERSISTENCE:
	    return "cache_persistence";
	case BOBOL_LOD_PRODUCER_STAGE_COVERAGE_PREVIEW:
	    return "coverage_preview";
	case BOBOL_LOD_PRODUCER_STAGE_ASSET_SERIALIZATION:
	    return "asset_serialization";
	case BOBOL_LOD_PRODUCER_STAGE_NONE:
	default: return "none";
    }
}

static std::string
telemetry_vec3(const SbVec3f &value)
{
    std::string result("[");
    for (int i = 0; i < 3; ++i) {
	if (i)
	    result.push_back(',');
	telemetry_append_json_real(result, value[i],
	    std::numeric_limits<float>::max_digits10);
    }
    result.push_back(']');
    return result;
}

static std::string
telemetry_proxy_reasons(const BObolLodProxyReasonStatus &reasons)
{
    TelemetryJsonObject result;
    result.integer("mask", reasons.reasonMask);
    result.integer("source_preparation",
	reasons.sourcePreparationOccurrenceCount);
    result.integer("visibility_planning",
	reasons.visibilityPlanningOccurrenceCount);
    result.integer("geometry_preparation",
	reasons.geometryPreparationOccurrenceCount);
    result.integer("renderer_preparation",
	reasons.rendererPreparationOccurrenceCount);
    result.integer("intentional_subpixel",
	reasons.intentionalSubpixelOccurrenceCount);
    result.integer("frame_budget", reasons.frameBudgetOccurrenceCount);
    result.integer("memory_budget", reasons.memoryBudgetOccurrenceCount);
    result.integer("terminal_failure",
	reasons.terminalFailureOccurrenceCount);
    result.integer("unclassified", reasons.unclassifiedOccurrenceCount);
    result.integer("temporary", reasons.temporaryStructuralOccurrenceCount());
    result.integer("budget_limited",
	reasons.budgetLimitedStructuralOccurrenceCount());
    return result.take();
}

static std::string
telemetry_episode(const BObolLodEpisodeStatus &episode)
{
    TelemetryJsonObject result;
    result.integer("elapsed_ms", episode.elapsedMilliseconds);
    result.boolean("first_proxy_reached", episode.firstProxyReached != FALSE);
    result.integer("first_proxy_ms", episode.firstProxyMilliseconds);
    result.boolean("first_mesh_reached", episode.firstMeshReached != FALSE);
    result.integer("first_mesh_ms", episode.firstMeshMilliseconds);
    result.integer("structural_proxy_baseline",
	episode.structuralProxyBaselineCount);
    result.boolean("half_structural_proxies_replaced",
	episode.halfStructuralProxiesReplaced != FALSE);
    result.integer("half_structural_proxy_replacement_ms",
	episode.halfStructuralProxyReplacementMilliseconds);
    result.boolean("stable_view_reached", episode.stableViewReached != FALSE);
    result.integer("stable_view_ms", episode.stableViewMilliseconds);
    return result.take();
}

static std::string
telemetry_host_work(const BObolHostWorkSnapshot &host)
{
    TelemetryJsonObject result;
    result.integer("revision", host.revision);
    result.integer("render_revision", host.renderRevision);
    result.integer("flags", host.flags);
    result.boolean("pump_pending", host.pumpPending() != FALSE);
    result.boolean("render_pending", host.renderPending() != FALSE);
    result.boolean("capacity_sample_requested",
	host.capacitySampleRequested() != FALSE);
    result.boolean("frame_claimed", host.frameClaimed() != FALSE);
    result.boolean("capacity_sample_claimed",
	host.capacitySampleClaimed() != FALSE);
    return result.take();
}

static std::string
telemetry_control(const BObolLodConvergenceStatus &status)
{
    using Refinement = BObolLodControlRefinement;
    static const TelemetryBitName factNames[] = {
	{static_cast<uint32_t>(Refinement::Fact::INTERACTION), "interaction"},
	{static_cast<uint32_t>(Refinement::Fact::INVENTORY), "inventory"},
	{static_cast<uint32_t>(Refinement::Fact::AVAILABILITY), "availability"},
	{static_cast<uint32_t>(Refinement::Fact::RESULT), "result"},
	{static_cast<uint32_t>(Refinement::Fact::PUBLICATION), "publication"},
	{static_cast<uint32_t>(Refinement::Fact::SUBMISSION), "submission"},
	{static_cast<uint32_t>(Refinement::Fact::SUBMISSION_RESCAN),
	    "submission_rescan"},
	{static_cast<uint32_t>(Refinement::Fact::SUBMISSION_DELTA),
	    "submission_delta"},
	{static_cast<uint32_t>(Refinement::Fact::QUALITY_PROBE), "quality_probe"},
	{static_cast<uint32_t>(Refinement::Fact::RETAINED_ALLOCATION),
	    "retained_allocation"},
	{static_cast<uint32_t>(Refinement::Fact::RETAINED_ALLOCATION_TRANSACTION),
	    "retained_allocation_transaction"},
	{static_cast<uint32_t>(Refinement::Fact::IMPORTANCE_CENSUS),
	    "importance_census"},
	{static_cast<uint32_t>(Refinement::Fact::RESIDENT_ADMISSION_RETRY),
	    "resident_admission_retry"},
	{static_cast<uint32_t>(Refinement::Fact::CAPACITY_ALLOCATION),
	    "capacity_allocation"},
	{static_cast<uint32_t>(Refinement::Fact::RESIDENT_GROWTH),
	    "resident_growth"},
	{static_cast<uint32_t>(Refinement::Fact::POINT_TRIANGLE_RECOVERY),
	    "point_triangle_recovery"},
	{static_cast<uint32_t>(Refinement::Fact::STRUCTURAL_FRONTIER),
	    "structural_frontier"},
	{static_cast<uint32_t>(Refinement::Fact::PRESENTATION_REPLAY),
	    "presentation_replay"},
	{static_cast<uint32_t>(Refinement::Fact::PRESENTATION_BARRIER),
	    "presentation_barrier"},
	{static_cast<uint32_t>(Refinement::Fact::CAPACITY_FRAME),
	    "capacity_frame"},
	{static_cast<uint32_t>(Refinement::Fact::POINT_ADMISSION_FRAME),
	    "point_admission_frame"},
	{static_cast<uint32_t>(Refinement::Fact::POINT_CALIBRATION),
	    "point_calibration"},
	{static_cast<uint32_t>(Refinement::Fact::CAPACITY_CALIBRATION),
	    "capacity_calibration"},
	{static_cast<uint32_t>(Refinement::Fact::HEADROOM_PROBE),
	    "headroom_probe"},
	{static_cast<uint32_t>(Refinement::Fact::HANDOFF), "handoff"},
	{static_cast<uint32_t>(Refinement::Fact::COMPACTION), "compaction"},
	{static_cast<uint32_t>(Refinement::Fact::CACHE_WRITE), "cache_write"},
	{static_cast<uint32_t>(Refinement::Fact::DEMAND_REFRESH),
	    "demand_refresh"},
	{static_cast<uint32_t>(Refinement::Fact::EXACT_PRESENTATION),
	    "exact_presentation"}
    };
    static const TelemetryBitName obligationNames[] = {
	{BOBOL_LOD_CONTROL_OBLIGATION_INTERACTION, "interaction"},
	{BOBOL_LOD_CONTROL_OBLIGATION_INVENTORY, "inventory"},
	{BOBOL_LOD_CONTROL_OBLIGATION_AVAILABILITY, "availability"},
	{BOBOL_LOD_CONTROL_OBLIGATION_PUBLICATION, "publication"},
	{BOBOL_LOD_CONTROL_OBLIGATION_PLANNING, "planning"},
	{BOBOL_LOD_CONTROL_OBLIGATION_PRESENTATION, "presentation"},
	{BOBOL_LOD_CONTROL_OBLIGATION_HANDOFF, "handoff"},
	{BOBOL_LOD_CONTROL_OBLIGATION_COMPACTION, "compaction"},
	{BOBOL_LOD_CONTROL_OBLIGATION_CACHE_WRITE, "cache_write"}
    };
    static const TelemetryBitName violationNames[] = {
	{BOBOL_LOD_CONTROL_VIOLATION_OWNERLESS_WORK, "ownerless_work"},
	{BOBOL_LOD_CONTROL_VIOLATION_TERMINAL_WITH_WORK,
	    "terminal_with_work"},
	{BOBOL_LOD_CONTROL_VIOLATION_INVALID_READINESS,
	    "invalid_readiness"},
	{BOBOL_LOD_CONTROL_VIOLATION_INVALID_OWNER, "invalid_owner"},
	{BOBOL_LOD_CONTROL_VIOLATION_UNWITNESSED_PRESENTATION,
	    "unwitnessed_presentation"},
	{BOBOL_LOD_CONTROL_VIOLATION_UNWITNESSED_CONSTRAINT,
	    "unwitnessed_constraint"},
	{BOBOL_LOD_CONTROL_VIOLATION_UNWITNESSED_PLANNING,
	    "unwitnessed_planning"},
	{BOBOL_LOD_CONTROL_VIOLATION_NONTERMINAL_WITHOUT_PROGRESS,
	    "nonterminal_without_progress"}
    };
    static const TelemetryBitName witnessNames[] = {
	{BOBOL_LOD_PRESENTATION_WITNESS_RENDER, "render"},
	{BOBOL_LOD_PRESENTATION_WITNESS_CONTROLLER_PUMP, "controller_pump"},
	{BOBOL_LOD_PRESENTATION_WITNESS_TIMER, "timer"},
	{BOBOL_LOD_PRESENTATION_WITNESS_INDEPENDENT_PRODUCER,
	    "independent_producer"},
	{BOBOL_LOD_PRESENTATION_WITNESS_CLAIMED_FRAME, "claimed_frame"}
    };
    static const TelemetryBitName constraintNames[] = {
	{BOBOL_LOD_CONSTRAINT_STABLE_BUDGET, "stable_budget"},
	{BOBOL_LOD_CONSTRAINT_PROGRESSIVE_CEILING, "progressive_ceiling"},
	{BOBOL_LOD_CONSTRAINT_SUBPIXEL_AGGREGATION,
	    "subpixel_aggregation"},
	{BOBOL_LOD_CONSTRAINT_STATIC_DEADLINE, "static_deadline"},
	{BOBOL_LOD_CONSTRAINT_MEMORY, "memory"},
	{BOBOL_LOD_CONSTRAINT_TERMINAL_PROXY, "terminal_proxy"}
    };
    TelemetryJsonObject result;
    result.integer("fact_mask", status.controlFactMask);
    result.raw("facts", telemetry_named_bits(status.controlFactMask,
	factNames));
    result.integer("obligation_mask", status.controlObligationMask);
    result.raw("obligations", telemetry_named_bits(
	status.controlObligationMask, obligationNames));
    result.integer("owner", status.controlOwner);
    result.text("owner_name", telemetry_control_owner_name(
	status.controlOwner));
    result.integer("violation_mask", status.controlViolationMask);
    result.raw("violations", telemetry_named_bits(
	status.controlViolationMask, violationNames));
    result.integer("presentation_witness_mask",
	status.controlPresentationWitnessMask);
    result.raw("presentation_witnesses", telemetry_named_bits(
	status.controlPresentationWitnessMask, witnessNames));
    result.integer("constraint_evidence_mask", status.constraintEvidenceMask);
    result.raw("constraint_evidence", telemetry_named_bits(
	status.constraintEvidenceMask, constraintNames));
    return result.take();
}

static std::string
telemetry_revisions(const BObolLodConvergenceStatus &status)
{
    TelemetryJsonObject result;
    result.integer("inventory", status.inventoryRevision);
    result.integer("availability", status.availabilityRevision);
    result.integer("visibility", status.visibilityRevision);
    result.integer("view", status.viewRevision);
    result.integer("policy", status.policyRevision);
    result.integer("episode", status.episodeRevision);
    result.integer("capacity", status.capacityRevision);
    result.integer("cad", status.cadRevision);
    result.integer("resident_demand", status.residentDemandRevision);
    result.integer("service_resident_admission",
	status.serviceResidentAdmissionRevision);
    result.integer("observed_resident_admission",
	status.observedResidentAdmissionRevision);
    return result.take();
}

static std::string
telemetry_progress(const BObolLodConvergenceStatus &status)
{
    TelemetryJsonObject result;
    const BObolLodProgressDisplayStatus display =
	status.progressDisplayStatus();
    result.integer("phase", status.phase);
    result.text("phase_name", telemetry_phase_name(status.phase));
    result.integer("outcome", status.outcome);
    result.text("outcome_name", telemetry_outcome_name(status.outcome));
    result.real("fraction", status.fraction);
    result.boolean("estimate_available",
	status.progressEstimateAvailable != FALSE);
    result.boolean("estimate_refinement_cycle_based",
	status.progressEstimateRefinementCycleBased != FALSE);
    result.real("estimated_fraction", status.estimatedFraction);
    result.integer("estimated_remaining_ms",
	status.estimatedRemainingMilliseconds);
    result.integer("estimated_remaining_refinement_cycles",
	status.estimatedRemainingRefinementCycles);
    result.boolean("terminal", status.terminal != FALSE);
    result.boolean("terminal_error", status.terminalError != FALSE);
    result.boolean("view_ready", status.viewReady != FALSE);
    result.boolean("has_lod_state", status.hasLodState != FALSE);
    result.boolean("background_pending", status.backgroundPending != FALSE);
    result.boolean("performance_limited",
	status.performanceLimited != FALSE);
    result.boolean("memory_limited", status.memoryLimited != FALSE);
    result.boolean("gpu_memory_pressure",
	status.gpuMemoryPressure != FALSE);
    result.integer("failed_source_count", status.failedSourceCount);
    TelemetryJsonObject displayData;
    displayData.boolean("visible", display.visible != FALSE);
    displayData.boolean("terminal_ready", display.terminalReady != FALSE);
    displayData.integer("class", display.publicationClass);
    displayData.text("class_name", telemetry_progress_display_class_name(
	display.publicationClass));
    result.raw("display", displayData.take());
    result.raw("episode", telemetry_episode(status.episode));
    return result.take();
}

static std::string
telemetry_population(const BObolLodConvergenceStatus &status)
{
    TelemetryJsonObject result;
    result.integer("expected_leaves", status.expectedLeafCount);
    result.integer("available_leaves", status.availableLeafCount);
    result.integer("visible_targets", status.visibleTargetCount);
    result.integer("active_payloads", status.activePayloadCount);
    result.integer("satisfied_payloads", status.satisfiedPayloadCount);
    result.integer("presented_subpixel_occurrences",
	status.presentedSubpixelOccurrenceCount);
    result.integer("presented_structural_boxes",
	status.presentedStructuralBoxCount);
    result.integer("terminal_proxy_occurrences",
	status.terminalProxyOccurrenceCount);
    result.integer("terminal_failure_occurrences",
	status.terminalOccurrenceFailureCount);
    result.integer("active_faces", status.activeFaces);
    result.integer("active_source_faces", status.activeSourceFaces);
    result.integer("source_mesh_occurrences",
	status.sourceMeshOccurrenceCount);
    result.integer("temporary_coverage_occurrences",
	status.temporaryCoverageOccurrenceCount);
    result.boolean("presented_primitive_count_valid",
	status.presentedPrimitiveCountValid != FALSE);
    result.integer("presented_primitives", status.presentedPrimitiveCount);
    result.raw("proxy_reasons", telemetry_proxy_reasons(status.proxyReasons));
    return result.take();
}

static std::string
telemetry_capacity_search(const BObolLodConvergenceStatus &status)
{
    TelemetryJsonObject result;
    result.integer("phase", status.capacitySearchPhase);
    result.integer("goal", status.capacitySearchGoal);
    result.integer("samples_remaining",
	status.capacitySearchSamplesRemaining);
    result.integer("measured_candidates",
	status.capacitySearchMeasuredCandidates);
    result.integer("total_measured_candidates",
	status.capacitySearchTotalMeasuredCandidates);
    result.integer("candidate_limit", status.capacitySearchCandidateLimit);
    result.integer("maximum_candidates",
	status.capacitySearchMaximumCandidates);
    result.integer("sample_limit", status.capacitySearchSampleLimit);
    result.integer("invalid_frame_attempts",
	status.capacitySearchInvalidFrameAttempts);
    result.integer("invalid_frame_attempt_limit",
	status.capacitySearchInvalidFrameAttemptLimit);
    result.integer("completed_units", status.capacitySearchCompletedUnits);
    result.integer("total_units", status.capacitySearchTotalUnits);
    return result.take();
}

static std::string
telemetry_allocation(const BObolLodConvergenceStatus &status)
{
    TelemetryJsonObject result;
    result.integer("current_plan_serial", status.currentAllocationPlanSerial);
    result.integer("committed_plan_serial",
	status.committedAllocationPlanSerial);
    result.integer("active_render_cost", status.activeRenderCost);
    result.integer("retained_render_cost", status.retainedRenderCost);
    result.integer("allocation_managed_render_cost",
	status.allocationManagedRenderCost);
    result.integer("allocation_unmanaged_render_cost",
	status.allocationUnmanagedRenderCost);
    result.boolean("presented_render_cost_valid",
	status.presentedRenderCostValid != FALSE);
    result.integer("presented_render_cost", status.presentedRenderCost);
    result.integer("render_cost_budget", status.renderCostBudget);
    result.real("interactive_calibrated_render_cost_per_second",
	status.interactiveCalibratedRenderCostPerSecond);
    result.integer("interactive_entry_render_cost_budget",
	status.interactiveEntryRenderCostBudget);
    result.integer("interactive_entry_progressive_cut_ceiling",
	status.interactiveEntryProgressiveCutCeiling);
    result.integer("selected_presentation_cost",
	status.selectedPresentationCost);
    result.integer("certified_presentation_budget",
	status.certifiedPresentationBudget);
    result.integer("pixel_demand_presentation_cost",
	status.pixelDemandPresentationCost);
    result.integer("requested_presentation_budget",
	status.requestedPresentationBudget);
    result.integer("maximum_marginal_presentation_budget",
	status.maximumMarginalPresentationBudget);
    result.integer("maximum_protected_presentation_budget",
	status.maximumProtectedPresentationBudget);
    result.integer("external_presentation_cost",
	status.allocationExternalPresentationCost);
    result.integer("resident_admission_revision",
	status.allocationResidentAdmissionRevision);
    result.real("point_proxy_pixel_threshold",
	status.allocationPointProxyPixelThreshold);
    result.boolean("protected_floor_allowed",
	status.allocationProtectedFloorAllowed != FALSE);
    result.integer("point_proxy_candidates", status.pointProxyCandidateCount);
    result.integer("reachable_point_proxy_candidates",
	status.reachablePointProxyCandidateCount);
    result.integer("selected_point_proxies", status.selectedPointProxyCount);
    result.integer("progressive_cut_ceiling", status.progressiveCutCeiling);
    result.integer("maximum_active_progressive_cut",
	status.maximumActiveProgressiveCut);
    result.boolean("maximum_nonaggregated_progressive_cut_known",
	status.maximumNonAggregatedProgressiveCutKnown != FALSE);
    result.integer("maximum_nonaggregated_progressive_cut",
	status.maximumNonAggregatedProgressiveCut);
    result.boolean("certificate_current",
	status.allocationCertificateCurrent != FALSE);
    result.boolean("cuts_applied", status.allocationCutsApplied != FALSE);
    result.integer("active_mismatch_count",
	status.activeAllocationMismatchCount);
    result.boolean("presentation_realized",
	status.allocationPresentationRealized != FALSE);
    result.boolean("frame_exact", status.presentationFrameExact != FALSE);
    result.boolean("point_proxy_protection_classified",
	status.pointProxyProtectionClassified != FALSE);
    result.integer("prominent_candidates", status.prominentCandidateCount);
    result.integer("prominent_floor_violations",
	status.prominentQualityFloorViolationCount);
    result.real("maximum_normalized_visual_error",
	status.maximumNormalizedVisualError);
    result.real("visual_importance_debt", status.visualImportanceDebt);
    return result.take();
}

static std::string
telemetry_work(const BObolLodConvergenceStatus &status)
{
    TelemetryJsonObject result;
    result.integer("active_generation", status.activeGeneration);
    result.integer("submission_source_index", status.submissionSourceIndex);
    result.integer("submission_entry_offset", status.submissionEntryOffset);
    result.integer("pending_tasks", status.pendingTasks);
    result.integer("in_flight", status.inFlight);
    result.integer("queued_results", status.queuedResults);
    result.integer("queued_cache_writes", status.queuedCacheWrites);
    result.integer("shared_producer_leases", status.sharedProducerLeases);
    result.integer("runnable_queued_tasks", status.runnableQueuedTasks);
    result.integer("dependency_blocked_tasks",
	status.dependencyBlockedTasks);
    result.integer("transient_memory_blocked_tasks",
	status.transientMemoryBlockedTasks);
    result.integer("cpu_admission_waiting_tasks",
	status.cpuAdmissionWaitingTasks);
    result.boolean("task_submission_capacity_blocked",
	status.taskSubmissionCapacityBlocked != FALSE);
    result.boolean("result_submission_capacity_blocked",
	status.resultSubmissionCapacityBlocked != FALSE);
    result.integer("producer_stage_mask", status.producerStageMask);
    result.integer("producer_stage", status.producerStage);
    result.text("producer_stage_name", telemetry_producer_stage_name(
	status.producerStage));
    result.raw("producer_stage_task_counts", telemetry_size_array(
	status.producerStageTaskCounts, BOBOL_LOD_PRODUCER_STAGE_COUNT));
    result.integer("producer_stage_task_count",
	status.producerStageTaskCount);
    result.integer("active_producer_count", status.activeProducerCount);
    result.integer("producer_stage_completed_units",
	status.producerStageCompletedUnits);
    result.integer("producer_stage_total_units",
	status.producerStageTotalUnits);
    result.integer("oldest_pending_task_age_us",
	status.oldestPendingTaskAgeMicroseconds);
    result.integer("maximum_producer_queue_wait_us",
	status.maximumProducerQueueWaitMicroseconds);
    result.integer("maximum_producer_elapsed_us",
	status.maximumProducerElapsedMicroseconds);
    result.integer("producer_stage_elapsed_us",
	status.producerStageElapsedMicroseconds);
    result.integer("active_producer_source_faces",
	status.activeProducerSourceFaceCount);
    result.integer("active_producer_source_vertices",
	status.activeProducerSourcePointCount);
    result.integer("active_producer_source_bytes",
	status.activeProducerSourceByteCount);
    return result.take();
}

static std::string
telemetry_renderer_preparation(const BObolLodConvergenceStatus &status)
{
    TelemetryJsonObject result;
    result.integer("target_signature",
	status.rendererPreparationTargetSignature);
    result.integer("total_units", status.rendererPreparationTotalUnits);
    result.integer("completed_units",
	status.rendererPreparationCompletedUnits);
    result.integer("remaining_units",
	status.rendererPreparationRemainingUnits);
    result.integer("reserved_bytes",
	status.rendererPreparationReservedBytes);
    result.integer("target_count", status.rendererPreparationTargetCount);
    result.integer("preparing_target_count",
	status.rendererPreparationPreparingTargetCount);
    result.integer("constrained_target_count",
	status.rendererPreparationConstrainedTargetCount);
    result.integer("failed_target_count",
	status.rendererPreparationFailedTargetCount);
    result.integer("invalid_target_count",
	status.rendererPreparationInvalidTargetCount);
    return result.take();
}

static std::string
telemetry_memory(const BObolLodConvergenceStatus &status)
{
    TelemetryJsonObject result;
    result.integer("resident_mesh_bytes", status.residentMeshBytes);
    result.integer("stable_resident_mesh_bytes",
	status.stableResidentMeshBytes);
    result.integer("reserved_resident_growth_bytes",
	status.reservedResidentMeshGrowthBytes);
    result.integer("resident_mesh_limit_bytes", status.residentMeshLimitBytes);
    result.integer("memory_limited_payload_count",
	status.memoryLimitedPayloadCount);
    result.integer("active_working_set_bytes", status.activeWorkingSetBytes);
    result.integer("peak_working_set_bytes", status.peakWorkingSetBytes);
    result.integer("resident_compaction_count",
	status.residentCompactionCount);
    result.integer("resident_compaction_plan_revision",
	status.residentCompactionPlanRevision);
    result.integer("resident_compaction_candidate_count",
	status.residentCompactionCandidateCount);
    result.boolean("resident_compaction_plan_current",
	status.residentCompactionPlanCurrent != FALSE);
    return result.take();
}

static std::string
telemetry_gpu(const BObolLodConvergenceStatus &status)
{
    TelemetryJsonObject result;
    result.integer("tracked_buffer_bytes", status.gpuTrackedBufferBytes);
    result.integer("ordinary_part_buffer_bytes",
	status.gpuOrdinaryPartBufferBytes);
    result.integer("progressive_cut_buffer_bytes",
	status.gpuProgressiveCutBufferBytes);
    result.integer("progressive_active_cut_bytes",
	status.gpuProgressiveActiveCutBytes);
    result.integer("batch_buffer_bytes", status.gpuBatchBufferBytes);
    result.integer("triangle_atlas_allocated_bytes",
	status.gpuTriangleAtlasAllocatedBytes);
    result.integer("triangle_atlas_live_bytes",
	status.gpuTriangleAtlasLiveBytes);
    result.integer("triangle_atlas_capacity_bytes",
	status.gpuTriangleAtlasConfiguredCapacityBytes);
    result.integer("triangle_atlas_part_count",
	status.gpuTriangleAtlasPartCount);
    result.integer("triangle_atlas_page_count",
	status.gpuTriangleAtlasPageCount);
    result.integer("pressure_proxy_count", status.gpuPressureProxyCount);
    result.integer("progressive_eviction_count",
	status.gpuProgressiveEvictionCount);
    result.integer("triangle_atlas_reclamation_count",
	status.gpuTriangleAtlasReclamationCount);
    result.integer("resource_sample_serial", status.gpuResourceSampleSerial);
    return result.take();
}

static std::string
telemetry_pending(const BObolLodControlTraceState &state)
{
    const BObolLodConvergenceStatus &status = state.convergence;
    TelemetryJsonObject result;
    result.boolean("refinement_frame", status.refinementFramePending != FALSE);
    result.boolean("budget_calibration",
	status.budgetCalibrationPending != FALSE);
    result.boolean("stable_presentation_handoff",
	status.stablePresentationHandoffPending != FALSE);
    result.boolean("point_proxy_calibration",
	status.pointProxyCalibrationPending != FALSE);
    result.boolean("point_proxy_admission_frame",
	status.pointProxyAdmissionFramePending != FALSE);
    result.boolean("stable_point_proxy_calibration",
	status.stablePointProxyCalibrationPending != FALSE);
    result.boolean("point_proxy_triangle_recovery",
	status.pointProxyTriangleRecoveryPending != FALSE);
    result.boolean("resident_growth_reallocation",
	status.residentGrowthReallocationPending != FALSE);
    result.boolean("publication_frame",
	status.publicationFramePending != FALSE);
    result.boolean("semantic_presentation_frame",
	status.semanticPresentationFramePending != FALSE);
    result.boolean("source_preparation",
	status.sourcePreparationPending != FALSE);
    result.integer("source_preparation_provider_count",
	status.sourcePreparationProviderCount);
    result.integer("source_preparation_completed_units",
	status.sourcePreparationCompletedUnits);
    result.integer("source_preparation_total_units",
	status.sourcePreparationTotalUnits);
    result.integer("refinement_cooldown_remaining_us",
	state.refinementCooldownRemainingMicroseconds);
    return result.take();
}

static std::string
telemetry_submission(const BObolLodControlTraceState &state)
{
    TelemetryJsonObject result;
    result.boolean("delta_active", state.submissionDeltaTargetCount != 0);
    result.integer("delta_target_count", state.submissionDeltaTargetCount);
    result.integer("delta_plan_count", state.submissionDeltaPlanCount);
    result.integer("delta_entry_count", state.submissionDeltaEntryCount);
    result.boolean("structural_repair_active",
	state.structuralRepairFrontierCount != 0);
    result.boolean("structural_terminal_proxy",
	state.structuralRepairTerminalProxy != FALSE);
    result.boolean("structural_point_relaxation",
	state.structuralRepairPointRelaxation != FALSE);
    result.integer("structural_frontier_count",
	state.structuralRepairFrontierCount);
    result.integer("structural_coverage_cost_reservation",
	state.structuralRepairCoverageCostReservation);
    result.boolean("pass_admitted_work",
	state.submissionPassAdmittedWork != FALSE);
    result.boolean("pass_cut_advanced",
	state.submissionPassCutAdvanced != FALSE);
    result.boolean("pass_refinement_pending",
	state.submissionPassRefinementPending != FALSE);
    result.boolean("pass_residency_pending",
	state.submissionPassResidencyPending != FALSE);
    result.boolean("pass_budget_blocked",
	state.submissionPassBudgetBlocked != FALSE);
    result.integer("pass_missing_mesh_budget_blocked",
	state.submissionPassMissingMeshBudgetBlockedCount);
    result.integer("last_visited_meshes",
	state.lastSubmissionVisitedMeshCount);
    result.integer("last_submitted_tasks",
	state.lastSubmissionSubmittedTaskCount);
    result.integer("last_updated_cuts",
	state.lastSubmissionUpdatedCutCount);
    result.integer("last_skipped_meshes",
	state.lastSubmissionSkippedMeshCount);
    return result.take();
}

static std::string
telemetry_signals(const BObolLodControlTraceState &state)
{
    const BObolLodConvergenceStatus &status = state.convergence;
    const bool producerIdle = status.pendingTasks == 0 &&
	status.inFlight == 0 && status.queuedResults == 0 &&
	status.sharedProducerLeases == 0 && status.activeProducerCount == 0;
    const bool rendererPrepared =
	status.rendererPreparationRemainingUnits == 0 &&
	status.rendererPreparationPreparingTargetCount == 0;
    const bool activeOverBudget = status.renderCostBudget > 0 &&
	status.activeRenderCost > status.renderCostBudget;
    const bool handoffWithoutActiveRenderOrProducer =
	status.stablePresentationHandoffPending != FALSE && producerIdle &&
	rendererPrepared && state.hostWork.renderPending() == FALSE &&
	state.hostWork.frameClaimed() == FALSE;
    const bool handoffWaitingForCooldown =
	handoffWithoutActiveRenderOrProducer &&
	state.refinementCooldownRemainingMicroseconds > 0;
    const bool handoffWithoutProgressRoute =
	handoffWithoutActiveRenderOrProducer &&
	state.hostWork.pumpPending() == FALSE && !handoffWaitingForCooldown;
    TelemetryJsonObject result;
    result.boolean("producer_idle", producerIdle);
    result.boolean("renderer_prepared", rendererPrepared);
    result.boolean("all_visible_targets_satisfied",
	status.visibleTargetCount > 0 &&
	status.satisfiedPayloadCount >= status.visibleTargetCount);
    result.boolean("all_visible_targets_have_source_mesh",
	status.visibleTargetCount > 0 &&
	status.sourceMeshOccurrenceCount >= status.visibleTargetCount);
    result.boolean("structural_boxes_present",
	status.presentedStructuralBoxCount > 0);
    result.boolean("active_cost_exceeds_budget", activeOverBudget);
    result.integer("active_cost_excess", activeOverBudget ?
	status.activeRenderCost - status.renderCostBudget : 0);
    result.boolean("allocation_presentation_cost_mismatch",
	status.currentAllocationPlanSerial != 0 &&
	status.activeRenderCost != status.selectedPresentationCost);
    result.boolean("handoff_waiting_for_cooldown",
	handoffWaitingForCooldown);
    result.boolean("handoff_without_render_or_producer",
	handoffWithoutActiveRenderOrProducer);
    result.boolean("handoff_without_progress_route",
	handoffWithoutProgressRoute);
    return result.take();
}

static std::string
telemetry_state(const BObolLodControlTraceState &state, bool sanitize)
{
    const BObolLodConvergenceStatus &status = state.convergence;
    TelemetryJsonObject result;
    result.integer("controller_view_revision", state.viewRevision);
    result.integer("controller_policy_revision", state.policyRevision);
    result.integer("render_completion_serial", state.renderCompletionSerial);
    result.boolean("interaction_active", state.interactionActive != FALSE);
    if (!sanitize)
	result.text("render_reason", state.renderReason.getString());
    else if (state.renderReason.getLength() > 0)
	result.boolean("render_reason_redacted", true);
    result.raw("host_work", telemetry_host_work(state.hostWork));
    result.raw("control", telemetry_control(status));
    result.raw("revisions", telemetry_revisions(status));
    result.raw("progress", telemetry_progress(status));
    result.raw("population", telemetry_population(status));
    result.raw("capacity_search", telemetry_capacity_search(status));
    result.raw("allocation", telemetry_allocation(status));
    result.raw("work", telemetry_work(status));
    result.raw("renderer_preparation",
	telemetry_renderer_preparation(status));
    result.raw("memory", telemetry_memory(status));
    result.raw("gpu", telemetry_gpu(status));
    result.raw("pending", telemetry_pending(state));
    result.raw("submission", telemetry_submission(state));
    result.raw("signals", telemetry_signals(state));
    result.integer("presentation_transaction_serial",
	status.presentationTransactionSerial);
    result.integer("presentation_required_render_serial",
	status.presentationRequiredRenderSerial);
    result.integer("presented_frame_serial", status.presentedFrameSerial);
    return result.take();
}

static std::string
telemetry_camera(const BObolViewController &controller)
{
    TelemetryJsonObject result;
    const SoCamera *camera = controller.getCamera();
    if (!camera) {
	result.text("type", "none");
	return result.take();
    }

    const SbVec3f position = camera->position.getValue();
    const float *rotation = camera->orientation.getValue().getValue();
    result.text("type", camera->isOfType(
	SoOrthographicCamera::getClassTypeId()) ? "orthographic" :
	camera->isOfType(SoPerspectiveCamera::getClassTypeId()) ?
	    "perspective" : "camera");
    result.raw("position", telemetry_vec3(position));
    std::string quaternion("[");
    for (int i = 0; i < 4; ++i) {
	if (i)
	    quaternion.push_back(',');
	telemetry_append_json_real(quaternion, rotation[i],
	    std::numeric_limits<float>::max_digits10);
    }
    quaternion.push_back(']');
    result.raw("orientation", quaternion);
    result.real("near_distance", camera->nearDistance.getValue());
    result.real("far_distance", camera->farDistance.getValue());
    result.real("focal_distance", camera->focalDistance.getValue());
    if (camera->isOfType(SoOrthographicCamera::getClassTypeId()))
	result.real("height", static_cast<const SoOrthographicCamera *>(
	    camera)->height.getValue());
    if (camera->isOfType(SoPerspectiveCamera::getClassTypeId()))
	result.real("height_angle", static_cast<const SoPerspectiveCamera *>(
	    camera)->heightAngle.getValue());
    const SbVec2s viewport =
	controller.getViewportRegion().getViewportSizePixels();
    std::string viewportJson = "[" + std::to_string(viewport[0]) + "," +
	std::to_string(viewport[1]) + "]";
    result.raw("viewport_pixels", viewportJson);
    return result.take();
}

class BObolLodTelemetrySink {
public:
    BObolLodTelemetrySink(const std::string &path, bool sanitized) :
	pathValue(path),
	sanitizedValue(sanitized),
	started(std::chrono::steady_clock::now())
    {
	this->file = std::fopen(path.c_str(), "wb");
	if (!this->file) {
	    bu_log("Unable to open BObol LoD telemetry file: %s\n",
		path.c_str());
	    return;
	}
	TelemetryJsonObject schema;
	schema.integer("version", bobol_lod_telemetry_schema_version);
	schema.text("format", "jsonl");
	schema.boolean("sanitized", sanitized);
	this->write(0, "schema", schema.take());
    }

    ~BObolLodTelemetrySink()
    {
	if (this->file)
	    std::fclose(this->file);
    }

    bool valid() const
    {
	return this->file != NULL &&
	    !this->failed.load(std::memory_order_relaxed);
    }

    bool sanitized() const
    {
	return this->sanitizedValue;
    }

    uint64_t allocateController()
    {
	std::lock_guard<std::mutex> guard(this->mutex);
	if (!this->file && !this->failed.load(std::memory_order_relaxed))
	    this->file = std::fopen(this->pathValue.c_str(), "ab");
	if (!this->file) {
	    this->failed.store(true, std::memory_order_relaxed);
	    return 0;
	}
	++this->activeControllers;
	return this->nextController++;
    }

    void releaseController()
    {
	std::lock_guard<std::mutex> guard(this->mutex);
	if (this->activeControllers)
	    --this->activeControllers;
	if (!this->activeControllers && this->file) {
	    std::fclose(this->file);
	    this->file = NULL;
	}
    }

    void write(uint64_t controller, const char *kind,
	const std::string &data) noexcept
    {
	try {
	    std::lock_guard<std::mutex> guard(this->mutex);
	    if (!this->file ||
		this->failed.load(std::memory_order_relaxed))
		return;
	    const auto elapsed = std::chrono::duration_cast<
		std::chrono::microseconds>(
		    std::chrono::steady_clock::now() - this->started).count();
	    const auto wall = std::chrono::duration_cast<
		std::chrono::microseconds>(
		    std::chrono::system_clock::now().time_since_epoch()).count();
	    TelemetryJsonObject record;
	    record.text("record", kind);
	    record.integer("schema_version", bobol_lod_telemetry_schema_version);
	    record.integer("sequence", this->nextSequence++);
	    record.integer("elapsed_us", elapsed);
	    record.integer("wall_time_unix_us", wall);
	    record.integer("controller", controller);
	    record.raw("data", data);
	    std::string line = record.take();
	    line.push_back('\n');
	    if (std::fwrite(line.data(), 1, line.size(), this->file) !=
		line.size() || std::fflush(this->file) != 0) {
		this->failed.store(true, std::memory_order_relaxed);
		bu_log("BObol LoD telemetry output failed; logging disabled\n");
	    }
	} catch (...) {
	    this->failed.store(true, std::memory_order_relaxed);
	}
    }

private:
    std::string pathValue;
    FILE *file = NULL;
    bool sanitizedValue = false;
    std::atomic_bool failed {false};
    std::mutex mutex;
    std::chrono::steady_clock::time_point started;
    uint64_t nextSequence = 1;
    uint64_t nextController = 1;
    size_t activeControllers = 0;
};

static std::shared_ptr<BObolLodTelemetrySink>
telemetry_acquire_sink(const std::string &path, bool sanitized)
{
    struct Registry {
	std::mutex mutex;
	std::unordered_map<std::string,
	    std::shared_ptr<BObolLodTelemetrySink>> sinks;
    };
    static Registry registry;

    std::lock_guard<std::mutex> guard(registry.mutex);
    const auto found = registry.sinks.find(path);
    if (found != registry.sinks.end()) {
	std::shared_ptr<BObolLodTelemetrySink> sink = found->second;
	if (sink->sanitized() != sanitized)
	    bu_log("BObol LoD telemetry sanitization setting changed after "
		   "the log was opened; retaining the original setting\n");
	return sink;
    }

    std::shared_ptr<BObolLodTelemetrySink> sink =
	std::make_shared<BObolLodTelemetrySink>(path, sanitized);
    if (!sink->valid())
	return std::shared_ptr<BObolLodTelemetrySink>();
    registry.sinks[path] = sink;
    return sink;
}

class TelemetryAliases {
public:
    std::string database(const struct db_i *dbip)
    {
	const auto found = this->databases.find(dbip);
	if (found != this->databases.end())
	    return found->second;
	const std::string alias = "db" +
	    std::to_string(this->nextDatabase++) + ".g";
	this->databases.emplace(dbip, alias);
	return alias;
    }

    std::string path(const struct db_i *dbip, const char *raw,
	bool finalCombination)
    {
	const std::string input = raw ? raw : "";
	std::string result;
	size_t begin = 0;
	while (begin < input.size()) {
	    while (begin < input.size() &&
		(input[begin] == '/' || input[begin] == '\\'))
		++begin;
	    if (begin >= input.size())
		break;
	    size_t end = input.find_first_of("/\\", begin);
	    if (end == std::string::npos)
		end = input.size();
	    const std::string component = input.substr(begin, end - begin);
	    size_t next = end;
	    while (next < input.size() &&
		(input[next] == '/' || input[next] == '\\'))
		++next;
	    bool combination = next < input.size() || finalCombination;
	    if (dbip) {
		struct directory *directory = db_lookup(
		    const_cast<struct db_i *>(dbip), component.c_str(),
		    LOOKUP_QUIET);
		if (directory != RT_DIR_NULL)
		    combination = (directory->d_flags & RT_DIR_COMB) != 0;
	    }
	    if (!result.empty())
		result.push_back('/');
	    result += combination ? this->combination(dbip, component) :
		this->solid(dbip, component);
	    begin = next;
	}
	return result;
    }

    std::string solid(const struct db_i *dbip, const std::string &name)
    {
	return this->objectAlias(dbip, name, false);
    }

private:
    std::string objectAlias(const struct db_i *dbip,
	const std::string &name, bool combination)
    {
	std::string key = std::to_string(
	    reinterpret_cast<uintptr_t>(dbip));
	key.push_back('\n');
	key += name;
	auto &objects = combination ? this->combinations : this->solids;
	const auto found = objects.find(key);
	if (found != objects.end())
	    return found->second;
	const std::string alias = combination ?
	    "c" + std::to_string(this->nextCombination++) + ".c" :
	    "s" + std::to_string(this->nextSolid++) + ".s";
	objects.emplace(std::move(key), alias);
	return alias;
    }

    std::string combination(const struct db_i *dbip,
	const std::string &name)
    {
	return this->objectAlias(dbip, name, true);
    }

    std::unordered_map<const struct db_i *, std::string> databases;
    std::unordered_map<std::string, std::string> combinations;
    std::unordered_map<std::string, std::string> solids;
    uint64_t nextDatabase = 1;
    uint64_t nextCombination = 1;
    uint64_t nextSolid = 1;
};

static std::string
telemetry_source_key(const SoBRLDatabaseSource *source)
{
    const uint64_t routing = source ? source->getCompactSourceRoutingId() : 0;
    if (routing)
	return "routing:" + std::to_string(routing);
    return "pointer:" + std::to_string(
	reinterpret_cast<uintptr_t>(source));
}

static const char *
telemetry_database_path(const struct db_i *dbip)
{
    return dbip && dbip->dbi_filename ? dbip->dbi_filename : "";
}

static const char *
telemetry_sanitized_object_type(const char *label)
{
    static const char *allowed[] = {
	"annot", "annotation", "arb8", "arbn", "ars", "binunif", "bot",
	"brep", "cline", "dsp", "ehy", "ell", "epa", "eto", "extrude",
	"grip", "half", "heart", "hf", "hyp", "metaball", "nmg", "part",
	"pipe", "pnts", "poly", "rec", "rcc", "rhc", "rpc", "sketch",
	"sph", "superell", "tgc", "tor", "vol", "wire"
    };
    for (const char *candidate : allowed) {
	if (label && BU_STR_EQUAL(label, candidate))
	    return candidate;
    }
    return "unknown";
}

static const char *
telemetry_sanitized_geometry_kind(const char *label)
{
    static const char *allowed[] = {
	"aabb", "annotation", "line", "mesh", "obb", "point", "proxy",
	"surface", "wire"
    };
    for (const char *candidate : allowed) {
	if (label && BU_STR_EQUAL(label, candidate))
	    return candidate;
    }
    return "unknown";
}

} /* namespace */

class BObolLodTelemetry::Impl {
public:
    struct SourceCursor {
	uint64_t id = 0;
	uint64_t inventoryRevision = 0;
	uint64_t visibilityRevision = 0;
	uint64_t populationEpoch = 0;
    };

    explicit Impl(std::shared_ptr<BObolLodTelemetrySink> output) :
	sink(std::move(output))
    {
	this->controllerId = this->sink->allocateController();
	if (!this->controllerId)
	    throw std::runtime_error("telemetry controller allocation failed");
    }

    ~Impl()
    {
	if (this->started) {
	    TelemetryJsonObject data;
	    data.integer("transition_count", this->transitionSerial);
	    this->sink->write(this->controllerId, "controller_end", data.take());
	}
	this->sink->releaseController();
    }

    void recordTransition(const BObolViewController &controller,
	const BObolLodControlTransitionRecord &transition)
    {
	if (!this->started) {
	    TelemetryJsonObject start;
	    start.boolean("sanitized", this->sink->sanitized());
	    start.raw("camera", telemetry_camera(controller));
	    this->sink->write(this->controllerId, "controller_start",
		start.take());
	    this->started = true;
	}

	TelemetryJsonObject data;
	data.integer("transition_serial", ++this->transitionSerial);
	data.text("event",
	    bobol_lod_control_transition_event_name(transition.event));
	data.raw("before", telemetry_state(transition.before,
	    this->sink->sanitized()));
	data.raw("after", telemetry_state(transition.after,
	    this->sink->sanitized()));
	data.raw("camera", telemetry_camera(controller));
	this->sink->write(this->controllerId, "transition", data.take());

	const uint64_t inventoryRevision =
	    transition.after.convergence.inventoryRevision;
	const bool inspectInventory = !this->inventoryObserved ||
	    inventoryRevision != this->lastInventoryRevision ||
	    transition.event == BOBOL_LOD_CONTROL_TRANSITION_INITIAL ||
	    transition.event == BOBOL_LOD_CONTROL_TRANSITION_INVENTORY ||
	    transition.event == BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT;
	if (inspectInventory)
	    this->recordInventory(controller, inventoryRevision);
    }

    void recordPublication(
	const BObolLodPublicationTelemetryRecord &publication)
    {
	TelemetryJsonObject data;
	data.integer("after_transition_serial", this->transitionSerial);
	data.integer("processed", publication.processed);
	data.integer("matched", publication.matched);
	data.integer("applied", publication.applied);
	data.integer("rejected", publication.rejected);
	data.integer("unmatched", publication.unmatched);
	data.raw("provider_status",
	    telemetry_provider_status_counts(publication));

	TelemetryJsonObject authentication;
	authentication.integer("publish",
	    publication.authenticationPublishCount);
	authentication.integer("terminal_failure",
	    publication.authenticationTerminalFailureCount);
	authentication.integer("retry_current_demand",
	    publication.authenticationRetryCount);
	authentication.integer("supersede",
	    publication.authenticationSupersedeCount);
	data.raw("authentication_disposition", authentication.take());

	TelemetryJsonObject mismatches;
	mismatches.integer("source_route",
	    publication.sourceRouteMismatchCount);
	mismatches.integer("source_population",
	    publication.sourcePopulationMismatchCount);
	mismatches.integer("demand", publication.demandMismatchCount);
	data.raw("authentication_mismatch", mismatches.take());

	TelemetryJsonObject source;
	source.integer("accepted", publication.sourceAcceptedCount);
	source.integer("unchanged", publication.sourceUnchangedCount);
	source.integer("rejected", publication.sourceRejectedCount);
	source.integer("retry_current_demand",
	    publication.sourceRetryCount);
	source.integer("update_retry_current_demand",
	    publication.updateRetryCount);
	data.raw("publication_disposition", source.take());

	data.integer("retry_source_entry_count",
	    publication.retrySourceEntryCount);
	data.boolean("retry_source_entries_truncated",
	    publication.retrySourceEntryCount >
		publication.retrySourceEntryIndices.size());
	data.raw("retry_source_entry_indices",
	    telemetry_u32_array(publication.retrySourceEntryIndices));
	data.integer("result_sample_count", publication.resultSampleCount);
	data.boolean("result_samples_truncated",
	    publication.resultSampleCount > publication.resultSamples.size());
	data.raw("result_samples",
	    telemetry_result_samples(publication.resultSamples));

	TelemetryJsonObject replay;
	replay.boolean("authentication",
	    publication.replayFromAuthentication);
	replay.boolean("source", publication.replayFromSource);
	replay.boolean("update_action", publication.replayFromUpdate);
	replay.boolean("retained_publication",
	    publication.replayFromRetainedPublication);
	replay.boolean("partial_refinement",
	    publication.replayFromPartialRefinement);
	replay.boolean("actionable_quality_debt",
	    publication.actionableQualityDebt);
	replay.boolean("requested", publication.replayRequested);
	data.raw("replay", replay.take());

	TelemetryJsonObject cursor;
	cursor.boolean("active_before",
	    publication.submissionActiveBeforeReplay);
	cursor.integer("source_index_before",
	    publication.submissionSourceIndexBeforeReplay);
	cursor.integer("entry_offset_before",
	    publication.submissionEntryOffsetBeforeReplay);
	cursor.boolean("active_after",
	    publication.submissionActiveAfterReplay);
	cursor.boolean("rescan_pending_after",
	    publication.rescanPendingAfterReplay);
	cursor.integer("source_index_after",
	    publication.submissionSourceIndexAfterReplay);
	cursor.integer("entry_offset_after",
	    publication.submissionEntryOffsetAfterReplay);
	data.raw("submission_cursor", cursor.take());

	this->sink->write(this->controllerId, "publication_results",
	    data.take());
    }

    void recordGap()
    {
	TelemetryJsonObject data;
	data.integer("after_transition_serial", this->transitionSerial);
	this->sink->write(this->controllerId, "transition_gap", data.take());
    }

private:
    void recordInventoryObject(const SoBRLDatabaseSource &source,
	uint64_t sourceId, size_t entryIndex,
	const BObolCompactOccurrence &occurrence)
    {
	const struct db_i *database = source.getDatabase();
	const BObolRealizedShapeSummary &summary = occurrence.summary;
	const BObolSourceMeshRequest *request =
	    occurrence.sourceMeshRequestValid ?
		&occurrence.sourceMeshRequest : NULL;
	const char *sourceType = request &&
	    request->sourceType.getLength() > 0 ?
		request->sourceType.getString() : summary.sourceType.getString();
	const uint64_t faces = request && request->faceCount ?
	    request->faceCount : summary.lodFaceCount ? summary.lodFaceCount :
		static_cast<uint64_t>(std::max(0, summary.triangleCount));
	const uint64_t vertices = request && request->pointCount ?
	    request->pointCount : summary.lodPointCount ? summary.lodPointCount :
		static_cast<uint64_t>(std::max(0, summary.pointCount));

	TelemetryJsonObject data;
	data.integer("source_id", sourceId);
	data.integer("entry_index", entryIndex);
	data.text("operation", "upsert");
	data.text("object_kind", "solid");
	if (this->sink->sanitized()) {
	    data.text("path", this->aliases.path(database,
		summary.path.getString(), false));
	    data.text("object", this->aliases.solid(database,
		summary.sourceName.getString()));
	    if (request) {
		data.text("mesh_asset_path", this->aliases.path(database,
		    request->meshAssetPath.getString(), false));
		data.text("mesh_asset", this->aliases.solid(database,
		    request->meshAssetName.getString()));
	    }
	} else {
	    data.text("path", summary.path.getString());
	    data.text("object", summary.sourceName.getString());
	    if (request) {
		data.text("mesh_asset_path",
		    request->meshAssetPath.getString());
		data.text("mesh_asset", request->meshAssetName.getString());
		data.integer("mesh_asset_content_hash",
		    request->meshAssetContentHash);
	    }
	}
	data.text("object_type", this->sink->sanitized() ?
	    telemetry_sanitized_object_type(sourceType) :
	    (sourceType && sourceType[0] ? sourceType : "unknown"));
	data.text("geometry_kind", this->sink->sanitized() ?
	    telemetry_sanitized_geometry_kind(summary.geometryKind.getString()) :
	    summary.geometryKind.getString());
	data.integer("faces", faces);
	data.integer("vertices", vertices);
	data.boolean("lod_backed", occurrence.lodBacked != FALSE);
	data.boolean("source_mesh_request",
	    occurrence.sourceMeshRequestValid != FALSE);
	data.boolean("visible", summary.visible != FALSE);
	this->sink->write(this->controllerId, "inventory_object", data.take());
    }

    void recordInventorySource(const SoBRLDatabaseSource &source,
	const SourceCursor &cursor, bool reset)
    {
	const struct db_i *database = source.getDatabase();
	TelemetryJsonObject data;
	data.integer("source_id", cursor.id);
	data.boolean("inventory_reset", reset);
	if (this->sink->sanitized()) {
	    data.text("database", this->aliases.database(database));
	    data.text("path", this->aliases.path(database,
		source.path.getValue().getString(), true));
	} else {
	    data.text("database", telemetry_database_path(database));
	    data.text("path", source.path.getValue().getString());
	    data.text("instance_key", source.instanceKey.getValue().getString());
	}
	data.text("object_kind", "combination");
	data.integer("inventory_revision", cursor.inventoryRevision);
	data.integer("population_epoch", cursor.populationEpoch);
	data.integer("occurrence_count", source.getCompactInstanceCount());
	data.integer("expected_occurrence_count",
	    source.getCompactExpectedInstanceCount());
	data.boolean("compact_population_complete",
	    source.hasCompleteCompactInstancePopulation() != FALSE);
	data.integer("part_count", source.getCompactPartCount());
	data.integer("selected_occurrence_count",
	    source.getCompactSelectedInstanceCount());
	data.integer("visibility_revision",
	    source.getDisplayMeshLodVisibilityRevision());
	data.integer("draw_mode", source.drawMode.getValue());
	data.integer("representation_mode", source.representationMode.getValue());
	data.boolean("visible", source.visible.getValue() != FALSE);
	BObolCompactSourceProfile profile;
	if (source.getCompactSourceProfile(profile)) {
	    TelemetryJsonObject profileData;
	    profileData.integer("occurrences", profile.occurrenceCount);
	    profileData.integer("unique_assets", profile.uniqueAssetCount);
	    profileData.integer("encoded_source_bytes",
		profile.encodedSourceBytes);
	    profileData.integer("largest_asset_bytes", profile.largestAssetBytes);
	    profileData.integer("reused_occurrences",
		profile.reusedOccurrenceCount);
	    data.raw("profile", profileData.take());
	}
	this->sink->write(this->controllerId, "inventory_source", data.take());
    }

    void recordInventory(const BObolViewController &controller,
	uint64_t inventoryRevision)
    {
	const std::vector<SoBRLDatabaseSource *> renderSources =
	    controller.getRenderDatabaseSources();
	std::unordered_set<std::string> liveSources;
	for (SoBRLDatabaseSource *source : renderSources) {
	    if (!source)
		continue;
	    const std::string key = telemetry_source_key(source);
	    liveSources.insert(key);
	    auto position = this->sources.find(key);
	    const bool newSource = position == this->sources.end();
	    if (newSource) {
		SourceCursor cursor;
		cursor.id = this->nextSourceId++;
		position = this->sources.emplace(key, cursor).first;
	    }
	    SourceCursor &cursor = position->second;
	    const uint64_t currentRevision =
		source->getDisplayMeshLodRevision();
	    const uint64_t currentVisibilityRevision =
		source->getDisplayMeshLodVisibilityRevision();
	    const uint64_t currentPopulation =
		source->getCompactPopulationEpoch();
	    if (!newSource && cursor.inventoryRevision == currentRevision &&
		cursor.visibilityRevision == currentVisibilityRevision &&
		cursor.populationEpoch == currentPopulation)
		continue;

	    std::vector<size_t> changedEntries;
	    SbBool coverageInvalidated = FALSE;
	    bool completeDelta = !newSource;
	    if (completeDelta && cursor.inventoryRevision != currentRevision) {
		std::vector<size_t> inventoryEntries;
		completeDelta = source->getDisplayMeshLodChangedEntries(
		    cursor.inventoryRevision, inventoryEntries,
		    &coverageInvalidated) != FALSE;
		changedEntries.insert(changedEntries.end(),
		    inventoryEntries.begin(), inventoryEntries.end());
	    }
	    if (completeDelta &&
		cursor.visibilityRevision != currentVisibilityRevision) {
		std::vector<size_t> visibilityEntries;
		completeDelta = source->getDisplayMeshLodVisibilityChangedEntries(
		    cursor.visibilityRevision, visibilityEntries) != FALSE;
		changedEntries.insert(changedEntries.end(),
		    visibilityEntries.begin(), visibilityEntries.end());
	    }
	    /* A coverage invalidation says that the scene-wide completeness proof
	     * changed; it does not invalidate the source-local entry identities.
	     * The delta journal still supplies every appended, replaced, or removed
	     * entry.  Treating each streamed append as a reset rewrites the growing
	     * prefix on every publication and turns telemetry into O(N^2) output. */
	    const bool reset = newSource || !completeDelta ||
		cursor.populationEpoch != currentPopulation;
	    if (reset) {
		changedEntries.clear();
		const int count = source->getCompactInstanceCount();
		changedEntries.reserve(static_cast<size_t>(std::max(0, count)));
		for (int index = 0; index < count; ++index)
		    changedEntries.push_back(static_cast<size_t>(index));
	    } else {
		std::sort(changedEntries.begin(), changedEntries.end());
		changedEntries.erase(std::unique(changedEntries.begin(),
		    changedEntries.end()), changedEntries.end());
	    }
	    cursor.inventoryRevision = currentRevision;
	    cursor.visibilityRevision = currentVisibilityRevision;
	    cursor.populationEpoch = currentPopulation;
	    this->recordInventorySource(*source, cursor, reset);
	    for (size_t entry : changedEntries) {
		if (entry > static_cast<size_t>(INT_MAX))
		    continue;
		BObolCompactOccurrence occurrence;
		if (source->getCompactOccurrence(
			static_cast<int>(entry), occurrence))
		    this->recordInventoryObject(*source, cursor.id, entry,
			occurrence);
		else {
		    TelemetryJsonObject data;
		    data.integer("source_id", cursor.id);
		    data.integer("entry_index", entry);
		    data.text("operation", "remove");
		    this->sink->write(this->controllerId, "inventory_object",
			data.take());
		}
	    }
	}

	std::vector<std::string> removed;
	for (const auto &source : this->sources) {
	    if (liveSources.find(source.first) == liveSources.end())
		removed.push_back(source.first);
	}
	for (const std::string &key : removed) {
	    const auto found = this->sources.find(key);
	    if (found == this->sources.end())
		continue;
	    TelemetryJsonObject data;
	    data.integer("source_id", found->second.id);
	    this->sink->write(this->controllerId, "inventory_source_removed",
		data.take());
	    this->sources.erase(found);
	}
	this->inventoryObserved = true;
	this->lastInventoryRevision = inventoryRevision;
    }

    std::shared_ptr<BObolLodTelemetrySink> sink;
    TelemetryAliases aliases;
    std::unordered_map<std::string, SourceCursor> sources;
    uint64_t controllerId = 0;
    uint64_t transitionSerial = 0;
    uint64_t nextSourceId = 1;
    uint64_t lastInventoryRevision = 0;
    bool inventoryObserved = false;
    bool started = false;
};

std::unique_ptr<BObolLodTelemetry>
BObolLodTelemetry::fromEnvironment(void) noexcept
{
    try {
	const char *path = std::getenv(bobol_lod_telemetry_file_environment);
	if (!path || !path[0])
	    return std::unique_ptr<BObolLodTelemetry>();
	const char *sanitizeSetting = std::getenv(
	    bobol_lod_telemetry_sanitize_environment);
	const bool sanitized = sanitizeSetting &&
	    bu_str_true(sanitizeSetting) != 0;
	std::shared_ptr<BObolLodTelemetrySink> sink =
	    telemetry_acquire_sink(path, sanitized);
	if (!sink)
	    return std::unique_ptr<BObolLodTelemetry>();
	return std::unique_ptr<BObolLodTelemetry>(new BObolLodTelemetry(
	    std::make_unique<Impl>(std::move(sink))));
    } catch (...) {
	bu_log("Unable to initialize BObol LoD telemetry; logging disabled\n");
	return std::unique_ptr<BObolLodTelemetry>();
    }
}

BObolLodTelemetry::BObolLodTelemetry(
    std::unique_ptr<BObolLodTelemetry::Impl> implementation) :
    impl(std::move(implementation))
{
}

BObolLodTelemetry::~BObolLodTelemetry(void) = default;

void
BObolLodTelemetry::recordTransition(const BObolViewController &controller,
    const BObolLodControlTransitionRecord &transition) noexcept
{
    try {
	if (this->impl)
	    this->impl->recordTransition(controller, transition);
    } catch (...) {
	if (this->impl)
	    this->impl->recordGap();
    }
}

void
BObolLodTelemetry::recordPublication(
    const BObolLodPublicationTelemetryRecord &publication) noexcept
{
    try {
	if (this->impl)
	    this->impl->recordPublication(publication);
    } catch (...) {
	if (this->impl)
	    this->impl->recordGap();
    }
}

void
BObolLodTelemetry::recordGap(void) noexcept
{
    try {
	if (this->impl)
	    this->impl->recordGap();
    } catch (...) {
    }
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
