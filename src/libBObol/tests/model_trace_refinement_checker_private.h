/*      M O D E L _ T R A C E _ R E F I N E M E N T _ C H E C K E R . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_TESTS_MODEL_TRACE_REFINEMENT_CHECKER_PRIVATE_H
#define LIBBOBOL_TESTS_MODEL_TRACE_REFINEMENT_CHECKER_PRIVATE_H

#include "model_conformance_harness_private.h"

#include "../lod_control_private.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

enum class BObolTraceRefinementFailureCode : uint8_t {
    NONE = 0,
    EMPTY_TRACE,
    TRACE_LOSS,
    INVALID_INITIAL,
    NONCONTIGUOUS_SERIAL,
    UNMAPPED_EVENT,
    DISCONTINUOUS_ENDPOINT,
    INVALID_STATE,
    CONTROL_PROJECTION_MISMATCH,
    EVENT_OWNER_MISMATCH,
    OBSERVER_MISMATCH,
    REVISION_REGRESSION,
    CLOCK_REGRESSION,
    BARRIER_IDENTITY_CHANGED,
    BARRIER_RETIRED_WITHOUT_FRAME,
    BARRIER_RETIRED_BEFORE_REQUIRED_FRAME,
    PUBLICATION_RETIRED_WITHOUT_BARRIER,
    PUBLICATION_BARRIER_SUPERSEDED,
    MISSING_PUBLICATION_PRESENTATION_FLOW,
    UNRETIRED_PUBLICATION_BARRIER
};

struct BObolTraceRefinementFailure {
    BObolTraceRefinementFailureCode code =
	BObolTraceRefinementFailureCode::NONE;
    size_t recordIndex = 0;
    uint64_t serial = 0;
    BObolLodControlTransitionEvent event =
	BOBOL_LOD_CONTROL_TRANSITION_UNNAMED;
    std::string detail;
};

struct BObolTraceRefinementStatistics {
    size_t publicationBarriers = 0;
    size_t completedPublicationBarriers = 0;
};

struct BObolTraceRefinementResult {
    bool valid = false;
    BObolTraceRefinementFailure failure;
    BObolTraceRefinementStatistics statistics;
};

/*
 * Offline observer for the production transition journal.  This checker is
 * deliberately test-side: it may reject an execution but must never select a
 * production successor.  Its independent fact projection is a refinement
 * oracle for ObolControlRefinement.tla rather than another control reducer.
 */
class BObolResultPresentationTraceChecker {
public:
    static BObolTraceRefinementResult check(
	const std::vector<BObolLodControlTransitionRecord> &records,
	uint64_t droppedRecordCount)
    {
	BObolTraceRefinementResult result;
	if (records.empty())
	    return failure(result, BObolTraceRefinementFailureCode::EMPTY_TRACE,
		0, NULL, "trace contains no initial endpoint");
	if (droppedRecordCount != 0) {
	    std::ostringstream detail;
	    detail << "journal dropped " << droppedRecordCount << " record(s)";
	    return failure(result, BObolTraceRefinementFailureCode::TRACE_LOSS,
		0, &records.front(), detail.str());
	}

	const BObolModelConformanceState initialBefore =
	    bobol_model_conformance_state(records.front().before);
	const BObolModelConformanceState initialAfter =
	    bobol_model_conformance_state(records.front().after);
	if (records.front().serial != FirstTraceSerial ||
		records.front().event != BOBOL_LOD_CONTROL_TRANSITION_INITIAL ||
		!bobol_model_conformance_state_equal(
		    initialBefore, initialAfter))
	    return failure(result,
		BObolTraceRefinementFailureCode::INVALID_INITIAL, 0,
		&records.front(),
		"first record is not one unchanged initial endpoint");

	bool publicationBarrierPending = false;
	uint64_t publicationTransaction = 0;
	uint64_t publicationView = 0;
	uint64_t publicationPolicy = 0;
	uint64_t publicationRequiredRender = 0;
	bool resultPublicationObserved = false;

	for (size_t index = 0; index < records.size(); ++index) {
	    const BObolLodControlTransitionRecord &record = records[index];
	    const BObolModelConformanceState before =
		bobol_model_conformance_state(record.before);
	    const BObolModelConformanceState after =
		bobol_model_conformance_state(record.after);

	    if (!knownEvent(record.event))
		return failure(result,
		    BObolTraceRefinementFailureCode::UNMAPPED_EVENT, index,
		    &record, "event is unnamed or outside the finite alphabet");
	    if (index > 0) {
		const BObolLodControlTransitionRecord &previous =
		    records[index - 1];
		if (previous.serial == std::numeric_limits<uint64_t>::max() ||
			record.serial != previous.serial + 1)
		    return failure(result,
			BObolTraceRefinementFailureCode::NONCONTIGUOUS_SERIAL,
			index, &record,
			"controller transition serial is not contiguous");
		if (!bobol_model_conformance_state_equal(
			bobol_model_conformance_state(previous.after), before))
		    return failure(result,
			BObolTraceRefinementFailureCode::DISCONTINUOUS_ENDPOINT,
			index, &record,
			"before endpoint differs from the preceding after endpoint");
		if (record.event == BOBOL_LOD_CONTROL_TRANSITION_INITIAL)
		    return failure(result,
			BObolTraceRefinementFailureCode::INVALID_INITIAL, index,
			&record, "initial appears after the first record");
	    }

	    const char *beforeViolation =
		bobol_model_conformance_state_violation(before);
	    const char *afterViolation =
		bobol_model_conformance_state_violation(after);
	    if (beforeViolation || afterViolation) {
		std::ostringstream detail;
		detail << (beforeViolation ? "before: " : "after: ")
		    << (beforeViolation ? beforeViolation : afterViolation)
		    << " (before_mask=" << before.controlViolationMask
		    << ", after_mask=" << after.controlViolationMask
		    << ", before={facts:" << before.controlFactMask
		    << ", obligations:" << before.controlObligationMask
		    << ", owner:" << before.controlOwner
		    << ", terminal:" << before.terminal << "}"
		    << ", after={facts:" << after.controlFactMask
		    << ", obligations:" << after.controlObligationMask
		    << ", owner:" << after.controlOwner
		    << ", terminal:" << after.terminal << "})";
		return failure(result,
		    BObolTraceRefinementFailureCode::INVALID_STATE, index,
		    &record, detail.str());
	    }

	    if (!observerMatches(record.before) ||
		    !observerMatches(record.after))
		return failure(result,
		    BObolTraceRefinementFailureCode::OBSERVER_MISMATCH, index,
		    &record,
		    "trace and convergence observers disagree on view or policy");

	    if (!controlProjectionMatches(before) ||
		    !controlProjectionMatches(after))
		return failure(result,
		    BObolTraceRefinementFailureCode::CONTROL_PROJECTION_MISMATCH,
		    index, &record,
		    "fact mask does not refine to the reported obligation and owner");

	    const bool ownerlessResultNotification = record.event ==
		BOBOL_LOD_CONTROL_TRANSITION_PUBLICATION &&
		before.controlOwner == BOBOL_LOD_CONTROL_OWNER_NONE &&
		after.controlOwner == BOBOL_LOD_CONTROL_OWNER_NONE;
	    const bool resultPublicationNotification =
		ownerlessResultNotification ||
		(record.event == BOBOL_LOD_CONTROL_TRANSITION_PUBLICATION &&
		 (factActive(before.controlFactMask, Refinement::Fact::RESULT) ||
		  factActive(after.controlFactMask, Refinement::Fact::RESULT)));
	    const int eventOwner = ownerForEvent(record.event);
	    if (eventOwner != UnmappedOwner) {
		/* A reducer may discharge its own obligation and expose the next
		 * owner, or introduce its obligation while another owner's work is
		 * still visible at the preceding endpoint.  Either endpoint is a
		 * valid ownership witness. */
		const bool endpointOwnsEvent =
		    eventOwner == before.controlOwner ||
		    eventOwner == after.controlOwner;
		if (!endpointOwnsEvent &&
			!ownerlessResultNotification) {
		    std::ostringstream detail;
		    detail << "owner event selects " << eventOwner
			<< " but neither endpoint owns it (before="
			<< before.controlOwner
			<< ", after=" << after.controlOwner << ")";
		    return failure(result,
			BObolTraceRefinementFailureCode::EVENT_OWNER_MISMATCH,
			index, &record, detail.str());
		}
	    }

	    if (regressesRevision(before, after))
		return failure(result,
		    BObolTraceRefinementFailureCode::REVISION_REGRESSION,
		    index, &record, "a semantic revision moved backward");
	    if (record.event != BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT &&
		    (before.viewEpoch != after.viewEpoch ||
		     before.policyEpoch != after.policyEpoch))
		return failure(result,
		    BObolTraceRefinementFailureCode::REVISION_REGRESSION,
		    index, &record,
		    "a non-input transition changed the demand identity");
	    if (regressesClock(before, after))
		return failure(result,
		    BObolTraceRefinementFailureCode::CLOCK_REGRESSION,
		    index, &record, "an ordering clock moved backward");

	    if (before.requiredRenderSerial != 0 &&
		    after.requiredRenderSerial != 0 &&
		    before.transactionIdentity == after.transactionIdentity &&
		    before.requiredRenderSerial != after.requiredRenderSerial)
		return failure(result,
		    BObolTraceRefinementFailureCode::BARRIER_IDENTITY_CHANGED,
		    index, &record,
		    "an active transaction moved its required frame identity");

	    const bool barrierRetired = before.requiredRenderSerial != 0 &&
		after.requiredRenderSerial == 0;
	    const bool publicationRetired = before.publicationFramePending &&
		!after.publicationFramePending;
	    const bool completedFrame =
		after.completedRenderSerial > before.completedRenderSerial;
	    if (barrierRetired &&
		    record.event != BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT &&
		    !completedFrame)
		return failure(result,
		    BObolTraceRefinementFailureCode::BARRIER_RETIRED_WITHOUT_FRAME,
		    index, &record,
		    "a non-cancellation transition retired a barrier without a frame");
	    if (barrierRetired && completedFrame &&
		    after.completedRenderSerial < before.requiredRenderSerial)
		return failure(result,
		    BObolTraceRefinementFailureCode::
			BARRIER_RETIRED_BEFORE_REQUIRED_FRAME,
		    index, &record,
		    "the completed frame did not reach the barrier identity");

	    if (!barrierRetired && resultPublicationNotification)
		resultPublicationObserved = true;
	    if (record.event == BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT)
		resultPublicationObserved = false;

	    const bool barrierArmed = before.requiredRenderSerial == 0 &&
		after.requiredRenderSerial != 0;
	    if (barrierArmed && resultPublicationObserved) {
		publicationBarrierPending = true;
		publicationTransaction = after.transactionIdentity;
		publicationView = after.viewEpoch;
		publicationPolicy = after.policyEpoch;
		publicationRequiredRender = after.requiredRenderSerial;
		resultPublicationObserved = false;
		++result.statistics.publicationBarriers;
	    }
	    if (publicationRetired && resultPublicationObserved &&
		    !publicationBarrierPending)
		return failure(result,
		    BObolTraceRefinementFailureCode::
			PUBLICATION_RETIRED_WITHOUT_BARRIER,
		    index, &record,
		    "result publication retired without an exact-frame barrier");

	    if (publicationBarrierPending) {
		const bool demandChanged = after.viewEpoch != publicationView ||
		    after.policyEpoch != publicationPolicy;
		if (demandChanged) {
		    if (record.event !=
			    BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT)
			return failure(result,
			    BObolTraceRefinementFailureCode::
				PUBLICATION_BARRIER_SUPERSEDED,
			    index, &record,
			    "publication barrier changed demand without input");
		    publicationBarrierPending = false;
		} else if (after.transactionIdentity != publicationTransaction) {
		    return failure(result,
			BObolTraceRefinementFailureCode::
			    PUBLICATION_BARRIER_SUPERSEDED,
			index, &record,
			"publication barrier was replaced before retirement");
		} else if (after.requiredRenderSerial == 0) {
		    if (!completedFrame ||
			    after.completedRenderSerial <
				publicationRequiredRender)
			return failure(result,
			    BObolTraceRefinementFailureCode::
				BARRIER_RETIRED_WITHOUT_FRAME,
			    index, &record,
			    "publication barrier ended without its exact frame");
		    ++result.statistics.completedPublicationBarriers;
		    publicationBarrierPending = false;
		}
	    }
	}

	if (publicationBarrierPending)
	    return failure(result,
		BObolTraceRefinementFailureCode::UNRETIRED_PUBLICATION_BARRIER,
		records.size() - 1, &records.back(),
		"trace ended with a result-publication barrier still pending");
	if (result.statistics.completedPublicationBarriers == 0)
	    return failure(result,
		BObolTraceRefinementFailureCode::
		    MISSING_PUBLICATION_PRESENTATION_FLOW,
		records.size() - 1, &records.back(),
		"trace contains no completed result-to-presentation transaction");

	result.valid = true;
	return result;
    }

private:
    using Refinement = BObolLodControlRefinement;
    static constexpr uint64_t FirstTraceSerial = 1;
    static constexpr int UnmappedOwner = -1;

    static BObolTraceRefinementResult failure(
	BObolTraceRefinementResult result,
	BObolTraceRefinementFailureCode code, size_t index,
	const BObolLodControlTransitionRecord *record,
	const std::string &detail)
    {
	result.valid = false;
	result.failure.code = code;
	result.failure.recordIndex = index;
	result.failure.detail = detail;
	if (record) {
	    result.failure.serial = record->serial;
	    result.failure.event = record->event;
	}
	return result;
    }

    static bool knownEvent(BObolLodControlTransitionEvent event)
    {
	return event >= BOBOL_LOD_CONTROL_TRANSITION_INITIAL &&
	    event <= BOBOL_LOD_CONTROL_TRANSITION_IDLE_SERVICE;
    }

    static int ownerForEvent(BObolLodControlTransitionEvent event)
    {
	switch (event) {
	    case BOBOL_LOD_CONTROL_TRANSITION_INTERACTION:
		return BOBOL_LOD_CONTROL_OWNER_INTERACTION;
	    case BOBOL_LOD_CONTROL_TRANSITION_INVENTORY:
		return BOBOL_LOD_CONTROL_OWNER_INVENTORY;
	    case BOBOL_LOD_CONTROL_TRANSITION_AVAILABILITY:
		return BOBOL_LOD_CONTROL_OWNER_AVAILABILITY;
	    case BOBOL_LOD_CONTROL_TRANSITION_PUBLICATION:
		return BOBOL_LOD_CONTROL_OWNER_PUBLICATION;
	    case BOBOL_LOD_CONTROL_TRANSITION_PLANNING:
		return BOBOL_LOD_CONTROL_OWNER_PLANNING;
	    case BOBOL_LOD_CONTROL_TRANSITION_PRESENTATION:
		return BOBOL_LOD_CONTROL_OWNER_PRESENTATION;
	    case BOBOL_LOD_CONTROL_TRANSITION_HANDOFF:
		return BOBOL_LOD_CONTROL_OWNER_HANDOFF;
	    case BOBOL_LOD_CONTROL_TRANSITION_COMPACTION:
		return BOBOL_LOD_CONTROL_OWNER_COMPACTION;
	    case BOBOL_LOD_CONTROL_TRANSITION_CACHE_WRITE:
		return BOBOL_LOD_CONTROL_OWNER_CACHE_WRITE;
	    case BOBOL_LOD_CONTROL_TRANSITION_IDLE_SERVICE:
		return BOBOL_LOD_CONTROL_OWNER_NONE;
	    default:
		return UnmappedOwner;
	}
    }

    static bool observerMatches(const BObolLodControlTraceState &state)
    {
	return state.viewRevision == state.convergence.viewRevision &&
	    state.policyRevision == state.convergence.policyRevision;
    }

    static bool factActive(uint32_t facts, Refinement::Fact fact)
    {
	return (facts & Refinement::bit(fact)) != 0;
    }

    static Refinement::Inputs controlInputs(uint32_t facts)
    {
	Refinement::Inputs inputs;
	inputs.interaction = factActive(facts, Refinement::Fact::INTERACTION);
	inputs.inventory = factActive(facts, Refinement::Fact::INVENTORY);
	inputs.availability = factActive(facts, Refinement::Fact::AVAILABILITY);
	inputs.result = factActive(facts, Refinement::Fact::RESULT);
	inputs.publication = factActive(facts, Refinement::Fact::PUBLICATION);
	inputs.submission = factActive(facts, Refinement::Fact::SUBMISSION);
	inputs.demandRefresh = factActive(
	    facts, Refinement::Fact::DEMAND_REFRESH);
	inputs.submissionRescan = factActive(
	    facts, Refinement::Fact::SUBMISSION_RESCAN);
	inputs.submissionDelta = factActive(
	    facts, Refinement::Fact::SUBMISSION_DELTA);
	inputs.qualityProbe = factActive(facts, Refinement::Fact::QUALITY_PROBE);
	inputs.retainedAllocation = factActive(
	    facts, Refinement::Fact::RETAINED_ALLOCATION);
	inputs.retainedAllocationTransaction = factActive(
	    facts, Refinement::Fact::RETAINED_ALLOCATION_TRANSACTION);
	inputs.importanceCensus = factActive(
	    facts, Refinement::Fact::IMPORTANCE_CENSUS);
	inputs.residentAdmissionRetry = factActive(
	    facts, Refinement::Fact::RESIDENT_ADMISSION_RETRY);
	inputs.capacityAllocation = factActive(
	    facts, Refinement::Fact::CAPACITY_ALLOCATION);
	inputs.residentGrowth = factActive(
	    facts, Refinement::Fact::RESIDENT_GROWTH);
	inputs.pointTriangleRecovery = factActive(
	    facts, Refinement::Fact::POINT_TRIANGLE_RECOVERY);
	inputs.structuralFrontier = factActive(
	    facts, Refinement::Fact::STRUCTURAL_FRONTIER);
	inputs.presentationReplay = factActive(
	    facts, Refinement::Fact::PRESENTATION_REPLAY);
	inputs.exactPresentation = factActive(
	    facts, Refinement::Fact::EXACT_PRESENTATION);
	inputs.presentationBarrier = factActive(
	    facts, Refinement::Fact::PRESENTATION_BARRIER);
	inputs.capacityFrame = factActive(
	    facts, Refinement::Fact::CAPACITY_FRAME);
	inputs.pointAdmissionFrame = factActive(
	    facts, Refinement::Fact::POINT_ADMISSION_FRAME);
	inputs.pointCalibration = factActive(
	    facts, Refinement::Fact::POINT_CALIBRATION);
	inputs.capacityCalibration = factActive(
	    facts, Refinement::Fact::CAPACITY_CALIBRATION);
	inputs.headroomProbe = factActive(
	    facts, Refinement::Fact::HEADROOM_PROBE);
	inputs.handoff = factActive(facts, Refinement::Fact::HANDOFF);
	inputs.compaction = factActive(facts, Refinement::Fact::COMPACTION);
	inputs.cacheWrite = factActive(facts, Refinement::Fact::CACHE_WRITE);
	return inputs;
    }

    static bool controlProjectionMatches(
	const BObolModelConformanceState &state)
    {
	const Refinement::Inputs inputs = controlInputs(state.controlFactMask);
	if (Refinement::factMask(inputs) != state.controlFactMask)
	    return false;
	const Refinement::Snapshot expected = Refinement::evaluate(inputs);
	return expected.obligations == state.controlObligationMask &&
	    static_cast<int>(expected.owner) == state.controlOwner;
    }

    static bool regressesRevision(const BObolModelConformanceState &before,
	const BObolModelConformanceState &after)
    {
	return after.inventoryEpoch < before.inventoryEpoch ||
	    after.availabilityEpoch < before.availabilityEpoch ||
	    after.visibilityEpoch < before.visibilityEpoch ||
	    after.viewEpoch < before.viewEpoch ||
	    after.policyEpoch < before.policyEpoch ||
	    after.capacityEpoch < before.capacityEpoch ||
	    after.cadEpoch < before.cadEpoch ||
	    after.residentDemandEpoch < before.residentDemandEpoch;
    }

    static bool regressesClock(const BObolModelConformanceState &before,
	const BObolModelConformanceState &after)
    {
	return after.transactionIdentity < before.transactionIdentity ||
	    after.completedRenderSerial < before.completedRenderSerial ||
	    after.presentedFrameIdentity < before.presentedFrameIdentity ||
	    after.hostWork.revision < before.hostWork.revision;
    }
};

#endif /* LIBBOBOL_TESTS_MODEL_TRACE_REFINEMENT_CHECKER_PRIVATE_H */
