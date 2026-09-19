/* M O D E L _ C O N F O R M A N C E _ H A R N E S S _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_TESTS_MODEL_CONFORMANCE_HARNESS_PRIVATE_H
#define LIBBOBOL_TESTS_MODEL_CONFORMANCE_HARNESS_PRIVATE_H

#include "../lod_result_authentication_private.h"

#include "BObol/BViewController.h"
#include "bu/log.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <utility>

enum BObolModelObservationDomain : uint32_t {
    BOBOL_MODEL_OBSERVES_NONE = 0,
    BOBOL_MODEL_OBSERVES_CONTROL = 1u << 0,
    BOBOL_MODEL_OBSERVES_RESULT = 1u << 1,
    BOBOL_MODEL_OBSERVES_TRANSACTION = 1u << 2,
    BOBOL_MODEL_OBSERVES_IDENTITY = 1u << 3,
    BOBOL_MODEL_OBSERVES_EVIDENCE = 1u << 4,
    BOBOL_MODEL_OBSERVES_SEMANTIC_PRESENTATION = 1u << 5
};

/*
 * Test-side projection shared by every conformance scenario.  Each field is
 * copied from a production reducer or public controller observation; the
 * harness never uses this state to schedule production work.  Domain bits
 * distinguish an intentionally absent concept from a meaningful zero value.
 */
struct BObolModelConformanceState {
    uint32_t domains = BOBOL_MODEL_OBSERVES_NONE;

    uint32_t controlFactMask = 0;
    uint32_t controlObligationMask = 0;
    int controlOwner = BOBOL_LOD_CONTROL_OWNER_NONE;
    uint32_t controlViolationMask = 0;
    uint32_t presentationWitnessMask = 0;
    bool terminal = true;
    bool terminalError = false;
    bool viewReady = true;
    uint32_t constraintEvidenceMask = 0;
    BObolHostWorkSnapshot hostWork;

    uint64_t inventoryEpoch = 0;
    uint64_t availabilityEpoch = 0;
    uint64_t visibilityEpoch = 0;
    uint64_t capacityEpoch = 0;
    uint64_t cadEpoch = 0;
    uint64_t residentDemandEpoch = 0;
    uint64_t allocationPlanIdentity = 0;
    uint64_t completedRenderSerial = 0;
    uint64_t presentedFrameIdentity = 0;
    bool publicationFramePending = false;
    bool refinementFramePending = false;
    bool hasLodState = false;

    uint64_t activeGeneration = 0;
    uint64_t sourceRoutingId = 0;
    uint64_t sourcePopulationEpoch = 0;
    uint64_t viewEpoch = 0;
    uint64_t policyEpoch = 0;
    uint8_t resultMismatchMask = 0;
    BObolLodResultDisposition resultDisposition =
	BObolLodResultDisposition::RETRY_CURRENT_DEMAND;

    bool transactionActive = false;
    bool transactionBarrierPending = false;
    uint64_t transactionIdentity = 0;
    uint64_t requiredRenderSerial = 0;

    bool semanticPresentationPending = false;
    bool semanticRequestRequired = false;
    bool semanticFramePending = false;

    uint64_t sourceIdentity = 0;
    uint64_t evidenceIdentity = 0;
    uint64_t evidenceValue = 0;
    bool evidenceValid = false;
};

inline BObolModelConformanceState
bobol_model_conformance_state(const BObolLodConvergenceStatus &status,
    const BObolHostWorkSnapshot &hostWork)
{
    BObolModelConformanceState state;
    state.domains = BOBOL_MODEL_OBSERVES_CONTROL |
	BOBOL_MODEL_OBSERVES_TRANSACTION;
    state.controlFactMask = status.controlFactMask;
    state.controlObligationMask = status.controlObligationMask;
    state.controlOwner = status.controlOwner;
    state.controlViolationMask = status.controlViolationMask;
    state.presentationWitnessMask =
	status.controlPresentationWitnessMask;
    state.terminal = status.terminal != FALSE;
    state.terminalError = status.terminalError != FALSE;
    state.viewReady = status.viewReady != FALSE;
    state.constraintEvidenceMask = status.constraintEvidenceMask;
    state.hostWork = hostWork;
    state.inventoryEpoch = status.inventoryRevision;
    state.availabilityEpoch = status.availabilityRevision;
    state.visibilityEpoch = status.visibilityRevision;
    state.capacityEpoch = status.capacityRevision;
    state.cadEpoch = status.cadRevision;
    state.residentDemandEpoch = status.residentDemandRevision;
    state.allocationPlanIdentity = status.currentAllocationPlanSerial;
    state.presentedFrameIdentity = status.presentedFrameSerial;
    state.publicationFramePending = status.publicationFramePending != FALSE;
    state.refinementFramePending = status.refinementFramePending != FALSE;
    state.hasLodState = status.hasLodState != FALSE;
    state.activeGeneration = status.activeGeneration;
    state.viewEpoch = status.viewRevision;
    state.policyEpoch = status.policyRevision;
    state.transactionActive = state.publicationFramePending ||
	status.presentationRequiredRenderSerial != 0;
    state.transactionBarrierPending =
	status.presentationRequiredRenderSerial != 0;
    state.transactionIdentity = status.presentationTransactionSerial;
    state.requiredRenderSerial = status.presentationRequiredRenderSerial;
    return state;
}

inline BObolModelConformanceState
bobol_model_conformance_state(const BObolLodControlTraceState &traceState)
{
    BObolModelConformanceState state = bobol_model_conformance_state(
	traceState.convergence, traceState.hostWork);
    state.completedRenderSerial = traceState.renderCompletionSerial;
    return state;
}

inline uint32_t
bobol_model_owner_obligation(int owner)
{
    switch (owner) {
	case BOBOL_LOD_CONTROL_OWNER_INTERACTION:
	    return BOBOL_LOD_CONTROL_OBLIGATION_INTERACTION;
	case BOBOL_LOD_CONTROL_OWNER_INVENTORY:
	    return BOBOL_LOD_CONTROL_OBLIGATION_INVENTORY;
	case BOBOL_LOD_CONTROL_OWNER_AVAILABILITY:
	    return BOBOL_LOD_CONTROL_OBLIGATION_AVAILABILITY;
	case BOBOL_LOD_CONTROL_OWNER_PUBLICATION:
	    return BOBOL_LOD_CONTROL_OBLIGATION_PUBLICATION;
	case BOBOL_LOD_CONTROL_OWNER_PLANNING:
	    return BOBOL_LOD_CONTROL_OBLIGATION_PLANNING;
	case BOBOL_LOD_CONTROL_OWNER_PRESENTATION:
	    return BOBOL_LOD_CONTROL_OBLIGATION_PRESENTATION;
	case BOBOL_LOD_CONTROL_OWNER_HANDOFF:
	    return BOBOL_LOD_CONTROL_OBLIGATION_HANDOFF;
	case BOBOL_LOD_CONTROL_OWNER_COMPACTION:
	    return BOBOL_LOD_CONTROL_OBLIGATION_COMPACTION;
	case BOBOL_LOD_CONTROL_OWNER_CACHE_WRITE:
	    return BOBOL_LOD_CONTROL_OBLIGATION_CACHE_WRITE;
	default:
	    return BOBOL_LOD_CONTROL_OBLIGATION_NONE;
    }
}

inline const char *
bobol_model_conformance_state_violation(
    const BObolModelConformanceState &state)
{
    if (state.domains & BOBOL_MODEL_OBSERVES_CONTROL) {
	const uint32_t background =
	    BOBOL_LOD_CONTROL_OBLIGATION_COMPACTION |
	    BOBOL_LOD_CONTROL_OBLIGATION_CACHE_WRITE;
	const uint32_t foreground =
	    state.controlObligationMask & ~background;
	if (state.controlViolationMask != 0)
	    return "production control projection reports a violation";
	if ((state.controlObligationMask == 0) !=
		(state.controlOwner == BOBOL_LOD_CONTROL_OWNER_NONE))
	    return "work and its selected owner disagree";
	const uint32_t owned = bobol_model_owner_obligation(
	    state.controlOwner);
	if (state.controlOwner != BOBOL_LOD_CONTROL_OWNER_NONE &&
		!(state.controlObligationMask & owned))
	    return "selected owner has no matching obligation";
	if (state.terminal && foreground != 0)
	    return "terminal observation retains foreground work";
	if (state.viewReady && (!state.terminal || state.terminalError))
	    return "ready observation is not a successful terminal state";
	if (state.controlOwner == BOBOL_LOD_CONTROL_OWNER_PRESENTATION &&
		state.presentationWitnessMask == 0)
	    return "presentation owner has no finite progress witness";
	if (state.hostWork.renderPending() &&
		state.hostWork.renderRevision == 0)
	    return "pending render has no transaction identity";
    }
    if (state.domains & BOBOL_MODEL_OBSERVES_RESULT) {
	const bool current = state.resultMismatchMask == 0;
	const bool superseded = state.resultDisposition ==
	    BObolLodResultDisposition::SUPERSEDE;
	if (current == superseded)
	    return "result disposition disagrees with exact identity";
	if (!state.sourceRoutingId || !state.sourcePopulationEpoch ||
		!state.viewEpoch || !state.policyEpoch)
	    return "compact result observation has an empty credential";
    }
    if (state.domains & BOBOL_MODEL_OBSERVES_TRANSACTION) {
	if (state.transactionActive &&
		(!state.transactionIdentity || !state.viewEpoch ||
		 !state.policyEpoch))
	    return "active transaction has an incomplete identity";
	if (state.transactionBarrierPending !=
		(state.requiredRenderSerial != 0))
	    return "frame barrier and required render identity disagree";
	if (state.transactionBarrierPending && !state.transactionActive)
	    return "frame barrier has no active transaction";
    }
    if (state.domains & BOBOL_MODEL_OBSERVES_SEMANTIC_PRESENTATION) {
	const bool owned = state.semanticRequestRequired ||
	    state.semanticFramePending;
	if (state.semanticPresentationPending != owned ||
		(state.semanticRequestRequired && state.semanticFramePending))
	    return "semantic presentation has no exclusive successor owner";
    }
    if ((state.domains & BOBOL_MODEL_OBSERVES_IDENTITY) &&
	    !state.sourceIdentity)
	return "source identity is empty";
    if ((state.domains & BOBOL_MODEL_OBSERVES_EVIDENCE) &&
	    state.evidenceValid && !state.evidenceIdentity &&
	    !state.evidenceValue)
	return "valid evidence has neither an identity nor a measured value";
    return NULL;
}

inline bool
bobol_model_conformance_state_equal(
    const BObolModelConformanceState &left,
    const BObolModelConformanceState &right)
{
    return left.domains == right.domains &&
	left.controlFactMask == right.controlFactMask &&
	left.controlObligationMask == right.controlObligationMask &&
	left.controlOwner == right.controlOwner &&
	left.controlViolationMask == right.controlViolationMask &&
	left.presentationWitnessMask == right.presentationWitnessMask &&
	left.terminal == right.terminal &&
	left.terminalError == right.terminalError &&
	left.viewReady == right.viewReady &&
	left.constraintEvidenceMask == right.constraintEvidenceMask &&
	left.hostWork.revision == right.hostWork.revision &&
	left.hostWork.renderRevision == right.hostWork.renderRevision &&
	left.hostWork.flags == right.hostWork.flags &&
	left.inventoryEpoch == right.inventoryEpoch &&
	left.availabilityEpoch == right.availabilityEpoch &&
	left.visibilityEpoch == right.visibilityEpoch &&
	left.capacityEpoch == right.capacityEpoch &&
	left.cadEpoch == right.cadEpoch &&
	left.residentDemandEpoch == right.residentDemandEpoch &&
	left.allocationPlanIdentity == right.allocationPlanIdentity &&
	left.completedRenderSerial == right.completedRenderSerial &&
	left.presentedFrameIdentity == right.presentedFrameIdentity &&
	left.publicationFramePending == right.publicationFramePending &&
	left.refinementFramePending == right.refinementFramePending &&
	left.hasLodState == right.hasLodState &&
	left.activeGeneration == right.activeGeneration &&
	left.sourceRoutingId == right.sourceRoutingId &&
	left.sourcePopulationEpoch == right.sourcePopulationEpoch &&
	left.viewEpoch == right.viewEpoch &&
	left.policyEpoch == right.policyEpoch &&
	left.resultMismatchMask == right.resultMismatchMask &&
	left.resultDisposition == right.resultDisposition &&
	left.transactionActive == right.transactionActive &&
	left.transactionBarrierPending == right.transactionBarrierPending &&
	left.transactionIdentity == right.transactionIdentity &&
	left.requiredRenderSerial == right.requiredRenderSerial &&
	left.semanticPresentationPending ==
	    right.semanticPresentationPending &&
	left.semanticRequestRequired == right.semanticRequestRequired &&
	left.semanticFramePending == right.semanticFramePending &&
	left.sourceIdentity == right.sourceIdentity &&
	left.evidenceIdentity == right.evidenceIdentity &&
	left.evidenceValue == right.evidenceValue &&
	left.evidenceValid == right.evidenceValid;
}

class BObolModelConformanceScenario {
public:
    typedef std::function<BObolModelConformanceState(void)> Observer;
    typedef std::function<void(void)> Mutation;
    typedef std::function<bool(const BObolModelConformanceState &,
	const BObolModelConformanceState &)> Expectation;

    BObolModelConformanceScenario(const char *modelName,
	const char *scenarioName, Observer observer) :
	model(modelName ? modelName : ""),
	name(scenarioName ? scenarioName : ""),
	observe(std::move(observer))
    {
    }

    bool checkInitial(const char *actionName = "Init")
    {
	const BObolModelConformanceState state = this->observe();
	return this->checkState(actionName, state);
    }

    bool step(const char *actionName, const Mutation &mutation,
	const Expectation &expectation = Expectation())
    {
	const BObolModelConformanceState before = this->observe();
	mutation();
	const BObolModelConformanceState after = this->observe();
	++this->stepCount;
	if (!this->checkState(actionName, after))
	    return false;
	if (expectation && !expectation(before, after)) {
	    this->report(actionName, "scenario postcondition failed");
	    return false;
	}
	return true;
    }

private:
    bool checkState(const char *actionName,
	const BObolModelConformanceState &state) const
    {
	if (!actionName || !actionName[0]) {
	    this->report("<unnamed>", "step has no action name");
	    return false;
	}
	const char *violation = bobol_model_conformance_state_violation(state);
	if (violation) {
	    this->report(actionName, violation);
	    return false;
	}
	return true;
    }

    void report(const char *actionName, const char *message) const
    {
	bu_log("FAIL: model=%s scenario=%s step=%zu action=%s: %s\n",
	    this->model.c_str(), this->name.c_str(), this->stepCount,
	    actionName ? actionName : "<null>", message);
    }

    std::string model;
    std::string name;
    Observer observe;
    size_t stepCount = 0;
};

#endif /* LIBBOBOL_TESTS_MODEL_CONFORMANCE_HARNESS_PRIVATE_H */
