/*            T E S T _ M O D E L _ C O N F O R M A N C E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "model_conformance_harness_private.h"
#include "model_trace_refinement_checker_private.h"

#include "exact_sample_identity_private.h"
#include "identity_counter_private.h"
#include "lod_control_private.h"
#include "lod_delivery_policy_private.h"
#include "BObol/BInit.h"
#include "BObol/BLodService.h"
#include "bu/app.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>

namespace {

using ConformanceState = BObolModelConformanceState;
using Scenario = BObolModelConformanceScenario;

constexpr uint64_t SourceRoute = 11;
constexpr uint64_t SourcePopulation = 13;
constexpr uint64_t ViewEpoch = 17;
constexpr uint64_t PolicyEpoch = 19;

unsigned int
set_bit_count(uint32_t value)
{
    unsigned int count = 0;
    while (value) {
	value &= value - 1;
	++count;
    }
    return count;
}

void
discharge_control_owner(BObolLodControlRefinement::Inputs &inputs)
{
    using Owner = BObolLodControlRefinement::Owner;
    switch (BObolLodControlRefinement::evaluate(inputs).owner) {
	case Owner::INTERACTION:
	    inputs.interaction = false;
	    break;
	case Owner::INVENTORY:
	    inputs.inventory = false;
	    break;
	case Owner::AVAILABILITY:
	    inputs.availability = false;
	    break;
	case Owner::PUBLICATION:
	    inputs.result = false;
	    inputs.publication = false;
	    break;
	case Owner::PLANNING:
	    inputs.submission = false;
	    inputs.demandRefresh = false;
	    inputs.submissionRescan = false;
	    inputs.submissionDelta = false;
	    inputs.qualityProbe = false;
	    inputs.retainedAllocation = false;
	    inputs.retainedAllocationTransaction = false;
	    inputs.importanceCensus = false;
	    inputs.residentAdmissionRetry = false;
	    inputs.capacityAllocation = false;
	    inputs.residentGrowth = false;
	    inputs.pointTriangleRecovery = false;
	    inputs.structuralFrontier = false;
	    break;
	case Owner::PRESENTATION:
	    inputs.presentationReplay = false;
	    inputs.exactPresentation = false;
	    inputs.presentationBarrier = false;
	    inputs.capacityFrame = false;
	    inputs.pointAdmissionFrame = false;
	    inputs.pointCalibration = false;
	    inputs.capacityCalibration = false;
	    inputs.headroomProbe = false;
	    break;
	case Owner::HANDOFF:
	    inputs.handoff = false;
	    break;
	case Owner::COMPACTION:
	    inputs.compaction = false;
	    break;
	case Owner::CACHE_WRITE:
	    inputs.cacheWrite = false;
	    break;
	case Owner::NONE:
	    break;
    }
}

int
control_refinement_scenario(void)
{
    using Refinement = BObolLodControlRefinement;
    Refinement::Inputs inputs;
    inputs.presentationReplay = true;
    inputs.interaction = true;
    inputs.result = true;
    inputs.inventory = true;
    inputs.availability = true;
    inputs.capacityAllocation = true;
    inputs.handoff = true;
    inputs.compaction = true;
    inputs.cacheWrite = true;

    const auto observe = [&inputs]() {
	ConformanceState state;
	state.domains = BOBOL_MODEL_OBSERVES_CONTROL;
	const Refinement::Snapshot control = Refinement::evaluate(inputs);
	state.controlFactMask = Refinement::factMask(inputs);
	state.controlObligationMask = control.obligations;
	state.controlOwner = static_cast<int>(control.owner);
	state.presentationWitnessMask = control.owner ==
	    Refinement::Owner::PRESENTATION ?
	    Refinement::bit(Refinement::PresentationWitness::RENDER) : 0;
	state.terminal = !control.foregroundPending();
	state.viewReady = state.terminal;
	state.controlViolationMask = Refinement::validate(control,
	    state.terminal, state.viewReady, false,
	    state.presentationWitnessMask != 0, false, false) |
	    Refinement::validateProducers(inputs);
	return state;
    };
    Scenario scenario("ObolControlRefinement", "owner-discharge", observe);
    if (!scenario.checkInitial() ||
	observe().controlOwner !=
	static_cast<int>(Refinement::Owner::PRESENTATION))
	return 1;

    const Refinement::Owner expectedOwners[] = {
	Refinement::Owner::INTERACTION,
	Refinement::Owner::PUBLICATION,
	Refinement::Owner::INVENTORY,
	Refinement::Owner::AVAILABILITY,
	Refinement::Owner::PLANNING,
	Refinement::Owner::HANDOFF,
	Refinement::Owner::COMPACTION,
	Refinement::Owner::CACHE_WRITE,
	Refinement::Owner::NONE
    };
    for (const Refinement::Owner expected : expectedOwners) {
	if (!scenario.step("DischargeOwner",
		[&inputs]() { discharge_control_owner(inputs); },
		[expected](const ConformanceState &before,
		    const ConformanceState &after) {
		    return set_bit_count(after.controlFactMask) <
			set_bit_count(before.controlFactMask) &&
			after.controlOwner == static_cast<int>(expected);
		}))
	    return 1;
    }
    return observe().controlFactMask == 0 ? 0 : 1;
}

int
public_observation_scenario(void)
{
    constexpr size_t WorkerCount = 1;

    BObolLodService service;
    BObolViewController controller;
    controller.clearRenderRequest();
    controller.clearProgressiveWorkPending();
    const auto observe = [&controller]() {
	BObolLodConvergenceStatus status;
	controller.getLodConvergenceStatus(status);
	return bobol_model_conformance_state(status,
	    controller.getHostWorkSnapshot());
    };
    Scenario scenario("BObolPublicObservation", "state-vector", observe);
    if (!scenario.checkInitial() ||
	!scenario.step("RequestPresentationRender", [&controller]() {
	    controller.requestPresentationRender("model-conformance");
	}, [](const ConformanceState &, const ConformanceState &after) {
	    return after.hostWork.renderPending() &&
		after.hostWork.renderRevision != 0;
	}) || !scenario.step("RetirePresentationRender", [&controller]() {
	    controller.clearRenderRequest();
	}, [](const ConformanceState &, const ConformanceState &after) {
	    return !after.hostWork.renderPending();
	}))
	return 1;

    if (!service.start(WorkerCount, FALSE)) {
	bu_log("FAIL: public conformance observation could not start service\n");
	return 1;
    }
    const bool passed = scenario.step("AttachService", [&]() {
	controller.setLodService(&service);
    }, [](const ConformanceState &, const ConformanceState &after) {
	return after.activeGeneration == 0;
    }) && scenario.step("BeginGeneration", [&controller]() {
	(void)controller.beginLodGeneration();
    }, [](const ConformanceState &, const ConformanceState &after) {
	return after.activeGeneration != 0;
    }) && scenario.step("DetachService", [&controller]() {
	controller.setLodService(NULL);
    }, [](const ConformanceState &, const ConformanceState &after) {
	return after.activeGeneration == 0;
    });
    if (controller.getLodService())
	controller.setLodService(NULL);
    service.stop();
    return passed ? 0 : 1;
}

struct ResultFixture {
    BObolLodRequest request;
    BObolLodResultAuthenticationContext current;
    int providerStatus = BOBOL_LOD_PROVIDER_READY;
    BObolLodResultAuthentication authentication;

    ResultFixture()
    {
	this->request.occurrenceKey = "model-conformance-result";
	this->request.sourceRoutingId = SourceRoute;
	this->request.sourcePopulationEpoch = SourcePopulation;
	this->request.viewRevision = ViewEpoch;
	this->request.policyRevision = PolicyEpoch;
	this->synchronize();
	this->evaluate();
    }

    void synchronize(void)
    {
	this->current.sourceRoutingId = this->request.sourceRoutingId;
	this->current.sourcePopulationEpoch =
	    this->request.sourcePopulationEpoch;
	this->current.viewRevision = this->request.viewRevision;
	this->current.policyRevision = this->request.policyRevision;
    }

    void evaluate(void)
    {
	this->authentication =
	    BObolLodResultAuthenticationContract::evaluate(
		this->request, this->providerStatus, this->current);
    }
};

int
result_authentication_scenario(void)
{
    ResultFixture fixture;
    const auto observe = [&fixture]() {
	ConformanceState state;
	state.domains = BOBOL_MODEL_OBSERVES_RESULT;
	state.sourceRoutingId = fixture.request.sourceRoutingId.value();
	state.sourcePopulationEpoch =
	    fixture.request.sourcePopulationEpoch.value();
	state.viewEpoch = fixture.request.viewRevision.value();
	state.policyEpoch = fixture.request.policyRevision.value();
	state.resultMismatchMask = fixture.authentication.identityMismatchMask;
	state.resultDisposition = fixture.authentication.disposition;
	return state;
    };
    Scenario scenario("ObolResultAuthentication", "late-result", observe);
    if (!scenario.checkInitial() || !scenario.step("ApplyCurrentResult",
	[&fixture]() {
	    fixture.providerStatus = BOBOL_LOD_PROVIDER_READY;
	    fixture.evaluate();
	}, [](const ConformanceState &, const ConformanceState &after) {
	    return after.resultDisposition == BObolLodResultDisposition::PUBLISH;
	}))
	return 1;

    if (!scenario.step("ChangeDemand", [&fixture]() {
	fixture.current.viewRevision = ViewEpoch + 1;
	fixture.providerStatus = BOBOL_LOD_PROVIDER_ERROR;
	fixture.evaluate();
    }) || !scenario.step("RejectSupersededResult", [&fixture]() {
	fixture.evaluate();
    }, [](const ConformanceState &, const ConformanceState &after) {
	return after.resultDisposition == BObolLodResultDisposition::SUPERSEDE;
    }))
	return 1;
    if (!scenario.step("ReplacePopulation", [&fixture]() {
	fixture.synchronize();
	fixture.current.sourcePopulationEpoch = SourcePopulation + 1;
	fixture.evaluate();
    }) || !scenario.step("ChangeRoute", [&fixture]() {
	fixture.synchronize();
	fixture.current.sourcePopulationEpoch = SourcePopulation;
	fixture.current.sourceRoutingId = SourceRoute + 1;
	fixture.evaluate();
    }))
	return 1;

    if (!scenario.step("RecordCurrentFailure", [&fixture]() {
	fixture.synchronize();
	fixture.providerStatus = BOBOL_LOD_PROVIDER_ERROR;
	fixture.evaluate();
    }, [](const ConformanceState &, const ConformanceState &after) {
	return after.resultDisposition ==
	    BObolLodResultDisposition::RECORD_TERMINAL_FAILURE;
    }) || !scenario.step("RetryCurrentResult", [&fixture]() {
	fixture.synchronize();
	fixture.providerStatus = BOBOL_LOD_PROVIDER_RUNNING;
	fixture.evaluate();
    }, [](const ConformanceState &, const ConformanceState &after) {
	return after.resultDisposition ==
	    BObolLodResultDisposition::RETRY_CURRENT_DEMAND;
    }))
	return 1;
    return 0;
}

int
identity_exhaustion_scenario(void)
{
    constexpr uint64_t FirstSourceIdentity = 31;
    constexpr uint64_t FirstSourceSerial = 1;
    constexpr uint64_t SecondSourceIdentity = 37;
    constexpr uint64_t SecondSourceSerial = 2;

    BObolExactSampleIdentity identity;
    BObolExactSampleIdentity::Entries samples = {
	{FirstSourceIdentity, FirstSourceSerial},
	{SecondSourceIdentity, SecondSourceSerial}
    };
    uint64_t evidence = identity.intern(samples);
    uint64_t accepted = 0;
    bool halted = false;
    uint64_t counter = std::numeric_limits<uint64_t>::max() - 1;

    const auto observe = [&]() {
	ConformanceState state;
	state.domains = BOBOL_MODEL_OBSERVES_IDENTITY |
	    BOBOL_MODEL_OBSERVES_EVIDENCE;
	state.sourceIdentity = FirstSourceIdentity;
	state.evidenceIdentity = evidence;
	state.evidenceValid = evidence != 0;
	return state;
    };
    Scenario scenario("ObolIdentityExhaustion",
	"identity-exhaustion", observe);
    const uint64_t initialEvidence = evidence;
    if (!scenario.checkInitial() ||
	!scenario.step("AuthenticateEvidence",
	    [&]() { accepted = evidence; },
	    [&](const ConformanceState &, const ConformanceState &after) {
		return accepted == after.evidenceIdentity;
	    }) ||
	!scenario.step("CaptureEvidence", [&]() {
	    std::reverse(samples.begin(), samples.end());
	    evidence = identity.intern(samples);
	}, [initialEvidence](const ConformanceState &,
	    const ConformanceState &after) {
	    return after.evidenceIdentity == initialEvidence;
	}))
	return 1;

    if (!scenario.step("RequestMutation", [&]() {
	++samples.front().sourceSerial;
    }, [initialEvidence](const ConformanceState &,
	const ConformanceState &after) {
	return after.evidenceIdentity == initialEvidence;
    }))
	return 1;

    if (!scenario.step("Advance", [&]() {
	evidence = identity.intern(samples);
	uint64_t successor = 0;
	if (bobol_identity_successor(counter, successor))
	    counter = successor;
    }, [initialEvidence, &counter](const ConformanceState &,
	const ConformanceState &after) {
	return after.evidenceIdentity != initialEvidence &&
	    counter == std::numeric_limits<uint64_t>::max();
    }))
	return 1;
    const uint64_t mutatedEvidence = evidence;

    if (!scenario.step("CaptureEvidence", [&]() {
	identity.invalidate();
	evidence = identity.intern(samples);
    }, [mutatedEvidence](const ConformanceState &,
	const ConformanceState &after) {
	return after.evidenceIdentity != mutatedEvidence;
    }) || !scenario.step("FailStop", [&]() {
	uint64_t successor = counter;
	halted = !bobol_identity_successor(counter, successor);
    }, [&](const ConformanceState &, const ConformanceState &) {
	return halted;
    }))
	return 1;
    return 0;
}

int
cad_timing_evidence_scenario(void)
{
    constexpr float pointThreshold = 2.0f;
    constexpr size_t renderCost = 300;
    constexpr uint64_t uploadBytes = 200;
    constexpr uint64_t initialNanoseconds = 80;
    constexpr uint64_t replayNanoseconds = 70;
    constexpr uint64_t structuralNanoseconds = 60;

    BObolLodRendererPerformanceEvidence evidence;
    const auto observe = [&]() {
	ConformanceState state;
	state.domains = BOBOL_MODEL_OBSERVES_EVIDENCE;
	state.evidenceValue = evidence.cadPresentationNanosecondsAt(
	    pointThreshold, renderCost);
	state.evidenceValid = state.evidenceValue != 0;
	return state;
    };
    Scenario scenario("ObolCadTimingEvidence", "retained-replay", observe);
    if (!scenario.checkInitial() ||
	!scenario.step("CompleteInitialCadFrame", [&]() {
	    evidence.noteCadPresentation(initialNanoseconds, pointThreshold,
		renderCost, false, uploadBytes);
	}, [](const ConformanceState &, const ConformanceState &after) {
	    return !after.evidenceValid;
	}) || !scenario.step("CompleteReusableCadReplay", [&]() {
	    evidence.noteCadPresentation(replayNanoseconds, pointThreshold,
		renderCost, false, uploadBytes);
	}, [](const ConformanceState &, const ConformanceState &after) {
	    return after.evidenceValid &&
		after.evidenceValue == replayNanoseconds;
	}) || !scenario.step("RecordNonCadFrame", [&]() {
	    evidence.noteStructuralPresentation(structuralNanoseconds,
		pointThreshold, renderCost);
	}, [](const ConformanceState &before, const ConformanceState &after) {
	    return after.evidenceValid && after.evidenceValue ==
		before.evidenceValue;
	}))
	return 1;
    return 0;
}

int
identity_evidence_scenario(void)
{
    if (identity_exhaustion_scenario())
	return 1;
    return cad_timing_evidence_scenario();
}

int
semantic_presentation_scenario(void)
{
    constexpr uint64_t InitialMutationBoundaryNanoseconds = 10;
    constexpr uint64_t SupersededFrameStartNanoseconds = 11;
    constexpr uint64_t ReplacementMutationBoundaryNanoseconds = 12;
    constexpr uint64_t CurrentFrameStartNanoseconds = 13;

    BObolLodExactPresentationFrame frame;
    uint64_t frameStartedNanoseconds = 0;
    const auto observe = [&]() {
	ConformanceState state;
	state.domains = BOBOL_MODEL_OBSERVES_SEMANTIC_PRESENTATION;
	state.semanticPresentationPending = frame.pending();
	state.semanticRequestRequired = frame.requestPending();
	state.semanticFramePending = frame.framePending();
	return state;
    };
    Scenario scenario("ObolSemanticPresentation", "semantic-frame", observe);
    if (!scenario.checkInitial() || !scenario.step("MutateStyle", [&]() {
	frame.require(InitialMutationBoundaryNanoseconds);
    }, [](const ConformanceState &, const ConformanceState &after) {
	return after.semanticRequestRequired;
    }) || !scenario.step("QueueFrame", [&]() {
	frame.noteFrameRequested();
    }, [](const ConformanceState &, const ConformanceState &after) {
	return after.semanticFramePending;
    }) || !scenario.step("BeginFrame", [&]() {
	frameStartedNanoseconds = SupersededFrameStartNanoseconds;
    }) || !scenario.step("MutateStyle", [&]() {
	frame.require(ReplacementMutationBoundaryNanoseconds);
    }, [](const ConformanceState &, const ConformanceState &after) {
	return after.semanticRequestRequired;
    }) || !scenario.step("CompleteStaleFrame", [&]() {
	if (!frame.confirm(frameStartedNanoseconds))
	    frame.noteRequestRetired();
    }, [](const ConformanceState &, const ConformanceState &after) {
	return after.semanticRequestRequired;
    }) || !scenario.step("QueueFrame", [&]() {
	frame.noteFrameRequested();
    }) || !scenario.step("BeginFrame", [&]() {
	frameStartedNanoseconds = CurrentFrameStartNanoseconds;
    }) || !scenario.step("CompleteCurrentFrame", [&]() {
	(void)frame.confirm(frameStartedNanoseconds);
    }, [](const ConformanceState &, const ConformanceState &after) {
	return !after.semanticPresentationPending;
    }))
	return 1;
    return 0;
}

int
presentation_barrier_scenario(void)
{
    constexpr uint64_t FirstCompletedRenderSerial = 10;
    constexpr uint64_t FirstRequiredRenderSerial = 11;
    constexpr uint64_t FirstViewEpoch = 3;
    constexpr uint64_t FirstPolicyEpoch = 7;
    constexpr uint64_t ReplacementCompletedRenderSerial = 20;
    constexpr uint64_t ReplacementCompletedFrameSerial = 21;
    constexpr uint64_t ReplacementViewEpoch = 4;
    constexpr uint64_t MismatchedViewEpoch = 5;
    constexpr uint64_t ReplacementPolicyEpoch = 8;

    BObolLodPresentationTransaction transaction;
    BObolLodPresentationTransaction::Completion completion;
    const auto observe = [&]() {
	ConformanceState state;
	state.domains = BOBOL_MODEL_OBSERVES_TRANSACTION;
	state.transactionActive = transaction.barrierPending() ||
	    transaction.publicationPending();
	state.transactionBarrierPending = transaction.barrierPending();
	state.transactionIdentity = state.transactionActive ?
	    transaction.sequence() : 0;
	state.requiredRenderSerial = transaction.requiredRenderSerial();
	state.viewEpoch = transaction.viewEpoch();
	state.policyEpoch = transaction.policyEpoch();
	return state;
    };
    Scenario scenario("BObolLodPresentationTransaction", "frame-barrier",
	observe);
    if (!scenario.checkInitial() || !scenario.step("ArmBarrier", [&]() {
	(void)transaction.arm(
	    BObolLodPresentationTransaction::REASON_CUT_PRESENTATION,
	    FirstCompletedRenderSerial, FirstViewEpoch, FirstPolicyEpoch);
    }, [](const ConformanceState &, const ConformanceState &after) {
	return after.transactionBarrierPending &&
	    after.requiredRenderSerial == FirstRequiredRenderSerial;
    }) || !scenario.step("FrameBelowBarrier", [&]() {
	    completion = transaction.complete(FirstCompletedRenderSerial,
		FirstViewEpoch, FirstPolicyEpoch, true);
	}, [&](const ConformanceState &, const ConformanceState &after) {
	    return !completion.retired && !completion.stale &&
		after.transactionBarrierPending;
	}) || !scenario.step("IncompleteFrame", [&]() {
	    completion = transaction.complete(FirstRequiredRenderSerial,
		FirstViewEpoch, FirstPolicyEpoch, false);
	}, [&](const ConformanceState &, const ConformanceState &after) {
	    return !completion.retired && after.transactionBarrierPending;
	}) || !scenario.step("CompleteBarrier", [&]() {
	    completion = transaction.complete(FirstRequiredRenderSerial,
		FirstViewEpoch, FirstPolicyEpoch, true);
	}, [&](const ConformanceState &, const ConformanceState &after) {
	    return completion.retired && !after.transactionActive;
	}))
	return 1;

    if (!scenario.step("ArmReplacementBarrier", [&]() {
	(void)transaction.arm(
	    BObolLodPresentationTransaction::REASON_RESULT_PUBLICATION,
	    ReplacementCompletedRenderSerial, ReplacementViewEpoch,
	    ReplacementPolicyEpoch);
	}) || !scenario.step("RejectEpochMismatch", [&]() {
	completion = transaction.complete(ReplacementCompletedFrameSerial,
	    MismatchedViewEpoch, ReplacementPolicyEpoch, true);
    }, [&](const ConformanceState &, const ConformanceState &after) {
	return completion.stale && !after.transactionActive;
    }))
	return 1;
    return 0;
}

int
presentation_transaction_scenario(void)
{
    if (semantic_presentation_scenario())
	return 1;
    return presentation_barrier_scenario();
}

static constexpr uint64_t TraceIdentity = 1;
static constexpr uint64_t TraceInitialSerial = 1;
static constexpr uint64_t TraceBarrierSerial = 2;
static constexpr uint64_t TraceCompletionSerial = 3;
static constexpr uint64_t TraceBarrierHostRevision = 1;
static constexpr uint64_t TraceCompletionHostRevision = 2;
static constexpr uint64_t TraceUnreachedRenderSerial = 2;
static constexpr uint64_t TraceReplacementIdentity = 2;
static constexpr size_t TraceBarrierRecordIndex = 1;
static constexpr size_t TraceCompletionRecordIndex = 2;

BObolLodControlTraceState
trace_publication_state(void)
{
    BObolLodControlTraceState state;
    state.viewRevision = TraceIdentity;
    state.policyRevision = TraceIdentity;
    state.convergence.viewRevision = state.viewRevision;
    state.convergence.policyRevision = state.policyRevision;
    state.convergence.phase = BOBOL_LOD_CONVERGENCE_REFINING;
    state.convergence.outcome = BOBOL_LOD_PRESENTATION_ACTIVE;
    state.convergence.controlFactMask =
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Fact::RESULT);
    state.convergence.controlObligationMask =
	BOBOL_LOD_CONTROL_OBLIGATION_PUBLICATION;
    state.convergence.controlOwner = BOBOL_LOD_CONTROL_OWNER_PUBLICATION;
    state.convergence.terminal = FALSE;
    state.convergence.viewReady = FALSE;
    state.convergence.hasLodState = TRUE;
    return state;
}

BObolLodControlTraceState
trace_barrier_state(const BObolLodControlTraceState &publication)
{
    BObolLodControlTraceState state = publication;
    state.convergence.controlFactMask =
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Fact::PUBLICATION) |
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Fact::PRESENTATION_BARRIER);
    state.convergence.controlObligationMask =
	BOBOL_LOD_CONTROL_OBLIGATION_PUBLICATION |
	BOBOL_LOD_CONTROL_OBLIGATION_PRESENTATION;
    state.convergence.presentationTransactionSerial = TraceIdentity;
    state.convergence.presentationRequiredRenderSerial = TraceIdentity;
    state.convergence.publicationFramePending = TRUE;
    state.convergence.refinementFramePending = TRUE;
    state.hostWork.revision = TraceBarrierHostRevision;
    state.hostWork.renderRevision = TraceIdentity;
    state.hostWork.flags = BOBOL_HOST_WORK_RENDER;
    return state;
}

BObolLodControlTraceState
trace_terminal_state(const BObolLodControlTraceState &barrier)
{
    BObolLodControlTraceState state = barrier;
    state.convergence.phase = BOBOL_LOD_CONVERGENCE_IDLE;
    state.convergence.outcome = BOBOL_LOD_PRESENTATION_READY;
    state.convergence.controlFactMask = 0;
    state.convergence.controlObligationMask = 0;
    state.convergence.controlOwner = BOBOL_LOD_CONTROL_OWNER_NONE;
    state.convergence.presentationRequiredRenderSerial = 0;
    state.convergence.publicationFramePending = FALSE;
    state.convergence.refinementFramePending = FALSE;
    state.convergence.terminal = TRUE;
    state.convergence.viewReady = TRUE;
    state.hostWork.revision = TraceCompletionHostRevision;
    state.hostWork.renderRevision = 0;
    state.hostWork.flags = BOBOL_HOST_WORK_NONE;
    state.renderCompletionSerial = TraceIdentity;
    return state;
}

BObolLodControlTransitionRecord
trace_record(uint64_t serial, BObolLodControlTransitionEvent event,
    const BObolLodControlTraceState &before,
    const BObolLodControlTraceState &after)
{
    BObolLodControlTransitionRecord record;
    record.serial = serial;
    record.event = event;
    record.before = before;
    record.after = after;
    return record;
}

int
trace_refinement_scenario(void)
{
    const BObolLodControlTraceState publication = trace_publication_state();
    const BObolLodControlTraceState barrier = trace_barrier_state(publication);
    const BObolLodControlTraceState terminal = trace_terminal_state(barrier);
    const std::vector<BObolLodControlTransitionRecord> valid = {
	trace_record(TraceInitialSerial, BOBOL_LOD_CONTROL_TRANSITION_INITIAL,
	    publication, publication),
	trace_record(TraceBarrierSerial, BOBOL_LOD_CONTROL_TRANSITION_PUBLICATION,
	    publication, barrier),
	trace_record(TraceCompletionSerial,
	    BOBOL_LOD_CONTROL_TRANSITION_PUBLICATION,
	    barrier, terminal)
    };

    const BObolTraceRefinementResult accepted =
	BObolResultPresentationTraceChecker::check(valid, 0);
    if (!accepted.valid || accepted.statistics.publicationBarriers != 1 ||
	accepted.statistics.completedPublicationBarriers != 1)
	return 1;

    const BObolTraceRefinementResult dropped =
	BObolResultPresentationTraceChecker::check(valid, TraceIdentity);
    if (dropped.valid || dropped.failure.code !=
	    BObolTraceRefinementFailureCode::TRACE_LOSS)
	return 1;

    const auto rejects = [&valid](
	BObolTraceRefinementFailureCode expected,
	const std::function<void(
	    std::vector<BObolLodControlTransitionRecord> &)> &mutate) {
	std::vector<BObolLodControlTransitionRecord> candidate = valid;
	mutate(candidate);
	const BObolTraceRefinementResult checked =
	    BObolResultPresentationTraceChecker::check(candidate, 0);
	return !checked.valid && checked.failure.code == expected;
    };

    if (!rejects(BObolTraceRefinementFailureCode::UNMAPPED_EVENT,
	    [](std::vector<BObolLodControlTransitionRecord> &records) {
		records[TraceBarrierRecordIndex].event =
		    BOBOL_LOD_CONTROL_TRANSITION_UNNAMED;
	    }) ||
	!rejects(BObolTraceRefinementFailureCode::NONCONTIGUOUS_SERIAL,
	    [](std::vector<BObolLodControlTransitionRecord> &records) {
		++records[TraceBarrierRecordIndex].serial;
	    }) ||
	!rejects(BObolTraceRefinementFailureCode::DISCONTINUOUS_ENDPOINT,
	    [](std::vector<BObolLodControlTransitionRecord> &records) {
		++records[TraceCompletionRecordIndex].before.hostWork.revision;
	    }) ||
	!rejects(BObolTraceRefinementFailureCode::EVENT_OWNER_MISMATCH,
	    [](std::vector<BObolLodControlTransitionRecord> &records) {
		records[TraceBarrierRecordIndex].event =
		    BOBOL_LOD_CONTROL_TRANSITION_PRESENTATION;
	    }) ||
	!rejects(
	    BObolTraceRefinementFailureCode::BARRIER_RETIRED_WITHOUT_FRAME,
	    [](std::vector<BObolLodControlTransitionRecord> &records) {
		records[TraceCompletionRecordIndex].after.renderCompletionSerial = 0;
	    }) ||
	!rejects(BObolTraceRefinementFailureCode::INVALID_INITIAL,
	    [](std::vector<BObolLodControlTransitionRecord> &records) {
		records.front().serial = 0;
	    }) ||
	!rejects(BObolTraceRefinementFailureCode::INVALID_STATE,
	    [](std::vector<BObolLodControlTransitionRecord> &records) {
		const uint32_t violation =
		    BOBOL_LOD_CONTROL_VIOLATION_INVALID_OWNER;
		records.front().before.convergence.controlViolationMask = violation;
		records.front().after.convergence.controlViolationMask = violation;
	    }) ||
	!rejects(
	    BObolTraceRefinementFailureCode::CONTROL_PROJECTION_MISMATCH,
	    [](std::vector<BObolLodControlTransitionRecord> &records) {
		const uint32_t compaction = BObolLodControlRefinement::bit(
		    BObolLodControlRefinement::Fact::COMPACTION);
		records.front().before.convergence.controlFactMask |= compaction;
		records.front().after.convergence.controlFactMask |= compaction;
	    }) ||
	!rejects(BObolTraceRefinementFailureCode::OBSERVER_MISMATCH,
	    [](std::vector<BObolLodControlTransitionRecord> &records) {
		++records.front().before.viewRevision;
		++records.front().after.viewRevision;
	    }) ||
	!rejects(BObolTraceRefinementFailureCode::REVISION_REGRESSION,
	    [](std::vector<BObolLodControlTransitionRecord> &records) {
		BObolLodControlTraceState &after =
		    records[TraceCompletionRecordIndex].after;
		after.viewRevision = 0;
		after.convergence.viewRevision = 0;
	    }) ||
	!rejects(BObolTraceRefinementFailureCode::CLOCK_REGRESSION,
	    [](std::vector<BObolLodControlTransitionRecord> &records) {
		records[TraceCompletionRecordIndex].after.hostWork.revision = 0;
	    }) ||
	!rejects(BObolTraceRefinementFailureCode::BARRIER_IDENTITY_CHANGED,
	    [](std::vector<BObolLodControlTransitionRecord> &records) {
		records[TraceCompletionRecordIndex].after.convergence.
		    presentationRequiredRenderSerial =
		    TraceUnreachedRenderSerial;
	    }) ||
	!rejects(BObolTraceRefinementFailureCode::
		BARRIER_RETIRED_BEFORE_REQUIRED_FRAME,
	    [](std::vector<BObolLodControlTransitionRecord> &records) {
		records[TraceBarrierRecordIndex].after.convergence.
		    presentationRequiredRenderSerial =
		    TraceUnreachedRenderSerial;
		records[TraceCompletionRecordIndex].before.convergence.
		    presentationRequiredRenderSerial =
		    TraceUnreachedRenderSerial;
	    }) ||
	!rejects(BObolTraceRefinementFailureCode::
		PUBLICATION_RETIRED_WITHOUT_BARRIER,
	    [](std::vector<BObolLodControlTransitionRecord> &records) {
		records[TraceBarrierRecordIndex].after.convergence.
		    presentationRequiredRenderSerial = 0;
		records[TraceCompletionRecordIndex].before.convergence.
		    presentationRequiredRenderSerial = 0;
	    }) ||
	!rejects(BObolTraceRefinementFailureCode::
		PUBLICATION_BARRIER_SUPERSEDED,
	    [](std::vector<BObolLodControlTransitionRecord> &records) {
		records[TraceCompletionRecordIndex].after.convergence.
		    presentationTransactionSerial = TraceReplacementIdentity;
	    }))
	return 1;

    const std::vector<BObolLodControlTransitionRecord> noFlow = {
	valid.front()
    };
    const BObolTraceRefinementResult missing =
	BObolResultPresentationTraceChecker::check(noFlow, 0);
    if (missing.valid || missing.failure.code !=
	    BObolTraceRefinementFailureCode::MISSING_PUBLICATION_PRESENTATION_FLOW)
	return 1;

    const std::vector<BObolLodControlTransitionRecord> unfinished(
	valid.begin(), valid.begin() + TraceCompletionRecordIndex);
    const BObolTraceRefinementResult pending =
	BObolResultPresentationTraceChecker::check(unfinished, 0);
    return !pending.valid && pending.failure.code ==
	BObolTraceRefinementFailureCode::UNRETIRED_PUBLICATION_BARRIER ? 0 : 1;
}

struct NamedScenario {
    const char *name;
    int (*run)(void);
};

const NamedScenario Scenarios[] = {
    {"control-refinement", control_refinement_scenario},
    {"public-observation", public_observation_scenario},
    {"result-authentication", result_authentication_scenario},
    {"identity-evidence", identity_evidence_scenario},
    {"presentation-transaction", presentation_transaction_scenario},
    {"trace-refinement", trace_refinement_scenario}
};

} // namespace

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    bobol_init(NULL);
    if (argc > 2) {
	bu_log("Usage: %s [scenario]\n", argv[0]);
	return 1;
    }
    for (const NamedScenario &scenario : Scenarios) {
	if (argc == 2 && std::strcmp(argv[1], scenario.name) != 0)
	    continue;
	const int result = scenario.run();
	if (result)
	    return result;
	if (argc == 2)
	    return 0;
    }
    if (argc == 2) {
	bu_log("Unknown conformance scenario: %s\n", argv[1]);
	return 1;
    }
    return 0;
}
