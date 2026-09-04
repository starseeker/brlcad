/*              T E S T _ L O D _ C O O R D I N A T O R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "exact_sample_identity_private.h"
#include "identity_counter_private.h"
#include "lod_coordinator_private.h"
#include "lod_control_private.h"
#include "presentation_preparation_private.h"
#include "retained_allocation_private.h"
#include "bu/app.h"

#include <cstdio>
#include <atomic>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>

static_assert(sizeof(BObolLodViewEpoch) == sizeof(uint64_t),
    "typed LoD epochs must remain zero-overhead");
static_assert(sizeof(BObolLodPolicyEpoch) == sizeof(uint64_t),
    "typed LoD epochs must remain zero-overhead");
static_assert(std::is_trivially_copyable<BObolLodViewEpoch>::value,
    "view epochs must remain trivially copyable");
static_assert(std::is_trivially_copyable<BObolLodPolicyEpoch>::value,
    "policy epochs must remain trivially copyable");
static_assert(!std::is_same<BObolLodViewEpoch, BObolLodPolicyEpoch>::value,
    "view and policy epochs must be distinct types");
static_assert(!std::is_convertible<BObolLodViewEpoch, uint64_t>::value,
    "epoch serialization must be explicit");
static_assert(std::is_trivially_copyable<BObolLodResourceEvidence>::value,
    "resource policy must remain an allocation-free value");
static_assert(std::is_trivially_copyable<BObolLodHeadroomEvidence>::value,
    "headroom policy must remain an allocation-free value");
static_assert(std::is_trivially_copyable<BObolLodCoveragePolicy>::value,
    "coverage policy must remain an allocation-free value");
static_assert(std::is_trivially_copyable<BObolLodConvergencePolicy>::value,
    "convergence policy must remain an allocation-free value");
static_assert(std::is_trivially_copyable<BObolLodCapacityEvidence>::value,
    "capacity evidence must remain an allocation-free value");
static_assert(std::is_trivially_copyable<BObolLodAdmissionCursor>::value,
    "admission cursors must remain allocation-free values");
static_assert(std::is_trivially_copyable<
    BObolLodAdmissionRevisionStamp>::value,
    "admission revision stamps must remain allocation-free values");
static_assert(!std::is_default_constructible<
    BObolLodAdmissionRevisionStamp>::value,
    "admission revision stamps must require all six domains");
static_assert(std::is_constructible<BObolLodAdmissionRevisionStamp,
    BObolLodInventoryEpoch, BObolLodAvailabilityEpoch,
    BObolLodVisibilityEpoch, BObolLodViewEpoch, BObolLodPolicyEpoch,
    BObolLodCapacityEpoch>::value,
    "admission revision stamps must accept exactly the canonical domains");
static_assert(!std::is_constructible<BObolLodAdmissionRevisionStamp,
    BObolLodInventoryEpoch, BObolLodAvailabilityEpoch,
    BObolLodViewEpoch, BObolLodPolicyEpoch,
    BObolLodCapacityEpoch>::value,
    "a visibility companion may not be omitted from an admission stamp");

static int
test_identity_exhaustion_contract(void)
{
    uint64_t successor = 17;
    if (!bobol_identity_successor<uint64_t>(0, successor) ||
	successor != 1)
	return 1;

    const uint64_t maximum = std::numeric_limits<uint64_t>::max();
    if (!bobol_identity_successor(maximum - 1, successor) ||
	successor != maximum)
	return 1;
    successor = 17;
    if (bobol_identity_successor(maximum, successor) || successor != 17)
	return 1;

    BObolLodViewEpoch epoch(maximum - 1);
    if (!epoch.canAdvance())
	return 1;
    epoch.advance();
    if (epoch.value() != maximum || epoch.canAdvance())
	return 1;

    uint64_t nextSequence = 0;
    if (bobol_identity_take(nextSequence) != 0 || nextSequence != 1)
	return 1;
    uint64_t nextNonzeroIdentity = 1;
    if (bobol_nonzero_identity_take(nextNonzeroIdentity) != 1 ||
	nextNonzeroIdentity != 2)
	return 1;

    std::atomic<uint64_t> atomicIdentity(41);
    if (bobol_atomic_identity_advance(atomicIdentity) != 42 ||
	atomicIdentity.load(std::memory_order_relaxed) != 42)
	return 1;
    if (bobol_atomic_nonzero_identity_take(atomicIdentity) != 42 ||
	atomicIdentity.load(std::memory_order_relaxed) != 43)
	return 1;

    uint32_t diagnostic = std::numeric_limits<uint32_t>::max();
    bobol_saturating_counter_advance(diagnostic);
    if (diagnostic != std::numeric_limits<uint32_t>::max())
	return 1;
    return 0;
}

static int
test_exact_sample_identity(void)
{
    BObolExactSampleIdentity identity;
    const BObolExactSampleIdentity::Entries first = {
	{17, 2}, {9, 4}
    };
    const uint64_t firstToken = identity.intern(first);
    if (!firstToken || identity.intern({{9, 4}, {17, 2}}) != firstToken)
	return 1;

    const uint64_t changedToken = identity.intern({{9, 5}, {17, 2}});
    if (changedToken == firstToken)
	return 1;
    const uint64_t returnedToken = identity.intern(first);
    if (returnedToken == firstToken || returnedToken == changedToken)
	return 1;

    identity.invalidate();
    if (identity.intern(first) == returnedToken)
	return 1;

    BObolExactSampleIdentity boundedIdentity(2);
    if (boundedIdentity.intern({{1, 1}}) != 1 ||
	boundedIdentity.intern({{1, 2}}) != 2 ||
	boundedIdentity.intern({{1, 2}}) != 2)
	return 1;
    return 0;
}

static BObolLodAdmissionRevisionStamp
admission_revision_stamp(uint64_t inventory = 1,
    uint64_t availability = 1, uint64_t visibility = 1,
    uint64_t view = 1, uint64_t policy = 1, uint64_t capacity = 1)
{
    return BObolLodAdmissionRevisionStamp(
	BObolLodInventoryEpoch(inventory),
	BObolLodAvailabilityEpoch(availability),
	BObolLodVisibilityEpoch(visibility), BObolLodViewEpoch(view),
	BObolLodPolicyEpoch(policy), BObolLodCapacityEpoch(capacity));
}

using BObolLodCapacityView = decltype(
    std::declval<BObolLodAdmissionEvidence &>().capacity());
static_assert(std::is_const<
    typename std::remove_reference<BObolLodCapacityView>::type>::value,
    "controllers must receive only const admission-evidence views");
static_assert(std::is_trivially_copyable<BObolLodViewDemandPolicy>::value,
    "view-demand policy must remain an allocation-free value");
static_assert(std::is_trivially_copyable<BObolLodInteractionSession>::value,
    "interaction sessions must remain allocation-free values");
static_assert(std::is_trivially_copyable<
    BObolLodInteractionStartCertificate>::value,
    "interaction start certificates must remain allocation-free values");
static_assert(std::is_trivially_copyable<BObolLodPresentationPolicy>::value,
    "presentation policy must remain an allocation-free value");
static_assert(std::is_trivially_copyable<
    BObolLodViewQualityHistory>::value,
    "view quality history must remain an allocation-free value");
static_assert(std::is_trivially_copyable<BObolLodQualityPolicy>::value,
    "quality policy must remain an allocation-free value");
static_assert(std::is_trivially_copyable<BObolLodDeliveryPolicy>::value,
    "delivery policy must remain an allocation-free value");
static_assert(std::is_trivially_copyable<
    BObolLodPresentationTransaction>::value,
    "presentation transactions must remain allocation-free values");
static_assert(!std::is_copy_constructible<BObolLodAvailabilityLedger>::value,
    "availability ledger must have one controller owner");
static_assert(std::is_trivially_copyable<BObolLodCompactionPolicy>::value,
    "compaction policy must remain an allocation-free value");
static_assert(std::is_trivially_copyable<BObolLodPointQualityPhase>::value,
    "point-quality phase ownership must remain an allocation-free value");
static_assert(std::is_trivially_copyable<
    BObolLodControlRefinement::Inputs>::value,
    "control refinement inputs must remain allocation-free values");
static_assert(std::is_trivially_copyable<
    BObolLodControlRefinement::Snapshot>::value,
    "control refinement snapshots must remain allocation-free values");
static_assert(std::is_trivially_copyable<
    BObolLodControlRefinement::PresentationProgress>::value,
    "presentation progress must remain an allocation-free value");
static_assert(std::is_trivially_copyable<
    BObolLodAdmissionPlanner::PointCalibrationProducerInputs>::value,
    "point-calibration producer inputs must remain allocation-free values");
static_assert(std::is_trivially_copyable<BObolLodSubmissionPass>::value,
    "submission-pass ownership must remain an allocation-free value");
static_assert(std::is_trivially_copyable<
    BObolLodDeadlineSafePresentation>::value,
    "deadline-safe presentation evidence must remain allocation-free");
static_assert(std::is_trivially_copyable<
    BObolLodRetainedQualityEvidence>::value,
    "retained quality evidence must remain allocation-free");
static_assert(std::is_trivially_copyable<
    BObolLodRendererPerformanceEvidence>::value,
    "renderer performance evidence must remain allocation-free");
static_assert(std::is_trivially_copyable<BObolLodSubmissionIntent>::value,
    "submission intent must remain an allocation-free value");
static_assert(std::is_trivially_copyable<
    BObolLodRetainedViewContinuity>::value,
    "retained-view continuity must remain an allocation-free value");
static_assert(std::is_trivially_copyable<
    BObolLodRetainedPassAnnotations>::value,
    "retained-pass annotations must remain an allocation-free value");
static_assert(std::is_trivially_copyable<
    BObolLodInterruptedPresentationReplay>::value,
    "presentation replay ownership must remain allocation-free");
static_assert(std::is_trivially_copyable<
    BObolLodExactPresentationFrame>::value,
    "exact-presentation ownership must remain allocation-free");
static_assert(std::is_trivially_copyable<BObolLodPointAdmissionFrame>::value,
    "point-admission frame ownership must remain allocation-free");
static_assert(std::is_trivially_copyable<BObolLodStructuralRepair>::value,
    "structural repair ownership must remain an allocation-free value");
static_assert(std::is_trivially_copyable<
    BObolLodPlanningObligations>::value,
    "planning obligations must remain an allocation-free value");
static_assert(std::is_trivially_copyable<
    BObolLodDiscreteTrialPermit>::value,
    "discrete trial permits must remain allocation-free values");
static_assert(std::is_trivially_copyable<
    BObolLodTerminalQualityReducer::Inputs>::value,
    "terminal-quality inputs must remain allocation-free values");
static_assert(std::is_trivially_copyable<
    BObolLodTerminalQualityReducer::Decision>::value,
    "terminal-quality decisions must remain allocation-free values");

static int
test_submission_pass(void)
{
    using Pass = BObolLodSubmissionPass;
    Pass pass;
    if (pass.state() != Pass::State::IDLE || pass.active() ||
	pass.rescanPending()) {
	std::fprintf(stderr, "FAIL: initial submission-pass state\n");
	return 1;
    }

    pass.activate();
    pass.setRescanPending(true);
    if (pass.state() != Pass::State::ACTIVE_RESCAN || !pass.active() ||
	!pass.rescanPending()) {
	std::fprintf(stderr, "FAIL: active submission rescan\n");
	return 1;
    }

    pass.deactivate();
    if (pass.state() != Pass::State::IDLE_RESCAN || pass.active() ||
	!pass.rescanPending()) {
	std::fprintf(stderr, "FAIL: paused submission rescan\n");
	return 1;
    }

    if (!pass.beginPendingRescan() || pass.beginPendingRescan()) {
	std::fprintf(stderr, "FAIL: pending submission rescan transition\n");
	return 1;
    }
    if (pass.state() != Pass::State::ACTIVE || !pass.active() ||
	pass.rescanPending()) {
	std::fprintf(stderr, "FAIL: consumed submission rescan\n");
	return 1;
    }

    pass.requestRescan();
    if (pass.state() != Pass::State::ACTIVE_RESCAN || !pass.active() ||
	!pass.rescanPending()) {
	std::fprintf(stderr, "FAIL: set submission rescan\n");
	return 1;
    }
    pass.clearRescan();
    if (pass.state() != Pass::State::ACTIVE || !pass.active() ||
	pass.rescanPending()) {
	std::fprintf(stderr, "FAIL: clear submission rescan through setter\n");
	return 1;
    }

    pass.setActive(false);
    pass.requestRescan();
    pass.beginFresh();
    if (pass.state() != Pass::State::ACTIVE || !pass.active() ||
	pass.rescanPending()) {
	std::fprintf(stderr,
	    "FAIL: fresh submission pass inherited predecessor state\n");
	return 1;
    }
    pass.requestRescan();
    pass.retire();
    if (pass.state() != Pass::State::IDLE || pass.active() ||
	pass.rescanPending()) {
	std::fprintf(stderr,
	    "FAIL: retired submission pass preserved transaction debt\n");
	return 1;
    }
    pass.beginFresh(false);
    if (pass.state() != Pass::State::IDLE || pass.active() ||
	pass.rescanPending()) {
	std::fprintf(stderr,
	    "FAIL: empty fresh submission pass did not remain idle\n");
	return 1;
    }
    return 0;
}

static int
test_deadline_safe_presentation(void)
{
    BObolLodDeadlineSafePresentation presentation;
    if (presentation.ceilingFor(3, 7) != -1) {
	std::fprintf(stderr, "FAIL: empty deadline-safe presentation\n");
	return 1;
    }

    presentation.remember(5, 3, 7);
    presentation.remember(2, 3, 7);
    if (presentation.ceilingFor(3, 7) != 5 ||
	presentation.ceilingFor(4, 7) != -1 ||
	presentation.ceilingFor(3, 8) != -1) {
	std::fprintf(stderr, "FAIL: deadline-safe presentation keying\n");
	return 1;
    }

    presentation.remember(4, 4, 8);
    if (presentation.ceilingFor(3, 7) != -1 ||
	presentation.ceilingFor(4, 8) != 4) {
	std::fprintf(stderr, "FAIL: deadline-safe presentation replacement\n");
	return 1;
    }

    presentation.remember(-1, 5, 9);
    if (presentation.ceilingFor(4, 8) != 4) {
	std::fprintf(stderr, "FAIL: invalid presentation replaced evidence\n");
	return 1;
    }

    presentation.reset();
    if (presentation.ceilingFor(4, 8) != -1) {
	std::fprintf(stderr, "FAIL: deadline-safe presentation reset\n");
	return 1;
    }
    return 0;
}

static int
test_retained_quality_evidence(void)
{
    static constexpr double TEST_EPSILON = 1.0e-12;
    BObolLodRetainedQualityEvidence evidence;
    if (std::isfinite(evidence.maximumNormalizedError()) ||
	std::isfinite(evidence.maximumProjectedErrorPixelsFor(2, 3, 1.0f))) {
	std::fprintf(stderr, "FAIL: empty retained quality evidence\n");
	return 1;
    }

    evidence.remember(0.25, 0.75, 2, 3, 1.0f);
    if (std::fabs(evidence.maximumNormalizedError() - 0.25) > TEST_EPSILON ||
	std::fabs(evidence.maximumProjectedErrorPixelsFor(2, 3, 1.0f) - 0.75) >
	    TEST_EPSILON ||
	std::isfinite(evidence.maximumProjectedErrorPixelsFor(4, 3, 1.0f)) ||
	std::isfinite(evidence.maximumProjectedErrorPixelsFor(2, 4, 1.0f)) ||
	std::isfinite(evidence.maximumProjectedErrorPixelsFor(2, 3, 1.1f)) ||
	std::isfinite(evidence.maximumProjectedErrorPixelsFor(2, 3,
		std::numeric_limits<float>::quiet_NaN()))) {
	std::fprintf(stderr, "FAIL: retained quality evidence keying\n");
	return 1;
    }

    evidence.reset();
    if (std::isfinite(evidence.maximumNormalizedError()) ||
	std::isfinite(evidence.maximumProjectedErrorPixelsFor(2, 3, 1.0f))) {
	std::fprintf(stderr, "FAIL: retained quality evidence reset\n");
	return 1;
    }
    return 0;
}

static int
test_renderer_performance_evidence(void)
{
    BObolLodRendererPerformanceEvidence evidence;
    if (evidence.lastGpuTimeNanoseconds() != 0 ||
	evidence.reusableCadPresentationAt(2.0f)) {
	std::fprintf(stderr, "FAIL: initial renderer performance evidence\n");
	return 1;
    }
    if (!BObolLodRendererPerformanceEvidence::pointThresholdMatches(
	    2.0f, 2.005f) ||
	BObolLodRendererPerformanceEvidence::pointThresholdMatches(
	    2.0f, 2.02f) ||
	BObolLodRendererPerformanceEvidence::pointThresholdMatches(
	    std::numeric_limits<float>::quiet_NaN(), 2.0f)) {
	std::fprintf(stderr, "FAIL: renderer point-threshold identity\n");
	return 1;
    }

    if (!evidence.acceptGpuMeasurement(2, 100) ||
	evidence.lastGpuTimeNanoseconds() != 100 ||
	evidence.acceptGpuMeasurement(2, 200) ||
	evidence.lastGpuTimeNanoseconds() != 100) {
	std::fprintf(stderr, "FAIL: renderer GPU sample identity\n");
	return 1;
    }

    evidence.noteCadPresentation(80, 2.0f, 300, false, 200);
    if (evidence.reusableCadPresentationAt(2.0f)) {
	std::fprintf(stderr, "FAIL: first retained upload marked reusable\n");
	return 1;
    }
    evidence.noteCadPresentation(70, 2.0f, 300, false, 200);
    if (!evidence.reusableCadPresentationAt(2.0f) ||
	evidence.cadPresentationNanosecondsAt(2.0f) != 70 ||
	evidence.cadPresentationNanosecondsAt(2.0f, 300) != 70 ||
	evidence.cadPresentationNanosecondsAt(2.0f, 301) != 0 ||
	evidence.cadPresentationNanosecondsAt(4.0f) != 0) {
	std::fprintf(stderr, "FAIL: unchanged retained upload not reusable\n");
	return 1;
    }
    evidence.noteCadPresentation(60, 4.0f, 400, true,
	std::optional<uint64_t>());
    if (!evidence.reusableCadPresentationAt(4.0f)) {
	std::fprintf(stderr, "FAIL: prepared replay not reusable\n");
	return 1;
    }
    evidence.noteCadPresentation(50, 4.0f, std::optional<size_t>(), false,
	std::optional<uint64_t>());
    if (evidence.reusableCadPresentationAt(4.0f)) {
	std::fprintf(stderr, "FAIL: unprepared unsampled frame reusable\n");
	return 1;
    }

    evidence.noteStructuralPresentation(90, 4.0f, 500);
    if (evidence.structuralPresentationNanosecondsAt(4.0f) != 90 ||
	evidence.structuralPresentationNanosecondsAt(4.0f, 500) != 90 ||
	evidence.structuralPresentationNanosecondsAt(4.0f, 501) != 0 ||
	evidence.structuralPresentationNanosecondsAt(2.0f) != 0 ||
	evidence.reusableCadPresentationAt(4.0f)) {
	std::fprintf(stderr, "FAIL: structural presentation evidence isolation\n");
	return 1;
    }

    evidence.reset();
    if (evidence.lastGpuTimeNanoseconds() != 0 ||
	evidence.reusableCadPresentationAt(4.0f) ||
	evidence.cadPresentationNanosecondsAt(4.0f) != 0 ||
	evidence.structuralPresentationNanosecondsAt(4.0f) != 0) {
	std::fprintf(stderr, "FAIL: renderer performance evidence reset\n");
	return 1;
    }
    return 0;
}

static int
test_submission_intent(void)
{
    using Intent = BObolLodSubmissionIntent;
    Intent intent;
    if (intent.mode() != Intent::Mode::ORDINARY ||
	intent.retainedAdmission() || !intent.refreshMissing() ||
	intent.resetExisting())
	return 1;

    intent.configure(false, true);
    intent.setRetainedAdmission(true);
    if (intent.mode() != Intent::Mode::RETAINED_ADMISSION ||
	!intent.retainedAdmission() || intent.refreshMissing() ||
	!intent.resetExisting() ||
	!intent.retainedAdmissionForPass(false, true) ||
	intent.retainedAdmissionForPass(false, false)) {
	std::fprintf(stderr, "FAIL: configured submission intent\n");
	return 1;
    }

    intent.reset();
    if (intent.mode() != Intent::Mode::ORDINARY ||
	intent.retainedAdmission() || !intent.refreshMissing() ||
	intent.resetExisting() ||
	!intent.retainedAdmissionForPass(true, false) ||
	intent.retainedAdmissionForPass(false, true)) {
	std::fprintf(stderr, "FAIL: reset submission intent\n");
	return 1;
    }
    return 0;
}

static int
test_retained_view_continuity(void)
{
    BObolLodRetainedViewContinuity continuity;
    if (continuity.retainOccurrenceCuts() ||
	continuity.startCapacityAtStatic() ||
	continuity.visibilityCensusDeferred())
	return 1;

    continuity.beginQuiet(true, false, true, true);
    continuity.deferVisibilityCensus();
    if (!continuity.retainOccurrenceCuts() ||
	!continuity.startCapacityAtStatic() ||
	!continuity.visibilityCensusDeferred())
	return 1;

    continuity.completeVisibilityCensus();
    if (!continuity.retainOccurrenceCuts() ||
	!continuity.startCapacityAtStatic() ||
	continuity.visibilityCensusDeferred()) {
	std::fprintf(stderr,
	    "FAIL: visibility census completion discarded retained handoff\n");
	return 1;
    }

    /* Zoom and perspective pose require new demand cuts, but neither may
     * restart a ready retained presentation at the preferred cadence. */
    continuity.beginQuiet(true, true, true, true);
    if (continuity.retainOccurrenceCuts() ||
	!continuity.startCapacityAtStatic()) {
	std::fprintf(stderr,
	    "FAIL: zoom handoff did not preserve its static capacity start\n");
	return 1;
    }
    continuity.beginQuiet(false, false, true, true);
    if (continuity.retainOccurrenceCuts() ||
	!continuity.startCapacityAtStatic()) {
	std::fprintf(stderr,
	    "FAIL: perspective handoff did not preserve its static capacity start\n");
	return 1;
    }

    continuity.beginQuiet(true, false, false, true);
    if (continuity.retainOccurrenceCuts() ||
	continuity.startCapacityAtStatic()) {
	std::fprintf(stderr,
	    "FAIL: non-ready view manufactured retained continuity\n");
	return 1;
    }
    continuity.beginQuiet(true, false, true, false);
    if (continuity.retainOccurrenceCuts() ||
	continuity.startCapacityAtStatic()) {
	std::fprintf(stderr,
	    "FAIL: empty retained view manufactured continuity\n");
	return 1;
    }

    continuity.beginQuiet(true, false, true, true);
    continuity.clearHandoff();
    if (continuity.retainOccurrenceCuts() ||
	continuity.startCapacityAtStatic()) {
	std::fprintf(stderr, "FAIL: retained handoff did not clear\n");
	return 1;
    }

    continuity.reset();
    return continuity.retainOccurrenceCuts() ||
	continuity.startCapacityAtStatic() ||
	continuity.visibilityCensusDeferred() ? 1 : 0;
}

static int
test_retained_pass_annotations(void)
{
    BObolLodRetainedPassAnnotations pass;
    if (pass.refinementPending() || pass.residencyPending() ||
	pass.cutAdvanced() || pass.budgetBlocked() || pass.admittedWork() ||
	pass.missingMeshBudgetBlockedCount() != 0)
	return 1;

    pass.noteRefinementPending();
    pass.noteResidencyPending();
    pass.noteCutAdvanced();
    pass.noteBudgetBlocked();
    pass.noteAdmittedWork();
    pass.addMissingMeshBudgetBlocked(7);
    pass.addMissingMeshBudgetBlocked(5);
    if (!pass.refinementPending() || !pass.residencyPending() ||
	!pass.cutAdvanced() || !pass.budgetBlocked() || !pass.admittedWork() ||
	pass.missingMeshBudgetBlockedCount() != 12)
	return 1;

    pass.retireRefinement();
    if (pass.refinementPending() || pass.cutAdvanced() ||
	pass.budgetBlocked() || !pass.residencyPending() ||
	!pass.admittedWork() || pass.missingMeshBudgetBlockedCount() != 12) {
	std::fprintf(stderr,
	    "FAIL: retained-pass refinement retirement changed other facts\n");
	return 1;
    }

    pass.reset();
    if (pass.refinementPending() || pass.residencyPending() ||
	pass.cutAdvanced() || pass.budgetBlocked() || pass.admittedWork() ||
	pass.missingMeshBudgetBlockedCount() != 0) {
	std::fprintf(stderr, "FAIL: retained-pass reset was incomplete\n");
	return 1;
    }
    return 0;
}

static int
test_interrupted_presentation_replay(void)
{
    using Replay = BObolLodInterruptedPresentationReplay;
    Replay replay;
    if (replay.pending() || replay.state() != Replay::State::IDLE)
	return 1;
    replay.arm();
    if (!replay.pending() ||
	replay.state() != Replay::State::AWAITING_FRAME)
	return 1;
    replay.retire();
    if (replay.pending() || replay.state() != Replay::State::IDLE)
	return 1;
    replay.arm();
    replay.reset();
    return replay.pending() ? 1 : 0;
}

static int
test_exact_presentation_frame(void)
{
    using Frame = BObolLodExactPresentationFrame;
    if (!Frame::recoveryRequired(true, true, false, false) ||
	Frame::recoveryRequired(false, true, false, false) ||
	Frame::recoveryRequired(true, false, false, false) ||
	Frame::recoveryRequired(true, true, true, false) ||
	Frame::recoveryRequired(true, true, false, true)) {
	std::fprintf(stderr,
	    "FAIL: exact-frame recovery did not isolate an ownerless stale presentation\n");
	return 1;
    }
    Frame frame;
    if (frame.pending() || frame.requestPending() || frame.framePending() ||
	frame.state() != Frame::State::CURRENT)
	return 1;
    frame.require(10);
    if (!frame.pending() || !frame.requestPending() || frame.framePending() ||
	frame.state() != Frame::State::REQUEST_REQUIRED)
	return 1;
    frame.noteFrameRequested();
    if (!frame.pending() || frame.requestPending() || !frame.framePending() ||
	frame.state() != Frame::State::AWAITING_FRAME)
	return 1;
    frame.noteRequestRetired();
    if (!frame.requestPending() || frame.framePending())
	return 1;
    frame.noteFrameRequested();
    /* A newer semantic mutation supersedes the queued frame. */
    frame.require(20);
    if (!frame.requestPending() || frame.framePending())
	return 1;
    /* A frame which began before the newer mutation cannot acknowledge it. */
    if (frame.confirm(19) || !frame.requestPending())
	return 1;
    if (!frame.confirm(21))
	return 1;
    if (frame.pending() || frame.state() != Frame::State::CURRENT)
	return 1;
    frame.require(30);
    frame.reset();
    return frame.pending() ? 1 : 0;
}

static int
test_point_admission_frame(void)
{
    using Frame = BObolLodPointAdmissionFrame;
    Frame frame;
    if (frame.pending() || frame.state() != Frame::State::IDLE)
	return 1;
    frame.request();
    if (!frame.pending() || frame.state() != Frame::State::AWAITING_FRAME)
	return 1;
    frame.setPending(false);
    if (frame.pending())
	return 1;
    frame.setPending(true);
    frame.retire();
    return frame.pending() ? 1 : 0;
}

static int
test_structural_repair(void)
{
    using Repair = BObolLodStructuralRepair;
    using State = Repair::PointRelaxationState;
    Repair repair;

    if (repair.active() || repair.pointRelaxationPending() ||
	repair.pointRelaxationState() != State::INACTIVE)
	return 1;

    repair.begin(5);
    repair.reserveCoverageCost(42);
    repair.reserveCoverageCost(99);
    if (!repair.active() || repair.frontierCount() != 5 ||
	repair.coverageCostReservation() != 42 ||
	repair.pointRelaxationPending() ||
	repair.pointRelaxationRemainingRank() != 0) {
	std::fprintf(stderr,
	    "FAIL: ordinary structural repair did not retain one reservation\n");
	return 1;
    }

    if (!repair.beginPointRelaxation(3, 7, 8.0f) ||
	repair.beginPointRelaxation(8, 7, 8.0f)) {
	std::fprintf(stderr,
	    "FAIL: point relaxation accepted an invalid initial batch\n");
	return 1;
    }
    repair.reserveCoverageCost(21);
    repair.notePointRelaxationPresented();
    if (!repair.active() || repair.frontierCount() != 3 ||
	repair.coverageCostReservation() != 21 ||
	repair.pointRelaxationState() != State::ADMISSION_PENDING ||
	repair.pointRelaxationRemainingRank() != 7 ||
	repair.pointRelaxationPresentationPending()) {
	std::fprintf(stderr,
	    "FAIL: point relaxation skipped its admission phase\n");
	return 1;
    }

    repair.completePointRelaxationAdmission();
    repair.completePointRelaxationAdmission();
    if (repair.active() || !repair.pointRelaxationPending() ||
	repair.pointRelaxationState() != State::PRESENTATION_PENDING ||
	!repair.pointRelaxationPresentationPending() ||
	repair.coverageCostReservation() != 0 ||
	repair.pointRelaxationRemainingRank() != 7) {
	std::fprintf(stderr,
	    "FAIL: point relaxation did not enter its batch frame phase\n");
	return 1;
    }

    repair.notePointRelaxationPresented();
    repair.notePointRelaxationPresented();
    if (!repair.pointRelaxationPending() ||
	repair.pointRelaxationState() != State::FINALIZATION_PENDING ||
	repair.pointRelaxationPresentationPending() ||
	repair.pointRelaxationRemainingRank() != 7) {
	std::fprintf(stderr,
	    "FAIL: point relaxation did not enter finalization exactly once\n");
	return 1;
    }

    if (repair.retryPointRelaxationAdmission(0, 3) ||
	repair.retryPointRelaxationAdmission(4, 3) ||
	repair.retryPointRelaxationAdmission(3, 7) ||
	repair.retryPointRelaxationAdmission(3, 8) ||
	repair.pointRelaxationState() != State::FINALIZATION_PENDING ||
	repair.pointRelaxationRemainingRank() != 7) {
	std::fprintf(stderr,
	    "FAIL: point relaxation accepted a non-decreasing retry\n");
	return 1;
    }
    if (!repair.retryPointRelaxationAdmission(3, 4) || !repair.active() ||
	repair.frontierCount() != 3 ||
	repair.pointRelaxationRemainingRank() != 4 ||
	repair.pointRelaxationState() != State::ADMISSION_PENDING) {
	std::fprintf(stderr,
	    "FAIL: point relaxation rejected a strictly decreasing retry\n");
	return 1;
    }
    repair.completePointRelaxationAdmission();
    repair.notePointRelaxationPresented();
    if (!repair.retryPointRelaxationAdmission(1, 1) ||
	repair.pointRelaxationRemainingRank() != 1) {
	std::fprintf(stderr,
	    "FAIL: point relaxation did not preserve its decreasing rank\n");
	return 1;
    }
    repair.completePointRelaxationAdmission();
    repair.notePointRelaxationPresented();

    repair.cancelPointRelaxation();
    if (repair.active() || repair.pointRelaxationPending() ||
	repair.pointRelaxationState() != State::INACTIVE ||
	std::fabs(repair.pointRelaxationTarget()) > 1.0e-6f ||
	repair.pointRelaxationRemainingRank() != 0)
	return 1;

    if (Repair::pointRelaxationDomainChanged(false, false, false) ||
	!Repair::pointRelaxationDomainChanged(true, false, false) ||
	!Repair::pointRelaxationDomainChanged(false, true, false) ||
	!Repair::pointRelaxationDomainChanged(false, false, true)) {
	std::fprintf(stderr,
	    "FAIL: point relaxation domain invalidation is inconsistent\n");
	return 1;
    }

    if (!repair.beginPointRelaxation(3, 5, 4.0f))
	return 1;
    repair.reserveCoverageCost(12);
    repair.cancelPointRelaxation();
    if (!repair.active() || repair.frontierCount() != 3 ||
	repair.coverageCostReservation() != 12 ||
	repair.pointRelaxationPending() ||
	repair.pointRelaxationRemainingRank() != 0) {
	std::fprintf(stderr,
	    "FAIL: cancelling relaxation discarded useful structural work\n");
	return 1;
    }

    repair.reset();
    return repair.active() || repair.pointRelaxationPending() ? 1 : 0;
}

static int
test_planning_obligations(void)
{
    BObolLodPlanningObligations work;
    if (work.importanceCensusPending() ||
	work.residentAdmissionRetryPending() ||
	work.exactVisibilityReallocationPending())
	return 1;
    work.requestImportanceCensus();
    work.setResidentAdmissionRetry(true);
    work.requestExactVisibilityReallocation();
    if (!work.importanceCensusPending() ||
	!work.residentAdmissionRetryPending() ||
	!work.exactVisibilityReallocationPending())
	return 1;
    if (!work.exactVisibilityReallocationReady(
	    false, false, false, false, false) ||
	work.exactVisibilityReallocationReady(
	    true, false, false, false, false) ||
	work.exactVisibilityReallocationReady(
	    false, true, false, false, false) ||
	work.exactVisibilityReallocationReady(
	    false, false, true, false, false) ||
	work.exactVisibilityReallocationReady(
	    false, false, false, true, false) ||
	work.exactVisibilityReallocationReady(
	    false, false, false, false, true)) {
	std::fprintf(stderr,
	    "FAIL: exact visibility reallocation bypassed a prerequisite owner\n");
	return 1;
    }
    work.retireImportanceCensus();
    work.retireExactVisibilityReallocation();
    if (work.importanceCensusPending() ||
	!work.residentAdmissionRetryPending() ||
	work.exactVisibilityReallocationPending())
	return 1;
    work.reset();
    return work.residentAdmissionRetryPending() ||
	work.exactVisibilityReallocationPending() ? 1 : 0;
}

static int
test_discrete_trial_permit(void)
{
    BObolLodDiscreteTrialPermit permit;
    if (permit.available() ||
	permit.state() != BObolLodDiscreteTrialPermit::State::UNAVAILABLE)
	return 1;
    permit.grant();
    if (!permit.available())
	return 1;
    permit.consume();
    if (permit.available())
	return 1;
    permit.setAvailable(true);
    permit.reset();
    return permit.available() ? 1 : 0;
}

static int
test_control_refinement(void)
{
    using Refinement = BObolLodControlRefinement;

    Refinement::Inputs inputs;
    Refinement::Snapshot snapshot = Refinement::evaluate(inputs);
    if (snapshot.obligations != 0 ||
	snapshot.owner != Refinement::Owner::NONE ||
	snapshot.foregroundPending() || snapshot.controlPending() ||
	Refinement::validate(snapshot, true, true, false,
	    true, false, true) != 0) {
	std::fprintf(stderr, "FAIL: empty control refinement\n");
	return 1;
    }

    inputs.compaction = true;
    inputs.cacheWrite = true;
    snapshot = Refinement::evaluate(inputs);
    if (snapshot.foregroundPending() ||
	snapshot.owner != Refinement::Owner::COMPACTION ||
	Refinement::validate(snapshot, true, true, false,
	    true, false, true) != 0) {
	std::fprintf(stderr, "FAIL: background control refinement\n");
	return 1;
    }

    inputs = Refinement::Inputs();
    inputs.result = true;
    inputs.submission = true;
    inputs.presentationReplay = true;
    snapshot = Refinement::evaluate(inputs);
    if (snapshot.owner != Refinement::Owner::PRESENTATION ||
	!snapshot.calibrationPending() || !snapshot.controlPending() ||
	!snapshot.foregroundPending()) {
	std::fprintf(stderr, "FAIL: interrupted presentation ownership\n");
	return 1;
    }
    const uint32_t terminalViolation =
	Refinement::validate(snapshot, true, true, false, true, false, true);
    if (!(terminalViolation & Refinement::bit(
	    Refinement::Violation::TERMINAL_WITH_WORK))) {
	std::fprintf(stderr, "FAIL: terminal work refinement violation\n");
	return 1;
    }

    inputs.presentationReplay = false;
    snapshot = Refinement::evaluate(inputs);
    if (snapshot.owner != Refinement::Owner::PUBLICATION) {
	std::fprintf(stderr, "FAIL: publication ownership precedence\n");
	return 1;
    }

    inputs = Refinement::Inputs();
    inputs.importanceCensus = true;
    if (!(Refinement::validateProducers(inputs) & Refinement::bit(
	    Refinement::Violation::UNWITNESSED_PLANNING))) {
	std::fprintf(stderr, "FAIL: ownerless importance-census producer\n");
	return 1;
    }
    inputs.submission = true;
    if (!(Refinement::validateProducers(inputs) & Refinement::bit(
	    Refinement::Violation::UNWITNESSED_PLANNING))) {
	std::fprintf(stderr,
	    "FAIL: unrelated submission witnessed importance census\n");
	return 1;
    }

    inputs.inventory = true;
    if (Refinement::validateProducers(inputs) != 0) {
	std::fprintf(stderr,
	    "FAIL: inventory coverage did not witness importance census\n");
	return 1;
    }

    inputs.inventory = false;
    inputs.demandRefresh = true;
    if (Refinement::validateProducers(inputs) != 0) {
	std::fprintf(stderr, "FAIL: witnessed importance-census producer\n");
	return 1;
    }

    inputs = Refinement::Inputs();
    inputs.submissionDelta = true;
    if (!(Refinement::validateProducers(inputs) & Refinement::bit(
	    Refinement::Violation::UNWITNESSED_PLANNING))) {
	std::fprintf(stderr, "FAIL: orphaned selective-submission scope\n");
	return 1;
    }
    inputs.submission = true;
    if (Refinement::validateProducers(inputs) != 0) {
	std::fprintf(stderr, "FAIL: active selective-submission cursor\n");
	return 1;
    }

    snapshot = Refinement::Snapshot();
    const uint32_t readinessViolation =
	Refinement::validate(snapshot, false, true, true, true, false, true);
    if (!(readinessViolation & Refinement::bit(
	    Refinement::Violation::INVALID_READINESS))) {
	std::fprintf(stderr, "FAIL: invalid readiness refinement violation\n");
	return 1;
    }

    const uint32_t stalledViolation = Refinement::validateLiveness(
	snapshot, false, true, false);
    if (!(stalledViolation & Refinement::bit(
	    Refinement::Violation::NONTERMINAL_WITHOUT_PROGRESS)) ||
	Refinement::validateLiveness(snapshot, true, true, false) != 0 ||
	Refinement::validateLiveness(snapshot, false, false, false) != 0 ||
	Refinement::validateLiveness(snapshot, false, true, true) != 0) {
	std::fprintf(stderr,
	    "FAIL: ownerless nonterminal refinement violation\n");
	return 1;
    }

    struct FieldExpectation {
	bool Refinement::Inputs::*field;
	Refinement::Fact fact;
	Refinement::Work work;
    };
    const FieldExpectation fields[] = {
	{&Refinement::Inputs::interaction, Refinement::Fact::INTERACTION,
	    Refinement::Work::INTERACTION},
	{&Refinement::Inputs::inventory, Refinement::Fact::INVENTORY,
	    Refinement::Work::INVENTORY},
	{&Refinement::Inputs::availability, Refinement::Fact::AVAILABILITY,
	    Refinement::Work::AVAILABILITY},
	{&Refinement::Inputs::result, Refinement::Fact::RESULT,
	    Refinement::Work::PUBLICATION},
	{&Refinement::Inputs::publication, Refinement::Fact::PUBLICATION,
	    Refinement::Work::PUBLICATION},
	{&Refinement::Inputs::submission, Refinement::Fact::SUBMISSION,
	    Refinement::Work::PLANNING},
	{&Refinement::Inputs::demandRefresh,
	    Refinement::Fact::DEMAND_REFRESH, Refinement::Work::PLANNING},
	{&Refinement::Inputs::submissionRescan,
	    Refinement::Fact::SUBMISSION_RESCAN, Refinement::Work::PLANNING},
	{&Refinement::Inputs::submissionDelta,
	    Refinement::Fact::SUBMISSION_DELTA, Refinement::Work::PLANNING},
	{&Refinement::Inputs::qualityProbe, Refinement::Fact::QUALITY_PROBE,
	    Refinement::Work::PLANNING},
	{&Refinement::Inputs::retainedAllocation,
	    Refinement::Fact::RETAINED_ALLOCATION, Refinement::Work::PLANNING},
	{&Refinement::Inputs::retainedAllocationTransaction,
	    Refinement::Fact::RETAINED_ALLOCATION_TRANSACTION,
	    Refinement::Work::PLANNING},
	{&Refinement::Inputs::importanceCensus,
	    Refinement::Fact::IMPORTANCE_CENSUS, Refinement::Work::PLANNING},
	{&Refinement::Inputs::residentAdmissionRetry,
	    Refinement::Fact::RESIDENT_ADMISSION_RETRY,
	    Refinement::Work::PLANNING},
	{&Refinement::Inputs::capacityAllocation,
	    Refinement::Fact::CAPACITY_ALLOCATION, Refinement::Work::PLANNING},
	{&Refinement::Inputs::residentGrowth, Refinement::Fact::RESIDENT_GROWTH,
	    Refinement::Work::PLANNING},
	{&Refinement::Inputs::pointTriangleRecovery,
	    Refinement::Fact::POINT_TRIANGLE_RECOVERY,
	    Refinement::Work::PLANNING},
	{&Refinement::Inputs::structuralFrontier,
	    Refinement::Fact::STRUCTURAL_FRONTIER, Refinement::Work::PLANNING},
	{&Refinement::Inputs::presentationReplay,
	    Refinement::Fact::PRESENTATION_REPLAY,
	    Refinement::Work::PRESENTATION},
	{&Refinement::Inputs::exactPresentation,
	    Refinement::Fact::EXACT_PRESENTATION,
	    Refinement::Work::PRESENTATION},
	{&Refinement::Inputs::presentationBarrier,
	    Refinement::Fact::PRESENTATION_BARRIER,
	    Refinement::Work::PRESENTATION},
	{&Refinement::Inputs::capacityFrame, Refinement::Fact::CAPACITY_FRAME,
	    Refinement::Work::PRESENTATION},
	{&Refinement::Inputs::pointAdmissionFrame,
	    Refinement::Fact::POINT_ADMISSION_FRAME,
	    Refinement::Work::PRESENTATION},
	{&Refinement::Inputs::pointCalibration,
	    Refinement::Fact::POINT_CALIBRATION,
	    Refinement::Work::PRESENTATION},
	{&Refinement::Inputs::capacityCalibration,
	    Refinement::Fact::CAPACITY_CALIBRATION,
	    Refinement::Work::PRESENTATION},
	{&Refinement::Inputs::headroomProbe, Refinement::Fact::HEADROOM_PROBE,
	    Refinement::Work::PRESENTATION},
	{&Refinement::Inputs::handoff, Refinement::Fact::HANDOFF,
	    Refinement::Work::HANDOFF},
	{&Refinement::Inputs::compaction, Refinement::Fact::COMPACTION,
	    Refinement::Work::COMPACTION},
	{&Refinement::Inputs::cacheWrite, Refinement::Fact::CACHE_WRITE,
	    Refinement::Work::CACHE_WRITE}
    };
    const uint32_t backgroundMask =
	Refinement::bit(Refinement::Work::COMPACTION) |
	Refinement::bit(Refinement::Work::CACHE_WRITE);

    /* Every concrete latch must enter its declared abstract obligation.  The
	 * ownership precedence primarily depends on the nine abstract work
	 * classes.  Capacity allocation and interrupted replay are deliberate
	 * exceptions, so include their exact facts with representative planning
	 * and presentation aliases below. */
    for (const FieldExpectation &field : fields) {
	inputs = Refinement::Inputs();
	inputs.*field.field = true;
	snapshot = Refinement::evaluate(inputs);
	if (Refinement::factMask(inputs) != Refinement::bit(field.fact) ||
	    !snapshot.has(field.work) || snapshot.obligations == 0 ||
	    snapshot.owner == Refinement::Owner::NONE ||
	    !(snapshot.obligations &
	      Refinement::ownerObligation(snapshot.owner))) {
	    std::fprintf(stderr, "FAIL: concrete control-refinement field\n");
	    return 1;
	}
    }

    bool Refinement::Inputs::*const representatives[] = {
	&Refinement::Inputs::interaction,
	&Refinement::Inputs::inventory,
	&Refinement::Inputs::availability,
	&Refinement::Inputs::publication,
	&Refinement::Inputs::submission,
	&Refinement::Inputs::capacityAllocation,
	&Refinement::Inputs::presentationReplay,
	&Refinement::Inputs::exactPresentation,
	&Refinement::Inputs::presentationBarrier,
	&Refinement::Inputs::handoff,
	&Refinement::Inputs::compaction,
	&Refinement::Inputs::cacheWrite
    };
    constexpr size_t representativeCount =
	sizeof(representatives) / sizeof(representatives[0]);
    constexpr uint64_t combinationCount =
	uint64_t(1) << representativeCount;
    for (uint64_t combination = 0;
	combination < combinationCount; ++combination) {
	inputs = Refinement::Inputs();
	for (size_t field = 0; field < representativeCount; ++field)
	    inputs.*representatives[field] =
		(combination & (uint64_t(1) << field)) != 0;
	snapshot = Refinement::evaluate(inputs);
	const bool foreground = (snapshot.obligations & ~backgroundMask) != 0;
	if ((snapshot.obligations == 0) !=
		(snapshot.owner == Refinement::Owner::NONE) ||
	    (snapshot.owner != Refinement::Owner::NONE &&
	     !(snapshot.obligations &
	       Refinement::ownerObligation(snapshot.owner))) ||
	    snapshot.foregroundPending() != foreground ||
	    Refinement::validate(snapshot, false, false, false,
		true, false, true) != 0) {
	    std::fprintf(stderr,
		"FAIL: exhaustive control refinement combination %llu\n",
		static_cast<unsigned long long>(combination));
	    return 1;
	}
    }

    snapshot = Refinement::Snapshot();
    snapshot.obligations = Refinement::bit(Refinement::Work::PLANNING);
    const uint32_t ownerlessViolation =
	Refinement::validate(snapshot, false, false, false, true, false, true);
    snapshot.owner = Refinement::Owner::COMPACTION;
    const uint32_t wrongOwnerViolation =
	Refinement::validate(snapshot, false, false, false, true, false, true);
    if (!(ownerlessViolation & Refinement::bit(
	    Refinement::Violation::OWNERLESS_WORK)) ||
	!(wrongOwnerViolation & Refinement::bit(
	    Refinement::Violation::INVALID_OWNER))) {
	std::fprintf(stderr, "FAIL: invalid control-owner detection\n");
	return 1;
    }

    inputs = Refinement::Inputs();
    inputs.pointCalibration = true;
    snapshot = Refinement::evaluate(inputs);
    const uint32_t unwitnessedPresentation =
	Refinement::validate(snapshot, false, false, false,
	    false, false, true);
    if (!(unwitnessedPresentation & Refinement::bit(
	    Refinement::Violation::UNWITNESSED_PRESENTATION)) ||
	Refinement::validate(snapshot, false, false, false,
	    true, false, true) != 0) {
	std::fprintf(stderr,
	    "FAIL: point presentation without a progress witness was accepted\n");
	return 1;
    }

    Refinement::PresentationProgress presentationProgress;
    if (presentationProgress.witnessed() ||
	presentationProgress.witnessMask() != 0) {
	std::fprintf(stderr, "FAIL: empty presentation progress witness\n");
	return 1;
    }
    struct PresentationWitnessExpectation {
	bool Refinement::PresentationProgress::*field;
	Refinement::PresentationWitness witness;
    };
    const PresentationWitnessExpectation presentationWitnesses[] = {
	{&Refinement::PresentationProgress::renderPending,
	    Refinement::PresentationWitness::RENDER},
	{&Refinement::PresentationProgress::controllerPumpPending,
	    Refinement::PresentationWitness::CONTROLLER_PUMP},
	{&Refinement::PresentationProgress::finiteTimerPending,
	    Refinement::PresentationWitness::TIMER},
	{&Refinement::PresentationProgress::independentProducerPending,
	    Refinement::PresentationWitness::INDEPENDENT_PRODUCER},
	{&Refinement::PresentationProgress::claimedFramePending,
	    Refinement::PresentationWitness::CLAIMED_FRAME}
    };
    for (const PresentationWitnessExpectation &expectation :
	    presentationWitnesses) {
	presentationProgress = Refinement::PresentationProgress();
	presentationProgress.*expectation.field = true;
	if (!presentationProgress.witnessed() ||
	    presentationProgress.witnessMask() !=
		Refinement::bit(expectation.witness) ||
	    Refinement::validate(snapshot, false, false, false,
		presentationProgress.witnessed(), false, true) != 0) {
	    std::fprintf(stderr,
		"FAIL: presentation progress witness classification\n");
	    return 1;
	}
    }

    snapshot = Refinement::Snapshot();
    const uint32_t unwitnessedConstraint = Refinement::validate(
	snapshot, true, true, false, true, true, false);
    if (!(unwitnessedConstraint & Refinement::bit(
	    Refinement::Violation::UNWITNESSED_CONSTRAINT)) ||
	Refinement::validate(
	    snapshot, true, true, false, true, true, true) != 0) {
	std::fprintf(stderr,
	    "FAIL: constrained terminal result without evidence was accepted\n");
	return 1;
    }

    inputs = Refinement::Inputs();
    inputs.capacityAllocation = true;
    snapshot = Refinement::evaluate(inputs);
    if (snapshot.owner != Refinement::Owner::PLANNING ||
	!snapshot.has(Refinement::Work::PLANNING) ||
	snapshot.has(Refinement::Work::PRESENTATION) ||
	Refinement::validate(snapshot, false, false, false,
	    false, false, true) != 0) {
	std::fprintf(stderr,
	    "FAIL: capacity allocation was not classified as planning work\n");
	return 1;
    }
    return 0;
}

static int
test_epochs(void)
{
    BObolLodViewEpoch view(1);
    BObolLodPolicyEpoch policy(9);

    view.advance();
    if (view.value() != 2 || policy.value() != 9) {
	std::fprintf(stderr, "FAIL: independent epoch advance\n");
	return 1;
    }

    view.reset();
    policy.set(17);
    if (view.value() != 0 || policy.value() != 17) {
	std::fprintf(stderr, "FAIL: epoch reset/set\n");
	return 1;
    }

    return 0;
}

static int
test_interaction_session(void)
{
    using Session = BObolLodInteractionSession;
    Session session;
    constexpr int64_t debounce = 150;

    if (session.active() || session.gestureActive() ||
	session.phase() != Session::Phase::QUIET ||
	session.releaseCutFloorActive() ||
	session.settleAfterRenderSerial() != 0 ||
	session.quietReady(1000, debounce)) {
	std::fprintf(stderr, "FAIL: initial interaction session\n");
	return 1;
    }

    session.beginGesture(1000);
    if (!session.hasNewCeilingFeedback(7)) {
	std::fprintf(stderr, "FAIL: initial interaction feedback ownership\n");
	return 1;
    }
    session.noteCeilingFeedback(7);
    if (session.hasNewCeilingFeedback(7) ||
	!session.hasNewCeilingFeedback(8)) {
	std::fprintf(stderr, "FAIL: interaction feedback was consumed twice\n");
	return 1;
    }
    session.observeCameraChange(1100, 7);
    if (!session.active() || !session.gestureActive() ||
	session.phase() != Session::Phase::GESTURE ||
	session.releaseCutFloorActive() ||
	session.settleAfterRenderSerial() != 8 ||
	session.noteMotionFrameCompleted(7, 1200) ||
	!session.noteMotionFrameCompleted(8, 1300) ||
	session.settleAfterRenderSerial() != 0 ||
	session.releaseExpiredMotionFrame(2000, debounce) ||
	session.quietReady(2000, debounce)) {
	std::fprintf(stderr, "FAIL: gesture interaction lifecycle\n");
	return 1;
    }

    if (!session.endGesture(3000) || session.endGesture(3000) ||
	!session.active() || session.gestureActive() ||
	session.phase() != Session::Phase::DEBOUNCING ||
	!session.releaseCutFloorActive()) {
	std::fprintf(stderr, "FAIL: gesture release transition\n");
	return 1;
    }
    if (session.quietReady(3149, debounce) ||
	!session.quietReady(3150, debounce)) {
	std::fprintf(stderr, "FAIL: quiet debounce boundary\n");
	return 1;
    }
    session.finishQuiet();
    if (session.active() || session.gestureActive() ||
	session.phase() != Session::Phase::QUIET ||
	session.releaseCutFloorActive() ||
	session.settleAfterRenderSerial() != 0) {
	std::fprintf(stderr, "FAIL: quiet interaction completion\n");
	return 1;
    }

    session.observeCameraChange(5000, UINT64_MAX - 1);
    if (!session.active() || session.gestureActive() ||
	session.phase() != Session::Phase::DEBOUNCING ||
	session.releaseCutFloorActive() ||
	session.settleAfterRenderSerial() != UINT64_MAX ||
	session.releaseExpiredMotionFrame(5149, debounce) ||
	!session.releaseExpiredMotionFrame(5150, debounce) ||
	!session.quietReady(5150, debounce)) {
	std::fprintf(stderr, "FAIL: unbracketed motion-frame gate\n");
	return 1;
    }

    /* A completed unbracketed motion frame starts the same debounce without
     * retaining a frame gate.  A later deadline observation must therefore
     * be able to hand control to the quiet transition directly. */
    session.finishQuiet();
    session.observeCameraChange(6000, 12);
    if (!session.noteMotionFrameCompleted(13, 6050) ||
	session.releaseExpiredMotionFrame(6200, debounce) ||
	session.quietReady(6199, debounce) ||
	!session.quietReady(6200, debounce)) {
	std::fprintf(stderr,
	    "FAIL: completed motion-frame debounce transition\n");
	return 1;
    }

    session.reset();
    if (session.active() || session.releaseCutFloorActive() ||
	session.settleAfterRenderSerial() != 0 ||
	!session.hasNewCeilingFeedback(7)) {
	std::fprintf(stderr, "FAIL: interaction session reset\n");
	return 1;
    }

    session.noteCeilingFeedback(7);
    session.resetCeilingFeedback();
    if (!session.hasNewCeilingFeedback(7))
	return 1;

    return 0;
}

static int
test_revision_contract(void)
{
    using Contract = BObolLodRevisionContract;
    const BObolLodAdmissionRevisionStamp stamp =
	admission_revision_stamp(3, 5, 7, 11, 13, 17);
    BObolLodCapacityEvidence::Inputs planningInputs;
    const BObolLodAdmissionPlan plan = BObolLodAdmissionPlanner::plan(
	BObolLodAdmissionEvidence(), BObolLodAdmissionCursor(), stamp,
	planningInputs);
    const BObolLodAdmissionPlan administrativePlan;
    if (!plan.certifiedFor(stamp) || !plan.nextCursor.matches(stamp) ||
	administrativePlan.certifiedFor(stamp)) {
	std::fprintf(stderr, "FAIL: complete/administrative stamp identity\n");
	return 1;
    }

    const BObolLodAdmissionRevisionDomain domains[] = {
	BObolLodAdmissionRevisionDomain::INVENTORY,
	BObolLodAdmissionRevisionDomain::AVAILABILITY,
	BObolLodAdmissionRevisionDomain::VISIBILITY,
	BObolLodAdmissionRevisionDomain::VIEW,
	BObolLodAdmissionRevisionDomain::POLICY,
	BObolLodAdmissionRevisionDomain::CAPACITY
    };
    for (const BObolLodAdmissionRevisionDomain domain : domains) {
	const BObolLodAdmissionRevisionStamp next =
	    Contract::advance(stamp, domain);
	const unsigned int changed =
	    (next.inventory() != stamp.inventory() ? 1u : 0u) +
	    (next.availability() != stamp.availability() ? 1u : 0u) +
	    (next.visibility() != stamp.visibility() ? 1u : 0u) +
	    (next.view() != stamp.view() ? 1u : 0u) +
	    (next.policy() != stamp.policy() ? 1u : 0u) +
	    (next.capacity() != stamp.capacity() ? 1u : 0u);
	if (changed != 1u) {
	    std::fprintf(stderr,
		"FAIL: revision transition changed %u domains\n", changed);
	    return 1;
	}
	if (Contract::planDisposition(stamp, next) !=
		Contract::PlanDisposition::STALE ||
	    plan.certifiedFor(next) || plan.nextCursor.matches(next)) {
	    std::fprintf(stderr,
		"FAIL: domain %u did not stale every asynchronous identity\n",
		static_cast<unsigned int>(domain));
	    return 1;
	}
    }

    const BObolLodAdmissionRevisionStamp policy =
	Contract::setPolicy(stamp, 19);
    if (policy.policy().value() != 19 ||
	policy.inventory() != stamp.inventory() ||
	policy.availability() != stamp.availability() ||
	policy.visibility() != stamp.visibility() ||
	policy.view() != stamp.view() || policy.capacity() != stamp.capacity()) {
	std::fprintf(stderr, "FAIL: policy synchronization changed evidence\n");
	return 1;
    }

    const BObolLodAdmissionRevisionStamp empty =
	BObolLodAdmissionRevisionStamp::administrative();
    if (Contract::planDisposition(empty, stamp) !=
	    Contract::PlanDisposition::ADMINISTRATIVE ||
	Contract::planDisposition(stamp, stamp) !=
	    Contract::PlanDisposition::CURRENT ||
	Contract::planDisposition(policy, stamp) !=
	    Contract::PlanDisposition::STALE) {
	std::fprintf(stderr, "FAIL: certified plan disposition\n");
	return 1;
    }
    return 0;
}

static BObolLodResourceEvidence::Decision
apply_resource_observation(BObolLodResourceEvidence &evidence,
    bool cpuPressure, bool gpuPressure, bool recoveryEnabled)
{
    BObolLodAdmissionEvidence allEvidence(evidence);
    const BObolLodAdmissionPlan plan =
	BObolLodAdmissionPlanner::planResourceObservation(allEvidence, BObolLodAdmissionCursor(),
	    cpuPressure, gpuPressure, recoveryEnabled);
    evidence = plan.nextEvidence.resources();
    return plan.resourceDecision;
}

static int
test_resource_policy(void)
{
    BObolLodResourceEvidence policy;
    BObolLodResourceEvidence::Decision decision =
	apply_resource_observation(policy, false, false, true);
    if (decision.changed || decision.queueRecovery || policy.anyPressure() ||
	policy.recoveryPending() || policy.revision() != 0) {
	std::fprintf(stderr, "FAIL: initial resource observation\n");
	return 1;
    }

    decision = apply_resource_observation(policy, true, false, true);
    if (!decision.changed || !decision.queueRecovery ||
	!policy.cpuPressure() || policy.gpuPressure() ||
	!policy.recoveryPending() || policy.revision() != 1 ||
	policy.recoveryHandledRevision() != 0) {
	std::fprintf(stderr, "FAIL: CPU pressure edge\n");
	return 1;
    }

    decision = apply_resource_observation(policy, true, false, true);
    if (decision.changed || decision.queueRecovery || policy.revision() != 1) {
	std::fprintf(stderr, "FAIL: duplicate pressure edge\n");
	return 1;
    }

    policy.markRecoveryHandled();
    if (policy.recoveryPending() || !policy.anyPressure() ||
	policy.recoveryHandledRevision() != policy.revision()) {
	std::fprintf(stderr, "FAIL: handled persistent pressure\n");
	return 1;
    }

    decision = apply_resource_observation(policy, true, true, true);
    if (!decision.changed || !decision.queueRecovery ||
	!policy.recoveryPending() || policy.revision() != 2) {
	std::fprintf(stderr, "FAIL: GPU pressure edge\n");
	return 1;
    }

    decision = apply_resource_observation(policy, false, false, true);
    if (!decision.changed || !decision.pressureCleared ||
	decision.queueRecovery || policy.anyPressure() ||
	policy.recoveryPending() || policy.revision() != 3 ||
	policy.recoveryHandledRevision() != 3) {
	std::fprintf(stderr, "FAIL: resource pressure clear\n");
	return 1;
    }

    (void)apply_resource_observation(policy, true, false, true);
    policy.resetForServiceChange();
    if (policy.anyPressure() || policy.recoveryPending() ||
	policy.recoveryHandledRevision() != policy.revision()) {
	std::fprintf(stderr, "FAIL: resource service reset\n");
	return 1;
    }

    decision = apply_resource_observation(policy, true, false, false);
    if (!decision.changed || decision.queueRecovery ||
	policy.recoveryPending()) {
	std::fprintf(stderr, "FAIL: disabled recovery pressure edge\n");
	return 1;
    }
    decision = apply_resource_observation(policy, true, false, true);
    if (decision.changed || !decision.queueRecovery ||
	!policy.recoveryPending()) {
	std::fprintf(stderr, "FAIL: late resource recovery enable\n");
	return 1;
    }
    policy.markRecoveryHandled();
    decision = apply_resource_observation(policy, true, false, true);
    if (decision.changed || decision.queueRecovery ||
	policy.recoveryPending()) {
	std::fprintf(stderr, "FAIL: handled resource recovery repeated\n");
	return 1;
    }
    return 0;
}

static bool
apply_headroom_retry_plan(BObolLodHeadroomEvidence &evidence,
    BObolLodViewEpoch viewEpoch, BObolLodPolicyEpoch policyEpoch,
    size_t activePopulationCost, size_t currentBudget)
{
    BObolLodAdmissionEvidence allEvidence(evidence);
    const BObolLodAdmissionPlan plan =
	BObolLodAdmissionPlanner::planHeadroomRetry(allEvidence, BObolLodAdmissionCursor(), viewEpoch,
	    policyEpoch, activePopulationCost, currentBudget);
    evidence = plan.nextEvidence.headroom();
    return plan.headroomAccepted;
}

static bool
apply_headroom_deferral_plan(BObolLodHeadroomEvidence &evidence,
    BObolLodViewEpoch viewEpoch, BObolLodPolicyEpoch policyEpoch,
    size_t activePopulationCost)
{
    BObolLodAdmissionEvidence allEvidence(evidence);
    const BObolLodAdmissionPlan plan =
	BObolLodAdmissionPlanner::planHeadroomTransientDeferral(allEvidence, BObolLodAdmissionCursor(),
	    viewEpoch, policyEpoch, activePopulationCost);
    evidence = plan.nextEvidence.headroom();
    return plan.headroomAccepted;
}

static bool
apply_headroom_consumption_plan(BObolLodHeadroomEvidence &evidence,
    BObolLodViewEpoch viewEpoch, BObolLodPolicyEpoch policyEpoch,
    size_t currentBudget, size_t activePopulationCost,
    long double demonstratedBudget, uint64_t elapsedNanoseconds,
    uint64_t targetNanoseconds, bool reusableSample)
{
    BObolLodAdmissionEvidence allEvidence(evidence);
    const BObolLodAdmissionPlan plan =
	BObolLodAdmissionPlanner::planHeadroomConsumption(allEvidence, BObolLodAdmissionCursor(),
	    viewEpoch, policyEpoch, currentBudget, activePopulationCost,
	    demonstratedBudget, elapsedNanoseconds, targetNanoseconds,
	    reusableSample);
    evidence = plan.nextEvidence.headroom();
    return plan.headroomAccepted;
}

static void
apply_headroom_rejection(BObolLodHeadroomEvidence &evidence,
    BObolLodViewEpoch viewEpoch, BObolLodPolicyEpoch policyEpoch,
    size_t activePopulationCost, size_t attemptedBudget)
{
    BObolLodAdmissionEvidence allEvidence(evidence);
    const BObolLodAdmissionPlan plan =
	BObolLodAdmissionPlanner::recordHeadroomRejection(allEvidence, BObolLodAdmissionCursor(),
	    viewEpoch, policyEpoch, activePopulationCost, attemptedBudget);
    evidence = plan.nextEvidence.headroom();
}

static int
test_headroom_policy(void)
{
    BObolLodHeadroomEvidence policy;
    const BObolLodViewEpoch view(4);
    const BObolLodPolicyEpoch lodPolicy(9);
    if (policy.retryPending() ||
	!apply_headroom_retry_plan(policy, view, lodPolicy, 1000, 100) ||
	!policy.retryPending() ||
	!policy.pendingMatches(view, lodPolicy, 1000) ||
	policy.pendingMatches(BObolLodViewEpoch(5), lodPolicy, 1000) ||
	policy.pendingMatches(view, lodPolicy, 1001) ||
	apply_headroom_retry_plan(policy, view, lodPolicy, 1000, 100) ||
	/* The measured allowance may legitimately grow between arming the replay
	 * and consuming its timing sample.  Witness identity is the immutable
	 * population and camera/policy epoch, not that evolving allowance. */
	!apply_headroom_consumption_plan(policy, view, lodPolicy, 101, 1000, 107.0L, 90, 100,
	    true) || policy.retryPending() ||
	apply_headroom_retry_plan(policy, view, lodPolicy, 1000, 101)) {
	std::fprintf(stderr, "FAIL: explicit headroom witness lifecycle\n");
	return 1;
    }
    /* Calibration may discover more capacity after the first witness.  Only
     * material growth above the monotonic high-water may revisit the same
     * population, which admits useful refinement without budget jitter. */
    if (apply_headroom_retry_plan(policy, view, lodPolicy, 1000, 126) ||
	!apply_headroom_retry_plan(policy, view, lodPolicy, 1000, 127) ||
	apply_headroom_retry_plan(policy, view, lodPolicy, 1000, 127) ||
	!apply_headroom_consumption_plan(policy, view, lodPolicy, 127, 1000, 140.0L, 90, 100,
	    true) ||
	apply_headroom_retry_plan(policy, view, lodPolicy, 1000, 158) ||
	!apply_headroom_retry_plan(policy, view, lodPolicy, 1000, 159) ||
	apply_headroom_consumption_plan(policy, view, lodPolicy, 159, 1000, 160.0L, 90, 100,
	    true) ||
	apply_headroom_retry_plan(policy, view, lodPolicy, 1000, 159)) {
	std::fprintf(stderr, "FAIL: monotonic headroom budget high-water\n");
	return 1;
    }
    /* A minimax redistribution changes active render cost without creating a
     * new renderer-capacity epoch.  It must not erase the high-water witness;
     * only material budget growth may re-arm in the same view/policy. */
    if (apply_headroom_retry_plan(policy, view, lodPolicy, 1001, 159) ||
	!apply_headroom_retry_plan(policy, view, lodPolicy, 1001, 200) ||
	!apply_headroom_consumption_plan(policy, view, lodPolicy, 200, 1001,
	    300.0L, 10, 100, true) ||
	!apply_headroom_retry_plan(policy, BObolLodViewEpoch(5), lodPolicy, 2000, 120) ||
	!apply_headroom_consumption_plan(policy, BObolLodViewEpoch(5), lodPolicy, 120, 2000,
	    200.0L, 10, 100, true) ||
	!apply_headroom_retry_plan(policy, BObolLodViewEpoch(5), BObolLodPolicyEpoch(10),
	    2000, 120) ||
	!apply_headroom_consumption_plan(policy, BObolLodViewEpoch(5),
	    BObolLodPolicyEpoch(10), 120, 2000, 200.0L, 10, 100, true)) {
	std::fprintf(stderr, "FAIL: headroom progress witnesses\n");
	return 1;
    }
    if (!apply_headroom_retry_plan(policy, BObolLodViewEpoch(6), lodPolicy, 3000, 200) ||
	apply_headroom_consumption_plan(policy, BObolLodViewEpoch(6), lodPolicy, 200, 3000,
	    210.0L, 10, 100, true) ||
	apply_headroom_retry_plan(policy, BObolLodViewEpoch(6), lodPolicy, 3000, 200)) {
	std::fprintf(stderr, "FAIL: exact threshold headroom admitted\n");
	return 1;
    }
    if (!apply_headroom_retry_plan(policy, BObolLodViewEpoch(7), lodPolicy, 3000, 200) ||
	apply_headroom_consumption_plan(policy, BObolLodViewEpoch(7), lodPolicy, 200, 3000,
	    300.0L, 120, 100, true)) {
	std::fprintf(stderr, "FAIL: slow headroom sample admitted\n");
	return 1;
    }
    if (!apply_headroom_retry_plan(policy, BObolLodViewEpoch(8), lodPolicy, 3000, 200) ||
	apply_headroom_consumption_plan(policy, BObolLodViewEpoch(8), lodPolicy, 200, 3000,
	    300.0L, 10, 100, false)) {
	std::fprintf(stderr, "FAIL: non-reusable headroom sample admitted\n");
	return 1;
    }
    /* Cold timing can leave an old conservative allocation stamped into the
     * occurrences after the scene allowance has already grown.  A reusable
     * frame which proves headroom over the drawn population must trigger one
     * recomputation even when it cannot saturate and thereby prove the full
     * (already larger) allowance. */
    if (!apply_headroom_retry_plan(policy, BObolLodViewEpoch(80), lodPolicy, 100, 500) ||
	!apply_headroom_consumption_plan(policy, BObolLodViewEpoch(80), lodPolicy, 500, 100,
	    250.0L, 50, 100, true) ||
	apply_headroom_retry_plan(policy, BObolLodViewEpoch(80), lodPolicy, 100, 500)) {
	std::fprintf(stderr, "FAIL: stale underallocated population witness\n");
	return 1;
    }
    const BObolLodViewEpoch transientView(81);
    if (!apply_headroom_retry_plan(policy, transientView, lodPolicy, 3001, 200) ||
	!apply_headroom_deferral_plan(policy, transientView, lodPolicy, 3001) ||
	!policy.retryPending() ||
	!apply_headroom_deferral_plan(policy, transientView, lodPolicy, 3001) ||
	apply_headroom_deferral_plan(policy, transientView, lodPolicy, 3001) ||
	!policy.retryPending() ||
	apply_headroom_deferral_plan(policy, BObolLodViewEpoch(82), lodPolicy,
	    3001) ||
	!apply_headroom_consumption_plan(policy, transientView, lodPolicy, 200, 3001,
	    300.0L, 10, 100, true) || policy.retryPending() ||
	apply_headroom_retry_plan(policy, transientView, lodPolicy, 3001, 200)) {
	std::fprintf(stderr, "FAIL: bounded transient headroom replay\n");
	return 1;
    }
    if (!apply_headroom_retry_plan(policy, BObolLodViewEpoch(9), lodPolicy, 3000, 200) ||
	apply_headroom_consumption_plan(policy, BObolLodViewEpoch(10), lodPolicy, 200, 3000,
	    300.0L, 10, 100, true) || policy.retryPending() ||
	!apply_headroom_retry_plan(policy, BObolLodViewEpoch(9), lodPolicy, 3000, 200)) {
	std::fprintf(stderr, "FAIL: stale headroom witness handling\n");
	return 1;
    }
    /* A deadline miss is a consumed negative witness, not a cancellation.
     * A recovery redistribution may not re-arm the failed allowance merely
     * because its active cost differs, but material capacity growth and a new
     * view epoch remain valid progress witnesses. */
    const BObolLodViewEpoch rejectedView(91);
    if (!apply_headroom_retry_plan(policy, rejectedView, lodPolicy, 4000, 300) ||
	!policy.retryPending()) {
	std::fprintf(stderr, "FAIL: rejected headroom witness setup\n");
	return 1;
    }
    apply_headroom_rejection(policy, rejectedView, lodPolicy, 4000, 320);
    if (policy.retryPending() ||
	apply_headroom_retry_plan(policy, rejectedView, lodPolicy, 3900, 400) ||
	!apply_headroom_retry_plan(policy, rejectedView, lodPolicy, 3900, 401)) {
	std::fprintf(stderr, "FAIL: rejected headroom witness reopened\n");
	return 1;
    }
    apply_headroom_rejection(policy, rejectedView, lodPolicy, 4000, 401);
    if (!apply_headroom_retry_plan(policy, BObolLodViewEpoch(92), lodPolicy, 4000, 300)) {
	std::fprintf(stderr, "FAIL: rejected headroom poisoned new view\n");
	return 1;
    }
    policy.cancelRetry();
    if (policy.retryPending() || policy.pendingMatches(
	    BObolLodViewEpoch(9), lodPolicy, 3000)) {
	std::fprintf(stderr, "FAIL: headroom witness cancellation\n");
	return 1;
    }
    return 0;
}

static int
test_coverage_policy(void)
{
    using Policy = BObolLodCoveragePolicy;
    Policy policy;
    if (!policy.active() || policy.sawBoundedSource() ||
	policy.coverageComplete() || policy.visibleCount() != 0 ||
	policy.coveredCount() != 0 || policy.hasCompleteVisibleCount() ||
	policy.required() || policy.effectiveActive() ||
	!policy.effectiveComplete() || !policy.compactionAllowed()) {
	std::fprintf(stderr, "FAIL: initial coverage policy\n");
	return 1;
    }
    policy.setRequired(true);
    if (!policy.required() || !policy.effectiveActive() ||
	policy.effectiveComplete() || policy.compactionAllowed() ||
	!policy.shouldResumeSubmission(false, false, false) ||
	policy.shouldResumeSubmission(true, false, false) ||
	policy.shouldResumeSubmission(false, true, false) ||
	policy.shouldResumeSubmission(false, false, true)) {
	std::fprintf(stderr, "FAIL: required coverage gate\n");
	return 1;
    }

    policy.markBoundedSource();
    policy.observe(4, 3);
    policy.observe(SIZE_MAX, SIZE_MAX);
    if (!policy.sawBoundedSource() || policy.visibleCount() != SIZE_MAX ||
	policy.coveredCount() != SIZE_MAX) {
	std::fprintf(stderr, "FAIL: bounded coverage accumulation\n");
	return 1;
    }
    Policy::Completion completion = policy.completeIfReady(false, false);
    if (completion.completed || !policy.active()) {
	std::fprintf(stderr, "FAIL: incomplete coverage pass retired\n");
	return 1;
    }
    completion = policy.completeIfReady(true, true);
    if (completion.completed || !policy.active()) {
	std::fprintf(stderr, "FAIL: rescan coverage pass retired\n");
	return 1;
    }

    policy.clearPassCounters();
    policy.observe(10, 7);
    policy.markBoundedSource();
    completion = policy.completeIfReady(true, false);
    if (!completion.completed || !completion.bounded ||
	!completion.missing || completion.visibleCount != 10 ||
	completion.coveredCount != 7 || policy.active() ||
	policy.coverageComplete() || !policy.hasCompleteVisibleCount() ||
	policy.completeVisibleCount() != 10 || policy.sawBoundedSource() ||
	policy.visibleCount() != 0 || policy.coveredCount() != 0 ||
	policy.shouldResumeSubmission(false, false, false)) {
	std::fprintf(stderr, "FAIL: missing bounded coverage proof\n");
	return 1;
    }
    if (policy.confirmPresentedCoverage(9) || policy.coverageComplete() ||
	!policy.confirmPresentedCoverage(10) || !policy.coverageComplete()) {
	std::fprintf(stderr, "FAIL: renderer-confirmed coverage proof\n");
	return 1;
    }
    policy.activate(false);
    policy.observe(3, 2);
    if (!policy.confirmPresentedCoverage(10) || policy.active() ||
	policy.visibleCount() != 0 || policy.coveredCount() != 0) {
	std::fprintf(stderr,
	    "FAIL: renderer proof did not supersede retry census\n");
	return 1;
    }

    /* A quality-only pass preserves the latest view denominator and its
	 * minimum-mesh proof while collecting a fresh set of counters. */
    policy.activate(false);
    policy.observe(5, 5);
    policy.noteDemandDeferred();
    completion = policy.completeIfReady(true, false);
    if (!completion.completed || completion.bounded || completion.missing ||
	!completion.demandDeferred || !policy.coverageComplete() ||
	policy.completeVisibleCount() != 5 ||
	!policy.demandCensusRequired()) {
	std::fprintf(stderr, "FAIL: complete unbounded coverage proof\n");
	return 1;
    }

    /* Structural coverage cannot authorize an empty sparse demand frontier.
     * A successor which is itself presentation-deferred retains the
     * obligation; the first complete ordinary demand census retires it. */
    policy.noteDemandDeferred();
    if (policy.completeDemandCensus() ||
	!policy.demandCensusRequired() ||
	!policy.completeDemandCensus() || policy.demandCensusRequired()) {
	std::fprintf(stderr, "FAIL: dense demand successor lifecycle\n");
	return 1;
    }

    policy.activate(true);
    if (!policy.active() || policy.coverageComplete() ||
	!policy.hasCompleteVisibleCount() || policy.demandCensusRequired()) {
	std::fprintf(stderr, "FAIL: view invalidation semantics\n");
	return 1;
    }

    /* A replacement view/inventory epoch must not inherit a partially
     * consumed pass.  The last completed denominator remains available until
     * its owner explicitly retires it, but new observations start at zero. */
    policy.observe(6, 6);
    policy.markBoundedSource();
    policy.activate(true);
    if (!policy.active() || policy.sawBoundedSource() ||
	policy.visibleCount() != 0 || policy.coveredCount() != 0 ||
	!policy.hasCompleteVisibleCount() ||
	policy.completeVisibleCount() != 5) {
	std::fprintf(stderr, "FAIL: invalidated coverage pass retained counters\n");
	return 1;
    }
    policy.observe(6, 6);
    completion = policy.completeIfReady(true, false);
    if (!completion.completed || completion.visibleCount != 6 ||
	completion.coveredCount != 6 || completion.missing ||
	policy.completeVisibleCount() != 6) {
	std::fprintf(stderr, "FAIL: replacement coverage pass denominator\n");
	return 1;
    }

    /* Non-invalidating activation is also used to resume a bounded logical
     * pass, so it must retain observations collected by earlier windows. */
    policy.activate(false);
    policy.observe(2, 2);
    policy.activate(false);
    policy.observe(3, 3);
    completion = policy.completeIfReady(true, false);
    if (!completion.completed || completion.visibleCount != 5 ||
	completion.coveredCount != 5 || completion.missing) {
	std::fprintf(stderr, "FAIL: resumed coverage pass lost counters\n");
	return 1;
    }

    policy.activate(true);
    policy.clearCompleteVisibleCount();
    if (policy.hasCompleteVisibleCount() ||
	policy.completeVisibleCount() != 0) {
	std::fprintf(stderr, "FAIL: visible denominator invalidation\n");
	return 1;
    }
    policy.setCompleteVisibleCount(17);
    if (!policy.hasCompleteVisibleCount() ||
	policy.completeVisibleCount() != 17) {
	std::fprintf(stderr, "FAIL: exact-delta visible denominator update\n");
	return 1;
    }
    if (!Policy::needsDeferredQualitySuccessor(true, false, true, false) ||
	Policy::needsDeferredQualitySuccessor(true, false, true, true) ||
	Policy::needsDeferredQualitySuccessor(true, true, true, false) ||
	Policy::needsDeferredQualitySuccessor(false, false, true, false) ||
	Policy::needsDeferredQualitySuccessor(true, false, false, false)) {
	std::fprintf(stderr,
	    "FAIL: owned quality work re-entered coverage successor\n");
	return 1;
    }
    policy.reset();
    if (policy.active() || policy.coverageComplete() ||
	policy.sawBoundedSource() || policy.hasCompleteVisibleCount() ||
	!policy.required() || policy.effectiveComplete() ||
	policy.compactionAllowed()) {
	std::fprintf(stderr, "FAIL: coverage reset\n");
	return 1;
    }
    policy.setRequired(false);
    if (policy.effectiveActive() || !policy.effectiveComplete() ||
	!policy.compactionAllowed()) {
	std::fprintf(stderr, "FAIL: optional coverage gate\n");
	return 1;
    }
    return 0;
}

static int
test_visibility_census(void)
{
    BObolLodVisibilityCensus census;
    static const BObolLodVisibilityCensus::SourceKey sourceA = 101;
    static const BObolLodVisibilityCensus::SourceKey sourceB = 202;

    census.begin(sourceA, 4);
    census.observe(sourceA, 4, 0, true);
    census.observe(sourceA, 4, 1, true);
    census.observe(sourceA, 4, 2, false);
    census.observe(sourceA, 4, 3, true);
    census.finish(sourceA);
    census.setCount(sourceB, 2);
    if (!census.complete(sourceA) || census.complete(sourceB) ||
	census.sourceCount(sourceA) != 3 ||
	census.sourceCount(sourceB) != 2 || census.total() != 5 ||
	!census.visible(sourceA, 0) || !census.visible(sourceA, 1) ||
	census.visible(sourceA, 2) || !census.visible(sourceA, 3) ||
	census.visible(sourceA, 4) || census.visible(sourceB, 0)) {
	std::fprintf(stderr, "FAIL: completed visibility census\n");
	return 1;
    }

    /* Exact deltas must revise, not erase, the complete all-source count.
	 * Replaying an identical delta is idempotent. */
    census.observe(sourceA, 4, 1, false);
    census.observe(sourceA, 4, 1, false);
    if (census.sourceCount(sourceA) != 2 || census.total() != 4 ||
	census.visible(sourceA, 1)) {
	std::fprintf(stderr, "FAIL: exact visibility removal census\n");
	return 1;
    }
    census.observe(sourceA, 4, 2, true);
    if (census.sourceCount(sourceA) != 3 || census.total() != 5 ||
	!census.visible(sourceA, 2)) {
	std::fprintf(stderr, "FAIL: exact visibility restoration census\n");
	return 1;
    }

    /* Appended dense indices remain bounded by one byte per entry and retain
	 * the already completed source count until a population rescan begins. */
    census.observe(sourceA, 150000, 149999, true);
    if (census.sourceCount(sourceA) != 4 || census.total() != 6) {
	std::fprintf(stderr, "FAIL: extended visibility census\n");
	return 1;
    }

    /* The live source-domain sweep creates placeholders for sources whose
     * bounded coverage has not reached them yet, is idempotent for an
     * unchanged domain, and retires stale source totals atomically. */
    census.beginSourceSetUpdate();
    const bool sourceAExisting = !census.retainSource(sourceA);
    const bool sourceBExisting = !census.retainSource(sourceB);
    const bool unchangedSources = !census.endSourceSetUpdate();
    if (!sourceAExisting || !sourceBExisting || !unchangedSources ||
	census.sourceEntryCount() != 2 || census.total() != 6) {
	std::fprintf(stderr, "FAIL: unchanged visibility source domain\n");
	return 1;
    }
    static const BObolLodVisibilityCensus::SourceKey sourceC = 303;
    census.beginSourceSetUpdate();
    const bool sourceCAdded = census.retainSource(sourceC);
    const bool retiredOldSources = census.endSourceSetUpdate();
    if (!sourceCAdded || !retiredOldSources ||
	census.sourceEntryCount() != 1 || census.total() != 0 ||
	census.sourceCount(sourceA) != 0) {
	std::fprintf(stderr, "FAIL: replaced visibility source domain\n");
	return 1;
    }

    census.begin(sourceA, 150000);
    if (census.complete(sourceA) || census.sourceCount(sourceA) != 0 ||
	census.total() != 0) {
	std::fprintf(stderr, "FAIL: invalidated visibility census\n");
	return 1;
    }
    census.clear();
    if (census.total() != 0) {
	std::fprintf(stderr, "FAIL: cleared visibility census\n");
	return 1;
    }

    /* A client may draw many independent top-level roots rather than one
     * compact hierarchy.  Exercise the fixed-ID index at production scale so
     * source-domain reconciliation cannot quietly return to linear lookup per
     * source. */
    static const size_t sourceScale = 50000;
    census.beginSourceSetUpdate();
    size_t addedSources = 0;
    for (size_t i = 0; i < sourceScale; ++i) {
	if (census.retainSource(static_cast<uint64_t>(i + 1)))
	    addedSources++;
    }
    if (census.endSourceSetUpdate() || addedSources != sourceScale ||
	census.sourceEntryCount() != sourceScale) {
	std::fprintf(stderr, "FAIL: large visibility source-domain creation\n");
	return 1;
    }
    census.beginSourceSetUpdate();
    addedSources = 0;
    for (size_t i = 0; i < sourceScale; i += 2) {
	if (census.retainSource(static_cast<uint64_t>(i + 1)))
	    addedSources++;
    }
    if (!census.endSourceSetUpdate() || addedSources != 0 ||
	census.sourceEntryCount() != sourceScale / 2) {
	std::fprintf(stderr, "FAIL: large visibility source-domain retirement\n");
	return 1;
    }
    return 0;
}

static int
test_source_profile_policy(void)
{
    using Policy = BObolLodSourceProfilePolicy;
    static const uint64_t mebibyte = 1024ULL * 1024ULL;
    if (!Policy::safeMeshFirstPreview(
	    true, 2249, 709, 12 * mebibyte, 2 * mebibyte, 709) ||
	!Policy::safeMeshFirstPreview(
	    true, 2248, 2242, 4 * mebibyte, 156088, 709) ||
	Policy::safeMeshFirstPreview(
	    false, 2249, 709, 12 * mebibyte, 2 * mebibyte, 709) ||
	Policy::safeMeshFirstPreview(
	    true, 50000, 709, 12 * mebibyte, 2 * mebibyte, 709) ||
	Policy::safeMeshFirstPreview(
	    true, 2249, 709, 128 * mebibyte, 2 * mebibyte, 709) ||
	Policy::safeMeshFirstPreview(
	    true, 2249, 709, 12 * mebibyte, 32 * mebibyte, 709) ||
	Policy::safeMeshFirstPreview(
	    true, 2249, 709, 12 * mebibyte, 2 * mebibyte, 3000)) {
	std::fprintf(stderr, "FAIL: source-profile mesh-first admission\n");
	return 1;
    }
    if (!Policy::admitTerminalMesh(
	    true, true, false, true, true, 4000, 0, 4000) ||
	!Policy::admitTerminalMesh(
	    true, true, false, true, true, 4000, 4000, 0) ||
	Policy::admitTerminalMesh(
	    false, true, false, true, true, 4000, 0, 8000) ||
	Policy::admitTerminalMesh(
	    true, false, false, true, true, 4000, 0, 8000) ||
	Policy::admitTerminalMesh(
	    true, true, true, true, true, 4000, 0, 8000) ||
	Policy::admitTerminalMesh(
	    true, true, false, false, true, 4000, 0, 8000) ||
	Policy::admitTerminalMesh(
	    true, true, false, true, false, 4000, 0, 8000) ||
	Policy::admitTerminalMesh(
	    true, true, false, true, true, 8001, 0, 8000) ||
	Policy::admitTerminalMesh(
	    true, true, false, true, true, 5000, 4000, 999)) {
	std::fprintf(stderr, "FAIL: terminal mesh budget admission\n");
	return 1;
    }
    return 0;
}

static int
test_projected_demand_cache(void)
{
    using Cache = BObolLodProjectedDemandCache;
    Cache cache;
    Cache::Source *source = cache.bind(101, 3, 7, 4);
    Cache::Evidence evidence;
    if (!source || Cache::lookup(source, 1, 11, 13, evidence)) {
	std::fprintf(stderr, "FAIL: initial projected-demand cache miss\n");
	return 1;
    }

    Cache::Evidence stored;
    stored.visible = true;
    stored.pixelDiameter = 42.0f;
    stored.pixelArea = 100.0f;
    stored.pixelPerimeter = 50.0f;
    stored.boundsContained = true;
    stored.presentationValid = true;
    stored.presentationVisible = true;
    stored.presentationContained = true;
    stored.presentationPixelWidth = 0.5f;
    stored.presentationPixelHeight = 0.75f;
    Cache::store(source, 1, 11, 13, stored);
    source = cache.bind(101, 3, 7, 8);
    if (!Cache::lookup(source, 1, 11, 13, evidence) ||
	std::fabs(evidence.pixelDiameter - 42.0f) > 1.0e-6f ||
	!evidence.presentationValid ||
	std::fabs(evidence.presentationPixelHeight - 0.75f) > 1.0e-6f) {
	std::fprintf(stderr,
	    "FAIL: same-camera projected-demand cache reuse\n");
	return 1;
    }
    /* Appending dense entries preserves prior evidence.  Exact geometry or
     * placement mutations invalidate only the changed entry. */
    if (Cache::lookup(source, 1, 12, 13, evidence) ||
	Cache::lookup(source, 1, 11, 14, evidence) ||
	Cache::lookup(source, 7, 11, 13, evidence)) {
	std::fprintf(stderr,
	    "FAIL: projected-demand entry revision validation\n");
	return 1;
    }
    /* A camera or compact-population epoch is a different proof even when
     * entry-local revisions happen to repeat. */
    source = cache.bind(101, 3, 8, 8);
    if (Cache::lookup(source, 1, 11, 13, evidence)) {
	std::fprintf(stderr, "FAIL: projected-demand camera invalidation\n");
	return 1;
    }
    Cache::store(source, 1, 11, 13, stored);
    source = cache.bind(101, 4, 8, 8);
    if (Cache::lookup(source, 1, 11, 13, evidence)) {
	std::fprintf(stderr,
	    "FAIL: projected-demand population invalidation\n");
	return 1;
    }
    cache.clear();
    source = cache.bind(101, 4, 8, 1);
    if (!source || Cache::lookup(source, 0, 1, 1, evidence)) {
	std::fprintf(stderr, "FAIL: projected-demand cache clear\n");
	return 1;
    }
    return 0;
}

static int
test_interaction_start_certificate(void)
{
    BObolLodInteractionStartCertificate certificate;
    BObolLodViewScaleSnapshot scale;
    scale.haveCamera = 1;
    scale.size = 8.0;

    if (certificate.captured() || certificate.exactScale() ||
	certificate.startedFromReadyView() || certificate.stableBudget() != 0 ||
	certificate.returnedToScale(scale)) {
	std::fprintf(stderr, "FAIL: empty interaction certificate\n");
	return 1;
    }

    certificate.capture(scale, false, true, 4096);
    if (!certificate.captured() || certificate.exactScale() ||
	!certificate.startedFromReadyView() ||
	certificate.stableBudget() != 4096 ||
	certificate.returnedToScale(scale)) {
	std::fprintf(stderr, "FAIL: interaction certificate without scale proof\n");
	return 1;
    }

    certificate.capture(scale, true, false, 8192);
    BObolLodViewScaleSnapshot changed = scale;
    changed.size = 4.0;
    if (!certificate.exactScale() ||
	certificate.startedFromReadyView() ||
	certificate.stableBudget() != 8192 ||
	!certificate.returnedToScale(scale) ||
	certificate.returnedToScale(changed)) {
	std::fprintf(stderr, "FAIL: exact interaction certificate\n");
	return 1;
    }

    certificate.reset();
    if (certificate.captured() || certificate.exactScale() ||
	certificate.startedFromReadyView() || certificate.stableBudget() != 0 ||
	certificate.returnedToScale(scale)) {
	std::fprintf(stderr, "FAIL: reset interaction certificate\n");
	return 1;
    }
    return 0;
}

static int
test_static_quality_policy(void)
{
    using Policy = BObolLodAdmissionPlanner;
    BObolLodStaticQualityTrial trial;
    if (trial.inProgress() || trial.blocksNewTrial() ||
	trial.capacityRejected() || trial.accepted() ||
	trial.sampledCeiling() != -1) {
	std::fprintf(stderr, "FAIL: static quality initial state\n");
	return 1;
    }
    trial.begin(16.0f);
    trial.noteSampledCeiling(7);

    const BObolLodAdmissionRevisionStamp stamp = admission_revision_stamp();
    if (!trial.probing() || !trial.inProgress() ||
	trial.capacityRejected() ||
	trial.sampledCeiling() != 7 ||
	std::fabs(trial.baselinePointProxyPixelThreshold() - 16.0f) >
	    1.0e-6f) {
	std::fprintf(stderr, "FAIL: static quality active trial\n");
	return 1;
    }
    BObolLodStaticQualityTrial::Constraint rejection;
    rejection.reason = BObolLodStaticQualityTrial::ConstraintReason::
	PRESENTATION_DEADLINE;
    rejection.revisionStamp = stamp;
    rejection.committedCeiling = 7;
    rejection.candidateCeiling = 8;
    rejection.committedCost = 800;
    rejection.candidateCost = 1200;
    rejection.allowedCost = 1000;
    if (!trial.reject(rejection) || !trial.inProgress() ||
	!trial.capacityRejected() || trial.accepted() ||
	trial.sampledCeiling() != -1 ||
	std::fabs(trial.baselinePointProxyPixelThreshold() - 16.0f) >
	    1.0e-6f) {
	std::fprintf(stderr, "FAIL: static quality reconciliation state\n");
	return 1;
    }
    if (!trial.completeReconciliation() ||
	trial.completeReconciliation()) {
	std::fprintf(stderr,
	    "FAIL: static quality reconciliation was not consumed once\n");
	return 1;
    }
    if (trial.inProgress() || !trial.capacityRejected() ||
	!trial.blocksNewTrial() ||
	!trial.usesStaticDeadline() ||
	trial.constraint().reason != rejection.reason ||
	!trial.rejectedFor(stamp)) {
	std::fprintf(stderr, "FAIL: static quality constrained state\n");
	return 1;
    }
    BObolLodAdmissionRevisionStamp changedStamp =
	BObolLodRevisionContract::advance(
	    stamp, BObolLodAdmissionRevisionDomain::VIEW);
    if (trial.rejectedFor(changedStamp)) {
	std::fprintf(stderr,
	    "FAIL: stale static quality constraint remained current\n");
	return 1;
    }
    changedStamp = BObolLodRevisionContract::advance(
	stamp, BObolLodAdmissionRevisionDomain::CAPACITY);
    if (!trial.rejectedFor(changedStamp)) {
	std::fprintf(stderr,
	    "FAIL: completed frame invalidated its capacity constraint\n");
	return 1;
    }
    trial.begin();
    if (trial.inProgress() || !trial.capacityRejected()) {
	std::fprintf(stderr, "FAIL: static quality constraint reopened\n");
	return 1;
    }

    trial.deactivate();
    if (trial.inProgress() || trial.blocksNewTrial() ||
	trial.capacityRejected() || trial.usesStaticDeadline()) {
	std::fprintf(stderr,
	    "FAIL: static quality constraint survived invalidation\n");
	return 1;
    }

    BObolLodStaticQualityTrial::Constraint predicted;
    predicted.reason = BObolLodStaticQualityTrial::ConstraintReason::
	PROTECTED_MINIMUM;
    predicted.revisionStamp = stamp;
    predicted.candidateCost = 200;
    predicted.allowedCost = 100;
    if (!trial.constrain(predicted) || trial.inProgress() ||
	!trial.capacityRejected() || !trial.blocksNewTrial() ||
	!trial.usesStaticDeadline()) {
	std::fprintf(stderr,
	    "FAIL: predicted static constraint was not terminal\n");
	return 1;
    }
    if (trial.constrain(predicted)) {
	std::fprintf(stderr,
	    "FAIL: predicted static constraint was overwritten\n");
	return 1;
    }

    trial.reset();
    trial.begin();
    predicted.committedCeiling = 4;
    predicted.candidateCeiling = 4;
    if (!trial.constrainPresented(predicted) || trial.inProgress() ||
	    !trial.rejectedFor(stamp)) {
	std::fprintf(stderr,
	    "FAIL: completed presentation constraint was not terminal\n");
	return 1;
    }
    if (!trial.retainsRendererCeilingFor(stamp, predicted.committedCeiling) ||
	trial.retainsRendererCeilingFor(stamp,
	    predicted.committedCeiling + 1)) {
	std::fprintf(stderr,
	    "FAIL: completed presentation constraint ceiling scope\n");
	return 1;
    }

    /* A protected-minimum constraint is the sole terminal case which may
     * retain its renderer guard after occurrence-local handoff. */
    BObolLodStaticQualityTrial::HandoffCompletion handoff =
	trial.completeOccurrenceHandoff(stamp, predicted.committedCeiling);
    if (handoff.completedRejectedConstraint ||
	handoff.releaseRendererCeiling ||
	handoff.measureCeilingFreeCandidate) {
	std::fprintf(stderr,
	    "FAIL: protected minimum released its renderer ceiling\n");
	return 1;
    }
    if (trial.constrainPresented(predicted)) {
	std::fprintf(stderr,
	    "FAIL: completed presentation constraint was overwritten\n");
	return 1;
    }

    /* A successful handoff behind a probing ceiling must expose exactly one
     * ceiling-free candidate without retiring the static measurement owner. */
    trial.reset();
    trial.begin();
    handoff = trial.completeOccurrenceHandoff(stamp, 5);
    if (handoff.completedRejectedConstraint ||
	!handoff.releaseRendererCeiling ||
	!handoff.measureCeilingFreeCandidate || !trial.probing()) {
	std::fprintf(stderr,
	    "FAIL: probing handoff did not expose a ceiling-free candidate\n");
	return 1;
    }

    /* A rejected candidate becomes a terminal constraint at this same edge,
     * and its ordinary renderer ceiling is then redundant. */
    trial.reset();
    trial.begin();
    rejection.revisionStamp = stamp;
    if (!trial.reject(rejection)) {
	std::fprintf(stderr,
	    "FAIL: rejected handoff setup\n");
	return 1;
    }
    handoff = trial.completeOccurrenceHandoff(stamp,
	rejection.committedCeiling);
    if (!handoff.completedRejectedConstraint ||
	!handoff.releaseRendererCeiling ||
	handoff.measureCeilingFreeCandidate || trial.inProgress() ||
	!trial.capacityRejected()) {
	std::fprintf(stderr,
	    "FAIL: rejected handoff did not complete reconciliation\n");
	return 1;
    }
    const BObolLodAdmissionRevisionStamp staleConstraintStamp =
	BObolLodRevisionContract::advance(
	    stamp, BObolLodAdmissionRevisionDomain::VIEW);
    if (trial.constrainedPresentationBudgetFor(stamp) !=
	    rejection.allowedCost ||
	trial.constrainedPresentationBudgetFor(staleConstraintStamp) != 0) {
	std::fprintf(stderr,
	    "FAIL: static constraint budget crossed its semantic epoch\n");
	return 1;
    }

    trial.reset();
    BObolLodStaticQualityTrial::Constraint invalidPredicted = predicted;
    invalidPredicted.candidateCost = 0;
    if (trial.constrain(invalidPredicted) || trial.blocksNewTrial()) {
	std::fprintf(stderr,
	    "FAIL: invalid predicted static constraint was accepted\n");
	return 1;
    }

    trial.reset();
    if (std::fabs(trial.baselinePointProxyPixelThreshold() - 1.0f) >
	    1.0e-6f) {
	std::fprintf(stderr, "FAIL: static quality baseline reset\n");
	return 1;
    }
    trial.begin();
    BObolLodStaticQualityTrial::Acceptance acceptance;
    acceptance.revisionStamp = stamp;
    acceptance.ceiling = 7;
    acceptance.nextFraction = 0.625f;
    acceptance.presentedCost = 900;
    acceptance.allowedCost = 1000;
    if (!trial.acceptFractionalCeiling(acceptance) || trial.inProgress() ||
	trial.capacityRejected() || !trial.accepted() ||
	!trial.blocksNewTrial() || !trial.usesStaticDeadline() ||
	trial.acceptance().presentedCost != 900) {
	std::fprintf(stderr, "FAIL: static quality accepted fractional state\n");
	return 1;
    }

    trial.reset();
    trial.restoreRendererCeiling(true);
    if (!trial.probing() || !trial.inProgress() ||
	trial.capacityRejected()) {
	std::fprintf(stderr, "FAIL: static quality renderer restoration\n");
	return 1;
    }
    trial.reset();

    /* A normal quiet recovery miss must not pre-reject the later static
	 * quality phase.  Only a frame attempted by an active static trial owns
	 * that terminal evidence. */
    if (Policy::rejectAfterInterruptedFrame(false, true, false, false) ||
	Policy::rejectAfterInterruptedFrame(true, true, false, true) ||
	Policy::rejectAfterInterruptedFrame(false, false, false, true) ||
	Policy::rejectAfterInterruptedFrame(false, true, true, true) ||
	!Policy::rejectAfterInterruptedFrame(false, true, false, true)) {
	std::fprintf(stderr, "FAIL: static quality rejection scope\n");
	return 1;
    }
    if (Policy::terminalCapacityMiss(true, true, false, true, false) ||
	Policy::terminalCapacityMiss(false, false, false, true, false) ||
	Policy::terminalCapacityMiss(false, true, true, true, false) ||
	Policy::terminalCapacityMiss(false, true, false, false, false) ||
	Policy::terminalCapacityMiss(false, true, false, true, true) ||
	!Policy::terminalCapacityMiss(false, true, false, true, false)) {
	std::fprintf(stderr,
	    "FAIL: changing population poisoned static capacity epoch\n");
	return 1;
    }
    if (Policy::ordinaryHeadroomAllowed(true) ||
	!Policy::ordinaryHeadroomAllowed(false) ||
	!Policy::marginalStaticCapacityAllowed(true, false, false) ||
	!Policy::marginalStaticCapacityAllowed(false, true, false) ||
	!Policy::marginalStaticCapacityAllowed(false, false, true) ||
	Policy::marginalStaticCapacityAllowed(false, false, false)) {
	std::fprintf(stderr,
	    "FAIL: accepted static quality returned to streaming budget\n");
	return 1;
    }
    if (!Policy::rejectedOnePixelTrialReleasesCalibration(
	    true, true, 1.0f) ||
	Policy::rejectedOnePixelTrialReleasesCalibration(
	    false, true, 1.0f) ||
	Policy::rejectedOnePixelTrialReleasesCalibration(
	    true, false, 1.0f) ||
	Policy::rejectedOnePixelTrialReleasesCalibration(
	    true, true, 2.0f) ||
	Policy::rejectedOnePixelTrialReleasesCalibration(
	    true, true, std::numeric_limits<float>::quiet_NaN())) {
	std::fprintf(stderr,
	    "FAIL: rejected static point trial release scope\n");
	return 1;
    }
    if (!Policy::actionableProgressiveQualityDebt(1, 1, 0, 20, 9) ||
	!Policy::actionableProgressiveQualityDebt(4, 3, 0, 9, -1) ||
	Policy::actionableProgressiveQualityDebt(4, 3, 1, 9, -1) ||
	Policy::actionableProgressiveQualityDebt(1, 1, 0, 9, 9) ||
	Policy::actionableProgressiveQualityDebt(1, 1, 0, 9, -1) ||
	!Policy::progressivePresentationQualityDebt(20, 9) ||
	Policy::progressivePresentationQualityDebt(9, 9) ||
	Policy::progressivePresentationQualityDebt(9, -1)) {
	std::fprintf(stderr,
	    "FAIL: static quality did not distinguish resident and presentation debt\n");
	return 1;
    }
    {
	using Reducer = BObolLodTerminalQualityReducer;
	using Owner = Reducer::Owner;
	const std::array<int, 2> ceilings = {{-1, 0}};
	const std::array<size_t, 3> payloadCounts = {{0, 1, 2}};
	for (unsigned int stable = 0; stable != 2; ++stable) {
	    for (unsigned int exact = 0; exact != 2; ++exact) {
		for (const int ceiling : ceilings) {
		    for (const size_t payloadCount : payloadCounts) {
			for (unsigned int combination = 0; combination != 256;
			     ++combination) {
			    const bool staticActive = (combination & 1u) != 0;
			    const bool staticRejected = (combination & 2u) != 0;
			    const bool pointRequired = (combination & 4u) != 0;
			    const bool capacityRequired = (combination & 8u) != 0;
			    const bool ceilingConstrained =
				(combination & 16u) != 0;
			    const bool allocationReplacesCeiling =
				(combination & 32u) != 0;
			    const bool allocationPlanningRequired =
				(combination & 64u) != 0;
			    const bool structuralFrontierPending =
				(combination & 128u) != 0;
			    Reducer::Inputs inputs;
			    inputs.stableTerminalContext = stable != 0;
			    inputs.exactCompletedFrame = exact != 0;
			    inputs.progressiveCeiling = ceiling;
			    inputs.progressivePayloadCount = payloadCount;
			    inputs.staticPhaseActive = staticActive;
			    inputs.staticPhaseRejected = staticRejected;
			    inputs.retainedCeilingConstraint =
				ceilingConstrained;
			    inputs.retainedAllocationReplacesCeiling =
				allocationReplacesCeiling;
			    inputs.structuralFrontierPending =
				structuralFrontierPending;
			    inputs.pointPresentationRequired = pointRequired;
			    inputs.allocationPlanningRequired =
				allocationPlanningRequired;
			    inputs.allocationCapacityRequired = capacityRequired;

			    const bool globalCeiling = ceiling >= 0 &&
				payloadCount != 1;
			    const bool effectiveConstraint = ceilingConstrained &&
				!allocationReplacesCeiling;
			    const bool ceilingDebt = globalCeiling &&
				!effectiveConstraint;
			    Owner expected = Owner::NONE;
			    if (stable && exact && !structuralFrontierPending) {
				if (ceilingDebt)
				    expected = Owner::CEILING_RECONCILIATION;
				else if (globalCeiling && effectiveConstraint)
				    expected = Owner::NONE;
				else if (pointRequired)
				    expected = Owner::POINT_QUALITY;
				else if (allocationPlanningRequired)
				    expected = Owner::ALLOCATION_PLANNING;
				else if (capacityRequired)
				    expected = Owner::ALLOCATION_CAPACITY;
				else if (staticActive && !staticRejected)
				    expected = Owner::STATIC_QUALITY;
			    }

			    const Reducer::Decision decision =
				Reducer::reduce(inputs);
			    if (decision.owner != expected ||
				decision.globalCeilingDebt != ceilingDebt) {
				std::fprintf(stderr,
				    "FAIL: terminal-quality successor ownership\n");
				return 1;
			    }
			}
		    }
		}
	    }
	}
	constexpr size_t SelectedCost = 500;
	constexpr size_t DemandCost = 510;
	constexpr size_t CertifiedBudget = 500;
	if (!Reducer::allocationDemandNeedsCapacity(2, true, true,
		SelectedCost, DemandCost, CertifiedBudget, false) ||
	    Reducer::allocationDemandNeedsCapacity(1, true, true,
		SelectedCost, DemandCost, CertifiedBudget, false) ||
	    Reducer::allocationDemandNeedsCapacity(2, false, true,
		SelectedCost, DemandCost, CertifiedBudget, false) ||
	    Reducer::allocationDemandNeedsCapacity(2, true, false,
		SelectedCost, DemandCost, CertifiedBudget, false) ||
	    Reducer::allocationDemandNeedsCapacity(2, true, true,
		DemandCost, DemandCost, DemandCost, false) ||
	    !Reducer::allocationDemandNeedsCapacity(2, true, true,
		DemandCost, DemandCost + 1, CertifiedBudget, false) ||
	    !Reducer::allocationDemandNeedsCapacity(2, true, true,
		DemandCost, DemandCost, CertifiedBudget, false) ||
	    Reducer::allocationDemandNeedsCapacity(2, true, true,
		SelectedCost, DemandCost, CertifiedBudget, true)) {
	    std::fprintf(stderr,
		"FAIL: terminal allocation demand did not acquire bounded capacity\n");
	    return 1;
	}
	if (!Reducer::allocationDemandNeedsPlanning(
		2, false, false, true, false) ||
	    !Reducer::allocationDemandNeedsPlanning(
		2, true, false, true, false) ||
	    Reducer::allocationDemandNeedsPlanning(
		1, false, false, true, false) ||
	    Reducer::allocationDemandNeedsPlanning(
		2, true, true, true, false) ||
	    Reducer::allocationDemandNeedsPlanning(
		2, false, false, false, false) ||
	    Reducer::allocationDemandNeedsPlanning(
		2, false, false, true, true)) {
	    std::fprintf(stderr,
		"FAIL: stale terminal allocation did not acquire planning\n");
	    return 1;
	}
	if (!Reducer::allocationConstraintPresented(true, true, true) ||
	    Reducer::allocationConstraintPresented(true, false, true) ||
	    Reducer::allocationConstraintPresented(true, true, false) ||
	    Reducer::allocationConstraintPresented(false, true, true)) {
	    std::fprintf(stderr,
		"FAIL: stale allocation retained terminal constraint evidence\n");
	    return 1;
	}
    }
    if (Policy::onePixelTrialRequired(1.0f) ||
	Policy::onePixelTrialRequired(1.01f) ||
	Policy::onePixelTrialRequired(
	    std::numeric_limits<float>::quiet_NaN()) ||
	!Policy::onePixelTrialRequired(1.02f) ||
	!Policy::onePixelTrialRequired(64.0f)) {
	std::fprintf(stderr, "FAIL: static one-pixel trial predicate\n");
	return 1;
    }
    if (std::fabs(Policy::maximumPointProxyPixelThreshold() - 64.0f) >
	    1.0e-6f) {
	std::fprintf(stderr, "FAIL: maximum point-proxy threshold contract\n");
	return 1;
    }
    if (!Policy::acceptSettledOnePixelFrame(false, true, false, false,
	    false, false, false, -1, 1.0f, true,
	    60000000ULL, 50000000ULL, 100000000ULL) ||
	Policy::acceptSettledOnePixelFrame(true, true, false, false,
	    false, false, false, -1, 1.0f, true,
	    60000000ULL, 50000000ULL, 100000000ULL) ||
	Policy::acceptSettledOnePixelFrame(false, true, false, false,
	    false, true, false, -1, 1.0f, true,
	    60000000ULL, 50000000ULL, 100000000ULL) ||
	Policy::acceptSettledOnePixelFrame(false, true, false, false,
	    false, false, true, -1, 1.0f, true,
	    60000000ULL, 50000000ULL, 100000000ULL) ||
	Policy::acceptSettledOnePixelFrame(false, true, false, false,
	    false, false, false, 0, 1.0f, true,
	    60000000ULL, 50000000ULL, 100000000ULL) ||
	Policy::acceptSettledOnePixelFrame(false, true, false, false,
	    false, false, false, -1, 2.0f, true,
	    60000000ULL, 50000000ULL, 100000000ULL) ||
	Policy::acceptSettledOnePixelFrame(false, true, false, false,
	    false, false, false, -1, 1.0f, false,
	    60000000ULL, 50000000ULL, 100000000ULL) ||
	Policy::acceptSettledOnePixelFrame(false, true, false, false,
	    false, false, false, -1, 1.0f, true,
	    101000000ULL, 50000000ULL, 100000000ULL)) {
	std::fprintf(stderr,
	    "FAIL: settled one-pixel frame did not enter static phase directly\n");
	return 1;
    }
    if (Policy::stagedProgressiveCeiling(-1, 10, 8.0f, 100) != -1 ||
	Policy::stagedProgressiveCeiling(3, 3, 8.0f, 100) != -1 ||
	Policy::stagedProgressiveCeiling(3, 10, 8.0f, 100) != -1 ||
	Policy::stagedProgressiveCeiling(3, 10, 1.0f, 1) != 4 ||
	Policy::stagedProgressiveCeiling(3, 10, 8.0f, 1) != 4 ||
	Policy::stagedProgressiveCeiling(9, 10, 8.0f, 1) != 10) {
	std::fprintf(stderr, "FAIL: static staged progressive ceiling\n");
	return 1;
    }
    if (Policy::staticProgressiveBaselineCeiling(-1, 10, 1) != 10 ||
	Policy::staticProgressiveBaselineCeiling(3, 10, 1) != -1 ||
	Policy::staticProgressiveBaselineCeiling(-1, -1, 1) != -1 ||
	Policy::staticProgressiveBaselineCeiling(-1, 10, 2) != -1) {
	std::fprintf(stderr,
	    "FAIL: static progressive baseline ceiling contract\n");
	return 1;
    }
    if (!Policy::staticOrdinalOverscanApplicable(1) ||
	Policy::staticOrdinalOverscanApplicable(0) ||
	Policy::staticOrdinalOverscanApplicable(2)) {
	std::fprintf(stderr,
	    "FAIL: static ordinal overscan escaped the single-payload domain\n");
	return 1;
    }
    if (!Policy::retainedPointAggregationApplicable(true, 1) ||
	Policy::retainedPointAggregationApplicable(false, 1) ||
	Policy::retainedPointAggregationApplicable(true, 0)) {
	std::fprintf(stderr,
	    "FAIL: inert retained point population admitted aggregation\n");
	return 1;
    }
    if (Policy::deadlinePointAggregationApplicable(false, false, 150000) ||
	Policy::deadlinePointAggregationApplicable(false, true, 1) ||
	!Policy::deadlinePointAggregationApplicable(false, true, 150000) ||
	!Policy::deadlinePointAggregationApplicable(true, false, 1)) {
	std::fprintf(stderr,
	    "FAIL: point deadline lost its structural-admission certificate\n");
	return 1;
    }
    if (!Policy::pointProxyThresholdInert(0, true, 0, 0) ||
	!Policy::pointProxyThresholdInert(0, true, 5, 5) ||
	Policy::pointProxyThresholdInert(1, true, 0, 0) ||
	Policy::pointProxyThresholdInert(0, false, 0, 0) ||
	Policy::pointProxyThresholdInert(0, true, 1, 0)) {
	std::fprintf(stderr,
	    "FAIL: structural population made the point threshold inert\n");
	return 1;
    }
    if (!Policy::structuralPointAggregationRequired(50000, 8192) ||
	Policy::structuralPointAggregationRequired(8, 512) ||
	Policy::structuralPointAggregationRequired(50000, SIZE_MAX)) {
	std::fprintf(stderr,
	    "FAIL: internal mesh chunks changed structural aggregation regime\n");
	return 1;
    }
    if (Policy::pointPolicyExhaustedForStructuralFrontier(
	    false, true, true, 0, 0) ||
	!Policy::pointPolicyExhaustedForStructuralFrontier(
	    true, true, false, 20, 0) ||
	!Policy::pointPolicyExhaustedForStructuralFrontier(
	    false, false, true, 20, 0) ||
	!Policy::pointPolicyExhaustedForStructuralFrontier(
	    false, true, true, 20, 20) ||
	Policy::pointPolicyExhaustedForStructuralFrontier(
	    false, true, true, 20, 19) ||
	Policy::pointPolicyExhaustedForStructuralFrontier(
	    false, true, false, 20, 20)) {
	std::fprintf(stderr,
	    "FAIL: structural repair did not consume the point successor exactly\n");
	return 1;
    }
    if (!Policy::pointSuccessorRequiredForStructuralFrontier(
	    true, 20, false, false) ||
	Policy::pointSuccessorRequiredForStructuralFrontier(
	    false, 20, false, false) ||
	Policy::pointSuccessorRequiredForStructuralFrontier(
	    true, 0, false, false) ||
	Policy::pointSuccessorRequiredForStructuralFrontier(
	    true, 20, true, false) ||
	Policy::pointSuccessorRequiredForStructuralFrontier(
	    true, 20, false, true)) {
	std::fprintf(stderr,
	    "FAIL: structural point successor ownership\n");
	return 1;
    }
    if (!Policy::structuralCapacityFrameApplicable(true, true, 1, 1) ||
	!Policy::structuralCapacityFrameApplicable(true, true, 10, 20) ||
	Policy::structuralCapacityFrameApplicable(false, true, 10, 20) ||
	Policy::structuralCapacityFrameApplicable(true, false, 10, 20) ||
	Policy::structuralCapacityFrameApplicable(true, true, 0, 20) ||
	Policy::structuralCapacityFrameApplicable(true, true, 10, 0)) {
	std::fprintf(stderr,
	    "FAIL: zero-cost structural capacity frame ownership\n");
	return 1;
    }
    if (Policy::structuralRepairPresentationDeadline(
	    false, 100, 400, 2) != 200 ||
	Policy::structuralRepairPresentationDeadline(
	    true, 100, 400, 2) != 400 ||
	Policy::structuralRepairPresentationDeadline(
	    true, 250, 400, 2) != 500 ||
	Policy::structuralRepairPresentationDeadline(
	    false, UINT64_MAX, 0, 2) != UINT64_MAX) {
	std::fprintf(stderr,
	    "FAIL: structural repair deadline profile\n");
	return 1;
    }
    if (Policy::structuralFirstWaveOccurrenceLimit(SIZE_MAX) != SIZE_MAX ||
	Policy::structuralFirstWaveOccurrenceLimit(50000) != 781 ||
	Policy::structuralFirstWaveOccurrenceLimit(640000) != 8192) {
	std::fprintf(stderr,
	    "FAIL: structural first-wave occurrence limit\n");
	return 1;
    }
    return 0;
}

static int
test_terminal_reconciliation_composition(void)
{
    using Presentation = BObolLodPresentationPolicy;
    using Terminal = BObolLodTerminalQualityReducer;

    constexpr size_t quietBudget = 100;
    constexpr size_t staticBudget = 200;
    constexpr size_t selectedCost = 150;
    constexpr int rendererCeiling = 5;

    Presentation presentation;
    presentation.armHandoff(false, quietBudget, quietBudget);

    Presentation::CompletedPassInputs pass;
    pass.completed = true;
    pass.retainedAllocationCompleted = true;
    pass.retainedAllocationCertified = true;
    pass.populationQuiescent = true;
    pass.presentationLimitsReconciled =
	Presentation::presentationLimitsReconciled(true, true,
	    selectedCost, quietBudget, quietBudget);
    Presentation::CompletedPassDecision decision =
	presentation.completePass(pass);
    if (!decision.requestLocalPresentationReduction ||
	decision.finishHandoff) {
	std::fprintf(stderr,
	    "FAIL: over-budget terminal handoff skipped local reduction\n");
	return 1;
    }

    /* The local representation cannot be reduced further, so production
     * spends its one distinct static population.  Its larger hard-deadline
     * allocation fits and must commit directly to a ceiling-free frame. */
    BObolLodStaticQualityTrial trial;
    trial.begin();
    presentation.armHandoff(false, quietBudget, staticBudget);
    pass.presentationLimitsReconciled =
	Presentation::presentationLimitsReconciled(true, true,
	    selectedCost, staticBudget, staticBudget);
    decision = presentation.completePass(pass);
    if (!decision.finishHandoff ||
	decision.requestLocalPresentationReduction ||
	presentation.handoffPending()) {
	std::fprintf(stderr,
	    "FAIL: fitting static handoff did not commit once\n");
	return 1;
    }

    const BObolLodAdmissionRevisionStamp stamp = admission_revision_stamp();
    const BObolLodStaticQualityTrial::HandoffCompletion completion =
	trial.completeOccurrenceHandoff(stamp, rendererCeiling);
    if (!completion.releaseRendererCeiling ||
	!completion.measureCeilingFreeCandidate || !trial.probing()) {
	std::fprintf(stderr,
	    "FAIL: fitting static allocation remained behind its old ceiling\n");
	return 1;
    }

    Terminal::Inputs terminalInputs;
    terminalInputs.stableTerminalContext = true;
    terminalInputs.exactCompletedFrame = true;
    terminalInputs.progressiveCeiling = -1;
    terminalInputs.progressivePayloadCount = 2;
    terminalInputs.staticPhaseActive = trial.inProgress();
    if (Terminal::reduce(terminalInputs).owner !=
	    Terminal::Owner::STATIC_QUALITY) {
	std::fprintf(stderr,
	    "FAIL: ceiling-free static candidate lost its terminal owner\n");
	return 1;
    }

    /* A satisfied candidate retires the static owner.  With the semantic
     * inputs unchanged, the reducer must then be a fixed point rather than
     * recreating allocation or reconciliation work. */
    trial.deactivate();
    terminalInputs.staticPhaseActive = trial.inProgress();
    if (Terminal::reduce(terminalInputs).owner != Terminal::Owner::NONE) {
	std::fprintf(stderr,
	    "FAIL: settled static candidate reopened terminal work\n");
	return 1;
    }
    return 0;
}

static int
test_delivery_policy(void)
{
    using Policy = BObolLodDeliveryPolicy;

    /* Speculative residency must not block an already-admitted visible
     * improvement. */
    if (Policy::visibleFirstCut(26, 27, 29, true) != 27) {
	std::fprintf(stderr,
	    "FAIL: delivery policy put speculative residency before presentation\n");
	return 1;
    }
    /* Small independent leaves are cheaper to publish in one batched wave;
     * splitting each one would multiply result/application traffic at 50k. */
    if (Policy::visibleFirstCut(26, 27, 29, false) != 29) {
	std::fprintf(stderr,
	    "FAIL: delivery policy split an inexpensive leaf update\n");
	return 1;
    }
    /* If presentation is budget-limited at the current cut, independent
     * resident prefetch is still required for future frames. */
    if (Policy::visibleFirstCut(26, 26, 29, true) != 29) {
	std::fprintf(stderr,
	    "FAIL: delivery policy suppressed resident-only prefetch\n");
	return 1;
    }
    if (Policy::visibleFirstCut(26, 29, 29, true) != 29 ||
	Policy::visibleFirstCut(26, 29, 27, true) != 27 ||
	Policy::visibleFirstCut(29, 26, 26, true) != 26) {
	std::fprintf(stderr,
	    "FAIL: delivery policy changed a non-speculative delivery\n");
	return 1;
    }

    /* The live allocation remains authoritative while its census exists. */
    if (Policy::authorizedPresentationCut(17, false, 0, 0, -1,
	    4, 7, 19) != 17) {
	std::fprintf(stderr,
	    "FAIL: delivery policy discarded the current allocation\n");
	return 1;
    }
    /* Publishing a new immutable generation may invalidate that census.  A
     * matching provider certificate still owns the already-admitted result. */
    if (Policy::authorizedPresentationCut(-1, true, 4, 7, 21,
	    4, 7, 19) != 21) {
	std::fprintf(stderr,
	    "FAIL: delivery policy discarded a current provider certificate\n");
	return 1;
    }
    if (Policy::authorizedPresentationCut(-1, true, 3, 7, 21,
	    4, 7, 19) >= 0 ||
	Policy::authorizedPresentationCut(-1, true, 4, 6, 21,
	    4, 7, 19) >= 0 ||
	Policy::authorizedPresentationCut(-1, true, 4, 7, 18,
	    4, 7, 19) >= 0 ||
	Policy::authorizedPresentationCut(-1, false, 4, 7, 21,
	    4, 7, 19) >= 0) {
	std::fprintf(stderr,
	    "FAIL: delivery policy admitted an uncertified provider advance\n");
	return 1;
    }
    return 0;
}

static int
test_convergence_policy(void)
{
    using Policy = BObolLodConvergencePolicy;
    Policy policy;
    const auto baseInput = []() {
	Policy::Inputs input;
	input.viewEpoch.set(1);
	input.policyEpoch.set(1);
	return input;
    };
    Policy::Inputs input = baseInput();
    Policy::Decision decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::IDLE ||
	decision.fraction < 0.999f || decision.fraction > 1.001f ||
	!decision.terminal || decision.terminalError ||
	decision.outcome != Policy::Outcome::READY || !decision.viewReady ||
	decision.backgroundPending || decision.visualPending ||
	decision.hasLodState) {
	std::fprintf(stderr, "FAIL: empty convergence state\n");
	return 1;
    }

    input = baseInput();
    input.enabled = false;
    input.expectedLeafCount = 10;
    input.availableLeafCount = 10;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::IDLE || !decision.terminal ||
	!decision.viewReady || decision.hasLodState || decision.visualPending) {
	std::fprintf(stderr, "FAIL: disabled LoD convergence state\n");
	return 1;
    }

    input.publicationPending = true;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::REFINING || decision.terminal ||
	decision.outcome != Policy::Outcome::ACTIVE || decision.viewReady ||
	!decision.visualPending || !decision.hasLodState) {
	std::fprintf(stderr,
	    "FAIL: explicit publication hidden by disabled automatic control\n");
	return 1;
    }

    input.publicationPending = false;
    input.controlPending = true;
    decision = policy.evaluate(input);
    if (decision.terminal || decision.viewReady || !decision.visualPending) {
	std::fprintf(stderr,
	    "FAIL: manual control debt hidden by disabled automatic control\n");
	return 1;
    }

    input = baseInput();
    input.expectedLeafCount = 10;
    input.availableLeafCount = 4;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::DISCOVERING ||
	decision.terminal || decision.outcome != Policy::Outcome::ACTIVE ||
	decision.viewReady || !decision.visualPending ||
	decision.fraction < 0.159f || decision.fraction > 0.161f) {
	std::fprintf(stderr, "FAIL: structural convergence phase\n");
	return 1;
    }

    input = baseInput();
    input.expectedLeafCount = 2;
    input.availableLeafCount = 2;
    input.visibleTargetCount = 2;
    input.activePayloadCount = 2;
    input.visibilityCensusComplete = true;
    decision = policy.evaluate(input);
    if (decision.terminal || decision.viewReady ||
	decision.outcome != Policy::Outcome::ACTIVE) {
	std::fprintf(stderr,
	    "FAIL: unresolved fractional allocation became terminal\n");
	return 1;
    }
    input.pixelDemandPresentationProven = true;
    decision = policy.evaluate(input);
    if (!decision.terminal || !decision.viewReady ||
	decision.outcome != Policy::Outcome::READY ||
	decision.visualPending) {
	std::fprintf(stderr,
	    "FAIL: exact pixel-demand allocation did not become terminal\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(2);
    input.structuralDiscovery = true;
    input.progressiveWorkPending = true;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::DISCOVERING ||
	decision.viewReady || !decision.visualPending ||
	!decision.backgroundPending || decision.fraction > 0.001f) {
	std::fprintf(stderr, "FAIL: provider discovery phase\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(3);
    input.sourcePreparationPending = true;
    input.visibleTargetCount = 10;
    input.activePayloadCount = 4;
    input.calibrationPending = true;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::PREPARING ||
	decision.viewReady || !decision.visualPending ||
	decision.fraction < 0.539f || decision.fraction > 0.541f) {
	std::fprintf(stderr, "FAIL: source preparation precedence\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(32);
    input.sourcePreparationPending = true;
    input.sourcePreparationCompletedUnits = 2;
    input.sourcePreparationTotalUnits = 8;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::PREPARING ||
	decision.viewReady || !decision.visualPending ||
	decision.fraction < 0.487f || decision.fraction > 0.488f) {
	std::fprintf(stderr, "FAIL: exact source preparation rank\n");
	return 1;
    }

    input.structuralDiscovery = true;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::PREPARING ||
	decision.fraction < 0.487f || decision.fraction > 0.488f) {
	std::fprintf(stderr,
	    "FAIL: derived discovery obscured source preparation rank\n");
	return 1;
    }
    input.structuralDiscovery = false;
    input.sourcePreparationCompletedUnits = 8;
    decision = policy.evaluate(input);
    if (decision.fraction < 0.749f || decision.fraction > 0.751f) {
	std::fprintf(stderr, "FAIL: completed source preparation rank\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(30);
    input.sourcePreparationPending = true;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::PREPARING ||
	decision.viewReady || decision.fraction < 0.399f ||
	decision.fraction > 0.401f) {
	std::fprintf(stderr, "FAIL: indeterminate source preparation\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(31);
    input.interactive = true;
    input.visibleTargetCount = 4;
    input.activePayloadCount = 4;
    input.satisfiedPayloadCount = 2;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::INTERACTIVE ||
	decision.viewReady || decision.fraction < 0.624f ||
	decision.fraction > 0.626f) {
	std::fprintf(stderr, "FAIL: interactive convergence phase\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(4);
    input.calibrationPending = true;
    input.visibleTargetCount = 10;
    input.activePayloadCount = 10;
    input.satisfiedPayloadCount = 5;
    input.visibilityCensusComplete = true;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::CALIBRATING ||
	decision.viewReady || decision.fraction < 0.844f ||
	decision.fraction > 0.846f) {
	std::fprintf(stderr, "FAIL: calibration convergence phase\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(42);
    input.calibrationPending = true;
    input.visibleTargetCount = 10;
    input.activePayloadCount = 10;
    input.satisfiedPayloadCount = 10;
    input.visibilityCensusComplete = true;
    input.capacitySearchTotalUnits = 20;
    decision = policy.evaluate(input);
    if (decision.fraction < 0.939f || decision.fraction > 0.941f) {
	std::fprintf(stderr, "FAIL: capacity-search progress origin\n");
	return 1;
    }
    input.capacitySearchCompletedUnits = 5;
    decision = policy.evaluate(input);
    if (decision.fraction < 0.951f || decision.fraction > 0.953f) {
	std::fprintf(stderr, "FAIL: capacity-search partial progress\n");
	return 1;
    }
    input.capacitySearchCompletedUnits = 20;
    decision = policy.evaluate(input);
    if (decision.fraction < 0.989f || decision.fraction > 0.991f) {
	std::fprintf(stderr, "FAIL: capacity-search completed progress\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(41);
    input.calibrationPending = true;
    input.availableLeafCount = 1000;
    input.visibleTargetCount = 10;
    input.activePayloadCount = 10;
    input.satisfiedPayloadCount = 10;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::CALIBRATING ||
	decision.fraction < 0.405f || decision.fraction > 0.407f) {
	std::fprintf(stderr, "FAIL: partial visibility progress denominator\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(5);
    input.resultPending = true;
    input.visibleTargetCount = 4;
    input.activePayloadCount = 4;
    input.satisfiedPayloadCount = 2;
    input.visibilityCensusComplete = true;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::REFINING ||
	decision.viewReady || !decision.visualPending ||
	decision.fraction < 0.844f || decision.fraction > 0.846f) {
	std::fprintf(stderr, "FAIL: refining convergence phase\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(53);
    input.controlPending = true;
    input.visibleTargetCount = 1;
    input.activePayloadCount = 1;
    input.satisfiedPayloadCount = 1;
    input.visibilityCensusComplete = true;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::REFINING || decision.terminal ||
	decision.viewReady || !decision.visualPending) {
	std::fprintf(stderr,
	    "FAIL: named control obligation was misreported as terminal\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(51);
    input.resultPending = true;
    input.visibleTargetCount = 1000;
    input.activePayloadCount = 980;
    input.satisfiedPayloadCount = 970;
    input.presentedSubpixelOccurrenceCount = 20;
    input.visibilityCensusComplete = true;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::REFINING ||
	decision.viewReady || decision.fraction < 0.937f ||
	decision.fraction > 0.939f) {
	std::fprintf(stderr, "FAIL: late refining work estimate\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(52);
    input.resultPending = true;
    input.visibleTargetCount = 150000;
    input.activePayloadCount = 44;
    input.satisfiedPayloadCount = 44;
    input.presentedSubpixelOccurrenceCount = 149479;
    input.presentedStructuralBoxCount = 511;
    input.visibilityCensusComplete = true;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::REFINING ||
	decision.viewReady || decision.fraction < 0.763f ||
	decision.fraction > 0.765f) {
	std::fprintf(stderr, "FAIL: expensive mesh tail work estimate\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(6);
    input.activePayloadCount = 1;
    input.queuedCacheWrites = 1;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::BACKGROUND ||
	!decision.terminal || decision.outcome != Policy::Outcome::READY ||
	!decision.viewReady || !decision.backgroundPending ||
	decision.fraction < 0.999f || decision.fraction > 1.001f) {
	std::fprintf(stderr, "FAIL: background convergence phase\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(61);
    input.visibleTargetCount = 1;
    input.activePayloadCount = 1;
    input.satisfiedPayloadCount = 1;
    input.inFlight = 1;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::BACKGROUND ||
	!decision.viewReady || !decision.backgroundPending ||
	decision.visualPending || decision.fraction < 0.999f ||
	decision.fraction > 1.001f) {
	std::fprintf(stderr, "FAIL: reusable cache background task phase\n");
	return 1;
    }

    input.satisfiedPayloadCount = 0;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::REFINING || decision.viewReady ||
	!decision.visualPending) {
	std::fprintf(stderr, "FAIL: unresolved task remains foreground\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(63);
    input.visibleTargetCount = 1;
    input.activePayloadCount = 1;
    input.satisfiedPayloadCount = 1;
    input.visibilityCensusComplete = true;
    input.resultPending = true;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::BACKGROUND ||
	!decision.terminal || !decision.viewReady ||
	!decision.backgroundPending || decision.visualPending) {
	std::fprintf(stderr,
	    "FAIL: reusable cache result reopened a resolved view\n");
	return 1;
    }

    input.satisfiedPayloadCount = 0;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::REFINING || decision.terminal ||
	decision.viewReady || !decision.visualPending) {
	std::fprintf(stderr,
	    "FAIL: result required by an unresolved view became background\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(62);
    input.visibilityCensusComplete = true;
    input.inFlight = 1;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::BACKGROUND ||
	!decision.viewReady || !decision.backgroundPending ||
	decision.visualPending || decision.fraction < 0.999f ||
	decision.fraction > 1.001f) {
	std::fprintf(stderr,
	    "FAIL: empty-view reusable task remained foreground\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(7);
    input.visibleTargetCount = 2;
    input.activePayloadCount = 2;
    input.satisfiedPayloadCount = 1;
    input.memoryLimitedPayloadCount = 1;
    decision = policy.evaluate(input);
    if (!decision.viewReady || !decision.memoryLimited ||
	!decision.performanceLimited ||
	!decision.terminal ||
	decision.outcome != Policy::Outcome::CONSTRAINED ||
	decision.phase != Policy::Phase::IDLE) {
	std::fprintf(stderr, "FAIL: limited terminal convergence\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(74);
    input.expectedLeafCount = 10;
    input.availableLeafCount = 10;
    input.visibleTargetCount = 10;
    input.visibilityCensusComplete = true;
    decision = policy.evaluate(input);
    if (decision.viewReady || decision.terminal ||
	!decision.visualPending || decision.phase != Policy::Phase::REFINING) {
	std::fprintf(stderr,
	    "FAIL: unrepresented nonempty view reported terminal readiness\n");
	return 1;
    }

    /* Pixel demand may exceed the measured renderer allowance.  A current,
     * fully presented allocation is then the terminal visual contract, but it
     * cannot hide a missing occurrence merely by certifying the represented
     * subset. */
    input = baseInput();
    input.viewEpoch.set(75);
    input.visibleTargetCount = 10;
    input.activePayloadCount = 8;
    input.satisfiedPayloadCount = 4;
    input.presentedSubpixelOccurrenceCount = 2;
    input.visibilityCensusComplete = true;
    input.stableBudgetLimited = true;
    input.constrainedPresentationProven = true;
    decision = policy.evaluate(input);
    if (!decision.viewReady || !decision.terminal ||
	!decision.performanceLimited ||
	decision.outcome != Policy::Outcome::CONSTRAINED ||
	std::fabs(decision.fraction - 1.0f) > 1.0e-6f) {
	std::fprintf(stderr,
	    "FAIL: presented capacity allocation did not terminate quality debt\n");
	return 1;
    }
    input.activePayloadCount = 7;
    decision = policy.evaluate(input);
    if (decision.viewReady || decision.terminal || !decision.visualPending ||
	decision.phase != Policy::Phase::REFINING) {
	std::fprintf(stderr,
	    "FAIL: capacity allocation concealed a missing occurrence\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(71);
    input.visibleTargetCount = 10;
    input.activePayloadCount = 6;
    input.satisfiedPayloadCount = 6;
    input.presentedSubpixelOccurrenceCount = 4;
    decision = policy.evaluate(input);
    if (!decision.viewReady || decision.visualPending ||
	decision.performanceLimited || decision.phase != Policy::Phase::IDLE) {
	std::fprintf(stderr,
	    "FAIL: subpixel occurrence terminal convergence\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(72);
    input.visibleTargetCount = 10;
    input.activePayloadCount = 6;
    input.satisfiedPayloadCount = 6;
    input.presentedSubpixelOccurrenceCount = 3;
    input.presentedStructuralBoxCount = 1;
    decision = policy.evaluate(input);
    if (decision.viewReady || !decision.visualPending ||
	decision.phase != Policy::Phase::REFINING) {
	std::fprintf(stderr,
	    "FAIL: unresolved structural presentation convergence\n");
	return 1;
    }

    input.terminalOccurrenceFailureCount = 1;
    decision = policy.evaluate(input);
    if (decision.viewReady || !decision.terminal ||
	!decision.terminalError || decision.visualPending ||
	decision.outcome != Policy::Outcome::ERROR ||
	decision.phase != Policy::Phase::CONVERROR) {
	std::fprintf(stderr,
	    "FAIL: terminal structural presentation error convergence\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(70);
    input.visibleTargetCount = 1;
    input.activePayloadCount = 1;
    input.satisfiedPayloadCount = 1;
    input.presentationLimited = true;
    decision = policy.evaluate(input);
    if (!decision.viewReady || decision.memoryLimited ||
	!decision.performanceLimited ||
	decision.phase != Policy::Phase::IDLE) {
	std::fprintf(stderr,
	    "FAIL: presentation-limited terminal convergence\n");
	return 1;
    }

    input.failedSourceCount = 1;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::CONVERROR ||
	decision.outcome != Policy::Outcome::ERROR || decision.viewReady ||
	!decision.terminal || !decision.terminalError) {
	std::fprintf(stderr, "FAIL: convergence error phase\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(73);
    input.expectedLeafCount = 10;
    input.availableLeafCount = 4;
    input.failedSourceCount = 1;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::CONVERROR ||
	decision.outcome != Policy::Outcome::ERROR || decision.viewReady ||
	!decision.terminal || !decision.terminalError ||
	decision.visualPending) {
	std::fprintf(stderr,
	    "FAIL: failed source left structural discovery live\n");
	return 1;
    }

    input.sourcePreparationPending = true;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::CONVERROR ||
	decision.outcome != Policy::Outcome::ERROR || decision.viewReady ||
	decision.terminal || decision.terminalError ||
	!decision.visualPending) {
	std::fprintf(stderr,
	    "FAIL: failed source hid remaining preparation work\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(8);
    input.calibrationPending = true;
    decision = policy.evaluate(input);
	if (decision.fraction < 0.399f || decision.fraction > 0.401f) {
	std::fprintf(stderr, "FAIL: convergence fraction seed\n");
	return 1;
    }
    input.calibrationPending = false;
    input.expectedLeafCount = 10;
    input.availableLeafCount = 1;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::DISCOVERING ||
	decision.fraction < 0.399f || decision.fraction > 0.401f ||
	policy.fractionFloor() < 0.399f ||
	policy.fractionFloor() > 0.401f) {
	std::fprintf(stderr, "FAIL: monotonic convergence fraction\n");
	return 1;
    }
    input.viewEpoch.set(9);
    decision = policy.evaluate(input);
    if (decision.fraction < 0.039f || decision.fraction > 0.041f) {
	std::fprintf(stderr, "FAIL: convergence epoch reset\n");
	return 1;
    }
    policy.resetFraction();
    if (policy.fractionFloor() < -0.001f ||
	policy.fractionFloor() > 0.001f) {
	std::fprintf(stderr, "FAIL: explicit convergence reset\n");
	return 1;
    }
    return 0;
}

static int
test_availability_scheduler(void)
{
    using Policy = BObolLodAvailabilityScheduler;
    using Successor = Policy::CompletedPassSuccessor;
    if (Policy::effectiveMicroseconds(8000, true, false) != 8000 ||
	Policy::effectiveMicroseconds(8000, true, true) != 4000 ||
	Policy::effectiveMicroseconds(8000, false, true) != 8000 ||
	Policy::effectiveMicroseconds(2500, true, false) != 2500) {
	std::fprintf(stderr, "FAIL: progressive availability scheduling\n");
	return 1;
    }
    enum SchedulerInput : unsigned int {
	PASS_COMPLETED = 0,
	RESIDENCY_DRAIN,
	RESIDENT_GROWTH,
	RESIDENT_REFINEMENT,
	RESIDENT_WORK,
	POINT_TRIANGLE_RECOVERY,
	POINT_CALIBRATION,
	CAPACITY_PLANNING,
	SCHEDULER_INPUT_COUNT
    };
    const unsigned int schedulerCombinationCount =
	1u << SCHEDULER_INPUT_COUNT;
    for (unsigned int combination = 0;
	    combination < schedulerCombinationCount; ++combination) {
	auto enabled = [combination](SchedulerInput input) {
	    return (combination & (1u << input)) != 0;
	};
	Successor expected = Successor::NONE;
	if (enabled(PASS_COMPLETED)) {
	    if (enabled(RESIDENCY_DRAIN))
		expected = Successor::COMPLETE_RESIDENCY_DRAIN;
	    else if (enabled(RESIDENT_GROWTH))
		expected = Successor::YIELD_TO_RESIDENT_GROWTH;
	    else if (enabled(RESIDENT_REFINEMENT))
		expected = enabled(RESIDENT_WORK) ?
		    Successor::AWAIT_RESIDENT_RESULT :
		    Successor::SUBMIT_RESIDENT_REQUESTS;
	    /* Capacity planning is the allocation owner for an already active
	     * point phase.  The point phase consumes its resulting presentation;
	     * it cannot pause the only transaction able to produce that frame. */
	    else if (enabled(CAPACITY_PLANNING))
		expected = Successor::CALIBRATE_CAPACITY;
	    else if (!enabled(POINT_TRIANGLE_RECOVERY) &&
		    enabled(POINT_CALIBRATION))
		expected = Successor::PRESENT_POINT_CALIBRATION;
	}
	const Successor actual = Policy::completedPassSuccessor(
	    enabled(PASS_COMPLETED), enabled(RESIDENCY_DRAIN),
	    enabled(RESIDENT_GROWTH), enabled(RESIDENT_REFINEMENT),
	    enabled(RESIDENT_WORK), enabled(POINT_TRIANGLE_RECOVERY),
	    enabled(POINT_CALIBRATION), enabled(CAPACITY_PLANNING));
	if (actual != expected) {
	    std::fprintf(stderr,
		"FAIL: completed pass successor combination=%u expected=%u "
		"actual=%u\n", combination,
		static_cast<unsigned int>(expected),
		static_cast<unsigned int>(actual));
	    return 1;
	}
    }
    return 0;
}

static BObolLodCapacityEvidence::Decision
apply_admission_plan(BObolLodCapacityEvidence &evidence,
	BObolLodAdmissionCursor &cursor,
    const BObolLodCapacityEvidence::Inputs &inputs)
{
    BObolLodAdmissionEvidence allEvidence(evidence);
    const BObolLodAdmissionPlan plan =
	BObolLodAdmissionPlanner::plan(allEvidence, cursor,
	    admission_revision_stamp(), inputs);
    evidence = plan.nextEvidence.capacity();
    cursor = plan.nextCursor;
    return plan.capacityDecision;
}

static void
complete_applied_allocation_plan(BObolLodCapacityEvidence &evidence,
    BObolLodAdmissionCursor &cursor,
    const BObolLodCapacityEvidence::CompletedAllocationInputs &inputs)
{
    BObolLodAdmissionEvidence allEvidence(evidence);
    const BObolLodAdmissionPlan plan =
	BObolLodAdmissionPlanner::completeAppliedAllocation(
	    allEvidence, cursor, inputs);
    evidence = plan.nextEvidence.capacity();
    cursor = plan.nextCursor;
}

static void
apply_capacity_evidence_action(BObolLodCapacityEvidence &evidence,
    BObolLodAdmissionCursor &cursor,
    BObolLodAdmissionPlanner::EvidenceAction action)
{
    const BObolLodAdmissionPlan plan =
	BObolLodAdmissionPlanner::applyEvidenceAction(
	    BObolLodAdmissionEvidence(evidence), cursor, action);
    evidence = plan.nextEvidence.capacity();
    cursor = plan.nextCursor;
}

static void
accept_static_presentation_constraint(BObolLodCapacityEvidence &evidence,
    BObolLodAdmissionCursor &cursor, size_t budget)
{
    const BObolLodAdmissionPlan plan =
	BObolLodAdmissionPlanner::acceptStaticPresentationConstraint(
	    BObolLodAdmissionEvidence(evidence), cursor, budget);
    evidence = plan.nextEvidence.capacity();
    cursor = plan.nextCursor;
}

static BObolLodCapacityEvidence::CalibrationDecision
finish_blocked_capacity_plan(BObolLodCapacityEvidence &evidence,
    BObolLodAdmissionCursor &cursor,
    const BObolLodCapacityEvidence::CalibrationInputs &inputs)
{
    BObolLodAdmissionEvidence allEvidence(evidence);
    const BObolLodAdmissionPlan plan =
	BObolLodAdmissionPlanner::finishBlockedCapacityPass(
	    allEvidence, cursor, inputs);
    evidence = plan.nextEvidence.capacity();
    cursor = plan.nextCursor;
    return plan.calibrationDecision;
}

static bool
resume_capacity_candidate_plan(BObolLodCapacityEvidence &evidence,
    BObolLodAdmissionCursor &cursor)
{
    BObolLodAdmissionEvidence allEvidence(evidence);
    const BObolLodAdmissionPlan plan =
	BObolLodAdmissionPlanner::resumeCapacityCandidateAllocation(
	    allEvidence, cursor);
    evidence = plan.nextEvidence.capacity();
    cursor = plan.nextCursor;
    return plan.transitionChanged;
}

static BObolLodCapacityEvidence::CompletedFrameDecision
complete_capacity_frame_plan(BObolLodCapacityEvidence &evidence,
    BObolLodAdmissionCursor &cursor)
{
    BObolLodAdmissionEvidence allEvidence(evidence);
    const BObolLodAdmissionPlan plan =
	BObolLodAdmissionPlanner::completeCapacityCalibrationFrame(
	    allEvidence, cursor);
    evidence = plan.nextEvidence.capacity();
    cursor = plan.nextCursor;
    return plan.completedFrameDecision;
}

static BObolLodStructuralAdmissionEvidence::Decision
apply_structural_admission_plan(
    BObolLodStructuralAdmissionEvidence &evidence,
    const BObolLodStructuralAdmissionEvidence::Inputs &inputs)
{
    BObolLodAdmissionEvidence allEvidence(evidence);
    const BObolLodAdmissionPlan plan =
	BObolLodAdmissionPlanner::planStructural(allEvidence, BObolLodAdmissionCursor(), inputs);
    evidence = plan.nextEvidence.structural();
    return plan.structuralDecision;
}

static BObolLodPointProxyEvidence::Decision
apply_point_interrupted_plan(BObolLodPointProxyEvidence &evidence,
    float currentThreshold, uint64_t renderNanoseconds, float targetFps)
{
    BObolLodAdmissionEvidence allEvidence(evidence);
    const BObolLodAdmissionPlan plan =
	BObolLodAdmissionPlanner::planPointInterrupted(allEvidence, BObolLodAdmissionCursor(),
	    BObolLodPointCalibrationGoal::RESPONSIVE,
	    currentThreshold, renderNanoseconds, targetFps);
    evidence = plan.nextEvidence.pointProxy();
    return plan.pointProxyDecision;
}

static BObolLodPointProxyEvidence::Decision
apply_point_completed_plan(BObolLodPointProxyEvidence &evidence,
    float currentThreshold, uint64_t renderNanoseconds, float targetFps,
    bool reusableSample, size_t unresolvedStructuralCount = 0,
    bool allowCoarsening = true)
{
    BObolLodAdmissionEvidence allEvidence(evidence);
    const BObolLodAdmissionPlan plan =
	BObolLodAdmissionPlanner::planPointCompleted(allEvidence, BObolLodAdmissionCursor(),
	    BObolLodPointCalibrationGoal::RESPONSIVE,
	    currentThreshold, renderNanoseconds, targetFps, reusableSample,
	    allowCoarsening,
	    unresolvedStructuralCount);
    evidence = plan.nextEvidence.pointProxy();
    return plan.pointProxyDecision;
}

static BObolLodPointProxyEvidence::Decision
apply_point_structural_distribution_plan(
    BObolLodPointProxyEvidence &evidence, float currentThreshold,
    const std::array<size_t, 7> &cumulativeCount, size_t visibleCount,
    size_t maximumUncollapsedCount)
{
    BObolLodAdmissionEvidence allEvidence(evidence);
    const BObolLodAdmissionPlan plan =
	BObolLodAdmissionPlanner::planPointStructuralDistribution(allEvidence, BObolLodAdmissionCursor(),
	    currentThreshold, cumulativeCount, visibleCount,
	    maximumUncollapsedCount);
    evidence = plan.nextEvidence.pointProxy();
    return plan.pointProxyDecision;
}

static BObolLodPointProxyEvidence::Decision
apply_point_structural_blocked_plan(BObolLodPointProxyEvidence &evidence,
    float currentThreshold, size_t unresolvedStructuralCount)
{
    BObolLodAdmissionEvidence allEvidence(evidence);
    const BObolLodAdmissionPlan plan =
	BObolLodAdmissionPlanner::planPointStructuralCoverageBlocked(allEvidence, BObolLodAdmissionCursor(),
	    currentThreshold, unresolvedStructuralCount);
    evidence = plan.nextEvidence.pointProxy();
    return plan.pointProxyDecision;
}

static int
test_point_structural_preload_threshold(void)
{
    using Planner = BObolLodAdmissionPlanner;
    struct ThresholdCase {
	float current;
	float candidate;
	float expected;
    };
    const ThresholdCase cases[] = {
	{1.0f, 1.25f, 1.25f},
	{1.5f, 1.75f, 1.75f},
	{1.0f, 2.0f, 2.0f},
	{4.0f, 3.5f, 2.0f},
	{2.0f, 1.99f, 1.0f},
	{8.0f, 4.0f, 4.0f},
	{1.0f, 128.0f, 64.0f},
	{4.0f, 0.25f, 1.0f}
    };
    for (const ThresholdCase &entry : cases) {
	const float actual = Planner::pointStructuralPreloadThreshold(
	    entry.current, entry.candidate);
	if (std::fabs(actual - entry.expected) > 0.001f) {
	    std::fprintf(stderr,
		"FAIL: point structural preload threshold "
		"current=%.3f candidate=%.3f expected=%.3f actual=%.3f\n",
		entry.current, entry.candidate, entry.expected, actual);
	    return 1;
	}
    }

    static const int thresholdQuantaPerPixel = 8;
    static const int maximumThresholdPixels = 64;
    const int maximumQuanta =
	thresholdQuantaPerPixel * maximumThresholdPixels;
    for (int currentQuanta = thresholdQuantaPerPixel;
	 currentQuanta <= maximumQuanta; ++currentQuanta) {
	const float current = static_cast<float>(currentQuanta) /
	    thresholdQuantaPerPixel;
	for (int candidateQuanta = thresholdQuantaPerPixel;
	     candidateQuanta <= maximumQuanta; ++candidateQuanta) {
	    if (candidateQuanta == currentQuanta)
		continue;
	    const float candidate = static_cast<float>(candidateQuanta) /
		thresholdQuantaPerPixel;
	    const float mapped = Planner::pointStructuralPreloadThreshold(
		current, candidate);
	    const bool progresses = candidate > current ? mapped > current :
		mapped < current;
	    if (!progresses) {
		std::fprintf(stderr,
		    "FAIL: changed point candidate did not change the "
		    "presented threshold current=%.3f candidate=%.3f "
		    "mapped=%.3f\n", current, candidate, mapped);
		return 1;
	    }
	}
    }
    if (!Planner::pointRelaxationFitsDeadline(
	    686, 1208, 13452000, 200000000) ||
	!Planner::pointRelaxationFitsDeadline(
	    6650, 7884, 72618000, 200000000) ||
	!Planner::pointRelaxationFitsDeadline(
	    14534, 19589, 102691000, 400000000) ||
	Planner::pointRelaxationFitsDeadline(
	    14534, 19589, 102691000, 200000000) ||
	!Planner::pointRelaxationFitsDeadline(
	    14534, 0, 102691000, 200000000) ||
	Planner::pointRelaxationFitsDeadline(
	    0, 1, 102691000, 200000000) ||
	Planner::pointRelaxationFitsDeadline(
	    1, 1, 0, 200000000) ||
	!Planner::pointRelaxationFitsDeadline(
	    1, 1, 102691000, 0)) {
	std::fprintf(stderr,
	    "FAIL: point relaxation preload deadline projection\n");
	return 1;
    }
    return 0;
}

static int
test_retained_allocation_request(void)
{
    using Request = BObolLodRetainedAllocationRequest;

    Request request;
    if (request.pending() || request.kind() != Request::Kind::NONE ||
	request.preservesBudget() || request.capacityCandidate() ||
	request.reconcilesPresentation() ||
	request.reconciliationBudget() != SIZE_MAX) {
	std::fprintf(stderr, "FAIL: retained allocation request default\n");
	return 1;
    }

    request.requestReallocation(true);
    if (!request.pending() || !request.preservesBudget() ||
	request.kind() != Request::Kind::PRESERVE_BUDGET) {
	std::fprintf(stderr, "FAIL: retained allocation preserve request\n");
	return 1;
    }

    request.requestReallocation(false);
    if (!request.pending() || request.preservesBudget() ||
	request.kind() != Request::Kind::RECOMPUTE_BUDGET) {
	std::fprintf(stderr, "FAIL: retained allocation recompute request\n");
	return 1;
    }

    request.requestCapacityCandidateAllocation();
    request.requestReallocation(false);
    request.requestPresentationReconciliation(120000);
    if (!request.pending() || !request.capacityCandidate() ||
	!request.preservesBudget() || request.reconcilesPresentation() ||
	request.kind() != Request::Kind::CAPACITY_CANDIDATE) {
	std::fprintf(stderr,
	    "FAIL: weaker request displaced capacity candidate\n");
	return 1;
    }
    request.retireCapacityCandidate();
    if (request.pending()) {
	std::fprintf(stderr, "FAIL: capacity candidate did not retire\n");
	return 1;
    }

    request.reset();
    request.requestPresentationReconciliation(120000);
    request.requestPresentationReconciliation(180000);
    request.requestReallocation(true);
    if (!request.pending() || request.preservesBudget() ||
	!request.reconcilesPresentation() ||
	request.kind() != Request::Kind::PRESENTATION_RECONCILIATION ||
	request.reconciliationBudget() != 120000) {
	std::fprintf(stderr,
	    "FAIL: weaker request displaced presentation reconciliation\n");
	return 1;
    }

    request.requestPresentationReconciliation(90000);
    if (request.reconciliationBudget() != 90000) {
	std::fprintf(stderr,
	    "FAIL: presentation reconciliation did not retain strict budget\n");
	return 1;
    }

    request.reset();
    request.requestPresentationReconciliation(0);
    request.requestPresentationReconciliation(SIZE_MAX);
    if (request.pending() || request.reconciliationBudget() != SIZE_MAX) {
	std::fprintf(stderr,
	    "FAIL: invalid presentation reconciliation became pending\n");
	return 1;
    }

    return 0;
}

static int
test_admission_capacity(void)
{
    using Policy = BObolLodCapacityEvidence;
    Policy policy;
    Policy::Inputs input;
    input.activeCost = 1234;
    input.minimumActiveCost = 321;
    BObolLodAdmissionEvidence initialEvidence(policy);
    BObolLodAdmissionCursor initialCursor;
    initialEvidence = BObolLodAdmissionPlanner::planResourceObservation(
	initialEvidence, initialCursor, true, false, true).nextEvidence;
    initialEvidence = BObolLodAdmissionPlanner::planHeadroomRetry(
	initialEvidence, initialCursor, BObolLodViewEpoch(2), BObolLodPolicyEpoch(3),
	1000, 100).nextEvidence;
    BObolLodStructuralAdmissionEvidence::Inputs structuralInputs;
    structuralInputs.viewRevision = 2;
    structuralInputs.policyRevision = 3;
    structuralInputs.frontierIdentity = 4;
    structuralInputs.unresolvedCount = 1;
    structuralInputs.activeCost = 100;
    structuralInputs.certifiedBudget = 200;
    structuralInputs.exactProjection = true;
    structuralInputs.pointPolicyExhausted = true;
    initialEvidence = BObolLodAdmissionPlanner::planStructural(
	initialEvidence, initialCursor, structuralInputs).nextEvidence;
    const BObolLodAdmissionPlan firstPlan =
	BObolLodAdmissionPlanner::plan(initialEvidence, initialCursor,
	    admission_revision_stamp(), input);
    const BObolLodAdmissionPlan repeatedPlan =
	BObolLodAdmissionPlanner::plan(initialEvidence, initialCursor,
	    admission_revision_stamp(), input);
    const BObolLodAdmissionPlan sameRevisionPlan =
	BObolLodAdmissionPlanner::plan(firstPlan.nextEvidence,
	    firstPlan.nextCursor, admission_revision_stamp(), input);
    const BObolLodAdmissionRevisionStamp nextViewStamp =
	admission_revision_stamp(1, 1, 1, 2, 1, 1);
    const BObolLodAdmissionPlan supersedingRevisionPlan =
	BObolLodAdmissionPlanner::plan(firstPlan.nextEvidence,
	    firstPlan.nextCursor, nextViewStamp, input);
    const BObolLodAdmissionPlan unrelatedTransition =
	BObolLodAdmissionPlanner::planResourceObservation(
	    firstPlan.nextEvidence, firstPlan.nextCursor, true, false, true);
    const BObolLodAdmissionPlan reallocationRequest =
	BObolLodAdmissionPlanner::requestRetainedReallocation(
	    firstPlan.nextEvidence, firstPlan.nextCursor);
    const BObolLodAdmissionPlan reallocationPlan =
	BObolLodAdmissionPlanner::plan(reallocationRequest.nextEvidence,
	    reallocationRequest.nextCursor, admission_revision_stamp(), input);
    if (initialCursor.initialized() || policy.currentBudget() != 50000 ||
	firstPlan.capacityDecision.initialized !=
	    repeatedPlan.capacityDecision.initialized ||
	firstPlan.capacityDecision.totalBudget !=
	    repeatedPlan.capacityDecision.totalBudget ||
	firstPlan.capacityDecision.refinementBudget !=
	    repeatedPlan.capacityDecision.refinementBudget ||
	firstPlan.nextEvidence.capacity().currentBudget() !=
	    repeatedPlan.nextEvidence.capacity().currentBudget() ||
	firstPlan.nextCursor.activeCost() !=
	    repeatedPlan.nextCursor.activeCost() ||
	!firstPlan.certifiedFor(admission_revision_stamp()) ||
	firstPlan.certifiedFor(nextViewStamp) ||
	sameRevisionPlan.capacityDecision.initialized ||
	!sameRevisionPlan.nextCursor.matches(admission_revision_stamp()) ||
	!supersedingRevisionPlan.capacityDecision.initialized ||
	!supersedingRevisionPlan.nextCursor.matches(nextViewStamp) ||
	!supersedingRevisionPlan.certifiedFor(nextViewStamp) ||
	unrelatedTransition.nextCursor.activeCost() !=
	    firstPlan.nextCursor.activeCost() ||
	unrelatedTransition.nextCursor.refinementRemaining() !=
	    firstPlan.nextCursor.refinementRemaining() ||
	unrelatedTransition.nextEvidence.capacity().currentBudget() !=
	    firstPlan.nextEvidence.capacity().currentBudget() ||
	reallocationRequest.nextCursor.initialized() ||
	!reallocationPlan.capacityDecision.retainedAdmission ||
	!reallocationPlan.nextCursor.retainedAdmission() ||
	!firstPlan.nextEvidence.resources().cpuPressure() ||
	firstPlan.nextEvidence.resources().revision() != 1 ||
	!firstPlan.nextEvidence.headroom().pendingMatches(
	    BObolLodViewEpoch(2), BObolLodPolicyEpoch(3), 1000) ||
	!BObolLodAdmissionPlanner::planStructural(
	    firstPlan.nextEvidence, firstPlan.nextCursor, structuralInputs).
	    structuralDecision.duplicateRejected) {
	std::fprintf(stderr,
	    "FAIL: admission planner was mutable, nondeterministic, or lost "
	    "unrelated evidence\n");
	return 1;
    }
    const BObolLodAdmissionPlan reconciliationPlan =
	BObolLodAdmissionPlanner::requestPresentationReconciliation(
	    firstPlan.nextEvidence, firstPlan.nextCursor, 1200);
    const BObolLodAdmissionPlan duplicateReconciliationPlan =
	BObolLodAdmissionPlanner::requestPresentationReconciliation(
	    reconciliationPlan.nextEvidence, reconciliationPlan.nextCursor,
	    1800);
    if (!reconciliationPlan.transitionChanged ||
	reconciliationPlan.nextCursor.initialized() ||
	reconciliationPlan.nextEvidence.capacity().currentBudget() != 1200 ||
	duplicateReconciliationPlan.transitionChanged ||
	duplicateReconciliationPlan.nextEvidence.capacity().currentBudget() !=
	    1200) {
	std::fprintf(stderr,
	    "FAIL: presentation reconciliation was not an atomic cursor/budget "
	    "transition\n");
	return 1;
    }
    policy = firstPlan.nextEvidence.capacity();
    BObolLodAdmissionCursor cursor = firstPlan.nextCursor;
    Policy::Decision decision = firstPlan.capacityDecision;
    if (!decision.initialized || decision.totalBudget != 50000 ||
	decision.refinementBudget != 48766 || decision.retainedAdmission ||
	!cursor.initialized() || policy.currentBudget() != 50000 ||
	cursor.activeCost() != 1234 ||
	cursor.minimumActiveCost() != 321) {
	std::fprintf(stderr, "FAIL: seed scene budget\n");
	return 1;
    }
    input.activeCost = 9999;
    input.minimumActiveCost = 8888;
    decision = apply_admission_plan(policy, cursor, input);
	if (decision.initialized || decision.totalBudget != 50000 ||
	cursor.activeCost() != 1234 ||
	cursor.minimumActiveCost() != 321) {
	std::fprintf(stderr, "FAIL: duplicate scene budget initialization\n");
	return 1;
    }
    cursor.consumeRefinement(12000);
    if (cursor.refinementRemaining() != 36766) {
	std::fprintf(stderr, "FAIL: refinement budget consumption\n");
	return 1;
    }

    cursor.reset();
    if (cursor.activeCost() != 0 ||
	cursor.minimumActiveCost() != 0) {
	std::fprintf(stderr, "FAIL: bounded-pass cost snapshot reset\n");
	return 1;
    }
    input.activeCost = 100000;
    input.minimumActiveCost = 0;
    input.targetFps = 20.0f;
    input.calibratedCostPerSecond = 10000000.0L;
    input.observedStableNanoseconds = 10000000ULL;
    decision = apply_admission_plan(policy, cursor, input);
    if (!decision.initialized || decision.totalBudget != 400000 ||
	decision.refinementBudget != 300000 ||
	decision.retainedAdmission) {
	std::fprintf(stderr, "FAIL: bounded fast-frame budget growth\n");
	return 1;
    }

    /* An exact reusable frame is stronger evidence than a linear throughput
     * estimate.  Many-instance dispatch/classification work is not encoded in
     * triangle cost, so the estimator may place a population below a budget
     * even though that exact population was just presented within deadline. */
    policy.reset();
    cursor.reset();
    input = Policy::Inputs();
    input.activeCost = 2070262;
    input.targetFps = 20.0f;
    input.calibratedCostPerSecond = 6300000.0L;
    input.observedStableNanoseconds = 37000000ULL;
    decision = apply_admission_plan(policy, cursor, input);
    if (!decision.initialized || decision.totalBudget != input.activeCost ||
	decision.refinementBudget != 0 || decision.retainedAdmission) {
	std::fprintf(stderr,
	    "FAIL: exact deadline-safe population was placed below budget "
	    "budget=%zu active=%zu admission=%d\n",
	    decision.totalBudget, input.activeCost,
	    decision.retainedAdmission ? 1 : 0);
	return 1;
    }

    /* A bounded static presentation uses the event-driven hard deadline
     * rather than the many-object growth ladder.  This is used both for a
     * single-occurrence quality step and for terminal structural coverage
     * repair; an endpoint abort supplies recovery if the prediction is
     * optimistic. */
    policy.reset();
    cursor.reset();
    input = Policy::Inputs();
    input.activeCost = 100000;
    input.targetFps = 10.0f;
    input.calibratedCostPerSecond = 10000000.0L;
    input.observedStableNanoseconds = 40000000ULL;
    input.hardDeadlinePresentation = true;
    decision = apply_admission_plan(policy, cursor, input);
    if (!decision.initialized || decision.totalBudget != 1000000 ||
	decision.refinementBudget != 900000 ||
	decision.retainedAdmission) {
	std::fprintf(stderr,
	    "FAIL: hard-deadline static presentation budget\n");
	return 1;
    }

    /* Replacing the last boxes in a many-occurrence static view is a terminal
     * coverage obligation, not ordinary 20 Hz quality refinement.  It must be
     * able to spend the complete calibrated hard-frame allowance in one pass
     * instead of repeatedly probing an already saturated preferred budget. */
    policy.reset();
    cursor.reset();
    input = Policy::Inputs();
    input.activeCost = 294000;
    input.targetFps = 10.0f;
    input.calibratedCostPerSecond = 6300000.0L;
    input.observedStableNanoseconds = 50000000ULL;
    input.hardDeadlinePresentation = true;
    decision = apply_admission_plan(policy, cursor, input);
    if (!decision.initialized || decision.totalBudget != 630000 ||
	decision.refinementBudget != 336000 ||
	decision.retainedAdmission) {
	std::fprintf(stderr,
	    "FAIL: many-occurrence terminal coverage repair budget\n");
	return 1;
    }

    /* The exact current population, rather than an older throughput EMA,
     * supplies bounded static headroom.  This is what lets terminal coverage
     * replace the last structural proxies after a large publication wave. */
    policy.reset();
    cursor.reset();
    input = Policy::Inputs();
    input.activeCost = 100000;
    input.targetFps = 10.0f;
    input.calibratedCostPerSecond = 1000000.0L;
    input.observedStableNanoseconds = 25000000ULL;
    input.hardDeadlinePresentation = true;
    decision = apply_admission_plan(policy, cursor, input);
    if (!decision.initialized || decision.totalBudget != 320000 ||
	decision.refinementBudget != 220000 ||
	decision.retainedAdmission) {
	std::fprintf(stderr,
	    "FAIL: exact static frame did not supply direct headroom "
	    "budget=%zu refinement=%zu admission=%d\n",
	    decision.totalBudget, decision.refinementBudget,
	    decision.retainedAdmission ? 1 : 0);
	return 1;
    }

    /* Once coverage proves a single visible occurrence, an empty retained
     * population may start at a useful medium-mesh allowance.  This removes
     * several warm-cache blocky frames without changing the conservative seed
     * for an unknown/many-object population. */
    policy.reset();
    cursor.reset();
    input = Policy::Inputs();
    input.coldSingleOccurrence = true;
    decision = apply_admission_plan(policy, cursor, input);
    if (!decision.initialized ||
	decision.totalBudget != Policy::singleOccurrenceBootstrapBudget() ||
	decision.refinementBudget !=
	    Policy::singleOccurrenceBootstrapBudget()) {
	std::fprintf(stderr, "FAIL: single-occurrence bootstrap budget\n");
	return 1;
    }

    /* A hard-deadline handoff may replace a stale, optimistic allowance, but
     * its renderer-only ceiling must not rewrite the retained occurrence
     * population.  An exact over-budget frame or an explicit recovery
     * request supplies that separate authority. */
    policy.reset();
    cursor.reset();
    input = Policy::Inputs();
    input.activeCost = 400000;
    input.targetFps = 20.0f;
    input.calibratedCostPerSecond = 2000000.0L;
    input.observedStableNanoseconds = 40000000ULL;
    input.stablePresentationHandoff = true;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget != 100000 ||
	decision.retainedAdmission ||
	decision.retainedAdmissionBudget != SIZE_MAX) {
	std::fprintf(stderr, "FAIL: deadline handoff budget reconciliation\n");
	return 1;
    }

    policy.raiseCurrentBudget(400000);
    policy.resetOverloadRecovery();
    policy.invalidateCalibration();
    input = Policy::Inputs();
    cursor.reset();
    input.activeCost = 400000;
    input.targetFps = 20.0f;
    input.calibratedCostPerSecond = 4000000.0L;
    input.observedStableNanoseconds = 70000000ULL;
    policy.noteDeadlineCapacityMiss(400000);
    decision = apply_admission_plan(policy, cursor, input);
    if (!decision.overloadRecovery || decision.totalBudget != 200000 ||
	!decision.retainedAdmission ||
	decision.retainedAdmissionBudget != 200000 ||
	!policy.overloadRecoveryPerformed() ||
	policy.overloadRecoveryActiveCost() != 400000) {
	std::fprintf(stderr,
	    "FAIL: stable overload recovery decision=%d budget=%zu "
	    "admission=%d admission_budget=%zu performed=%d active=%zu\n",
	    decision.overloadRecovery ? 1 : 0, decision.totalBudget,
	    decision.retainedAdmission ? 1 : 0,
	    decision.retainedAdmissionBudget,
	    policy.overloadRecoveryPerformed() ? 1 : 0,
	    policy.overloadRecoveryActiveCost());
	return 1;
    }

    /* The quiet preference is not a hard failure threshold.  A static frame
     * below the independently enforced presentation deadline keeps its
     * proven population across a pose-only interaction; reducing it here can
     * expose structural boxes which were absent before the rotation. */
    policy.reset();
    cursor.reset();
    policy.raiseCurrentBudget(400000);
    input = Policy::Inputs();
    input.activeCost = 400000;
    input.targetFps = 20.0f;
    input.calibratedCostPerSecond = 4000000.0L;
    input.observedStableNanoseconds = 70000000ULL;
    input.hardPresentationDeadlineNanoseconds = 100000000ULL;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.overloadRecovery || decision.totalBudget != 400000 ||
	decision.retainedAdmission) {
	std::fprintf(stderr,
	    "FAIL: hard-deadline-safe static population was coarsened "
	    "recovery=%d budget=%zu admission=%d\n",
	    decision.overloadRecovery ? 1 : 0, decision.totalBudget,
	    decision.retainedAdmission ? 1 : 0);
	return 1;
    }

    /* A many-occurrence fixed cost is not described by the triangle-rate
	 * calibration.  A typed deadline miss must therefore authorize the direct
	 * duration-based recovery even when that estimator predicts the active cut
	 * should fit. */
    policy.reset();
    cursor.reset();
    policy.raiseCurrentBudget(400000);
    cursor.reset();
    input = Policy::Inputs();
    input.activeCost = 400000;
    input.targetFps = 20.0f;
    input.calibratedCostPerSecond = 20000000.0L;
    input.observedStableNanoseconds = 80000000ULL;
    policy.noteDeadlineCapacityMiss(400000);
    decision = apply_admission_plan(policy, cursor, input);
    if (!decision.overloadRecovery || decision.totalBudget != 200000 ||
	!decision.retainedAdmission) {
	std::fprintf(stderr,
	    "FAIL: repeated exact deadline misses did not override the "
	    "throughput model recovery=%d budget=%zu admission=%d\n",
	    decision.overloadRecovery ? 1 : 0, decision.totalBudget,
	    decision.retainedAdmission ? 1 : 0);
	return 1;
    }
    cursor.consumeRetainedAdmission(75000);
    if (cursor.retainedAdmissionRemaining() != 125000) {
	std::fprintf(stderr, "FAIL: retained admission consumption\n");
	return 1;
    }

    /* Minimum mesh coverage is pre-reserved.  Bounded windows carry only the
     * upgrade allowance so every slice cannot charge the same floor and
     * normalize the whole scene before refinement restarts. */
    policy.reset();
    cursor.reset();
    policy.requestRetainedRecovery(200000);
    input = Policy::Inputs();
    input.activeCost = 400000;
    input.minimumActiveCost = 125000;
    input.targetFps = 20.0f;
    input.calibratedCostPerSecond = 4000000.0L;
    decision = apply_admission_plan(policy, cursor, input);
    if (!decision.retainedAdmission ||
	decision.totalBudget != 200000 ||
	decision.retainedAdmissionBudget != 75000) {
	std::fprintf(stderr,
	    "FAIL: retained admission did not reserve coverage floor "
	    "budget=%zu upgrade=%zu\n",
	    decision.totalBudget, decision.retainedAdmissionBudget);
	return 1;
    }

    /* Returning a temporary point-proxy cut to one pixel explicitly
     * compacts reducible PoP detail before whole visible objects may be
     * aggregated.  This recovery must not depend on a slow-frame probe: the
     * point cut which hid the population can itself be fast. */
    policy.reset();
    cursor.reset();
    policy.requestRetainedRecovery(120000);
    if (!policy.retainedRecoveryCeilingActive()) {
	std::fprintf(stderr, "FAIL: measured recovery ceiling provenance\n");
	return 1;
    }
    cursor.reset();
    input = Policy::Inputs();
    input.activeCost = 400000;
    input.targetFps = 20.0f;
    input.calibratedCostPerSecond = 20000000.0L;
    input.observedStableNanoseconds = 10000000ULL;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget != 120000 ||
	!decision.retainedAdmission ||
	decision.retainedAdmissionBudget != 120000) {
	std::fprintf(stderr,
	    "FAIL: point-to-triangle recovery budget=%zu admission=%d "
	    "admission_budget=%zu\n",
	    decision.totalBudget, decision.retainedAdmission ? 1 : 0,
	    decision.retainedAdmissionBudget);
	return 1;
    }

    /* The measured recovery ceiling belongs to the current view/policy
     * epoch, not just the first admission pass.  A cheap coarse frame must
     * not immediately re-admit the rejected discrete cut. */
    cursor.reset();
    input.activeCost = 100000;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget != 120000) {
	std::fprintf(stderr,
	    "FAIL: retained recovery ceiling did not persist budget=%zu\n",
	    decision.totalBudget);
	return 1;
    }
    policy.raiseCurrentBudget(400000);
    cursor.reset();
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget != 120000) {
	std::fprintf(stderr,
	    "FAIL: retained recovery ceiling did not bound raise budget=%zu\n",
	    decision.totalBudget);
	return 1;
    }
    if (policy.confirmRetainedRecoveryPresentation(false, cursor) ||
	!policy.retainedRecoveryCeilingActive()) {
	std::fprintf(stderr,
	    "FAIL: unconfirmed one-pixel recovery retired its ceiling\n");
	return 1;
    }

    /* A persistent recovery ceiling protects richer triangle prefixes, but
     * cannot strand visible structural boxes.  Terminal coverage repair uses
     * the separately deadline-bounded static allowance for minimum meshes and
     * leaves the ceiling installed for subsequent detail allocation. */
    cursor.reset();
    input = Policy::Inputs();
    input.activeCost = 100000;
    input.targetFps = 10.0f;
    input.calibratedCostPerSecond = 10000000.0L;
    input.observedStableNanoseconds = 40000000ULL;
    input.hardDeadlinePresentation = true;
    input.structuralCoverageRepair = true;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget != 1000000 ||
	decision.refinementBudget != 900000 ||
	decision.retainedAdmission ||
	!policy.retainedRecoveryCeilingActive()) {
	std::fprintf(stderr,
	    "FAIL: structural coverage remained recovery limited "
	    "budget=%zu refinement=%zu admission=%d ceiling=%d\n",
	    decision.totalBudget, decision.refinementBudget,
	    decision.retainedAdmission ? 1 : 0,
	    policy.retainedRecoveryCeilingActive() ? 1 : 0);
	return 1;
    }
    if (!policy.confirmRetainedRecoveryPresentation(true, cursor) ||
	policy.retainedRecoveryCeilingActive() ||
	cursor.initialized()) {
	std::fprintf(stderr,
	    "FAIL: confirmed one-pixel recovery retained policy state\n");
	return 1;
    }
    policy.raiseCurrentBudget(400000);
    cursor.reset();
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget <= 120000) {
	std::fprintf(stderr,
	    "FAIL: retired recovery ceiling still constrained budget=%zu\n",
	    decision.totalBudget);
	return 1;
    }

    /* A handoff normalization is intentionally one-pass.  Its successor may
     * grow from the newly measured coherent baseline. */
    policy.reset();
    cursor.reset();
    policy.requestRetainedNormalization(120000);
    if (policy.retainedRecoveryCeilingActive()) {
	std::fprintf(stderr, "FAIL: normalization acquired recovery ceiling\n");
	return 1;
    }
    cursor.reset();
    input.activeCost = 400000;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget != 120000 || !decision.retainedAdmission) {
	std::fprintf(stderr,
	    "FAIL: retained normalization budget=%zu admission=%d\n",
	    decision.totalBudget, decision.retainedAdmission ? 1 : 0);
	return 1;
    }
    cursor.reset();
    input.activeCost = 100000;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget <= 120000) {
	std::fprintf(stderr,
	    "FAIL: retained normalization incorrectly persisted budget=%zu\n",
	    decision.totalBudget);
	return 1;
    }

    /* A changed quiet camera pose does not imply overload: it may leave the
     * same total cost affordable while moving projected importance between
     * occurrences.  Authorize one reallocation at the retained budget, then
     * consume the authority so an unchanged view cannot churn its cuts. */
    policy.reset();
    cursor.reset();
    policy.raiseCurrentBudget(200000);
    policy.requestRetainedReallocation();
    input = Policy::Inputs();
    input.activeCost = 150000;
    input.minimumActiveCost = 50000;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget != 200000 ||
	!decision.retainedAdmission ||
	decision.retainedAdmissionBudget != 150000) {
	std::fprintf(stderr,
	    "FAIL: retained importance reallocation budget=%zu "
	    "admission=%d upgrade=%zu\n",
	    decision.totalBudget, decision.retainedAdmission ? 1 : 0,
	    decision.retainedAdmissionBudget);
	return 1;
    }
    cursor.reset();
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.retainedAdmission ||
	decision.retainedAdmissionBudget != SIZE_MAX) {
	std::fprintf(stderr,
	    "FAIL: retained importance reallocation was not one-shot\n");
	return 1;
    }

    /* A global presentation ceiling can make the exact submitted frame look
     * cheaper than the richer retained population.  Its hard-deadline proof
     * must nevertheless force one full allocation at the demonstrated
     * allowance, overriding both a stale low seed and a soft protected floor. */
    policy.reset();
    cursor.reset();
    policy.raiseCurrentBudget(200000);
    policy.setRetainedQualityFloorBudget(cursor,
	300000, 0x2222222222222222ULL, 80000, 50000);
    policy.requestPresentationReconciliation(120000);
    policy.requestPresentationReconciliation(180000);
    if (policy.currentBudget() != 120000) {
	std::fprintf(stderr,
	    "FAIL: weaker handoff request replaced canonical budget=%zu\n",
	    policy.currentBudget());
	return 1;
    }
    input = Policy::Inputs();
    input.activeCost = 80000;
    input.minimumActiveCost = 50000;
    input.stablePresentationHandoff = true;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget != 120000 || !decision.retainedAdmission ||
	decision.retainedAdmissionBudget != 70000 ||
	policy.currentBudget() != 120000) {
	std::fprintf(stderr,
	    "FAIL: presentation reconciliation budget=%zu admission=%d "
	    "upgrade=%zu current=%zu\n",
	    decision.totalBudget, decision.retainedAdmission ? 1 : 0,
	    decision.retainedAdmissionBudget, policy.currentBudget());
	return 1;
    }

    /* An exact-frame structural coverage allowance bypasses the soft
     * preferred-cadence estimator for one repair pass, without turning that
     * pass into a retained-cut reallocation or a permanent budget floor. */
    policy.reset();
    cursor.reset();
    policy.raiseCurrentBudget(200000);
    if (policy.requestCoverageCompletion(300000, 450000) != 450000) {
	std::fprintf(stderr,
	    "FAIL: structural coverage budget request was rejected\n");
	return 1;
    }
    input = Policy::Inputs();
    /* Presentation handoff may legitimately change the retained accounting
     * baseline between certification and consumption.  Preserve the proven
     * 150k marginal headroom instead of replaying the stale 450k absolute
     * total. */
    input.activeCost = 330000;
    input.minimumActiveCost = 100000;
    input.structuralCoverageRepair = true;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget != 480000 || decision.retainedAdmission ||
	decision.refinementBudget != 150000) {
	std::fprintf(stderr,
	    "FAIL: structural coverage transaction budget=%zu retained=%d "
	    "refinement=%zu\n", decision.totalBudget,
	    decision.retainedAdmission ? 1 : 0, decision.refinementBudget);
	return 1;
    }
    cursor.reset();
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget == 480000) {
	std::fprintf(stderr,
	    "FAIL: structural coverage transaction was not one-shot\n");
	return 1;
    }

    /* A hard-deadline recovery budget is a one-pass safe handoff, while the
     * failed attempted population is a persistent unsafe bound for that
     * deadline class.  A preferred-cadence miss must not reject a prominent
     * floor which is still eligible for the longer static deadline. */
    policy.reset();
    cursor.reset();
    policy.raiseCurrentBudget(500000);
    policy.setRetainedQualityFloorBudget(cursor,
	480000, 0x4444444444444444ULL, 400000, 100000);
    policy.noteDeadlineCapacityMiss(500000);
    if (policy.deadlineCapacityCeiling() != 499999 ||
	policy.staticDeadlineCapacityCeiling() != SIZE_MAX ||
	policy.currentBudget() != 499999 ||
	!policy.retainedQualityFloorActive() ||
	policy.retainedQualityFloorRejected()) {
	std::fprintf(stderr,
	    "FAIL: steady deadline capacity bound erased static quality\n");
	return 1;
    }

    /* The preferred-cadence miss must not erase the independent static-view
     * quality domain.  A miss at that longer deadline constrains both goals. */
	policy.noteDeadlineCapacityMiss(600000, true);
	if (policy.staticDeadlineCapacityCeiling() != 599999 ||
	policy.deadlineCapacityCeiling() != 499999 ||
	!policy.retainedQualityFloorActive()) {
	std::fprintf(stderr,
	    "FAIL: deadline capacity evidence mixed steady/static domains\n");
	return 1;
	}
	policy.setRetainedQualityFloorBudget(cursor,
	    610000, 0x5555555555555555ULL, 400000, 100000);
	if (policy.retainedQualityFloorActive() ||
	    !policy.retainedQualityFloorRejected()) {
	std::fprintf(stderr,
	    "FAIL: static deadline miss did not reject richer floor\n");
	return 1;
    }
    /* The deadline ceiling is already a strict successor below the failed
     * population.  Marginal recovery may apply its independent throughput
     * margin before this cap, but must not compound the margin here. */
    if (BObolLodAdmissionPlanner::capBudgetAtDeadlineCeiling(
	    600000, policy.deadlineCapacityCeiling()) != 499999 ||
	BObolLodAdmissionPlanner::capBudgetAtDeadlineCeiling(
	    400000, SIZE_MAX) != 400000) {
	std::fprintf(stderr,
	    "FAIL: static marginal budget compounded deadline safety margin\n");
	return 1;
    }
    policy.raiseCurrentBudget(600000);
    if (policy.currentBudget() != 499999) {
	std::fprintf(stderr,
	    "FAIL: ordinary budget raise exceeded deadline capacity bound\n");
	return 1;
    }
    policy.requestPresentationReconciliation(400000);
    input = Policy::Inputs();
    input.activeCost = 350000;
    input.minimumActiveCost = 100000;
    input.stablePresentationHandoff = true;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget != 400000 || !decision.retainedAdmission) {
	std::fprintf(stderr,
	    "FAIL: deadline capacity bound blocked safe reconciliation\n");
	return 1;
    }
    cursor.reset();
    input.stablePresentationHandoff = false;
    input.activeCost = 400000;
    input.calibratedCostPerSecond = 20000000.0L;
    input.targetFps = 20.0f;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget > 499999) {
	std::fprintf(stderr,
	    "FAIL: calibration restored a known failed population\n");
	return 1;
    }
    cursor.reset();
    if (policy.requestCoverageCompletion(400000, 600000) != 600000) {
	std::fprintf(stderr,
	    "FAIL: exact structural coverage was serialized by an old allocation\n");
	return 1;
    }
    input = Policy::Inputs();
    input.activeCost = 400000;
    input.minimumActiveCost = 100000;
    input.structuralCoverageRepair = true;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget != 600000 ||
	decision.refinementBudget != 200000) {
	std::fprintf(stderr,
	    "FAIL: exact structural coverage did not use certified batch budget\n");
	return 1;
    }
    policy.clearDeadlineCapacityCeiling();
    policy.raiseCurrentBudget(600000);
    if (policy.currentBudget() != 600000) {
	std::fprintf(stderr,
	    "FAIL: capacity epoch did not release deadline bound\n");
	return 1;
    }

    /* A hard miss rejects exactly the attempted candidate.  It must not
     * apply an arbitrary percentage margin to the persistent search bound:
     * an adjacent discrete PoP population can occupy that narrow interval
     * and still fit comfortably inside the renderer deadline. */
    policy.reset();
    constexpr size_t failedDiscreteCandidate = 100;
    constexpr size_t adjacentDiscretePopulation = 96;
    policy.noteDeadlineCapacityMiss(failedDiscreteCandidate, true);
    if (policy.staticDeadlineCapacityCeiling() !=
	    failedDiscreteCandidate - 1 ||
	policy.staticDeadlineCapacityCeiling() < adjacentDiscretePopulation) {
	std::fprintf(stderr,
	    "FAIL: deadline evidence skipped an untested discrete population "
	    "(ceiling=%zu)\n",
	    policy.staticDeadlineCapacityCeiling());
	return 1;
    }

    /* A late reusable frame may prove a larger allowance after the original
     * allocation became terminal.  Reallocate the whole population using the
     * newly calibrated capacity rather than preserving the smaller allowance
     * merely because the current cut is already under it. */
    policy.reset();
    cursor.reset();
    policy.raiseCurrentBudget(200000);
    policy.requestRetainedReallocation(false);
    input = Policy::Inputs();
    input.activeCost = 150000;
    input.minimumActiveCost = 50000;
    input.targetFps = 20.0f;
    input.calibratedCostPerSecond = 8000000.0L;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget != 400000 || !decision.retainedAdmission ||
	decision.retainedAdmissionBudget != 350000) {
	std::fprintf(stderr,
	    "FAIL: demonstrated headroom did not reallocate grown budget\n");
	return 1;
    }

    /* A prominent-quality allocation may deliberately use capacity between
     * the preferred stable cadence and the hard quality-frame deadline.  It
     * must remain a fixed point across ordinary 20 Hz overload evaluation;
     * an actual hard-deadline miss explicitly rejects the floor. */
    policy.reset();
    cursor.reset();
    policy.raiseCurrentBudget(200000);
    policy.requestRetainedReallocation();
    input = Policy::Inputs();
    input.activeCost = 150000;
    input.minimumActiveCost = 50000;
    decision = apply_admission_plan(policy, cursor, input);
    policy.setRetainedQualityFloorBudget(cursor,
	300000, 0x1111111111111111ULL, 150000, 50000);
    if (policy.currentBudget() != 300000 ||
	cursor.retainedAdmissionRemaining() != 250000) {
	std::fprintf(stderr, "FAIL: retained quality-floor admission setup\n");
	return 1;
    }
    cursor.reset();
    policy.resetOverloadRecovery();
    input.activeCost = 300000;
    input.targetFps = 20.0f;
    input.calibratedCostPerSecond = 4000000.0L;
    input.observedStableNanoseconds = 80000000ULL;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget != 300000 || decision.retainedAdmission) {
	std::fprintf(stderr,
	    "FAIL: retained quality floor was not a stable fixed point\n");
	return 1;
    }
    if (policy.noteRetainedQualityFloorMiss() ||
	policy.noteRetainedQualityFloorMiss() ||
	policy.retainedQualityFloorRejected()) {
	std::fprintf(stderr,
	    "FAIL: transient deadline miss rejected quality floor\n");
	return 1;
    }
    policy.setRetainedQualityFloorBudget(cursor,
	300000, 0x2222222222222222ULL, 300000, 50000);
    if (policy.noteRetainedQualityFloorMiss() ||
	policy.retainedQualityFloorRejected()) {
	std::fprintf(stderr,
	    "FAIL: distinct quality-floor populations shared misses\n");
	return 1;
    }
    if (policy.noteRetainedQualityFloorMet(
	    false, 0x2222222222222222ULL, 300000) ||
	policy.retainedQualityFloorMissCount() != 1 ||
	policy.noteRetainedQualityFloorMet(
	    true, 0x1111111111111111ULL, 300000) ||
	policy.retainedQualityFloorMissCount() != 1 ||
	policy.noteRetainedQualityFloorMet(
	    true, 0x2222222222222222ULL, 299999) ||
	policy.retainedQualityFloorMissCount() != 1 ||
	!policy.noteRetainedQualityFloorMet(
	    true, 0x2222222222222222ULL, 300000) ||
	policy.retainedQualityFloorMissCount() != 0) {
	std::fprintf(stderr,
	    "FAIL: unproven presentation reset quality-floor misses\n");
	return 1;
    }
    if (policy.noteRetainedQualityFloorMiss() ||
	policy.noteRetainedQualityFloorMiss() ||
	!policy.noteRetainedQualityFloorMiss()) {
	std::fprintf(stderr,
	    "FAIL: repeated deadline misses did not reject quality floor\n");
	return 1;
    }
    cursor.reset();
    policy.resetOverloadRecovery();
    policy.noteDeadlineCapacityMiss(300000);
    decision = apply_admission_plan(policy, cursor, input);
    if (!policy.retainedQualityFloorRejected() ||
	decision.totalBudget >= 300000 || !decision.retainedAdmission) {
	std::fprintf(stderr,
	    "FAIL: hard-deadline rejection retained the quality floor\n");
	return 1;
    }

    /* A richer static candidate and the protected floor are different
     * populations.  Rejecting the richer candidate must preserve the floor;
     * only a hard miss at the floor's exact cost may retire it. */
    policy.reset();
    cursor.reset();
    policy.raiseCurrentBudget(200000);
    policy.requestRetainedReallocation();
    input = Policy::Inputs();
    input.activeCost = 150000;
    input.minimumActiveCost = 50000;
    decision = apply_admission_plan(policy, cursor, input);
    policy.setRetainedQualityFloorBudget(cursor,
	300000, 0x3333333333333333ULL, 150000, 50000);
	if (policy.capacitySearchMinimumBudget(
		50000, 300000, 0x3333333333333333ULL) != 300000 ||
	    policy.capacitySearchMinimumBudget(
		50000, 300000, 0x9999999999999999ULL) != 50000) {
	    std::fprintf(stderr,
		"FAIL: protected floor did not bound capacity search domain\n");
	    return 1;
	}
	policy.noteDeadlineCapacityMiss(300001, true);
    if (!policy.retainedQualityFloorActive() ||
	policy.retainedQualityFloorRejected()) {
	std::fprintf(stderr,
	    "FAIL: richer static miss retired protected floor\n");
	return 1;
    }
    policy.noteDeadlineCapacityMiss(300000, true);
    if (policy.retainedQualityFloorActive() ||
	!policy.retainedQualityFloorRejected()) {
	std::fprintf(stderr,
	    "FAIL: exact protected-floor miss did not retire floor\n");
	return 1;
    }
    cursor.reset();
    policy.resetOverloadRecovery();
    policy.noteDeadlineCapacityMiss(300000);
    input.activeCost = 300000;
    input.targetFps = 20.0f;
    input.calibratedCostPerSecond = 4000000.0L;
    input.observedStableNanoseconds = 80000000ULL;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget >= 300000 || !decision.retainedAdmission) {
	std::fprintf(stderr,
	    "FAIL: explicit static deadline recovery restored quality floor\n");
	return 1;
    }

    /* A calibrated successor to a minimax allocation must recompute that
     * allocation using the newly demonstrated allowance.  Reusing the old
     * per-occurrence allocation stamps makes the richer budget inert and can
     * leave a warm Hubble view indefinitely coarse. */
    policy.reset();
    cursor.reset();
    policy.raiseCurrentBudget(200000);
    policy.requestRetainedReallocation();
    input = Policy::Inputs();
    input.activeCost = 150000;
    input.minimumActiveCost = 50000;
    decision = apply_admission_plan(policy, cursor, input);
    Policy::CalibrationInputs retainedCalibration;
    retainedCalibration.activeCost = 150000;
    retainedCalibration.observedNanoseconds = 10000000ULL;
    retainedCalibration.passAdmittedWork = true;
    retainedCalibration.allocationCutsApplied = true;
    Policy::CalibrationDecision retainedProbe =
	policy.finishBlockedPass(retainedCalibration);
    Policy::CompletedFrameDecision retainedFrame =
	policy.completeCalibrationFrame(cursor);
    input.targetFps = 20.0f;
    input.calibratedCostPerSecond = 8000000.0L;
    decision = apply_admission_plan(policy, cursor, input);
    if (!retainedProbe.requestFrame ||
	!retainedFrame.restartSubmission ||
	decision.totalBudget != 400000 || !decision.retainedAdmission ||
	decision.retainedAdmissionBudget != 350000) {
	std::fprintf(stderr,
	    "FAIL: calibrated minimax successor did not reallocate its grown "
	    "budget (frame=%d budget=%zu admission=%d upgrade=%zu)\n",
	    retainedFrame.restartSubmission ? 1 : 0, decision.totalBudget,
	    decision.retainedAdmission ? 1 : 0,
	    decision.retainedAdmissionBudget);
	return 1;
    }

    policy.reset();
    cursor.reset();
    input = Policy::Inputs();
    input.interactive = true;
    input.activeCost = 100000;
    input.targetFps = 60.0f;
    input.lastRenderNanoseconds = 4000000ULL;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget < 332000 || decision.totalBudget > 334000 ||
	decision.refinementBudget < 232000 ||
	decision.refinementBudget > 234000) {
	std::fprintf(stderr, "FAIL: interactive retained-cut bootstrap\n");
	return 1;
    }

    policy.reset();
    cursor.reset();
    input = Policy::Inputs();
    input.interactive = true;
    input.scaleQualityProbe = true;
    input.activeCost = 200000;
    input.targetFps = 10.0f;
    input.lastRenderNanoseconds = 40000000ULL;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget != 800000 ||
	decision.refinementBudget != 600000) {
	std::fprintf(stderr, "FAIL: scale-quality discrete probe\n");
	return 1;
    }

    policy.reset();
    cursor.reset();
    input = Policy::Inputs();
    input.activeCost = 100000;
    input.targetFps = 20.0f;
    input.calibratedCostPerSecond = 1000000.0L;
    input.stablePresentationHandoff = true;
    input.stablePresentationCostFloor = 90000;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget != 90000 ||
	decision.retainedAdmission ||
	decision.retainedAdmissionBudget != SIZE_MAX) {
	std::fprintf(stderr, "FAIL: proven presentation cost floor\n");
	return 1;
    }

    cursor.reset();
    input = Policy::Inputs();
    input.forceTerminal = true;
    input.activeCost = 123;
    decision = apply_admission_plan(policy, cursor, input);
    if (decision.totalBudget != SIZE_MAX ||
	decision.refinementBudget != SIZE_MAX ||
	decision.retainedAdmission) {
	std::fprintf(stderr, "FAIL: unbounded terminal policy\n");
	return 1;
    }

    policy.reduceCurrentBudget(90000);
    if (policy.currentBudget() != 90000) {
	std::fprintf(stderr, "FAIL: deadline budget backoff\n");
	return 1;
    }
    policy.raiseCurrentBudget(120000);
    policy.raiseCurrentBudget(100000);
    if (policy.currentBudget() != 120000) {
	std::fprintf(stderr, "FAIL: stable budget restoration\n");
	return 1;
    }
    policy.resetOverloadRecovery();
    if (policy.overloadRecoveryPerformed() ||
	policy.overloadRecoveryActiveCost() != 0) {
	std::fprintf(stderr, "FAIL: overload witness reset\n");
	return 1;
    }

    policy.reset();
    cursor.reset();
    Policy::CalibrationInputs calibrationInput;
    calibrationInput.activeCost = 100000;
    calibrationInput.observedNanoseconds = 10000000ULL;
    Policy::CalibrationDecision calibration =
	finish_blocked_capacity_plan(policy, cursor, calibrationInput);
    if (calibration.candidateReallocation || calibration.searchActive ||
	calibration.requestFrame || calibration.sampleFrame ||
	!calibration.restartSubmission || policy.capacityTransactionPending() ||
	policy.stableBudgetLimited() || cursor.initialized()) {
	std::fprintf(stderr,
	    "FAIL: unchanged ordinary pass did not transfer to allocation\n");
	return 1;
    }
    Policy::Inputs reallocationInputs;
    reallocationInputs.activeCost = 100000;
    reallocationInputs.minimumActiveCost = 25000;
    decision = apply_admission_plan(policy, cursor, reallocationInputs);
    if (!decision.retainedAdmission ||
	decision.retainedAdmissionBudget != 25000) {
	std::fprintf(stderr,
	    "FAIL: ordinary quality debt bypassed retained allocation\n");
	return 1;
    }

    policy.reset();
    cursor.reset();
    (void)apply_admission_plan(policy, cursor, reallocationInputs);
    calibrationInput.passAdmittedWork = true;
	calibration = finish_blocked_capacity_plan(
	    policy, cursor, calibrationInput);
    if (!calibration.requestFrame || calibration.sampleFrame ||
	calibration.restartSubmission || !policy.capacityTransactionPending() ||
	policy.stableBudgetLimited()) {
	std::fprintf(stderr, "FAIL: changed ordinary pass frame barrier\n");
	return 1;
    }
	Policy::CompletedFrameDecision completed =
	    complete_capacity_frame_plan(policy, cursor);
    if (!completed.restartSubmission || policy.capacityTransactionPending()) {
	std::fprintf(stderr, "FAIL: changed ordinary pass allocation transfer\n");
	return 1;
    }
	decision = apply_admission_plan(policy, cursor, reallocationInputs);
    if (!decision.retainedAdmission) {
	std::fprintf(stderr,
	    "FAIL: presented ordinary population did not enter allocation\n");
	return 1;
    }

    policy.reset();
    cursor.reset();
    policy.requestPopulationBarrier();
    if (!policy.capacityTransactionPending() ||
	!policy.presentationFramePending() || policy.stableBudgetLimited()) {
	std::fprintf(stderr, "FAIL: explicit presentation frame barrier\n");
	return 1;
    }
    completed = policy.retireUnmeasurableCalibrationFrame(cursor);
    if (completed.requestSampleFrame ||
	!completed.restartSubmission || policy.capacityTransactionPending()) {
	std::fprintf(stderr, "FAIL: presentation barrier retirement\n");
	return 1;
    }
	decision = apply_admission_plan(policy, cursor, reallocationInputs);
    if (decision.retainedAdmission) {
	std::fprintf(stderr,
	    "FAIL: coverage barrier manufactured quality reallocation\n");
	return 1;
    }
    policy.requestPopulationBarrier();
    policy.invalidateCalibration();
    if (policy.capacityTransactionPending() || policy.stableBudgetLimited() ||
	policy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::INACTIVE) {
	std::fprintf(stderr, "FAIL: full capacity measurement reset\n");
	return 1;
    }

    /* The production capacity seam maps a changed retained allocation to one
     * frozen candidate and exactly three completed-frame samples.  It must
     * choose a distinct successor directly, without the old extra O(scene)
     * scan between unchanged calibration frames. */
    policy.reset();
    cursor.reset();
    policy.reduceCurrentBudget(200);
    BObolLodCapacitySearchKey searchKey;
    searchKey.inventory.set(1);
    searchKey.availability.set(2);
    searchKey.view.set(3);
    searchKey.policy.set(4);
    searchKey.preferredTargetNanoseconds = 50000000ULL;
    searchKey.maximumTargetNanoseconds = 50000000ULL;
    searchKey.maximumCandidateBudget = 800;
    calibrationInput = Policy::CalibrationInputs();
    calibrationInput.activeCost = 150;
    calibrationInput.passAdmittedWork = true;
    calibrationInput.allocationCutsApplied = true;
    calibrationInput.boundedSearch = true;
    calibrationInput.searchKey = searchKey;
	calibrationInput.candidateBudget = 200;
    calibrationInput.maximumCandidateBudget =
	searchKey.maximumCandidateBudget;
    /* A handoff request may already be queued when capacity measurement
     * begins.  A later search candidate must supersede it. */
    policy.requestPresentationReconciliation(180);
    calibration = policy.finishBlockedPass(calibrationInput);
	const uint64_t capacityCandidateProgressUnits =
	    2 + BObolLodCapacitySearchCertificate::sampleLimit() *
		(BObolLodCapacitySearchCertificate::invalidSampleLimit() + 1);
    if (!calibration.requestFrame || !calibration.sampleFrame ||
	policy.capacitySearch().candidateBudget() != 200 ||
	policy.capacitySearch().maximumCandidateCount() != 4 ||
	policy.capacitySearch().progressTotalUnits() !=
	    policy.capacitySearch().maximumCandidateCount() *
		capacityCandidateProgressUnits ||
	policy.capacitySearch().progressCompletedUnits() >=
	    policy.capacitySearch().progressTotalUnits()) {
	std::fprintf(stderr, "FAIL: bounded capacity candidate barrier\n");
	return 1;
    }
    Policy::CompletedFrameInputs frameInput;
    frameInput.searchKey = searchKey;
    frameInput.candidateBudget = 200;
    frameInput.presentedCost = 180;
    frameInput.observedNanoseconds = 25000000ULL;
    frameInput.validSample = true;
    const size_t activeSearchBudget = policy.currentBudget();
    policy.requestPresentationReconciliation(170);
    if (policy.currentBudget() != activeSearchBudget) {
	std::fprintf(stderr,
	    "FAIL: handoff request replaced an active capacity candidate\n");
	return 1;
    }
    for (unsigned int i = 0;
	 i < BObolLodCapacitySearchCertificate::sampleLimit(); ++i) {
	completed = policy.completeCapacitySearchFrame(cursor, frameInput);
	if (i + 1 < BObolLodCapacitySearchCertificate::sampleLimit() &&
	    (!completed.requestSampleFrame ||
	     completed.capacityCandidateChanged ||
	     completed.restartSubmission)) {
	    std::fprintf(stderr,
		"FAIL: bounded capacity sample series index=%u\n", i);
	    return 1;
	}
    }
    if (!completed.restartSubmission ||
	!completed.capacityCandidateChanged ||
	completed.requestSampleFrame ||
	policy.currentBudget() <= 200 ||
	policy.capacitySearch().measuredCandidateCount() != 1 ||
	policy.capacitySearch().progressCompletedUnits() !=
	    capacityCandidateProgressUnits ||
	policy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::ALLOCATING ||
	policy.presentationFramePending() ||
	BObolLodAdmissionPlanner::presentationPausesSubmission(
	    false, false, policy.capacityAllocationPending(),
	    policy.presentationFramePending(), true) ||
	cursor.initialized()) {
	std::fprintf(stderr,
	    "FAIL: bounded capacity successor budget=%zu measured=%u\n",
	    policy.currentBudget(),
	    policy.capacitySearch().measuredCandidateCount());
	return 1;
    }

    const size_t frozenCandidateBudget =
	policy.capacitySearch().candidateBudget();
    completed = policy.completeCapacitySearchFrame(cursor, frameInput);
    if (completed.requestCeilingFreeFrame ||
	completed.capacityCandidateChanged || completed.requestSampleFrame ||
	completed.restartSubmission ||
	policy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::ALLOCATING ||
	policy.capacitySearch().candidateBudget() != frozenCandidateBudget ||
	policy.presentationFramePending()) {
	std::fprintf(stderr,
	    "FAIL: completed frame advanced an unapplied capacity candidate\n");
	return 1;
    }

    /* A completed safety handoff may consume the retained request which was
     * carrying this exact ALLOCATING candidate.  The certificate remains
     * active, but cannot execute itself: recreating its typed producer must
     * restore the candidate budget without manufacturing a presentation
     * frame. */
    policy.reduceCurrentBudget(1);
    if (!resume_capacity_candidate_plan(policy, cursor) ||
	policy.currentBudget() != frozenCandidateBudget ||
	policy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::ALLOCATING ||
	policy.presentationFramePending() || cursor.initialized()) {
	std::fprintf(stderr,
	    "FAIL: safety handoff did not resume capacity allocation\n");
	return 1;
    }

    Policy interruptedAllocationPolicy = policy;
    BObolLodAdmissionCursor interruptedAllocationCursor = cursor;
    const size_t interruptedCandidate =
	interruptedAllocationPolicy.capacitySearch().candidateBudget();
    const unsigned int interruptedMeasured =
	interruptedAllocationPolicy.capacitySearch().measuredCandidateCount();
    Policy::DeadlineMissInputs activeSearchMiss;
    activeSearchMiss.attemptedBudget = interruptedCandidate;
    Policy::CompletedFrameDecision interruptedDecision =
	interruptedAllocationPolicy.noteDeadlineCapacityMiss(
	    activeSearchMiss, &interruptedAllocationCursor);
    if (!interruptedDecision.capacityCandidateChanged ||
	!interruptedDecision.restartSubmission ||
	interruptedDecision.searchResult !=
	    BObolLodCapacitySearchCertificate::Result::REALLOCATE ||
	interruptedAllocationPolicy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::ALLOCATING ||
	interruptedAllocationPolicy.capacitySearch().safeBudget() >=
	    interruptedCandidate ||
	interruptedAllocationPolicy.capacitySearch().measuredCandidateCount() !=
	    interruptedMeasured + 1 ||
	interruptedAllocationPolicy.currentBudget() !=
	    interruptedAllocationPolicy.capacitySearch().candidateBudget() ||
	interruptedAllocationPolicy.currentBudget() <=
	    interruptedAllocationPolicy.capacitySearch().safeBudget() ||
	interruptedAllocationPolicy.currentBudget() >= interruptedCandidate ||
	interruptedAllocationCursor.initialized()) {
	std::fprintf(stderr,
	    "FAIL: interrupted candidate allocation did not narrow search\n");
	return 1;
    }

    /* A complete retained allocation may encounter its first hard deadline
     * before it supplies an exact timing frame.  That miss must create the
     * same bounded certificate as an ordinary completed sample; a scalar
     * five-percent budget decrement would require an unbounded number of
     * whole-scene allocation passes as candidate size grows. */
    Policy initialDeadlinePolicy;
    BObolLodAdmissionCursor initialDeadlineCursor;
    Policy::DeadlineMissInputs initialDeadlineMiss;
    initialDeadlineMiss.attemptedBudget = 400;
    initialDeadlineMiss.searchKey = searchKey;
    initialDeadlineMiss.searchKey.minimumBudget = 100;
    initialDeadlineMiss.candidateBudget = 400;
	Policy::CompletedFrameDecision initialDeadlineDecision =
	    initialDeadlinePolicy.noteDeadlineCapacityMiss(
		initialDeadlineMiss, &initialDeadlineCursor);
    if (!initialDeadlineDecision.capacityCandidateChanged ||
	!initialDeadlineDecision.restartSubmission ||
	initialDeadlineDecision.searchResult !=
	    BObolLodCapacitySearchCertificate::Result::REALLOCATE ||
	initialDeadlinePolicy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::ALLOCATING ||
	initialDeadlinePolicy.capacitySearch().candidateBudget() >= 400 ||
	initialDeadlinePolicy.capacitySearch().candidateBudget() < 100 ||
	initialDeadlinePolicy.capacitySearch().measuredCandidateCount() != 1 ||
	initialDeadlinePolicy.currentBudget() !=
	    initialDeadlinePolicy.capacitySearch().candidateBudget() ||
	initialDeadlineCursor.initialized()) {
	std::fprintf(stderr,
	    "FAIL: first interrupted allocation did not enter bounded search\n");
	return 1;
    }

    /* A completed changed-population frame may leave its generic barrier in
     * place at the same instant a hard deadline selects a distinct capacity
     * candidate.  The completed frame must retire the older barrier and
     * preserve the exact ALLOCATING candidate as the only successor.  If the
     * barrier remains pending, presentationPausesSubmission() prevents that
     * successor from ever applying. */
    Policy barrierCandidatePolicy = policy;
    BObolLodAdmissionCursor barrierCandidateCursor;
    Policy::CalibrationInputs barrierInput;
    barrierInput.passAdmittedWork = true;
    Policy::CalibrationDecision barrierDecision =
	barrierCandidatePolicy.finishBlockedPass(barrierInput);
    if (!barrierDecision.requestFrame ||
	!barrierCandidatePolicy.calibrationFramePending() ||
	!barrierCandidatePolicy.presentationFramePending()) {
	std::fprintf(stderr,
	    "FAIL: capacity candidate barrier setup did not request a frame\n");
	return 1;
    }
    const size_t barrierCandidateBudget =
	barrierCandidatePolicy.capacitySearch().candidateBudget();
    completed = barrierCandidatePolicy.completeCalibrationFrame(
	barrierCandidateCursor);
    if (!completed.restartSubmission ||
	!completed.capacityCandidateChanged ||
	barrierCandidatePolicy.calibrationFramePending() ||
	barrierCandidatePolicy.presentationFramePending() ||
	barrierCandidatePolicy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::ALLOCATING ||
	barrierCandidatePolicy.currentBudget() != barrierCandidateBudget ||
	BObolLodAdmissionPlanner::presentationPausesSubmission(
	    false, false, barrierCandidatePolicy.capacityAllocationPending(),
	    barrierCandidatePolicy.presentationFramePending(), true) ||
	barrierCandidateCursor.initialized()) {
	std::fprintf(stderr,
	    "FAIL: completed barrier did not release capacity candidate "
	    "budget=%zu expected=%zu\n",
	    barrierCandidatePolicy.currentBudget(), barrierCandidateBudget);
	return 1;
    }

	const size_t searchCandidateBudget =
	    policy.capacitySearch().candidateBudget();
	/* Triangle recovery may leave a conservative prefix ceiling while a
	 * bounded capacity search selects its next occurrence population.  The
	 * candidate is the stronger owner and must be applied verbatim; otherwise
	 * the following exact frame describes a stale population and restarts the
	 * same search. */
	policy.requestRetainedRecovery(searchCandidateBudget / 2);
	Policy::Inputs successorInputs;
    successorInputs.activeCost = 180;
    successorInputs.minimumActiveCost = 100;
    successorInputs.targetFps = 20.0f;
    successorInputs.calibratedCostPerSecond = 1000000000.0L;
	decision = apply_admission_plan(policy, cursor, successorInputs);
	if (!decision.retainedAdmission ||
	    decision.totalBudget != searchCandidateBudget ||
	    !policy.retainedRecoveryCeilingActive()) {
	    std::fprintf(stderr,
		"FAIL: recovery ceiling displaced bounded search candidate "
		"candidate=%zu applied=%zu retained=%d ceiling=%d\n",
		searchCandidateBudget, decision.totalBudget,
		decision.retainedAdmission ? 1 : 0,
		policy.retainedRecoveryCeilingActive() ? 1 : 0);
	    return 1;
	}

    /* The successor allocation may be fully applied while an older global
     * ceiling still protects the framebuffer.  That is sufficient to bind
     * the frozen candidate, but not to consume a sample.  The ceiling-free
     * successor frame owns the first measurement. */
    calibrationInput.activeCost = searchCandidateBudget;
    calibrationInput.passAdmittedWork = false;
    calibrationInput.allocationPresentationRealized = false;
    calibrationInput.candidateBudget = searchCandidateBudget;
    calibration = policy.finishBlockedPass(calibrationInput);
    if (!calibration.requestFrame || !calibration.sampleFrame ||
	calibration.restartSubmission ||
	policy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::PRESENTING ||
	policy.capacitySearch().samplesRemaining() !=
	    0 || !policy.presentationFramePending()) {
	std::fprintf(stderr,
	    "FAIL: applied hidden capacity candidate consumed a sample\n");
	return 1;
    }

    /* A presentation-hidden candidate is complete current work, not invalid
     * work.  Preserve its search and request one ceiling-free frame rather
     * than reallocating the same immutable population indefinitely. */
    Policy::CompletedFrameInputs hiddenFrameInput;
    hiddenFrameInput.searchKey = searchKey;
    hiddenFrameInput.candidateBudget = searchCandidateBudget;
    hiddenFrameInput.candidateState = Policy::CompletedFrameInputs::
	CandidateState::PRESENTATION_PENDING;
    Policy presentationGuardPolicy = policy;
    BObolLodAdmissionCursor presentationGuardCursor = cursor;
    completed = presentationGuardPolicy.completeCapacitySearchFrame(
	presentationGuardCursor, hiddenFrameInput);
    if (completed.requestSampleFrame ||
	!completed.requestCeilingFreeFrame ||
	completed.capacityCandidateChanged ||
	completed.restartSubmission ||
	completed.searchResult !=
	    BObolLodCapacitySearchCertificate::Result::PRESENTATION_REQUIRED ||
	!presentationGuardPolicy.capacityTransactionPending() ||
	!presentationGuardPolicy.presentationFramePending() ||
	!BObolLodAdmissionPlanner::presentationPausesSubmission(
	    false, false, presentationGuardPolicy.capacityAllocationPending(),
	    presentationGuardPolicy.presentationFramePending(), true) ||
	presentationGuardPolicy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::PRESENTING ||
	presentationGuardCursor.initialized()) {
	std::fprintf(stderr,
	    "FAIL: capacity presentation guard did not preserve candidate\n");
	return 1;
    }

    /* Removing the ceiling does not itself consume evidence.  The following
     * exact current frame owns the first sample and advances the preserved
     * candidate monotonically into its finite measurement series. */
    Policy::CompletedFrameInputs exactFrameInput = hiddenFrameInput;
    exactFrameInput.candidateState = Policy::CompletedFrameInputs::
	CandidateState::CURRENT;
    exactFrameInput.presentedCost = searchCandidateBudget;
    exactFrameInput.observedNanoseconds = 25000000ULL;
    exactFrameInput.validSample = true;
    completed = presentationGuardPolicy.completeCapacitySearchFrame(
	presentationGuardCursor, exactFrameInput);
    if (!completed.requestSampleFrame ||
	completed.requestCeilingFreeFrame ||
	completed.capacityCandidateChanged || completed.restartSubmission ||
	presentationGuardPolicy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::MEASURING ||
	presentationGuardPolicy.capacitySearch().samplesRemaining() !=
	    BObolLodCapacitySearchCertificate::sampleLimit() - 1 ||
	!presentationGuardPolicy.capacityTransactionPending() ||
	!presentationGuardPolicy.presentationFramePending()) {
	std::fprintf(stderr,
	    "FAIL: ceiling-free capacity frame did not begin measurement\n");
	return 1;
    }

    /* A resident population change invalidates the candidate; it is not a
     * presentation retry.  When the resident-growth producer owns the next
     * allocation, retire the search without racing it or retaining a frame
     * latch.  Without that producer, request one immediate reallocation. */
    Policy deferredInvalidationPolicy = policy;
    BObolLodAdmissionCursor deferredInvalidationCursor = cursor;
    Policy::CompletedFrameInputs invalidatedFrameInput = hiddenFrameInput;
    invalidatedFrameInput.candidateState = Policy::CompletedFrameInputs::
	CandidateState::REALLOCATION_REQUIRED;
    invalidatedFrameInput.reallocationProducerPending = true;
    completed = deferredInvalidationPolicy.completeCapacitySearchFrame(
	deferredInvalidationCursor, invalidatedFrameInput);
    if (completed.requestSampleFrame || completed.restartSubmission ||
	!completed.capacityCandidateChanged ||
	completed.searchResult !=
	    BObolLodCapacitySearchCertificate::Result::STALE_POPULATION ||
	deferredInvalidationPolicy.capacityTransactionPending() ||
	deferredInvalidationCursor.initialized()) {
	std::fprintf(stderr,
	    "FAIL: resident-growth invalidation retained capacity ownership\n");
	return 1;
    }
    Policy immediateInvalidationPolicy = policy;
    BObolLodAdmissionCursor immediateInvalidationCursor = cursor;
    invalidatedFrameInput.reallocationProducerPending = false;
    completed = immediateInvalidationPolicy.completeCapacitySearchFrame(
	immediateInvalidationCursor, invalidatedFrameInput);
    if (completed.requestSampleFrame || !completed.restartSubmission ||
	!completed.capacityCandidateChanged ||
	completed.searchResult !=
	    BObolLodCapacitySearchCertificate::Result::STALE_POPULATION ||
	immediateInvalidationPolicy.capacityTransactionPending() ||
	immediateInvalidationCursor.initialized()) {
	std::fprintf(stderr,
	    "FAIL: unowned invalid capacity population did not reallocate\n");
	return 1;
    }

    /* No ordinary admission policy may change a candidate between its
     * completed-frame samples.  In particular, a stricter 20 Hz throughput
     * estimate must not coarsen a candidate being certified against the
     * static deadline. */
    if (!cursor.initialized()) {
	std::fprintf(stderr,
	    "FAIL: capacity measurement test lost production cursor state\n");
	return 1;
    }
    Policy::Inputs measurementInputs;
    measurementInputs.activeCost = searchCandidateBudget;
    measurementInputs.minimumActiveCost = 100;
    measurementInputs.targetFps = 20.0f;
    measurementInputs.calibratedCostPerSecond = 1000.0L;
    decision = apply_admission_plan(policy, cursor, measurementInputs);
    if (!decision.initialized ||
	decision.totalBudget != searchCandidateBudget ||
	decision.refinementBudget != 0 || decision.retainedAdmission ||
	policy.currentBudget() != searchCandidateBudget) {
	std::fprintf(stderr,
	    "FAIL: capacity measurement did not own population "
	    "candidate=%zu applied=%zu refinement=%zu retained=%d\n",
	    searchCandidateBudget, decision.totalBudget,
	    decision.refinementBudget, decision.retainedAdmission ? 1 : 0);
	return 1;
    }

    /* The production deadline seam must carry an aborted measurement into
     * the same certificate.  Updating the persistent ceiling may not make
     * finishBlockedPass interpret the successor as a new key. */
    static constexpr size_t deadlineCandidateBudget = 400;
    policy.reset();
    cursor.reset();
    policy.reduceCurrentBudget(deadlineCandidateBudget);
    searchKey = BObolLodCapacitySearchKey();
    searchKey.inventory.set(1);
    searchKey.availability.set(2);
    searchKey.view.set(3);
    searchKey.policy.set(4);
    searchKey.preferredTargetNanoseconds = 50000000ULL;
    searchKey.maximumTargetNanoseconds = 50000000ULL;
    searchKey.maximumCandidateBudget = 800;
    calibrationInput = Policy::CalibrationInputs();
    calibrationInput.activeCost = deadlineCandidateBudget;
    calibrationInput.allocationCutsApplied = true;
    calibrationInput.boundedSearch = true;
    calibrationInput.allocationPresentationRealized = true;
    calibrationInput.searchKey = searchKey;
    calibrationInput.candidateBudget = deadlineCandidateBudget;
    calibrationInput.maximumCandidateBudget =
	searchKey.maximumCandidateBudget;
    calibrationInput.observedNanoseconds = 10000000ULL;
    calibration = policy.finishBlockedPass(calibrationInput);
    if (!calibration.requestFrame ||
	policy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::MEASURING) {
	std::fprintf(stderr,
	    "FAIL: deadline regression did not begin capacity measurement\n");
	return 1;
    }
    policy.noteDeadlineCapacityMiss(
	deadlineCandidateBudget, false, &cursor);
    if (policy.capacitySearch().phase() ==
	    BObolLodCapacitySearchCertificate::Phase::INACTIVE ||
	policy.capacitySearch().measuredCandidateCount() != 1 ||
	policy.currentBudget() >= deadlineCandidateBudget ||
	cursor.initialized()) {
	std::fprintf(stderr,
	    "FAIL: deadline seam reset active capacity search "
	    "phase=%u measured=%u budget=%zu\n",
	    static_cast<unsigned int>(policy.capacitySearch().phase()),
	    policy.capacitySearch().measuredCandidateCount(),
	    policy.currentBudget());
	return 1;
    }
    searchKey.preferredBudgetCeiling =
	policy.deadlineCapacityCeiling();
    calibrationInput.searchKey = searchKey;
    calibrationInput.candidateBudget = policy.currentBudget();
    calibrationInput.activeCost = policy.currentBudget();
    calibrationInput.allocationPresentationRealized = false;
    calibration = policy.finishBlockedPass(calibrationInput);
    if (!calibration.requestFrame || calibration.restartSubmission ||
	policy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::PRESENTING ||
	policy.capacitySearch().measuredCandidateCount() != 1) {
	std::fprintf(stderr,
	    "FAIL: deadline successor reopened capacity certificate\n");
	return 1;
    }

    /* Applying the richer allocation selected by an active search is the
     * candidate's presentation boundary, not a capacity invalidation.  The
     * exact successor frame must retain the frozen search and its static
     * deadline.  This failed in the 5k zoom case: the richer population drew
     * once and ordinary steady planning immediately coarsened it again. */
    Policy activeCandidatePolicy;
    BObolLodAdmissionCursor activeCandidateCursor;
    BObolLodCapacitySearchKey activeCandidateKey;
    activeCandidateKey.inventory.set(21);
    activeCandidateKey.availability.set(22);
    activeCandidateKey.view.set(23);
    activeCandidateKey.policy.set(24);
    activeCandidateKey.preferredTargetNanoseconds = 50000000ULL;
    activeCandidateKey.maximumTargetNanoseconds = 400000000ULL;
    activeCandidateKey.startAtStatic = true;
    activeCandidateKey.minimumBudget = 100;
    activeCandidateKey.maximumCandidateBudget = 1000;
    Policy::CalibrationInputs activeCandidateInput;
    activeCandidateInput.activeCost = 100;
    activeCandidateInput.observedNanoseconds = 10000000ULL;
    activeCandidateInput.allocationCutsApplied = true;
    activeCandidateInput.boundedSearch = true;
    activeCandidateInput.allocationPresentationRealized = true;
    activeCandidateInput.searchKey = activeCandidateKey;
    activeCandidateInput.candidateBudget = 100;
    activeCandidateInput.maximumCandidateBudget = 1000;
    Policy::CalibrationDecision activeCandidateDecision =
	activeCandidatePolicy.finishBlockedPass(activeCandidateInput);
    if (!activeCandidateDecision.requestFrame ||
	activeCandidateDecision.searchTerminal) {
	std::fprintf(stderr,
	    "FAIL: active-candidate preservation setup did not start search\n");
	return 1;
    }
    Policy::CompletedFrameInputs activeCandidateFrame;
    activeCandidateFrame.searchKey = activeCandidateKey;
    activeCandidateFrame.candidateBudget = 100;
    activeCandidateFrame.presentedCost = 100;
    activeCandidateFrame.observedNanoseconds = 10000000ULL;
    activeCandidateFrame.validSample = true;
    Policy::CompletedFrameDecision activeCandidateCompletion;
    for (unsigned int i = 1;
	 i < BObolLodCapacitySearchCertificate::sampleLimit(); ++i)
	activeCandidateCompletion = activeCandidatePolicy.
	    completeCapacitySearchFrame(activeCandidateCursor,
		activeCandidateFrame);
    const size_t richerCandidate =
	activeCandidatePolicy.capacitySearch().candidateBudget();
    if (!activeCandidateCompletion.capacityCandidateChanged ||
	activeCandidatePolicy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::ALLOCATING ||
	richerCandidate <= 100 ||
	activeCandidatePolicy.currentBudget() != richerCandidate) {
	std::fprintf(stderr,
	    "FAIL: active-candidate preservation setup did not select richer "
	    "allocation phase=%u candidate=%zu budget=%zu\n",
	    static_cast<unsigned int>(
		activeCandidatePolicy.capacitySearch().phase()),
	    richerCandidate, activeCandidatePolicy.currentBudget());
	return 1;
    }
    const BObolLodAdmissionRevisionStamp activeCandidateStamp(
	activeCandidateKey.inventory, activeCandidateKey.availability,
	BObolLodVisibilityEpoch(1), activeCandidateKey.view,
	activeCandidateKey.policy, BObolLodCapacityEpoch(99));
    Policy::CompletedAllocationInputs activeCandidateAllocation;
    activeCandidateAllocation.revisionStamp = activeCandidateStamp;
    activeCandidateAllocation.requestedSceneBudget = richerCandidate;
    activeCandidateAllocation.certifiedPresentationBudget = richerCandidate;
    activeCandidateAllocation.allocationCertificateCurrent = true;
    activeCandidateAllocation.allocationCutsApplied = true;
    complete_applied_allocation_plan(activeCandidatePolicy,
	activeCandidateCursor, activeCandidateAllocation);
    if (activeCandidatePolicy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::ALLOCATING ||
	activeCandidatePolicy.capacitySearch().candidateBudget() !=
	    richerCandidate ||
	activeCandidatePolicy.stableBudgetLimited()) {
	std::fprintf(stderr,
	    "FAIL: applied active capacity candidate erased frozen search "
	    "phase=%u candidate=%zu\n",
	    static_cast<unsigned int>(
		activeCandidatePolicy.capacitySearch().phase()),
	    activeCandidatePolicy.capacitySearch().candidateBudget());
	return 1;
    }

    /* Applying the allocation selected by a terminal capacity certificate is
     * an ordinary changed-cut pass mechanically, but it is not a new
     * capacity problem.  Retiring that pass must preserve the terminal
     * witness.  Clearing it reproduced a quiet-view loop which repeatedly
     * searched between the same two Lucy wire populations. */
    Policy terminalPolicy;
    BObolLodAdmissionCursor terminalCursor;
    BObolLodCapacitySearchKey terminalKey;
    terminalKey.inventory.set(11);
    terminalKey.availability.set(12);
    terminalKey.view.set(13);
    terminalKey.policy.set(14);
    terminalKey.preferredTargetNanoseconds = 50000000ULL;
    terminalKey.maximumTargetNanoseconds = 50000000ULL;
    terminalKey.minimumBudget = 100;
    terminalKey.maximumCandidateBudget = 100;
    Policy::CalibrationInputs terminalInput;
    terminalInput.activeCost = 100;
    terminalInput.observedNanoseconds = 25000000ULL;
    terminalInput.allocationCutsApplied = true;
    terminalInput.boundedSearch = true;
    terminalInput.allocationPresentationRealized = true;
    terminalInput.searchKey = terminalKey;
    terminalInput.candidateBudget = 100;
    terminalInput.maximumCandidateBudget = 100;
    Policy::CalibrationDecision terminalDecision =
	terminalPolicy.finishBlockedPass(terminalInput);
    if (!terminalDecision.requestFrame || terminalDecision.searchTerminal) {
	std::fprintf(stderr,
	    "FAIL: capacity terminal-retirement setup did not start\n");
	return 1;
    }
    Policy::CompletedFrameInputs terminalFrame;
    terminalFrame.searchKey = terminalKey;
    terminalFrame.candidateBudget = 100;
    terminalFrame.presentedCost = 100;
    terminalFrame.observedNanoseconds = 25000000ULL;
    terminalFrame.validSample = true;
    for (unsigned int i = 1;
	 i < BObolLodCapacitySearchCertificate::sampleLimit(); ++i)
	(void)terminalPolicy.completeCapacitySearchFrame(
	    terminalCursor, terminalFrame);
    const BObolLodAdmissionRevisionStamp terminalStamp(
	terminalKey.inventory, terminalKey.availability,
	BObolLodVisibilityEpoch(1), terminalKey.view, terminalKey.policy,
	BObolLodCapacityEpoch(1));
    if (terminalPolicy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::TERMINAL ||
	!terminalPolicy.stableBudgetLimited() ||
	!terminalPolicy.terminalCertificateCovers(terminalStamp, 100)) {
	std::fprintf(stderr,
	    "FAIL: capacity terminal-retirement setup did not certify\n");
	return 1;
    }
    const BObolLodAdmissionRevisionStamp differentTerminalStamp =
	BObolLodRevisionContract::advance(
	    terminalStamp, BObolLodAdmissionRevisionDomain::VIEW);
    const BObolLodAdmissionRevisionStamp visibilitySuccessorStamp =
	BObolLodRevisionContract::advance(
	    terminalStamp, BObolLodAdmissionRevisionDomain::VISIBILITY);
    if (terminalPolicy.terminalCertificateCovers(
	    differentTerminalStamp, 100) ||
	!terminalPolicy.terminalCertificateCovers(
	    visibilitySuccessorStamp, 100) ||
	terminalPolicy.terminalCertificateCovers(terminalStamp, 99)) {
	std::fprintf(stderr,
	    "FAIL: terminal capacity proof crossed a different renderer problem "
	    "or rejected an exact visibility successor\n");
	return 1;
    }
    Policy::CompletedAllocationInputs terminalAllocation;
    terminalAllocation.revisionStamp = terminalStamp;
    terminalAllocation.requestedSceneBudget = 100;
    terminalAllocation.certifiedPresentationBudget = 100;
    terminalAllocation.allocationCertificateCurrent = true;
    terminalAllocation.allocationCutsApplied = true;
    complete_applied_allocation_plan(
	terminalPolicy, terminalCursor, terminalAllocation);
    terminalDecision = terminalPolicy.finishBlockedPass(terminalInput);
    if (!terminalDecision.searchTerminal || terminalDecision.requestFrame ||
	terminalDecision.restartSubmission ||
	terminalDecision.searchResult !=
	    BObolLodCapacitySearchCertificate::Result::CERTIFIED ||
	terminalPolicy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::TERMINAL ||
	!terminalPolicy.stableBudgetLimited()) {
	std::fprintf(stderr,
	    "FAIL: completed terminal allocation rearmed capacity search\n");
	return 1;
    }

    /* A completed static-quality reconciliation is the terminal authority for
     * the longer event-driven deadline.  It must atomically replace an older
     * steady-capacity certificate; retaining the latter made the next plan
     * fall from the reconciled population back to its cheapest predecessor. */
    Policy staticConstraintPolicy = terminalPolicy;
    BObolLodAdmissionCursor staticConstraintCursor = terminalCursor;
    accept_static_presentation_constraint(staticConstraintPolicy,
	staticConstraintCursor, 150);
    if (staticConstraintPolicy.currentBudget() != 150 ||
	staticConstraintPolicy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::INACTIVE ||
	staticConstraintPolicy.stableBudgetLimited() ||
	staticConstraintCursor.initialized()) {
	std::fprintf(stderr,
	    "FAIL: static presentation constraint retained older capacity authority\n");
	return 1;
    }
    Policy::Inputs staticConstraintInputs;
    staticConstraintInputs.activeCost = 150;
    staticConstraintInputs.minimumActiveCost = 100;
    staticConstraintInputs.targetFps = 2.5f;
    staticConstraintInputs.calibratedCostPerSecond = 100.0L;
    staticConstraintInputs.observedStableNanoseconds = 300000000ULL;
    staticConstraintInputs.hardPresentationDeadlineNanoseconds =
	400000000ULL;
    staticConstraintInputs.hardDeadlinePresentation = true;
    staticConstraintInputs.stablePresentationCostFloor = 150;
    const Policy::Decision staticConstraintDecision = apply_admission_plan(
	staticConstraintPolicy, staticConstraintCursor,
	staticConstraintInputs);
    if (staticConstraintDecision.totalBudget < 150) {
	std::fprintf(stderr,
	    "FAIL: static presentation constraint was undercut budget=%zu\n",
	    staticConstraintDecision.totalBudget);
	return 1;
    }

    /* The terminal certificate and its capacity guard are one proof.  An
     * explicit semantic invalidation may discard both, but must never leave a
     * terminal static deadline paired with unrestricted ordinary admission. */
    Policy invalidatedTerminalPolicy = terminalPolicy;
    BObolLodAdmissionCursor invalidatedTerminalCursor = terminalCursor;
    apply_capacity_evidence_action(invalidatedTerminalPolicy,
	invalidatedTerminalCursor,
	BObolLodAdmissionPlanner::EvidenceAction::CLEAR_CAPACITY_LIMIT);
    if (invalidatedTerminalPolicy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::INACTIVE ||
	invalidatedTerminalPolicy.stableBudgetLimited()) {
	std::fprintf(stderr,
	    "FAIL: terminal capacity invalidation split certificate and guard\n");
	return 1;
    }

    Policy::Inputs postTerminalInputs;
    postTerminalInputs.activeCost = 100;
    postTerminalInputs.minimumActiveCost = 100;
    postTerminalInputs.targetFps = 20.0f;
    postTerminalInputs.calibratedCostPerSecond = 100000.0L;
    postTerminalInputs.observedStableNanoseconds = 25000000ULL;
    const Policy::Decision postTerminalDecision = apply_admission_plan(
	terminalPolicy, terminalCursor, postTerminalInputs);
    if (postTerminalDecision.totalBudget != 100) {
	std::fprintf(stderr,
	    "FAIL: ordinary planning exceeded terminal capacity budget "
	    "budget=%zu\n", postTerminalDecision.totalBudget);
	return 1;
    }
    terminalAllocation.requestedSceneBudget = 99;
    complete_applied_allocation_plan(
	terminalPolicy, terminalCursor, terminalAllocation);
    if (terminalPolicy.capacitySearch().phase() !=
	    BObolLodCapacitySearchCertificate::Phase::INACTIVE ||
	terminalPolicy.stableBudgetLimited()) {
	std::fprintf(stderr,
	    "FAIL: mismatched completed allocation retained terminal evidence\n");
	return 1;
    }

    /* A reversible motion ceiling does not authorize destructive retained
     * coarsening after an orthographic pose-only interaction.  Preserve the
     * currently useful population until a completed miss records a strict
     * ceiling; that ceiling must then remain authoritative. */
    policy.reset();
    cursor.reset();
    policy.requestRetainedReallocation(false);
    Policy::Inputs poseInputs;
    poseInputs.activeCost = 200000;
    poseInputs.minimumActiveCost = 100;
    poseInputs.targetFps = 20.0f;
    poseInputs.calibratedCostPerSecond = 1000000.0L;
    poseInputs.preserveActivePopulation = true;
    decision = apply_admission_plan(policy, cursor, poseInputs);
    if (decision.totalBudget < poseInputs.activeCost ||
	!decision.retainedAdmission) {
	std::fprintf(stderr,
	    "FAIL: pose continuity discarded an unrefuted population "
	    "active=%zu budget=%zu retained=%d\n",
	    poseInputs.activeCost, decision.totalBudget,
	    decision.retainedAdmission ? 1 : 0);
	return 1;
    }

    policy.reset();
    cursor.reset();
    policy.noteDeadlineCapacityMiss(100000, true);
    policy.requestRetainedReallocation(false);
    decision = apply_admission_plan(policy, cursor, poseInputs);
    if (decision.totalBudget >= 100000) {
	std::fprintf(stderr,
	    "FAIL: pose continuity overrode a completed deadline miss "
	    "budget=%zu\n", decision.totalBudget);
	return 1;
    }
    return 0;
}

static int
test_quality_policy(void)
{
    using Policy = BObolLodQualityPolicy;
    using DeadlineSuccessor = Policy::DeadlineSuccessor;
    const auto near = [](float left, float right) {
	return std::fabs(left - right) <= 0.0001f;
    };
    if (Policy::deadlineSuccessor(true, true) !=
	    DeadlineSuccessor::RETRY_TRANSACTION ||
	Policy::deadlineSuccessor(true, false) !=
	    DeadlineSuccessor::RETRY_TRANSACTION ||
	Policy::deadlineSuccessor(false, true) !=
	    DeadlineSuccessor::CONTINUE_POPULATION ||
	Policy::deadlineSuccessor(false, false) !=
	    DeadlineSuccessor::RECOVER_PRESENTATION) {
	std::fprintf(stderr, "FAIL: presentation deadline successor ownership\n");
	return 1;
    }
    if (BObolLodAdmissionPlanner::canonicalSceneCostAtCadCeiling(
	    300, 250, 100) != 150 ||
	BObolLodAdmissionPlanner::canonicalSceneCostAtCadCeiling(
	    200, 250, 100) != 100 ||
	BObolLodAdmissionPlanner::canonicalSceneCostAtCadCeiling(
	    SIZE_MAX, 1, SIZE_MAX) != SIZE_MAX) {
	std::fprintf(stderr,
	    "FAIL: canonical renderer-ceiling allocation currency\n");
	return 1;
    }
    if (BObolLodAdmissionPlanner::capacitySearchPresentedCost(
	    true, 1983470, 1980818) != 1980818 ||
	BObolLodAdmissionPlanner::capacitySearchPresentedCost(
	    false, 1983470, 1980818) != 1983470 ||
	BObolLodAdmissionPlanner::capacitySearchPresentedCost(
	    true, 1983470, 0) != 1983470) {
	std::fprintf(stderr,
	    "FAIL: capacity search did not use its certified population cost\n");
	return 1;
    }
    if (!near(Policy::interactivePixelError(0, 60.0f), 4.0f) ||
	!near(Policy::interactivePixelError(10000000ULL, 60.0f), 1.0f) ||
	!near(Policy::interactivePixelError(10000000ULL, 0.0f), 1.0f) ||
	!near(Policy::interactivePixelError(UINT64_MAX, 60.0f), 64.0f)) {
	std::fprintf(stderr, "FAIL: interactive pixel-error policy\n");
	return 1;
    }
    const float error =
	Policy::interactivePixelError(66666668ULL, 60.0f);
    if (error < 1.99f || error > 2.01f ||
	Policy::progressiveCeiling(21, error) != 20) {
	std::fprintf(stderr, "FAIL: measured pixel-error scaling\n");
	return 1;
    }

    if (!near(Policy::staticFrameTargetFps(
	    20.0f, 100000000ULL), 10.0f) ||
	!near(Policy::staticFrameTargetFps(
	    5.0f, 100000000ULL), 5.0f) ||
	!near(Policy::staticFrameTargetFps(20.0f, 0), 20.0f)) {
	std::fprintf(stderr, "FAIL: static frame deadline cadence\n");
	return 1;
    }
    if (Policy::staticPresentationRenderCostLimit(
	    1000, 50000000ULL, 100000000ULL) != 1900 ||
	Policy::staticPresentationRenderCostLimit(
	    1000, 95000000ULL, 100000000ULL) != 1000 ||
	Policy::staticPresentationRenderCostLimit(
	    1000, 100000001ULL, 100000000ULL) != 0 ||
	Policy::staticPresentationRenderCostLimit(
	    0, 50000000ULL, 100000000ULL) != 0) {
	std::fprintf(stderr, "FAIL: static presentation render-cost limit\n");
	return 1;
    }
    if (Policy::staticLocalMinimumRetryBudget(
	    1000, 1200, 5000, 100000000ULL, 400000000ULL) != 3800 ||
	Policy::staticLocalMinimumRetryBudget(
	    1000, 4200, 5000, 100000000ULL, 400000000ULL) != 4200 ||
	Policy::staticLocalMinimumRetryBudget(
	    1000, 4200, 4300, 100000000ULL, 400000000ULL) != 4200 ||
	Policy::staticLocalMinimumRetryBudget(
	    1000, 0, 5000, 100000000ULL, 400000000ULL) != 0) {
	std::fprintf(stderr, "FAIL: static local-minimum retry budget\n");
	return 1;
    }
    if (Policy::pixelErrorRenderCostFloor(1000, 1.0f, 0.5f) != 4000 ||
	Policy::pixelErrorRenderCostFloor(1000, 0.75f, 0.5f) != 2250 ||
	Policy::pixelErrorRenderCostFloor(1001, 1.0f, 0.75f) != 1780 ||
	Policy::pixelErrorRenderCostFloor(1000, 0.5f, 0.25f) != 4000 ||
	Policy::pixelErrorRenderCostFloor(1000, 1.0f, 1.0f) != 1000 ||
	Policy::pixelErrorRenderCostFloor(0, 1.0f, 0.5f) != 0) {
	std::fprintf(stderr, "FAIL: pixel-error render-cost floor\n");
	return 1;
    }
    if (Policy::incrementalSceneCostBudget(1000, 3000, 9000) != 7000 ||
	Policy::incrementalSceneCostBudget(1000, 3000, 3000) != 1000 ||
	Policy::incrementalSceneCostBudget(1000, 3000, 2000) != 1000 ||
	Policy::incrementalSceneCostBudget(0, 3000, 9000) != 0 ||
	Policy::incrementalSceneCostBudget(
	    SIZE_MAX - 5, 1, SIZE_MAX) != SIZE_MAX) {
	std::fprintf(stderr, "FAIL: incremental scene-cost budget\n");
	return 1;
    }

    for (unsigned int mask = 0; mask < 8; ++mask) {
	const bool stableContext = (mask & 1u) != 0;
	const bool exactFrame = (mask & 2u) != 0;
	const uint64_t reusableNanoseconds = (mask & 4u) ? 1 : 0;
	const bool expected = stableContext && exactFrame &&
	    reusableNanoseconds == 0;
	if (Policy::staticQualityTimingReplayRequired(stableContext,
		exactFrame, reusableNanoseconds) != expected) {
	    std::fprintf(stderr,
		"FAIL: static-quality timing replay ownership\n");
	    return 1;
	}
    }

    /* Exact completed frames turn a coarse interaction carry-over into the
     * finest measured intermediate target, then admit the raster-stable
     * subpixel tiers whose inverse-square work estimate fits both time and
     * scene-cost headroom.  Below one pixel, the first pass may move directly
     * to 0.5, take the intermediate 0.75 rung, or remain at the fast one-
     * pixel target.  A separate completed 0.5-pixel witness is required
     * before the raster-stable 0.25 tier.  The witness deliberately includes
     * a completed first-use frame: upload/preparation cost only makes it
     * conservative. */
    if (!near(Policy::stablePixelError(
	    1.0f, 100, 400, 10000000ULL, 20.0f, true, true), 0.5f) ||
	!near(Policy::stablePixelError(
	    1.0f, 100, SIZE_MAX, 10000000ULL, 20.0f, true, true), 0.5f) ||
	!near(Policy::stablePixelError(
	    0.75f, 100, 225, 10000000ULL, 20.0f, true, true), 0.5f) ||
	!near(Policy::stablePixelError(
	    0.5f, 100, 400, 10000000ULL, 20.0f, true, true), 0.25f) ||
	!near(Policy::stablePixelError(
	    0.5f, 101, 400, 10000000ULL, 20.0f, true, true), 0.5f) ||
	!near(Policy::stablePixelError(
	    0.25f, 100, 10000, 1000000ULL, 20.0f, true, true), 0.25f) ||
	!near(Policy::stablePixelError(
	    1.0f, 101, 400, 10000000ULL, 20.0f, true, true), 0.75f) ||
	!near(Policy::stablePixelError(
	    1.0f, 100, 200, 20000000ULL, 20.0f, true, true), 0.75f) ||
	!near(Policy::stablePixelError(
	    1.0f, 113, 200, 29000000ULL, 20.0f, true, true), 1.0f) ||
	!near(Policy::stablePixelError(
	    1.0f, 110, 200, 29000000ULL,
	    Policy::staticFrameTargetFps(20.0f, 100000000ULL),
	    true, true), 0.75f) ||
	!near(Policy::stablePixelError(
	    1.0f, 100, 400, 10000000ULL, 20.0f, false, true), 1.0f) ||
	!near(Policy::stablePixelError(
	    1.0f, 100, 400, 10000000ULL, 20.0f, true, false), 1.0f) ||
	!near(Policy::stablePixelError(
	    4.0f, 100, 400, 10000000ULL, 20.0f, true, true), 2.0f) ||
	!near(Policy::stablePixelError(
	    4.0f, 100, 100, 20000000ULL, 20.0f, true, true), 4.0f)) {
	std::fprintf(stderr, "FAIL: stable subpixel headroom policy\n");
	return 1;
    }


    if (!near(Policy::pointProxyThreshold(
	    0.0f, 10000000ULL, 60.0f), 1.0f) ||
	!near(Policy::pointProxyThreshold(
	    70.0f, 10000000ULL, 60.0f), 64.0f) ||
	!near(Policy::pointProxyThreshold(
	    std::numeric_limits<float>::quiet_NaN(),
	    10000000ULL, 60.0f), 1.0f)) {
	std::fprintf(stderr, "FAIL: point-proxy threshold bounds\n");
	return 1;
    }

    if (Policy::deadlineRecoveryCostLimit(
	    1000, 400000000ULL, 400000000ULL) != 950 ||
	Policy::deadlineRecoveryCostLimit(
	    1000, 800000000ULL, 400000000ULL) != 475 ||
	Policy::deadlineRecoveryCostLimit(
	    0, 800000000ULL, 400000000ULL) != 0) {
	std::fprintf(stderr, "FAIL: deadline cost correction\n");
	return 1;
    }

    if (Policy::quietPresentationDeadline(
	    400000000ULL, 200000000ULL, true, true) != 400000000ULL ||
	Policy::quietPresentationDeadline(
	    400000000ULL, 200000000ULL, true, false) != 200000000ULL ||
	Policy::quietPresentationDeadline(
	    400000000ULL, 200000000ULL, false, false) != 400000000ULL) {
	std::fprintf(stderr, "FAIL: quiet deadline owner precedence\n");
	return 1;
    }

    {
	using PreparationState = Obol::CadPresentationPreparationState;
	Obol::CadPresentationPreparationSnapshot empty;
	Obol::CadPresentationPreparationSnapshot started;
	started.target.kind =
	    Obol::CadPresentationPreparationKind::RetainedIndirect;
	started.target.obligationRevision = 7;
	started.target.planRevision = 11;
	started.state = PreparationState::Preparing;
	started.totalUnits = 20;
	if (BObolCadPreparationPolicy::classify(empty, started) !=
		BOBOL_CAD_PREPARATION_STARTED) {
	    std::fprintf(stderr, "FAIL: preparation start classification\n");
	    return 1;
	}
	Obol::CadPresentationPreparationSnapshot advanced = started;
	advanced.completedUnits = 5;
	if (BObolCadPreparationPolicy::classify(started, advanced) !=
		BOBOL_CAD_PREPARATION_ADVANCED ||
	    BObolCadPreparationPolicy::classify(advanced, advanced) !=
		BOBOL_CAD_PREPARATION_NONE) {
	    std::fprintf(stderr, "FAIL: strict preparation progress witness\n");
	    return 1;
	}
	Obol::CadPresentationPreparationSnapshot completed = advanced;
	completed.state = PreparationState::Complete;
	completed.completedUnits = completed.totalUnits;
	if (BObolCadPreparationPolicy::classify(advanced, completed) !=
		BOBOL_CAD_PREPARATION_COMPLETED ||
	    BObolCadPreparationPolicy::classify(completed, completed) !=
		BOBOL_CAD_PREPARATION_NONE) {
	    std::fprintf(stderr, "FAIL: preparation completion edge\n");
	    return 1;
	}
	Obol::CadPresentationPreparationSnapshot regressed = advanced;
	regressed.completedUnits = 4;
	if (BObolCadPreparationPolicy::classify(advanced, regressed) !=
		BOBOL_CAD_PREPARATION_FAILED) {
	    std::fprintf(stderr, "FAIL: preparation regression rejection\n");
	    return 1;
	}
	Obol::CadPresentationPreparationSnapshot malformedStart = started;
	malformedStart.completedUnits = malformedStart.totalUnits + 1;
	Obol::CadPresentationPreparationSnapshot malformedConstraint = advanced;
	malformedConstraint.state = PreparationState::Constrained;
	malformedConstraint.totalUnits++;
	Obol::CadPresentationPreparationSnapshot malformedCompletion = advanced;
	malformedCompletion.state = PreparationState::Complete;
	if (BObolCadPreparationPolicy::classify(empty, malformedStart) !=
		BOBOL_CAD_PREPARATION_FAILED ||
	    BObolCadPreparationPolicy::classify(
		advanced, malformedConstraint) !=
		BOBOL_CAD_PREPARATION_FAILED ||
	    BObolCadPreparationPolicy::classify(
		advanced, malformedCompletion) !=
		BOBOL_CAD_PREPARATION_FAILED) {
	    std::fprintf(stderr,
		"FAIL: preparation certificate validation precedence\n");
	    return 1;
	}
    }

    if (!Policy::retryTransientPresentation(
	    false, 1, false, true, false, false, false) ||
	!Policy::retryTransientPresentation(
	    false, 1, false, false, true, false, false) ||
	!Policy::retryTransientPresentation(
	    false, 1, false, false, false, true, false) ||
	!Policy::retryTransientPresentation(
	    false, 1, false, false, false, false, true) ||
	Policy::retryTransientPresentation(
	    true, 1, false, true, true, true, true) ||
	Policy::retryTransientPresentation(
	    false, 2, false, true, true, true, true) ||
	!Policy::retryTransientPresentation(
	    true, 4, true, false, false, false, false)) {
	std::fprintf(stderr,
	    "FAIL: certified preparation deadline retry classification\n");
	return 1;
    }
    const float proxy =
	Policy::pointProxyThreshold(4.0f, 66666668ULL, 60.0f);
    if (proxy < 7.99f || proxy > 8.01f) {
	std::fprintf(stderr, "FAIL: cumulative point-proxy correction\n");
	return 1;
    }

    BObolLodPointProxyEvidence pointCalibration;
    if (!BObolLodAdmissionPlanner::adaptivePointAggregationAllowed(
	    false, false) ||
	!BObolLodAdmissionPlanner::adaptivePointAggregationAllowed(
	    false, true) ||
	!BObolLodAdmissionPlanner::adaptivePointAggregationAllowed(
	    true, true) ||
	BObolLodAdmissionPlanner::adaptivePointAggregationAllowed(
	    true, false)) {
	std::fprintf(stderr,
	    "FAIL: compact-scene static quality did not gate point aggregation\n");
	return 1;
    }
    const auto directStaticQuality = apply_point_completed_plan(
	pointCalibration, 1.0f, 100000000ULL, 20.0f, true, 0, false);
    if (directStaticQuality.changed ||
	std::fabs(directStaticQuality.threshold - 1.0f) > 0.001f) {
	std::fprintf(stderr,
	    "FAIL: preferred cadence coarsened a compact scene before its hard trial\n");
	return 1;
    }
    const auto rejectedStaticQuality = apply_point_completed_plan(
	pointCalibration, 1.0f, 100000000ULL, 20.0f, true, 0, true);
    if (!rejectedStaticQuality.changed ||
	rejectedStaticQuality.threshold <= 1.01f) {
	std::fprintf(stderr,
	    "FAIL: hard static rejection did not release compact-scene aggregation\n");
	return 1;
    }
    pointCalibration.reset();
    if (!BObolLodAdmissionPlanner::pointDeadlineRequiresPopulationAggregation(1000, 1000, 4, 3, 900) ||
	!BObolLodAdmissionPlanner::pointDeadlineRequiresPopulationAggregation(2000, 1000, 0, 0, 900) ||
	!BObolLodAdmissionPlanner::pointDeadlineRequiresPopulationAggregation(2000, 1000, -1, -1, 900) ||
	!BObolLodAdmissionPlanner::pointDeadlineRequiresPopulationAggregation(2000, 1000, 4, 0, 900) ||
	BObolLodAdmissionPlanner::pointDeadlineRequiresPopulationAggregation(2000, 1000, 4, 1, 900) ||
	BObolLodAdmissionPlanner::pointDeadlineRequiresPopulationAggregation(2000, 1000, 4, 0, 1100)) {
	std::fprintf(stderr,
	    "FAIL: point population deadline floor classification\n");
	return 1;
    }
    if (!Policy::quietPresentationIrreducible(
	    false, true, -1, 2000, 1000, 0, false) ||
	!Policy::quietPresentationIrreducible(
	    false, true, 0, 2000, 1000, 2, false) ||
	!Policy::quietPresentationIrreducible(
	    false, true, -1, 1000, 1000, 2, false) ||
	Policy::quietPresentationIrreducible(
	    true, true, -1, 2000, 1000, 0, false) ||
	Policy::quietPresentationIrreducible(
	    false, false, -1, 2000, 1000, 0, false) ||
	Policy::quietPresentationIrreducible(
	    false, true, -1, 2000, 1000, 0, true) ||
	Policy::quietPresentationIrreducible(
	    false, true, -1, 2000, 1000, 2, false)) {
	std::fprintf(stderr,
	    "FAIL: quiet irreducible presentation classification\n");
	return 1;
    }
    if (BObolLodAdmissionPlanner::pointAggregationApplicable(0) ||
	BObolLodAdmissionPlanner::pointAggregationApplicable(1) ||
	!BObolLodAdmissionPlanner::pointAggregationApplicable(2)) {
	std::fprintf(stderr,
	    "FAIL: point aggregation admitted a single prominent mesh\n");
	return 1;
    }
    /* Both an unbracketed camera invalidation and bracketed interaction entry
     * must preserve an already presented aggregate-point population while the
     * replacement visibility census is transiently empty. */
    if (BObolLodAdmissionPlanner::pointAggregationApplicableAcrossCameraInvalidation(false, 1.0f) ||
	!BObolLodAdmissionPlanner::pointAggregationApplicableAcrossCameraInvalidation(true, 1.0f) ||
	!BObolLodAdmissionPlanner::pointAggregationApplicableAcrossCameraInvalidation(false, 64.0f) ||
	!BObolLodAdmissionPlanner::pointStreamingPopulationWorkPending(false, false, false, 1) ||
	!BObolLodAdmissionPlanner::pointStreamingPopulationWorkPending(true, false, false, 0) ||
	BObolLodAdmissionPlanner::pointStreamingPopulationWorkPending(false, false, false, 0)) {
	std::fprintf(stderr,
	    "FAIL: point aggregation discarded camera/streaming capacity proof\n");
	return 1;
    }
    /* Subpixel points are the intended terminal representation at the
     * ordinary one-pixel cut.  Requiring this count to be zero retains the
     * recovery budget ceiling forever on Hubble-like many-part scenes. */
    if (!BObolLodAdmissionPlanner::onePixelPresentationReady(
	    true, true, true, 700, 0, -1, 12, 1.0f, 1000) ||
	!BObolLodAdmissionPlanner::onePixelPresentationReady(
	    true, true, true, 700, 0, 12, 12, 1.0f, 1000) ||
	BObolLodAdmissionPlanner::onePixelPresentationReady(
	    true, true, true, 700, 1, -1, 12, 1.0f, 1000) ||
	BObolLodAdmissionPlanner::onePixelPresentationReady(
	    true, true, true, 700, 0, 4, 5, 1.0f, 1000) ||
	BObolLodAdmissionPlanner::onePixelPresentationReady(
	    true, true, true, 700, 0, -1, 12, 1.02f, 1000) ||
	BObolLodAdmissionPlanner::onePixelPresentationReady(
	    true, false, true, 700, 0, -1, 12, 1.0f, 1000)) {
	std::fprintf(stderr,
	    "FAIL: terminal one-pixel point population classification\n");
	return 1;
    }
    if (!BObolLodAdmissionPlanner::capacitySamplePopulationReady(
	    true, true, true, 0, true) ||
	BObolLodAdmissionPlanner::capacitySamplePopulationReady(
	    true, true, true, 1, true) ||
	BObolLodAdmissionPlanner::capacitySamplePopulationReady(
	    true, true, true, 0, false)) {
	std::fprintf(stderr,
	    "FAIL: exact capacity-sample population classification\n");
	return 1;
    }
    if (!BObolLodAdmissionPlanner::
	    capacitySampleYieldsToStructuralFrontier(true, true, true, 1) ||
	BObolLodAdmissionPlanner::
	    capacitySampleYieldsToStructuralFrontier(false, true, true, 1) ||
	BObolLodAdmissionPlanner::
	    capacitySampleYieldsToStructuralFrontier(true, false, true, 1) ||
	BObolLodAdmissionPlanner::
	    capacitySampleYieldsToStructuralFrontier(true, true, false, 1) ||
	BObolLodAdmissionPlanner::
	    capacitySampleYieldsToStructuralFrontier(true, true, true, 0)) {
	std::fprintf(stderr,
	    "FAIL: structural frontier did not release impossible capacity sample\n");
	return 1;
    }
    if (!BObolLodAdmissionPlanner::pointBlocksSourceAdmission(
	    true, false) ||
	BObolLodAdmissionPlanner::pointBlocksSourceAdmission(
	    false, true) ||
	!BObolLodAdmissionPlanner::pointBlocksSourceAdmission(
	    true, true)) {
	std::fprintf(stderr,
	    "FAIL: stable point calibration blocked its source prerequisite\n");
	return 1;
    }
    if (!BObolLodAdmissionPlanner::structuralFallbackAdmissionDeferred(
	    true, true, true, 12, 0, false, false) ||
	BObolLodAdmissionPlanner::structuralFallbackAdmissionDeferred(
	    false, true, true, 12, 0, false, false) ||
	BObolLodAdmissionPlanner::structuralFallbackAdmissionDeferred(
	    true, false, true, 12, 0, false, false) ||
	BObolLodAdmissionPlanner::structuralFallbackAdmissionDeferred(
	    true, true, false, 12, 0, false, false) ||
	BObolLodAdmissionPlanner::structuralFallbackAdmissionDeferred(
	    true, true, true, 12, 12, false, false) ||
	BObolLodAdmissionPlanner::structuralFallbackAdmissionDeferred(
	    true, true, true, 12, 0, true, false) ||
	BObolLodAdmissionPlanner::structuralFallbackAdmissionDeferred(
	    true, true, true, 12, 0, false, true)) {
	std::fprintf(stderr,
	    "FAIL: structural fallback admission escaped its exclusive owner\n");
	return 1;
    }
    if (!BObolLodAdmissionPlanner::structuralCoverageOnly(
	    false, true, false, false, false) ||
	!BObolLodAdmissionPlanner::structuralCoverageOnly(
	    false, false, false, false, true) ||
	BObolLodAdmissionPlanner::structuralCoverageOnly(
	    true, true, false, false, false) ||
	BObolLodAdmissionPlanner::structuralCoverageOnly(
	    true, false, false, false, true) ||
	BObolLodAdmissionPlanner::structuralCoverageOnly(
	    false, true, true, false, false) ||
	BObolLodAdmissionPlanner::structuralCoverageOnly(
	    false, true, false, true, false)) {
	std::fprintf(stderr,
	    "FAIL: structural-only coverage escaped its controller domain\n");
	return 1;
    }
    if (!BObolLodAdmissionPlanner::presentationPausesSubmission(
	    true, false, false, false, false) ||
	!BObolLodAdmissionPlanner::presentationPausesSubmission(
	    false, false, false, true, false) ||
	!BObolLodAdmissionPlanner::presentationPausesSubmission(
	    false, true, false, false, true) ||
	BObolLodAdmissionPlanner::presentationPausesSubmission(
	    false, true, false, false, false) ||
	BObolLodAdmissionPlanner::presentationPausesSubmission(
	    false, false, false, false, true) ||
	BObolLodAdmissionPlanner::presentationPausesSubmission(
	    false, true, true, false, true)) {
	std::fprintf(stderr,
	    "FAIL: presentation measurement did not own submission pause\n");
	return 1;
    }

    /* Exercise the composition used by the completed-frame callback, not
     * only its two predicates in isolation.  A stable CAD calibration owns
     * the presentation and pauses the otherwise active submission cursor;
     * that same cursor may therefore never stand in for the successor frame.
     */
    BObolLodAdmissionPlanner::PointCalibrationProducerInputs producerInputs;
    producerInputs.submissionPending = true;
    producerInputs.stableCalibrationPending = true;
    producerInputs.stablePresentationAvailable = true;
    if (!BObolLodAdmissionPlanner::presentationPausesSubmission(
	    false, true, false, false, true) ||
	BObolLodAdmissionPlanner::pointProducerOwnsCalibrationFrame(
	    producerInputs)) {
	std::fprintf(stderr,
	    "FAIL: stable point calibration formed a submission closed wait\n");
	return 1;
    }
    /* A newly armed one-pixel frame must not be replaced by another coarse
     * structural seed during the completion callback which armed it. */
    if (!BObolLodAdmissionPlanner::maySeedPointStructuralDistribution(false, false, false) ||
	BObolLodAdmissionPlanner::maySeedPointStructuralDistribution(true, false, false) ||
	BObolLodAdmissionPlanner::maySeedPointStructuralDistribution(false, true, false) ||
	BObolLodAdmissionPlanner::maySeedPointStructuralDistribution(false, false, true)) {
	std::fprintf(stderr,
	    "FAIL: point presentation barrier permits structural reseeding\n");
	return 1;
    }
    /* Admission is deliberately paused until the changed point population
     * has rendered.  That paused cursor cannot also be the producer which is
     * expected to publish the frame, or both transitions wait forever. */

    producerInputs =
	BObolLodAdmissionPlanner::PointCalibrationProducerInputs();
    producerInputs.submissionPending = true;
    producerInputs.stableCalibrationPending = true;
    producerInputs.stablePresentationAvailable = true;
    const bool pausedSubmissionOwnsFrame =
	BObolLodAdmissionPlanner::pointProducerOwnsCalibrationFrame(
	    producerInputs);
    producerInputs.capacityAllocationPending = true;
    const bool capacityAllocationSubmissionOwnsFrame =
	BObolLodAdmissionPlanner::pointProducerOwnsCalibrationFrame(
	    producerInputs);
    producerInputs.capacityAllocationPending = false;
    producerInputs.stableCalibrationPending = false;
    const bool activeSubmissionOwnsFrame =
	BObolLodAdmissionPlanner::pointProducerOwnsCalibrationFrame(
	    producerInputs);
    producerInputs =
	BObolLodAdmissionPlanner::PointCalibrationProducerInputs();
    producerInputs.stableCalibrationPending = true;
    producerInputs.stablePresentationAvailable = true;
    producerInputs.providerPending = true;
    const bool providerOwnsFrame =
	BObolLodAdmissionPlanner::pointProducerOwnsCalibrationFrame(
	    producerInputs);
    producerInputs.providerPending = false;
    producerInputs.servicePending = true;
    const bool serviceOwnsFrame =
	BObolLodAdmissionPlanner::pointProducerOwnsCalibrationFrame(
	    producerInputs);
    producerInputs.servicePending = false;
    producerInputs.publicationAwaitingFrameRequest = true;
    const bool publicationOwnsFrame =
	BObolLodAdmissionPlanner::pointProducerOwnsCalibrationFrame(
	    producerInputs);
    if (pausedSubmissionOwnsFrame ||
	!capacityAllocationSubmissionOwnsFrame ||
	!activeSubmissionOwnsFrame ||
	!providerOwnsFrame || !serviceOwnsFrame || !publicationOwnsFrame) {
	std::fprintf(stderr,
	    "FAIL: point-calibration producer ownership permits a closed wait\n");
	return 1;
    }
    for (unsigned int combination = 0; combination < 512; ++combination) {
	producerInputs =
	    BObolLodAdmissionPlanner::PointCalibrationProducerInputs();
	producerInputs.submissionPending = (combination & (1u << 0)) != 0;
	producerInputs.discoveryCalibrationPending =
	    (combination & (1u << 1)) != 0;
	producerInputs.stableCalibrationPending =
	    (combination & (1u << 2)) != 0;
	producerInputs.capacityAllocationPending =
	    (combination & (1u << 3)) != 0;
	producerInputs.capacitySamplePending =
	    (combination & (1u << 4)) != 0;
	producerInputs.stablePresentationAvailable =
	    (combination & (1u << 5)) != 0;
	producerInputs.providerPending = (combination & (1u << 6)) != 0;
	producerInputs.servicePending = (combination & (1u << 7)) != 0;
	producerInputs.publicationAwaitingFrameRequest =
	    (combination & (1u << 8)) != 0;
	const bool submissionPaused =
	    producerInputs.capacitySamplePending ||
	    producerInputs.discoveryCalibrationPending ||
	    (producerInputs.stableCalibrationPending &&
	     producerInputs.stablePresentationAvailable &&
	     !producerInputs.capacityAllocationPending);
	const bool expectedProducer =
	    producerInputs.publicationAwaitingFrameRequest ||
	    producerInputs.providerPending || producerInputs.servicePending ||
	    (producerInputs.submissionPending && !submissionPaused);
	if (BObolLodAdmissionPlanner::pointProducerOwnsCalibrationFrame(
		producerInputs) != expectedProducer) {
	    std::fprintf(stderr,
		"FAIL: point-calibration producer combination %u\n",
		combination);
	    return 1;
	}
    }
    if (!BObolLodAdmissionPlanner::pointCalibrationOwnsCompletedFrame(
	    true, false) ||
	BObolLodAdmissionPlanner::pointCalibrationOwnsCompletedFrame(
	    true, true) ||
	BObolLodAdmissionPlanner::pointCalibrationOwnsCompletedFrame(
	    false, false)) {
	std::fprintf(stderr,
	    "FAIL: point calibration consumed a capacity-owned frame\n");
	return 1;
    }
    float threshold = 4.0f;
    bool settled = false;
    for (size_t iteration = 0; iteration < 16; ++iteration) {
	BObolLodPointProxyEvidence::Decision decision;
	if (threshold < 6.0f) {
	    decision = apply_point_interrupted_plan(pointCalibration,
		threshold, 66666668ULL, 60.0f);
	} else {
	    decision = apply_point_completed_plan(pointCalibration,
		threshold, 10000000ULL, 60.0f, true);
	}
	threshold = decision.threshold;
	if (!decision.changed) {
	    settled = true;
	    break;
	}
    }
    if (!settled || threshold < 6.0f || threshold > 6.5f) {
	std::fprintf(stderr,
	    "FAIL: point-proxy bracket did not converge (threshold=%g)\n",
	    threshold);
	return 1;
    }

    /* The responsive and event-driven static contracts intentionally use
     * different deadlines.  A frame may miss 20 FPS while remaining a safe
     * static endpoint, so the responsive unsafe bracket must not force the
     * static search toward a coarser cut. */
    BObolLodAdmissionEvidence pointGoalEvidence;
    const BObolLodAdmissionPlan responsiveMiss =
	BObolLodAdmissionPlanner::planPointInterrupted(
	    pointGoalEvidence, BObolLodAdmissionCursor(),
	    BObolLodPointCalibrationGoal::RESPONSIVE,
	    32.0f, 112000000ULL, 20.0f);
    pointGoalEvidence = responsiveMiss.nextEvidence;
    const BObolLodAdmissionPlan staticHeadroom =
	BObolLodAdmissionPlanner::planPointCompleted(
	    pointGoalEvidence, BObolLodAdmissionCursor(),
	    BObolLodPointCalibrationGoal::STATIC,
	    32.0f, 112000000ULL, 2.5f, true, true, 0);
    if (!responsiveMiss.pointProxyDecision.changed ||
	responsiveMiss.pointProxyDecision.threshold <= 32.0f ||
	!staticHeadroom.pointProxyDecision.changed ||
	staticHeadroom.pointProxyDecision.threshold >= 32.0f) {
	std::fprintf(stderr,
	    "FAIL: responsive point evidence contaminated static quality\n");
	return 1;
    }

    /* The terminal reducer runs after the transient point-phase latch has
     * retired.  Its goal must still select the event-driven deadline rather
     * than falling back to the ordinary quiet cadence. */
    if (std::fabs(BObolLodAdmissionPlanner::pointCalibrationTargetFps(
	    BObolLodPointCalibrationGoal::RESPONSIVE, 10.0f, 2.5f) -
	    10.0f) > 0.001f ||
	std::fabs(BObolLodAdmissionPlanner::pointCalibrationTargetFps(
	    BObolLodPointCalibrationGoal::STATIC, 10.0f, 2.5f) -
	    2.5f) > 0.001f ||
	BObolLodAdmissionPlanner::pointCalibrationPresentationDeadline(
	    BObolLodPointCalibrationGoal::RESPONSIVE,
	    100000000ULL, 400000000ULL) != 100000000ULL ||
	BObolLodAdmissionPlanner::pointCalibrationPresentationDeadline(
	    BObolLodPointCalibrationGoal::STATIC,
	    100000000ULL, 400000000ULL) != 400000000ULL) {
	std::fprintf(stderr,
	    "FAIL: point evidence goal lost its deadline binding\n");
	return 1;
    }

    /* Once a coarse safe side exists, an overloaded frame's measured timing
     * correction may skip redundant geometric probes but may never cross
     * that bound.  This is the 150k-software-renderer case: several nearby
     * thresholds can classify exactly the same occurrence population. */
    pointCalibration.reset();
    constexpr float coarseSafeThreshold = 64.0f;
    (void)apply_point_completed_plan(pointCalibration,
	coarseSafeThreshold, 350000000ULL, 2.5f, true);
    const auto measuredCorrection = apply_point_interrupted_plan(
	pointCalibration, 32.0f, 308000000ULL, 10.0f);
    if (!measuredCorrection.changed ||
	measuredCorrection.threshold < 56.0f ||
	measuredCorrection.threshold > 56.5f) {
	std::fprintf(stderr,
	    "FAIL: point overload ignored measured bracket correction "
	    "(threshold=%g)\n", measuredCorrection.threshold);
	return 1;
    }
    const auto boundedCorrection = apply_point_interrupted_plan(
	pointCalibration, measuredCorrection.threshold,
	270000000ULL, 10.0f);
    if (!boundedCorrection.changed ||
	!near(boundedCorrection.threshold, coarseSafeThreshold)) {
	std::fprintf(stderr,
	    "FAIL: point overload did not reach proven coarse bound "
	    "(threshold=%g)\n", boundedCorrection.threshold);
	return 1;
    }

    /* A hard endpoint abort is not a completed-frame timing sample.  Even if
     * it falls inside the latter's five-percent jitter deadband, retrying the
     * unchanged irreducible population would be a no-progress render loop. */
    pointCalibration.reset();
    const auto nearDeadlineAbort = apply_point_interrupted_plan(pointCalibration,
	1.0f, 101000000ULL, 10.0f);
    if (!nearDeadlineAbort.changed ||
	!near(nearDeadlineAbort.threshold, 1.25f)) {
	std::fprintf(stderr,
	    "FAIL: hard point deadline retained an unchanged population\n");
	return 1;
    }
    pointCalibration.reset();
    const auto preparationOnly = apply_point_completed_plan(pointCalibration,
	8.0f, 10000000ULL, 60.0f, false);
    if (preparationOnly.changed || preparationOnly.continueRelaxation ||
	!near(preparationOnly.threshold, 8.0f)) {
	std::fprintf(stderr,
	    "FAIL: non-reusable point-proxy preparation changed quality\n");
	return 1;
    }

    /* A structural repair which exhausts the scene allowance must aggregate
     * its low-importance tail rather than spending one fast frame per mesh
     * wave.  The population-sized step is deterministic and shares the
     * timing controller's safe/unsafe bracket. */
    pointCalibration.reset();
    const std::array<size_t, 7> projectedStructural = {
	100, 800, 3000, 7000, 9000, 9900, 10000
    };
    const auto structuralSeed = apply_point_structural_distribution_plan(pointCalibration,
	1.0f, projectedStructural, 10000, 1500);
    if (!structuralSeed.changed ||
	!near(structuralSeed.threshold, 16.0f)) {
	std::fprintf(stderr,
	    "FAIL: projected structural census did not bound first mesh wave\n");
	return 1;
    }
    const auto retainedStructuralSeed =
	apply_point_structural_distribution_plan(pointCalibration,
	    structuralSeed.threshold, projectedStructural, 10000, 1500);
    if (retainedStructuralSeed.changed ||
	!near(retainedStructuralSeed.threshold, 16.0f)) {
	std::fprintf(stderr,
	    "FAIL: projected structural seed advanced without new evidence\n");
	return 1;
    }

    pointCalibration.reset();
    auto structuralBlocked = apply_point_structural_blocked_plan(pointCalibration,
	1.0f, 5000);
    if (!structuralBlocked.changed ||
	!near(structuralBlocked.threshold, 2.0f)) {
	std::fprintf(stderr,
	    "FAIL: large structural frontier did not aggregate promptly\n");
	return 1;
    }
    const auto incompleteStructuralFrame = apply_point_completed_plan(pointCalibration,
	structuralBlocked.threshold, 10000000ULL, 60.0f, true, 4000);
    if (incompleteStructuralFrame.changed ||
	!near(incompleteStructuralFrame.threshold, 2.0f)) {
	std::fprintf(stderr,
	    "FAIL: incomplete structural frame polluted safe bracket\n");
	return 1;
    }
    structuralBlocked = apply_point_structural_blocked_plan(pointCalibration,
	structuralBlocked.threshold, 400);
    if (!structuralBlocked.changed ||
	!near(structuralBlocked.threshold, 2.5f)) {
	std::fprintf(stderr,
	    "FAIL: small structural frontier did not advance deterministically\n");
	return 1;
    }

    pointCalibration.reset();
    threshold = apply_point_structural_blocked_plan(pointCalibration,
	1.0f, 5000).threshold;
    const auto firstStructuralSafe = apply_point_completed_plan(pointCalibration,
	threshold, 10000000ULL, 60.0f, true, 0);
    if (!firstStructuralSafe.changed ||
	!near(firstStructuralSafe.threshold,
	    static_cast<float>(std::sqrt(2.0)))) {
	std::fprintf(stderr,
	    "FAIL: structural safe bracket did not begin quality recovery\n");
	return 1;
    }
    threshold = firstStructuralSafe.threshold;
    settled = false;
    for (size_t iteration = 0; iteration < 16; ++iteration) {
	BObolLodPointProxyEvidence::Decision decision;
	if (threshold < 1.6f)
	    decision = apply_point_structural_blocked_plan(pointCalibration,
		threshold, 100);
	else
	    decision = apply_point_completed_plan(pointCalibration,
		threshold, 10000000ULL, 60.0f, true, 0);
	threshold = decision.threshold;
	if (!decision.changed) {
	    settled = true;
	    break;
	}
    }
    if (!settled || threshold < 1.6f || threshold > 1.7f) {
	std::fprintf(stderr,
	    "FAIL: structural point bracket cycled (threshold=%g)\n",
	    threshold);
	return 1;
    }

    /* Losing the initiating pressure-probe latch while a large scene streams
     * results must not strand a coarse point cut.  Reusable headroom alone is
     * the durable continuation witness. */
    pointCalibration.reset();
    const auto recoveredCoarse = apply_point_completed_plan(pointCalibration,
	8.0f, 10000000ULL, 60.0f, true);
    if (!recoveredCoarse.changed || !recoveredCoarse.continueRelaxation ||
	!(recoveredCoarse.threshold < 8.0f)) {
	std::fprintf(stderr,
	    "FAIL: reusable coarse point cut did not recover quality\n");
	return 1;
    }
    /* A coarse structural classifier cut cannot become terminal merely
     * because cache/result scheduling omitted its stable replay.  The exact
     * threshold must acquire a reusable timing witness, and a changed
     * successor must acquire its own distinct witness. */
    pointCalibration.reset();
    if (!BObolLodAdmissionPlanner::pointTerminalReplayRequired(
	    BObolLodAdmissionEvidence(pointCalibration),
	    BObolLodPointCalibrationGoal::RESPONSIVE, 8.0f)) {
	std::fprintf(stderr,
	    "FAIL: unsampled coarse point cut accepted as terminal\n");
	return 1;
    }
    (void)apply_point_completed_plan(pointCalibration,
	8.0f, 45000000ULL, 20.0f, false);
    if (!BObolLodAdmissionPlanner::pointTerminalReplayRequired(
	    BObolLodAdmissionEvidence(pointCalibration),
	    BObolLodPointCalibrationGoal::RESPONSIVE, 8.0f)) {
	std::fprintf(stderr,
	    "FAIL: transient point frame settled the terminal witness\n");
	return 1;
    }
    const auto settledCoarse = apply_point_completed_plan(pointCalibration,
	8.0f, 45000000ULL, 20.0f, true);
    if (settledCoarse.changed ||
	BObolLodAdmissionPlanner::pointTerminalReplayRequired(
	    BObolLodAdmissionEvidence(pointCalibration),
	    BObolLodPointCalibrationGoal::RESPONSIVE, 8.0f)) {
	std::fprintf(stderr,
	    "FAIL: reusable point frame did not settle its exact threshold\n");
	return 1;
    }
    pointCalibration.reset();
    const auto finerCoarse = apply_point_completed_plan(pointCalibration,
	8.0f, 10000000ULL, 20.0f, true);
    if (!finerCoarse.changed ||
	!BObolLodAdmissionPlanner::pointTerminalReplayRequired(
	    BObolLodAdmissionEvidence(pointCalibration),
	    BObolLodPointCalibrationGoal::RESPONSIVE,
	    finerCoarse.threshold)) {
	std::fprintf(stderr,
	    "FAIL: changed point cut inherited its predecessor's witness\n");
	return 1;
    }
    pointCalibration.reset();
    const auto saturatedCoarse = apply_point_interrupted_plan(
	pointCalibration, 64.0f, 100000000ULL, 20.0f);
    if (saturatedCoarse.changed ||
	BObolLodAdmissionPlanner::pointTerminalReplayRequired(
	    BObolLodAdmissionEvidence(pointCalibration),
	    BObolLodPointCalibrationGoal::RESPONSIVE, 64.0f)) {
	std::fprintf(stderr,
	    "FAIL: saturated point cut retained an endless replay\n");
	return 1;
    }
    if (BObolLodAdmissionPlanner::pointRequiresReusableConfirmation(1.0f) ||
	BObolLodAdmissionPlanner::pointRequiresReusableConfirmation(
	    std::numeric_limits<float>::quiet_NaN()) ||
	BObolLodAdmissionPlanner::pointRequiresReusableConfirmation(8.0f, 1) ||
	!BObolLodAdmissionPlanner::pointRequiresReusableConfirmation(1.02f) ||
	!BObolLodAdmissionPlanner::pointRequiresReusableConfirmation(8.0f)) {
	std::fprintf(stderr,
	    "FAIL: point confirmation did not yield to structural repair\n");
	return 1;
    }
    pointCalibration.reset();
    if (!BObolLodAdmissionPlanner::pointTerminalReplayRequired(
	    BObolLodAdmissionEvidence(pointCalibration),
	    BObolLodPointCalibrationGoal::RESPONSIVE, 8.0f)) {
	std::fprintf(stderr,
	    "FAIL: structural-limit test lacked an unsettled cut\n");
	return 1;
    }
    const BObolLodAdmissionPlan structuralLimitPlan =
	BObolLodAdmissionPlanner::settlePointAtStructuralLimit(
	    BObolLodAdmissionEvidence(pointCalibration),
	    BObolLodAdmissionCursor(),
	    BObolLodPointCalibrationGoal::RESPONSIVE, 8.0f);
    pointCalibration = structuralLimitPlan.nextEvidence.pointProxy();
    if (BObolLodAdmissionPlanner::pointTerminalReplayRequired(
	    BObolLodAdmissionEvidence(pointCalibration),
	    BObolLodPointCalibrationGoal::RESPONSIVE, 8.0f)) {
	std::fprintf(stderr,
	    "FAIL: structural capacity limit did not settle point cut\n");
	return 1;
    }
    const std::array<size_t, 7> unchangedStructuralCensus = {
	100, 200, 400, 700, 900, 1000, 1000
    };
    const auto unchangedSeed = apply_point_structural_distribution_plan(
	pointCalibration, 8.0f, unchangedStructuralCensus, 1000, 900);
    if (unchangedSeed.changed ||
	BObolLodAdmissionPlanner::pointTerminalReplayRequired(
	    BObolLodAdmissionEvidence(pointCalibration),
	    BObolLodPointCalibrationGoal::RESPONSIVE, 8.0f)) {
	std::fprintf(stderr,
	    "FAIL: unchanged structural census erased terminal point witness\n");
	return 1;
    }
    if (!BObolLodAdmissionPlanner::pointCalibrationYieldsToStructuralRepair(true, true, 147, false) ||
	BObolLodAdmissionPlanner::pointCalibrationYieldsToStructuralRepair(false, true, 147, false) ||
	BObolLodAdmissionPlanner::pointCalibrationYieldsToStructuralRepair(true, false, 147, false) ||
	BObolLodAdmissionPlanner::pointCalibrationYieldsToStructuralRepair(true, true, 0, false) ||
	BObolLodAdmissionPlanner::pointCalibrationYieldsToStructuralRepair(true, true, 147, true)) {
	std::fprintf(stderr,
	    "FAIL: point calibration did not transfer a quiet structural frontier\n");
	return 1;
    }

    /* An exact structural population with no useful point successor gets one
     * hard-deadline-capacity admission, not an unlimited series of
     * preferred-cadence retries.  A changed population is a new transaction;
     * an inexact frontier or one still owned by point calibration is not a
     * structural transaction. */
    BObolLodStructuralAdmissionEvidence structuralAdmissionPolicy;
    const std::array<size_t, 7> structuralCapacityHistogram = {
	10, 20, 30, 40, 50, 60, 80
    };
    if (BObolLodAdmissionPlanner::unaggregatableStructuralCount(
	    structuralCapacityHistogram, 100) != 20) {
	std::fprintf(stderr,
	    "FAIL: structural maximum-threshold census classification\n");
	return 1;
    }
    if (BObolLodAdmissionPlanner::structuralPerOccurrenceReservation(
	    450000, 300000, 20) != 7500 ||
	BObolLodAdmissionPlanner::structuralPerOccurrenceReservation(
	    300000, 300000, 20) != 0 ||
	BObolLodAdmissionPlanner::structuralPerOccurrenceReservation(
	    300019, 300000, 20) != 0 ||
	BObolLodAdmissionPlanner::structuralPerOccurrenceReservation(
	    450000, 300000, 0) != 0) {
	std::fprintf(stderr,
	    "FAIL: structural frontier allowance was not divided exactly\n");
	return 1;
    }
    BObolLodStructuralAdmissionEvidence::Inputs structuralInputs;
    structuralInputs.viewRevision = 7;
    structuralInputs.policyRevision = 11;
    structuralInputs.frontierIdentity = 13;
    structuralInputs.unresolvedCount = 20;
    structuralInputs.activeCost = 300000;
    structuralInputs.currentBudget = 320000;
    structuralInputs.certifiedBudget = 450000;
    structuralInputs.exactProjection = true;
    structuralInputs.pointPolicyExhausted = true;
    auto structuralAdmission =
	apply_structural_admission_plan(structuralAdmissionPolicy, structuralInputs);
    if (!structuralAdmission.ownsFrontier ||
	!structuralAdmission.requestAdmission ||
	structuralAdmission.duplicateRejected ||
	structuralAdmission.budget != 450000) {
	std::fprintf(stderr,
	    "FAIL: irreducible structural frontier did not receive one trial\n");
	return 1;
    }
    structuralAdmission = apply_structural_admission_plan(structuralAdmissionPolicy, structuralInputs);
    if (!structuralAdmission.ownsFrontier ||
	structuralAdmission.requestAdmission ||
	!structuralAdmission.duplicateRejected) {
	std::fprintf(stderr,
	    "FAIL: unchanged structural frontier was rearmed\n");
	return 1;
    }
    structuralInputs.unresolvedCount = 12;
    structuralAdmission = apply_structural_admission_plan(structuralAdmissionPolicy, structuralInputs);
    if (!structuralAdmission.requestAdmission ||
	structuralAdmission.duplicateRejected) {
	std::fprintf(stderr,
	    "FAIL: changed structural population did not open a new transaction\n");
	return 1;
    }
    structuralInputs.frontierIdentity = 14;
    structuralAdmission = apply_structural_admission_plan(
	structuralAdmissionPolicy, structuralInputs);
    if (!structuralAdmission.requestAdmission ||
	structuralAdmission.duplicateRejected) {
	std::fprintf(stderr,
	    "FAIL: changed structural frontier identity did not open a new transaction\n");
	return 1;
    }
    structuralInputs.exactProjection = false;
    structuralAdmission = apply_structural_admission_plan(structuralAdmissionPolicy, structuralInputs);
    if (structuralAdmission.ownsFrontier ||
	structuralAdmission.requestAdmission) {
	std::fprintf(stderr,
	    "FAIL: inexact structural census owned admission\n");
	return 1;
    }
    structuralInputs.exactProjection = true;
    structuralInputs.pointPolicyExhausted = false;
    structuralAdmission = apply_structural_admission_plan(structuralAdmissionPolicy, structuralInputs);
    if (structuralAdmission.ownsFrontier ||
	structuralAdmission.requestAdmission) {
	std::fprintf(stderr,
	    "FAIL: structural frontier bypassed an available point successor\n");
	return 1;
    }

    /* A safe coarse small-part cut must not collapse all reducible mesh
     * prefixes.  That old coupling made Hubble alternate forever between a
     * one-pixel population and a roughly two-pixel population, starving its
     * large solar arrays of the detail budget. */
    if (BObolLodAdmissionPlanner::shouldRecoverTriangleDetail(true, false, true, false) ||
	BObolLodAdmissionPlanner::shouldRecoverTriangleDetail(false, true, true, false) ||
	BObolLodAdmissionPlanner::shouldRecoverTriangleDetail(true, true, true, true) ||
	!BObolLodAdmissionPlanner::shouldRecoverTriangleDetail(true, true, false, false) ||
	!BObolLodAdmissionPlanner::shouldRecoverTriangleDetail(true, true, true, false)) {
	std::fprintf(stderr,
	    "FAIL: point/triangle recovery ignores protected quality floor\n");
	return 1;
    }

    /* A settled event-driven mesh population goes straight to one pixel
     * under the hard static deadline.  Producer activity, resource pressure,
     * a structural fallback population, inexact/preparation frames, active
     * input, and an already rejected trial must all retain the conservative
     * convergence path. */
    if (!BObolLodAdmissionPlanner::
	    startOnePixelTrialFromSettledPointFrame(false, true, false,
		false, false, false, 4.0f, 60000000ULL, 50000000ULL,
		100000000ULL) ||
	BObolLodAdmissionPlanner::
	    startOnePixelTrialFromSettledPointFrame(true, true, false,
		false, false, false, 4.0f, 60000000ULL, 50000000ULL,
		100000000ULL) ||
	BObolLodAdmissionPlanner::
	    startOnePixelTrialFromSettledPointFrame(false, false, false,
		false, false, false, 4.0f, 60000000ULL, 50000000ULL,
		100000000ULL) ||
	BObolLodAdmissionPlanner::
	    startOnePixelTrialFromSettledPointFrame(false, true, true,
		false, false, false, 4.0f, 60000000ULL, 50000000ULL,
		100000000ULL) ||
	BObolLodAdmissionPlanner::
	    startOnePixelTrialFromSettledPointFrame(false, true, false,
		true, false, false, 4.0f, 60000000ULL, 50000000ULL,
		100000000ULL) ||
	BObolLodAdmissionPlanner::
	    startOnePixelTrialFromSettledPointFrame(false, true, false,
		false, true, false, 4.0f, 60000000ULL, 50000000ULL,
		100000000ULL) ||
	BObolLodAdmissionPlanner::
	    startOnePixelTrialFromSettledPointFrame(false, true, false,
		false, false, true, 4.0f, 60000000ULL, 50000000ULL,
		100000000ULL) ||
	BObolLodAdmissionPlanner::
	    startOnePixelTrialFromSettledPointFrame(false, true, false,
		false, false, false, 1.0f, 60000000ULL, 50000000ULL,
		100000000ULL) ||
	BObolLodAdmissionPlanner::
	    startOnePixelTrialFromSettledPointFrame(false, true, false,
		false, false, false, 4.0f, 110000000ULL, 50000000ULL,
		100000000ULL)) {
	std::fprintf(stderr,
	    "FAIL: settled point population static one-pixel trial policy\n");
	return 1;
    }

    if (std::fabs(BObolLodAdmissionPlanner::triangleRecoveryPointThreshold(16.0f, true, 4) - 16.0f) >
		0.001f ||
	std::fabs(BObolLodAdmissionPlanner::triangleRecoveryPointThreshold(16.0f, false, 0) - 16.0f) >
		0.001f ||
	std::fabs(BObolLodAdmissionPlanner::triangleRecoveryPointThreshold(16.0f, true, 0) - 1.0f) >
		0.001f) {
	std::fprintf(stderr,
	    "FAIL: triangle recovery exposed structural fallback boxes\n");
	return 1;
    }

    if (!BObolLodAdmissionPlanner::retainAggregatedPresentation(
	    -1, 16.0f, 100, true) ||
	!BObolLodAdmissionPlanner::retainAggregatedPresentation(
	    4, 16.0f, 100, false) ||
	BObolLodAdmissionPlanner::retainAggregatedPresentation(
	    -1, 16.0f, 100, false) ||
	BObolLodAdmissionPlanner::retainAggregatedPresentation(
	    -1, 1.0f, 100, true)) {
	std::fprintf(stderr,
	    "FAIL: static quality exposed structural fallback boxes\n");
	return 1;
    }

    if (Policy::progressiveCeiling(-1, 4.0f) != -1 ||
	Policy::progressiveCeiling(10, 1.0f) != -1 ||
	Policy::progressiveCeiling(10, 3.99f) != 8 ||
	Policy::progressiveCeiling(10, 2.0f, 7) != 6 ||
	Policy::progressiveCeiling(3, 64.0f) != 0 ||
	Policy::progressiveCeiling(
	    3, std::numeric_limits<float>::infinity()) != 0 ||
	Policy::progressiveCeiling(
	    3, std::numeric_limits<float>::quiet_NaN()) != -1) {
	std::fprintf(stderr, "FAIL: progressive ceiling policy\n");
	return 1;
    }

    const Policy::StableCapacityEvidence exactStable =
	Policy::stableCapacityEvidence(1000, 50000000ULL, 20.0f, true);
    const Policy::StableCapacityEvidence fasterStable =
	Policy::stableCapacityEvidence(1000, 40000000ULL, 20.0f, true);
    if (!exactStable.proven ||
	std::fabs(exactStable.renderCostPerSecond - 20000.0L) > 0.001L ||
	!fasterStable.proven ||
	std::fabs(fasterStable.renderCostPerSecond - 25000.0L) > 0.001L ||
	Policy::stableCapacityEvidence(
	    1000, 50000001ULL, 20.0f, true).proven ||
	Policy::stableCapacityEvidence(
	    1000, 40000000ULL, 20.0f, false).proven ||
	Policy::stableCapacityEvidence(
	    0, 40000000ULL, 20.0f, true).proven ||
	Policy::stableCapacityEvidence(
	    1000, 0, 20.0f, true).proven ||
	Policy::stableCapacityEvidence(
	    1000, 40000000ULL,
	    std::numeric_limits<float>::infinity(), true).proven ||
	Policy::stableCapacityEvidence(
	    1000, 40000000ULL,
	    std::numeric_limits<float>::quiet_NaN(), true).proven) {
	std::fprintf(stderr, "FAIL: reusable quality-to-stable evidence\n");
	return 1;
    }

    const size_t lucyTrial =
	Policy::discreteTrialOverBudgetAllowance(222780, 431301);
    const size_t overflowTrial =
	Policy::discreteTrialOverBudgetAllowance(
	    SIZE_MAX / 2, SIZE_MAX / 2);
    if (lucyTrial != 431301 ||
	Policy::discreteTrialOverBudgetAllowance(100, 150) != 150 ||
	Policy::discreteTrialOverBudgetAllowance(100, 500) != 0 ||
	Policy::discreteTrialOverBudgetAllowance(0, 100) != 0 ||
	Policy::discreteTrialOverBudgetAllowance(100, 0) != 0 ||
	Policy::discreteTrialOverBudgetAllowance(101, 100) != 0 ||
	Policy::discreteTrialOverBudgetAllowance(100, SIZE_MAX) != 0 ||
	overflowTrial == 0 || overflowTrial > SIZE_MAX / 2) {
	std::fprintf(stderr, "FAIL: bounded discrete population trial policy "
	    "(lucy=%zu overflow=%zu)\n", lucyTrial, overflowTrial);
	return 1;
    }
    return 0;
}

static int
test_view_demand_policy(void)
{
    using Policy = BObolLodViewDemandPolicy;
    Policy policy;
    if (Policy::staticHistoryDeadlineNanoseconds(0) != 0 ||
	Policy::staticHistoryDeadlineNanoseconds(100000000ULL) !=
	120000000ULL ||
	Policy::staticHistoryDeadlineNanoseconds(UINT64_MAX) != UINT64_MAX) {
	std::fprintf(stderr, "FAIL: bounded static-history deadline policy\n");
	return 1;
    }
    if (policy.demandRefreshActive() || policy.viewScaleChanging() ||
	policy.interactionScaleChanged() || policy.qualityProbeActive() ||
	policy.qualityProbePending() || policy.qualityProbePresented() ||
	policy.qualityBudgetActive() || policy.qualityCeilingLimit() != -1 ||
	policy.qualityCeilingFailedWork() != 0 ||
	policy.scaleChangingInteraction(true)) {
	std::fprintf(stderr, "FAIL: initial view-demand state\n");
	return 1;
    }

    Policy::CameraChangeDecision camera =
	policy.observeCameraChange(true, 10000000ULL);
    policy.refreshForViewRevision(true);
    policy.beginCameraInteraction(true, true);
    if (camera.retainPriorQualityCeiling ||
	!policy.viewScaleChanging() ||
	!policy.interactionScaleChanged() ||
	!policy.demandRefreshActive() ||
	!policy.scaleChangingInteraction(true) ||
	!policy.noteMotionFrameSettled() ||
	!policy.qualityProbePending()) {
	std::fprintf(stderr, "FAIL: scale-demand camera edge\n");
	return 1;
    }

    Policy::QualityProbeInputs probeInput;
    probeInput.activeMaximum = 6;
    probeInput.motionFramePending = true;
    Policy::QualityProbeDecision probe =
	policy.beginQualityProbe(probeInput);
    if (probe.consumed || probe.begin ||
	!policy.qualityProbePending()) {
	std::fprintf(stderr, "FAIL: quality probe bypassed frame barrier\n");
	return 1;
    }
    probeInput.motionFramePending = false;
    probe = policy.beginQualityProbe(probeInput);
    if (!probe.consumed || !probe.begin ||
	probe.progressiveCeiling != 7 ||
	!policy.qualityProbeActive() ||
	policy.qualityProbePending() ||
	policy.qualityProbePresented() ||
	!policy.qualityBudgetActive()) {
	std::fprintf(stderr, "FAIL: bounded zoom-quality probe\n");
	return 1;
    }

    policy.noteFramePresented(7, true, 100000000ULL);
    camera = policy.observeCameraChange(true, 100000000ULL);
    if (!camera.retainPriorQualityCeiling ||
	policy.qualityProbeActive() || policy.qualityProbePending() ||
	policy.qualityProbePresented() ||
	!policy.qualityBudgetActive() || policy.qualityFloor() != 7 ||
	!policy.viewScaleChanging()) {
	std::fprintf(stderr, "FAIL: responsive quality-floor retention\n");
	return 1;
    }
    policy.noteFramePresented(7, true, 100000001ULL);
    camera = policy.observeCameraChange(true, 100000001ULL);
    if (camera.retainPriorQualityCeiling) {
	std::fprintf(stderr, "FAIL: slow quality floor retained\n");
	return 1;
    }

    policy.noteQualityMiss(6, 1000);
    if (policy.qualityFloor() != 6 || policy.qualityCeilingLimit() != 6 ||
	policy.qualityCeilingFailedWork() != 1000) {
	std::fprintf(stderr, "FAIL: slow quality floor correction\n");
	return 1;
    }

    policy.beginCameraInteraction(false, true);
    if (!policy.rearmAfterResidentGrowth(true) ||
	policy.qualityProbeActive() || !policy.qualityProbePending() ||
	policy.qualityProbePresented() || !policy.qualityBudgetActive()) {
	std::fprintf(stderr, "FAIL: resident-growth quality rearm\n");
	return 1;
    }
    probeInput.activeMaximum = 15;
    probeInput.presentationCeiling = 6;
    probe = policy.beginQualityProbe(probeInput);
    if (!probe.consumed || probe.begin ||
	probe.progressiveCeiling != 7 || policy.qualityProbeActive()) {
	std::fprintf(stderr,
	    "FAIL: known-bad quality population was reprobed\n");
	return 1;
    }

    /* A cut is known-bad only for the population which missed its deadline.
     * Spatial culling after a zoom can make the same immutable prefix much
     * cheaper, in which case one deadline-bounded retry is valid. */
    if (!policy.rearmAfterResidentGrowth(true)) {
	std::fprintf(stderr, "FAIL: cheaper-view quality rearm\n");
	return 1;
    }
    probeInput.presentedWork = 400;
    probe = policy.beginQualityProbe(probeInput);
    if (!probe.consumed || !probe.begin ||
	probe.progressiveCeiling != 7 ||
	policy.qualityCeilingLimit() != -1 ||
	policy.qualityCeilingFailedWork() != 0) {
	std::fprintf(stderr, "FAIL: cheaper visible population stayed capped\n");
	return 1;
    }

    if (!policy.finishQuietInteraction() ||
	policy.qualityCeilingLimit() != -1 ||
	policy.qualityCeilingFailedWork() != 0) {
	std::fprintf(stderr, "FAIL: quiet epoch retained motion miss\n");
	return 1;
    }
    policy.beginGesture(true);
    (void)policy.observeCameraChange(true, 10000000ULL);
    policy.beginCameraInteraction(true, true);
    policy.seedQualityFloor(9);
    if (!policy.noteMotionFrameSettled()) {
	std::fprintf(stderr, "FAIL: quality-floor probe was not armed\n");
	return 1;
    }
    probeInput = Policy::QualityProbeInputs();
    probeInput.activeMaximum = 12;
    probeInput.presentationCeiling = 5;
    probe = policy.beginQualityProbe(probeInput);
    if (!probe.begin || probe.progressiveCeiling != 9) {
	std::fprintf(stderr, "FAIL: proven quality floor walked from motion cut\n");
	return 1;
    }
    policy.noteFramePresented(9, true, 90000000ULL);
    if (!policy.rearmAfterQualityFrame(
	    true, 11, 9, true, 90000000ULL) ||
	!policy.qualityProbePending() || policy.qualityProbeActive()) {
	std::fprintf(stderr, "FAIL: resident successor quality rearm\n");
	return 1;
    }
    probeInput.presentationCeiling = 9;
    probe = policy.beginQualityProbe(probeInput);
    if (!probe.begin || probe.progressiveCeiling != 10 ||
	policy.rearmAfterQualityFrame(
	    true, 11, 10, true, 100000001ULL)) {
	std::fprintf(stderr, "FAIL: quality successor deadline bound\n");
	return 1;
    }

    /* Resident growth rearms exactly one measured presentation step.  It may
     * not reopen the complete hidden prefix: a missed software frame would
     * then back down through every ordinal and later repaint edges could
     * repeat that same staircase. */
    if (!policy.rearmAfterResidentGrowth(true)) {
	std::fprintf(stderr, "FAIL: measured quality probe rearm\n");
	return 1;
    }
    probe = policy.beginQualityProbe(probeInput);
    if (!probe.begin || probe.progressiveCeiling != 10) {
	std::fprintf(stderr, "FAIL: resident quality probe skipped witness\n");
	return 1;
    }
    policy.endGesture();
    if (policy.qualityProbeActive() || policy.qualityProbePending() ||
	policy.qualityProbePresented() || !policy.qualityBudgetActive()) {
	std::fprintf(stderr, "FAIL: gesture-end view-demand state\n");
	return 1;
    }

    /* A bracketed rotation may begin before an unbracketed wheel epoch's
     * quiet debounce.  The earlier scale change remains a stable-convergence
     * obligation, but it must not classify the rotation frame or a late
     * resident suffix as another scale-quality probe. */
    policy.beginGesture(false);
    camera = policy.observeCameraChange(false, 10000000ULL);
    policy.beginCameraInteraction(false, false);
    if (camera.retainPriorQualityCeiling ||
	!policy.interactionScaleChanged() ||
	policy.scaleChangingInteraction(true) ||
	policy.rearmAfterResidentGrowth(true)) {
	std::fprintf(stderr, "FAIL: mixed zoom/pose demand isolation\n");
	return 1;
    }
    if (!policy.finishQuietInteraction() ||
	policy.viewScaleChanging() || policy.interactionScaleChanged() ||
	policy.qualityProbeActive() || policy.qualityProbePending() ||
	policy.qualityProbePresented() || policy.qualityBudgetActive()) {
	std::fprintf(stderr, "FAIL: quiet scale handoff\n");
	return 1;
    }

    policy.reset();
    (void)policy.observeCameraChange(true, 10000000ULL);
    policy.beginCameraInteraction(true, true);
    if (policy.finishQuietInteraction(true)) {
	std::fprintf(stderr, "FAIL: scale round trip required retarget\n");
	return 1;
    }

    policy.refreshForPolicyRevision(true, false);
    if (!policy.demandRefreshActive()) {
	std::fprintf(stderr, "FAIL: physical-demand policy continuation\n");
	return 1;
    }

    Policy::DemandPassInputs demandPassInputs;
    if (!policy.demandPassRequired(demandPassInputs)) {
	std::fprintf(stderr,
	    "FAIL: pending demand refresh did not require a pass\n");
	return 1;
    }
    bool Policy::DemandPassInputs::*const passOwners[] = {
	&Policy::DemandPassInputs::submissionActive,
	&Policy::DemandPassInputs::rescanPending,
	&Policy::DemandPassInputs::selectivePass,
	&Policy::DemandPassInputs::structuralRepair,
	&Policy::DemandPassInputs::retainedAllocation,
	&Policy::DemandPassInputs::strongerOwnerPending
    };
    for (bool Policy::DemandPassInputs::*owner : passOwners) {
	demandPassInputs = Policy::DemandPassInputs();
	demandPassInputs.*owner = true;
	if (!policy.demandPassRequired(demandPassInputs))
	    continue;
	std::fprintf(stderr,
	    "FAIL: demand-refresh pass restart ownership\n");
	return 1;
    }

    policy.clearDemandRefresh();
    policy.requestDemandRefresh();
    if (!policy.demandRefreshActive() ||
	!policy.demandPassRequired(Policy::DemandPassInputs())) {
	std::fprintf(stderr,
	    "FAIL: explicit result-demand refresh did not arm a pass\n");
	return 1;
    }
    if (policy.completeDemandRefresh(false, false, false, false, false) ||
	policy.completeDemandRefresh(true, true, false, false, false) ||
	policy.completeDemandRefresh(true, false, true, false, false) ||
	policy.completeDemandRefresh(true, false, false, true, false) ||
	policy.completeDemandRefresh(true, false, false, false, true) ||
	!policy.demandRefreshActive() ||
	!policy.completeDemandRefresh(true, false, false, false, false) ||
	policy.demandRefreshActive() ||
	policy.demandPassRequired(Policy::DemandPassInputs())) {
	std::fprintf(stderr, "FAIL: dense demand-refresh completion proof\n");
	return 1;
    }
    policy.beginGesture(true);
    camera = policy.observeCameraChange(false, 10000000ULL);
    policy.refreshForViewRevision(true);
    policy.beginCameraInteraction(true, false);
    if (camera.retainPriorQualityCeiling ||
	policy.demandRefreshActive() || policy.viewScaleChanging() ||
	policy.interactionScaleChanged() ||
	policy.scaleChangingInteraction(true) ||
	policy.rearmAfterResidentGrowth(true)) {
	std::fprintf(stderr, "FAIL: pose-only demand isolation\n");
	return 1;
    }

    policy.reset();
    (void)policy.observeCameraChange(true, 0);
    policy.beginCameraInteraction(true, true);
    (void)policy.noteMotionFrameSettled();
    probeInput = Policy::QualityProbeInputs();
    probe = policy.beginQualityProbe(probeInput);
    if (!probe.consumed || probe.begin || policy.qualityProbePending()) {
	std::fprintf(stderr, "FAIL: unavailable population probe retirement\n");
	return 1;
    }
    policy.reset();
    if (policy.demandRefreshActive() || policy.viewScaleChanging() ||
	policy.interactionScaleChanged() || policy.qualityBudgetActive() ||
	policy.qualityCeilingLimit() != -1) {
	std::fprintf(stderr, "FAIL: view-demand reset\n");
	return 1;
    }

    policy.beginGesture(true);
    (void)policy.observeCameraChange(true, 0);
    policy.beginCameraInteraction(true, true);
    (void)policy.noteMotionFrameSettled();
    probeInput = Policy::QualityProbeInputs();
    probeInput.activeMaximum = 31;
    probeInput.presentationCeiling = 20;
    probe = policy.beginQualityProbe(probeInput);
    if (!probe.begin || probe.progressiveCeiling != 21) {
	std::fprintf(stderr,
	    "FAIL: quality probe retained legacy 16-cut ceiling\n");
	return 1;
    }
    return 0;
}

static int
test_quiet_successor_reducer(void)
{
    using Reducer = BObolLodQuietSuccessorReducer;
    Reducer::Inputs inputs;
    inputs.retainedMeshPayloads = true;
    inputs.exactView.valid = true;
    inputs.exactView.targetPixelError = 0.5f;
    inputs.exactView.progressiveCeiling = 17;
    inputs.exactView.progressiveNextFraction = 0.625f;
    inputs.exactView.pointProxyPixelThreshold = 1.25f;
    inputs.exactView.presentedRenderCost = 7000;
    inputs.exactView.provenRenderCostCapacity = 9000;

    /* Intermediate publication order may alter the last motion controls, but
     * it cannot alter the certified quiet successor. */
    for (int ceiling = -1; ceiling <= 9; ++ceiling) {
	for (int pixel = 1; pixel <= 16; pixel *= 2) {
	    inputs.current.valid = true;
	    inputs.current.targetPixelError = static_cast<float>(pixel);
	    inputs.current.progressiveCeiling = ceiling;
	    inputs.current.progressiveNextFraction =
		ceiling < 0 ? 0.0f : static_cast<float>(ceiling % 4) / 4.0f;
	    inputs.current.pointProxyPixelThreshold =
		static_cast<float>(pixel);
	    const Reducer::Decision decision = Reducer::reduce(inputs);
	    if (decision.source != Reducer::Source::EXACT_VIEW ||
		decision.needsHandoff() || decision.clearPresentationLimits ||
		std::fabs(decision.target.targetPixelError - 0.5f) > 1.0e-6f ||
		decision.target.progressiveCeiling != 17 ||
		std::fabs(decision.target.progressiveNextFraction - 0.625f) >
		    1.0e-6f ||
		std::fabs(decision.target.pointProxyPixelThreshold - 1.25f) >
		    1.0e-6f ||
		decision.target.presentedRenderCost != 7000 ||
		decision.target.provenRenderCostCapacity != 9000) {
		std::fprintf(stderr,
		    "FAIL: quiet successor depended on transient publication\n");
		return 1;
	    }
	}
    }

    /* With no certificate, motion ceilings survive only as handoff safety
     * controls.  They cannot retain the motion pixel-error demand. */
    inputs.exactView = Reducer::Target();
    inputs.current.targetPixelError = 8.0f;
    inputs.current.progressiveCeiling = 4;
    inputs.current.progressiveNextFraction = 0.5f;
    inputs.current.pointProxyPixelThreshold = 3.0f;
    Reducer::Decision decision = Reducer::reduce(inputs);
    if (decision.source != Reducer::Source::STABLE_DEMAND ||
	decision.handoff != Reducer::Handoff::ALLOCATION || !decision.apply ||
	std::fabs(decision.target.targetPixelError - 1.0f) > 1.0e-6f ||
	decision.target.progressiveCeiling != 4 ||
	std::fabs(decision.target.progressiveNextFraction - 0.5f) > 1.0e-6f ||
	std::fabs(decision.target.pointProxyPixelThreshold - 3.0f) > 1.0e-6f) {
	std::fprintf(stderr, "FAIL: quiet stable-demand handoff\n");
	return 1;
    }

    /* Stable controls survive an orthographic zoom, but the old occurrence
     * cuts do not.  The first quiet allocation must therefore start from the
     * stable threshold without claiming pose-only presentation reuse. */
    inputs.priorStable.valid = true;
    inputs.priorStable.targetPixelError = 0.5f;
    inputs.priorStable.progressiveCeiling = 8;
    inputs.priorStable.progressiveNextFraction = 0.625f;
    inputs.priorStable.pointProxyPixelThreshold = 1.25f;
    inputs.priorStableCutsReusable = false;
    decision = Reducer::reduce(inputs);
    if (decision.source != Reducer::Source::PRIOR_STABLE ||
	decision.handoff != Reducer::Handoff::ALLOCATION ||
	std::fabs(decision.target.targetPixelError - 0.5f) > 1.0e-6f ||
	decision.target.progressiveCeiling != 8 ||
	std::fabs(decision.target.pointProxyPixelThreshold - 1.25f) > 1.0e-6f) {
	std::fprintf(stderr, "FAIL: quiet zoom lost stable controls\n");
	return 1;
    }

    /* A completed current-scale frame is stronger than the pre-zoom control
     * vector, regardless of which one was recorded first. */
    inputs.provenScale.valid = true;
    inputs.provenScale.targetPixelError = 0.75f;
    inputs.provenScale.progressiveCeiling = 11;
    inputs.provenScale.pointProxyPixelThreshold = 1.0f;
    decision = Reducer::reduce(inputs);
    if (decision.source != Reducer::Source::PROVEN_SCALE ||
	decision.target.progressiveCeiling != 11 ||
	std::fabs(decision.target.targetPixelError - 0.5f) > 1.0e-6f ||
	std::fabs(decision.target.pointProxyPixelThreshold - 1.25f) >
	    1.0e-6f) {
	std::fprintf(stderr, "FAIL: current-scale proof lost certificate precedence\n");
	return 1;
    }

    inputs.retainedMeshPayloads = false;
    decision = Reducer::reduce(inputs);
    if (decision.source != Reducer::Source::STABLE_DEMAND ||
	decision.needsHandoff() || !decision.clearPresentationLimits ||
	decision.target.progressiveCeiling != -1 ||
	std::fabs(decision.target.pointProxyPixelThreshold - 1.0f) > 1.0e-6f) {
	std::fprintf(stderr, "FAIL: quiet payload-free successor\n");
	return 1;
    }
    return 0;
}

static int
test_presentation_policy(void)
{
    using Policy = BObolLodPresentationPolicy;
    if (Policy::presentationLimitsReconciled(
	    false, true, 900, 1000) ||
	Policy::presentationLimitsReconciled(
	    true, false, 900, 1000) ||
	Policy::presentationLimitsReconciled(
	    true, true, 1001, 1000) ||
	Policy::presentationLimitsReconciled(
	    true, true, 900, 1000, 800) ||
	!Policy::presentationLimitsReconciled(
	    true, true, 1000, 1000)) {
	std::fprintf(stderr,
	    "FAIL: retained presentation reconciliation budget predicate\n");
	return 1;
    }

    Policy policy;
    Policy::Population population;
    population.available = true;
    population.sceneDomainRevision = 3;

    Policy::RestoreInputs input;
    input.orthographic = true;
    input.retainedMeshPayloads = true;
    input.viewEpoch.set(11);
    input.population = population;
    input.currentTargetPixelError = 2.0f;
    input.currentProgressiveCeiling = 4;
    input.currentProgressiveNextFraction = 0.25f;
    input.currentPointProxyPixelThreshold = 3.0f;

    policy.capturePrior(0.5f, 8, 0.625f, 2.0f, population,
	BObolLodViewEpoch(10));
    Policy::RestoreDecision restore = policy.beginQuiet(input);
	if (!restore.apply || !restore.restoredPriorStable ||
	std::fabs(restore.targetPixelError - 0.5f) > 0.0001f ||
	restore.progressiveCeiling != 8 ||
	std::fabs(restore.progressiveNextFraction - 0.625f) > 0.0001f ||
	std::fabs(restore.pointProxyPixelThreshold - 2.0f) > 0.0001f ||
	restore.handoff != Policy::RestoreDecision::Handoff::PRESENTATION ||
	!policy.handoffPending() || !policy.handoffPresentationPending() ||
	policy.handoffCostFloor() != 0) {
	std::fprintf(stderr,
	    "FAIL: quiet pose restore did not require current-pose proof\n");
	return 1;
    }

    policy.reset();
    Policy::QuietInputs exactInputs;
    exactInputs.presentation = input;
    exactInputs.exactView.valid = true;
    exactInputs.exactView.targetPixelError = 0.25f;
    exactInputs.exactView.progressiveCeiling = 19;
    exactInputs.exactView.progressiveNextFraction = 0.75f;
    exactInputs.exactView.pointProxyPixelThreshold = 1.0f;
    exactInputs.exactView.presentedRenderCost = 8000;
    exactInputs.exactView.provenRenderCostCapacity = 10000;
    restore = policy.beginQuiet(exactInputs);
    if (!restore.apply || !restore.restoredExactView ||
	restore.restoredPriorStable || restore.restoredProvenQuality ||
	restore.needsHandoff() || policy.handoffPending() ||
	std::fabs(restore.targetPixelError - 0.25f) > 0.0001f ||
	restore.progressiveCeiling != 19 ||
	std::fabs(restore.progressiveNextFraction - 0.75f) > 0.0001f ||
	restore.provenRenderCostFloor != 8000 ||
	restore.provenRenderCostCapacity != 10000) {
	std::fprintf(stderr, "FAIL: exact-view quiet successor\n");
	return 1;
    }

    /* A late interaction-quiet edge may arrive after an offline caller has
     * entered unbounded terminal convergence.  It must clear presentation
     * limits without restoring history or creating a finite allocation
     * handoff that the terminal allocator cannot satisfy. */
    policy.armHandoff(true, 4000, 5000);
    Policy::QuietInputs terminalInputs = exactInputs;
    terminalInputs.unboundedTerminal = true;
    terminalInputs.presentation.currentTargetPixelError = 0.25f;
    restore = policy.beginQuiet(terminalInputs);
    if (!restore.apply || !restore.clearPresentationLimits ||
	restore.restoredExactView || restore.restoredPriorStable ||
	restore.restoredProvenQuality || restore.needsHandoff() ||
	policy.handoffPending() ||
	std::fabs(restore.targetPixelError - 0.25f) > 0.0001f ||
	restore.progressiveCeiling != -1 ||
	std::fabs(restore.progressiveNextFraction) > 0.0001f ||
	std::fabs(restore.pointProxyPixelThreshold - 1.0f) > 0.0001f) {
	std::fprintf(stderr, "FAIL: terminal quiet transition restored capacity policy\n");
	return 1;
    }

    policy.reset();
    policy.capturePrior(0.5f, 8, 0.625f, 2.0f, population,
	BObolLodViewEpoch(10));
    input.population.sceneDomainRevision++;
    restore = policy.beginQuiet(input);
    if (!restore.apply || !restore.restoredPriorStable ||
	restore.handoff != Policy::RestoreDecision::Handoff::ALLOCATION ||
	!policy.handoffPending() ||
	policy.handoffCostFloor() != 0 ||
	std::fabs(restore.targetPixelError - 0.5f) > 0.0001f ||
	restore.progressiveCeiling != input.currentProgressiveCeiling ||
	std::fabs(restore.progressiveNextFraction -
	    input.currentProgressiveNextFraction) > 0.0001f ||
	std::fabs(restore.pointProxyPixelThreshold - 2.0f) > 0.0001f) {
	std::fprintf(stderr,
	    "FAIL: population growth lost stable demand controls\n");
	return 1;
    }
    Policy::CompletedPassInputs completed;
    completed.completed = true;
    completed.populationQuiescent = true;

    /* Structural source records may remain behind an exact subpixel aggregate
     * without leaving any repairable fallback in the framebuffer.  Only the
     * latter is structural-frontier debt. */
    if (Policy::nonterminalStructuralFrontier(true, 0, 0) ||
	!Policy::nonterminalStructuralFrontier(true, 1, 0) ||
	Policy::nonterminalStructuralFrontier(false, 1, 0) ||
	Policy::nonterminalStructuralFrontier(true, 1, 1)) {
	std::fprintf(stderr,
	    "FAIL: nonterminal structural frontier classification\n");
	return 1;
    }
    if (Policy::presentationOnlyFrameAdvancesPlanning(
	    false, true, true, true) ||
	Policy::presentationOnlyFrameAdvancesPlanning(
	    true, false, false, false) ||
	!Policy::presentationOnlyFrameAdvancesPlanning(
	    true, true, false, false) ||
	!Policy::presentationOnlyFrameAdvancesPlanning(
	    true, false, true, false) ||
	!Policy::presentationOnlyFrameAdvancesPlanning(
	    true, false, false, true)) {
	std::fprintf(stderr,
	    "FAIL: presentation-only planning-edge classification\n");
	return 1;
    }

    /* A targeted source/delta pass and an allocation completed while worker
     * results are still enlarging the population are not handoff proofs. */
    Policy::CompletedPassDecision completion = policy.completePass(completed);
    if (completion.finishHandoff ||
	!completion.requestRetainedAllocation ||
	!policy.handoffPending()) {
	std::fprintf(stderr,
	    "FAIL: non-allocation pass did not request scene allocation\n");
	return 1;
    }

    /* A retained allocation cannot repair records which exist only as
     * structural fallbacks.  Leave the handoff armed but yield its planning
     * edge until the structural owner publishes those occurrences. */
    completed.structuralFrontierPending = true;
    completed.retainedRefinementPending = true;
    completion = policy.completePass(completed);
    if (completion.finishHandoff || completion.requestRetainedAllocation ||
	completion.requestLocalPresentationReduction ||
	!completion.retireRetainedObservation ||
	!policy.handoffPending()) {
	std::fprintf(stderr,
	    "FAIL: allocation handoff preempted structural frontier repair\n");
	return 1;
    }
    /* Structural ownership is stronger than a stale capacity/budget
     * annotation from the ordinary pass.  The exact frontier has its own
     * bounded successor and must be able to retire that observation even when
     * no handoff is active. */
    policy.reset();
    completed = Policy::CompletedPassInputs();
    completed.completed = true;
    completed.structuralFrontierPending = true;
    completed.capacityTransactionPending = true;
    completed.retainedRefinementPending = true;
    completed.retainedRefinementBudgetBlocked = true;
    completion = policy.completePass(completed);
    if (completion.finishHandoff || completion.requestRetainedAllocation ||
	completion.requestLocalPresentationReduction ||
	!completion.retireRetainedObservation || policy.handoffPending()) {
	std::fprintf(stderr,
	    "FAIL: structural owner retained stale pass annotations\n");
	return 1;
    }
    policy.armHandoff(false);
    completed.structuralFrontierPending = false;
    completed.retainedRefinementPending = false;
    completed.retainedRefinementBudgetBlocked = false;
    completed.capacityTransactionPending = false;
    completed.retainedAllocationCompleted = true;
    completed.populationQuiescent = false;
    completion = policy.completePass(completed);
    if (completion.finishHandoff || completion.requestRetainedAllocation ||
	!policy.handoffPending()) {
	std::fprintf(stderr,
	    "FAIL: growing population completed presentation handoff\n");
	return 1;
    }
    completed.populationQuiescent = true;
    completed.changedCut = true;
    if (policy.completePass(completed).finishHandoff) {
	std::fprintf(stderr, "FAIL: changed cut completed handoff\n");
	return 1;
    }
    completed.changedCut = false;
    completed.retainedRefinementPending = true;
    completed.retainedAllocationCertified = true;
    completed.presentationLimitsReconciled = true;
    completion = policy.completePass(completed);
    if (!completion.finishHandoff || !completion.requestRetainedRescan ||
	completion.retireRetainedObservation || policy.handoffPending()) {
	std::fprintf(stderr, "FAIL: unchanged pass did not complete handoff\n");
	return 1;
    }

    /* A quiet render-deadline correction has not presented its constrained
     * cut yet.  An unchanged planning pass cannot remove that ceiling; only a
     * completed presentation supplies the missing proof. */
    policy.reset();
    policy.armHandoff(true);
    completed = Policy::CompletedPassInputs();
    completed.completed = true;
    completed.retainedAllocationCompleted = true;
    completed.retainedAllocationCertified = true;
    completed.presentationLimitsReconciled = true;
    completed.populationQuiescent = true;
    completion = policy.completePass(completed);
    if (completion.finishHandoff || !policy.handoffPending() ||
	!policy.handoffPresentationPending()) {
	std::fprintf(stderr, "FAIL: unpresented deadline handoff completed\n");
	return 1;
    }
    if (!policy.noteFramePresented(700)) {
	std::fprintf(stderr, "FAIL: deadline frame did not release handoff barrier\n");
	return 1;
    }
    if (policy.handoffCostFloor() != 700) {
	std::fprintf(stderr,
	    "FAIL: deadline frame did not retain proven cost floor\n");
	return 1;
    }
    completed.retainedAllocationCertified = false;
    completed.presentationLimitsReconciled = false;
    completion = policy.completePass(completed);
    if (completion.finishHandoff ||
	!completion.requestRetainedAllocation || !policy.handoffPending()) {
	std::fprintf(stderr,
	    "FAIL: unproved deadline handoff did not request allocation\n");
	return 1;
    }
    completed.retainedAllocationCertified = true;
    completion = policy.completePass(completed);
    if (completion.finishHandoff ||
	!completion.requestLocalPresentationReduction ||
	!policy.handoffPending()) {
	std::fprintf(stderr,
	    "FAIL: over-budget deadline handoff became a terminal ceiling\n");
	return 1;
    }
    completed.presentationLimitsReconciled = true;
    completion = policy.completePass(completed);
    if (!completion.finishHandoff || policy.handoffPending() ||
	policy.handoffPresentationPending() ||
	completion.requestRetainedRescan) {
	std::fprintf(stderr, "FAIL: certified deadline handoff not completed\n");
	return 1;
    }

    /* Once a complete importance allocation fits the constrained frame's
     * measured budget, the global ceiling is redundant even when staging did
     * not need to alter a cut. */
    policy.reset();
    policy.armHandoff(true);
    if (!policy.noteFramePresented(700, 1200) ||
	policy.handoffReconciliationBudget() != 1200 ||
	policy.allocationReconciliationBudget(false) != 1200 ||
	policy.allocationReconciliationBudget(true) != 0) {
	std::fprintf(stderr,
	    "FAIL: deadline frame did not retain reconciliation budget\n");
	return 1;
    }
    completed = Policy::CompletedPassInputs();
    completed.completed = true;
    completed.retainedAllocationCompleted = true;
    completed.retainedAllocationCertified = true;
    completed.populationQuiescent = true;
    completed.presentationLimitsReconciled = true;
    completion = policy.completePass(completed);
    if (!completion.finishHandoff ||
	policy.handoffPending() || policy.handoffReconciliationBudget() != 0) {
	std::fprintf(stderr,
	    "FAIL: bounded occurrence reconciliation retained global ceiling\n");
	return 1;
    }

    /* A fast renderer-limited fallback cannot overrule the safe retained
     * population budget computed from the full frame which just missed its
     * deadline.  The two samples describe different cuts. */
    policy.reset();
    policy.armHandoff(true, 0, 800);
    if (!policy.noteFramePresented(700, 1200) ||
	policy.handoffReconciliationBudget() != 800) {
	std::fprintf(stderr,
	    "FAIL: coarse fallback enlarged deadline reconciliation budget\n");
	return 1;
    }
    completed = Policy::CompletedPassInputs();
    completed.completed = true;
    completed.retainedAllocationCompleted = true;
    completed.retainedAllocationCertified = true;
    completed.populationQuiescent = true;
    completed.presentationLimitsReconciled = true;
    if (!policy.completePass(completed).finishHandoff ||
	policy.handoffPending()) {
	std::fprintf(stderr,
	    "FAIL: limited deadline reconciliation did not complete\n");
	return 1;
    }

    /* A corrected renderer frame which completes just beyond its endpoint
	* has no new safe timing sample.  It must nevertheless preserve the
	* recovery limit derived from the interrupted richer frame, or the
	* handoff will remove its ceiling and retry that known miss forever. */
    policy.reset();
    policy.armHandoff(true, 0, 800);
    if (!policy.noteFramePresented(0, 0) ||
	policy.handoffPresentationPending() ||
	policy.handoffReconciliationBudget() != 800) {
	std::fprintf(stderr,
	    "FAIL: timer-edge fallback discarded deadline recovery limit\n");
	return 1;
    }
    completed = Policy::CompletedPassInputs();
    completed.completed = true;
    completed.retainedAllocationCompleted = true;
    completed.retainedAllocationCertified = true;
    completed.populationQuiescent = true;
    completion = policy.completePass(completed);
    if (!completion.requestLocalPresentationReduction ||
	completion.finishHandoff || !policy.handoffPending()) {
	std::fprintf(stderr,
	    "FAIL: timer-edge fallback accepted an unreconciled population\n");
	return 1;
    }

    /* Changing some occurrence cuts is not proof that the complete hidden
	* population fits the constrained frame's deadline.  Without an explicit
	* reconciliation witness, request an occurrence-local reduction; a global
	* ceiling is never a terminal quality policy. */
    policy.reset();
    policy.armHandoff(true);
    completed.presentationLimitsReconciled = false;
    if (!policy.noteFramePresented()) {
	std::fprintf(stderr, "FAIL: changed deadline handoff frame proof\n");
	return 1;
    }
    completed.changedCut = true;
    if (policy.completePass(completed).finishHandoff) {
	std::fprintf(stderr, "FAIL: changed deadline cut completed handoff\n");
	return 1;
    }
    completed.changedCut = false;
    completion = policy.completePass(completed);
    if (completion.finishHandoff ||
	!completion.requestLocalPresentationReduction ||
	!policy.handoffPending()) {
	std::fprintf(stderr,
	    "FAIL: unproven deadline handoff retained a terminal renderer ceiling\n");
	return 1;
    }

    completed.presentationLimitsReconciled = true;
    completion = policy.completePass(completed);
    if (!completion.finishHandoff ||
	policy.handoffPending()) {
	std::fprintf(stderr,
	    "FAIL: reconciled deadline handoff retained renderer ceiling\n");
	return 1;
    }

    policy.reset();
    policy.armHandoff(true);
    if (!policy.noteFramePresented(700, 1200)) {
	std::fprintf(stderr, "FAIL: changed reconciled deadline frame proof\n");
	return 1;
    }
    completed.changedCut = true;
    completed.presentationLimitsReconciled = true;
    if (policy.completePass(completed).finishHandoff) {
	std::fprintf(stderr,
	    "FAIL: changed reconciled deadline cut completed handoff early\n");
	return 1;
    }
    completed.changedCut = false;
    completion = policy.completePass(completed);
    if (!completion.finishHandoff ||
	policy.handoffPending()) {
	std::fprintf(stderr,
	    "FAIL: proven changed deadline reconciliation retained ceiling\n");
	return 1;
    }
    completed.presentationLimitsReconciled = false;

    /* A fully applied allocation above its certified budget is a typed
	* presentation-reconciliation successor, not an invitation to run the
	* same generic capacity allocation again. */
    policy.reset();
    if (!policy.claimOverBudgetAllocation(true, true, 1201, 1200, 700) ||
	!policy.handoffPending() || policy.handoffPresentationPending() ||
	policy.handoffCostFloor() != 700 ||
	policy.handoffReconciliationBudget() != 1200) {
	std::fprintf(stderr,
	    "FAIL: over-budget allocation did not establish handoff ownership\n");
	return 1;
    }
    completed = Policy::CompletedPassInputs();
    completed.completed = true;
    completed.retainedAllocationCompleted = true;
    completed.retainedAllocationCertified = true;
    completed.populationQuiescent = true;
    completion = policy.completePass(completed);
    if (!completion.requestLocalPresentationReduction ||
	completion.finishHandoff || !policy.handoffPending()) {
	std::fprintf(stderr,
	    "FAIL: claimed over-budget allocation did not request reduction\n");
	return 1;
    }
    policy.reset();
    if (policy.claimOverBudgetAllocation(true, true, 1200, 1200, 700) ||
	policy.claimOverBudgetAllocation(false, true, 1201, 1200, 700) ||
	policy.claimOverBudgetAllocation(true, false, 1201, 1200, 700) ||
	policy.handoffPending()) {
	std::fprintf(stderr,
	    "FAIL: fitting or incomplete allocation claimed handoff ownership\n");
	return 1;
    }

    policy.reset();
    policy.armHandoff(false);
    completed.presentationLimitsReconciled = true;
    completed.changedCut = true;
    completed.retainedRefinementPending = true;
    completed.retainedRefinementBudgetBlocked = true;
    if (policy.completePass(completed).finishHandoff) {
	std::fprintf(stderr, "FAIL: changed budget handoff completed early\n");
	return 1;
    }
    completed.changedCut = false;
    completion = policy.completePass(completed);
    if (!completion.finishHandoff ||
	completion.requestRetainedRescan || policy.handoffPending()) {
	std::fprintf(stderr,
	    "FAIL: reconciled budget handoff retained global ceiling\n");
	return 1;
    }
    completed.retainedRefinementPending = false;
    completed.retainedRefinementBudgetBlocked = false;

    /* A static staircase has already presented its guarded cut.  Its exact
	* predicted work allowance must be available to the occurrence allocator
	* immediately rather than being replaced by the broader scene budget. */
    policy.reset();
    policy.armHandoff(false, 700, 1200);
    if (policy.handoffCostFloor() != 700 ||
	policy.handoffReconciliationBudget() != 1200 ||
	policy.handoffPresentationPending()) {
	std::fprintf(stderr,
	    "FAIL: presented static handoff lost reconciliation budget\n");
	return 1;
    }

    /* A terminal capacity certificate is measured only after the global
     * presentation ceiling is absent.  Its exact population supersedes an
     * older, lower reconciliation request; retaining that request would
     * restart the allocator below the newly certified fixed point. */
    policy.acceptCapacityCertificate();
    if (policy.handoffPending() || policy.handoffCostFloor() != 0 ||
	policy.handoffReconciliationBudget() != 0) {
	std::fprintf(stderr,
	    "FAIL: capacity certificate retained stale presentation handoff\n");
	return 1;
    }

    policy.armHandoff(false, 700, 1200);
    completed = Policy::CompletedPassInputs();
    completed.completed = true;
    completed.retainedAllocationCompleted = true;
    completed.retainedAllocationCertified = true;
    completed.populationQuiescent = true;
    completed.presentationLimitsReconciled = true;
    completion = policy.completePass(completed);
    if (!completion.finishHandoff || policy.handoffPending() ||
	policy.handoffReconciliationBudget() != 0) {
	std::fprintf(stderr,
	    "FAIL: presented static reconciliation did not complete\n");
	return 1;
    }

    /* A richer retained target with no handoff, cut, budget barrier, or
     * rescan is terminal quality state, not self-sustaining background work. */
    completed = Policy::CompletedPassInputs();
    completed.completed = true;
    completed.retainedRefinementPending = true;
    completion = policy.completePass(completed);
    if (!completion.retireRetainedObservation ||
	completion.finishHandoff || completion.requestRetainedRescan) {
	std::fprintf(stderr, "FAIL: unwitnessed retained observation not retired\n");
	return 1;
    }
    completed.retainedRefinementBudgetBlocked = true;
    completion = policy.completePass(completed);
    if (completion.retireRetainedObservation) {
	std::fprintf(stderr, "FAIL: budget-owned retained observation retired\n");
	return 1;
    }
    completed.retainedRefinementBudgetBlocked = false;
    completed.changedCut = true;
    completion = policy.completePass(completed);
    if (completion.retireRetainedObservation) {
	std::fprintf(stderr, "FAIL: frame-owned retained observation retired\n");
	return 1;
    }

    /* A quiet handoff may discover that the desired retained cut is beyond
     * the calibrated renderer budget.  The constrained cut is then terminal:
     * finish the handoff without manufacturing an identical rescan. */
    policy.reset();
    policy.armHandoff(false);
    completed = Policy::CompletedPassInputs();
    completed.completed = true;
    completed.retainedAllocationCompleted = true;
    completed.retainedAllocationCertified = true;
    completed.populationQuiescent = true;
    completed.retainedRefinementPending = true;
    completed.retainedRefinementBudgetBlocked = true;
    completed.presentationLimitsReconciled = true;
    completion = policy.completePass(completed);
    if (!completion.finishHandoff || completion.requestRetainedRescan ||
	completion.retireRetainedObservation || policy.handoffPending()) {
	std::fprintf(stderr,
	    "FAIL: budget-limited handoff scheduled another refinement pass\n");
	return 1;
    }

    /* Successor ownership is selected before generic capacity calibration can
     * install a reallocation request.  An allocation handoff owns a clean
     * completed pass while it waits for either a certificate or quiescence;
     * concrete cut/rescan work retains precedence. */
    policy.reset();
    policy.armHandoff(false);
    completed = Policy::CompletedPassInputs();
    completed.completed = true;
    if (policy.completedPassSelection(completed, true, false, false, -1).owner !=
	    Policy::CompletedPassOwner::ALLOCATION_HANDOFF) {
	std::fprintf(stderr,
	    "FAIL: allocation handoff did not claim clean completed pass\n");
	return 1;
    }
    completed.changedCut = true;
    if (policy.completedPassSelection(completed, true, false, false, -1).owner !=
	    Policy::CompletedPassOwner::CAPACITY) {
	std::fprintf(stderr,
	    "FAIL: allocation handoff hid changed-cut successor\n");
	return 1;
    }
    completed.changedCut = false;
    completed.capacityTransactionPending = true;
    if (policy.completedPassSelection(completed, true, false, false, -1).owner !=
	    Policy::CompletedPassOwner::CAPACITY) {
	std::fprintf(stderr,
	    "FAIL: allocation handoff hid frame-rescan successor\n");
	return 1;
    }
    completed.completed = false;
    completed.capacityTransactionPending = false;

    if (policy.completedPassSelection(completed, true, false, false, -1).owner !=
	    Policy::CompletedPassOwner::NONE) {
	std::fprintf(stderr,
	    "FAIL: incomplete pass selected a successor owner\n");
	return 1;
    }
    completed.completed = true;
    policy.reset();
    if (policy.completedPassSelection(completed, true, false, false, -1).owner !=
	    Policy::CompletedPassOwner::CAPACITY) {
	std::fprintf(stderr,
	    "FAIL: inactive handoff claimed completed pass\n");
	return 1;
    }
    if (policy.completedPassSelection(completed, false, false, false, -1).owner !=
	    Policy::CompletedPassOwner::NONE) {
	std::fprintf(stderr,
	    "FAIL: ineligible capacity path claimed completed pass\n");
	return 1;
    }
    completed.retainedRefinementPending = true;
    if (!policy.completedPassSelection(completed, false, false, false, -1).
	    retainedObservationOwns()) {
	std::fprintf(stderr,
	    "FAIL: retained observation had no completed-pass owner\n");
	return 1;
    }
    completed.retainedRefinementPending = false;
    policy.armHandoff(true);
    if (policy.completedPassSelection(completed, true, false, false, -1).owner !=
	    Policy::CompletedPassOwner::PRESENTATION_HANDOFF) {
	std::fprintf(stderr,
	    "FAIL: presentation handoff did not claim clean completed pass\n");
	return 1;
    }

    /* Exhaust the immutable selector's concrete Boolean domain.  This is the
     * executable refinement map for ObolCompletedPassOwnership.tla: ceiling
     * removal is part of owner selection, never a later cancellation of an
     * already executed capacity effect. */
    enum SelectionBit : unsigned int {
	PASS_COMPLETED = 0,
	SUBMISSION_PENDING,
	RESCAN_AFTER_FRAME,
	CHANGED_CUT,
	ALLOCATION_COMPLETED,
	ALLOCATION_CERTIFIED,
	POPULATION_QUIESCENT,
	LIMITS_RECONCILED,
	CAPACITY_ELIGIBLE,
	CAPACITY_SAMPLE_PENDING,
	CAPACITY_PRESENTATION_PENDING,
	STRUCTURAL_FRONTIER_PENDING,
	SELECTION_BIT_COUNT
    };
    const unsigned int selectionCombinationCount =
	1u << SELECTION_BIT_COUNT;
    constexpr unsigned int HandoffStateCount = 3;
    constexpr unsigned int PresentationHandoffState = 1;
    constexpr unsigned int AllocationHandoffState = 2;
    for (unsigned int handoffState = 0;
	    handoffState < HandoffStateCount; ++handoffState) {
	for (int ceiling : {-1, 0}) {
	    for (unsigned int combination = 0;
		    combination < selectionCombinationCount; ++combination) {
		auto enabled = [combination](SelectionBit bit) {
		    return (combination & (1u << bit)) != 0;
		};
		Policy candidate;
		if (handoffState == PresentationHandoffState)
		    candidate.armHandoff(true);
		else if (handoffState == AllocationHandoffState)
		    candidate.armHandoff(false);
		Policy::CompletedPassInputs inputs;
		inputs.completed = enabled(PASS_COMPLETED);
		inputs.submissionPending = enabled(SUBMISSION_PENDING);
		inputs.capacityTransactionPending = enabled(RESCAN_AFTER_FRAME);
		inputs.structuralFrontierPending =
		    enabled(STRUCTURAL_FRONTIER_PENDING);
		inputs.changedCut = enabled(CHANGED_CUT);
		inputs.retainedAllocationCompleted =
		    enabled(ALLOCATION_COMPLETED);
		inputs.retainedAllocationCertified =
		    enabled(ALLOCATION_CERTIFIED);
		inputs.populationQuiescent = enabled(POPULATION_QUIESCENT);
		inputs.presentationLimitsReconciled =
		    enabled(LIMITS_RECONCILED);
		const bool capacityEligible = enabled(CAPACITY_ELIGIBLE);
		const bool capacitySamplePending =
		    enabled(CAPACITY_SAMPLE_PENDING);
		const bool capacityPresentationPending =
		    enabled(CAPACITY_PRESENTATION_PENDING);
		const bool consumeAnnotations =
		    !inputs.structuralFrontierPending &&
		    capacityPresentationPending && ceiling >= 0 &&
		    (inputs.capacityTransactionPending || inputs.changedCut) &&
		    inputs.completed &&
		    !inputs.submissionPending &&
		    handoffState == AllocationHandoffState &&
		    inputs.retainedAllocationCompleted &&
		    inputs.retainedAllocationCertified &&
		    inputs.populationQuiescent &&
		    inputs.presentationLimitsReconciled;
		const bool deferredCapacitySuccessor =
		    consumeAnnotations && capacitySamplePending;
		const bool cleanPass =
		    (!inputs.capacityTransactionPending && !inputs.changedCut) ||
		    consumeAnnotations;
		Policy::CompletedPassOwner expectedOwner =
		    Policy::CompletedPassOwner::NONE;
		if (inputs.completed && !inputs.submissionPending)
		    expectedOwner = capacityEligible ?
			Policy::CompletedPassOwner::CAPACITY :
			Policy::CompletedPassOwner::NONE;
		if (inputs.completed && !inputs.submissionPending &&
			inputs.structuralFrontierPending)
		    expectedOwner =
			Policy::CompletedPassOwner::STRUCTURAL_FRONTIER;
		else if (inputs.completed && !inputs.submissionPending &&
			cleanPass && handoffState != 0 &&
			(!capacitySamplePending || consumeAnnotations))
		    expectedOwner = handoffState == PresentationHandoffState ?
			    Policy::CompletedPassOwner::PRESENTATION_HANDOFF :
			    Policy::CompletedPassOwner::ALLOCATION_HANDOFF;
		const Policy::CompletedPassSelection selection =
		    candidate.completedPassSelection(inputs, capacityEligible,
			capacitySamplePending, capacityPresentationPending,
			ceiling);
		if (selection.owner != expectedOwner ||
			selection.consumePassAnnotations != consumeAnnotations ||
			selection.deferredCapacitySuccessor !=
			    deferredCapacitySuccessor ||
			selection.capacityOwns() !=
			    (expectedOwner ==
			     Policy::CompletedPassOwner::CAPACITY) ||
			selection.structuralFrontierOwns() !=
			    (expectedOwner ==
			     Policy::CompletedPassOwner::STRUCTURAL_FRONTIER) ||
			selection.retainedObservationOwns() !=
			    (expectedOwner ==
			     Policy::CompletedPassOwner::RETAINED_OBSERVATION) ||
			selection.handoffOwns() !=
			    (expectedOwner ==
				 Policy::CompletedPassOwner::PRESENTATION_HANDOFF ||
			     expectedOwner ==
				 Policy::CompletedPassOwner::ALLOCATION_HANDOFF)) {
		    std::fprintf(stderr,
			"FAIL: completed-pass selector combination=%u "
			"handoff=%u ceiling=%d\n",
			combination, handoffState, ceiling);
		    return 1;
		}
	    }
	}
    }

    for (size_t payloadCount : {size_t(0), size_t(1), size_t(2),
	    size_t(50000)}) {
	for (bool deferred : {false, true}) {
	    for (bool finished : {false, true}) {
		const bool expected = deferred && finished && payloadCount > 1;
		if (Policy::capacityProducerRequiredAfterHandoff(
			deferred, finished, payloadCount) != expected) {
		    std::fprintf(stderr,
			"FAIL: deferred capacity producer boundary payloads=%zu "
			"deferred=%d finished=%d\n", payloadCount,
			deferred ? 1 : 0, finished ? 1 : 0);
		    return 1;
		}
	    }
	}
    }

    /* An exact capacity sample cannot be taken through a renderer-wide
     * ceiling.  A complete handoff allocation removes that ceiling first,
     * while incomplete or unproved populations retain normal rescan order. */
    constexpr unsigned int CapacityPendingInputCount = 3;
    constexpr unsigned int CapacityPendingCombinationCount =
	1u << CapacityPendingInputCount;
    for (unsigned int combination = 0;
	    combination < CapacityPendingCombinationCount; ++combination) {
	const bool searchAwaiting = (combination & 1u) != 0;
	const bool allocationCurrent = (combination & 2u) != 0;
	const bool presentationRealized = (combination & 4u) != 0;
	const bool expected = searchAwaiting ||
	    (allocationCurrent && !presentationRealized);
	if (Policy::capacitySamplePending(searchAwaiting, allocationCurrent,
		presentationRealized) != expected) {
	    std::fprintf(stderr,
		"FAIL: predicted capacity sample combination=%u\n",
		combination);
	    return 1;
	}
    }
    policy.reset();
    policy.armHandoff(false);
    completed = Policy::CompletedPassInputs();
    completed.completed = true;
    completed.capacityTransactionPending = true;
    completed.retainedAllocationCompleted = true;
    completed.retainedAllocationCertified = true;
    completed.populationQuiescent = true;
    completed.presentationLimitsReconciled = true;
    if (!policy.capacitySampleRequiresCeilingFreeHandoff(
	    completed, true, 4)) {
	std::fprintf(stderr,
	    "FAIL: exact capacity sample did not yield to handoff\n");
	return 1;
    }
    completed.changedCut = true;
    if (!policy.capacitySampleRequiresCeilingFreeHandoff(
	    completed, true, 4)) {
	std::fprintf(stderr,
	    "FAIL: changed certified allocation retained capacity ceiling\n");
	return 1;
    }
    completed.changedCut = false;
    completed.capacityTransactionPending = false;
    completed.changedCut = true;
    if (!policy.capacitySampleRequiresCeilingFreeHandoff(
	    completed, true, 4)) {
	std::fprintf(stderr,
	    "FAIL: changed capacity allocation did not yield to handoff\n");
	return 1;
    }
    completed.changedCut = false;
    completed.capacityTransactionPending = true;
    completed.presentationLimitsReconciled = false;
    if (policy.capacitySampleRequiresCeilingFreeHandoff(
	    completed, true, 4)) {
	std::fprintf(stderr,
	    "FAIL: unproved capacity sample released handoff ceiling\n");
	return 1;
    }
    completed.presentationLimitsReconciled = true;
    if (policy.capacitySampleRequiresCeilingFreeHandoff(
	    completed, false, 4) ||
	policy.capacitySampleRequiresCeilingFreeHandoff(
	    completed, true, -1)) {
	std::fprintf(stderr,
	    "FAIL: absent capacity barrier yielded presentation order\n");
	return 1;
    }

    policy.reset();
    input.population = population;
    input.scaleChanged = true;
    input.viewEpoch.set(11);
    policy.capturePrior(0.5f, 8, 0.625f, 1.25f, population,
	BObolLodViewEpoch(10));
    restore = policy.beginQuiet(input);
    if (!restore.apply || !restore.restoredPriorStable ||
	restore.handoff != Policy::RestoreDecision::Handoff::ALLOCATION ||
	std::fabs(restore.targetPixelError - 0.5f) > 0.0001f ||
	restore.progressiveCeiling != 8 ||
	std::fabs(restore.pointProxyPixelThreshold - 1.25f) > 0.0001f ||
	policy.handoffPresentationPending() || policy.handoffCostFloor() != 0) {
	std::fprintf(stderr,
	    "FAIL: quiet zoom did not restore stable controls for reallocation\n");
	return 1;
    }

    policy.reset();
    input.population = population;
    input.scaleChanged = true;
    input.viewEpoch.set(11);
    input.currentProgressiveCeiling = 6;
    input.currentPointProxyPixelThreshold = 2.0f;

    /* The current-scale frame supplies the safe renderer presentation and
     * cost floor.  It must not replace the stable semantic fidelity controls
     * captured at gesture start. */
    policy.capturePrior(0.25f, 18, 0.75f, 1.0f, population,
	BObolLodViewEpoch(10));
    policy.noteStableQuality(0.75f, 7, 0.375f, 1.5f, population,
	BObolLodViewEpoch(11), 7000);
    restore = policy.beginQuiet(input);
    if (!restore.apply || !restore.restoredProvenQuality ||
	std::fabs(restore.targetPixelError - 0.25f) > 0.0001f ||
	restore.handoff != Policy::RestoreDecision::Handoff::ALLOCATION ||
	restore.progressiveCeiling != 7 ||
	std::fabs(restore.progressiveNextFraction - 0.375f) > 0.0001f ||
	std::fabs(restore.pointProxyPixelThreshold - 1.0f) > 0.0001f ||
	restore.provenRenderCostFloor != 7000 ||
	policy.handoffCostFloor() != 7000 ||
	policy.handoffReconciliationBudget() != 7000 ||
	policy.allocationReconciliationBudget(false) != 7000 ||
	policy.allocationReconciliationBudget(true) != 0) {
	std::fprintf(stderr, "FAIL: stable zoom-quality proof restore\n");
	return 1;
    }
    completed = Policy::CompletedPassInputs();
    completed.completed = true;
    completed.retainedAllocationCompleted = true;
    completed.retainedAllocationCertified = true;
    completed.presentationLimitsReconciled = true;
    completed.populationQuiescent = true;
    completed.submissionPending = true;
    if (policy.completePass(completed).finishHandoff ||
	policy.handoffCostFloor() != 7000) {
	std::fprintf(stderr, "FAIL: pending submission released quality proof\n");
	return 1;
    }
    completed.submissionPending = false;
    completion = policy.completePass(completed);
    if (!completion.finishHandoff || policy.handoffCostFloor() != 0 ||
	policy.handoffReconciliationBudget() != 0) {
	std::fprintf(stderr, "FAIL: quality proof handoff completion\n");
	return 1;
    }

    policy.reset();
    policy.noteStableQuality(0.75f, 7, 0.375f, 1.5f, population,
	BObolLodViewEpoch(11), 7000);
    input.viewEpoch.set(12);
    restore = policy.beginQuiet(input);
    if (!restore.apply || restore.restoredProvenQuality ||
	restore.handoff != Policy::RestoreDecision::Handoff::ALLOCATION ||
	policy.handoffCostFloor() != 0 ||
	std::fabs(restore.targetPixelError - 1.0f) > 0.0001f) {
	std::fprintf(stderr, "FAIL: stale view reused quality proof\n");
	return 1;
    }

    policy.reset();
    input.retainedMeshPayloads = false;
    restore = policy.beginQuiet(input);
    if (!restore.apply || !restore.clearPresentationLimits ||
	std::fabs(restore.targetPixelError - 1.0f) > 0.0001f ||
	restore.progressiveCeiling != -1 ||
	std::fabs(restore.progressiveNextFraction) > 0.0001f ||
	std::fabs(restore.pointProxyPixelThreshold - 1.0f) > 0.0001f ||
	policy.handoffPending()) {
	std::fprintf(stderr, "FAIL: payload-free presentation reset\n");
	return 1;
    }
    return 0;
}

static int
test_completed_pass_composed_lifecycle(void)
{
    using Policy = BObolLodPresentationPolicy;
    using Scheduler = BObolLodAvailabilityScheduler;
    constexpr uint64_t TraceCount = 512;
    constexpr size_t StepBound = 20;
    constexpr uint64_t RandomMultiplier = UINT64_C(6364136223846793005);
    constexpr uint64_t RandomIncrement = UINT64_C(1442695040888963407);
    constexpr unsigned int HandoffStateCount = 3;
    constexpr unsigned int PresentationHandoffState = 1;
    constexpr unsigned int AllocationHandoffState = 2;
    constexpr unsigned int CapacityPopulationCount = 2;
    constexpr size_t NonzeroPresentationCost = 1;

    for (uint64_t trace = 1; trace <= TraceCount; ++trace) {
	uint64_t randomState = trace;
	auto random = [&randomState]() {
	    randomState = randomState * RandomMultiplier + RandomIncrement;
	    return randomState;
	};
	auto randomBool = [&random]() { return (random() & 1u) != 0; };

	Policy policy;
	const unsigned int initialHandoff =
	    static_cast<unsigned int>(random() % HandoffStateCount);
	if (initialHandoff == PresentationHandoffState)
	    policy.armHandoff(true);
	else if (initialHandoff == AllocationHandoffState)
	    policy.armHandoff(false);

	BObolLodRetainedPassAnnotations annotations;
	if (randomBool())
	    annotations.noteRefinementPending();
	if (randomBool())
	    annotations.noteResidencyPending();
	if (randomBool())
	    annotations.noteCutAdvanced();
	if (randomBool())
	    annotations.noteBudgetBlocked();
	if (randomBool())
	    annotations.noteAdmittedWork();

	bool residencyDrainActive = randomBool();
	bool residentGrowthPending = randomBool();
	bool residentWorkPending = randomBool();
	bool coveragePending = randomBool();
	bool capacityTransactionPending = randomBool();
	bool capacityEligible = randomBool() || capacityTransactionPending ||
	    annotations.cutAdvanced() || annotations.budgetBlocked();
	bool structuralFrontierPending = randomBool();
    bool allocationCertified = randomBool();
    const bool capacitySamplePending = randomBool();
    const bool capacityPresentationPending =
	capacitySamplePending && randomBool();
	const int progressiveCeiling = randomBool() ? 0 : -1;
	unsigned int capacityPopulation =
	    static_cast<unsigned int>(random() % CapacityPopulationCount);
	const unsigned int capacityCandidate =
	    static_cast<unsigned int>(random() % CapacityPopulationCount);
	uint64_t capacityRevision = 0;
	bool terminal = false;

	for (size_t step = 0; step < StepBound && !terminal; ++step) {
	    const Scheduler::CompletedPassSuccessor availabilityOwner =
		Scheduler::completedPassSuccessor(true,
		    residencyDrainActive, residentGrowthPending,
		    annotations.residencyPending(), residentWorkPending,
		    false, false, annotations.budgetBlocked());
	    if (availabilityOwner ==
		    Scheduler::CompletedPassSuccessor::COMPLETE_RESIDENCY_DRAIN) {
		residencyDrainActive = false;
		annotations.reset();
		continue;
	    }
	    if (availabilityOwner ==
		    Scheduler::CompletedPassSuccessor::YIELD_TO_RESIDENT_GROWTH) {
		const Scheduler::ResidentGrowthSuccessor successor =
		    Scheduler::residentGrowthSuccessor(!coveragePending);
		residentGrowthPending = false;
		coveragePending = successor ==
		    Scheduler::ResidentGrowthSuccessor::RETRY_COVERAGE;
		annotations.reset();
		continue;
	    }
	    if (availabilityOwner ==
		    Scheduler::CompletedPassSuccessor::AWAIT_RESIDENT_RESULT) {
		residentWorkPending = false;
		annotations.clearResidencyPending();
		residentGrowthPending = true;
		continue;
	    }
	    if (availabilityOwner ==
		    Scheduler::CompletedPassSuccessor::SUBMIT_RESIDENT_REQUESTS) {
		annotations.clearResidencyPending();
		residentGrowthPending = true;
		annotations.reset();
		continue;
	    }
	    if (coveragePending) {
		coveragePending = false;
		annotations.reset();
		continue;
	    }
	    if (availabilityOwner ==
		    Scheduler::CompletedPassSuccessor::CALIBRATE_CAPACITY)
		capacityEligible = true;

	    Policy::CompletedPassInputs inputs;
	    inputs.completed = true;
	    inputs.capacityTransactionPending = capacityTransactionPending;
	    inputs.structuralFrontierPending = structuralFrontierPending;
	    inputs.changedCut = annotations.cutAdvanced();
	    inputs.retainedAllocationCompleted = allocationCertified;
	    inputs.retainedAllocationCertified = allocationCertified;
	    inputs.populationQuiescent = true;
	    inputs.presentationLimitsReconciled = allocationCertified;
	    inputs.retainedRefinementPending =
		annotations.refinementPending();
	    inputs.retainedRefinementBudgetBlocked =
		annotations.budgetBlocked();
	    const Policy::CompletedPassSelection selection =
		policy.completedPassSelection(inputs, capacityEligible,
		    capacitySamplePending, capacityPresentationPending,
		    progressiveCeiling);
	    if (selection.consumePassAnnotations) {
		inputs.capacityTransactionPending = false;
		inputs.changedCut = false;
	    }

	    const uint64_t previousRevision = capacityRevision;
	    const unsigned int previousPopulation = capacityPopulation;
	    unsigned int effectCount = 0;
	    bool successorPass = false;
	    switch (selection.owner) {
		case Policy::CompletedPassOwner::STRUCTURAL_FRONTIER:
		    ++effectCount;
		    successorPass = true;
		    structuralFrontierPending = false;
		    capacityTransactionPending = false;
		    annotations.reset();
		    break;
		case Policy::CompletedPassOwner::RETAINED_OBSERVATION:
		    ++effectCount;
		    annotations.clearRefinementPending();
		    break;
		case Policy::CompletedPassOwner::CAPACITY:
		    ++effectCount;
		    successorPass = true;
		    if (inputs.capacityTransactionPending || inputs.changedCut) {
			capacityTransactionPending = false;
			annotations.reset();
		    } else if (capacityPopulation != capacityCandidate) {
			capacityPopulation = capacityCandidate;
			++capacityRevision;
			annotations.reset();
		    } else {
			capacityEligible = false;
			annotations.reset();
		    }
		    break;
		case Policy::CompletedPassOwner::PRESENTATION_HANDOFF:
		    ++effectCount;
		    successorPass = true;
		    if (!policy.noteFramePresented(NonzeroPresentationCost,
			    NonzeroPresentationCost)) {
			std::fprintf(stderr,
			    "FAIL: completed-pass presentation transition "
			    "trace=%llu step=%zu\n",
			    static_cast<unsigned long long>(trace), step);
			return 1;
		    }
		    annotations.reset();
		    capacityTransactionPending = false;
		    break;
		case Policy::CompletedPassOwner::ALLOCATION_HANDOFF: {
		    ++effectCount;
		    successorPass = true;
		    const Policy::CompletedPassDecision decision =
			policy.completePass(inputs);
		    if (!allocationCertified) {
			if (!decision.requestRetainedAllocation) {
			    std::fprintf(stderr,
				"FAIL: completed-pass allocation transition "
				"trace=%llu step=%zu\n",
				static_cast<unsigned long long>(trace), step);
			    return 1;
			}
			allocationCertified = true;
		    } else if (!decision.finishHandoff) {
			std::fprintf(stderr,
			    "FAIL: completed-pass finish transition "
			    "trace=%llu step=%zu\n",
			    static_cast<unsigned long long>(trace), step);
			return 1;
		    }
		    annotations.reset();
		    capacityTransactionPending = false;
		    break;
		}
		case Policy::CompletedPassOwner::NONE:
		    terminal = true;
		    break;
	    }

	    if (effectCount > 1 ||
		    (capacityRevision != previousRevision) !=
			(capacityPopulation != previousPopulation) ||
		    (capacityRevision != previousRevision &&
		     selection.owner !=
			 Policy::CompletedPassOwner::CAPACITY) ||
		    (successorPass &&
		     (annotations.refinementPending() ||
		      annotations.residencyPending() ||
		      annotations.cutAdvanced() ||
		      annotations.budgetBlocked() ||
		      annotations.admittedWork()))) {
		std::fprintf(stderr,
		    "FAIL: completed-pass composed invariant trace=%llu "
		    "step=%zu owner=%u effects=%u\n",
		    static_cast<unsigned long long>(trace), step,
		    static_cast<unsigned int>(selection.owner), effectCount);
		return 1;
	    }
	}
	if (!terminal) {
	    std::fprintf(stderr,
		"FAIL: completed-pass trace did not terminate trace=%llu\n",
		static_cast<unsigned long long>(trace));
	    return 1;
	}
    }
    return 0;
}

static int
test_view_quality_history(void)
{
    using History = BObolLodViewQualityHistory;
    History history;
    History::RememberInputs input;
    input.view.haveCamera = 1;
    input.view.width = 1024;
    input.view.height = 768;
    input.view.size = 42.0;
    input.view.lodScale = 1.0;
    input.view.curveScale = 1.0;
    input.view.pointScale = 1.0;
    input.view.botThreshold = 1;
    input.view.viewVolumeMatrix[0] = 1.0f;
    input.view.viewVolumeMatrix[5] = 1.0f;
    input.view.viewVolumeMatrix[10] = 1.0f;
    input.view.viewVolumeMatrix[15] = 1.0f;
    input.domainRevision = 7;
    input.sceneAvailable = true;
    input.quality.targetPixelError = 0.5f;
    input.quality.progressiveCeiling = 12;
    input.quality.progressiveNextFraction = 0.625f;
    input.quality.pointProxyPixelThreshold = 1.5f;
    input.quality.presentedRenderCost = 24000;
    input.exactCompletedFrame = true;
    input.terminalPresentationComplete = true;
    input.producersSettled = true;
    input.presentationDeadlineMet = true;

    if (!history.remember(input) || history.size() != 1 ||
	history.rememberCount() != 1 || history.recallCount() != 0) {
	std::fprintf(stderr, "FAIL: exact view quality was not remembered\n");
	return 1;
    }

    History::RecallInputs recall;
    recall.view = input.view;
    recall.domainRevision = input.domainRevision;
    recall.sceneAvailable = true;
    History::Quality quality;
    if (!history.recall(recall, quality) ||
	std::fabs(quality.targetPixelError - 0.5f) > 0.0001f ||
	quality.progressiveCeiling != 12 ||
	std::fabs(quality.progressiveNextFraction - 0.625f) > 0.0001f ||
	std::fabs(quality.pointProxyPixelThreshold - 1.5f) > 0.0001f ||
	quality.presentedRenderCost != 24000 ||
	quality.provenRenderCostCapacity != 24000 ||
	history.recallCount() != 1) {
	std::fprintf(stderr, "FAIL: exact view quality recall\n");
	return 1;
    }

    History::RecallInputs mismatch = recall;
    mismatch.view.width++;
    if (history.recall(mismatch, quality)) {
	std::fprintf(stderr, "FAIL: viewport mismatch recalled quality\n");
	return 1;
    }
    mismatch = recall;
    mismatch.view.viewVolumeMatrix[12] = 1.0f;
    if (history.recall(mismatch, quality)) {
	std::fprintf(stderr, "FAIL: camera mismatch recalled quality\n");
	return 1;
    }
    mismatch = recall;
    mismatch.view.viewVolumeMatrix[1] = -0.0f;
    if (!history.recall(mismatch, quality)) {
	std::fprintf(stderr, "FAIL: signed-zero camera missed exact history\n");
	return 1;
    }
    mismatch = recall;
    mismatch.view.size = std::numeric_limits<double>::quiet_NaN();
    if (history.recall(mismatch, quality)) {
	std::fprintf(stderr, "FAIL: nonfinite camera recalled quality\n");
	return 1;
    }
    mismatch = recall;
    mismatch.domainRevision++;
    if (history.recall(mismatch, quality)) {
	std::fprintf(stderr, "FAIL: scene domain mismatch recalled quality\n");
	return 1;
    }
    mismatch = recall;
    mismatch.resourcePressure = true;
    if (history.recall(mismatch, quality)) {
	std::fprintf(stderr, "FAIL: resource pressure recalled quality\n");
	return 1;
    }

    /* A preliminary one-pixel frame on a revisit must not erase the richer
     * proof which that same exact view established earlier. */
    History::RememberInputs coarse = input;
    coarse.quality.targetPixelError = 1.0f;
    coarse.quality.presentedRenderCost = 6000;
    if (history.remember(coarse) ||
	!history.recall(recall, quality) ||
	std::fabs(quality.targetPixelError - 0.5f) > 0.0001f) {
	std::fprintf(stderr, "FAIL: preliminary frame erased rich history\n");
	return 1;
    }

    /* Rendering more total work proves throughput, not visual superiority.
     * A costly coarse distribution may raise the capacity seed without
     * replacing the exact view's better fidelity controls. */
    History::RememberInputs costlyCoarse = input;
    costlyCoarse.quality.targetPixelError = 1.0f;
    costlyCoarse.quality.progressiveCeiling = 2;
    costlyCoarse.quality.pointProxyPixelThreshold = 8.0f;
    costlyCoarse.quality.presentedRenderCost = 100000;
    if (history.remember(costlyCoarse) ||
	!history.recall(recall, quality) ||
	std::fabs(quality.targetPixelError - 0.5f) > 0.0001f ||
	quality.progressiveCeiling != 12 ||
	quality.presentedRenderCost != 24000 ||
	quality.provenRenderCostCapacity != 100000) {
	std::fprintf(stderr,
	    "FAIL: renderer capacity overwrote exact-view fidelity\n");
	return 1;
    }

    History::RememberInputs richer = input;
    richer.quality.targetPixelError = 0.25f;
    richer.quality.progressiveCeiling = -1;
    richer.quality.progressiveNextFraction = 0.0f;
    richer.quality.pointProxyPixelThreshold = 1.0f;
    richer.quality.presentedRenderCost = 48000;
    if (!history.remember(richer) ||
	!history.recall(recall, quality) ||
	std::fabs(quality.targetPixelError - 0.25f) > 0.0001f ||
	quality.progressiveCeiling != -1 ||
	quality.provenRenderCostCapacity != 100000) {
	std::fprintf(stderr, "FAIL: richer exact proof did not replace history\n");
	return 1;
    }

    /* The retained allocator's projected-error witness may enrich an
     * otherwise identical fidelity snapshot without using triangle count as
     * a visual metric. */
    History::RememberInputs measured = richer;
    measured.quality.maximumProjectedErrorPixels = 2.0;
    measured.quality.presentedRenderCost = 36000;
    if (!history.remember(measured) ||
	!history.recall(recall, quality) ||
	!std::isfinite(quality.maximumProjectedErrorPixels) ||
	std::fabs(quality.maximumProjectedErrorPixels - 2.0) > 1.0e-9 ||
	quality.progressiveCeiling != -1 ||
	quality.provenRenderCostCapacity != 100000) {
	std::fprintf(stderr,
	    "FAIL: measured screen-error proof was not retained\n");
	return 1;
    }

    /* Once both frames carry directly comparable pixel-error evidence, the
     * lower observed bound wins even when the global control vectors alone
     * are incomparable. */
    History::RememberInputs measuredBetter = measured;
    measuredBetter.quality.targetPixelError = 0.5f;
    measuredBetter.quality.progressiveCeiling = 10;
    measuredBetter.quality.progressiveNextFraction = 0.5f;
    measuredBetter.quality.pointProxyPixelThreshold = 1.25f;
    measuredBetter.quality.maximumProjectedErrorPixels = 0.75;
    measuredBetter.quality.presentedRenderCost = 40000;
    if (!history.remember(measuredBetter) ||
	!history.recall(recall, quality) ||
	std::fabs(quality.maximumProjectedErrorPixels - 0.75) > 1.0e-9 ||
	quality.progressiveCeiling != 10 ||
	std::fabs(quality.progressiveNextFraction - 0.5f) > 0.0001f ||
	quality.presentedRenderCost != 40000 ||
	quality.provenRenderCostCapacity != 100000) {
	std::fprintf(stderr,
	    "FAIL: lower measured projected error was not preferred\n");
	return 1;
    }

    History::RememberInputs invalid = input;
    invalid.exactCompletedFrame = false;
    if (history.remember(invalid)) {
	std::fprintf(stderr, "FAIL: inexact frame entered history\n");
	return 1;
    }
    invalid = input;
    invalid.resourcePressure = true;
    if (history.remember(invalid)) {
	std::fprintf(stderr, "FAIL: pressure-limited frame entered history\n");
	return 1;
    }
    invalid = input;
    invalid.terminalPresentationComplete = false;
    if (history.remember(invalid)) {
	std::fprintf(stderr,
	    "FAIL: structural/failure presentation entered history\n");
	return 1;
    }
    invalid = input;
    invalid.producersSettled = false;
    if (history.remember(invalid)) {
	std::fprintf(stderr,
	    "FAIL: frame with an unsettled geometry producer entered history\n");
	return 1;
    }
    invalid = input;
    invalid.presentationDeadlineMet = false;
    if (history.remember(invalid)) {
	std::fprintf(stderr, "FAIL: missed frame entered history\n");
	return 1;
    }
    invalid = input;
    invalid.quality.maximumProjectedErrorPixels =
	std::numeric_limits<double>::quiet_NaN();
    if (history.remember(invalid)) {
	std::fprintf(stderr, "FAIL: invalid pixel-error proof entered history\n");
	return 1;
    }
    invalid = input;
    invalid.quality.progressiveCeiling =
	BOBOL_MESH_LOD_CUT_COUNT_MAX;
    if (history.remember(invalid)) {
	std::fprintf(stderr, "FAIL: invalid PoP ceiling entered history\n");
	return 1;
    }

    /* Quality is multidimensional.  Without a comparable measured error
     * bound, neither a lower PoP ceiling nor a larger point-proxy threshold
     * may replace an established fidelity proof. */
    History::RememberInputs lowerCost = input;
    lowerCost.quality.progressiveCeiling = 3;
    lowerCost.quality.pointProxyPixelThreshold = 8.0f;
    lowerCost.quality.presentedRenderCost = 12000;
    if (history.remember(lowerCost) || !history.recall(recall, quality) ||
	quality.presentedRenderCost != 40000 ||
	quality.progressiveCeiling != 10 ||
	std::fabs(quality.pointProxyPixelThreshold - 1.25f) > 0.0001f ||
	quality.provenRenderCostCapacity != 100000) {
	std::fprintf(stderr, "FAIL: lower-cost frame erased richer proof\n");
	return 1;
    }

    /* Fixed-size move-to-front behavior must be deterministic: touching the
     * oldest entry protects it while the next insertion evicts the true LRU. */
    history.reset();
    for (size_t i = 0; i < History::capacity(); ++i) {
	History::RememberInputs entry = input;
	entry.view.width = static_cast<int32_t>(100 + i);
	if (!history.remember(entry)) {
	    std::fprintf(stderr, "FAIL: view history fill\n");
	    return 1;
	}
    }
    History::RecallInputs oldest = recall;
    oldest.view.width = 100;
    if (!history.recall(oldest, quality)) {
	std::fprintf(stderr, "FAIL: view history LRU promotion\n");
	return 1;
    }
    History::RememberInputs next = input;
    next.view.width = 999;
    if (!history.remember(next) ||
	history.size() != History::capacity()) {
	std::fprintf(stderr, "FAIL: bounded view history insertion\n");
	return 1;
    }
    History::RecallInputs evicted = recall;
    evicted.view.width = 101;
    if (history.recall(evicted, quality) ||
	!history.recall(oldest, quality)) {
	std::fprintf(stderr, "FAIL: view history LRU eviction\n");
	return 1;
    }

    history.reset();
    if (history.size() != 0 || history.rememberCount() != 0 ||
	history.recallCount() != 0 || history.recall(recall, quality)) {
	std::fprintf(stderr, "FAIL: view history reset\n");
	return 1;
    }
    return 0;
}

static int
test_availability_and_publication(void)
{
    using Producer = BObolLodProducerPolicy;
    if (Producer::canProduceGeometry(false, false, false, false) ||
	!Producer::canProduceGeometry(true, false, false, false) ||
	Producer::canProduceGeometry(true, true, false, false) ||
	!Producer::canProduceGeometry(true, true, true, false) ||
	!Producer::canProduceGeometry(true, true, false, true) ||
	Producer::ownsFutureFrame(true, true, false, false, false) ||
	!Producer::ownsFutureFrame(true, true, false, false, true) ||
	Producer::ownsFutureFrame(false, false, false, false, false)) {
	std::fprintf(stderr,
	    "FAIL: geometry producer continuation/paused-frame contract\n");
	return 1;
    }

    BObolLodAvailabilityLedger availability;
    if (availability.deferInventoryDelta(
	    true, true, false, false, false, 100) ||
	availability.inventoryFirstPendingMicroseconds() != 0) {
	std::fprintf(stderr,
	    "FAIL: first source inventory must submit immediately\n");
	return 1;
    }
    if (!availability.deferInventoryDelta(
	    true, true, false, true, false, 100) ||
	availability.inventoryFirstPendingMicroseconds() != 100 ||
	!availability.deferInventoryDelta(
	    true, true, false, true, false, 250099) ||
	availability.deferInventoryDelta(
	    true, true, false, true, false, 250100)) {
	std::fprintf(stderr,
	    "FAIL: stable source-inventory delta coalescing bound\n");
	return 1;
    }
    availability.commitInventoryDelta();
    if (!availability.deferInventoryDelta(
	    true, true, false, true, true, 500000) ||
	availability.deferInventoryDelta(
	    true, true, false, true, true, 600000)) {
	std::fprintf(stderr,
	    "FAIL: interactive source-inventory delta coalescing bound\n");
	return 1;
    }
    if (availability.deferInventoryDelta(
	    true, false, false, true, false, 700000) ||
	availability.deferInventoryDelta(
	    true, true, true, true, false, 700000) ||
	availability.deferInventoryDelta(
	    false, true, false, true, false, 700000) ||
	availability.inventoryFirstPendingMicroseconds() != 0) {
	std::fprintf(stderr,
	    "FAIL: source-inventory semantic/final edge bypass\n");
	return 1;
    }

    availability.setProviderPendingCount(7);
    availability.noteResultsReady(1234);
    availability.noteResultsReady(5678);
    if (availability.providerPendingCount() != 7 ||
	!availability.resultsPending() ||
	availability.firstResultReadyMicroseconds() != 1234) {
	std::fprintf(stderr, "FAIL: availability producer/result ledger\n");
	return 1;
    }
    availability.resetResultQueue();
    if (availability.resultsPending() ||
	availability.firstResultReadyMicroseconds() != 0) {
	std::fprintf(stderr, "FAIL: availability result reset\n");
	return 1;
    }

    using Policy = BObolLodPresentationTransaction;
    if (Policy::presentationIntervalMicroseconds(0, false) != 50000 ||
	Policy::presentationIntervalMicroseconds(20000000, false) != 50000 ||
	Policy::presentationIntervalMicroseconds(200000000, false) != 250000 ||
	Policy::presentationIntervalMicroseconds(0, true) != 33333 ||
	Policy::presentationIntervalMicroseconds(60000000, true) != 100000 ||
	Policy::refinementCooldownMicroseconds(20000000, false) != 0 ||
	Policy::refinementCooldownMicroseconds(40000000, false) != 50000 ||
	Policy::refinementCooldownMicroseconds(250000000, false) != 250000 ||
	Policy::refinementCooldownMicroseconds(3000000000ULL, false) !=
	    2000000 ||
	Policy::refinementCooldownMicroseconds(250000000, true) != 0) {
	std::fprintf(stderr, "FAIL: adaptive presentation cadence/cooldown\n");
	return 1;
    }
    Policy policy;
    Policy::Inputs input;
    input.nowMicroseconds = 100;
    Policy::Decision decision = policy.decide(input);
    if (decision.keepPumpAlive || decision.requestFrame ||
	policy.publicationPending() ||
	policy.publicationFramePending() ||
	policy.publicationAwaitingFrameRequest() ||
	policy.unpresentedCount() != 0) {
	std::fprintf(stderr, "FAIL: initial publication policy\n");
	return 1;
    }

    policy.noteApplied(2, 100, 3, 7);
    policy.noteApplied(3, 110, 3, 7);
    input.nowMicroseconds = 50099;
    decision = policy.decide(input);
    if (!decision.keepPumpAlive || decision.requestFrame ||
	!policy.publicationPending() || policy.publicationFramePending() ||
	!policy.publicationAwaitingFrameRequest() ||
	policy.unpresentedCount() != 5 ||
	policy.firstUnpresentedMicroseconds() != 100) {
	std::fprintf(stderr, "FAIL: bounded publication batch\n");
	return 1;
    }
    input.nowMicroseconds = 50100;
    decision = policy.decide(input);
    if (decision.keepPumpAlive || !decision.requestFrame ||
	!policy.publicationFramePending() ||
	policy.publicationAwaitingFrameRequest()) {
	std::fprintf(stderr, "FAIL: stable publication deadline\n");
	return 1;
    }
    /* A batch timer which has not yet requested its frame can own a future
     * point-calibration presentation.  The requested publication frame is
     * that presentation obligation, not another producer which may justify
     * dropping its render edge. */
    BObolLodAdmissionPlanner::PointCalibrationProducerInputs
	publicationProducerInputs;
    publicationProducerInputs.publicationAwaitingFrameRequest =
	policy.publicationAwaitingFrameRequest();
    if (BObolLodAdmissionPlanner::pointProducerOwnsCalibrationFrame(
	    publicationProducerInputs)) {
	std::fprintf(stderr,
	    "FAIL: requested publication frame impersonates a future producer\n");
	return 1;
    }
    decision = policy.decide(input);
    if (decision.keepPumpAlive || decision.requestFrame ||
	!policy.publicationFramePending()) {
	std::fprintf(stderr, "FAIL: duplicate publication frame\n");
	return 1;
    }
    (void)policy.complete(1, 3, 7, true);
    if (policy.publicationPending() ||
	policy.publicationFramePending() ||
	policy.publicationAwaitingFrameRequest() ||
	policy.unpresentedCount() != 0 ||
	policy.firstUnpresentedMicroseconds() != 0) {
	std::fprintf(stderr, "FAIL: publication completion\n");
	return 1;
    }

    policy.noteApplied(1, 200, 3, 7);
    input = Policy::Inputs();
    input.nowMicroseconds = 200;
    input.firstUseful = true;
    decision = policy.decide(input);
    if (!decision.requestFrame || !policy.publicationFramePending()) {
	std::fprintf(stderr, "FAIL: first-useful publication\n");
	return 1;
    }
    policy.reset();
    policy.noteApplied(1, 300, 3, 7);
    input = Policy::Inputs();
    input.nowMicroseconds = 300;
    input.streamIdle = true;
    decision = policy.decide(input);
    if (!decision.requestFrame || !policy.publicationFramePending()) {
	std::fprintf(stderr, "FAIL: idle-stream publication\n");
	return 1;
    }

    policy.reset();
    policy.noteApplied(1, 1000, 3, 7);
    input = Policy::Inputs();
    input.interactive = true;
    input.nowMicroseconds = 34332;
    decision = policy.decide(input);
    if (decision.requestFrame) {
	std::fprintf(stderr, "FAIL: early interactive publication\n");
	return 1;
    }
    input.nowMicroseconds = 34333;
    decision = policy.decide(input);
    if (!decision.requestFrame) {
	std::fprintf(stderr, "FAIL: interactive publication deadline\n");
	return 1;
    }

    policy.reset();
    policy.noteApplied(1, 1000, 3, 7);
    input = Policy::Inputs();
    input.observedRenderNanoseconds = 100000000ULL;
    input.nowMicroseconds = 200999;
    decision = policy.decide(input);
    if (decision.requestFrame) {
	std::fprintf(stderr, "FAIL: adaptive publication fired early\n");
	return 1;
    }
    input.nowMicroseconds = 201000;
    decision = policy.decide(input);
    if (!decision.requestFrame) {
	std::fprintf(stderr, "FAIL: adaptive publication deadline\n");
	return 1;
    }

    policy.reset();
    policy.noteApplied(SIZE_MAX - 1, 500, 3, 7);
    policy.noteApplied(10, 600, 3, 7);
    input = Policy::Inputs();
    input.nowMicroseconds = 499;
    decision = policy.decide(input);
    if (!decision.requestFrame || policy.unpresentedCount() != SIZE_MAX) {
	std::fprintf(stderr, "FAIL: publication saturation/backward clock\n");
	return 1;
    }
    policy.reset();
    if (policy.publicationPending() ||
	policy.publicationFramePending()) {
	std::fprintf(stderr, "FAIL: publication reset\n");
	return 1;
    }
    return 0;
}

static int
test_availability_growth(void)
{
    BObolLodAvailabilityLedger availability;
    if (!BObolLodAvailabilityScheduler::needsIndependentResidentGrowth(
	    true, false) ||
	BObolLodAvailabilityScheduler::needsIndependentResidentGrowth(
	    true, true) ||
	BObolLodAvailabilityScheduler::needsIndependentResidentGrowth(
	    false, false) ||
	BObolLodAvailabilityScheduler::needsIndependentResidentGrowth(
	    false, true)) {
	std::fprintf(stderr,
	    "FAIL: capacity candidate did not retain suffix-result ownership\n");
	return 1;
    }
    if (!BObolLodAvailabilityScheduler::
	    retainedPublicationRequiresDemandReplay(true, true, false) ||
	BObolLodAvailabilityScheduler::
	    retainedPublicationRequiresDemandReplay(false, true, false) ||
	BObolLodAvailabilityScheduler::
	    retainedPublicationRequiresDemandReplay(true, false, false) ||
	BObolLodAvailabilityScheduler::
	    retainedPublicationRequiresDemandReplay(true, true, true)) {
	std::fprintf(stderr,
	    "FAIL: retained publication lost its quality successor\n");
	return 1;
    }
    if (!BObolLodAvailabilityScheduler::allocationPopulationSettled(
	    true, true, true, false, false) ||
	BObolLodAvailabilityScheduler::allocationPopulationSettled(
	    false, true, true, false, false) ||
	BObolLodAvailabilityScheduler::allocationPopulationSettled(
	    true, false, true, false, false) ||
	BObolLodAvailabilityScheduler::allocationPopulationSettled(
	    true, true, false, false, false) ||
	BObolLodAvailabilityScheduler::allocationPopulationSettled(
	    true, true, true, true, false) ||
	BObolLodAvailabilityScheduler::allocationPopulationSettled(
	    true, true, true, false, true)) {
	std::fprintf(stderr,
	    "FAIL: resident-growth allocation population boundary\n");
	return 1;
    }
    if (!BObolLodAvailabilityScheduler::structuralRepairMayOwn(false, false) ||
	BObolLodAvailabilityScheduler::structuralRepairMayOwn(true, false) ||
	BObolLodAvailabilityScheduler::structuralRepairMayOwn(false, true) ||
	BObolLodAvailabilityScheduler::structuralRepairMayOwn(true, true)) {
	std::fprintf(stderr,
	    "FAIL: resident growth did not exclude structural-repair ownership\n");
	return 1;
    }
    using AvailabilitySuccessor =
	BObolLodAvailabilityScheduler::CompletedPassSuccessor;
    if (BObolLodAvailabilityScheduler::residentGrowthOwnsSuccessor(
	    AvailabilitySuccessor::NONE) ||
	BObolLodAvailabilityScheduler::residentGrowthOwnsSuccessor(
	    AvailabilitySuccessor::CALIBRATE_CAPACITY) ||
	!BObolLodAvailabilityScheduler::residentGrowthOwnsSuccessor(
	    AvailabilitySuccessor::COMPLETE_RESIDENCY_DRAIN) ||
	!BObolLodAvailabilityScheduler::residentGrowthOwnsSuccessor(
	    AvailabilitySuccessor::YIELD_TO_RESIDENT_GROWTH)) {
	std::fprintf(stderr,
	    "FAIL: resident-growth successor ownership mismatch\n");
	return 1;
    }
    if (BObolLodAvailabilityScheduler::residentRefinementOwnsSuccessor(
	    AvailabilitySuccessor::NONE) ||
	BObolLodAvailabilityScheduler::residentRefinementOwnsSuccessor(
	    AvailabilitySuccessor::CALIBRATE_CAPACITY) ||
	!BObolLodAvailabilityScheduler::residentRefinementOwnsSuccessor(
	    AvailabilitySuccessor::AWAIT_RESIDENT_RESULT) ||
	!BObolLodAvailabilityScheduler::residentRefinementOwnsSuccessor(
	    AvailabilitySuccessor::SUBMIT_RESIDENT_REQUESTS)) {
	std::fprintf(stderr,
	    "FAIL: resident-refinement successor ownership mismatch\n");
	return 1;
    }
    using ResidentGrowthSuccessor =
	BObolLodAvailabilityScheduler::ResidentGrowthSuccessor;
    if (BObolLodAvailabilityScheduler::residentGrowthSuccessor(false) !=
	    ResidentGrowthSuccessor::RETRY_COVERAGE ||
	BObolLodAvailabilityScheduler::residentGrowthSuccessor(true) !=
	    ResidentGrowthSuccessor::REALLOCATE) {
	std::fprintf(stderr,
	    "FAIL: resident-growth coverage handoff mismatch\n");
	return 1;
    }
    if (!BObolLodAvailabilityScheduler::
	    presentationCutDowngradeAllowed(false, false, true, false, true) ||
	!BObolLodAvailabilityScheduler::
	    presentationCutDowngradeAllowed(false, false, false, true, false) ||
	BObolLodAvailabilityScheduler::
	    presentationCutDowngradeAllowed(false, false, true, true, false) ||
	BObolLodAvailabilityScheduler::
	    presentationCutDowngradeAllowed(true, false, false, true, true) ||
	BObolLodAvailabilityScheduler::
	    presentationCutDowngradeAllowed(false, true, false, true, true) ||
	BObolLodAvailabilityScheduler::
	    presentationCutDowngradeAllowed(false, false, false, false, false)) {
	std::fprintf(stderr,
	    "FAIL: resident-growth presentation-cut ownership boundary\n");
	return 1;
    }
    if (!BObolLodAvailabilityScheduler::
	    residentPrefetchPastAllocationAllowed(true, false) ||
	!BObolLodAvailabilityScheduler::
	    residentPrefetchPastAllocationAllowed(false, true) ||
	!BObolLodAvailabilityScheduler::
	    residentPrefetchPastAllocationAllowed(true, true) ||
	BObolLodAvailabilityScheduler::
	    residentPrefetchPastAllocationAllowed(false, false)) {
	std::fprintf(stderr,
	    "FAIL: resident prefetch retained a stale allocation boundary\n");
	return 1;
    }
    if (!BObolLodAvailabilityScheduler::canRetainPresentation(
	    true, true, true, true, true) ||
	BObolLodAvailabilityScheduler::canRetainPresentation(
	    false, true, true, true, true) ||
	BObolLodAvailabilityScheduler::canRetainPresentation(
	    true, false, true, true, true) ||
	BObolLodAvailabilityScheduler::canRetainPresentation(
	    true, true, false, true, true) ||
	BObolLodAvailabilityScheduler::canRetainPresentation(
	    true, true, true, false, true) ||
	BObolLodAvailabilityScheduler::canRetainPresentation(
	    true, true, true, true, false)) {
	std::fprintf(stderr,
	    "FAIL: resident-growth presentation-retention contract\n");
	return 1;
    }
    if (availability.residentGrowthPending() ||
	availability.residencyDrainRequired() ||
	availability.residencyDrainActive() ||
	availability.beginResidencyDrainIfReady(true, true, true) ||
	availability.consumeResidentGrowthIfReady(true, true, true)) {
	std::fprintf(stderr, "FAIL: initial resident-growth ledger\n");
	return 1;
    }

    availability.noteRicherPrefixAvailable();
    availability.noteRicherPrefixAvailable();

    if (!availability.residentGrowthPending() ||
	!availability.residencyDrainRequired() ||
	availability.consumeResidentGrowthIfReady(false, true, true) ||
	availability.consumeResidentGrowthIfReady(true, false, true) ||
	availability.consumeResidentGrowthIfReady(true, true, false) ||
	!availability.residentGrowthPending()) {
	std::fprintf(stderr,
	    "FAIL: resident growth did not coalesce behind readiness\n");
	return 1;
    }
    if (!availability.transferResidentGrowthToCapacityCandidate() ||
	availability.residentGrowthPending() ||
	availability.transferResidentGrowthToCapacityCandidate()) {
	std::fprintf(stderr,
	    "FAIL: capacity candidate did not subsume pending resident growth\n");
	return 1;
    }

    availability.noteRicherPrefixAvailable();
    if (availability.beginResidencyDrainIfReady(false, true, true) ||
	availability.beginResidencyDrainIfReady(true, false, true) ||
	availability.beginResidencyDrainIfReady(true, true, false) ||
	!availability.beginResidencyDrainIfReady(true, true, true) ||
	availability.residencyDrainRequired() ||
	!availability.residencyDrainActive() ||
	availability.beginResidencyDrainIfReady(true, true, true)) {
	std::fprintf(stderr,
	    "FAIL: resident growth residency drain readiness\n");
	return 1;
    }
    if (!availability.completeResidencyDrain() ||
	availability.residencyDrainActive() ||
	!availability.consumeResidentGrowthIfReady(true, true, true) ||
	availability.residentGrowthPending() ||
	availability.consumeResidentGrowthIfReady(true, true, true)) {
	std::fprintf(stderr,
	    "FAIL: resident growth allocation edge was not one-shot\n");
	return 1;
    }

    availability.noteRicherPrefixAvailable();
    if (!availability.beginResidencyDrainIfReady(true, true, true)) {
	std::fprintf(stderr, "FAIL: resident-growth second drain\n");
	return 1;
    }
    availability.noteRicherPrefixAvailable();
    if (!availability.residencyDrainRequired() ||
	availability.consumeResidentGrowthIfReady(true, true, true) ||
	!availability.completeResidencyDrain() ||
	!availability.beginResidencyDrainIfReady(true, true, true) ||
	!availability.completeResidencyDrain() ||
	!availability.consumeResidentGrowthIfReady(true, true, true)) {
	std::fprintf(stderr,
	    "FAIL: resident growth during drain was not coalesced\n");
	return 1;
    }
    availability.noteRicherPrefixAvailable();
    if (!availability.beginResidencyDrainIfReady(true, true, true)) {
	std::fprintf(stderr, "FAIL: resident-growth interrupt setup\n");
	return 1;
    }
    availability.interruptResidencyDrain();
    if (availability.residencyDrainActive() ||
	!availability.residencyDrainRequired() ||
	!availability.beginResidencyDrainIfReady(true, true, true) ||
	!availability.completeResidencyDrain() ||
	!availability.consumeResidentGrowthIfReady(true, true, true)) {
	std::fprintf(stderr, "FAIL: resident-growth drain interruption\n");
	return 1;
    }
    availability.noteRicherPrefixAvailable();
    availability.resetResidentGrowth();
    if (availability.residentGrowthPending() ||
	availability.residencyDrainRequired() ||
	availability.residencyDrainActive()) {
	std::fprintf(stderr, "FAIL: resident-growth reset\n");
	return 1;
    }
    return 0;
}

static int
test_compaction_policy(void)
{
    using Policy = BObolLodCompactionPolicy;
    Policy policy;
    const auto readyInput = []() {
	Policy::Inputs input;
	input.automatic = true;
	input.coverageRequired = true;
	input.coverageComplete = true;
	input.serviceAvailable = true;
	input.residentDataAvailable = true;
	input.nowMicroseconds = 150;
	return input;
    };
    Policy::Decision decision = policy.decide(readyInput());
    Policy::Inputs input;
    if (decision.admission != Policy::Admission::INACTIVE ||
	decision.keepPumpAlive || decision.retiredRequest ||
	policy.pending() || policy.planning()) {
	std::fprintf(stderr, "FAIL: initial compaction policy\n");
	return 1;
    }

    policy.requestAfter(100, 50);
    if (!policy.pending() || policy.planning() ||
	policy.deadlineMicroseconds() != 150) {
	std::fprintf(stderr, "FAIL: delayed compaction request\n");
	return 1;
    }

    input = readyInput();
    input.residentDataAvailable = false;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::INACTIVE ||
	decision.keepPumpAlive || !decision.retiredRequest ||
	policy.pending() || policy.planning() ||
	policy.deadlineMicroseconds() != 0) {
	std::fprintf(stderr, "FAIL: empty resident store retained compaction\n");
	return 1;
    }

    policy.requestAfter(100, 50);
    policy.requestAfter(125, 500);
    if (policy.deadlineMicroseconds() != 150) {
	std::fprintf(stderr, "FAIL: repeated compaction request slid deadline\n");
	return 1;
    }
    input = readyInput();
    input.coverageComplete = false;
    input.coverageProgressPending = true;
    input.nowMicroseconds = 200;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::DEFER ||
	!decision.keepPumpAlive || decision.retiredRequest ||
	!policy.pending()) {
	std::fprintf(stderr, "FAIL: active coverage not deferred\n");
	return 1;
    }
    input.coverageProgressPending = false;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::INACTIVE ||
	decision.keepPumpAlive || !decision.retiredRequest ||
	policy.pending() || policy.deadlineMicroseconds() != 0) {
	std::fprintf(stderr, "FAIL: unwitnessed coverage request retained\n");
	return 1;
    }

    /* With mesh LoD disabled there is no mesh-coverage prerequisite.  An
	 * empty demand pass must be allowed to retire service-initialization
	 * state instead of becoming an impossible background latch. */
    policy.requestAfter(100, 50);
    input = readyInput();
    input.coverageRequired = false;
    input.coverageComplete = false;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::PLAN ||
	decision.keepPumpAlive || decision.retiredRequest) {
	std::fprintf(stderr, "FAIL: optional coverage blocked compaction\n");
	return 1;
    }
    input = readyInput();
    input.nowMicroseconds = 149;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::DEFER ||
	!decision.keepPumpAlive) {
	std::fprintf(stderr, "FAIL: compaction bypassed quiet deadline\n");
	return 1;
    }
    input = readyInput();
    input.realizationPending = true;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::DEFER ||
	!decision.keepPumpAlive) {
	std::fprintf(stderr, "FAIL: compaction bypassed realization\n");
	return 1;
    }
    input = readyInput();
    input.submissionPending = true;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::DEFER ||
	!decision.keepPumpAlive) {
	std::fprintf(stderr, "FAIL: compaction bypassed submission\n");
	return 1;
    }
    input = readyInput();
    input.settlingPending = true;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::DEFER ||
	!decision.keepPumpAlive) {
	std::fprintf(stderr, "FAIL: compaction bypassed visual settling\n");
	return 1;
    }
    input = readyInput();
    input.automatic = false;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::INACTIVE ||
	decision.keepPumpAlive) {
	std::fprintf(stderr, "FAIL: disabled compaction kept pump alive\n");
	return 1;
    }
    input = readyInput();
    input.interactive = true;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::INACTIVE ||
	decision.keepPumpAlive) {
	std::fprintf(stderr, "FAIL: interactive compaction admitted\n");
	return 1;
    }
    input = readyInput();
    input.serviceAvailable = false;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::INACTIVE ||
	decision.keepPumpAlive) {
	std::fprintf(stderr, "FAIL: service-free compaction admitted\n");
	return 1;
    }
    decision = policy.decide(readyInput());
    if (decision.admission != Policy::Admission::PLAN ||
	decision.keepPumpAlive || policy.continues(17)) {
	std::fprintf(stderr, "FAIL: ready compaction not admitted\n");
	return 1;
    }

    policy.finishPlanning(false, 17, 151);
    if (!policy.pending() || !policy.planning() ||
	!policy.continues(17) || policy.continues(18) ||
	policy.demandRevision() != 17 ||
	policy.deadlineMicroseconds() != 151) {
	std::fprintf(stderr, "FAIL: partial compaction continuation\n");
	return 1;
    }
    policy.finishPlanning(true, 17, 152);
    if (policy.pending() || policy.planning() ||
	policy.demandRevision() != 0 ||
	policy.deadlineMicroseconds() != 0) {
	std::fprintf(stderr, "FAIL: completed compaction retained state\n");
	return 1;
    }

    policy.requestAfter(200, 100);
    policy.requestImmediate(250);
    if (policy.deadlineMicroseconds() != 250) {
	std::fprintf(stderr, "FAIL: immediate compaction deadline\n");
	return 1;
    }
    policy.requestImmediate(275);
    if (policy.deadlineMicroseconds() != 250) {
	std::fprintf(stderr, "FAIL: immediate compaction delayed request\n");
	return 1;
    }
    policy.resetForServiceChange(false, 300, 50);
    if (policy.pending() || policy.planning() ||
	policy.deadlineMicroseconds() != 0) {
	std::fprintf(stderr, "FAIL: compaction service removal\n");
	return 1;
    }
    policy.resetForServiceChange(true, INT64_MAX - 10, 100);
    if (!policy.pending() || policy.planning() ||
	policy.deadlineMicroseconds() != INT64_MAX) {
	std::fprintf(stderr, "FAIL: compaction deadline saturation\n");
	return 1;
    }
    policy.retire();
    if (policy.pending() || policy.planning() ||
	policy.demandRevision() != 0 ||
	policy.deadlineMicroseconds() != 0) {
	std::fprintf(stderr, "FAIL: policy disable retained compaction state\n");
	return 1;
    }
    return 0;
}

/*
 * A deterministic composed-lifecycle model for the coordinator policies.
 *
 * The direct tests above prove each allocation-free value at its boundaries;
 * this runner proves that their witnesses compose.  It deliberately models
 * the owner-thread contracts rather than sleeping or launching workers: every
 * seed is an exact replay token and FakeClock is the only source of time.
 * Random perturbations include repeated partial coverage checkpoints,
 * provider failure during interaction, cancellation under memory pressure,
 * delayed publication, view invalidation, and bounded compaction.
 */
class BObolLodSeededSequence {
public:
    enum Witness : uint32_t {
	WITNESS_NONE = 0,
	WITNESS_SERVICE_TASK = 1u << 0,
	WITNESS_COVERAGE_CURSOR = 1u << 1,
	WITNESS_PUBLICATION_DEADLINE = 1u << 2,
	WITNESS_PRESENTATION_FRAME = 1u << 3,
	WITNESS_INTERACTION_INPUT = 1u << 4,
	WITNESS_QUIET_TIMER = 1u << 5,
	WITNESS_COMPACTION_DEADLINE = 1u << 6,
	WITNESS_COMPACTION_CURSOR = 1u << 7
    };

    explicit BObolLodSeededSequence(uint64_t seedValue) :
	seed(seedValue ? seedValue : 1), randomState(seedValue ? seedValue : 1)
    {
	coverage.reset();
    }

    uint32_t random(void)
    {
	/* xorshift64* has no platform library state and therefore preserves the
	 * seed as a stable failure reproducer. */
	uint64_t value = randomState;
	value ^= value >> 12;
	value ^= value << 25;
	value ^= value >> 27;
	randomState = value;
	return static_cast<uint32_t>(
	    (value * UINT64_C(2685821657736338717)) >> 32);
    }

    void perturb(unsigned int action)
    {
	switch (action % 13) {
	    case 0: this->startSource(); break;
	    case 1: this->coverageCheckpoint(); break;
	    case 2: this->completeProvider(); break;
	    case 3: this->startInteraction(); break;
	    case 4: this->endGesture(); break;
	    case 5: this->togglePressure(); break;
	    case 6: this->cancelGeneration(); break;
	    case 7: this->failProvider(); break;
	    case 8: this->applyResult(); break;
	    case 9: this->advanceClock(); break;
	    case 10: this->requestCompaction(); break;
	    case 11: this->invalidateView(); break;
	    case 12: this->changePolicy(); break;
	}
    }

    bool progressOne(void)
    {
	if (gestureActive) {
	    this->endGesture();
	    return true;
	}
	if (interactive || quietTimer) {
	    nowMicroseconds += 200000;
	    interactive = false;
	    quietTimer = false;
	    policyEpoch.advance();
	    return true;
	}
	if (servicePending) {
	    this->completeProvider();
	    return true;
	}
	if (coverageCursor) {
	    this->coverageCheckpoint();
	    return true;
	}
	if (publication.publicationFramePending()) {
            (void)publication.complete(1, viewEpoch.value(),
                policyEpoch.value(), true);
	    return true;
	}
	if (publication.publicationAwaitingFrameRequest()) {
	    nowMicroseconds += 300000;
	    BObolLodPresentationTransaction::Inputs inputs;
	    inputs.nowMicroseconds = nowMicroseconds;
	    inputs.interactive = interactive;
	    inputs.streamIdle = !servicePending;
	    (void)publication.decide(inputs);
	    return true;
	}
	if (compactionCursor) {
	    compaction.finishPlanning(true, residentDemandRevision,
		nowMicroseconds);
	    compactionCursor = false;
	    if (resources.recoveryPending())
		resources.markRecoveryHandled();
	    return true;
	}
	if (compaction.pending() || resources.recoveryPending()) {
	    BObolLodCompactionPolicy::Inputs inputs =
		this->compactionInputs();
	    if (nowMicroseconds < compaction.deadlineMicroseconds())
		nowMicroseconds = compaction.deadlineMicroseconds();
	    inputs.nowMicroseconds = nowMicroseconds;
	    const BObolLodCompactionPolicy::Decision decision =
		compaction.decide(inputs);
	    if (decision.admission ==
		    BObolLodCompactionPolicy::Admission::PLAN) {
		residentDemandRevision++;
		if (!residentDemandRevision)
		    residentDemandRevision++;
		compaction.finishPlanning(false, residentDemandRevision,
		    nowMicroseconds);
		compactionCursor = true;
	    } else if (decision.retiredRequest &&
		resources.recoveryPending()) {
		resources.markRecoveryHandled();
	    } else if (decision.admission ==
		BObolLodCompactionPolicy::Admission::DEFER) {
		/* Its prerequisite must be one of the earlier named witnesses. */
		return false;
	    }
	    return true;
	}
	return false;
    }

    uint32_t witnesses(void) const
    {
	uint32_t result = WITNESS_NONE;
	if (servicePending)
	    result |= WITNESS_SERVICE_TASK;
	if (coverageCursor)
	    result |= WITNESS_COVERAGE_CURSOR;
	if (publication.publicationAwaitingFrameRequest())
	    result |= WITNESS_PUBLICATION_DEADLINE;
	if (publication.publicationFramePending())
	    result |= WITNESS_PRESENTATION_FRAME;
	if (gestureActive)
	    result |= WITNESS_INTERACTION_INPUT;
	else if (interactive || quietTimer)
	    result |= WITNESS_QUIET_TIMER;
	if (compactionCursor)
	    result |= WITNESS_COMPACTION_CURSOR;
	else if (compaction.pending() || resources.recoveryPending())
	    result |= WITNESS_COMPACTION_DEADLINE;
	return result;
    }

    bool terminal(void)
    {
	const BObolLodConvergencePolicy::Decision decision =
	    this->convergenceDecision();
	return decision.terminal;
    }

    bool valid(void)
    {
	const BObolLodConvergencePolicy::Decision decision =
	    this->convergenceDecision();
	const uint32_t foregroundWitnesses = this->witnesses() &
	    ~(WITNESS_COMPACTION_DEADLINE | WITNESS_COMPACTION_CURSOR);
	if (!decision.terminal && foregroundWitnesses == WITNESS_NONE)
	    return false;
	if (decision.terminal && foregroundWitnesses != WITNESS_NONE)
	    return false;
	const bool successfulTerminal = decision.terminal &&
	    (decision.outcome == BObolLodConvergencePolicy::Outcome::READY ||
	     decision.outcome ==
		BObolLodConvergencePolicy::Outcome::CONSTRAINED);
	if (decision.viewReady != successfulTerminal)
	    return false;
	if (decision.terminalError !=
	    (decision.terminal && decision.outcome ==
		BObolLodConvergencePolicy::Outcome::ERROR))
	    return false;
	return true;
    }

    uint64_t replaySeed(void) const { return seed; }
    BObolLodConvergencePolicy::Phase phase(void)
    {
	return this->convergenceDecision().phase;
    }
    bool error(void)
    {
	return this->convergenceDecision().outcome ==
	    BObolLodConvergencePolicy::Outcome::ERROR;
    }
    uint32_t invariantMask(void) const { return 0; }

private:
    void startSource(void)
    {
	if (sourceActive || servicePending)
	    return;
	providerFailed = false;
	providerCancelled = false;
	sourceActive = true;
	coverage.setRequired(true);
	servicePending = true;
	visibleCount = 7 + random() % 57;
	coveredCount = 0;
	coverage.reset();
	coverage.activate(true);
	coverageCursor = true;
	viewEpoch.advance();
    }

    void coverageCheckpoint(void)
    {
	if (!coverageCursor)
	    return;
	const size_t remaining = visibleCount > coveredCount ?
	    visibleCount - coveredCount : 0;
	const size_t quantum = std::min<size_t>(
	    remaining, 1 + random() % 11);
	coverage.observe(quantum, quantum);
	coveredCount += quantum;
	if (quantum)
	    publication.noteApplied(quantum, nowMicroseconds,
		viewEpoch.value(), policyEpoch.value());
	if (coveredCount >= visibleCount) {
	    (void)coverage.completeIfReady(true, false);
	    coverageCursor = false;
	}
	this->publicationDecision(!coveredCount || coveredCount == quantum,
	    !servicePending && !coverageCursor);
    }

    void completeProvider(void)
    {
	if (!servicePending)
	    return;
	servicePending = false;
	this->publicationDecision(false, true);
    }

    void startInteraction(void)
    {
	if (gestureActive)
	    return;
	interactive = true;
	gestureActive = true;
	quietTimer = false;
    }

    void endGesture(void)
    {
	if (!gestureActive)
	    return;
	gestureActive = false;
	interactive = true;
	quietTimer = true;
    }

    void togglePressure(void)
    {
	const bool cpu = !resources.cpuPressure();
	const bool gpu = cpu && ((random() & 1u) != 0);
	const BObolLodResourceEvidence::Decision decision =
	    apply_resource_observation(resources, cpu, gpu, true);
	if (decision.queueRecovery)
	    compaction.requestImmediate(nowMicroseconds);
    }

    void cancelGeneration(void)
    {
	if (!sourceActive && !servicePending)
	    return;
	servicePending = false;
	sourceActive = false;
	coverage.setRequired(false);
	providerCancelled = true;
	providerFailed = false;
	coverageCursor = false;
	coverage.reset();
	publication.reset();
    }

    void failProvider(void)
    {
	if (!servicePending)
	    return;
	servicePending = false;
	sourceActive = false;
	coverage.setRequired(false);
	providerFailed = true;
	providerCancelled = false;
	coverageCursor = false;
	coverage.reset();
	publication.reset();
    }

    void applyResult(void)
    {
	if (!sourceActive || providerFailed || providerCancelled)
	    return;
	publication.noteApplied(1 + random() % 8, nowMicroseconds,
	    viewEpoch.value(), policyEpoch.value());
	this->publicationDecision(false, !servicePending);
    }

    void advanceClock(void)
    {
	nowMicroseconds += 1 + random() % 400000;
	if (quietTimer && nowMicroseconds > 150000) {
	    interactive = false;
	    quietTimer = false;
	    policyEpoch.advance();
	}
	this->publicationDecision(false, !servicePending);
    }

    void requestCompaction(void)
    {
	compaction.requestAfter(nowMicroseconds, 50000);
    }

    void invalidateView(void)
    {
	viewEpoch.advance();
	if (sourceActive) {
	    coverage.reset();
	    coverage.activate(true);
	    coverageCursor = true;
	    coveredCount = 0;
	}
    }

    void changePolicy(void)
    {
	policyEpoch.advance();
    }

    void publicationDecision(bool firstUseful, bool streamIdle)
    {
	BObolLodPresentationTransaction::Inputs inputs;
	inputs.nowMicroseconds = nowMicroseconds;
	inputs.interactive = interactive;
	inputs.firstUseful = firstUseful;
	inputs.streamIdle = streamIdle;
	(void)publication.decide(inputs);
    }

    BObolLodCompactionPolicy::Inputs compactionInputs(void) const
    {
	BObolLodCompactionPolicy::Inputs inputs;
	inputs.automatic = true;
	inputs.interactive = interactive;
	inputs.coverageRequired = sourceActive;
	inputs.coverageComplete =
	    !sourceActive || coverage.coverageComplete();
	inputs.coverageProgressPending = coverageCursor;
	inputs.settlingPending =
	    publication.publicationPending() || quietTimer;
	inputs.realizationPending = servicePending;
	inputs.submissionPending = coverageCursor;
	inputs.serviceAvailable = true;
	inputs.nowMicroseconds = nowMicroseconds;
	return inputs;
    }

    BObolLodConvergencePolicy::Decision convergenceDecision(void)
    {
	BObolLodConvergencePolicy::Inputs inputs;
	inputs.viewEpoch = viewEpoch;
	inputs.policyEpoch = policyEpoch;
	inputs.expectedLeafCount = sourceActive || providerFailed ?
	    visibleCount : 0;
	inputs.availableLeafCount = sourceActive ? coveredCount : 0;
	inputs.visibleTargetCount = sourceActive ? visibleCount : 0;
	inputs.activePayloadCount = sourceActive ? coveredCount : 0;
	inputs.satisfiedPayloadCount = inputs.activePayloadCount;
	inputs.pendingTasks = servicePending ? 1 : 0;
	inputs.submissionPending = coverageCursor;
	inputs.publicationPending = publication.publicationPending();
	inputs.interactive = interactive;
	inputs.compactionPending =
	    compaction.pending() || compactionCursor ||
	    resources.recoveryPending();
	inputs.structuralDiscovery = servicePending;
	inputs.failedSourceCount = providerFailed ? 1 : 0;
	inputs.gpuMemoryPressure = resources.gpuPressure();
	return convergence.evaluate(inputs);
    }
    const uint64_t seed;
    uint64_t randomState;
    int64_t nowMicroseconds = 1;
    BObolLodViewEpoch viewEpoch {1};
    BObolLodPolicyEpoch policyEpoch {1};
    uint64_t residentDemandRevision = 0;
    size_t visibleCount = 0;
    size_t coveredCount = 0;
    bool sourceActive = false;
    bool servicePending = false;
    bool coverageCursor = false;
    bool providerFailed = false;
    bool providerCancelled = false;
    bool interactive = false;
    bool gestureActive = false;
    bool quietTimer = false;
    bool compactionCursor = false;
    BObolLodResourceEvidence resources;
    BObolLodCoveragePolicy coverage;
    BObolLodPresentationTransaction publication;
    BObolLodCompactionPolicy compaction;
    BObolLodConvergencePolicy convergence;
};

static int
test_seeded_composed_lifecycle(void)
{
    static const size_t seedCount = 512;
    static const size_t perturbationsPerSeed = 96;
    static const size_t dischargeBound = 128;
    const auto discharge = [](BObolLodSeededSequence &sequence,
	    const char *scenario) {
	for (size_t step = 0;
		step < dischargeBound && !sequence.terminal(); ++step) {
	    if (sequence.progressOne() && sequence.valid())
		continue;
	    std::fprintf(stderr,
		"FAIL: composed lifecycle %s seed=%llu step=%zu "
		"phase=%u witnesses=0x%x invariant=0x%x\n",
		scenario,
		static_cast<unsigned long long>(sequence.replaySeed()),
		step, static_cast<unsigned int>(sequence.phase()),
		sequence.witnesses(), sequence.invariantMask());
	    return 1;
	}
	if (sequence.terminal())
	    return 0;
	std::fprintf(stderr,
	    "FAIL: composed lifecycle %s exceeded discharge bound seed=%llu "
	    "phase=%u witnesses=0x%x\n", scenario,
	    static_cast<unsigned long long>(sequence.replaySeed()),
	    static_cast<unsigned int>(sequence.phase()),
	    sequence.witnesses());
	return 1;
    };

    /* Keep the release-critical compositions explicit even though the seeded
     * corpus also reaches them.  The action numbers are part of this fake
     * service's stable replay vocabulary. */
    {
	BObolLodSeededSequence repeatedCheckpoints(UINT64_C(0x43504f494e54));
	repeatedCheckpoints.perturb(0);  /* source/service start */
	repeatedCheckpoints.perturb(1);
	repeatedCheckpoints.perturb(1);
	repeatedCheckpoints.perturb(1);
	repeatedCheckpoints.perturb(2);  /* provider stream idle */
	if (!repeatedCheckpoints.valid() ||
	    discharge(repeatedCheckpoints, "repeated-checkpoints"))
	    return 1;
    }
    {
	BObolLodSeededSequence failureDuringInput(UINT64_C(0x4641494c5552));
	failureDuringInput.perturb(0); /* source/service start */
	failureDuringInput.perturb(3); /* interaction start */
	failureDuringInput.perturb(7); /* terminal provider failure */
	if (!failureDuringInput.valid() || !failureDuringInput.error() ||
	    discharge(failureDuringInput, "failure-during-interaction")) {
	    std::fprintf(stderr,
		"FAIL: provider failure during interaction did not converge\n");
	    return 1;
	}
    }
    {
	BObolLodSeededSequence pressureCancellation(UINT64_C(0x505245535355));
	pressureCancellation.perturb(0); /* source/service start */
	pressureCancellation.perturb(5); /* memory pressure edge */
	pressureCancellation.perturb(6); /* generation cancellation */
	if (!pressureCancellation.valid() ||
	    discharge(pressureCancellation, "pressure-during-cancellation"))
	    return 1;
    }

    for (uint64_t seed = 1; seed <= seedCount; ++seed) {
	BObolLodSeededSequence sequence(seed);
	for (size_t step = 0; step < perturbationsPerSeed; ++step) {
	    sequence.perturb(sequence.random());
	    const size_t opportunistic = sequence.random() % 3;
	    for (size_t i = 0; i < opportunistic; ++i)
		(void)sequence.progressOne();
	    if (!sequence.valid()) {
		std::fprintf(stderr,
		    "FAIL: composed lifecycle seed=%llu step=%zu "
		    "phase=%u witnesses=0x%x invariant=0x%x\n",
		    static_cast<unsigned long long>(sequence.replaySeed()),
		    step, static_cast<unsigned int>(sequence.phase()),
		    sequence.witnesses(), sequence.invariantMask());
		return 1;
	    }
	}
	if (discharge(sequence, "seeded"))
	    return 1;
    }
    return 0;
}

static int
test_presentation_transaction(void)
{
    BObolLodPresentationTransaction transaction;
    if (transaction.barrierPending() ||
	transaction.publicationPending() ||
	transaction.requiredRenderSerial() != 0) {
	std::fprintf(stderr,
	    "FAIL: presentation transaction does not initialize empty\n");
	return 1;
    }

    if (!transaction.arm(
	    BObolLodPresentationTransaction::REASON_CUT_PRESENTATION,
	    10, 3, 7) || !transaction.barrierPending() ||
	transaction.requiredRenderSerial() != 11 ||
	transaction.viewEpoch() != 3 || transaction.policyEpoch() != 7) {
	std::fprintf(stderr,
	    "FAIL: frame barrier did not capture its transaction token\n");
	return 1;
    }
    const uint64_t sequence = transaction.sequence();
    transaction.noteApplied(2, 100, 3, 7);
    if (transaction.arm(
	    BObolLodPresentationTransaction::REASON_RESULT_PUBLICATION,
	    12, 3, 7) || transaction.requiredRenderSerial() != 11 ||
	transaction.sequence() != sequence ||
	!transaction.publicationPending() ||
	(transaction.reasons() &
	 (BObolLodPresentationTransaction::REASON_CUT_PRESENTATION |
	  BObolLodPresentationTransaction::REASON_RESULT_PUBLICATION)) == 0) {
	std::fprintf(stderr,
	    "FAIL: shared frame reasons moved an existing barrier forward\n");
	return 1;
    }

    BObolLodPresentationTransaction::Completion completion =
	transaction.complete(10, 3, 7, true);
    if (completion.retired || completion.stale ||
	!transaction.barrierPending() ||
	!transaction.publicationPending()) {
	std::fprintf(stderr,
	    "FAIL: an early frame retired its obligation\n");
	return 1;
    }
    completion = transaction.complete(11, 3, 7, false);
    if (completion.retired || completion.stale ||
	!transaction.barrierPending() ||
	!transaction.publicationPending()) {
	std::fprintf(stderr,
	    "FAIL: an incomplete frame retired its presentation obligation\n");
	return 1;
    }

    completion = transaction.complete(11, 3, 7, true);
    if (!completion.retired || completion.stale ||
	transaction.barrierPending() || transaction.publicationPending() ||
	(completion.reasons &
	 BObolLodPresentationTransaction::REASON_RESULT_PUBLICATION) == 0) {
	std::fprintf(stderr,
	    "FAIL: matching frame did not retire the presentation transaction\n");
	return 1;
    }


    transaction.noteApplied(1, 200, 4, 8);
    (void)transaction.arm(BObolLodPresentationTransaction::REASON_HANDOFF,
	20, 4, 8);
    completion = transaction.complete(21, 5, 8, false);
    if (completion.retired || !completion.stale ||
	transaction.barrierPending() || transaction.publicationPending()) {
	std::fprintf(stderr,
	    "FAIL: stale camera epoch acknowledged a presentation transaction\n");
	return 1;
    }


    (void)transaction.arm(
	BObolLodPresentationTransaction::REASON_CALIBRATION,
	UINT64_MAX - 1, 9, 11);
    if (transaction.requiredRenderSerial() != UINT64_MAX ||
	!transaction.complete(UINT64_MAX, 9, 11, true).retired) {
	std::fprintf(stderr,
	    "FAIL: presentation serial final successor is not deterministic\n");
	return 1;
    }

    transaction.noteApplied(1, 300, 12, 13);
    completion = transaction.complete(1, 12, 13, true);
    if (completion.retired || completion.stale ||
	transaction.publicationPending()) {
	std::fprintf(stderr,
	    "FAIL: ordinary CAD frame did not publish an unbarriered batch\n");
	return 1;
    }
    return 0;
}

static BObolLodCapacitySearchKey
capacity_search_key(size_t demand, uint64_t view = 3,
    uint64_t maximumNanoseconds = 50000000ULL)
{
    BObolLodCapacitySearchKey key;
    key.inventory.set(1);
    key.availability.set(2);
    key.view.set(view);
    key.policy.set(4);
    key.preferredTargetNanoseconds = 50000000ULL;
    key.maximumTargetNanoseconds = maximumNanoseconds;
    key.maximumCandidateBudget = demand;
    return key;
}

static int
test_capacity_search_certificate(void)
{
    using Certificate = BObolLodCapacitySearchCertificate;
    using Result = Certificate::Result;
    static constexpr size_t demand = 800;
    static constexpr size_t initialCandidate = 400;

    const BObolLodCapacitySearchKey pointProblem =
	capacity_search_key(demand);
    BObolLodCapacitySearchKey boxProblem = pointProblem;
    boxProblem.pointProxyPixelThreshold = 8.0f;
    if (pointProblem.sameProblem(boxProblem)) {
	std::fprintf(stderr,
	    "FAIL: capacity certificate crossed aggregate classification\n");
	return 1;
    }

    /* The current pass's optional allowances do not bound later allocator
     * requests.  Capacity search owns a sequence of complete scene-budget
     * allocations, so its finite endpoint is resident pixel demand. */
    BObolRetainedAllocationResult allocationBoundary;
    allocationBoundary.pixelDemandPresentationCost = demand;
    allocationBoundary.maximumMarginalBudget = 500;
    if (allocationBoundary.maximumCapacitySearchBudget() != demand) {
	std::fprintf(stderr,
	    "FAIL: current marginal allowance truncated capacity domain\n");
	return 1;
    }
    allocationBoundary.maximumProtectedBudget = 650;
    if (allocationBoundary.maximumCapacitySearchBudget() != demand) {
	std::fprintf(stderr,
	    "FAIL: current protected allowance truncated capacity domain\n");
	return 1;
    }
    allocationBoundary.pixelDemandPresentationCost = 600;
    if (allocationBoundary.maximumCapacitySearchBudget() != 600) {
	std::fprintf(stderr,
	    "FAIL: capacity search exceeded view demand\n");
	return 1;
    }

    /* Exhaust the numeric classifier boundary.  The certificate need not
     * identify every integer budget: each budget denotes one member of the
     * bounded, complete-allocation candidate set.  It must never repeat a
     * member, cross the true safe/unsafe boundary, or exceed its finite
     * measurement rank. */
    for (size_t trueCapacity = 0; trueCapacity <= demand;
	 trueCapacity += 25) {
	Certificate certificate;
	Certificate::Observation observation;
	observation.key = capacity_search_key(demand);
	observation.candidateBudget = initialCandidate;
	observation.knownSafeBudget = 0;
	std::array<size_t, Certificate::candidateLimit()> visited {{}};
	unsigned int visitedCount = 0;
	Certificate::Decision decision;
	for (unsigned int step = 0;
	     step < Certificate::candidateLimit() *
		Certificate::sampleLimit() + 1; ++step) {
	    observation.presentedCost = observation.candidateBudget;
	    observation.validSample = true;
	    observation.observedNanoseconds =
		observation.candidateBudget <= trueCapacity ?
		    certificate.activeTargetNanoseconds() / 2 :
		    certificate.activeTargetNanoseconds() * 2;
	    decision = certificate.observe(observation);
	    if (decision.requestsFrame())
		continue;
	    if (decision.requestsReallocation()) {
		for (unsigned int i = 0; i < visitedCount; ++i) {
		    if (visited[i] == decision.budget) {
			std::fprintf(stderr,
			    "FAIL: capacity search repeated candidate %zu\n",
			    decision.budget);
			return 1;
		    }
		}
		visited[visitedCount++] = decision.budget;
		observation.candidateBudget = decision.budget;
		continue;
	    }
	    break;
	}
	if (decision.result != Result::CERTIFIED ||
	    !decision.terminal() || decision.safeBudget > trueCapacity ||
	    decision.unsafeBudget <= trueCapacity ||
	    decision.measuredCandidates > Certificate::candidateLimit() ||
	    decision.totalMeasuredCandidates >
		Certificate::totalCandidateLimit()) {
	    std::fprintf(stderr,
		"FAIL: capacity search terminal true=%zu result=%u safe=%zu "
		"unsafe=%zu measured=%u\n", trueCapacity,
		static_cast<unsigned int>(decision.result),
		decision.safeBudget, decision.unsafeBudget,
		decision.measuredCandidates);
	    return 1;
	}
    }

    /* The event-driven static allowance is the second goal of the same
     * finite certificate, not an independently re-armed trial.  It inherits
     * the population proven at the preferred cadence, may remeasure a budget
     * under the longer deadline, and can advance only once. */
    {
	static constexpr size_t steadyCapacity = 250;
	static constexpr size_t staticCapacity = 650;
	Certificate certificate;
	Certificate::Observation observation;
	observation.key = capacity_search_key(demand, 3, 400000000ULL);
	observation.key.preferredBudgetCeiling = steadyCapacity;
	observation.key.maximumBudgetCeiling = staticCapacity;
	observation.candidateBudget = steadyCapacity;
	observation.knownSafeBudget = 0;
	bool observedStaticGoal = false;
	bool observedStaticCandidateAbovePreferred = false;
	Certificate::Decision decision;
	for (unsigned int step = 0;
	     step < 2 * Certificate::candidateLimit() *
		Certificate::sampleLimit() + 2; ++step) {
	    const bool staticGoal = certificate.goal() ==
		Certificate::Goal::STATIC;
	    observedStaticGoal = observedStaticGoal || staticGoal;
	    observedStaticCandidateAbovePreferred =
		observedStaticCandidateAbovePreferred ||
		(staticGoal && observation.candidateBudget > steadyCapacity);
	    const size_t trueCapacity = staticGoal ?
		staticCapacity : steadyCapacity;
	    observation.presentedCost = observation.candidateBudget;
	    observation.validSample = true;
	    observation.observedNanoseconds =
		observation.candidateBudget <= trueCapacity ?
		    certificate.activeTargetNanoseconds() / 2 :
		    certificate.activeTargetNanoseconds() * 2;
	    decision = certificate.observe(observation);
	    if (decision.requestsReallocation())
		observation.candidateBudget = decision.budget;
	    if (decision.terminal())
		break;
	}
	if (!observedStaticGoal || !observedStaticCandidateAbovePreferred ||
	    !decision.terminal() ||
	    decision.result != Result::CERTIFIED ||
	    certificate.goal() != Certificate::Goal::STATIC ||
	    decision.safeBudget > staticCapacity ||
	    decision.unsafeBudget <= staticCapacity ||
	    decision.totalMeasuredCandidates >
		Certificate::totalCandidateLimit()) {
	    std::fprintf(stderr,
		"FAIL: bounded steady-to-static capacity progression "
		"goal=%u result=%u safe=%zu unsafe=%zu\n",
		static_cast<unsigned int>(certificate.goal()),
		static_cast<unsigned int>(decision.result),
		decision.safeBudget, decision.unsafeBudget);
	    return 1;
	}
	const Certificate::Decision repeated = certificate.observe(observation);
	if (!repeated.terminal() || repeated.safeBudget != decision.safeBudget ||
	    certificate.goal() != Certificate::Goal::STATIC) {
	    std::fprintf(stderr,
		"FAIL: terminal static capacity certificate reopened\n");
	    return 1;
	}
    }

    /* A throughput miss proposes the next candidate directly.  A blind
     * midpoint is safe but throws away all three completed measurements and
     * makes a large discrete mesh visibly walk through coarse populations. */
    {
	Certificate certificate;
	Certificate::Observation observation;
	observation.key = capacity_search_key(demand);
	observation.candidateBudget = initialCandidate;
	observation.presentedCost = initialCandidate;
	observation.validSample = true;
	observation.observedNanoseconds =
	    observation.key.preferredTargetNanoseconds * 2;
	Certificate::Decision decision;
	for (unsigned int i = 0; i < Certificate::sampleLimit(); ++i)
	    decision = certificate.observe(observation);
	const size_t expectedProposal = initialCandidate * 2 / 5;
	if (!decision.requestsReallocation() ||
	    decision.budget != expectedProposal) {
	    std::fprintf(stderr,
		"FAIL: capacity miss ignored measured proposal: got=%zu "
		"expected=%zu\n", decision.budget, expectedProposal);
	    return 1;
	}
    }

    /* An interrupted capacity frame is a typed unsafe observation of the
     * active candidate.  Repeated hard misses must consume the same bounded
     * certificate rather than changing a ceiling, rekeying, and restarting
     * the search forever. */
    {
	Certificate certificate;
	Certificate::Observation deadlineObservation;
	deadlineObservation.key = capacity_search_key(demand);
	deadlineObservation.candidateBudget = initialCandidate;
	deadlineObservation.presentedCost = initialCandidate;
	deadlineObservation.validSample = true;
	deadlineObservation.observedNanoseconds =
	    deadlineObservation.key.preferredTargetNanoseconds / 2;
	Certificate::Decision deadlineDecision =
	    certificate.observe(deadlineObservation);
	unsigned int priorMeasured = 0;
	for (unsigned int miss = 0;
	     miss < Certificate::candidateLimit() &&
		 !deadlineDecision.terminal(); ++miss) {
	    const size_t candidate = certificate.candidateBudget();
	    const size_t reduction = std::max<size_t>(1, candidate / 20);
	    deadlineObservation.key.preferredBudgetCeiling = std::max(
		certificate.safeBudget(), candidate - reduction);
	    deadlineDecision = certificate.observeDeadlineMiss(
		deadlineObservation.key);
	    if (certificate.measuredCandidateCount() <= priorMeasured) {
		std::fprintf(stderr,
		    "FAIL: deadline miss restarted capacity certificate\n");
		return 1;
	    }
	    priorMeasured = certificate.measuredCandidateCount();
	    if (!deadlineDecision.requestsReallocation())
		continue;
	    deadlineObservation.candidateBudget = deadlineDecision.budget;
	    deadlineObservation.presentedCost = deadlineDecision.budget;
	    deadlineObservation.validSample = false;
	    deadlineObservation.observedNanoseconds = 0;
	    deadlineDecision = certificate.prepare(deadlineObservation);
	    if (!deadlineDecision.requestsFrame() ||
		certificate.phase() != Certificate::Phase::PRESENTING ||
		certificate.measuredCandidateCount() != priorMeasured) {
		std::fprintf(stderr,
		    "FAIL: deadline successor rekeyed capacity search\n");
		return 1;
	    }
	    /* The successor must become presentation-current before another hard
	     * deadline can classify it.  prepare(), unlike observe(), binds that
	     * exact population without consuming an ordinary timing sample. */
	    deadlineObservation.validSample = true;
	    deadlineObservation.presentedCost = deadlineDecision.budget;
	    deadlineDecision = certificate.prepare(deadlineObservation);
	    if (!deadlineDecision.requestsFrame() ||
		certificate.phase() != Certificate::Phase::MEASURING) {
		std::fprintf(stderr,
		    "FAIL: exact deadline successor did not enter measurement\n");
		return 1;
	    }
	}
	if (!deadlineDecision.terminal() ||
	    deadlineDecision.result != Result::CERTIFIED ||
	    certificate.measuredCandidateCount() >
		Certificate::candidateLimit()) {
	    std::fprintf(stderr,
		"FAIL: repeated deadline misses did not terminate "
		"result=%u measured=%u\n",
		static_cast<unsigned int>(deadlineDecision.result),
		certificate.measuredCandidateCount());
	    return 1;
	}
    }

    /* Once a safe population exists, a rejected richer population permits one
     * midpoint bridge on the final goal.  This recovers a wide unused bracket
     * without exposing an open-ended binary search to the framebuffer. */
    {
	Certificate certificate;
	Certificate::Observation bounded;
	bounded.key = capacity_search_key(demand);
	bounded.candidateBudget = 200;
	bounded.knownSafeBudget = 200;
	bounded.presentedCost = 200;
	bounded.validSample = true;
	bounded.observedNanoseconds =
	    bounded.key.preferredTargetNanoseconds / 2;
	Certificate::Decision boundedDecision;
	for (unsigned int i = 0; i < Certificate::sampleLimit(); ++i)
	    boundedDecision = certificate.observe(bounded);
	if (!boundedDecision.requestsReallocation() ||
	    boundedDecision.budget <= 200) {
	    std::fprintf(stderr,
		"FAIL: safe capacity population did not request enrichment\n");
	    return 1;
	}
	bounded.candidateBudget = boundedDecision.budget;
	bounded.presentedCost = boundedDecision.budget;
	bounded.observedNanoseconds =
	    bounded.key.preferredTargetNanoseconds * 2;
	for (unsigned int i = 0; i < Certificate::sampleLimit(); ++i)
	    boundedDecision = certificate.observe(bounded);
	const size_t rejectedBudget = bounded.candidateBudget;
	const size_t expectedBridge = 200 + (rejectedBudget - 200) / 2;
	if (!boundedDecision.requestsReallocation() ||
	    boundedDecision.budget != expectedBridge) {
	    std::fprintf(stderr,
		"FAIL: rejected richer population did not request one midpoint "
		"bridge result=%u budget=%zu expected=%zu safe=%zu unsafe=%zu "
		"measured=%u\n",
		static_cast<unsigned int>(boundedDecision.result),
		boundedDecision.budget, expectedBridge,
		boundedDecision.safeBudget, boundedDecision.unsafeBudget,
		boundedDecision.measuredCandidates);
	    return 1;
	}
	bounded.candidateBudget = boundedDecision.budget;
	bounded.presentedCost = boundedDecision.budget;
	bounded.observedNanoseconds =
	    bounded.key.preferredTargetNanoseconds / 2;
	for (unsigned int i = 0; i < Certificate::sampleLimit(); ++i)
	    boundedDecision = certificate.observe(bounded);
	if (!boundedDecision.terminal() ||
	    boundedDecision.result != Result::CERTIFIED ||
	    boundedDecision.safeBudget != expectedBridge ||
	    certificate.measuredCandidateCount() != 3) {
	    std::fprintf(stderr,
		"FAIL: final midpoint bridge did not settle proven visual\n");
	    return 1;
	}
    }

    /* A ready retained-view handoff uses the event-driven deadline as its only
     * goal.  A population which already meets that deadline is a safe lower
     * bound: the search may enrich it, but must never replace it with a
     * coarser calibration candidate.  Zoom differs from pose in demand
     * planning, not in this capacity-start contract. */
    {
	Certificate certificate;
	Certificate::Observation observation;
	observation.key = capacity_search_key(demand, 3, 400000000ULL);
	observation.key.startAtStatic = true;
	BObolLodCapacitySearchKey steadyKey = observation.key;
	steadyKey.startAtStatic = false;
	if (observation.key.sameProblem(steadyKey)) {
	    std::fprintf(stderr,
		"FAIL: static-start choice is absent from capacity problem key\n");
	    return 1;
	}
	observation.candidateBudget = initialCandidate;
	observation.knownSafeBudget = initialCandidate;
	observation.presentedCost = initialCandidate;
	observation.validSample = true;
	observation.observedNanoseconds =
	    observation.key.maximumTargetNanoseconds / 2;
	Certificate::Decision decision;
	for (unsigned int step = 0;
	     step < Certificate::candidateLimit() * Certificate::sampleLimit() + 1;
	     ++step) {
	    decision = certificate.observe(observation);
	    if (decision.requestsReallocation()) {
		if (decision.budget < initialCandidate) {
		    std::fprintf(stderr,
			"FAIL: retained capacity search coarsened a safe "
			"population\n");
		    return 1;
		}
		observation.candidateBudget = decision.budget;
		observation.presentedCost = decision.budget;
	    }
	    if (decision.terminal())
		break;
	}
	if (!decision.terminal() || certificate.goal() != Certificate::Goal::STATIC ||
	    certificate.activeTargetNanoseconds() !=
		observation.key.maximumTargetNanoseconds ||
	    decision.safeBudget < initialCandidate) {
	    std::fprintf(stderr,
		"FAIL: retained capacity search lost its static lower bound\n");
	    return 1;
	}
    }

    /* The allocator cannot express a candidate below the complete visible
     * minimum.  Once that floor misses the preferred cadence, the search must
     * advance to the static goal without requesting an impossible smaller
     * allocation or repeating the floor. */
    {
	Certificate certificate;
	Certificate::Observation observation;
	observation.key = capacity_search_key(demand, 3, 400000000ULL);
	observation.key.minimumBudget = 300;
	observation.candidateBudget = 400;
	Certificate::Decision decision;
	bool observedMinimum = false;
	bool observedStaticGoal = false;
	for (unsigned int step = 0;
	     step < 2 * Certificate::candidateLimit() *
		Certificate::sampleLimit() + 2; ++step) {
	    observation.presentedCost = observation.candidateBudget;
	    observation.validSample = true;
	    const bool staticGoal = certificate.goal() ==
		Certificate::Goal::STATIC;
	    observedStaticGoal = observedStaticGoal || staticGoal;
	    observation.observedNanoseconds = staticGoal ?
		certificate.activeTargetNanoseconds() / 2 :
		observation.key.maximumTargetNanoseconds * 2;
	    decision = certificate.observe(observation);
	    if (decision.requestsReallocation()) {
		if (decision.budget < observation.key.minimumBudget) {
		    std::fprintf(stderr,
			"FAIL: capacity search crossed allocation floor\n");
		    return 1;
		}
		observedMinimum = observedMinimum ||
		    decision.budget == observation.key.minimumBudget;
		observation.candidateBudget = decision.budget;
	    }
	    if (decision.terminal())
		break;
	}
	if (!observedMinimum || !observedStaticGoal || !decision.terminal() ||
	    certificate.goal() != Certificate::Goal::STATIC ||
	    decision.safeBudget < observation.key.minimumBudget) {
	    std::fprintf(stderr,
		"FAIL: irreducible minimum did not advance to static goal\n");
	    return 1;
	}
    }

    /* A completed steady sample also supplies evidence for the ordered,
     * longer static deadline.  If a population misses only the steady target,
     * carry that exact visual into the static goal instead of coarsening and
     * reconstructing it. */
    {
	Certificate certificate;
	Certificate::Observation handoff;
	handoff.key = capacity_search_key(demand, 3, 400000000ULL);
	handoff.candidateBudget = initialCandidate;
	handoff.presentedCost = initialCandidate;
	handoff.validSample = true;
	handoff.observedNanoseconds = 100000000ULL;
	Certificate::Decision handoffDecision;
	for (unsigned int i = 0; i < Certificate::sampleLimit(); ++i)
	    handoffDecision = certificate.observe(handoff);
	if (certificate.goal() != Certificate::Goal::STATIC ||
	    certificate.safeBudget() < initialCandidate ||
	    (handoffDecision.requestsReallocation() &&
	     handoffDecision.budget <= initialCandidate)) {
	    std::fprintf(stderr,
		"FAIL: static goal discarded cross-deadline frame evidence\n");
	    return 1;
	}
    }

    /* A semantic input edge starts a new certificate.  The capacity-output
     * epoch is intentionally absent from the key, but camera/input epochs
     * and the complete demand are not. */
    Certificate resetCertificate;
    Certificate::Observation observation;
    observation.key = capacity_search_key(demand);
    observation.candidateBudget = 200;
    observation.presentedCost = 200;
    observation.validSample = true;
    observation.observedNanoseconds = 10000000ULL;
    Certificate::Decision decision = resetCertificate.observe(observation);
    if (!decision.requestsFrame() || decision.samplesRemaining != 2)
	return 1;
    observation.key = capacity_search_key(demand, 9);
    decision = resetCertificate.observe(observation);
    if (!decision.requestsFrame() || decision.samplesRemaining != 2 ||
	resetCertificate.measuredCandidateCount() != 0) {
	std::fprintf(stderr,
	    "FAIL: capacity search did not reset on semantic revision\n");
	return 1;
    }

    /* An unchanged candidate cannot silently change its population, and
     * invalid presentations cannot manufacture an unbounded repaint loop. */
    Certificate staleCertificate;
    observation = Certificate::Observation();
    observation.key = capacity_search_key(demand);
    observation.candidateBudget = 200;
    observation.presentedCost = 150;
    observation.validSample = true;
    observation.observedNanoseconds = 10000000ULL;
    (void)staleCertificate.observe(observation);
    observation.presentedCost = 151;
    decision = staleCertificate.observe(observation);
    if (decision.result != Result::STALE_POPULATION || !decision.terminal()) {
	std::fprintf(stderr,
	    "FAIL: capacity search accepted a changing candidate population\n");
	return 1;
    }

    /* Numeric scene budgets are not discrete population identities.  A mesh
     * cut staircase can map a wide interval of budgets to exactly the same
     * occurrence cuts.  Once that population has three measurements, a
     * successor budget selecting the same exact identity must reuse the result
     * immediately instead of repainting the unchanged view three more times. */
    {
	Certificate equivalentCertificate;
	Certificate::Observation equivalent;
	equivalent.key = capacity_search_key(demand);
	equivalent.candidateBudget = initialCandidate;
	equivalent.presentedCost = 200;
	equivalent.populationDigest = 0x4c554359ULL;
	equivalent.populationIdentity = 0x1d3f7a91ULL;
	equivalent.validSample = true;
	equivalent.observedNanoseconds =
	    equivalent.key.preferredTargetNanoseconds / 2;
	Certificate::Decision equivalentDecision;
	for (unsigned int i = 0; i < Certificate::sampleLimit(); ++i)
	    equivalentDecision = equivalentCertificate.observe(equivalent);
	if (!equivalentDecision.requestsReallocation()) {
	    std::fprintf(stderr,
		"FAIL: equivalent-population test did not request a successor\n");
	    return 1;
	}
	equivalent.candidateBudget = equivalentDecision.budget;
	equivalent.populationDigest ^= 0x55aa55aaULL;
	const size_t repeatedPopulationBudget = equivalent.candidateBudget;
	equivalentDecision = equivalentCertificate.observe(equivalent);
	if (equivalentDecision.requestsFrame() ||
	    equivalentDecision.samplesRemaining != 0 ||
	    equivalentCertificate.measuredCandidateCount() != 2 ||
	    (equivalentDecision.requestsReallocation() &&
	     equivalentDecision.budget <= repeatedPopulationBudget + 1)) {
	    std::fprintf(stderr,
		"FAIL: equivalent capacity population did not make bounded "
		"progress\n");
	    return 1;
	}
    }

    /* A compact digest is diagnostic only.  Equal digests may not authorize
     * a different exact population identity. */
    {
	Certificate collisionCertificate;
	Certificate::Observation collision;
	collision.key = capacity_search_key(demand);
	collision.candidateBudget = initialCandidate;
	collision.presentedCost = 200;
	collision.populationDigest = 73;
	collision.populationIdentity = 101;
	collision.validSample = true;
	collision.observedNanoseconds =
	    collision.key.preferredTargetNanoseconds / 2;
	(void)collisionCertificate.observe(collision);
	collision.populationIdentity = 102;
	const Certificate::Decision collisionDecision =
	    collisionCertificate.observe(collision);
	if (collisionDecision.result != Result::STALE_POPULATION) {
	    std::fprintf(stderr,
		"FAIL: equal digest authorized a different exact population\n");
	    return 1;
	}
    }

    /* A completed retained allocation knows the exact numeric interval which
     * selects its discrete PoP population.  Classifying that population must
     * classify the whole interval: otherwise a 50k-object scene can spend an
     * O(N) allocation pass proving that candidate+1 selects unchanged cuts. */
    {
	Certificate safeIntervalCertificate;
	Certificate::Observation interval;
	interval.key = capacity_search_key(demand);
	interval.key.startAtStatic = true;
	interval.candidateBudget = initialCandidate;
	interval.presentedCost = 200;
	interval.populationMinimumBudget = 200;
	interval.nextDistinctPopulationBudget = 700;
	interval.nextDistinctPopulationBudgetKnown = true;
	interval.validSample = true;
	interval.observedNanoseconds =
	    interval.key.maximumTargetNanoseconds / 2;
	Certificate::Decision intervalDecision;
	for (unsigned int i = 0; i < Certificate::sampleLimit(); ++i)
	    intervalDecision = safeIntervalCertificate.observe(interval);
	if (!intervalDecision.requestsReallocation() ||
	    intervalDecision.safeBudget != 699 ||
	    intervalDecision.budget < 700) {
	    std::fprintf(stderr,
		"FAIL: safe population interval requested a redundant allocation "
		"(safe=%zu next=%zu result=%u)\n",
		intervalDecision.safeBudget, intervalDecision.budget,
		static_cast<unsigned int>(intervalDecision.result));
	    return 1;
	}

	Certificate unsafeIntervalCertificate;
	interval.observedNanoseconds =
	    interval.key.maximumTargetNanoseconds * 2;
	for (unsigned int i = 0; i < Certificate::sampleLimit(); ++i)
	    intervalDecision = unsafeIntervalCertificate.observe(interval);
	if (!intervalDecision.requestsReallocation() ||
	    intervalDecision.unsafeBudget != 200 ||
	    intervalDecision.budget >= 200) {
	    std::fprintf(stderr,
		"FAIL: unsafe population interval requested a redundant allocation "
		"(unsafe=%zu next=%zu result=%u)\n",
		intervalDecision.unsafeBudget, intervalDecision.budget,
		static_cast<unsigned int>(intervalDecision.result));
	    return 1;
	}

	Certificate endpointCertificate;
	interval.nextDistinctPopulationBudget = 0;
	interval.observedNanoseconds =
	    interval.key.maximumTargetNanoseconds / 2;
	for (unsigned int i = 0; i < Certificate::sampleLimit(); ++i)
	    intervalDecision = endpointCertificate.observe(interval);
	if (!intervalDecision.terminal() ||
	    intervalDecision.result != Result::CERTIFIED ||
	    intervalDecision.safeBudget != demand) {
	    std::fprintf(stderr,
		"FAIL: richest safe population did not terminate at demand "
		"(safe=%zu result=%u)\n", intervalDecision.safeBudget,
		static_cast<unsigned int>(intervalDecision.result));
	    return 1;
	}
    }

    /* Allocation and rendering use different cost currencies.  In
     * particular, OSMesa expands an indexed retained allocation into a
     * submitted position stream.  The allocation barrier may therefore not
     * bind the completed-frame population identity; the first rendered frame
     * owns that value. */
    Certificate currencyCertificate;
    observation = Certificate::Observation();
    observation.key = capacity_search_key(demand);
    observation.candidateBudget = 200;
    observation.presentedCost = 200;
    observation.validSample = true;
    observation.observedNanoseconds = 100000000ULL;
    for (unsigned int i = 0; i < Certificate::sampleLimit(); ++i)
	decision = currencyCertificate.observe(observation);
    if (!decision.requestsReallocation()) {
	std::fprintf(stderr,
	    "FAIL: capacity currency test did not request a successor\n");
	return 1;
    }
    observation.candidateBudget = decision.budget;
    observation.presentedCost = decision.budget;
    observation.validSample = false;
    observation.observedNanoseconds = 0;
    decision = currencyCertificate.prepare(observation);
    if (!decision.requestsFrame() ||
	currencyCertificate.phase() != Certificate::Phase::PRESENTING) {
	std::fprintf(stderr,
	    "FAIL: allocation barrier consumed a capacity measurement\n");
	return 1;
    }
    observation.presentedCost = decision.budget * 2;
    observation.validSample = true;
    observation.observedNanoseconds = 10000000ULL;
    decision = currencyCertificate.observe(observation);
    if (!decision.requestsFrame() ||
	decision.result == Result::STALE_POPULATION ||
	decision.samplesRemaining != Certificate::sampleLimit() - 1) {
	std::fprintf(stderr,
	    "FAIL: allocation cost was mistaken for rendered population cost\n");
	return 1;
    }

    Certificate invalidCertificate;
    observation = Certificate::Observation();
    observation.key = capacity_search_key(demand);
    observation.candidateBudget = 200;
	observation.presentedCost = 200;
	observation.observedNanoseconds = 10000000ULL;
	observation.validSample = true;
	decision = invalidCertificate.observe(observation);
    observation.presentedCost = 0;
    observation.observedNanoseconds = 0;
    observation.validSample = false;
	uint64_t previousInvalidProgress =
	    invalidCertificate.progressCompletedUnits();
    for (unsigned int i = 0; i <= Certificate::invalidSampleLimit(); ++i) {
	decision = invalidCertificate.observe(observation);
	const uint64_t currentInvalidProgress =
	    invalidCertificate.progressCompletedUnits();
	if (currentInvalidProgress <= previousInvalidProgress) {
	    std::fprintf(stderr,
		"FAIL: bounded invalid measurement did not advance progress\n");
	    return 1;
	}
	previousInvalidProgress = currentInvalidProgress;
    }
    if (decision.result != Result::UNMEASURABLE || !decision.terminal()) {
	std::fprintf(stderr,
	    "FAIL: capacity search invalid samples did not terminate\n");
	return 1;
    }

    /* Planning may inspect a candidate while its presentation producer is
     * pending without consuming the retry allowance.  Completed renderer
     * attempts which still carry no exact CAD work are different: they must
     * retire the search after a finite number of requests. */
    Certificate presentationCertificate;
    observation = Certificate::Observation();
    observation.key = capacity_search_key(demand);
    observation.candidateBudget = 200;
    for (unsigned int i = 0; i <= Certificate::invalidSampleLimit(); ++i) {
	decision = presentationCertificate.prepare(observation);
	if (!decision.requestsFrame() || decision.terminal()) {
	    std::fprintf(stderr,
		"FAIL: capacity planning consumed presentation retries\n");
	    return 1;
	}
    }
	previousInvalidProgress =
	    presentationCertificate.progressCompletedUnits();
    for (unsigned int i = 0; i <= Certificate::invalidSampleLimit(); ++i) {
	decision = presentationCertificate.observe(observation);
	const uint64_t currentInvalidProgress =
	    presentationCertificate.progressCompletedUnits();
	if (currentInvalidProgress <= previousInvalidProgress) {
	    std::fprintf(stderr,
		"FAIL: bounded inexact presentation did not advance progress\n");
	    return 1;
	}
	previousInvalidProgress = currentInvalidProgress;
    }
    if (decision.result != Result::UNMEASURABLE || !decision.terminal()) {
	std::fprintf(stderr,
	    "FAIL: inexact capacity presentations did not terminate\n");
	return 1;
    }
    return 0;
}

static int
test_point_quality_phase(void)
{
    BObolLodPointQualityPhase phase;


    if (phase.pending() || phase.adaptiveCalibrationPending() ||
	phase.staticCalibrationPending() ||
	phase.handoffConfirmationPending())
	return 1;

    phase.requestCalibration();
    if (!phase.presentationPending() ||
	!phase.adaptiveCalibrationPending() ||
	phase.staticCalibrationPending() ||
	phase.handoffConfirmationPending() ||
	phase.triangleRecoveryPending())
	return 1;

    phase.requestStaticCalibration();
    phase.requestCalibration();
    if (!phase.presentationPending() ||
	phase.adaptiveCalibrationPending() ||
	!phase.staticCalibrationPending() ||
	phase.handoffConfirmationPending()) {
	std::fprintf(stderr,
	    "FAIL: ordinary point calibration displaced static deadline ownership\n");
	return 1;
    }

    phase.requestHandoffConfirmation();
    if (!phase.presentationPending() ||
	phase.adaptiveCalibrationPending() ||
	phase.staticCalibrationPending() ||
	!phase.handoffConfirmationPending()) {
	std::fprintf(stderr,
	    "FAIL: point handoff confirmation retained adaptive census authority\n");
	return 1;
    }

    phase.beginTriangleRecovery();
    phase.requestCalibration();
    phase.requestHandoffConfirmation();
    phase.completeCalibration();
    if (phase.presentationPending() || !phase.triangleRecoveryPending()) {
	std::fprintf(stderr,
	    "FAIL: point calibration displaced active triangle recovery\n");
	return 1;
    }
    if (!phase.requiresRecoveryPresentation(true, true, false, false) ||
	phase.requiresRecoveryPresentation(false, true, false, false) ||
	phase.requiresRecoveryPresentation(true, false, false, false) ||
	phase.requiresRecoveryPresentation(true, true, true, false) ||
	phase.requiresRecoveryPresentation(true, true, false, true)) {
	std::fprintf(stderr,
	    "FAIL: triangle recovery lost its changed-pass presentation owner\n");
	return 1;
    }
    if (!phase.recoveryAdmissionPending(true, false) ||
	phase.recoveryAdmissionPending(true, true) ||
	phase.recoveryAdmissionPending(false, false)) {
	std::fprintf(stderr,
	    "FAIL: triangle recovery confused residual quality debt with pending admission\n");
	return 1;
    }

    phase.completeTriangleRecovery(42);
    if (phase.pending() || !phase.triangleRecoverySaturatedBy(42) ||
	phase.triangleRecoverySaturatedBy(0) ||
	phase.triangleRecoverySaturatedBy(41))
	return 1;
    phase.requestCalibration();
    phase.completeCalibration();
    if (phase.pending() || !phase.triangleRecoverySaturatedBy(42))
	return 1;
    phase.reset();
    if (phase.triangleRecoverySaturatedBy(42))
	return 1;
    return 0;
}

int
main(int argc, char **argv)
{
    (void)argc;
    bu_setprogname(argv[0]);
    if (test_identity_exhaustion_contract())
	return 1;
    if (test_exact_sample_identity())
	return 1;
    if (test_submission_pass())
	return 1;
    if (test_deadline_safe_presentation())
	return 1;
    if (test_retained_quality_evidence())
	return 1;
    if (test_renderer_performance_evidence())
	return 1;
    if (test_submission_intent())
	return 1;
    if (test_retained_view_continuity())
	return 1;
    if (test_retained_pass_annotations())
	return 1;
    if (test_interrupted_presentation_replay())
	return 1;
    if (test_exact_presentation_frame())
	return 1;
    if (test_point_admission_frame())
	return 1;
    if (test_structural_repair())
	return 1;
    if (test_planning_obligations())
	return 1;
    if (test_discrete_trial_permit())
	return 1;
    if (test_control_refinement())
	return 1;
    if (test_epochs())
	return 1;
    if (test_interaction_session())
	return 1;
    if (test_revision_contract())
	return 1;
    if (test_resource_policy())
	return 1;
    if (test_headroom_policy())
	return 1;
    if (test_coverage_policy())
	return 1;
    if (test_visibility_census())
	return 1;
    if (test_source_profile_policy())
	return 1;
    if (test_projected_demand_cache())
	return 1;
    if (test_interaction_start_certificate())
	return 1;
    if (test_static_quality_policy())
	return 1;
    if (test_terminal_reconciliation_composition())
	return 1;
    if (test_delivery_policy())
	return 1;
    if (test_convergence_policy())
	return 1;
    if (test_availability_scheduler())
	return 1;
    if (test_point_structural_preload_threshold())
	return 1;
    if (test_retained_allocation_request())
	return 1;
    if (test_admission_capacity())
	return 1;
    if (test_capacity_search_certificate())
	return 1;
    if (test_point_quality_phase())
	return 1;
    if (test_quality_policy())
	return 1;
    if (test_view_demand_policy())
	return 1;
    if (test_quiet_successor_reducer())
	return 1;
    if (test_presentation_policy())
	return 1;
    if (test_completed_pass_composed_lifecycle())
	return 1;
    if (test_view_quality_history())
	return 1;
    if (test_availability_and_publication())
	return 1;
    if (test_availability_growth())
	return 1;
    if (test_compaction_policy())
	return 1;
    if (test_presentation_transaction())
	return 1;
    if (test_seeded_composed_lifecycle())
	return 1;
    return 0;
}
