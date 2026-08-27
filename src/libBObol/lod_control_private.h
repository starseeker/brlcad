/*                L O D _ C O N T R O L _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_CONTROL_PRIVATE_H
#define LIBBOBOL_LOD_CONTROL_PRIVATE_H

#include <cstddef>
#include <cstdint>
#include <limits>

/*
 * One submission cursor may be active while a complete successor rescan is
 * owed.  Discovery can also pause the cursor while preserving that rescan
 * until a later inventory revision supplies more entries.  These four states
 * were formerly represented by two freely written booleans in three
 * translation units.
 */
class BObolLodSubmissionPass {
public:
    enum class State : uint8_t {
	IDLE = 0,
	ACTIVE,
	IDLE_RESCAN,
	ACTIVE_RESCAN
    };

    void activate(void)
    {
	this->stateValue = this->rescanPending() ?
	    State::ACTIVE_RESCAN : State::ACTIVE;
    }

    void deactivate(void)
    {
	this->stateValue = this->rescanPending() ?
	    State::IDLE_RESCAN : State::IDLE;
    }

    void setActive(bool active)
    {
	if (active)
	    this->activate();
	else
	    this->deactivate();
    }

    void requestRescan(void)
    {
	this->stateValue = this->active() ?
	    State::ACTIVE_RESCAN : State::IDLE_RESCAN;
    }

    void clearRescan(void)
    {
	this->stateValue = this->active() ? State::ACTIVE : State::IDLE;
    }

    void setRescanPending(bool pending)
    {
	if (pending)
	    this->requestRescan();
	else
	    this->clearRescan();
    }

    void reset(void)
    {
	this->stateValue = State::IDLE;
    }

    bool active(void) const
    {
	return this->stateValue == State::ACTIVE ||
	    this->stateValue == State::ACTIVE_RESCAN;
    }

    bool rescanPending(void) const
    {
	return this->stateValue == State::IDLE_RESCAN ||
	    this->stateValue == State::ACTIVE_RESCAN;
    }

    State state(void) const
    {
	return this->stateValue;
    }

private:
    State stateValue = State::IDLE;
};

/* Immutable-in-meaning configuration retained while a bounded submission
 * cursor is active.  Mode, refresh behavior, and reset intent used to be
 * independent fields, allowing a resumed cursor to inherit only part of the
 * request which created it. */
class BObolLodSubmissionIntent {
public:
    enum class Mode : uint8_t {
	ORDINARY = 0,
	RETAINED_ADMISSION
    };

    void configure(bool refreshMissing, bool resetExisting)
    {
	this->refreshMissingValue = refreshMissing;
	this->resetExistingValue = resetExisting;
    }

    void setRetainedAdmission(bool enabled)
    {
	this->modeValue = enabled ?
	    Mode::RETAINED_ADMISSION : Mode::ORDINARY;
    }

    void reset(void)
    {
	this->modeValue = Mode::ORDINARY;
	this->refreshMissingValue = true;
	this->resetExistingValue = false;
    }

    bool retainedAdmission(void) const
    {
	return this->modeValue == Mode::RETAINED_ADMISSION;
    }

    bool refreshMissing(void) const
    {
	return this->refreshMissingValue;
    }

    bool resetExisting(void) const
    {
	return this->resetExistingValue;
    }

    Mode mode(void) const
    {
	return this->modeValue;
    }

private:
    Mode modeValue = Mode::ORDINARY;
    bool refreshMissingValue = true;
    bool resetExistingValue = false;
};

/* Pose-only continuity proof for one interaction handoff.  Cut retention and
 * its deferred exact-visibility census are related but not interchangeable:
 * motion may keep a proven presentation immediately, while quiet readiness
 * still waits for a current-camera census. */
class BObolLodPoseContinuity {
public:
    void setRetainOccurrenceCuts(bool retain)
    {
	this->retainOccurrenceCutsValue = retain;
    }

    void clearRetainOccurrenceCuts(void)
    {
	this->retainOccurrenceCutsValue = false;
    }

    bool retainOccurrenceCuts(void) const
    {
	return this->retainOccurrenceCutsValue;
    }

    void deferVisibilityCensus(void)
    {
	this->visibilityCensusDeferredValue = true;
    }

    void completeVisibilityCensus(void)
    {
	this->visibilityCensusDeferredValue = false;
    }

    bool visibilityCensusDeferred(void) const
    {
	return this->visibilityCensusDeferredValue;
    }

    void reset(void)
    {
	this->retainOccurrenceCutsValue = false;
	this->visibilityCensusDeferredValue = false;
    }

private:
    bool retainOccurrenceCutsValue = false;
    bool visibilityCensusDeferredValue = false;
};

/* Complete annotations produced by one bounded retained-admission pass.
 * These facts may coexist, but they share one pass lifetime.  Keeping their
 * mutation private prevents a successor pass from inheriting a convenient
 * subset of the preceding pass's cut, residency, or budget evidence. */
class BObolLodRetainedPassAnnotations {
public:
    void reset(void)
    {
	this->refinementPendingValue = false;
	this->residencyPendingValue = false;
	this->cutAdvancedValue = false;
	this->budgetBlockedValue = false;
	this->admittedWorkValue = false;
	this->missingMeshBudgetBlockedCountValue = 0;
    }

    void retireRefinement(void)
    {
	this->refinementPendingValue = false;
	this->cutAdvancedValue = false;
	this->budgetBlockedValue = false;
    }

    void noteRefinementPending(void) { this->refinementPendingValue = true; }
    void clearRefinementPending(void)
    {
	this->refinementPendingValue = false;
    }
    bool refinementPending(void) const { return this->refinementPendingValue; }

    void noteResidencyPending(void) { this->residencyPendingValue = true; }
    void clearResidencyPending(void) { this->residencyPendingValue = false; }
    bool residencyPending(void) const { return this->residencyPendingValue; }

    void noteCutAdvanced(void) { this->cutAdvancedValue = true; }
    bool cutAdvanced(void) const { return this->cutAdvancedValue; }

    void noteBudgetBlocked(void) { this->budgetBlockedValue = true; }
    bool budgetBlocked(void) const { return this->budgetBlockedValue; }

    void noteAdmittedWork(void) { this->admittedWorkValue = true; }
    void clearAdmittedWork(void) { this->admittedWorkValue = false; }
    bool admittedWork(void) const { return this->admittedWorkValue; }

    void addMissingMeshBudgetBlocked(size_t count)
    {
	const size_t maximum = (std::numeric_limits<size_t>::max)();
	this->missingMeshBudgetBlockedCountValue =
	    count > maximum - this->missingMeshBudgetBlockedCountValue ?
		maximum : this->missingMeshBudgetBlockedCountValue + count;
    }
    void clearMissingMeshBudgetBlocked(void)
    {
	this->missingMeshBudgetBlockedCountValue = 0;
    }
    size_t missingMeshBudgetBlockedCount(void) const
    {
	return this->missingMeshBudgetBlockedCountValue;
    }

private:
    bool refinementPendingValue = false;
    bool residencyPendingValue = false;
    bool cutAdvancedValue = false;
    bool budgetBlockedValue = false;
    bool admittedWorkValue = false;
    size_t missingMeshBudgetBlockedCountValue = 0;
};

/* Exclusive owner-thread gate for one interrupted retained traversal.  The
 * replay is either absent or waiting for its exact successor frame; source
 * publication cannot infer or clear this state indirectly. */
class BObolLodInterruptedPresentationReplay {
public:
    enum class State : uint8_t {
	IDLE = 0,
	AWAITING_FRAME
    };

    void arm(void) { this->stateValue = State::AWAITING_FRAME; }
    void retire(void) { this->stateValue = State::IDLE; }
    void reset(void) { this->retire(); }
    bool pending(void) const
    {
	return this->stateValue == State::AWAITING_FRAME;
    }
    State state(void) const { return this->stateValue; }

private:
    State stateValue = State::IDLE;
};

/* One classifier-changing point threshold must be presented before source
 * admission may consume it.  This is a frame obligation, not a free boolean
 * or a stable point-quality calibration phase. */
class BObolLodPointAdmissionFrame {
public:
    enum class State : uint8_t {
	IDLE = 0,
	AWAITING_FRAME
    };

    void request(void) { this->stateValue = State::AWAITING_FRAME; }
    void setPending(bool pending)
    {
	this->stateValue = pending ? State::AWAITING_FRAME : State::IDLE;
    }
    void retire(void) { this->stateValue = State::IDLE; }
    void reset(void) { this->retire(); }
    bool pending(void) const
    {
	return this->stateValue == State::AWAITING_FRAME;
    }
    State state(void) const { return this->stateValue; }

private:
    State stateValue = State::IDLE;
};

/* Concurrent one-shot planning obligations owned by the current evidence
 * epoch.  These are work requests, not outcome flags: each must be explicitly
 * requested and retired by its sole effect owner. */
class BObolLodPlanningObligations {
private:
    enum class Work : uint8_t {
	RETAINED_IMPORTANCE_CENSUS = 1u << 0,
	RESIDENT_ADMISSION_RETRY = 1u << 1
    };

    static constexpr uint8_t bit(Work work)
    {
	return static_cast<uint8_t>(work);
    }

    void request(Work work)
    {
	this->workValue |= bit(work);
    }

    void retire(Work work)
    {
	this->workValue &= static_cast<uint8_t>(~bit(work));
    }

    void setPending(Work work, bool pending)
    {
	if (pending)
	    this->request(work);
	else
	    this->retire(work);
    }

    void clearAll(void)
    {
	this->workValue = 0;
    }

    bool pending(Work work) const
    {
	return (this->workValue & bit(work)) != 0;
    }

public:
    void reset(void)
    {
	this->clearAll();
    }

    void requestImportanceCensus(void)
    {
	this->request(Work::RETAINED_IMPORTANCE_CENSUS);
    }

    void retireImportanceCensus(void)
    {
	this->retire(Work::RETAINED_IMPORTANCE_CENSUS);
    }

    void setImportanceCensus(bool pending)
    {
	this->setPending(Work::RETAINED_IMPORTANCE_CENSUS, pending);
    }

    bool importanceCensusPending(void) const
    {
	return this->pending(Work::RETAINED_IMPORTANCE_CENSUS);
    }

    void setResidentAdmissionRetry(bool pending)
    {
	this->setPending(Work::RESIDENT_ADMISSION_RETRY, pending);
    }

    void retireResidentAdmissionRetry(void)
    {
	this->retire(Work::RESIDENT_ADMISSION_RETRY);
    }

    bool residentAdmissionRetryPending(void) const
    {
	return this->pending(Work::RESIDENT_ADMISSION_RETRY);
    }

private:
    uint8_t workValue = 0;
};

/* Scene-wide permission to try one otherwise-unaffordable discrete
 * population.  Granting and consuming the permit are explicit so multiple
 * time-sliced source actions cannot each infer their own allowance. */
class BObolLodDiscreteTrialPermit {
public:
    enum class State : uint8_t {
	UNAVAILABLE = 0,
	AVAILABLE
    };

    void grant(void) { this->stateValue = State::AVAILABLE; }
    void revoke(void) { this->stateValue = State::UNAVAILABLE; }
    void consume(void) { this->revoke(); }
    void reset(void) { this->revoke(); }
    void setAvailable(bool available)
    {
	this->stateValue = available ? State::AVAILABLE : State::UNAVAILABLE;
    }
    bool available(void) const
    {
	return this->stateValue == State::AVAILABLE;
    }
    State state(void) const { return this->stateValue; }

private:
    State stateValue = State::UNAVAILABLE;
};

/*
 * Allocation-free refinement map from the remaining concrete controller
 * latches to the finite work ledger in ObolProgressivePipeline.tla.  This is
 * deliberately a projection, not another scheduler: effect ownership moves
 * behind the reducer in deletion-producing slices, while readiness and trace
 * diagnostics use this single mapping immediately.
 */
class BObolLodControlRefinement {
public:
    enum class Work : uint32_t {
	INTERACTION = 1u << 0,
	INVENTORY = 1u << 1,
	AVAILABILITY = 1u << 2,
	PUBLICATION = 1u << 3,
	PLANNING = 1u << 4,
	PRESENTATION = 1u << 5,
	HANDOFF = 1u << 6,
	COMPACTION = 1u << 7,
	CACHE_WRITE = 1u << 8
    };

    enum class Owner : uint8_t {
	NONE = 0,
	INTERACTION,
	INVENTORY,
	AVAILABILITY,
	PUBLICATION,
	PLANNING,
	PRESENTATION,
	HANDOFF,
	COMPACTION,
	CACHE_WRITE
    };

    enum class Violation : uint32_t {
	OWNERLESS_WORK = 1u << 0,
	TERMINAL_WITH_WORK = 1u << 1,
	INVALID_READINESS = 1u << 2,
	INVALID_OWNER = 1u << 3
    };

    struct Inputs {
	bool interaction = false;
	bool inventory = false;
	bool availability = false;
	bool result = false;
	bool publication = false;
	bool submission = false;
	bool submissionRescan = false;
	bool retainedAllocation = false;
	bool residentGrowth = false;
	bool pointTriangleRecovery = false;
	bool structuralFrontier = false;
	bool presentationReplay = false;
	bool presentationBarrier = false;
	bool capacityFrame = false;
	bool pointAdmissionFrame = false;
	bool pointCalibration = false;
	bool capacityCalibration = false;
	bool headroomProbe = false;
	bool handoff = false;
	bool compaction = false;
	bool cacheWrite = false;
    };

    struct Snapshot {
	uint32_t obligations = 0;
	Owner owner = Owner::NONE;

	bool has(Work obligation) const
	{
	    return (this->obligations & bit(obligation)) != 0;
	}

	bool calibrationPending(void) const
	{
	    return this->has(Work::PRESENTATION);
	}

	bool controlPending(void) const
	{
	    return this->has(Work::PLANNING) || this->has(Work::HANDOFF);
	}

	bool foregroundPending(void) const
	{
	    const uint32_t background =
		bit(Work::COMPACTION) | bit(Work::CACHE_WRITE);
	    return (this->obligations & ~background) != 0;
	}

	bool refinementFramePending(void) const
	{
	    return this->has(Work::PRESENTATION);
	}
    };

    static constexpr uint32_t bit(Work work)
    {
	return static_cast<uint32_t>(work);
    }

    static constexpr uint32_t bit(Violation violation)
    {
	return static_cast<uint32_t>(violation);
    }

    static constexpr uint32_t ownerObligation(Owner owner)
    {
	return owner == Owner::INTERACTION ? bit(Work::INTERACTION) :
	    owner == Owner::INVENTORY ? bit(Work::INVENTORY) :
	    owner == Owner::AVAILABILITY ? bit(Work::AVAILABILITY) :
	    owner == Owner::PUBLICATION ? bit(Work::PUBLICATION) :
	    owner == Owner::PLANNING ? bit(Work::PLANNING) :
	    owner == Owner::PRESENTATION ? bit(Work::PRESENTATION) :
	    owner == Owner::HANDOFF ? bit(Work::HANDOFF) :
	    owner == Owner::COMPACTION ? bit(Work::COMPACTION) :
	    owner == Owner::CACHE_WRITE ? bit(Work::CACHE_WRITE) : 0;
    }

    static Snapshot evaluate(const Inputs &inputs)
    {
	Snapshot snapshot;
	if (inputs.interaction)
	    snapshot.obligations |= bit(Work::INTERACTION);
	if (inputs.inventory)
	    snapshot.obligations |= bit(Work::INVENTORY);
	if (inputs.availability)
	    snapshot.obligations |= bit(Work::AVAILABILITY);
	if (inputs.result || inputs.publication)
	    snapshot.obligations |= bit(Work::PUBLICATION);
	if (inputs.submission || inputs.submissionRescan ||
	    inputs.retainedAllocation || inputs.residentGrowth ||
	    inputs.pointTriangleRecovery || inputs.structuralFrontier)
	    snapshot.obligations |= bit(Work::PLANNING);
	if (inputs.presentationReplay || inputs.presentationBarrier ||
	    inputs.capacityFrame || inputs.pointAdmissionFrame ||
	    inputs.pointCalibration || inputs.capacityCalibration ||
	    inputs.headroomProbe)
	    snapshot.obligations |= bit(Work::PRESENTATION);
	if (inputs.handoff)
	    snapshot.obligations |= bit(Work::HANDOFF);
	if (inputs.compaction)
	    snapshot.obligations |= bit(Work::COMPACTION);
	if (inputs.cacheWrite)
	    snapshot.obligations |= bit(Work::CACHE_WRITE);

	/* A closed interrupted presentation transaction blocks every owner-thread
	 * scene mutation until its exact target commits or returns typed capacity
	 * evidence.  The remaining order mirrors the current bounded pump and is
	 * the explicit seam along which those effects will be migrated. */
	if (inputs.presentationReplay) {
	    snapshot.owner = Owner::PRESENTATION;
	} else if (inputs.interaction) {
	    snapshot.owner = Owner::INTERACTION;
	} else if (inputs.result || inputs.publication) {
	    snapshot.owner = Owner::PUBLICATION;
	} else if (inputs.inventory) {
	    snapshot.owner = Owner::INVENTORY;
	} else if (inputs.availability) {
	    snapshot.owner = Owner::AVAILABILITY;
	} else if (snapshot.has(Work::PRESENTATION)) {
	    snapshot.owner = Owner::PRESENTATION;
	} else if (snapshot.has(Work::PLANNING)) {
	    snapshot.owner = Owner::PLANNING;
	} else if (inputs.handoff) {
	    snapshot.owner = Owner::HANDOFF;
	} else if (inputs.compaction) {
	    snapshot.owner = Owner::COMPACTION;
	} else if (inputs.cacheWrite) {
	    snapshot.owner = Owner::CACHE_WRITE;
	}
	return snapshot;
    }

    static uint32_t validate(const Snapshot &snapshot, bool terminal,
	bool viewReady, bool terminalError)
    {
	uint32_t violations = 0;
	if (snapshot.obligations != 0 && snapshot.owner == Owner::NONE)
	    violations |= bit(Violation::OWNERLESS_WORK);
	const uint32_t ownedObligation = ownerObligation(snapshot.owner);
	if (snapshot.owner != Owner::NONE &&
	    !(snapshot.obligations & ownedObligation))
	    violations |= bit(Violation::INVALID_OWNER);
	if (terminal && snapshot.foregroundPending())
	    violations |= bit(Violation::TERMINAL_WITH_WORK);
	if (viewReady && (!terminal || terminalError))
	    violations |= bit(Violation::INVALID_READINESS);
	return violations;
    }
};

#endif /* LIBBOBOL_LOD_CONTROL_PRIVATE_H */
