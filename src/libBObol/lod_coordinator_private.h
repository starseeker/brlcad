/*          L O D _ C O O R D I N A T O R _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_COORDINATOR_PRIVATE_H
#define LIBBOBOL_LOD_COORDINATOR_PRIVATE_H

#include "common.h"
#include "BObol/BLodIdentifiers.h"

#include <array>
#include <cstddef>
#include <cstdint>

/*
 * Use the same strong epoch domains as requests and results.  Keeping a
 * second coordinator-only implementation previously allowed the state
 * machine and service boundary to drift while appearing type safe.
 */
using BObolLodViewEpoch = BObolViewEpoch;
using BObolLodPolicyEpoch = BObolPolicyEpoch;
using BObolLodInventoryEpoch = BObolInventoryEpoch;
using BObolLodSourceRoutingId = BObolSourceRoutingId;

/*
 * Owner-thread phase tracker.  It stores no strings, allocates no memory, and
 * is never consulted by per-occurrence planning.  Each phase retains its most
 * recent progress witness so deterministic tests and diagnostics can prove
 * whether repeated pumps are making progress instead of merely observing a
 * collection of Boolean flags.
 */
class BObolLodPhaseTracker {
public:
    enum class Phase : uint8_t {
	FALLBACK = 0,
	COVERAGE,
	INTERACTIVE,
	SETTLING,
	STABLE,
	COMPACTING,
	COUNT
    };

    /*
     * Allocation-free phase decision inputs.  The view coordinator projects
     * its detailed counters and latches into this compact contract; tests can
     * then exercise the state graph without constructing a GL view or relying
     * on wall-clock timing.
     */
    struct Inputs {
	bool interactive = false;
	bool compacting = false;
	bool coverageActive = false;
	bool coverageComplete = false;
	bool generationActive = false;
	bool settlingWork = false;
    };

    struct Witness {
	uint64_t sequence = 0;
	uint64_t viewEpoch = 0;
	uint64_t policyEpoch = 0;
	uint64_t renderSerial = 0;
	uint64_t activeGeneration = 0;
	uint64_t residentDemandRevision = 0;
	size_t visibleCount = 0;
	size_t completedCount = 0;
	size_t pendingCount = 0;
    };

    BObolLodPhaseTracker(void) = default;

    static Phase phaseFor(const Inputs &inputs)
    {
	if (inputs.interactive)
	    return Phase::INTERACTIVE;
	if (inputs.compacting)
	    return Phase::COMPACTING;
	if (inputs.coverageActive || !inputs.coverageComplete)
	    return inputs.generationActive ? Phase::COVERAGE :
		Phase::FALLBACK;
	if (inputs.settlingWork)
	    return Phase::SETTLING;
	return Phase::STABLE;
    }

    Phase currentPhase(void) const
    {
	return this->phase;
    }

    uint64_t phaseTransitionSerial(void) const
    {
	return this->transitionSerial;
    }

    const Witness &phaseWitness(Phase requestedPhase) const
    {
	const size_t index = static_cast<size_t>(requestedPhase);
	return this->witnesses[index < this->witnesses.size() ?
	    index : static_cast<size_t>(Phase::FALLBACK)];
    }

    void enterFallback(const Witness &witness)
    {
	this->transitionTo(Phase::FALLBACK, witness);
    }

    void enterCoverage(const Witness &witness)
    {
	this->transitionTo(Phase::COVERAGE, witness);
    }

    void enterInteractive(const Witness &witness)
    {
	this->transitionTo(Phase::INTERACTIVE, witness);
    }

    void enterSettling(const Witness &witness)
    {
	this->transitionTo(Phase::SETTLING, witness);
    }

    void enterStable(const Witness &witness)
    {
	this->transitionTo(Phase::STABLE, witness);
    }

    void enterCompacting(const Witness &witness)
    {
	this->transitionTo(Phase::COMPACTING, witness);
    }

private:
    static bool sameProgress(const Witness &lhs, const Witness &rhs)
    {
	return lhs.viewEpoch == rhs.viewEpoch &&
	    lhs.policyEpoch == rhs.policyEpoch &&
	    lhs.renderSerial == rhs.renderSerial &&
	    lhs.activeGeneration == rhs.activeGeneration &&
	    lhs.residentDemandRevision == rhs.residentDemandRevision &&
	    lhs.visibleCount == rhs.visibleCount &&
	    lhs.completedCount == rhs.completedCount &&
	    lhs.pendingCount == rhs.pendingCount;
    }

    void transitionTo(Phase nextPhase, const Witness &nextWitness)
    {
	const size_t index = static_cast<size_t>(nextPhase);
	if (index >= this->witnesses.size())
	    return;

	if (this->phase != nextPhase) {
	    this->phase = nextPhase;
	    this->transitionSerial++;
	    if (this->transitionSerial == 0)
		this->transitionSerial++;
	}

	Witness &recorded = this->witnesses[index];
	if (!recorded.sequence ||
	    !BObolLodPhaseTracker::sameProgress(recorded, nextWitness)) {
	    recorded = nextWitness;
	    this->progressSerial++;
	    if (this->progressSerial == 0)
		this->progressSerial++;
	    recorded.sequence = this->progressSerial;
	}
    }

    Phase phase = Phase::FALLBACK;
    uint64_t transitionSerial = 0;
    uint64_t progressSerial = 0;
    std::array<Witness, static_cast<size_t>(Phase::COUNT)> witnesses {};
};

#endif /* LIBBOBOL_LOD_COORDINATOR_PRIVATE_H */
