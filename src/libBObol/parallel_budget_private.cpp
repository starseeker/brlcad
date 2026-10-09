/*           P A R A L L E L _ B U D G E T _ P R I V A T E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "parallel_budget_private.h"

#include "bu/parallel.h"

#include <algorithm>
#include <condition_variable>
#include <cstdint>
#include <mutex>
#include <vector>

namespace {

struct BObolParallelBudgetState {
    BObolParallelBudgetState(void)
    {
	const int available = bu_avail_cpus();
	/* The display service itself is capped at eight workers.  Applying the
	 * same ceiling here bounds both explicitly larger services and nested
	 * cache helpers without changing the established default on ordinary
	 * workstations. */
	limit = std::min<size_t>(8,
	    static_cast<size_t>(std::max(1, available)));
    }

    std::mutex mutex;
    std::condition_variable available;
    size_t limit = 1;
    size_t outerSlots = 0;
    size_t helperSlots = 0;
    struct OuterWaiter {
	int priority = 0;
	uint64_t ticket = 0;
    };
    std::vector<OuterWaiter *> outerWaiters;
    uint64_t nextOuterTicket = 1;
};

bool
parallel_budget_waiter_is_next(const BObolParallelBudgetState &state,
    const BObolParallelBudgetState::OuterWaiter *candidate)
{
    if (!candidate)
	return false;
    for (const BObolParallelBudgetState::OuterWaiter *other :
	 state.outerWaiters) {
	if (!other || other == candidate)
	    continue;
	if (other->priority > candidate->priority ||
	    (other->priority == candidate->priority &&
	     other->ticket < candidate->ticket))
	    return false;
    }
    return true;
}

BObolParallelBudgetState &
parallel_budget_state(void)
{
    /* Process-lifetime storage avoids static-destruction ordering between
     * independently owned services and cache handles during DLL shutdown. */
    static BObolParallelBudgetState *state =
	new BObolParallelBudgetState();
    return *state;
}

} /* namespace */

BObolParallelBudgetLease::~BObolParallelBudgetLease(void)
{
    this->release();
}

void
BObolParallelBudgetLease::acquireOuter(int priority)
{
    if (this->slotCount)
	return;
    BObolParallelBudgetState &state = parallel_budget_state();
    std::unique_lock<std::mutex> lock(state.mutex);
    BObolParallelBudgetState::OuterWaiter waiter;
    waiter.priority = priority;
    waiter.ticket = state.nextOuterTicket++;
    if (!state.nextOuterTicket)
	state.nextOuterTicket = 1;
    state.outerWaiters.push_back(&waiter);
    try {
	state.available.wait(lock, [&state, &waiter]() {
	    return state.outerSlots + state.helperSlots < state.limit &&
		parallel_budget_waiter_is_next(state, &waiter);
	});
    } catch (...) {
	const auto failed = std::find(state.outerWaiters.begin(),
	    state.outerWaiters.end(), &waiter);
	if (failed != state.outerWaiters.end())
	    state.outerWaiters.erase(failed);
	state.available.notify_all();
	throw;
    }
    const auto found = std::find(state.outerWaiters.begin(),
	state.outerWaiters.end(), &waiter);
    if (found != state.outerWaiters.end())
	state.outerWaiters.erase(found);
    state.outerSlots++;
    this->slotCount = 1;
    this->role = Role::OUTER;
}

size_t
BObolParallelBudgetLease::tryAcquireHelpers(size_t maximum)
{
    if (!maximum || this->slotCount)
	return 0;
    BObolParallelBudgetState &state = parallel_budget_state();
    std::lock_guard<std::mutex> lock(state.mutex);
    const size_t occupied = state.outerSlots + state.helperSlots;
    const size_t available = occupied < state.limit ?
	state.limit - occupied : 0;

    /* Keep two top-level lanes available when service work exists.  One lane
     * can construct the visually dominant asset while the other obtains a
     * fast first-mesh result; helpers may consume all capacity for a direct,
     * non-service cache operation. */
    const size_t desiredOuterSlots =
	(state.outerSlots || !state.outerWaiters.empty()) ?
	    std::min<size_t>(2, state.limit) : 0;
    const size_t reservedOuterSlots = desiredOuterSlots > state.outerSlots ?
	desiredOuterSlots - state.outerSlots : 0;
    const size_t borrowable = available > reservedOuterSlots ?
	available - reservedOuterSlots : 0;
    this->slotCount = std::min(maximum, borrowable);
    if (!this->slotCount)
	return 0;
    state.helperSlots += this->slotCount;
    this->role = Role::HELPER;
    return this->slotCount;
}

void
BObolParallelBudgetLease::release(void)
{
    if (!this->slotCount)
	return;
    BObolParallelBudgetState &state = parallel_budget_state();
    {
	std::lock_guard<std::mutex> lock(state.mutex);
	if (this->role == Role::OUTER) {
	    state.outerSlots = this->slotCount >= state.outerSlots ?
		0 : state.outerSlots - this->slotCount;
	} else if (this->role == Role::HELPER) {
	    state.helperSlots = this->slotCount >= state.helperSlots ?
		0 : state.helperSlots - this->slotCount;
	}
	this->slotCount = 0;
	this->role = Role::NONE;
    }
    state.available.notify_all();
}

size_t
bobol_parallel_budget_limit(void)
{
    return parallel_budget_state().limit;
}

size_t
bobol_parallel_budget_waiting_outer_count(void)
{
    BObolParallelBudgetState &state = parallel_budget_state();
    std::lock_guard<std::mutex> lock(state.mutex);
    return state.outerWaiters.size();
}
