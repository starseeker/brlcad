/*             P A R A L L E L _ B U D G E T _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_PARALLEL_BUDGET_PRIVATE_H
#define LIBBOBOL_PARALLEL_BUDGET_PRIVATE_H

#include "common.h"

#include <cstddef>

/*
 * Process-wide CPU admission for background LoD work.
 *
 * A service callback owns one outer slot while it can execute CPU work.
 * Synchronous helpers may opportunistically borrow the remaining slots, but
 * never wait for them: falling back to the caller prevents nested parallelism
 * from deadlocking.  Helper admission leaves room for a second outer producer
 * whenever the machine has enough capacity, preserving the cold-preview
 * hero/quick-win lanes.
 */
class BObolParallelBudgetLease {
public:
    BObolParallelBudgetLease(void) = default;
    ~BObolParallelBudgetLease(void);

    BObolParallelBudgetLease(const BObolParallelBudgetLease &) = delete;
    BObolParallelBudgetLease &operator=(
	const BObolParallelBudgetLease &) = delete;

    /* Wait for one top-level background-work slot.  Larger priorities pass
     * already-waiting lower-priority work; equal priorities remain FIFO.
     * Running work is never preempted. */
    void acquireOuter(int priority = 0);

    /* Borrow up to maximum helper slots without waiting. */
    size_t tryAcquireHelpers(size_t maximum);

    void release(void);
    size_t size(void) const { return this->slotCount; }

private:
    enum class Role {
	NONE,
	OUTER,
	HELPER
    };

    size_t slotCount = 0;
    Role role = Role::NONE;
};

size_t bobol_parallel_budget_limit(void);

/* Lock-consistent test/diagnostic observation of queued outer leases. */
size_t bobol_parallel_budget_waiting_outer_count(void);

#endif /* LIBBOBOL_PARALLEL_BUDGET_PRIVATE_H */
