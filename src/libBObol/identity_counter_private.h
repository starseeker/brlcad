/*       I D E N T I T Y _ C O U N T E R _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_IDENTITY_COUNTER_PRIVATE_H
#define LIBBOBOL_IDENTITY_COUNTER_PRIVATE_H

#include "common.h"

#include <atomic>
#include <exception>
#include <limits>
#include <type_traits>

/* Identity exhaustion must never reuse an earlier value.  There is no safe
 * local recovery once an externally visible asynchronous credential is out of
 * values: continuing could authenticate arbitrary retired work.  Termination
 * is therefore the deliberate fail-stop boundary.  The pure successor helper
 * makes the machine-integer edge directly testable without terminating a
 * test process. */
template <typename UInt>
constexpr bool
bobol_identity_successor(UInt current, UInt &successor) noexcept
{
    static_assert(std::is_integral<UInt>::value &&
	std::is_unsigned<UInt>::value,
	"identity counters require an unsigned integral representation");
    if (current == std::numeric_limits<UInt>::max())
	return false;
    successor = static_cast<UInt>(current + 1);
    return true;
}

[[noreturn]] inline void
bobol_identity_exhausted(void) noexcept
{
    std::terminate();
}

template <typename UInt>
UInt
bobol_identity_successor_or_terminate(UInt current)
{
    UInt successor = current;
    if (!bobol_identity_successor(current, successor))
	bobol_identity_exhausted();
    return successor;
}

template <typename UInt>
void
bobol_identity_advance(UInt &identity)
{
    identity = bobol_identity_successor_or_terminate(identity);
}

/* Allocate the current value and advance the source.  UINT_MAX is reserved as
 * the exhausted next-value state and is never issued by this operation. */
template <typename UInt>
UInt
bobol_identity_take(UInt &nextIdentity)
{
    const UInt identity = nextIdentity;
    bobol_identity_advance(nextIdentity);
    return identity;
}

template <typename UInt>
UInt
bobol_nonzero_identity_take(UInt &nextIdentity)
{
    if (nextIdentity == 0)
	bobol_identity_exhausted();
    return bobol_identity_take(nextIdentity);
}

template <typename UInt>
UInt
bobol_atomic_identity_advance(std::atomic<UInt> &identity,
    std::memory_order successOrder = std::memory_order_relaxed,
    std::memory_order failureOrder = std::memory_order_relaxed)
{
    UInt observed = identity.load(std::memory_order_relaxed);
    for (;;) {
	const UInt successor =
	    bobol_identity_successor_or_terminate(observed);
	if (identity.compare_exchange_weak(observed, successor,
		successOrder, failureOrder))
	    return successor;
    }
}

template <typename UInt>
UInt
bobol_atomic_nonzero_identity_take(std::atomic<UInt> &nextIdentity,
    std::memory_order successOrder = std::memory_order_relaxed,
    std::memory_order failureOrder = std::memory_order_relaxed)
{
    UInt observed = nextIdentity.load(std::memory_order_relaxed);
    for (;;) {
	if (observed == 0)
	    bobol_identity_exhausted();
	const UInt successor =
	    bobol_identity_successor_or_terminate(observed);
	if (nextIdentity.compare_exchange_weak(observed, successor,
		successOrder, failureOrder))
	    return observed;
    }
}

/* Diagnostic totals do not authenticate work and need not fail-stop at their
 * representation boundary.  Saturation preserves their only useful promise:
 * once nonzero, they never misleadingly return to zero. */
template <typename UInt>
void
bobol_saturating_counter_advance(UInt &counter) noexcept
{
    static_assert(std::is_integral<UInt>::value &&
	std::is_unsigned<UInt>::value,
	"diagnostic counters require an unsigned integral representation");
    if (counter != std::numeric_limits<UInt>::max())
	++counter;
}

#endif /* LIBBOBOL_IDENTITY_COUNTER_PRIVATE_H */
