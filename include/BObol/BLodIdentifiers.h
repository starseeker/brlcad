/*              B L O D I D E N T I F I E R S . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BLodIdentifiers.h */

#ifndef BOBOL_BLODIDENTIFIERS_H
#define BOBOL_BLODIDENTIFIERS_H

#include "BObol/BDefines.h"

#include <exception>
#include <limits>
#include <stdint.h>
#include <type_traits>

/**
 * Allocation-free strong value used for LoD identities and epochs.
 *
 * Raw integers may be assigned at API/serialization boundaries, but values
 * of different semantic types do not implicitly convert to one another.
 * Inner loops retain the exact layout and comparison cost of uint64_t.
 */
template <typename Tag>
class BObolLodStrongUInt64 {
public:
    constexpr BObolLodStrongUInt64(void) : rawValue(0) {}
    explicit constexpr BObolLodStrongUInt64(uint64_t value) :
	rawValue(value) {}

    BObolLodStrongUInt64 &operator=(uint64_t value)
    {
	this->rawValue = value;
	return *this;
    }

    template <typename OtherTag>
    BObolLodStrongUInt64 &operator=(
	const BObolLodStrongUInt64<OtherTag> &) = delete;

    explicit constexpr operator uint64_t(void) const
    {
	return this->rawValue;
    }

    constexpr uint64_t value(void) const
    {
	return this->rawValue;
    }

    void clear(void)
    {
	this->rawValue = 0;
    }

    void reset(void)
    {
	this->clear();
    }

    void set(uint64_t value)
    {
	this->rawValue = value;
    }

    void advance(void)
    {
	if (!this->canAdvance())
	    std::terminate();
	++this->rawValue;
    }

    constexpr bool canAdvance(void) const
    {
	return this->rawValue != std::numeric_limits<uint64_t>::max();
    }

    BObolLodStrongUInt64 &operator++(void)
    {
	this->advance();
	return *this;
    }

    BObolLodStrongUInt64 operator++(int)
    {
	BObolLodStrongUInt64 prior(*this);
	++(*this);
	return prior;
    }

    friend constexpr bool operator==(BObolLodStrongUInt64 lhs,
	    BObolLodStrongUInt64 rhs)
    {
	return lhs.rawValue == rhs.rawValue;
    }

    friend constexpr bool operator!=(BObolLodStrongUInt64 lhs,
	    BObolLodStrongUInt64 rhs)
    {
	return !(lhs == rhs);
    }

    friend constexpr bool operator<(BObolLodStrongUInt64 lhs,
	    BObolLodStrongUInt64 rhs)
    {
	return lhs.rawValue < rhs.rawValue;
    }

    friend constexpr bool operator>(BObolLodStrongUInt64 lhs,
	    BObolLodStrongUInt64 rhs)
    {
	return rhs < lhs;
    }

    friend constexpr bool operator<=(BObolLodStrongUInt64 lhs,
	    BObolLodStrongUInt64 rhs)
    {
	return !(rhs < lhs);
    }

    friend constexpr bool operator>=(BObolLodStrongUInt64 lhs,
	    BObolLodStrongUInt64 rhs)
    {
	return !(lhs < rhs);
    }

    /*
     * Integer comparisons are useful at serialization and sentinel
     * boundaries, but they must not make the strong value itself implicitly
     * convertible to an integer.
     */
    friend constexpr bool operator==(BObolLodStrongUInt64 lhs, uint64_t rhs)
    {
	return lhs.rawValue == rhs;
    }

    friend constexpr bool operator==(uint64_t lhs, BObolLodStrongUInt64 rhs)
    {
	return lhs == rhs.rawValue;
    }

    friend constexpr bool operator!=(BObolLodStrongUInt64 lhs, uint64_t rhs)
    {
	return !(lhs == rhs);
    }

    friend constexpr bool operator!=(uint64_t lhs, BObolLodStrongUInt64 rhs)
    {
	return !(lhs == rhs);
    }

    friend constexpr bool operator<(BObolLodStrongUInt64 lhs, uint64_t rhs)
    {
	return lhs.rawValue < rhs;
    }

    friend constexpr bool operator<(uint64_t lhs, BObolLodStrongUInt64 rhs)
    {
	return lhs < rhs.rawValue;
    }

    friend constexpr bool operator>(BObolLodStrongUInt64 lhs, uint64_t rhs)
    {
	return rhs < lhs;
    }

    friend constexpr bool operator>(uint64_t lhs, BObolLodStrongUInt64 rhs)
    {
	return rhs < lhs;
    }

    friend constexpr bool operator<=(BObolLodStrongUInt64 lhs, uint64_t rhs)
    {
	return !(rhs < lhs);
    }

    friend constexpr bool operator<=(uint64_t lhs, BObolLodStrongUInt64 rhs)
    {
	return !(rhs < lhs);
    }

    friend constexpr bool operator>=(BObolLodStrongUInt64 lhs, uint64_t rhs)
    {
	return !(lhs < rhs);
    }

    friend constexpr bool operator>=(uint64_t lhs, BObolLodStrongUInt64 rhs)
    {
	return !(lhs < rhs);
    }

    template <typename OtherTag>
    friend bool operator==(BObolLodStrongUInt64,
	    BObolLodStrongUInt64<OtherTag>) = delete;
    template <typename OtherTag>
    friend bool operator!=(BObolLodStrongUInt64,
	    BObolLodStrongUInt64<OtherTag>) = delete;
    template <typename OtherTag>
    friend bool operator<(BObolLodStrongUInt64,
	    BObolLodStrongUInt64<OtherTag>) = delete;
    template <typename OtherTag>
    friend bool operator>(BObolLodStrongUInt64,
	    BObolLodStrongUInt64<OtherTag>) = delete;
    template <typename OtherTag>
    friend bool operator<=(BObolLodStrongUInt64,
	    BObolLodStrongUInt64<OtherTag>) = delete;
    template <typename OtherTag>
    friend bool operator>=(BObolLodStrongUInt64,
	    BObolLodStrongUInt64<OtherTag>) = delete;

private:
    uint64_t rawValue;
};

struct BObolDatabaseEpochTag;
struct BObolSourceEpochTag;
struct BObolViewEpochTag;
struct BObolPolicyEpochTag;
struct BObolInventoryEpochTag;
struct BObolSourceRoutingIdTag;
struct BObolSourcePopulationEpochTag;
struct BObolOccurrenceIdTag;
struct BObolAssetIdTag;

using BObolDatabaseEpoch = BObolLodStrongUInt64<BObolDatabaseEpochTag>;
using BObolSourceEpoch = BObolLodStrongUInt64<BObolSourceEpochTag>;
using BObolViewEpoch = BObolLodStrongUInt64<BObolViewEpochTag>;
using BObolPolicyEpoch = BObolLodStrongUInt64<BObolPolicyEpochTag>;
using BObolInventoryEpoch = BObolLodStrongUInt64<BObolInventoryEpochTag>;
using BObolSourceRoutingId =
    BObolLodStrongUInt64<BObolSourceRoutingIdTag>;
using BObolSourcePopulationEpoch =
    BObolLodStrongUInt64<BObolSourcePopulationEpochTag>;
using BObolOccurrenceId = BObolLodStrongUInt64<BObolOccurrenceIdTag>;
using BObolAssetId = BObolLodStrongUInt64<BObolAssetIdTag>;

static_assert(sizeof(BObolSourceEpoch) == sizeof(uint64_t),
    "strong LoD identities must remain zero-overhead");
static_assert(std::is_trivially_copyable<BObolSourceEpoch>::value,
    "strong LoD identities must remain trivially copyable");
static_assert(!std::is_convertible<BObolSourceEpoch, uint64_t>::value,
    "strong LoD identities must not decay implicitly to integers");
static_assert(!std::is_assignable<BObolViewEpoch &, BObolPolicyEpoch>::value,
    "different LoD epoch domains must not be assignable");
static_assert(!std::is_convertible<BObolViewEpoch, BObolPolicyEpoch>::value,
    "different LoD epoch domains must not be convertible");
static_assert(sizeof(BObolSourcePopulationEpoch) == sizeof(uint64_t),
    "compact population epochs must remain fixed-width");

#endif /* BOBOL_BLODIDENTIFIERS_H */
