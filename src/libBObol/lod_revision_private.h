/*          L O D _ R E V I S I O N _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_REVISION_PRIVATE_H
#define LIBBOBOL_LOD_REVISION_PRIVATE_H

#include "common.h"
#include "BObol/BLodIdentifiers.h"

#include <cstdint>
#include <type_traits>

/* Use the same strong epoch domains as requests and results.  Keeping a
 * second coordinator-only implementation previously allowed the state
 * machine and service boundary to drift while appearing type safe. */
using BObolLodViewEpoch = BObolViewEpoch;
using BObolLodPolicyEpoch = BObolPolicyEpoch;
using BObolLodInventoryEpoch = BObolInventoryEpoch;
using BObolLodSourceRoutingId = BObolSourceRoutingId;

struct BObolLodAvailabilityEpochTag;
struct BObolLodVisibilityEpochTag;
struct BObolLodCapacityEpochTag;
using BObolLodAvailabilityEpoch =
    BObolLodStrongUInt64<BObolLodAvailabilityEpochTag>;
using BObolLodVisibilityEpoch =
    BObolLodStrongUInt64<BObolLodVisibilityEpochTag>;
using BObolLodCapacityEpoch =
    BObolLodStrongUInt64<BObolLodCapacityEpochTag>;

enum class BObolLodAdmissionRevisionDomain : uint8_t {
    INVENTORY = 0,
    AVAILABILITY,
    VISIBILITY,
    VIEW,
    POLICY,
    CAPACITY
};

/* Exact inputs which certify one admission plan.  Cardinality, cache warmth,
 * renderer backend, and model identity influence the ledgers behind these
 * epochs; they do not add control modes or appear in the cursor itself. */
class BObolLodAdmissionRevisionStamp {
public:
    constexpr BObolLodAdmissionRevisionStamp(
	BObolLodInventoryEpoch inventory,
	BObolLodAvailabilityEpoch availability,
	BObolLodVisibilityEpoch visibility, BObolLodViewEpoch view,
	BObolLodPolicyEpoch policy, BObolLodCapacityEpoch capacity) :
	inventoryValue(inventory), availabilityValue(availability),
	visibilityValue(visibility), viewValue(view), policyValue(policy),
	capacityValue(capacity)
    {
    }

    /* A zero stamp is an explicit administrative marker, never a partially
     * assembled certificate. */
    static constexpr BObolLodAdmissionRevisionStamp administrative(void)
    {
	return BObolLodAdmissionRevisionStamp(BObolLodInventoryEpoch(0),
	    BObolLodAvailabilityEpoch(0), BObolLodVisibilityEpoch(0),
	    BObolLodViewEpoch(0), BObolLodPolicyEpoch(0),
	    BObolLodCapacityEpoch(0));
    }

    constexpr BObolLodInventoryEpoch inventory(void) const
    {
	return this->inventoryValue;
    }

    constexpr BObolLodAvailabilityEpoch availability(void) const
    {
	return this->availabilityValue;
    }

    constexpr BObolLodVisibilityEpoch visibility(void) const
    {
	return this->visibilityValue;
    }

    constexpr BObolLodViewEpoch view(void) const
    {
	return this->viewValue;
    }

    constexpr BObolLodPolicyEpoch policy(void) const
    {
	return this->policyValue;
    }

    constexpr BObolLodCapacityEpoch capacity(void) const
    {
	return this->capacityValue;
    }

    bool same(const BObolLodAdmissionRevisionStamp &other) const
    {
	return this->inventoryValue == other.inventoryValue &&
	    this->availabilityValue == other.availabilityValue &&
	    this->visibilityValue == other.visibilityValue &&
	    this->viewValue == other.viewValue &&
	    this->policyValue == other.policyValue &&
	    this->capacityValue == other.capacityValue;
    }

    /* A renderer-capacity measurement is reusable across a visibility-only
     * reallocation and across the capacity revision produced by that very
     * measurement.  Every other domain changes the measured problem.  Keep
     * this intentional projection here rather than open-coding partial stamp
     * comparisons at certificate consumers. */
    bool sameRendererCapacityProblem(
	const BObolLodAdmissionRevisionStamp &other) const
    {
	return this->inventoryValue == other.inventoryValue &&
	    this->availabilityValue == other.availabilityValue &&
	    this->viewValue == other.viewValue &&
	    this->policyValue == other.policyValue;
    }

    bool empty(void) const
    {
	return this->inventoryValue == 0 && this->availabilityValue == 0 &&
	    this->visibilityValue == 0 && this->viewValue == 0 &&
	    this->policyValue == 0 && this->capacityValue == 0;
    }

private:
    BObolLodInventoryEpoch inventoryValue;
    BObolLodAvailabilityEpoch availabilityValue;
    BObolLodVisibilityEpoch visibilityValue;
    BObolLodViewEpoch viewValue;
    BObolLodPolicyEpoch policyValue;
    BObolLodCapacityEpoch capacityValue;
};

static_assert(!std::is_default_constructible<
	BObolLodAdmissionRevisionStamp>::value,
    "admission revision stamps must name all six domains");
static_assert(std::is_trivially_copyable<
	BObolLodAdmissionRevisionStamp>::value,
    "admission revision stamps must remain allocation-free values");

/* Executable refinement boundary for the six-domain formal revision tuple.
 * Production revision owners use this pure value transition, and focused
 * tests exercise it without constructing a Coin controller. */
class BObolLodRevisionContract {
public:
    enum class PlanDisposition : uint8_t {
	ADMINISTRATIVE = 0,
	CURRENT,
	STALE
    };

    static BObolLodAdmissionRevisionStamp advance(
	const BObolLodAdmissionRevisionStamp &current,
	BObolLodAdmissionRevisionDomain domain)
    {
	BObolLodInventoryEpoch inventory = current.inventory();
	BObolLodAvailabilityEpoch availability = current.availability();
	BObolLodVisibilityEpoch visibility = current.visibility();
	BObolLodViewEpoch view = current.view();
	BObolLodPolicyEpoch policy = current.policy();
	BObolLodCapacityEpoch capacity = current.capacity();
	switch (domain) {
	    case BObolLodAdmissionRevisionDomain::INVENTORY:
		inventory.advance();
		break;
	    case BObolLodAdmissionRevisionDomain::AVAILABILITY:
		availability.advance();
		break;
	    case BObolLodAdmissionRevisionDomain::VISIBILITY:
		visibility.advance();
		break;
	    case BObolLodAdmissionRevisionDomain::VIEW:
		view.advance();
		break;
	    case BObolLodAdmissionRevisionDomain::POLICY:
		policy.advance();
		break;
	    case BObolLodAdmissionRevisionDomain::CAPACITY:
		capacity.advance();
		break;
	}
	return BObolLodAdmissionRevisionStamp(inventory, availability,
	    visibility, view, policy, capacity);
    }

    static BObolLodAdmissionRevisionStamp setPolicy(
	const BObolLodAdmissionRevisionStamp &current, uint64_t revision)
    {
	return BObolLodAdmissionRevisionStamp(current.inventory(),
	    current.availability(), current.visibility(), current.view(),
	    BObolLodPolicyEpoch(revision), current.capacity());
    }

    static PlanDisposition planDisposition(
	const BObolLodAdmissionRevisionStamp &plan,
	const BObolLodAdmissionRevisionStamp &current)
    {
	if (plan.empty())
	    return PlanDisposition::ADMINISTRATIVE;
	return plan.same(current) ? PlanDisposition::CURRENT :
	    PlanDisposition::STALE;
    }
};

#endif /* LIBBOBOL_LOD_REVISION_PRIVATE_H */
