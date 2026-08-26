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
struct BObolLodCapacityEpochTag;
using BObolLodAvailabilityEpoch =
    BObolLodStrongUInt64<BObolLodAvailabilityEpochTag>;
using BObolLodCapacityEpoch =
    BObolLodStrongUInt64<BObolLodCapacityEpochTag>;

enum class BObolLodAdmissionRevisionDomain : uint8_t {
    INVENTORY = 0,
    AVAILABILITY,
    VIEW,
    POLICY,
    CAPACITY
};

/* Exact inputs which certify one admission plan.  Cardinality, cache warmth,
 * renderer backend, and model identity influence the ledgers behind these
 * epochs; they do not add control modes or appear in the cursor itself. */
struct BObolLodAdmissionRevisionStamp {
    BObolLodInventoryEpoch inventory;
    BObolLodAvailabilityEpoch availability;
    BObolLodViewEpoch view;
    BObolLodPolicyEpoch policy;
    BObolLodCapacityEpoch capacity;

    bool same(const BObolLodAdmissionRevisionStamp &other) const
    {
	return this->inventory == other.inventory &&
	    this->availability == other.availability &&
	    this->view == other.view &&
	    this->policy == other.policy &&
	    this->capacity == other.capacity;
    }

    bool empty(void) const
    {
	return this->inventory == 0 && this->availability == 0 &&
	    this->view == 0 && this->policy == 0 && this->capacity == 0;
    }
};

static_assert(std::is_trivially_copyable<
	BObolLodAdmissionRevisionStamp>::value,
    "admission revision stamps must remain allocation-free values");

/* Executable refinement boundary for the five-domain formal revision tuple.
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
	BObolLodAdmissionRevisionStamp next = current;
	switch (domain) {
	    case BObolLodAdmissionRevisionDomain::INVENTORY:
		next.inventory.advance();
		break;
	    case BObolLodAdmissionRevisionDomain::AVAILABILITY:
		next.availability.advance();
		break;
	    case BObolLodAdmissionRevisionDomain::VIEW:
		next.view.advance();
		break;
	    case BObolLodAdmissionRevisionDomain::POLICY:
		next.policy.advance();
		break;
	    case BObolLodAdmissionRevisionDomain::CAPACITY:
		next.capacity.advance();
		break;
	}
	return next;
    }

    static BObolLodAdmissionRevisionStamp setPolicy(
	const BObolLodAdmissionRevisionStamp &current, uint64_t revision)
    {
	BObolLodAdmissionRevisionStamp next = current;
	next.policy.set(revision);
	return next;
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
