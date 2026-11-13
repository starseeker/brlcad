/*    L O D _ R E S U L T _ A U T H E N T I C A T I O N _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_RESULT_AUTHENTICATION_PRIVATE_H
#define LIBBOBOL_LOD_RESULT_AUTHENTICATION_PRIVATE_H

#include "BObol/BLodRealization.h"

#include <cstdint>
#include <type_traits>

enum class BObolLodResultIdentityDomain : uint8_t {
    SOURCE_ROUTE = 0,
    SOURCE_POPULATION,
    DEMAND
};

enum class BObolLodResultDisposition : uint8_t {
    PUBLISH = 0,
    RECORD_TERMINAL_FAILURE,
    RETRY_CURRENT_DEMAND,
    SUPERSEDE
};

struct BObolLodResultAuthenticationContext {
    BObolSourceRoutingId sourceRoutingId;
    BObolSourcePopulationEpoch sourcePopulationEpoch;
    BObolViewEpoch viewRevision;
    BObolPolicyEpoch policyRevision;
};

struct BObolLodResultAuthentication {
    uint8_t identityMismatchMask = 0;
    BObolLodResultDisposition disposition =
	BObolLodResultDisposition::RETRY_CURRENT_DEMAND;

    bool identityCurrent(void) const
    {
	return this->identityMismatchMask == 0;
    }

    bool mismatches(BObolLodResultIdentityDomain domain) const
    {
	return (this->identityMismatchMask & mismatchBit(domain)) != 0;
    }

    static constexpr uint8_t mismatchBit(
	BObolLodResultIdentityDomain domain)
    {
	return static_cast<uint8_t>(1u << static_cast<uint8_t>(domain));
    }
};

/**
 * Allocation-free refinement of ObolResultAuthentication.tla.
 *
 * Compact result identity has three independent domains: source route,
 * compact population, and view/policy demand.  Completion status is examined
 * only after all three are current, so a terminal error for obsolete work can
 * never become failure evidence for its successor.
 */
class BObolLodResultAuthenticationContract {
public:
    static BObolLodResultAuthentication evaluate(
	const BObolLodRequest &request,
	int providerStatus,
	const BObolLodResultAuthenticationContext &current)
    {
	BObolLodResultAuthentication authentication;
	const bool compactResult = request.occurrenceKey.getLength() > 0;
	if (compactResult &&
	    (request.sourceRoutingId == 0 || current.sourceRoutingId == 0 ||
	     request.sourceRoutingId != current.sourceRoutingId))
	    authentication.identityMismatchMask |=
		BObolLodResultAuthentication::mismatchBit(
		    BObolLodResultIdentityDomain::SOURCE_ROUTE);
	if (compactResult &&
	    (request.sourcePopulationEpoch == 0 ||
	     current.sourcePopulationEpoch == 0 ||
	     request.sourcePopulationEpoch != current.sourcePopulationEpoch))
	    authentication.identityMismatchMask |=
		BObolLodResultAuthentication::mismatchBit(
		    BObolLodResultIdentityDomain::SOURCE_POPULATION);
	if (request.viewRevision != current.viewRevision ||
	    request.policyRevision != current.policyRevision)
	    authentication.identityMismatchMask |=
		BObolLodResultAuthentication::mismatchBit(
		    BObolLodResultIdentityDomain::DEMAND);

	if (!authentication.identityCurrent()) {
	    authentication.disposition = BObolLodResultDisposition::SUPERSEDE;
	    return authentication;
	}

	switch (providerStatus) {
	    case BOBOL_LOD_PROVIDER_READY:
		authentication.disposition =
		    BObolLodResultDisposition::PUBLISH;
		break;
	    case BOBOL_LOD_PROVIDER_CACHE_MISS:
	    case BOBOL_LOD_PROVIDER_STALE:
	    case BOBOL_LOD_PROVIDER_ERROR:
		authentication.disposition =
		    BObolLodResultDisposition::RECORD_TERMINAL_FAILURE;
		break;
	    default:
		authentication.disposition =
		    BObolLodResultDisposition::RETRY_CURRENT_DEMAND;
		break;
	}
	return authentication;
    }
};

static_assert(std::is_trivially_copyable<
	BObolLodResultAuthenticationContext>::value,
    "result authentication context must remain allocation free");
static_assert(std::is_trivially_copyable<
	BObolLodResultAuthentication>::value,
    "result authentication must remain allocation free");

#endif /* LIBBOBOL_LOD_RESULT_AUTHENTICATION_PRIVATE_H */
