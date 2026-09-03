/*       E X A C T _ S A M P L E _ I D E N T I T Y _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_EXACT_SAMPLE_IDENTITY_PRIVATE_H
#define LIBBOBOL_EXACT_SAMPLE_IDENTITY_PRIVATE_H

#include "identity_counter_private.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

/*
 * Assign a non-reused token to the exact unordered set of source samples.
 * Only the current set is retained: returning from A through B to A is a new
 * observation and therefore receives a new token.  Source identities must
 * themselves be non-reused; an object address does not satisfy that contract.
 */
class BObolExactSampleIdentity {
public:
    struct Entry {
	uint64_t sourceIdentity = 0;
	uint64_t sourceSerial = 0;

	bool operator==(const Entry &other) const noexcept
	{
	    return this->sourceIdentity == other.sourceIdentity &&
		this->sourceSerial == other.sourceSerial;
	}
    };

    using Entries = std::vector<Entry>;

    explicit BObolExactSampleIdentity(
	uint64_t maximumIdentity = std::numeric_limits<uint64_t>::max()) :
	maximumIdentityValue(maximumIdentity)
    {
	if (!this->maximumIdentityValue)
	    bobol_identity_exhausted();
    }

    uint64_t intern(Entries entries)
    {
	std::sort(entries.begin(), entries.end(),
	    [](const Entry &left, const Entry &right) {
		return left.sourceIdentity < right.sourceIdentity ||
		    (left.sourceIdentity == right.sourceIdentity &&
		     left.sourceSerial < right.sourceSerial);
	    });
	if (this->currentIdentity && entries == this->currentEntries)
	    return this->currentIdentity;
	if (this->nextIdentity > this->maximumIdentityValue)
	    bobol_identity_exhausted();
	this->currentIdentity =
	    bobol_nonzero_identity_take(this->nextIdentity);
	this->currentEntries = std::move(entries);
	return this->currentIdentity;
    }

    void invalidate() noexcept
    {
	this->currentEntries.clear();
	this->currentIdentity = 0;
    }

private:
    const uint64_t maximumIdentityValue;
    uint64_t nextIdentity = 1;
    uint64_t currentIdentity = 0;
    Entries currentEntries;
};

#endif /* LIBBOBOL_EXACT_SAMPLE_IDENTITY_PRIVATE_H */
