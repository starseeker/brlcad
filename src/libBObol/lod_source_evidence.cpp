/*          L O D _ S O U R C E _ E V I D E N C E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"
#include "lod_source_evidence_private.h"

#include <algorithm>
#include <utility>

namespace {

constexpr size_t PERFORMANCE_EPOCH_RETAINED_DIVISOR = 2;

}

bool
BObolLodSourceSnapshot::sameIdentity(const BObolLodSourceSnapshot &other) const
{
    return this->database == other.database &&
	this->routingId == other.routingId &&
	this->databaseId == other.databaseId &&
	this->path == other.path &&
	this->drawMode == other.drawMode &&
	this->representationMode == other.representationMode &&
	this->visible == other.visible &&
	this->lodBotThreshold == other.lodBotThreshold &&
	this->sourceRevision == other.sourceRevision &&
	this->inputsRevision == other.inputsRevision;
}

BObolLodSourceSnapshotSet::BObolLodSourceSnapshotSet(Sources sources)
{
    if (!sources.empty())
	this->sourceValues = std::make_shared<const Sources>(std::move(sources));
}

const BObolLodSourceSnapshotSet::Sources &
BObolLodSourceSnapshotSet::sources(void) const
{
    static const Sources emptySources;
    return this->sourceValues ? *this->sourceValues : emptySources;
}

bool
BObolLodSourceSnapshotSet::sameIdentities(const Sources &other) const
{
    const Sources &current = this->sources();
    if (current.size() != other.size())
	return false;
    for (size_t i = 0; i < current.size(); ++i) {
	if (!current[i].sameIdentity(other[i]))
	    return false;
    }
    return true;
}

bool
BObolLodSourceSnapshotSet::sameInventories(const Sources &other) const
{
    if (!this->sameIdentities(other))
	return false;
    const Sources &current = this->sources();
    for (size_t i = 0; i < current.size(); ++i) {
	if (current[i].inventoryRevision != other[i].inventoryRevision)
	    return false;
    }
    return true;
}

bool
BObolLodSourceSnapshotSet::sameVisibilityInputs(const Sources &other) const
{
    if (!this->sameIdentities(other))
	return false;
    const Sources &current = this->sources();
    for (size_t i = 0; i < current.size(); ++i) {
	if (current[i].visibilityRevision != other[i].visibilityRevision)
	    return false;
    }
    return true;
}

bool
BObolLodSourceSnapshotSet::samePlanningInputs(const Sources &other) const
{
    return this->sameInventories(other) && this->sameVisibilityInputs(other);
}

bool
BObolLodSourceSnapshotSet::startsNewPerformanceEpoch(
    const Sources &previous) const
{
    const Sources &current = this->sources();
    if (current.empty() || previous.empty())
	return current.empty() != previous.empty();

    size_t unchangedCount = 0;
    std::vector<bool> matched(previous.size(), false);
    for (const BObolLodSourceSnapshot &source : current) {
	for (size_t i = 0; i < previous.size(); ++i) {
	    if (!matched[i] && source.sameIdentity(previous[i])) {
		matched[i] = true;
		unchangedCount++;
		break;
	    }
	}
    }

    const size_t largerCount = std::max(current.size(), previous.size());
    /* Preserve calibration only when at least half of the larger source set
     * survives unchanged.  Exact render cost still adapts ordinary complexity
     * changes; this threshold detects a different workload composition. */
    const size_t minimumUnchangedCount =
	largerCount / PERFORMANCE_EPOCH_RETAINED_DIVISOR +
	largerCount % PERFORMANCE_EPOCH_RETAINED_DIVISOR;
    return unchangedCount < minimumUnchangedCount;
}

std::optional<BObolLodAdmissionRevisionDomain>
BObolLodSourceEvidence::observe(const BObolLodSourceSnapshotSet &current)
{
    const bool sameInventories =
	current.sameInventories(this->observedValue.sources());
    if (sameInventories &&
	current.sameVisibilityInputs(this->observedValue.sources()))
	return std::nullopt;
    const BObolLodAdmissionRevisionDomain domain =
	sameInventories ?
	BObolLodAdmissionRevisionDomain::VISIBILITY :
	BObolLodAdmissionRevisionDomain::INVENTORY;
    this->observedValue = current;
    return domain;
}

bool
BObolLodSourceEvidence::consume(const BObolLodSourceSnapshotSet &current)
{
    if (!current.samePlanningInputs(this->observedValue.sources()))
	return false;
    this->submittedValue = this->observedValue;
    return true;
}
