/*          L O D _ S O U R C E _ E V I D E N C E _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_SOURCE_EVIDENCE_PRIVATE_H
#define LIBBOBOL_LOD_SOURCE_EVIDENCE_PRIVATE_H

#include "common.h"
#include "lod_revision_private.h"

#include <Inventor/SbString.h>

#include <memory>
#include <optional>
#include <vector>

class SoBRLDatabaseSource;
struct db_i;

struct BObolLodSourceSnapshot {
    SoBRLDatabaseSource *source = NULL;
    struct db_i *database = NULL;
    BObolLodSourceRoutingId routingId;
    BObolLodInventoryEpoch inventoryRevision;
    uint64_t visibilityRevision = 0;
    SbString databaseId;
    SbString path;
    int drawMode = 0;
    int representationMode = 0;
    SbBool visible = FALSE;
    int lodBotThreshold = 0;
    uint32_t sourceRevision = 0;
    uint32_t inputsRevision = 0;

    bool sameIdentity(const BObolLodSourceSnapshot &other) const;
};

/* Source-contract signatures contain no occurrence records or geometry.
 * Observed and submitted evidence share this immutable storage once the
 * bounded submission catches up; intermediate observations replace only the
 * observed handle, preserving the submitted baseline for sparse deltas. */
class BObolLodSourceSnapshotSet {
public:
    using Sources = std::vector<BObolLodSourceSnapshot>;

    BObolLodSourceSnapshotSet() = default;
    explicit BObolLodSourceSnapshotSet(Sources sources);

    const Sources &sources(void) const;
    bool empty(void) const { return this->sources().empty(); }
    bool sameIdentities(const Sources &other) const;
    bool sameInventories(const Sources &other) const;
    bool sameVisibilityInputs(const Sources &other) const;
    bool samePlanningInputs(const Sources &other) const;

private:
    std::shared_ptr<const Sources> sourceValues;
};

class BObolLodSourceEvidence {
public:
    std::optional<BObolLodAdmissionRevisionDomain> observe(
	const BObolLodSourceSnapshotSet &current);

    /* A reentrant source notification may supersede the snapshot whose
     * submission just returned.  Such a pass cannot consume newer evidence. */
    bool consume(const BObolLodSourceSnapshotSet &current);

    const BObolLodSourceSnapshotSet::Sources &submitted(void) const
    {
	return this->submittedValue.sources();
    }

    const BObolLodSourceSnapshotSet::Sources &observed(void) const
    {
	return this->observedValue.sources();
    }

    void reset(void)
    {
	this->observedValue = BObolLodSourceSnapshotSet();
	this->submittedValue = BObolLodSourceSnapshotSet();
    }

private:
    BObolLodSourceSnapshotSet observedValue;
    BObolLodSourceSnapshotSet submittedValue;
};

#endif /* LIBBOBOL_LOD_SOURCE_EVIDENCE_PRIVATE_H */
