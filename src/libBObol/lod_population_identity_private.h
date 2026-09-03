/*      L O D _ P O P U L A T I O N _ I D E N T I T Y _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_POPULATION_IDENTITY_PRIVATE_H
#define LIBBOBOL_LOD_POPULATION_IDENTITY_PRIVATE_H

#include <cstddef>
#include <cstdint>
#include <vector>

/*
 * Exact identity of one retained occurrence representation.  The payload
 * addresses are valid only inside the allocation's already checked semantic
 * epoch; the remaining fields distinguish every representation choice made
 * for those immutable payloads.
 */
struct BObolLodPopulationMember {
    uintptr_t payload = 0;
    int cut = -1;
    int drawMode = 0;
    bool pointProxy = false;

    bool operator==(const BObolLodPopulationMember &other) const
    {
	return this->payload == other.payload && this->cut == other.cut &&
	    this->drawMode == other.drawMode &&
	    this->pointProxy == other.pointProxy;
    }
};

struct BObolLodPopulationIdentity {
    size_t fixedPresentationCost = 0;
    std::vector<BObolLodPopulationMember> members;

    bool operator==(const BObolLodPopulationIdentity &other) const
    {
	return this->fixedPresentationCost == other.fixedPresentationCost &&
	    this->members == other.members;
    }
};

#endif /* LIBBOBOL_LOD_POPULATION_IDENTITY_PRIVATE_H */
