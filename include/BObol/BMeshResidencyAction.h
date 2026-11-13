/*        B M E S H R E S I D E N C Y A C T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BMeshResidencyAction.h */

#ifndef BOBOL_BMESHRESIDENCYACTION_H
#define BOBOL_BMESHRESIDENCYACTION_H

#include "BObol/BDefines.h"

#include <Inventor/SbBasic.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/actions/SoSubAction.h>

#include <stddef.h>
#include <vector>

class SoBRLMeshShape;

class BOBOL_EXPORT SoBRLMeshResidencyAction : public SoAction {
    typedef SoAction inherited;

    SO_ACTION_HEADER(SoBRLMeshResidencyAction);

public:
    SoBRLMeshResidencyAction(void);
    virtual ~SoBRLMeshResidencyAction(void);
    static void initClass(void);

    void setMaxResidentMeshBytes(size_t maxBytes);
    size_t getMaxResidentMeshBytes(void) const;
    void setEvictDisplayPayloads(SbBool enabled);
    SbBool isEvictDisplayPayloadsEnabled(void) const;

    unsigned int getVisitedMeshCount(void) const;
    unsigned int getEvictedFullDetailMeshCount(void) const;
    unsigned int getEvictedDisplayMeshCount(void) const;
    unsigned int getSkippedDisplayMeshCount(void) const;
    size_t getInitialResidentMeshBytes(void) const;
    size_t getFinalResidentMeshBytes(void) const;
    size_t getFreedFullDetailBytes(void) const;
    size_t getFreedDisplayBytes(void) const;
    size_t getFreedResidentMeshBytes(void) const;

protected:
    virtual void beginTraversal(SoNode *node);

private:
    struct Entry {
	SoBRLMeshShape *shape;
	size_t residentBytes;
	size_t fullDetailBytes;
	size_t displayBytes;
    };

    static void nodeAction(SoAction *action, SoNode *node);
    static void meshShapeAction(SoAction *action, SoNode *node);
    void resetResults(void);
    void evictToBudget(void);
    void recomputeFinalResidentBytes(void);

    size_t maxResidentMeshBytes;
    SbBool evictDisplayPayloads;
    std::vector<Entry> entries;
    unsigned int visitedMeshCount;
    unsigned int evictedFullDetailMeshCount;
    unsigned int evictedDisplayMeshCount;
    unsigned int skippedDisplayMeshCount;
    size_t initialResidentMeshBytes;
    size_t finalResidentMeshBytes;
    size_t freedFullDetailBytes;
    size_t freedDisplayBytes;
};

#endif /* BOBOL_BMESHRESIDENCYACTION_H */
