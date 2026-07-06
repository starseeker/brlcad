/*              L O D _ U P D A T E _ A C T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/lod_update_action.h */

#ifndef BRLOBOL_LOD_UPDATE_ACTION_H
#define BRLOBOL_LOD_UPDATE_ACTION_H

#include "brlobol/lod_realization.h"

#include <Inventor/SbString.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/actions/SoSubAction.h>

#include <stddef.h>
#include <vector>

class BRLObolLodService;
class BRLObolViewLodState;

class BRLOBOL_EXPORT SoBRLLodUpdateAction : public SoAction {
    typedef SoAction inherited;

    SO_ACTION_HEADER(SoBRLLodUpdateAction);

public:
    SoBRLLodUpdateAction(void);
    virtual ~SoBRLLodUpdateAction(void);
    static void initClass(void);

    void clearResults(void);
    void addResult(const BRLObolLodResult &result);
    void setResults(const std::vector<BRLObolLodResult> &results);
    size_t drainService(BRLObolLodService &service, size_t maxResults = 0);
    void setViewLodState(BRLObolViewLodState *viewState);
    BRLObolViewLodState *getViewLodState(void) const;

    size_t getResultCount(void) const;
    unsigned int getMatchedResultCount(void) const;
    unsigned int getAppliedResultCount(void) const;
    unsigned int getRejectedResultCount(void) const;
    unsigned int getUnmatchedResultCount(void) const;
    const SbString &getDiagnostics(void) const;

protected:
    virtual void beginTraversal(SoNode *node);

private:
    static void nodeAction(SoAction *action, SoNode *node);
    static void databaseSourceAction(SoAction *action, SoNode *node);
    static void meshShapeAction(SoAction *action, SoNode *node);
    void appendDiagnostic(const BRLObolLodResult &result,
	const char *message);
    void finalizeUnmatchedDiagnostics(void);

    std::vector<BRLObolLodResult> results;
    std::vector<SbBool> matched;
    BRLObolViewLodState *viewState;
    unsigned int matchedResultCount;
    unsigned int appliedResultCount;
    unsigned int rejectedResultCount;
    unsigned int unmatchedResultCount;
    SbString diagnostics;
};

#endif /* BRLOBOL_LOD_UPDATE_ACTION_H */
