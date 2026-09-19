/*            B L O D U P D A T E A C T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BLodUpdateAction.h */

#ifndef BOBOL_BLODUPDATEACTION_H
#define BOBOL_BLODUPDATEACTION_H

#include "BObol/BLodRealization.h"

#include <Inventor/SbString.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/actions/SoSubAction.h>

#include <stddef.h>
#include <vector>

class BObolLodService;
class BObolViewLodState;

class BOBOL_EXPORT SoBRLLodUpdateAction : public SoAction {
    typedef SoAction inherited;

    SO_ACTION_HEADER(SoBRLLodUpdateAction);

public:
    SoBRLLodUpdateAction(void);
    virtual ~SoBRLLodUpdateAction(void);
    static void initClass(void);

    void clearResults(void);
    void addResult(const BObolLodResult &result);
    void addResult(BObolLodResult &&result);
    void setResults(const std::vector<BObolLodResult> &results);
    size_t drainService(BObolLodService &service, size_t maxResults = 0);
    void setViewLodState(BObolViewLodState *viewState);
    BObolViewLodState *getViewLodState(void) const;

    size_t getResultCount(void) const;
    unsigned int getMatchedResultCount(void) const;
    unsigned int getAppliedResultCount(void) const;
    unsigned int getRejectedResultCount(void) const;
    unsigned int getCurrentDemandRetryResultCount(void) const;
    unsigned int getUnmatchedResultCount(void) const;
    const SbString &getDiagnostics(void) const;

protected:
    virtual void beginTraversal(SoNode *node);

private:
    static void nodeAction(SoAction *action, SoNode *node);
    static void databaseSourceAction(SoAction *action, SoNode *node);
    static void meshShapeAction(SoAction *action, SoNode *node);
    void appendDiagnostic(const BObolLodResult &result,
	const char *message);
    void finalizeUnmatchedDiagnostics(void);

    std::vector<BObolLodResult> results;
    std::vector<SbBool> matched;
    BObolViewLodState *viewState;
    unsigned int matchedResultCount;
    unsigned int appliedResultCount;
    unsigned int rejectedResultCount;
    unsigned int currentDemandRetryResultCount;
    unsigned int unmatchedResultCount;
    SbString diagnostics;
};

#endif /* BOBOL_BLODUPDATEACTION_H */
