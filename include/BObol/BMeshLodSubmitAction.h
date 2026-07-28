/*        B M E S H L O D S U B M I T A C T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BMeshLodSubmitAction.h */

#ifndef BOBOL_BMESHLODSUBMITACTION_H
#define BOBOL_BMESHLODSUBMITACTION_H

#include "BObol/BDefines.h"
#include "BObol/BLodRealization.h"
#include "bv/view.h"

#include <Inventor/SbString.h>
#include <Inventor/SbViewVolume.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/actions/SoSubAction.h>

#include <stdint.h>
#include <stddef.h>
#include <vector>

class BObolLodService;
class BObolViewLodState;
struct db_i;

class BOBOL_EXPORT SoBRLMeshLodSubmitAction : public SoAction {
    typedef SoAction inherited;

    SO_ACTION_HEADER(SoBRLMeshLodSubmitAction);

public:
    SoBRLMeshLodSubmitAction(void);
    virtual ~SoBRLMeshLodSubmitAction(void);
    static void initClass(void);

    void setService(BObolLodService *service);
    BObolLodService *getService(void) const;
    void setDatabase(struct db_i *dbip, const char *databaseId = NULL,
	    uint64_t databaseRevision = 0);
    void setViewInfo(const struct bv_view_info *info);
    const struct bv_view_info &getViewInfo(void) const;
    void setViewVolume(const SbViewVolume *volume,
	float targetPixelError = 1.0f);
    void setGeneration(uint64_t generation);
    uint64_t getGeneration(void) const;
    void setRevisions(uint64_t viewRevision, uint64_t policyRevision);
    void setProvider(const char *providerId, const char *providerVersion);
    void setQualityTier(int qualityTier);
    void setRefreshMissing(SbBool refreshMissing);
    void setReset(int reset);
    void setForcedLevel(int level);
    void clearForcedLevel(void);
    SbBool hasForcedLevel(void) const;
    int getForcedLevel(void) const;
    void setRequireLodBacked(SbBool requireLodBacked);
    SbBool getRequireLodBacked(void) const;
    void setAllowLevelDowngrade(SbBool allow);
    SbBool getAllowLevelDowngrade(void) const;
    void setAllowRetainedRefinement(SbBool allow);
    SbBool getAllowRetainedRefinement(void) const;
    /* Bound aggregate upward PoP growth selected by this traversal.  The
     * budget applies equally to an already-resident cut change and to a
     * provider request for the next missing level. */
    void setRefinementFaceBudget(size_t additionalFaces);
    size_t getRefinementFaceBudget(void) const;
    size_t getRefinementFaceBudgetUsed(void) const;
    unsigned int getRefinementBudgetBlockedCount(void) const;
    void setViewLodState(BObolViewLodState *viewState);
    const BObolViewLodState *getViewLodState(void) const;
    void setCompactEntryRange(size_t first, size_t count);
    /* Pin one deterministic compact-entry order across bounded service
     * windows.  Without this plan a many-leaf source is projected and sorted
     * again for every result-queue window, and a camera change can reorder a
     * partially consumed cursor so entries are skipped. */
    void setCompactEntryPlan(const std::vector<size_t> &entryIndices);
    void getCompactEntryPlan(std::vector<size_t> &entryIndices) const;
    /* Provider queue capacity limits tasks, not cheap scene decisions.
     * Retargeting resident cuts and binding shared assets may scan the whole
     * pinned plan without consuming a result slot. */
    void setSubmissionTaskLimit(size_t taskCount);
    /* When an already-resident scene exceeds its calibrated total budget,
     * re-admit retained occurrences in pinned priority order.  Occurrences
     * whose minimum cut cannot fit fall back to their structural proxy. */
    void setRetainedSceneFaceBudget(size_t totalFaces);
    size_t getRetainedSceneFaceBudgetUsed(void) const;
    size_t getCompactEntryNext(void) const;
    size_t getCompactEntryTotal(void) const;
    SbBool hasDeferredCompactEntries(void) const;

    unsigned int getVisitedMeshCount(void) const;
    unsigned int getSubmittedTaskCount(void) const;
    unsigned int getUpdatedCutCount(void) const;
    unsigned int getPendingRetainedRefinementCount(void) const;
    unsigned int getSkippedMeshCount(void) const;
    const SbString &getDiagnostics(void) const;

protected:
    virtual void beginTraversal(SoNode *node);

private:
    static void nodeAction(SoAction *action, SoNode *node);
    static void databaseSourceAction(SoAction *action, SoNode *node);
    static void meshShapeAction(SoAction *action, SoNode *node);
    void appendDiagnostic(const SbString &target, const char *message);
    SbBool reserveRefinementFaces(
	const BObolLodProgressiveMeshPtr &progressiveMesh,
	int activeLevel, int nextLevel);
    /* Reserve a conservative first-cut population before its hierarchy has
     * been opened by a worker.  This closes the all-box zero-face blind spot:
     * thousands of independent warm-cache requests must not each interpret
     * an aggregate scene budget as their own private allowance. */
    SbBool reserveInitialFaces(uint64_t sourceFaces,
	size_t &providerFaceAllowance);
    SbBool reserveInitialFaceCost(size_t faceCount);

    BObolLodService *service;
    struct db_i *dbip;
    SbString databaseId;
    uint64_t databaseRevision;
    struct bv_view_info view;
    SbViewVolume viewVolume;
    SbBool useViewVolume;
    float targetPixelError;
    uint64_t generation;
    uint64_t viewRevision;
    uint64_t policyRevision;
    SbString providerId;
    SbString providerVersion;
    int qualityTier;
    SbBool refreshMissing;
    int reset;
    SbBool useForcedLevel;
    int forcedLevel;
    SbBool requireLodBacked;
    SbBool allowLevelDowngrade;
    SbBool allowRetainedRefinement;
    size_t refinementFaceBudget;
    size_t refinementFaceBudgetUsed;
    unsigned int refinementBudgetBlockedCount;
    BObolViewLodState *viewState;
    size_t compactEntryFirst;
    size_t compactEntryLimit;
    size_t compactEntryNext;
    size_t compactEntryTotal;
    std::vector<size_t> compactEntryPlan;
    SbBool compactEntryPlanSupplied;
    size_t submissionTaskLimit;
    size_t retainedSceneFaceBudget;
    size_t retainedSceneFaceBudgetUsed;
    SbBool deferredCompactEntries;
    unsigned int visitedMeshCount;
    unsigned int submittedTaskCount;
    unsigned int updatedCutCount;
    unsigned int pendingRetainedRefinementCount;
    unsigned int skippedMeshCount;
    unsigned int diagnosticCount;
    unsigned int suppressedDiagnosticCount;
    SbString diagnostics;
};

#endif /* BOBOL_BMESHLODSUBMITACTION_H */
