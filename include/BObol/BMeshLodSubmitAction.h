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
#include <string>
#include <unordered_set>
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
    /* Permit a zoom-in pass to request a missing pixel-demanded cache suffix
     * even when the aggregate render-cost budget cannot expose it yet.
     * Worker working-set and resident-memory admission remain authoritative;
     * this only separates residency from presentation cost. */
    void setAllowResidentPrefetch(SbBool allow);
    SbBool getAllowResidentPrefetch(void) const;
    /* Bound upward occurrence-cut admission to the renderer's current
     * presentation ceiling.  Existing richer cuts are retained and hidden by
     * the renderer; a negative value disables the bound. */
    void setRefinementLevelCeiling(int level);
    int getRefinementLevelCeiling(void) const;
    /* Source representation generation (currently adaptive BREP
     * tessellation) is quiet-view work.  It is independent of inexpensive
     * resident PoP cut changes, which may continue during zoom interaction. */
    void setAllowRepresentationRefinement(SbBool allow);
    SbBool getAllowRepresentationRefinement(void) const;
    /* A coverage pass may exceed the learned triangle budget by the minimum
     * drawable prefix of each visible mesh.  This is a presentation floor,
     * not refinement: renderer-side point aggregation is responsible for
     * keeping that floor interactive without reverting available meshes to
     * boxes. */
    void setPreserveMeshCoverage(SbBool preserve);
    SbBool getPreserveMeshCoverage(void) const;
    /* Bound aggregate upward PoP growth selected by this traversal.  The
     * budget applies equally to an already-resident cut change and to a
     * provider request for the next missing level. */
    void setRefinementCostBudget(size_t additionalCost);
    size_t getRefinementCostBudget(void) const;
    size_t getRefinementCostBudgetUsed(void) const;
    unsigned int getRefinementBudgetBlockedCount(void) const;
    /* Permit exactly one next-population transition whose marginal cost may
     * exceed the remaining ordinary refinement budget by at most this amount.
     * The complete marginal cost is still reported as budget used.  A
     * controller owns scene-global uniqueness across bounded actions. */
    void setOneOverBudgetRefinementLimit(size_t excessCost);
    size_t getOneOverBudgetRefinementLimit(void) const;
    SbBool getOneOverBudgetRefinementUsed(void) const;
    /* Limit one finite-budget visit to one populated PoP transition.  This is
     * used by perceptually ordered refinement frontiers so the score for one
     * marginal step cannot accidentally grant the same occurrence every
     * remaining level and starve other visible features. */
    void setTransitionLimitedRefinement(SbBool limited);
    SbBool getTransitionLimitedRefinement(void) const;
    void setViewLodState(BObolViewLodState *viewState);
    const BObolViewLodState *getViewLodState(void) const;
    void setCompactEntryRange(size_t first, size_t count);
    /* Pin one deterministic compact-entry order across bounded service
     * windows.  Without this plan a many-leaf source is projected and sorted
     * again for every result-queue window, and a camera change can reorder a
     * partially consumed cursor so entries are skipped. */
    void setCompactEntryPlan(const std::vector<size_t> &entryIndices);
    /* Borrow a controller-owned plan for this synchronous action traversal.
     * The caller must retain the vector through apply().  This avoids
     * allocating and copying tens of thousands of indices for every bounded
     * UI-thread window. */
    void setCompactEntryPlanView(
	const std::vector<size_t> *entryIndices);
    void getCompactEntryPlan(std::vector<size_t> &entryIndices) const;
    /* Provider queue capacity limits tasks, not cheap scene decisions.
     * Retargeting resident cuts and binding shared assets may scan the whole
     * pinned plan without consuming a result slot. */
    void setSubmissionTaskLimit(size_t taskCount);
    /* Bound compact occurrence planning on an interactive caller thread.
     * The current entry is atomic; zero disables the deadline. */
    void setSubmissionTimeLimit(uint64_t microseconds);
    /* When an already-resident scene exceeds its calibrated total budget,
     * retarget retained occurrences to minimum prefixes in pinned priority
     * order.  The budget limits refinement above the mesh-coverage floor. */
    void setRetainedSceneCostBudget(size_t totalCost);
    size_t getRetainedSceneCostBudgetUsed(void) const;
    size_t getCompactEntryNext(void) const;
    size_t getCompactEntryTotal(void) const;
    SbBool hasDeferredCompactEntries(void) const;

    unsigned int getVisitedMeshCount(void) const;
    unsigned int getSubmittedTaskCount(void) const;
    unsigned int getUpdatedCutCount(void) const;
    unsigned int getPendingRetainedRefinementCount(void) const;
    unsigned int getSkippedMeshCount(void) const;
    /* Projected compact occurrences visited by this bounded window, and the
     * subset which already had a drawable payload when visited.  Controllers
     * use these counters to finish a coverage-only pass before spending the
     * scene budget on richer prefixes. */
    size_t getVisibleMeshCount(void) const;
    size_t getCoveredVisibleMeshCount(void) const;
    const SbString &getDiagnostics(void) const;

protected:
    virtual void beginTraversal(SoNode *node);

private:
    static void nodeAction(SoAction *action, SoNode *node);
    static void databaseSourceAction(SoAction *action, SoNode *node);
    static void meshShapeAction(SoAction *action, SoNode *node);
    void appendDiagnostic(const SbString &target, const char *message);
    SbBool reserveRefinementCost(
	const BObolLodProgressiveMeshPtr &progressiveMesh,
	int activeLevel, int nextLevel, int drawMode);
    int reserveRefinementLevel(
	const BObolLodProgressiveMeshPtr &progressiveMesh,
	int activeLevel, int preferredLevel, int drawMode);
    /* Reserve a conservative first-cut population before its hierarchy has
     * been opened by a worker.  This closes the all-box zero-face blind spot:
     * thousands of independent warm-cache requests must not each interpret
     * an aggregate scene budget as their own private allowance. */
    SbBool reserveInitialCost(uint64_t sourceFaces, uint64_t sourcePoints,
	int drawMode, size_t &providerFaceAllowance);
    SbBool reserveInitialCost(const BObolLodCounts &counts, int drawMode);

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
    SbBool allowResidentPrefetch;
    int refinementLevelCeiling;
    SbBool allowRepresentationRefinement;
    SbBool preserveMeshCoverage;
    size_t refinementCostBudget;
    size_t refinementCostBudgetUsed;
    unsigned int refinementBudgetBlockedCount;
    size_t oneOverBudgetRefinementLimit;
    SbBool oneOverBudgetRefinementUsed;
    SbBool transitionLimitedRefinement;
    BObolViewLodState *viewState;
    size_t compactEntryFirst;
    size_t compactEntryLimit;
    size_t compactEntryNext;
    size_t compactEntryTotal;
    std::vector<size_t> compactEntryPlan;
    const std::vector<size_t> *compactEntryPlanView;
    SbBool compactEntryPlanSupplied;
    size_t submissionTaskLimit;
    uint64_t submissionTimeLimitMicroseconds;
    size_t retainedSceneCostBudget;
    size_t retainedSceneCostBudgetUsed;
    /* Occurrences charged by the scene-wide retained recovery pass.  The
     * bounded provider window must not charge them a second time. */
    std::unordered_set<std::string> retainedRecoveredOccurrences;
    SbBool deferredCompactEntries;
    unsigned int visitedMeshCount;
    unsigned int submittedTaskCount;
    unsigned int updatedCutCount;
    unsigned int pendingRetainedRefinementCount;
    unsigned int skippedMeshCount;
    size_t visibleMeshCount;
    size_t coveredVisibleMeshCount;
    unsigned int diagnosticCount;
    unsigned int suppressedDiagnosticCount;
    SbString diagnostics;
};

#endif /* BOBOL_BMESHLODSUBMITACTION_H */
