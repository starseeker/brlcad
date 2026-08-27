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
#include <utility>
#include <vector>

class BObolLodService;
class BObolLodProjectedDemandCache;
class BObolViewController;
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
    /* Mirror SoCADAssembly's current aggregate-point threshold during source
     * admission.  An unselected structural leaf already below this projected
     * size is a complete pixel-appropriate presentation; opening/building its
     * PoP asset cannot improve the current frame. */
    void setPointProxyPixelThreshold(float pixels);
    /* During an initial/current-inventory coverage census, the standing leaf
     * proxy is already useful data.  Count it as coverage without opening a
     * PoP hierarchy; the following quality pass promotes view-significant
     * leaves under the aggregate scene budget. */
    void setStructuralCoverageOnly(SbBool coverageOnly);
    /* Permit visible, non-subpixel BoTs from a complete small-scene profile
     * to request their terminal mesh directly when its exact render cost fits
     * the aggregate scene allowance.  Service working-set, result, resident,
     * and frame governors remain authoritative. */
    void setAllowTerminalMeshAdmission(SbBool allow);
    SbBool getAllowTerminalMeshAdmission(void) const;
    /* A completed renderer frame may prove that one or more predicted
     * subpixel structural proxies did not actually collapse.  During the
     * resulting repair pass, require those otherwise eligible fallbacks to
     * enter the mesh provider path. */
    void setStructuralPresentationRepair(SbBool repair);
    /* Divide one exact-frame structural repair allowance across the complete
     * unresolved frontier.  The controller supplies a per-occurrence share;
     * workers still report the exact minimum PoP population and the completed
     * framebuffer remains the capacity authority. */
    void setStructuralCoverageCostReservation(size_t cost);
    /* Exactly one effective selected occurrence may bypass a subpixel proxy
     * so Coin geometry and edit manipulators can be prepared.  Bulk
     * selection only restyles valid point/box/mesh presentations and must not
     * create scene-wide geometry demand.  Callers may saturate counts above
     * one.  The scene owner must supply its global saturated count; an unset
     * count deliberately does not promote geometry because a source-local
     * population cannot prove that the scene has only one selection. */
    void setSelectedOccurrenceCount(size_t count);
    size_t getSelectedOccurrenceCount(void) const;
    void setGeneration(uint64_t generation);
    uint64_t getGeneration(void) const;
    void setRevisions(uint64_t viewRevision, uint64_t policyRevision);
    void setProvider(const char *providerId, const char *providerVersion);
    void setQualityTier(int qualityTier);
    void setRefreshMissing(SbBool refreshMissing);
    void setReset(SbBool newResetExisting);
    void setForcedCut(int cut);
    void clearForcedCut(void);
    SbBool hasForcedCut(void) const;
    int getForcedCut(void) const;
    void setRequireLodBacked(SbBool requireLodBacked);
    SbBool getRequireLodBacked(void) const;
    void setAllowCutDowngrade(SbBool allow);
    SbBool getAllowCutDowngrade(void) const;
    void setAllowRetainedRefinement(SbBool allow);
    SbBool getAllowRetainedRefinement(void) const;
    /* Cut hysteresis is a presentation stabilizer for an actively changing
     * view, not part of the physical pixel-error demand.  Quiet/final passes
     * must leave it disabled so requestedCut remains the first producer cut
     * which actually satisfies targetPixelError. */
    void setCutHysteresisEnabled(SbBool enabled);
    SbBool getCutHysteresisEnabled(void) const;
    /* Permit a pass to request a missing suffix up to its current allocated
     * presentation cut.  Worker working-set and resident-memory admission
     * remain authoritative. */
    void setAllowResidentPrefetch(SbBool allow);
    SbBool getAllowResidentPrefetch(void) const;
    /* Active scale interaction may prefetch one bounded transition beyond a
     * stale/conservative allocation so zoom does not magnify a fixed coarse
     * prefix.  Quiet views leave this disabled: unconstrained physical demand
     * remains quality debt, not permission to fill memory with cuts the scene
     * allocator cannot present. */
    void setAllowResidentPrefetchPastAllocation(SbBool allow);
    SbBool getAllowResidentPrefetchPastAllocation(void) const;
    /* Bound upward occurrence-cut admission to the renderer's current
     * presentation ceiling.  Existing richer cuts are retained and hidden by
     * the renderer; a negative value disables the bound. */
    void setRefinementCutCeiling(int cut);
    int getRefinementCutCeiling(void) const;
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
     * provider request for the next missing cut. */
    void setRefinementCostBudget(size_t additionalCost);
    size_t getRefinementCostBudget(void) const;
    size_t getRefinementCostBudgetUsed(void) const;
    /* Override the ordinary many-object first-prefix quantum for a provider
     * created by this action.  Controllers use this only when one uncovered
     * visible occurrence owns the complete scene allowance. */
    void setInitialProviderCostBudget(size_t cost);
    size_t getInitialProviderCostBudget(void) const;
    unsigned int getRefinementBudgetBlockedCount(void) const;
    /* Visible occurrences whose structural proxy could not be replaced by a
     * first mesh because the scene allowance was exhausted.  This excludes
     * richer-cut debt on occurrences which already own a mesh. */
    unsigned int getMissingMeshBudgetBlockedCount(void) const;
    /* Split retained minimax observations into a deliberately coarser
     * scene-quality ceiling and failure to reach the cut allocated by that
     * ceiling.  Explicit recovery terminates at the former but must retry the
     * latter. */
    unsigned int getRetainedQualityLimitedCount(void) const;
    unsigned int getRetainedAdmissionBlockedCount(void) const;
    /* Permit exactly one next-population transition whose marginal cost may
     * exceed the remaining ordinary refinement budget by at most this amount.
     * The complete marginal cost is still reported as budget used.  A
     * controller owns scene-global uniqueness across bounded actions. */
    void setOneOverBudgetRefinementLimit(size_t excessCost);
    size_t getOneOverBudgetRefinementLimit(void) const;
    SbBool getOneOverBudgetRefinementUsed(void) const;
    /* Limit one finite-budget visit to one populated PoP transition.  Stable
     * retained refinement normally jumps to the richest budget-fitting cut;
     * this stricter pacing is reserved for motion-time quality probes which
     * require a completed-frame deadline between transitions. */
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
    /* Large-scene recovery already guarantees one minimum drawable mesh per
     * retained occurrence.  This allowance is only for PoP detail above that
     * floor and is therefore safely carried across bounded planning windows
     * without charging the coverage population repeatedly. */
    void setRetainedSceneUpgradeCostBudget(size_t additionalCost);
    size_t getRetainedSceneUpgradeCostBudgetUsed(void) const;
    /* A retained importance pass uses one scene-wide normalized screen-error
     * ceiling.  Each occurrence selects the cheapest resident PoP cut which
     * meets that ceiling; the aggregate upgrade budget remains the hard
     * discrete-cost safety check.  Non-finite values disable the ceiling. */
    void setRetainedSceneMaximumNormalizedError(double error);
    double getRetainedSceneMaximumNormalizedError(void) const;
    size_t getCompactEntryNext(void) const;
    size_t getCompactEntryTotal(void) const;
    SbBool hasDeferredCompactEntries(void) const;

    unsigned int getVisitedMeshCount(void) const;
    unsigned int getSubmittedTaskCount(void) const;
    unsigned int getUpdatedCutCount(void) const;
    unsigned int getPendingRetainedRefinementCount(void) const;
    /* Retained allocation found a pixel-demanded cut whose immutable PoP
     * suffix is not resident yet.  Unlike a budget-limited richer-cut
     * observation, this is an actionable provider-work witness. */
    unsigned int getPendingResidentRefinementCount(void) const;
    unsigned int getSkippedMeshCount(void) const;
    /* Projected compact occurrences visited by this bounded window, and the
     * subset which already had a drawable payload when visited.  Controllers
     * use these counters to finish a coverage-only pass before spending the
     * scene budget on richer prefixes. */
    size_t getVisibleMeshCount(void) const;
    size_t getCoveredVisibleMeshCount(void) const;
    /* Compact-entry visibility observations made by this bounded action.
     * The entry index is source-stable and the boolean says whether that
     * occurrence was an on-screen mesh LoD target for this view.  A
     * controller uses these observations to update a completed visibility
     * census after an exact source delta without rescanning every unchanged
     * occurrence. */
    const std::vector<std::pair<size_t, SbBool>> &
	getCompactEntryVisibilityObservations(void) const;
    const SbString &getDiagnostics(void) const;

protected:
    virtual void beginTraversal(SoNode *node);

private:
    friend class BObolViewController;

    /* Controller-owned, presentation-thread-only cache.  It is deliberately
     * private API: projected evidence is an implementation detail of one view
     * and must not become source-global mutable state. */
    void setProjectedDemandCache(BObolLodProjectedDemandCache *cache);
    static void nodeAction(SoAction *action, SoNode *node);
    static void databaseSourceAction(SoAction *action, SoNode *node);
    static void meshShapeAction(SoAction *action, SoNode *node);
    void appendDiagnostic(const SbString &target, const char *message);
    SbBool reserveRefinementCost(
	const BObolLodProgressiveMeshPtr &progressiveMesh,
	const std::vector<uint32_t> &chunkIds, int activeCut, int nextCut,
	int drawMode, SbBool hasNormals);
    int reserveRefinementCut(
	const BObolLodProgressiveMeshPtr &progressiveMesh,
	const std::vector<uint32_t> &chunkIds, int activeCut, int preferredCut,
	int drawMode, SbBool hasNormals);
    int admitAllocatedRefinementCut(
	const BObolLodProgressiveMeshPtr &progressiveMesh,
	const std::vector<uint32_t> &chunkIds, int activeCut, int preferredCut,
	int drawMode, SbBool hasNormals, SbBool allocationCoversCut);
    /* Reserve a conservative first-cut population before its hierarchy has
     * been opened by a worker.  This closes the all-box zero-face blind spot:
     * thousands of independent warm-cache requests must not each interpret
     * an aggregate scene budget as their own private allowance. */
    SbBool reserveInitialCost(uint64_t sourceFaces, uint64_t sourcePoints,
	int drawMode, size_t &providerCostAllowance);
    SbBool reserveInitialCost(const BObolLodCounts &counts, int drawMode);

    BObolLodService *service;
    struct db_i *dbip;
    SbString databaseId;
    uint64_t databaseRevision;
    struct bv_view_info view;
    SbViewVolume viewVolume;
    SbBool useViewVolume;
    float targetPixelError;
    float pointProxyPixelThreshold;
    SbBool structuralCoverageOnly;
    SbBool allowTerminalMeshAdmission;
    SbBool structuralPresentationRepair;
    size_t structuralCoverageCostReservation;
    size_t selectedOccurrenceCount;
    uint64_t generation;
    uint64_t viewRevision;
    uint64_t policyRevision;
    SbString providerId;
    SbString providerVersion;
    int qualityTier;
    SbBool refreshMissing;
    SbBool resetExisting;
    SbBool useForcedCut;
    int forcedCut;
    SbBool requireLodBacked;
    SbBool allowCutDowngrade;
    SbBool allowRetainedRefinement;
    SbBool cutHysteresisEnabled;
    SbBool allowResidentPrefetch;
    SbBool allowResidentPrefetchPastAllocation;
    int refinementCutCeiling;
    SbBool allowRepresentationRefinement;
    SbBool preserveMeshCoverage;
    size_t refinementCostBudget;
    size_t refinementCostBudgetUsed;
    size_t initialProviderCostBudget;
    unsigned int refinementBudgetBlockedCount;
    unsigned int missingMeshBudgetBlockedCount;
    unsigned int retainedQualityLimitedCount;
    unsigned int retainedAdmissionBlockedCount;
    size_t oneOverBudgetRefinementLimit;
    SbBool oneOverBudgetRefinementUsed;
    SbBool transitionLimitedRefinement;
    BObolViewLodState *viewState;
    BObolLodProjectedDemandCache *projectedDemandCache;
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
    size_t retainedSceneUpgradeCostBudget;
    size_t retainedSceneUpgradeCostBudgetUsed;
    double retainedSceneMaximumNormalizedError;
    /* Occurrences charged by the scene-wide retained recovery pass.  The
     * bounded provider window must not charge them a second time. */
    std::unordered_set<std::string> retainedRecoveredOccurrences;
    SbBool deferredCompactEntries;
    unsigned int visitedMeshCount;
    unsigned int submittedTaskCount;
    unsigned int updatedCutCount;
    unsigned int pendingRetainedRefinementCount;
    unsigned int pendingResidentRefinementCount;
    unsigned int skippedMeshCount;
    size_t visibleMeshCount;
    size_t coveredVisibleMeshCount;
    std::vector<std::pair<size_t, SbBool>>
	compactEntryVisibilityObservations;
    unsigned int diagnosticCount;
    unsigned int suppressedDiagnosticCount;
    SbString diagnostics;
};

#endif /* BOBOL_BMESHLODSUBMITACTION_H */
