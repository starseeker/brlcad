/*                    B V I E W L O D . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BViewLod.h */

#ifndef BOBOL_BVIEWLOD_H
#define BOBOL_BVIEWLOD_H

#include "BObol/BDefines.h"
#include "BObol/BLodRealization.h"

#include <Inventor/SbBasic.h>
#include <Inventor/SbBox.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/elements/SoElement.h>
#include <Inventor/elements/SoSubElement.h>
#include <Inventor/nodes/SoGroup.h>

#include <array>
#include <memory>
#include <stddef.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class SoAction;
class SoCallbackAction;
class SoGetBoundingBoxAction;
class SoGLRenderAction;
class SoNode;
class SoPickAction;
class SoState;
class SoBRLMeshShape;
class SoBRLDatabaseSource;
class SoCADAssembly;

/**
 * View-local active LoD bindings for shared Obol geometry.
 *
 * Shared database objects keep their full-detail/source identity in the scene
 * graph.  Each view controller owns one BObolViewLodState that records the
 * active display payloads selected for that view.  Coin actions reach the
 * state through SoBRLViewLodElement rather than by mutating shared mesh nodes.
 */
class BOBOL_EXPORT BObolViewLodState
{
public:
    /** Renderer-neutral aggregate of the most recent complete-frame CAD GPU
     * resource snapshots.  Assembly identity is deduplicated before summing. */
    struct BOBOL_EXPORT CadGpuResourceStatus {
        size_t trackedBufferBytes = 0;
        size_t ordinaryPartBufferBytes = 0;
        size_t progressiveCutBufferBytes = 0;
        size_t progressiveActiveCutBytes = 0;
        size_t batchBufferBytes = 0;
        size_t triangleAtlasAllocatedBytes = 0;
        size_t triangleAtlasLiveBytes = 0;
        size_t triangleAtlasConfiguredCapacityBytes = 0;
        size_t triangleAtlasPartCount = 0;
        size_t triangleAtlasPageCount = 0;
        uint64_t ordinaryPartFullUploadBytes = 0;
        uint64_t ordinaryPartSuffixUploadBytes = 0;
        uint64_t ordinaryPartGpuCopyBytes = 0;
        uint64_t ordinaryPartLineageReuseCount = 0;
        uint64_t ordinaryPartLineageReplacementCount = 0;
        uint64_t triangleAtlasFullUploadBytes = 0;
        uint64_t triangleAtlasSuffixUploadBytes = 0;
        uint64_t triangleAtlasLineageReuseCount = 0;
        size_t pressureProxyCount = 0;
        uint64_t progressiveEvictionCount = 0;
        uint64_t triangleAtlasReclamationCount = 0;
        uint64_t sampleSerial = 0;
        SbBool memoryPressure = FALSE;
    };

    enum NormalStyle {
	NORMAL_AUTHORED = 0,
	NORMAL_FLAT = 1,
	NORMAL_SMOOTH = 2
    };

    struct BOBOL_EXPORT MeshPayload {
	BObolLodMeshPayload mesh;
	BObolLodProgressiveMeshPtr progressiveMesh;
	SbString sourcePath;
	SbString sourceName;
	SbString sourceIdentity;
	SbString cacheIdentity;
	SbString cacheKey;
	uint64_t sourceContentHash;
	int resultKind;
	int qualityTier;
	int providerStatus;
	int activeCut;
	int residentCut;
	int requestedCut;
	std::vector<uint32_t> requiredChunks;
	uint64_t residentAdmissionRevision;
	uint64_t viewRevision;
	uint64_t policyRevision;
	BObolLodCounts counts;
	SbBox3f bounds;
	SbBool hasSnappedPoints;
	SbBool hasNormals;
	SbBool shadedCullBackfaces;
	SbBool memoryLimited;
	SbString diagnostic;

	MeshPayload(void);
	SbBool isValid(void) const;
	size_t estimateBytes(void) const;
	int getTriangleCount(void) const;
	SbBool getTriangleVertexIndices(int triangleIndex,
					int &indexA,
					int &indexB,
					int &indexC) const;
	SbBool getTriangle(int triangleIndex,
			   SbVec3f &a,
			   SbVec3f &b,
			   SbVec3f &c) const;
    };
    typedef std::shared_ptr<MeshPayload> MeshPayloadPtr;

    struct BOBOL_EXPORT ProxyPayload {
	BObolLodProxy proxy;
	SbString sourcePath;
	SbString sourceName;
	SbString sourceIdentity;
	SbString cacheIdentity;
	SbString cacheKey;
	int resultKind;
	int qualityTier;
	int providerStatus;
	uint64_t viewRevision;
	uint64_t policyRevision;
	BObolLodCounts counts;
	SbBox3f bounds;
	SbString diagnostic;

	ProxyPayload(void);
	SbBool isValid(void) const;
	size_t estimateBytes(void) const;
    };
    typedef std::shared_ptr<ProxyPayload> ProxyPayloadPtr;

    struct BOBOL_EXPORT CadPayload {
	typedef std::array<BObolLodCounts,
	    BOBOL_MESH_LOD_CUT_COUNT_MAX> ProjectedCutCounts;

	BObolLodMeshPayload mesh;
	BObolLodProgressiveMeshPtr progressiveMesh;
	std::shared_ptr<const Obol::PartGeometry> preparedCadGeometry;
	uint64_t preparedCadGeometryRevision;
	BObolLodProxy proxy;
	SbString sourcePath;
	SbString sourceName;
	SbString sourceIdentity;
	SbString sourceInstanceKey;
	SbString sourceBindingKey;
	/* Validated source-local compact index, or UINT32_MAX for a noncompact
	 * binding.  This avoids repeated string hashing in 50k/150k retained
	 * allocation passes without making an index the semantic identity. */
	uint32_t sourceEntryIndex;
	SbString cacheIdentity;
	SbString cacheKey;
	uint64_t sourceContentHash;
	uint64_t databaseRevision;
	uint64_t sourceRevision;
	int resultKind;
	int qualityTier;
	int providerStatus;
	int drawMode;
	int activeCut;
	int residentCut;
	int requestedCut;
	std::vector<uint32_t> requiredChunks;
	/* Exact scene-budget allocation for this occurrence.  requestedCut is
	 * unconstrained view demand; allocatedCut is the richest cut this
	 * view/policy/mode may present inside the measured aggregate allowance.
	 * Keeping the result on the authoritative occurrence makes bounded GUI
	 * windows and asynchronous resident growth apply one immutable allocation
	 * rather than reconstructing it from a scalar ceiling. */
	int allocatedCut;
	int allocationDrawMode;
	uint64_t allocationViewRevision;
	uint64_t allocationPolicyRevision;
	/* Nonzero values belong to an owner-state allocation transaction.  The
	 * fields above are visible to consumers only after this serial becomes
	 * the state's active serial; serial zero is an immediate single-payload
	 * update used by small/offline paths. */
	uint64_t allocationPlanSerial;
	float projectedPixelDiameter;
	float projectedPixelArea;
	float projectedPixelPerimeter;
	SbBool projectedBoundsContained;
	float targetPixelError;
	/* A bounded demand window already classifies the spatial ranges visible
	 * for this occurrence.  Retain that immutable per-cut population only
	 * when the occurrence straddles the frustum; the quiet scene allocator
	 * can then consume the exact census without traversing the mesh hierarchy
	 * again on the owner thread.  Fully contained and unclustered occurrences
	 * deliberately keep this null and use ordinary whole-prefix counts. */
	std::shared_ptr<const ProjectedCutCounts> projectedCutCounts;
	uint64_t projectedCutCountsViewRevision;
	uint64_t projectedCutCountsPolicyRevision;
	uint64_t projectedCutCountsMeshRevision;
	uint64_t residentAdmissionRevision;
	uint64_t viewRevision;
	uint64_t policyRevision;
	uint8_t visualEmphasis;
	BObolLodCounts counts;
	SbBox3f bounds;
	SbBool hasSnappedPoints;
	SbBool hasNormals;
	SbBool shadedCullBackfaces;
	SbBool memoryLimited;
	SbString diagnostic;
	/* Non-owning provenance for O(1) validation of metadata-only hot-path
	 * updates.  A CadPayload is never published independently of its view
	 * state; raw payload handles have the same lifetime contract. */
	const BObolViewLodState *ownerState;

	CadPayload(void);
	SbBool isValid(void) const;
	size_t estimateBytes(void) const;
    };
    typedef std::shared_ptr<CadPayload> CadPayloadPtr;

    BObolViewLodState(void);
    ~BObolViewLodState(void);

    void clear(void);
    SbBool applyMeshResult(const SoBRLMeshShape *shape,
			   const BObolLodResult &result);
    SbBool applyProxyResult(const SoBRLMeshShape *shape,
			    const BObolLodResult &result);
    SbBool applyDisplayResult(const SoBRLMeshShape *shape,
			      const BObolLodResult &result);
    SbBool applySourceResult(const SoBRLDatabaseSource *source,
			     const BObolLodResult &result);
    SbBool consumeDisplayResult(const SoBRLMeshShape *shape,
	BObolLodResult &result);
    SbBool consumeSourceResult(const SoBRLDatabaseSource *source,
	BObolLodResult &result);
    const MeshPayload *findMesh(const SoBRLMeshShape *shape) const;
    const MeshPayload *findMeshForResult(const BObolLodResult &result) const;
    const ProxyPayload *findProxy(const SoBRLMeshShape *shape) const;
    const ProxyPayload *findProxyForResult(const BObolLodResult &result) const;
    const CadPayload *findCad(const SoBRLDatabaseSource *source) const;
    void findCadPayloads(const SoBRLDatabaseSource *source,
	std::vector<const CadPayload *> &payloads) const;
    /* Hot-path variant for scene planning and telemetry.  Compact
     * occurrence identity, rather than container iteration order, provides
     * determinism to those consumers, so avoid sorting thousands of payload
     * strings on every camera epoch. */
    void findCadPayloadsUnordered(const SoBRLDatabaseSource *source,
	std::vector<const CadPayload *> &payloads) const;
    const CadPayload *findCadForOccurrence(
	const SoBRLDatabaseSource *source,
	const SbString &occurrenceKey) const;
    /* Return a retained progressive asset representative for direct
     * occurrence binding.  Asset residency outlives an individual display
     * occurrence, so off-frustum admission does not force a cache reload. */
    const CadPayload *findCadForAsset(
	const SoBRLDatabaseSource *source,
	const SbString &assetPath) const;
    const CadPayload *findCadForResult(const BObolLodResult &result) const;
    const CadPayload *findCadForResult(const SoBRLDatabaseSource *source,
	const BObolLodResult &result) const;
    /* A terminal provider failure is coverage by the source occurrence's
     * structural fallback, not an invitation to resubmit the identical work
     * every quiet-frame pump.  Failures are demand scoped: a new camera,
     * policy, source revision, or requested cut may try again. */
    SbBool hasCadOccurrenceTerminalFailure(
	const SoBRLDatabaseSource *source,
	const BObolLodRequest &request) const;
    /** Number of occurrence-scoped provider failures retained for the
     * current demand epochs.  This is diagnostic/error state, not a mesh
     * payload and not permission to call an ordinary structural box
     * converged. */
    size_t cadOccurrenceTerminalFailureCount(void) const;
    size_t cadOccurrenceTerminalFailureCountForSource(
	const SoBRLDatabaseSource *source) const;
    /* Change only a view-local PoP cut.  No source/cache work is performed;
     * the call succeeds only when the retained progressive asset already
     * contains the requested prefix. */
    SbBool retargetMeshPayload(const MeshPayload *payload, int activeCut,
	int requestedCut, uint64_t viewRevision, uint64_t policyRevision);
    /** Retarget one retained CAD occurrence to an already resident cut and
     * record the exact view demand which selected it.  Keeping the projection
     * with the payload lets later bounded recovery retain perceptual ordering
     * without reprojecting or sorting an entire large scene on the owner
     * thread. */
    SbBool retargetCadPayload(const CadPayload *payload, int activeCut,
	const BObolLodRequest &demand);
    /* Record a metadata-only aggregate render-budget decision.  This neither
     * changes the active presentation nor journals an occurrence mutation;
     * later bounded windows consume it when its epochs and draw mode match. */
    SbBool setCadAllocatedCut(const CadPayload *payload, int allocatedCut,
	uint64_t viewRevision, uint64_t policyRevision, int drawMode);
    /** Reserve an unpublished allocation serial for a resumable plan. */
    uint64_t beginCadAllocationPlan(void);
    /** Stage one allocation inside an unpublished resumable plan.  Staged
     * metadata is ignored by presentation consumers until the matching plan
     * is committed. */
    SbBool stageCadAllocatedCut(const CadPayload *payload, int allocatedCut,
	uint64_t viewRevision, uint64_t policyRevision, int drawMode,
	uint64_t planSerial);
    /** Atomically expose all payload fields staged with @p planSerial. */
    SbBool commitCadAllocationPlan(uint64_t planSerial,
	uint64_t cadRevision, uint64_t residentDemandRevision);
    /** Serial against which consumers validate staged allocation metadata. */
    uint64_t activeCadAllocationPlan(void) const;
    /* Remove one view-local display binding while retaining its shared asset.
     * The source occurrence's structural fallback becomes visible again.
     * Used by scene-budget/frustum admission when an insignificant occurrence
     * should cost zero triangles rather than its minimum populated PoP cut. */
    SbBool removeCadPayload(const CadPayload *payload);
    SbBool removeMeshPayload(const MeshPayload *payload);
    size_t bindingCount(void) const;
    size_t payloadCount(void) const;
    size_t meshPayloadCount(void) const;
    size_t proxyPayloadCount(int proxyKind = BOBOL_LOD_PROXY_NONE) const;
    size_t cadPayloadCount(void) const;
    /** Return the number of retained occurrence bindings owned by one source
     * without copying or walking the source's payload table. */
    size_t cadPayloadCountForSource(
	const SoBRLDatabaseSource *source) const;
    size_t cadMeshPayloadCount(void) const;
    size_t cadMeshPayloadCountForSource(
	const SoBRLDatabaseSource *source) const;
    /** Copy the exact occurrence identities whose retained mesh cut has not
     * yet reached its current view target.  Maintained incrementally so a
     * stable refinement pass can visit only its work frontier rather than
     * rediscovering it by projecting every source leaf. */
    void unsatisfiedCadOccurrenceKeys(
        const SoBRLDatabaseSource *source,
        uint64_t residentAdmissionRevision,
        std::vector<SbString> &occurrenceKeys) const;
    /** Copy only resident-memory denials made actionable by a newer service
     * admission epoch.  This is the capacity-recovery frontier: ordinary
     * scene-budget quality debt is deliberately excluded so one reclaimed
     * asset cannot restart a complete large-scene demand pass. */
    void retriableMemoryLimitedCadOccurrenceKeys(
        const SoBRLDatabaseSource *source,
        uint64_t residentAdmissionRevision,
        std::vector<SbString> &occurrenceKeys) const;
    size_t cadProxyPayloadCount(int proxyKind = BOBOL_LOD_PROXY_NONE) const;
    /** Count current display payloads and those which have reached their
     * recorded view target.  This is telemetry only and does not mutate or
     * sort the hot occurrence maps. */
    void convergencePayloadCounts(size_t &active,
	size_t &satisfied, size_t &memoryLimited) const;
    /* O(1) for the compact CAD path.  Legacy shape payloads are included for
     * completeness and are expected to remain a small population. */
    size_t memoryLimitedPayloadCount(void) const;
    /* Estimated triangles submitted by the active view-local cuts.  Shared
     * arrays are counted once per displayed occurrence because render cost
     * follows instances, not storage aliases. */
    size_t activeFaceCount(void) const;
    /* Multi-cost scheduler population in shaded-triangle equivalents. */
    size_t activeRenderCost(void) const;
    /** O(1) CAD-only portion of activeRenderCost(). */
    size_t activeCadRenderCost(void) const;
    /* Irreducible retained PoP population in the same weighted units as
     * activeRenderCost().  This is the cost after every currently displayed
     * progressive occurrence has been retargeted to its minimum drawable
     * prefix; it lets the controller distinguish reducible triangle detail
     * from a genuine per-visible-object floor before aggregating objects into
     * point proxies. */
    size_t minimumActiveRenderCost(void) const;
    /** O(1) CAD-only portion of minimumActiveRenderCost(). */
    size_t minimumActiveCadRenderCost(void) const;
    /* Actual retained CAD primitive submissions from the most recently
     * completed render.  Returns FALSE when any active CAD presentation did
     * not publish an exact completed-frame work record. */
    SbBool lastCadPresentedPrimitiveCount(size_t &primitives) const;
    /* Exact submitted work translated into the same weighted units as
     * activeRenderCost and scene admission.  Unlike a triangle-ratio
     * estimate, this preserves shaded and wire population costs. */
    SbBool lastCadPresentedRenderCost(size_t &cost) const;
    /** TRUE when every retained CAD assembly completed all drawing channels
     * for the latest traversal.  Deadline-stopped point/wire-only buffers are
     * coherent GL results, but are not complete presentation frames. */
    SbBool lastCadPresentationFrameExact(void) const;
    /** Exact presentation classification for the newest completed CAD frame.
     * Mesh occurrences, subpixel point surrogates, and ordinary structural
     * boxes are deliberately separate: only the first two can satisfy a
     * successful current-view coverage proof. */
    SbBool lastCadPresentationOccurrenceCoverage(
	size_t &subpixelOccurrences, size_t &structuralBoxes) const;
    /* Latest completed asynchronous GPU timer aggregate.  The serial changes
     * only when at least one retained CAD context publishes a newer result. */
    SbBool lastCadGpuMeasurement(size_t &faces,
	uint64_t &nanoseconds, uint64_t &serial,
	float &pointProxyPixelThreshold) const;
    /** TRUE when this view currently owns at least one retained CAD
     * presentation.  A completed zero-work CAD frame is not evidence for
     * calibrating the still-active occurrence population. */
    SbBool hasCadPresentationAssemblies(void) const;
    /** Aggregate monotonic token for retained CAD draw attempts.  Comparing
     * this value around one render distinguishes deadline spent preparing a
     * resumable plan from deadline spent submitting actual CAD geometry. */
    uint64_t cadPresentationExecutionSerial(void) const;
    /** Aggregate monotonic token for retained CAD preparation work.  A
     * deadline abort which changes this token is not a steady draw-capacity
     * sample. */
    uint64_t cadPresentationPreparationSerial(void) const;
    /** Snapshot the retained CAD assemblies immediately before one endpoint
     * presentation traversal.  The matching refresh call can then reject a
     * stale exact work record when a deadline stops traversal between CAD
     * assemblies or before an assembly reaches its renderer. */
    void beginCadPresentationFrame(void) const;
    /** Sample every unique active CAD assembly once after a complete frame.
     * The individual presentation, timer, replay, and resource queries below
     * then read one retained aggregate instead of independently walking the
     * source-presentation map. */
    void refreshCadPresentationFrameStatus(void) const;
    /** Aggregate exact Obol-owned buffer accounting from unique active CAD
     * presentations.  Once the completed frame has been sampled this query
     * is O(1), never O(parts) or O(source presentations). */
    SbBool cadGpuResourceStatus(CadGpuResourceStatus &status) const;
    /* TRUE only when every active retained CAD presentation used a reusable
     * indirect command record or flattened transformed atlas in the most
     * recent frame. */
    SbBool lastCadPresentationUsedPreparedReplay(void) const;
    int maximumActiveProgressiveCut(void) const;
    /** For exactly one retained progressive CAD occurrence, return the
     * richest active cut whose cached submitted population fits primitives.
     * The query is O(cuts), performs no realization or allocation, and
     * returns -1 when the scene is not the single-occurrence case. */
    int singleCadProgressiveCutWithinPrimitiveCount(
	size_t primitives) const;
    /* Apply an O(1)-per-assembly render-only ceiling while the precise
     * occurrence allocator catches up with an interactive view. */
    void setCadPresentationProgressiveCutCeiling(int cut) const;
    /* Collapse eligible micro-geometry into one point batch using the current
     * view's measured screen-error tolerance.  One pixel is the stable,
     * pixel-exact setting. */
    void setCadPresentationPointProxyPixelThreshold(float pixels) const;
    /* Permit bounded reuse of the last exact camera-dependent CAD submission
     * during rotation/translation.  Zoom and quiet frames must disable it. */
    void setCadPresentationCameraMotionFrameReuse(SbBool enabled) const;
    size_t estimateDisplayMeshBytes(void) const;
    /* Report the richest already-resident, view-requested prefix of each
     * retained asset.  Presentation admission may draw a coarser activeCut
     * to meet the frame budget, but that must not discard useful data which
     * still fits the resident-memory budget.  A lower request after zoom-out
     * and explicit memory-pressure recovery remain valid compaction edges. */
    void residentMeshDemands(
	std::vector<BObolLodResidentDemand> &demands) const;
    uint64_t residentMeshDemandRevision(void) const;
    /* Adopt a compacted immutable asset generation and journal one exact
     * occurrence per affected source assembly.  Since every occurrence of
     * the asset references the same retained part, one source-local
     * publication swaps its geometry without scanning all instances. */
    size_t applyResidentMeshCompaction(
	const BObolLodResidentCompaction &result);
    uint64_t cadRevision(void) const;
    /* Compact presentations consume occurrence-local changes without
     * rebuilding every leaf.  fullResync is TRUE when the authoritative
     * source state, rather than the returned keys, must be scanned. */
    void cadOccurrenceChangesSince(const SoBRLDatabaseSource *source,
	uint64_t revision, std::vector<SbString> &occurrenceKeys,
	SbBool &fullResync) const;
    void acknowledgeCadOccurrenceChanges(
	const SoBRLDatabaseSource *source, uint64_t revision) const;
    void noteResidentMeshesChanged(const char *reason);
    void setNormalStyle(NormalStyle style, float creaseAngleDegrees = 60.0f);
    NormalStyle getNormalStyle(void) const;
    float getNormalCreaseAngle(void) const;
    /* Drop only mesh/full-detail display data.  Coarse proxy bindings remain
     * resident so memory pressure degrades detail without emptying a scene. */
    size_t evictDisplayMeshPayloads(unsigned int *evictedMeshCount = NULL);
    size_t evictDisplayMeshes(unsigned int *evictedMeshCount = NULL);

    /* Presentation nodes are view-owned, not payload-owned.  A stable node
     * lets an active frame remain drawable while individual LoD occurrences
     * are replaced in place. */
    SoCADAssembly *findCadPresentation(const SoBRLDatabaseSource *source,
	SbString *contentKey = NULL) const;
    void setCadPresentation(const SoBRLDatabaseSource *source,
	SoCADAssembly *assembly, const SbString &contentKey = SbString("")) const;

private:
    struct CadPresentation {
	CadPresentation(void) : assembly(NULL), contentKey(""), sourceRoutingId(0) {}
	SoCADAssembly *assembly;
	SbString contentKey;
	/* Presentation nodes are reusable across view epochs, but never across
	 * database-source object lifetimes.  Semantic keys deliberately survive
	 * erase/redraw for payload reuse, so pointer/revision equality is not a
	 * sufficient lifetime test. */
	uint64_t sourceRoutingId;
    };

    void clearCadPresentations(void) const;
    SbBool applyMeshResultInternal(const SoBRLMeshShape *shape,
	BObolLodResult &result, SbBool consume);
    SbBool applyProxyResultInternal(const SoBRLMeshShape *shape,
	BObolLodResult &result, SbBool consume);
    SbBool applySourceResultInternal(const SoBRLDatabaseSource *source,
	BObolLodResult &result, SbBool consume);
    std::unordered_map<std::string, MeshPayloadPtr> meshBindings;
    std::unordered_map<std::string, ProxyPayloadPtr> proxyBindings;
    /* One authoritative payload per source/occurrence.  cadBindings is only
     * for source-wide legacy shape/path/name lookups; compact occurrences
     * resolve directly through cadSourceBindings and never allocate a second
     * string-keyed alias. */
    std::unordered_map<std::string,
	std::unordered_map<std::string, CadPayloadPtr> > cadSourceBindings;
    std::unordered_map<std::string,
	std::unordered_map<std::string, CadPayloadPtr> > cadAssetBindings;
    std::unordered_map<std::string,
	std::unordered_set<CadPayload *> > cadPayloadsByAssetKey;
    std::unordered_map<std::string, CadPayloadPtr> cadBindings;
    struct CadOccurrenceFailure {
	uint64_t databaseRevision = 0;
	uint64_t sourceRevision = 0;
	uint64_t sourceContentHash = 0;
	uint64_t viewRevision = 0;
	uint64_t policyRevision = 0;
	int requestedCut = -1;
	int drawMode = BOBOL_LOD_DRAW_UNKNOWN;
	int qualityTier = BOBOL_LOD_QUALITY_METADATA;
	int providerStatus = BOBOL_LOD_PROVIDER_UNKNOWN;
    };
    std::unordered_map<std::string,
	std::unordered_map<std::string, CadOccurrenceFailure> >
	cadOccurrenceFailures;
    struct CadOccurrenceChange {
	uint64_t revision;
	SbString occurrenceKey;
    };
    mutable std::unordered_map<std::string,
	std::vector<CadOccurrenceChange> > cadOccurrenceChanges;
    uint64_t cadFullResyncRevision;
    uint64_t cadBindingsRevision;
    NormalStyle normalStyle;
    float normalCreaseAngle;
    mutable int cadPresentationProgressiveLodCeiling;
    mutable float cadPresentationPointProxyPixelThreshold;
    mutable SbBool cadPresentationCameraMotionFrameReuse;
    mutable std::unordered_map<std::string, CadPresentation> cadPresentations;
    /* Multiple semantic source keys may share one retained assembly.  Track
     * that uniqueness at mutation time so completed-frame policy never has
     * to allocate a temporary set or scan duplicate bindings. */
    mutable std::unordered_map<const SoCADAssembly *, size_t>
	cadPresentationAssemblyUseCounts;
    mutable std::unordered_map<const SoCADAssembly *, uint64_t>
	cadPresentationFrameStartExecutionSerials;
    mutable SbBool cadPresentationFrameObservationArmed;
    mutable SbBool cadPresentationFrameStatusValid;
    mutable SbBool cadLastPresentedPrimitiveCountValid;
    mutable size_t cadLastPresentedPrimitiveCount;
    mutable SbBool cadLastPresentedRenderCostValid;
    mutable size_t cadLastPresentedRenderCost;
    mutable SbBool cadLastPresentationFrameExact;
    mutable size_t cadLastSubpixelProxyCount;
    mutable size_t cadLastUncollapsedStructuralProxyCount;
    mutable SbBool cadLastGpuMeasurementValid;
    mutable size_t cadLastGpuFaces;
    mutable uint64_t cadLastGpuNanoseconds;
    mutable uint64_t cadLastGpuSerial;
    mutable float cadLastGpuPointProxyPixelThreshold;
    mutable SbBool cadLastPreparedReplay;
    mutable SbBool cadGpuResourceStatusValid;
    mutable CadGpuResourceStatus cadGpuResourceStatusValue;
    /* Authoritative CAD occurrence telemetry is maintained at mutation
     * points.  These values are read on the presentation/UI thread several
     * times per frame; scanning and locking every resident progressive mesh
     * made nominally cheap HUD and scene-budget queries O(scene size). */
    size_t cadValidPayloadCount;
    size_t cadMeshPayloadCountValue;
    std::unordered_map<std::string, size_t> cadMeshPayloadCountsBySource;
    std::unordered_map<std::string,
	std::unordered_set<std::string>>
	cadUnsatisfiedOccurrencesBySource;
    size_t cadProxyPayloadCountValue;
    size_t cadSatisfiedMeshPayloadCount;
    size_t cadMemoryLimitedMeshPayloadCount;
    size_t cadActiveFaceCount;
    size_t cadActiveRenderCost;
    size_t cadMinimumActiveRenderCost;
    size_t cadDisplayMeshBytes;
    size_t cadProxyKindCounts[5];
    size_t cadProgressiveCutCounts[BOBOL_MESH_LOD_CUT_COUNT_MAX];
    /*
     * Retain the resident view-demand aggregate at the same mutation
     * points as the other payload telemetry.  Reconstructing this table from
     * tens of thousands of occurrences during a paint/quiet compaction made
     * an otherwise unrelated selection event stall the GUI thread.
     */
    struct CadResidentDemandState {
	size_t cutCounts[BOBOL_MESH_LOD_CUT_COUNT_MAX] = {};
	size_t channelCounts[4] = {};
	std::unordered_map<uint32_t, size_t> chunkCounts;
	size_t demandIndex = 0;
	int maximumCut = -1;
	unsigned int channelMask = 0;
    };
    std::unordered_map<std::string, CadResidentDemandState>
	cadResidentDemandStates;
    std::vector<BObolLodResidentDemand> cadResidentDemands;
    uint64_t cadResidentDemandRevision;
    uint64_t cadActiveAllocationPlanSerial;
    uint64_t cadNextAllocationPlanSerial;
    void clearCadPayloadMetrics(void);
    void addCadPayloadMetrics(const CadPayload *payload);
    void removeCadPayloadMetrics(const CadPayload *payload);
    void addCadResidentDemand(const CadPayload *payload);
    void removeCadResidentDemand(const CadPayload *payload);
    void noteCadOccurrenceChanged(const std::string &sourceBindingKey,
	const SbString &occurrenceKey);
};

class BOBOL_EXPORT SoBRLViewLodElement : public SoElement
{
    typedef SoElement inherited;

    SO_ELEMENT_HEADER(SoBRLViewLodElement);

public:
    static void initClass(void);

    virtual void init(SoState *state);
    virtual void push(SoState *state);
    virtual SbBool matches(const SoElement *element) const;
    virtual SoElement *copyMatchInfo(void) const;

    static void set(SoState *state,
		    SoNode *node,
		    const BObolViewLodState *viewState);
    static const BObolViewLodState *get(SoState *state);

protected:
    virtual ~SoBRLViewLodElement(void);

private:
    const BObolViewLodState *viewState;
};

class BOBOL_EXPORT SoBRLViewLodGroup : public SoGroup
{
    typedef SoGroup inherited;

    SO_NODE_HEADER(SoBRLViewLodGroup);

public:
    SoBRLViewLodGroup(void);
    static void initClass(void);

    void setViewLodState(BObolViewLodState *viewState);
    BObolViewLodState *getViewLodState(void) const;
    void setSoftwareWireMode(int mode);
    int getSoftwareWireMode(void) const;

    virtual void doAction(SoAction *action) override;
    virtual void GLRender(SoGLRenderAction *action) override;
    virtual void callback(SoCallbackAction *action) override;
    virtual void getBoundingBox(SoGetBoundingBoxAction *action) override;
    virtual void pick(SoPickAction *action) override;

protected:
    virtual ~SoBRLViewLodGroup(void);

private:
    SbBool pushViewState(SoAction *action);
    void popViewState(SoAction *action, SbBool pushed);

    BObolViewLodState *viewState;
    int softwareWireMode;
};

BOBOL_EXPORT const BObolViewLodState::MeshPayload *
bobol_view_lod_mesh_for_action(SoAction *action,
				 const SoBRLMeshShape *shape);

BOBOL_EXPORT const BObolViewLodState::ProxyPayload *
bobol_view_lod_proxy_for_action(SoAction *action,
				  const SoBRLMeshShape *shape);

BOBOL_EXPORT const BObolViewLodState::CadPayload *
bobol_view_lod_cad_for_action(SoAction *action,
				const SoBRLDatabaseSource *source);

BOBOL_EXPORT const BObolViewLodState *
bobol_view_lod_state_for_action(SoAction *action);

#endif /* BOBOL_BVIEWLOD_H */
