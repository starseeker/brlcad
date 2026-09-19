/*            B S C E N E C O N T R O L L E R . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BSceneController.h */

#ifndef BOBOL_BSCENECONTROLLER_H
#define BOBOL_BSCENECONTROLLER_H

#include "BObol/BDefines.h"
#include "BObol/BDatabaseSource.h"

#include <Inventor/SbBasic.h>
#include <Inventor/SbBox.h>
#include <Inventor/SbColor.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec3f.h>

#include <stdint.h>
#include <memory>
#include <vector>

class SoNode;
class SoGroup;
class SoBRLDatabaseSource;
class BObolRealizationRepository;
struct BObolSourceRealizationEffects;
struct BObolAuxiliaryLineSetDisplayState;
struct BObolDatabaseSourceDisplayPatch;
struct BObolDrawMetadataRecord;
struct BObolExternalAnnotation;
struct BObolExternalLineSet;
struct BObolExternalPointSet;
struct BObolExternalTriangleMesh;
struct db_i;
struct BObolDatabaseSourceSummary;
struct BObolRealizedMaterialSummary;
struct BObolRealizedShapeSummary;
struct BObolSceneBoundsSummary;
struct BObolSceneDisplaySummary;
struct BObolSceneMaterialSummary;
struct BObolSceneTreeSummary;

struct BOBOL_EXPORT BObolSceneSummary {
    BObolSceneSummary(void);

    SbBool valid;
    SbBool hasRoot;
    SbBool rootIsGroup;
    int rootChildCount;
    int databaseSourceCount;
    int nonDatabaseRootChildCount;
    uint64_t structuralRevision;
    uint64_t frameRevision;
    unsigned int lastVisitedSourceCount;
    unsigned int lastRealizedSourceCount;
    unsigned int lastFailedSourceCount;
    SbString lastDiagnostics;
};

/* Metadata for the target group of a source publication. The scene commits
 * it with the source and membership before either node notifies observers. */
struct BOBOL_EXPORT BObolSceneGroupPublishState {
    BObolSceneGroupPublishState(void);

    const char *intentPath;
    int drawMode;
    int fallbackDrawMode;
    SbBool overlayIntent;
    uint32_t revalidationRevision;
    SbBool visible;
    SbBool selected;
    SbBool highlighted;
    int lineStyle;
    int lineWidth;
    float transparency;
    SbBool colorOverride;
    SbColor color;
    SbBool materialColorValid;
    SbColor materialColor;
    uint32_t materialRevision;
};

struct BOBOL_EXPORT BObolSceneSourcePresentationTarget {
    SbString sourceInstanceKey;
    BObolDatabaseSourcePresentationPatch presentation;
    /* When valid, the target must still be the exact compact-presentation
     * revision accepted before an earlier publication notified observers. */
    BObolSourcePresentationStamp expectedStamp;
};

struct BOBOL_EXPORT BObolSceneGroupPresentationTarget {
    SbString groupPath;
    BObolDatabaseSourceDisplayPatch display;
};

struct BOBOL_EXPORT BObolScenePresentationTransaction {
    std::vector<BObolSceneSourcePresentationTarget> sources;
    std::vector<BObolSceneGroupPresentationTarget> groups;
};

/* Scene edges retired by one structural publication.  Descendant targets may
 * be included; the controller reduces them beneath the selected ancestor. */
struct BOBOL_EXPORT BObolSceneRemovalTransaction {
    std::vector<SbString> sourceInstanceKeys;
    std::vector<SbString> groupPaths;
};

class BOBOL_EXPORT BObolSceneController {
public:
    BObolSceneController(void);
    /**
     * Create a controller over an application-owned Obol scene root.
     *
     * The controller retains the root with normal Obol reference counting, but
     * it does not own an authoritative hierarchy separate from Obol. After
     * attachment, publish database-source and named-group membership changes
     * through this controller so its indexes, repository ownership and scene
     * revisions remain one transaction. Direct edits to shape data and compact
     * source journals retain their documented source-level notification paths.
     */
    explicit BObolSceneController(SoNode *root);
    ~BObolSceneController(void);

    /**
     * Replace the retained scene root used by subsequent realization passes.
     */
    void setSceneRoot(SoNode *root);
    SoNode *getSceneRoot(void) const;
    uint64_t getStructuralRevision(void) const;
    uint64_t getFrameRevision(void) const;
    SbBool getSceneSummary(BObolSceneSummary &summary) const;

    /* Each changed source advances its scenes' frames before observers.
     * A callback that changes this scene's hierarchy or repository stops the
     * enclosing traversals; a subsequent call visits the current scene.
     * Completed progress survives interruption, and mutation batches retain
     * their explicit revision coalescing semantics. */
    SbBool realizePending(void);
    /* Realize only the source accepted by this key and worker stamp. A stale,
     * reconfigured or same-key replacement target is left untouched. */
    SbBool realizeDatabaseSourceInstance(const char *sourceInstanceKey,
	const BObolSourceRealizationStamp &stamp);
    /* Atomically transfer a prepared mesh-LoD reader and its bounds to the
     * exact source accepted before preparation.  On rejection the caller
     * retains ownership of lod.  This private runtime resource does not
     * advance scene presentation revisions. */
    SbBool adoptDatabaseSourceInstanceMeshLod(
	const char *sourceInstanceKey,
	const BObolSourceRealizationStamp &stamp,
	struct BObolMeshLod *lod,
	const SbVec3f &bmin,
	const SbVec3f &bmax);
    void shareRealizationRepository(BObolSceneController *source);
    /* Repository mutations stop an active realization pass. Retired cache
     * values stay alive until repository bookkeeping is complete, so their
     * destruction callbacks observe the committed state and later edits win. */
    void clearRealizationRepository(void);
    void invalidateRealizationViewVariants(void);
    void renameRealizationObject(const char *oldName, const char *newName);
    void beginSceneMutationBatch(size_t expectedDatabaseSources = 0,
	size_t expectedGroups = 0);
    void endSceneMutationBatch(void);

    SoGroup *findGroup(const char *groupPath) const;
    /* Publish the complete missing hierarchy, indexes and revisions before
     * observers. Return the group currently at the requested path afterward,
     * or null if a callback removed or renamed it. */
    SoGroup *ensureGroup(const char *groupPath);
    int setGroupDrawIntent(const char *groupPath,
	const char *intentPath,
	int drawMode,
	int fallbackDrawMode,
	SbBool overlayIntent,
	uint32_t revalidationRevision);
    int setGroupDisplayState(const char *groupPath,
	SbBool visible,
	SbBool selected,
	SbBool highlighted,
	int lineStyle,
	int lineWidth,
	float transparency,
	SbBool colorOverride,
	const SbColor &color,
	SbBool materialColorValid,
	const SbColor &materialColor,
	uint32_t materialRevision);
    /* Publish the name, descendant paths, group indexes and scene revisions
     * together before observers. Preserve source state and later callback edits. */
    int renameGroup(const char *groupPath, const char *newLeafName);
    /* Publish complete subtree membership, indexes, repository ownership and
     * revisions before observers; retain nodes still reachable through other edges. */
    int appendChildToGroup(const char *groupPath, SoNode *child);
    int removeChildFromGroup(const char *groupPath, SoNode *child);
    /* Apply each removal or clear as one complete child publication. Counts
     * describe this operation even when an observer subsequently changes it. */
    int eraseGroupSubpath(const char *parentGroupPath,
	const char *subpath);
    int removeGroup(const char *groupPath);
    int clearGroup(const char *groupPath);
    int getGroupChildCount(const char *groupPath) const;
    int getGroupDescendantGroupCount(const char *groupPath) const;
    int getGroupDatabaseSourceCount(const char *groupPath) const;
    SbBool getDatabaseSourceBounds(SbBox3f &bounds,
	SbBool padForAutoview) const;
    SbBool getSceneSubtreeBounds(const char *nodePath,
	SbBool includeOverlays,
	SbBox3f &bounds) const;

    SoNode *findShape(const char *shapePath) const;
    SoGroup *findShapeParent(const char *shapePath) const;
    /* Publish shape membership, child paths and scene revisions together before
     * observers. Preserve shape data, shared geometry and later callback edits. */
    int moveShapeToGroup(const char *shapePath, const char *groupPath);
    int removeShape(const char *shapePath);
    /* Shape-state setters prepare requested fields and commit their frame effect
     * before observers. Later callback edits remain current. */
    int setShapeDrawState(const char *shapePath,
	int drawMode,
	SbBool databaseIntent,
	SbBool overlayIntent,
	SbBool hudIntent);
    int setShapeDisplayState(const char *shapePath,
	SbBool visible,
	SbBool selected,
	SbBool highlighted,
	int lineStyle,
	int lineWidth,
	float transparency,
	SbBool colorOverride,
	const SbColor &color,
	SbBool materialColorValid,
	const SbColor &materialColor,
	uint32_t materialRevision);
    int setShapeSourceState(const char *shapePath,
	const char *ownerSourcePath,
	uint32_t ownerSourceRevision,
	uint32_t ownerInputsRevision,
	uint32_t ownerViewRevision,
	uint32_t ownerRealizedRevision,
	uint32_t ownerRealizedSourceRevision,
	uint32_t ownerRealizedInputsRevision,
	uint32_t ownerRealizedViewRevision,
	int ownerRealizationStatus,
	const char *ownerRealizationDiagnostic,
	const char *ownerRealizationIdentity,
	SbBool ownerSourceStale,
	uint32_t ownerStaleReason);
    int setShapePlacementState(const char *shapePath,
	SbBool drawMatrixValid,
	const SbMatrix &drawMatrix,
	SbBool drawCenterValid,
	const SbVec3f &drawCenter,
	SbBool drawSizeValid,
	float drawSize);
    int publishDatabaseSourceAuxiliaryLineSet(const char *sourcePath,
	const char *name,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	const BObolAuxiliaryLineSetDisplayState *displayState = NULL);
    int publishDatabaseSourceAuxiliarySourceLineSet(
	const char *sourcePath,
	const char *auxiliarySourcePath,
	const char *displayName,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	const BObolAuxiliaryLineSetDisplayState *displayState = NULL);
    int publishDatabaseSourceInstanceAuxiliaryLineSet(
	const char *sourceInstanceKey,
	const char *name,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	const BObolAuxiliaryLineSetDisplayState *displayState = NULL);
    int publishDatabaseSourceInstanceAuxiliarySourceLineSet(
	const char *sourceInstanceKey,
	const char *auxiliarySourcePath,
	const char *displayName,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	const BObolAuxiliaryLineSetDisplayState *displayState = NULL);
    /* External primary geometry is source-local.  The database source's
     * placement state carries local-to-scene transforms separately.
     */
    int publishDatabaseSourceExternalLineSet(const char *sourcePath,
	const BObolExternalLineSet &lineSet);
    int publishDatabaseSourceInstanceExternalLineSet(
	const char *sourceInstanceKey,
	const BObolExternalLineSet &lineSet);
    int publishDatabaseSourceInstancePrimitiveWireframe(
	const char *sourceInstanceKey,
	struct rt_db_internal *intern,
	const struct bg_tess_tol *ttol = NULL,
	const struct bn_tol *tol = NULL);
    int publishDatabaseSourceExternalPointSet(const char *sourcePath,
	const BObolExternalPointSet &pointSet);
    int publishDatabaseSourceInstanceExternalPointSet(
	const char *sourceInstanceKey,
	const BObolExternalPointSet &pointSet);
    int publishDatabaseSourceExternalTriangleMesh(const char *sourcePath,
	const BObolExternalTriangleMesh &triangleMesh);
    int publishDatabaseSourceInstanceExternalTriangleMesh(
	const char *sourceInstanceKey,
	const BObolExternalTriangleMesh &triangleMesh);
    int publishDatabaseSourceExternalAnnotation(const char *sourcePath,
	const BObolExternalAnnotation &annotation);
    int publishDatabaseSourceInstanceExternalAnnotation(
	const char *sourceInstanceKey,
	const BObolExternalAnnotation &annotation);
    int clearDatabaseSourceExternalPrimaryGeometry(const char *sourcePath);
    int clearDatabaseSourceInstanceExternalPrimaryGeometry(
	const char *sourceInstanceKey);
    int clearDatabaseSourceAuxiliaryShapes(const char *sourcePath);
    int clearDatabaseSourceInstanceAuxiliaryShapes(
	const char *sourceInstanceKey);

    SoBRLDatabaseSource *getDatabaseSource(int index) const;
    int getDatabaseSourceCount(void) const;
    /* Constant-time owner lookup for asynchronous compact-LoD results.
     * The routing id is an in-process lifetime token, not scene or cache
     * identity.  Rebuilds occur only after structural scene mutations. */
    SoBRLDatabaseSource *findDatabaseSourceRoutingId(uint64_t routingId) const;
    SoBRLDatabaseSource *findDatabaseSource(const char *sourcePath) const;
    SoBRLDatabaseSource *findDatabaseSourceInstance(
	const char *sourceInstanceKey) const;
    int replaceDatabaseSource(const char *sourcePath,
	struct db_i *database,
	int drawMode,
	uint32_t sourceRevision);
    int replaceDatabaseSourceInstance(const char *sourceInstanceKey,
	const char *sourcePath,
	struct db_i *database,
	int drawMode,
	uint32_t sourceRevision);
    int replaceDatabaseSourceInstanceRepresentation(
	const char *sourceInstanceKey,
	const char *sourcePath,
	const char *sourceRepresentationKey,
	int sourceRepresentationMode,
	struct db_i *database,
	int drawMode,
	uint32_t sourceRevision);
    int publishDatabaseSourceInstance(
	const BObolDatabaseSourcePublishState &state);
    /* A non-null group state requires a database, a nonempty targetGroupPath
     * and a SoBRLSceneGroup target. Invalid requests return -1 without publication. */
    int publishDatabaseSourceInstance(
	const BObolDatabaseSourcePublishState &state,
	const BObolSceneGroupPublishState *groupState);
    /* Publish an existing instance under state.sourceInstanceKey. Identity,
     * source state, hierarchy, indexes and revisions commit before observers.
     * The destination key must be unowned or already identify that instance. */
    int publishDatabaseSourceInstance(
	const char *currentSourceInstanceKey,
	const BObolDatabaseSourcePublishState &state,
	const BObolSceneGroupPublishState *groupState);
    /* Share the immutable resident occurrences of the specified sources with
     * the target and retire their selected scene edges in one publication.
     * Preparation failure preserves every source. Callback edits, replacement
     * and removal of the target remain current. Returns retired edge count. */
    int subsumeDatabaseSourceInstances(
	const char *targetSourceInstanceKey,
	const char *const *sourceInstanceKeys,
	size_t sourceInstanceCount);
    /* Publish a completed detached realization through its attached scene.
     * The stamped source, terminal fields, child/index effects and frame
     * revision commit before observers. A stale or replaced target is ignored.
     * publicationCommitted distinguishes preparation failure from a throwing
     * observer after the complete commit. */
    int adoptDatabaseSourceInstanceRealization(
	const char *sourceInstanceKey,
	const BObolSourceRealizationStamp &stamp,
	SoBRLDatabaseSource *detached,
	SbBool authoritativeStreamDrained = FALSE,
	const std::shared_ptr<BObolCompactOccurrenceStream> &stagedSourceStream =
	    std::shared_ptr<BObolCompactOccurrenceStream>(),
	SbBool *publicationCommitted = nullptr);
    /* Merge one producer batch into the exact stamped scene source. Source
     * changes and all affected frame revisions commit before field observers.
     * batchCompleted distinguishes a merge failure from a throwing observer
     * after the entire requested batch was processed. */
    int mergeDatabaseSourceInstanceCompactOccurrences(
	const char *sourceInstanceKey,
	const BObolSourceRealizationStamp &stamp,
	const std::vector<BObolCompactOccurrence> &occurrences,
	SbBool authoritativeGeometry = FALSE,
	SbBool *batchCompleted = nullptr,
	size_t reserveCapacity = 0);
    /* Publish immutable discovery facts for one exact source/producer epoch.
     * The expected count is certified without allocating occurrence storage;
     * a supplied profile must describe that same population. This private
     * planning state does not advance scene revisions or notify fields. */
    int certifyDatabaseSourceInstanceCompactStream(
	const char *sourceInstanceKey,
	const BObolSourceRealizationStamp &stamp,
	size_t expectedCount,
	const BObolCompactSourceProfile *profile = nullptr);
    /* Replace the exact stamped source with one complete compact snapshot.
     * The registry, optional producer-certified bounds/profile, terminal
     * realization state, primary children and scene revisions commit before
     * observers. Existing auxiliary children are retained. A stale or
     * replaced target is ignored. publicationCommitted distinguishes a
     * preparation failure from a throwing observer after the complete commit. */
    int publishDatabaseSourceInstanceCompactSnapshot(
	const char *sourceInstanceKey,
	const BObolSourceRealizationStamp &stamp,
	const std::vector<BObolCompactOccurrence> &occurrences,
	const SbBox3f *certifiedBounds = nullptr,
	const BObolCompactSourceProfile *profile = nullptr,
	SbBool *publicationCommitted = nullptr);
    int renameDatabaseSource(const char *sourcePath,
	const char *newSourcePath,
	uint32_t sourceRevision);
    /* Rename publishes source metadata, affected indexes and revisions before
     * observers. Return -1 if replacing a conflict would remove this source
     * or leave another occurrence of the conflicting instance reachable. */
    int renameDatabaseSourceInstance(const char *sourceInstanceKey,
	const char *newSourceInstanceKey,
	const char *newSourcePath,
	uint32_t sourceRevision);
    /* Compose one database object rename across cache residency, every source
     * path, compact occurrence semantics, group paths and scene indexes.  A
     * key listed in pathDerivedSourceInstanceKeys is retargeted with its
     * source path; unlisted keys remain stable.  All changes and revisions
     * commit before observers. */
    int renameDatabaseObject(const char *oldObjectPath,
	const char *newObjectPath,
	const std::vector<SbString> &pathDerivedSourceInstanceKeys,
	uint32_t sourceRevision = 0);
    int setDatabaseSourceState(const char *sourcePath,
	SbBool sourceRevisionValid,
	uint32_t sourceRevision,
	uint32_t inputsRevision,
	SbBool visible,
	SbBool selected,
	SbBool highlighted,
	int lineStyle,
	int lineWidth,
	float transparency,
	SbBool colorOverride,
	const SbColor &color,
	SbBool materialColorValid,
	const SbColor &materialColor,
	uint32_t materialRevision);
    int setDatabaseSourceInstanceState(const char *sourceInstanceKey,
	SbBool sourceRevisionValid,
	uint32_t sourceRevision,
	uint32_t inputsRevision,
	SbBool visible,
	SbBool selected,
	SbBool highlighted,
	int lineStyle,
	int lineWidth,
	float transparency,
	SbBool colorOverride,
	const SbColor &color,
	SbBool materialColorValid,
	const SbColor &materialColor,
	uint32_t materialRevision);
    int setDatabaseSourceDisplayPatch(const char *sourcePath,
	const BObolDatabaseSourceDisplayPatch &patch);
    int setDatabaseSourceInstanceDisplayPatch(const char *sourceInstanceKey,
	const BObolDatabaseSourceDisplayPatch &patch);
    /* Snapshot every target before publication.  Each changed source or group
     * commits as one complete presentation record before observers.  A valid
     * expectedStamp also rejects work superseded before this call.  Callback
     * edits, removals or replacements supersede pending work for that target;
     * allocation or observer failure leaves a complete prefix. */
    int applyPresentationTransaction(
	const BObolScenePresentationTransaction &transaction);
    int setDatabaseSourceDisplayName(const char *sourcePath,
	const char *displayName);
    int setDatabaseSourceInstanceDisplayName(const char *sourceInstanceKey,
	const char *displayName);
    int setDatabaseSourceDrawMode(const char *sourcePath,
	int drawMode);
    int setDatabaseSourceInstanceDrawMode(const char *sourceInstanceKey,
	int drawMode);
    int setDatabaseSourceInstanceRepresentation(
	const char *sourceInstanceKey,
	const char *sourceRepresentationKey,
	int sourceRepresentationMode);
    int setDatabaseSourceMaterialPolicy(const char *sourcePath,
	int materialPolicy);
    int setDatabaseSourceInstanceMaterialPolicy(const char *sourceInstanceKey,
	int materialPolicy);
    /* Publish each changed source and its frame revision before observers.
     * Failures retain the completed prefix; later callback edits or source
     * replacements supersede the remaining original targets. */
    int setDatabaseSourcesEvaluatedRegionForPath(SoBRLDatabaseSource *const *sources,
	size_t count, const char *path, SbBool evaluated);
    /* A path applies metadata to matching compact occurrences; NULL updates
     * aggregate source metadata. Noncompact sources retain aggregate semantics.
     * Return -1 for invalid/missing targets, 0 for no change, or a positive count. */
    int applyDatabaseSourceInstanceDrawMetadata(const char *sourceInstanceKey,
	const BObolDrawMetadataRecord &metadata, const char *path = nullptr);
    int refreshDatabaseSourceInstanceMaterialColorFromDatabase(
	const char *sourceInstanceKey,
	uint32_t materialRevision,
	struct db_i *overrideDbip = NULL);
    /* Full-scene refresh retains the first source's database as its fallback. */
    int refreshDatabaseSourceMaterialColorsFromDatabase(
	uint32_t materialRevision,
	struct db_i *overrideDbip = NULL);
    /* Retain original owners and publish each changed source's frame effect
     * before observers. Later edits/replacements supersede pending targets.
     * Without an override, each target uses its own database. */
    int refreshDatabaseSourcesMaterialColorsFromDatabase(SoBRLDatabaseSource *const *sources,
	size_t count, uint32_t materialRevision, struct db_i *overrideDbip = NULL);
    int setDatabaseSourcePlacementState(const char *sourcePath,
	SbBool drawMatrixValid,
	const SbMatrix &drawMatrix,
	SbBool drawCenterValid,
	const SbVec3f &drawCenter,
	SbBool drawSizeValid,
	float drawSize);
    int setDatabaseSourceInstancePlacementState(const char *sourceInstanceKey,
	SbBool drawMatrixValid,
	const SbMatrix &drawMatrix,
	SbBool drawCenterValid,
	const SbVec3f &drawCenter,
	SbBool drawSizeValid,
	float drawSize);
    int setDatabaseSourceInstanceHierarchyState(const char *sourceInstanceKey,
	const char *parentInstanceKey,
	uint32_t occurrenceIndex,
	int booleanOperation);
    int setDatabaseSourceBoundsState(const char *sourcePath,
	SbBool boundsValid,
	const SbVec3f &boundsMin,
	const SbVec3f &boundsMax,
	SbBool boundsExact = FALSE);
    /* publicationCommitted becomes true after bounds and affected frame
     * revisions commit, before observers. */
    int setDatabaseSourceInstanceBoundsState(const char *sourceInstanceKey,
	SbBool boundsValid,
	const SbVec3f &boundsMin,
	const SbVec3f &boundsMax,
	SbBool boundsExact = FALSE,
	SbBool *publicationCommitted = nullptr);
    int markDatabaseSourceStale(const char *sourcePath,
	uint32_t staleReason);
    int markDatabaseSourceInstanceStale(const char *sourceInstanceKey,
	uint32_t staleReason);
    int refreshDatabaseSourceInstanceObject(const char *sourceInstanceKey,
	const char *objectPath, uint32_t sourceRevision = 0);
    int setDatabaseSourceRealizationState(const char *sourcePath,
	int realizationStatus,
	uint32_t realizedSourceRevision,
	uint32_t realizedInputsRevision,
	uint32_t staleReason,
	const char *diagnostic = NULL);
    int setDatabaseSourceInstanceRealizationState(
	const char *sourceInstanceKey,
	int realizationStatus,
	uint32_t realizedSourceRevision,
	uint32_t realizedInputsRevision,
	uint32_t staleReason,
	const char *diagnostic = NULL);
    /* Publish realization status and role ownership through one frame effect.
     * Both values commit before observers and later callback edits remain current. */
    int setDatabaseSourceInstanceRealizationState(
	const char *sourceInstanceKey,
	int realizationStatus,
	uint32_t realizedSourceRevision,
	uint32_t realizedInputsRevision,
	uint32_t staleReason,
	const char *diagnostic,
	int roleFlags);
    int setDatabaseSourceRealizationRoleFlags(const char *sourcePath,
	int roleFlags);
    int setDatabaseSourceInstanceRealizationRoleFlags(
	const char *sourceInstanceKey,
	int roleFlags);
    int setDatabaseSourceRealizationViewPolicy(const char *sourcePath,
	SbBool viewDependent,
	SbBool csgLodEnabled,
	SbBool meshLodEnabled,
	float viewScale,
	float lodScale,
	int viewWidth,
	int viewHeight,
	uint32_t botThreshold,
	float curveScale,
	float pointScale);
    int setDatabaseSourceInstanceRealizationViewPolicy(
	const char *sourceInstanceKey,
	SbBool viewDependent,
	SbBool csgLodEnabled,
	SbBool meshLodEnabled,
	float viewScale,
	float lodScale,
	int viewWidth,
	int viewHeight,
	uint32_t botThreshold,
	float curveScale,
	float pointScale);
    int moveDatabaseSourceToGroup(const char *sourcePath,
	const char *groupPath);
    /* Move membership without reconfiguring the source. Publish the complete
     * hierarchy, indexes and scene effects before observers; reject cycles
     * through existing or newly created destination groups. */
    int moveDatabaseSourceInstanceToGroup(const char *sourceInstanceKey,
	const char *groupPath);
    /* Publish source-edge removal with complete descendant indexes, repository
     * ownership and scene effects before observers. Path removal keeps its
     * selected node even when another source has the same instance key. */
    int removeDatabaseSource(const char *sourcePath);
    int removeDatabaseSourceInstance(const char *sourceInstanceKey);
    /* Commit all selected source and group edges, indexes, repository
     * ownership and scene revisions before delivering any observer. */
    int applyRemovalTransaction(
	const BObolSceneRemovalTransaction &transaction);
    /* Clear non-auxiliary source edges together. Preserve ordinary nodes and
     * auxiliary source subtrees; return the number of removed parent edges. */
    int clearDatabaseSources(void);
    SbBool getDatabaseSourceSummary(int index,
	BObolDatabaseSourceSummary &summary) const;
    SbBool getDatabaseSourceSummaryForPath(const char *sourcePath,
	BObolDatabaseSourceSummary &summary) const;
    int getDatabaseSourceInstanceCountForPath(const char *sourcePath) const;
    SbBool getDatabaseSourceInstanceSummaryForPath(const char *sourcePath,
	int instanceIndex, BObolDatabaseSourceSummary &summary) const;
    SbBool getDatabaseSourceSummaryForInstance(const char *sourceInstanceKey,
	BObolDatabaseSourceSummary &summary) const;
    int getRealizedShapeSummaryCount(void) const;
    SbBool getRealizedShapeSummary(int index,
	BObolRealizedShapeSummary &summary) const;
    int getRealizedMaterialSummaryCount(void) const;
    SbBool getRealizedMaterialSummary(int index,
	BObolRealizedMaterialSummary &summary) const;
    SbBool getRealizedMaterialProperty(int materialIndex, int propertyIndex,
	SbString &groupOut, SbString &nameOut, SbString &valueOut) const;
    int getSceneTreeSummaryCount(void) const;
    SbBool getSceneTreeSummary(int index,
	BObolSceneTreeSummary &summary) const;
    SbBool getSceneTreeSummaryForPath(const char *nodePath,
	BObolSceneTreeSummary &summary) const;
    /** Resolve a compact occurrence without creating a per-leaf scene node. */
    SbBool getCompactSceneTreeSummaryForPath(const char *nodePath,
	SbBool includeDescendants,
	BObolSceneTreeSummary &summary) const;
    SbBool getSceneChildTreeSummary(const char *nodePath,
	int childIndex,
	BObolSceneTreeSummary &summary) const;
    int getSceneDisplaySummaryCount(void) const;
    SbBool getSceneDisplaySummary(int index,
	BObolSceneDisplaySummary &summary) const;
    int getSceneMaterialSummaryCount(void) const;
    SbBool getSceneMaterialSummary(int index,
	BObolSceneMaterialSummary &summary) const;
    int getSceneBoundsSummaryCount(void) const;
    SbBool getSceneBoundsSummary(int index,
	BObolSceneBoundsSummary &summary) const;

    unsigned int getLastVisitedSourceCount(void) const;
    unsigned int getLastRealizedSourceCount(void) const;
    unsigned int getLastFailedSourceCount(void) const;
    const SbString &getLastDiagnostics(void) const;

private:
    friend class BObolViewController;
    SbBool realizePending(BObolSourceRealizationEffects *effects);
    SbBool realizeSubtree(SoNode *root,
	BObolSourceRealizationEffects *effects,
	const char *compactRootInstanceKey);
    struct Impl;
    class SourceIndexPublication;
    class SourcePublication;
    class GroupPublication;
    class ShapePublication;
    class HierarchyPublication;
    class ChildPublication;
    class SourceChildEffects;
    class RootOwnership;
    class PeerEffects;
    class RootPublication;
    class MembershipPublication;
    static std::unique_ptr<BObolSourceRealizationEffects> prepareRealizationEffects(
	BObolSceneController &scene, BObolSourceRealizationEffects *downstream);

    template <typename Configure>
    int publishShapeState(const char *shapePath, Configure configure);

    int publishDatabaseSourceInstanceImpl(
	const char *currentSourceInstanceKey,
	SbBool requireExisting,
	const BObolDatabaseSourcePublishState &state,
	const BObolSceneGroupPublishState *groupState);

    int removeGroupSubpath(SoGroup *parent, const char *subpath);

    void advanceFrameRevision(void);
    void advanceStructuralRevision(SbBool stopRealizationTraversal = TRUE);
    static bool databaseSourcePublicationAccepted(SoBRLDatabaseSource *source, void *context);
    void clearDatabaseSourceIndex(void) const;
    void rebuildDatabaseSourceIndex(void) const;
    SoGroup *findIndexedGroup(const char *groupPath) const;
    SoBRLDatabaseSource *findIndexedDatabaseSource(
	const char *sourcePath) const;
    SbString databaseSourceInstanceKeyForPath(
	const char *sourcePath) const;
    SoBRLDatabaseSource *findIndexedDatabaseSourceInstance(
	const char *sourceInstanceKey) const;
    SoGroup *findIndexedDatabaseSourceInstanceParent(
	const char *sourceInstanceKey) const;
    SbBool databaseSourceSummaryForSource(SoBRLDatabaseSource *source,
	BObolDatabaseSourceSummary &summary) const;

    std::unique_ptr<Impl> d;
};

#endif /* BOBOL_BSCENECONTROLLER_H */
