/*             B D A T A B A S E S O U R C E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BDatabaseSource.h */

#ifndef BOBOL_BDATABASESOURCE_H
#define BOBOL_BDATABASESOURCE_H

#include "BObol/BDefines.h"
#include "BObol/BSourceMeshRequest.h"
#include "BObol/BViewLod.h"

#include <stdint.h>

#include <Inventor/SbBox.h>
#include <Inventor/SbColor.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFColor.h>
#include <Inventor/fields/SoSFEnum.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFInt32.h>
#include <Inventor/fields/SoSFMatrix.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFUInt32.h>
#include <Inventor/fields/SoSFVec3f.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Obol/cad/CadProgressive.h>

#include <memory>
#include <array>
#include <vector>

class SoBRLVListShape;
class SoBRLMeshShape;
class SoBRLMaterialObject;
class SoBRLCadAssembly;
class SoBRLCadRenderBatch;
class SoBRLExportAction;
class SoBRLMeasureAction;
class SoBRLSnapAction;
struct BObolViewPickRecord;
class SoFieldSensor;
class SoCallbackAction;
class SoGetBoundingBoxAction;
class SoGLRenderAction;
class SoRayPickAction;
class SoSensor;
struct db_i;
struct BObolMeshLod;
struct db_full_path;
struct rt_db_internal;
struct bg_tess_tol;
struct bn_tol;
struct BObolDatabaseSourceRealizationCache;
struct BObolCompactInstanceIndex;
struct BObolCompactInstanceEntry;
struct BObolCadBatchBuildState;
namespace Obol { struct PartGeometry; }

/** Build immutable structural box geometry.  Ordinary AABBs share a unit
 * wire primitive through geometryTransform; oriented bounds retain their
 * object-coordinate corners so Obol can classify and render them as OBBs. */
BOBOL_EXPORT std::shared_ptr<const Obol::PartGeometry>
bobol_cad_structural_bounds_geometry(
    const SbBox3f &bounds,
    SbMatrix &geometryTransform,
    const SbVec3f orientedBounds[8] = nullptr);

BOBOL_EXPORT SbBool bobol_database_source_fullpath_material_color(
	struct db_i *dbip,
	const struct db_full_path *pathp,
	SbColor &color);
BOBOL_EXPORT SbBool bobol_database_source_path_material_color(
	struct db_i *dbip,
	const char *path,
	SbColor &color);

/** Kinds of in-scene light realized from a database "light"-shader region. */
enum BObolSceneLightKind {
    BOBOL_SCENE_LIGHT_POINT = 0,
    BOBOL_SCENE_LIGHT_SPOT = 1,
    BOBOL_SCENE_LIGHT_DIRECTIONAL = 2
};

/** Explicit semantic path matching for compact occurrence updates. */
enum BObolCompactPathMatch {
    BOBOL_COMPACT_PATH_EXACT = 0,
    BOBOL_COMPACT_PATH_SUBTREE,
    BOBOL_COMPACT_PATH_OBJECT
};

/** A light source derived from a database region whose optical shader is
 * "light".  Positions/directions are world-space (the realize walk's
 * accumulated transform is already applied), so consumers need no further
 * transform.  Mirrors the parameters of liboptical's sh_light.c. */
struct BOBOL_EXPORT BObolSceneLightRealization {
    BObolSceneLightRealization(void);

    int kind;             /**< BObolSceneLightKind */
    SbVec3f position;     /**< world position (point/spot) */
    SbVec3f direction;    /**< world aim direction (spot/directional) */
    SbColor color;
    float intensity;      /**< normalized [0,1] */
    float coneAngleDeg;   /**< full beam angle in degrees (spot) */
    SbString name;        /**< region name (for idempotent rebuilds) */
};

struct BOBOL_EXPORT BObolDatabaseSourceSummary {
    BObolDatabaseSourceSummary(void);

    SbBool valid;
    SbString path;
    SbString instanceKey;
    SbString parentInstanceKey;
    uint32_t occurrenceIndex;
    int booleanOperation;
    SbString displayName;
    SbBool hasParent;
    int drawTreeDepth;
    SbString parentGroupPath;
    SbString representationKey;
    int representationMode;
    int drawMode;
    uint32_t sourceRevision;
    uint32_t inputsRevision;
    uint32_t viewRevision;
    uint32_t realizedRevision;
    uint32_t realizedSourceRevision;
    uint32_t realizedInputsRevision;
    uint32_t realizedViewRevision;
    int realizationStatus;
    SbString realizationDiagnostic;
    SbString realizationIdentity;
    int realizationRoleFlags;
    SbBool realizationViewDependent;
    SbBool realizationCsgLodEnabled;
    SbBool realizationMeshLodEnabled;
    float realizationViewScale;
    float realizationLodScale;
    int realizationViewWidth;
    int realizationViewHeight;
    uint32_t realizationBotThreshold;
    float realizationCurveScale;
    float realizationPointScale;
    SbBool visible;
    SbBool selected;
    SbBool highlighted;
    int lineStyle;
    int lineWidth;
    float transparency;
    SbBool materialColorValid;
    SbColor materialColor;
    uint32_t materialRevision;
    int materialPolicy;
    SbBool databaseMetadataValid;
    int databaseRegionId;
    int databaseAirCode;
    int databaseMaterialId;
    int databaseLos;
    SbBool databaseMaterialColorValid;
    SbColor databaseMaterialColor;
    SbString databaseMaterialShader;
    SbBool colorOverride;
    SbColor color;
    SbColor selectedColor;
    SbColor highlightedColor;
    SbColor ghostedColor;
    SbBool drawMatrixValid;
    SbMatrix drawMatrix;
    SbBool drawCenterValid;
    SbVec3f drawCenter;
    SbBool drawSizeValid;
    float drawSize;
    SbBool sourceBoundsValid;
    SbBool sourceBoundsExact;
    SbBox3f sourceBounds;
    /* At least one retained compact payload was produced by a primitive's
     * view-dependent CSG LoD callback.  Such geometry must not seed a cache
     * for a different view policy; immutable mesh/PoP payloads remain safe. */
    SbBool hasViewDependentCsgGeometry;
    SbBool stale;
    uint32_t staleReason;
    int realizedShapeCount;
    int realizedMeshCount;
    int realizedMaterialObjectCount;
};

struct BOBOL_EXPORT BObolCompactInstanceHandle {
    BObolCompactInstanceHandle(void);

    uint64_t sourceNodeId;
    uint64_t instanceWord0;
    uint64_t instanceWord1;

    SbBool isValid(void) const;
};

struct BOBOL_EXPORT BObolCompactInstanceSummary {
    BObolCompactInstanceSummary(void);

    SbBool valid;
    BObolCompactInstanceHandle handle;
    SbString path;
    SbString sourceName;
    SbString sourceInstanceKey;
    /* Current resident presentation tier ("aabb", "obb", "surface",
     * "line", ...).  This describes the compact source fallback, not a
     * view-local LoD payload that may temporarily supersede it. */
    SbString geometryKind;
    /* Canonical source asset for progressive mesh residency.  Empty values
     * mean the occurrence's sourceName/path are already canonical. */
    SbString meshAssetPath;
    SbString meshAssetName;
    SbBox3f meshAssetBounds;
    uint64_t sourceContentHash;
    uint64_t sourceFaceCount;
    uint64_t sourcePointCount;
    SbMatrix localToSource;
    SbBox3f localBounds;
    uint64_t geometryIdentity;
    uint64_t geometryRevision;
    uint64_t appearanceRevision;
    uint64_t placementRevision;
    uint64_t visibilityRevision;
    uint64_t selectionRevision;
    uint32_t occurrenceIndex;
    int booleanOperation;
    int regionId;
    int airCode;
    int materialId;
    int los;
    SbBool materialColorValid;
    SbColor materialColor;
    SbString materialShader;
    /* Effective compact-instance presentation after selected/highlighted
     * state has been applied. */
    SbBool appearanceColorValid;
    SbColor appearanceColor;
    int lineStyle;
    int lineWidth;
    float transparency;
    SbBool wireGeometry;
    SbBool pointGeometry;
    SbBool meshGeometry;
    SbBool lodBacked;
    SbBool sourceMeshRequestValid;
    SbBool visible;
    SbBool selectable;
    SbBool selected;
    SbBool highlighted;
};

/* Provider-only tail of a compact source request.  This is materialized only
 * after retained/shared/structural fast paths have proved that an asynchronous
 * cache task is required. */
struct BOBOL_EXPORT BObolCompactLodProviderSummary {
    BObolCompactLodProviderSummary(void);

    /* One-shot strong handoff from a live source's bounded staging stream.
     * Materializing this provider-only tail means task submission has already
     * passed all cheap resident/shared/admission checks. */
    std::shared_ptr<const BObolStagedSourceMesh> stagedSource;
    SbBool lodAvailable;
    int lodActiveCut;
    uint32_t lodFaceCount;
    uint32_t lodPointCount;
    SbBool lodHasNormals;
};

/* Minimal immutable occurrence snapshot consumed by the view LoD planner.
 * Keep selection/material/edit metadata out of this hot path: a many-leaf
 * camera epoch should copy only what is needed to cull, prioritize, and
 * construct a progressive asset request. */
struct BOBOL_EXPORT BObolCompactLodInstanceSummary {
    BObolCompactLodInstanceSummary(void);

    SbBool valid;
    SbString path;
    SbString sourceName;
    SbString sourceInstanceKey;
    SbString meshAssetPath;
    SbString meshAssetName;
    SbBox3f meshAssetBounds;
    uint64_t sourceContentHash;
    uint64_t sourceFaceCount;
    uint64_t sourcePointCount;
    /* Fixed-size representation facts needed before the planner knows
     * whether provider work is necessary.  Keeping them here avoids copying
     * the legacy, string-heavy source request for every retained leaf. */
    SbBool brepSource;
    double meshAssetTessellationAbsTol;
    double meshAssetTessellationRelTol;
    double meshAssetTessellationNormTol;
    SbMatrix localToSource;
    SbBox3f localBounds;
    /* The currently authored/displayed part may be a normalized structural
     * proxy while localBounds/localToSource describe the canonical mesh
     * asset.  Point-coverage decisions must use the former, since that is the
     * geometry SoCADAssembly actually classifies. */
    SbMatrix presentationLocalToSource;
    SbBox3f presentationLocalBounds;
    std::array<SbVec3f, 8> presentationCorners;
    SbBool presentationCornersValid;
    SbBool meshGeometry;
    SbBool lodBacked;
    SbBool sourceMeshRequestValid;
    SbBool visible;
    SbBool selected;
    SbBool highlighted;
};

/* Allocation-light subset used while ranking every leaf in a view.  Full
 * object/path/cache request strings are materialized only for the bounded
 * entries that the submit action will actually process. */
struct BOBOL_EXPORT BObolCompactLodPlanningSummary {
    BObolCompactLodPlanningSummary(void);

    SbBool valid;
    SbString sourceInstanceKey;
    /* Projection inputs change only with geometry/placement.  View-local
     * planners use these fixed-width revisions to reuse projected evidence
     * across policy-only passes without hashing strings or comparing
     * matrices for every occurrence. */
    uint64_t geometryRevision;
    uint64_t placementRevision;
    uint64_t sourceContentHash;
    uint64_t sourceFaceCount;
    uint64_t sourcePointCount;
    SbBool botSource;
    SbBool brepSource;
    double meshAssetTessellationAbsTol;
    double meshAssetTessellationRelTol;
    double meshAssetTessellationNormTol;
    SbMatrix localToSource;
    SbBox3f localBounds;
    SbMatrix presentationLocalToSource;
    SbBox3f presentationLocalBounds;
    std::array<SbVec3f, 8> presentationCorners;
    SbBool presentationCornersValid;
    SbBool meshGeometry;
    SbBool lodBacked;
    SbBool sourceMeshRequestValid;
    /* The authored compact part contains a producer-certified progressive
     * WireRep.  It is retargeted directly for wire drawing.  If a separate
     * source-mesh request exists, shaded and hybrid modes still require that
     * provider channel. */
    SbBool residentProgressiveGeometry;
    SbBool visible;
    SbBool selected;
    SbBool highlighted;
};

/* On-demand detail for an authored compact progressive part.  Keeping the
 * fixed cut tables out of BObolCompactLodPlanningSummary avoids copying them
 * while scanning ordinary 50k/150k occurrence populations. */
struct BOBOL_EXPORT BObolCompactResidentProgressiveSummary {
    BObolCompactResidentProgressiveSummary(void);

    SbBool valid;
    SbBool wire;
    int minimumCut;
    int residentCut;
    std::array<size_t, Obol::ProgressiveCutLimit> primitiveCounts;
    std::array<float, Obol::ProgressiveCutLimit> normalizedErrors;
};

struct BOBOL_EXPORT BObolAuxiliaryLineSetDisplayState {
    BObolAuxiliaryLineSetDisplayState(void);

    SbBool valid;
    int drawMode;
    SbBool visible;
    SbBool highlighted;
    int lineStyle;
    int lineWidth;
    float transparency;
    SbBool materialColorValid;
    SbColor materialColor;
    uint32_t materialRevision;
};

/*
 * External primary geometry published to a database source is source-local.
 * Callers must not pre-apply the source drawMatrix to these point arrays.
 * Bridges from legacy instance/world-space producers are responsible for
 * normalizing points back to source-local coordinates before publication.
 */
struct BOBOL_EXPORT BObolExternalLineSet {
    BObolExternalLineSet(void);

    const SbVec3f *points;
    const int32_t *commands;
    const double *precisePoints;
    int count;
    const char *sourceType;
    const char *geometryKind;
};

struct BOBOL_EXPORT BObolExternalPointSet {
    BObolExternalPointSet(void);

    const SbVec3f *points;
    const double *precisePoints;
    int count;
    const char *sourceType;
    const char *geometryKind;
};

struct BOBOL_EXPORT BObolExternalTriangleMesh {
    BObolExternalTriangleMesh(void);

    const SbVec3f *points;
    int pointCount;
    const int32_t *indices;
    int indexCount;
    /* Optional triangle-corner normals.  When present normalCount must match
     * indexCount (three normals per triangle). */
    const SbVec3f *normals;
    int normalCount;
    const char *sourceType;
    const char *geometryKind;
    SbBool lodBacked;
};

struct BOBOL_EXPORT BObolExternalAnnotationSegment {
    enum SegmentKind {
	SEGMENT_NONE = 0,
	SEGMENT_LINE = 1,
	SEGMENT_TEXT = 2
    };

    BObolExternalAnnotationSegment(void);

    int kind;
    int lineStart;
    int lineEnd;
    int textRefPoint;
    const char *text;
};

struct BOBOL_EXPORT BObolExternalAnnotation {
    BObolExternalAnnotation(void);

    SbVec3f basePoint;
    const SbVec3f *linePoints;
    const int32_t *lineCommands;
    const double *preciseLinePoints;
    int linePointCount;
    const SbVec3f *annotationPoints;
    const double *preciseAnnotationPoints;
    int annotationPointCount;
    const BObolExternalAnnotationSegment *segments;
    int segmentCount;
    const char *sourceType;
    const char *geometryKind;
};

struct BOBOL_EXPORT BObolDatabaseSourceDisplayPatch {
    BObolDatabaseSourceDisplayPatch(void);

    SbBool visibleValid;
    SbBool visible;
    SbBool selectedValid;
    SbBool selected;
    SbBool highlightedValid;
    SbBool highlighted;
    SbBool lineStyleValid;
    int lineStyle;
    SbBool lineWidthValid;
    int lineWidth;
    SbBool transparencyValid;
    float transparency;
    SbBool colorOverrideValid;
    SbBool colorOverride;
    SbBool colorValid;
    SbColor color;
    SbBool selectedColorValid;
    SbColor selectedColor;
    SbBool highlightedColorValid;
    SbColor highlightedColor;
    SbBool ghostedColorValid;
    SbColor ghostedColor;
};

struct BOBOL_EXPORT BObolRealizedShapeSummary {
    enum ShapeKind {
	SHAPE_UNKNOWN = 0,
	SHAPE_VLIST = 1,
	SHAPE_MESH = 2
    };

    BObolRealizedShapeSummary(void);

    SbBool valid;
    int shapeKind;
    SbString path;
    SbString sourceName;
    SbString sourceType;
    uint32_t sourceId;
    int ownerSourceIndex;
    SbString ownerSourcePath;
    SbString ownerSourceInstanceKey;
    int ownerDrawMode;
    uint32_t ownerSourceRevision;
    uint32_t ownerInputsRevision;
    uint32_t ownerViewRevision;
    uint32_t ownerRealizedRevision;
    uint32_t ownerRealizedSourceRevision;
    uint32_t ownerRealizedInputsRevision;
    uint32_t ownerRealizedViewRevision;
    int ownerRealizationStatus;
    SbString ownerRealizationDiagnostic;
    SbString ownerRealizationIdentity;
    SbBool ownerSourceStale;
    uint32_t ownerStaleReason;
    SbString displayName;
    SbString geometryName;
    SbString cacheIdentity;
    SbString sourceIdentity;
    SbBool databaseIntent;
    SbBool overlayIntent;
    SbBool hudIntent;
    SbBool localSource;
    SbBool sharedSource;
    SbBool nonDatabaseSource;
    int drawMode;
    SbString recordRole;
    SbString geometryKind;
    int regionId;
    int airCode;
    int materialId;
    int los;
    SbBool materialColorValid;
    SbColor materialColor;
    SbString materialShader;
    uint32_t materialRevision;
    SbBool visible;
    SbBool selectable;
    SbBool selected;
    SbBool highlighted;
    SbBool ghosted;
    SbBool hiddenLine;
    SbBool editEmphasis;
    int lineStyle;
    int lineWidth;
    float transparency;
    SbString editIntentId;
    SbString editIntentRole;
    uint32_t lodPolicy;
    SbBool colorOverride;
    SbColor color;
    SbBool lodAvailable;
    int lodActiveCut;
    uint64_t lodFaceCount;
    uint64_t lodPointCount;
    uint64_t lodOriginalPointCount;
    uint64_t lodNormalCount;
    SbBool lodHasSnappedPoints;
    SbBool lodHasNormals;
    SbVec3f lodBoundsMin;
    SbVec3f lodBoundsMax;
    int pointCount;
    int commandCount;
    int segmentCount;
    int pointPrimitiveCount;
    int triangleCount;
    int indexCount;
    SbBool boundsValid;
    SbBox3f bounds;
};

/** Immutable geometry and semantic state for one compact CAD occurrence. */
struct BOBOL_EXPORT BObolCompactOccurrence {
    BObolCompactOccurrence(void);

    std::shared_ptr<const Obol::PartGeometry> geometry;
    BObolRealizedShapeSummary summary;
    /* Maps geometry coordinates into the occurrence's object-local frame. */
    SbMatrix geometryTransform;
    SbMatrix localTransform;
    SbBool viewDependentCsgGeometry;
    SbBool lodBacked;
    SbBool sourceMeshRequestValid;
    BObolSourceMeshRequest sourceMeshRequest;
    uint32_t occurrenceIndex;
    int booleanOperation;
};

/* Immutable aggregate metadata for one complete compact-source discovery.
 * It is mesh-buffer-free: consumers may use it to bound planning work, but it
 * never authorizes allocation or drawing by itself. */
struct BOBOL_EXPORT BObolCompactSourceProfile {
    BObolCompactSourceProfile(void);

    SbBool valid;
    uint64_t occurrenceCount;
    uint64_t uniqueAssetCount;
    uint64_t encodedSourceBytes;
    uint64_t largestAssetBytes;
    uint64_t reusedOccurrenceCount;
};

/** Lightweight, mesh-buffer-free persistence record for one authoritative
 * compact occurrence.  A progressive producer may be drained by the scene
 * owner before detached realization completes, so the producer retains this
 * small immutable subset instead of a second occurrence registry. */
struct BOBOL_EXPORT BObolCompactManifestOccurrence {
    BObolCompactManifestOccurrence(void);

    SbString path;
    SbString sourceName;
    SbMatrix localTransform;
    SbBox3f bounds;
    /* Optional object-coordinate PCA bounds for a retained aggregate proxy.
     * The AABB above remains the authoritative coverage contract. */
    SbBool orientedBoundsValid;
    std::array<SbVec3f, 8> orientedBounds;
    int booleanOperation;
    uint32_t occurrenceIndex;
    SbBool sourceMeshRequestValid;
    SbString sourceType;
    SbString meshAssetPath;
    SbString meshAssetName;
    uint64_t meshAssetContentHash;
    double meshAssetTessellationAbsTol;
    double meshAssetTessellationRelTol;
    double meshAssetTessellationNormTol;
    SbBox3f meshAssetBounds;
    SbMatrix meshAssetTransform;
    uint64_t sourceFaceCount;
    uint64_t sourcePointCount;
    int regionId;
    int airCode;
    int materialId;
    int los;
    SbBool materialColorValid;
    SbColor materialColor;
    SbString materialShader;
};

/** Thread-safe hand-off of completed compact occurrences from a realization
 * worker thread to the progressive pump thread.  The worker pushes copies of
 * finished per-leaf occurrences as they are tessellated/plotted; the pump drains
 * batches and appends them to the live compact root so geometry streams in
 * incrementally instead of appearing all at once when the walk completes.
 *
 * BObolCompactOccurrence::geometry is a shared_ptr<const PartGeometry> (immutable
 * pointee), so pushing a copy across threads is a race-free refcount bump; all
 * other fields are value types. */
struct BOBOL_EXPORT BObolCompactOccurrenceStream {
    BObolCompactOccurrenceStream(void);
    ~BObolCompactOccurrenceStream(void);
    BObolCompactOccurrenceStream(
	const BObolCompactOccurrenceStream &) = delete;
    BObolCompactOccurrenceStream &operator=(
	const BObolCompactOccurrenceStream &) = delete;

    void push(const BObolCompactOccurrence &occurrence);
    void push(BObolCompactOccurrence &&occurrence);
    /* Publish a producer-local batch with one queue lock and one amortized
     * capacity change.  Ordering within the batch is preserved; ordering
     * between concurrent producers is intentionally unspecified, just as it
     * is for repeated push() calls. */
    void pushBatch(std::vector<BObolCompactOccurrence> &&occurrences);
    /* Priority occurrences describe the current whole draw target rather than
     * one leaf.  The newest undrained occurrence supersedes its predecessor;
     * drain it before an already queued leaf backlog so a completed aggregate
     * extent can frame a very large cold draw immediately. */
    void pushPriority(const BObolCompactOccurrence &occurrence);
    /* Retain only the cache-persistent metadata for an authoritative leaf.
     * Repeated records with the same semantic path update the prior entry, as
     * happens when transformed-reuse proof enriches an already visible box.
     * sealManifest succeeds only for one complete semantic population.  A
     * mixed population may contain structural records for non-mesh leaves;
     * warm consumers independently decide whether that seed is terminal.
     * takeManifest is unavailable until sealing and transfers the one-shot
     * journal to its persistence callback without a second large copy. */
    void recordManifestOccurrence(const BObolCompactOccurrence &occurrence);
    bool sealManifest(size_t expectedCount);
    bool takeManifest(
	std::vector<BObolCompactManifestOccurrence> &occurrences);
    /* Retain a bounded LRU window of full cold imports long enough for the
     * first view-selected LoD task to reuse them.  Occurrences hold weak
     * references; the live source may retain this stream across terminal
     * realization adoption until a provider claims the matching import. */
    SbBool retainStagedSource(
	const std::shared_ptr<const BObolStagedSourceMesh> &source);
    /* Transfer one retained import into a provider.  The returned shared
     * pointer is the new owner and the stream's staging-budget charge ends.
     * A weak request which is absent, stale, or already claimed returns an
     * empty pointer. */
    std::shared_ptr<const BObolStagedSourceMesh> claimStagedSource(
	const std::weak_ptr<const BObolStagedSourceMesh> &source);
    size_t stagedSourceByteCount(void);
    size_t drain(std::vector<BObolCompactOccurrence> &out, size_t cap);
    size_t size(void);
    void setExpectedCount(size_t count);
    size_t getExpectedCount(void) const;
    /** Publish an exact finite denominator for the source-representation work
     * which follows structural coverage.  Completion is monotonic for the
     * lifetime of this stream and counts resolved producer items whether they
     * produce LoD-backed or terminal geometry. */
    void setPreparationWorkCount(size_t count);
    void notePreparationWorkCompleted(void);
    void getPreparationWorkCount(size_t &completed, size_t &total) const;
    void setSourceProfile(const BObolCompactSourceProfile &profile);
    SbBool getSourceProfile(BObolCompactSourceProfile &profile) const;
    /* A warm manifest may describe the complete semantic population while
     * some non-mesh occurrences still need terminal geometry regeneration.
     * Keep that census proof distinct from representation readiness: it
     * authorizes selective replay, not declaring the source realized. */
    void setWarmCensusComplete(bool complete);
    bool hasWarmCensusComplete(void) const;
    void setWarmCoverageComplete(bool complete);
    bool hasWarmCoverageComplete(void) const;
    void recordWarmTerminalOccurrence(
	const BObolCompactManifestOccurrence &occurrence);
    bool takeWarmTerminalOccurrences(
	std::vector<BObolCompactManifestOccurrence> &occurrences);
    /* Publish/query the source-local extent of the draw target.  An early
     * value may be a conservative, monotonically growing overview.  Once
     * hasCoverageBoundsComplete() is true the value is immutable and
     * semantically exact. */
    void setCoverageBounds(const SbBox3f &bounds);
    bool getCoverageBounds(SbBox3f &bounds);
    /*
     * True once a semantically exact whole-target bound has been queued in the priority
     * lane.  Progressive autoview must not chase the append-only union of
     * partial leaf coverage: doing so recenters the camera on every streamed
     * batch and makes an otherwise monotonic cold draw appear to flicker.
     */
    void setCoverageBoundsComplete(bool complete);
    bool hasCoverageBoundsComplete(void) const;
    /** True after a complete whole-target coverage occurrence has been
     * consumed by the owner thread.  The extent may be conservative when
     * hasCoverageBoundsComplete() is false, but it already covers every leaf
     * and is therefore safe for the one-shot progressive autoview. */
    bool hasCoverageBoundsDrained(void);
    void requestCancel(void);
    bool isCancelled(void) const;

private:
    struct Impl;
    std::unique_ptr<Impl> d;
};

struct BOBOL_EXPORT BObolRealizedMaterialSummary {
    BObolRealizedMaterialSummary(void);

    SbBool valid;
    SbString sourcePath;
    SbString sourceName;
    SbString sourceType;
    uint32_t sourceId;
    int ownerSourceIndex;
    SbString ownerSourcePath;
    SbString ownerSourceInstanceKey;
    int ownerDrawMode;
    uint32_t ownerSourceRevision;
    uint32_t ownerInputsRevision;
    uint32_t ownerViewRevision;
    uint32_t ownerRealizedRevision;
    uint32_t ownerRealizedSourceRevision;
    uint32_t ownerRealizedInputsRevision;
    uint32_t ownerRealizedViewRevision;
    int ownerRealizationStatus;
    SbString ownerRealizationDiagnostic;
    SbString ownerRealizationIdentity;
    SbBool ownerSourceStale;
    uint32_t ownerStaleReason;
    SbString materialName;
    SbString parentName;
    SbString materialSource;
    int propertyCount;
};

struct BOBOL_EXPORT BObolSceneTreeSummary {
    enum NodeKind {
	NODE_UNKNOWN = 0,
	NODE_GROUP = 1,
	NODE_DATABASE_SOURCE = 2,
	NODE_VLIST_SHAPE = 3,
	NODE_MESH_SHAPE = 4,
	NODE_MATERIAL_OBJECT = 5,
	NODE_OTHER = 6
    };

    BObolSceneTreeSummary(void);

    SbBool valid;
    int nodeKind;
    SbBool isGroup;
    SbBool isShape;
    SbBool isDatabaseSource;
    SbBool isMaterialObject;
    SbBool hasParent;
    int drawTreeDepth;
    int childCount;
    int ownerSourceIndex;
    SbString ownerSourcePath;
    SbString ownerSourceInstanceKey;
    SbString path;
    SbString sourceName;
    SbString sourceType;
    uint32_t sourceId;
    SbString displayName;
    SbString geometryName;
};

struct BOBOL_EXPORT BObolSceneDisplaySummary {
    BObolSceneDisplaySummary(void);

    SbBool valid;
    int nodeKind;
    SbBool isDatabaseSource;
    SbBool hasDrawIntent;
    SbString intentPath;
    int intentDrawMode;
    SbBool visible;
    SbBool selected;
    SbBool highlighted;
    int lineStyle;
    int lineWidth;
    double transparency;
    int drawMode;
    SbBool materialValid;
    SbColor materialColor;
    uint32_t materialRevision;
    SbBool drawMatrixValid;
    SbMatrix drawMatrix;
    SbBool drawCenterValid;
    SbVec3f drawCenter;
    SbBool drawSizeValid;
    float drawSize;
    int ownerSourceIndex;
    SbString ownerSourcePath;
    SbString ownerSourceInstanceKey;
    SbString path;
};

struct BOBOL_EXPORT BObolSceneMaterialSummary {
    BObolSceneMaterialSummary(void);

    SbBool valid;
    int nodeKind;
    SbBool materialValid;
    uint32_t materialRevision;
    SbColor materialColor;
    int ownerSourceIndex;
    SbString ownerSourcePath;
    SbString ownerSourceInstanceKey;
    SbString path;
};

struct BOBOL_EXPORT BObolSceneBoundsSummary {
    BObolSceneBoundsSummary(void);

    SbBool valid;
    int nodeKind;
    SbBool boundsValid;
    SbBox3f bounds;
    int ownerSourceIndex;
    SbString ownerSourcePath;
    SbString ownerSourceInstanceKey;
    SbString path;
};

class BOBOL_EXPORT SoBRLDatabaseSource : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLDatabaseSource);

public:
    enum DrawMode {
	WIREFRAME = 0,
	SHADED = 1
    };

    enum RepresentationMode {
	REPRESENTATION_DEFAULT = -1,
	REPRESENTATION_WIRE = 0,
	REPRESENTATION_SHADED_BOTS = 1,
	REPRESENTATION_SHADED = 2,
	REPRESENTATION_EVAL_WIRE = 3,
	REPRESENTATION_HIDDEN_LINE = 4,
	REPRESENTATION_EVAL_POINTS = 5
    };

    enum RealizationStatus {
	UNREALIZED = 0,
	REALIZED = 1,
	FAILED = 2
    };

    enum StaleReason {
	STALE_NONE = 0,
	STALE_SOURCE = 1,
	STALE_INPUTS = 2,
	STALE_VIEW = 4,
	STALE_DATABASE = 8,
	STALE_DRAW = 16,
	STALE_TESSELLATION = 32
    };

    enum RealizationRoleFlag {
	REALIZATION_ROLE_NONE = 0,
	REALIZATION_ROLE_CSG = 1,
	REALIZATION_ROLE_MESH = 2,
	REALIZATION_ROLE_EXTERNAL = 4
    };

    enum MaterialPolicy {
	MATERIAL_INHERIT = 0,
	MATERIAL_DATABASE = 1
    };

    enum BooleanOperation {
	BOOLEAN_UNION = 0,
	BOOLEAN_SUBTRACT = 1,
	BOOLEAN_INTERSECT = 2
    };

    SoSFString instanceKey;
    SoSFString path;
    SoSFString parentInstanceKey;
    SoSFUInt32 occurrenceIndex;
    SoSFInt32 booleanOperation;
    SoSFString displayName;
    SoSFString representationKey;
    SoSFInt32 representationMode;
    SoSFBool auxiliarySource;
    SoSFEnum drawMode;
    SoSFBool visible;
    SoSFBool selected;
    SoSFBool highlighted;
    SoSFInt32 lineStyle;
    SoSFInt32 lineWidth;
    SoSFFloat transparency;
    SoSFBool materialColorValid;
    SoSFColor materialColor;
    SoSFUInt32 materialRevision;
    SoSFEnum materialPolicy;
    SoSFBool databaseMetadataValid;
    SoSFInt32 databaseRegionId;
    SoSFInt32 databaseAirCode;
    SoSFInt32 databaseMaterialId;
    SoSFInt32 databaseLos;
    SoSFBool databaseMaterialColorValid;
    SoSFColor databaseMaterialColor;
    SoSFString databaseMaterialShader;
    SoSFBool colorOverride;
    SoSFColor color;
    SoSFColor selectedColor;
    SoSFColor highlightedColor;
    SoSFColor ghostedColor;
    /* Local-to-scene placement metadata for this source instance.  Primary
     * geometry remains source-local; render traversal and summaries apply
     * this transform explicitly instead of baking it into point/index arrays.
     */
    SoSFBool drawMatrixValid;
    SoSFMatrix drawMatrix;
    SoSFBool drawCenterValid;
    SoSFVec3f drawCenter;
    SoSFBool drawSizeValid;
    SoSFFloat drawSize;
    SoSFBool sourceBoundsValid;
    /* TRUE only when sourceBoundsMin/Max describe the complete immutable
     * target extent rather than the union of a partially streamed registry. */
    SoSFBool sourceBoundsExact;
    SoSFVec3f sourceBoundsMin;
    SoSFVec3f sourceBoundsMax;
    SoSFFloat tessellationAbsTol;
    SoSFFloat tessellationRelTol;
    SoSFFloat tessellationNormTol;
    SoSFUInt32 lodBotThreshold;
    SoSFUInt32 sourceRevision;
    SoSFUInt32 inputsRevision;
    SoSFUInt32 viewRevision;
    SoSFUInt32 realizedRevision;
    SoSFUInt32 realizedSourceRevision;
    SoSFUInt32 realizedInputsRevision;
    SoSFUInt32 realizedViewRevision;
    SoSFEnum realizationStatus;
    SoSFString realizationDiagnostic;
    SoSFString realizationIdentity;
    SoSFInt32 realizationRoleFlags;
    SoSFBool realizationViewDependent;
    SoSFBool realizationCsgLodEnabled;
    SoSFBool realizationMeshLodEnabled;
    SoSFFloat realizationViewScale;
    SoSFFloat realizationLodScale;
    SoSFInt32 realizationViewWidth;
    SoSFInt32 realizationViewHeight;
    SoSFUInt32 realizationBotThreshold;
    SoSFFloat realizationCurveScale;
    SoSFFloat realizationPointScale;
    SoSFBool stale;
    SoSFUInt32 staleReason;

    SoBRLDatabaseSource(void);
    static void initClass(void);

    /**
     * Set the BRL-CAD database used for realization.
     *
     * The caller owns the db_i pointer.  SoBRLDatabaseSource stores a borrowed
     * pointer and does not close or otherwise manage the database lifetime.
     * The database must remain valid while this source can be realized.
     */
    void setDatabase(struct db_i *dbip);
    struct db_i *getDatabase(void) const;
    void setMeshLod(struct BObolMeshLod *lod);
    struct BObolMeshLod *getMeshLod(void) const;
    void clearMeshLod(void);
    int setMeshLodBounds(const SbVec3f &bmin, const SbVec3f &bmax);
    SbBool getMeshLodBounds(SbVec3f &bmin, SbVec3f &bmax) const;
    void clearMeshLodBounds(void);
    void configureDatabaseSource(const char *sourcePath,
	struct db_i *database,
	int mode,
	uint32_t revision);
    void configureDatabaseSourceInstance(const char *sourceInstanceKey,
	const char *sourcePath,
	struct db_i *database,
	int mode,
	uint32_t revision);
    void configureDatabaseSourceInstanceRepresentation(
	const char *sourceInstanceKey,
	const char *sourcePath,
	const char *sourceRepresentationKey,
	int sourceRepresentationMode,
	struct db_i *database,
	int mode,
	uint32_t revision);
    /* Capture the small, immutable field configuration needed by detached
     * realization without copying database geometry.  This is safe to call on
     * the scene owner thread and returns a ref-counted unattached source with
     * field sensors disabled.  After launch it is worker-exclusive until the
     * realization job publishes a terminal state. */
    SoBRLDatabaseSource *createDetachedRealizationTemplate(void) const;
    /* Populate a detached template on a realization worker.  File-backed
     * inputs use an independent read handle and therefore require the caller
     * to reject the complete stream when its captured source revision is no
     * longer current; non-reopenable inputs use an immutable closure
     * snapshot.  The caller owns the returned database and any nonempty
     * temporary snapshot path. */
    SbBool initializeDetachedRealizationDatabase(
	struct db_i *sourceDatabase,
	struct db_i **databaseOut,
	SbString *snapshotPathOut = NULL);
    /* Synchronous immutable-snapshot form for consumers without a
     * revision-gated stream.  The caller must unref the returned source,
     * db_close the returned database, and delete snapshotPathOut when it is
     * nonempty. */
    SoBRLDatabaseSource *createDetachedRealizationSource(
	struct db_i **databaseOut, SbString *snapshotPathOut = NULL) const;
    /* Transfer a completed detached compact registry into this live source.
     * Call only on the live source's owner thread. */
    int adoptDetachedCompactRealization(SoBRLDatabaseSource *detached,
	SbBool authoritativeStreamDrained = FALSE,
	const std::shared_ptr<BObolCompactOccurrenceStream> &stagedSourceStream =
	    std::shared_ptr<BObolCompactOccurrenceStream>());
    int retargetDatabaseSource(const char *sourcePath,
	uint32_t revision);
    int retargetDatabaseSourceInstance(const char *sourceInstanceKey,
	const char *sourcePath,
	uint32_t revision);
    int setRepresentationState(const char *sourceRepresentationKey,
	int sourceRepresentationMode);
    int setHierarchyState(const char *sourceParentInstanceKey,
	uint32_t sourceOccurrenceIndex,
	int sourceBooleanOperation);

    void markStale(void);
    void markStale(uint32_t reason);
    int setDrawModeState(int drawMode);
    /* Resolve the authoritative representation and legacy draw-mode fields
     * to the BOBOL_LOD_DRAW_* value used by realization, worker preparation,
     * and presentation.  All three stages must use this same contract or the
     * scene owner is forced to derive a replacement geometry after handoff. */
    int getEffectiveLodDrawMode(void) const;
    /* Return whether database realization must preserve triangle-mesh
     * identity.  REALIZATION_ROLE_* fields are orchestration hints and may
     * transiently describe an external publication; view-managed wire BoTs
     * still require the mesh path because their edges come from a PoP
     * triangle prefix.  Deferred workers and synchronous realization must
     * consult this one geometry contract. */
    SbBool usesMeshRealization(void) const;
    int setDisplayNameState(const char *name);
    int setMaterialPolicyState(int materialPolicy);
    int setRealizationState(int realizationStatus,
	uint32_t realizedSourceRevision,
	uint32_t realizedInputsRevision,
	uint32_t staleReason,
	const char *diagnostic = NULL);
    int setRealizationRoleFlags(int roleFlags);
    int setRealizationViewPolicy(SbBool viewDependent,
	SbBool csgLodEnabled,
	SbBool meshLodEnabled,
	float viewScale,
	float lodScale,
	int viewWidth,
	int viewHeight,
	uint32_t botThreshold,
	float curveScale,
	float pointScale);
    int setDisplayState(SbBool sourceRevisionValid,
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
    int setDatabaseMetadataState(SbBool metadataValid,
	int regionId,
	int airCode,
	int materialId,
	int los,
	SbBool metadataMaterialColorValid,
	const SbColor &metadataMaterialColor,
	const SbString &metadataMaterialShader);
    int refreshMaterialColorFromDatabase(uint32_t materialRevision,
	struct db_i *overrideDbip = NULL);
    int applyDisplayPatch(const BObolDatabaseSourceDisplayPatch &patch);
    /* Update local-to-scene placement metadata.  This does not mutate primary
     * geometry coordinates; database realization and external publication
     * must keep drawable arrays in source-local space.
     */
    int setPlacementState(SbBool drawMatrixValid,
	const SbMatrix &drawMatrix,
	SbBool drawCenterValid,
	const SbVec3f &drawCenter,
	SbBool drawSizeValid,
	float drawSize);
    int setSourceBoundsState(SbBool boundsValid,
	const SbVec3f &boundsMin,
	const SbVec3f &boundsMax,
	SbBool boundsExact = FALSE);
    int setSourceBoundsExactState(SbBool boundsExact);
    void clearSourceBounds(void);
    SbBool getSourceBounds(SbBox3f &bounds) const;
    SbBool getEffectiveSourceBounds(SbBox3f &bounds) const;
    SbBool hasExactSourceBounds(void) const;
    SbBool needsRealization(void) const;
    SbBool realizePrototypeWireframe(void);
    SbBool realizeDatabaseWireframe(void);
    SbBool realizeDatabaseMesh(void);
    /* Streaming variants: publish each completed per-leaf occurrence to the
     * given stream as the walk produces it, so a consumer can adopt geometry
     * incrementally.  Passing NULL is identical to the zero-argument form. */
    SbBool realizeDatabaseWireframe(BObolCompactOccurrenceStream *stream);
    SbBool realizeDatabaseMesh(BObolCompactOccurrenceStream *stream);
    int clearRealizedGeometry(SbBool preserveAuxiliary = TRUE);
    int clearExternalPrimaryGeometry(void);
    /* Publish source-local primary geometry for externally supplied source
     * realizations.  These APIs preserve source placement separately via
     * setPlacementState/drawMatrix.
     */
    int publishExternalLineSet(const BObolExternalLineSet &lineSet);
    int publishExternalPointSet(const BObolExternalPointSet &pointSet);
    int publishExternalTriangleMesh(
	const BObolExternalTriangleMesh &triangleMesh);
    int publishExternalAnnotation(
	const BObolExternalAnnotation &annotation);
    int publishPrimitiveWireframe(
	struct rt_db_internal *intern,
	const struct bg_tess_tol *ttol = NULL,
	const struct bn_tol *tol = NULL);
    SoBRLVListShape *getRealizedShape(void) const;
    SoBRLVListShape *getRealizedShape(int index) const;
    int getRealizedShapeCount(void) const;
    SbBool hasRealizedWireGeometry(void) const;
    SoBRLVListShape *findAuxiliaryVListShape(const char *name) const;
    SoBRLDatabaseSource *findAuxiliarySource(const char *sourcePath) const;
    int setAuxiliaryLineSet(const char *name,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	const BObolAuxiliaryLineSetDisplayState *displayState = NULL);
    int setAuxiliarySourceLineSet(const char *sourcePath,
	const char *auxDisplayName,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	const BObolAuxiliaryLineSetDisplayState *displayState = NULL);
    int clearAuxiliaryShapes(void);
    SoBRLMeshShape *getRealizedMesh(void) const;
    SoBRLMeshShape *getRealizedMesh(int index) const;
    int getRealizedMeshCount(void) const;
    SbBool hasRealizedMeshGeometry(void) const;
    /** True when the compact presentation contains any geometry whose
     * displayed fidelity is selected by the view.  This includes both
     * asynchronous source-mesh requests and immutable resident progressive
     * parts; it is the controller scheduling predicate. */
    SbBool hasDisplayLodTargets(void) const;
    /** Number of current view-managed compact occurrences.  The count is
     * maintained with the compact index so querying a 150k-object source is
     * O(1). */
    size_t getDisplayLodTargetCount(void) const;
    /** True when at least one compact part owns immutable producer-certified
     * cuts which the view retargets without a provider request. */
    SbBool hasDisplayResidentProgressiveGeometry(void) const;
    size_t getDisplayResidentProgressiveGeometryCount(void) const;
    /** True when compact occurrences can refine their current display
     * geometry from source-backed mesh requests validated for the source's
     * current source/inputs epoch.  Streamed contracts may be usable before
     * the source-wide detached realization reaches REALIZED. */
    SbBool hasDisplayMeshLodRequests(void) const;
    /** Number of compact occurrences carrying a current source-mesh LoD
     * request contract.  Returns zero when that contract is stale. */
    size_t getDisplayMeshLodRequestCount(void) const;
    /** Monotonic source-local display revision used to invalidate view demand
     * when compact occurrences are added or replaced. */
    uint64_t getDisplayMeshLodRevision(void) const;
    /** Return the compact entries whose view-LoD demand changed after
     * @p revision.  FALSE means the bounded delta history cannot prove a
     * complete answer and the caller must rescan the source.  When non-NULL,
     * @p coverageInvalidated reports whether any returned mutation changed
     * the occurrence population and therefore invalidated a scene-wide
     * coverage proof.  An exact visibility or in-place data mutation leaves
     * it FALSE.  The journal is source-owned and non-consuming, so
     * independent views may advance at different rates. */
    SbBool getDisplayMeshLodChangedEntries(uint64_t revision,
	std::vector<size_t> &entryIndices,
	SbBool *coverageInvalidated = NULL) const;
    /** Monotonic revision for exact occurrence-visibility changes.  Visibility
     * is a view-planning input, not immutable mesh inventory, so it has a
     * separate source-local journal. */
    uint64_t getDisplayMeshLodVisibilityRevision(void) const;
    /** Return occurrence indices whose effective presentation visibility
     * changed after @p revision.  FALSE requires one authoritative source
     * rescan because the bounded non-consuming history is no longer complete. */
    SbBool getDisplayMeshLodVisibilityChangedEntries(uint64_t revision,
	std::vector<size_t> &entryIndices) const;
    SoBRLMaterialObject *getRealizedMaterialObject(void) const;
    SoBRLMaterialObject *getRealizedMaterialObject(int index) const;
    int getRealizedMaterialObjectCount(void) const;
    int setCompactOccurrence(const BObolCompactOccurrence &occurrence);
    int setCompactOccurrenceRegistry(
	const std::vector<BObolCompactOccurrence> &occurrences);
    /* Incrementally merge occurrences into the live compact registry without
     * replacing the whole index.  Each incoming occurrence upgrades the leaf's
     * existing occurrence in place when it carries higher-tier drawing data
     * (AABB < OBB < realized geometry, ranked by geometryKind) and is appended
     * only when that leaf has no occurrence yet.  Used to stream per-leaf
     * geometry onto a standing coarse (box) root during progressive
     * realization so boxes are replaced, not left lingering. */
    int mergeCompactOccurrences(
	const std::vector<BObolCompactOccurrence> &occurrences,
	SbBool authoritativeGeometry = FALSE);
    /** Reserve storage and certify the exact leaf count for the active
     * streaming epoch.  The count excludes any temporary whole-target
     * overview; reaching it permits that overview to retire after the leaf
     * frontier has been presented. */
    void reserveCompactOccurrenceCapacity(size_t expectedCount);
    /** Expected leaf occurrence population for the active streaming epoch.
     * This is a progress denominator, not a residency requirement: a
     * view-dependent terminal state may intentionally retain mesh payloads
     * for only a visible subset. */
    size_t getCompactExpectedInstanceCount(void) const;
    SbBool getCompactSourceProfile(BObolCompactSourceProfile &profile) const;
    void setCompactSourceProfile(const BObolCompactSourceProfile &profile);
    /** Resolve a stable compact occurrence identity to its source-local entry
     * index in expected O(1) time.  Used by retained view work frontiers; it
     * does not expose or pin the underlying registry storage. */
    SbBool getCompactInstanceIndex(
	const char *occurrenceKey, size_t &entryIndex) const;
    /**
     * Share the resident immutable compact payloads of @p source with this
     * source, translating occurrence placement into this source's coordinate
     * frame.  Existing equal-or-richer occurrences remain authoritative.
     *
     * This is used when a broader draw intent subsumes a narrower one: the
     * broader source can display every useful resident result immediately
     * while its own progressive realization requests only missing data.
     */
    int adoptCompactOccurrencesFrom(const SoBRLDatabaseSource *source);
    /**
     * Apply an explicit database rename to retained compact occurrence
     * semantics before realization.  Geometry and occurrence handles remain
     * unchanged; the next authoritative index can therefore transfer runtime
     * state by exact path rather than ambiguous leaf-name/order heuristics.
     */
    int retargetCompactOccurrencePaths(const char *oldPath,
	const char *newPath);
    SbBool hasCompactInstanceIndex(void) const;
    SbBool isCompactOccurrenceRegistry(void) const;
    int getCompactInstanceCount(void) const;
    /* Number of occurrences this source will submit in its selected
     * presentation set.  This is distinct from semantic GED membership. */
    int getCompactSelectedInstanceCount(void) const;
    /* Unique immutable payloads shared by compact occurrences. */
    int getCompactPartCount(void) const;
    SbBool getCompactInstanceHandle(int index,
	BObolCompactInstanceHandle &handle) const;
    SbBool getCompactOccurrence(int index,
	BObolCompactOccurrence &occurrence) const;
    /**
     * Append selectable compact occurrences whose conservative projected
     * bounds overlap an NDC rectangle.  This is the scalable object/window
     * selection path: it uses the same retained occurrence bounds as culling
     * and never scans triangles or performs per-sample raytraces.
     *
     * Returns -1 when this source has no compact occurrence registry;
     * otherwise returns the number of appended records.  Rectangle
     * coordinates use OpenGL NDC (-1..1, lower-left origin).
     */
    int queryCompactRectangle(const SbMatrix &parentToWorld,
	const SbMatrix &viewProjection,
	float minimumX, float minimumY,
	float maximumX, float maximumY,
	std::vector<BObolViewPickRecord> &records) const;
    /** Append this source's conservative whole-target bound when it overlaps
     * an NDC rectangle.  Returns -1 when an occurrence-level registry is
     * available so callers do not publish both root and leaf identities. */
    int querySourceRectangle(const SbMatrix &parentToWorld,
	const SbMatrix &viewProjection,
	float minimumX, float minimumY,
	float maximumX, float maximumY,
	std::vector<BObolViewPickRecord> &records) const;
    /* Flatten visible compact wire instances for non-raster consumers such as
     * vector export without constructing per-occurrence Coin shape nodes. */
    SbBool copyCompactWireGeometry(std::vector<SbVec3f> &points,
	std::vector<int32_t> &commands) const;
    /* Copy one compact occurrence into source coordinates for a transient
     * selection/edit presentation.  The compact index remains authoritative
     * and retains no Coin shape for this operation. */
    SbBool copyCompactInstanceEditGeometry(
	const BObolCompactInstanceHandle &handle,
	std::vector<SbVec3f> &points,
	std::vector<int32_t> &commands,
	BObolCompactInstanceSummary &summary) const;
    SbBool getCompactInstanceSummary(
	const BObolCompactInstanceHandle &handle,
	BObolCompactInstanceSummary &summary) const;
    /* Direct, allocation-light view-planning lookup.  Unlike the general
     * handle API this is intentionally scoped to one immutable compact-index
     * epoch and avoids a second instance-id hash lookup. */
    SbBool getCompactLodInstanceSummary(int index,
	BObolCompactLodInstanceSummary &summary) const;
    SbBool getCompactLodProviderSummary(int index,
	BObolCompactLodProviderSummary &summary) const;
    SbBool getCompactLodPlanningSummary(int index,
	BObolCompactLodPlanningSummary &summary) const;
    SbBool getCompactResidentProgressiveSummary(int index,
	BObolCompactResidentProgressiveSummary &summary) const;
    /* Constant-time occurrence lookup for scene-wide retained-LoD
     * admission.  This avoids rescanning every structural leaf when only
     * the already-displayed mesh bindings need to be coarsened. */
    SbBool getCompactLodPlanningSummaryForKey(const char *occurrenceKey,
	BObolCompactLodPlanningSummary &summary) const;
    SbBool isCompactInstanceHandleValid(
	const BObolCompactInstanceHandle &handle) const;
    /* Constant-time occurrence lookup used by per-view LoD result routing. */
    SbBool hasCompactInstanceKey(const char *occurrenceKey) const;
    /* Stable in-process owner token copied into compact LoD requests. */
    uint64_t getCompactSourceRoutingId(void) const;
    /* Lifetime of the current dense compact-entry index.  Appends preserve
     * this epoch; replacement/reinstallation advances it.  Together with the
     * routing id it lets a view-local cache reject recycled entry indices. */
    uint64_t getCompactPopulationEpoch(void) const;
    /**
     * Resolve the first compact occurrence matching a semantic path.
     *
     * This is the bounded, allocation-free lookup used when an interactive
     * operation must promote one occurrence from a compact source.  It does
     * not create a Coin node or change compact presentation state.
     * visibleOnly skips occurrences suppressed by the presentation frontier.
     */
    SbBool getCompactInstanceForPath(const char *path,
	SbBool includeDescendants,
	SbBool visibleOnly,
	BObolCompactInstanceHandle &handle,
	BObolCompactInstanceSummary &summary) const;
    int getCompactInstanceCountForPath(const char *path,
	SbBool includeDescendants = TRUE) const;
    SbBool getCompactInstanceBoundsForPath(const char *path,
	SbBool includeDescendants, SbBox3f &bounds) const;
    int setCompactInstanceDisplayStateForPath(const char *path,
	SbBool includeDescendants,
	int visibleValid, SbBool visible,
	int selectedValid, SbBool selected,
	int highlightedValid, SbBool highlighted);
    int setCompactInstanceDisplayStateForPathMatch(const char *path,
	BObolCompactPathMatch match,
	int visibleValid, SbBool visible,
	int selectedValid, SbBool selected,
	int highlightedValid, SbBool highlighted);
    int setCompactInstanceTransparencyForPathMatch(const char *path,
	BObolCompactPathMatch match, float nextTransparency);
    /* Retained semantic presentation overrides.  Unlike the immediate
     * display-state helpers above, these rules are also applied to compact
     * occurrences that stream into the source later. */
    int setCompactInstanceVisibilityOverrideForPathMatch(const char *path,
	BObolCompactPathMatch match, SbBool visible);
    int setCompactInstanceHighlightOverrideForPathMatch(const char *path,
	BObolCompactPathMatch match, SbBool highlighted);
    int setCompactInstanceTransparencyOverrideForPathMatch(const char *path,
	BObolCompactPathMatch match, float nextTransparency);
    int clearCompactInstanceHighlightOverrides(void);
    /**
     * Restrict this source's visible compact occurrences to the union of the
     * supplied semantic path subtrees.  Geometry, instances, PoP residency,
     * and LoD request ownership remain on this source; only its presentation
     * mask changes.  The rule is retained and applied to occurrences streamed
     * in later.  An empty frontier hides every compact occurrence.
     */
    int setCompactInstanceVisibilityFrontier(
	const std::vector<SbString> &paths);
    /**
     * Apply a compact, retained set of subtree visibility overrides.
     * Occurrences are authored-visible by default; each path changes the
     * state of its complete subtree, and a deeper path takes precedence over
     * an ancestor.  This represents "draw root except these paths" without
     * expanding one erased path into every visible sibling.  The rule set is
     * also applied to compact occurrences streamed in later.
     */
    int setCompactInstanceVisibilityOverrides(
	const std::vector<SbString> &paths,
	const std::vector<SbBool> &states);
    int clearCompactInstanceVisibilityFrontier(void);
    SbBool hasCompactInstanceVisibilityFrontier(void) const;
    /**
     * Replace the semantic selection mask for compact occurrences.  The mask
     * is retained and applied to occurrences streamed in later.  Paths select
     * their complete semantic subtree; an empty list clears compact selection.
     */
    int syncCompactInstanceSelectedPaths(const std::vector<SbString> &paths);
    /**
     * Apply an incremental semantic selection change.  Removed paths are
     * cleared before added paths are asserted, allowing one commit to replace
     * an ancestor selection with a descendant (or the reverse) without a
     * full compact-index pass.  The resulting frontier is retained for
     * occurrences streamed in later.
     */
    int applyCompactInstanceSelectionDelta(
	const std::vector<SbString> &addedPaths,
	const std::vector<SbString> &removedPaths);
    int setCompactInstanceSelectableForPath(const char *path,
	SbBool includeDescendants, SbBool selectable);
    int setCompactInstanceRegionIdForPath(const char *path,
	SbBool includeDescendants, int regionId);
    int setCompactInstanceRegionMetadataForPath(const char *path,
	SbBool includeDescendants, int regionId, int airCode,
	int materialId, int los);
    int setCompactInstanceMetadataForPath(const char *path,
	SbBool includeDescendants, int regionId, int airCode,
	int materialId, int los, SbBool materialColorValid,
	const SbColor &materialColor, const SbString &materialShader);
    int setCompactSubtractLineStyle(int lineStyle);
    int refreshCompactObjectGeometry(const char *objectPath,
	uint32_t sourceRevision = 0);
    int prepareCompiledAssembly(void);
    SoBRLCadAssembly *compactViewLodAssembly(
	const std::vector<const BObolViewLodState::CadPayload *> &payloads,
	const BObolViewLodState *viewState)
	const;
    SoBRLCadAssembly *currentCompactViewLodAssembly(
	const BObolViewLodState *viewState) const;
    /* Diagnostic invariant: an occurrence with an active view-local LoD
     * payload must not still present its structural AABB/OBB base part. */
    int getCompactViewLodSupersededFallbackCount(
	const BObolViewLodState *viewState,
	std::vector<SbString> *paths = NULL) const;
    /* Number of structural AABB/OBB leaf parts that are still the active
     * presentation, whether or not a richer payload has already arrived.
     * Unlike the superseded-only invariant above, this reports legitimate
     * cold coverage as well as terminal budget starvation. */
    int getCompactViewLodActiveFallbackCount(
	const BObolViewLodState *viewState,
	std::vector<SbString> *paths = NULL) const;
    SbBool hasCompiledAssembly(void) const;
    int getCompiledAssemblyPartCount(void) const;
    int getCompiledAssemblyInstanceCount(void) const;
    SbBool getSummary(BObolDatabaseSourceSummary &summary) const;
    int getRealizedShapeSummaryCount(void) const;
    SbBool getRealizedShapeSummary(int index,
	BObolRealizedShapeSummary &summary) const;
    int getRealizedMaterialSummaryCount(void) const;
    SbBool getRealizedMaterialSummary(int index,
	BObolRealizedMaterialSummary &summary) const;
    SbBool getRealizedMaterialProperty(int materialIndex, int propertyIndex,
	SbString &groupOut, SbString &nameOut, SbString &valueOut) const;
    int getRealizedTreeSummaryCount(void) const;
    SbBool getRealizedTreeSummary(int index,
	BObolSceneTreeSummary &summary) const;
    int getRealizedDisplaySummaryCount(void) const;
    SbBool getRealizedDisplaySummary(int index,
	BObolSceneDisplaySummary &summary) const;
    int getRealizedSceneMaterialSummaryCount(void) const;
    SbBool getRealizedSceneMaterialSummary(int index,
	BObolSceneMaterialSummary &summary) const;
    int getRealizedBoundsSummaryCount(void) const;
    SbBool getRealizedBoundsSummary(int index,
	BObolSceneBoundsSummary &summary) const;

protected:
    virtual ~SoBRLDatabaseSource(void);
    void GLRender(SoGLRenderAction *action) override;
    void GLRenderBelowPath(SoGLRenderAction *action) override;
    void callback(SoCallbackAction *action) override;
    void getBoundingBox(SoGetBoundingBoxAction *action) override;
    void rayPick(SoRayPickAction *action) override;

private:
    friend class SoBRLCadRenderBatch;
    friend class SoBRLExportAction;
    friend class SoBRLMeasureAction;
    friend class SoBRLSnapAction;
    friend SbBool bobol_database_source_realize_wireframe_with_cache(
	    SoBRLDatabaseSource *source,
	    BObolDatabaseSourceRealizationCache *cache);
    friend SbBool bobol_database_source_realize_mesh_with_cache(
	    SoBRLDatabaseSource *source,
	    BObolDatabaseSourceRealizationCache *cache);
    friend int bobol_database_source_realize_wireframe_compact_with_cache(
	    SoBRLDatabaseSource *source,
	    BObolDatabaseSourceRealizationCache *cache,
	    BObolCompactOccurrenceStream *stream);
    friend int bobol_database_source_realize_mesh_compact_with_cache(
	    SoBRLDatabaseSource *source,
	    BObolDatabaseSourceRealizationCache *cache,
	    BObolCompactOccurrenceStream *stream);
    friend void bobol_database_source_seed_realization_cache(
	    SoBRLDatabaseSource *source,
	    BObolDatabaseSourceRealizationCache *cache);

    static void fieldSensorCB(void *data, SoSensor *sensor);
    void attachFieldSensors(void);
    void detachFieldSensors(void);
    void syncRealizedShapeOwnerState(void);
    void syncCompactInstanceDisplayState(void);
    void rebuildCompactInstanceDisplayState(SbBool syncSourceState);
    int reapplyCompactInstanceVisibilityFrontier(size_t firstEntry = 0);
    int reapplyCompactInstanceSelectedPaths(size_t firstEntry = 0);
    int reapplyCompactInstancePresentationOverrides(size_t firstEntry = 0);
    void syncCompactInstancePlacementState(void);
    void retireCompactOverviewAfterPresentation(SoBRLCadAssembly *assembly);
    void seedCompactRealizationCache(
	BObolDatabaseSourceRealizationCache *cache) const;
    void clearCompiledAssembly(void);
    void ensureCompiledAssemblyChild(void);
    void markCompiledAssemblyDirty(void);
    void markCadBatchDirty(void);
    void markCadBatchDirty(const std::vector<size_t> &entryIndices);
    SbBool getCadBatchChangedEntries(uint64_t revision,
	std::vector<size_t> &entryIndices) const;
    void markDisplayMeshLodDirty(void);
    void markDisplayMeshLodDirty(
	const std::vector<size_t> &entryIndices,
	SbBool coverageInvalidated = FALSE);
    void markDisplayMeshLodVisibilityDirty(
	const std::vector<size_t> &entryIndices);
    uint64_t cadBatchRevisionGet(void) const;
    void clearCompactInstanceIndex(void);
    void discardCompactInstanceHistory(void);
    void installCompactInstanceIndex(
	struct BObolCompactInstanceIndex *index,
	SbBool occurrenceRegistry);
    int syncCompiledAssembly(void);
    uint64_t cadBatchStructureSignature(void) const;
    uint64_t cadBatchStyleSignature(void) const;
    uint64_t cadBatchSemanticSignature(void) const;
    int appendCadRenderBatch(struct BObolCadBatchBuildState *state,
	SbBool includeGeometry, SbBool includeSemantics);
    const struct BObolCompactInstanceEntry *findCompactInstanceEntry(
	const BObolCompactInstanceHandle &handle) const;
    int exportCompactInstances(SoBRLExportAction *action,
	const SbMatrix &parentToWorld);
    int measureCompactInstances(SoBRLMeasureAction *action,
	const SbMatrix &parentToWorld);
    int snapCompactInstances(SoBRLSnapAction *action,
	const SbMatrix &parentToWorld);

    struct Impl;
    std::unique_ptr<Impl> d;
};

/**
 * Refresh effective database material colors for a retained source set.
 * Combination state and shared path prefixes are resolved once per sweep;
 * callers should prefer this over invoking the single-source database lookup
 * repeatedly.
 */
BOBOL_EXPORT int
bobol_database_sources_refresh_material_colors(
    SoBRLDatabaseSource *const *sources,
    size_t sourceCount,
    uint32_t materialRevision,
    struct db_i *dbip);

#endif /* BOBOL_BDATABASESOURCE_H */
