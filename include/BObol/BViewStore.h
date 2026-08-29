/*                  B V I E W S T O R E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BViewStore.h
 *
 * Obol-owned view feature, polygon, and selection stores.
 *
 * These classes are the C++ backing-store target for neutral RT/GED view
 * adapters.  They deliberately expose BObol/Inventor-friendly value types
 * instead of mirroring the transitional BSG C API.
 */

#ifndef BOBOL_BVIEWSTORE_H
#define BOBOL_BVIEWSTORE_H

#include "BObol/BDefines.h"

#include "bg/polygon_types.h"
#include "bn/tol.h"
#include "bu/color.h"
#include "vmath.h"

#include <Inventor/SbBasic.h>
#include <Inventor/SbColor.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec2f.h>
#include <Inventor/SbVec3f.h>

#include <stddef.h>
#include <stdint.h>
#include <vector>

class BObolFeatureStore;
class BObolPolygonStore;
class BObolSelectionStore;
class BObolViewController;
class SoBRLDatabaseSource;
class SoNode;
struct BObolCompactInstanceHandle;
struct BObolCompactInstanceSummary;
struct bg_line_layer_builder;
struct db_i;
struct directory;

struct BOBOL_EXPORT BObolFeatureHandle {
    uint64_t id;
    uint64_t revision;

    BObolFeatureHandle(void);
    BObolFeatureHandle(uint64_t featureId, uint64_t featureRevision);
    SbBool isValid(void) const;
};

struct BOBOL_EXPORT BObolPolygonHandle {
    uint64_t id;
    uint64_t revision;

    BObolPolygonHandle(void);
    BObolPolygonHandle(uint64_t polygonId, uint64_t polygonRevision);
    SbBool isValid(void) const;
};

BOBOL_EXPORT SbBool
operator==(const BObolFeatureHandle &a, const BObolFeatureHandle &b);

BOBOL_EXPORT SbBool
operator!=(const BObolFeatureHandle &a, const BObolFeatureHandle &b);

BOBOL_EXPORT SbBool
operator==(const BObolPolygonHandle &a, const BObolPolygonHandle &b);

BOBOL_EXPORT SbBool
operator!=(const BObolPolygonHandle &a, const BObolPolygonHandle &b);

enum class BObolFeatureKind {
    Unknown = 0,
    Lines,
    IndexedLines,
    Points,
    Labels,
    Arrow,
    Axes,
    LineLayer,
    EditPreview,
    IndexedFaceSet,
    PolygonOverlay,
    HudLabel,
    CustomNode
};

enum class BObolFeatureScope {
    Shared = 0,
    Local = 1
};

static const unsigned int BOBOL_FEATURE_SCOPE_SHARED = 0x01u;
static const unsigned int BOBOL_FEATURE_SCOPE_LOCAL = 0x02u;
static const unsigned int BOBOL_FEATURE_SCOPE_ALL =
    BOBOL_FEATURE_SCOPE_SHARED | BOBOL_FEATURE_SCOPE_LOCAL;

static const int BOBOL_SELECTION_ALL = -1;
static const int BOBOL_SELECTION_SELECTED_PATH = 0;
static const int BOBOL_SELECTION_HIGHLIGHTED_REF = 1;

enum class BObolLineCommand {
    Move = 0,
    Draw = 1,
    Point = 2
};

enum class BObolCommandResultStatus {
    None = 0,
    Accepted,
    Updated,
    Removed,
    Failed
};

enum class BObolPolygonType {
    General = 0,
    Circle = 1,
    Ellipse = 2,
    Rectangle = 3,
    Square = 4
};

enum class BObolPolygonUpdate {
    Default = 0,
    PropsOnly = 1,
    PointSelect = 2,
    PointSelectClear = 3,
    PointMove = 4,
    PointAppend = 5,
    PointDelete = 6
};

enum BObolPolygonFillFlag {
    BOBOL_POLYGON_FILL_NONE = 0,
    BOBOL_POLYGON_FILL_HATCH = 1,
    BOBOL_POLYGON_FILL_MESH = 2
};

struct BOBOL_EXPORT BObolFeatureStyle {
    SbBool hasVisible;
    SbBool visible;
    SbBool hasSelectable;
    SbBool selectable;
    SbBool hasColor;
    SbColor color;
    SbBool hasLineWidth;
    int lineWidth;
    SbBool hasLineStyle;
    int lineStyle;
    SbBool hasTransparency;
    float transparency;
    SbBool hasArrow;
    SbBool arrow;
    SbBool hasArrowTip;
    float arrowTipLength;
    float arrowTipWidth;

    /* When TRUE the feature is a screen-locked HUD overlay: its line geometry
     * is expressed in pixel coordinates (origin bottom-left, 0..width/0..height)
     * and is rendered through an SoHUDKit rather than the model/view pipeline.
     * Used by the faceplate GUI so its lines stay fixed on screen like the
     * HUD text labels published through the view-feature batch adapter. */
    SbBool hud;

    BObolFeatureStyle(void);
};

struct BOBOL_EXPORT BObolCommandResult {
    BObolCommandResultStatus status;
    BObolFeatureHandle feature;
    BObolPolygonHandle polygon;
    SbString command;
    SbString diagnostic;

    BObolCommandResult(void);
};

typedef void (*BObolCommandResultCallback)(
	const BObolCommandResult &result,
	void *userData);

struct BOBOL_EXPORT BObolFeatureOwner {
    const void *ownerToken;
    SbString ownerId;
    SbString ownerRole;
    uint64_t generation;
    BObolCommandResultCallback resultCallback;
    void *callbackUserData;

    BObolFeatureOwner(void);
};

enum class BObolOverlayRole {
    None,
    Model,
    Screen,
    XRay
};

enum class BObolOverlayClass {
    None,
    Faceplate,
    EditHandle,
    Measure,
    SelectionRubberBand,
    SnapGuide,
    CommandResult,
    Diagnostic,
    TclOverlay,
    PolygonEdit,
    SketchEdit,
    UserAnnotation
};

enum class BObolOverlayLifecycle {
    None,
    Persistent,
    PerFrame,
    PerCommand,
    PerTool,
    PerView,
    SharedViewSet,
    AutoRemoveOnSource
};

enum class BObolOverlayOrder {
    Model,
    Screen,
    XRay,
    PostTransparent
};

struct BOBOL_EXPORT BObolOverlayInfo {
    SbBool isOverlay;
    const void *ownerToken;
    BObolOverlayRole role;
    BObolOverlayClass overlayClass;
    BObolOverlayLifecycle lifecycle;
    BObolOverlayOrder order;
    int sortOrder;
    SbString sourcePath;

    BObolOverlayInfo(void);
};

struct BOBOL_EXPORT BObolLabel {
    SbString text;
    SbVec3f point;
    SbBool hasColor;
    SbColor color;
    SbBool hasLeader;
    SbVec3f target;
    int anchor;
    SbBool arrow;
    float fontSize;
    uint32_t sourceId;

    BObolLabel(void);
};

struct BOBOL_EXPORT BObolLineLayer {
    SbString name;
    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    BObolFeatureStyle style;

    BObolLineLayer(void);
};

struct BOBOL_EXPORT BObolFeatureMetadata {
    SbString key;
    SbString value;

    BObolFeatureMetadata(void);
};

struct BOBOL_EXPORT BObolFeaturePrimitiveMetadata {
    int32_t primitiveIndex;
    std::vector<BObolFeatureMetadata> metadata;

    BObolFeaturePrimitiveMetadata(void);
};

struct BOBOL_EXPORT BObolFeaturePrimitivePick {
    BObolFeatureHandle handle;
    SbString featureName;
    int32_t primitiveIndex;
    std::vector<BObolFeatureMetadata> metadata;

    BObolFeaturePrimitivePick(void);
};

struct BOBOL_EXPORT BObolFeatureSummary {
    SbBool exists;
    SbBool visible;
    SbBool realized;
    BObolFeatureKind kind;
    BObolFeatureScope scope;
    size_t pointCount;
    size_t commandCount;
    size_t childCount;
    size_t metadataCount;
    size_t primitiveMetadataCount;
    size_t selectedPrimitiveCount;
    size_t highlightedPrimitiveCount;
    BObolFeatureOwner owner;
    BObolOverlayInfo overlay;

    BObolFeatureSummary(void);
};

struct BOBOL_EXPORT BObolFeatureRecord {
    BObolFeatureHandle handle;
    SbString name;
    BObolFeatureKind kind;
    BObolFeatureScope scope;
    BObolFeatureStyle style;
    BObolFeatureOwner owner;
    BObolOverlayInfo overlay;
    SbBool realized;
    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    std::vector<int32_t> indices;
    std::vector<SbVec3f> normals;
    std::vector<BObolLabel> labels;
    std::vector<SbVec3f> axesCenters;
    float halfAxesSize;
    std::vector<BObolLineLayer> layers;
    std::vector<BObolFeatureMetadata> metadata;
    std::vector<BObolFeaturePrimitiveMetadata> primitiveMetadata;
    std::vector<int32_t> selectedPrimitives;
    std::vector<int32_t> highlightedPrimitives;
    SbString identity;
    SbString editIntentId;
    SbString editIntentRole;
    uint32_t sourceRevision;
    uint32_t inputsRevision;

    BObolFeatureRecord(void);
};

typedef int (*BObolFeatureRecordCallback)(
	const BObolFeatureRecord &record,
	void *userData);

/** Allocation-free retained-node visitor for owner-thread scene queries. */
typedef int (*BObolFeatureNodeCallback)(
	BObolFeatureHandle handle,
	SoNode *node,
	void *userData);

struct BOBOL_EXPORT BObolPolygonVisual {
    SbColor edgeColor;
    SbColor fillColor;
    SbBool fill;
    unsigned int fillFlags;
    SbVec2f fillSlope;
    float fillSpacing;
    float viewZ;

    BObolPolygonVisual(void);
};

struct BOBOL_EXPORT BObolPolygonRecord {
    BObolPolygonHandle handle;
    SbString name;
    BObolFeatureScope scope;
    BObolPolygonType type;
    SbBool selected;
    SbBool fill;
    unsigned int fillFlags;
    SbVec2f fillSlope;
    float fillSpacing;
    SbColor fillColor;
    SbColor edgeColor;
    long currentContour;
    long currentPoint;
    SbBool firstContourOpen;
    size_t contourCount;
    size_t pointCount;
    SbVec3f originPoint;
    plane_t viewPlane;
    float viewZ;
    SbString sketchName;
    void *userData;

    BObolPolygonRecord(void);
};

typedef int (*BObolPolygonRecordCallback)(
	const BObolPolygonRecord &record,
	void *userData);

struct BOBOL_EXPORT BObolSelectionRecord {
    SbString path;
    BObolFeatureHandle feature;
    BObolFeatureOwner owner;
    int kind;
    int primitiveIndex;
    double hitDistance;

    BObolSelectionRecord(void);
};

class BOBOL_EXPORT BObolFeatureStore {
public:
    BObolFeatureStore(void);
    explicit BObolFeatureStore(BObolViewController *controller);
    ~BObolFeatureStore(void);

    void setController(BObolViewController *controller);
    BObolViewController *controller(void) const;

    /**
     * Stable identity for references issued by this store instance.
     *
     * The identity changes when a controller/store is replaced.  Callers
     * must resolve it on the store owner thread before using an object id;
     * it is not a pointer and does not retain the store.
     */
    uint64_t referenceGeneration(void) const;

    /**
     * Monotonic revision of realized presentation content.
     *
     * Unlike a controller render-request serial, this advances for every
     * actual retained-content mutation even when several mutations are
     * coalesced into one frame.  Shared-scene hosts use it to wake every view
     * which renders this store's scene root.
     */
    uint64_t presentationRevision(void) const;

    void clear(void);
    BObolFeatureHandle find(const SbString &name,
	unsigned int scopeMask = BOBOL_FEATURE_SCOPE_ALL) const;
    BObolFeatureHandle findOwned(const SbString &name,
	unsigned int scopeMask,
	const BObolFeatureOwner *owner) const;
    SbBool exists(const SbString &name,
	unsigned int scopeMask = BOBOL_FEATURE_SCOPE_ALL) const;
    SbBool existsOwned(const SbString &name,
	unsigned int scopeMask,
	const BObolFeatureOwner *owner) const;
    SbBool remove(BObolFeatureHandle handle);
    SbBool remove(const SbString &name);
    SbBool removeOwned(const SbString &name,
	unsigned int scopeMask,
	const BObolFeatureOwner *owner);
    size_t removeScope(unsigned int scopeMask,
	const BObolFeatureOwner *owner = NULL);
    size_t removePrefix(const SbString &prefix);
    size_t removePrefix(const SbString &prefix,
	unsigned int scopeMask,
	const BObolFeatureOwner *owner);
    size_t countPrefix(const SbString &prefix,
	unsigned int scopeMask = BOBOL_FEATURE_SCOPE_ALL,
	const BObolFeatureOwner *owner = NULL) const;
    void markCommandOwnerGeneration(const BObolFeatureOwner &owner);
    SbBool commandOwnerGenerationCurrent(
	const BObolFeatureOwner &owner) const;

    BObolFeatureHandle publishLineSet(const SbString &name,
	BObolFeatureScope scope,
	const std::vector<SbVec3f> &points,
	const std::vector<int32_t> &commands,
	const BObolFeatureStyle *style = NULL,
	const BObolFeatureOwner *owner = NULL);
    BObolFeatureHandle publishIndexedLineSet(const SbString &name,
	BObolFeatureScope scope,
	const std::vector<SbVec3f> &points,
	const std::vector<int32_t> &indices,
	const BObolFeatureStyle *style = NULL,
	const BObolFeatureOwner *owner = NULL);
    BObolFeatureHandle publishPointSet(const SbString &name,
	BObolFeatureScope scope,
	const std::vector<SbVec3f> &points,
	const BObolFeatureStyle *style = NULL,
	const BObolFeatureOwner *owner = NULL);
    BObolFeatureHandle publishLabels(const SbString &name,
	BObolFeatureScope scope,
	const std::vector<BObolLabel> &labels,
	const BObolFeatureStyle *style = NULL,
	const BObolFeatureOwner *owner = NULL);
    BObolFeatureHandle publishHudLabels(const SbString &name,
	BObolFeatureScope scope,
	const std::vector<BObolLabel> &labels,
	const BObolFeatureStyle *style = NULL,
	const BObolFeatureOwner *owner = NULL);
    BObolFeatureHandle publishArrow(const SbString &name,
	BObolFeatureScope scope,
	const std::vector<SbVec3f> &points,
	const BObolFeatureStyle *style = NULL,
	const BObolFeatureOwner *owner = NULL);
    BObolFeatureHandle publishAxes(const SbString &name,
	BObolFeatureScope scope,
	const std::vector<SbVec3f> &centers,
	float halfAxesSize,
	const BObolFeatureStyle *style = NULL,
	const BObolFeatureOwner *owner = NULL);
    BObolFeatureHandle publishLineLayers(const SbString &name,
	BObolFeatureScope scope,
	const std::vector<BObolLineLayer> &layers,
	const BObolFeatureStyle *style = NULL,
	const BObolFeatureOwner *owner = NULL);
    BObolFeatureHandle publishLineLayerBuilder(const SbString &name,
	BObolFeatureScope scope,
	const struct bg_line_layer_builder *builder,
	const BObolFeatureStyle *style = NULL,
	const BObolFeatureOwner *owner = NULL);
    BObolFeatureHandle publishIndexedFaceSet(const SbString &name,
	BObolFeatureScope scope,
	const std::vector<SbVec3f> &points,
	const std::vector<SbVec3f> &normals,
	const std::vector<int32_t> &indices,
	const BObolFeatureStyle *style = NULL,
	const BObolFeatureOwner *owner = NULL);

    /**
     * Patch indexed-face vertices without rebuilding or recopying topology.
     * Intended for retained edit presenters whose authoritative geometry is
     * owned elsewhere.  Each entry in pointIndices selects the corresponding
     * entry in points; all indices are validated before any mutation.
     */
    SbBool updateIndexedFaceSetPoints(BObolFeatureHandle handle,
	const std::vector<int32_t> &pointIndices,
	const std::vector<SbVec3f> &points);
    BObolFeatureHandle publishCustomNode(const SbString &name,
	BObolFeatureScope scope,
	SoNode *node,
	const BObolFeatureStyle *style = NULL,
	const BObolFeatureOwner *owner = NULL);
    BObolFeatureHandle publishEditPreview(const SbString &name,
	const SbString &identity,
	const SbString &editIntentId,
	const SbString &editIntentRole,
	const std::vector<SbVec3f> &points,
	const std::vector<int32_t> &commands,
	uint32_t sourceRevision = 0,
	uint32_t inputsRevision = 0,
	const BObolFeatureOwner *owner = NULL);
    /**
     * Materialize one compact occurrence as a transient local edit preview.
     *
     * Geometry is copied into source coordinates; identity, effective
     * appearance, placement revisions, and LoD provenance remain available
     * through compactEditBinding().  Neither operation changes or retains a
     * render node in the compact database source.
     */
    BObolFeatureHandle promoteCompactInstanceForEdit(
	const SbString &name,
	const SoBRLDatabaseSource &source,
	const BObolCompactInstanceHandle &instance,
	const SbString &editIntentId,
	const SbString &editIntentRole,
	const BObolFeatureOwner *owner = NULL);
    SbBool demoteCompactInstanceFromEdit(
	BObolFeatureHandle preview,
	const BObolCompactInstanceHandle &instance);
    SbBool compactEditBinding(BObolFeatureHandle preview,
	BObolCompactInstanceHandle &instanceOut,
	BObolCompactInstanceSummary &summaryOut) const;
    SbBool replaceEditPreviewGeometry(BObolFeatureHandle handle,
	const SbString &identity,
	const std::vector<SbVec3f> &points,
	const std::vector<int32_t> &commands,
	uint32_t sourceRevision = 0,
	uint32_t inputsRevision = 0);
    SbBool touch(BObolFeatureHandle handle);

    SbBool appendLinePoint(BObolFeatureHandle handle, const SbVec3f &point);
    SbBool replaceLabels(BObolFeatureHandle handle,
	const std::vector<BObolLabel> &labels);
    SbBool clearGeometry(BObolFeatureHandle handle);
    SbBool points(BObolFeatureHandle handle,
	std::vector<SbVec3f> &pointsOut) const;
    SbBool commands(BObolFeatureHandle handle,
	std::vector<int32_t> &commandsOut) const;
    SbBool lineCommandAt(BObolFeatureHandle handle,
	size_t index,
	int32_t &commandOut) const;
    SbBool labels(BObolFeatureHandle handle,
	std::vector<BObolLabel> &labelsOut) const;
    SbBool axesCenters(BObolFeatureHandle handle,
	std::vector<SbVec3f> &centersOut,
	float *halfAxesSizeOut = NULL) const;
    SbBool indices(BObolFeatureHandle handle,
	std::vector<int32_t> &indicesOut) const;
    SbBool normals(BObolFeatureHandle handle,
	std::vector<SbVec3f> &normalsOut) const;

    SbBool applyStyle(BObolFeatureHandle handle,
	const BObolFeatureStyle &style,
	SbBool recursive = FALSE);
    SbBool style(BObolFeatureHandle handle,
	BObolFeatureStyle &styleOut) const;
    SbBool setVisible(BObolFeatureHandle handle, SbBool visible);
    SbBool setColor(BObolFeatureHandle handle, const SbColor &color);
    SbBool setLineWidth(BObolFeatureHandle handle, int lineWidth);
    SbBool arrowTip(BObolFeatureHandle handle,
	float &tipLength,
	float &tipWidth) const;
    SbBool setArrowTip(BObolFeatureHandle handle,
	float tipLength,
	float tipWidth);
    SbBool setOverlayInfo(BObolFeatureHandle handle,
	const BObolOverlayInfo &overlay);
    SbBool clearOverlayInfo(BObolFeatureHandle handle);
    SbBool overlayInfo(BObolFeatureHandle handle,
	BObolOverlayInfo &overlayOut) const;
    SbBool replaceMetadata(BObolFeatureHandle handle,
	const std::vector<BObolFeatureMetadata> &metadata);
    SbBool metadata(BObolFeatureHandle handle,
	std::vector<BObolFeatureMetadata> &metadataOut) const;
    SbBool replacePrimitiveMetadata(BObolFeatureHandle handle,
	int32_t primitiveIndex,
	const std::vector<BObolFeatureMetadata> &metadata);
    SbBool primitiveMetadata(BObolFeatureHandle handle,
	int32_t primitiveIndex,
	std::vector<BObolFeatureMetadata> &metadataOut) const;
    SbBool resolvePrimitivePick(const SbString &name,
	int32_t primitiveIndex,
	BObolFeaturePrimitivePick &pickOut,
	unsigned int scopeMask = BOBOL_FEATURE_SCOPE_ALL,
	const BObolFeatureOwner *owner = NULL) const;
    SbBool replaceSelectedPrimitives(BObolFeatureHandle handle,
	const std::vector<int32_t> &primitives);
    SbBool replaceHighlightedPrimitives(BObolFeatureHandle handle,
	const std::vector<int32_t> &primitives);
    SbBool selectedPrimitives(BObolFeatureHandle handle,
	std::vector<int32_t> &primitivesOut) const;
    SbBool highlightedPrimitives(BObolFeatureHandle handle,
	std::vector<int32_t> &primitivesOut) const;

    SbBool realize(BObolFeatureHandle handle, SbBool recursive = FALSE);
    SbBool summary(const SbString &name,
	BObolFeatureSummary &summaryOut,
	unsigned int scopeMask = BOBOL_FEATURE_SCOPE_ALL) const;
    /**
     * Summarize one exact retained feature.  Unlike the name-based overload,
     * this preserves the identity selected by find()/findOwned() when
     * multiple producers use the same display name.
     */
    SbBool summary(BObolFeatureHandle handle,
	BObolFeatureSummary &summaryOut) const;
    SbBool summaryOwned(const SbString &name,
	BObolFeatureSummary &summaryOut,
	unsigned int scopeMask,
	const BObolFeatureOwner *owner) const;
    SbBool record(BObolFeatureHandle handle,
	BObolFeatureRecord &recordOut) const;
    void visitRecords(BObolFeatureRecordCallback callback,
	void *userData = NULL,
	unsigned int scopeMask = BOBOL_FEATURE_SCOPE_ALL,
	const BObolFeatureOwner *owner = NULL) const;
    /** Visit realized feature roots without copying their geometry records. */
    void visitNodes(BObolFeatureNodeCallback callback,
	void *userData = NULL) const;
    SoNode *node(BObolFeatureHandle handle) const;

private:
    BObolFeatureStore(const BObolFeatureStore &);
    BObolFeatureStore &operator=(const BObolFeatureStore &);

    struct Impl;
    Impl *impl;
};

class BOBOL_EXPORT BObolPolygonStore {
public:
    BObolPolygonStore(void);
    explicit BObolPolygonStore(BObolViewController *controller);
    ~BObolPolygonStore(void);

    void setController(BObolViewController *controller);
    BObolViewController *controller(void) const;

    /** See BObolFeatureStore::referenceGeneration(). */
    uint64_t referenceGeneration(void) const;

    /** See BObolFeatureStore::presentationRevision(). */
    uint64_t presentationRevision(void) const;

    /**
     * Return the store-owned name for a live polygon handle.
     *
     * The returned pointer remains valid until the polygon is renamed or
     * removed, or until the store is destroyed.  The call must be made on
     * the store owner thread.
     */
    const char *name(BObolPolygonHandle handle) const;

    void clear(void);
    BObolPolygonHandle create(const SbString &name,
	BObolFeatureScope scope,
	BObolPolygonType type,
	const SbVec3f &originPoint,
	const fastf_t *viewPlane = NULL,
	float viewZ = 0.0f);
    BObolPolygonHandle find(const SbString &name,
	unsigned int scopeMask = BOBOL_FEATURE_SCOPE_ALL) const;
    BObolPolygonHandle selectAtModelPoint(const SbVec3f &point) const;
    BObolPolygonHandle duplicate(BObolPolygonHandle handle,
	const SbString &newName);

    SbBool update(BObolPolygonHandle handle,
	BObolPolygonUpdate update);
    SbBool updateScreenPoint(BObolPolygonHandle handle,
	int x,
	int y,
	BObolPolygonUpdate update);
    SbBool updateModelPoint(BObolPolygonHandle handle,
	const SbVec3f &point,
	BObolPolygonUpdate update);
    SbBool move(BObolPolygonHandle handle,
	const SbVec3f &currentPoint,
	const SbVec3f &previousPoint);
    SbBool rename(BObolPolygonHandle handle, const SbString &newName);
    SbBool remove(BObolPolygonHandle handle);
    size_t removeScope(unsigned int scopeMask);

    SbBool record(BObolPolygonHandle handle,
	BObolPolygonRecord &recordOut) const;
    void visitRecords(BObolPolygonRecordCallback callback,
	void *userData = NULL) const;

    SbBool setCurrent(BObolPolygonHandle handle,
	long contour,
	long point);
    SbBool setSelected(BObolPolygonHandle handle, SbBool selected);
    SbBool clearSelection(void);
    SbBool setContourOpen(BObolPolygonHandle handle,
	long contour,
	SbBool open);
    SbBool setAllContoursOpen(BObolPolygonHandle handle,
	SbBool open);
    SbBool clearSelectedPoint(BObolPolygonHandle handle);
    SbBool clearAllPointSelections(void);

    SbBool setVisible(BObolPolygonHandle handle, SbBool visible);
    SbBool isVisible(BObolPolygonHandle handle, SbBool &visibleOut) const;
    SbBool setVisual(BObolPolygonHandle handle,
	const BObolPolygonVisual &visual);
    SbBool visual(BObolPolygonHandle handle,
	BObolPolygonVisual &visualOut) const;
    SbBool setEdgeColor(BObolPolygonHandle handle,
	const SbColor &edgeColor);
    SbBool edgeColor(BObolPolygonHandle handle,
	SbColor &edgeColorOut) const;
    SbBool setFill(BObolPolygonHandle handle,
	SbBool fill,
	const SbVec2f &slope,
	float spacing);
    SbBool setFillFlags(BObolPolygonHandle handle,
	unsigned int fillFlags);
    SbBool setFillColor(BObolPolygonHandle handle,
	const SbColor &fillColor);
    SbBool fillColor(BObolPolygonHandle handle,
	SbColor &fillColorOut) const;

    SbBool setGeometry(BObolPolygonHandle handle,
	const struct bg_polygon *polygon);
    const struct bg_polygon *geometry(BObolPolygonHandle handle) const;
    SbBool copyGeometry(BObolPolygonHandle handle,
	struct bg_polygon *polygonOut) const;
    SbBool setSketchName(BObolPolygonHandle handle,
	const SbString &sketchName);
    const char *sketchName(BObolPolygonHandle handle) const;
    SbBool setUserData(BObolPolygonHandle handle, void *userData);
    void *userData(BObolPolygonHandle handle) const;

    double area(BObolPolygonHandle handle, double viewScale) const;
    SbBool overlaps(BObolPolygonHandle a,
	BObolPolygonHandle b,
	const struct bn_tol &tol,
	double viewScale) const;
    SbBool csg(BObolPolygonHandle target,
	BObolPolygonHandle stencil,
	enum bg_polygon_boolean_op op);

    BObolPolygonHandle importSketch(const SbString &name,
	BObolFeatureScope scope,
	struct db_i *dbip,
	struct directory *dp);
    SbBool exportSketch(BObolPolygonHandle handle,
	struct db_i *dbip,
	const SbString &name) const;
    SbBool updateSketch(BObolPolygonHandle handle,
	struct db_i *dbip,
	const SbString &name) const;

    size_t snapCount(BObolPolygonHandle exclude = BObolPolygonHandle()) const;
    /** Exclude one actively edited polygon from snap queries.  An invalid
     * handle clears the exclusion. */
    SbBool setSnapExclude(BObolPolygonHandle handle);
    BObolPolygonHandle snapExclude(void) const;
    SoNode *node(BObolPolygonHandle handle) const;

private:
    BObolPolygonStore(const BObolPolygonStore &);
    BObolPolygonStore &operator=(const BObolPolygonStore &);

    struct Impl;
    Impl *impl;
};

class BOBOL_EXPORT BObolSelectionStore {
public:
    BObolSelectionStore(void);
    ~BObolSelectionStore(void);

    void clear(const BObolFeatureOwner *owner = NULL,
	int kind = BOBOL_SELECTION_ALL);
    size_t count(const BObolFeatureOwner *owner = NULL,
	int kind = BOBOL_SELECTION_ALL) const;
    SbBool containsPath(const SbString &path,
	int kind = BOBOL_SELECTION_SELECTED_PATH,
	const BObolFeatureOwner *owner = NULL) const;
    SbBool addPath(const SbString &path,
	int kind = BOBOL_SELECTION_SELECTED_PATH,
	const BObolFeatureOwner *owner = NULL);
    SbBool setPath(const SbString &path,
	int kind = BOBOL_SELECTION_SELECTED_PATH,
	const BObolFeatureOwner *owner = NULL);
    SbBool removePath(const SbString &path,
	int kind = BOBOL_SELECTION_ALL,
	const BObolFeatureOwner *owner = NULL);
    /** Apply a path delta with one linear store pass.  This is the scalable
     * publication path for large semantic selections; repeated addPath calls
     * would otherwise perform a growing linear duplicate scan. */
    SbBool applyPathDelta(const std::vector<SbString> &addedPaths,
	const std::vector<SbString> &removedPaths,
	int kind = BOBOL_SELECTION_SELECTED_PATH,
	const BObolFeatureOwner *owner = NULL);
    const BObolSelectionRecord *record(size_t index) const;
    SbBool addRecord(const BObolSelectionRecord &record);
    SbBool setRecords(const std::vector<BObolSelectionRecord> &records);
    void visitPaths(int (*callback)(const SbString &path, void *userData),
	void *userData = NULL,
	const BObolFeatureOwner *owner = NULL,
	int kind = BOBOL_SELECTION_ALL) const;
    SbBool applyPickResults(const std::vector<BObolSelectionRecord> &records,
	void (*selectedPathCallback)(const SbString &path, void *userData) = NULL,
	void *userData = NULL,
	const BObolFeatureOwner *owner = NULL);

private:
    BObolSelectionStore(const BObolSelectionStore &);
    BObolSelectionStore &operator=(const BObolSelectionStore &);

    struct Impl;
    Impl *impl;
};

#endif /* BOBOL_BVIEWSTORE_H */
