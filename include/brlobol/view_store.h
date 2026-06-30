/*                    V I E W _ S T O R E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/view_store.h
 *
 * Obol-owned view feature, polygon, and selection stores.
 *
 * These classes are the C++ backing-store target for neutral RT/GED view
 * adapters.  They deliberately expose BRLObol/Inventor-friendly value types
 * instead of mirroring the transitional BSG C API.
 */

#ifndef BRLOBOL_VIEW_STORE_H
#define BRLOBOL_VIEW_STORE_H

#include "brlobol/defines.h"

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

class BRLObolFeatureStore;
class BRLObolPolygonStore;
class BRLObolSelectionStore;
class BRLObolViewController;
class SoNode;
struct bg_line_layer_builder;
struct db_i;
struct directory;

struct BRLOBOL_EXPORT BRLObolFeatureHandle {
    uint64_t id;
    uint64_t revision;

    BRLObolFeatureHandle(void);
    BRLObolFeatureHandle(uint64_t featureId, uint64_t featureRevision);
    SbBool isValid(void) const;
};

struct BRLOBOL_EXPORT BRLObolPolygonHandle {
    uint64_t id;
    uint64_t revision;

    BRLObolPolygonHandle(void);
    BRLObolPolygonHandle(uint64_t polygonId, uint64_t polygonRevision);
    SbBool isValid(void) const;
};

BRLOBOL_EXPORT SbBool
operator==(const BRLObolFeatureHandle &a, const BRLObolFeatureHandle &b);

BRLOBOL_EXPORT SbBool
operator!=(const BRLObolFeatureHandle &a, const BRLObolFeatureHandle &b);

BRLOBOL_EXPORT SbBool
operator==(const BRLObolPolygonHandle &a, const BRLObolPolygonHandle &b);

BRLOBOL_EXPORT SbBool
operator!=(const BRLObolPolygonHandle &a, const BRLObolPolygonHandle &b);

enum class BRLObolFeatureKind {
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
    HudLabel
};

enum class BRLObolFeatureScope {
    Shared = 0,
    Local = 1
};

static const unsigned int BRLOBOL_FEATURE_SCOPE_SHARED = 0x01u;
static const unsigned int BRLOBOL_FEATURE_SCOPE_LOCAL = 0x02u;
static const unsigned int BRLOBOL_FEATURE_SCOPE_ALL =
    BRLOBOL_FEATURE_SCOPE_SHARED | BRLOBOL_FEATURE_SCOPE_LOCAL;

enum class BRLObolLineCommand {
    Move = 0,
    Draw = 1,
    Point = 2
};

enum class BRLObolCommandResultStatus {
    None = 0,
    Accepted,
    Updated,
    Removed,
    Failed
};

enum class BRLObolPolygonType {
    General = 0,
    Circle = 1,
    Ellipse = 2,
    Rectangle = 3,
    Square = 4
};

enum class BRLObolPolygonUpdate {
    Default = 0,
    PropsOnly = 1,
    PointSelect = 2,
    PointSelectClear = 3,
    PointMove = 4,
    PointAppend = 5
};

struct BRLOBOL_EXPORT BRLObolFeatureStyle {
    SbBool hasVisible;
    SbBool visible;
    SbBool hasColor;
    SbColor color;
    SbBool hasLineWidth;
    int lineWidth;
    SbBool hasLineStyle;
    int lineStyle;
    SbBool hasArrow;
    SbBool arrow;
    SbBool hasArrowTip;
    float arrowTipLength;
    float arrowTipWidth;

    BRLObolFeatureStyle(void);
};

struct BRLOBOL_EXPORT BRLObolCommandResult {
    BRLObolCommandResultStatus status;
    BRLObolFeatureHandle feature;
    BRLObolPolygonHandle polygon;
    SbString command;
    SbString diagnostic;

    BRLObolCommandResult(void);
};

typedef void (*BRLObolCommandResultCallback)(
	const BRLObolCommandResult &result,
	void *userData);

struct BRLOBOL_EXPORT BRLObolFeatureOwner {
    const void *ownerToken;
    SbString ownerId;
    SbString ownerRole;
    BRLObolCommandResultCallback resultCallback;
    void *callbackUserData;

    BRLObolFeatureOwner(void);
};

struct BRLOBOL_EXPORT BRLObolEditPreviewCallbacks {
    void *previewContext;
    uint64_t (*revisionCallback)(void *previewContext);
    int (*updateCallback)(void *previewContext);
    int (*pickCallback)(void *previewContext, int x, int y, void *pickOut);

    BRLObolEditPreviewCallbacks(void);
};

struct BRLOBOL_EXPORT BRLObolLabel {
    SbString text;
    SbVec3f point;
    SbBool hasColor;
    SbColor color;
    SbBool hasLeader;
    SbVec3f target;
    int anchor;
    SbBool arrow;

    BRLObolLabel(void);
};

struct BRLOBOL_EXPORT BRLObolLineLayer {
    SbString name;
    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    BRLObolFeatureStyle style;

    BRLObolLineLayer(void);
};

struct BRLOBOL_EXPORT BRLObolFeatureSummary {
    SbBool exists;
    SbBool visible;
    SbBool realized;
    BRLObolFeatureKind kind;
    BRLObolFeatureScope scope;
    size_t pointCount;
    size_t commandCount;
    size_t childCount;
    BRLObolFeatureOwner owner;

    BRLObolFeatureSummary(void);
};

struct BRLOBOL_EXPORT BRLObolPolygonVisual {
    SbColor edgeColor;
    SbColor fillColor;
    SbBool fill;
    SbVec2f fillSlope;
    float fillSpacing;
    float viewZ;

    BRLObolPolygonVisual(void);
};

struct BRLOBOL_EXPORT BRLObolPolygonRecord {
    BRLObolPolygonHandle handle;
    SbString name;
    BRLObolFeatureScope scope;
    BRLObolPolygonType type;
    SbBool fill;
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
    void *userData;

    BRLObolPolygonRecord(void);
};

typedef int (*BRLObolPolygonRecordCallback)(
	const BRLObolPolygonRecord &record,
	void *userData);

struct BRLOBOL_EXPORT BRLObolSelectionRecord {
    SbString path;
    BRLObolFeatureHandle feature;
    int primitiveIndex;
    double hitDistance;

    BRLObolSelectionRecord(void);
};

class BRLOBOL_EXPORT BRLObolFeatureStore {
public:
    BRLObolFeatureStore(void);
    explicit BRLObolFeatureStore(BRLObolViewController *controller);
    ~BRLObolFeatureStore(void);

    void setController(BRLObolViewController *controller);
    BRLObolViewController *controller(void) const;

    void clear(void);
    BRLObolFeatureHandle find(const SbString &name,
	unsigned int scopeMask = BRLOBOL_FEATURE_SCOPE_ALL) const;
    SbBool exists(const SbString &name,
	unsigned int scopeMask = BRLOBOL_FEATURE_SCOPE_ALL) const;
    SbBool remove(const SbString &name);
    size_t removePrefix(const SbString &prefix);

    BRLObolFeatureHandle publishLineSet(const SbString &name,
	BRLObolFeatureScope scope,
	const std::vector<SbVec3f> &points,
	const std::vector<int32_t> &commands,
	const BRLObolFeatureStyle *style = NULL,
	const BRLObolFeatureOwner *owner = NULL);
    BRLObolFeatureHandle publishIndexedLineSet(const SbString &name,
	BRLObolFeatureScope scope,
	const std::vector<SbVec3f> &points,
	const std::vector<int32_t> &indices,
	const BRLObolFeatureStyle *style = NULL,
	const BRLObolFeatureOwner *owner = NULL);
    BRLObolFeatureHandle publishPointSet(const SbString &name,
	BRLObolFeatureScope scope,
	const std::vector<SbVec3f> &points,
	const BRLObolFeatureStyle *style = NULL,
	const BRLObolFeatureOwner *owner = NULL);
    BRLObolFeatureHandle publishLabels(const SbString &name,
	BRLObolFeatureScope scope,
	const std::vector<BRLObolLabel> &labels,
	const BRLObolFeatureStyle *style = NULL,
	const BRLObolFeatureOwner *owner = NULL);
    BRLObolFeatureHandle publishArrow(const SbString &name,
	BRLObolFeatureScope scope,
	const std::vector<SbVec3f> &points,
	const BRLObolFeatureStyle *style = NULL,
	const BRLObolFeatureOwner *owner = NULL);
    BRLObolFeatureHandle publishAxes(const SbString &name,
	BRLObolFeatureScope scope,
	const std::vector<SbVec3f> &centers,
	float halfAxesSize,
	const BRLObolFeatureStyle *style = NULL,
	const BRLObolFeatureOwner *owner = NULL);
    BRLObolFeatureHandle publishLineLayers(const SbString &name,
	BRLObolFeatureScope scope,
	const std::vector<BRLObolLineLayer> &layers,
	const BRLObolFeatureStyle *style = NULL,
	const BRLObolFeatureOwner *owner = NULL);
    BRLObolFeatureHandle publishLineLayerBuilder(const SbString &name,
	BRLObolFeatureScope scope,
	const struct bg_line_layer_builder *builder,
	const BRLObolFeatureStyle *style = NULL,
	const BRLObolFeatureOwner *owner = NULL);
    BRLObolFeatureHandle publishIndexedFaceSet(const SbString &name,
	BRLObolFeatureScope scope,
	const std::vector<SbVec3f> &points,
	const std::vector<SbVec3f> &normals,
	const std::vector<int32_t> &indices,
	const BRLObolFeatureStyle *style = NULL,
	const BRLObolFeatureOwner *owner = NULL);
    BRLObolFeatureHandle publishEditPreview(const SbString &name,
	const SbString &identity,
	const SbString &editIntentId,
	const SbString &editIntentRole,
	const std::vector<SbVec3f> &points,
	const std::vector<int32_t> &commands,
	uint32_t sourceRevision = 0,
	uint32_t inputsRevision = 0,
	const BRLObolEditPreviewCallbacks *callbacks = NULL,
	const BRLObolFeatureOwner *owner = NULL);
    uint64_t editPreviewRevision(BRLObolFeatureHandle handle) const;
    int updateEditPreview(BRLObolFeatureHandle handle);
    int pickEditPreview(BRLObolFeatureHandle handle,
	int x,
	int y,
	void *pickOut) const;

    SbBool appendLinePoint(BRLObolFeatureHandle handle, const SbVec3f &point);
    SbBool replaceLabels(BRLObolFeatureHandle handle,
	const std::vector<BRLObolLabel> &labels);
    SbBool clearGeometry(BRLObolFeatureHandle handle);
    SbBool points(BRLObolFeatureHandle handle,
	std::vector<SbVec3f> &pointsOut) const;
    SbBool commands(BRLObolFeatureHandle handle,
	std::vector<int32_t> &commandsOut) const;
    SbBool lineCommandAt(BRLObolFeatureHandle handle,
	size_t index,
	int32_t &commandOut) const;
    SbBool labels(BRLObolFeatureHandle handle,
	std::vector<BRLObolLabel> &labelsOut) const;
    SbBool axesCenters(BRLObolFeatureHandle handle,
	std::vector<SbVec3f> &centersOut,
	float *halfAxesSizeOut = NULL) const;

    SbBool applyStyle(BRLObolFeatureHandle handle,
	const BRLObolFeatureStyle &style,
	SbBool recursive = FALSE);
    SbBool style(BRLObolFeatureHandle handle,
	BRLObolFeatureStyle &styleOut) const;
    SbBool setVisible(BRLObolFeatureHandle handle, SbBool visible);
    SbBool setColor(BRLObolFeatureHandle handle, const SbColor &color);
    SbBool setLineWidth(BRLObolFeatureHandle handle, int lineWidth);
    SbBool arrowTip(BRLObolFeatureHandle handle,
	float &tipLength,
	float &tipWidth) const;
    SbBool setArrowTip(BRLObolFeatureHandle handle,
	float tipLength,
	float tipWidth);

    SbBool realize(BRLObolFeatureHandle handle, SbBool recursive = FALSE);
    SbBool summary(const SbString &name,
	BRLObolFeatureSummary &summaryOut,
	unsigned int scopeMask = BRLOBOL_FEATURE_SCOPE_ALL) const;
    SoNode *node(BRLObolFeatureHandle handle) const;

private:
    BRLObolFeatureStore(const BRLObolFeatureStore &);
    BRLObolFeatureStore &operator=(const BRLObolFeatureStore &);

    struct Impl;
    Impl *impl;
};

class BRLOBOL_EXPORT BRLObolPolygonStore {
public:
    BRLObolPolygonStore(void);
    explicit BRLObolPolygonStore(BRLObolViewController *controller);
    ~BRLObolPolygonStore(void);

    void setController(BRLObolViewController *controller);
    BRLObolViewController *controller(void) const;

    void clear(void);
    BRLObolPolygonHandle create(const SbString &name,
	BRLObolFeatureScope scope,
	BRLObolPolygonType type,
	const SbVec3f &originPoint,
	const fastf_t *viewPlane = NULL,
	float viewZ = 0.0f);
    BRLObolPolygonHandle find(const SbString &name,
	unsigned int scopeMask = BRLOBOL_FEATURE_SCOPE_ALL) const;
    BRLObolPolygonHandle selectAtModelPoint(const SbVec3f &point) const;
    BRLObolPolygonHandle duplicate(BRLObolPolygonHandle handle,
	const SbString &newName);

    SbBool update(BRLObolPolygonHandle handle,
	BRLObolPolygonUpdate update);
    SbBool updateScreenPoint(BRLObolPolygonHandle handle,
	int x,
	int y,
	BRLObolPolygonUpdate update);
    SbBool move(BRLObolPolygonHandle handle,
	const SbVec3f &currentPoint,
	const SbVec3f &previousPoint);
    SbBool rename(BRLObolPolygonHandle handle, const SbString &newName);
    SbBool remove(BRLObolPolygonHandle handle);

    SbBool record(BRLObolPolygonHandle handle,
	BRLObolPolygonRecord &recordOut) const;
    void visitRecords(BRLObolPolygonRecordCallback callback,
	void *userData = NULL) const;

    SbBool setCurrent(BRLObolPolygonHandle handle,
	long contour,
	long point);
    SbBool setContourOpen(BRLObolPolygonHandle handle,
	long contour,
	SbBool open);
    SbBool setAllContoursOpen(BRLObolPolygonHandle handle,
	SbBool open);
    SbBool clearSelectedPoint(BRLObolPolygonHandle handle);
    SbBool clearAllPointSelections(void);

    SbBool setVisual(BRLObolPolygonHandle handle,
	const BRLObolPolygonVisual &visual);
    SbBool visual(BRLObolPolygonHandle handle,
	BRLObolPolygonVisual &visualOut) const;
    SbBool setEdgeColor(BRLObolPolygonHandle handle,
	const SbColor &edgeColor);
    SbBool edgeColor(BRLObolPolygonHandle handle,
	SbColor &edgeColorOut) const;
    SbBool setFill(BRLObolPolygonHandle handle,
	SbBool fill,
	const SbVec2f &slope,
	float spacing);
    SbBool setFillColor(BRLObolPolygonHandle handle,
	const SbColor &fillColor);
    SbBool fillColor(BRLObolPolygonHandle handle,
	SbColor &fillColorOut) const;

    SbBool setGeometry(BRLObolPolygonHandle handle,
	const struct bg_polygon *polygon);
    const struct bg_polygon *geometry(BRLObolPolygonHandle handle) const;
    SbBool copyGeometry(BRLObolPolygonHandle handle,
	struct bg_polygon *polygonOut) const;
    SbBool setUserData(BRLObolPolygonHandle handle, void *userData);
    void *userData(BRLObolPolygonHandle handle) const;

    double area(BRLObolPolygonHandle handle, double viewScale) const;
    SbBool overlaps(BRLObolPolygonHandle a,
	BRLObolPolygonHandle b,
	const struct bn_tol &tol,
	double viewScale) const;
    SbBool csg(BRLObolPolygonHandle target,
	BRLObolPolygonHandle stencil,
	bg_clip_t op);

    BRLObolPolygonHandle importSketch(const SbString &name,
	BRLObolFeatureScope scope,
	struct db_i *dbip,
	struct directory *dp);
    SbBool exportSketch(BRLObolPolygonHandle handle,
	struct db_i *dbip,
	const SbString &name) const;

    size_t snapCount(BRLObolPolygonHandle exclude = BRLObolPolygonHandle()) const;
    SbBool setSnapExclude(BRLObolPolygonHandle handle);
    BRLObolPolygonHandle snapExclude(void) const;
    SoNode *node(BRLObolPolygonHandle handle) const;

private:
    BRLObolPolygonStore(const BRLObolPolygonStore &);
    BRLObolPolygonStore &operator=(const BRLObolPolygonStore &);

    struct Impl;
    Impl *impl;
};

class BRLOBOL_EXPORT BRLObolSelectionStore {
public:
    BRLObolSelectionStore(void);
    ~BRLObolSelectionStore(void);

    void clear(void);
    size_t count(void) const;
    SbBool containsPath(const SbString &path) const;
    SbBool addPath(const SbString &path);
    SbBool setPath(const SbString &path);
    SbBool removePath(const SbString &path);
    const BRLObolSelectionRecord *record(size_t index) const;
    SbBool addRecord(const BRLObolSelectionRecord &record);
    SbBool setRecords(const std::vector<BRLObolSelectionRecord> &records);
    void visitPaths(int (*callback)(const SbString &path, void *userData),
	void *userData = NULL) const;
    SbBool applyPickResults(const std::vector<BRLObolSelectionRecord> &records,
	void (*selectedPathCallback)(const SbString &path, void *userData) = NULL,
	void *userData = NULL);

private:
    BRLObolSelectionStore(const BRLObolSelectionStore &);
    BRLObolSelectionStore &operator=(const BRLObolSelectionStore &);

    struct Impl;
    Impl *impl;
};

#endif /* BRLOBOL_VIEW_STORE_H */
