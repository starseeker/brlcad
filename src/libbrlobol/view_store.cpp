/*                  V I E W _ S T O R E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file view_store.cpp */

#include "common.h"

#include "brlobol/axes.h"
#include "brlobol/edit_preview.h"
#include "brlobol/line_layer_overlay.h"
#include "brlobol/lod_realization.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/scene_group.h"
#include "brlobol/view_controller.h"
#include "brlobol/view_store.h"
#include "brlobol/vlist_shape.h"

#include "bg/line_layer.h"
#include "bg/plane.h"
#include "bg/polygon.h"
#include "bn/tol.h"
#include "rt/primitives/sketch.h"
#include "bu/malloc.h"

#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoFont.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoText2.h>
#include <Inventor/nodes/SoTranslation.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>

static std::string
store_string(const SbString &s)
{
    const char *str = s.getString();
    return str ? std::string(str) : std::string();
}

static std::string
store_owner_key(const BRLObolFeatureOwner *owner)
{
    if (!owner)
	return std::string();

    if (owner->ownerToken) {
	char buf[64] = {0};
	snprintf(buf, sizeof(buf), "T:%p", owner->ownerToken);
	return std::string(buf);
    }

    const char *id = owner->ownerId.getString();
    if (id && id[0])
	return std::string("I:") + id;

    return std::string();
}

static std::string
store_key(BRLObolFeatureScope scope,
	const SbString &name,
	const BRLObolFeatureOwner *owner = NULL)
{
    if (scope != BRLObolFeatureScope::Local)
	return std::string("S:") + store_string(name);

    return std::string("L:") + store_owner_key(owner) + ":" +
	store_string(name);
}

static SbBool
store_owner_matches(const BRLObolFeatureOwner &recordOwner,
	const BRLObolFeatureOwner *queryOwner)
{
    if (!queryOwner)
	return TRUE;

    if (queryOwner->ownerToken)
	return recordOwner.ownerToken == queryOwner->ownerToken ? TRUE : FALSE;

    const char *queryId = queryOwner->ownerId.getString();
    if (queryId && queryId[0]) {
	const char *recordId = recordOwner.ownerId.getString();
	return recordId && strcmp(recordId, queryId) == 0 ? TRUE : FALSE;
    }

    const char *queryRole = queryOwner->ownerRole.getString();
    if (queryRole && queryRole[0]) {
	const char *recordRole = recordOwner.ownerRole.getString();
	return recordRole && strcmp(recordRole, queryRole) == 0 ? TRUE : FALSE;
    }

    return TRUE;
}

static SbBool
store_selection_kind_matches(int recordKind, int queryKind)
{
    return queryKind == BRLOBOL_SELECTION_ALL || recordKind == queryKind ?
	TRUE : FALSE;
}

static int
store_selection_record_kind(int kind)
{
    return kind == BRLOBOL_SELECTION_ALL ?
	BRLOBOL_SELECTION_SELECTED_PATH : kind;
}

static unsigned int
store_scope_bit(BRLObolFeatureScope scope)
{
    return scope == BRLObolFeatureScope::Local ?
	BRLOBOL_FEATURE_SCOPE_LOCAL : BRLOBOL_FEATURE_SCOPE_SHARED;
}

static SbVec3f
store_vec3(const point_t p)
{
    return SbVec3f(static_cast<float>(p[X]),
	    static_cast<float>(p[Y]),
	    static_cast<float>(p[Z]));
}

static void
store_point(point_t p, const SbVec3f &v)
{
    VSET(p, v[0], v[1], v[2]);
}

static int32_t
store_shape_command(int32_t command)
{
    if (command == static_cast<int32_t>(BRLObolLineCommand::Point))
	return SoBRLVListShape::POINT;
    if (command == static_cast<int32_t>(BRLObolLineCommand::Draw))
	return SoBRLVListShape::DRAW;
    return SoBRLVListShape::MOVE;
}

static std::vector<int32_t>
store_normalized_line_commands(const std::vector<SbVec3f> &points,
	const std::vector<int32_t> &commands)
{
    std::vector<int32_t> normalized = commands;

    if (normalized.size() != points.size()) {
	normalized.assign(points.size(), static_cast<int32_t>(
		BRLObolLineCommand::Draw));
	if (!normalized.empty())
	    normalized[0] = static_cast<int32_t>(BRLObolLineCommand::Move);
    }

    return normalized;
}

static std::vector<int32_t>
store_shape_commands(const std::vector<int32_t> &commands)
{
    std::vector<int32_t> shapeCommands;
    shapeCommands.reserve(commands.size());
    for (size_t i = 0; i < commands.size(); i++)
	shapeCommands.push_back(store_shape_command(commands[i]));
    return shapeCommands;
}

static void
store_apply_vlist_style(SoBRLVListShape *shape,
	const BRLObolFeatureStyle &style)
{
    if (!shape)
	return;
    if (style.hasVisible)
	shape->visible = style.visible;
    if (style.hasColor) {
	shape->colorOverride = TRUE;
	shape->color = style.color;
    }
    if (style.hasLineWidth)
	shape->lineWidth = style.lineWidth;
    if (style.hasLineStyle)
	shape->lineStyle = style.lineStyle;
}

static void
store_apply_mesh_style(SoBRLMeshShape *shape,
	const BRLObolFeatureStyle &style)
{
    if (!shape)
	return;
    if (style.hasVisible)
	shape->visible = style.visible;
    if (style.hasColor) {
	shape->colorOverride = TRUE;
	shape->color = style.color;
    }
    if (style.hasLineWidth)
	shape->lineWidth = style.lineWidth;
    if (style.hasLineStyle)
	shape->lineStyle = style.lineStyle;
}

static void
store_apply_mesh_color(SoBRLMeshShape *shape, const SbColor &color)
{
    if (!shape)
	return;
    shape->colorOverride = TRUE;
    shape->color = color;
}

static void
store_sbcolor_to_bu(const SbColor &src, struct bu_color *dst)
{
    if (!dst)
	return;

    dst->buc_rgb[RED] = std::max(0.0f, std::min(1.0f, src[0]));
    dst->buc_rgb[GRN] = std::max(0.0f, std::min(1.0f, src[1]));
    dst->buc_rgb[BLU] = std::max(0.0f, std::min(1.0f, src[2]));
    dst->buc_rgb[ALP] = 0.0;
}

static SbColor
store_bu_to_sbcolor(const struct bu_color &src)
{
    return SbColor(static_cast<float>(src.buc_rgb[RED]),
	    static_cast<float>(src.buc_rgb[GRN]),
	    static_cast<float>(src.buc_rgb[BLU]));
}

static SoGroup *
store_controller_root_group(BRLObolViewController *controller)
{
    if (!controller)
	return NULL;

    SoNode *root = controller->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return NULL;

    return static_cast<SoGroup *>(root);
}

static void
store_detach_node(BRLObolViewController *controller, SoNode *node)
{
    SoGroup *group = store_controller_root_group(controller);
    if (!group || !node)
	return;

    for (int i = 0; i < group->getNumChildren(); i++) {
	if (group->getChild(i) == node) {
	    group->removeChild(i);
	    return;
	}
    }
}

static void
store_attach_node(BRLObolViewController *controller, SoNode *node)
{
    SoGroup *group = store_controller_root_group(controller);
    if (group && node)
	group->addChild(node);
}

static void
store_release_node(BRLObolViewController *controller, SoNode *node)
{
    if (!node)
	return;

    store_detach_node(controller, node);
    node->unref();
}

static void
store_set_node(BRLObolViewController *controller, SoNode *&slot, SoNode *node)
{
    if (slot == node)
	return;

    if (slot)
	store_release_node(controller, slot);

    slot = node;
    if (slot) {
	slot->ref();
	store_attach_node(controller, slot);
    }
}

BRLObolFeatureHandle::BRLObolFeatureHandle(void) : id(0), revision(0)
{
}

BRLObolFeatureHandle::BRLObolFeatureHandle(uint64_t featureId,
	uint64_t featureRevision) : id(featureId), revision(featureRevision)
{
}

SbBool
BRLObolFeatureHandle::isValid(void) const
{
    return id != 0 && revision != 0 ? TRUE : FALSE;
}

BRLObolPolygonHandle::BRLObolPolygonHandle(void) : id(0), revision(0)
{
}

BRLObolPolygonHandle::BRLObolPolygonHandle(uint64_t polygonId,
	uint64_t polygonRevision) : id(polygonId), revision(polygonRevision)
{
}

SbBool
BRLObolPolygonHandle::isValid(void) const
{
    return id != 0 && revision != 0 ? TRUE : FALSE;
}

SbBool
operator==(const BRLObolFeatureHandle &a, const BRLObolFeatureHandle &b)
{
    return a.id == b.id && a.revision == b.revision ? TRUE : FALSE;
}

SbBool
operator!=(const BRLObolFeatureHandle &a, const BRLObolFeatureHandle &b)
{
    return !(a == b);
}

SbBool
operator==(const BRLObolPolygonHandle &a, const BRLObolPolygonHandle &b)
{
    return a.id == b.id && a.revision == b.revision ? TRUE : FALSE;
}

SbBool
operator!=(const BRLObolPolygonHandle &a, const BRLObolPolygonHandle &b)
{
    return !(a == b);
}

BRLObolFeatureStyle::BRLObolFeatureStyle(void) :
    hasVisible(FALSE),
    visible(TRUE),
    hasColor(FALSE),
    color(1.0f, 1.0f, 1.0f),
    hasLineWidth(FALSE),
    lineWidth(1),
    hasLineStyle(FALSE),
    lineStyle(0),
    hasArrow(FALSE),
    arrow(FALSE),
    hasArrowTip(FALSE),
    arrowTipLength(0.0f),
    arrowTipWidth(0.0f)
{
}

BRLObolCommandResult::BRLObolCommandResult(void) :
    status(BRLObolCommandResultStatus::None),
    feature(),
    polygon(),
    command(""),
    diagnostic("")
{
}

BRLObolFeatureOwner::BRLObolFeatureOwner(void) :
    ownerToken(NULL),
    ownerId(""),
    ownerRole(""),
    resultCallback(NULL),
    callbackUserData(NULL)
{
}

BRLObolOverlayInfo::BRLObolOverlayInfo(void) :
    isOverlay(FALSE),
    ownerToken(NULL),
    role(BRLObolOverlayRole::None),
    overlayClass(BRLObolOverlayClass::None),
    lifecycle(BRLObolOverlayLifecycle::None),
    order(BRLObolOverlayOrder::Model),
    sortOrder(0),
    sourcePath("")
{
}

BRLObolEditPreviewCallbacks::BRLObolEditPreviewCallbacks(void) :
    previewContext(NULL),
    revisionCallback(NULL),
    updateCallback(NULL),
    pickCallback(NULL)
{
}

BRLObolLabel::BRLObolLabel(void) :
    text(""),
    point(0.0f, 0.0f, 0.0f),
    hasColor(FALSE),
    color(1.0f, 1.0f, 1.0f),
    hasLeader(FALSE),
    target(0.0f, 0.0f, 0.0f),
    anchor(0),
    arrow(FALSE),
    fontSize(20.0f)
{
}

BRLObolLineLayer::BRLObolLineLayer(void) :
    name(""),
    points(),
    commands(),
    style()
{
}

BRLObolFeatureSummary::BRLObolFeatureSummary(void) :
    exists(FALSE),
    visible(FALSE),
    realized(FALSE),
    kind(BRLObolFeatureKind::Unknown),
    scope(BRLObolFeatureScope::Shared),
    pointCount(0),
    commandCount(0),
    childCount(0),
    owner(),
    overlay()
{
}

BRLObolFeatureRecord::BRLObolFeatureRecord(void) :
    handle(),
    name(""),
    kind(BRLObolFeatureKind::Unknown),
    scope(BRLObolFeatureScope::Shared),
    style(),
    owner(),
    overlay(),
    realized(FALSE),
    points(),
    commands(),
    indices(),
    normals(),
    labels(),
    axesCenters(),
    halfAxesSize(1.0f),
    layers()
{
}

BRLObolPolygonVisual::BRLObolPolygonVisual(void) :
    edgeColor(1.0f, 1.0f, 0.0f),
    fillColor(0.0f, 0.0f, 1.0f),
    fill(FALSE),
    fillFlags(BRLOBOL_POLYGON_FILL_NONE),
    fillSlope(1.0f, 0.0f),
    fillSpacing(1.0f),
    viewZ(0.0f)
{
}

BRLObolPolygonRecord::BRLObolPolygonRecord(void) :
    handle(),
    name(""),
    scope(BRLObolFeatureScope::Shared),
    type(BRLObolPolygonType::General),
    fill(FALSE),
    fillFlags(BRLOBOL_POLYGON_FILL_NONE),
    fillSlope(1.0f, 0.0f),
    fillSpacing(1.0f),
    fillColor(0.0f, 0.0f, 1.0f),
    edgeColor(1.0f, 1.0f, 0.0f),
    currentContour(-1),
    currentPoint(-1),
    firstContourOpen(FALSE),
    contourCount(0),
    pointCount(0),
    originPoint(0.0f, 0.0f, 0.0f),
    viewZ(0.0f),
    userData(NULL)
{
    HSET(this->viewPlane, 0.0, 0.0, 1.0, 0.0);
}

BRLObolSelectionRecord::BRLObolSelectionRecord(void) :
    path(""),
    feature(),
    owner(),
    kind(BRLOBOL_SELECTION_SELECTED_PATH),
    primitiveIndex(-1),
    hitDistance(0.0)
{
}

struct BRLObolFeatureStoreRecord {
    uint64_t id;
    uint64_t revision;
    SbString name;
    BRLObolFeatureKind kind;
    BRLObolFeatureScope scope;
    BRLObolFeatureStyle style;
    BRLObolFeatureOwner owner;
    BRLObolOverlayInfo overlay;
    BRLObolEditPreviewCallbacks previewCallbacks;
    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    std::vector<int32_t> indices;
    std::vector<SbVec3f> normals;
    std::vector<BRLObolLabel> labels;
    std::vector<SbVec3f> axesCenters;
    float halfAxesSize;
    std::vector<BRLObolLineLayer> layers;
    SbString identity;
    SbString editIntentId;
    SbString editIntentRole;
    uint32_t sourceRevision;
    uint32_t inputsRevision;
    SoNode *node;

    BRLObolFeatureStoreRecord(void) :
	id(0),
	revision(0),
	name(""),
	kind(BRLObolFeatureKind::Unknown),
	scope(BRLObolFeatureScope::Shared),
	style(),
	owner(),
	overlay(),
	previewCallbacks(),
	points(),
	commands(),
	indices(),
	normals(),
	labels(),
	axesCenters(),
	halfAxesSize(1.0f),
	layers(),
	identity(""),
	editIntentId(""),
	editIntentRole(""),
	sourceRevision(0),
	inputsRevision(0),
	node(NULL)
    {
    }
};

struct BRLObolFeatureStore::Impl {
    BRLObolViewController *controller;
    uint64_t nextId;
    std::map<uint64_t, BRLObolFeatureStoreRecord *> records;
    std::map<std::string, uint64_t> names;

    Impl(void) : controller(NULL), nextId(1), records(), names()
    {
    }

    ~Impl(void)
    {
	clear();
    }

    void clear(void)
    {
	for (std::map<uint64_t, BRLObolFeatureStoreRecord *>::iterator it =
		records.begin(); it != records.end(); ++it) {
	    if (it->second) {
		store_release_node(controller, it->second->node);
		delete it->second;
	    }
	}
	records.clear();
	names.clear();
    }

    BRLObolFeatureStoreRecord *record(BRLObolFeatureHandle handle) const
    {
	std::map<uint64_t, BRLObolFeatureStoreRecord *>::const_iterator it =
	    records.find(handle.id);
	if (it == records.end() || !it->second)
	    return NULL;
	return it->second;
    }

    BRLObolFeatureStoreRecord *recordByName(
	    const SbString &name,
	    unsigned int scopeMask,
	    const BRLObolFeatureOwner *owner = NULL) const
    {
	const std::string cleanName = store_string(name);
	if (cleanName.empty())
	    return NULL;

	for (std::map<uint64_t, BRLObolFeatureStoreRecord *>::const_iterator it =
		records.begin(); it != records.end(); ++it) {
	    BRLObolFeatureStoreRecord *rec = it->second;
	    if (!rec)
		continue;
	    if (store_string(rec->name) != cleanName)
		continue;
	    if (!(store_scope_bit(rec->scope) & scopeMask))
		continue;
	    if (rec->scope == BRLObolFeatureScope::Local &&
		    !store_owner_matches(rec->owner, owner))
		continue;
	    return rec;
	}
	return NULL;
    }

    BRLObolFeatureStoreRecord *upsert(const SbString &name,
	    BRLObolFeatureScope scope,
	    BRLObolFeatureKind kind,
	    const BRLObolFeatureStyle *style,
	    const BRLObolFeatureOwner *owner)
    {
	if (store_string(name).empty())
	    return NULL;
	const std::string key = store_key(scope, name, owner);

	uint64_t id = 0;
	std::map<std::string, uint64_t>::iterator nit = names.find(key);
	if (nit != names.end())
	    id = nit->second;

	BRLObolFeatureStoreRecord *rec = NULL;
	if (id) {
	    std::map<uint64_t, BRLObolFeatureStoreRecord *>::iterator rit =
		records.find(id);
	    if (rit != records.end())
		rec = rit->second;
	}

	if (!rec) {
	    rec = new BRLObolFeatureStoreRecord;
	    rec->id = nextId++;
	    rec->revision = 1;
	    rec->name = name;
	    rec->scope = scope;
	    records[rec->id] = rec;
	    names[key] = rec->id;
	} else {
	    rec->revision++;
	}

	rec->kind = kind;
	if (style)
	    rec->style = *style;
	if (owner)
	    rec->owner = *owner;
	return rec;
    }

    BRLObolFeatureHandle handle(const BRLObolFeatureStoreRecord *rec) const
    {
	return rec ? BRLObolFeatureHandle(rec->id, rec->revision) :
	    BRLObolFeatureHandle();
    }

    void notify(const BRLObolFeatureStoreRecord *rec,
	    BRLObolCommandResultStatus status,
	    const char *command,
	    const char *diagnostic = NULL) const
    {
	if (!rec || !rec->owner.resultCallback)
	    return;

	BRLObolCommandResult result;
	result.status = status;
	result.feature = handle(rec);
	result.command = command ? command : "";
	result.diagnostic = diagnostic ? diagnostic : "";
	rec->owner.resultCallback(result, rec->owner.callbackUserData);
    }

    void setNode(BRLObolFeatureStoreRecord *rec, SoNode *node)
    {
	if (!rec)
	    return;
	store_set_node(controller, rec->node, node);
	if (controller && node)
	    controller->requestRender("view-feature-store");
    }
};

static SoBRLVListShape *
store_vlist_node(const BRLObolFeatureStoreRecord &rec)
{
    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->sourcePath = rec.identity.getLength() > 0 ? rec.identity : rec.name;
    shape->sourceName = rec.name;
    switch (rec.kind) {
	case BRLObolFeatureKind::Points:
	    shape->sourceType = "point-set";
	    shape->geometryKind = "point";
	    break;
	case BRLObolFeatureKind::Arrow:
	    shape->sourceType = "arrow";
	    shape->geometryKind = "line";
	    break;
	case BRLObolFeatureKind::IndexedLines:
	    shape->sourceType = "indexed-line-set";
	    shape->geometryKind = "line";
	    break;
	case BRLObolFeatureKind::EditPreview:
	    shape->sourceType = "edit-preview";
	    shape->geometryKind = "line";
	    shape->editEmphasis = TRUE;
	    shape->editIntentId = rec.editIntentId;
	    shape->editIntentRole = rec.editIntentRole;
	    break;
	default:
	    shape->sourceType = "line-set";
	    shape->geometryKind = "line";
	    break;
    }
    shape->displayName = rec.name;
    shape->geometryName = rec.name;
    shape->sourceIdentity = shape->sourcePath.getValue();
    shape->cacheIdentity = shape->sourcePath.getValue();
    shape->databaseIntent = FALSE;
    shape->overlayIntent = TRUE;
    shape->hudIntent = FALSE;
    shape->localSource = rec.scope == BRLObolFeatureScope::Local ? TRUE : FALSE;
    shape->sharedSource = rec.scope == BRLObolFeatureScope::Shared ? TRUE : FALSE;
    shape->nonDatabaseSource = TRUE;
    shape->drawMode = BRLOBOL_LOD_DRAW_DIAGNOSTIC;
    shape->recordRole = "view-feature";
    shape->sourceId = rec.sourceRevision;
    store_apply_vlist_style(shape, rec.style);

    std::vector<int32_t> commands =
	store_normalized_line_commands(rec.points, rec.commands);
    std::vector<int32_t> shapeCommands = store_shape_commands(commands);

    if (!rec.points.empty())
	shape->setLineSet(&rec.points[0], &shapeCommands[0],
		static_cast<int>(rec.points.size()));
    return shape;
}

static SoBRLEditPreview *
store_edit_preview_node(const BRLObolFeatureStoreRecord &rec)
{
    SoBRLEditPreview *preview = new SoBRLEditPreview;
    preview->previewId = rec.name;
    preview->setEditIntent(rec.editIntentId, rec.editIntentRole);
    preview->sourceRevision = rec.sourceRevision;
    preview->inputsRevision = rec.inputsRevision;
    if (!rec.points.empty()) {
	std::vector<int32_t> shapeCommands =
	    store_shape_commands(rec.commands);
	preview->setLineSet(
		rec.identity.getLength() > 0 ? rec.identity : rec.name,
		&rec.points[0],
		&shapeCommands[0],
		static_cast<int>(rec.points.size()));
    }
    return preview;
}

static SoBRLSceneGroup *
store_feature_group_node(const BRLObolFeatureStoreRecord &rec)
{
    SoBRLSceneGroup *group = new SoBRLSceneGroup;
    const SbString path =
	rec.identity.getLength() > 0 ? rec.identity : rec.name;
    group->groupPath = path;
    group->drawIntentValid = TRUE;
    group->drawIntentPath = path;
    group->drawMode = BRLOBOL_LOD_DRAW_DIAGNOSTIC;
    group->fallbackDrawMode = BRLOBOL_LOD_DRAW_WIRE;
    group->overlayIntent = TRUE;
    if (rec.style.hasVisible)
	group->visible = rec.style.visible;
    if (rec.style.hasLineStyle)
	group->lineStyle = rec.style.lineStyle;
    if (rec.style.hasLineWidth)
	group->lineWidth = rec.style.lineWidth;
    if (rec.style.hasColor) {
	group->colorOverride = TRUE;
	group->color = rec.style.color;
    }
    return group;
}

static void
store_append_face_triangles(const std::vector<int32_t> &face,
	size_t pointCount,
	std::vector<int32_t> &triangles)
{
    if (face.size() < 3)
	return;

    std::vector<int32_t> cleanFace = face;
    if (cleanFace.size() > 3 && cleanFace.front() == cleanFace.back())
	cleanFace.pop_back();
    if (cleanFace.size() < 3)
	return;

    for (size_t i = 0; i < cleanFace.size(); i++) {
	if (cleanFace[i] < 0 ||
		static_cast<size_t>(cleanFace[i]) >= pointCount)
	    return;
    }

    for (size_t i = 1; i + 1 < cleanFace.size(); i++) {
	triangles.push_back(cleanFace[0]);
	triangles.push_back(cleanFace[i]);
	triangles.push_back(cleanFace[i + 1]);
    }
}

static std::vector<int32_t>
store_indexed_faces_to_triangles(const std::vector<SbVec3f> &points,
	const std::vector<int32_t> &indices)
{
    std::vector<int32_t> triangles;
    if (points.empty() || indices.empty())
	return triangles;

    const SbBool hasSeparators =
	std::find_if(indices.begin(), indices.end(),
		[](int32_t idx) { return idx < 0; }) != indices.end() ?
	TRUE : FALSE;

    if (!hasSeparators && indices.size() % 3 == 0) {
	for (size_t i = 0; i < indices.size(); i += 3) {
	    std::vector<int32_t> face;
	    face.push_back(indices[i]);
	    face.push_back(indices[i + 1]);
	    face.push_back(indices[i + 2]);
	    store_append_face_triangles(face, points.size(), triangles);
	}
	return triangles;
    }

    std::vector<int32_t> face;
    for (size_t i = 0; i < indices.size(); i++) {
	if (indices[i] < 0) {
	    store_append_face_triangles(face, points.size(), triangles);
	    face.clear();
	} else {
	    face.push_back(indices[i]);
	}
    }
    store_append_face_triangles(face, points.size(), triangles);
    return triangles;
}

static SoNode *
store_indexed_face_node(const BRLObolFeatureStoreRecord &rec)
{
    SoBRLMeshShape *shape = new SoBRLMeshShape;
    shape->sourcePath = rec.identity.getLength() > 0 ? rec.identity : rec.name;
    shape->sourceName = rec.name;
    shape->sourceType = "indexed-face-set";
    shape->displayName = rec.name;
    shape->geometryName = rec.name;
    shape->sourceIdentity = shape->sourcePath.getValue();
    shape->cacheIdentity = shape->sourcePath.getValue();
    shape->databaseIntent = FALSE;
    shape->overlayIntent = TRUE;
    shape->hudIntent = FALSE;
    shape->localSource = rec.scope == BRLObolFeatureScope::Local ? TRUE : FALSE;
    shape->sharedSource = rec.scope == BRLObolFeatureScope::Shared ? TRUE : FALSE;
    shape->nonDatabaseSource = TRUE;
    shape->drawMode = BRLOBOL_LOD_DRAW_SHADED;
    shape->recordRole = "view-feature";
    shape->geometryKind = "surface";
    shape->sourceId = rec.sourceRevision;
    store_apply_mesh_style(shape, rec.style);

    std::vector<int32_t> triangles =
	store_indexed_faces_to_triangles(rec.points, rec.indices);
    if (!rec.points.empty() && !triangles.empty())
	shape->setIndexedTriangles(&rec.points[0],
		static_cast<int>(rec.points.size()),
		&triangles[0],
		static_cast<int>(triangles.size()));
    return shape;
}

static SoNode *
store_line_layers_node(const BRLObolFeatureStoreRecord &rec)
{
    SoBRLSceneGroup *sep = store_feature_group_node(rec);
    for (size_t i = 0; i < rec.layers.size(); i++) {
	BRLObolFeatureStoreRecord layerRec = rec;
	layerRec.name = rec.layers[i].name.getLength() > 0 ?
	    rec.layers[i].name : rec.name;
	layerRec.identity = layerRec.name;
	layerRec.points = rec.layers[i].points;
	layerRec.commands = rec.layers[i].commands;
	layerRec.style = rec.layers[i].style;
	layerRec.kind = BRLObolFeatureKind::Lines;
	SoBRLVListShape *shape = store_vlist_node(layerRec);
	if (shape)
	    sep->addChild(shape);
    }
    return sep;
}

static SoNode *
store_axes_node(const BRLObolFeatureStoreRecord &rec)
{
    SoBRLSceneGroup *sep = store_feature_group_node(rec);
    const float size = rec.halfAxesSize > 0.0f ? rec.halfAxesSize : 1.0f;

    for (size_t i = 0; i < rec.axesCenters.size(); i++) {
	SoBRLAxes *axes = new SoBRLAxes;
	axes->overlayId = rec.name;
	axes->origin = rec.axesCenters[i];
	axes->size = size;
	if (rec.style.hasVisible)
	    axes->visible = rec.style.visible;
	SoBRLVListShape *shape = axes->rebuildGeometry();
	store_apply_vlist_style(shape, rec.style);
	sep->addChild(axes);
    }

    return sep;
}

static SoNode *
store_label_node(const BRLObolFeatureStoreRecord &rec)
{
    SoBRLSceneGroup *sep = store_feature_group_node(rec);
    const SbColor fallbackColor = rec.style.hasColor ?
	rec.style.color : SbColor(1.0f, 1.0f, 1.0f);

    for (size_t i = 0; i < rec.labels.size(); i++) {
	const BRLObolLabel &label = rec.labels[i];
	const SbColor color = label.hasColor ? label.color : fallbackColor;

	if (label.hasLeader) {
	    BRLObolFeatureStoreRecord leader = rec;
	    leader.kind = BRLObolFeatureKind::Lines;
	    leader.points.clear();
	    leader.commands.clear();
	    leader.points.push_back(label.target);
	    leader.points.push_back(label.point);
	    leader.commands.push_back(static_cast<int32_t>(
		    BRLObolLineCommand::Move));
	    leader.commands.push_back(static_cast<int32_t>(
		    BRLObolLineCommand::Draw));
	    leader.style.hasColor = TRUE;
	    leader.style.color = color;
	    SoBRLVListShape *shape = store_vlist_node(leader);
	    shape->sourceType = "label-leader";
	    shape->geometryKind = "line";
	    sep->addChild(shape);
	}

	if (label.text.getLength() == 0)
	    continue;

	SoSeparator *textSep = new SoSeparator;
	SoTranslation *translation = new SoTranslation;
	translation->translation = label.point;
	textSep->addChild(translation);

	SoBaseColor *baseColor = new SoBaseColor;
	baseColor->rgb = color;
	textSep->addChild(baseColor);

	SoFont *font = new SoFont;
	font->size = label.fontSize > 0.0f ? label.fontSize : 20.0f;
	textSep->addChild(font);

	SoText2 *text = new SoText2;
	text->string.set1Value(0, label.text);
	text->justification = label.anchor == 2 ? SoText2::RIGHT :
	    label.anchor == 1 ? SoText2::CENTER : SoText2::LEFT;
	text->depthTest = FALSE;
	textSep->addChild(text);

	sep->addChild(textSep);
    }
    return sep;
}

static SoNode *
store_node_for_feature(const BRLObolFeatureStoreRecord &rec)
{
    switch (rec.kind) {
	case BRLObolFeatureKind::Labels:
	case BRLObolFeatureKind::HudLabel:
	    return store_label_node(rec);
	case BRLObolFeatureKind::Axes:
	    return store_axes_node(rec);
	case BRLObolFeatureKind::LineLayer:
	    return store_line_layers_node(rec);
	case BRLObolFeatureKind::IndexedFaceSet:
	    return store_indexed_face_node(rec);
	case BRLObolFeatureKind::EditPreview:
	    return store_edit_preview_node(rec);
	default:
	    return store_vlist_node(rec);
    }
}

BRLObolFeatureStore::BRLObolFeatureStore(void) : impl(new Impl)
{
}

BRLObolFeatureStore::BRLObolFeatureStore(BRLObolViewController *controller) :
    impl(new Impl)
{
    this->impl->controller = controller;
}

BRLObolFeatureStore::~BRLObolFeatureStore(void)
{
    delete this->impl;
    this->impl = NULL;
}

void
BRLObolFeatureStore::setController(BRLObolViewController *controller)
{
    this->impl->controller = controller;
}

BRLObolViewController *
BRLObolFeatureStore::controller(void) const
{
    return this->impl->controller;
}

void
BRLObolFeatureStore::clear(void)
{
    this->impl->clear();
}

BRLObolFeatureHandle
BRLObolFeatureStore::find(const SbString &name, unsigned int scopeMask) const
{
    return this->impl->handle(this->impl->recordByName(name, scopeMask));
}

BRLObolFeatureHandle
BRLObolFeatureStore::findOwned(const SbString &name,
	unsigned int scopeMask,
	const BRLObolFeatureOwner *owner) const
{
    return this->impl->handle(this->impl->recordByName(name, scopeMask,
	    owner));
}

SbBool
BRLObolFeatureStore::exists(const SbString &name, unsigned int scopeMask) const
{
    return this->find(name, scopeMask).isValid();
}

SbBool
BRLObolFeatureStore::existsOwned(const SbString &name,
	unsigned int scopeMask,
	const BRLObolFeatureOwner *owner) const
{
    return this->findOwned(name, scopeMask, owner).isValid();
}

SbBool
BRLObolFeatureStore::remove(BRLObolFeatureHandle handle)
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    const std::string key = store_key(rec->scope, rec->name, &rec->owner);
    this->impl->names.erase(key);
    this->impl->notify(rec, BRLObolCommandResultStatus::Removed, "remove");
    store_release_node(this->impl->controller, rec->node);
    this->impl->records.erase(rec->id);
    delete rec;
    if (this->impl->controller)
	this->impl->controller->requestRender("view-feature-remove");
    return TRUE;
}

SbBool
BRLObolFeatureStore::remove(const SbString &name)
{
    return this->remove(this->find(name));
}

SbBool
BRLObolFeatureStore::removeOwned(const SbString &name,
	unsigned int scopeMask,
	const BRLObolFeatureOwner *owner)
{
    return this->remove(this->findOwned(name, scopeMask, owner));
}

size_t
BRLObolFeatureStore::removePrefix(const SbString &prefix)
{
    return this->removePrefix(prefix, BRLOBOL_FEATURE_SCOPE_ALL, NULL);
}

size_t
BRLObolFeatureStore::removePrefix(const SbString &prefix,
	unsigned int scopeMask,
	const BRLObolFeatureOwner *owner)
{
    const std::string p = store_string(prefix);
    if (p.empty())
	return 0;

    std::vector<BRLObolFeatureHandle> handles;
    for (std::map<uint64_t, BRLObolFeatureStoreRecord *>::const_iterator it =
	    this->impl->records.begin(); it != this->impl->records.end(); ++it) {
	if (!it->second)
	    continue;
	if (!(store_scope_bit(it->second->scope) & scopeMask))
	    continue;
	if (it->second->scope == BRLObolFeatureScope::Local &&
		!store_owner_matches(it->second->owner, owner))
	    continue;
	const std::string name = store_string(it->second->name);
	if (name.compare(0, p.size(), p) == 0)
	    handles.push_back(this->impl->handle(it->second));
    }

    size_t removed = 0;
    for (size_t i = 0; i < handles.size(); i++)
	if (this->remove(handles[i]))
	    removed++;
    return removed;
}

BRLObolFeatureHandle
BRLObolFeatureStore::publishLineSet(const SbString &name,
	BRLObolFeatureScope scope,
	const std::vector<SbVec3f> &points,
	const std::vector<int32_t> &commands,
	const BRLObolFeatureStyle *style,
	const BRLObolFeatureOwner *owner)
{
    BRLObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
	    BRLObolFeatureKind::Lines, style, owner);
    if (!rec)
	return BRLObolFeatureHandle();

    rec->points = points;
    rec->commands = commands;
    SoNode *node = store_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BRLObolCommandResultStatus::Updated, "publishLineSet");
    return this->impl->handle(rec);
}

BRLObolFeatureHandle
BRLObolFeatureStore::publishIndexedLineSet(const SbString &name,
	BRLObolFeatureScope scope,
	const std::vector<SbVec3f> &points,
	const std::vector<int32_t> &indices,
	const BRLObolFeatureStyle *style,
	const BRLObolFeatureOwner *owner)
{
    std::vector<SbVec3f> linePoints;
    std::vector<int32_t> commands;
    linePoints.reserve(indices.size());
    commands.reserve(indices.size());
    for (size_t i = 0; i < indices.size(); i++) {
	int32_t idx = indices[i];
	if (idx < 0 || static_cast<size_t>(idx) >= points.size())
	    continue;
	linePoints.push_back(points[static_cast<size_t>(idx)]);
	commands.push_back(linePoints.size() % 2 == 1 ?
		static_cast<int32_t>(BRLObolLineCommand::Move) :
		static_cast<int32_t>(BRLObolLineCommand::Draw));
    }

    BRLObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
	    BRLObolFeatureKind::IndexedLines, style, owner);
    if (!rec)
	return BRLObolFeatureHandle();

    rec->points = linePoints;
    rec->commands = commands;
    rec->indices = indices;
    SoNode *node = store_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BRLObolCommandResultStatus::Updated,
	    "publishIndexedLineSet");
    return this->impl->handle(rec);
}

BRLObolFeatureHandle
BRLObolFeatureStore::publishPointSet(const SbString &name,
	BRLObolFeatureScope scope,
	const std::vector<SbVec3f> &points,
	const BRLObolFeatureStyle *style,
	const BRLObolFeatureOwner *owner)
{
    std::vector<int32_t> commands(points.size(),
	    static_cast<int32_t>(BRLObolLineCommand::Point));
    BRLObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
	    BRLObolFeatureKind::Points, style, owner);
    if (!rec)
	return BRLObolFeatureHandle();
    rec->points = points;
    rec->commands = commands;
    SoNode *node = store_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BRLObolCommandResultStatus::Updated,
	    "publishPointSet");
    return this->impl->handle(rec);
}

BRLObolFeatureHandle
BRLObolFeatureStore::publishLabels(const SbString &name,
	BRLObolFeatureScope scope,
	const std::vector<BRLObolLabel> &labels,
	const BRLObolFeatureStyle *style,
	const BRLObolFeatureOwner *owner)
{
    BRLObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
	    BRLObolFeatureKind::Labels, style, owner);
    if (!rec)
	return BRLObolFeatureHandle();
    rec->labels = labels;
    SoNode *node = store_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BRLObolCommandResultStatus::Updated,
	    "publishLabels");
    return this->impl->handle(rec);
}

BRLObolFeatureHandle
BRLObolFeatureStore::publishArrow(const SbString &name,
	BRLObolFeatureScope scope,
	const std::vector<SbVec3f> &points,
	const BRLObolFeatureStyle *style,
	const BRLObolFeatureOwner *owner)
{
    BRLObolFeatureStyle arrowStyle = style ? *style : BRLObolFeatureStyle();
    arrowStyle.hasArrow = TRUE;
    arrowStyle.arrow = TRUE;

    BRLObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
	    BRLObolFeatureKind::Arrow, &arrowStyle, owner);
    if (!rec)
	return BRLObolFeatureHandle();

    rec->points = points;
    rec->commands.clear();
    SoNode *node = store_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BRLObolCommandResultStatus::Updated,
	    "publishArrow");
    return this->impl->handle(rec);
}

BRLObolFeatureHandle
BRLObolFeatureStore::publishAxes(const SbString &name,
	BRLObolFeatureScope scope,
	const std::vector<SbVec3f> &centers,
	float halfAxesSize,
	const BRLObolFeatureStyle *style,
	const BRLObolFeatureOwner *owner)
{
    BRLObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
	    BRLObolFeatureKind::Axes, style, owner);
    if (!rec)
	return BRLObolFeatureHandle();
    rec->axesCenters = centers;
    rec->halfAxesSize = halfAxesSize;
    SoNode *node = store_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BRLObolCommandResultStatus::Updated,
	    "publishAxes");
    return this->impl->handle(rec);
}

BRLObolFeatureHandle
BRLObolFeatureStore::publishLineLayers(const SbString &name,
	BRLObolFeatureScope scope,
	const std::vector<BRLObolLineLayer> &layers,
	const BRLObolFeatureStyle *style,
	const BRLObolFeatureOwner *owner)
{
    BRLObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
	    BRLObolFeatureKind::LineLayer, style, owner);
    if (!rec)
	return BRLObolFeatureHandle();
    rec->layers = layers;
    rec->points.clear();
    rec->commands.clear();
    for (size_t i = 0; i < layers.size(); i++) {
	rec->points.insert(rec->points.end(), layers[i].points.begin(),
		layers[i].points.end());
	rec->commands.insert(rec->commands.end(), layers[i].commands.begin(),
		layers[i].commands.end());
    }
    SoNode *node = store_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BRLObolCommandResultStatus::Updated,
	    "publishLineLayers");
    return this->impl->handle(rec);
}

BRLObolFeatureHandle
BRLObolFeatureStore::publishLineLayerBuilder(const SbString &name,
	BRLObolFeatureScope scope,
	const struct bg_line_layer_builder *builder,
	const BRLObolFeatureStyle *style,
	const BRLObolFeatureOwner *owner)
{
    if (!builder)
	return BRLObolFeatureHandle();

    BRLObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
	    BRLObolFeatureKind::LineLayer, style, owner);
    if (!rec)
	return BRLObolFeatureHandle();

    SoBRLLineLayerOverlay *overlay = new SoBRLLineLayerOverlay;
    overlay->overlayId = name;
    overlay->sourceId = static_cast<uint32_t>(rec->revision);
    overlay->selectable = TRUE;
    overlay->rebuildGeometry(builder);
    this->impl->setNode(rec, overlay);
    this->impl->notify(rec, BRLObolCommandResultStatus::Updated,
	    "publishLineLayerBuilder");
    return this->impl->handle(rec);
}

BRLObolFeatureHandle
BRLObolFeatureStore::publishIndexedFaceSet(const SbString &name,
	BRLObolFeatureScope scope,
	const std::vector<SbVec3f> &points,
	const std::vector<SbVec3f> &normals,
	const std::vector<int32_t> &indices,
	const BRLObolFeatureStyle *style,
	const BRLObolFeatureOwner *owner)
{
    BRLObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
	    BRLObolFeatureKind::IndexedFaceSet, style, owner);
    if (!rec)
	return BRLObolFeatureHandle();

    rec->points = points;
    rec->normals = normals;
    rec->indices = indices;
    rec->commands.clear();
    SoNode *node = store_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BRLObolCommandResultStatus::Updated,
	    "publishIndexedFaceSet");
    return this->impl->handle(rec);
}

BRLObolFeatureHandle
BRLObolFeatureStore::publishEditPreview(const SbString &name,
	const SbString &identity,
	const SbString &editIntentId,
	const SbString &editIntentRole,
	const std::vector<SbVec3f> &points,
	const std::vector<int32_t> &commands,
	uint32_t sourceRevision,
	uint32_t inputsRevision,
	const BRLObolEditPreviewCallbacks *callbacks,
	const BRLObolFeatureOwner *owner)
{
    BRLObolFeatureStoreRecord *rec = this->impl->upsert(name,
	    BRLObolFeatureScope::Local, BRLObolFeatureKind::EditPreview,
	    NULL, owner);
    if (!rec)
	return BRLObolFeatureHandle();

    rec->identity = identity;
    rec->editIntentId = editIntentId.getLength() > 0 ? editIntentId : name;
    rec->editIntentRole = editIntentRole.getLength() > 0 ?
	editIntentRole : SbString("preview");
    rec->points = points;
    rec->commands = store_normalized_line_commands(points, commands);
    rec->sourceRevision = sourceRevision ? sourceRevision :
	static_cast<uint32_t>(rec->revision);
    rec->inputsRevision = inputsRevision ? inputsRevision :
	static_cast<uint32_t>(rec->revision);
    if (callbacks)
	rec->previewCallbacks = *callbacks;

    SoNode *node = store_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BRLObolCommandResultStatus::Updated,
	    "publishEditPreview");
    return this->impl->handle(rec);
}

SbBool
BRLObolFeatureStore::replaceEditPreviewGeometry(
	BRLObolFeatureHandle handle,
	const SbString &identity,
	const std::vector<SbVec3f> &points,
	const std::vector<int32_t> &commands,
	uint32_t sourceRevision,
	uint32_t inputsRevision)
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec || rec->kind != BRLObolFeatureKind::EditPreview)
	return FALSE;

    if (identity.getLength() > 0)
	rec->identity = identity;
    rec->points = points;
    rec->commands = store_normalized_line_commands(points, commands);
    rec->revision++;
    rec->sourceRevision = sourceRevision ? sourceRevision :
	static_cast<uint32_t>(rec->revision);
    rec->inputsRevision = inputsRevision ? inputsRevision :
	static_cast<uint32_t>(rec->revision);

    SoNode *node = store_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BRLObolCommandResultStatus::Updated,
	    "replaceEditPreviewGeometry");
    return TRUE;
}

uint64_t
BRLObolFeatureStore::editPreviewRevision(BRLObolFeatureHandle handle) const
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec || rec->kind != BRLObolFeatureKind::EditPreview)
	return 0;
    if (rec->previewCallbacks.revisionCallback)
	return rec->previewCallbacks.revisionCallback(
		rec->previewCallbacks.previewContext);
    return rec->sourceRevision;
}

int
BRLObolFeatureStore::updateEditPreview(BRLObolFeatureHandle handle)
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec || rec->kind != BRLObolFeatureKind::EditPreview ||
	    !rec->previewCallbacks.updateCallback)
	return -1;
    return rec->previewCallbacks.updateCallback(
	    rec->previewCallbacks.previewContext);
}

int
BRLObolFeatureStore::pickEditPreview(BRLObolFeatureHandle handle,
	int x,
	int y,
	void *pickOut) const
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec || rec->kind != BRLObolFeatureKind::EditPreview ||
	    !rec->previewCallbacks.pickCallback)
	return -1;
    return rec->previewCallbacks.pickCallback(
	    rec->previewCallbacks.previewContext, x, y, pickOut);
}

SbBool
BRLObolFeatureStore::appendLinePoint(BRLObolFeatureHandle handle,
	const SbVec3f &point)
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    rec->points.push_back(point);
    rec->commands.push_back(static_cast<int32_t>(BRLObolLineCommand::Draw));
    rec->revision++;
    SoNode *node = store_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BRLObolCommandResultStatus::Updated,
	    "appendLinePoint");
    return TRUE;
}

SbBool
BRLObolFeatureStore::replaceLabels(BRLObolFeatureHandle handle,
	const std::vector<BRLObolLabel> &labels)
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    rec->labels = labels;
    rec->revision++;
    SoNode *node = store_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BRLObolCommandResultStatus::Updated,
	    "replaceLabels");
    return TRUE;
}

SbBool
BRLObolFeatureStore::clearGeometry(BRLObolFeatureHandle handle)
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    rec->points.clear();
    rec->commands.clear();
    rec->indices.clear();
    rec->normals.clear();
    rec->labels.clear();
    rec->axesCenters.clear();
    rec->layers.clear();
    rec->revision++;
    this->impl->setNode(rec, NULL);
    this->impl->notify(rec, BRLObolCommandResultStatus::Updated,
	    "clearGeometry");
    return TRUE;
}

SbBool
BRLObolFeatureStore::points(BRLObolFeatureHandle handle,
	std::vector<SbVec3f> &pointsOut) const
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    pointsOut = rec->points;
    return TRUE;
}

SbBool
BRLObolFeatureStore::commands(BRLObolFeatureHandle handle,
	std::vector<int32_t> &commandsOut) const
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    commandsOut = rec->commands;
    return TRUE;
}

SbBool
BRLObolFeatureStore::lineCommandAt(BRLObolFeatureHandle handle,
	size_t index,
	int32_t &commandOut) const
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec || index >= rec->commands.size())
	return FALSE;
    commandOut = rec->commands[index];
    return TRUE;
}

SbBool
BRLObolFeatureStore::labels(BRLObolFeatureHandle handle,
	std::vector<BRLObolLabel> &labelsOut) const
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    labelsOut = rec->labels;
    return TRUE;
}

SbBool
BRLObolFeatureStore::axesCenters(BRLObolFeatureHandle handle,
	std::vector<SbVec3f> &centersOut,
	float *halfAxesSizeOut) const
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    centersOut = rec->axesCenters;
    if (halfAxesSizeOut)
	*halfAxesSizeOut = rec->halfAxesSize;
    return TRUE;
}

SbBool
BRLObolFeatureStore::indices(BRLObolFeatureHandle handle,
	std::vector<int32_t> &indicesOut) const
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    indicesOut = rec->indices;
    return TRUE;
}

SbBool
BRLObolFeatureStore::normals(BRLObolFeatureHandle handle,
	std::vector<SbVec3f> &normalsOut) const
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    normalsOut = rec->normals;
    return TRUE;
}

SbBool
BRLObolFeatureStore::applyStyle(BRLObolFeatureHandle handle,
	const BRLObolFeatureStyle &style,
	SbBool UNUSED(recursive))
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    if (style.hasVisible) {
	rec->style.hasVisible = TRUE;
	rec->style.visible = style.visible;
    }
    if (style.hasColor) {
	rec->style.hasColor = TRUE;
	rec->style.color = style.color;
    }
    if (style.hasLineWidth) {
	rec->style.hasLineWidth = TRUE;
	rec->style.lineWidth = style.lineWidth;
    }
    if (style.hasLineStyle) {
	rec->style.hasLineStyle = TRUE;
	rec->style.lineStyle = style.lineStyle;
    }
    if (style.hasArrow) {
	rec->style.hasArrow = TRUE;
	rec->style.arrow = style.arrow;
    }
    if (style.hasArrowTip) {
	rec->style.hasArrowTip = TRUE;
	rec->style.arrowTipLength = style.arrowTipLength;
	rec->style.arrowTipWidth = style.arrowTipWidth;
    }

    rec->revision++;
    SoNode *node = store_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BRLObolCommandResultStatus::Updated,
	    "applyStyle");
    return TRUE;
}

SbBool
BRLObolFeatureStore::style(BRLObolFeatureHandle handle,
	BRLObolFeatureStyle &styleOut) const
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    styleOut = rec->style;
    return TRUE;
}

SbBool
BRLObolFeatureStore::setVisible(BRLObolFeatureHandle handle, SbBool visible)
{
    BRLObolFeatureStyle style;
    style.hasVisible = TRUE;
    style.visible = visible;
    return this->applyStyle(handle, style);
}

SbBool
BRLObolFeatureStore::setColor(BRLObolFeatureHandle handle,
	const SbColor &color)
{
    BRLObolFeatureStyle style;
    style.hasColor = TRUE;
    style.color = color;
    return this->applyStyle(handle, style);
}

SbBool
BRLObolFeatureStore::setLineWidth(BRLObolFeatureHandle handle,
	int lineWidth)
{
    BRLObolFeatureStyle style;
    style.hasLineWidth = TRUE;
    style.lineWidth = lineWidth;
    return this->applyStyle(handle, style);
}

SbBool
BRLObolFeatureStore::arrowTip(BRLObolFeatureHandle handle,
	float &tipLength,
	float &tipWidth) const
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    tipLength = rec->style.arrowTipLength;
    tipWidth = rec->style.arrowTipWidth;
    return TRUE;
}

SbBool
BRLObolFeatureStore::setArrowTip(BRLObolFeatureHandle handle,
	float tipLength,
	float tipWidth)
{
    BRLObolFeatureStyle style;
    style.hasArrowTip = TRUE;
    style.arrowTipLength = tipLength;
    style.arrowTipWidth = tipWidth;
    return this->applyStyle(handle, style);
}

SbBool
BRLObolFeatureStore::setOverlayInfo(BRLObolFeatureHandle handle,
	const BRLObolOverlayInfo &overlay)
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    rec->overlay = overlay;
    rec->revision++;
    SoNode *node = store_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BRLObolCommandResultStatus::Updated,
	    "setOverlayInfo");
    return TRUE;
}

SbBool
BRLObolFeatureStore::clearOverlayInfo(BRLObolFeatureHandle handle)
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    rec->overlay = BRLObolOverlayInfo();
    rec->revision++;
    SoNode *node = store_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BRLObolCommandResultStatus::Updated,
	    "clearOverlayInfo");
    return TRUE;
}

SbBool
BRLObolFeatureStore::overlayInfo(BRLObolFeatureHandle handle,
	BRLObolOverlayInfo &overlayOut) const
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    overlayOut = rec->overlay;
    return TRUE;
}

SbBool
BRLObolFeatureStore::realize(BRLObolFeatureHandle handle,
	SbBool UNUSED(recursive))
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    SoNode *node = store_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    return TRUE;
}

SbBool
BRLObolFeatureStore::summary(const SbString &name,
	BRLObolFeatureSummary &summaryOut,
	unsigned int scopeMask) const
{
    return this->summaryOwned(name, summaryOut, scopeMask, NULL);
}

SbBool
BRLObolFeatureStore::summaryOwned(const SbString &name,
	BRLObolFeatureSummary &summaryOut,
	unsigned int scopeMask,
	const BRLObolFeatureOwner *owner) const
{
    summaryOut = BRLObolFeatureSummary();
    BRLObolFeatureStoreRecord *rec = this->impl->recordByName(name,
	    scopeMask, owner);
    if (!rec)
	return TRUE;

    summaryOut.exists = TRUE;
    summaryOut.visible = rec->style.hasVisible ? rec->style.visible : TRUE;
    summaryOut.realized = rec->node ? TRUE : FALSE;
    summaryOut.kind = rec->kind;
    summaryOut.scope = rec->scope;
    summaryOut.pointCount = rec->points.size();
    summaryOut.commandCount = rec->commands.size();
    summaryOut.owner = rec->owner;
    summaryOut.overlay = rec->overlay;
    if (rec->node && rec->node->isOfType(SoGroup::getClassTypeId()))
	summaryOut.childCount =
	    static_cast<size_t>(static_cast<SoGroup *>(rec->node)->getNumChildren());
    else
	summaryOut.childCount = rec->node ? 1 : 0;
    return TRUE;
}

SbBool
BRLObolFeatureStore::record(BRLObolFeatureHandle handle,
	BRLObolFeatureRecord &recordOut) const
{
    recordOut = BRLObolFeatureRecord();
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    recordOut.handle = this->impl->handle(rec);
    recordOut.name = rec->name;
    recordOut.kind = rec->kind;
    recordOut.scope = rec->scope;
    recordOut.style = rec->style;
    recordOut.owner = rec->owner;
    recordOut.overlay = rec->overlay;
    recordOut.realized = rec->node ? TRUE : FALSE;
    recordOut.points = rec->points;
    recordOut.commands = rec->commands;
    recordOut.indices = rec->indices;
    recordOut.normals = rec->normals;
    recordOut.labels = rec->labels;
    recordOut.axesCenters = rec->axesCenters;
    recordOut.halfAxesSize = rec->halfAxesSize;
    recordOut.layers = rec->layers;
    return TRUE;
}

void
BRLObolFeatureStore::visitRecords(BRLObolFeatureRecordCallback callback,
	void *userData,
	unsigned int scopeMask,
	const BRLObolFeatureOwner *owner) const
{
    if (!callback)
	return;

    for (std::map<uint64_t, BRLObolFeatureStoreRecord *>::const_iterator it =
	    this->impl->records.begin(); it != this->impl->records.end(); ++it) {
	BRLObolFeatureStoreRecord *rec = it->second;
	if (!rec)
	    continue;
	if (!(store_scope_bit(rec->scope) & scopeMask))
	    continue;
	if (rec->scope == BRLObolFeatureScope::Local &&
		!store_owner_matches(rec->owner, owner))
	    continue;

	BRLObolFeatureRecord record;
	if (!this->record(this->impl->handle(rec), record))
	    continue;
	if (!callback(record, userData))
	    return;
    }
}

SoNode *
BRLObolFeatureStore::node(BRLObolFeatureHandle handle) const
{
    BRLObolFeatureStoreRecord *rec = this->impl->record(handle);
    return rec ? rec->node : NULL;
}

struct BRLObolPolygonStoreRecord {
    uint64_t id;
    uint64_t revision;
    SbString name;
    BRLObolFeatureScope scope;
    BRLObolPolygonType type;
    BRLObolPolygonVisual visual;
    long currentContour;
    long currentPoint;
    SbVec3f originPoint;
    plane_t viewPlane;
    void *userData;
    struct bg_polygon polygon;
    SoNode *node;

    BRLObolPolygonStoreRecord(void) :
	id(0),
	revision(0),
	name(""),
	scope(BRLObolFeatureScope::Shared),
	type(BRLObolPolygonType::General),
	visual(),
	currentContour(-1),
	currentPoint(-1),
	originPoint(0.0f, 0.0f, 0.0f),
	userData(NULL),
	polygon(),
	node(NULL)
    {
	HSET(viewPlane, 0.0, 0.0, 1.0, 0.0);
    }
};

static SoNode *
store_polygon_node(const BRLObolPolygonStoreRecord &rec);

struct BRLObolPolygonStore::Impl {
    BRLObolViewController *controller;
    uint64_t nextId;
    std::map<uint64_t, BRLObolPolygonStoreRecord *> records;
    std::map<std::string, uint64_t> names;
    BRLObolPolygonHandle snapExclude;

    Impl(void) : controller(NULL), nextId(1), records(), names(), snapExclude()
    {
    }

    ~Impl(void)
    {
	clear();
    }

    void clear(void)
    {
	for (std::map<uint64_t, BRLObolPolygonStoreRecord *>::iterator it =
		records.begin(); it != records.end(); ++it) {
	    if (it->second) {
		store_release_node(controller, it->second->node);
		bg_polygon_free(&it->second->polygon);
		delete it->second;
	    }
	}
	records.clear();
	names.clear();
	snapExclude = BRLObolPolygonHandle();
    }

    BRLObolPolygonStoreRecord *record(BRLObolPolygonHandle handle) const
    {
	std::map<uint64_t, BRLObolPolygonStoreRecord *>::const_iterator it =
	    records.find(handle.id);
	if (it == records.end() || !it->second)
	    return NULL;
	return it->second;
    }

    BRLObolPolygonStoreRecord *recordByName(const SbString &name, unsigned int scopeMask) const
    {
	const std::string cleanName = store_string(name);
	if (cleanName.empty())
	    return NULL;
	for (std::map<std::string, uint64_t>::const_iterator it =
		names.begin(); it != names.end(); ++it) {
	    const std::string &key = it->first;
	    if (key.size() < 3 || key.substr(2) != cleanName)
		continue;
	    std::map<uint64_t, BRLObolPolygonStoreRecord *>::const_iterator rit =
		records.find(it->second);
	    if (rit == records.end() || !rit->second)
		continue;
	    if (!(store_scope_bit(rit->second->scope) & scopeMask))
		continue;
	    return rit->second;
	}
	return NULL;
    }

    BRLObolPolygonHandle handle(const BRLObolPolygonStoreRecord *rec) const
    {
	return rec ? BRLObolPolygonHandle(rec->id, rec->revision) :
	    BRLObolPolygonHandle();
    }

    void setNode(BRLObolPolygonStoreRecord *rec, SoNode *node)
    {
	if (!rec)
	    return;
	store_set_node(controller, rec->node, node);
	if (controller && node)
	    controller->requestRender("view-polygon-store");
    }

    void realize(BRLObolPolygonStoreRecord *rec)
    {
	if (!rec)
	    return;
	SoNode *node = store_polygon_node(*rec);
	setNode(rec, node);
    }
};

static void
store_polygon_init_one_point(struct bg_polygon *poly, const SbVec3f &point)
{
    if (!poly)
	return;
    bg_polygon_free(poly);
    poly->num_contours = 1;
    poly->hole = (int *)bu_calloc(1, sizeof(int), "BRLObol polygon hole");
    poly->contour = (struct bg_poly_contour *)bu_calloc(1,
	    sizeof(struct bg_poly_contour), "BRLObol polygon contour");
    poly->contour[0].num_points = 1;
    poly->contour[0].open = 1;
    poly->contour[0].point = (point_t *)bu_calloc(1, sizeof(point_t),
	    "BRLObol polygon point");
    store_point(poly->contour[0].point[0], point);
}

static size_t
store_polygon_point_count(const struct bg_polygon &poly)
{
    size_t count = 0;
    for (size_t i = 0; i < poly.num_contours; i++)
	count += poly.contour[i].num_points;
    return count;
}

static int
store_polygon_type_to_rt(BRLObolPolygonType type)
{
    switch (type) {
	case BRLObolPolygonType::Circle:
	    return RT_SKETCH_POLYGON_CIRCLE;
	case BRLObolPolygonType::Ellipse:
	    return RT_SKETCH_POLYGON_ELLIPSE;
	case BRLObolPolygonType::Rectangle:
	    return RT_SKETCH_POLYGON_RECTANGLE;
	case BRLObolPolygonType::Square:
	    return RT_SKETCH_POLYGON_SQUARE;
	default:
	    return RT_SKETCH_POLYGON_GENERAL;
    }
}

static BRLObolPolygonType
store_polygon_type_from_rt(int type)
{
    switch (type) {
	case RT_SKETCH_POLYGON_CIRCLE:
	    return BRLObolPolygonType::Circle;
	case RT_SKETCH_POLYGON_ELLIPSE:
	    return BRLObolPolygonType::Ellipse;
	case RT_SKETCH_POLYGON_RECTANGLE:
	    return BRLObolPolygonType::Rectangle;
	case RT_SKETCH_POLYGON_SQUARE:
	    return BRLObolPolygonType::Square;
	default:
	    return BRLObolPolygonType::General;
    }
}

static unsigned int
store_polygon_valid_fill_flags(unsigned int flags)
{
    return flags & (BRLOBOL_POLYGON_FILL_HATCH | BRLOBOL_POLYGON_FILL_MESH);
}

static unsigned int
store_polygon_fill_flags(const BRLObolPolygonVisual &visual)
{
    const unsigned int flags =
	store_polygon_valid_fill_flags(visual.fillFlags);
    if (flags)
	return flags;
    return visual.fill ? BRLOBOL_POLYGON_FILL_HATCH :
	BRLOBOL_POLYGON_FILL_NONE;
}

static void
store_polygon_set_fill_flags(BRLObolPolygonVisual &visual, unsigned int flags)
{
    visual.fillFlags = store_polygon_valid_fill_flags(flags);
    visual.fill = (visual.fillFlags & BRLOBOL_POLYGON_FILL_HATCH) ?
	TRUE : FALSE;
}

static SbVec3f
store_polygon_origin(const struct bg_polygon &poly, const point_t fallback)
{
    if (poly.num_contours > 0 && poly.contour && poly.contour[0].num_points > 0 &&
	    poly.contour[0].point)
	return store_vec3(poly.contour[0].point[0]);
    return store_vec3(fallback);
}

static void
store_polygon_zplane(plane_t dst, const BRLObolPolygonStoreRecord *rec)
{
    if (!rec) {
	HSET(dst, 0.0, 0.0, 1.0, 0.0);
	return;
    }

    HMOVE(dst, rec->viewPlane);
    dst[3] += rec->visual.viewZ;
}

static void
store_polygon_plane_uv(fastf_t *u, fastf_t *v, const fastf_t *plane,
	const SbVec3f &point)
{
    if (u)
	*u = 0.0;
    if (v)
	*v = 0.0;
    if (!u || !v || !plane)
	return;

    plane_t local_plane;
    point_t model_point;
    HMOVE(local_plane, plane);
    store_point(model_point, point);
    (void)bg_plane_closest_pt(u, v, &local_plane, &model_point);
}

static SbVec3f
store_polygon_plane_point(const fastf_t *plane, fastf_t u, fastf_t v)
{
    if (!plane)
	return SbVec3f(0.0f, 0.0f, 0.0f);

    plane_t local_plane;
    point_t model_point = VINIT_ZERO;
    HMOVE(local_plane, plane);
    (void)bg_plane_pt_at(&model_point, &local_plane, u, v);
    return store_vec3(model_point);
}

static SbVec3f
store_polygon_project_to_plane(const fastf_t *plane, const SbVec3f &point)
{
    fastf_t u = 0.0;
    fastf_t v = 0.0;
    store_polygon_plane_uv(&u, &v, plane, point);
    return store_polygon_plane_point(plane, u, v);
}

static SbVec3f
store_polygon_project_to_zplane(const BRLObolPolygonStoreRecord *rec,
	const SbVec3f &point)
{
    plane_t zplane;
    store_polygon_zplane(zplane, rec);
    return store_polygon_project_to_plane(zplane, point);
}

static void
store_polygon_canonical_hatch_slope(vect2d_t out, const SbVec2f &slope)
{
    V2SET(out, static_cast<fastf_t>(slope[0]),
	    static_cast<fastf_t>(slope[1]));
    if (MAG2SQ(out) < SMALL_FASTF)
	V2SET(out, 1.0, 0.0);
    V2UNITIZE(out);

    /* Hatch slopes are unoriented line families: d and -d are equivalent.
     * Canonicalize to one half-plane so small view/input changes do not flip
     * the generated perpendicular stepping direction. */
    if (out[X] < 0.0 || (NEAR_ZERO(out[X], SMALL_FASTF) && out[Y] < 0.0)) {
	out[X] = -out[X];
	out[Y] = -out[Y];
    }
}

static struct bg_polygon *
store_polygon_hatch_segments(
	const struct bg_polygon *poly,
	const plane_t *vp,
	const SbVec2f &slope,
	float spacing)
{
    if (!poly || !vp || poly->num_contours < 1 ||
	    poly->contour[0].num_points < 3 ||
	    fabs(static_cast<fastf_t>(spacing)) < BN_TOL_DIST)
	return NULL;

    vect2d_t line_slope;
    store_polygon_canonical_hatch_slope(line_slope, slope);
    const fastf_t line_spacing = fabs(static_cast<fastf_t>(spacing));

    struct bg_polygon poly_2d = BG_POLYGON_NULL;
    vect2d_t b2d_min = {MAX_FASTF, MAX_FASTF};
    vect2d_t b2d_max = {-MAX_FASTF, -MAX_FASTF};

    poly_2d.num_contours = poly->num_contours;
    poly_2d.hole = (int *)bu_calloc(poly->num_contours, sizeof(int),
	    "BRLObol hatch mask holes");
    poly_2d.contour = (struct bg_poly_contour *)bu_calloc(
	    poly->num_contours, sizeof(struct bg_poly_contour),
	    "BRLObol hatch mask contours");
    for (size_t i = 0; i < poly->num_contours; i++) {
	poly_2d.hole[i] = poly->hole[i];
	poly_2d.contour[i].num_points = poly->contour[i].num_points;
	poly_2d.contour[i].point = (point_t *)bu_calloc(
		poly->contour[i].num_points, sizeof(point_t),
		"BRLObol hatch mask points");
	for (size_t j = 0; j < poly->contour[i].num_points; j++) {
	    vect2d_t p2d;
	    point_t src_point;
	    VMOVE(src_point, poly->contour[i].point[j]);
	    bg_plane_closest_pt(&p2d[0], &p2d[1], const_cast<plane_t *>(vp),
		    &src_point);
	    VSET(poly_2d.contour[i].point[j], p2d[0], p2d[1], 0.0);
	    V2MINMAX(b2d_min, b2d_max, p2d);
	}
    }

    vect2d_t bcenter, lseg, per;
    bcenter[0] = (b2d_max[0] - b2d_min[0]) * 0.5 + b2d_min[0];
    bcenter[1] = (b2d_max[1] - b2d_min[1]) * 0.5 + b2d_min[1];
    const fastf_t ldiag = DIST_PNT2_PNT2(b2d_max, b2d_min);
    V2MOVE(lseg, line_slope);

    const int dir_step_cnt = static_cast<int>(0.5 * ldiag /
	    line_spacing + 1.0);
    if (dir_step_cnt < 2) {
	bg_polygon_free(&poly_2d);
	return NULL;
    }

    const int step_cnt = 2 * dir_step_cnt - 1;
    struct bg_polygon poly_lines = BG_POLYGON_NULL;
    poly_lines.num_contours = static_cast<size_t>(step_cnt);
    poly_lines.hole = (int *)bu_calloc(poly_lines.num_contours,
	    sizeof(int), "BRLObol hatch line holes");
    poly_lines.contour = (struct bg_poly_contour *)bu_calloc(
	    poly_lines.num_contours, sizeof(struct bg_poly_contour),
	    "BRLObol hatch line contours");

    struct bg_poly_contour *c = &poly_lines.contour[0];
    vect2d_t p2d1, p2d2;
    c->num_points = 2;
    c->open = 1;
    c->point = (point_t *)bu_calloc(c->num_points, sizeof(point_t),
	    "BRLObol hatch line points");
    V2JOIN1(p2d1, bcenter, ldiag * 0.51, lseg);
    V2JOIN1(p2d2, bcenter, -ldiag * 0.51, lseg);
    VSET(c->point[0], p2d1[0], p2d1[1], 0.0);
    VSET(c->point[1], p2d2[0], p2d2[1], 0.0);

    V2SET(per, -lseg[1], lseg[0]);
    V2UNITIZE(per);
    for (int i = 1; i < dir_step_cnt + 1; i++) {
	c = &poly_lines.contour[i];
	c->num_points = 2;
	c->open = 1;
	c->point = (point_t *)bu_calloc(c->num_points, sizeof(point_t),
		"BRLObol hatch line points");
	V2JOIN2(p2d1, bcenter, line_spacing * i, per, ldiag * 0.51, lseg);
	V2JOIN2(p2d2, bcenter, line_spacing * i, per, -ldiag * 0.51, lseg);
	VSET(c->point[0], p2d1[0], p2d1[1], 0.0);
	VSET(c->point[1], p2d2[0], p2d2[1], 0.0);
    }

    V2SET(per, lseg[1], -lseg[0]);
    V2UNITIZE(per);
    for (int i = 1 + dir_step_cnt; i < step_cnt; i++) {
	c = &poly_lines.contour[i];
	c->num_points = 2;
	c->open = 1;
	c->point = (point_t *)bu_calloc(c->num_points, sizeof(point_t),
		"BRLObol hatch line points");
	const fastf_t offset = line_spacing *
	    (static_cast<fastf_t>(i) - static_cast<fastf_t>(dir_step_cnt));
	V2JOIN2(p2d1, bcenter, offset, per, ldiag * 0.51, lseg);
	V2JOIN2(p2d2, bcenter, offset, per, -ldiag * 0.51, lseg);
	VSET(c->point[0], p2d1[0], p2d1[1], 0.0);
	VSET(c->point[1], p2d2[0], p2d2[1], 0.0);
    }

    struct bg_polygon *fpoly = bg_clip_polygon(bg_Intersection,
	    &poly_lines, &poly_2d, CLIPPER_MAX, NULL);
    if (!fpoly || !fpoly->num_contours) {
	bg_polygon_free(&poly_lines);
	bg_polygon_free(&poly_2d);
	if (fpoly) {
	    bg_polygon_free(fpoly);
	    BU_PUT(fpoly, struct bg_polygon);
	}
	return NULL;
    }

    struct bg_polygon *poly_fill = NULL;
    BU_GET(poly_fill, struct bg_polygon);
    *poly_fill = BG_POLYGON_NULL;
    poly_fill->num_contours = fpoly->num_contours;
    poly_fill->hole = (int *)bu_calloc(fpoly->num_contours, sizeof(int),
	    "BRLObol hatch output holes");
    poly_fill->contour = (struct bg_poly_contour *)bu_calloc(
	    fpoly->num_contours, sizeof(struct bg_poly_contour),
	    "BRLObol hatch output contours");
    for (size_t i = 0; i < fpoly->num_contours; i++) {
	poly_fill->hole[i] = fpoly->hole[i];
	poly_fill->contour[i].open = 1;
	poly_fill->contour[i].num_points = fpoly->contour[i].num_points;
	poly_fill->contour[i].point = (point_t *)bu_calloc(
		fpoly->contour[i].num_points, sizeof(point_t),
		"BRLObol hatch output points");
	for (size_t j = 0; j < fpoly->contour[i].num_points; j++) {
	    bg_plane_pt_at(&poly_fill->contour[i].point[j],
		    const_cast<plane_t *>(vp),
		    fpoly->contour[i].point[j][0],
		    fpoly->contour[i].point[j][1]);
	}
    }

    bg_polygon_free(&poly_lines);
    bg_polygon_free(&poly_2d);
    bg_polygon_free(fpoly);
    BU_PUT(fpoly, struct bg_polygon);

    return poly_fill;
}

static SoNode *
store_polygon_mesh_fill_node(const BRLObolPolygonStoreRecord &rec)
{
    if (!(store_polygon_fill_flags(rec.visual) & BRLOBOL_POLYGON_FILL_MESH) ||
	    rec.polygon.num_contours == 0 ||
	    !rec.polygon.contour || rec.polygon.contour[0].num_points < 3 ||
	    rec.polygon.contour[0].open)
	return NULL;

    for (size_t i = 0; i < rec.polygon.num_contours; i++) {
	if (rec.polygon.contour[i].open)
	    return NULL;
    }

    std::vector<SbVec3f> points;
    std::vector<int32_t> indices;

    if (rec.polygon.num_contours == 1) {
	const struct bg_poly_contour &contour = rec.polygon.contour[0];
	points.reserve(contour.num_points);
	for (size_t i = 0; i < contour.num_points; i++)
	    points.push_back(store_vec3(contour.point[i]));

	point_t center;
	vect_t normal;
	plane_t plane;
	if (bg_fit_plane(&center, &normal, contour.num_points, contour.point) ||
		bg_plane_pt_nrml(&plane, center, normal))
	    return NULL;

	point2d_t *projected = (point2d_t *)bu_calloc(contour.num_points,
		sizeof(point2d_t), "BRLObol projected polygon fill points");
	for (size_t i = 0; i < contour.num_points; i++)
	    bg_plane_closest_pt(&projected[i][0], &projected[i][1],
		    &plane, &contour.point[i]);

	int *faces = NULL;
	int numFaces = 0;
	int ret = bg_poly_triangulate(&faces, &numFaces, NULL, NULL, NULL, 0,
		projected, contour.num_points, TRI_EAR_CLIPPING);
	bu_free(projected, "BRLObol projected polygon fill points");
	if (ret || numFaces <= 0 || !faces) {
	    if (faces)
		bu_free(faces, "BRLObol polygon fill faces");
	    return NULL;
	}

	indices.reserve(static_cast<size_t>(numFaces) * 3);
	for (int i = 0; i < numFaces * 3; i++)
	    indices.push_back(static_cast<int32_t>(faces[i]));
	bu_free(faces, "BRLObol polygon fill faces");
    } else {
	struct bg_polygon poly = BG_POLYGON_NULL;
	bg_polygon_cpy(&poly, const_cast<struct bg_polygon *>(&rec.polygon));

	int *faces = NULL;
	int numFaces = 0;
	point_t *outPts = NULL;
	int numOutPts = 0;
	int ret = bg_polygon_triangulate(&faces, &numFaces, &outPts, &numOutPts,
		&poly, TRI_EAR_CLIPPING);
	bg_polygon_free(&poly);

	if (ret || numFaces <= 0 || numOutPts <= 0 || !faces || !outPts) {
	    if (faces)
		bu_free(faces, "BRLObol polygon fill faces");
	    if (outPts)
		bu_free(outPts, "BRLObol polygon fill points");
	    return NULL;
	}

	points.reserve(static_cast<size_t>(numOutPts));
	for (int i = 0; i < numOutPts; i++)
	    points.push_back(store_vec3(outPts[i]));

	indices.reserve(static_cast<size_t>(numFaces) * 3);
	for (int i = 0; i < numFaces * 3; i++)
	    indices.push_back(static_cast<int32_t>(faces[i]));

	bu_free(faces, "BRLObol polygon fill faces");
	bu_free(outPts, "BRLObol polygon fill points");
    }

    if (points.empty() || indices.empty())
	return NULL;

    SoBRLMeshShape *shape = new SoBRLMeshShape;
    shape->sourcePath = rec.name;
    shape->sourceName = rec.name;
    shape->sourceType = "view-polygon-mesh-fill";
    shape->displayName = rec.name;
    shape->geometryName = rec.name;
    shape->sourceIdentity = rec.name;
    shape->cacheIdentity = rec.name;
    shape->databaseIntent = FALSE;
    shape->overlayIntent = TRUE;
    shape->hudIntent = FALSE;
    shape->localSource = rec.scope == BRLObolFeatureScope::Local ? TRUE : FALSE;
    shape->sharedSource = rec.scope == BRLObolFeatureScope::Shared ? TRUE : FALSE;
    shape->nonDatabaseSource = TRUE;
    shape->drawMode = BRLOBOL_LOD_DRAW_SHADED;
    shape->recordRole = "view-polygon";
    shape->geometryKind = "surface";
    shape->sourceId = static_cast<uint32_t>(rec.revision);
    shape->transparency = 0.55f;
    store_apply_mesh_color(shape, rec.visual.fillColor);

    if (!points.empty() && !indices.empty())
	shape->setIndexedTriangles(&points[0], static_cast<int>(points.size()),
		&indices[0], static_cast<int>(indices.size()));
    return shape;
}

static SoNode *
store_polygon_hatch_fill_node(const BRLObolPolygonStoreRecord &rec)
{
    if (!(store_polygon_fill_flags(rec.visual) & BRLOBOL_POLYGON_FILL_HATCH) ||
	    rec.polygon.num_contours == 0 ||
	    !rec.polygon.contour || rec.polygon.contour[0].num_points < 3 ||
	    rec.polygon.contour[0].open)
	return NULL;

    plane_t zplane;
    store_polygon_zplane(zplane, &rec);
    struct bg_polygon *hatch = store_polygon_hatch_segments(&rec.polygon,
	    &zplane, rec.visual.fillSlope, rec.visual.fillSpacing);
    if (!hatch)
	return NULL;

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    for (size_t i = 0; i < hatch->num_contours; i++) {
	const struct bg_poly_contour &contour = hatch->contour[i];
	if (!contour.num_points || !contour.point)
	    continue;
	for (size_t j = 0; j < contour.num_points; j++) {
	    points.push_back(store_vec3(contour.point[j]));
	    commands.push_back(j == 0 ?
		    static_cast<int32_t>(BRLObolLineCommand::Move) :
		    static_cast<int32_t>(BRLObolLineCommand::Draw));
	}
    }

    bg_polygon_free(hatch);
    BU_PUT(hatch, struct bg_polygon);

    if (points.empty())
	return NULL;

    BRLObolFeatureStoreRecord feature;
    const std::string hatchName = store_string(rec.name) + ":hatch";
    feature.name = hatchName.c_str();
    feature.scope = rec.scope;
    feature.kind = BRLObolFeatureKind::PolygonOverlay;
    feature.points = points;
    feature.commands = commands;
    feature.style.hasVisible = TRUE;
    feature.style.visible = TRUE;
    feature.style.hasColor = TRUE;
    feature.style.color = rec.visual.fillColor;
    feature.style.hasLineWidth = TRUE;
    feature.style.lineWidth = 1;
    SoBRLVListShape *shape = store_vlist_node(feature);
    shape->sourcePath = rec.name;
    shape->sourceType = "view-polygon-hatch-fill";
    shape->displayName = rec.name;
    shape->geometryName = rec.name;
    shape->sourceIdentity = rec.name;
    shape->cacheIdentity = rec.name;
    shape->recordRole = "view-polygon";
    shape->geometryKind = "line";
    shape->sourceId = static_cast<uint32_t>(rec.revision);
    return shape;
}

static SoNode *
store_polygon_node(const BRLObolPolygonStoreRecord &rec)
{
    std::vector<SbVec3f> points;
    std::vector<double> precisePoints;
    std::vector<int32_t> commands;
    for (size_t i = 0; i < rec.polygon.num_contours; i++) {
	const struct bg_poly_contour &contour = rec.polygon.contour[i];
	if (!contour.num_points || !contour.point)
	    continue;
	for (size_t j = 0; j < contour.num_points; j++) {
	    points.push_back(store_vec3(contour.point[j]));
	    precisePoints.push_back(contour.point[j][X]);
	    precisePoints.push_back(contour.point[j][Y]);
	    precisePoints.push_back(contour.point[j][Z]);
	    commands.push_back(j == 0 ?
		    static_cast<int32_t>(BRLObolLineCommand::Move) :
		    static_cast<int32_t>(BRLObolLineCommand::Draw));
	}
	if (!contour.open && contour.num_points > 1) {
	    points.push_back(store_vec3(contour.point[0]));
	    precisePoints.push_back(contour.point[0][X]);
	    precisePoints.push_back(contour.point[0][Y]);
	    precisePoints.push_back(contour.point[0][Z]);
	    commands.push_back(static_cast<int32_t>(BRLObolLineCommand::Draw));
	}
    }

    BRLObolFeatureStoreRecord feature;
    feature.name = rec.name;
    feature.scope = rec.scope;
    feature.kind = BRLObolFeatureKind::PolygonOverlay;
    feature.points = points;
    feature.commands = commands;
    feature.style.hasVisible = TRUE;
    feature.style.visible = TRUE;
    feature.style.hasColor = TRUE;
    feature.style.color = rec.visual.edgeColor;
    feature.style.hasLineWidth = TRUE;
    feature.style.lineWidth = 1;
    SoBRLVListShape *shape = store_vlist_node(feature);
    shape->setPrecisePoints(precisePoints.empty() ? NULL :
	    precisePoints.data(), static_cast<int>(points.size()));
    shape->sourceType = "view-polygon-edge";
    shape->geometryKind = "line";

    SoSeparator *sep = new SoSeparator;
    SoNode *meshFillNode = store_polygon_mesh_fill_node(rec);
    if (meshFillNode)
	sep->addChild(meshFillNode);
    SoNode *hatchFillNode = store_polygon_hatch_fill_node(rec);
    if (hatchFillNode)
	sep->addChild(hatchFillNode);
    sep->addChild(shape);
    return sep;
}

static void
store_polygon_append_point(BRLObolPolygonStoreRecord *rec,
	const SbVec3f &point)
{
    if (!rec)
	return;

    SbVec3f model_point = store_polygon_project_to_plane(rec->viewPlane,
	    point);
    if (rec->polygon.num_contours == 0)
	store_polygon_init_one_point(&rec->polygon, model_point);
    else {
	size_t contourIndex = rec->currentContour >= 0 ?
	    static_cast<size_t>(rec->currentContour) : 0;
	if (contourIndex >= rec->polygon.num_contours)
	    contourIndex = rec->polygon.num_contours - 1;

	struct bg_poly_contour *contour = &rec->polygon.contour[contourIndex];
	const size_t oldCount = contour->num_points;
	contour->point = (point_t *)bu_realloc(contour->point,
		(oldCount + 1) * sizeof(point_t),
		"BRLObol polygon append point");
	store_point(contour->point[oldCount], model_point);
	contour->num_points = oldCount + 1;
    }
}

static void
store_polygon_set_rectangle(BRLObolPolygonStoreRecord *rec,
	const SbVec3f &corner,
	SbBool square)
{
    if (!rec)
	return;

    plane_t zplane;
    store_polygon_zplane(zplane, rec);

    fastf_t pfx = 0.0;
    fastf_t pfy = 0.0;
    fastf_t fx = 0.0;
    fastf_t fy = 0.0;
    store_polygon_plane_uv(&pfx, &pfy, zplane, rec->originPoint);
    store_polygon_plane_uv(&fx, &fy, zplane, corner);
    if (square) {
	const fastf_t dx = fx - pfx;
	const fastf_t dy = fy - pfy;
	const fastf_t side = std::max(std::fabs(dx), std::fabs(dy));
	fx = dx < 0.0 ? pfx - side : pfx + side;
	fy = dy < 0.0 ? pfy - side : pfy + side;
    }

    bg_polygon_free(&rec->polygon);
    rec->polygon.num_contours = 1;
    rec->polygon.hole = (int *)bu_calloc(1, sizeof(int),
	    "BRLObol rectangle hole");
    rec->polygon.contour = (struct bg_poly_contour *)bu_calloc(1,
	    sizeof(struct bg_poly_contour), "BRLObol rectangle contour");
    rec->polygon.contour[0].num_points = 4;
    rec->polygon.contour[0].open = 0;
    rec->polygon.contour[0].point = (point_t *)bu_calloc(4,
	    sizeof(point_t), "BRLObol rectangle points");
    store_point(rec->polygon.contour[0].point[0],
	    store_polygon_plane_point(zplane, pfx, pfy));
    store_point(rec->polygon.contour[0].point[1],
	    store_polygon_plane_point(zplane, pfx, fy));
    store_point(rec->polygon.contour[0].point[2],
	    store_polygon_plane_point(zplane, fx, fy));
    store_point(rec->polygon.contour[0].point[3],
	    store_polygon_plane_point(zplane, fx, pfy));
}

static void
store_polygon_set_ellipse(BRLObolPolygonStoreRecord *rec,
	const SbVec3f &corner,
	SbBool circle)
{
    if (!rec)
	return;

    const int nsegs = 64;
    plane_t zplane;
    store_polygon_zplane(zplane, rec);

    fastf_t pfx = 0.0;
    fastf_t pfy = 0.0;
    fastf_t fx = 0.0;
    fastf_t fy = 0.0;
    store_polygon_plane_uv(&pfx, &pfy, zplane, rec->originPoint);
    store_polygon_plane_uv(&fx, &fy, zplane, corner);
    fastf_t rx = fx - pfx;
    fastf_t ry = fy - pfy;
    if (circle) {
	const fastf_t r = std::sqrt(rx * rx + ry * ry);
	rx = r;
	ry = r;
    }

    bg_polygon_free(&rec->polygon);
    rec->polygon.num_contours = 1;
    rec->polygon.hole = (int *)bu_calloc(1, sizeof(int),
	    "BRLObol ellipse hole");
    rec->polygon.contour = (struct bg_poly_contour *)bu_calloc(1,
	    sizeof(struct bg_poly_contour), "BRLObol ellipse contour");
    rec->polygon.contour[0].num_points = nsegs;
    rec->polygon.contour[0].open = 0;
    rec->polygon.contour[0].point = (point_t *)bu_calloc(nsegs,
	    sizeof(point_t), "BRLObol ellipse points");
    for (int i = 0; i < nsegs; i++) {
	const double twoPi = 6.283185307179586476925286766559;
	const float a = static_cast<float>(twoPi *
		static_cast<double>(i) / static_cast<double>(nsegs));
	store_point(rec->polygon.contour[0].point[i],
		store_polygon_plane_point(zplane,
		    pfx + std::cos(a) * rx,
		    pfy + std::sin(a) * ry));
    }
}

BRLObolPolygonStore::BRLObolPolygonStore(void) : impl(new Impl)
{
}

BRLObolPolygonStore::BRLObolPolygonStore(BRLObolViewController *controller) :
    impl(new Impl)
{
    this->impl->controller = controller;
}

BRLObolPolygonStore::~BRLObolPolygonStore(void)
{
    delete this->impl;
    this->impl = NULL;
}

void
BRLObolPolygonStore::setController(BRLObolViewController *controller)
{
    this->impl->controller = controller;
}

BRLObolViewController *
BRLObolPolygonStore::controller(void) const
{
    return this->impl->controller;
}

void
BRLObolPolygonStore::clear(void)
{
    this->impl->clear();
}

BRLObolPolygonHandle
BRLObolPolygonStore::create(const SbString &name,
	BRLObolFeatureScope scope,
	BRLObolPolygonType type,
	const SbVec3f &originPoint,
	const fastf_t *viewPlane,
	float viewZ)
{
    if (store_string(name).empty())
	return BRLObolPolygonHandle();

    const std::string key = store_key(scope, name);
    if (this->impl->names.find(key) != this->impl->names.end())
	return BRLObolPolygonHandle();

    BRLObolPolygonStoreRecord *rec = new BRLObolPolygonStoreRecord;
    rec->id = this->impl->nextId++;
    rec->revision = 1;
    rec->name = name;
    rec->scope = scope;
    rec->type = type;
    rec->originPoint = originPoint;
    rec->visual.viewZ = viewZ;
    if (viewPlane)
	HMOVE(rec->viewPlane, viewPlane);
    else
	HSET(rec->viewPlane, 0.0, 0.0, 1.0, originPoint[2]);
    store_polygon_init_one_point(&rec->polygon, originPoint);
    this->impl->records[rec->id] = rec;
    this->impl->names[key] = rec->id;
    this->impl->realize(rec);
    return this->impl->handle(rec);
}

BRLObolPolygonHandle
BRLObolPolygonStore::find(const SbString &name, unsigned int scopeMask) const
{
    return this->impl->handle(this->impl->recordByName(name, scopeMask));
}

BRLObolPolygonHandle
BRLObolPolygonStore::selectAtModelPoint(const SbVec3f &point) const
{
    double best = INFINITY;
    const BRLObolPolygonStoreRecord *bestRec = NULL;
    for (std::map<uint64_t, BRLObolPolygonStoreRecord *>::const_iterator it =
	    this->impl->records.begin(); it != this->impl->records.end(); ++it) {
	BRLObolPolygonStoreRecord *rec = it->second;
	if (!rec)
	    continue;
	for (size_t i = 0; i < rec->polygon.num_contours; i++) {
	    const struct bg_poly_contour &contour = rec->polygon.contour[i];
	    for (size_t j = 0; j < contour.num_points; j++) {
		SbVec3f p = store_vec3(contour.point[j]);
		const double dx = static_cast<double>(p[0] - point[0]);
		const double dy = static_cast<double>(p[1] - point[1]);
		const double dz = static_cast<double>(p[2] - point[2]);
		const double d = dx * dx + dy * dy + dz * dz;
		if (d < best) {
		    best = d;
		    bestRec = rec;
		}
	    }
	}
    }
    return this->impl->handle(bestRec);
}

BRLObolPolygonHandle
BRLObolPolygonStore::duplicate(BRLObolPolygonHandle handle,
	const SbString &newName)
{
    BRLObolPolygonStoreRecord *src = this->impl->record(handle);
    if (!src || store_string(newName).empty())
	return BRLObolPolygonHandle();

    BRLObolPolygonHandle dstHandle = this->create(newName, src->scope,
	    src->type, src->originPoint, src->viewPlane, src->visual.viewZ);
    BRLObolPolygonStoreRecord *dst = this->impl->record(dstHandle);
    if (!dst)
	return BRLObolPolygonHandle();

    bg_polygon_free(&dst->polygon);
    bg_polygon_cpy(&dst->polygon, &src->polygon);
    dst->visual = src->visual;
    dst->currentContour = src->currentContour;
    dst->currentPoint = src->currentPoint;
    dst->revision++;
    this->impl->realize(dst);
    return this->impl->handle(dst);
}

SbBool
BRLObolPolygonStore::update(BRLObolPolygonHandle handle,
	BRLObolPolygonUpdate update)
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    if (update == BRLObolPolygonUpdate::PointSelectClear) {
	rec->currentContour = -1;
	rec->currentPoint = -1;
    }
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BRLObolPolygonStore::updateScreenPoint(BRLObolPolygonHandle handle,
	int x,
	int y,
	BRLObolPolygonUpdate update)
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    return this->updateModelPoint(handle,
	    SbVec3f(static_cast<float>(x), static_cast<float>(y),
		rec->originPoint[2]), update);
}

SbBool
BRLObolPolygonStore::updateModelPoint(BRLObolPolygonHandle handle,
	const SbVec3f &point,
	BRLObolPolygonUpdate update)
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    if (update == BRLObolPolygonUpdate::PointAppend) {
	store_polygon_append_point(rec, point);
	rec->currentPoint = -1;
    } else if (update == BRLObolPolygonUpdate::PointMove) {
	if (rec->currentContour < 0 || rec->currentPoint < 0)
	    return FALSE;
	struct bg_poly_contour &contour =
	    rec->polygon.contour[rec->currentContour];
	if (static_cast<size_t>(rec->currentPoint) >= contour.num_points)
	    return FALSE;
	store_point(contour.point[rec->currentPoint],
		store_polygon_project_to_zplane(rec, point));
    } else if (update == BRLObolPolygonUpdate::PointSelect) {
	SbVec3f model_point = store_polygon_project_to_zplane(rec, point);
	const long selected_contour = rec->currentContour;
	double best = INFINITY;
	rec->currentContour = -1;
	rec->currentPoint = -1;
	size_t start = 0;
	size_t end = rec->polygon.num_contours;
	if (selected_contour >= 0 &&
		static_cast<size_t>(selected_contour) < end) {
	    start = static_cast<size_t>(selected_contour);
	    end = start + 1;
	}
	for (size_t i = start; i < end; i++) {
	    const struct bg_poly_contour &contour = rec->polygon.contour[i];
	    for (size_t j = 0; j < contour.num_points; j++) {
		SbVec3f p = store_vec3(contour.point[j]);
		const double dx = static_cast<double>(p[0] - model_point[0]);
		const double dy = static_cast<double>(p[1] - model_point[1]);
		const double dz = static_cast<double>(p[2] - model_point[2]);
		const double d = dx * dx + dy * dy + dz * dz;
		if (d < best) {
		    best = d;
		    rec->currentContour = static_cast<long>(i);
		    rec->currentPoint = static_cast<long>(j);
		}
	    }
	}
    } else if (rec->type == BRLObolPolygonType::Circle) {
	store_polygon_set_ellipse(rec, point, TRUE);
    } else if (rec->type == BRLObolPolygonType::Ellipse) {
	store_polygon_set_ellipse(rec, point, FALSE);
    } else if (rec->type == BRLObolPolygonType::Square) {
	store_polygon_set_rectangle(rec, point, TRUE);
    } else if (rec->type == BRLObolPolygonType::Rectangle) {
	store_polygon_set_rectangle(rec, point, FALSE);
    }

    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BRLObolPolygonStore::move(BRLObolPolygonHandle handle,
	const SbVec3f &currentPoint,
	const SbVec3f &previousPoint)
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    SbVec3f current = store_polygon_project_to_zplane(rec, currentPoint);
    SbVec3f previous = store_polygon_project_to_zplane(rec, previousPoint);
    SbVec3f delta = current - previous;
    for (size_t i = 0; i < rec->polygon.num_contours; i++) {
	struct bg_poly_contour &contour = rec->polygon.contour[i];
	for (size_t j = 0; j < contour.num_points; j++) {
	    contour.point[j][X] += delta[0];
	    contour.point[j][Y] += delta[1];
	    contour.point[j][Z] += delta[2];
	}
    }
    rec->originPoint += delta;
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BRLObolPolygonStore::rename(BRLObolPolygonHandle handle,
	const SbString &newName)
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec || store_string(newName).empty())
	return FALSE;

    const std::string newKey = store_key(rec->scope, newName);
    std::map<std::string, uint64_t>::const_iterator existing =
	this->impl->names.find(newKey);
    if (existing != this->impl->names.end() && existing->second != rec->id)
	return FALSE;

    this->impl->names.erase(store_key(rec->scope, rec->name));
    rec->name = newName;
    rec->revision++;
    this->impl->names[newKey] = rec->id;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BRLObolPolygonStore::remove(BRLObolPolygonHandle handle)
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    this->impl->names.erase(store_key(rec->scope, rec->name));
    store_release_node(this->impl->controller, rec->node);
    bg_polygon_free(&rec->polygon);
    this->impl->records.erase(rec->id);
    delete rec;
    return TRUE;
}

SbBool
BRLObolPolygonStore::record(BRLObolPolygonHandle handle,
	BRLObolPolygonRecord &recordOut) const
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    recordOut = BRLObolPolygonRecord();
    recordOut.handle = this->impl->handle(rec);
    recordOut.name = rec->name;
    recordOut.scope = rec->scope;
    recordOut.type = rec->type;
    recordOut.fillFlags = store_polygon_fill_flags(rec->visual);
    recordOut.fill = (recordOut.fillFlags & BRLOBOL_POLYGON_FILL_HATCH) ?
	TRUE : FALSE;
    recordOut.fillSlope = rec->visual.fillSlope;
    recordOut.fillSpacing = rec->visual.fillSpacing;
    recordOut.fillColor = rec->visual.fillColor;
    recordOut.edgeColor = rec->visual.edgeColor;
    recordOut.currentContour = rec->currentContour;
    recordOut.currentPoint = rec->currentPoint;
    recordOut.contourCount = rec->polygon.num_contours;
    recordOut.pointCount = store_polygon_point_count(rec->polygon);
    recordOut.originPoint = rec->originPoint;
    HMOVE(recordOut.viewPlane, rec->viewPlane);
    recordOut.viewZ = rec->visual.viewZ;
    recordOut.userData = rec->userData;
    if (rec->polygon.num_contours > 0)
	recordOut.firstContourOpen = rec->polygon.contour[0].open ? TRUE : FALSE;
    return TRUE;
}

void
BRLObolPolygonStore::visitRecords(BRLObolPolygonRecordCallback callback,
	void *userData) const
{
    if (!callback)
	return;
    for (std::map<uint64_t, BRLObolPolygonStoreRecord *>::const_iterator it =
	    this->impl->records.begin(); it != this->impl->records.end(); ++it) {
	BRLObolPolygonRecord rec;
	if (this->record(this->impl->handle(it->second), rec) &&
		!callback(rec, userData))
	    return;
    }
}

SbBool
BRLObolPolygonStore::setCurrent(BRLObolPolygonHandle handle,
	long contour,
	long point)
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    if (contour >= static_cast<long>(rec->polygon.num_contours))
	return FALSE;
    if (contour >= 0 && point >= 0) {
	const struct bg_poly_contour &c = rec->polygon.contour[contour];
	if (point >= static_cast<long>(c.num_points))
	    return FALSE;
    }
    rec->currentContour = contour;
    rec->currentPoint = point;
    rec->revision++;
    return TRUE;
}

SbBool
BRLObolPolygonStore::setContourOpen(BRLObolPolygonHandle handle,
	long contour,
	SbBool open)
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec || contour < 0 ||
	    contour >= static_cast<long>(rec->polygon.num_contours))
	return FALSE;
    rec->polygon.contour[contour].open = open ? 1 : 0;
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BRLObolPolygonStore::setAllContoursOpen(BRLObolPolygonHandle handle,
	SbBool open)
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    for (size_t i = 0; i < rec->polygon.num_contours; i++)
	rec->polygon.contour[i].open = open ? 1 : 0;
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BRLObolPolygonStore::clearSelectedPoint(BRLObolPolygonHandle handle)
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    rec->currentContour = -1;
    rec->currentPoint = -1;
    rec->revision++;
    return TRUE;
}

SbBool
BRLObolPolygonStore::clearAllPointSelections(void)
{
    for (std::map<uint64_t, BRLObolPolygonStoreRecord *>::iterator it =
	    this->impl->records.begin(); it != this->impl->records.end(); ++it) {
	if (it->second) {
	    it->second->currentContour = -1;
	    it->second->currentPoint = -1;
	    it->second->revision++;
	}
    }
    return TRUE;
}

SbBool
BRLObolPolygonStore::setVisual(BRLObolPolygonHandle handle,
	const BRLObolPolygonVisual &visual)
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    rec->visual = visual;
    store_polygon_set_fill_flags(rec->visual,
	    store_polygon_fill_flags(rec->visual));
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BRLObolPolygonStore::visual(BRLObolPolygonHandle handle,
	BRLObolPolygonVisual &visualOut) const
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    visualOut = rec->visual;
    return TRUE;
}

SbBool
BRLObolPolygonStore::setEdgeColor(BRLObolPolygonHandle handle,
	const SbColor &edgeColor)
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    rec->visual.edgeColor = edgeColor;
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BRLObolPolygonStore::edgeColor(BRLObolPolygonHandle handle,
	SbColor &edgeColorOut) const
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    edgeColorOut = rec->visual.edgeColor;
    return TRUE;
}

SbBool
BRLObolPolygonStore::setFill(BRLObolPolygonHandle handle,
	SbBool fill,
	const SbVec2f &slope,
	float spacing)
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    unsigned int flags = store_polygon_fill_flags(rec->visual);
    if (fill)
	flags |= BRLOBOL_POLYGON_FILL_HATCH;
    else
	flags &= ~BRLOBOL_POLYGON_FILL_HATCH;
    store_polygon_set_fill_flags(rec->visual, flags);
    rec->visual.fillSlope = slope;
    rec->visual.fillSpacing = spacing;
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BRLObolPolygonStore::setFillFlags(BRLObolPolygonHandle handle,
	unsigned int fillFlags)
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    store_polygon_set_fill_flags(rec->visual, fillFlags);
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BRLObolPolygonStore::setFillColor(BRLObolPolygonHandle handle,
	const SbColor &fillColor)
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    rec->visual.fillColor = fillColor;
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BRLObolPolygonStore::fillColor(BRLObolPolygonHandle handle,
	SbColor &fillColorOut) const
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    fillColorOut = rec->visual.fillColor;
    return TRUE;
}

SbBool
BRLObolPolygonStore::setGeometry(BRLObolPolygonHandle handle,
	const struct bg_polygon *polygon)
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec || !polygon)
	return FALSE;
    bg_polygon_free(&rec->polygon);
    bg_polygon_cpy(&rec->polygon, const_cast<struct bg_polygon *>(polygon));
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

const struct bg_polygon *
BRLObolPolygonStore::geometry(BRLObolPolygonHandle handle) const
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    return rec ? &rec->polygon : NULL;
}

SbBool
BRLObolPolygonStore::copyGeometry(BRLObolPolygonHandle handle,
	struct bg_polygon *polygonOut) const
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec || !polygonOut)
	return FALSE;
    bg_polygon_cpy(polygonOut, &rec->polygon);
    return TRUE;
}

SbBool
BRLObolPolygonStore::setUserData(BRLObolPolygonHandle handle, void *userData)
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    rec->userData = userData;
    rec->revision++;
    return TRUE;
}

void *
BRLObolPolygonStore::userData(BRLObolPolygonHandle handle) const
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    return rec ? rec->userData : NULL;
}

double
BRLObolPolygonStore::area(BRLObolPolygonHandle handle, double viewScale) const
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return 0.0;
    return static_cast<double>(bg_find_polygon_area(&rec->polygon,
	    CLIPPER_MAX, &rec->viewPlane, viewScale));
}

SbBool
BRLObolPolygonStore::overlaps(BRLObolPolygonHandle a,
	BRLObolPolygonHandle b,
	const struct bn_tol &tol,
	double viewScale) const
{
    BRLObolPolygonStoreRecord *ra = this->impl->record(a);
    BRLObolPolygonStoreRecord *rb = this->impl->record(b);
    if (!ra || !rb)
	return FALSE;
    return bg_polygons_overlap(&ra->polygon, &rb->polygon, &ra->viewPlane,
	    &tol, viewScale) ? TRUE : FALSE;
}

SbBool
BRLObolPolygonStore::csg(BRLObolPolygonHandle target,
	BRLObolPolygonHandle stencil,
	bg_clip_t op)
{
    BRLObolPolygonStoreRecord *rt = this->impl->record(target);
    BRLObolPolygonStoreRecord *rs = this->impl->record(stencil);
    if (!rt || !rs)
	return FALSE;
    if (op == bg_None || rs->polygon.num_contours == 0)
	return FALSE;

    if (rt->polygon.num_contours == 0) {
	if (op != bg_Union)
	    return FALSE;
	bg_polygon_cpy(&rt->polygon, &rs->polygon);
	rt->type = rs->type;
	rt->currentContour = rs->currentContour;
	rt->currentPoint = rs->currentPoint;
	rt->originPoint = rs->originPoint;
	HMOVE(rt->viewPlane, rs->viewPlane);
	rt->visual.viewZ = rs->visual.viewZ;
	rt->revision++;
	this->impl->realize(rt);
	return TRUE;
    }

    const struct bn_tol tol = BN_TOL_INIT_TOL;
    if (!bg_polygons_overlap(&rt->polygon, &rs->polygon, &rt->viewPlane,
	    &tol, CLIPPER_MAX))
	return FALSE;

    struct bg_polygon *result = bg_clip_polygon(op, &rt->polygon,
	    &rs->polygon, CLIPPER_MAX, &rt->viewPlane);
    if (!result)
	return FALSE;

    bg_polygon_free(&rt->polygon);
    bg_polygon_cpy(&rt->polygon, result);
    bg_polygon_free(result);
    bu_free(result, "BRLObol polygon CSG result");
    rt->type = BRLObolPolygonType::General;
    rt->revision++;
    this->impl->realize(rt);
    return TRUE;
}

BRLObolPolygonHandle
BRLObolPolygonStore::importSketch(const SbString &name,
	BRLObolFeatureScope scope,
	struct db_i *dbip,
	struct directory *dp)
{
    if (store_string(name).empty() || !dbip || !dp)
	return BRLObolPolygonHandle();

    const std::string key = store_key(scope, name);
    if (this->impl->names.find(key) != this->impl->names.end())
	return BRLObolPolygonHandle();

    struct rt_sketch_polygon_data data;
    rt_sketch_polygon_data_init(&data);
    if (db_sketch_to_polygon_data(&data, store_string(name).c_str(),
		dbip, dp) != 0) {
	rt_sketch_polygon_data_free(&data);
	return BRLObolPolygonHandle();
    }

    if (data.polygon.num_contours == 0 ||
	    store_polygon_point_count(data.polygon) == 0) {
	rt_sketch_polygon_data_free(&data);
	return BRLObolPolygonHandle();
    }

    BRLObolPolygonStoreRecord *rec = new BRLObolPolygonStoreRecord;
    rec->id = this->impl->nextId++;
    rec->revision = 1;
    rec->name = name;
    rec->scope = scope;
    rec->type = store_polygon_type_from_rt(data.type);
    rec->originPoint = store_polygon_origin(data.polygon, data.origin_point);
    HMOVE(rec->viewPlane, data.vp);
    store_polygon_set_fill_flags(rec->visual,
	    data.fill_flag ? BRLOBOL_POLYGON_FILL_HATCH :
	    BRLOBOL_POLYGON_FILL_NONE);
    rec->visual.fillSlope = SbVec2f(static_cast<float>(data.fill_dir[0]),
	    static_cast<float>(data.fill_dir[1]));
    rec->visual.fillSpacing = static_cast<float>(data.fill_delta);
    rec->visual.fillColor = store_bu_to_sbcolor(data.fill_color);
    if (data.have_edge_color)
	rec->visual.edgeColor = store_bu_to_sbcolor(data.edge_color);
    rec->visual.viewZ = static_cast<float>(data.vZ);
    bg_polygon_cpy(&rec->polygon, &data.polygon);

    rt_sketch_polygon_data_free(&data);

    this->impl->records[rec->id] = rec;
    this->impl->names[key] = rec->id;
    this->impl->realize(rec);
    return this->impl->handle(rec);
}

SbBool
BRLObolPolygonStore::exportSketch(BRLObolPolygonHandle handle,
	struct db_i *dbip,
	const SbString &name) const
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec || !dbip || store_string(name).empty())
	return FALSE;

    struct rt_sketch_polygon_data data;
    rt_sketch_polygon_data_init(&data);
    data.type = store_polygon_type_to_rt(rec->type);
    data.fill_flag =
	(store_polygon_fill_flags(rec->visual) & BRLOBOL_POLYGON_FILL_HATCH) ?
	1 : 0;
    V2SET(data.fill_dir, rec->visual.fillSlope[0], rec->visual.fillSlope[1]);
    data.fill_delta = rec->visual.fillSpacing;
    store_sbcolor_to_bu(rec->visual.fillColor, &data.fill_color);
    store_point(data.origin_point, rec->originPoint);
    HMOVE(data.vp, rec->viewPlane);
    data.vZ = rec->visual.viewZ;
    data.have_edge_color = 1;
    store_sbcolor_to_bu(rec->visual.edgeColor, &data.edge_color);
    bg_polygon_cpy(&data.polygon, &rec->polygon);

    struct directory *dp = db_sketch_polygon_data_to_sketch(dbip,
	    store_string(name).c_str(), &data);
    rt_sketch_polygon_data_free(&data);
    return dp ? TRUE : FALSE;
}

size_t
BRLObolPolygonStore::snapCount(BRLObolPolygonHandle exclude) const
{
    size_t count = 0;
    for (std::map<uint64_t, BRLObolPolygonStoreRecord *>::const_iterator it =
	    this->impl->records.begin(); it != this->impl->records.end(); ++it) {
	if (!it->second)
	    continue;
	BRLObolPolygonHandle handle = this->impl->handle(it->second);
	if (exclude.isValid() && handle.id == exclude.id)
	    continue;
	count++;
    }
    return count;
}

SbBool
BRLObolPolygonStore::setSnapExclude(BRLObolPolygonHandle handle)
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    this->impl->snapExclude = this->impl->handle(rec);
    return TRUE;
}

BRLObolPolygonHandle
BRLObolPolygonStore::snapExclude(void) const
{
    return this->impl->snapExclude;
}

SoNode *
BRLObolPolygonStore::node(BRLObolPolygonHandle handle) const
{
    BRLObolPolygonStoreRecord *rec = this->impl->record(handle);
    return rec ? rec->node : NULL;
}

struct BRLObolSelectionStore::Impl {
    std::vector<BRLObolSelectionRecord> records;
};

BRLObolSelectionStore::BRLObolSelectionStore(void) : impl(new Impl)
{
}

BRLObolSelectionStore::~BRLObolSelectionStore(void)
{
    delete this->impl;
    this->impl = NULL;
}

void
BRLObolSelectionStore::clear(const BRLObolFeatureOwner *owner, int kind)
{
    if (!owner && kind == BRLOBOL_SELECTION_ALL) {
	this->impl->records.clear();
	return;
    }

    this->impl->records.erase(std::remove_if(
		this->impl->records.begin(),
		this->impl->records.end(),
		[owner, kind](const BRLObolSelectionRecord &rec) {
		    return store_owner_matches(rec.owner, owner) &&
			store_selection_kind_matches(rec.kind, kind);
		}),
	    this->impl->records.end());
}

size_t
BRLObolSelectionStore::count(const BRLObolFeatureOwner *owner, int kind) const
{
    if (!owner && kind == BRLOBOL_SELECTION_ALL)
	return this->impl->records.size();

    size_t cnt = 0;
    for (size_t i = 0; i < this->impl->records.size(); i++) {
	const BRLObolSelectionRecord &rec = this->impl->records[i];
	if (store_owner_matches(rec.owner, owner) &&
		store_selection_kind_matches(rec.kind, kind))
	    cnt++;
    }
    return cnt;
}

SbBool
BRLObolSelectionStore::containsPath(const SbString &path, int kind,
	const BRLObolFeatureOwner *owner) const
{
    const std::string target = store_string(path);
    if (target.empty())
	return FALSE;
    for (size_t i = 0; i < this->impl->records.size(); i++) {
	const BRLObolSelectionRecord &rec = this->impl->records[i];
	if (store_string(rec.path) == target &&
		store_owner_matches(rec.owner, owner) &&
		store_selection_kind_matches(rec.kind, kind))
	    return TRUE;
    }
    return FALSE;
}

SbBool
BRLObolSelectionStore::addPath(const SbString &path, int kind,
	const BRLObolFeatureOwner *owner)
{
    if (store_string(path).empty())
	return FALSE;
    int recordKind = store_selection_record_kind(kind);
    if (this->containsPath(path, recordKind, owner))
	return TRUE;
    BRLObolSelectionRecord rec;
    rec.path = path;
    rec.kind = recordKind;
    if (owner)
	rec.owner = *owner;
    this->impl->records.push_back(rec);
    return TRUE;
}

SbBool
BRLObolSelectionStore::setPath(const SbString &path, int kind,
	const BRLObolFeatureOwner *owner)
{
    this->clear(owner, kind);
    if (store_string(path).empty())
	return TRUE;
    return this->addPath(path, kind, owner);
}

SbBool
BRLObolSelectionStore::removePath(const SbString &path, int kind,
	const BRLObolFeatureOwner *owner)
{
    const std::string target = store_string(path);
    for (std::vector<BRLObolSelectionRecord>::iterator it =
	    this->impl->records.begin(); it != this->impl->records.end(); ++it) {
	if (store_string(it->path) == target &&
		store_owner_matches(it->owner, owner) &&
		store_selection_kind_matches(it->kind, kind)) {
	    this->impl->records.erase(it);
	    return TRUE;
	}
    }
    return FALSE;
}

const BRLObolSelectionRecord *
BRLObolSelectionStore::record(size_t index) const
{
    return index < this->impl->records.size() ?
	&this->impl->records[index] : NULL;
}

SbBool
BRLObolSelectionStore::addRecord(const BRLObolSelectionRecord &record)
{
    if (record.path.getLength() == 0)
	return FALSE;
    BRLObolSelectionRecord rec = record;
    if (rec.kind == BRLOBOL_SELECTION_ALL)
	rec.kind = BRLOBOL_SELECTION_SELECTED_PATH;
    if (this->containsPath(rec.path, rec.kind, &rec.owner))
	return TRUE;
    this->impl->records.push_back(rec);
    return TRUE;
}

SbBool
BRLObolSelectionStore::setRecords(
	const std::vector<BRLObolSelectionRecord> &records)
{
    this->impl->records.clear();
    for (size_t i = 0; i < records.size(); i++)
	this->addRecord(records[i]);
    return TRUE;
}

void
BRLObolSelectionStore::visitPaths(
	int (*callback)(const SbString &path, void *userData),
	void *userData,
	const BRLObolFeatureOwner *owner,
	int kind) const
{
    if (!callback)
	return;
    for (size_t i = 0; i < this->impl->records.size(); i++) {
	const BRLObolSelectionRecord &rec = this->impl->records[i];
	if (!store_owner_matches(rec.owner, owner) ||
		!store_selection_kind_matches(rec.kind, kind))
	    continue;
	if (!callback(rec.path, userData))
	    return;
    }
}

SbBool
BRLObolSelectionStore::applyPickResults(
	const std::vector<BRLObolSelectionRecord> &records,
	void (*selectedPathCallback)(const SbString &path, void *userData),
	void *userData,
	const BRLObolFeatureOwner *owner)
{
    this->clear(owner, BRLOBOL_SELECTION_ALL);
    for (size_t i = 0; i < records.size(); i++) {
	BRLObolSelectionRecord rec = records[i];
	if (owner)
	    rec.owner = *owner;
	this->addRecord(rec);
	if (selectedPathCallback &&
		store_selection_record_kind(rec.kind) ==
		BRLOBOL_SELECTION_SELECTED_PATH)
	    selectedPathCallback(rec.path, userData);
    }
    return TRUE;
}
