/*                  V I E W _ S T O R E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file view_store.cpp */

#include "common.h"

#include "bu/str.h"

#include "BObol/BAxes.h"
#include "BObol/BEditPreview.h"
#include "BObol/BHUDLabelOverlay.h"
#include "BObol/BLineLayerOverlay.h"
#include "BObol/BLodRealization.h"
#include "BObol/BMeshShape.h"
#include "BObol/BSceneGroup.h"
#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"
#include "BObol/BVListShape.h"

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
#include <Inventor/annex/HUD/nodekits/SoHUDKit.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>

static std::atomic<uint64_t> store_reference_generation_counter(1);

static uint64_t
store_reference_generation_next(void)
{
    uint64_t generation =
	store_reference_generation_counter.fetch_add(1,
		std::memory_order_relaxed);

    /* Zero is reserved for null references.  Exhaustion is not realistic,
     * but keeping the invariant here makes malformed/null checks exact. */
    while (!generation)
	generation = store_reference_generation_counter.fetch_add(1,
		std::memory_order_relaxed);
    return generation;
}

static std::string
store_string(const SbString &s)
{
    const char *str = s.getString();
    return str ? std::string(str) : std::string();
}

static std::string
store_owner_key(const BObolFeatureOwner *owner)
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
store_owner_generation_key(const BObolFeatureOwner *owner)
{
    if (!owner)
	return std::string();

    if (owner->ownerToken) {
	char buf[64] = {0};
	snprintf(buf, sizeof(buf), "T:%p", owner->ownerToken);
	return std::string(buf);
    }

    const char *id = owner->ownerId.getString();
    const char *role = owner->ownerRole.getString();
    if ((id && id[0]) || (role && role[0]))
	return std::string("I:") + (id ? id : "") + "|R:" +
	       (role ? role : "");

    return std::string();
}

static std::string
store_key(BObolFeatureScope scope,
	  const SbString &name,
	  const BObolFeatureOwner *owner = NULL)
{
    if (scope != BObolFeatureScope::Local)
	return std::string("S:") + store_string(name);

    return std::string("L:") + store_owner_key(owner) + ":" +
	   store_string(name);
}

static SbBool
store_owner_matches(const BObolFeatureOwner &recordOwner,
		    const BObolFeatureOwner *queryOwner)
{
    if (!queryOwner)
	return TRUE;

    if (queryOwner->ownerToken)
	return recordOwner.ownerToken == queryOwner->ownerToken ? TRUE : FALSE;

    const char *queryId = queryOwner->ownerId.getString();
    if (queryId && queryId[0]) {
	const char *recordId = recordOwner.ownerId.getString();
	return recordId && bu_strcmp(recordId, queryId) == 0 ? TRUE : FALSE;
    }

    const char *queryRole = queryOwner->ownerRole.getString();
    if (queryRole && queryRole[0]) {
	const char *recordRole = recordOwner.ownerRole.getString();
	return recordRole && bu_strcmp(recordRole, queryRole) == 0 ? TRUE : FALSE;
    }

    return TRUE;
}

static SbBool
store_selection_kind_matches(int recordKind, int queryKind)
{
    return queryKind == BOBOL_SELECTION_ALL || recordKind == queryKind ?
	   TRUE : FALSE;
}

static int
store_selection_record_kind(int kind)
{
    return kind == BOBOL_SELECTION_ALL ?
	   BOBOL_SELECTION_SELECTED_PATH : kind;
}

static unsigned int
store_scope_bit(BObolFeatureScope scope)
{
    return scope == BObolFeatureScope::Local ?
	   BOBOL_FEATURE_SCOPE_LOCAL : BOBOL_FEATURE_SCOPE_SHARED;
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
    if (command == static_cast<int32_t>(BObolLineCommand::Point))
	return SoBRLVListShape::POINT;
    if (command == static_cast<int32_t>(BObolLineCommand::Draw))
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
			      BObolLineCommand::Draw));
	if (!normalized.empty())
	    normalized[0] = static_cast<int32_t>(BObolLineCommand::Move);
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

static int32_t
store_line_command_from_bg(int command, size_t pointIndex)
{
    switch (command) {
	case BG_GEOMETRY_LINE_MOVE:
	    return static_cast<int32_t>(BObolLineCommand::Move);
	case BG_GEOMETRY_LINE_DRAW:
	    return static_cast<int32_t>(BObolLineCommand::Draw);
	case BG_GEOMETRY_POINT_DRAW:
	    return static_cast<int32_t>(BObolLineCommand::Point);
	default:
	    break;
    }
    return pointIndex ? static_cast<int32_t>(BObolLineCommand::Draw) :
	   static_cast<int32_t>(BObolLineCommand::Move);
}

static std::vector<BObolLineLayer>
store_line_layers_from_builder(const SbString &featureName,
			       const struct bg_line_layer_builder *builder)
{
    std::vector<BObolLineLayer> layers;
    const size_t layerCount = bg_line_layer_builder_layer_count(builder);
    layers.reserve(layerCount);

    for (size_t i = 0; i < layerCount; i++) {
	const struct bg_line_layer *bgLayer =
	    bg_line_layer_builder_layer_at(builder, i);
	const size_t pointCount = bg_line_layer_point_count(bgLayer);
	if (!bgLayer || !pointCount)
	    continue;

	BObolLineLayer layer;
	const char *layerName = bg_line_layer_name(bgLayer);
	std::string fullLayerName = store_string(featureName);
	if (layerName && layerName[0]) {
	    fullLayerName += "/";
	    fullLayerName += layerName;
	}
	layer.name = fullLayerName.c_str();
	layer.points.reserve(pointCount);
	layer.commands.reserve(pointCount);

	unsigned char r = 255;
	unsigned char g = 255;
	unsigned char b = 255;
	if (bg_line_layer_color(bgLayer, &r, &g, &b)) {
	    layer.style.hasColor = TRUE;
	    layer.style.color = SbColor(
				    static_cast<float>(r) / 255.0f,
				    static_cast<float>(g) / 255.0f,
				    static_cast<float>(b) / 255.0f);
	}

	const point_t *points = bg_line_layer_points(bgLayer);
	const int *commands = bg_line_layer_commands(bgLayer);
	for (size_t j = 0; j < pointCount; j++) {
	    if (points) {
		layer.points.push_back(SbVec3f(
					   static_cast<float>(points[j][0]),
					   static_cast<float>(points[j][1]),
					   static_cast<float>(points[j][2])));
	    }
	    const int command = commands ? commands[j] : -1;
	    layer.commands.push_back(store_line_command_from_bg(command, j));
	}

	if (!layer.points.empty())
	    layers.push_back(layer);
    }

    return layers;
}

static void
store_apply_vlist_style(SoBRLVListShape *shape,
			const BObolFeatureStyle &style)
{
    if (!shape)
	return;
    if (style.hasVisible)
	shape->visible = style.visible;
    if (style.hasSelectable)
	shape->selectable = style.selectable;
    if (style.hasColor) {
	shape->colorOverride = TRUE;
	shape->color = style.color;
    }
    if (style.hasLineWidth)
	shape->lineWidth = style.lineWidth;
    if (style.hasLineStyle)
	shape->lineStyle = style.lineStyle;
    if (style.hasTransparency)
	shape->transparency = style.transparency;
}

static void
store_apply_mesh_style(SoBRLMeshShape *shape,
		       const BObolFeatureStyle &style)
{
    if (!shape)
	return;
    if (style.hasVisible)
	shape->visible = style.visible;
    if (style.hasSelectable)
	shape->selectable = style.selectable;
    if (style.hasColor) {
	shape->colorOverride = TRUE;
	shape->color = style.color;
    }
    if (style.hasLineWidth)
	shape->lineWidth = style.lineWidth;
    if (style.hasLineStyle)
	shape->lineStyle = style.lineStyle;
    if (style.hasTransparency)
	shape->transparency = style.transparency;
}

static void
store_apply_mesh_color(SoBRLMeshShape *shape, const SbColor &color)
{
    if (!shape)
	return;
    shape->colorOverride = TRUE;
    shape->color = color;
}

static BObolFeatureStyle
store_merge_feature_style(const BObolFeatureStyle &base,
			  const BObolFeatureStyle &overrideStyle)
{
    BObolFeatureStyle out = base;
    if (overrideStyle.hasVisible) {
	out.hasVisible = TRUE;
	out.visible = overrideStyle.visible;
    }
    if (overrideStyle.hasSelectable) {
	out.hasSelectable = TRUE;
	out.selectable = overrideStyle.selectable;
    }
    if (overrideStyle.hasColor) {
	out.hasColor = TRUE;
	out.color = overrideStyle.color;
    }
    if (overrideStyle.hasLineWidth) {
	out.hasLineWidth = TRUE;
	out.lineWidth = overrideStyle.lineWidth;
    }
    if (overrideStyle.hasLineStyle) {
	out.hasLineStyle = TRUE;
	out.lineStyle = overrideStyle.lineStyle;
    }
    if (overrideStyle.hasTransparency) {
	out.hasTransparency = TRUE;
	out.transparency = overrideStyle.transparency;
    }
    if (overrideStyle.hasArrow) {
	out.hasArrow = TRUE;
	out.arrow = overrideStyle.arrow;
    }
    if (overrideStyle.hasArrowTip) {
	out.hasArrowTip = TRUE;
	out.arrowTipLength = overrideStyle.arrowTipLength;
	out.arrowTipWidth = overrideStyle.arrowTipWidth;
    }
    /* HUD is a whole-feature trait; a per-layer override can only turn it on. */
    if (overrideStyle.hud)
	out.hud = TRUE;
    return out;
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
store_controller_root_group(BObolViewController *controller)
{
    if (!controller)
	return NULL;

    SoNode *root = controller->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return NULL;

    return static_cast<SoGroup *>(root);
}

static void
store_detach_node(BObolViewController *controller, SoNode *node)
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
store_attach_node(BObolViewController *controller, SoNode *node)
{
    SoGroup *group = store_controller_root_group(controller);
    if (group && node)
	group->addChild(node);
}

static void
store_release_node(BObolViewController *controller, SoNode *node)
{
    if (!node)
	return;

    store_detach_node(controller, node);
    node->unref();
}

static void
store_set_node(BObolViewController *controller, SoNode *&slot, SoNode *node)
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

BObolFeatureHandle::BObolFeatureHandle(void) : id(0), revision(0)
{
}

BObolFeatureHandle::BObolFeatureHandle(uint64_t featureId,
	uint64_t featureRevision) : id(featureId), revision(featureRevision)
{
}

SbBool
BObolFeatureHandle::isValid(void) const
{
    return id != 0 && revision != 0 ? TRUE : FALSE;
}

BObolPolygonHandle::BObolPolygonHandle(void) : id(0), revision(0)
{
}

BObolPolygonHandle::BObolPolygonHandle(uint64_t polygonId,
	uint64_t polygonRevision) : id(polygonId), revision(polygonRevision)
{
}

SbBool
BObolPolygonHandle::isValid(void) const
{
    return id != 0 && revision != 0 ? TRUE : FALSE;
}

SbBool
operator==(const BObolFeatureHandle &a, const BObolFeatureHandle &b)
{
    return a.id == b.id && a.revision == b.revision ? TRUE : FALSE;
}

SbBool
operator!=(const BObolFeatureHandle &a, const BObolFeatureHandle &b)
{
    return !(a == b);
}

SbBool
operator==(const BObolPolygonHandle &a, const BObolPolygonHandle &b)
{
    return a.id == b.id && a.revision == b.revision ? TRUE : FALSE;
}

SbBool
operator!=(const BObolPolygonHandle &a, const BObolPolygonHandle &b)
{
    return !(a == b);
}

BObolFeatureStyle::BObolFeatureStyle(void) :
    hasVisible(FALSE),
    visible(TRUE),
    hasSelectable(FALSE),
    selectable(TRUE),
    hasColor(FALSE),
    color(1.0f, 1.0f, 1.0f),
    hasLineWidth(FALSE),
    lineWidth(1),
    hasLineStyle(FALSE),
    lineStyle(0),
    hasTransparency(FALSE),
    transparency(0.0f),
    hasArrow(FALSE),
    arrow(FALSE),
    hasArrowTip(FALSE),
    arrowTipLength(0.0f),
    arrowTipWidth(0.0f),
    hud(FALSE)
{
}

BObolCommandResult::BObolCommandResult(void) :
    status(BObolCommandResultStatus::None),
    feature(),
    polygon(),
    command(""),
    diagnostic("")
{
}

BObolFeatureOwner::BObolFeatureOwner(void) :
    ownerToken(NULL),
    ownerId(""),
    ownerRole(""),
    generation(0),
    resultCallback(NULL),
    callbackUserData(NULL)
{
}

BObolOverlayInfo::BObolOverlayInfo(void) :
    isOverlay(FALSE),
    ownerToken(NULL),
    role(BObolOverlayRole::None),
    overlayClass(BObolOverlayClass::None),
    lifecycle(BObolOverlayLifecycle::None),
    order(BObolOverlayOrder::Model),
    sortOrder(0),
    sourcePath("")
{
}

BObolLabel::BObolLabel(void) :
    text(""),
    point(0.0f, 0.0f, 0.0f),
    hasColor(FALSE),
    color(1.0f, 1.0f, 1.0f),
    hasLeader(FALSE),
    target(0.0f, 0.0f, 0.0f),
    anchor(0),
    arrow(FALSE),
    fontSize(20.0f),
    sourceId(0)
{
}

BObolLineLayer::BObolLineLayer(void) :
    name(""),
    points(),
    commands(),
    style()
{
}

BObolFeatureMetadata::BObolFeatureMetadata(void) :
    key(""),
    value("")
{
}

BObolFeaturePrimitiveMetadata::BObolFeaturePrimitiveMetadata(void) :
    primitiveIndex(-1),
    metadata()
{
}

BObolFeaturePrimitivePick::BObolFeaturePrimitivePick(void) :
    handle(),
    featureName(""),
    primitiveIndex(-1),
    metadata()
{
}

BObolFeatureSummary::BObolFeatureSummary(void) :
    exists(FALSE),
    visible(FALSE),
    realized(FALSE),
    kind(BObolFeatureKind::Unknown),
    scope(BObolFeatureScope::Shared),
    pointCount(0),
    commandCount(0),
    childCount(0),
    metadataCount(0),
    primitiveMetadataCount(0),
    selectedPrimitiveCount(0),
    highlightedPrimitiveCount(0),
    owner(),
    overlay()
{
}

BObolFeatureRecord::BObolFeatureRecord(void) :
    handle(),
    name(""),
    kind(BObolFeatureKind::Unknown),
    scope(BObolFeatureScope::Shared),
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
    layers(),
    metadata(),
    primitiveMetadata(),
    selectedPrimitives(),
    highlightedPrimitives(),
    identity(""),
    editIntentId(""),
    editIntentRole(""),
    sourceRevision(0),
    inputsRevision(0)
{
}

BObolPolygonVisual::BObolPolygonVisual(void) :
    edgeColor(1.0f, 1.0f, 0.0f),
    fillColor(0.0f, 0.0f, 1.0f),
    fill(FALSE),
    fillFlags(BOBOL_POLYGON_FILL_NONE),
    fillSlope(1.0f, 0.0f),
    fillSpacing(1.0f),
    viewZ(0.0f)
{
}

BObolPolygonRecord::BObolPolygonRecord(void) :
    handle(),
    name(""),
    scope(BObolFeatureScope::Shared),
    type(BObolPolygonType::General),
    fill(FALSE),
    fillFlags(BOBOL_POLYGON_FILL_NONE),
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

BObolSelectionRecord::BObolSelectionRecord(void) :
    path(""),
    feature(),
    owner(),
    kind(BOBOL_SELECTION_SELECTED_PATH),
    primitiveIndex(-1),
    hitDistance(0.0)
{
}

struct BObolFeatureStoreRecord {
    uint64_t id;
    uint64_t revision;
    SbString name;
    BObolFeatureKind kind;
    BObolFeatureScope scope;
    BObolFeatureStyle style;
    BObolFeatureOwner owner;
    BObolOverlayInfo overlay;
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
    SbBool compactEdit;
    BObolCompactInstanceSummary compactSummary;
    SoNode *node;

    BObolFeatureStoreRecord(void) :
	id(0),
	revision(0),
	name(""),
	kind(BObolFeatureKind::Unknown),
	scope(BObolFeatureScope::Shared),
	style(),
	owner(),
	overlay(),
	points(),
	commands(),
	indices(),
	normals(),
	labels(),
	axesCenters(),
	halfAxesSize(1.0f),
	layers(),
	metadata(),
	primitiveMetadata(),
	selectedPrimitives(),
	highlightedPrimitives(),
	identity(""),
	editIntentId(""),
	editIntentRole(""),
	sourceRevision(0),
	inputsRevision(0),
	compactEdit(FALSE),
	compactSummary(),
	node(NULL)
    {
    }
};

static void
store_primitive_metadata_for_record(const BObolFeatureStoreRecord *rec,
				    int32_t primitiveIndex,
				    std::vector<BObolFeatureMetadata> &metadataOut)
{
    metadataOut.clear();
    if (!rec || primitiveIndex < 0)
	return;

    for (std::vector<BObolFeaturePrimitiveMetadata>::const_iterator it =
	     rec->primitiveMetadata.begin();
	 it != rec->primitiveMetadata.end(); ++it) {
	if (it->primitiveIndex != primitiveIndex)
	    continue;
	metadataOut = it->metadata;
	return;
    }
}

struct BObolFeatureStore::Impl {
    BObolViewController *controller;
    uint64_t referenceGeneration;
    uint64_t nextId;
    std::map<uint64_t, BObolFeatureStoreRecord *> records;
    std::map<std::string, uint64_t> names;
    std::map<std::string, uint64_t> ownerGenerations;

    Impl(void) : controller(NULL),
	referenceGeneration(store_reference_generation_next()), nextId(1),
	records(), names(), ownerGenerations()
    {
    }

    ~Impl(void)
    {
	clear();
    }

    void clear(void)
    {
	for (std::map<uint64_t, BObolFeatureStoreRecord *>::iterator it =
		 records.begin(); it != records.end(); ++it) {
	    if (it->second) {
		store_release_node(controller, it->second->node);
		delete it->second;
	    }
	}
	records.clear();
	names.clear();
	ownerGenerations.clear();
    }

    BObolFeatureStoreRecord *record(BObolFeatureHandle handle) const
    {
	std::map<uint64_t, BObolFeatureStoreRecord *>::const_iterator it =
	    records.find(handle.id);
	if (it == records.end() || !it->second)
	    return NULL;
	return it->second;
    }

    BObolFeatureStoreRecord *recordByName(
	const SbString &name,
	unsigned int scopeMask,
	const BObolFeatureOwner *owner = NULL) const
    {
	const std::string cleanName = store_string(name);
	if (cleanName.empty())
	    return NULL;

	for (std::map<uint64_t, BObolFeatureStoreRecord *>::const_iterator it =
		 records.begin(); it != records.end(); ++it) {
	    BObolFeatureStoreRecord *rec = it->second;
	    if (!rec)
		continue;
	    if (store_string(rec->name) != cleanName)
		continue;
	    if (!(store_scope_bit(rec->scope) & scopeMask))
		continue;
	    if (owner && !store_owner_matches(rec->owner, owner))
		continue;
	    return rec;
	}
	return NULL;
    }

    BObolFeatureStoreRecord *upsert(const SbString &name,
				      BObolFeatureScope scope,
				      BObolFeatureKind kind,
				      const BObolFeatureStyle *style,
				      const BObolFeatureOwner *owner)
    {
	if (store_string(name).empty())
	    return NULL;
	const std::string key = store_key(scope, name, owner);

	uint64_t id = 0;
	std::map<std::string, uint64_t>::iterator nit = names.find(key);
	if (nit != names.end())
	    id = nit->second;

	BObolFeatureStoreRecord *rec = NULL;
	if (id) {
	    std::map<uint64_t, BObolFeatureStoreRecord *>::iterator rit =
		records.find(id);
	    if (rit != records.end())
		rec = rit->second;
	}

	if (!rec) {
	    rec = new BObolFeatureStoreRecord;
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
	rec->metadata.clear();
	rec->primitiveMetadata.clear();
	rec->selectedPrimitives.clear();
	rec->highlightedPrimitives.clear();
	rec->compactEdit = FALSE;
	rec->compactSummary = BObolCompactInstanceSummary();
	if (style)
	    rec->style = *style;
	if (owner)
	    rec->owner = *owner;
	return rec;
    }

    BObolFeatureHandle handle(const BObolFeatureStoreRecord *rec) const
    {
	return rec ? BObolFeatureHandle(rec->id, rec->revision) :
	       BObolFeatureHandle();
    }

    void notify(const BObolFeatureStoreRecord *rec,
		BObolCommandResultStatus status,
		const char *command,
		const char *diagnostic = NULL) const
    {
	if (!rec || !rec->owner.resultCallback)
	    return;

	BObolCommandResult result;
	result.status = status;
	result.feature = handle(rec);
	result.command = command ? command : "";
	result.diagnostic = diagnostic ? diagnostic : "";
	rec->owner.resultCallback(result, rec->owner.callbackUserData);
    }

    void setNode(BObolFeatureStoreRecord *rec, SoNode *node)
    {
	if (!rec)
	    return;
	if (rec->node == node)
	    return;
	store_set_node(controller, rec->node, node);
	if (controller)
	    controller->requestRender("view-feature-store");
    }

    void markOwnerGeneration(const BObolFeatureOwner &owner)
    {
	if (owner.generation == 0)
	    return;

	const std::string key = store_owner_generation_key(&owner);
	if (key.empty())
	    return;

	uint64_t &generation = ownerGenerations[key];
	if (owner.generation > generation)
	    generation = owner.generation;
    }

    SbBool ownerGenerationCurrent(const BObolFeatureOwner &owner) const
    {
	if (owner.generation == 0)
	    return TRUE;

	const std::string key = store_owner_generation_key(&owner);
	if (key.empty())
	    return TRUE;

	std::map<std::string, uint64_t>::const_iterator it =
	    ownerGenerations.find(key);
	return it == ownerGenerations.end() ||
	       owner.generation >= it->second ? TRUE : FALSE;
    }
};

static void
store_apply_vlist_primitive_field(SoMFInt32 &field,
				  const std::vector<int32_t> &primitives)
{
    field.setNum(0);
    if (!primitives.empty())
	field.setValues(0, static_cast<int>(primitives.size()),
			&primitives[0]);
}

static void
store_apply_vlist_primitives(SoBRLVListShape *shape,
			     const BObolFeatureStoreRecord &rec)
{
    if (!shape)
	return;

    store_apply_vlist_primitive_field(shape->selectedPrimitive,
				      rec.selectedPrimitives);
    store_apply_vlist_primitive_field(shape->highlightedPrimitive,
				      rec.highlightedPrimitives);
}

static size_t
store_line_segment_primitive_count(const std::vector<int32_t> &commands)
{
    size_t count = 0;
    SbBool haveLast = FALSE;
    for (size_t i = 0; i < commands.size(); i++) {
	switch (commands[i]) {
	    case static_cast<int32_t>(BObolLineCommand::Move):
		haveLast = TRUE;
		break;
	    case static_cast<int32_t>(BObolLineCommand::Draw):
		if (haveLast)
		    count++;
		haveLast = TRUE;
		break;
	    default:
		break;
	}
    }
    return count;
}

static std::vector<int32_t>
store_primitive_subset_for_layer(const std::vector<int32_t> &primitives,
				 size_t firstPrimitive,
				 size_t primitiveCount)
{
    std::vector<int32_t> out;
    const size_t end = firstPrimitive + primitiveCount;
    for (size_t i = 0; i < primitives.size(); i++) {
	if (primitives[i] < 0)
	    continue;
	const size_t primitive = static_cast<size_t>(primitives[i]);
	if (primitive >= firstPrimitive && primitive < end)
	    out.push_back(static_cast<int32_t>(primitive - firstPrimitive));
    }
    return out;
}

static SoBRLVListShape *
store_vlist_node(const BObolFeatureStoreRecord &rec)
{
    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->sourcePath = rec.identity.getLength() > 0 ? rec.identity : rec.name;
    shape->sourceName = rec.name;
    switch (rec.kind) {
	case BObolFeatureKind::Points:
	    shape->sourceType = "point-set";
	    shape->geometryKind = "point";
	    break;
	case BObolFeatureKind::Arrow:
	    shape->sourceType = "arrow";
	    shape->geometryKind = "line";
	    break;
	case BObolFeatureKind::IndexedLines:
	    shape->sourceType = "indexed-line-set";
	    shape->geometryKind = "line";
	    break;
	case BObolFeatureKind::EditPreview:
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
    /* HUD features are screen-locked: their pixel-space geometry is projected
     * by an SoHUDKit (see store_hud_wrap_if_needed).  Non-HUD line features keep
     * the normal model/view pipeline. */
    shape->hudIntent = rec.style.hud ? TRUE : FALSE;
    shape->localSource = rec.scope == BObolFeatureScope::Local ? TRUE : FALSE;
    shape->sharedSource = rec.scope == BObolFeatureScope::Shared ? TRUE : FALSE;
    shape->nonDatabaseSource = TRUE;
    shape->drawMode = BOBOL_LOD_DRAW_DIAGNOSTIC;
    shape->recordRole = "view-feature";
    shape->sourceId = rec.sourceRevision;
    store_apply_vlist_style(shape, rec.style);
    store_apply_vlist_primitives(shape, rec);

    std::vector<int32_t> commands =
	store_normalized_line_commands(rec.points, rec.commands);
    std::vector<int32_t> shapeCommands = store_shape_commands(commands);

    if (!rec.points.empty())
	shape->setLineSet(&rec.points[0], &shapeCommands[0],
			  static_cast<int>(rec.points.size()));
    return shape;
}

/* Wrap a HUD (screen-locked) line shape in an SoHUDKit so its pixel-space
 * geometry is projected to the screen, matching the HUD text labels and the
 * grid overlay (see grid.cpp).  A non-HUD shape is returned unchanged so the
 * normal model/view pipeline still applies. */
static SoNode *
store_hud_wrap_if_needed(SoBRLVListShape *shape)
{
    if (!shape || !shape->hudIntent.getValue())
	return shape;
    SoHUDKit *hud = new SoHUDKit;
    hud->addWidget(shape);
    return hud;
}

static SoBRLEditPreview *
store_edit_preview_node(const BObolFeatureStoreRecord &rec)
{
    SoBRLEditPreview *preview = new SoBRLEditPreview;
    preview->previewId = rec.name;
    preview->setEditIntent(rec.editIntentId, rec.editIntentRole);
    preview->sourceRevision = rec.sourceRevision;
    preview->inputsRevision = rec.inputsRevision;
    if (!rec.points.empty()) {
	std::vector<int32_t> shapeCommands =
	    store_shape_commands(rec.commands);
	SoBRLVListShape *shape = preview->setLineSet(
	    rec.identity.getLength() > 0 ? rec.identity : rec.name,
	    &rec.points[0],
	    &shapeCommands[0],
	    static_cast<int>(rec.points.size()));
	store_apply_vlist_style(shape, rec.style);
    }
    return preview;
}

static SoBRLSceneGroup *
store_feature_group_node(const BObolFeatureStoreRecord &rec)
{
    SoBRLSceneGroup *group = new SoBRLSceneGroup;
    const SbString path =
	rec.identity.getLength() > 0 ? rec.identity : rec.name;
    group->groupPath = path;
    group->drawIntentValid = TRUE;
    group->drawIntentPath = path;
    group->drawMode = BOBOL_LOD_DRAW_DIAGNOSTIC;
    group->fallbackDrawMode = BOBOL_LOD_DRAW_WIRE;
    group->overlayIntent = TRUE;
    if (rec.style.hasVisible)
	group->visible = rec.style.visible;
    if (rec.style.hasLineStyle)
	group->lineStyle = rec.style.lineStyle;
    if (rec.style.hasLineWidth)
	group->lineWidth = rec.style.lineWidth;
    if (rec.style.hasTransparency)
	group->transparency = rec.style.transparency;
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
    [](int32_t idx) {
	return idx < 0;
    }) != indices.end() ?
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
store_indexed_face_node(const BObolFeatureStoreRecord &rec)
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
    shape->localSource = rec.scope == BObolFeatureScope::Local ? TRUE : FALSE;
    shape->sharedSource = rec.scope == BObolFeatureScope::Shared ? TRUE : FALSE;
    shape->nonDatabaseSource = TRUE;
    shape->drawMode = BOBOL_LOD_DRAW_SHADED;
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
store_line_layers_node(const BObolFeatureStoreRecord &rec)
{
    SoBRLSceneGroup *sep = store_feature_group_node(rec);
    size_t primitiveOffset = 0;
    for (size_t i = 0; i < rec.layers.size(); i++) {
	BObolFeatureStoreRecord layerRec = rec;
	layerRec.name = rec.layers[i].name.getLength() > 0 ?
			rec.layers[i].name : rec.name;
	layerRec.identity = layerRec.name;
	layerRec.points = rec.layers[i].points;
	layerRec.commands = rec.layers[i].commands;
	layerRec.style = store_merge_feature_style(rec.style,
			 rec.layers[i].style);
	const size_t primitiveCount =
	    store_line_segment_primitive_count(layerRec.commands);
	layerRec.selectedPrimitives = store_primitive_subset_for_layer(
					  rec.selectedPrimitives, primitiveOffset, primitiveCount);
	layerRec.highlightedPrimitives = store_primitive_subset_for_layer(
					     rec.highlightedPrimitives, primitiveOffset, primitiveCount);
	primitiveOffset += primitiveCount;
	layerRec.kind = BObolFeatureKind::Lines;
	SoNode *layerNode = store_hud_wrap_if_needed(store_vlist_node(layerRec));
	if (layerNode)
	    sep->addChild(layerNode);
    }
    return sep;
}

static SoNode *
store_axes_node(const BObolFeatureStoreRecord &rec)
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
store_label_node(const BObolFeatureStoreRecord &rec)
{
    SoBRLSceneGroup *sep = store_feature_group_node(rec);
    const SbColor fallbackColor = rec.style.hasColor ?
				  rec.style.color : SbColor(1.0f, 1.0f, 1.0f);

    for (size_t i = 0; i < rec.labels.size(); i++) {
	const BObolLabel &label = rec.labels[i];
	const SbColor color = label.hasColor ? label.color : fallbackColor;

	if (label.hasLeader) {
	    BObolFeatureStoreRecord leader = rec;
	    leader.kind = BObolFeatureKind::Lines;
	    leader.points.clear();
	    leader.commands.clear();
	    leader.points.push_back(label.target);
	    leader.points.push_back(label.point);
	    leader.commands.push_back(static_cast<int32_t>(
					  BObolLineCommand::Move));
	    leader.commands.push_back(static_cast<int32_t>(
					  BObolLineCommand::Draw));
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
store_hud_label_node(const BObolFeatureStoreRecord &rec)
{
    SoBRLSceneGroup *sep = store_feature_group_node(rec);
    const SbColor fallbackColor = rec.style.hasColor ?
				  rec.style.color : SbColor(1.0f, 1.0f, 1.0f);
    const SbBool visible = rec.style.hasVisible ? rec.style.visible : TRUE;

    for (size_t i = 0; i < rec.labels.size(); i++) {
	const BObolLabel &label = rec.labels[i];
	if (label.text.getLength() == 0)
	    continue;

	const SbColor color = label.hasColor ? label.color : fallbackColor;
	SoBRLHUDLabelOverlay *overlay = new SoBRLHUDLabelOverlay;
	overlay->labelId = rec.name;
	overlay->sourceId = label.sourceId;
	overlay->text = label.text;
	overlay->position = SbVec2f(label.point[0], label.point[1]);
	overlay->color = color;
	overlay->fontSize = label.fontSize > 0.0f ? label.fontSize : 12.0f;
	overlay->visible = visible;
	overlay->rebuildGeometry();
	sep->addChild(overlay);
    }
    return sep;
}

static SoNode *
store_node_for_feature(const BObolFeatureStoreRecord &rec)
{
    switch (rec.kind) {
	case BObolFeatureKind::Labels:
	    return store_label_node(rec);
	case BObolFeatureKind::HudLabel:
	    return store_hud_label_node(rec);
	case BObolFeatureKind::Axes:
	    return store_axes_node(rec);
	case BObolFeatureKind::LineLayer:
	    return store_line_layers_node(rec);
	case BObolFeatureKind::IndexedFaceSet:
	    return store_indexed_face_node(rec);
	case BObolFeatureKind::EditPreview:
	    return store_edit_preview_node(rec);
	case BObolFeatureKind::CustomNode:
	    return rec.node;
	default:
	    return store_hud_wrap_if_needed(store_vlist_node(rec));
    }
}

static SoNode *
store_rebuild_node_for_feature(const BObolFeatureStoreRecord &rec)
{
    if (rec.kind == BObolFeatureKind::CustomNode)
	return rec.node;
    return store_node_for_feature(rec);
}

BObolFeatureStore::BObolFeatureStore(void) : impl(new Impl)
{
}

BObolFeatureStore::BObolFeatureStore(BObolViewController *controller) :
    impl(new Impl)
{
    this->impl->controller = controller;
}

BObolFeatureStore::~BObolFeatureStore(void)
{
    delete this->impl;
    this->impl = NULL;
}

void
BObolFeatureStore::setController(BObolViewController *controller)
{
    this->impl->controller = controller;
}

BObolViewController *
BObolFeatureStore::controller(void) const
{
    return this->impl->controller;
}

uint64_t
BObolFeatureStore::referenceGeneration(void) const
{
    return this->impl->referenceGeneration;
}

void
BObolFeatureStore::clear(void)
{
    this->impl->clear();
}

BObolFeatureHandle
BObolFeatureStore::find(const SbString &name, unsigned int scopeMask) const
{
    return this->impl->handle(this->impl->recordByName(name, scopeMask));
}

BObolFeatureHandle
BObolFeatureStore::findOwned(const SbString &name,
			       unsigned int scopeMask,
			       const BObolFeatureOwner *owner) const
{
    return this->impl->handle(this->impl->recordByName(name, scopeMask,
			      owner));
}

SbBool
BObolFeatureStore::exists(const SbString &name, unsigned int scopeMask) const
{
    return this->find(name, scopeMask).isValid();
}

SbBool
BObolFeatureStore::existsOwned(const SbString &name,
				 unsigned int scopeMask,
				 const BObolFeatureOwner *owner) const
{
    return this->findOwned(name, scopeMask, owner).isValid();
}

SbBool
BObolFeatureStore::remove(BObolFeatureHandle handle)
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    const std::string key = store_key(rec->scope, rec->name, &rec->owner);
    this->impl->names.erase(key);
    this->impl->notify(rec, BObolCommandResultStatus::Removed, "remove");
    store_release_node(this->impl->controller, rec->node);
    this->impl->records.erase(rec->id);
    delete rec;
    if (this->impl->controller)
	this->impl->controller->requestRender("view-feature-remove");
    return TRUE;
}

SbBool
BObolFeatureStore::remove(const SbString &name)
{
    return this->remove(this->find(name));
}

SbBool
BObolFeatureStore::removeOwned(const SbString &name,
				 unsigned int scopeMask,
				 const BObolFeatureOwner *owner)
{
    return this->remove(this->findOwned(name, scopeMask, owner));
}

size_t
BObolFeatureStore::removeScope(unsigned int scopeMask,
				 const BObolFeatureOwner *owner)
{
    std::vector<BObolFeatureHandle> handles;
    for (std::map<uint64_t, BObolFeatureStoreRecord *>::const_iterator it =
	     this->impl->records.begin(); it != this->impl->records.end(); ++it) {
	if (!it->second)
	    continue;
	if (!(store_scope_bit(it->second->scope) & scopeMask))
	    continue;
	if (owner && !store_owner_matches(it->second->owner, owner))
	    continue;
	handles.push_back(this->impl->handle(it->second));
    }

    size_t removed = 0;
    for (size_t i = 0; i < handles.size(); i++)
	if (this->remove(handles[i]))
	    removed++;
    return removed;
}

size_t
BObolFeatureStore::removePrefix(const SbString &prefix)
{
    return this->removePrefix(prefix, BOBOL_FEATURE_SCOPE_ALL, NULL);
}

size_t
BObolFeatureStore::removePrefix(const SbString &prefix,
				  unsigned int scopeMask,
				  const BObolFeatureOwner *owner)
{
    const std::string p = store_string(prefix);
    if (p.empty())
	return 0;

    std::vector<BObolFeatureHandle> handles;
    for (std::map<uint64_t, BObolFeatureStoreRecord *>::const_iterator it =
	     this->impl->records.begin(); it != this->impl->records.end(); ++it) {
	if (!it->second)
	    continue;
	if (!(store_scope_bit(it->second->scope) & scopeMask))
	    continue;
	if (owner && !store_owner_matches(it->second->owner, owner))
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

size_t
BObolFeatureStore::countPrefix(const SbString &prefix,
			       unsigned int scopeMask,
			       const BObolFeatureOwner *owner) const
{
    const std::string p = store_string(prefix);
    if (p.empty())
	return 0;

    size_t count = 0;
    for (std::map<uint64_t, BObolFeatureStoreRecord *>::const_iterator it =
	 this->impl->records.begin(); it != this->impl->records.end(); ++it) {
	if (!it->second ||
	    !(store_scope_bit(it->second->scope) & scopeMask) ||
	    (owner && !store_owner_matches(it->second->owner, owner)))
	    continue;
	const std::string name = store_string(it->second->name);
	if (name.compare(0, p.size(), p) == 0)
	    count++;
    }
    return count;
}

void
BObolFeatureStore::markCommandOwnerGeneration(
    const BObolFeatureOwner &owner)
{
    this->impl->markOwnerGeneration(owner);
}

SbBool
BObolFeatureStore::commandOwnerGenerationCurrent(
    const BObolFeatureOwner &owner) const
{
    return this->impl->ownerGenerationCurrent(owner);
}

BObolFeatureHandle
BObolFeatureStore::publishLineSet(const SbString &name,
				    BObolFeatureScope scope,
				    const std::vector<SbVec3f> &points,
				    const std::vector<int32_t> &commands,
				    const BObolFeatureStyle *style,
				    const BObolFeatureOwner *owner)
{
    BObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
				     BObolFeatureKind::Lines, style, owner);
    if (!rec)
	return BObolFeatureHandle();

    rec->points = points;
    rec->commands = commands;
    SoNode *node = store_rebuild_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated, "publishLineSet");
    return this->impl->handle(rec);
}

BObolFeatureHandle
BObolFeatureStore::publishIndexedLineSet(const SbString &name,
	BObolFeatureScope scope,
	const std::vector<SbVec3f> &points,
	const std::vector<int32_t> &indices,
	const BObolFeatureStyle *style,
	const BObolFeatureOwner *owner)
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
			   static_cast<int32_t>(BObolLineCommand::Move) :
			   static_cast<int32_t>(BObolLineCommand::Draw));
    }

    BObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
				     BObolFeatureKind::IndexedLines, style, owner);
    if (!rec)
	return BObolFeatureHandle();

    rec->points = linePoints;
    rec->commands = commands;
    rec->indices = indices;
    SoNode *node = store_rebuild_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "publishIndexedLineSet");
    return this->impl->handle(rec);
}

BObolFeatureHandle
BObolFeatureStore::publishPointSet(const SbString &name,
				     BObolFeatureScope scope,
				     const std::vector<SbVec3f> &points,
				     const BObolFeatureStyle *style,
				     const BObolFeatureOwner *owner)
{
    std::vector<int32_t> commands(points.size(),
				  static_cast<int32_t>(BObolLineCommand::Point));
    BObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
				     BObolFeatureKind::Points, style, owner);
    if (!rec)
	return BObolFeatureHandle();
    rec->points = points;
    rec->commands = commands;
    SoNode *node = store_rebuild_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "publishPointSet");
    return this->impl->handle(rec);
}

BObolFeatureHandle
BObolFeatureStore::publishLabels(const SbString &name,
				   BObolFeatureScope scope,
				   const std::vector<BObolLabel> &labels,
				   const BObolFeatureStyle *style,
				   const BObolFeatureOwner *owner)
{
    BObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
				     BObolFeatureKind::Labels, style, owner);
    if (!rec)
	return BObolFeatureHandle();
    rec->labels = labels;
    SoNode *node = store_rebuild_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "publishLabels");
    return this->impl->handle(rec);
}

BObolFeatureHandle
BObolFeatureStore::publishHudLabels(const SbString &name,
				      BObolFeatureScope scope,
				      const std::vector<BObolLabel> &labels,
				      const BObolFeatureStyle *style,
				      const BObolFeatureOwner *owner)
{
    BObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
				     BObolFeatureKind::HudLabel, style, owner);
    if (!rec)
	return BObolFeatureHandle();
    rec->labels = labels;
    SoNode *node = store_rebuild_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "publishHudLabels");
    return this->impl->handle(rec);
}

BObolFeatureHandle
BObolFeatureStore::publishArrow(const SbString &name,
				  BObolFeatureScope scope,
				  const std::vector<SbVec3f> &points,
				  const BObolFeatureStyle *style,
				  const BObolFeatureOwner *owner)
{
    BObolFeatureStyle arrowStyle = style ? *style : BObolFeatureStyle();
    arrowStyle.hasArrow = TRUE;
    arrowStyle.arrow = TRUE;

    BObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
				     BObolFeatureKind::Arrow, &arrowStyle, owner);
    if (!rec)
	return BObolFeatureHandle();

    rec->points = points;
    rec->commands.clear();
    SoNode *node = store_rebuild_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "publishArrow");
    return this->impl->handle(rec);
}

BObolFeatureHandle
BObolFeatureStore::publishAxes(const SbString &name,
				 BObolFeatureScope scope,
				 const std::vector<SbVec3f> &centers,
				 float halfAxesSize,
				 const BObolFeatureStyle *style,
				 const BObolFeatureOwner *owner)
{
    BObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
				     BObolFeatureKind::Axes, style, owner);
    if (!rec)
	return BObolFeatureHandle();
    rec->axesCenters = centers;
    rec->halfAxesSize = halfAxesSize;
    SoNode *node = store_rebuild_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "publishAxes");
    return this->impl->handle(rec);
}

BObolFeatureHandle
BObolFeatureStore::publishLineLayers(const SbString &name,
				       BObolFeatureScope scope,
				       const std::vector<BObolLineLayer> &layers,
				       const BObolFeatureStyle *style,
				       const BObolFeatureOwner *owner)
{
    BObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
				     BObolFeatureKind::LineLayer, style, owner);
    if (!rec)
	return BObolFeatureHandle();
    rec->layers = layers;
    rec->points.clear();
    rec->commands.clear();
    for (size_t i = 0; i < layers.size(); i++) {
	rec->points.insert(rec->points.end(), layers[i].points.begin(),
			   layers[i].points.end());
	rec->commands.insert(rec->commands.end(), layers[i].commands.begin(),
			     layers[i].commands.end());
    }
    SoNode *node = store_rebuild_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "publishLineLayers");
    return this->impl->handle(rec);
}

BObolFeatureHandle
BObolFeatureStore::publishLineLayerBuilder(const SbString &name,
	BObolFeatureScope scope,
	const struct bg_line_layer_builder *builder,
	const BObolFeatureStyle *style,
	const BObolFeatureOwner *owner)
{
    if (!builder)
	return BObolFeatureHandle();

    BObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
				     BObolFeatureKind::LineLayer, style, owner);
    if (!rec)
	return BObolFeatureHandle();

    rec->layers = store_line_layers_from_builder(name, builder);
    rec->points.clear();
    rec->commands.clear();
    for (size_t i = 0; i < rec->layers.size(); i++) {
	rec->points.insert(rec->points.end(), rec->layers[i].points.begin(),
			   rec->layers[i].points.end());
	rec->commands.insert(rec->commands.end(),
			     rec->layers[i].commands.begin(), rec->layers[i].commands.end());
    }

    SoBRLLineLayerOverlay *overlay = new SoBRLLineLayerOverlay;
    overlay->overlayId = name;
    overlay->sourceId = static_cast<uint32_t>(rec->revision);
    overlay->selectable = rec->style.hasSelectable ?
			  rec->style.selectable : TRUE;
    overlay->rebuildGeometry(builder);
    this->impl->setNode(rec, overlay);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "publishLineLayerBuilder");
    return this->impl->handle(rec);
}

BObolFeatureHandle
BObolFeatureStore::publishIndexedFaceSet(const SbString &name,
	BObolFeatureScope scope,
	const std::vector<SbVec3f> &points,
	const std::vector<SbVec3f> &normals,
	const std::vector<int32_t> &indices,
	const BObolFeatureStyle *style,
	const BObolFeatureOwner *owner)
{
    BObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
				     BObolFeatureKind::IndexedFaceSet, style, owner);
    if (!rec)
	return BObolFeatureHandle();

    rec->points = points;
    rec->normals = normals;
    rec->indices = indices;
    rec->commands.clear();
    SoNode *node = store_rebuild_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "publishIndexedFaceSet");
    return this->impl->handle(rec);
}

BObolFeatureHandle
BObolFeatureStore::publishCustomNode(const SbString &name,
				       BObolFeatureScope scope,
				       SoNode *node,
				       const BObolFeatureStyle *style,
				       const BObolFeatureOwner *owner)
{
    if (!node)
	return BObolFeatureHandle();

    BObolFeatureStoreRecord *rec = this->impl->upsert(name, scope,
				     BObolFeatureKind::CustomNode, style, owner);
    if (!rec)
	return BObolFeatureHandle();

    rec->points.clear();
    rec->commands.clear();
    rec->indices.clear();
    rec->normals.clear();
    rec->labels.clear();
    rec->axesCenters.clear();
    rec->layers.clear();
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "publishCustomNode");
    return this->impl->handle(rec);
}

BObolFeatureHandle
BObolFeatureStore::publishEditPreview(const SbString &name,
					const SbString &identity,
					const SbString &editIntentId,
					const SbString &editIntentRole,
					const std::vector<SbVec3f> &points,
					const std::vector<int32_t> &commands,
					uint32_t sourceRevision,
					uint32_t inputsRevision,
					const BObolFeatureOwner *owner)
{
    BObolFeatureStoreRecord *rec = this->impl->upsert(name,
				     BObolFeatureScope::Local, BObolFeatureKind::EditPreview,
				     NULL, owner);
    if (!rec)
	return BObolFeatureHandle();

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
    SoNode *node = store_rebuild_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "publishEditPreview");
    return this->impl->handle(rec);
}

static SbBool
store_compact_handle_equal(const BObolCompactInstanceHandle &a,
	const BObolCompactInstanceHandle &b)
{
    return a.sourceNodeId == b.sourceNodeId &&
	a.instanceWord0 == b.instanceWord0 &&
	a.instanceWord1 == b.instanceWord1 ? TRUE : FALSE;
}

BObolFeatureHandle
BObolFeatureStore::promoteCompactInstanceForEdit(
    const SbString &name,
    const SoBRLDatabaseSource &source,
    const BObolCompactInstanceHandle &instance,
    const SbString &editIntentId,
    const SbString &editIntentRole,
    const BObolFeatureOwner *owner)
{
    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    BObolCompactInstanceSummary summary;
    if (!source.copyCompactInstanceEditGeometry(instance, points, commands,
	summary))
	return BObolFeatureHandle();

    BObolFeatureStyle style;
    style.hasVisible = TRUE;
    style.visible = summary.visible;
    style.hasSelectable = TRUE;
    style.selectable = summary.selectable;
    style.hasColor = summary.appearanceColorValid;
    style.color = summary.appearanceColor;
    style.hasLineWidth = TRUE;
    style.lineWidth = summary.lineWidth;
    style.hasLineStyle = TRUE;
    style.lineStyle = summary.lineStyle;
    style.hasTransparency = TRUE;
    style.transparency = summary.transparency;

    BObolFeatureStoreRecord *rec = this->impl->upsert(name,
	BObolFeatureScope::Local, BObolFeatureKind::EditPreview, &style, owner);
    if (!rec)
	return BObolFeatureHandle();

    rec->identity = summary.sourceInstanceKey.getLength() > 0 ?
	summary.sourceInstanceKey : summary.path;
    rec->editIntentId = editIntentId.getLength() > 0 ? editIntentId : name;
    rec->editIntentRole = editIntentRole.getLength() > 0 ?
	editIntentRole : SbString("compact-instance");
    rec->points = points;
    rec->commands = store_normalized_line_commands(points, commands);
    rec->sourceRevision = source.sourceRevision.getValue();
    rec->inputsRevision = source.inputsRevision.getValue();
    rec->compactEdit = TRUE;
    rec->compactSummary = summary;

    SoNode *node = store_rebuild_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Accepted,
	"promoteCompactInstanceForEdit");
    return this->impl->handle(rec);
}

SbBool
BObolFeatureStore::demoteCompactInstanceFromEdit(
    BObolFeatureHandle preview,
    const BObolCompactInstanceHandle &instance)
{
    BObolFeatureStoreRecord *rec = this->impl->record(preview);
    if (!rec || rec->kind != BObolFeatureKind::EditPreview ||
	!rec->compactEdit ||
	!store_compact_handle_equal(rec->compactSummary.handle, instance))
	return FALSE;
    return this->remove(this->impl->handle(rec));
}

SbBool
BObolFeatureStore::compactEditBinding(
    BObolFeatureHandle preview,
    BObolCompactInstanceHandle &instanceOut,
    BObolCompactInstanceSummary &summaryOut) const
{
    instanceOut = BObolCompactInstanceHandle();
    summaryOut = BObolCompactInstanceSummary();
    BObolFeatureStoreRecord *rec = this->impl->record(preview);
    if (!rec || rec->kind != BObolFeatureKind::EditPreview ||
	!rec->compactEdit)
	return FALSE;
    instanceOut = rec->compactSummary.handle;
    summaryOut = rec->compactSummary;
    return TRUE;
}

SbBool
BObolFeatureStore::replaceEditPreviewGeometry(
    BObolFeatureHandle handle,
    const SbString &identity,
    const std::vector<SbVec3f> &points,
    const std::vector<int32_t> &commands,
    uint32_t sourceRevision,
    uint32_t inputsRevision)
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec || rec->kind != BObolFeatureKind::EditPreview)
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

    SoNode *node = store_rebuild_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "replaceEditPreviewGeometry");
    return TRUE;
}

SbBool
BObolFeatureStore::touch(BObolFeatureHandle handle)
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    rec->revision++;
    if (rec->node)
	rec->node->touch();
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
	"touchFeature");
    return TRUE;
}

SbBool
BObolFeatureStore::appendLinePoint(BObolFeatureHandle handle,
				     const SbVec3f &point)
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec || rec->kind == BObolFeatureKind::CustomNode)
	return FALSE;

    rec->points.push_back(point);
    rec->commands.push_back(static_cast<int32_t>(BObolLineCommand::Draw));
    rec->revision++;
    SoNode *node = store_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "appendLinePoint");
    return TRUE;
}

SbBool
BObolFeatureStore::replaceLabels(BObolFeatureHandle handle,
				   const std::vector<BObolLabel> &labels)
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec || rec->kind == BObolFeatureKind::CustomNode)
	return FALSE;

    rec->labels = labels;
    rec->revision++;
    SoNode *node = store_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "replaceLabels");
    return TRUE;
}

SbBool
BObolFeatureStore::clearGeometry(BObolFeatureHandle handle)
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
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
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "clearGeometry");
    return TRUE;
}

SbBool
BObolFeatureStore::points(BObolFeatureHandle handle,
			    std::vector<SbVec3f> &pointsOut) const
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    pointsOut = rec->points;
    return TRUE;
}

SbBool
BObolFeatureStore::commands(BObolFeatureHandle handle,
			      std::vector<int32_t> &commandsOut) const
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    commandsOut = rec->commands;
    return TRUE;
}

SbBool
BObolFeatureStore::lineCommandAt(BObolFeatureHandle handle,
				   size_t index,
				   int32_t &commandOut) const
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec || index >= rec->commands.size())
	return FALSE;
    commandOut = rec->commands[index];
    return TRUE;
}

SbBool
BObolFeatureStore::labels(BObolFeatureHandle handle,
			    std::vector<BObolLabel> &labelsOut) const
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    labelsOut = rec->labels;
    return TRUE;
}

SbBool
BObolFeatureStore::axesCenters(BObolFeatureHandle handle,
				 std::vector<SbVec3f> &centersOut,
				 float *halfAxesSizeOut) const
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    centersOut = rec->axesCenters;
    if (halfAxesSizeOut)
	*halfAxesSizeOut = rec->halfAxesSize;
    return TRUE;
}

SbBool
BObolFeatureStore::indices(BObolFeatureHandle handle,
			     std::vector<int32_t> &indicesOut) const
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    indicesOut = rec->indices;
    return TRUE;
}

SbBool
BObolFeatureStore::normals(BObolFeatureHandle handle,
			     std::vector<SbVec3f> &normalsOut) const
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    normalsOut = rec->normals;
    return TRUE;
}

SbBool
BObolFeatureStore::applyStyle(BObolFeatureHandle handle,
				const BObolFeatureStyle &style,
				SbBool UNUSED(recursive))
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    if (style.hasVisible) {
	rec->style.hasVisible = TRUE;
	rec->style.visible = style.visible;
    }
    if (style.hasSelectable) {
	rec->style.hasSelectable = TRUE;
	rec->style.selectable = style.selectable;
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
    if (style.hasTransparency) {
	rec->style.hasTransparency = TRUE;
	rec->style.transparency = style.transparency;
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
    SoNode *node = store_rebuild_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "applyStyle");
    return TRUE;
}

SbBool
BObolFeatureStore::style(BObolFeatureHandle handle,
			   BObolFeatureStyle &styleOut) const
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    styleOut = rec->style;
    return TRUE;
}

SbBool
BObolFeatureStore::setVisible(BObolFeatureHandle handle, SbBool visible)
{
    BObolFeatureStyle style;
    style.hasVisible = TRUE;
    style.visible = visible;
    return this->applyStyle(handle, style);
}

SbBool
BObolFeatureStore::setColor(BObolFeatureHandle handle,
			      const SbColor &color)
{
    BObolFeatureStyle style;
    style.hasColor = TRUE;
    style.color = color;
    return this->applyStyle(handle, style);
}

SbBool
BObolFeatureStore::setLineWidth(BObolFeatureHandle handle,
				  int lineWidth)
{
    BObolFeatureStyle style;
    style.hasLineWidth = TRUE;
    style.lineWidth = lineWidth;
    return this->applyStyle(handle, style);
}

SbBool
BObolFeatureStore::arrowTip(BObolFeatureHandle handle,
			      float &tipLength,
			      float &tipWidth) const
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    tipLength = rec->style.arrowTipLength;
    tipWidth = rec->style.arrowTipWidth;
    return TRUE;
}

SbBool
BObolFeatureStore::setArrowTip(BObolFeatureHandle handle,
				 float tipLength,
				 float tipWidth)
{
    BObolFeatureStyle style;
    style.hasArrowTip = TRUE;
    style.arrowTipLength = tipLength;
    style.arrowTipWidth = tipWidth;
    return this->applyStyle(handle, style);
}

SbBool
BObolFeatureStore::setOverlayInfo(BObolFeatureHandle handle,
				    const BObolOverlayInfo &overlay)
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    rec->overlay = overlay;
    rec->revision++;
    SoNode *node = store_rebuild_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "setOverlayInfo");
    return TRUE;
}

SbBool
BObolFeatureStore::clearOverlayInfo(BObolFeatureHandle handle)
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    rec->overlay = BObolOverlayInfo();
    rec->revision++;
    SoNode *node = store_rebuild_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "clearOverlayInfo");
    return TRUE;
}

SbBool
BObolFeatureStore::overlayInfo(BObolFeatureHandle handle,
				 BObolOverlayInfo &overlayOut) const
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    overlayOut = rec->overlay;
    return TRUE;
}

SbBool
BObolFeatureStore::replaceMetadata(BObolFeatureHandle handle,
				     const std::vector<BObolFeatureMetadata> &metadata)
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    rec->metadata = metadata;
    rec->revision++;
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "replaceMetadata");
    return TRUE;
}

SbBool
BObolFeatureStore::metadata(BObolFeatureHandle handle,
			      std::vector<BObolFeatureMetadata> &metadataOut) const
{
    metadataOut.clear();
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    metadataOut = rec->metadata;
    return TRUE;
}

SbBool
BObolFeatureStore::replacePrimitiveMetadata(BObolFeatureHandle handle,
	int32_t primitiveIndex,
	const std::vector<BObolFeatureMetadata> &metadata)
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec || primitiveIndex < 0)
	return FALSE;

    for (std::vector<BObolFeaturePrimitiveMetadata>::iterator it =
	     rec->primitiveMetadata.begin();
	 it != rec->primitiveMetadata.end(); ++it) {
	if (it->primitiveIndex != primitiveIndex)
	    continue;
	if (metadata.empty())
	    rec->primitiveMetadata.erase(it);
	else
	    it->metadata = metadata;
	rec->revision++;
	this->impl->notify(rec, BObolCommandResultStatus::Updated,
			   "replacePrimitiveMetadata");
	return TRUE;
    }

    if (!metadata.empty()) {
	BObolFeaturePrimitiveMetadata item;
	item.primitiveIndex = primitiveIndex;
	item.metadata = metadata;
	rec->primitiveMetadata.push_back(item);
	rec->revision++;
	this->impl->notify(rec, BObolCommandResultStatus::Updated,
			   "replacePrimitiveMetadata");
    }
    return TRUE;
}

SbBool
BObolFeatureStore::primitiveMetadata(BObolFeatureHandle handle,
				       int32_t primitiveIndex,
				       std::vector<BObolFeatureMetadata> &metadataOut) const
{
    metadataOut.clear();
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec || primitiveIndex < 0)
	return FALSE;

    for (std::vector<BObolFeaturePrimitiveMetadata>::const_iterator it =
	     rec->primitiveMetadata.begin();
	 it != rec->primitiveMetadata.end(); ++it) {
	if (it->primitiveIndex != primitiveIndex)
	    continue;
	metadataOut = it->metadata;
	return TRUE;
    }

    return TRUE;
}

SbBool
BObolFeatureStore::resolvePrimitivePick(const SbString &name,
	int32_t primitiveIndex,
	BObolFeaturePrimitivePick &pickOut,
	unsigned int scopeMask,
	const BObolFeatureOwner *owner) const
{
    pickOut = BObolFeaturePrimitivePick();
    if (primitiveIndex < 0)
	return FALSE;

    BObolFeatureStoreRecord *rec = this->impl->recordByName(name,
				     scopeMask, owner);
    if (rec) {
	pickOut.handle = BObolFeatureHandle(rec->id, rec->revision);
	pickOut.featureName = rec->name;
	pickOut.primitiveIndex = primitiveIndex;
	store_primitive_metadata_for_record(rec, primitiveIndex,
					    pickOut.metadata);
	return TRUE;
    }

    const std::string cleanName = store_string(name);
    if (cleanName.empty())
	return FALSE;

    for (std::map<uint64_t, BObolFeatureStoreRecord *>::const_iterator it =
	     this->impl->records.begin(); it != this->impl->records.end();
	 ++it) {
	rec = it->second;
	if (!rec || rec->kind != BObolFeatureKind::LineLayer)
	    continue;
	if (!(store_scope_bit(rec->scope) & scopeMask))
	    continue;
	if (owner && !store_owner_matches(rec->owner, owner))
	    continue;

	size_t primitiveOffset = 0;
	for (size_t i = 0; i < rec->layers.size(); i++) {
	    const BObolLineLayer &layer = rec->layers[i];
	    const SbString layerName = layer.name.getLength() > 0 ?
				       layer.name : rec->name;
	    const size_t primitiveCount =
		store_line_segment_primitive_count(layer.commands);
	    if (store_string(layerName) != cleanName) {
		primitiveOffset += primitiveCount;
		continue;
	    }
	    if (static_cast<size_t>(primitiveIndex) >= primitiveCount)
		return FALSE;

	    const int32_t resolvedPrimitive = static_cast<int32_t>(
						  primitiveOffset + static_cast<size_t>(primitiveIndex));
	    pickOut.handle = BObolFeatureHandle(rec->id, rec->revision);
	    pickOut.featureName = rec->name;
	    pickOut.primitiveIndex = resolvedPrimitive;
	    store_primitive_metadata_for_record(rec, resolvedPrimitive,
						pickOut.metadata);
	    return TRUE;
	}
    }

    return FALSE;
}

SbBool
BObolFeatureStore::replaceSelectedPrimitives(BObolFeatureHandle handle,
	const std::vector<int32_t> &primitives)
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    rec->selectedPrimitives = primitives;
    rec->revision++;
    SoNode *node = store_rebuild_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "replaceSelectedPrimitives");
    return TRUE;
}

SbBool
BObolFeatureStore::replaceHighlightedPrimitives(BObolFeatureHandle handle,
	const std::vector<int32_t> &primitives)
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    rec->highlightedPrimitives = primitives;
    rec->revision++;
    SoNode *node = store_rebuild_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    this->impl->notify(rec, BObolCommandResultStatus::Updated,
		       "replaceHighlightedPrimitives");
    return TRUE;
}

SbBool
BObolFeatureStore::selectedPrimitives(BObolFeatureHandle handle,
					std::vector<int32_t> &primitivesOut) const
{
    primitivesOut.clear();
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    primitivesOut = rec->selectedPrimitives;
    return TRUE;
}

SbBool
BObolFeatureStore::highlightedPrimitives(BObolFeatureHandle handle,
	std::vector<int32_t> &primitivesOut) const
{
    primitivesOut.clear();
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    primitivesOut = rec->highlightedPrimitives;
    return TRUE;
}

SbBool
BObolFeatureStore::realize(BObolFeatureHandle handle,
			     SbBool UNUSED(recursive))
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    SoNode *node = store_rebuild_node_for_feature(*rec);
    this->impl->setNode(rec, node);
    return TRUE;
}

SbBool
BObolFeatureStore::summary(const SbString &name,
			     BObolFeatureSummary &summaryOut,
			     unsigned int scopeMask) const
{
    return this->summaryOwned(name, summaryOut, scopeMask, NULL);
}

SbBool
BObolFeatureStore::summaryOwned(const SbString &name,
				  BObolFeatureSummary &summaryOut,
				  unsigned int scopeMask,
				  const BObolFeatureOwner *owner) const
{
    summaryOut = BObolFeatureSummary();
    BObolFeatureStoreRecord *rec = this->impl->recordByName(name,
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
    summaryOut.metadataCount = rec->metadata.size();
    summaryOut.primitiveMetadataCount = rec->primitiveMetadata.size();
    summaryOut.selectedPrimitiveCount = rec->selectedPrimitives.size();
    summaryOut.highlightedPrimitiveCount = rec->highlightedPrimitives.size();
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
BObolFeatureStore::record(BObolFeatureHandle handle,
			    BObolFeatureRecord &recordOut) const
{
    recordOut = BObolFeatureRecord();
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
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
    recordOut.metadata = rec->metadata;
    recordOut.primitiveMetadata = rec->primitiveMetadata;
    recordOut.selectedPrimitives = rec->selectedPrimitives;
    recordOut.highlightedPrimitives = rec->highlightedPrimitives;
    recordOut.identity = rec->identity;
    recordOut.editIntentId = rec->editIntentId;
    recordOut.editIntentRole = rec->editIntentRole;
    recordOut.sourceRevision = rec->sourceRevision;
    recordOut.inputsRevision = rec->inputsRevision;
    return TRUE;
}

void
BObolFeatureStore::visitRecords(BObolFeatureRecordCallback callback,
				  void *userData,
				  unsigned int scopeMask,
				  const BObolFeatureOwner *owner) const
{
    if (!callback)
	return;

    for (std::map<uint64_t, BObolFeatureStoreRecord *>::const_iterator it =
	     this->impl->records.begin(); it != this->impl->records.end(); ++it) {
	BObolFeatureStoreRecord *rec = it->second;
	if (!rec)
	    continue;
	if (!(store_scope_bit(rec->scope) & scopeMask))
	    continue;
	if (rec->scope == BObolFeatureScope::Local &&
	    !store_owner_matches(rec->owner, owner))
	    continue;

	BObolFeatureRecord record;
	if (!this->record(this->impl->handle(rec), record))
	    continue;
	if (!callback(record, userData))
	    return;
    }
}

SoNode *
BObolFeatureStore::node(BObolFeatureHandle handle) const
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    return rec ? rec->node : NULL;
}

struct BObolPolygonStoreRecord {
    uint64_t id;
    uint64_t revision;
    SbString name;
    BObolFeatureScope scope;
    BObolPolygonType type;
    SbBool visible;
    BObolPolygonVisual visual;
    long currentContour;
    long currentPoint;
    SbVec3f originPoint;
    plane_t viewPlane;
    void *userData;
    struct bg_polygon polygon;
    SoNode *node;

    BObolPolygonStoreRecord(void) :
	id(0),
	revision(0),
	name(""),
	scope(BObolFeatureScope::Shared),
	type(BObolPolygonType::General),
	visible(TRUE),
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
store_polygon_node(const BObolPolygonStoreRecord &rec);

struct BObolPolygonStore::Impl {
    BObolViewController *controller;
    uint64_t referenceGeneration;
    uint64_t nextId;
    std::map<uint64_t, BObolPolygonStoreRecord *> records;
    std::map<std::string, uint64_t> names;
    BObolPolygonHandle snapExclude;

    Impl(void) : controller(NULL),
	referenceGeneration(store_reference_generation_next()), nextId(1),
	records(), names(), snapExclude()
    {
    }

    ~Impl(void)
    {
	clear();
    }

    void clear(void)
    {
	for (std::map<uint64_t, BObolPolygonStoreRecord *>::iterator it =
		 records.begin(); it != records.end(); ++it) {
	    if (it->second) {
		store_release_node(controller, it->second->node);
		bg_polygon_clear(&it->second->polygon);
		delete it->second;
	    }
	}
	records.clear();
	names.clear();
	snapExclude = BObolPolygonHandle();
    }

    BObolPolygonStoreRecord *record(BObolPolygonHandle handle) const
    {
	std::map<uint64_t, BObolPolygonStoreRecord *>::const_iterator it =
	    records.find(handle.id);
	if (it == records.end() || !it->second)
	    return NULL;
	return it->second;
    }

    BObolPolygonStoreRecord *recordByName(const SbString &name, unsigned int scopeMask) const
    {
	const std::string cleanName = store_string(name);
	if (cleanName.empty())
	    return NULL;

	for (std::map<uint64_t, BObolPolygonStoreRecord *>::const_iterator it =
		 records.begin(); it != records.end(); ++it) {
	    BObolPolygonStoreRecord *rec = it->second;
	    if (!rec)
		continue;
	    if (store_string(rec->name) != cleanName)
		continue;
	    if (!(store_scope_bit(rec->scope) & scopeMask))
		continue;
	    return rec;
	}
	return NULL;
    }

    BObolPolygonHandle handle(const BObolPolygonStoreRecord *rec) const
    {
	return rec ? BObolPolygonHandle(rec->id, rec->revision) :
	       BObolPolygonHandle();
    }

    void setNode(BObolPolygonStoreRecord *rec, SoNode *node)
    {
	if (!rec)
	    return;
	if (rec->node == node)
	    return;
	store_set_node(controller, rec->node, node);
	if (controller)
	    controller->requestRender("view-polygon-store");
    }

    void realize(BObolPolygonStoreRecord *rec)
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
    bg_polygon_clear(poly);
    poly->num_contours = 1;
    poly->hole = (int *)bu_calloc(1, sizeof(int), "BObol polygon hole");
    poly->contour = (struct bg_poly_contour *)bu_calloc(1,
		    sizeof(struct bg_poly_contour), "BObol polygon contour");
    poly->contour[0].num_points = 1;
    poly->contour[0].open = 1;
    poly->contour[0].point = (point_t *)bu_calloc(1, sizeof(point_t),
			     "BObol polygon point");
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
store_polygon_type_to_rt(BObolPolygonType type)
{
    switch (type) {
	case BObolPolygonType::Circle:
	    return RT_SKETCH_POLYGON_CIRCLE;
	case BObolPolygonType::Ellipse:
	    return RT_SKETCH_POLYGON_ELLIPSE;
	case BObolPolygonType::Rectangle:
	    return RT_SKETCH_POLYGON_RECTANGLE;
	case BObolPolygonType::Square:
	    return RT_SKETCH_POLYGON_SQUARE;
	default:
	    return RT_SKETCH_POLYGON_GENERAL;
    }
}

static BObolPolygonType
store_polygon_type_from_rt(int type)
{
    switch (type) {
	case RT_SKETCH_POLYGON_CIRCLE:
	    return BObolPolygonType::Circle;
	case RT_SKETCH_POLYGON_ELLIPSE:
	    return BObolPolygonType::Ellipse;
	case RT_SKETCH_POLYGON_RECTANGLE:
	    return BObolPolygonType::Rectangle;
	case RT_SKETCH_POLYGON_SQUARE:
	    return BObolPolygonType::Square;
	default:
	    return BObolPolygonType::General;
    }
}

static unsigned int
store_polygon_valid_fill_flags(unsigned int flags)
{
    return flags & (BOBOL_POLYGON_FILL_HATCH | BOBOL_POLYGON_FILL_MESH);
}

static unsigned int
store_polygon_fill_flags(const BObolPolygonVisual &visual)
{
    const unsigned int flags =
	store_polygon_valid_fill_flags(visual.fillFlags);
    if (flags)
	return flags;
    return visual.fill ? BOBOL_POLYGON_FILL_HATCH :
	   BOBOL_POLYGON_FILL_NONE;
}

static void
store_polygon_set_fill_flags(BObolPolygonVisual &visual, unsigned int flags)
{
    visual.fillFlags = store_polygon_valid_fill_flags(flags);
    visual.fill = (visual.fillFlags & BOBOL_POLYGON_FILL_HATCH) ?
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
store_polygon_zplane(plane_t dst, const BObolPolygonStoreRecord *rec)
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
store_polygon_project_to_zplane(const BObolPolygonStoreRecord *rec,
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
    vect2d_t line_slope;
    store_polygon_canonical_hatch_slope(line_slope, slope);
    struct bg_polygon *poly_fill = NULL;
    BU_GET(poly_fill, struct bg_polygon);
    bg_polygon_init(poly_fill);
    if (bg_polygon_hatch(poly_fill, poly, *vp, line_slope,
	    static_cast<fastf_t>(spacing)) || !poly_fill->num_contours) {
	bg_polygon_clear(poly_fill);
	BU_PUT(poly_fill, struct bg_polygon);
	return NULL;
    }
    return poly_fill;
}

static SoNode *
store_polygon_mesh_fill_node(const BObolPolygonStoreRecord &rec)
{
    if (!(store_polygon_fill_flags(rec.visual) & BOBOL_POLYGON_FILL_MESH) ||
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
			       sizeof(point2d_t), "BObol projected polygon fill points");
	for (size_t i = 0; i < contour.num_points; i++)
	    bg_plane_closest_pt(&projected[i][0], &projected[i][1],
				&plane, &contour.point[i]);

	int *faces = NULL;
	int numFaces = 0;
	int ret = bg_poly_triangulate(&faces, &numFaces, NULL, NULL, NULL, 0,
				      projected, contour.num_points, TRI_EAR_CLIPPING);
	bu_free(projected, "BObol projected polygon fill points");
	if (ret || numFaces <= 0 || !faces) {
	    if (faces)
		bu_free(faces, "BObol polygon fill faces");
	    return NULL;
	}

	indices.reserve(static_cast<size_t>(numFaces) * 3);
	for (int i = 0; i < numFaces * 3; i++)
	    indices.push_back(static_cast<int32_t>(faces[i]));
	bu_free(faces, "BObol polygon fill faces");
    } else {
	struct bg_polygon poly = BG_POLYGON_INIT_ZERO;
	(void)bg_polygon_copy(&poly, &rec.polygon);

	int *faces = NULL;
	int numFaces = 0;
	point_t *outPts = NULL;
	int numOutPts = 0;
	int ret = bg_polygon_triangulate(&faces, &numFaces, &outPts, &numOutPts,
					 &poly, TRI_EAR_CLIPPING);
	bg_polygon_clear(&poly);

	if (ret || numFaces <= 0 || numOutPts <= 0 || !faces || !outPts) {
	    if (faces)
		bu_free(faces, "BObol polygon fill faces");
	    if (outPts)
		bu_free(outPts, "BObol polygon fill points");
	    return NULL;
	}

	points.reserve(static_cast<size_t>(numOutPts));
	for (int i = 0; i < numOutPts; i++)
	    points.push_back(store_vec3(outPts[i]));

	indices.reserve(static_cast<size_t>(numFaces) * 3);
	for (int i = 0; i < numFaces * 3; i++)
	    indices.push_back(static_cast<int32_t>(faces[i]));

	bu_free(faces, "BObol polygon fill faces");
	bu_free(outPts, "BObol polygon fill points");
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
    shape->localSource = rec.scope == BObolFeatureScope::Local ? TRUE : FALSE;
    shape->sharedSource = rec.scope == BObolFeatureScope::Shared ? TRUE : FALSE;
    shape->nonDatabaseSource = TRUE;
    shape->drawMode = BOBOL_LOD_DRAW_SHADED;
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
store_polygon_hatch_fill_node(const BObolPolygonStoreRecord &rec)
{
    if (!(store_polygon_fill_flags(rec.visual) & BOBOL_POLYGON_FILL_HATCH) ||
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
			       static_cast<int32_t>(BObolLineCommand::Move) :
			       static_cast<int32_t>(BObolLineCommand::Draw));
	}
    }

    bg_polygon_clear(hatch);
    BU_PUT(hatch, struct bg_polygon);

    if (points.empty())
	return NULL;

    BObolFeatureStoreRecord feature;
    const std::string hatchName = store_string(rec.name) + ":hatch";
    feature.name = hatchName.c_str();
    feature.scope = rec.scope;
    feature.kind = BObolFeatureKind::PolygonOverlay;
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
store_polygon_node(const BObolPolygonStoreRecord &rec)
{
    if (!rec.visible)
	return new SoSeparator;

    std::vector<SbVec3f> points;
    std::vector<double> precisePoints;
    std::vector<int32_t> commands;
    for (size_t i = 0; i < rec.polygon.num_contours; i++) {
	const struct bg_poly_contour &contour = rec.polygon.contour[i];
	if (!contour.num_points || !contour.point)
	    continue;
	if (contour.num_points == 1) {
	    points.push_back(store_vec3(contour.point[0]));
	    precisePoints.push_back(contour.point[0][X]);
	    precisePoints.push_back(contour.point[0][Y]);
	    precisePoints.push_back(contour.point[0][Z]);
	    commands.push_back(static_cast<int32_t>(BObolLineCommand::Move));
	    continue;
	}

	/* Realize polygon outlines as independent edge pairs.  A repeated first
	 * vertex at the end of a chained line strip is not reliably preserved by
	 * all hosted rendering paths, and losing it drops the closure edge after
	 * interactive updates.  Explicit pairs are also a better match for
	 * per-edge picking and styling. */
	const bool closed = !contour.open ||
	    rec.type != BObolPolygonType::General;
	const size_t edgeCount = closed ? contour.num_points :
	    contour.num_points - 1;
	for (size_t j = 0; j < edgeCount; j++) {
	    const size_t next = (j + 1) % contour.num_points;
	    for (int endpoint = 0; endpoint < 2; endpoint++) {
		const size_t pointIndex = endpoint ? next : j;
		points.push_back(store_vec3(contour.point[pointIndex]));
		precisePoints.push_back(contour.point[pointIndex][X]);
		precisePoints.push_back(contour.point[pointIndex][Y]);
		precisePoints.push_back(contour.point[pointIndex][Z]);
		commands.push_back(endpoint ?
		    static_cast<int32_t>(BObolLineCommand::Draw) :
		    static_cast<int32_t>(BObolLineCommand::Move));
	    }
	}
    }

    BObolFeatureStoreRecord feature;
    feature.name = rec.name;
    feature.scope = rec.scope;
    feature.kind = BObolFeatureKind::PolygonOverlay;
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
store_polygon_append_point(BObolPolygonStoreRecord *rec,
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

	point_t appended;
	store_point(appended, model_point);
	(void)bg_polygon_append_point(&rec->polygon, contourIndex, appended);
    }
}

static void
store_polygon_set_rectangle(BObolPolygonStoreRecord *rec,
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

    point_t first_corner, opposite_corner;
    store_point(first_corner, store_polygon_plane_point(zplane, pfx, pfy));
    store_point(opposite_corner, store_polygon_plane_point(zplane, fx, fy));
    (void)bg_polygon_make_rectangle(&rec->polygon, zplane, first_corner,
	opposite_corner, square ? 1 : 0);
}

static void
store_polygon_set_ellipse(BObolPolygonStoreRecord *rec,
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

    point_t center, radius_point;
    store_point(center, store_polygon_plane_point(zplane, pfx, pfy));
    store_point(radius_point, store_polygon_plane_point(zplane,
	fx, fy));
    (void)bg_polygon_make_ellipse(&rec->polygon, zplane, center,
	radius_point, circle ? 1 : 0, nsegs);
}

BObolPolygonStore::BObolPolygonStore(void) : impl(new Impl)
{
}

BObolPolygonStore::BObolPolygonStore(BObolViewController *controller) :
    impl(new Impl)
{
    this->impl->controller = controller;
}

BObolPolygonStore::~BObolPolygonStore(void)
{
    delete this->impl;
    this->impl = NULL;
}

void
BObolPolygonStore::setController(BObolViewController *controller)
{
    this->impl->controller = controller;
}

BObolViewController *
BObolPolygonStore::controller(void) const
{
    return this->impl->controller;
}

uint64_t
BObolPolygonStore::referenceGeneration(void) const
{
    return this->impl->referenceGeneration;
}

const char *
BObolPolygonStore::name(BObolPolygonHandle handle) const
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    return rec ? rec->name.getString() : NULL;
}

void
BObolPolygonStore::clear(void)
{
    this->impl->clear();
}

BObolPolygonHandle
BObolPolygonStore::create(const SbString &name,
			    BObolFeatureScope scope,
			    BObolPolygonType type,
			    const SbVec3f &originPoint,
			    const fastf_t *viewPlane,
			    float viewZ)
{
    if (store_string(name).empty())
	return BObolPolygonHandle();

    const std::string key = store_key(scope, name);
    if (this->impl->names.find(key) != this->impl->names.end())
	return BObolPolygonHandle();

    BObolPolygonStoreRecord *rec = new BObolPolygonStoreRecord;
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

BObolPolygonHandle
BObolPolygonStore::find(const SbString &name, unsigned int scopeMask) const
{
    return this->impl->handle(this->impl->recordByName(name, scopeMask));
}

BObolPolygonHandle
BObolPolygonStore::selectAtModelPoint(const SbVec3f &point) const
{
    double best = INFINITY;
    const BObolPolygonStoreRecord *bestRec = NULL;
    for (std::map<uint64_t, BObolPolygonStoreRecord *>::const_iterator it =
	     this->impl->records.begin(); it != this->impl->records.end(); ++it) {
	BObolPolygonStoreRecord *rec = it->second;
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

BObolPolygonHandle
BObolPolygonStore::duplicate(BObolPolygonHandle handle,
			       const SbString &newName)
{
    BObolPolygonStoreRecord *src = this->impl->record(handle);
    if (!src || store_string(newName).empty())
	return BObolPolygonHandle();

    BObolPolygonHandle dstHandle = this->create(newName, src->scope,
				     src->type, src->originPoint, src->viewPlane, src->visual.viewZ);
    BObolPolygonStoreRecord *dst = this->impl->record(dstHandle);
    if (!dst)
	return BObolPolygonHandle();

    (void)bg_polygon_copy(&dst->polygon, &src->polygon);
    dst->visual = src->visual;
    dst->currentContour = src->currentContour;
    dst->currentPoint = src->currentPoint;
    dst->revision++;
    this->impl->realize(dst);
    return this->impl->handle(dst);
}

SbBool
BObolPolygonStore::update(BObolPolygonHandle handle,
			    BObolPolygonUpdate update)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    if (update == BObolPolygonUpdate::PointSelectClear) {
	rec->currentContour = -1;
	rec->currentPoint = -1;
    }
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BObolPolygonStore::updateScreenPoint(BObolPolygonHandle handle,
				       int x,
				       int y,
				       BObolPolygonUpdate update)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    return this->updateModelPoint(handle,
				  SbVec3f(static_cast<float>(x), static_cast<float>(y),
					  rec->originPoint[2]), update);
}

SbBool
BObolPolygonStore::updateModelPoint(BObolPolygonHandle handle,
				      const SbVec3f &point,
				      BObolPolygonUpdate update)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    if (update == BObolPolygonUpdate::PointAppend) {
	store_polygon_append_point(rec, point);
	rec->currentPoint = -1;
    } else if (update == BObolPolygonUpdate::PointMove) {
	if (rec->currentContour < 0 || rec->currentPoint < 0)
	    return FALSE;
	struct bg_poly_contour &contour =
		rec->polygon.contour[rec->currentContour];
	if (static_cast<size_t>(rec->currentPoint) >= contour.num_points)
	    return FALSE;
	store_point(contour.point[rec->currentPoint],
		    store_polygon_project_to_zplane(rec, point));
    } else if (update == BObolPolygonUpdate::PointSelect) {
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
    } else if (rec->type == BObolPolygonType::Circle) {
	store_polygon_set_ellipse(rec, point, TRUE);
    } else if (rec->type == BObolPolygonType::Ellipse) {
	store_polygon_set_ellipse(rec, point, FALSE);
    } else if (rec->type == BObolPolygonType::Square) {
	store_polygon_set_rectangle(rec, point, TRUE);
    } else if (rec->type == BObolPolygonType::Rectangle) {
	store_polygon_set_rectangle(rec, point, FALSE);
    }

    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BObolPolygonStore::move(BObolPolygonHandle handle,
			  const SbVec3f &currentPoint,
			  const SbVec3f &previousPoint)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    SbVec3f current = store_polygon_project_to_zplane(rec, currentPoint);
    SbVec3f previous = store_polygon_project_to_zplane(rec, previousPoint);
    SbVec3f delta = current - previous;
    vect_t translation = {delta[0], delta[1], delta[2]};
    (void)bg_polygon_translate(&rec->polygon, translation);
    rec->originPoint += delta;
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BObolPolygonStore::rename(BObolPolygonHandle handle,
			    const SbString &newName)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
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
BObolPolygonStore::remove(BObolPolygonHandle handle)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    this->impl->names.erase(store_key(rec->scope, rec->name));
    store_release_node(this->impl->controller, rec->node);
    bg_polygon_clear(&rec->polygon);
    this->impl->records.erase(rec->id);
    delete rec;
    if (this->impl->controller)
	this->impl->controller->requestRender("view-polygon-remove");
    return TRUE;
}

size_t
BObolPolygonStore::removeScope(unsigned int scopeMask)
{
    std::vector<BObolPolygonHandle> handles;
    for (std::map<uint64_t, BObolPolygonStoreRecord *>::const_iterator it =
	     this->impl->records.begin(); it != this->impl->records.end(); ++it) {
	if (!it->second)
	    continue;
	if (!(store_scope_bit(it->second->scope) & scopeMask))
	    continue;
	handles.push_back(this->impl->handle(it->second));
    }

    size_t removed = 0;
    for (size_t i = 0; i < handles.size(); i++)
	if (this->remove(handles[i]))
	    removed++;
    return removed;
}

SbBool
BObolPolygonStore::record(BObolPolygonHandle handle,
			    BObolPolygonRecord &recordOut) const
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    recordOut = BObolPolygonRecord();
    recordOut.handle = this->impl->handle(rec);
    recordOut.name = rec->name;
    recordOut.scope = rec->scope;
    recordOut.type = rec->type;
    recordOut.fillFlags = store_polygon_fill_flags(rec->visual);
    recordOut.fill = (recordOut.fillFlags & BOBOL_POLYGON_FILL_HATCH) ?
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
BObolPolygonStore::visitRecords(BObolPolygonRecordCallback callback,
				  void *userData) const
{
    if (!callback)
	return;
    for (std::map<uint64_t, BObolPolygonStoreRecord *>::const_iterator it =
	     this->impl->records.begin(); it != this->impl->records.end(); ++it) {
	BObolPolygonRecord rec;
	if (this->record(this->impl->handle(it->second), rec) &&
	    !callback(rec, userData))
	    return;
    }
}

SbBool
BObolPolygonStore::setCurrent(BObolPolygonHandle handle,
				long contour,
				long point)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
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
BObolPolygonStore::setContourOpen(BObolPolygonHandle handle,
				    long contour,
				    SbBool open)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec || contour < 0 ||
	contour >= static_cast<long>(rec->polygon.num_contours))
	return FALSE;
    if (bg_polygon_contour_open_set(&rec->polygon,
	static_cast<size_t>(contour), open ? 1 : 0))
	return FALSE;
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BObolPolygonStore::setAllContoursOpen(BObolPolygonHandle handle,
					SbBool open)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    if (bg_polygon_contours_open_set(&rec->polygon, open ? 1 : 0))
	return FALSE;
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BObolPolygonStore::clearSelectedPoint(BObolPolygonHandle handle)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    rec->currentContour = -1;
    rec->currentPoint = -1;
    rec->revision++;
    return TRUE;
}

SbBool
BObolPolygonStore::clearAllPointSelections(void)
{
    for (std::map<uint64_t, BObolPolygonStoreRecord *>::iterator it =
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
BObolPolygonStore::setVisible(BObolPolygonHandle handle, SbBool visible)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    rec->visible = visible;
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BObolPolygonStore::isVisible(BObolPolygonHandle handle,
			       SbBool &visibleOut) const
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    visibleOut = rec->visible;
    return TRUE;
}

SbBool
BObolPolygonStore::setVisual(BObolPolygonHandle handle,
			       const BObolPolygonVisual &visual)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
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
BObolPolygonStore::visual(BObolPolygonHandle handle,
			    BObolPolygonVisual &visualOut) const
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    visualOut = rec->visual;
    return TRUE;
}

SbBool
BObolPolygonStore::setEdgeColor(BObolPolygonHandle handle,
				  const SbColor &edgeColor)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    rec->visual.edgeColor = edgeColor;
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BObolPolygonStore::edgeColor(BObolPolygonHandle handle,
			       SbColor &edgeColorOut) const
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    edgeColorOut = rec->visual.edgeColor;
    return TRUE;
}

SbBool
BObolPolygonStore::setFill(BObolPolygonHandle handle,
			     SbBool fill,
			     const SbVec2f &slope,
			     float spacing)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    unsigned int flags = store_polygon_fill_flags(rec->visual);
    if (fill)
	flags |= BOBOL_POLYGON_FILL_HATCH;
    else
	flags &= ~BOBOL_POLYGON_FILL_HATCH;
    store_polygon_set_fill_flags(rec->visual, flags);
    rec->visual.fillSlope = slope;
    rec->visual.fillSpacing = spacing;
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BObolPolygonStore::setFillFlags(BObolPolygonHandle handle,
				  unsigned int fillFlags)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    store_polygon_set_fill_flags(rec->visual, fillFlags);
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BObolPolygonStore::setFillColor(BObolPolygonHandle handle,
				  const SbColor &fillColor)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    rec->visual.fillColor = fillColor;
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BObolPolygonStore::fillColor(BObolPolygonHandle handle,
			       SbColor &fillColorOut) const
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    fillColorOut = rec->visual.fillColor;
    return TRUE;
}

SbBool
BObolPolygonStore::setGeometry(BObolPolygonHandle handle,
				 const struct bg_polygon *polygon)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec || !polygon)
	return FALSE;
    (void)bg_polygon_copy(&rec->polygon, polygon);
    rec->revision++;
    this->impl->realize(rec);
    return TRUE;
}

const struct bg_polygon *
BObolPolygonStore::geometry(BObolPolygonHandle handle) const {
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    return rec ? &rec->polygon : NULL;
}

SbBool
BObolPolygonStore::copyGeometry(BObolPolygonHandle handle,
				  struct bg_polygon *polygonOut) const
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec || !polygonOut)
	return FALSE;
    return bg_polygon_copy(polygonOut, &rec->polygon) == 0 ? TRUE : FALSE;
}

SbBool
BObolPolygonStore::setUserData(BObolPolygonHandle handle, void *userData)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    rec->userData = userData;
    rec->revision++;
    return TRUE;
}

void *
BObolPolygonStore::userData(BObolPolygonHandle handle) const
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    return rec ? rec->userData : NULL;
}

double
BObolPolygonStore::area(BObolPolygonHandle handle, double viewScale) const
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return 0.0;
    return static_cast<double>(bg_polygon_area(&rec->polygon,
			       CLIPPER_MAX, rec->viewPlane, viewScale));
}

SbBool
BObolPolygonStore::overlaps(BObolPolygonHandle a,
			      BObolPolygonHandle b,
			      const struct bn_tol &tol,
			      double viewScale) const
{
    BObolPolygonStoreRecord *ra = this->impl->record(a);
    BObolPolygonStoreRecord *rb = this->impl->record(b);
    if (!ra || !rb)
	return FALSE;
    return bg_polygon_overlaps(&ra->polygon, &rb->polygon, ra->viewPlane,
			       &tol, viewScale) ? TRUE : FALSE;
}

SbBool
BObolPolygonStore::csg(BObolPolygonHandle target,
			 BObolPolygonHandle stencil,
			 enum bg_polygon_boolean_op op)
{
    BObolPolygonStoreRecord *rt = this->impl->record(target);
    BObolPolygonStoreRecord *rs = this->impl->record(stencil);
    if (!rt || !rs)
	return FALSE;
    if (op == BG_POLYGON_BOOLEAN_NONE || rs->polygon.num_contours == 0)
	return FALSE;

    if (rt->polygon.num_contours == 0) {
	if (op != BG_POLYGON_BOOLEAN_UNION)
	    return FALSE;
	if (bg_polygon_copy(&rt->polygon, &rs->polygon))
	    return FALSE;
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
    if (!bg_polygon_overlaps(&rt->polygon, &rs->polygon, rt->viewPlane,
			     &tol, CLIPPER_MAX))
	return FALSE;

    struct bg_polygon result = BG_POLYGON_INIT_ZERO;
    if (bg_polygon_boolean(&result, op, &rt->polygon,
	    &rs->polygon, CLIPPER_MAX, rt->viewPlane))
	return FALSE;

    (void)bg_polygon_move(&rt->polygon, &result);
    rt->type = BObolPolygonType::General;
    rt->revision++;
    this->impl->realize(rt);
    return TRUE;
}

BObolPolygonHandle
BObolPolygonStore::importSketch(const SbString &name,
				  BObolFeatureScope scope,
				  struct db_i *dbip,
				  struct directory *dp)
{
    if (store_string(name).empty() || !dbip || !dp)
	return BObolPolygonHandle();

    const std::string key = store_key(scope, name);
    if (this->impl->names.find(key) != this->impl->names.end())
	return BObolPolygonHandle();

    struct rt_sketch_polygon_data data;
    rt_sketch_polygon_data_init(&data);
    if (db_sketch_to_polygon_data(&data, store_string(name).c_str(),
				  dbip, dp) != 0) {
	rt_sketch_polygon_data_free(&data);
	return BObolPolygonHandle();
    }

    if (data.polygon.num_contours == 0 ||
	store_polygon_point_count(data.polygon) == 0) {
	rt_sketch_polygon_data_free(&data);
	return BObolPolygonHandle();
    }

    BObolPolygonStoreRecord *rec = new BObolPolygonStoreRecord;
    rec->id = this->impl->nextId++;
    rec->revision = 1;
    rec->name = name;
    rec->scope = scope;
    rec->type = store_polygon_type_from_rt(data.type);
    rec->originPoint = store_polygon_origin(data.polygon, data.origin_point);
    HMOVE(rec->viewPlane, data.vp);
    store_polygon_set_fill_flags(rec->visual,
				 data.fill_flag ? BOBOL_POLYGON_FILL_HATCH :
				 BOBOL_POLYGON_FILL_NONE);
    rec->visual.fillSlope = SbVec2f(static_cast<float>(data.fill_dir[0]),
				    static_cast<float>(data.fill_dir[1]));
    rec->visual.fillSpacing = static_cast<float>(data.fill_delta);
    rec->visual.fillColor = store_bu_to_sbcolor(data.fill_color);
    if (data.have_edge_color)
	rec->visual.edgeColor = store_bu_to_sbcolor(data.edge_color);
    rec->visual.viewZ = static_cast<float>(data.vZ);
    (void)bg_polygon_copy(&rec->polygon, &data.polygon);

    rt_sketch_polygon_data_free(&data);

    this->impl->records[rec->id] = rec;
    this->impl->names[key] = rec->id;
    this->impl->realize(rec);
    return this->impl->handle(rec);
}

SbBool
BObolPolygonStore::exportSketch(BObolPolygonHandle handle,
				  struct db_i *dbip,
				  const SbString &name) const
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec || !dbip || store_string(name).empty())
	return FALSE;

    struct rt_sketch_polygon_data data;
    rt_sketch_polygon_data_init(&data);
    data.type = store_polygon_type_to_rt(rec->type);
    data.fill_flag =
	(store_polygon_fill_flags(rec->visual) & BOBOL_POLYGON_FILL_HATCH) ?
	1 : 0;
    V2SET(data.fill_dir, rec->visual.fillSlope[0], rec->visual.fillSlope[1]);
    data.fill_delta = rec->visual.fillSpacing;
    store_sbcolor_to_bu(rec->visual.fillColor, &data.fill_color);
    store_point(data.origin_point, rec->originPoint);
    HMOVE(data.vp, rec->viewPlane);
    data.vZ = rec->visual.viewZ;
    data.have_edge_color = 1;
    store_sbcolor_to_bu(rec->visual.edgeColor, &data.edge_color);
    (void)bg_polygon_copy(&data.polygon, &rec->polygon);

    struct directory *dp = db_sketch_polygon_data_to_sketch(dbip,
			   store_string(name).c_str(), &data);
    rt_sketch_polygon_data_free(&data);
    return dp ? TRUE : FALSE;
}

size_t
BObolPolygonStore::snapCount(BObolPolygonHandle exclude) const
{
    size_t count = 0;
    for (std::map<uint64_t, BObolPolygonStoreRecord *>::const_iterator it =
	     this->impl->records.begin(); it != this->impl->records.end(); ++it) {
	if (!it->second)
	    continue;
	BObolPolygonHandle handle = this->impl->handle(it->second);
	if (exclude.isValid() && handle.id == exclude.id)
	    continue;
	count++;
    }
    return count;
}

SbBool
BObolPolygonStore::setSnapExclude(BObolPolygonHandle handle)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    this->impl->snapExclude = this->impl->handle(rec);
    return TRUE;
}

BObolPolygonHandle
BObolPolygonStore::snapExclude(void) const
{
    return this->impl->snapExclude;
}

SoNode *
BObolPolygonStore::node(BObolPolygonHandle handle) const
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    return rec ? rec->node : NULL;
}

struct BObolSelectionStore::Impl {
    std::vector<BObolSelectionRecord> records;
};

BObolSelectionStore::BObolSelectionStore(void) : impl(new Impl)
{
}

BObolSelectionStore::~BObolSelectionStore(void)
{
    delete this->impl;
    this->impl = NULL;
}

void
BObolSelectionStore::clear(const BObolFeatureOwner *owner, int kind)
{
    if (!owner && kind == BOBOL_SELECTION_ALL) {
	this->impl->records.clear();
	return;
    }

    this->impl->records.erase(std::remove_if(
				  this->impl->records.begin(),
				  this->impl->records.end(),
    [owner, kind](const BObolSelectionRecord &rec) {
	return store_owner_matches(rec.owner, owner) &&
	       store_selection_kind_matches(rec.kind, kind);
    }),
    this->impl->records.end());
}

size_t
BObolSelectionStore::count(const BObolFeatureOwner *owner, int kind) const
{
    if (!owner && kind == BOBOL_SELECTION_ALL)
	return this->impl->records.size();

    size_t cnt = 0;
    for (size_t i = 0; i < this->impl->records.size(); i++) {
	const BObolSelectionRecord &rec = this->impl->records[i];
	if (store_owner_matches(rec.owner, owner) &&
	    store_selection_kind_matches(rec.kind, kind))
	    cnt++;
    }
    return cnt;
}

SbBool
BObolSelectionStore::containsPath(const SbString &path, int kind,
				    const BObolFeatureOwner *owner) const
{
    const std::string target = store_string(path);
    if (target.empty())
	return FALSE;
    for (size_t i = 0; i < this->impl->records.size(); i++) {
	const BObolSelectionRecord &rec = this->impl->records[i];
	if (store_string(rec.path) == target &&
	    store_owner_matches(rec.owner, owner) &&
	    store_selection_kind_matches(rec.kind, kind))
	    return TRUE;
    }
    return FALSE;
}

SbBool
BObolSelectionStore::addPath(const SbString &path, int kind,
			       const BObolFeatureOwner *owner)
{
    if (store_string(path).empty())
	return FALSE;
    int recordKind = store_selection_record_kind(kind);
    if (this->containsPath(path, recordKind, owner))
	return TRUE;
    BObolSelectionRecord rec;
    rec.path = path;
    rec.kind = recordKind;
    if (owner)
	rec.owner = *owner;
    this->impl->records.push_back(rec);
    return TRUE;
}

SbBool
BObolSelectionStore::setPath(const SbString &path, int kind,
			       const BObolFeatureOwner *owner)
{
    this->clear(owner, kind);
    if (store_string(path).empty())
	return TRUE;
    return this->addPath(path, kind, owner);
}

SbBool
BObolSelectionStore::removePath(const SbString &path, int kind,
				  const BObolFeatureOwner *owner)
{
    const std::string target = store_string(path);
    for (std::vector<BObolSelectionRecord>::iterator it =
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

const BObolSelectionRecord *
BObolSelectionStore::record(size_t index) const
{
    return index < this->impl->records.size() ?
	   &this->impl->records[index] : NULL;
}

SbBool
BObolSelectionStore::addRecord(const BObolSelectionRecord &record)
{
    if (record.path.getLength() == 0)
	return FALSE;
    BObolSelectionRecord rec = record;
    if (rec.kind == BOBOL_SELECTION_ALL)
	rec.kind = BOBOL_SELECTION_SELECTED_PATH;
    if (this->containsPath(rec.path, rec.kind, &rec.owner))
	return TRUE;
    this->impl->records.push_back(rec);
    return TRUE;
}

SbBool
BObolSelectionStore::setRecords(
    const std::vector<BObolSelectionRecord> &records)
{
    this->impl->records.clear();
    for (size_t i = 0; i < records.size(); i++)
	this->addRecord(records[i]);
    return TRUE;
}

void
BObolSelectionStore::visitPaths(
    int (*callback)(const SbString &path, void *userData),
    void *userData,
    const BObolFeatureOwner *owner,
    int kind) const
{
    if (!callback)
	return;
    for (size_t i = 0; i < this->impl->records.size(); i++) {
	const BObolSelectionRecord &rec = this->impl->records[i];
	if (!store_owner_matches(rec.owner, owner) ||
	    !store_selection_kind_matches(rec.kind, kind))
	    continue;
	if (!callback(rec.path, userData))
	    return;
    }
}

SbBool
BObolSelectionStore::applyPickResults(
    const std::vector<BObolSelectionRecord> &records,
    void (*selectedPathCallback)(const SbString &path, void *userData),
    void *userData,
    const BObolFeatureOwner *owner)
{
    this->clear(owner, BOBOL_SELECTION_ALL);
    for (size_t i = 0; i < records.size(); i++) {
	BObolSelectionRecord rec = records[i];
	if (owner)
	    rec.owner = *owner;
	this->addRecord(rec);
	if (selectedPathCallback &&
	    store_selection_record_kind(rec.kind) ==
	    BOBOL_SELECTION_SELECTED_PATH)
	    selectedPathCallback(rec.path, userData);
    }
    return TRUE;
}
