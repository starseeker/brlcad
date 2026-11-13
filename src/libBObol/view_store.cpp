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
#include "identity_counter_private.h"
#include "view_controller_private.h"

#include "bg/line_layer.h"
#include "bg/plane.h"
#include "bg/polygon.h"
#include "bn/tol.h"
#include "rt/primitives/sketch.h"
#include "bu/malloc.h"

#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoDepthBuffer.h>
#include <Inventor/nodes/SoFont.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoText2.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/annex/HUD/nodekits/SoHUDKit.h>
#include <Inventor/fields/SoMField.h>
#include <Inventor/misc/SoChildList.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

static std::atomic<uint64_t> store_reference_generation_counter(1);
static constexpr const char *store_feature_controller_detach_reason =
    "view-feature-controller-detach";
static constexpr const char *store_feature_controller_attach_reason =
    "view-feature-controller-attach";
static constexpr const char *store_feature_publication_reason =
    "view-feature-store";
static constexpr const char *store_feature_removal_reason =
    "view-feature-remove";
static constexpr const char *store_feature_clear_reason =
    "view-feature-clear";
static constexpr const char *store_feature_publication_command =
    "applyPublication";
static constexpr const char *store_feature_custom_node_command =
    "publishCustomNode";
static constexpr const char *store_feature_indexed_points_reason =
    "view-feature-indexed-face-points";
static constexpr const char *store_feature_selected_primitives_reason =
    "view-feature-selected-primitives";
static constexpr const char *store_feature_highlighted_primitives_reason =
    "view-feature-highlighted-primitives";
static constexpr const char *store_feature_remove_command = "remove";
static constexpr const char *store_feature_clear_command = "clear";

static uint64_t
store_reference_generation_next(void)
{
    return bobol_atomic_nonzero_identity_take(
	store_reference_generation_counter);
}

static void
store_revision_advance(uint64_t &revision)
{
    bobol_identity_advance(revision);
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
store_overlay_equal(const BObolOverlayInfo &a, const BObolOverlayInfo &b)
{
    return a.isOverlay == b.isOverlay &&
	a.ownerToken == b.ownerToken &&
	a.role == b.role &&
	a.overlayClass == b.overlayClass &&
	a.lifecycle == b.lifecycle &&
	a.order == b.order &&
	a.sortOrder == b.sortOrder &&
	a.sourcePath == b.sourcePath ? TRUE : FALSE;
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
store_detach_node(SoGroup *group, SoNode *node)
{
    if (!group || !node)
	return;

    const int index = group->findChild(node);
    if (index >= 0)
	group->removeChild(index);
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

BObolFeaturePublication::BObolFeaturePublication(void) :
    action(BObolFeaturePublicationAction::Replace),
    name(""),
    kind(BObolFeatureKind::Unknown),
    scope(BObolFeatureScope::Shared),
    style(),
    owner(),
    overlay(),
    points(),
    commands(),
    labels(),
    selectedPrimitives(),
    highlightedPrimitives(),
    identity(""),
    editIntentId(""),
    editIntentRole(""),
    sourceRevision(0),
    inputsRevision(0),
    customNode(NULL)
{
}

BObolFeatureStorePublication::BObolFeatureStorePublication(void) :
    store(NULL),
    features()
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
    selected(FALSE),
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
    sketchName(""),
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
    SoGroup *attachmentRoot;

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
	node(NULL),
	attachmentRoot(NULL)
    {
    }
};

static SoNode *store_rebuild_node_for_feature(
    const BObolFeatureStoreRecord &rec);

static bool
store_feature_overlay_less(const BObolFeatureStoreRecord *a,
	const BObolFeatureStoreRecord *b)
{
    if (a->overlay.order != b->overlay.order)
	return static_cast<int>(a->overlay.order) <
	    static_cast<int>(b->overlay.order);
    if (a->overlay.sortOrder != b->overlay.sortOrder)
	return a->overlay.sortOrder < b->overlay.sortOrder;
    return a->id < b->id;
}

static void
store_append_unique_node(std::vector<SoNode *> &nodes, SoNode *node)
{
    if (node && std::find(nodes.begin(), nodes.end(), node) == nodes.end())
	nodes.push_back(node);
}

static SoGroup *
store_feature_attachment_root(BObolViewController *controller,
	const BObolFeatureStoreRecord *rec)
{
    if (!controller)
	return NULL;

    /* Screen overlays belong outside the CAD render batch.  Besides giving
     * them deterministic last-pass ordering, this keeps retained HUD nodes
     * visible when a compact CAD batch bypasses ordinary source traversal. */
    if (rec && rec->overlay.isOverlay &&
	(rec->overlay.role == BObolOverlayRole::Screen ||
	 rec->overlay.order == BObolOverlayOrder::Screen))
	return controller->getFramebufferOverlayRoot();

    return store_controller_root_group(controller);
}

static void
store_feature_release_node(BObolFeatureStoreRecord *rec)
{
    if (!rec || !rec->node)
	return;

    store_detach_node(rec->attachmentRoot, rec->node);
    rec->node->unref();
    rec->node = NULL;
    rec->attachmentRoot = NULL;
}

/* Rebuilt typed features are immutable presentation values.  Compare their
 * field trees before replacing the retained node so periodic HUD/faceplate
 * synchronization is idempotent.  Group children are not Coin fields and
 * therefore need an explicit recursive comparison.  Custom nodes are opaque
 * plugin-owned objects and are compared by identity only below. */
static const SoNode *
store_feature_first_node_of_type(const SoNode *node, SoType type)
{
    if (!node)
	return NULL;
    if (node->isOfType(type))
	return node;

    const SoChildList *children = node->getChildren();
    if (!children)
	return NULL;
    for (int i = 0; i < children->getLength(); i++) {
	const SoNode *found = store_feature_first_node_of_type(
	    (*children)[i], type);
	if (found)
	    return found;
    }
    return NULL;
}

static SbBool
store_feature_nodes_equal(const SoNode *a, const SoNode *b)
{
    if (a == b)
	return TRUE;
    if (!a || !b || a->getTypeId() != b->getTypeId())
	return FALSE;

    /* Feature-store HUD kits wrap exactly one SoBRLVListShape.  The kit owns
     * renderer-maintained viewport/camera state, so comparing its fields (or
     * its complete nodekit child graph) makes identical line geometry look
     * different after the first render or any widget resize.  Compare the
     * wrapped immutable feature instead and ignore runtime projection state. */
    if (a->isOfType(SoHUDKit::getClassTypeId())) {
	return store_feature_nodes_equal(
	    store_feature_first_node_of_type(a,
		SoBRLVListShape::getClassTypeId()),
	    store_feature_first_node_of_type(b,
		SoBRLVListShape::getClassTypeId()));
    }
    if (!a->fieldsAreEqual(b))
	return FALSE;

    const SbBool aGroup = a->isOfType(SoGroup::getClassTypeId());
    const SbBool bGroup = b->isOfType(SoGroup::getClassTypeId());
    if (aGroup != bGroup)
	return FALSE;
    if (!aGroup)
	return TRUE;

    const SoGroup *ga = static_cast<const SoGroup *>(a);
    const SoGroup *gb = static_cast<const SoGroup *>(b);

    /* These BRL-CAD adapters expose the complete semantic value in fields and
     * rebuild renderer/helper children from those fields.  Comparing those
     * implementation children would make an identical HUD/grid/ADC publish
     * look different merely because nodekit internals have fresh identities. */
    if (a->isOfType(SoBRLHUDLabelOverlay::getClassTypeId()))
	return TRUE;
    if (ga->getNumChildren() != gb->getNumChildren())
	return FALSE;
    for (int i = 0; i < ga->getNumChildren(); i++) {
	if (!store_feature_nodes_equal(ga->getChild(i), gb->getChild(i)))
	    return FALSE;
    }
    return TRUE;
}

static bool
store_feature_styles_equal(const BObolFeatureStyle &a,
	const BObolFeatureStyle &b)
{
    return a.hasVisible == b.hasVisible && a.visible == b.visible &&
	a.hasSelectable == b.hasSelectable &&
	a.selectable == b.selectable && a.hasColor == b.hasColor &&
	a.color == b.color && a.hasLineWidth == b.hasLineWidth &&
	a.lineWidth == b.lineWidth && a.hasLineStyle == b.hasLineStyle &&
	a.lineStyle == b.lineStyle &&
	a.hasTransparency == b.hasTransparency &&
	std::memcmp(&a.transparency, &b.transparency,
	    sizeof(a.transparency)) == 0 && a.hasArrow == b.hasArrow &&
	a.arrow == b.arrow && a.hasArrowTip == b.hasArrowTip &&
	std::memcmp(&a.arrowTipLength, &b.arrowTipLength,
	    sizeof(a.arrowTipLength)) == 0 &&
	std::memcmp(&a.arrowTipWidth, &b.arrowTipWidth,
	    sizeof(a.arrowTipWidth)) == 0 && a.hud == b.hud;
}

static bool
store_feature_owners_equal(const BObolFeatureOwner &a,
	const BObolFeatureOwner &b)
{
    return a.ownerToken == b.ownerToken && a.ownerId == b.ownerId &&
	a.ownerRole == b.ownerRole && a.generation == b.generation &&
	a.resultCallback == b.resultCallback &&
	a.callbackUserData == b.callbackUserData;
}

static bool
store_feature_overlays_equal(const BObolOverlayInfo &a,
	const BObolOverlayInfo &b)
{
    return a.isOverlay == b.isOverlay && a.ownerToken == b.ownerToken &&
	a.role == b.role && a.overlayClass == b.overlayClass &&
	a.lifecycle == b.lifecycle && a.order == b.order &&
	a.sortOrder == b.sortOrder && a.sourcePath == b.sourcePath;
}

static bool
store_feature_publications_equal(const BObolFeatureStoreRecord &current,
	const BObolFeatureStoreRecord &candidate, const SoNode *node)
{
    return current.kind == candidate.kind &&
	current.scope == candidate.scope &&
	store_feature_styles_equal(current.style, candidate.style) &&
	store_feature_owners_equal(current.owner, candidate.owner) &&
	store_feature_overlays_equal(current.overlay, candidate.overlay) &&
	current.metadata.empty() && current.primitiveMetadata.empty() &&
	current.selectedPrimitives == candidate.selectedPrimitives &&
	current.highlightedPrimitives == candidate.highlightedPrimitives &&
	!current.compactEdit &&
	current.identity == candidate.identity &&
	current.editIntentId == candidate.editIntentId &&
	current.editIntentRole == candidate.editIntentRole &&
	current.sourceRevision == candidate.sourceRevision &&
	current.inputsRevision == candidate.inputsRevision &&
	(candidate.kind == BObolFeatureKind::CustomNode ?
	 current.node == node : store_feature_nodes_equal(current.node, node));
}

static bool
store_custom_node_primitives_equal(const BObolFeaturePublication &entry)
{
    if (!entry.customNode || !entry.customNode->isOfType(
	    SoBRLMeshShape::getClassTypeId()))
	return entry.selectedPrimitives.empty() &&
	    entry.highlightedPrimitives.empty();

    const SoBRLMeshShape *mesh =
	static_cast<const SoBRLMeshShape *>(entry.customNode);
    if (mesh->selectedPrimitive.getNum() !=
	    static_cast<int>(entry.selectedPrimitives.size()) ||
	mesh->highlightedPrimitive.getNum() !=
	    static_cast<int>(entry.highlightedPrimitives.size()))
	return false;
    for (int i = 0; i < mesh->selectedPrimitive.getNum(); ++i)
	if (mesh->selectedPrimitive[i] !=
	    entry.selectedPrimitives[static_cast<size_t>(i)])
	    return false;
    for (int i = 0; i < mesh->highlightedPrimitive.getNum(); ++i)
	if (mesh->highlightedPrimitive[i] !=
	    entry.highlightedPrimitives[static_cast<size_t>(i)])
	    return false;
    return true;
}

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
    uint64_t presentationRevision;
    uint64_t nextId;
    std::map<uint64_t, BObolFeatureStoreRecord *> records;
    std::map<std::string, uint64_t> names;
    std::map<std::string, uint64_t> ownerGenerations;

    Impl(void) : controller(NULL),
	referenceGeneration(store_reference_generation_next()),
	presentationRevision(0), nextId(1),
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
		notify(it->second, BObolCommandResultStatus::Removed, "clear");
		store_feature_release_node(it->second);
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

    BObolFeatureStoreRecord *recordByKey(const SbString &name,
	BObolFeatureScope scope, const BObolFeatureOwner *owner) const
    {
	const std::string key = store_key(scope, name, owner);
	std::map<std::string, uint64_t>::const_iterator nit = names.find(key);
	if (nit == names.end())
	    return NULL;
	std::map<uint64_t, BObolFeatureStoreRecord *>::const_iterator rit =
	    records.find(nit->second);
	return rit == records.end() ? NULL : rit->second;
    }

    static void configureUpsert(BObolFeatureStoreRecord &rec,
	BObolFeatureKind kind, const BObolFeatureStyle *style,
	const BObolFeatureOwner *owner)
    {
	rec.kind = kind;
	rec.metadata.clear();
	rec.primitiveMetadata.clear();
	rec.selectedPrimitives.clear();
	rec.highlightedPrimitives.clear();
	rec.compactEdit = FALSE;
	rec.compactSummary = BObolCompactInstanceSummary();
	if (style)
	    rec.style = *style;
	if (owner)
	    rec.owner = *owner;
    }

    std::unique_ptr<BObolFeatureStoreRecord> prepareExistingUpsert(
	const SbString &name, BObolFeatureScope scope, BObolFeatureKind kind,
	const BObolFeatureStyle *style, const BObolFeatureOwner *owner) const
    {
	BObolFeatureStoreRecord *current = recordByKey(name, scope, owner);
	if (!current)
	    return std::unique_ptr<BObolFeatureStoreRecord>();
	auto candidate = std::make_unique<BObolFeatureStoreRecord>(*current);
	store_revision_advance(candidate->revision);
	configureUpsert(*candidate, kind, style, owner);
	return candidate;
    }

    BObolFeatureHandle handle(const BObolFeatureStoreRecord *rec) const
    {
	return rec ? BObolFeatureHandle(rec->id, rec->revision) :
	       BObolFeatureHandle();
    }

    struct FeaturePlacement {
	BObolFeatureStoreRecord *record;
	SoNode *node;
	SoGroup *root;
    };

    std::vector<FeaturePlacement> candidateFeaturePlacements(
	BObolFeatureStoreRecord &candidate)
    {
	std::vector<FeaturePlacement> placements;
	placements.reserve(records.size() + 1);
	bool replacedExisting = false;
	for (const auto &entry : records) {
	    const bool replace = entry.first == candidate.id;
	    BObolFeatureStoreRecord *rec = replace ? &candidate : entry.second;
	    replacedExisting = replacedExisting || replace;
	    if (rec && rec->node)
		placements.push_back({rec, rec->node, rec->attachmentRoot});
	}
	if (!replacedExisting && candidate.node)
	    placements.push_back({&candidate, candidate.node,
		candidate.attachmentRoot});
	return placements;
    }

    std::vector<FeaturePlacement> controllerFeaturePlacements(
	BObolViewController *nextController)
    {
	std::vector<FeaturePlacement> placements;
	placements.reserve(records.size());
	for (const auto &entry : records) {
	    BObolFeatureStoreRecord *rec = entry.second;
	    if (rec && rec->node) {
		placements.push_back({rec, rec->node,
		    store_feature_attachment_root(nextController, rec)});
	    }
	}
	return placements;
    }

    static bool rootOrderChanged(SoGroup *root,
	const std::vector<SoNode *> &next)
    {
	if (next.size() != static_cast<size_t>(root->getNumChildren()))
	    return true;
	for (size_t i = 0; i < next.size(); ++i)
	    if (root->getChild(static_cast<int>(i)) != next[i])
		return true;
	return false;
    }

    static void notifyRootReplacements(
	std::vector<std::unique_ptr<SoChildList::Replacement>> &replacements,
	std::exception_ptr &failure)
    {
	for (auto &replacement : replacements) {
	    try {
		replacement->notify();
	    } catch (...) {
		if (!failure)
		    failure = std::current_exception();
	    }
	}
    }

    struct PreparedFeaturePresentation {
	BObolViewController *controller = NULL;
	std::optional<BObolLodControlTransitionScope> transition;
	std::optional<BObolPreparedRenderRequest> request;
    };

    static void prepareControllerPresentation(BObolViewController *target,
	const char *reason, PreparedFeaturePresentation &publication)
    {
	if (!target)
	    return;
	publication.controller = target;
	publication.transition.emplace(target);
	publication.request.emplace(target->prepareRenderRequest(reason,
	    BObolViewController::RenderRequestIntent::PRESENTATION));
    }

    static void commitControllerPresentation(
	PreparedFeaturePresentation &publication)
    {
	if (publication.request)
	    publication.controller->commitRenderRequest(*publication.request);
    }

    static void notifyControllerPresentation(
	const PreparedFeaturePresentation &publication,
	std::exception_ptr &failure)
    {
	if (!publication.request)
	    return;
	try {
	    publication.controller->notifyRenderRequest(*publication.request);
	} catch (...) {
	    if (!failure)
		failure = std::current_exception();
	}
    }

    void prepareFeaturePresentation(SbBool changed,
	PreparedFeaturePresentation &publication,
	const char *reason = store_feature_publication_reason) const
    {
	if (!changed)
	    return;
	prepareControllerPresentation(controller,
	    reason ? reason : store_feature_publication_reason, publication);
    }

    void commitFeaturePresentation(SbBool changed,
	PreparedFeaturePresentation &publication)
    {
	if (!changed)
	    return;
	store_revision_advance(presentationRevision);
	commitControllerPresentation(publication);
    }

    void notifyFeaturePresentation(
	const PreparedFeaturePresentation &publication,
	std::exception_ptr &failure) const
    {
	notifyControllerPresentation(publication, failure);
    }

    std::vector<SoNode *> featureRootOrder(SoGroup *root,
	const std::vector<FeaturePlacement> &placements,
	SoNode *replacedNode = NULL, SoNode *replacementNode = NULL) const
    {
	std::unordered_set<SoNode *> desiredNodes;
	desiredNodes.reserve(placements.size());
	for (const FeaturePlacement &placement : placements)
	    if (placement.node && placement.root == root)
		desiredNodes.insert(placement.node);

	std::unordered_set<SoNode *> managedNodes;
	managedNodes.reserve(records.size());
	for (const auto &entry : records) {
	    const BObolFeatureStoreRecord *rec = entry.second;
	    if (rec && rec->node && rec->attachmentRoot == root)
		managedNodes.insert(rec->node);
	}

	/* Preserve every unrelated child in place. A replacement keeps its old
	 * position when the record held the node's last attachment. */
	std::vector<SoNode *> next;
	next.reserve(static_cast<size_t>(root->getNumChildren()) +
	    placements.size());
	int replacedIndex = -1;
	for (int i = 0; i < root->getNumChildren(); ++i) {
	    SoNode *child = root->getChild(i);
	    if (managedNodes.count(child) && !desiredNodes.count(child)) {
		if (child == replacedNode)
		    replacedIndex = static_cast<int>(next.size());
		continue;
	    }
	    next.push_back(child);
	}

	if (replacementNode && desiredNodes.count(replacementNode) &&
	    std::find(next.begin(), next.end(), replacementNode) == next.end() &&
	    replacedIndex >= 0)
	    next.insert(next.begin() + replacedIndex, replacementNode);
	for (const FeaturePlacement &placement : placements)
	    if (placement.root == root)
		store_append_unique_node(next, placement.node);

	/* Overlay order is part of presentation. Shared custom-node records retain
	 * one occurrence, and unrelated root children keep their relative order. */
	std::vector<BObolFeatureStoreRecord *> overlays;
	overlays.reserve(placements.size());
	for (const FeaturePlacement &placement : placements)
	    if (placement.node && placement.root == root &&
		placement.record->overlay.isOverlay)
		overlays.push_back(placement.record);
	std::sort(overlays.begin(), overlays.end(), store_feature_overlay_less);
	std::vector<SoNode *> overlayNodes;
	overlayNodes.reserve(overlays.size());
	for (const BObolFeatureStoreRecord *rec : overlays)
	    store_append_unique_node(overlayNodes, rec->node);
	if (!overlayNodes.empty()) {
	    next.erase(std::remove_if(next.begin(), next.end(),
		[&overlayNodes](SoNode *node) {
		    return std::find(overlayNodes.begin(), overlayNodes.end(), node) !=
			overlayNodes.end();
		}), next.end());
	    next.insert(next.end(), overlayNodes.begin(), overlayNodes.end());
	}
	return next;
    }

    void reorderOverlayNodes(SoGroup *root)
    {
	if (!root)
	    return;

	const std::vector<FeaturePlacement> placements =
	    controllerFeaturePlacements(controller);
	const std::vector<SoNode *> next = featureRootOrder(root, placements);
	if (!rootOrderChanged(root, next))
	    return;

	auto replacement = root->getChildren()->prepareReplacement(next);
	replacement->commit();
	replacement->notify();
    }

    void migrateController(BObolViewController *nextController)
    {
	BObolViewController *previousController = controller;
	if (previousController == nextController)
	    return;

	std::vector<FeaturePlacement> placements =
	    controllerFeaturePlacements(nextController);
	std::vector<SoGroup *> roots;
	roots.reserve(placements.size() + placements.size());
	auto appendRoot = [&roots](SoGroup *root) {
	    if (root && std::find(roots.begin(), roots.end(), root) == roots.end())
		roots.push_back(root);
	};
	for (const auto &entry : records) {
	    const BObolFeatureStoreRecord *rec = entry.second;
	    if (rec && rec->node)
		appendRoot(rec->attachmentRoot);
	}
	for (const FeaturePlacement &placement : placements)
	    appendRoot(placement.root);

	std::vector<std::unique_ptr<SoChildList::Replacement>> replacements;
	replacements.reserve(roots.size());
	for (SoGroup *root : roots) {
	    std::vector<SoNode *> next = featureRootOrder(root, placements);
	    if (rootOrderChanged(root, next))
		replacements.push_back(root->getChildren()->prepareReplacement(next));
	}
	PreparedFeaturePresentation detachPresentation;
	PreparedFeaturePresentation attachPresentation;
	if (!roots.empty()) {
	    prepareControllerPresentation(previousController,
		store_feature_controller_detach_reason, detachPresentation);
	    prepareControllerPresentation(nextController,
		store_feature_controller_attach_reason, attachPresentation);
	}

	for (auto &replacement : replacements)
	    replacement->commit();
	controller = nextController;
	for (FeaturePlacement &placement : placements)
	    placement.record->attachmentRoot = placement.root;
	if (roots.empty())
	    return;
	store_revision_advance(presentationRevision);
	commitControllerPresentation(detachPresentation);
	commitControllerPresentation(attachPresentation);

	std::exception_ptr failure;
	notifyRootReplacements(replacements, failure);
	notifyControllerPresentation(detachPresentation, failure);
	notifyControllerPresentation(attachPresentation, failure);
	if (failure)
	    std::rethrow_exception(failure);
    }

    SbBool publishExistingRecord(
	std::unique_ptr<BObolFeatureStoreRecord> candidate, SoNode *node,
	SbBool rollbackUnchangedPublish = FALSE)
    {
	if (!candidate)
	    return FALSE;
	std::map<uint64_t, BObolFeatureStoreRecord *>::iterator currentEntry =
	    records.find(candidate->id);
	if (currentEntry == records.end() || !currentEntry->second)
	    return FALSE;
	BObolFeatureStoreRecord *current = currentEntry->second;
	SoNode *currentNode = current->node;
	SoGroup *currentRoot = current->attachmentRoot;
	SoGroup *desiredRoot = store_feature_attachment_root(controller,
	    candidate.get());
	bool retainedCandidateNode = node && node != currentNode;
	if (retainedCandidateNode)
	    node->ref();
	try {
	    if (candidate->kind != BObolFeatureKind::CustomNode && currentNode &&
		node && currentRoot == desiredRoot &&
		store_feature_nodes_equal(currentNode, node)) {
		node->unref();
		retainedCandidateNode = false;
		node = currentNode;
		if (rollbackUnchangedPublish && candidate->revision > 1)
		    candidate->revision--;
	    }
	    candidate->node = node;
	    candidate->attachmentRoot = node ? desiredRoot : NULL;
	    const std::vector<FeaturePlacement> placements =
		candidateFeaturePlacements(*candidate);

	    std::vector<SoGroup *> roots;
	    if (currentRoot)
		roots.push_back(currentRoot);
	    if (desiredRoot && std::find(roots.begin(), roots.end(), desiredRoot) ==
		    roots.end())
		roots.push_back(desiredRoot);
	    std::vector<std::unique_ptr<SoChildList::Replacement>> replacements;
	    replacements.reserve(roots.size());
	    for (SoGroup *root : roots) {
		std::vector<SoNode *> next = featureRootOrder(root, placements,
		    currentNode, node);
		if (rootOrderChanged(root, next))
		    replacements.push_back(
			root->getChildren()->prepareReplacement(next));
	    }
	    const SbBool presentationChanged = node != currentNode ||
		desiredRoot != currentRoot || !replacements.empty();
	    PreparedFeaturePresentation presentation;
	    prepareFeaturePresentation(presentationChanged, presentation);

	    for (auto &replacement : replacements)
		replacement->commit();
	    BObolFeatureStoreRecord *installed = candidate.release();
	    currentEntry->second = installed;
	    if (currentNode && currentNode != node)
		currentNode->unref();
	    delete current;
	    commitFeaturePresentation(presentationChanged, presentation);

	    std::exception_ptr failure;
	    notifyRootReplacements(replacements, failure);
	    notifyFeaturePresentation(presentation, failure);
	    if (failure)
		std::rethrow_exception(failure);
	    return presentationChanged;
	} catch (...) {
	    if (candidate && retainedCandidateNode)
		node->unref();
	    throw;
	}
    }

    struct RecordPublication {
	BObolFeatureStoreRecord *record;
	SbBool nodeChanged;

	RecordPublication(BObolFeatureStoreRecord *published = NULL,
	    SbBool changed = FALSE) : record(published), nodeChanged(changed)
	{
	}
    };

    RecordPublication publishNewRecord(
	std::unique_ptr<BObolFeatureStoreRecord> candidate, SoNode *node,
	const std::string &key, uint64_t followingId)
    {
	if (!candidate)
	    return RecordPublication();
	const uint64_t id = candidate->id;
	const bool retainedNode = node != NULL;
	if (retainedNode)
	    node->ref();
	try {
	    SoGroup *desiredRoot = store_feature_attachment_root(controller,
		candidate.get());
	    candidate->node = node;
	    candidate->attachmentRoot = node ? desiredRoot : NULL;
	    const std::vector<FeaturePlacement> placements =
		candidateFeaturePlacements(*candidate);

	    std::unique_ptr<SoChildList::Replacement> rootReplacement;
	    if (desiredRoot && node) {
		const std::vector<SoNode *> next = featureRootOrder(desiredRoot,
		    placements);
		if (rootOrderChanged(desiredRoot, next))
		    rootReplacement =
			desiredRoot->getChildren()->prepareReplacement(next);
	    }

	    std::map<uint64_t, BObolFeatureStoreRecord *> preparedRecords;
	    const auto preparedRecord = preparedRecords.emplace(id,
		candidate.get());
	    auto recordEntry = preparedRecords.extract(preparedRecord.first);
	    std::map<std::string, uint64_t> preparedNames;
	    const auto preparedName = preparedNames.emplace(key, id);
	    auto nameEntry = preparedNames.extract(preparedName.first);

	    const SbBool presentationChanged = node ? TRUE : FALSE;
	    PreparedFeaturePresentation presentation;
	    prepareFeaturePresentation(presentationChanged, presentation);

	    records.insert(std::move(recordEntry));
	    names.insert(std::move(nameEntry));
	    if (rootReplacement)
		rootReplacement->commit();
	    nextId = followingId;
	    candidate.release();
	    commitFeaturePresentation(presentationChanged, presentation);

	    std::exception_ptr failure;
	    if (rootReplacement) {
		try {
		    rootReplacement->notify();
		} catch (...) {
		    failure = std::current_exception();
		}
	    }
	    notifyFeaturePresentation(presentation, failure);
	    if (failure)
		std::rethrow_exception(failure);
	    return RecordPublication(record(BObolFeatureHandle(id, 0)),
		presentationChanged);
	} catch (...) {
	    if (candidate && retainedNode)
		node->unref();
	    throw;
	}
    }

    template <typename Configure, typename BuildNode>
    RecordPublication publishExistingEdit(BObolFeatureStoreRecord *current,
	Configure configure, BuildNode buildNode,
	SbBool rollbackUnchangedPublish = FALSE)
    {
	if (!current)
	    return RecordPublication();
	auto candidate = std::make_unique<BObolFeatureStoreRecord>(*current);
	configure(*candidate);
	SoNode *node = buildNode(*candidate);
	const uint64_t id = candidate->id;
	const SbBool changed = publishExistingRecord(std::move(candidate), node,
	    rollbackUnchangedPublish);
	return RecordPublication(record(BObolFeatureHandle(id, 0)), changed);
    }

    struct PreparedFeatureResult {
	BObolFeatureOwner owner;
	BObolCommandResult result;
	bool active = false;

	void notify(std::exception_ptr &failure) const
	{
	    if (!active)
		return;
	    try {
		owner.resultCallback(result, owner.callbackUserData);
	    } catch (...) {
		if (!failure)
		    failure = std::current_exception();
	    }
	}
    };

    PreparedFeatureResult prepareFeatureResult(
	const BObolFeatureStoreRecord &candidate,
	BObolCommandResultStatus status, const char *command) const
    {
	PreparedFeatureResult notification;
	if (!candidate.owner.resultCallback)
	    return notification;
	notification.owner = candidate.owner;
	notification.result.status = status;
	notification.result.feature = handle(&candidate);
	notification.result.command = command ? command : "";
	notification.active = true;
	return notification;
    }

    struct PreparedFeatureChange {
	std::string key;
	BObolFeatureStoreRecord *current = NULL;
	std::unique_ptr<BObolFeatureStoreRecord> candidate;
	SoNode *node = NULL;
	bool retainedNode = false;
	bool remove = false;
	bool changed = false;
	std::map<uint64_t, BObolFeatureStoreRecord *>::node_type recordEntry;
	std::map<std::string, uint64_t>::node_type nameEntry;
    };

    struct PreparedFeaturePublication {
	explicit PreparedFeaturePublication(Impl &target) :
	    store(target), followingId(target.nextId)
	{
	}

	~PreparedFeaturePublication()
	{
	    if (committed)
		return;
	    for (PreparedFeatureChange &change : changes) {
		if (change.retainedNode && change.node)
		    change.node->unref();
	    }
	}

	Impl &store;
	std::vector<PreparedFeatureChange> changes;
	std::vector<PreparedFeatureResult> notifications;
	std::vector<std::unique_ptr<SoChildList::Replacement>> replacements;
	PreparedFeaturePresentation presentation;
	uint64_t followingId;
	bool clearOwnerGenerations = false;
	bool changed = false;
	bool committed = false;
    };

    template <typename Configure>
    SbBool publishExistingDataEdit(BObolFeatureStoreRecord *current,
	Configure configure, const char *command)
    {
	if (!current)
	    return FALSE;
	std::map<uint64_t, BObolFeatureStoreRecord *>::iterator currentEntry =
	    records.find(current->id);
	if (currentEntry == records.end() || currentEntry->second != current)
	    return FALSE;

	auto candidate = std::make_unique<BObolFeatureStoreRecord>(*current);
	configure(*candidate);
	const PreparedFeatureResult notification = prepareFeatureResult(*candidate,
	    BObolCommandResultStatus::Updated, command);
	BObolFeatureStoreRecord *installed = candidate.release();
	currentEntry->second = installed;
	delete current;
	std::exception_ptr failure;
	notification.notify(failure);
	if (failure)
	    std::rethrow_exception(failure);
	return TRUE;
    }

    using FeatureFieldValue = std::pair<SoMField *, const SoMField *>;

    SbBool publishExistingFieldEdit(
	std::unique_ptr<BObolFeatureStoreRecord> candidate,
	const std::vector<FeatureFieldValue> &fields,
	SoSFUInt32 &sourceId, const char *reason, const char *command)
    {
	if (!candidate || !candidate->node)
	    return FALSE;
	std::map<uint64_t, BObolFeatureStoreRecord *>::iterator currentEntry =
	    records.find(candidate->id);
	if (currentEntry == records.end() || !currentEntry->second)
	    return FALSE;
	BObolFeatureStoreRecord *current = currentEntry->second;
	if (current->node != candidate->node)
	    return FALSE;

	std::vector<std::unique_ptr<SoMField::ValueReplacement>> replacements;
	replacements.reserve(fields.size());
	for (const FeatureFieldValue &field : fields) {
	    if (!field.first || !field.second)
		return FALSE;
	    std::unique_ptr<SoMField::ValueReplacement> replacement =
		field.first->prepareValueReplacement(*field.second);
	    if (replacement)
		replacements.push_back(std::move(replacement));
	}

	PreparedFeaturePresentation presentation;
	prepareFeaturePresentation(TRUE, presentation, reason);
	const PreparedFeatureResult notification = prepareFeatureResult(*candidate,
	    BObolCommandResultStatus::Updated, command);
	const SbBool sourceNotifications = sourceId.enableNotify(FALSE);
	SoNode *node = candidate->node;
	node->ref();

	BObolFeatureStoreRecord *installed = candidate.release();
	currentEntry->second = installed;
	for (auto &replacement : replacements)
	    replacement->commit();
	sourceId = static_cast<uint32_t>(installed->revision);
	commitFeaturePresentation(TRUE, presentation);
	delete current;
	sourceId.enableNotify(sourceNotifications);

	std::exception_ptr failure;
	for (auto &replacement : replacements) {
	    try {
		replacement->notify();
	    } catch (...) {
		if (!failure)
		    failure = std::current_exception();
	    }
	}
	if (sourceNotifications) {
	    try {
		sourceId.touch();
	    } catch (...) {
		if (!failure)
		    failure = std::current_exception();
	    }
	}
	notifyFeaturePresentation(presentation, failure);
	notification.notify(failure);
	node->unref();
	if (failure)
	    std::rethrow_exception(failure);
	return TRUE;
    }

    template <typename Configure, typename BuildNode>
    RecordPublication publishUpsert(const SbString &name,
	BObolFeatureScope scope, BObolFeatureKind kind,
	const BObolFeatureStyle *style, const BObolFeatureOwner *owner,
	Configure configure, BuildNode buildNode,
	SbBool rollbackUnchangedPublish = FALSE)
    {
	std::unique_ptr<BObolFeatureStoreRecord> candidate =
	    prepareExistingUpsert(name, scope, kind, style, owner);
	if (candidate) {
	    configure(*candidate);
	    SoNode *node = buildNode(*candidate);
	    const uint64_t id = candidate->id;
	    const SbBool changed = publishExistingRecord(std::move(candidate),
		node, rollbackUnchangedPublish);
	    return RecordPublication(record(BObolFeatureHandle(id, 0)), changed);
	}

	if (store_string(name).empty())
	    return RecordPublication();
	const std::string key = store_key(scope, name, owner);
	uint64_t followingId = nextId;
	auto inserted = std::make_unique<BObolFeatureStoreRecord>();
	inserted->id = bobol_nonzero_identity_take(followingId);
	inserted->revision = 1;
	inserted->name = name;
	inserted->scope = scope;
	configureUpsert(*inserted, kind, style, owner);
	configure(*inserted);
	SoNode *node = buildNode(*inserted);
	return publishNewRecord(std::move(inserted), node, key, followingId);
    }

    BObolFeatureHandle publishCustomNode(const SbString &name,
	BObolFeatureScope scope, SoNode *node,
	const BObolFeatureStyle *style, const BObolFeatureOwner *owner,
	const BObolOverlayInfo *overlay)
    {
	if (!node)
	    return BObolFeatureHandle();

	BObolFeatureStoreRecord *current = recordByKey(name, scope, owner);
	BObolFeaturePublication entry;
	entry.name = name;
	entry.kind = BObolFeatureKind::CustomNode;
	entry.scope = scope;
	entry.customNode = node;
	if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	    const SoBRLMeshShape *mesh =
		static_cast<const SoBRLMeshShape *>(node);
	    entry.selectedPrimitives.reserve(
		static_cast<size_t>(mesh->selectedPrimitive.getNum()));
	    for (int i = 0; i < mesh->selectedPrimitive.getNum(); ++i)
		entry.selectedPrimitives.push_back(mesh->selectedPrimitive[i]);
	    entry.highlightedPrimitives.reserve(
		static_cast<size_t>(mesh->highlightedPrimitive.getNum()));
	    for (int i = 0; i < mesh->highlightedPrimitive.getNum(); ++i)
		entry.highlightedPrimitives.push_back(
		    mesh->highlightedPrimitive[i]);
	}
	if (current) {
	    entry.style = current->style;
	    entry.owner = current->owner;
	    entry.overlay = current->overlay;
	}
	if (style)
	    entry.style = *style;
	if (owner)
	    entry.owner = *owner;
	if (overlay)
	    entry.overlay = *overlay;

	std::vector<BObolFeaturePublication> publication;
	publication.reserve(1);
	publication.push_back(entry);
	if (!applyPublication(publication, FALSE,
		store_feature_publication_reason, store_feature_remove_command,
		store_feature_custom_node_command))
	    return BObolFeatureHandle();
	return handle(recordByKey(name, scope, &entry.owner));
    }

    SbBool preparePublication(
	const std::vector<BObolFeaturePublication> &publication,
	SbBool clearOwnerGenerations,
	const char *presentationReason,
	const char *removeCommand,
	const char *updateCommand,
	PreparedFeaturePublication &prepared)
    {
	if (&prepared.store != this)
	    return FALSE;
	prepared.clearOwnerGenerations = clearOwnerGenerations;
	if (publication.empty())
	    return TRUE;

	std::vector<std::string> keys;
	keys.reserve(publication.size());
	for (const BObolFeaturePublication &entry : publication) {
	    if (store_string(entry.name).empty())
		return FALSE;
	    if (entry.scope != BObolFeatureScope::Shared &&
		entry.scope != BObolFeatureScope::Local)
		return FALSE;
	    if (entry.action != BObolFeaturePublicationAction::Replace &&
		entry.action != BObolFeaturePublicationAction::Remove)
		return FALSE;
	    if (entry.action == BObolFeaturePublicationAction::Replace &&
		entry.kind != BObolFeatureKind::Lines &&
		entry.kind != BObolFeatureKind::HudLabel &&
		entry.kind != BObolFeatureKind::EditPreview &&
		entry.kind != BObolFeatureKind::CustomNode)
		return FALSE;
	    if (entry.action == BObolFeaturePublicationAction::Replace &&
		entry.kind == BObolFeatureKind::CustomNode && !entry.customNode)
		return FALSE;
	    if (entry.action == BObolFeaturePublicationAction::Replace &&
		entry.kind == BObolFeatureKind::CustomNode &&
		!store_custom_node_primitives_equal(entry))
		return FALSE;
	    const std::string key = store_key(entry.scope, entry.name,
		&entry.owner);
	    if (std::find(keys.begin(), keys.end(), key) != keys.end())
		return FALSE;
	    keys.push_back(key);
	}

	prepared.changes.reserve(publication.size());
	prepared.notifications.reserve(publication.size());
	for (size_t i = 0; i < publication.size(); ++i) {
		const BObolFeaturePublication &entry = publication[i];
		prepared.changes.emplace_back();
		PreparedFeatureChange &change = prepared.changes.back();
		change.key = keys[i];
		change.current = recordByKey(entry.name, entry.scope,
		    &entry.owner);
		change.remove =
		    entry.action == BObolFeaturePublicationAction::Remove;
		if (change.remove) {
		    if (!change.current)
			continue;
		    change.changed = true;
		    if (change.current->owner.resultCallback) {
			prepared.notifications.push_back(prepareFeatureResult(
			    *change.current, BObolCommandResultStatus::Removed,
			    removeCommand));
		    }
		    continue;
		}

		change.candidate = std::make_unique<BObolFeatureStoreRecord>();
		BObolFeatureStoreRecord &candidate = *change.candidate;
		if (change.current) {
		    candidate.id = change.current->id;
		    candidate.revision = change.current->revision;
		    store_revision_advance(candidate.revision);
		} else {
		    candidate.id = bobol_nonzero_identity_take(
			prepared.followingId);
		    candidate.revision = 1;
		}
		candidate.name = entry.name;
		candidate.scope = entry.scope;
		configureUpsert(candidate, entry.kind, &entry.style,
		    &entry.owner);
		candidate.overlay = entry.overlay;
		candidate.points = entry.points;
		candidate.commands = entry.kind == BObolFeatureKind::EditPreview ?
		    store_normalized_line_commands(entry.points, entry.commands) :
		    entry.commands;
		candidate.labels = entry.labels;
		candidate.selectedPrimitives = entry.selectedPrimitives;
		candidate.highlightedPrimitives = entry.highlightedPrimitives;
		if (entry.kind == BObolFeatureKind::EditPreview) {
		    candidate.identity = entry.identity;
		    candidate.editIntentId = entry.editIntentId.getLength() > 0 ?
			entry.editIntentId : entry.name;
		    candidate.editIntentRole = entry.editIntentRole.getLength() > 0 ?
			entry.editIntentRole : SbString("preview");
		    candidate.sourceRevision = entry.sourceRevision ?
			entry.sourceRevision : static_cast<uint32_t>(candidate.revision);
		    candidate.inputsRevision = entry.inputsRevision ?
			entry.inputsRevision : static_cast<uint32_t>(candidate.revision);
		}
		change.node = entry.kind == BObolFeatureKind::CustomNode ?
		    entry.customNode : store_rebuild_node_for_feature(candidate);
		if (change.current && store_feature_publications_equal(
			*change.current, candidate, change.node)) {
		    if (change.node != change.current->node) {
			change.node->ref();
			change.node->unref();
		    }
		    change.node = NULL;
		    change.candidate.reset();
		    continue;
		}

		change.changed = true;
		if (change.node != (change.current ? change.current->node : NULL)) {
		    change.node->ref();
		    change.retainedNode = true;
		}
		candidate.node = change.node;
		candidate.attachmentRoot = store_feature_attachment_root(
		    controller, &candidate);
		if (!change.current) {
		    std::map<uint64_t, BObolFeatureStoreRecord *> preparedRecords;
		    const auto preparedRecord = preparedRecords.emplace(
			candidate.id, &candidate);
		    change.recordEntry = preparedRecords.extract(
			preparedRecord.first);
		    std::map<std::string, uint64_t> preparedNames;
		    const auto preparedName = preparedNames.emplace(change.key,
			candidate.id);
		    change.nameEntry = preparedNames.extract(preparedName.first);
		}
		if (candidate.owner.resultCallback) {
		    prepared.notifications.push_back(prepareFeatureResult(candidate,
			BObolCommandResultStatus::Updated, updateCommand));
		}
	}

	prepared.changed = std::any_of(prepared.changes.begin(),
	    prepared.changes.end(), [](const PreparedFeatureChange &change) {
		return change.changed;
	    });
	if (!prepared.changed)
		return TRUE;

	std::vector<FeaturePlacement> placements;
	placements.reserve(records.size() + prepared.changes.size());
	for (const auto &recordEntry : records) {
	    BObolFeatureStoreRecord *record = recordEntry.second;
	    const PreparedFeatureChange *replacement = NULL;
	    for (const PreparedFeatureChange &change : prepared.changes) {
		if (change.changed && change.current == record) {
		    replacement = &change;
		    break;
		}
	    }
	    if (replacement) {
		if (!replacement->remove && replacement->candidate &&
		    replacement->candidate->node)
		    placements.push_back({replacement->candidate.get(),
			replacement->candidate->node,
			replacement->candidate->attachmentRoot});
	    } else if (record && record->node) {
		placements.push_back({record, record->node,
		    record->attachmentRoot});
	    }
	}
	for (const PreparedFeatureChange &change : prepared.changes) {
	    if (change.changed && !change.current && change.candidate &&
		change.candidate->node) {
		placements.push_back({change.candidate.get(),
		    change.candidate->node,
		    change.candidate->attachmentRoot});
	    }
	}

	std::vector<SoGroup *> roots;
	roots.reserve(prepared.changes.size() * 2);
	const auto appendRoot = [&roots](SoGroup *root) {
	    if (root && std::find(roots.begin(), roots.end(), root) ==
		    roots.end())
		roots.push_back(root);
	};
	for (const PreparedFeatureChange &change : prepared.changes) {
	    if (!change.changed)
		continue;
	    if (change.current)
		appendRoot(change.current->attachmentRoot);
	    if (change.candidate)
		appendRoot(change.candidate->attachmentRoot);
	}

	prepared.replacements.reserve(roots.size());
	for (SoGroup *root : roots) {
	    const std::vector<SoNode *> next = featureRootOrder(root,
		placements);
	    if (rootOrderChanged(root, next))
		prepared.replacements.push_back(
		    root->getChildren()->prepareReplacement(next));
	}
	prepareFeaturePresentation(TRUE, prepared.presentation,
	    presentationReason);
	return TRUE;
    }

    void commitPublication(PreparedFeaturePublication &prepared) noexcept
    {
	if (prepared.committed)
	    return;
	if (prepared.changed) {
	    for (auto &replacement : prepared.replacements)
		replacement->commit();
	    for (PreparedFeatureChange &change : prepared.changes) {
		if (!change.changed)
		    continue;
		if (change.remove) {
		    names.erase(change.key);
		    records.erase(change.current->id);
		} else if (change.current) {
		    records.find(change.current->id)->second =
			change.candidate.get();
		} else {
		    records.insert(std::move(change.recordEntry));
		    names.insert(std::move(change.nameEntry));
		}
	    }
	    nextId = prepared.followingId;
	    for (PreparedFeatureChange &change : prepared.changes) {
		if (!change.changed)
		    continue;
		BObolFeatureStoreRecord *installed = change.candidate.get();
		if (!change.remove)
		    change.candidate.release();
		if (change.current) {
		    if (installed &&
			change.current->node == installed->node)
			change.current->node = NULL;
		    else if (change.current->node)
			change.current->node->unref();
		    delete change.current;
		    change.current = NULL;
		}
		change.retainedNode = false;
	    }
	    commitFeaturePresentation(TRUE, prepared.presentation);
	}
	if (prepared.clearOwnerGenerations)
	    ownerGenerations.clear();
	prepared.committed = true;
    }

    void notifyPublication(PreparedFeaturePublication &prepared,
	std::exception_ptr &failure) const
    {
	if (!prepared.committed || !prepared.changed)
	    return;
	notifyRootReplacements(prepared.replacements, failure);
	notifyFeaturePresentation(prepared.presentation, failure);
	for (const PreparedFeatureResult &notification : prepared.notifications)
	    notification.notify(failure);
    }

    SbBool applyPublication(
	const std::vector<BObolFeaturePublication> &publication,
	SbBool clearOwnerGenerations = FALSE,
	const char *presentationReason = store_feature_publication_reason,
	const char *removeCommand = store_feature_remove_command,
	const char *updateCommand = store_feature_publication_command)
    {
	PreparedFeaturePublication prepared(*this);
	if (!preparePublication(publication, clearOwnerGenerations,
		presentationReason, removeCommand, updateCommand, prepared))
	    return FALSE;
	commitPublication(prepared);
	std::exception_ptr failure;
	notifyPublication(prepared, failure);
	if (failure)
	    std::rethrow_exception(failure);
	return TRUE;
    }

    static BObolFeaturePublication removalEntry(
	const BObolFeatureStoreRecord &record)
    {
	BObolFeaturePublication entry;
	entry.action = BObolFeaturePublicationAction::Remove;
	entry.name = record.name;
	entry.scope = record.scope;
	entry.owner = record.owner;
	return entry;
    }

    std::vector<BObolFeaturePublication> prepareRemovalPublication(
	unsigned int scopeMask, const BObolFeatureOwner *owner,
	const std::string *prefix = NULL) const
    {
	std::vector<BObolFeaturePublication> publication;
	for (const auto &recordEntry : records) {
	    const BObolFeatureStoreRecord *record = recordEntry.second;
	    if (!record || !(store_scope_bit(record->scope) & scopeMask) ||
		(owner && !store_owner_matches(record->owner, owner)))
		continue;
	    if (prefix) {
		const std::string name = store_string(record->name);
		if (name.compare(0, prefix->size(), *prefix) != 0)
		    continue;
	    }
	    publication.push_back(removalEntry(*record));
	}
	return publication;
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

    void requestPresentation(const char *reason)
    {
	store_revision_advance(presentationRevision);
	if (controller)
	    controller->requestPresentationRender(reason);
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

    /* SoHUDKit disables the raw GL depth test around its widgets, but a
     * renderer may subsequently reapply Coin's still-current depth element.
     * Express the HUD contract in the retained graph as well: screen lines
     * neither test nor write model depth, consistently across System GL and
     * OSMesa.  The separator scopes that override to this widget. */
    SoSeparator *widget = new SoSeparator;
    SoDepthBuffer *depth = new SoDepthBuffer;
    depth->test = FALSE;
    depth->write = FALSE;
    widget->addChild(depth);
    widget->addChild(shape);

    SoHUDKit *hud = new SoHUDKit;
    hud->addWidget(widget);
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
    this->impl->migrateController(controller);
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

uint64_t
BObolFeatureStore::presentationRevision(void) const
{
    return this->impl->presentationRevision;
}

SbBool
BObolFeatureStore::applyPublication(
    const std::vector<BObolFeaturePublication> &publication)
{
    return this->impl->applyPublication(publication);
}

SbBool
BObolFeatureStore::applyCoordinatedPublications(
    const std::vector<BObolFeatureStorePublication> &publications,
    BObolFeaturePublicationCommit committed, void *context)
{
    std::vector<BObolFeatureStore *> stores;
    stores.reserve(publications.size());
    for (const BObolFeatureStorePublication &publication : publications) {
	if (!publication.store ||
	    std::find(stores.begin(), stores.end(), publication.store) !=
		stores.end())
	    return FALSE;
	stores.push_back(publication.store);
    }

    std::vector<std::unique_ptr<Impl::PreparedFeaturePublication>> prepared;
    prepared.reserve(publications.size());
    bool changed = false;
    for (const BObolFeatureStorePublication &publication : publications) {
	auto candidate = std::make_unique<Impl::PreparedFeaturePublication>(
	    *publication.store->impl);
	if (!publication.store->impl->preparePublication(publication.features,
		FALSE, store_feature_publication_reason,
		store_feature_remove_command, store_feature_publication_command,
		*candidate))
	    return FALSE;
	changed = changed || candidate->changed;
	prepared.push_back(std::move(candidate));
    }

    if (changed && committed)
	committed(context);
    for (const auto &publication : prepared)
	publication->store.commitPublication(*publication);

    std::exception_ptr failure;
    for (const auto &publication : prepared)
	publication->store.notifyPublication(*publication, failure);
    if (failure)
	std::rethrow_exception(failure);
    return TRUE;
}

void
BObolFeatureStore::clear(void)
{
    const std::vector<BObolFeaturePublication> publication =
	this->impl->prepareRemovalPublication(BOBOL_FEATURE_SCOPE_ALL, NULL);
    (void)this->impl->applyPublication(publication, TRUE,
	store_feature_clear_reason, store_feature_clear_command);
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

    std::vector<BObolFeaturePublication> publication;
    publication.reserve(1);
    publication.push_back(this->impl->removalEntry(*rec));
    return this->impl->applyPublication(publication, FALSE,
	store_feature_removal_reason, store_feature_remove_command);
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
    const std::vector<BObolFeaturePublication> publication =
	this->impl->prepareRemovalPublication(scopeMask, owner);
    const size_t removed = publication.size();
    if (removed)
	(void)this->impl->applyPublication(publication, FALSE,
	    store_feature_removal_reason, store_feature_remove_command);
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

    const std::vector<BObolFeaturePublication> publication =
	this->impl->prepareRemovalPublication(scopeMask, owner, &p);
    const size_t removed = publication.size();
    if (removed)
	(void)this->impl->applyPublication(publication, FALSE,
	    store_feature_removal_reason, store_feature_remove_command);
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
    const auto publication = this->impl->publishUpsert(name, scope,
	BObolFeatureKind::Lines, style, owner,
	[&points, &commands](BObolFeatureStoreRecord &rec) {
	    rec.points = points;
	    rec.commands = commands;
	},
	store_rebuild_node_for_feature, TRUE);
    if (!publication.record)
	return BObolFeatureHandle();
    this->impl->notify(publication.record, BObolCommandResultStatus::Updated,
	"publishLineSet");
    return this->impl->handle(publication.record);
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

    const auto publication = this->impl->publishUpsert(name, scope,
	BObolFeatureKind::IndexedLines, style, owner,
	[&linePoints, &commands, &indices](BObolFeatureStoreRecord &rec) {
	    rec.points = linePoints;
	    rec.commands = commands;
	    rec.indices = indices;
	},
	store_rebuild_node_for_feature, TRUE);
    if (!publication.record)
	return BObolFeatureHandle();
    this->impl->notify(publication.record, BObolCommandResultStatus::Updated,
	"publishIndexedLineSet");
    return this->impl->handle(publication.record);
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
    const auto publication = this->impl->publishUpsert(name, scope,
	BObolFeatureKind::Points, style, owner,
	[&points, &commands](BObolFeatureStoreRecord &rec) {
	    rec.points = points;
	    rec.commands = commands;
	},
	store_rebuild_node_for_feature, TRUE);
    if (!publication.record)
	return BObolFeatureHandle();
    this->impl->notify(publication.record, BObolCommandResultStatus::Updated,
	"publishPointSet");
    return this->impl->handle(publication.record);
}

BObolFeatureHandle
BObolFeatureStore::publishLabels(const SbString &name,
				   BObolFeatureScope scope,
				   const std::vector<BObolLabel> &labels,
				   const BObolFeatureStyle *style,
				   const BObolFeatureOwner *owner)
{
    const auto publication = this->impl->publishUpsert(name, scope,
	BObolFeatureKind::Labels, style, owner,
	[&labels](BObolFeatureStoreRecord &rec) { rec.labels = labels; },
	store_rebuild_node_for_feature, TRUE);
    if (!publication.record)
	return BObolFeatureHandle();
    this->impl->notify(publication.record, BObolCommandResultStatus::Updated,
	"publishLabels");
    return this->impl->handle(publication.record);
}

BObolFeatureHandle
BObolFeatureStore::publishHudLabels(const SbString &name,
				      BObolFeatureScope scope,
				      const std::vector<BObolLabel> &labels,
				      const BObolFeatureStyle *style,
				      const BObolFeatureOwner *owner)
{
    const auto publication = this->impl->publishUpsert(name, scope,
	BObolFeatureKind::HudLabel, style, owner,
	[&labels](BObolFeatureStoreRecord &rec) { rec.labels = labels; },
	store_rebuild_node_for_feature, TRUE);
    if (!publication.record)
	return BObolFeatureHandle();
    this->impl->notify(publication.record, BObolCommandResultStatus::Updated,
	"publishHudLabels");
    return this->impl->handle(publication.record);
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

    const auto publication = this->impl->publishUpsert(name, scope,
	BObolFeatureKind::Arrow, &arrowStyle, owner,
	[&points](BObolFeatureStoreRecord &rec) {
	    rec.points = points;
	    rec.commands.clear();
	},
	store_rebuild_node_for_feature, TRUE);
    if (!publication.record)
	return BObolFeatureHandle();
    this->impl->notify(publication.record, BObolCommandResultStatus::Updated,
	"publishArrow");
    return this->impl->handle(publication.record);
}

BObolFeatureHandle
BObolFeatureStore::publishAxes(const SbString &name,
				 BObolFeatureScope scope,
				 const std::vector<SbVec3f> &centers,
				 float halfAxesSize,
				 const BObolFeatureStyle *style,
				 const BObolFeatureOwner *owner)
{
    const auto publication = this->impl->publishUpsert(name, scope,
	BObolFeatureKind::Axes, style, owner,
	[&centers, halfAxesSize](BObolFeatureStoreRecord &rec) {
	    rec.axesCenters = centers;
	    rec.halfAxesSize = halfAxesSize;
	},
	store_rebuild_node_for_feature, TRUE);
    if (!publication.record)
	return BObolFeatureHandle();
    this->impl->notify(publication.record, BObolCommandResultStatus::Updated,
	"publishAxes");
    return this->impl->handle(publication.record);
}

BObolFeatureHandle
BObolFeatureStore::publishLineLayers(const SbString &name,
				       BObolFeatureScope scope,
				       const std::vector<BObolLineLayer> &layers,
				       const BObolFeatureStyle *style,
				       const BObolFeatureOwner *owner)
{
    const auto publication = this->impl->publishUpsert(name, scope,
	BObolFeatureKind::LineLayer, style, owner,
	[&layers](BObolFeatureStoreRecord &rec) {
	    rec.layers = layers;
	    rec.points.clear();
	    rec.commands.clear();
	    for (const BObolLineLayer &layer : rec.layers) {
		rec.points.insert(rec.points.end(), layer.points.begin(),
		    layer.points.end());
		rec.commands.insert(rec.commands.end(), layer.commands.begin(),
		    layer.commands.end());
	    }
	},
	store_rebuild_node_for_feature, TRUE);
    if (!publication.record)
	return BObolFeatureHandle();
    this->impl->notify(publication.record, BObolCommandResultStatus::Updated,
	"publishLineLayers");
    return this->impl->handle(publication.record);
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

    const auto publication = this->impl->publishUpsert(name, scope,
	BObolFeatureKind::LineLayer, style, owner,
	[&name, builder](BObolFeatureStoreRecord &rec) {
	    rec.layers = store_line_layers_from_builder(name, builder);
	    rec.points.clear();
	    rec.commands.clear();
	    for (const BObolLineLayer &layer : rec.layers) {
		rec.points.insert(rec.points.end(), layer.points.begin(),
		    layer.points.end());
		rec.commands.insert(rec.commands.end(), layer.commands.begin(),
		    layer.commands.end());
	    }
	},
	[&name, builder](const BObolFeatureStoreRecord &rec) -> SoNode * {
	    SoBRLLineLayerOverlay *overlay = new SoBRLLineLayerOverlay;
	    overlay->overlayId = name;
	    overlay->sourceId = static_cast<uint32_t>(rec.revision);
	    overlay->selectable = rec.style.hasSelectable ?
		rec.style.selectable : TRUE;
	    overlay->rebuildGeometry(builder);
	    return overlay;
	});
    if (!publication.record)
	return BObolFeatureHandle();
    this->impl->notify(publication.record, BObolCommandResultStatus::Updated,
	"publishLineLayerBuilder");
    return this->impl->handle(publication.record);
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
    const auto publication = this->impl->publishUpsert(name, scope,
	BObolFeatureKind::IndexedFaceSet, style, owner,
	[&points, &normals, &indices](BObolFeatureStoreRecord &rec) {
	    rec.points = points;
	    rec.normals = normals;
	    rec.indices = indices;
	    rec.commands.clear();
	},
	store_rebuild_node_for_feature, TRUE);
    if (!publication.record)
	return BObolFeatureHandle();
    this->impl->notify(publication.record, BObolCommandResultStatus::Updated,
	"publishIndexedFaceSet");
    return this->impl->handle(publication.record);
}

SbBool
BObolFeatureStore::updateIndexedFaceSetPoints(
    BObolFeatureHandle handle,
    const std::vector<int32_t> &pointIndices,
    const std::vector<SbVec3f> &points)
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec || rec->kind != BObolFeatureKind::IndexedFaceSet ||
	pointIndices.empty() || pointIndices.size() != points.size() ||
	!rec->node ||
	!rec->node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return FALSE;

    /* Validate the complete patch before changing either the record or its
     * retained node.  A rejected edit presentation update must be atomic. */
    for (size_t i = 0; i < pointIndices.size(); i++) {
	if (pointIndices[i] < 0 ||
	    static_cast<size_t>(pointIndices[i]) >= rec->points.size())
	    return FALSE;
    }

    SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(rec->node);
    if (shape->point.getNum() != static_cast<int>(rec->points.size()))
	return FALSE;

    auto candidate = std::make_unique<BObolFeatureStoreRecord>(*rec);
    for (size_t i = 0; i < pointIndices.size(); i++)
	candidate->points[static_cast<size_t>(pointIndices[i])] = points[i];

    /* A vertex edit invalidates supplied normals.  Feature meshes normally
     * omit them and let SoBRLMeshShape derive face normals, but clearing here
     * keeps the generic patch API correct for callers that did supply them. */
    candidate->normals.clear();
    store_revision_advance(candidate->revision);

    SoMFVec3f nextPoints;
    nextPoints.setValues(0, static_cast<int>(candidate->points.size()),
	candidate->points.data());
    SoMFVec3f nextNormals;
    const std::vector<BObolFeatureStore::Impl::FeatureFieldValue> fields = {
	{&shape->point, &nextPoints}, {&shape->normal, &nextNormals}};
    return this->impl->publishExistingFieldEdit(std::move(candidate), fields,
	shape->sourceId, store_feature_indexed_points_reason,
	"updateIndexedFaceSetPoints");
}

BObolFeatureHandle
BObolFeatureStore::publishCustomNode(const SbString &name,
				       BObolFeatureScope scope,
				       SoNode *node,
				       const BObolFeatureStyle *style,
				       const BObolFeatureOwner *owner)
{
    return this->impl->publishCustomNode(name, scope, node, style, owner, NULL);
}

BObolFeatureHandle
BObolFeatureStore::publishCustomNode(const SbString &name,
				       BObolFeatureScope scope,
				       SoNode *node,
				       const BObolOverlayInfo &overlay,
				       const BObolFeatureStyle *style,
				       const BObolFeatureOwner *owner)
{
    return this->impl->publishCustomNode(name, scope, node, style, owner,
	&overlay);
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
    const auto publication = this->impl->publishUpsert(name,
	BObolFeatureScope::Local, BObolFeatureKind::EditPreview, NULL, owner,
	[&](BObolFeatureStoreRecord &rec) {
	    rec.identity = identity;
	    rec.editIntentId = editIntentId.getLength() > 0 ?
		editIntentId : name;
	    rec.editIntentRole = editIntentRole.getLength() > 0 ?
		editIntentRole : SbString("preview");
	    rec.points = points;
	    rec.commands = store_normalized_line_commands(points, commands);
	    rec.sourceRevision = sourceRevision ? sourceRevision :
		static_cast<uint32_t>(rec.revision);
	    rec.inputsRevision = inputsRevision ? inputsRevision :
		static_cast<uint32_t>(rec.revision);
	},
	store_rebuild_node_for_feature);
    if (!publication.record)
	return BObolFeatureHandle();
    this->impl->notify(publication.record, BObolCommandResultStatus::Updated,
	"publishEditPreview");
    return this->impl->handle(publication.record);
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

    const auto publication = this->impl->publishUpsert(name,
	BObolFeatureScope::Local, BObolFeatureKind::EditPreview, &style, owner,
	[&](BObolFeatureStoreRecord &rec) {
	    rec.identity = summary.sourceInstanceKey.getLength() > 0 ?
		summary.sourceInstanceKey : summary.path;
	    rec.editIntentId = editIntentId.getLength() > 0 ?
		editIntentId : name;
	    rec.editIntentRole = editIntentRole.getLength() > 0 ?
		editIntentRole : SbString("compact-instance");
	    rec.points = points;
	    rec.commands = store_normalized_line_commands(points, commands);
	    rec.sourceRevision = source.sourceRevision.getValue();
	    rec.inputsRevision = source.inputsRevision.getValue();
	    rec.compactEdit = TRUE;
	    rec.compactSummary = summary;
	},
	store_rebuild_node_for_feature);
    if (!publication.record)
	return BObolFeatureHandle();

    this->impl->notify(publication.record, BObolCommandResultStatus::Accepted,
	"promoteCompactInstanceForEdit");
    return this->impl->handle(publication.record);
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

    const auto publication = this->impl->publishExistingEdit(rec,
	[&](BObolFeatureStoreRecord &candidate) {
	    if (identity.getLength() > 0)
		candidate.identity = identity;
	    candidate.points = points;
	    candidate.commands = store_normalized_line_commands(points, commands);
	    store_revision_advance(candidate.revision);
	    candidate.sourceRevision = sourceRevision ? sourceRevision :
		static_cast<uint32_t>(candidate.revision);
	    candidate.inputsRevision = inputsRevision ? inputsRevision :
		static_cast<uint32_t>(candidate.revision);
	},
	store_rebuild_node_for_feature);
    this->impl->notify(publication.record, BObolCommandResultStatus::Updated,
		       "replaceEditPreviewGeometry");
    return TRUE;
}

SbBool
BObolFeatureStore::appendLinePoint(BObolFeatureHandle handle,
				     const SbVec3f &point)
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec || rec->kind == BObolFeatureKind::CustomNode)
	return FALSE;

    const auto publication = this->impl->publishExistingEdit(rec,
	[&point](BObolFeatureStoreRecord &candidate) {
	    candidate.points.push_back(point);
	    candidate.commands.push_back(
		static_cast<int32_t>(BObolLineCommand::Draw));
	    store_revision_advance(candidate.revision);
	},
	store_node_for_feature);
    this->impl->notify(publication.record, BObolCommandResultStatus::Updated,
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

    const auto publication = this->impl->publishExistingEdit(rec,
	[&labels](BObolFeatureStoreRecord &candidate) {
	    candidate.labels = labels;
	    store_revision_advance(candidate.revision);
	},
	store_node_for_feature);
    this->impl->notify(publication.record, BObolCommandResultStatus::Updated,
		       "replaceLabels");
    return TRUE;
}

SbBool
BObolFeatureStore::clearGeometry(BObolFeatureHandle handle)
{
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

    const auto publication = this->impl->publishExistingEdit(rec,
	[](BObolFeatureStoreRecord &candidate) {
	    candidate.points.clear();
	    candidate.commands.clear();
	    candidate.indices.clear();
	    candidate.normals.clear();
	    candidate.labels.clear();
	    candidate.axesCenters.clear();
	    candidate.layers.clear();
	    store_revision_advance(candidate.revision);
	},
	[](const BObolFeatureStoreRecord &) -> SoNode * { return NULL; });
    this->impl->notify(publication.record, BObolCommandResultStatus::Updated,
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

    const auto publication = this->impl->publishExistingEdit(rec,
	[&style](BObolFeatureStoreRecord &candidate) {
	    if (style.hasVisible) {
		candidate.style.hasVisible = TRUE;
		candidate.style.visible = style.visible;
	    }
	    if (style.hasSelectable) {
		candidate.style.hasSelectable = TRUE;
		candidate.style.selectable = style.selectable;
	    }
	    if (style.hasColor) {
		candidate.style.hasColor = TRUE;
		candidate.style.color = style.color;
	    }
	    if (style.hasLineWidth) {
		candidate.style.hasLineWidth = TRUE;
		candidate.style.lineWidth = style.lineWidth;
	    }
	    if (style.hasLineStyle) {
		candidate.style.hasLineStyle = TRUE;
		candidate.style.lineStyle = style.lineStyle;
	    }
	    if (style.hasTransparency) {
		candidate.style.hasTransparency = TRUE;
		candidate.style.transparency = style.transparency;
	    }
	    if (style.hasArrow) {
		candidate.style.hasArrow = TRUE;
		candidate.style.arrow = style.arrow;
	    }
	    if (style.hasArrowTip) {
		candidate.style.hasArrowTip = TRUE;
		candidate.style.arrowTipLength = style.arrowTipLength;
		candidate.style.arrowTipWidth = style.arrowTipWidth;
	    }
	    store_revision_advance(candidate.revision);
	},
	store_rebuild_node_for_feature);
    this->impl->notify(publication.record, BObolCommandResultStatus::Updated,
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
    if (store_overlay_equal(rec->overlay, overlay))
	return TRUE;

    auto candidate = std::make_unique<BObolFeatureStoreRecord>(*rec);
    candidate->overlay = overlay;
    store_revision_advance(candidate->revision);
    SoNode *node = store_rebuild_node_for_feature(*candidate);
    const uint64_t id = candidate->id;
    const SbBool nodeChanged = this->impl->publishExistingRecord(
	std::move(candidate), node);
    rec = this->impl->record(BObolFeatureHandle(id, 0));
    if (!nodeChanged)
	this->impl->requestPresentation("view-feature-overlay-order");
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

    const BObolOverlayInfo empty;
    if (store_overlay_equal(rec->overlay, empty))
	return TRUE;

    auto candidate = std::make_unique<BObolFeatureStoreRecord>(*rec);
    candidate->overlay = empty;
    store_revision_advance(candidate->revision);
    SoNode *node = store_rebuild_node_for_feature(*candidate);
    const uint64_t id = candidate->id;
    const SbBool nodeChanged = this->impl->publishExistingRecord(
	std::move(candidate), node);
    rec = this->impl->record(BObolFeatureHandle(id, 0));
    if (!nodeChanged)
	this->impl->requestPresentation("view-feature-overlay-order");
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

    return this->impl->publishExistingDataEdit(rec,
	[&metadata](BObolFeatureStoreRecord &candidate) {
	    candidate.metadata = metadata;
	    store_revision_advance(candidate.revision);
	}, "replaceMetadata");
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

    for (std::vector<BObolFeaturePrimitiveMetadata>::const_iterator it =
	     rec->primitiveMetadata.begin();
	 it != rec->primitiveMetadata.end(); ++it) {
	if (it->primitiveIndex != primitiveIndex)
	    continue;
	return this->impl->publishExistingDataEdit(rec,
	    [primitiveIndex, &metadata](BObolFeatureStoreRecord &candidate) {
		auto current = std::find_if(candidate.primitiveMetadata.begin(),
		    candidate.primitiveMetadata.end(),
		    [primitiveIndex](const BObolFeaturePrimitiveMetadata &item) {
			return item.primitiveIndex == primitiveIndex;
		    });
		if (metadata.empty())
		    candidate.primitiveMetadata.erase(current);
		else
		    current->metadata = metadata;
		store_revision_advance(candidate.revision);
	    }, "replacePrimitiveMetadata");
    }

    if (!metadata.empty()) {
	return this->impl->publishExistingDataEdit(rec,
	    [primitiveIndex, &metadata](BObolFeatureStoreRecord &candidate) {
		BObolFeaturePrimitiveMetadata item;
		item.primitiveIndex = primitiveIndex;
		item.metadata = metadata;
		candidate.primitiveMetadata.push_back(item);
		store_revision_advance(candidate.revision);
	    }, "replacePrimitiveMetadata");
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

    if (rec->node &&
	rec->node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(rec->node);
	auto candidate = std::make_unique<BObolFeatureStoreRecord>(*rec);
	candidate->selectedPrimitives = primitives;
	store_revision_advance(candidate->revision);
	SoMFInt32 next;
	if (!primitives.empty())
	    next.setValues(0, static_cast<int>(primitives.size()),
		primitives.data());
	const std::vector<BObolFeatureStore::Impl::FeatureFieldValue> fields = {
	    {&shape->selectedPrimitive, &next}};
	return this->impl->publishExistingFieldEdit(std::move(candidate), fields,
	    shape->sourceId, store_feature_selected_primitives_reason,
	    "replaceSelectedPrimitives");
    }
    const auto publication = this->impl->publishExistingEdit(rec,
	[&primitives](BObolFeatureStoreRecord &candidate) {
	    candidate.selectedPrimitives = primitives;
	    store_revision_advance(candidate.revision);
	},
	store_rebuild_node_for_feature);
    this->impl->notify(publication.record, BObolCommandResultStatus::Updated,
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

    if (rec->node &&
	rec->node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(rec->node);
	auto candidate = std::make_unique<BObolFeatureStoreRecord>(*rec);
	candidate->highlightedPrimitives = primitives;
	store_revision_advance(candidate->revision);
	SoMFInt32 next;
	if (!primitives.empty())
	    next.setValues(0, static_cast<int>(primitives.size()),
		primitives.data());
	const std::vector<BObolFeatureStore::Impl::FeatureFieldValue> fields = {
	    {&shape->highlightedPrimitive, &next}};
	return this->impl->publishExistingFieldEdit(std::move(candidate), fields,
	    shape->sourceId, store_feature_highlighted_primitives_reason,
	    "replaceHighlightedPrimitives");
    }
    const auto publication = this->impl->publishExistingEdit(rec,
	[&primitives](BObolFeatureStoreRecord &candidate) {
	    candidate.highlightedPrimitives = primitives;
	    store_revision_advance(candidate.revision);
	},
	store_rebuild_node_for_feature);
    this->impl->notify(publication.record, BObolCommandResultStatus::Updated,
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
    this->impl->publishExistingEdit(rec,
	[](BObolFeatureStoreRecord &) {},
	store_rebuild_node_for_feature);
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
BObolFeatureStore::summary(BObolFeatureHandle handle,
			     BObolFeatureSummary &summaryOut) const
{
    summaryOut = BObolFeatureSummary();
    BObolFeatureStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;

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

void
BObolFeatureStore::visitNodes(BObolFeatureNodeCallback callback,
	void *userData) const
{
    if (!callback)
	return;
    for (std::map<uint64_t, BObolFeatureStoreRecord *>::const_iterator it =
	    this->impl->records.begin(); it != this->impl->records.end(); ++it) {
	BObolFeatureStoreRecord *rec = it->second;
	if (!rec || !rec->node)
	    continue;
	if (!callback(this->impl->handle(rec), rec->node, userData))
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
    SbBool selected;
    SbBool visible;
    BObolPolygonVisual visual;
    long currentContour;
    long currentPoint;
    SbVec3f originPoint;
    plane_t viewPlane;
    SbString sketchName;
    void *userData;
    struct bg_polygon polygon;
    SoNode *node;

    BObolPolygonStoreRecord(void) :
	id(0),
	revision(0),
	name(""),
	scope(BObolFeatureScope::Shared),
	type(BObolPolygonType::General),
	selected(FALSE),
	visible(TRUE),
	visual(),
	currentContour(-1),
	currentPoint(-1),
	originPoint(0.0f, 0.0f, 0.0f),
	sketchName(""),
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
    uint64_t presentationRevision;
    uint64_t nextId;
    std::map<uint64_t, BObolPolygonStoreRecord *> records;
    std::map<std::string, uint64_t> names;
    BObolPolygonHandle snapExclude;

    Impl(void) : controller(NULL),
	referenceGeneration(store_reference_generation_next()),
	presentationRevision(0), nextId(1),
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
	requestPresentation("view-polygon-store");
    }

    void requestPresentation(const char *reason)
    {
	store_revision_advance(presentationRevision);
	if (controller)
	    controller->requestPresentationRender(reason);
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
    shape->selected = rec.selected;

    SoBRLVListShape *handles = NULL;
    if (rec.type == BObolPolygonType::General &&
	(rec.selected || (rec.currentContour >= 0 && rec.currentPoint >= 0))) {
	std::vector<SbVec3f> handlePoints;
	std::vector<double> preciseHandlePoints;
	std::vector<int32_t> handleCommands;
	std::vector<int> colorValid;
	std::vector<SbColor> colors;
	std::vector<int> scaleValid;
	std::vector<float> scales;
	std::vector<int> normalValid;
	std::vector<SbVec3f> normals;
	for (size_t i = 0; i < rec.polygon.num_contours; i++) {
	    const struct bg_poly_contour &contour = rec.polygon.contour[i];
	    for (size_t j = 0; j < contour.num_points; j++) {
		const bool active = rec.currentContour == static_cast<long>(i) &&
		    rec.currentPoint == static_cast<long>(j);
		handlePoints.push_back(store_vec3(contour.point[j]));
		preciseHandlePoints.push_back(contour.point[j][X]);
		preciseHandlePoints.push_back(contour.point[j][Y]);
		preciseHandlePoints.push_back(contour.point[j][Z]);
		handleCommands.push_back(SoBRLVListShape::POINT);
		colorValid.push_back(1);
		colors.push_back(active ? SbColor(1.0f, 1.0f, 0.0f) :
		    SbColor(1.0f, 1.0f, 1.0f));
		scaleValid.push_back(1);
		scales.push_back(active ? 8.0f : 4.0f);
		normalValid.push_back(0);
		normals.push_back(SbVec3f(0.0f, 0.0f, 1.0f));
	    }
	}
	if (!handlePoints.empty()) {
	    handles = new SoBRLVListShape;
	    handles->setLineSet(handlePoints.data(), handleCommands.data(),
		static_cast<int>(handlePoints.size()));
	    handles->setPrecisePoints(preciseHandlePoints.data(),
		static_cast<int>(handlePoints.size()));
	    handles->setPointAttributes(colorValid.data(), colors.data(),
		scaleValid.data(), scales.data(), normalValid.data(),
		normals.data(), static_cast<int>(handlePoints.size()));
	    handles->sourcePath = rec.name;
	    handles->sourceName = rec.name;
	    handles->sourceType = "view-polygon-handle";
	    handles->displayName = rec.name;
	    handles->geometryName = rec.name;
	    handles->sourceIdentity = rec.name;
	    handles->cacheIdentity = rec.name;
	    handles->overlayIntent = TRUE;
	    handles->nonDatabaseSource = TRUE;
	    handles->localSource = rec.scope == BObolFeatureScope::Local;
	    handles->sharedSource = rec.scope == BObolFeatureScope::Shared;
	    handles->recordRole = "view-polygon";
	    handles->geometryKind = "point";
	    handles->sourceId = static_cast<uint32_t>(rec.revision);
	}
    }

    SoSeparator *sep = new SoSeparator;
    SoNode *meshFillNode = store_polygon_mesh_fill_node(rec);
    if (meshFillNode)
	sep->addChild(meshFillNode);
    SoNode *hatchFillNode = store_polygon_hatch_fill_node(rec);
    if (hatchFillNode)
	sep->addChild(hatchFillNode);
    sep->addChild(shape);
    if (handles)
	sep->addChild(handles);
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

uint64_t
BObolPolygonStore::presentationRevision(void) const
{
    return this->impl->presentationRevision;
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
    const bool changed = !this->impl->records.empty();
    this->impl->clear();
    if (changed)
	this->impl->requestPresentation("view-polygon-clear");
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
    rec->id = bobol_nonzero_identity_take(this->impl->nextId);
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
    store_revision_advance(dst->revision);
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
    } else if (update == BObolPolygonUpdate::PointDelete) {
	if (rec->type != BObolPolygonType::General ||
	    rec->currentContour < 0 || rec->currentPoint < 0 ||
	    static_cast<size_t>(rec->currentContour) >=
		rec->polygon.num_contours)
	    return FALSE;
	struct bg_poly_contour &contour =
	    rec->polygon.contour[rec->currentContour];
	if (static_cast<size_t>(rec->currentPoint) >= contour.num_points)
	    return FALSE;
	/* Keep every interactive contour valid: an open contour represents at
	 * least one line segment, and a closed contour at least one triangle. */
	const size_t minimum = contour.open ? 2 : 3;
	if (contour.num_points <= minimum)
	    return FALSE;
	const size_t oldPoint = static_cast<size_t>(rec->currentPoint);
	if (bg_polygon_remove_point(&rec->polygon,
		static_cast<size_t>(rec->currentContour), oldPoint))
	    return FALSE;
	const struct bg_poly_contour &updated =
	    rec->polygon.contour[rec->currentContour];
	rec->currentPoint = static_cast<long>(
	    oldPoint < updated.num_points ? oldPoint : updated.num_points - 1);
    }
    store_revision_advance(rec->revision);
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

    store_revision_advance(rec->revision);
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
    store_revision_advance(rec->revision);
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
    store_revision_advance(rec->revision);
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

    if (this->impl->snapExclude.isValid() &&
	this->impl->snapExclude.id == handle.id)
	this->impl->snapExclude = BObolPolygonHandle();

    this->impl->names.erase(store_key(rec->scope, rec->name));
    store_release_node(this->impl->controller, rec->node);
    bg_polygon_clear(&rec->polygon);
    this->impl->records.erase(rec->id);
    delete rec;
    this->impl->requestPresentation("view-polygon-remove");
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
    recordOut.selected = rec->selected;
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
    recordOut.sketchName = rec->sketchName;
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
    const SbBool presentationChanged = rec->currentPoint >= 0 || point >= 0;
    rec->currentContour = contour;
    rec->currentPoint = point;
    store_revision_advance(rec->revision);
    if (presentationChanged)
	this->impl->realize(rec);
    return TRUE;
}

SbBool
BObolPolygonStore::setSelected(BObolPolygonHandle handle, SbBool selected)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    const SbBool normalized = selected ? TRUE : FALSE;
    if (rec->selected == normalized)
	return TRUE;
    rec->selected = normalized;
    store_revision_advance(rec->revision);
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BObolPolygonStore::clearSelection(void)
{
    SbBool changed = FALSE;
    for (std::map<uint64_t, BObolPolygonStoreRecord *>::iterator it =
	    this->impl->records.begin(); it != this->impl->records.end(); ++it) {
	BObolPolygonStoreRecord *rec = it->second;
	if (!rec || !rec->selected)
	    continue;
	rec->selected = FALSE;
	store_revision_advance(rec->revision);
	this->impl->realize(rec);
	changed = TRUE;
    }
    return changed;
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
    store_revision_advance(rec->revision);
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
    store_revision_advance(rec->revision);
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BObolPolygonStore::clearSelectedPoint(BObolPolygonHandle handle)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    if (rec->currentContour < 0 && rec->currentPoint < 0)
	return TRUE;
    rec->currentContour = -1;
    rec->currentPoint = -1;
    store_revision_advance(rec->revision);
    this->impl->realize(rec);
    return TRUE;
}

SbBool
BObolPolygonStore::clearAllPointSelections(void)
{
    for (std::map<uint64_t, BObolPolygonStoreRecord *>::iterator it =
	     this->impl->records.begin(); it != this->impl->records.end(); ++it) {
	if (it->second) {
	    if (it->second->currentContour < 0 &&
		it->second->currentPoint < 0)
		continue;
	    it->second->currentContour = -1;
	    it->second->currentPoint = -1;
	    store_revision_advance(it->second->revision);
	    this->impl->realize(it->second);
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
    store_revision_advance(rec->revision);
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
    store_revision_advance(rec->revision);
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
    store_revision_advance(rec->revision);
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
    store_revision_advance(rec->revision);
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
    store_revision_advance(rec->revision);
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
    store_revision_advance(rec->revision);
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
    store_revision_advance(rec->revision);
    this->impl->realize(rec);
    return TRUE;
}

const struct bg_polygon *
BObolPolygonStore::geometry(BObolPolygonHandle handle) const {
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    return rec ? &rec->polygon : NULL;
}

SbBool
BObolPolygonStore::setSketchName(BObolPolygonHandle handle,
	const SbString &sketchName)
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    if (!rec)
	return FALSE;
    rec->sketchName = sketchName;
    store_revision_advance(rec->revision);
    return TRUE;
}

const char *
BObolPolygonStore::sketchName(BObolPolygonHandle handle) const
{
    BObolPolygonStoreRecord *rec = this->impl->record(handle);
    return rec ? rec->sketchName.getString() : NULL;
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
    store_revision_advance(rec->revision);
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
	store_revision_advance(rt->revision);
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
    store_revision_advance(rt->revision);
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
    rec->id = bobol_nonzero_identity_take(this->impl->nextId);
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
    rec->sketchName = dp->d_namep;
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

SbBool
BObolPolygonStore::updateSketch(BObolPolygonHandle handle,
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

    struct directory *dp = db_sketch_polygon_data_update_sketch(dbip,
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
    if (!handle.isValid()) {
	this->impl->snapExclude = BObolPolygonHandle();
	return TRUE;
    }
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

SbBool
BObolSelectionStore::applyPathDelta(
    const std::vector<SbString> &addedPaths,
    const std::vector<SbString> &removedPaths,
    int kind,
    const BObolFeatureOwner *owner)
{
    const int recordKind = store_selection_record_kind(kind);
    std::unordered_set<std::string> removed;
    removed.reserve(removedPaths.size());
    for (const SbString &path : removedPaths) {
	const std::string value = store_string(path);
	if (!value.empty())
	    removed.insert(value);
    }

    std::unordered_set<std::string> existing;
    existing.reserve(this->impl->records.size() + addedPaths.size());
    std::vector<BObolSelectionRecord> next;
    next.reserve(this->impl->records.size() + addedPaths.size());
    for (const BObolSelectionRecord &record : this->impl->records) {
	const bool target = store_owner_matches(record.owner, owner) &&
	    store_selection_kind_matches(record.kind, recordKind);
	const std::string path = store_string(record.path);
	if (target && removed.find(path) != removed.end())
	    continue;
	if (target && !path.empty())
	    existing.insert(path);
	next.push_back(record);
    }

    for (const SbString &pathValue : addedPaths) {
	const std::string path = store_string(pathValue);
	if (path.empty() || !existing.insert(path).second)
	    continue;
	BObolSelectionRecord record;
	record.path = pathValue;
	record.kind = recordKind;
	if (owner)
	    record.owner = *owner;
	next.push_back(record);
    }
    this->impl->records.swap(next);
    return TRUE;
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
