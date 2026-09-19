/*                 L I N E _ L A Y E R _ O V E R L A Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BLineLayerOverlay.h"
#include "BObol/BLodRealization.h"
#include "BObol/BVListShape.h"

#include "bg/line_layer.h"

#include <Inventor/misc/SoChildList.h>
#include <Inventor/nodes/SoDepthBuffer.h>
#include <Inventor/tools/SbModernUtils.h>

#include <limits>
#include <utility>
#include <vector>

static constexpr unsigned char DEFAULT_COLOR_COMPONENT =
    std::numeric_limits<unsigned char>::max();
static constexpr float COLOR_COMPONENT_SCALE =
    1.0f / static_cast<float>(DEFAULT_COLOR_COMPONENT);

SO_NODE_SOURCE(SoBRLLineLayerOverlay);

SoBRLLineLayerOverlay::SoBRLLineLayerOverlay(void)
{
    SO_NODE_CONSTRUCTOR(SoBRLLineLayerOverlay);

    SO_NODE_ADD_FIELD(overlayId, ("overlay::line-layer"));
    SO_NODE_ADD_FIELD(sourceId, (0));
    SO_NODE_ADD_FIELD(visible, (TRUE));
    SO_NODE_ADD_FIELD(selectable, (TRUE));
    SO_NODE_ADD_FIELD(depthTest, (TRUE));
}

SoBRLLineLayerOverlay::~SoBRLLineLayerOverlay(void)
{
}

void
SoBRLLineLayerOverlay::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLLineLayerOverlay, SoSeparator, "Separator");
}

static int
convert_line_layer_command(int command)
{
    switch (command) {
	case BG_GEOMETRY_LINE_MOVE:
	    return SoBRLVListShape::MOVE;
	case BG_GEOMETRY_LINE_DRAW:
	    return SoBRLVListShape::DRAW;
	case BG_GEOMETRY_POINT_DRAW:
	    return SoBRLVListShape::POINT;
	default:
	    break;
    }
    return -1;
}

static SbString
line_layer_path(const SbString &overlayId, const char *layerName)
{
    SbString path = overlayId;
    if (path.getLength() > 0 && layerName && layerName[0]) {
	path += "/";
	path += layerName;
    }
    return path;
}

int
SoBRLLineLayerOverlay::rebuildGeometry(const struct bg_line_layer_builder *builder)
{
    if (!this->visible.getValue() || !builder) {
	if (!this->getNumChildren())
	    return 0;
	auto replacement = this->getChildren()->prepareReplacement({});
	replacement->commit();
	replacement->notify();
	return 0;
    }

    const size_t layerCount = bg_line_layer_builder_layer_count(builder);
    std::vector<SbModernUtils::SoNodeRef> owners;
    std::vector<SoNode *> publishedChildren;
    owners.reserve(layerCount + 1);
    publishedChildren.reserve(layerCount + 1);

    SbModernUtils::SoNodeRef depthOwner(new SoDepthBuffer);
    auto *depth = static_cast<SoDepthBuffer *>(depthOwner.get());
    depth->test = this->depthTest.getValue();
    depth->write = this->depthTest.getValue();
    depth->function = SoDepthBuffer::LEQUAL;
    publishedChildren.push_back(depth);
    owners.push_back(std::move(depthOwner));

    int realized = 0;
    for (size_t layerIndex = 0; layerIndex < layerCount; layerIndex++) {
	const struct bg_line_layer *layer =
	    bg_line_layer_builder_layer_at(builder, layerIndex);
	if (!layer)
	    continue;
	const size_t pointCount = bg_line_layer_point_count(layer);
	if (pointCount == 0)
	    continue;

	std::vector<SbVec3f> points;
	std::vector<int32_t> commands;
	points.reserve(pointCount);
	commands.reserve(pointCount);

	for (size_t i = 0; i < pointCount; i++) {
	    point_t p;
	    int command = -1;
	    if (!bg_line_layer_point_at(layer, i, p) ||
		!bg_line_layer_command_at(layer, i, &command))
		continue;

	    int shapeCommand = convert_line_layer_command(command);
	    if (shapeCommand < 0)
		continue;

	    points.push_back(SbVec3f(static_cast<float>(p[X]),
				     static_cast<float>(p[Y]),
				     static_cast<float>(p[Z])));
	    commands.push_back(shapeCommand);
	}

	if (points.empty() || points.size() != commands.size())
	    continue;

	unsigned char r = DEFAULT_COLOR_COMPONENT;
	unsigned char g = DEFAULT_COLOR_COMPONENT;
	unsigned char b = DEFAULT_COLOR_COMPONENT;
	(void)bg_line_layer_color(layer, &r, &g, &b);

	const char *layerName = bg_line_layer_name(layer);
	SbModernUtils::SoNodeRef shapeOwner(new SoBRLVListShape);
	auto *shape = static_cast<SoBRLVListShape *>(shapeOwner.get());
	shape->sourcePath = line_layer_path(this->overlayId.getValue(), layerName);
	shape->sourceName = layerName ? layerName : "";
	shape->sourceType = "line-layer";
	shape->sourceId = this->sourceId.getValue();
	shape->displayName = layerName ? SbString(layerName) :
			     this->overlayId.getValue();
	shape->geometryName = layerName ? SbString(layerName) :
			      SbString("line-layer");
	shape->sourceIdentity = shape->sourcePath.getValue();
	shape->cacheIdentity = shape->sourcePath.getValue();
	shape->databaseIntent = FALSE;
	shape->overlayIntent = TRUE;
	shape->hudIntent = FALSE;
	shape->localSource = TRUE;
	shape->sharedSource = FALSE;
	shape->nonDatabaseSource = TRUE;
	shape->drawMode = BOBOL_LOD_DRAW_DIAGNOSTIC;
	shape->recordRole = "overlay";
	shape->geometryKind = "line";
	shape->visible = TRUE;
	shape->selectable = this->selectable.getValue();
	shape->colorOverride = TRUE;
	shape->color = SbColor(
			   static_cast<float>(r) * COLOR_COMPONENT_SCALE,
			   static_cast<float>(g) * COLOR_COMPONENT_SCALE,
			   static_cast<float>(b) * COLOR_COMPONENT_SCALE);
	shape->setLineSet(points.data(), commands.data(),
			  static_cast<int>(points.size()));
	publishedChildren.push_back(shape);
	owners.push_back(std::move(shapeOwner));
	realized++;
    }

    auto replacement = this->getChildren()->prepareReplacement(publishedChildren);
    replacement->commit();
    replacement->notify();
    return realized;
}

SoBRLVListShape *
SoBRLLineLayerOverlay::getLayerShape(int index) const
{
    if (index < 0)
	return NULL;

    int seen = 0;
    for (int i = 0; i < this->getNumChildren(); i++) {
	SoNode *node = this->getChild(i);
	if (!node || !node->isOfType(SoBRLVListShape::getClassTypeId()))
	    continue;
	if (seen == index)
	    return static_cast<SoBRLVListShape *>(node);
	seen++;
    }
    return NULL;
}

int
SoBRLLineLayerOverlay::getLayerShapeCount(void) const
{
    int ret = 0;
    for (int i = 0; i < this->getNumChildren(); i++) {
	SoNode *node = this->getChild(i);
	if (node && node->isOfType(SoBRLVListShape::getClassTypeId()))
	    ret++;
    }
    return ret;
}

size_t
SoBRLLineLayerOverlay::getPointCount(void) const
{
    size_t ret = 0;
    for (int i = 0; i < this->getLayerShapeCount(); i++) {
	SoBRLVListShape *shape = this->getLayerShape(i);
	if (shape)
	    ret += static_cast<size_t>(shape->point.getNum());
    }
    return ret;
}
