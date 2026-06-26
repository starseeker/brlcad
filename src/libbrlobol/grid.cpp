/*                         G R I D . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/grid.h"
#include "brlobol/lod_realization.h"
#include "brlobol/vlist_shape.h"

#include <algorithm>
#include <vector>

SO_NODE_SOURCE(SoBRLGrid);

SoBRLGrid::SoBRLGrid(void)
{
    SO_NODE_CONSTRUCTOR(SoBRLGrid);

    SO_NODE_ADD_FIELD(overlayId, ("overlay::grid"));
    SO_NODE_ADD_FIELD(center, (0.0f, 0.0f, 0.0f));
    SO_NODE_ADD_FIELD(spacing, (1.0f));
    SO_NODE_ADD_FIELD(divisions, (4));
    SO_NODE_ADD_FIELD(visible, (TRUE));
}

SoBRLGrid::~SoBRLGrid(void)
{
}

void
SoBRLGrid::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLGrid, SoSeparator, "Separator");
}

SoBRLVListShape *
SoBRLGrid::rebuildGeometry(void)
{
    this->removeAllChildren();
    if (!this->visible.getValue())
	return NULL;

    int divs = std::max(1, this->divisions.getValue());
    float step = this->spacing.getValue();
    if (step <= 0.0f)
	step = 1.0f;

    SbVec3f c = this->center.getValue();
    float extent = step * static_cast<float>(divs);
    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    points.reserve(static_cast<size_t>((divs * 2 + 1) * 4));
    commands.reserve(static_cast<size_t>((divs * 2 + 1) * 4));

    for (int i = -divs; i <= divs; i++) {
	float v = static_cast<float>(i) * step;

	points.push_back(SbVec3f(c[0] - extent, c[1] + v, c[2]));
	commands.push_back(SoBRLVListShape::MOVE);
	points.push_back(SbVec3f(c[0] + extent, c[1] + v, c[2]));
	commands.push_back(SoBRLVListShape::DRAW);

	points.push_back(SbVec3f(c[0] + v, c[1] - extent, c[2]));
	commands.push_back(SoBRLVListShape::MOVE);
	points.push_back(SbVec3f(c[0] + v, c[1] + extent, c[2]));
	commands.push_back(SoBRLVListShape::DRAW);
    }

    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->sourcePath = this->overlayId.getValue();
    shape->displayName = this->overlayId.getValue();
    shape->geometryName = "grid";
    shape->sourceIdentity = this->overlayId.getValue();
    shape->cacheIdentity = this->overlayId.getValue();
    shape->databaseIntent = FALSE;
    shape->overlayIntent = TRUE;
    shape->hudIntent = FALSE;
    shape->localSource = TRUE;
    shape->sharedSource = FALSE;
    shape->nonDatabaseSource = TRUE;
    shape->drawMode = BRLOBOL_LOD_DRAW_DIAGNOSTIC;
    shape->recordRole = "overlay";
    shape->geometryKind = "line";
    shape->sourceId = static_cast<uint32_t>(divs);
    shape->setLineSet(points.data(), commands.data(), static_cast<int>(points.size()));
    this->addChild(shape);
    return shape;
}

SoBRLVListShape *
SoBRLGrid::getGeometryShape(void) const
{
    for (int i = 0; i < this->getNumChildren(); i++) {
	SoNode *node = this->getChild(i);
	if (node && node->isOfType(SoBRLVListShape::getClassTypeId()))
	    return static_cast<SoBRLVListShape *>(node);
    }
    return NULL;
}
