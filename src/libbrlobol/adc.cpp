/*                         A D C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/adc.h"
#include "brlobol/lod_realization.h"
#include "brlobol/vlist_shape.h"

#include <cmath>

SO_NODE_SOURCE(SoBRLADC);

SoBRLADC::SoBRLADC(void)
{
    SO_NODE_CONSTRUCTOR(SoBRLADC);

    SO_NODE_ADD_FIELD(overlayId, ("overlay::adc"));
    SO_NODE_ADD_FIELD(center, (0.0f, 0.0f, 0.0f));
    SO_NODE_ADD_FIELD(angleDegrees, (0.0f));
    SO_NODE_ADD_FIELD(distance, (1.0f));
    SO_NODE_ADD_FIELD(crosshairSize, (1.0f));
    SO_NODE_ADD_FIELD(tickSize, (0.5f));
    SO_NODE_ADD_FIELD(visible, (TRUE));
}

SoBRLADC::~SoBRLADC(void)
{
}

void
SoBRLADC::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLADC, SoSeparator, "Separator");
}

SoBRLVListShape *
SoBRLADC::rebuildGeometry(void)
{
    this->removeAllChildren();
    if (!this->visible.getValue())
	return NULL;

    float d = this->distance.getValue();
    if (d < 0.0f)
	d = 0.0f;

    float cross = this->crosshairSize.getValue();
    if (cross <= 0.0f)
	cross = 1.0f;

    float tick = this->tickSize.getValue();
    if (tick <= 0.0f)
	tick = cross * 0.5f;

    const float radians = this->angleDegrees.getValue() *
	static_cast<float>(M_PI / 180.0);
    SbVec3f c = this->center.getValue();
    SbVec3f dir(std::cos(radians), std::sin(radians), 0.0f);
    SbVec3f perp(-dir[1], dir[0], 0.0f);
    SbVec3f end = c + dir * d;

    SbVec3f points[8] = {
	SbVec3f(c[0] - cross, c[1], c[2]),
	SbVec3f(c[0] + cross, c[1], c[2]),
	SbVec3f(c[0], c[1] - cross, c[2]),
	SbVec3f(c[0], c[1] + cross, c[2]),
	c,
	end,
	end - perp * tick,
	end + perp * tick
    };
    int32_t commands[8] = {
	SoBRLVListShape::MOVE, SoBRLVListShape::DRAW,
	SoBRLVListShape::MOVE, SoBRLVListShape::DRAW,
	SoBRLVListShape::MOVE, SoBRLVListShape::DRAW,
	SoBRLVListShape::MOVE, SoBRLVListShape::DRAW
    };

    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->sourcePath = this->overlayId.getValue();
    shape->displayName = this->overlayId.getValue();
    shape->geometryName = "adc";
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
    shape->sourceId = static_cast<uint32_t>(d);
    shape->setLineSet(points, commands, 8);
    this->addChild(shape);
    return shape;
}

SoBRLVListShape *
SoBRLADC::getGeometryShape(void) const
{
    for (int i = 0; i < this->getNumChildren(); i++) {
	SoNode *node = this->getChild(i);
	if (node && node->isOfType(SoBRLVListShape::getClassTypeId()))
	    return static_cast<SoBRLVListShape *>(node);
    }
    return NULL;
}
