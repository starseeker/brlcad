/*                         A D C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/str.h"

#include "brlobol/adc.h"
#include "brlobol/lod_realization.h"
#include "brlobol/vlist_shape.h"

#include <algorithm>
#include <cmath>
#include <cstring>

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
    SO_NODE_ADD_FIELD(lineColor, (1.0f, 1.0f, 1.0f));
    SO_NODE_ADD_FIELD(tickColor, (1.0f, 1.0f, 1.0f));
    SO_NODE_ADD_FIELD(lineWidth, (1));
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

    SbVec3f linePoints[6] = {
	SbVec3f(c[0] - cross, c[1], c[2]),
	SbVec3f(c[0] + cross, c[1], c[2]),
	SbVec3f(c[0], c[1] - cross, c[2]),
	SbVec3f(c[0], c[1] + cross, c[2]),
	c,
	end
    };
    int32_t lineCommands[6] = {
	SoBRLVListShape::MOVE, SoBRLVListShape::DRAW,
	SoBRLVListShape::MOVE, SoBRLVListShape::DRAW,
	SoBRLVListShape::MOVE, SoBRLVListShape::DRAW
    };
    SbVec3f tickPoints[2] = {
	end - perp * tick,
	end + perp * tick
    };
	int32_t tickCommands[2] = {
	SoBRLVListShape::MOVE, SoBRLVListShape::DRAW
    };

    const int effectiveLineWidth = std::max(1, this->lineWidth.getValue());
    SoBRLVListShape *lineShape = new SoBRLVListShape;
    lineShape->sourcePath = this->overlayId.getValue();
    lineShape->displayName = this->overlayId.getValue();
    lineShape->geometryName = "adc-line";
    lineShape->sourceIdentity = this->overlayId.getValue();
    lineShape->cacheIdentity = this->overlayId.getValue();
    lineShape->databaseIntent = FALSE;
    lineShape->overlayIntent = TRUE;
    lineShape->hudIntent = FALSE;
    lineShape->localSource = TRUE;
    lineShape->sharedSource = FALSE;
    lineShape->nonDatabaseSource = TRUE;
    lineShape->drawMode = BRLOBOL_LOD_DRAW_DIAGNOSTIC;
    lineShape->recordRole = "overlay";
    lineShape->geometryKind = "line";
    lineShape->sourceId = static_cast<uint32_t>(d);
    lineShape->color = this->lineColor.getValue();
    lineShape->lineWidth = effectiveLineWidth;
    lineShape->setLineSet(linePoints, lineCommands, 6);
    this->addChild(lineShape);

    SoBRLVListShape *tickShape = new SoBRLVListShape;
    tickShape->sourcePath = this->overlayId.getValue();
    tickShape->displayName = this->overlayId.getValue();
    tickShape->geometryName = "adc-tick";
    tickShape->sourceIdentity = this->overlayId.getValue();
    tickShape->cacheIdentity = this->overlayId.getValue();
    tickShape->databaseIntent = FALSE;
    tickShape->overlayIntent = TRUE;
    tickShape->hudIntent = FALSE;
    tickShape->localSource = TRUE;
    tickShape->sharedSource = FALSE;
    tickShape->nonDatabaseSource = TRUE;
    tickShape->drawMode = BRLOBOL_LOD_DRAW_DIAGNOSTIC;
    tickShape->recordRole = "overlay";
    tickShape->geometryKind = "line";
    tickShape->sourceId = static_cast<uint32_t>(d);
    tickShape->color = this->tickColor.getValue();
    tickShape->lineWidth = effectiveLineWidth;
    tickShape->setLineSet(tickPoints, tickCommands, 2);
    this->addChild(tickShape);
    return lineShape;
}

SoBRLVListShape *
SoBRLADC::getGeometryShape(void) const
{
    for (int i = 0; i < this->getNumChildren(); i++) {
	SoNode *node = this->getChild(i);
        if (node && node->isOfType(SoBRLVListShape::getClassTypeId()) &&
	    bu_strcmp(static_cast<SoBRLVListShape *>(node)->geometryName.getValue().getString(),
		"adc-line") == 0)
	    return static_cast<SoBRLVListShape *>(node);
    }
    return NULL;
}

SoBRLVListShape *
SoBRLADC::getTickGeometryShape(void) const
{
    for (int i = 0; i < this->getNumChildren(); i++) {
	SoNode *node = this->getChild(i);
	if (node && node->isOfType(SoBRLVListShape::getClassTypeId()) &&
	    bu_strcmp(static_cast<SoBRLVListShape *>(node)->geometryName.getValue().getString(),
		"adc-tick") == 0)
	    return static_cast<SoBRLVListShape *>(node);
    }
    return NULL;
}
