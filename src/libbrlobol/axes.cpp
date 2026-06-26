/*                         A X E S . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/axes.h"
#include "brlobol/lod_realization.h"
#include "brlobol/vlist_shape.h"

SO_NODE_SOURCE(SoBRLAxes);

SoBRLAxes::SoBRLAxes(void)
{
    SO_NODE_CONSTRUCTOR(SoBRLAxes);

    SO_NODE_ADD_FIELD(overlayId, ("overlay::axes"));
    SO_NODE_ADD_FIELD(origin, (0.0f, 0.0f, 0.0f));
    SO_NODE_ADD_FIELD(size, (1.0f));
    SO_NODE_ADD_FIELD(visible, (TRUE));
}

SoBRLAxes::~SoBRLAxes(void)
{
}

void
SoBRLAxes::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLAxes, SoSeparator, "Separator");
}

SoBRLVListShape *
SoBRLAxes::rebuildGeometry(void)
{
    this->removeAllChildren();
    if (!this->visible.getValue())
	return NULL;

    float s = this->size.getValue();
    if (s <= 0.0f)
	s = 1.0f;

    SbVec3f o = this->origin.getValue();
    SbVec3f points[6] = {
	o, SbVec3f(o[0] + s, o[1], o[2]),
	o, SbVec3f(o[0], o[1] + s, o[2]),
	o, SbVec3f(o[0], o[1], o[2] + s)
    };
    int32_t commands[6] = {
	SoBRLVListShape::MOVE, SoBRLVListShape::DRAW,
	SoBRLVListShape::MOVE, SoBRLVListShape::DRAW,
	SoBRLVListShape::MOVE, SoBRLVListShape::DRAW
    };

    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->sourcePath = this->overlayId.getValue();
    shape->displayName = this->overlayId.getValue();
    shape->geometryName = "axes";
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
    shape->sourceId = static_cast<uint32_t>(s);
    shape->setLineSet(points, commands, 6);
    this->addChild(shape);
    return shape;
}

SoBRLVListShape *
SoBRLAxes::getGeometryShape(void) const
{
    for (int i = 0; i < this->getNumChildren(); i++) {
	SoNode *node = this->getChild(i);
	if (node && node->isOfType(SoBRLVListShape::getClassTypeId()))
	    return static_cast<SoBRLVListShape *>(node);
    }
    return NULL;
}
