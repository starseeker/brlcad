/*                         A X E S . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bv.h"

#include "BObol/BAxes.h"
#include "BObol/BLodRealization.h"
#include "BObol/BVListShape.h"
#include "scalar_publication_private.h"

#include <Inventor/misc/SoChildList.h>
#include <Inventor/tools/SbModernUtils.h>

#include <vector>

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
    if (!this->visible.getValue()) {
	if (!this->getNumChildren())
	    return NULL;
	auto replacement = this->getChildren()->prepareReplacement({});
	replacement->commit();
	replacement->notify();
	return NULL;
    }

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
    SbModernUtils::SoNodeRef owner(shape);
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
    shape->drawMode = BOBOL_LOD_DRAW_DIAGNOSTIC;
    shape->recordRole = "overlay";
    shape->geometryKind = "line";
    shape->sourceId = static_cast<uint32_t>(s);
    shape->setLineSet(points, commands, 6);
    auto replacement = this->getChildren()->prepareReplacement({shape});
    replacement->commit();
    replacement->notify();
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

int
bobol_axes_configure_from_view(SoBRLAxes *axes,
	const struct bv_axes_state *state)
{
    if (!axes || !state)
	return 0;

    constexpr float DEFAULT_AXES_SIZE = 1.0f;
    SbModernUtils::SoNodeRef candidateOwner(new SoBRLAxes);
    auto *candidate = static_cast<SoBRLAxes *>(candidateOwner.get());
    copy_publication_scalar_fields(*candidate, *axes);
    candidate->origin = SbVec3f(
	static_cast<float>(state->axes_pos[X]),
	static_cast<float>(state->axes_pos[Y]),
	static_cast<float>(state->axes_pos[Z]));
    candidate->size = static_cast<float>(
	state->axes_size > SMALL_FASTF ? state->axes_size : DEFAULT_AXES_SIZE);
    candidate->visible = state->draw ? TRUE : FALSE;
    candidate->rebuildGeometry();

    publish_scalar_fields_and_children(*axes, *candidate);
    return state->draw ? 1 : 0;
}
