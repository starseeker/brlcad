/*                         A D C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/str.h"
#include "bv.h"

#include "BObol/BADC.h"
#include "BObol/BLodRealization.h"
#include "BObol/BVListShape.h"
#include "scalar_publication_private.h"

#include <Inventor/misc/SoChildList.h>
#include <Inventor/tools/SbModernUtils.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

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
    if (!this->visible.getValue()) {
	if (!this->getNumChildren())
	    return NULL;
	auto replacement = this->getChildren()->prepareReplacement({});
	replacement->commit();
	replacement->notify();
	return NULL;
    }

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
    SbModernUtils::SoNodeRef lineOwner(lineShape);
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
    lineShape->drawMode = BOBOL_LOD_DRAW_DIAGNOSTIC;
    lineShape->recordRole = "overlay";
    lineShape->geometryKind = "line";
    lineShape->sourceId = static_cast<uint32_t>(d);
    lineShape->color = this->lineColor.getValue();
	lineShape->lineWidth = effectiveLineWidth;
    lineShape->setLineSet(linePoints, lineCommands, 6);

    SoBRLVListShape *tickShape = new SoBRLVListShape;
    SbModernUtils::SoNodeRef tickOwner(tickShape);
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
    tickShape->drawMode = BOBOL_LOD_DRAW_DIAGNOSTIC;
    tickShape->recordRole = "overlay";
    tickShape->geometryKind = "line";
    tickShape->sourceId = static_cast<uint32_t>(d);
    tickShape->color = this->tickColor.getValue();
    tickShape->lineWidth = effectiveLineWidth;
    tickShape->setLineSet(tickPoints, tickCommands, 2);
    auto replacement = this->getChildren()->prepareReplacement(
	std::vector<SoNode *>{lineShape, tickShape});
    replacement->commit();
    replacement->notify();
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

int
bobol_adc_configure_from_view(SoBRLADC *adc,
	const struct bv_adc_state *state)
{
    if (!adc || !state)
	return 0;

    constexpr float DEFAULT_ADC_DISTANCE = 1.0f;
    constexpr int DEFAULT_LINE_WIDTH = 1;
    constexpr float COLOR_CHANNEL_MAX = 255.0f;
    SbModernUtils::SoNodeRef candidateOwner(new SoBRLADC);
    auto *candidate = static_cast<SoBRLADC *>(candidateOwner.get());
    copy_publication_scalar_fields(*candidate, *adc);
    candidate->center = SbVec3f(
	static_cast<float>(state->pos_model[X]),
	static_cast<float>(state->pos_model[Y]),
	static_cast<float>(state->pos_model[Z]));
    candidate->angleDegrees = static_cast<float>(state->a1);
    candidate->distance = static_cast<float>(
	state->dst > SMALL_FASTF ? state->dst : DEFAULT_ADC_DISTANCE);
    candidate->lineColor = SbColor(
	static_cast<float>(state->line_color[0]) / COLOR_CHANNEL_MAX,
	static_cast<float>(state->line_color[1]) / COLOR_CHANNEL_MAX,
	static_cast<float>(state->line_color[2]) / COLOR_CHANNEL_MAX);
    candidate->tickColor = SbColor(
	static_cast<float>(state->tick_color[0]) / COLOR_CHANNEL_MAX,
	static_cast<float>(state->tick_color[1]) / COLOR_CHANNEL_MAX,
	static_cast<float>(state->tick_color[2]) / COLOR_CHANNEL_MAX);
    candidate->lineWidth = state->line_width > 0 ?
	state->line_width : DEFAULT_LINE_WIDTH;
    candidate->visible = state->draw ? TRUE : FALSE;
    candidate->rebuildGeometry();

    publish_scalar_fields_and_children(*adc, *candidate);
    return state->draw ? 1 : 0;
}
