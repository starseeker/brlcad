/*             H U D _ L A B E L _ O V E R L A Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BHUDLabelOverlay.h"

#include <Inventor/annex/HUD/nodekits/SoHUDKit.h>
#include <Inventor/annex/HUD/nodes/SoHUDLabel.h>

SO_NODE_SOURCE(SoBRLHUDLabelOverlay);

SoBRLHUDLabelOverlay::SoBRLHUDLabelOverlay(void) :
    hudLabel(NULL)
{
    SO_NODE_CONSTRUCTOR(SoBRLHUDLabelOverlay);

    SO_NODE_ADD_FIELD(labelId, ("hud::label"));
    SO_NODE_ADD_FIELD(sourceId, (0));
    SO_NODE_ADD_FIELD(text, (""));
    SO_NODE_ADD_FIELD(position, (0.0f, 0.0f));
    SO_NODE_ADD_FIELD(color, (1.0f, 1.0f, 1.0f));
    SO_NODE_ADD_FIELD(fontSize, (12.0f));
    SO_NODE_ADD_FIELD(visible, (TRUE));
}

SoBRLHUDLabelOverlay::~SoBRLHUDLabelOverlay(void)
{
}

void
SoBRLHUDLabelOverlay::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLHUDLabelOverlay, SoSeparator, "Separator");
}

SoHUDLabel *
SoBRLHUDLabelOverlay::rebuildGeometry(void)
{
    this->hudLabel = NULL;
    this->removeAllChildren();
    if (!this->visible.getValue() || this->text.getValue().getLength() == 0)
	return NULL;

    SoHUDKit *hud = new SoHUDKit;
    SoHUDLabel *label = new SoHUDLabel;
    label->position = this->position.getValue();
    label->string.set1Value(0, this->text.getValue());
    label->color = this->color.getValue();

    float size = this->fontSize.getValue();
    if (size <= 0.0f)
	size = 12.0f;
    label->fontSize = size;

    hud->addWidget(label);
    this->addChild(hud);
    this->hudLabel = label;
    return label;
}

SoHUDKit *
SoBRLHUDLabelOverlay::getHUDKit(void) const
{
    for (int i = 0; i < this->getNumChildren(); i++) {
	SoNode *node = this->getChild(i);
	if (node && node->isOfType(SoHUDKit::getClassTypeId()))
	    return static_cast<SoHUDKit *>(node);
    }
    return NULL;
}

SoHUDLabel *
SoBRLHUDLabelOverlay::getHUDLabel(void) const
{
    if (!this->getHUDKit())
	return NULL;
    return this->hudLabel;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
