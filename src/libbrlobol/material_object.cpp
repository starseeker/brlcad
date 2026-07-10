/*             M A T E R I A L _ O B J E C T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/material_object.h"

#include <string.h>

SO_NODE_SOURCE(SoBRLMaterialObject);

SoBRLMaterialObject::SoBRLMaterialObject(void)
{
    SO_NODE_CONSTRUCTOR(SoBRLMaterialObject);

    SO_NODE_ADD_FIELD(sourcePath, (""));
    SO_NODE_ADD_FIELD(sourceName, (""));
    SO_NODE_ADD_FIELD(sourceType, ("material"));
    SO_NODE_ADD_FIELD(sourceId, (0));
    SO_NODE_ADD_FIELD(materialName, (""));
    SO_NODE_ADD_FIELD(parentName, (""));
    SO_NODE_ADD_FIELD(materialSource, (""));
    SO_NODE_ADD_FIELD(propertyGroup, (""));
    SO_NODE_ADD_FIELD(propertyName, (""));
    SO_NODE_ADD_FIELD(propertyValue, (""));
    this->clearProperties();
}

SoBRLMaterialObject::~SoBRLMaterialObject(void)
{
}

void
SoBRLMaterialObject::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLMaterialObject, SoNode, "Node");
}

int
SoBRLMaterialObject::getPropertyCount(void) const
{
    return this->propertyName.getNum();
}

void
SoBRLMaterialObject::clearProperties(void)
{
    this->propertyGroup.setNum(0);
    this->propertyName.setNum(0);
    this->propertyValue.setNum(0);
}

void
SoBRLMaterialObject::addProperty(const char *group, const char *name,
				 const char *value)
{
    const int index = this->getPropertyCount();
    this->propertyGroup.set1Value(index, group ? group : "");
    this->propertyName.set1Value(index, name ? name : "");
    this->propertyValue.set1Value(index, value ? value : "");
}

SbBool
SoBRLMaterialObject::getProperty(int index, SbString &groupOut,
				 SbString &nameOut, SbString &valueOut) const
{
    if (index < 0 ||
	index >= this->propertyGroup.getNum() ||
	index >= this->propertyName.getNum() ||
	index >= this->propertyValue.getNum())
	return FALSE;

    groupOut = this->propertyGroup[index];
    nameOut = this->propertyName[index];
    valueOut = this->propertyValue[index];
    return TRUE;
}

SbBool
SoBRLMaterialObject::findProperty(const char *group, const char *name,
				  SbString &valueOut) const
{
    if (!group || !name)
	return FALSE;

    for (int i = 0; i < this->getPropertyCount(); i++) {
	SbString propertyGroupValue;
	SbString propertyNameValue;
	SbString propertyValueValue;
	if (!this->getProperty(i, propertyGroupValue, propertyNameValue,
			       propertyValueValue))
	    continue;
	if (strcmp(propertyGroupValue.getString(), group) == 0 &&
	    strcmp(propertyNameValue.getString(), name) == 0) {
	    valueOut = propertyValueValue;
	    return TRUE;
	}
    }

    return FALSE;
}
