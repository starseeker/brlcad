/*             B M A T E R I A L O B J E C T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BMaterialObject.h */

#ifndef BOBOL_BMATERIALOBJECT_H
#define BOBOL_BMATERIALOBJECT_H

#include "BObol/BDefines.h"

#include <stdint.h>

#include <Inventor/fields/SoMFString.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFUInt32.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoSubNode.h>

class BOBOL_EXPORT SoBRLMaterialObject : public SoNode {
    typedef SoNode inherited;

    SO_NODE_HEADER(SoBRLMaterialObject);

public:
    SoSFString sourcePath;
    SoSFString sourceName;
    SoSFString sourceType;
    SoSFUInt32 sourceId;
    SoSFString materialName;
    SoSFString parentName;
    SoSFString materialSource;
    SoMFString propertyGroup;
    SoMFString propertyName;
    SoMFString propertyValue;

    SoBRLMaterialObject(void);
    static void initClass(void);

    int getPropertyCount(void) const;
    void clearProperties(void);
    void addProperty(const char *group, const char *name, const char *value);
    SbBool getProperty(int index, SbString &groupOut, SbString &nameOut,
	    SbString &valueOut) const;
    SbBool findProperty(const char *group, const char *name,
	    SbString &valueOut) const;

protected:
    virtual ~SoBRLMaterialObject(void);
};

#endif /* BOBOL_BMATERIALOBJECT_H */
