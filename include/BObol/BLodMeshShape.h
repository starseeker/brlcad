/*               B L O D M E S H S H A P E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BLodMeshShape.h */

#ifndef BOBOL_BLODMESHSHAPE_H
#define BOBOL_BLODMESHSHAPE_H

#include "BObol/BMeshShape.h"

class BOBOL_EXPORT SoBRLLodMeshShape : public SoBRLMeshShape {
    typedef SoBRLMeshShape inherited;

    SO_NODE_HEADER(SoBRLLodMeshShape);

public:
    SoBRLLodMeshShape(void);
    static void initClass(void);

protected:
    virtual ~SoBRLLodMeshShape(void);
};

#endif /* BOBOL_BLODMESHSHAPE_H */
