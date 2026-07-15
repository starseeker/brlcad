/*                 L O D _ M E S H _ S H A P E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/lod_mesh_shape.h */

#ifndef BRLOBOL_LOD_MESH_SHAPE_H
#define BRLOBOL_LOD_MESH_SHAPE_H

#include "brlobol/mesh_shape.h"

class BRLOBOL_EXPORT SoBRLLodMeshShape : public SoBRLMeshShape {
    typedef SoBRLMeshShape inherited;

    SO_NODE_HEADER(SoBRLLodMeshShape);

public:
    SoBRLLodMeshShape(void);
    static void initClass(void);

protected:
    virtual ~SoBRLLodMeshShape(void);
};

#endif /* BRLOBOL_LOD_MESH_SHAPE_H */
