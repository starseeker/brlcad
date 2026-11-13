/*              L O D _ M E S H _ S H A P E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BLodMeshShape.h"

SO_NODE_SOURCE(SoBRLLodMeshShape);

SoBRLLodMeshShape::SoBRLLodMeshShape(void)
{
    SO_NODE_CONSTRUCTOR(SoBRLLodMeshShape);
    this->setLodBackedMesh(TRUE);
    this->setLodPreserveFullDetail(FALSE);
}

SoBRLLodMeshShape::~SoBRLLodMeshShape(void)
{
}

void
SoBRLLodMeshShape::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLLodMeshShape, SoBRLMeshShape, "SoBRLMeshShape");
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
