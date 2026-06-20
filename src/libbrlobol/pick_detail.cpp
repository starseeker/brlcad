/*                   P I C K _ D E T A I L . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/pick_detail.h"

SO_DETAIL_SOURCE(SoBRLPickDetail);

SoBRLPickDetail::SoBRLPickDetail(void) :
    dbpath(""),
    sourceName(""),
    sourceType(""),
    materialShader(""),
    modelPoint(0.0f, 0.0f, 0.0f),
    materialColor(1.0f, 1.0f, 1.0f),
    sourceId(0),
    regionId(0),
    airCode(0),
    materialId(0),
    los(0),
    materialColorValid(FALSE),
    primitiveKind(UNKNOWN),
    primitiveIndex(-1)
{
    this->faceVertexIndex[0] = -1;
    this->faceVertexIndex[1] = -1;
    this->faceVertexIndex[2] = -1;
    this->nearestFaceEdgeSlot = -1;
    this->nearestFaceEdgeVertexIndex[0] = -1;
    this->nearestFaceEdgeVertexIndex[1] = -1;
    this->nearestFaceVertexSlot = -1;
    this->nearestFaceVertexIndex = -1;
}

SoBRLPickDetail::SoBRLPickDetail(const SoBRLPickDetail &other) :
    SoDetail(other),
    dbpath(other.dbpath),
    sourceName(other.sourceName),
    sourceType(other.sourceType),
    materialShader(other.materialShader),
    modelPoint(other.modelPoint),
    materialColor(other.materialColor),
    sourceId(other.sourceId),
    regionId(other.regionId),
    airCode(other.airCode),
    materialId(other.materialId),
    los(other.los),
    materialColorValid(other.materialColorValid),
    primitiveKind(other.primitiveKind),
    primitiveIndex(other.primitiveIndex)
{
    this->faceVertexIndex[0] = other.faceVertexIndex[0];
    this->faceVertexIndex[1] = other.faceVertexIndex[1];
    this->faceVertexIndex[2] = other.faceVertexIndex[2];
    this->nearestFaceEdgeSlot = other.nearestFaceEdgeSlot;
    this->nearestFaceEdgeVertexIndex[0] = other.nearestFaceEdgeVertexIndex[0];
    this->nearestFaceEdgeVertexIndex[1] = other.nearestFaceEdgeVertexIndex[1];
    this->nearestFaceVertexSlot = other.nearestFaceVertexSlot;
    this->nearestFaceVertexIndex = other.nearestFaceVertexIndex;
}

SoBRLPickDetail::~SoBRLPickDetail(void)
{
}

void
SoBRLPickDetail::initClass(void)
{
    SO_DETAIL_INIT_CLASS(SoBRLPickDetail, SoDetail);
}

SoDetail *
SoBRLPickDetail::copy(void) const
{
    return new SoBRLPickDetail(*this);
}

void
SoBRLPickDetail::setPath(const SbString &path)
{
    this->dbpath = path;
}

const SbString &
SoBRLPickDetail::getPath(void) const
{
    return this->dbpath;
}

void
SoBRLPickDetail::setSourceName(const SbString &name)
{
    this->sourceName = name;
}

const SbString &
SoBRLPickDetail::getSourceName(void) const
{
    return this->sourceName;
}

void
SoBRLPickDetail::setSourceType(const SbString &type)
{
    this->sourceType = type;
}

const SbString &
SoBRLPickDetail::getSourceType(void) const
{
    return this->sourceType;
}

void
SoBRLPickDetail::setSourceId(uint32_t id)
{
    this->sourceId = id;
}

uint32_t
SoBRLPickDetail::getSourceId(void) const
{
    return this->sourceId;
}

void
SoBRLPickDetail::setRegionId(int id)
{
    this->regionId = id;
}

int
SoBRLPickDetail::getRegionId(void) const
{
    return this->regionId;
}

void
SoBRLPickDetail::setAirCode(int air)
{
    this->airCode = air;
}

int
SoBRLPickDetail::getAirCode(void) const
{
    return this->airCode;
}

void
SoBRLPickDetail::setMaterialId(int material)
{
    this->materialId = material;
}

int
SoBRLPickDetail::getMaterialId(void) const
{
    return this->materialId;
}

void
SoBRLPickDetail::setLos(int value)
{
    this->los = value;
}

int
SoBRLPickDetail::getLos(void) const
{
    return this->los;
}

void
SoBRLPickDetail::setMaterialColor(SbBool valid, const SbColor &color)
{
    this->materialColorValid = valid;
    this->materialColor = color;
}

SbBool
SoBRLPickDetail::hasMaterialColor(void) const
{
    return this->materialColorValid;
}

const SbColor &
SoBRLPickDetail::getMaterialColor(void) const
{
    return this->materialColor;
}

void
SoBRLPickDetail::setMaterialShader(const SbString &shader)
{
    this->materialShader = shader;
}

const SbString &
SoBRLPickDetail::getMaterialShader(void) const
{
    return this->materialShader;
}

void
SoBRLPickDetail::setPrimitive(PrimitiveKind kind, int index)
{
    this->primitiveKind = kind;
    this->primitiveIndex = index;
}

SoBRLPickDetail::PrimitiveKind
SoBRLPickDetail::getPrimitiveKind(void) const
{
    return this->primitiveKind;
}

int
SoBRLPickDetail::getPrimitiveIndex(void) const
{
    return this->primitiveIndex;
}

void
SoBRLPickDetail::setFaceVertexIndices(int indexA, int indexB, int indexC)
{
    this->faceVertexIndex[0] = indexA;
    this->faceVertexIndex[1] = indexB;
    this->faceVertexIndex[2] = indexC;
}

int
SoBRLPickDetail::getFaceVertexIndex(int vertexSlot) const
{
    if (vertexSlot < 0 || vertexSlot >= 3)
	return -1;
    return this->faceVertexIndex[vertexSlot];
}

int
SoBRLPickDetail::getFaceVertexIndexA(void) const
{
    return this->getFaceVertexIndex(0);
}

int
SoBRLPickDetail::getFaceVertexIndexB(void) const
{
    return this->getFaceVertexIndex(1);
}

int
SoBRLPickDetail::getFaceVertexIndexC(void) const
{
    return this->getFaceVertexIndex(2);
}

void
SoBRLPickDetail::setNearestFaceEdge(int edgeSlot, int indexA, int indexB)
{
    this->nearestFaceEdgeSlot = edgeSlot;
    this->nearestFaceEdgeVertexIndex[0] = indexA;
    this->nearestFaceEdgeVertexIndex[1] = indexB;
}

int
SoBRLPickDetail::getNearestFaceEdgeSlot(void) const
{
    return this->nearestFaceEdgeSlot;
}

int
SoBRLPickDetail::getNearestFaceEdgeVertexIndexA(void) const
{
    return this->nearestFaceEdgeVertexIndex[0];
}

int
SoBRLPickDetail::getNearestFaceEdgeVertexIndexB(void) const
{
    return this->nearestFaceEdgeVertexIndex[1];
}

void
SoBRLPickDetail::setNearestFaceVertex(int vertexSlot, int vertexIndex)
{
    this->nearestFaceVertexSlot = vertexSlot;
    this->nearestFaceVertexIndex = vertexIndex;
}

int
SoBRLPickDetail::getNearestFaceVertexSlot(void) const
{
    return this->nearestFaceVertexSlot;
}

int
SoBRLPickDetail::getNearestFaceVertexIndex(void) const
{
    return this->nearestFaceVertexIndex;
}

void
SoBRLPickDetail::setModelPoint(const SbVec3f &point)
{
    this->modelPoint = point;
}

const SbVec3f &
SoBRLPickDetail::getModelPoint(void) const
{
    return this->modelPoint;
}
