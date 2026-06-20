/*                   P I C K _ D E T A I L . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/pick_detail.h */

#ifndef BRLOBOL_PICK_DETAIL_H
#define BRLOBOL_PICK_DETAIL_H

#include "brlobol/defines.h"

#include <Inventor/SbColor.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/details/SoDetail.h>
#include <Inventor/details/SoSubDetail.h>

class BRLOBOL_EXPORT SoBRLPickDetail : public SoDetail {
    typedef SoDetail inherited;

    SO_DETAIL_HEADER(SoBRLPickDetail);

public:
    enum PrimitiveKind {
	UNKNOWN = 0,
	LINE_SEGMENT = 1,
	POINT = 2,
	FACE = 3
    };

    SoBRLPickDetail(void);
    SoBRLPickDetail(const SoBRLPickDetail &other);
    virtual ~SoBRLPickDetail(void);

    static void initClass(void);
    virtual SoDetail *copy(void) const;

    void setPath(const SbString &path);
    const SbString &getPath(void) const;

    void setSourceName(const SbString &name);
    const SbString &getSourceName(void) const;

    void setSourceType(const SbString &type);
    const SbString &getSourceType(void) const;

    void setSourceId(uint32_t id);
    uint32_t getSourceId(void) const;

    void setRegionId(int id);
    int getRegionId(void) const;

    void setAirCode(int air);
    int getAirCode(void) const;

    void setMaterialId(int material);
    int getMaterialId(void) const;

    void setLos(int los);
    int getLos(void) const;

    void setMaterialColor(SbBool valid, const SbColor &color);
    SbBool hasMaterialColor(void) const;
    const SbColor &getMaterialColor(void) const;

    void setMaterialShader(const SbString &shader);
    const SbString &getMaterialShader(void) const;

    void setPrimitive(PrimitiveKind kind, int index);
    PrimitiveKind getPrimitiveKind(void) const;
    int getPrimitiveIndex(void) const;

    void setFaceVertexIndices(int indexA, int indexB, int indexC);
    int getFaceVertexIndex(int vertexSlot) const;
    int getFaceVertexIndexA(void) const;
    int getFaceVertexIndexB(void) const;
    int getFaceVertexIndexC(void) const;

    void setNearestFaceEdge(int edgeSlot, int indexA, int indexB);
    int getNearestFaceEdgeSlot(void) const;
    int getNearestFaceEdgeVertexIndexA(void) const;
    int getNearestFaceEdgeVertexIndexB(void) const;

    void setNearestFaceVertex(int vertexSlot, int vertexIndex);
    int getNearestFaceVertexSlot(void) const;
    int getNearestFaceVertexIndex(void) const;

    void setModelPoint(const SbVec3f &point);
    const SbVec3f &getModelPoint(void) const;

private:
    SbString dbpath;
    SbString sourceName;
    SbString sourceType;
    SbString materialShader;
    SbVec3f modelPoint;
    SbColor materialColor;
    uint32_t sourceId;
    int regionId;
    int airCode;
    int materialId;
    int los;
    SbBool materialColorValid;
    PrimitiveKind primitiveKind;
    int primitiveIndex;
    int faceVertexIndex[3];
    int nearestFaceEdgeSlot;
    int nearestFaceEdgeVertexIndex[2];
    int nearestFaceVertexSlot;
    int nearestFaceVertexIndex;
};

#endif /* BRLOBOL_PICK_DETAIL_H */
