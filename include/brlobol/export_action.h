/*                  E X P O R T _ A C T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/export_action.h */

#ifndef BRLOBOL_EXPORT_ACTION_H
#define BRLOBOL_EXPORT_ACTION_H

#include "brlobol/defines.h"
#include "brlobol/lod_realization.h"
#include "brlobol/source_mesh_request.h"

#include <Inventor/SbBox.h>
#include <Inventor/SbColor.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/actions/SoSubAction.h>

#include <stdint.h>
#include <vector>

class BRLObolLodService;
class SoBRLMeshShape;
struct db_i;

class BRLOBOL_EXPORT SoBRLExportAction : public SoAction {
    typedef SoAction inherited;

    SO_ACTION_HEADER(SoBRLExportAction);

public:
    enum GeometryPolicy {
	FULL_DETAIL = 0,
	DISPLAY_LEVEL = 1
    };

    struct LineRecord {
	SbString path;
	SbString sourceName;
	SbString sourceType;
	uint32_t sourceId;
	int regionId;
	int airCode;
	int materialId;
	int los;
	int materialColorValid;
	SbColor materialColor;
	SbString materialShader;
	int primitiveIndex;
	int selected;
	int highlighted;
	int ghosted;
	int hiddenLine;
	int editEmphasis;
	uint32_t lodPolicy;
	int colorOverride;
	SbColor color;
	SbVec3f a;
	SbVec3f b;
    };

    struct PointRecord {
	SbString path;
	SbString sourceName;
	SbString sourceType;
	uint32_t sourceId;
	int regionId;
	int airCode;
	int materialId;
	int los;
	int materialColorValid;
	SbColor materialColor;
	SbString materialShader;
	int primitiveIndex;
	int selected;
	int highlighted;
	int ghosted;
	int hiddenLine;
	int editEmphasis;
	uint32_t lodPolicy;
	int colorOverride;
	SbColor color;
	SbVec3f point;
    };

    struct TriangleRecord {
	SbString path;
	SbString sourceName;
	SbString sourceType;
	uint32_t sourceId;
	int regionId;
	int airCode;
	int materialId;
	int los;
	int materialColorValid;
	SbColor materialColor;
	SbString materialShader;
	int primitiveIndex;
	int vertexIndexA;
	int vertexIndexB;
	int vertexIndexC;
	int selected;
	int highlighted;
	int ghosted;
	int hiddenLine;
	int editEmphasis;
	uint32_t lodPolicy;
	int lodAvailable;
	int lodActiveLevel;
	uint32_t lodFaceCount;
	uint32_t lodPointCount;
	uint32_t lodOriginalPointCount;
	uint32_t lodNormalCount;
	int lodHasSnappedPoints;
	int lodHasNormals;
	SbVec3f lodBoundsMin;
	SbVec3f lodBoundsMax;
	int colorOverride;
	SbColor color;
	SbVec3f a;
	SbVec3f b;
	SbVec3f c;
    };

    SoBRLExportAction(void);
    virtual ~SoBRLExportAction(void);
    static void initClass(void);

    int getLineCount(void) const;
    const LineRecord &getLine(int index) const;
    int getPointCount(void) const;
    const PointRecord &getPoint(int index) const;
    int getTriangleCount(void) const;
    const TriangleRecord &getTriangle(int index) const;
    const SbBox3f &getBounds(void) const;
    void setGeometryPolicy(GeometryPolicy policy);
    GeometryPolicy getGeometryPolicy(void) const;
    unsigned int getSkippedLodDisplayMeshCount(void) const;
    int getSourceBackedFullDetailRequestCount(void) const;
    const BRLObolSourceMeshRequest &getSourceBackedFullDetailRequest(int index) const;
    SbBool makeSourceBackedFullDetailLodRequest(int index,
	    BRLObolLodRequest &request,
	    const BRLObolLodRequest *templateRequest = 0) const;
    SbBool appendSourceBackedFullDetailResult(
	    const BRLObolSourceMeshRequest &sourceRequest,
	    const BRLObolLodResult &result);
    int submitSourceBackedFullDetailRequests(BRLObolLodService *service,
	    uint64_t generation, struct db_i *dbip,
	    const BRLObolLodRequest *templateRequest = 0,
	    uint64_t maxFullDetailFaceCount = 0,
	    uint64_t maxFullDetailPointCount = 0) const;
    int consumeSourceBackedFullDetailResults(
	    const std::vector<BRLObolLodResult> &results,
	    const BRLObolLodRequest *templateRequest = 0);

protected:
    virtual void beginTraversal(SoNode *node);

private:
    static void nodeAction(SoAction *action, SoNode *node);
    static void vlistShapeAction(SoAction *action, SoNode *node);
    static void meshShapeAction(SoAction *action, SoNode *node);

    void resetResults(void);
    void appendSourceBackedFullDetailRequest(const SoBRLMeshShape *shape,
	    const SbMatrix &localToWorld);
    void appendLine(const SbString &path, const SbString &sourceName,
	    const SbString &sourceType, uint32_t sourceId,
	    int regionId, int airCode, int materialId, int los,
	    int materialColorValid, const SbColor &materialColor,
	    const SbString &materialShader, int primitiveIndex,
	    int selected, int highlighted, int ghosted,
	    int hiddenLine, int editEmphasis, uint32_t lodPolicy,
	    int colorOverride, const SbColor &color,
	    const SbVec3f &a, const SbVec3f &b);
    void appendPoint(const SbString &path, const SbString &sourceName,
	    const SbString &sourceType, uint32_t sourceId,
	    int regionId, int airCode, int materialId, int los,
	    int materialColorValid, const SbColor &materialColor,
	    const SbString &materialShader, int primitiveIndex,
	    int selected, int highlighted, int ghosted,
	    int hiddenLine, int editEmphasis, uint32_t lodPolicy,
	    int colorOverride, const SbColor &color,
	    const SbVec3f &point);
    void appendTriangle(const SbString &path, const SbString &sourceName,
	    const SbString &sourceType, uint32_t sourceId,
	    int regionId, int airCode, int materialId, int los,
	    int materialColorValid, const SbColor &materialColor,
	    const SbString &materialShader, int primitiveIndex,
	    int vertexIndexA, int vertexIndexB, int vertexIndexC,
	    int selected, int highlighted, int ghosted,
	    int hiddenLine, int editEmphasis, uint32_t lodPolicy,
	    int lodAvailable, int lodActiveLevel, uint32_t lodFaceCount,
	    uint32_t lodPointCount, uint32_t lodOriginalPointCount,
	    uint32_t lodNormalCount, int lodHasSnappedPoints,
	    int lodHasNormals, const SbVec3f &lodBoundsMin,
	    const SbVec3f &lodBoundsMax,
	    int colorOverride, const SbColor &color,
	    const SbVec3f &a, const SbVec3f &b, const SbVec3f &c);

    std::vector<LineRecord> lines;
    std::vector<PointRecord> points;
    std::vector<TriangleRecord> triangles;
    std::vector<BRLObolSourceMeshRequest> sourceBackedFullDetailRequests;
    SbBox3f bounds;
    GeometryPolicy geometryPolicy;
    unsigned int skippedLodDisplayMeshCount;
};

#endif /* BRLOBOL_EXPORT_ACTION_H */
