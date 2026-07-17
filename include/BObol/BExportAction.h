/*               B E X P O R T A C T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BExportAction.h */

#ifndef BOBOL_BEXPORTACTION_H
#define BOBOL_BEXPORTACTION_H

#include "BObol/BDefines.h"
#include "BObol/BLodRealization.h"
#include "BObol/BSourceMeshRequest.h"

#include <Inventor/SbBox.h>
#include <Inventor/SbColor.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/actions/SoSubAction.h>

#include <stddef.h>
#include <stdint.h>
#include <vector>

class BObolLodService;
class SoBRLDatabaseSource;
class SoBRLMeshShape;
struct db_i;

class BOBOL_EXPORT SoBRLExportAction : public SoAction {
    typedef SoAction inherited;

    SO_ACTION_HEADER(SoBRLExportAction);

public:
    enum GeometryPolicy {
	FULL_DETAIL = 0,
	DISPLAY_LEVEL = 1
    };

    enum ObjectQueryFlag {
	QUERY_VISIBLE_ONLY = 0x01u,
	QUERY_DATABASE_OBJECTS = 0x02u,
	QUERY_VIEW_OBJECTS = 0x04u,
	QUERY_LOCAL_ONLY = 0x08u
    };

    struct LineRecord {
	size_t sequence;
	SbString path;
	SbString sourceName;
	SbString sourceType;
	uint32_t sourceId;
	SbString displayName;
	SbString geometryName;
	SbString cacheIdentity;
	SbString sourceIdentity;
	uint64_t cacheIdentityValue;
	uint64_t sourceIdentityValue;
	int databaseIntent;
	int overlayIntent;
	int hudIntent;
	int localSource;
	int sharedSource;
	int nonDatabaseSource;
	int drawMode;
	SbString recordRole;
	SbString geometryKind;
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
	SbString editIntentId;
	SbString editIntentRole;
	uint32_t lodPolicy;
	int colorOverride;
	SbColor color;
	SbVec3f a;
	SbVec3f b;
    };

    struct PointRecord {
	size_t sequence;
	SbString path;
	SbString sourceName;
	SbString sourceType;
	uint32_t sourceId;
	SbString displayName;
	SbString geometryName;
	SbString cacheIdentity;
	SbString sourceIdentity;
	uint64_t cacheIdentityValue;
	uint64_t sourceIdentityValue;
	int databaseIntent;
	int overlayIntent;
	int hudIntent;
	int localSource;
	int sharedSource;
	int nonDatabaseSource;
	int drawMode;
	SbString recordRole;
	SbString geometryKind;
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
	SbString editIntentId;
	SbString editIntentRole;
	uint32_t lodPolicy;
	int colorOverride;
	SbColor color;
	int pointColorValid;
	SbColor pointColor;
	int pointScaleValid;
	float pointScale;
	int pointNormalValid;
	SbVec3f pointNormal;
	SbVec3f point;
    };

    struct TriangleRecord {
	size_t sequence;
	SbString path;
	SbString sourceName;
	SbString sourceType;
	uint32_t sourceId;
	SbString displayName;
	SbString geometryName;
	SbString cacheIdentity;
	SbString sourceIdentity;
	uint64_t cacheIdentityValue;
	uint64_t sourceIdentityValue;
	int databaseIntent;
	int overlayIntent;
	int hudIntent;
	int localSource;
	int sharedSource;
	int nonDatabaseSource;
	int drawMode;
	SbString recordRole;
	SbString geometryKind;
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
	SbString editIntentId;
	SbString editIntentRole;
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

    struct ObjectRecord {
	size_t sequence;
	SbString path;
	SbString sourceName;
	SbString sourceType;
	uint32_t sourceId;
	SbString displayName;
	SbString geometryName;
	SbString cacheIdentity;
	SbString sourceIdentity;
	uint64_t cacheIdentityValue;
	uint64_t sourceIdentityValue;
	int databaseIntent;
	int overlayIntent;
	int hudIntent;
	int localSource;
	int sharedSource;
	int nonDatabaseSource;
	int drawMode;
	SbString recordRole;
	SbString geometryKind;
	int selected;
	int highlighted;
	int visible;
	int colorOverride;
	SbColor color;
	SbBox3f bounds;
	std::vector<int> lineIndices;
	std::vector<int> pointIndices;
	std::vector<int> triangleIndices;
    };

    enum ObjectLineCommand {
	LINE_MOVE = 0,
	LINE_DRAW = 1
    };

    struct ObjectLineSummary {
	int valid;
	size_t pointCount;
	size_t segmentCount;
	SbString cacheIdentity;
	SbString sourceIdentity;
	uint64_t cacheIdentityValue;
	uint64_t sourceIdentityValue;
    };

    struct ObjectPointSummary {
	int valid;
	size_t pointCount;
	SbString cacheIdentity;
	SbString sourceIdentity;
	uint64_t cacheIdentityValue;
	uint64_t sourceIdentityValue;
    };

    struct ObjectSurfaceSummary {
	int valid;
	size_t pointCount;
	size_t normalCount;
	size_t indexCount;
	size_t faceCount;
	int normalsPerIndex;
	int materialColorValid;
	SbColor materialColor;
	int materialDrawMode;
	int materialHighlighted;
	SbString cacheIdentity;
	SbString sourceIdentity;
	uint64_t cacheIdentityValue;
	uint64_t sourceIdentityValue;
    };

    SoBRLExportAction(void);
    virtual ~SoBRLExportAction(void);
    static void initClass(void);
    static uint64_t identityValue(const char *identity);
    static uint64_t identityValue(const SbString &identity);

    int getLineCount(void) const;
    const LineRecord &getLine(int index) const;
    int getPointCount(void) const;
    const PointRecord &getPoint(int index) const;
    int getTriangleCount(void) const;
    const TriangleRecord &getTriangle(int index) const;
    int getObjectRecordCount(void) const;
    const ObjectRecord &getObjectRecord(int index) const;
    int collectObjectRecords(std::vector<ObjectRecord> &records,
	    unsigned int queryFlags = 0,
	    const char *glob = 0,
	    int drawMode = -1) const;
    SbBool getObjectRecordLineSummary(const ObjectRecord &record,
	    ObjectLineSummary &summary) const;
    SbBool getObjectRecordLinePoint(const ObjectRecord &record,
	    size_t index, SbVec3f &point) const;
    SbBool getObjectRecordLineCommand(const ObjectRecord &record,
	    size_t index, int &command) const;
    SbBool getObjectRecordPointSummary(const ObjectRecord &record,
	    ObjectPointSummary &summary) const;
    SbBool getObjectRecordPoint(const ObjectRecord &record,
	    size_t index, SbVec3f &point) const;
    SbBool getObjectRecordSurfaceSummary(const ObjectRecord &record,
	    ObjectSurfaceSummary &summary) const;
    SbBool getObjectRecordSurfaceDetail(const ObjectRecord &record,
	    std::vector<SbVec3f> *surfacePoints,
	    std::vector<int> *surfaceIndices) const;
    SbBool getObjectRecordSurfacePoint(const ObjectRecord &record,
	    size_t index, SbVec3f &point) const;
    SbBool getObjectRecordSurfaceIndex(const ObjectRecord &record,
	    size_t index, int &out) const;
    const SbBox3f &getBounds(void) const;
    void setRecordStorageEnabled(SbBool enabled);
    SbBool isRecordStorageEnabled(void) const;
    void setGeometryPolicy(GeometryPolicy policy);
    GeometryPolicy getGeometryPolicy(void) const;
    unsigned int getSkippedLodDisplayMeshCount(void) const;
    int getSourceBackedFullDetailRequestCount(void) const;
    const BObolSourceMeshRequest &getSourceBackedFullDetailRequest(int index) const;
    SbBool makeSourceBackedFullDetailLodRequest(int index,
	    BObolLodRequest &request,
	    const BObolLodRequest *templateRequest = 0) const;
    SbBool appendSourceBackedFullDetailResult(
	    const BObolSourceMeshRequest &sourceRequest,
	    const BObolLodResult &result);
    int submitSourceBackedFullDetailRequests(BObolLodService *service,
	    uint64_t generation, struct db_i *dbip,
	    const BObolLodRequest *templateRequest = 0,
	    uint64_t maxFullDetailFaceCount = 0,
	    uint64_t maxFullDetailPointCount = 0) const;
    int consumeSourceBackedFullDetailResults(
	    const std::vector<BObolLodResult> &results,
	    const BObolLodRequest *templateRequest = 0);

protected:
    virtual void beginTraversal(SoNode *node);

private:
    friend class SoBRLDatabaseSource;
    static void nodeAction(SoAction *action, SoNode *node);
    static void databaseSourceAction(SoAction *action, SoNode *node);
    static void vlistShapeAction(SoAction *action, SoNode *node);
    static void meshShapeAction(SoAction *action, SoNode *node);

    void resetResults(void);
    void invalidateObjectRecords(void);
    void rebuildObjectRecords(void) const;
    void appendSourceBackedFullDetailRequest(const SoBRLMeshShape *shape,
	    const SbMatrix &localToWorld);
    void appendLineSummary(const SbVec3f &a, const SbVec3f &b);
    void appendPointSummary(int pointScaleValid, float pointScale,
	    const SbVec3f &point);
    void appendTriangleSummary(const SbVec3f &a, const SbVec3f &b,
	    const SbVec3f &c);
    void appendLine(const SbString &path, const SbString &sourceName,
	    const SbString &sourceType, uint32_t sourceId,
	    int regionId, int airCode, int materialId, int los,
	    int materialColorValid, const SbColor &materialColor,
	    const SbString &materialShader, int primitiveIndex,
	    int selected, int highlighted, int ghosted,
	    int hiddenLine, int editEmphasis,
	    const SbString &editIntentId,
	    const SbString &editIntentRole,
	    uint32_t lodPolicy,
	    int colorOverride, const SbColor &color,
	    const SbVec3f &a, const SbVec3f &b);
    void appendPoint(const SbString &path, const SbString &sourceName,
	    const SbString &sourceType, uint32_t sourceId,
	    int regionId, int airCode, int materialId, int los,
	    int materialColorValid, const SbColor &materialColor,
	    const SbString &materialShader, int primitiveIndex,
	    int selected, int highlighted, int ghosted,
	    int hiddenLine, int editEmphasis,
	    const SbString &editIntentId,
	    const SbString &editIntentRole,
	    uint32_t lodPolicy,
	    int colorOverride, const SbColor &color,
	    int pointColorValid, const SbColor &pointColor,
	    int pointScaleValid, float pointScale,
	    int pointNormalValid, const SbVec3f &pointNormal,
	    const SbVec3f &point);
    void appendTriangle(const SbString &path, const SbString &sourceName,
	    const SbString &sourceType, uint32_t sourceId,
	    int regionId, int airCode, int materialId, int los,
	    int materialColorValid, const SbColor &materialColor,
	    const SbString &materialShader, int primitiveIndex,
	    int vertexIndexA, int vertexIndexB, int vertexIndexC,
	    int selected, int highlighted, int ghosted,
	    int hiddenLine, int editEmphasis,
	    const SbString &editIntentId,
	    const SbString &editIntentRole,
	    uint32_t lodPolicy,
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
    mutable std::vector<ObjectRecord> objects;
    std::vector<BObolSourceMeshRequest> sourceBackedFullDetailRequests;
    SbBox3f bounds;
    int lineCount;
    int pointCount;
    int triangleCount;
    GeometryPolicy geometryPolicy;
    SbBool recordStorageEnabled;
    unsigned int skippedLodDisplayMeshCount;
    size_t recordSequence;
    mutable SbBool objectRecordsDirty;
};

#endif /* BOBOL_BEXPORTACTION_H */
