/*                  Q G O B O L E X P O R T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file qtcad/QgObolExport.h */

#ifndef QGOBOL_EXPORT_H
#define QGOBOL_EXPORT_H

#include "qtcad/defines.h"

#include <Inventor/SbColor.h>
#include <Inventor/SbVec3f.h>

#include <stddef.h>
#include <stdint.h>
#include <string>
#include <vector>

class QgView;

struct QTCAD_EXPORT QgObolExportTriangleRecord {
    QgObolExportTriangleRecord(void);

    std::string path;
    std::string sourceName;
    std::string sourceType;
    std::string displayName;
    std::string geometryName;
    std::string cacheIdentity;
    std::string sourceIdentity;
    uint64_t cacheIdentityValue;
    uint64_t sourceIdentityValue;
    std::string editIntentId;
    std::string editIntentRole;
    uint32_t sourceId;
    bool databaseIntent;
    bool overlayIntent;
    bool hudIntent;
    bool localSource;
    bool sharedSource;
    bool nonDatabaseSource;
    int drawMode;
    std::string recordRole;
    std::string geometryKind;
    int primitiveIndex;
    int vertexIndexA;
    int vertexIndexB;
    int vertexIndexC;
    SbVec3f a;
    SbVec3f b;
    SbVec3f c;
};

struct QTCAD_EXPORT QgObolExportGeometryRecord {
    QgObolExportGeometryRecord(void);

    int lineCount;
    int pointCount;
    int triangleCount;
    int submittedSourceRequestCount;
    SbVec3f boundsMin;
    SbVec3f boundsMax;
    bool boundsValid;
    bool sourceFullDetailPending;
    std::vector<QgObolExportTriangleRecord> triangles;
};

enum QgObolExportObjectQueryFlag {
    QG_OBOL_EXPORT_QUERY_VISIBLE_ONLY = 0x01u,
    QG_OBOL_EXPORT_QUERY_DATABASE_OBJECTS = 0x02u,
    QG_OBOL_EXPORT_QUERY_VIEW_OBJECTS = 0x04u,
    QG_OBOL_EXPORT_QUERY_LOCAL_ONLY = 0x08u
};

enum QgObolExportGeometryPolicy {
    QG_OBOL_EXPORT_FULL_DETAIL = 0,
    QG_OBOL_EXPORT_DISPLAY_LEVEL = 1
};

struct QTCAD_EXPORT QgObolExportObjectLineSummary {
    QgObolExportObjectLineSummary(void);

    bool valid;
    size_t pointCount;
    size_t segmentCount;
    std::string cacheIdentity;
    std::string sourceIdentity;
    uint64_t cacheIdentityValue;
    uint64_t sourceIdentityValue;
};

struct QTCAD_EXPORT QgObolExportObjectPointSummary {
    QgObolExportObjectPointSummary(void);

    bool valid;
    size_t pointCount;
    std::string cacheIdentity;
    std::string sourceIdentity;
    uint64_t cacheIdentityValue;
    uint64_t sourceIdentityValue;
};

struct QTCAD_EXPORT QgObolExportObjectSurfaceSummary {
    QgObolExportObjectSurfaceSummary(void);

    bool valid;
    size_t pointCount;
    size_t normalCount;
    size_t indexCount;
    size_t faceCount;
    int normalsPerIndex;
    bool materialColorValid;
    SbColor materialColor;
    int materialDrawMode;
    bool materialHighlighted;
    std::string cacheIdentity;
    std::string sourceIdentity;
    uint64_t cacheIdentityValue;
    uint64_t sourceIdentityValue;
};

struct QTCAD_EXPORT QgObolExportObjectRecord {
    QgObolExportObjectRecord(void);

    std::string path;
    std::string sourceName;
    std::string sourceType;
    std::string displayName;
    std::string geometryName;
    std::string cacheIdentity;
    std::string sourceIdentity;
    uint64_t cacheIdentityValue;
    uint64_t sourceIdentityValue;
    uint32_t sourceId;
    bool databaseIntent;
    bool overlayIntent;
    bool hudIntent;
    bool localSource;
    bool sharedSource;
    bool nonDatabaseSource;
    int drawMode;
    std::string recordRole;
    std::string geometryKind;
    bool selected;
    bool highlighted;
    bool visible;
    bool colorOverride;
    SbColor color;
    bool boundsValid;
    SbVec3f boundsMin;
    SbVec3f boundsMax;
    size_t lineCount;
    size_t pointPrimitiveCount;
    size_t triangleCount;
    QgObolExportObjectLineSummary lineSummary;
    QgObolExportObjectPointSummary pointSummary;
    QgObolExportObjectSurfaceSummary surfaceSummary;
    std::vector<SbVec3f> linePoints;
    std::vector<int> lineCommands;
    std::vector<SbVec3f> points;
    std::vector<SbVec3f> surfacePoints;
    std::vector<int> surfaceIndices;
};

struct QTCAD_EXPORT QgObolExportObjectQuery {
    QgObolExportObjectQuery(void);

    unsigned int flags;
    std::string glob;
    int drawMode;
    int geometryPolicy;
};

struct QTCAD_EXPORT QgObolExportObjectQueryResult {
    QgObolExportObjectQueryResult(void);

    int lineCount;
    int pointCount;
    int triangleCount;
    int submittedSourceRequestCount;
    bool sourceFullDetailPending;
    std::vector<QgObolExportObjectRecord> objects;
};

/**
 * Export Obol geometry using exact full-detail mesh policy.  LoD-backed
 * meshes with no resident full-detail payload consume matching ready
 * source-backed LoD results, or submit nonblocking source-backed full-detail
 * requests when a database source and LoD service are available.
 */
QTCAD_EXPORT int qg_obol_export_geometry_full_detail(QgView *display,
	QgObolExportGeometryRecord &record);

/**
 * Export grouped Obol object records using GED-compatible query semantics.
 * This is the qtcad-facing bridge for querying Obol/libbrlobol-owned record
 * data without exposing scene implementation details.
 */
QTCAD_EXPORT int qg_obol_export_object_records(QgView *display,
	const QgObolExportObjectQuery &query,
	QgObolExportObjectQueryResult &record);

#endif /* QGOBOL_EXPORT_H */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
