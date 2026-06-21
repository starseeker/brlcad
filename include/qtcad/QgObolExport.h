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

#include <Inventor/SbVec3f.h>

#include <stdint.h>
#include <string>
#include <vector>

class QgView;

struct QTCAD_EXPORT QgObolExportTriangleRecord {
    QgObolExportTriangleRecord(void);

    std::string path;
    std::string sourceName;
    std::string sourceType;
    std::string editIntentId;
    std::string editIntentRole;
    uint32_t sourceId;
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

/**
 * Export Obol geometry using exact full-detail mesh policy.  LoD-backed
 * meshes with no resident full-detail payload consume matching ready
 * source-backed LoD results, or submit nonblocking source-backed full-detail
 * requests when a database source and LoD service are available.
 */
QTCAD_EXPORT int qg_obol_export_geometry_full_detail(QgView *display,
	QgObolExportGeometryRecord &record);

#endif /* QGOBOL_EXPORT_H */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
