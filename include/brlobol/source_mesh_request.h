/*              S O U R C E _ M E S H _ R E Q U E S T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/source_mesh_request.h */

#ifndef BRLOBOL_SOURCE_MESH_REQUEST_H
#define BRLOBOL_SOURCE_MESH_REQUEST_H

#include "brlobol/defines.h"

#include <Inventor/SbBox.h>
#include <Inventor/SbColor.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec3f.h>

#include <stdint.h>

struct BRLOBOL_EXPORT BRLObolSourceMeshRequest {
    SbString path;
    SbString sourceName;
    SbString sourceType;
    uint32_t sourceId;
    uint64_t faceCount;
    uint64_t pointCount;
    SbBox3f bounds;
    SbMatrix localToWorld;
    int queryBoundsValid;
    SbBox3f queryBounds;
    int queryRayValid;
    SbVec3f queryRayOrigin;
    SbVec3f queryRayDirection;
    int queryToleranceValid;
    float queryTolerance;
    int regionId;
    int airCode;
    int materialId;
    int los;
    int materialColorValid;
    SbColor materialColor;
    SbString materialShader;
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

    BRLObolSourceMeshRequest(void)
    {
	clear();
    }

    void clear(void)
    {
	path = "";
	sourceName = "";
	sourceType = "";
	sourceId = 0;
	faceCount = 0;
	pointCount = 0;
	bounds.makeEmpty();
	localToWorld.makeIdentity();
	queryBoundsValid = 0;
	queryBounds.makeEmpty();
	queryRayValid = 0;
	queryRayOrigin = SbVec3f(0.0f, 0.0f, 0.0f);
	queryRayDirection = SbVec3f(0.0f, 0.0f, -1.0f);
	queryToleranceValid = 0;
	queryTolerance = 0.0f;
	regionId = 0;
	airCode = 0;
	materialId = 0;
	los = 0;
	materialColorValid = 0;
	materialColor = SbColor(0.0f, 0.0f, 0.0f);
	materialShader = "";
	selected = 0;
	highlighted = 0;
	ghosted = 0;
	hiddenLine = 0;
	editEmphasis = 0;
	lodPolicy = 0;
	lodAvailable = 0;
	lodActiveLevel = -1;
	lodFaceCount = 0;
	lodPointCount = 0;
	lodOriginalPointCount = 0;
	lodNormalCount = 0;
	lodHasSnappedPoints = 0;
	lodHasNormals = 0;
	lodBoundsMin = SbVec3f(0.0f, 0.0f, 0.0f);
	lodBoundsMax = SbVec3f(0.0f, 0.0f, 0.0f);
	colorOverride = 0;
	color = SbColor(0.0f, 0.0f, 0.0f);
    }
};

#endif /* BRLOBOL_SOURCE_MESH_REQUEST_H */
