/*          B S O U R C E M E S H R E Q U E S T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BSourceMeshRequest.h */

#ifndef BOBOL_BSOURCEMESHREQUEST_H
#define BOBOL_BSOURCEMESHREQUEST_H

#include "BObol/BDefines.h"

#include <Inventor/SbBox.h>
#include <Inventor/SbColor.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec3f.h>

#include <stdint.h>

struct BOBOL_EXPORT BObolSourceMeshRequest {
    SbString path;
    SbString sourceName;
    SbString sourceType;
    uint32_t sourceId;
    SbString displayName;
    SbString geometryName;
    SbString cacheIdentity;
    SbString sourceIdentity;
    SbString ownerSourceInstanceKey;
    int databaseIntent;
    int overlayIntent;
    int hudIntent;
    int localSource;
    int sharedSource;
    int nonDatabaseSource;
    int drawMode;
    SbString recordRole;
    SbString geometryKind;
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

    BObolSourceMeshRequest(void)
    {
	clear();
    }

    void clear(void)
    {
	path = "";
	sourceName = "";
	sourceType = "";
	sourceId = 0;
	displayName = "";
	geometryName = "";
	cacheIdentity = "";
	sourceIdentity = "";
	ownerSourceInstanceKey = "";
	databaseIntent = 0;
	overlayIntent = 0;
	hudIntent = 0;
	localSource = 0;
	sharedSource = 0;
	nonDatabaseSource = 0;
	drawMode = 0;
	recordRole = "";
	geometryKind = "";
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
	editIntentId = "";
	editIntentRole = "";
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

#endif /* BOBOL_BSOURCEMESHREQUEST_H */
