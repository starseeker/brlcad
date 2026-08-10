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
#include "vmath.h"

#include <Inventor/SbBox.h>
#include <Inventor/SbColor.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec3f.h>

#include <memory>
#include <stddef.h>
#include <stdint.h>

struct rt_bot_internal;

/* A short-lived, immutable source import retained by the cold occurrence
 * stream.  Requests carry only a weak reference: the standing scene index
 * must never pin full source arrays after PoP data is resident.  The opaque
 * owner releases the rt_db_internal when the bounded stream lease expires. */
struct BOBOL_EXPORT BObolStagedSourceMesh {
    std::shared_ptr<void> owner;
    const struct rt_bot_internal *bot;
    /* Generic immutable triangle source.  These arrays are owned by owner
     * and permit analytic/BREP providers to use the same bounded cold-path
     * handoff as BoTs without retaining database internals or callbacks. */
    const point_t *points;
    const vect_t *normals;
    const int *faces;
    size_t pointCount;
    size_t faceCount;
    uint64_t contentKey;
    int shadedCullBackfaces;
    SbString assetName;
    uint32_t sourceRevision;
    size_t byteCount;

    BObolStagedSourceMesh(void) :
	owner(), bot(NULL), points(NULL), normals(NULL), faces(NULL),
	pointCount(0), faceCount(0), contentKey(0),
	shadedCullBackfaces(0), assetName(""), sourceRevision(0),
	byteCount(0)
    {
    }

    bool isValid(void) const
    {
	return owner && assetName.getLength() > 0 &&
	    (bot || (points && faces && pointCount && faceCount));
    }
};

struct BOBOL_EXPORT BObolSourceMeshRequest {
    SbString path;
    SbString sourceName;
    SbString sourceType;
    uint32_t sourceId;
    SbString displayName;
    SbString geometryName;
    SbString cacheIdentity;
    SbString sourceIdentity;
    /* The occurrence identity above is what picking/export reports.  These
     * fields identify the canonical database mesh whose arrays back the
     * request.  They differ after rigid transformed-copy reuse (for example,
     * an xpush result): all occurrences retain their own semantic identity
     * while loading one shared progressive asset. */
    SbString meshAssetPath;
    SbString meshAssetName;
    SbBox3f meshAssetBounds;
    /* Immutable representation payload identity.  Zero retains the named
     * authored-mesh lookup; nonzero reopens an exact BREP/tolerance variant
     * even if another view has since made a different variant the latest. */
    uint64_t meshAssetContentHash;
    /* Source-space tessellation contract for generated triangle assets.
     * Zero values identify authored meshes.  BREP requests retain the
     * canonical band here so a stable close-up may derive a finer immutable
     * variant without making the standing scene geometry camera-dependent. */
    double meshAssetTessellationAbsTol;
    double meshAssetTessellationRelTol;
    double meshAssetTessellationNormTol;
    /* Optional cold-path handoff.  This is deliberately transient and is
     * never part of persistent draw/cache identity. */
    std::weak_ptr<const BObolStagedSourceMesh> stagedSource;
    /* Maps canonical mesh-asset coordinates into this occurrence's
     * object-local coordinates.  This is independent of any temporary
     * display-proxy transform. */
    SbMatrix meshAssetTransform;
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
    float transparency;

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
	meshAssetPath = "";
	meshAssetName = "";
	meshAssetBounds.makeEmpty();
	meshAssetContentHash = 0;
	meshAssetTessellationAbsTol = 0.0;
	meshAssetTessellationRelTol = 0.0;
	meshAssetTessellationNormTol = 0.0;
	stagedSource.reset();
	meshAssetTransform.makeIdentity();
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
	transparency = 0.0f;
    }
};

#endif /* BOBOL_BSOURCEMESHREQUEST_H */
