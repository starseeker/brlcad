/*                    M E S H _ S H A P E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/mesh_shape.h */

#ifndef BRLOBOL_MESH_SHAPE_H
#define BRLOBOL_MESH_SHAPE_H

#include "brlobol/defines.h"
#include "brlobol/source_mesh_request.h"

#include <Inventor/fields/SoMFInt32.h>
#include <Inventor/fields/SoMFBool.h>
#include <Inventor/fields/SoMFString.h>
#include <Inventor/fields/SoMFVec3f.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFColor.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFInt32.h>
#include <Inventor/fields/SoSFMatrix.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFUInt32.h>
#include <Inventor/fields/SoSFVec3f.h>
#include <Inventor/SbBox.h>
#include <Inventor/nodes/SoShape.h>

#include <stddef.h>
#include <stdint.h>
#include <vector>

class SoPickedPoint;
class SoGLRenderAction;
class SoRayPickAction;
struct BRLObolLodRequest;
struct BRLObolLodResult;

class BRLOBOL_EXPORT SoBRLMeshShape : public SoShape {
    typedef SoShape inherited;

    SO_NODE_HEADER(SoBRLMeshShape);

public:
    enum PickGeometryPolicy {
	PICK_DISPLAY_LEVEL = 0,
	PICK_FULL_DETAIL = 1
    };

    SoMFVec3f point;
    SoMFInt32 coordIndex;
    SoSFString sourcePath;
    SoSFString sourceName;
    SoSFString sourceType;
    SoSFUInt32 sourceId;
    SoSFString displayName;
    SoSFString geometryName;
    SoSFString cacheIdentity;
    SoSFString sourceIdentity;
    SoSFString ownerSourcePath;
    SoSFString ownerSourceInstanceKey;
    SoSFUInt32 ownerSourceRevision;
    SoSFUInt32 ownerInputsRevision;
    SoSFUInt32 ownerViewRevision;
    SoSFUInt32 ownerRealizedRevision;
    SoSFUInt32 ownerRealizedSourceRevision;
    SoSFUInt32 ownerRealizedInputsRevision;
    SoSFUInt32 ownerRealizedViewRevision;
    SoSFInt32 ownerRealizationStatus;
    SoSFString ownerRealizationDiagnostic;
    SoSFString ownerRealizationIdentity;
    SoSFBool ownerSourceStale;
    SoSFUInt32 ownerStaleReason;
    SoSFBool databaseIntent;
    SoSFBool overlayIntent;
    SoSFBool hudIntent;
    SoSFBool localSource;
    SoSFBool sharedSource;
    SoSFBool nonDatabaseSource;
    SoSFInt32 drawMode;
    SoSFString recordRole;
    SoSFString geometryKind;
    SoSFInt32 regionId;
    SoSFInt32 airCode;
    SoSFInt32 materialId;
    SoSFInt32 los;
    SoSFBool materialColorValid;
    SoSFColor materialColor;
    SoSFString materialShader;
    SoSFUInt32 materialRevision;
    SoSFBool drawMatrixValid;
    SoSFMatrix drawMatrix;
    SoSFBool drawCenterValid;
    SoSFVec3f drawCenter;
    SoSFBool drawSizeValid;
    SoSFFloat drawSize;
    SoSFBool colorOverride;
    SoSFColor color;
    SoSFColor selectedColor;
    SoSFColor highlightedColor;
    SoSFColor ghostedColor;
    SoSFBool visible;
    SoSFBool selectable;
    SoSFBool selected;
    SoSFBool highlighted;
    SoSFBool ghosted;
    SoSFInt32 lineStyle;
    SoSFInt32 lineWidth;
    SoSFFloat transparency;
    SoSFBool hiddenLine;
    SoSFBool editEmphasis;
    SoSFString editIntentId;
    SoSFString editIntentRole;
    SoSFUInt32 lodPolicy;
    SoSFBool lodAvailable;
    SoSFInt32 lodActiveLevel;
    SoSFUInt32 lodFaceCount;
    SoSFUInt32 lodPointCount;
    SoSFUInt32 lodOriginalPointCount;
    SoSFUInt32 lodNormalCount;
    SoSFBool lodHasSnappedPoints;
    SoSFBool lodHasNormals;
    SoSFVec3f lodBoundsMin;
    SoSFVec3f lodBoundsMax;
    SoSFBool lodStagedAvailable;
    SoSFInt32 lodResultKind;
    SoSFInt32 lodQualityTier;
    SoSFInt32 lodProviderStatus;
    SoSFString lodCacheKey;
    SoSFString lodProviderId;
    SoSFString lodProviderVersion;
    SoSFString lodDiagnostic;
    SoMFString lodDependencyPath;
    SoMFString lodDependencyName;
    SoMFString lodDependencySourceRevision;
    SoMFString lodDependencySourceHash;
    SoMFInt32 lodDependencyQualityTier;
    SoMFBool lodDependencyOptional;
    SoMFString lodAttributeName;
    SoMFString lodAttributeValue;
    SoSFInt32 lodProxyKind;
    SoSFVec3f lodProxyCenter;
    SoSFVec3f lodProxyAxisX;
    SoSFVec3f lodProxyAxisY;
    SoSFVec3f lodProxyAxisZ;
    SoSFVec3f lodProxyHalfExtents;
    SoSFFloat lodProxyRadius;
    SoSFBool lodBacked;
    SoSFBool lodDisplayActive;
    SoSFBool lodPreserveFullDetail;
    SoSFInt32 pickGeometryPolicy;
    SoMFInt32 selectedPrimitive;
    SoMFInt32 highlightedPrimitive;

    SoBRLMeshShape(void);
    static void initClass(void);

    void setIndexedTriangles(const SbVec3f *points, int pointCount,
	    const int32_t *indices, int indexCount);
    int getTriangleCount(void) const;
    SbBool getTriangle(int triangleIndex, SbVec3f &a, SbVec3f &b, SbVec3f &c) const;
    SbBool getTriangleVertexIndices(int triangleIndex, int &indexA, int &indexB, int &indexC) const;
    SbBool hasFullDetailMesh(void) const;
    int getFullDetailTriangleCount(void) const;
    SbBool getFullDetailTriangle(int triangleIndex, SbVec3f &a, SbVec3f &b, SbVec3f &c) const;
    SbBool getFullDetailTriangleVertexIndices(int triangleIndex, int &indexA, int &indexB, int &indexC) const;
    SbBool makeSourceMeshRequest(BRLObolSourceMeshRequest &request) const;
    SbBool needsSourceBackedFullDetail(void) const;
    size_t estimateDisplayMeshBytes(void) const;
    size_t estimateFullDetailMeshBytes(void) const;
    size_t estimateResidentMeshBytes(void) const;
    size_t evictFullDetailMesh(void);
    size_t evictActiveDisplayMesh(void);
    SbBool isPrimitiveSelected(int primitiveIndex) const;
    SbBool isPrimitiveHighlighted(int primitiveIndex) const;
    void setLodBackedMesh(SbBool enabled);
    SbBool isLodBackedMesh(void) const;
    void setLodPreserveFullDetail(SbBool enabled);
    SbBool isLodPreserveFullDetailEnabled(void) const;
    void setPickGeometryPolicy(PickGeometryPolicy policy);
    PickGeometryPolicy getPickGeometryPolicy(void) const;
    SbBool isLodDisplayActive(void) const;
    void makeLodRequest(BRLObolLodRequest &request,
	    const char *databaseId,
	    uint64_t databaseRevision,
	    uint64_t viewRevision,
	    uint64_t policyRevision,
	    int requestDrawMode,
	    const char *providerId,
	    const char *providerVersion,
	    int qualityTier) const;
    SbBool applyStagedLodResult(const BRLObolLodResult &result,
	    const BRLObolLodRequest *expectedRequest = NULL);
    void clearStagedLodResult(void);

    void GLRender(SoGLRenderAction *action) override;
    void computeBBox(SoAction *action, SbBox3f &box, SbVec3f &center) override;

protected:
    virtual ~SoBRLMeshShape(void);
    void generatePrimitives(SoAction *action) override;
    SoDetail *createTriangleDetail(SoRayPickAction *action,
	    const SoPrimitiveVertex *v1,
	    const SoPrimitiveVertex *v2,
	    const SoPrimitiveVertex *v3,
	    SoPickedPoint *pp) override;

private:
    void setIndexedTriangleFields(const SbVec3f *points, int pointCount,
	    const int32_t *indices, int indexCount);
    void updateSourceMeshMetricsFromFields(void);
    void updateSourceMeshMetricsFromFullDetail(void);
    void captureFullDetailMesh(void);
    void restoreFullDetailMesh(void);
    void clearFullDetailMesh(void);

    std::vector<SbVec3f> fullDetailPoint;
    std::vector<int32_t> fullDetailCoordIndex;
    uint64_t sourceMeshFaceCount;
    uint64_t sourceMeshPointCount;
    SbBox3f sourceMeshBounds;
    SbBool sourceMeshMetricsValid;
};

#endif /* BRLOBOL_MESH_SHAPE_H */
