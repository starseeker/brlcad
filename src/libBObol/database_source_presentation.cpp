/* D A T A B A S E _ S O U R C E _ P R E S E N T A T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file database_source_presentation.cpp
 *
 * View-local CAD presentation planning and publication for database sources.
 */

#include "common.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BLodRealization.h"
#include "BObol/BPickDetail.h"
#include "BObol/BViewLod.h"
#include "cad_assembly_private.h"
#include "cad_normals_private.h"
#include "cad_publication_private.h"
#include "compact_presentation_staging_private.h"
#include "transaction_fault_private.h"
#include "compact_occurrence_registry_private.h"
#include "database_source_private.h"
#include "database_source_presentation_private.h"

#include "bu/datetime.h"
#include "bu/log.h"
#include "bu/str.h"

#include <Inventor/actions/SoAction.h>
#include <Obol/cad/CadGeometryValidation.h>
#include <Obol/cad/CadProjectedProxy.h>

#include <algorithm>
#include <array>
#include <inttypes.h>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

static int
cad_proxy_corners(const BObolLodProxy &proxy, SbVec3f corners[8])
{
    if (!proxy.isValid())
	return 0;

    if (proxy.kind == BOBOL_LOD_PROXY_AABB) {
	const SbVec3f bmin = proxy.bounds.getMin();
	const SbVec3f bmax = proxy.bounds.getMax();
	for (size_t corner = 0; corner < 8; ++corner)
	    corners[corner].setValue(
		(corner & 1u) ? bmax[0] : bmin[0],
		(corner & 2u) ? bmax[1] : bmin[1],
		(corner & 4u) ? bmax[2] : bmin[2]);
	return 1;
    }

    if (proxy.kind == BOBOL_LOD_PROXY_OBB) {
	SbVec3f ax = proxy.axisX;
	SbVec3f ay = proxy.axisY;
	SbVec3f az = proxy.axisZ;
	if (ax.length() > 0.0f)
	    ax.normalize();
	if (ay.length() > 0.0f)
	    ay.normalize();
	if (az.length() > 0.0f)
	    az.normalize();
	ax *= proxy.halfExtents[0];
	ay *= proxy.halfExtents[1];
	az *= proxy.halfExtents[2];
	for (size_t corner = 0; corner < 8; ++corner)
	    corners[corner] = proxy.center +
		((corner & 1u) ? ax : -ax) +
		((corner & 2u) ? ay : -ay) +
		((corner & 4u) ? az : -az);
	return 1;
    }

    return 0;
}

int
bobol_cad_wire_geometry_from_binary_corners(const SbVec3f corners[8],
			       Obol::PartGeometryBuilder &geometry)
{
    static const int edges[12][2] = {
	{0, 1}, {2, 3}, {4, 5}, {6, 7},
	{0, 2}, {1, 3}, {4, 6}, {5, 7},
	{0, 4}, {1, 5}, {2, 6}, {3, 7}
    };

    Obol::WireRep wire;
    wire.bounds.makeEmpty();
    wire.segmentPoints.reserve(24);
    wire.segmentIds.reserve(12);
    for (size_t i = 0; i < 12; i++) {
	const SbVec3f &a = corners[edges[i][0]];
	const SbVec3f &b = corners[edges[i][1]];
	wire.segmentPoints.push_back(a);
	wire.segmentPoints.push_back(b);
	wire.segmentIds.push_back(static_cast<uint32_t>(i));
	wire.bounds.extendBy(a);
	wire.bounds.extendBy(b);
    }
    if (wire.bounds.isEmpty())
	return 0;
    geometry.wire = wire;
    return 1;
}

std::shared_ptr<const Obol::PartGeometry>
bobol_cad_structural_bounds_geometry(const SbBox3f &bounds,
	SbMatrix &geometryTransform, const SbVec3f orientedBounds[8])
{
    geometryTransform = SbMatrix::identity();
    if (bounds.isEmpty())
	return std::shared_ptr<const Obol::PartGeometry>();

    if (orientedBounds) {
	Obol::PartGeometryBuilder geometry;
	if (bobol_cad_wire_geometry_from_binary_corners(
		orientedBounds, geometry)) {
	    geometry.aggregateProxyCorners = std::array<SbVec3f, 8>{
		orientedBounds[0], orientedBounds[1],
		orientedBounds[2], orientedBounds[3],
		orientedBounds[4], orientedBounds[5],
		orientedBounds[6], orientedBounds[7]};
	    geometry.subpixelProxyEligible = true;
	    geometry.structuralProxy = true;
	    /* A cached OBB is an optional presentation optimization.  Float
	     * conversion or an older cache producer may leave its corners outside
	     * Obol's current orthogonality/conservativeness tolerance.  Reject it
	     * before the reporting publication boundary and use the authoritative
	     * AABB below.  Re-reporting the same optional miss on every retained
	     * scene rebuild both obscures real producer errors and adds avoidable
	     * owner-thread work. */
	    if (Obol::cadValidatePartGeometry(geometry)) {
		std::shared_ptr<const Obol::PartGeometry> oriented =
		    bobol_cad_build_geometry(
			std::move(geometry), "oriented structural bounds");
		if (oriented)
		    return oriented;
	    }
	}
    }

    const SbVec3f minimum = bounds.getMin();
    const SbVec3f maximum = bounds.getMax();
    const SbVec3f extent = maximum - minimum;
    if (extent[0] > SMALL_FASTF && extent[1] > SMALL_FASTF &&
	extent[2] > SMALL_FASTF) {
	static const std::shared_ptr<const Obol::PartGeometry> unitAabb = []() {
	    const SbVec3f corners[8] = {
		SbVec3f(0.0f, 0.0f, 0.0f), SbVec3f(1.0f, 0.0f, 0.0f),
		SbVec3f(0.0f, 1.0f, 0.0f), SbVec3f(1.0f, 1.0f, 0.0f),
		SbVec3f(0.0f, 0.0f, 1.0f), SbVec3f(1.0f, 0.0f, 1.0f),
		SbVec3f(0.0f, 1.0f, 1.0f), SbVec3f(1.0f, 1.0f, 1.0f)
	    };
	    Obol::PartGeometryBuilder geometry;
	    if (!bobol_cad_wire_geometry_from_binary_corners(corners, geometry))
		return std::shared_ptr<const Obol::PartGeometry>();
	    geometry.subpixelProxyEligible = true;
	    geometry.structuralProxy = true;
	    return bobol_cad_build_geometry(
		std::move(geometry), "unit structural bounds");
	}();
	if (!unitAabb)
	    return std::shared_ptr<const Obol::PartGeometry>();
	geometryTransform.setScale(extent);
	SbMatrix translation;
	translation.setTranslate(minimum);
	geometryTransform.multRight(translation);
	return unitAabb;
    }

    const SbVec3f corners[8] = {
	SbVec3f(minimum[0], minimum[1], minimum[2]),
	SbVec3f(maximum[0], minimum[1], minimum[2]),
	SbVec3f(minimum[0], maximum[1], minimum[2]),
	SbVec3f(maximum[0], maximum[1], minimum[2]),
	SbVec3f(minimum[0], minimum[1], maximum[2]),
	SbVec3f(maximum[0], minimum[1], maximum[2]),
	SbVec3f(minimum[0], maximum[1], maximum[2]),
	SbVec3f(maximum[0], maximum[1], maximum[2])
    };
    Obol::PartGeometryBuilder geometry;
    if (!bobol_cad_wire_geometry_from_binary_corners(corners, geometry))
	return std::shared_ptr<const Obol::PartGeometry>();
    geometry.subpixelProxyEligible = true;
    geometry.structuralProxy = true;
    return bobol_cad_build_geometry(
	std::move(geometry), "structural bounds");
}

static int
cad_shaded_geometry_from_corners(const SbVec3f corners[8],
	Obol::PartGeometryBuilder &geometry)
{
    static const uint32_t triangleIndices[36] = {
	0, 4, 6, 0, 6, 2,
	1, 3, 7, 1, 7, 5,
	0, 1, 5, 0, 5, 4,
	2, 6, 7, 2, 7, 3,
	0, 2, 3, 0, 3, 1,
	4, 5, 7, 4, 7, 6
    };
    Obol::TriMesh mesh;
    mesh.positions.assign(corners, corners + 8);
    mesh.indices.assign(triangleIndices, triangleIndices + 36);
    mesh.bounds.makeEmpty();
    for (const SbVec3f &corner : mesh.positions)
	mesh.bounds.extendBy(corner);
    if (mesh.bounds.isEmpty())
	return 0;
    geometry.shaded = std::move(mesh);
    geometry.shadedCullBackfaces = false;
    return 1;
}

static int
cad_proxy_part_geometry(const BObolLodProxy &proxy,
			int wire, int shaded,
			Obol::PartGeometryBuilder &geometry)
{
    SbVec3f corners[8];
    if ((!wire && !shaded) || !cad_proxy_corners(proxy, corners))
	return 0;
    if ((wire && !bobol_cad_wire_geometry_from_binary_corners(
		    corners, geometry)) ||
	(shaded && !cad_shaded_geometry_from_corners(corners, geometry)))
	return 0;
    /* An AABB is temporary structural coverage.  An OBB is a terminal
     * renderer representation of the retained CAD object and must therefore
     * retain its semantic selection/picking identity without participating
     * in structural-proxy lifecycle decisions. */
    geometry.subpixelProxyEligible = true;
    geometry.structuralProxy = proxy.kind == BOBOL_LOD_PROXY_AABB;
    if (proxy.kind == BOBOL_LOD_PROXY_OBB)
	geometry.aggregateProxyCorners = std::array<SbVec3f, 8>{
	    corners[0], corners[1], corners[2], corners[3],
	    corners[4], corners[5], corners[6], corners[7]};
    return 1;
}

static int
cad_mesh_payload_part_geometry(const BObolLodMeshPayload &payload,
			       int wire,
			       int shaded,
			       int shadedCullBackfaces,
			       Obol::PartGeometryBuilder &geometry,
			       BObolViewLodState::NormalStyle normalStyle =
				   BObolViewLodState::NORMAL_AUTHORED,
			       float normalCreaseAngle = 60.0f)
{
    if (!payload.isValid() || (!wire && !shaded))
	return 0;

    SbBox3f bounds;
    bounds.makeEmpty();
    for (size_t i = 0; i < payload.points.size(); i++)
	bounds.extendBy(payload.points[i]);
    if (bounds.isEmpty())
	return 0;

    if (shaded) {
	Obol::TriMesh mesh;
	mesh.positions = payload.points;
	mesh.bounds = bounds;
	mesh.indices.reserve(payload.coordIndex.size());
	for (size_t i = 0; i < payload.coordIndex.size(); i++) {
	    const int idx = payload.coordIndex[i];
	    if (idx < 0 || static_cast<size_t>(idx) >= payload.points.size())
		return 0;
	    mesh.indices.push_back(static_cast<uint32_t>(idx));
	}
	if (mesh.indices.empty())
	    return 0;
	std::vector<int32_t> cornerIndices;
	cornerIndices.reserve(payload.coordIndex.size());
	for (const int index : payload.coordIndex)
	    cornerIndices.push_back(index);
	std::vector<SbVec3f> cornerNormals;
	if (normalStyle == BObolViewLodState::NORMAL_FLAT) {
	    generate_flat_triangle_normals(cornerNormals, mesh.positions,
		cornerIndices);
	} else if (normalStyle == BObolViewLodState::NORMAL_SMOOTH) {
	    generate_smooth_triangle_normals(cornerNormals, mesh.positions,
		cornerIndices, normalCreaseAngle);
	} else if (payload.normals.size() == payload.coordIndex.size()) {
	    cornerNormals = payload.normals;
	}
	/* An indexed scan mesh without authored normals must stay indexed.
	 * Expanding flat corner normals here can turn N shared vertices into
	 * three vertices per face (tens of millions for Lucy), defeating PoP
	 * residency and creating multi-gigabyte transient maps.  Obol derives a
	 * flat face normal in the fragment shader, or per triangle in its
	 * immediate fallback, when this array remains empty. */
	if (!cornerNormals.empty()) {
	    sanitize_triangle_normals(cornerNormals, mesh.positions,
		cornerIndices);
	    if (!canonicalize_corner_normal_mesh(mesh, cornerNormals))
		return 0;
	}
	geometry.shaded = mesh;
	geometry.shadedCullBackfaces = shadedCullBackfaces != 0;
    }

    if (wire) {
	Obol::WireRep wireRep;
	wireRep.bounds = bounds;
	uint32_t edgeId = 0;
	wireRep.segmentPoints.reserve(payload.coordIndex.size() * 2);
	wireRep.segmentIds.reserve(payload.coordIndex.size());
	for (size_t i = 0; i + 2 < payload.coordIndex.size(); i += 3) {
	    int tri[3] = {
		payload.coordIndex[i],
		payload.coordIndex[i + 1],
		payload.coordIndex[i + 2]
	    };
	    for (int e = 0; e < 3; e++) {
		int a = tri[e];
		int b = tri[(e + 1) % 3];
		if (a < 0 || b < 0 ||
		    static_cast<size_t>(a) >= payload.points.size() ||
		    static_cast<size_t>(b) >= payload.points.size())
		    return 0;
		wireRep.segmentPoints.push_back(
		    payload.points[static_cast<size_t>(a)]);
		wireRep.segmentPoints.push_back(
		    payload.points[static_cast<size_t>(b)]);
		wireRep.segmentIds.push_back(edgeId++);
	    }
	}
	if (wireRep.segmentPoints.empty())
	    return 0;
	geometry.wire = wireRep;
    }

    /*
     * This is view-LoD presentation data with conservative full-mesh bounds.
     * When that entire extent is below one pixel, Obol may keep the retained
     * mesh and picking identity while drawing it through the aggregate point
     * channel.  Promotion is then a camera-local mask change, not a reload or
     * scene traversal.
     */
    geometry.subpixelProxyEligible = true;
    return 1;
}

static int
cad_progressive_mesh_part_geometry(
    const BObolLodProgressiveMeshPtr &progressive,
    int wire,
    int shaded,
    int shadedCullBackfaces,
    Obol::PartGeometryBuilder &geometry,
    BObolViewLodState::NormalStyle normalStyle =
	BObolViewLodState::NORMAL_AUTHORED,
    float normalCreaseAngle = 60.0f)
{
    if (!progressive || !progressive->isValid())
	return 0;
    const int residentCut = progressive->residentCut();
    /* residentCut is the highest population prefix loaded from storage,
     * not necessarily the finest coordinate cut the retained exact arrays
     * can draw.  Once all points and faces needed by later PoP cuts are
     * resident, those cuts differ only in the renderer's coordinate snap and
     * require no further I/O.  Advertise that drawable frontier to Obol;
     * otherwise it clamps a correctly requested cut back to the last
     * population-changing cut, leaving small
     * meshes visibly blocky forever. */
    int drawableCut = residentCut;
    for (int cut = residentCut + 1;
	 cut <= progressive->maximumCut(); ++cut) {
	if (!progressive->canDrawCut(cut))
	    break;
	drawableCut = cut;
    }
    BObolLodMeshPayload payload;
    if (!progressive->copyCut(payload, residentCut) ||
	!cad_mesh_payload_part_geometry(payload, wire, shaded,
	    shadedCullBackfaces, geometry, normalStyle, normalCreaseAngle))
	return 0;

    const int minimumCut = progressive->minimumCut();
    const SbVec3f quantizationMinimum =
	progressive->quantizationMinimum();
    const SbVec3f quantizationMaximum =
	progressive->quantizationMaximum();
    const SbBox3f conservativeBounds(
	quantizationMinimum, quantizationMaximum);
    if (geometry.shaded.has_value()) {
	Obol::TriMesh &mesh = *geometry.shaded;
	/*
	 * A resident PoP prefix does not necessarily touch the asset's full
	 * extent.  Its quantization domain does, and is therefore the stable
	 * conservative bound required by view culling and the aggregate
	 * subpixel proxy.  Prefix-derived bounds caused objects to cross the
	 * proxy threshold as data arrived, repeatedly rebuilding large retained
	 * submissions at an unchanged camera.
	 */
	if (!conservativeBounds.isEmpty())
	    mesh.bounds = conservativeBounds;
	mesh.progressiveMinimumCut =
	    static_cast<uint8_t>(std::max(0, minimumCut));
	mesh.progressiveResidentCut =
	    static_cast<uint8_t>(std::max(0, drawableCut));
	mesh.progressiveQuantizationMinimum = quantizationMinimum;
	mesh.progressiveQuantizationMaximum = quantizationMaximum;
	mesh.progressiveCuts.resize(
	    static_cast<size_t>(progressive->maximumCut()) + 1);
	for (int cut = 0; cut < BOBOL_MESH_LOD_CUT_COUNT_MAX; ++cut) {
	    if (static_cast<size_t>(cut) >= mesh.progressiveCuts.size())
		break;
	    const size_t count = progressive->faceCount(cut) * 3;
	    mesh.progressiveCuts[cut].indexCount = static_cast<uint32_t>(
		std::min(count,
		    static_cast<size_t>(std::numeric_limits<uint32_t>::max())));
	    BObolMeshLodCutInfo cutInfo = {};
	    if (!progressive->cutInfo(cut, &cutInfo))
		return 0;
	    mesh.progressiveCuts[cut].quantization = {
		cutInfo.quantization_bits[X], cutInfo.quantization_bits[Y],
		cutInfo.quantization_bits[Z]
	    };
	}
	/* Record the vertex extent of each cumulative index prefix once while
	 * constructing the geometry.  Recomputing it from every visible part on
	 * every render is quadratic bookkeeping for large distinct-mesh scenes.
	 * Derive it from the finalized Obol indices rather than the source point
	 * counts because normal canonicalization may split shared vertices. */
	size_t scannedIndexCount = 0;
	uint32_t maximumIndex = 0;
	for (size_t cut = 0; cut < mesh.progressiveCuts.size(); ++cut) {
	    const size_t indexCount = std::min<size_t>(
		mesh.progressiveCuts[cut].indexCount, mesh.indices.size());
	    for (; scannedIndexCount < indexCount; ++scannedIndexCount)
		maximumIndex = std::max(maximumIndex,
		    mesh.indices[scannedIndexCount]);
	    const size_t positionCount = indexCount ?
		static_cast<size_t>(maximumIndex) + 1 : 0;
	    mesh.progressiveCuts[cut].positionCount = static_cast<uint32_t>(
		std::min(positionCount,
		    static_cast<size_t>(std::numeric_limits<uint32_t>::max())));
	}
    }
    if (geometry.wire.has_value()) {
	Obol::WireRep &wireRep = *geometry.wire;
	if (!conservativeBounds.isEmpty())
	    wireRep.bounds = conservativeBounds;
	wireRep.progressiveMinimumCut =
	    static_cast<uint8_t>(std::max(0, minimumCut));
	wireRep.progressiveResidentCut =
	    static_cast<uint8_t>(std::max(0, drawableCut));
	wireRep.progressiveQuantizationMinimum = quantizationMinimum;
	wireRep.progressiveQuantizationMaximum = quantizationMaximum;
	wireRep.progressiveCuts.resize(
	    static_cast<size_t>(progressive->maximumCut()) + 1);
	for (size_t cut = 0; cut < wireRep.progressiveCuts.size(); ++cut) {
	    const size_t count = progressive->faceCount(cut) * 3;
	    wireRep.progressiveCuts[cut].segmentFirst = 0;
	    wireRep.progressiveCuts[cut].segmentCount = static_cast<uint32_t>(
		std::min(count,
		    static_cast<size_t>(std::numeric_limits<uint32_t>::max())));
	    BObolMeshLodCutInfo cutInfo = {};
	    if (!progressive->cutInfo(static_cast<int>(cut), &cutInfo))
		return 0;
	    wireRep.progressiveCuts[cut].quantization = {
		cutInfo.quantization_bits[X], cutInfo.quantization_bits[Y],
		cutInfo.quantization_bits[Z]
	    };
	}
    }
    return 1;
}

static bool
cad_progressive_geometry_can_draw(const Obol::PartGeometry *geometry,
				  int wire, int shaded, int activeCut)
{
    if (!geometry || activeCut < 0)
	return false;
    if (wire) {
	if (!geometry->wire.has_value() ||
	    !geometry->wire->isProgressive() ||
	    activeCut >= static_cast<int>(geometry->wire->progressiveCuts.size()) ||
	    geometry->wire->progressiveResidentCut < activeCut)
	    return false;
    }
    if (shaded) {
	if (!geometry->shaded.has_value() ||
	    !geometry->shaded->isProgressive() ||
	    activeCut >= static_cast<int>(geometry->shaded->progressiveCuts.size()) ||
	    geometry->shaded->progressiveResidentCut < activeCut)
	    return false;
    }
    return wire || shaded;
}

static std::string
cad_view_lod_assembly_key(const SoBRLDatabaseSource *source,
			  const BObolViewLodState::CadPayload *payload,
			  const SbMatrix &matrix)
{
    std::string key = source ? source->path.getValue().getString() : "";
    key += "|";
    key += source ? source->instanceKey.getValue().getString() : "";
    key += "|";
    key += payload ? payload->cacheKey.getString() : "";
    char buf[256] = {0};
    snprintf(buf, sizeof(buf), "|%d|%d|%llu|%llu|%u|%u|%d|%d|%d|%d|%d|%d|%d|%g|%d",
	     payload ? payload->resultKind : 0,
	     payload ? payload->qualityTier : 0,
	     payload ? (unsigned long long)payload->viewRevision : 0ULL,
	     payload ? (unsigned long long)payload->policyRevision : 0ULL,
	     source ? source->sourceRevision.getValue() : 0,
	     source ? source->inputsRevision.getValue() : 0,
	     source ? source->drawMode.getValue() : 0,
	     source ? source->visible.getValue() : 0,
	     source ? source->selected.getValue() : 0,
	     source ? source->highlighted.getValue() : 0,
	     source ? source->colorOverride.getValue() : 0,
	     source ? source->materialColorValid.getValue() : 0,
	     source ? source->lineWidth.getValue() : 0,
	     source ? source->transparency.getValue() : 0.0f,
	     source ? source->drawMatrixValid.getValue() : 0);
    key += buf;
    const SbColor selectedColor = source ?
	source->selectedColor.getValue() : SbColor(1.0f, 1.0f, 1.0f);
    const SbColor highlightedColor = source ?
	source->highlightedColor.getValue() : SbColor(1.0f, 1.0f, 0.0f);
    const SbColor ghostedColor = source ?
	source->ghostedColor.getValue() : SbColor(0.55f, 0.55f, 0.55f);
    for (int i = 0; i < 3; i++) {
	char cbuf[128] = {0};
	snprintf(cbuf, sizeof(cbuf), "|sc=%.9g|hc=%.9g|gc=%.9g",
		 selectedColor[i], highlightedColor[i], ghostedColor[i]);
	key += cbuf;
    }
    const float *m = matrix[0];
    for (int i = 0; i < 16; i++) {
	char mbuf[64] = {0};
	snprintf(mbuf, sizeof(mbuf), "|%.9g", m[i]);
	key += mbuf;
    }
    return key;
}

Obol::CadDrawMode
cad_presentation_draw_mode(int sourceDrawMode, size_t shadedCount,
    size_t wireCount)
{
    if (sourceDrawMode == BOBOL_LOD_DRAW_HIDDEN_LINE)
	return Obol::CadDrawMode::HiddenLine;
    if (shadedCount > 0 && wireCount > 0)
	return Obol::CadDrawMode::ShadedWithEdges;
    if (shadedCount > 0)
	return Obol::CadDrawMode::Shaded;
    return Obol::CadDrawMode::Wireframe;
}

static SoBRLCadAssembly *
cad_view_lod_assembly(const SoBRLDatabaseSource *source,
		      const BObolViewLodState::CadPayload *payload,
		      const BObolViewLodState *viewState)
{
    if (!source || !payload || !viewState || !payload->isValid() ||
	!source->visible.getValue())
	return NULL;

    const SbMatrix matrix = cad_instance_matrix(source, SbMatrix::identity());
    const BObolViewLodState::NormalStyle normalStyle =
	viewState->getNormalStyle();
    const float normalCreaseAngle = viewState->getNormalCreaseAngle();
    std::string assemblyKey =
	cad_view_lod_assembly_key(source, payload, matrix);
    assemblyKey += "|normal-style=";
    assemblyKey += std::to_string(static_cast<int>(normalStyle));
    assemblyKey += "|normal-crease=";
    assemblyKey += std::to_string(normalCreaseAngle);
    SbString cachedKey;
    SoBRLCadAssembly *assembly = dynamic_cast<SoBRLCadAssembly *>(
	viewState->findCadPresentation(source, &cachedKey));
    if (assembly &&
	bu_strcmp(cachedKey.getString(), assemblyKey.c_str()) == 0)
	return assembly;

    Obol::PartGeometryBuilder geometry;
    std::shared_ptr<const Obol::PartGeometry> preparedGeometry;
    int shadedCount = 0;
    int wireCount = 0;
    SoBRLPickDetail::PrimitiveKind primitiveKind =
	SoBRLPickDetail::LINE_SEGMENT;
    const int sourceDrawMode = source_record_draw_mode(source);
    const bool hiddenLine = sourceDrawMode == BOBOL_LOD_DRAW_HIDDEN_LINE;
    const int wire =
	(sourceDrawMode == BOBOL_LOD_DRAW_WIRE || hiddenLine ||
	 payload->drawMode == BOBOL_LOD_DRAW_WIRE) ? 1 : 0;
    const int shaded = hiddenLine || !wire ? 1 : 0;
    if (payload->resultKind == BOBOL_LOD_RESULT_AABB ||
	payload->resultKind == BOBOL_LOD_RESULT_PROXY) {
	if (!cad_proxy_part_geometry(payload->proxy, wire, shaded, geometry))
	    return NULL;
	wireCount = geometry.wire ? 1 : 0;
	shadedCount = geometry.shaded ? 1 : 0;
	if (shadedCount)
	    primitiveKind = SoBRLPickDetail::FACE;
    } else {
	if (payload->preparedCadGeometry && !payload->progressiveMesh) {
	    preparedGeometry = payload->preparedCadGeometry;
	} else if (payload->progressiveMesh &&
	    normalStyle == BObolViewLodState::NORMAL_AUTHORED &&
	    payload->preparedCadGeometry &&
	    payload->preparedCadGeometryRevision ==
		payload->progressiveMesh->revision() &&
	    cad_progressive_geometry_can_draw(
		payload->preparedCadGeometry.get(), wire, shaded,
		payload->activeCut))
	    preparedGeometry = payload->preparedCadGeometry;
	if (!preparedGeometry &&
	    !cad_mesh_payload_part_geometry(payload->mesh, wire, shaded,
		payload->shadedCullBackfaces,
		geometry, normalStyle, normalCreaseAngle))
	    return NULL;
	const bool hasWire = preparedGeometry ?
	    preparedGeometry->wire.has_value() : geometry.wire.has_value();
	const bool hasShaded = preparedGeometry ?
	    preparedGeometry->shaded.has_value() : geometry.shaded.has_value();
	if (hasWire)
	    wireCount = 1;
	if (hasShaded) {
	    shadedCount = 1;
	    primitiveKind = SoBRLPickDetail::FACE;
	}
    }

    std::string partKey = "view-lod:";
    partKey += source->path.getValue().getString();
    partKey += ":";
    partKey += payload->cacheKey.getString();
    partKey += ":";
    partKey += payload->resultKind == BOBOL_LOD_RESULT_AABB ? "aabb" :
	       (payload->resultKind == BOBOL_LOD_RESULT_PROXY ? "proxy" :
		(wireCount ? "wire-mesh" : "mesh"));
    Obol::PartId partId = Obol::CadIdBuilder::partId(partKey);
    Obol::PartUpdate part;
    part.part = partId;
    const std::shared_ptr<const Obol::PartGeometry> sharedGeometry =
	preparedGeometry ? preparedGeometry :
	bobol_cad_build_geometry(
	    std::move(geometry), "single-view LoD staging");
    if (!bobol_cad_admit_geometry(sharedGeometry, part.geometry,
	    "single-view LoD staging"))
	return NULL;
    std::vector<Obol::PartUpdate> parts(1, std::move(part));

    std::string instanceKey =
	cad_instance_key(source, source->path.getValue().getString(), 0);
    instanceKey += ":view-lod:";
    instanceKey += payload->cacheKey.getString();
    Obol::InstanceId instanceId =
	Obol::CadIdBuilder::instanceId(instanceKey);
    Obol::InstanceRecord record;
    record.part = partId;
    record.localToRoot = matrix;
    record.parent = cad_source_parent_instance(source);
    record.childName = cad_source_leaf_name(source);
    record.occurrenceIndex = source->occurrenceIndex.getValue();
    record.boolOp = cad_source_boolean_operation(source);
    record.lodStructuralProxy =
	payload->resultKind == BOBOL_LOD_RESULT_AABB ||
	(payload->resultKind == BOBOL_LOD_RESULT_PROXY &&
	 payload->proxy.kind == BOBOL_LOD_PROXY_AABB);
    record.style = cad_source_style(source);
    Obol::InstanceUpdate update;
    update.instance = instanceId;
    update.record = record;
    std::vector<Obol::InstanceUpdate> instances(1, update);
    if (!assembly) {
	assembly = new SoBRLCadAssembly;
	viewState->setCadPresentation(source, assembly);
    }
    if (!bobol_cad_replace_scene(assembly, parts, instances,
	    "single-view LoD replacement"))
	return NULL;
    assembly->clearSemanticMap();
    assembly->setInstanceSemantic(instanceId,
	cad_source_semantic(source, primitiveKind));

    assembly->setPresentationDrawMode(cad_presentation_draw_mode(
	source_record_draw_mode(source), shadedCount, wireCount));
    if (payload->resultKind == BOBOL_LOD_RESULT_AABB ||
	payload->resultKind == BOBOL_LOD_RESULT_PROXY)
	assembly->setPresentationPickMode(Obol::CadPickMode::Bounds);
    else
	assembly->setPresentationPickMode(Obol::CadPickMode::Automatic);
    viewState->setCadPresentation(source, assembly, assemblyKey.c_str());
    return assembly;
}

static bool
cad_lod_path_equal(const SbString &a, const SbString &b)
{
    const char *ap = a.getString();
    const char *bp = b.getString();
    if (!ap || !bp)
	return false;
    if (ap[0] == '/')
	ap++;
    if (bp[0] == '/')
	bp++;
    return bu_strcmp(ap, bp) == 0;
}

static const BObolViewLodState::CadPayload *
cad_compact_payload_for_entry(
    const std::vector<const BObolViewLodState::CadPayload *> &payloads,
    const BObolCompactInstanceEntry &entry)
{
    const SbString &occurrenceKey = compact_instance_identity(entry);
    for (const BObolViewLodState::CadPayload *payload : payloads) {
	if (payload && payload->sourceInstanceKey.getLength() > 0 &&
	    bu_strcmp(payload->sourceInstanceKey.getString(),
		occurrenceKey.getString()) == 0)
	    return payload;
    }
    /* Source-wide legacy results deliberately have no occurrence key.  They
     * may still bind by path, but an occurrence-specific result never aliases
     * a sibling merely because both paths are textual matches. */
    for (const BObolViewLodState::CadPayload *payload : payloads) {
	if (payload && payload->sourceInstanceKey.getLength() == 0 &&
	    cad_lod_path_equal(payload->sourcePath, entry.semantic.path))
	    return payload;
    }
    return NULL;
}

SoBRLCadAssembly *
SoBRLDatabaseSource::currentCompactViewLodAssembly(
    const BObolViewLodState *viewState) const
{
    if (!this->d->compactIndex || !viewState)
	return NULL;
    SoBRLCadAssembly *assembly = dynamic_cast<SoBRLCadAssembly *>(
	viewState->findCadPresentation(this));
    if (!assembly || !assembly->compactPresentationInitialized ||
	assembly->compactPresentationSourceRoutingId != this->d->routingId ||
	assembly->compactPresentationPopulationEpoch !=
	    this->d->compactPopulationEpoch ||
	assembly->compactPresentationIndex != this->d->compactIndex ||
	assembly->compactPresentationSourceRevision !=
	    this->sourceRevision.getValue() ||
	assembly->compactPresentationInputsRevision !=
	    this->inputsRevision.getValue() ||
	assembly->compactPresentationCadBatchRevision !=
	    this->cadBatchRevisionGet() ||
	assembly->compactPresentationDrawMode != source_record_draw_mode(this) ||
	assembly->compactPresentationPayloadRevision !=
	    viewState->cadRevision())
	return NULL;
    return assembly;
}

int
SoBRLDatabaseSource::getCompactViewLodSupersededFallbackCount(
    const BObolViewLodState *viewState, std::vector<SbString> *paths) const
{
    if (paths)
	paths->clear();
    if (!this->d->compactIndex || !viewState)
	return 0;

    SoBRLCadAssembly *assembly = dynamic_cast<SoBRLCadAssembly *>(
	viewState->findCadPresentation(this));
    if (!assembly || !assembly->compactPresentationInitialized ||
	assembly->compactPresentationSourceRoutingId != this->d->routingId ||
	assembly->compactPresentationPopulationEpoch !=
	    this->d->compactPopulationEpoch ||
	assembly->compactPresentationIndex != this->d->compactIndex)
	return 0;

    std::vector<const BObolViewLodState::CadPayload *> payloads;
    viewState->findCadPayloads(this, payloads);
    std::unordered_set<std::string_view> payloadKeys;
    payloadKeys.reserve(payloads.size());
    for (const BObolViewLodState::CadPayload *payload : payloads)
	if (payload && payload->isValid() &&
	    payload->sourceInstanceKey.getLength() > 0)
	    payloadKeys.emplace(payload->sourceInstanceKey.getString(),
		static_cast<size_t>(
		    payload->sourceInstanceKey.getLength()));

    int count = 0;
    for (const BObolCompactInstanceEntry &entry :
	 this->d->compactIndex->entries) {
	/* Geometry kind alone is not presentation semantics: authored AABB-like
	 * structure and the root overview remain useful scene content, whereas
	 * only an explicitly marked temporary LoD leaf fallback is an obligation
	 * that may request further refinement. */
	if (!entry.lodBacked || !entry.geometry ||
	    !entry.geometry->structuralProxy)
	    continue;
	const char *kind = entry.shapeSummary.geometryKind.getString();
	if (!BU_STR_EQUAL(kind, "aabb") && !BU_STR_EQUAL(kind, "obb"))
	    continue;
	const SbString &occurrenceKey = compact_instance_identity(entry);
	const std::string_view occurrenceView(
	    occurrenceKey.getString(),
	    static_cast<size_t>(occurrenceKey.getLength()));
	if (payloadKeys.find(occurrenceView) == payloadKeys.end())
	    continue;
	const auto presentation =
	    assembly->compactInstancePresentations.find(entry.instance);
	if (presentation == assembly->compactInstancePresentations.end() ||
	    assembly->compactSpatiallyCulledInstances.find(entry.instance) !=
		assembly->compactSpatiallyCulledInstances.end() ||
	    presentation->second.activePart != entry.part)
	    continue;
	count++;
	if (paths)
	    paths->push_back(entry.semantic.path);
    }
    return count;
}

int
SoBRLDatabaseSource::getCompactViewLodActiveFallbackCount(
    const BObolViewLodState *viewState, std::vector<SbString> *paths) const
{
    if (paths)
	paths->clear();
    if (!this->d->compactIndex || !viewState)
	return 0;

    SoBRLCadAssembly *assembly = dynamic_cast<SoBRLCadAssembly *>(
	viewState->findCadPresentation(this));
    if (!assembly || !assembly->compactPresentationInitialized ||
	assembly->compactPresentationSourceRoutingId != this->d->routingId ||
	assembly->compactPresentationPopulationEpoch !=
	    this->d->compactPopulationEpoch ||
	assembly->compactPresentationIndex != this->d->compactIndex)
	return 0;

    int count = 0;
    for (const BObolCompactInstanceEntry &entry :
	 this->d->compactIndex->entries) {
	/* See getCompactViewLodSupersededFallbackCount: only temporary LoD
	 * fallback geometry drives convergence, not arbitrary AABB-like scene
	 * content. */
	if (!entry.lodBacked || !entry.geometry ||
	    !entry.geometry->structuralProxy)
	    continue;
	const char *kind = entry.shapeSummary.geometryKind.getString();
	if (!BU_STR_EQUAL(kind, "aabb") && !BU_STR_EQUAL(kind, "obb"))
	    continue;
	const auto presentation =
	    assembly->compactInstancePresentations.find(entry.instance);
	if (presentation == assembly->compactInstancePresentations.end() ||
	    assembly->compactSpatiallyCulledInstances.find(entry.instance) !=
		assembly->compactSpatiallyCulledInstances.end() ||
	    presentation->second.activePart != entry.part)
	    continue;
	count++;
	if (paths)
	    paths->push_back(entry.semantic.path);
    }
    return count;
}

/* A detached producer can finish its index before a view has committed that
 * index to its retained assembly.  Keep the temporary extent through that
 * gap: otherwise the producer may hide the only visible object one frame
 * before the replacement boxes/meshes exist in the rendered plan. */
void
SoBRLDatabaseSource::retireCompactOverviewAfterPresentation(
	SoBRLCadAssembly *assembly)
{
	if (!assembly || !this->d->compactIndex ||
	this->d->compactOverviewState !=
	    BObolCompactOccurrenceRegistryState::OverviewState::RetirementPending)
	return;
    BObolCompactInstanceIndex &index = *this->d->compactIndex;

    bool replacementPresented = false;
    for (const BObolCompactInstanceEntry &entry : index.entries) {
	if (BU_STR_EQUAL(entry.shapeSummary.recordRole.getString(),
		"lod-overview") || !entry.visible)
	    continue;
	const auto presentation =
	    assembly->compactInstancePresentations.find(entry.instance);
	if (presentation == assembly->compactInstancePresentations.end() ||
	    assembly->isInstanceHidden(entry.instance))
	    continue;
	/* The source record deliberately keeps its AABB while view-specific
	 * PoP geometry lives under a private payload part.  Either a payload
	 * rebind or a committed leaf part proves that the view now has useful
	 * per-leaf coverage.  The latter deliberately includes structural leaf
	 * boxes: the whole-target overview is only an early extent cue and must
	 * retire once the complete leaf frontier is presented. */
	const bool payloadGeometryPresented =
	    !presentation->second.payloadKey.empty() &&
	    presentation->second.activePart != entry.part;
	const bool leafGeometryPresented = entry.geometry &&
	    presentation->second.activePart == entry.part;
	if (!payloadGeometryPresented && !leafGeometryPresented)
	    continue;
	replacementPresented = true;
	break;
    }
    if (!replacementPresented)
	return;

    std::vector<size_t> retiredEntries;
    for (size_t i = 0; i < index.entries.size(); ++i) {
	BObolCompactInstanceEntry &entry = index.entries[i];
	if (!BU_STR_EQUAL(entry.shapeSummary.recordRole.getString(),
		"lod-overview") || !entry.visible)
	    continue;
	entry.authoredVisible = FALSE;
	entry.visible = FALSE;
	entry.visibilityRevision =
	    compact_next_revision(entry.visibilityRevision);
	compact_sync_shape_summary_state(entry);
	index.hiddenInstances.push_back(entry.instance);
	retiredEntries.push_back(i);
    }
    if (retiredEntries.empty())
	return;

    std::sort(index.hiddenInstances.begin(), index.hiddenInstances.end());
    index.hiddenInstances.erase(std::unique(index.hiddenInstances.begin(),
	index.hiddenInstances.end()), index.hiddenInstances.end());
    /* This presentation has already atomically installed its replacement
     * records, so its visibility update may retire the overview in the same
     * frame.  Other views will build the same replacement records before
     * observing the shared hidden state. */
    assembly->setHiddenInstances(index.hiddenInstances);
    this->d->compactOverviewState =
	BObolCompactOccurrenceRegistryState::OverviewState::Retired;
    this->markCompiledAssemblyDirty();
    this->markCadBatchDirty(retiredEntries);
}

static Obol::PartId
cad_compact_shared_lod_part_id(
    const BObolViewLodState::CadPayload &payload, int drawMode)
{
    std::string key(payload.progressiveMesh ?
	"compact-progressive-asset:" : "compact-terminal-asset:");
    key += payload.cacheKey.getString();
    key += ':';
    key += std::to_string(drawMode);
    return Obol::CadIdBuilder::partId(key);
}

template <typename Geometry>
static uint8_t
cad_compact_geometry_channels(const Geometry *geometry)
{
    if (!geometry)
	return 0;
    return (geometry->wire.has_value() ? 1u : 0u) |
	(geometry->shaded.has_value() ? 2u : 0u) |
	(geometry->points.has_value() ? 4u : 0u);
}

static Obol::PartId
cad_compact_layer_part_id(
    const BObolViewLodState::CadPayload &payload,
    const BObolLodPresentationLayer &layer, int drawMode)
{
    std::string key("compact-layer-part:");
    key += payload.cacheKey.getString();
    key += ':';
    key += layer.layerKey.getString();
    key += ':';
    key += std::to_string(drawMode);
    return Obol::CadIdBuilder::partId(key);
}

static Obol::InstanceId
cad_compact_layer_instance_id(Obol::InstanceId base,
    const BObolLodPresentationLayer &layer)
{
    std::string key("compact-layer-instance:");
    key += std::to_string(base.w0);
    key += ':';
    key += std::to_string(base.w1);
    key += ':';
    key += layer.layerKey.getString();
    return Obol::CadIdBuilder::instanceId(key);
}

static Obol::PartId
cad_compact_payload_part_id(
    const BObolViewLodState::CadPayload &payload,
    const SbString &occurrenceKey, const std::string &payloadKey,
    int drawMode)
{
    if (!payload.presentationLayers.empty())
	return cad_compact_layer_part_id(
	    payload, payload.presentationLayers.front(), drawMode);
    if ((payload.progressiveMesh && payload.progressiveMesh->isValid()) ||
	(payload.resultKind == BOBOL_LOD_RESULT_FULL_DETAIL &&
	 payload.preparedCadGeometry))
	return cad_compact_shared_lod_part_id(payload, drawMode);

    std::string partKey("compact-lod:");
    partKey.reserve(partKey.size() +
	static_cast<size_t>(occurrenceKey.getLength()) +
	payloadKey.size() + 2u);
    partKey += occurrenceKey.getString();
    partKey += ':';
    partKey += payloadKey;
    return Obol::CadIdBuilder::partId(partKey);
}

static void
cad_compact_append_layer_identity(
    std::string &key,
    const std::vector<BObolLodPresentationLayer> &layers)
{
    key += ":layers=";
    key += std::to_string(layers.size());
    for (const BObolLodPresentationLayer &layer : layers) {
	const char *layerKey = layer.layerKey.getString();
	const size_t length = static_cast<size_t>(layer.layerKey.getLength());
	key += ':';
	key += std::to_string(length);
	key += '#';
	key.append(layerKey, length);
	key += ':';
	key += std::to_string(layer.geometryRevision);
	key += ':';
	key += std::to_string(layer.activeCut);
	key += ':';
	key += layer.coverage ? '1' : '0';
    }
}

static void
cad_compact_adjust_channel_count(uint8_t channels, bool add,
    size_t &wireCount, size_t &shadedCount)
{
    const bool wire = channels & (1u | 4u);
    const bool shaded = channels & 2u;
    if (wire)
	wireCount = add ? wireCount + 1 : (wireCount ? wireCount - 1 : 0);
    if (shaded)
	shadedCount = add ? shadedCount + 1 :
	    (shadedCount ? shadedCount - 1 : 0);
}

static int
cad_compact_default_progressive_cut(
    const BObolCompactInstanceEntry &entry)
{
    if (!bobol_compact_geometry_is_resident_progressive(entry.geometry))
	return -1;
    const Obol::WireRep &wire = *entry.geometry->wire;
    return static_cast<int>(wire.progressiveMinimumCut);
}

SoBRLCadAssembly *
SoBRLDatabaseSource::compactViewLodAssembly(
    const std::vector<const BObolViewLodState::CadPayload *> &payloads,
    const BObolViewLodState *viewState)
    const
{
    const bool presentationTiming =
	getenv("BOBOL_COMPACT_PRESENTATION_TIMING") != NULL;
    const int64_t presentationStarted =
	presentationTiming ? bu_gettime() : 0;
    if (!this->d->compactIndex || !viewState)
	return NULL;

    SoBRLCadAssembly *assembly = dynamic_cast<SoBRLCadAssembly *>(
	viewState->findCadPresentation(this));
    const bool newAssembly = assembly == NULL;
    if (!assembly) {
	assembly = new SoBRLCadAssembly;
    }
    if (!assembly->reserveCompactPresentationCapacity(std::max(
	    this->d->compactExpectedInstanceCount,
	    this->d->compactIndex->entries.size()))) {
	if (newAssembly) {
	    assembly->ref();
	    assembly->unref();
	    return NULL;
	}
	return assembly;
    }
    if (newAssembly)
	viewState->setCadPresentation(this, assembly);

    const int sourceDrawMode = source_record_draw_mode(this);
    const uint64_t payloadRevision = viewState->cadRevision();
    const BObolViewLodState::NormalStyle normalStyle =
	viewState->getNormalStyle();
    const float normalCreaseAngle = viewState->getNormalCreaseAngle();
    if (assembly->compactPresentationInitialized &&
	assembly->compactPresentationSourceRoutingId == this->d->routingId &&
	assembly->compactPresentationPopulationEpoch ==
	    this->d->compactPopulationEpoch &&
	assembly->compactPresentationIndex == this->d->compactIndex &&
	assembly->compactPresentationSourceRevision ==
	    this->sourceRevision.getValue() &&
	assembly->compactPresentationInputsRevision ==
	    this->inputsRevision.getValue() &&
	assembly->compactPresentationCadBatchRevision ==
	    this->cadBatchRevisionGet() &&
	assembly->compactPresentationDrawMode == sourceDrawMode &&
	assembly->compactPresentationPayloadRevision == payloadRevision) {
	const_cast<SoBRLDatabaseSource *>(this)->
	    retireCompactOverviewAfterPresentation(assembly);
	return assembly;
    }
    const bool hiddenLine = sourceDrawMode == BOBOL_LOD_DRAW_HIDDEN_LINE;
    bool reset = !assembly->compactPresentationInitialized ||
	assembly->compactPresentationSourceRoutingId != this->d->routingId ||
	assembly->compactPresentationPopulationEpoch !=
	    this->d->compactPopulationEpoch ||
	assembly->compactPresentationIndex != this->d->compactIndex ||
	assembly->compactPresentationSourceRevision !=
	    this->sourceRevision.getValue() ||
	assembly->compactPresentationInputsRevision !=
	    this->inputsRevision.getValue() ||
	assembly->compactPresentationDrawMode != sourceDrawMode;

    /*
     * Streaming a large root appends/replaces a bounded occurrence batch in
     * the same compact index.  Clearing and reinserting the whole compiled
     * assembly on every batch made population O(n^2), invalidated every part
     * generation, and forced the renderer to re-snap every prior PoP cut.
     * Preserve the assembly when every old occurrence still exists and route
     * only new or geometry-changed records through the normal upsert path.
     * A removal or index replacement remains an authoritative full reset.
     */
    std::vector<size_t> structuralEntryIndices;
    const bool cadBatchChanged = !reset &&
	assembly->compactPresentationCadBatchRevision !=
	    this->cadBatchRevisionGet();
    SbBool sparseCadBatchUpdate = FALSE;
    if (cadBatchChanged)
	sparseCadBatchUpdate = this->getCadBatchChangedEntries(
	    assembly->compactPresentationCadBatchRevision,
	    structuralEntryIndices);
    if (cadBatchChanged && !sparseCadBatchUpdate) {
	structuralEntryIndices.reserve(
	    this->d->compactIndex->entries.size() >
		assembly->compactInstancePresentations.size() ?
	    this->d->compactIndex->entries.size() -
		assembly->compactInstancePresentations.size() : 0);
	size_t retainedPresentations = 0;
	for (size_t i = 0; i < this->d->compactIndex->entries.size(); i++) {
	    const BObolCompactInstanceEntry &entry =
		this->d->compactIndex->entries[i];
	    const auto found =
		assembly->compactInstancePresentations.find(entry.instance);
	    if (found == assembly->compactInstancePresentations.end()) {
		structuralEntryIndices.push_back(i);
	    } else {
		retainedPresentations++;
		const SoBRLCadAssembly::CompactInstancePresentation &prior =
		    found->second;
		if (prior.geometryRevision != entry.geometryRevision ||
		    prior.activePart != entry.part ||
		    prior.appearanceRevision != entry.appearanceRevision ||
		    prior.placementRevision != entry.placementRevision ||
		    prior.visibilityRevision != entry.visibilityRevision ||
		    prior.selectionRevision != entry.selectionRevision)
		    structuralEntryIndices.push_back(i);
	    }
	}
	if (retainedPresentations !=
		assembly->compactInstancePresentations.size()) {
	    reset = true;
	    structuralEntryIndices.clear();
	}
    } else if (sparseCadBatchUpdate) {
	/*
	 * The source mutation supplied an exact, non-consuming occurrence
	 * journal.  Validate only those records; other occurrences are known to
	 * match the presentation acknowledged at the prior CAD-batch revision.
	 */
	for (auto it = structuralEntryIndices.begin();
	     it != structuralEntryIndices.end();) {
	    if (*it >= this->d->compactIndex->entries.size()) {
		reset = true;
		structuralEntryIndices.clear();
		break;
	    }
	    const BObolCompactInstanceEntry &entry =
		this->d->compactIndex->entries[*it];
	    const auto prior =
		assembly->compactInstancePresentations.find(entry.instance);
	    if (prior == assembly->compactInstancePresentations.end()) {
		++it;
		continue;
	    }
	    const SoBRLCadAssembly::CompactInstancePresentation &presented =
		prior->second;
	    if (presented.geometryRevision == entry.geometryRevision &&
		presented.activePart == entry.part &&
		presented.appearanceRevision == entry.appearanceRevision &&
		presented.placementRevision == entry.placementRevision &&
		presented.visibilityRevision == entry.visibilityRevision &&
		presented.selectionRevision == entry.selectionRevision)
		it = structuralEntryIndices.erase(it);
	    else
		++it;
	}
    }

    /* A camera epoch usually changes only the active PoP cut on a bounded
     * wave of occurrences.  Preserve the compiled assembly and route those
     * changes directly by occurrence key.  Source structure/style revisions
     * still take the authoritative full-scan path. */
    const bool payloadOnlyUpdate = !reset &&
	assembly->compactPresentationCadBatchRevision ==
	    this->cadBatchRevisionGet();
    const bool incrementalStructureUpdate =
	!reset && !payloadOnlyUpdate;
    std::vector<SbString> changedOccurrenceKeys;
    SbBool fullPayloadResync = TRUE;
    if (!reset)
	viewState->cadOccurrenceChangesSince(
	    this, assembly->compactPresentationPayloadRevision,
	    changedOccurrenceKeys, fullPayloadResync);
    const bool sparsePayloadUpdate = !reset && !fullPayloadResync;
    if (sparsePayloadUpdate && changedOccurrenceKeys.empty()) {
	if (incrementalStructureUpdate && !structuralEntryIndices.empty()) {
	    /* Continue below with just the compact-index delta. */
	} else {
	    if (incrementalStructureUpdate)
		assembly->compactPresentationCadBatchRevision =
		    this->cadBatchRevisionGet();
	    assembly->compactPresentationPayloadRevision = payloadRevision;
	    viewState->acknowledgeCadOccurrenceChanges(this, payloadRevision);
	    return assembly;
	}
    }

    std::vector<const BObolViewLodState::CadPayload *>
	authoritativePayloads;
    const std::vector<const BObolViewLodState::CadPayload *> *
	payloadInput = &payloads;
    if (!sparsePayloadUpdate && payloads.empty()) {
	viewState->findCadPayloadsUnordered(this, authoritativePayloads);
	payloadInput = &authoritativePayloads;
    }

    /* Both collections may contain thousands of records.  Route modern
     * occurrence-specific payloads by key instead of scanning every payload
     * for every compact leaf.  Retain the small path scan only for legacy
     * source-wide payloads with no occurrence identity. */
    std::unordered_map<std::string_view,
	const BObolViewLodState::CadPayload *> payloadByInstance;
    std::vector<const BObolViewLodState::CadPayload *> sourceWidePayloads;
    payloadByInstance.reserve(payloadInput->size());
    for (const BObolViewLodState::CadPayload *payload : *payloadInput) {
	/*
	 * Payload input comes from the validated view-state occurrence table or
	 * from the caller's result lookup.  It cannot become invalid in place:
	 * replacement and eviction publish/erase whole shared payload objects.
	 */
	if (!payload)
	    continue;
	if (payload->sourceInstanceKey.getLength() > 0)
	    payloadByInstance[std::string_view(
		payload->sourceInstanceKey.getString(),
		static_cast<size_t>(
		    payload->sourceInstanceKey.getLength()))] = payload;
	else
	    sourceWidePayloads.push_back(payload);
    }

    std::vector<size_t> entryIndices;
    if (sparsePayloadUpdate || incrementalStructureUpdate) {
	entryIndices = structuralEntryIndices;
	entryIndices.reserve(entryIndices.size() +
	    changedOccurrenceKeys.size());
	for (const SbString &occurrenceKey : changedOccurrenceKeys) {
	    const auto found =
		this->d->compactIndex->entryIndexByKey.find(
		    occurrenceKey.getString());
	    if (found != this->d->compactIndex->entryIndexByKey.end())
		entryIndices.push_back(found->second);
	}
	std::sort(entryIndices.begin(), entryIndices.end());
	entryIndices.erase(
	    std::unique(entryIndices.begin(), entryIndices.end()),
	entryIndices.end());
    } else {
	entryIndices.resize(this->d->compactIndex->entries.size());
	std::iota(entryIndices.begin(), entryIndices.end(), 0);
    }
    const int64_t presentationIndexed =
	presentationTiming ? bu_gettime() : 0;

    BObolCadPresentationMutation mutation;
    mutation.reset = reset;
    BObolCompactPresentationStaging presentationStaging(assembly, reset);
    auto &baseParts = mutation.baseParts;
    auto &lodSharedParts = mutation.lodParts;
    auto &changedInstanceIndices = mutation.changedInstanceIndices;
    auto &instanceStyles = mutation.styles;
    auto &lodCutUpdates = mutation.cuts;
    auto &layerInstanceUpdates = mutation.layerInstances;
    auto &layerInstancesToRemove = mutation.removedLayerInstances;
    auto &layerSemantics = mutation.layerSemantics;
    auto &layerSemanticsToRemove = mutation.removedLayerSemantics;
    std::unordered_set<Obol::PartId, std::hash<Obol::PartId>>
	sharedLodPartsUpdated;
    auto &partsToRemove = mutation.removedParts;
    bool &instanceSetsChanged = mutation.instanceSetsChanged;
    instanceSetsChanged = reset;
    std::vector<std::pair<Obol::PartId, uint8_t>> basePartChannels;

    if (reset) {
	baseParts.reserve(this->d->compactIndex->parts.size());
	basePartChannels.reserve(this->d->compactIndex->parts.size());
	for (const BObolCompactPartReference &partRef :
	     this->d->compactIndex->parts) {
	    if (!partRef.geometry)
		continue;
	    Obol::PartUpdate part;
	    part.part = partRef.part;
	    if (!bobol_cad_admit_geometry(partRef.geometry, part.geometry,
		    "compact base-part staging"))
		return assembly;
	    basePartChannels.emplace_back(part.part,
		cad_compact_geometry_channels(part.geometry.get()));
	    baseParts.push_back(std::move(part));
	}
	for (const auto &partChannel : basePartChannels)
	    presentationStaging.setPartChannel(
		partChannel.first, partChannel.second);
    } else if (!structuralEntryIndices.empty()) {
	std::unordered_set<Obol::PartId, std::hash<Obol::PartId>>
	    addedBaseParts;
	addedBaseParts.reserve(structuralEntryIndices.size());
	for (const size_t i : structuralEntryIndices) {
	    if (i >= this->d->compactIndex->entries.size())
		continue;
	    const BObolCompactInstanceEntry &entry =
		this->d->compactIndex->entries[i];
	    const auto prior =
		assembly->compactInstancePresentations.find(entry.instance);
	    /* A producer may keep the numeric geometry revision while changing
	     * the immutable base part (for example, when a temporary LoD box
	     * resolves to an authored line with byte-identical coordinates).  That
	     * replacement still needs the new base part registered before its
	     * instance record can be republished. */
	    if (prior != assembly->compactInstancePresentations.end() &&
		prior->second.geometryRevision == entry.geometryRevision &&
		prior->second.activePart == entry.part)
		continue;
	    if (!entry.geometry || !addedBaseParts.insert(entry.part).second)
		continue;
	    Obol::PartUpdate part;
	    part.part = entry.part;
	    if (!bobol_cad_admit_geometry(entry.geometry, part.geometry,
		    "incremental compact base-part staging"))
		return assembly;
	    baseParts.push_back(std::move(part));
	    basePartChannels.emplace_back(entry.part,
		cad_compact_geometry_channels(entry.geometry.get()));
	}
	for (const auto &partChannel : basePartChannels)
	    presentationStaging.setPartChannel(
		partChannel.first, partChannel.second);
    }

    const bool incrementalUpdate = !reset;
    size_t wireCount = incrementalUpdate ?
	assembly->compactWirePresentationCount : 0;
    size_t shadedCount = incrementalUpdate ?
	assembly->compactShadedPresentationCount : 0;
    for (const size_t i : entryIndices) {
	const BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[i];
	const SbString &occurrenceKey = compact_instance_identity(entry);
	const BObolViewLodState::CadPayload *payload = NULL;
	if (sparsePayloadUpdate) {
	    payload = viewState->findCadForOccurrence(
		this, occurrenceKey);
	} else {
	    const auto payloadIt = payloadByInstance.find(
		std::string_view(occurrenceKey.getString(),
		    static_cast<size_t>(occurrenceKey.getLength())));
	    payload =
		payloadIt != payloadByInstance.end() ? payloadIt->second :
		cad_compact_payload_for_entry(sourceWidePayloads, entry);
	}
	/* A snapshot can only be created by successful Obol admission.  Reject
	 * missing payload storage before changing retained presentation state. */
	if (payload) {
	    for (const BObolLodPresentationLayer &layer :
		 payload->presentationLayers) {
		if (!layer.geometry) {
		    payload = NULL;
		    break;
		}
	    }
	}
	const bool spatiallyCulled = payload && payload->progressiveMesh &&
	    payload->progressiveMesh->hasSpatialClusters() &&
	    payload->requiredChunks.empty();
	const bool wasSpatiallyCulled =
	    presentationStaging.isSpatiallyCulled(entry.instance);
	presentationStaging.setSpatiallyCulled(
	    entry.instance, spatiallyCulled);
	if (spatiallyCulled != wasSpatiallyCulled)
	    instanceSetsChanged = true;

	const auto existingPresentation =
	    presentationStaging.findPresentation(entry.instance);
	const bool presentationWasPresent =
	    existingPresentation != NULL;
	SoBRLCadAssembly::CompactInstancePresentation &presentation =
	    presentationStaging.editPresentation(entry.instance);
	const uint8_t previousChannels = presentation.channels;
	const Obol::PartId previousPart = presentation.activePart;
	if (reset) {
	    presentation.activePart = entry.part;
	    presentation.channels = presentationStaging.partChannel(entry.part);
	    presentation.activeCut = -1;
	}

	std::string payloadKey;
	Obol::PartId desiredPart = entry.part;
	uint8_t desiredChannels = presentationStaging.partChannel(entry.part);
	int desiredActiveCut = viewState->residentCadProgressiveCut(
	    this, i <= static_cast<size_t>(UINT32_MAX) ?
		static_cast<uint32_t>(i) : UINT32_MAX, occurrenceKey,
	    entry.geometryRevision);
	if (desiredActiveCut < 0)
	    desiredActiveCut = cad_compact_default_progressive_cut(entry);
	/* A direct progressive cut reuses the authored immutable part.  It is
	 * view state, but not a cache/service payload and therefore must not
	 * manufacture a second part identity or geometry publication. */
	bool desiredGeometryValid = desiredActiveCut >= 0;
	bool desiredLodStructuralProxy = entry.geometry &&
	    entry.geometry->structuralProxy && entry.lodBacked;
	if (payload) {
	    payloadKey.reserve(
		static_cast<size_t>(payload->cacheKey.getLength()) + 128u);
	    payloadKey.assign(payload->cacheKey.getString(),
		static_cast<size_t>(payload->cacheKey.getLength()));
	    payloadKey += ':';
	    payloadKey += std::to_string(payload->resultKind);
	    payloadKey += ':';
	    payloadKey += std::to_string(payload->qualityTier);
	    payloadKey += ':';
	    if (payload->progressiveMesh) {
		payloadKey += "resident-revision=";
		payloadKey += std::to_string(
		    payload->progressiveMesh->revision());
		desiredActiveCut = payload->activeCut;
	    } else {
		payloadKey += "cut=";
		payloadKey += std::to_string(payload->activeCut);
	    }
	    payloadKey += ':';
	    payloadKey += std::to_string(sourceDrawMode);
	    payloadKey += ":normal-style=";
	    payloadKey += std::to_string(static_cast<int>(normalStyle));
	    payloadKey += ":normal-crease=";
	    payloadKey += std::to_string(normalCreaseAngle);
	    const bool layeredPresentation =
		!payload->presentationLayers.empty();
	    if (layeredPresentation)
		cad_compact_append_layer_identity(
		    payloadKey, payload->presentationLayers);
	    const Obol::PartId expectedPayloadPart =
		cad_compact_payload_part_id(
		    *payload, occurrenceKey, payloadKey, sourceDrawMode);
	    /* A payload key describes immutable data, not the retained instance
	     * binding.  A failed or interrupted sparse rebind may have advanced the
	     * key while leaving activePart on the source box.  Trust the fast path
	     * only when both pieces of the presentation transaction agree. */
	    if (payloadKey == presentation.payloadKey &&
		presentation.activePart == expectedPayloadPart) {
		desiredPart = presentation.activePart;
		desiredChannels = presentation.channels;
		desiredGeometryValid = true;
	    } else {
		Obol::PartGeometryBuilder geometry;
		const bool progressive =
		    payload->progressiveMesh &&
		    payload->progressiveMesh->isValid();
		const bool sharedTerminal = !progressive &&
		    payload->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL &&
		    payload->preparedCadGeometry;
		const bool sharedLayer = layeredPresentation;
		bool reusedSharedPart = false;
		desiredPart = expectedPayloadPart;
		if (sharedLayer) {
		    const BObolLodPresentationLayer &baseLayer =
			payload->presentationLayers.front();
		    desiredActiveCut = baseLayer.activeCut;
		}
		const int wire =
		    sourceDrawMode == BOBOL_LOD_DRAW_WIRE || hiddenLine;
		/* SHADED_BOTS creates source-mesh requests only for BoTs.
		 * Those payloads are shaded just like mode 2 payloads. */
		const int shaded =
		    sourceDrawMode == BOBOL_LOD_DRAW_SHADED_BOTS ||
		    sourceDrawMode == BOBOL_LOD_DRAW_SHADED ||
		    hiddenLine;
		std::shared_ptr<const Obol::PartGeometry> preparedGeometry;
		if (sharedLayer) {
		    preparedGeometry =
			payload->presentationLayers.front().geometry;
		} else if (!progressive && payload->preparedCadGeometry) {
		    preparedGeometry = payload->preparedCadGeometry;
		} else if (progressive &&
		    normalStyle == BObolViewLodState::NORMAL_AUTHORED &&
		    payload->preparedCadGeometry &&
		    payload->preparedCadGeometryRevision ==
			payload->progressiveMesh->revision() &&
		    cad_progressive_geometry_can_draw(
			payload->preparedCadGeometry.get(), wire, shaded,
			desiredActiveCut)) {
		    preparedGeometry = payload->preparedCadGeometry;
		}

		/*
		 * A result without a newly prepared renderer generation need not
		 * force a full CPU/GPU realization while the view is asking for a
		 * coarser prefix: keep the older cumulative part when it covers the
		 * active cut.  A different prepared geometry pointer is an
		 * authoritative immutable-generation handoff, however.  Its ordinal
		 * cut may match the old part while its spatial page population does
		 * not; retaining the old part in that case leaves visible holes until
		 * an unrelated richer cut happens to replace it.
		 */
		const Obol::PartGeometry *presentedGeometry =
		    assembly->partGeometry(desiredPart);
		const bool preparedGenerationChanged = preparedGeometry &&
		    preparedGeometry.get() != presentedGeometry;
		if (progressive && !preparedGenerationChanged &&
		    presentation.activePart == desiredPart &&
		    cad_progressive_geometry_can_draw(
			presentedGeometry, wire, shaded,
			desiredActiveCut)) {
		    desiredGeometryValid = true;
		    reusedSharedPart = true;
		    desiredChannels =
			presentationStaging.partChannel(desiredPart);
		    payloadKey = presentation.payloadKey;
		} else if ((progressive || sharedTerminal || sharedLayer) &&
		    sharedLodPartsUpdated.find(desiredPart) !=
			sharedLodPartsUpdated.end()) {
		    desiredGeometryValid = true;
		    reusedSharedPart = true;
		    desiredChannels =
			presentationStaging.partChannel(desiredPart);
		} else if (payload->resultKind == BOBOL_LOD_RESULT_AABB ||
		    payload->resultKind == BOBOL_LOD_RESULT_PROXY) {
		    desiredGeometryValid =
			cad_proxy_part_geometry(
			    payload->proxy, wire, shaded, geometry) != 0;
		} else if (preparedGeometry) {
		    desiredGeometryValid = true;
		} else {
		    desiredGeometryValid = progressive ?
			cad_progressive_mesh_part_geometry(
			    payload->progressiveMesh, wire, shaded,
			    payload->shadedCullBackfaces, geometry,
			    normalStyle, normalCreaseAngle) != 0 :
			cad_mesh_payload_part_geometry(payload->mesh,
			    wire, shaded, payload->shadedCullBackfaces,
			    geometry, normalStyle, normalCreaseAngle) != 0;
		}
		if (desiredGeometryValid && !reusedSharedPart) {
		    desiredChannels = preparedGeometry ?
			cad_compact_geometry_channels(preparedGeometry.get()) :
			cad_compact_geometry_channels(&geometry);
		    const bool uniqueSharedPart =
			progressive || sharedTerminal || sharedLayer;
		    const bool newlyRegistered = uniqueSharedPart &&
			sharedLodPartsUpdated.insert(desiredPart).second;
		    if (!uniqueSharedPart || newlyRegistered) {
			bool admitted = false;
			if (preparedGeometry) {
			    Obol::PartUpdate part;
			    part.part = desiredPart;
			    admitted = bobol_cad_admit_geometry(preparedGeometry,
				part.geometry,
				"prepared compact LoD staging");
			    /*
			     * Every immutable PoP generation carries the same
			     * full quantization-domain bounds.  Let Obol validate
			     * and retain occurrence bounds/BVH state while merely
			     * swapping the renderer geometry handle.
			     */
			    part.preservesBounds = true;
			    if (admitted)
				lodSharedParts.push_back(std::move(part));
			} else {
			    /* Publish completed private geometry by move so proxy
			     * and non-progressive fallback publication stay zero-copy
			     * at the presentation boundary. */
			    Obol::PartUpdate part;
			    part.part = desiredPart;
			    const auto sharedGeometry =
				bobol_cad_build_geometry(std::move(geometry),
				    "generated compact LoD staging");
			    admitted = bobol_cad_admit_geometry(sharedGeometry,
				part.geometry,
				"generated compact LoD staging");
			    if (admitted)
				lodSharedParts.push_back(std::move(part));
			}
			if (!admitted) {
			    desiredGeometryValid = false;
			    if (newlyRegistered)
				sharedLodPartsUpdated.erase(desiredPart);
			}
		    }
		    if (desiredGeometryValid) {
			presentationStaging.setPartChannel(
			    desiredPart, desiredChannels);
			presentationStaging.addLodPart(desiredPart);
		    }
		}
	    }
	}
	if (!desiredGeometryValid)
	{
	    /* A rejected or stale LoD payload must fall back to the authored
	     * structural part as one coherent record.  Keeping the hashed LoD
	     * part id after geometry conversion failed leaves the instance
	     * referencing a part that was never inserted and makes it vanish. */
	    payloadKey.clear();
	    desiredPart = entry.part;
	    desiredChannels = presentationStaging.partChannel(entry.part);
	    desiredActiveCut = -1;
	    desiredLodStructuralProxy = entry.geometry &&
		entry.geometry->structuralProxy && entry.lodBacked;
	}
	else if (payload) {
	    desiredLodStructuralProxy =
		payload->resultKind == BOBOL_LOD_RESULT_AABB ||
		(payload->resultKind == BOBOL_LOD_RESULT_PROXY &&
		 payload->proxy.kind == BOBOL_LOD_PROXY_AABB);
	}

	/* Structural proxies and source meshes generally do not share a local
	 * coordinate system.  Cold coverage may use one normalized unit box
	 * transformed by geometryTransform, while the cached PoP arrays remain
	 * in the BoT asset's original coordinates.  Swapping only the part id
	 * therefore draws a valid mesh with the proxy transform—often far
	 * outside the view or collapsed to an apparently empty wireframe.
	 * Choose the transform with the chosen representation and update it in
	 * the same instance mutation as the part. */
	const bool payloadUsesSourceMeshCoordinates = payload &&
	    (payload->resultKind == BOBOL_LOD_RESULT_MESH ||
	     payload->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL) &&
	    entry.sourceMeshRequestValid;
	const SbMatrix desiredLocalToRoot =
	    payloadUsesSourceMeshCoordinates ?
		compact_mesh_asset_matrix(this, entry) : entry.localToSource;

	/* A live cold result may contain several independently drawable renderer
	 * layers for this one logical CAD occurrence.  The first layer occupies the
	 * authored instance record; later layers use private deterministic instance
	 * IDs which inherit the same transform, style, and semantic identity. */
	std::vector<SoBRLCadAssembly::CompactLayerPresentation> desiredLayers;
	const std::vector<SoBRLCadAssembly::CompactLayerPresentation> *
	    previousLayers = presentationStaging.findLayers(entry.instance);
	std::unordered_map<std::string_view,
	    const SoBRLCadAssembly::CompactLayerPresentation *> previousByKey;
	if (previousLayers) {
	    previousByKey.reserve(previousLayers->size());
	    for (const SoBRLCadAssembly::CompactLayerPresentation &layer :
		 *previousLayers)
		previousByKey.emplace(layer.layerKey, &layer);
	}
	if (payload && desiredGeometryValid &&
	    !payload->presentationLayers.empty()) {
	    desiredLayers.reserve(payload->presentationLayers.size());
	    for (size_t layerIndex = 0;
		 layerIndex < payload->presentationLayers.size(); ++layerIndex) {
		const BObolLodPresentationLayer &sourceLayer =
		    payload->presentationLayers[layerIndex];
		SoBRLCadAssembly::CompactLayerPresentation layer;
		layer.layerKey = sourceLayer.layerKey.getString();
		layer.part = cad_compact_layer_part_id(
		    *payload, sourceLayer, sourceDrawMode);
		layer.instance = layerIndex == 0 ? entry.instance :
		    cad_compact_layer_instance_id(entry.instance, sourceLayer);
		layer.channels = cad_compact_geometry_channels(
		    sourceLayer.geometry.get());
		layer.activeCut = sourceLayer.activeCut;
		layer.coverage = sourceLayer.coverage ? true : false;
		layer.geometryRevision = sourceLayer.geometryRevision;

		const auto prior = previousByKey.find(layer.layerKey);
		const SoBRLCadAssembly::CompactLayerPresentation *oldLayer =
		    prior == previousByKey.end() ? NULL : prior->second;
		const bool geometryChanged = !oldLayer ||
		    oldLayer->part != layer.part ||
		    oldLayer->geometryRevision != layer.geometryRevision ||
		    assembly->partGeometry(layer.part) !=
			sourceLayer.geometry.get();
		if (geometryChanged &&
		    sharedLodPartsUpdated.insert(layer.part).second) {
		    Obol::PartUpdate part;
		    part.part = layer.part;
		    if (!bobol_cad_admit_geometry(sourceLayer.geometry,
			    part.geometry, "spatial LoD layer staging")) {
			sharedLodPartsUpdated.erase(layer.part);
			continue;
		    }
		    part.preservesBounds = true;
		    lodSharedParts.push_back(std::move(part));
		}
		presentationStaging.setPartChannel(layer.part, layer.channels);
		presentationStaging.addLodPart(layer.part);

		if (layerIndex > 0) {
		    /* Auxiliary presentation bookkeeping describes the desired logical
		     * layer, not proof that its physical Obol instance survived a preceding
		     * reset, cull, or rejected publication.  A sparse cut/style update may
		     * target only a live record; recover a missing layer with the complete
		     * idempotent upsert already used for a changed record. */
		    const bool retainedInstanceMissing =
			!assembly->getInstanceRecord(layer.instance).has_value();
		    const bool recordChanged = !oldLayer ||
			retainedInstanceMissing || oldLayer->part != layer.part ||
			presentation.placementRevision != entry.placementRevision;
		    if (recordChanged) {
			Obol::InstanceUpdate update;
			update.instance = layer.instance;
			update.record = this->d->compactIndex->instances[i].record;
			update.record.part = layer.part;
			update.record.localToRoot = desiredLocalToRoot;
			update.record.style = entry.style;
			update.record.lodStructuralProxy = false;
			update.record.lodCut = layer.activeCut >= 0 ?
			    static_cast<uint8_t>(std::min<int>(
				Obol::ProgressiveCutLimit - 1,
				layer.activeCut)) : 255;
			layerInstanceUpdates.push_back(std::move(update));
		    } else if (presentation.appearanceRevision !=
			    entry.appearanceRevision ||
			presentation.selectionRevision !=
			    entry.selectionRevision) {
			Obol::InstanceStyleUpdate update;
			update.instance = layer.instance;
			update.style = entry.style;
			instanceStyles.push_back(std::move(update));
		    } else if (oldLayer->activeCut != layer.activeCut) {
			Obol::InstanceLodUpdate update;
			update.instance = layer.instance;
			update.lodCut = layer.activeCut >= 0 ?
			    static_cast<uint8_t>(std::min<int>(
				Obol::ProgressiveCutLimit - 1,
				layer.activeCut)) : 255;
			lodCutUpdates.push_back(update);
		    }
		    if (!oldLayer || incrementalStructureUpdate) {
			SoBRLCadAssembly::InstanceSemantic semantic =
			    entry.semantic;
			semantic.sourceInstanceKey = occurrenceKey;
			layerSemantics.push_back(
			    std::make_pair(layer.instance, semantic));
		    }
		}
		desiredLayers.push_back(std::move(layer));
	    }
	}

	/* Recompute only the auxiliary physical-instance contributions for this
	 * logical occurrence.  The base contribution remains handled by the
	 * ordinary compact presentation accounting below. */
	if (previousLayers) {
	    for (size_t layerIndex = 1;
		 layerIndex < previousLayers->size(); ++layerIndex) {
		const SoBRLCadAssembly::CompactLayerPresentation &oldLayer =
		    (*previousLayers)[layerIndex];
		if (incrementalUpdate)
		    cad_compact_adjust_channel_count(
			oldLayer.channels, false, wireCount, shadedCount);
		presentationStaging.removePartReference(oldLayer.part);
		partsToRemove.insert(oldLayer.part);
		const auto replacement = std::find_if(
		    desiredLayers.begin(), desiredLayers.end(),
		    [&oldLayer](const auto &candidate) {
			return candidate.instance == oldLayer.instance;
		    });
		if (replacement == desiredLayers.end()) {
		    layerInstancesToRemove.push_back(oldLayer.instance);
		    layerSemanticsToRemove.push_back(oldLayer.instance);
		}
	    }
	}
	for (size_t layerIndex = 1; layerIndex < desiredLayers.size();
	     ++layerIndex) {
	    const SoBRLCadAssembly::CompactLayerPresentation &layer =
		desiredLayers[layerIndex];
	    presentationStaging.addPartReference(layer.part);
	    cad_compact_adjust_channel_count(
		layer.channels, true, wireCount, shadedCount);
	}
	const bool desiredLayerPresentation = !desiredLayers.empty();
	presentationStaging.setLayers(entry.instance, std::move(desiredLayers));
	if ((previousLayers != NULL) !=
		desiredLayerPresentation ||
	    !layerInstanceUpdates.empty() || !layerInstancesToRemove.empty())
	    instanceSetsChanged = true;

	const bool partChanged = presentation.activePart != desiredPart ||
	    presentation.payloadKey != payloadKey;
	const bool cutChanged =
	    presentation.activeCut != desiredActiveCut;
	const bool appearanceChanged = !reset &&
	    presentation.appearanceRevision != entry.appearanceRevision;
	const bool selectionChanged = !reset &&
	    presentation.selectionRevision != entry.selectionRevision;
	/* Replacing immutable arrays behind the same part is wholly described by
	 * the shared-part geometry journal.  Re-publishing its unchanged instance
	 * record would mix a same-part update into an otherwise append-only leaf
	 * batch and force the complete growing frame plan to be rebuilt. */
	const bool overviewGeometryOnly = !partChanged &&
	    BU_STR_EQUAL(entry.shapeSummary.recordRole.getString(),
		"lod-overview");
	const bool retainedInstanceMissing = !reset &&
	    !assembly->getInstanceRecord(entry.instance).has_value();
	const bool recordChanged = !reset &&
	    (retainedInstanceMissing || partChanged ||
	     (!overviewGeometryOnly &&
		presentation.geometryRevision != entry.geometryRevision) ||
	     presentation.placementRevision != entry.placementRevision ||
	     presentation.lodStructuralProxy !=
		desiredLodStructuralProxy);
	if (partChanged) {
	    if (presentation.activePart.isValid() &&
		presentation.activePart != entry.part &&
		presentation.activePart != desiredPart)
		partsToRemove.insert(presentation.activePart);
	    presentation.payloadKey = payloadKey;
	    presentation.activePart = desiredPart;
	    presentation.channels = desiredChannels;
	}
	if (!presentationWasPresent || partChanged) {
	    if (presentationWasPresent)
		presentationStaging.removePartReference(previousPart);
	    presentationStaging.addPartReference(desiredPart);
	}
	if (reset && i < this->d->compactIndex->instances.size()) {
	    Obol::InstanceUpdate &update =
		this->d->compactIndex->instances[i];
	    update.record.part = desiredPart;
	    update.record.localToRoot = desiredLocalToRoot;
	    update.record.style = entry.style;
	    update.record.lodStructuralProxy =
		desiredLodStructuralProxy;
	    update.record.lodCut = desiredActiveCut >= 0 ?
		static_cast<uint8_t>(std::min<int>(
		    Obol::ProgressiveCutLimit - 1, desiredActiveCut)) : 255;
	} else if (recordChanged &&
	    i < this->d->compactIndex->instances.size()) {
	    Obol::InstanceUpdate &update =
		this->d->compactIndex->instances[i];
	    update.record.part = desiredPart;
	    update.record.localToRoot = desiredLocalToRoot;
	    /* The retained index record is structural storage and may predate a
	     * selection/highlight delta.  A geometry or placement change publishes
	     * a complete record, so carry the entry's authoritative effective style
	     * in that same atomic update.  Otherwise the full record can overwrite a
	     * newer selected style and the following selection revision is marked as
	     * consumed without ever reaching the renderer. */
	    update.record.style = entry.style;
	    update.record.lodStructuralProxy =
		desiredLodStructuralProxy;
	    update.record.lodCut = desiredActiveCut >= 0 ?
		static_cast<uint8_t>(std::min<int>(
		    Obol::ProgressiveCutLimit - 1, desiredActiveCut)) : 255;
	    changedInstanceIndices.push_back(i);
	} else if (appearanceChanged || selectionChanged) {
	    Obol::InstanceStyleUpdate update;
	    update.instance = entry.instance;
	    update.style = entry.style;
	    instanceStyles.push_back(update);
	} else if (cutChanged) {
	    Obol::InstanceLodUpdate update;
	    update.instance = entry.instance;
	    update.lodCut = desiredActiveCut >= 0 ?
		static_cast<uint8_t>(std::min<int>(
		    Obol::ProgressiveCutLimit - 1, desiredActiveCut)) : 255;
	    lodCutUpdates.push_back(update);
	}
	presentation.activeCut = desiredActiveCut;
	presentation.lodStructuralProxy = desiredLodStructuralProxy;

	if (incrementalUpdate) {
	    const bool previousWire = previousChannels & (1u | 4u);
	    const bool currentWire = presentation.channels & (1u | 4u);
	    const bool previousShaded = previousChannels & 2u;
	    const bool currentShaded = presentation.channels & 2u;
	    if (previousWire != currentWire)
		wireCount = currentWire ? wireCount + 1 :
		    (wireCount > 0 ? wireCount - 1 : 0);
	    if (previousShaded != currentShaded)
		shadedCount = currentShaded ? shadedCount + 1 :
		    (shadedCount > 0 ? shadedCount - 1 : 0);
	    } else {
	    if (presentation.channels & (1u | 4u))
		wireCount++;
	    if (presentation.channels & 2u)
		shadedCount++;
	}

	if (!reset &&
	    (presentation.visibilityRevision != entry.visibilityRevision ||
	     presentation.selectionRevision != entry.selectionRevision)) {
	    instanceSetsChanged = true;
	}
	presentation.geometryRevision = entry.geometryRevision;
	presentation.appearanceRevision = entry.appearanceRevision;
	presentation.placementRevision = entry.placementRevision;
	presentation.visibilityRevision = entry.visibilityRevision;
	presentation.selectionRevision = entry.selectionRevision;
    }
    const int64_t presentationPlanned =
	presentationTiming ? bu_gettime() : 0;

    /* Progressive assets are shared by every occurrence of the same source
     * geometry.  One occurrence changing back to its base part must not
     * remove a retained part that another occurrence still references. */
    if (!partsToRemove.empty()) {
	for (auto part = partsToRemove.begin();
	     part != partsToRemove.end();) {
	    if (presentationStaging.activePartReferences(*part) > 0)
		part = partsToRemove.erase(part);
	    else
		++part;
	}
    }

    const Obol::CadDrawMode presentationDrawMode =
	cad_presentation_draw_mode(hiddenLine ?
	    BOBOL_LOD_DRAW_HIDDEN_LINE : BOBOL_LOD_DRAW_UNKNOWN,
	    shadedCount, wireCount);
    const bool drawModeChanged =
	assembly->presentationDrawMode() != presentationDrawMode;
    mutation.drawModeChanged = drawModeChanged;
    const bool publicationRequired = mutation.publicationRequired();

    const char *presentationDebug =
	getenv("BOBOL_COMPACT_PRESENTATION_DEBUG");
    const bool presentationDebugEnabled = presentationDebug &&
	presentationDebug[0] && !BU_STR_EQUAL(presentationDebug, "0");
    std::vector<Obol::PartId> removedParts(
	partsToRemove.begin(), partsToRemove.end());
    std::sort(removedParts.begin(), removedParts.end());
    for (const Obol::PartId part : removedParts)
	presentationStaging.removePart(part);

    Obol::CadSceneMutation geometryMutation;
    if (!reset) {
	geometryMutation.parts.reserve(baseParts.size() + lodSharedParts.size());
	geometryMutation.parts.insert(geometryMutation.parts.end(),
	    baseParts.begin(), baseParts.end());
	geometryMutation.parts.insert(geometryMutation.parts.end(),
	    lodSharedParts.begin(), lodSharedParts.end());
    }
    Obol::CadSceneMutation layerMutation;
    layerMutation.instances = std::move(layerInstanceUpdates);
    layerMutation.styles = std::move(instanceStyles);
    layerMutation.cuts = std::move(lodCutUpdates);
    layerMutation.removedInstances = std::move(layerInstancesToRemove);
    layerMutation.removedParts = std::move(removedParts);

    if (publicationRequired) {
	static constexpr size_t CadInstancePublicationBatchSize = 512;
	const std::vector<Obol::InstanceUpdate> &compactInstances =
	    this->d->compactIndex->instances;
	if (!bobol_cad_validate_mutation(assembly, geometryMutation,
		"compact geometry preflight"))
	    return assembly;
	if (reset) {
	    if (!bobol_cad_validate_instances(compactInstances,
		    "compact reset preflight"))
		return assembly;
	} else {
	    for (const size_t instanceIndex : changedInstanceIndices) {
		if (instanceIndex >= compactInstances.size()) {
		    bu_log("libBObol: rejected compact scene preflight: "
			"instance index %zu exceeds population %zu\n",
			instanceIndex, compactInstances.size());
		    return assembly;
		}
		if (!bobol_cad_validate_instance(
			compactInstances[instanceIndex], instanceIndex,
			"compact incremental preflight")) {
		    if (presentationDebugEnabled &&
			instanceIndex < this->d->compactIndex->entries.size()) {
			const BObolCompactInstanceEntry &entry =
			    this->d->compactIndex->entries[instanceIndex];
			const SbMatrix &matrix = compactInstances[instanceIndex].
			    record.localToRoot;
			bu_log("libBObol: invalid compact transform path=%s "
			       "source=%s asset=%s\n",
			    entry.semantic.path.getString(),
			    entry.semantic.sourceName.getString(),
			    entry.sourceMeshRequest.meshAssetName.getString());
			for (int row = 0; row < 4; ++row)
			    bu_log("libBObol: transform[%d] %.9g %.9g %.9g %.9g\n",
				row, matrix[row][0], matrix[row][1],
				matrix[row][2], matrix[row][3]);
		    }
		    return assembly;
		}
	    }
	}
	/* Complete auxiliary upserts are staged whenever a remembered physical
	 * layer is absent.  Consequently every sparse style/cut target either
	 * exists in the current assembly or is introduced by layerMutation, so the
	 * whole journal can be rejected before any live state changes. */
	if (!bobol_cad_validate_mutation(assembly, layerMutation,
		"compact layer preflight"))
	    return assembly;
	/* Everything above is private staging.  A resource denial must not open
	 * the outer update scope: even an empty scope has a notification edge. */
	if (bobol_transaction_fault_requested(
		BObolTransactionFaultPoint::RETAINED_SCENE_COMMIT) ||
	    bobol_transaction_fault_requested(
		BObolTransactionFaultPoint::PRESENTATION_COMMIT))
	    return assembly;

	SoCADAssembly::UpdateScope updateScope = assembly->batchUpdate();
	if (reset) {
	    std::vector<Obol::PartUpdate> replacementParts;
	    replacementParts.reserve(baseParts.size() + lodSharedParts.size());
	    replacementParts.insert(replacementParts.end(), baseParts.begin(),
		baseParts.end());
	    replacementParts.insert(replacementParts.end(), lodSharedParts.begin(),
		lodSharedParts.end());
	    if (!bobol_cad_replace_scene(assembly, replacementParts,
		    compactInstances, "compact scene replacement"))
		return assembly;
	    assembly->clearSemanticMap();
	    for (const BObolCompactInstanceEntry &entry :
		 this->d->compactIndex->entries) {
		SoBRLCadAssembly::InstanceSemantic semantic = entry.semantic;
		semantic.sourceInstanceKey = compact_instance_identity(entry);
		assembly->setInstanceSemantic(entry.instance, semantic);
	    }
	} else if (!bobol_cad_publish_mutation(assembly, geometryMutation,
		"compact geometry publication")) {
	    return assembly;
	}
	/* Do not form a second full 150k-instance population merely to publish a
	 * sparse representation change.  The retained compact index already owns
	 * the authoritative records.  Small bounded batches keep transient memory
	 * proportional to publication work rather than scene population. */
	std::vector<Obol::InstanceUpdate> instanceUpdates;
	instanceUpdates.reserve(CadInstancePublicationBatchSize);
	bool instancePublicationValid = true;
	auto publishInstanceUpdates = [&assembly, &compactInstances, &instanceUpdates,
		&instancePublicationValid](size_t instanceIndex) {
	    if (!instancePublicationValid)
		return;
	    instanceUpdates.push_back(compactInstances[instanceIndex]);
	    if (instanceUpdates.size() == CadInstancePublicationBatchSize) {
		Obol::CadSceneMutation batch;
		batch.instances.swap(instanceUpdates);
		instancePublicationValid = bobol_cad_publish_mutation(assembly,
		    batch, "compact instance publication");
		instanceUpdates.reserve(CadInstancePublicationBatchSize);
	    }
	};
	if (!reset) {
	    for (const size_t i : changedInstanceIndices)
		publishInstanceUpdates(i);
	}
	if (instancePublicationValid && !instanceUpdates.empty()) {
	    Obol::CadSceneMutation batch;
	    batch.instances.swap(instanceUpdates);
	    instancePublicationValid = bobol_cad_publish_mutation(assembly,
		batch, "compact instance publication");
	}
	if (!instancePublicationValid)
	    return assembly;
	if (!bobol_cad_publish_mutation(assembly, layerMutation,
		"compact layer publication"))
	    return assembly;
	for (const auto &semantic : layerSemantics)
	    assembly->setInstanceSemantic(semantic.first, semantic.second);
	for (const Obol::InstanceId instance : layerSemanticsToRemove)
	    assembly->semantics.erase(instance);
	if (!reset) {
	    for (const size_t i : structuralEntryIndices) {
		if (i >= this->d->compactIndex->entries.size())
		    continue;
		const BObolCompactInstanceEntry &entry =
		    this->d->compactIndex->entries[i];
		SoBRLCadAssembly::InstanceSemantic semantic = entry.semantic;
		semantic.sourceInstanceKey = compact_instance_identity(entry);
		assembly->setInstanceSemantic(entry.instance, semantic);
	    }
	}
	std::vector<Obol::InstanceId> hiddenInstances;
	std::vector<Obol::InstanceId> selectedInstances;
	std::vector<Obol::InstanceId> unpickableInstances;
	if (instanceSetsChanged) {
	    hiddenInstances = this->d->compactIndex->hiddenInstances;
	    presentationStaging.appendSpatiallyCulled(hiddenInstances);
	    std::sort(hiddenInstances.begin(), hiddenInstances.end());
	    hiddenInstances.erase(std::unique(hiddenInstances.begin(),
		hiddenInstances.end()), hiddenInstances.end());
	    hiddenInstances =
		presentationStaging.expandInstanceSet(hiddenInstances);
	    selectedInstances = presentationStaging.expandInstanceSet(
		this->d->compactIndex->selectedInstances);
	    unpickableInstances = presentationStaging.expandInstanceSet(
		this->d->compactIndex->unpickableInstances);
	}
	presentationStaging.commit();
	if (instanceSetsChanged) {
	    assembly->setHiddenInstances(hiddenInstances);
	    assembly->setSelectedInstances(selectedInstances);
	    assembly->setUnpickableInstances(unpickableInstances);
	}
	assembly->setPresentationDrawMode(presentationDrawMode);
	assembly->setPresentationPickMode(Obol::CadPickMode::Automatic);
	updateScope.finish();
    } else {
	presentationStaging.commit();
    }

    /*
     * Compact presentation is intentionally sparse, so an intended active
     * part and the retained Obol instance must never drift silently.  Keep
     * this full invariant audit behind an opt-in diagnostic: it is useful for
     * validating streaming box-to-mesh waves, but an O(N) scan does not
     * belong in the normal publication path for very large scenes.
     */
    if (presentationDebugEnabled) {
	size_t missingInstances = 0;
	size_t mismatchedParts = 0;
	size_t retainedStructural = 0;
	size_t retainedLodStructural = 0;
	size_t structuralLodBacked = 0;
	size_t structuralSourceRequests = 0;
	size_t structuralEntryGeometry = 0;
	size_t structuralGeometryPointerMismatch = 0;
	std::map<std::string, size_t> structuralKinds;
	std::map<std::string, size_t> structuralTypes;
	for (const auto &item : assembly->compactInstancePresentations) {
	    const std::optional<Obol::InstanceRecord> record =
		assembly->getInstanceRecord(item.first);
	    if (!record) {
		missingInstances++;
		continue;
	    }
	    if (!(record->part == item.second.activePart))
		mismatchedParts++;
	    const Obol::PartGeometry *geometry =
		assembly->partGeometry(record->part);
	    if (geometry && geometry->structuralProxy) {
		retainedStructural++;
		if (record->lodStructuralProxy)
		    retainedLodStructural++;
		const auto entryIndex =
		    this->d->compactIndex->entryIndex.find(item.first);
		if (entryIndex !=
			this->d->compactIndex->entryIndex.end() &&
		    entryIndex->second <
			this->d->compactIndex->entries.size()) {
		    const BObolCompactInstanceEntry &entry =
			this->d->compactIndex->entries[entryIndex->second];
		    if (entry.geometry && entry.geometry->structuralProxy)
			structuralEntryGeometry++;
		    if (!entry.geometry ||
			entry.geometry.get() != geometry)
			structuralGeometryPointerMismatch++;
		    if (entry.lodBacked)
			structuralLodBacked++;
		    if (entry.sourceMeshRequestValid)
			structuralSourceRequests++;
		    structuralKinds[
			entry.shapeSummary.geometryKind.getString()]++;
		    structuralTypes[
			entry.shapeSummary.sourceType.getString()]++;
		}
	    }
	}
	bu_log("BObol compact presentation invariant presentations=%zu "
	       "retained_instances=%zu retained_structural=%zu "
	       "retained_lod_structural=%zu "
	       "structural_lod=%zu structural_requests=%zu "
	       "structural_entry_geometry=%zu structural_pointer_mismatch=%zu "
	       "missing_instances=%zu mismatched_parts=%zu "
	       "updates=%zu lod_updates=%zu base_parts=%zu lod_parts=%zu "
	       "removed_parts=%zu reset=%d\n",
	       assembly->compactInstancePresentations.size(),
	       assembly->instanceCount(), retainedStructural,
	       retainedLodStructural,
	       structuralLodBacked, structuralSourceRequests,
	       structuralEntryGeometry, structuralGeometryPointerMismatch,
	       missingInstances, mismatchedParts,
	       reset ? this->d->compactIndex->instances.size() :
	       changedInstanceIndices.size(),
	       layerMutation.cuts.size(), baseParts.size(), lodSharedParts.size(),
	       partsToRemove.size(), reset ? 1 : 0);
	if (retainedStructural) {
	    for (const auto &kind : structuralKinds)
		bu_log("  structural kind=%s count=%zu\n",
		    kind.first.c_str(), kind.second);
	    for (const auto &type : structuralTypes)
		bu_log("  structural type=%s count=%zu\n",
		    type.first.c_str(), type.second);
	}
    }

    assembly->compactPresentationInitialized = TRUE;
    assembly->compactPresentationIndex = this->d->compactIndex;
    assembly->compactPresentationSourceRoutingId = this->d->routingId;
    assembly->compactPresentationPopulationEpoch =
	this->d->compactPopulationEpoch;
    assembly->compactPresentationSourceRevision =
	this->sourceRevision.getValue();
    assembly->compactPresentationInputsRevision =
	this->inputsRevision.getValue();
    assembly->compactPresentationCadBatchRevision =
	this->cadBatchRevisionGet();
    assembly->compactPresentationPayloadRevision = payloadRevision;
    assembly->compactPresentationDrawMode = sourceDrawMode;
    assembly->compactWirePresentationCount = wireCount;
    assembly->compactShadedPresentationCount = shadedCount;
    viewState->acknowledgeCadOccurrenceChanges(this, payloadRevision);
    /* A structural leaf box is useful per-object coverage.  Once its retained
     * presentation is committed, retireCompactOverviewAfterPresentation can
     * safely replace the whole-target extent cue. */
    const_cast<SoBRLDatabaseSource *>(this)->
	retireCompactOverviewAfterPresentation(assembly);
    if (presentationTiming) {
	const int64_t presentationCompleted = bu_gettime();
	const int64_t total =
	    presentationCompleted - presentationStarted;
	if (total >= 10000) {
	    bu_log("BObol compact presentation total=%" PRId64 "us "
		   "index=%" PRId64 "us plan=%" PRId64 "us "
		   "publish=%" PRId64 "us reset=%d full_resync=%d "
		   "sparse_payload=%d cad_batch_changed=%d "
		   "changed_keys=%zu structural_entries=%zu entries=%zu "
		   "payload_input=%zu shared_parts=%zu "
		   "instance_updates=%zu lod_updates=%zu\n",
		   total,
		   presentationIndexed - presentationStarted,
		   presentationPlanned - presentationIndexed,
		   presentationCompleted - presentationPlanned,
		   reset ? 1 : 0, fullPayloadResync ? 1 : 0,
		   sparsePayloadUpdate ? 1 : 0,
		   cadBatchChanged ? 1 : 0,
		   changedOccurrenceKeys.size(),
		   structuralEntryIndices.size(), entryIndices.size(),
		   payloadInput->size(), lodSharedParts.size(),
		   reset ? this->d->compactIndex->instances.size() :
		   changedInstanceIndices.size(),
		   layerMutation.cuts.size());
	}
    }
    return assembly;
}

SoBRLCadAssembly *
cad_view_lod_assembly_for_action(SoAction *action,
				 const SoBRLDatabaseSource *source)
{
    const BObolViewLodState *viewState =
	bobol_view_lod_state_for_action(action);
    if (source && (source->isCompactOccurrenceRegistry() ||
	    source->hasDisplayResidentProgressiveGeometry())) {
	if (viewState) {
	    SoBRLCadAssembly *current =
		source->currentCompactViewLodAssembly(viewState);
	    if (current)
		return current;
	    std::vector<const BObolViewLodState::CadPayload *> payloads;
	    return source->compactViewLodAssembly(payloads, viewState);
	}
	/* Aggregate sources require an exact source binding.  The legacy
	 * path/name lookup may otherwise attach a sibling source's payload. */
	return NULL;
    }
    const BObolViewLodState::CadPayload *payload =
	bobol_view_lod_cad_for_action(action, source);
    return cad_view_lod_assembly(source, payload, viewState);
}
