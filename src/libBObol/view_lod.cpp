/*                  V I E W _ L O D . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/str.h"

#include "BObol/BExportAction.h"
#include "BObol/BDatabaseSource.h"
#include "BObol/BMeasureAction.h"
#include "BObol/BMeshLodSubmitAction.h"
#include "BObol/BMeshShape.h"
#include "BObol/BSnapAction.h"
#include "BObol/BViewLod.h"

#include <Obol/cad/SoCADAssembly.h>
#include <Obol/cad/SoCADViewState.h>

#include <Inventor/actions/SoAction.h>
#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/actions/SoPickAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/misc/SoState.h>
#include <Inventor/nodes/SoNode.h>

#include <algorithm>
#include <cmath>
#include <string.h>
#include <utility>
#include <unordered_set>

SO_ELEMENT_SOURCE(SoBRLViewLodElement);
SO_NODE_SOURCE(SoBRLViewLodGroup);

static std::string
view_lod_string_key(const char *prefix, const SbString &value)
{
    if (!prefix || value.getLength() == 0)
	return std::string();

    std::string key(prefix);
    key += ':';
    key += value.getString();
    return key;
}

static std::string
view_lod_shape_primary_key(const SoBRLMeshShape *shape)
{
    if (!shape)
	return std::string();

    std::string key = view_lod_string_key("identity",
					  shape->sourceIdentity.getValue());
    if (!key.empty())
	return key;

    key = view_lod_string_key("path", shape->sourcePath.getValue());
    if (!key.empty())
	return key;

    key = view_lod_string_key("name", shape->sourceName.getValue());
    if (!key.empty())
	return key;

    char buffer[64] = {0};
    snprintf(buffer, sizeof(buffer), "node:%p", static_cast<const void *>(shape));
    return std::string(buffer);
}

static void
view_lod_append_unique_key(std::vector<std::string> &keys,
			   const std::string &key)
{
    if (key.empty())
	return;

    for (size_t i = 0; i < keys.size(); i++) {
	if (keys[i] == key)
	    return;
    }

    keys.push_back(key);
}

static void
view_lod_shape_keys(std::vector<std::string> &keys,
		    const SoBRLMeshShape *shape)
{
    if (!shape)
	return;

    view_lod_append_unique_key(keys, view_lod_string_key("identity",
			       shape->sourceIdentity.getValue()));
    view_lod_append_unique_key(keys, view_lod_string_key("path",
			       shape->sourcePath.getValue()));
    if (shape->sourcePath.getValue().getLength() > 1 &&
	shape->sourcePath.getValue().getString()[0] == '/')
	view_lod_append_unique_key(keys, std::string("path:") +
				   (shape->sourcePath.getValue().getString() + 1));
    view_lod_append_unique_key(keys, view_lod_string_key("name",
			       shape->sourceName.getValue()));
    view_lod_append_unique_key(keys, view_lod_shape_primary_key(shape));
}

static void
view_lod_result_keys(std::vector<std::string> &keys,
		     const BObolLodResult &result)
{
    /* A compact occurrence may share both object path and object name with
     * siblings.  Its immutable identity is consequently the only valid
     * display binding key. */
    if (result.request.occurrenceKey.getLength() > 0) {
	view_lod_append_unique_key(keys, view_lod_string_key("occurrence",
			       result.request.occurrenceKey));
	return;
    }
    view_lod_append_unique_key(keys, view_lod_string_key("path",
			       result.request.objectPath));
    if (result.request.objectPath.getLength() > 1 &&
	result.request.objectPath.getString()[0] == '/')
	view_lod_append_unique_key(keys, std::string("path:") +
				   (result.request.objectPath.getString() + 1));
    view_lod_append_unique_key(keys, view_lod_string_key("name",
			       result.request.objectName));
}

static const char *
view_lod_leaf_name_from_path(const SbString &path)
{
    const char *p = path.getString();
    if (!p || !p[0])
	return "";
    const char *slash = strrchr(p, '/');
    return (slash && slash[1]) ? slash + 1 : p;
}

static bool
view_lod_paths_equal(const SbString &a, const SbString &b)
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

static std::string
view_lod_source_primary_key(const SoBRLDatabaseSource *source)
{
    if (!source)
	return std::string();

    std::string key = view_lod_string_key("instance",
					  source->instanceKey.getValue());
    if (!key.empty())
	return key;

    key = view_lod_string_key("realization",
			      source->realizationIdentity.getValue());
    if (!key.empty())
	return key;

    key = view_lod_string_key("path", source->path.getValue());
    if (!key.empty())
	return key;

    char buffer[64] = {0};
    snprintf(buffer, sizeof(buffer), "source:%p",
	     static_cast<const void *>(source));
    return std::string(buffer);
}

static void
view_lod_source_keys(std::vector<std::string> &keys,
		     const SoBRLDatabaseSource *source)
{
    if (!source)
	return;

    view_lod_append_unique_key(keys, view_lod_string_key("instance",
			       source->instanceKey.getValue()));
    view_lod_append_unique_key(keys, view_lod_string_key("realization",
			       source->realizationIdentity.getValue()));
    view_lod_append_unique_key(keys, view_lod_string_key("path",
			       source->path.getValue()));
    if (source->path.getValue().getLength() > 1 &&
	source->path.getValue().getString()[0] == '/')
	view_lod_append_unique_key(keys, std::string("path:") +
				   (source->path.getValue().getString() + 1));
    view_lod_append_unique_key(keys, view_lod_string_key("name",
			       SbString(view_lod_leaf_name_from_path(
					    source->path.getValue()))));
    view_lod_append_unique_key(keys, view_lod_source_primary_key(source));
}

BObolViewLodState::MeshPayload::MeshPayload(void) :
    sourceContentHash(0),
    resultKind(BOBOL_LOD_RESULT_NONE),
    qualityTier(BOBOL_LOD_QUALITY_METADATA),
    providerStatus(BOBOL_LOD_PROVIDER_UNKNOWN),
    activeLevel(-1),
    residentLevel(-1),
    requestedLevel(-1),
    viewRevision(0),
    policyRevision(0),
    hasSnappedPoints(FALSE),
    hasNormals(FALSE),
    shadedCullBackfaces(FALSE)
{
}

SbBool
BObolViewLodState::MeshPayload::isValid(void) const
{
    return this->mesh.isValid() ||
	(this->progressiveMesh && this->progressiveMesh->isValid());
}

size_t
BObolViewLodState::MeshPayload::estimateBytes(void) const
{
    return sizeof(*this) + (this->progressiveMesh ? 0 :
	    this->mesh.points.size() * sizeof(SbVec3f) +
	   this->mesh.normals.size() * sizeof(SbVec3f) +
	   this->mesh.coordIndex.size() * sizeof(int32_t) +
	   this->mesh.faceIndex.size() * sizeof(int32_t) +
	   this->mesh.vertexIndex.size() * sizeof(int32_t));
}

int
BObolViewLodState::MeshPayload::getTriangleCount(void) const
{
    return static_cast<int>(this->mesh.coordIndex.size() / 3);
}

SbBool
BObolViewLodState::MeshPayload::getTriangleVertexIndices(int triangleIndex,
	int &indexA,
	int &indexB,
	int &indexC) const
{
    indexA = -1;
    indexB = -1;
    indexC = -1;
    if (triangleIndex < 0)
	return FALSE;

    size_t base = static_cast<size_t>(triangleIndex) * 3;
    if (base + 2 >= this->mesh.coordIndex.size())
	return FALSE;

    indexA = this->mesh.coordIndex[base];
    indexB = this->mesh.coordIndex[base + 1];
    indexC = this->mesh.coordIndex[base + 2];
    if (indexA < 0 || indexB < 0 || indexC < 0 ||
	static_cast<size_t>(indexA) >= this->mesh.points.size() ||
	static_cast<size_t>(indexB) >= this->mesh.points.size() ||
	static_cast<size_t>(indexC) >= this->mesh.points.size()) {
	indexA = -1;
	indexB = -1;
	indexC = -1;
	return FALSE;
    }

    return TRUE;
}

SbBool
BObolViewLodState::MeshPayload::getTriangle(int triangleIndex,
	SbVec3f &a,
	SbVec3f &b,
	SbVec3f &c) const
{
    int ia = -1;
    int ib = -1;
    int ic = -1;
    if (!this->getTriangleVertexIndices(triangleIndex, ia, ib, ic))
	return FALSE;

    a = this->mesh.points[static_cast<size_t>(ia)];
    b = this->mesh.points[static_cast<size_t>(ib)];
    c = this->mesh.points[static_cast<size_t>(ic)];
    return TRUE;
}

BObolViewLodState::ProxyPayload::ProxyPayload(void) :
    resultKind(BOBOL_LOD_RESULT_NONE),
    qualityTier(BOBOL_LOD_QUALITY_METADATA),
    providerStatus(BOBOL_LOD_PROVIDER_UNKNOWN),
    viewRevision(0),
    policyRevision(0),
    diagnostic("")
{
}

SbBool
BObolViewLodState::ProxyPayload::isValid(void) const
{
    return this->proxy.isValid();
}

size_t
BObolViewLodState::ProxyPayload::estimateBytes(void) const
{
    return sizeof(*this);
}

BObolViewLodState::CadPayload::CadPayload(void) :
    sourceContentHash(0),
    resultKind(BOBOL_LOD_RESULT_NONE),
    qualityTier(BOBOL_LOD_QUALITY_METADATA),
    providerStatus(BOBOL_LOD_PROVIDER_UNKNOWN),
    drawMode(BOBOL_LOD_DRAW_UNKNOWN),
    activeLevel(-1),
    residentLevel(-1),
    requestedLevel(-1),
    viewRevision(0),
    policyRevision(0),
    hasSnappedPoints(FALSE),
    hasNormals(FALSE),
    shadedCullBackfaces(FALSE)
{
}

SbBool
BObolViewLodState::CadPayload::isValid(void) const
{
    if (this->providerStatus != BOBOL_LOD_PROVIDER_READY)
	return FALSE;
    if (this->resultKind == BOBOL_LOD_RESULT_MESH ||
	this->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL)
	return this->mesh.isValid() ||
	    (this->progressiveMesh && this->progressiveMesh->isValid());
    if (this->resultKind == BOBOL_LOD_RESULT_AABB ||
	this->resultKind == BOBOL_LOD_RESULT_PROXY)
	return this->proxy.isValid();
    return FALSE;
}

size_t
BObolViewLodState::CadPayload::estimateBytes(void) const
{
    return (this->progressiveMesh ? 0 :
	    this->mesh.points.size() * sizeof(SbVec3f) +
	   this->mesh.normals.size() * sizeof(SbVec3f) +
	   this->mesh.coordIndex.size() * sizeof(int32_t) +
	   this->mesh.faceIndex.size() * sizeof(int32_t) +
	   this->mesh.vertexIndex.size() * sizeof(int32_t)) +
	   sizeof(*this);
}

/* One compact occurrence has one active display representation.  Retaining
 * every completed stage both wastes memory and makes an unordered binding map
 * choose AABB versus OBB versus mesh by accident. */
static int
view_lod_cad_payload_rank(const BObolViewLodState::CadPayload &payload)
{
    int rank = payload.qualityTier * 16;
    switch (payload.resultKind) {
	case BOBOL_LOD_RESULT_FULL_DETAIL:
	    return 512 + rank;
	case BOBOL_LOD_RESULT_MESH:
	    return 448 + rank;
	case BOBOL_LOD_RESULT_PROXY:
	    if (payload.proxy.kind == BOBOL_LOD_PROXY_OBB)
		return 320 + rank;
	    return 256 + rank;
	case BOBOL_LOD_RESULT_AABB:
	    return 192 + rank;
	default:
	    return rank;
    }
}

static SbBool
view_lod_cad_payload_targets(const BObolViewLodState::CadPayload &payload,
	const std::string &sourceBindingKey, const SbString &occurrenceKey)
{
    return bu_strcmp(payload.sourceBindingKey.getString(),
		     sourceBindingKey.c_str()) == 0 &&
	bu_strcmp(payload.sourceInstanceKey.getString(),
		  occurrenceKey.getString()) == 0 ? TRUE : FALSE;
}

static std::string
view_lod_cad_result_binding_key(const std::string &sourceBindingKey,
	const std::string &resultKey)
{
    std::string key = "source:";
    key += sourceBindingKey;
    key += '|';
    key += resultKey;
    return key;
}

static std::string
view_lod_cad_occurrence_key(const SbString &occurrenceKey)
{
    return occurrenceKey.getLength() > 0 ?
	std::string(occurrenceKey.getString()) : std::string();
}

static void
view_lod_remove_superseded_cad_payloads(
	std::unordered_map<std::string, BObolViewLodState::CadPayloadPtr>
	&bindings, const std::string &sourceBindingKey,
	const SbString &occurrenceKey,
	const BObolViewLodState::CadPayload &candidate)
{
    const int candidateRank = view_lod_cad_payload_rank(candidate);
    for (auto binding = bindings.begin(); binding != bindings.end();) {
	const BObolViewLodState::CadPayloadPtr &payload = binding->second;
	if (payload && view_lod_cad_payload_targets(*payload,
		sourceBindingKey, occurrenceKey) &&
	    view_lod_cad_payload_rank(*payload) <= candidateRank) {
	    binding = bindings.erase(binding);
	    continue;
	}
	++binding;
    }
}

BObolViewLodState::BObolViewLodState(void) :
    cadFullResyncRevision(1),
    cadBindingsRevision(1),
    normalStyle(NORMAL_AUTHORED),
    normalCreaseAngle(60.0f),
    cadPresentationProgressiveLodCeiling(-1)
{
}

BObolViewLodState::~BObolViewLodState(void)
{
    this->clearCadPresentations();
}

void
BObolViewLodState::clear(void)
{
    this->clearCadPresentations();
    this->meshBindings.clear();
    this->proxyBindings.clear();
    this->cadSourceBindings.clear();
    this->cadAssetBindings.clear();
    this->cadBindings.clear();
    this->cadOccurrenceChanges.clear();
    this->noteResidentMeshesChanged();
}

SoCADAssembly *
BObolViewLodState::findCadPresentation(
    const SoBRLDatabaseSource *source,
    SbString *contentKey) const
{
    if (!source)
	return NULL;
    const std::string key = view_lod_source_primary_key(source);
    const auto found = this->cadPresentations.find(key);
    if (found == this->cadPresentations.end())
	return NULL;
    if (contentKey)
	*contentKey = found->second.contentKey;
    return found->second.assembly;
}

void
BObolViewLodState::setCadPresentation(
    const SoBRLDatabaseSource *source,
    SoCADAssembly *assembly,
    const SbString &contentKey) const
{
    if (!source)
	return;
    const std::string key = view_lod_source_primary_key(source);
    CadPresentation &presentation = this->cadPresentations[key];
    if (presentation.assembly != assembly) {
	if (assembly)
	    assembly->ref();
	if (presentation.assembly)
	    presentation.assembly->unref();
	presentation.assembly = assembly;
    }
    presentation.contentKey = contentKey;
    if (presentation.assembly)
	presentation.assembly->progressiveLodCeiling =
	    this->cadPresentationProgressiveLodCeiling;
    if (!presentation.assembly)
	this->cadPresentations.erase(key);
}

void
BObolViewLodState::clearCadPresentations(void) const
{
    for (auto &binding : this->cadPresentations) {
	if (binding.second.assembly)
	    binding.second.assembly->unref();
	binding.second.assembly = NULL;
    }
    this->cadPresentations.clear();
}

void
BObolViewLodState::setNormalStyle(NormalStyle style, float creaseAngleDegrees)
{
    if (style < NORMAL_AUTHORED || style > NORMAL_SMOOTH)
	style = NORMAL_AUTHORED;
    if (!std::isfinite(creaseAngleDegrees))
	creaseAngleDegrees = 60.0f;
    creaseAngleDegrees = std::max(0.0f,
	std::min(180.0f, creaseAngleDegrees));
    if (this->normalStyle == style &&
	std::fabs(this->normalCreaseAngle - creaseAngleDegrees) <= 1.0e-6f)
	return;
    this->normalStyle = style;
    this->normalCreaseAngle = creaseAngleDegrees;
    /* This is a presentation policy change.  Retain all view-local PoP
     * payloads and invalidate only their renderer assemblies. */
    this->clearCadPresentations();
}

BObolViewLodState::NormalStyle
BObolViewLodState::getNormalStyle(void) const
{
    return this->normalStyle;
}

float
BObolViewLodState::getNormalCreaseAngle(void) const
{
    return this->normalCreaseAngle;
}

SbBool
BObolViewLodState::applyMeshResult(const SoBRLMeshShape *shape,
				     const BObolLodResult &result)
{
    BObolLodResult copy = result;
    return this->applyMeshResultInternal(shape, copy, FALSE);
}

SbBool
BObolViewLodState::applyMeshResultInternal(const SoBRLMeshShape *shape,
	BObolLodResult &result, SbBool consume)
{
    if (!shape ||
	result.resultKind != BOBOL_LOD_RESULT_MESH ||
	result.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	(!result.mesh.isValid() &&
	 (!result.progressiveMesh || !result.progressiveMesh->isValid())))
	return FALSE;

    MeshPayloadPtr payload(new MeshPayload);
    payload->progressiveMesh = result.progressiveMesh;
    if (consume)
	payload->mesh = std::move(result.mesh);
    else
	payload->mesh = result.mesh;
    payload->sourcePath = shape->sourcePath.getValue();
    payload->sourceName = shape->sourceName.getValue();
    payload->sourceIdentity = shape->sourceIdentity.getValue();
    payload->cacheIdentity = result.cacheKey.value;
    payload->cacheKey = result.geometry.cacheKey.isValid() ?
	result.geometry.cacheKey.value : result.cacheKey.value;
    payload->sourceContentHash = result.request.sourceContentHash;
    payload->resultKind = result.resultKind;
    payload->qualityTier = result.qualityTier;
    payload->providerStatus = result.providerStatus;
    payload->activeLevel = result.geometry.activeLevel;
    payload->residentLevel = result.residentLevel;
    payload->requestedLevel = result.request.requestedLevel;
    payload->viewRevision = result.request.viewRevision;
    payload->policyRevision = result.request.policyRevision;
    payload->counts = result.counts;
    payload->bounds = result.bounds;
    payload->hasSnappedPoints = result.hasSnappedPoints;
    payload->hasNormals = result.hasNormals;
    payload->shadedCullBackfaces = result.shadedCullBackfaces;
    payload->diagnostic = result.diagnostic;

    std::vector<std::string> keys;
    view_lod_shape_keys(keys, shape);
    view_lod_result_keys(keys, result);
    for (size_t i = 0; i < keys.size(); i++)
	this->meshBindings[keys[i]] = payload;

    return TRUE;
}

SbBool
BObolViewLodState::applyProxyResult(const SoBRLMeshShape *shape,
				      const BObolLodResult &result)
{
    BObolLodResult copy = result;
    return this->applyProxyResultInternal(shape, copy, FALSE);
}

SbBool
BObolViewLodState::applyProxyResultInternal(const SoBRLMeshShape *shape,
	BObolLodResult &result, SbBool consume)
{
    if (!shape ||
	(result.resultKind != BOBOL_LOD_RESULT_AABB &&
	 result.resultKind != BOBOL_LOD_RESULT_PROXY) ||
	result.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	!result.proxy.isValid())
	return FALSE;

    ProxyPayloadPtr payload(new ProxyPayload);
    if (consume)
	payload->proxy = std::move(result.proxy);
    else
	payload->proxy = result.proxy;
    payload->sourcePath = shape->sourcePath.getValue();
    payload->sourceName = shape->sourceName.getValue();
    payload->sourceIdentity = shape->sourceIdentity.getValue();
    payload->cacheIdentity = shape->cacheIdentity.getValue();
    payload->cacheKey = result.cacheKey.value;
    payload->resultKind = result.resultKind;
    payload->qualityTier = result.qualityTier;
    payload->providerStatus = result.providerStatus;
    payload->viewRevision = result.request.viewRevision;
    payload->policyRevision = result.request.policyRevision;
    payload->counts = result.counts;
    payload->bounds = result.bounds;
    payload->diagnostic = result.diagnostic;

    std::vector<std::string> keys;
    view_lod_shape_keys(keys, shape);
    view_lod_result_keys(keys, result);
    for (size_t i = 0; i < keys.size(); i++)
	this->proxyBindings[keys[i]] = payload;

    return TRUE;
}

SbBool
BObolViewLodState::applyDisplayResult(const SoBRLMeshShape *shape,
					const BObolLodResult &result)
{
    if (result.resultKind == BOBOL_LOD_RESULT_MESH ||
	result.resultKind == BOBOL_LOD_RESULT_FULL_DETAIL)
	return this->applyMeshResult(shape, result);

    return this->applyProxyResult(shape, result);
}

SbBool
BObolViewLodState::consumeDisplayResult(const SoBRLMeshShape *shape,
	BObolLodResult &result)
{
    if (result.resultKind == BOBOL_LOD_RESULT_MESH ||
	result.resultKind == BOBOL_LOD_RESULT_FULL_DETAIL)
	return this->applyMeshResultInternal(shape, result, TRUE);
    return this->applyProxyResultInternal(shape, result, TRUE);
}

SbBool
BObolViewLodState::applySourceResult(
    const SoBRLDatabaseSource *source,
    const BObolLodResult &result)
{
    BObolLodResult copy = result;
    return this->applySourceResultInternal(source, copy, FALSE);
}

SbBool
BObolViewLodState::consumeSourceResult(
    const SoBRLDatabaseSource *source,
    BObolLodResult &result)
{
    return this->applySourceResultInternal(source, result, TRUE);
}

SbBool
BObolViewLodState::applySourceResultInternal(
    const SoBRLDatabaseSource *source,
    BObolLodResult &result,
    SbBool consume)
{
    if (!source || result.providerStatus != BOBOL_LOD_PROVIDER_READY)
	return FALSE;

    if ((result.resultKind == BOBOL_LOD_RESULT_MESH ||
	 result.resultKind == BOBOL_LOD_RESULT_FULL_DETAIL) &&
	!result.mesh.isValid() &&
	(!result.progressiveMesh || !result.progressiveMesh->isValid()))
	return FALSE;
    if ((result.resultKind == BOBOL_LOD_RESULT_AABB ||
	 result.resultKind == BOBOL_LOD_RESULT_PROXY) &&
	!result.proxy.isValid())
	return FALSE;
    if (result.resultKind != BOBOL_LOD_RESULT_MESH &&
	result.resultKind != BOBOL_LOD_RESULT_FULL_DETAIL &&
	result.resultKind != BOBOL_LOD_RESULT_AABB &&
	result.resultKind != BOBOL_LOD_RESULT_PROXY)
	return FALSE;

    CadPayloadPtr payload(new CadPayload);
    payload->progressiveMesh = result.progressiveMesh;
    if (consume) {
	payload->mesh = std::move(result.mesh);
	payload->proxy = std::move(result.proxy);
    } else {
	payload->mesh = result.mesh;
	payload->proxy = result.proxy;
    }
    payload->sourcePath = result.request.objectPath.getLength() > 0 ?
	result.request.objectPath : source->path.getValue();
    payload->sourceName = result.request.objectName.getLength() > 0 ?
	result.request.objectName :
	SbString(view_lod_leaf_name_from_path(payload->sourcePath));
    payload->sourceIdentity =
	source->realizationIdentity.getValue().getLength() > 0 ?
	source->realizationIdentity.getValue() : source->path.getValue();
    /* This field intentionally means a compact occurrence key, not the
     * source-node key.  An empty value is a source-wide legacy binding that
     * may fall back to a path; a nonempty value must never alias a sibling. */
    payload->sourceInstanceKey = result.request.occurrenceKey;
    const std::string sourceBindingKey = view_lod_source_primary_key(source);
    payload->sourceBindingKey = sourceBindingKey.c_str();
    payload->cacheIdentity = result.cacheKey.value;
    payload->cacheKey = result.geometry.cacheKey.isValid() ?
	result.geometry.cacheKey.value : result.cacheKey.value;
    payload->sourceContentHash = result.request.sourceContentHash;
    payload->resultKind = result.resultKind;
    payload->qualityTier = result.qualityTier;
    payload->providerStatus = result.providerStatus;
    payload->drawMode = result.request.drawMode;
    payload->activeLevel = result.geometry.activeLevel;
    payload->residentLevel = result.residentLevel;
    payload->requestedLevel = result.request.requestedLevel;
    payload->viewRevision = result.request.viewRevision;
    payload->policyRevision = result.request.policyRevision;
    payload->counts = result.counts;
    payload->bounds = result.bounds;
    payload->hasSnappedPoints = result.hasSnappedPoints;
    payload->hasNormals = result.hasNormals;
    payload->shadedCullBackfaces = result.shadedCullBackfaces;
    payload->diagnostic = result.diagnostic;

    const std::string occurrenceKey =
	view_lod_cad_occurrence_key(payload->sourceInstanceKey);
    std::unordered_map<std::string, CadPayloadPtr> &sourcePayloads =
	this->cadSourceBindings[sourceBindingKey];
    const auto current = sourcePayloads.find(occurrenceKey);
    if (current != sourcePayloads.end() && current->second &&
	view_lod_cad_payload_rank(*current->second) >
	    view_lod_cad_payload_rank(*payload))
	return TRUE;

    /* A compact occurrence has exactly one alias and one authoritative
     * source/occurrence slot.  Replacing it by direct key is O(1).  Retain
     * the broader alias cleanup only for legacy source-wide results. */
    if (payload->sourceInstanceKey.getLength() == 0)
	view_lod_remove_superseded_cad_payloads(this->cadBindings,
	    sourceBindingKey, payload->sourceInstanceKey, *payload);

    std::vector<std::string> resultKeys;
    view_lod_result_keys(resultKeys, result);
    for (const std::string &resultKey : resultKeys) {
	const std::string key =
	    view_lod_cad_result_binding_key(sourceBindingKey, resultKey);
	this->cadBindings[key] = payload;
    }
    sourcePayloads[occurrenceKey] = payload;
    if (payload->progressiveMesh && payload->progressiveMesh->isValid() &&
	payload->sourcePath.getLength() > 0) {
	std::unordered_map<std::string, CadPayloadPtr> &sourceAssets =
	    this->cadAssetBindings[sourceBindingKey];
	const std::string assetKey = payload->sourcePath.getString();
	const auto resident = sourceAssets.find(assetKey);
	if (resident == sourceAssets.end() || !resident->second ||
	    bu_strcmp(resident->second->cacheKey.getString(),
		payload->cacheKey.getString()) != 0 ||
	    resident->second->residentLevel <= payload->residentLevel)
	    sourceAssets[assetKey] = payload;
    }

    /* Preserve the source-wide lookup for a root result.  Entry results must
     * coexist, so binding them under the root keys would make the last
     * completed occurrence replace every sibling. */
	if (result.request.occurrenceKey.getLength() == 0 &&
	(result.request.objectPath.getLength() == 0 ||
	view_lod_paths_equal(source->path.getValue(),
	    result.request.objectPath))) {
	std::vector<std::string> sourceKeys;
	view_lod_source_keys(sourceKeys, source);
	for (const std::string &key : sourceKeys)
	    this->cadBindings[key] = payload;
    }

    if (payload->sourceInstanceKey.getLength() > 0)
	this->noteCadOccurrenceChanged(
	    sourceBindingKey, payload->sourceInstanceKey);
    else
	this->noteResidentMeshesChanged();
    return TRUE;
}

const BObolViewLodState::MeshPayload *
BObolViewLodState::findMesh(const SoBRLMeshShape *shape) const
{
    if (this->meshBindings.empty())
	return NULL;

    std::vector<std::string> keys;
    view_lod_shape_keys(keys, shape);
    for (size_t i = 0; i < keys.size(); i++) {
	std::unordered_map<std::string, MeshPayloadPtr>::const_iterator it =
	    this->meshBindings.find(keys[i]);
	if (it != this->meshBindings.end() && it->second &&
	    it->second->isValid())
	    return it->second.get();
    }

    return NULL;
}

const BObolViewLodState::MeshPayload *
BObolViewLodState::findMeshForResult(
    const BObolLodResult &result) const
{
    if (this->meshBindings.empty())
	return NULL;

    std::vector<std::string> keys;
    view_lod_result_keys(keys, result);
    for (size_t i = 0; i < keys.size(); i++) {
	std::unordered_map<std::string, MeshPayloadPtr>::const_iterator it =
	    this->meshBindings.find(keys[i]);
	if (it != this->meshBindings.end() && it->second &&
	    it->second->isValid())
	    return it->second.get();
    }

    return NULL;
}

const BObolViewLodState::ProxyPayload *
BObolViewLodState::findProxy(const SoBRLMeshShape *shape) const
{
    if (this->proxyBindings.empty())
	return NULL;

    std::vector<std::string> keys;
    view_lod_shape_keys(keys, shape);
    for (size_t i = 0; i < keys.size(); i++) {
	std::unordered_map<std::string, ProxyPayloadPtr>::const_iterator it =
	    this->proxyBindings.find(keys[i]);
	if (it != this->proxyBindings.end() && it->second &&
	    it->second->isValid())
	    return it->second.get();
    }

    return NULL;
}

const BObolViewLodState::ProxyPayload *
BObolViewLodState::findProxyForResult(
    const BObolLodResult &result) const
{
    if (this->proxyBindings.empty())
	return NULL;

    std::vector<std::string> keys;
    view_lod_result_keys(keys, result);
    for (size_t i = 0; i < keys.size(); i++) {
	std::unordered_map<std::string, ProxyPayloadPtr>::const_iterator it =
	    this->proxyBindings.find(keys[i]);
	if (it != this->proxyBindings.end() && it->second &&
	    it->second->isValid())
	    return it->second.get();
    }

    return NULL;
}

const BObolViewLodState::CadPayload *
BObolViewLodState::findCad(const SoBRLDatabaseSource *source) const
{
    if (this->cadBindings.empty())
	return NULL;

    std::vector<std::string> keys;
    view_lod_source_keys(keys, source);
    for (size_t i = 0; i < keys.size(); i++) {
	std::unordered_map<std::string, CadPayloadPtr>::const_iterator it =
	    this->cadBindings.find(keys[i]);
	if (it != this->cadBindings.end() && it->second &&
	    it->second->isValid())
	    return it->second.get();
    }

    return NULL;
}

void
BObolViewLodState::findCadPayloadsUnordered(
    const SoBRLDatabaseSource *source,
    std::vector<const CadPayload *> &payloads) const
{
    payloads.clear();
    if (!source || this->cadSourceBindings.empty())
	return;

    const std::string sourceKey = view_lod_source_primary_key(source);
    const auto sourcePayloads = this->cadSourceBindings.find(sourceKey);
    if (sourcePayloads == this->cadSourceBindings.end())
	return;

    payloads.reserve(sourcePayloads->second.size());
    for (const auto &binding : sourcePayloads->second) {
	const CadPayloadPtr &payload = binding.second;
	if (payload && payload->isValid())
	    payloads.push_back(payload.get());
    }
}

void
BObolViewLodState::findCadPayloads(
    const SoBRLDatabaseSource *source,
    std::vector<const CadPayload *> &payloads) const
{
    this->findCadPayloadsUnordered(source, payloads);
    std::sort(payloads.begin(), payloads.end(),
	[](const CadPayload *a, const CadPayload *b) {
	    const int instanceOrder = bu_strcmp(a->sourceInstanceKey.getString(),
		b->sourceInstanceKey.getString());
	    if (instanceOrder != 0)
		return instanceOrder < 0;
	    const int pathOrder = bu_strcmp(a->sourcePath.getString(),
		b->sourcePath.getString());
	    if (pathOrder != 0)
		return pathOrder < 0;
	    return bu_strcmp(a->cacheKey.getString(),
		b->cacheKey.getString()) < 0;
	});
}

const BObolViewLodState::CadPayload *
BObolViewLodState::findCadForOccurrence(
    const SoBRLDatabaseSource *source,
    const SbString &occurrenceKey) const
{
    if (!source || occurrenceKey.getLength() == 0 ||
	this->cadSourceBindings.empty())
	return NULL;
    const auto sourcePayloads = this->cadSourceBindings.find(
	view_lod_source_primary_key(source));
    if (sourcePayloads == this->cadSourceBindings.end())
	return NULL;
    const auto payload =
	sourcePayloads->second.find(occurrenceKey.getString());
    return payload != sourcePayloads->second.end() && payload->second &&
	payload->second->isValid() ? payload->second.get() : NULL;
}

const BObolViewLodState::CadPayload *
BObolViewLodState::findCadForAsset(
    const SoBRLDatabaseSource *source,
    const SbString &assetPath) const
{
    if (!source || assetPath.getLength() == 0 ||
	this->cadAssetBindings.empty())
	return NULL;
    const auto sourceAssets = this->cadAssetBindings.find(
	view_lod_source_primary_key(source));
    if (sourceAssets == this->cadAssetBindings.end())
	return NULL;
    const auto payload = sourceAssets->second.find(assetPath.getString());
    return payload != sourceAssets->second.end() && payload->second &&
	payload->second->progressiveMesh &&
	payload->second->progressiveMesh->isValid() ?
	payload->second.get() : NULL;
}

const BObolViewLodState::CadPayload *
BObolViewLodState::findCadForResult(
    const SoBRLDatabaseSource *source,
    const BObolLodResult &result) const
{
    if (!source || this->cadSourceBindings.empty())
	return NULL;

    const std::string sourceKey = view_lod_source_primary_key(source);
    const auto sourcePayloads = this->cadSourceBindings.find(sourceKey);
    if (sourcePayloads == this->cadSourceBindings.end())
	return NULL;
    const std::string occurrenceKey =
	view_lod_cad_occurrence_key(result.request.occurrenceKey);
    const auto payload = sourcePayloads->second.find(occurrenceKey);
    if (payload == sourcePayloads->second.end() || !payload->second ||
	!payload->second->isValid())
	return NULL;

    if (result.request.occurrenceKey.getLength() > 0)
	return payload->second.get();
    if (view_lod_paths_equal(payload->second->sourcePath,
	    result.request.objectPath) ||
	(payload->second->sourceName.getLength() > 0 &&
	 result.request.objectName.getLength() > 0 &&
	 bu_strcmp(payload->second->sourceName.getString(),
	    result.request.objectName.getString()) == 0))
	return payload->second.get();
    return NULL;
}

const BObolViewLodState::CadPayload *
BObolViewLodState::findCadForResult(
    const BObolLodResult &result) const
{
    if (this->cadSourceBindings.empty())
	return NULL;

    for (const auto &source : this->cadSourceBindings) {
	for (const auto &binding : source.second) {
	    const CadPayloadPtr &payload = binding.second;
	    if (!payload || !payload->isValid())
		continue;
	if (result.request.occurrenceKey.getLength() > 0) {
	    if (bu_strcmp(payload->sourceInstanceKey.getString(),
		result.request.occurrenceKey.getString()) == 0)
		return payload.get();
	    continue;
	}
	if (view_lod_paths_equal(payload->sourcePath,
		result.request.objectPath) ||
	    (payload->sourceName.getLength() > 0 &&
	     result.request.objectName.getLength() > 0 &&
	     bu_strcmp(payload->sourceName.getString(),
		result.request.objectName.getString()) == 0))
	    return payload.get();
	}
    }

    return NULL;
}

SbBool
BObolViewLodState::retargetMeshPayload(
    const BObolViewLodState::MeshPayload *target,
    int activeLevel,
    int requestedLevel,
    uint64_t viewRevision,
    uint64_t policyRevision)
{
    if (!target || !target->progressiveMesh ||
	!target->progressiveMesh->canDrawLevel(activeLevel) ||
	requestedLevel < 0)
	return FALSE;

    MeshPayloadPtr payload;
    for (const auto &binding : this->meshBindings) {
	if (binding.second && binding.second.get() == target) {
	    payload = binding.second;
	    break;
	}
    }
    if (!payload)
	return FALSE;

    payload->activeLevel = activeLevel;
    payload->requestedLevel = requestedLevel;
    payload->residentLevel = payload->progressiveMesh->residentLevel();
    payload->viewRevision = viewRevision;
    payload->policyRevision = policyRevision;
    payload->counts.faceCount =
	payload->progressiveMesh->faceCount(activeLevel);
    payload->counts.pointCount =
	payload->progressiveMesh->pointCount(activeLevel);
    payload->counts.originalPointCount = payload->counts.pointCount;
    payload->counts.normalCount = payload->hasNormals ?
	payload->counts.faceCount * 3 : 0;
    return TRUE;
}

SbBool
BObolViewLodState::retargetCadPayload(
    const BObolViewLodState::CadPayload *target,
    int activeLevel,
    int requestedLevel,
    uint64_t viewRevision,
    uint64_t policyRevision)
{
    if (!target || !target->progressiveMesh ||
	!target->progressiveMesh->canDrawLevel(activeLevel) ||
	requestedLevel < 0)
	return FALSE;

    CadPayloadPtr payload;
    const auto source = this->cadSourceBindings.find(
	target->sourceBindingKey.getString());
    if (source != this->cadSourceBindings.end()) {
	const auto occurrence = source->second.find(
	    view_lod_cad_occurrence_key(target->sourceInstanceKey));
	if (occurrence != source->second.end() &&
	    occurrence->second.get() == target)
	    payload = occurrence->second;
    }
    if (!payload)
	return FALSE;

    payload->activeLevel = activeLevel;
    payload->requestedLevel = requestedLevel;
    payload->residentLevel = payload->progressiveMesh->residentLevel();
    payload->viewRevision = viewRevision;
    payload->policyRevision = policyRevision;
    payload->counts.faceCount =
	payload->progressiveMesh->faceCount(activeLevel);
    payload->counts.pointCount =
	payload->progressiveMesh->pointCount(activeLevel);
    payload->counts.originalPointCount = payload->counts.pointCount;
    payload->counts.normalCount = payload->hasNormals ?
	payload->counts.faceCount * 3 : 0;
    this->noteCadOccurrenceChanged(
	payload->sourceBindingKey.getString(), payload->sourceInstanceKey);
    return TRUE;
}

SbBool
BObolViewLodState::removeCadPayload(
    const BObolViewLodState::CadPayload *target)
{
    if (!target)
	return FALSE;
    SbBool removed = FALSE;
    const std::string sourceKey = target->sourceBindingKey.getString();
    const std::string occurrenceKey =
	view_lod_cad_occurrence_key(target->sourceInstanceKey);
    auto source = this->cadSourceBindings.find(sourceKey);
    if (source != this->cadSourceBindings.end()) {
	auto occurrence = source->second.find(occurrenceKey);
	if (occurrence != source->second.end() &&
	    occurrence->second.get() == target) {
	    source->second.erase(occurrence);
	    if (source->second.empty())
		this->cadSourceBindings.erase(source);
	    removed = TRUE;
	}
    }

    if (target->sourceInstanceKey.getLength() > 0) {
	const std::string resultKey = std::string("occurrence:") +
	    target->sourceInstanceKey.getString();
	const std::string aliasKey =
	    view_lod_cad_result_binding_key(sourceKey, resultKey);
	auto alias = this->cadBindings.find(aliasKey);
	if (alias != this->cadBindings.end() &&
	    alias->second.get() == target)
	    this->cadBindings.erase(alias);
    } else {
	for (auto binding = this->cadBindings.begin();
	     binding != this->cadBindings.end();) {
	    if (binding->second && binding->second.get() == target)
		binding = this->cadBindings.erase(binding);
	    else
		++binding;
	}
    }
    if (removed) {
	if (target->sourceInstanceKey.getLength() > 0)
	    this->noteCadOccurrenceChanged(
		sourceKey, target->sourceInstanceKey);
	else {
	    this->noteResidentMeshesChanged();
	    this->clearCadPresentations();
	}
    }
    return removed;
}

SbBool
BObolViewLodState::removeMeshPayload(
    const BObolViewLodState::MeshPayload *target)
{
    if (!target)
	return FALSE;
    SbBool removed = FALSE;
    for (auto binding = this->meshBindings.begin();
	 binding != this->meshBindings.end();) {
	if (binding->second && binding->second.get() == target) {
	    binding = this->meshBindings.erase(binding);
	    removed = TRUE;
	} else {
	    ++binding;
	}
    }
    return removed;
}

size_t
BObolViewLodState::bindingCount(void) const
{
    return this->meshBindings.size() + this->proxyBindings.size() +
	   this->cadBindings.size();
}

static std::vector<BObolViewLodState::MeshPayloadPtr>
view_lod_unique_payloads(
    const std::unordered_map<std::string,
    BObolViewLodState::MeshPayloadPtr> &bindings)
{
    std::vector<BObolViewLodState::MeshPayloadPtr> payloads;
    std::unordered_set<const BObolViewLodState::MeshPayload *> seen;
    seen.reserve(bindings.size());
    for (std::unordered_map<std::string,
	 BObolViewLodState::MeshPayloadPtr>::const_iterator it =
	     bindings.begin(); it != bindings.end(); ++it) {
	if (!it->second)
	    continue;
	if (seen.insert(it->second.get()).second)
	    payloads.push_back(it->second);
    }

    return payloads;
}

static std::vector<BObolViewLodState::ProxyPayloadPtr>
view_lod_unique_proxy_payloads(
    const std::unordered_map<std::string,
    BObolViewLodState::ProxyPayloadPtr> &bindings)
{
    std::vector<BObolViewLodState::ProxyPayloadPtr> payloads;
    std::unordered_set<const BObolViewLodState::ProxyPayload *> seen;
    seen.reserve(bindings.size());
    for (std::unordered_map<std::string,
	 BObolViewLodState::ProxyPayloadPtr>::const_iterator it =
	     bindings.begin(); it != bindings.end(); ++it) {
	if (!it->second)
	    continue;
	if (seen.insert(it->second.get()).second)
	    payloads.push_back(it->second);
    }

    return payloads;
}

size_t
BObolViewLodState::payloadCount(void) const
{
    return view_lod_unique_payloads(this->meshBindings).size() +
	   view_lod_unique_proxy_payloads(this->proxyBindings).size() +
	   this->cadPayloadCount();
}

size_t
BObolViewLodState::meshPayloadCount(void) const
{
    return view_lod_unique_payloads(this->meshBindings).size();
}

size_t
BObolViewLodState::proxyPayloadCount(int proxyKind) const
{
    std::vector<ProxyPayloadPtr> payloads =
	view_lod_unique_proxy_payloads(this->proxyBindings);
    size_t count = 0;
    for (size_t i = 0; i < payloads.size(); i++) {
	if (!payloads[i] || !payloads[i]->isValid())
	    continue;
	if (proxyKind == BOBOL_LOD_PROXY_NONE ||
	    payloads[i]->proxy.kind == proxyKind)
	    count++;
    }

    return count;
}

size_t
BObolViewLodState::cadPayloadCount(void) const
{
    size_t count = 0;
    for (const auto &source : this->cadSourceBindings)
	for (const auto &binding : source.second)
	    if (binding.second && binding.second->isValid())
		count++;
    return count;
}

size_t
BObolViewLodState::cadMeshPayloadCount(void) const
{
    size_t count = 0;
    for (const auto &source : this->cadSourceBindings)
	for (const auto &binding : source.second) {
	    const CadPayloadPtr &payload = binding.second;
	    if (payload && payload->isValid() &&
		(payload->resultKind == BOBOL_LOD_RESULT_MESH ||
		 payload->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL))
		count++;
	}
    return count;
}

size_t
BObolViewLodState::cadProxyPayloadCount(int proxyKind) const
{
    size_t count = 0;
    for (const auto &source : this->cadSourceBindings)
	for (const auto &binding : source.second) {
	    const CadPayloadPtr &payload = binding.second;
	    if (!payload || !payload->isValid() ||
		(payload->resultKind != BOBOL_LOD_RESULT_AABB &&
		 payload->resultKind != BOBOL_LOD_RESULT_PROXY) ||
		(proxyKind != BOBOL_LOD_PROXY_NONE &&
		 payload->proxy.kind != proxyKind))
		continue;
	    count++;
	}
    return count;
}

static size_t
view_lod_saturating_add(size_t lhs, size_t rhs)
{
    return rhs > SIZE_MAX - lhs ? SIZE_MAX : lhs + rhs;
}

size_t
BObolViewLodState::activeFaceCount(void) const
{
    size_t faces = 0;
    const std::vector<MeshPayloadPtr> meshPayloads =
	view_lod_unique_payloads(this->meshBindings);
    for (const MeshPayloadPtr &payload : meshPayloads) {
	if (!payload || !payload->isValid() ||
	    payload->resultKind != BOBOL_LOD_RESULT_MESH)
	    continue;
	faces = view_lod_saturating_add(faces, payload->counts.faceCount);
    }

    for (const auto &source : this->cadSourceBindings)
	for (const auto &binding : source.second) {
	    const CadPayloadPtr &payload = binding.second;
	    if (!payload || !payload->isValid() ||
		(payload->resultKind != BOBOL_LOD_RESULT_MESH &&
		 payload->resultKind != BOBOL_LOD_RESULT_FULL_DETAIL))
		continue;
	    faces = view_lod_saturating_add(faces,
		payload->counts.faceCount);
	}
    return faces;
}

int
BObolViewLodState::maximumActiveProgressiveLevel(void) const
{
    int maximum = -1;
    for (const auto &source : this->cadSourceBindings)
	for (const auto &binding : source.second) {
	    const CadPayloadPtr &payload = binding.second;
	    if (!payload || !payload->isValid() ||
		!payload->progressiveMesh || payload->activeLevel < 0)
		continue;
	    maximum = std::max(maximum, payload->activeLevel);
	}
    return maximum;
}

void
BObolViewLodState::setCadPresentationProgressiveLodCeiling(int level) const
{
    level = level < 0 ? -1 : std::min(15, level);
    this->cadPresentationProgressiveLodCeiling = level;
    for (const auto &presentation : this->cadPresentations)
	if (presentation.second.assembly)
	    presentation.second.assembly->progressiveLodCeiling = level;
}

size_t
BObolViewLodState::estimateDisplayMeshBytes(void) const
{
    std::vector<MeshPayloadPtr> payloads =
	view_lod_unique_payloads(this->meshBindings);
    size_t bytes = 0;
    for (size_t i = 0; i < payloads.size(); i++)
	bytes += payloads[i]->estimateBytes();
    std::vector<ProxyPayloadPtr> proxyPayloads =
	view_lod_unique_proxy_payloads(this->proxyBindings);
    for (size_t i = 0; i < proxyPayloads.size(); i++)
	bytes += proxyPayloads[i]->estimateBytes();
    for (const auto &source : this->cadSourceBindings)
	for (const auto &binding : source.second)
	    if (binding.second && binding.second->isValid())
		bytes += binding.second->estimateBytes();

    return bytes;
}

void
BObolViewLodState::residentMeshDemands(
    std::vector<BObolLodResidentDemand> &demands) const
{
    demands.clear();
    std::unordered_map<std::string, int> maximumLevels;

    std::vector<MeshPayloadPtr> meshPayloads =
	view_lod_unique_payloads(this->meshBindings);
    for (const MeshPayloadPtr &payload : meshPayloads) {
	if (!payload || !payload->progressiveMesh ||
	    !payload->progressiveMesh->isValid() ||
	    payload->cacheKey.getLength() == 0 ||
	    payload->requestedLevel < 0)
	    continue;
	int &level = maximumLevels[payload->cacheKey.getString()];
	level = std::max(level, payload->requestedLevel);
    }

    for (const auto &source : this->cadSourceBindings)
	for (const auto &binding : source.second) {
	    const CadPayloadPtr &payload = binding.second;
	    if (!payload || !payload->progressiveMesh ||
		!payload->progressiveMesh->isValid() ||
		payload->cacheKey.getLength() == 0 ||
		payload->requestedLevel < 0)
		continue;
	    int &level = maximumLevels[payload->cacheKey.getString()];
	    level = std::max(level, payload->requestedLevel);
	}

    demands.reserve(maximumLevels.size());
    for (const auto &entry : maximumLevels) {
	BObolLodResidentDemand demand;
	demand.assetKey = entry.first.c_str();
	demand.level = entry.second;
	demands.push_back(demand);
    }
    std::sort(demands.begin(), demands.end(),
	[](const BObolLodResidentDemand &a,
	   const BObolLodResidentDemand &b) {
	    return bu_strcmp(a.assetKey.getString(),
		b.assetKey.getString()) < 0;
	});
}

uint64_t
BObolViewLodState::cadRevision(void) const
{
    return this->cadBindingsRevision;
}

void
BObolViewLodState::cadOccurrenceChangesSince(
    const SoBRLDatabaseSource *source, uint64_t revision,
    std::vector<SbString> &occurrenceKeys, SbBool &fullResync) const
{
    occurrenceKeys.clear();
    fullResync = revision < this->cadFullResyncRevision ? TRUE : FALSE;
    if (fullResync || !source)
	return;

    const std::string sourceKey = view_lod_source_primary_key(source);
    const auto changes = this->cadOccurrenceChanges.find(sourceKey);
    if (changes == this->cadOccurrenceChanges.end())
	return;

    std::unordered_set<std::string> seen;
    for (const CadOccurrenceChange &change : changes->second) {
	if (change.revision <= revision ||
	    change.occurrenceKey.getLength() == 0)
	    continue;
	if (seen.insert(change.occurrenceKey.getString()).second)
	    occurrenceKeys.push_back(change.occurrenceKey);
    }
}

void
BObolViewLodState::acknowledgeCadOccurrenceChanges(
    const SoBRLDatabaseSource *source, uint64_t revision) const
{
    if (!source)
	return;
    const std::string sourceKey = view_lod_source_primary_key(source);
    auto changes = this->cadOccurrenceChanges.find(sourceKey);
    if (changes == this->cadOccurrenceChanges.end())
	return;
    std::vector<CadOccurrenceChange> &entries = changes->second;
    entries.erase(std::remove_if(entries.begin(), entries.end(),
	[revision](const CadOccurrenceChange &change) {
	    return change.revision <= revision;
	}), entries.end());
    if (entries.empty())
	this->cadOccurrenceChanges.erase(changes);
}

void
BObolViewLodState::noteCadOccurrenceChanged(
    const std::string &sourceBindingKey, const SbString &occurrenceKey)
{
    if (sourceBindingKey.empty() || occurrenceKey.getLength() == 0) {
	this->noteResidentMeshesChanged();
	return;
    }
    this->cadBindingsRevision++;
    if (!this->cadBindingsRevision) {
	this->cadBindingsRevision = 1;
	this->cadFullResyncRevision = 1;
	this->cadOccurrenceChanges.clear();
	return;
    }
    CadOccurrenceChange change;
    change.revision = this->cadBindingsRevision;
    change.occurrenceKey = occurrenceKey;
    this->cadOccurrenceChanges[sourceBindingKey].push_back(
	std::move(change));
}

void
BObolViewLodState::noteResidentMeshesChanged(void)
{
    this->cadBindingsRevision++;
    if (!this->cadBindingsRevision)
	this->cadBindingsRevision++;
    this->cadFullResyncRevision = this->cadBindingsRevision;
    this->cadOccurrenceChanges.clear();
}

size_t
BObolViewLodState::evictDisplayMeshes(unsigned int *evictedMeshCount)
{
    std::vector<MeshPayloadPtr> payloads =
	view_lod_unique_payloads(this->meshBindings);
    size_t bytes = 0;
    for (size_t i = 0; i < payloads.size(); i++)
	bytes += payloads[i]->estimateBytes();
    std::vector<ProxyPayloadPtr> proxyPayloads =
	view_lod_unique_proxy_payloads(this->proxyBindings);
    for (size_t i = 0; i < proxyPayloads.size(); i++)
	bytes += proxyPayloads[i]->estimateBytes();
    size_t cadPayloadCount = 0;
    for (const auto &source : this->cadSourceBindings)
	for (const auto &binding : source.second) {
	    if (!binding.second || !binding.second->isValid())
		continue;
	    bytes += binding.second->estimateBytes();
	    cadPayloadCount++;
	}

    if (evictedMeshCount)
	*evictedMeshCount = static_cast<unsigned int>(
	    payloads.size() + proxyPayloads.size() + cadPayloadCount);
    this->meshBindings.clear();
    this->proxyBindings.clear();
    this->cadSourceBindings.clear();
    this->cadAssetBindings.clear();
    this->cadBindings.clear();
    this->noteResidentMeshesChanged();
    return bytes;
}

size_t
BObolViewLodState::evictDisplayMeshPayloads(
    unsigned int *evictedMeshCount)
{
    std::vector<MeshPayloadPtr> meshPayloads =
	view_lod_unique_payloads(this->meshBindings);
    std::vector<CadPayloadPtr> cadMeshPayloads;
    for (const auto &source : this->cadSourceBindings)
	for (const auto &binding : source.second) {
	    const CadPayloadPtr &payload = binding.second;
	    if (payload && payload->isValid() &&
		(payload->resultKind == BOBOL_LOD_RESULT_MESH ||
		 payload->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL))
		cadMeshPayloads.push_back(payload);
	}

    size_t bytes = 0;
    for (const MeshPayloadPtr &payload : meshPayloads)
	bytes += payload->estimateBytes();
    for (const CadPayloadPtr &payload : cadMeshPayloads)
	bytes += payload->estimateBytes();

    if (evictedMeshCount)
	*evictedMeshCount = static_cast<unsigned int>(
	    meshPayloads.size() + cadMeshPayloads.size());
    this->meshBindings.clear();
    this->cadAssetBindings.clear();
    for (auto source = this->cadSourceBindings.begin();
	 source != this->cadSourceBindings.end();) {
	for (auto binding = source->second.begin();
	     binding != source->second.end();) {
	    const CadPayloadPtr &payload = binding->second;
	    if (payload &&
		(payload->resultKind == BOBOL_LOD_RESULT_MESH ||
		 payload->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL))
		binding = source->second.erase(binding);
	    else
		++binding;
	}
	if (source->second.empty())
	    source = this->cadSourceBindings.erase(source);
	else
	    ++source;
    }
    for (auto binding = this->cadBindings.begin();
	 binding != this->cadBindings.end();) {
	const CadPayloadPtr &payload = binding->second;
	if (payload &&
	    (payload->resultKind == BOBOL_LOD_RESULT_MESH ||
	     payload->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL)) {
	    binding = this->cadBindings.erase(binding);
	    continue;
	}
	++binding;
    }
    if (!cadMeshPayloads.empty())
	this->noteResidentMeshesChanged();

    return bytes;
}

SoBRLViewLodElement::~SoBRLViewLodElement(void)
{
}

void
SoBRLViewLodElement::initClass(void)
{
    SO_ELEMENT_INIT_CLASS(SoBRLViewLodElement, SoElement);
}

void
SoBRLViewLodElement::init(SoState *state)
{
    inherited::init(state);
    this->viewState = NULL;
}

void
SoBRLViewLodElement::push(SoState *state)
{
    const SoBRLViewLodElement *prev =
	static_cast<const SoBRLViewLodElement *>(
	    this->getNextInStack());
    this->viewState = prev ? prev->viewState : NULL;
    inherited::push(state);
}

SbBool
SoBRLViewLodElement::matches(const SoElement *element) const
{
    const SoBRLViewLodElement *other =
	static_cast<const SoBRLViewLodElement *>(element);
    return other && this->viewState == other->viewState ? TRUE : FALSE;
}

SoElement *
SoBRLViewLodElement::copyMatchInfo(void) const
{
    SoBRLViewLodElement *element = new SoBRLViewLodElement;
    element->viewState = this->viewState;
    return element;
}

void
SoBRLViewLodElement::set(SoState *state,
			 SoNode *UNUSED(node),
			 const BObolViewLodState *newViewState)
{
    if (!state || !state->isElementEnabled(SoBRLViewLodElement::getClassStackIndex()))
	return;

    SoBRLViewLodElement *element =
	static_cast<SoBRLViewLodElement *>(
	    SoElement::getElement(state,
				  SoBRLViewLodElement::getClassStackIndex()));
    element->viewState = newViewState;
}

const BObolViewLodState *
SoBRLViewLodElement::get(SoState *state)
{
    if (!state || !state->isElementEnabled(SoBRLViewLodElement::getClassStackIndex()))
	return NULL;

    const SoBRLViewLodElement *element =
	static_cast<const SoBRLViewLodElement *>(
	    SoElement::getConstElement(state,
				       SoBRLViewLodElement::getClassStackIndex()));
    return element ? element->viewState : NULL;
}

SoBRLViewLodGroup::SoBRLViewLodGroup(void) :
    viewState(NULL),
    softwareWireMode(SoCADViewState::SOFTWARE_WIRE_AUTO)
{
    SO_NODE_CONSTRUCTOR(SoBRLViewLodGroup);
}

SoBRLViewLodGroup::~SoBRLViewLodGroup(void)
{
}

void
SoBRLViewLodGroup::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLViewLodGroup, SoGroup, "Group");
    SO_ENABLE(SoGLRenderAction, SoBRLViewLodElement);
    SO_ENABLE(SoGetBoundingBoxAction, SoBRLViewLodElement);
    SO_ENABLE(SoRayPickAction, SoBRLViewLodElement);
}

void
SoBRLViewLodGroup::setViewLodState(BObolViewLodState *newViewState)
{
    this->viewState = newViewState;
}

BObolViewLodState *
SoBRLViewLodGroup::getViewLodState(void) const
{
    return this->viewState;
}

void
SoBRLViewLodGroup::setSoftwareWireMode(int mode)
{
    if (mode < SoCADViewState::SOFTWARE_WIRE_AUTO ||
	mode > SoCADViewState::SOFTWARE_WIRE_FAST)
	mode = SoCADViewState::SOFTWARE_WIRE_AUTO;
    this->softwareWireMode = mode;
}

int
SoBRLViewLodGroup::getSoftwareWireMode(void) const
{
    return this->softwareWireMode;
}

SbBool
SoBRLViewLodGroup::pushViewState(SoAction *action)
{
    if (!action || !this->viewState)
	return FALSE;

    SoState *state = action->getState();
    if (!state || !state->isElementEnabled(SoBRLViewLodElement::getClassStackIndex()))
	return FALSE;

    state->push();
    SoBRLViewLodElement::set(state, this, this->viewState);
    if (state->isElementEnabled(
	    SoCADViewStateElement::getClassStackIndex())) {
	Obol::CadViewState cadState = SoCADViewStateElement::get(state);
	cadState.softwareWireMode =
	    static_cast<Obol::CadSoftwareWireMode>(this->softwareWireMode);
	SoCADViewStateElement::set(state, cadState);
    }
    return TRUE;
}

void
SoBRLViewLodGroup::popViewState(SoAction *action, SbBool pushed)
{
    if (!pushed || !action)
	return;

    SoState *state = action->getState();
    if (state)
	state->pop();
}

void
SoBRLViewLodGroup::doAction(SoAction *action)
{
    const SbBool pushed = this->pushViewState(action);
    inherited::doAction(action);
    this->popViewState(action, pushed);
}

void
SoBRLViewLodGroup::GLRender(SoGLRenderAction *action)
{
    const SbBool pushed = this->pushViewState(action);
    inherited::GLRender(action);
    this->popViewState(action, pushed);
}

void
SoBRLViewLodGroup::callback(SoCallbackAction *action)
{
    const SbBool pushed = this->pushViewState(action);
    inherited::callback(action);
    this->popViewState(action, pushed);
}

void
SoBRLViewLodGroup::getBoundingBox(SoGetBoundingBoxAction *action)
{
    const SbBool pushed = this->pushViewState(action);
    inherited::getBoundingBox(action);
    this->popViewState(action, pushed);
}

void
SoBRLViewLodGroup::pick(SoPickAction *action)
{
    const SbBool pushed = this->pushViewState(action);
    inherited::pick(action);
    this->popViewState(action, pushed);
}

const BObolViewLodState::MeshPayload *
bobol_view_lod_mesh_for_action(SoAction *action,
				 const SoBRLMeshShape *shape)
{
    if (!action || !shape)
	return NULL;

    const BObolViewLodState *viewState =
	SoBRLViewLodElement::get(action->getState());
    return viewState ? viewState->findMesh(shape) : NULL;
}

const BObolViewLodState *
bobol_view_lod_state_for_action(SoAction *action)
{
    return action ? SoBRLViewLodElement::get(action->getState()) : NULL;
}

const BObolViewLodState::ProxyPayload *
bobol_view_lod_proxy_for_action(SoAction *action,
				  const SoBRLMeshShape *shape)
{
    if (!action || !shape)
	return NULL;

    const BObolViewLodState *viewState =
	SoBRLViewLodElement::get(action->getState());
    return viewState ? viewState->findProxy(shape) : NULL;
}

const BObolViewLodState::CadPayload *
bobol_view_lod_cad_for_action(SoAction *action,
				const SoBRLDatabaseSource *source)
{
    if (!action || !source)
	return NULL;

    const BObolViewLodState *viewState =
	SoBRLViewLodElement::get(action->getState());
    return viewState ? viewState->findCad(source) : NULL;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
