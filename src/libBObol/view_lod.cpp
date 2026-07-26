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

static SbBool
view_lod_cad_payload_can_replace(
	const std::unordered_map<std::string,
	BObolViewLodState::CadPayloadPtr> &bindings,
	const std::string &sourceBindingKey, const SbString &occurrenceKey,
	const BObolViewLodState::CadPayload &candidate)
{
    const int candidateRank = view_lod_cad_payload_rank(candidate);
    for (const auto &binding : bindings) {
	const BObolViewLodState::CadPayloadPtr &payload = binding.second;
	if (!payload || !view_lod_cad_payload_targets(*payload,
		sourceBindingKey, occurrenceKey))
	    continue;
	if (view_lod_cad_payload_rank(*payload) > candidateRank)
	    return FALSE;
    }
    return TRUE;
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
    cadBindingsRevision(1),
    normalStyle(NORMAL_AUTHORED),
    normalCreaseAngle(60.0f)
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
    this->cadBindings.clear();
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

    if (!view_lod_cad_payload_can_replace(this->cadBindings,
	    sourceBindingKey, payload->sourceInstanceKey, *payload))
	return TRUE;

    view_lod_remove_superseded_cad_payloads(this->cadBindings,
	sourceBindingKey, payload->sourceInstanceKey, *payload);

    std::vector<std::string> resultKeys;
    view_lod_result_keys(resultKeys, result);
    for (const std::string &resultKey : resultKeys) {
	std::string key = "source:";
	key += sourceBindingKey;
	key += '|';
	key += resultKey;
	this->cadBindings[key] = payload;
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
BObolViewLodState::findCadPayloads(
    const SoBRLDatabaseSource *source,
    std::vector<const CadPayload *> &payloads) const
{
    payloads.clear();
    if (!source || this->cadBindings.empty())
	return;

    const std::string sourceKey = view_lod_source_primary_key(source);
    for (const auto &binding : this->cadBindings) {
	const CadPayloadPtr &payload = binding.second;
	if (!payload || !payload->isValid() ||
	    bu_strcmp(payload->sourceBindingKey.getString(),
		sourceKey.c_str()) != 0)
	    continue;
	if (std::find(payloads.begin(), payloads.end(), payload.get()) ==
	    payloads.end())
	    payloads.push_back(payload.get());
    }
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
BObolViewLodState::findCadForResult(
    const BObolLodResult &result) const
{
    if (this->cadBindings.empty())
	return NULL;

    std::vector<const CadPayload *> seen;
    for (const auto &binding : this->cadBindings) {
	const CadPayloadPtr &payload = binding.second;
	if (!payload || !payload->isValid() ||
	    std::find(seen.begin(), seen.end(), payload.get()) != seen.end())
	    continue;
	seen.push_back(payload.get());
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

    return NULL;
}

SbBool
BObolViewLodState::retargetMeshPayload(
    const BObolViewLodState::MeshPayload *target,
    int activeLevel,
    uint64_t viewRevision,
    uint64_t policyRevision)
{
    if (!target || !target->progressiveMesh ||
	!target->progressiveMesh->canDrawLevel(activeLevel))
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
    payload->requestedLevel = activeLevel;
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
    uint64_t viewRevision,
    uint64_t policyRevision)
{
    if (!target || !target->progressiveMesh ||
	!target->progressiveMesh->canDrawLevel(activeLevel))
	return FALSE;

    CadPayloadPtr payload;
    for (const auto &binding : this->cadBindings) {
	if (binding.second && binding.second.get() == target) {
	    payload = binding.second;
	    break;
	}
    }
    if (!payload)
	return FALSE;

    payload->activeLevel = activeLevel;
    payload->requestedLevel = activeLevel;
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
    this->noteResidentMeshesChanged();
    return TRUE;
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
    for (std::unordered_map<std::string,
	 BObolViewLodState::MeshPayloadPtr>::const_iterator it =
	     bindings.begin(); it != bindings.end(); ++it) {
	if (!it->second)
	    continue;
	if (std::find(payloads.begin(), payloads.end(), it->second) ==
	    payloads.end())
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
    for (std::unordered_map<std::string,
	 BObolViewLodState::ProxyPayloadPtr>::const_iterator it =
	     bindings.begin(); it != bindings.end(); ++it) {
	if (!it->second)
	    continue;
	if (std::find(payloads.begin(), payloads.end(), it->second) ==
	    payloads.end())
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
    std::vector<CadPayloadPtr> payloads;
    for (std::unordered_map<std::string,
	 CadPayloadPtr>::const_iterator it =
	     this->cadBindings.begin(); it != this->cadBindings.end(); ++it) {
	if (!it->second || !it->second->isValid())
	    continue;
	if (std::find(payloads.begin(), payloads.end(), it->second) ==
	    payloads.end())
	    payloads.push_back(it->second);
    }
    return payloads.size();
}

size_t
BObolViewLodState::cadMeshPayloadCount(void) const
{
    std::vector<CadPayloadPtr> payloads;
    for (const auto &binding : this->cadBindings) {
	const CadPayloadPtr &payload = binding.second;
	if (!payload || !payload->isValid() ||
	    (payload->resultKind != BOBOL_LOD_RESULT_MESH &&
	     payload->resultKind != BOBOL_LOD_RESULT_FULL_DETAIL))
	    continue;
	if (std::find(payloads.begin(), payloads.end(), payload) ==
	    payloads.end())
	    payloads.push_back(payload);
    }
    return payloads.size();
}

size_t
BObolViewLodState::cadProxyPayloadCount(int proxyKind) const
{
    std::vector<CadPayloadPtr> payloads;
    for (const auto &binding : this->cadBindings) {
	const CadPayloadPtr &payload = binding.second;
	if (!payload || !payload->isValid() ||
	    (payload->resultKind != BOBOL_LOD_RESULT_AABB &&
	     payload->resultKind != BOBOL_LOD_RESULT_PROXY) ||
	    (proxyKind != BOBOL_LOD_PROXY_NONE &&
	     payload->proxy.kind != proxyKind))
	    continue;
	if (std::find(payloads.begin(), payloads.end(), payload) ==
	    payloads.end())
	    payloads.push_back(payload);
    }
    return payloads.size();
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
    std::vector<CadPayloadPtr> cadPayloads;
    for (std::unordered_map<std::string,
	 CadPayloadPtr>::const_iterator it =
	     this->cadBindings.begin(); it != this->cadBindings.end(); ++it) {
	if (!it->second || !it->second->isValid())
	    continue;
	if (std::find(cadPayloads.begin(), cadPayloads.end(), it->second) ==
	    cadPayloads.end())
	    cadPayloads.push_back(it->second);
    }
    for (size_t i = 0; i < cadPayloads.size(); i++)
	bytes += cadPayloads[i]->estimateBytes();

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
	    payload->cacheKey.getLength() == 0 || payload->activeLevel < 0)
	    continue;
	int &level = maximumLevels[payload->cacheKey.getString()];
	level = std::max(level, payload->activeLevel);
    }

    std::vector<CadPayloadPtr> cadPayloads;
    for (const auto &binding : this->cadBindings) {
	const CadPayloadPtr &payload = binding.second;
	if (!payload || !payload->progressiveMesh ||
	    !payload->progressiveMesh->isValid() ||
	    payload->cacheKey.getLength() == 0 || payload->activeLevel < 0 ||
	    std::find(cadPayloads.begin(), cadPayloads.end(), payload) !=
		cadPayloads.end())
	    continue;
	cadPayloads.push_back(payload);
	int &level = maximumLevels[payload->cacheKey.getString()];
	level = std::max(level, payload->activeLevel);
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
BObolViewLodState::noteResidentMeshesChanged(void)
{
    this->cadBindingsRevision++;
    if (!this->cadBindingsRevision)
	this->cadBindingsRevision++;
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
    std::vector<CadPayloadPtr> cadPayloads;
    for (std::unordered_map<std::string,
	 CadPayloadPtr>::const_iterator it =
	     this->cadBindings.begin(); it != this->cadBindings.end(); ++it) {
	if (!it->second || !it->second->isValid())
	    continue;
	if (std::find(cadPayloads.begin(), cadPayloads.end(), it->second) ==
	    cadPayloads.end())
	    cadPayloads.push_back(it->second);
    }
    for (size_t i = 0; i < cadPayloads.size(); i++)
	bytes += cadPayloads[i]->estimateBytes();

    if (evictedMeshCount)
	*evictedMeshCount = static_cast<unsigned int>(
	    payloads.size() + proxyPayloads.size() + cadPayloads.size());
    this->meshBindings.clear();
    this->proxyBindings.clear();
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
    for (const auto &binding : this->cadBindings) {
	const CadPayloadPtr &payload = binding.second;
	if (!payload || !payload->isValid() ||
	    (payload->resultKind != BOBOL_LOD_RESULT_MESH &&
	     payload->resultKind != BOBOL_LOD_RESULT_FULL_DETAIL) ||
	    std::find(cadMeshPayloads.begin(), cadMeshPayloads.end(), payload) !=
		cadMeshPayloads.end())
	    continue;
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
