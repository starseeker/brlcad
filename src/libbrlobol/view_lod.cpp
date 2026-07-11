/*                  V I E W _ L O D . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/export_action.h"
#include "brlobol/database_source.h"
#include "brlobol/measure_action.h"
#include "brlobol/mesh_lod_submit_action.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/snap_action.h"
#include "brlobol/view_lod.h"

#include <obol/cad/SoCADAssembly.h>

#include <Inventor/actions/SoAction.h>
#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/actions/SoPickAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/misc/SoState.h>
#include <Inventor/nodes/SoNode.h>

#include <algorithm>
#include <string.h>

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
		     const BRLObolLodResult &result)
{
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

BRLObolViewLodState::MeshPayload::MeshPayload(void) :
    resultKind(BRLOBOL_LOD_RESULT_NONE),
    qualityTier(BRLOBOL_LOD_QUALITY_METADATA),
    providerStatus(BRLOBOL_LOD_PROVIDER_UNKNOWN),
    activeLevel(-1),
    viewRevision(0),
    policyRevision(0),
    hasSnappedPoints(FALSE),
    hasNormals(FALSE)
{
}

SbBool
BRLObolViewLodState::MeshPayload::isValid(void) const
{
    return this->mesh.isValid();
}

size_t
BRLObolViewLodState::MeshPayload::estimateBytes(void) const
{
    return this->mesh.points.size() * sizeof(SbVec3f) +
	   this->mesh.normals.size() * sizeof(SbVec3f) +
	   this->mesh.coordIndex.size() * sizeof(int32_t) +
	   this->mesh.faceIndex.size() * sizeof(int32_t) +
	   this->mesh.vertexIndex.size() * sizeof(int32_t);
}

int
BRLObolViewLodState::MeshPayload::getTriangleCount(void) const
{
    return static_cast<int>(this->mesh.coordIndex.size() / 3);
}

SbBool
BRLObolViewLodState::MeshPayload::getTriangleVertexIndices(int triangleIndex,
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
BRLObolViewLodState::MeshPayload::getTriangle(int triangleIndex,
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

BRLObolViewLodState::ProxyPayload::ProxyPayload(void) :
    resultKind(BRLOBOL_LOD_RESULT_NONE),
    qualityTier(BRLOBOL_LOD_QUALITY_METADATA),
    providerStatus(BRLOBOL_LOD_PROVIDER_UNKNOWN),
    viewRevision(0),
    policyRevision(0),
    diagnostic("")
{
}

SbBool
BRLObolViewLodState::ProxyPayload::isValid(void) const
{
    return this->proxy.isValid();
}

size_t
BRLObolViewLodState::ProxyPayload::estimateBytes(void) const
{
    return sizeof(*this);
}

BRLObolViewLodState::CadPayload::CadPayload(void) :
    resultKind(BRLOBOL_LOD_RESULT_NONE),
    qualityTier(BRLOBOL_LOD_QUALITY_METADATA),
    providerStatus(BRLOBOL_LOD_PROVIDER_UNKNOWN),
    drawMode(BRLOBOL_LOD_DRAW_UNKNOWN),
    viewRevision(0),
    policyRevision(0),
    hasSnappedPoints(FALSE),
    hasNormals(FALSE),
    assembly(NULL),
    assemblyKey("")
{
}

BRLObolViewLodState::CadPayload::~CadPayload(void)
{
    this->clearAssembly();
}

void
BRLObolViewLodState::CadPayload::clearAssembly(void) const
{
    if (this->assembly)
	this->assembly->unref();
    this->assembly = NULL;
    this->assemblyKey = "";
}

SbBool
BRLObolViewLodState::CadPayload::isValid(void) const
{
    if (this->providerStatus != BRLOBOL_LOD_PROVIDER_READY)
	return FALSE;
    if (this->resultKind == BRLOBOL_LOD_RESULT_MESH ||
	this->resultKind == BRLOBOL_LOD_RESULT_FULL_DETAIL)
	return this->mesh.isValid();
    if (this->resultKind == BRLOBOL_LOD_RESULT_AABB ||
	this->resultKind == BRLOBOL_LOD_RESULT_PROXY)
	return this->proxy.isValid();
    return FALSE;
}

size_t
BRLObolViewLodState::CadPayload::estimateBytes(void) const
{
    return this->mesh.points.size() * sizeof(SbVec3f) +
	   this->mesh.normals.size() * sizeof(SbVec3f) +
	   this->mesh.coordIndex.size() * sizeof(int32_t) +
	   this->mesh.faceIndex.size() * sizeof(int32_t) +
	   this->mesh.vertexIndex.size() * sizeof(int32_t) +
	   sizeof(*this);
}

BRLObolViewLodState::BRLObolViewLodState(void)
{
}

BRLObolViewLodState::~BRLObolViewLodState(void)
{
}

void
BRLObolViewLodState::clear(void)
{
    this->meshBindings.clear();
    this->proxyBindings.clear();
    this->cadBindings.clear();
}

SbBool
BRLObolViewLodState::applyMeshResult(const SoBRLMeshShape *shape,
				     const BRLObolLodResult &result)
{
    if (!shape ||
	result.resultKind != BRLOBOL_LOD_RESULT_MESH ||
	result.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	!result.mesh.isValid())
	return FALSE;

    MeshPayloadPtr payload(new MeshPayload);
    payload->mesh = result.mesh;
    payload->sourcePath = shape->sourcePath.getValue();
    payload->sourceName = shape->sourceName.getValue();
    payload->sourceIdentity = shape->sourceIdentity.getValue();
    payload->cacheIdentity = shape->cacheIdentity.getValue();
    payload->cacheKey = result.cacheKey.value;
    payload->resultKind = result.resultKind;
    payload->qualityTier = result.qualityTier;
    payload->providerStatus = result.providerStatus;
    payload->activeLevel = result.geometry.activeLevel;
    payload->viewRevision = result.request.viewRevision;
    payload->policyRevision = result.request.policyRevision;
    payload->counts = result.counts;
    payload->bounds = result.bounds;
    payload->hasSnappedPoints = result.hasSnappedPoints;
    payload->hasNormals = result.hasNormals;
    payload->diagnostic = result.diagnostic;

    std::vector<std::string> keys;
    view_lod_shape_keys(keys, shape);
    view_lod_result_keys(keys, result);
    for (size_t i = 0; i < keys.size(); i++)
	this->meshBindings[keys[i]] = payload;

    return TRUE;
}

SbBool
BRLObolViewLodState::applyProxyResult(const SoBRLMeshShape *shape,
				      const BRLObolLodResult &result)
{
    if (!shape ||
	(result.resultKind != BRLOBOL_LOD_RESULT_AABB &&
	 result.resultKind != BRLOBOL_LOD_RESULT_PROXY) ||
	result.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	!result.proxy.isValid())
	return FALSE;

    ProxyPayloadPtr payload(new ProxyPayload);
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
BRLObolViewLodState::applyDisplayResult(const SoBRLMeshShape *shape,
					const BRLObolLodResult &result)
{
    if (result.resultKind == BRLOBOL_LOD_RESULT_MESH ||
	result.resultKind == BRLOBOL_LOD_RESULT_FULL_DETAIL)
	return this->applyMeshResult(shape, result);

    return this->applyProxyResult(shape, result);
}

SbBool
BRLObolViewLodState::applySourceResult(
    const SoBRLDatabaseSource *source,
    const BRLObolLodResult &result)
{
    if (!source || result.providerStatus != BRLOBOL_LOD_PROVIDER_READY)
	return FALSE;

    if ((result.resultKind == BRLOBOL_LOD_RESULT_MESH ||
	 result.resultKind == BRLOBOL_LOD_RESULT_FULL_DETAIL) &&
	!result.mesh.isValid())
	return FALSE;
    if ((result.resultKind == BRLOBOL_LOD_RESULT_AABB ||
	 result.resultKind == BRLOBOL_LOD_RESULT_PROXY) &&
	!result.proxy.isValid())
	return FALSE;
    if (result.resultKind != BRLOBOL_LOD_RESULT_MESH &&
	result.resultKind != BRLOBOL_LOD_RESULT_FULL_DETAIL &&
	result.resultKind != BRLOBOL_LOD_RESULT_AABB &&
	result.resultKind != BRLOBOL_LOD_RESULT_PROXY)
	return FALSE;

    CadPayloadPtr payload(new CadPayload);
    payload->mesh = result.mesh;
    payload->proxy = result.proxy;
    payload->sourcePath = source->path.getValue();
    payload->sourceName = view_lod_leaf_name_from_path(source->path.getValue());
    payload->sourceIdentity =
	source->realizationIdentity.getValue().getLength() > 0 ?
	source->realizationIdentity.getValue() : source->path.getValue();
    payload->sourceInstanceKey = source->instanceKey.getValue();
    payload->cacheIdentity = result.cacheKey.value;
    payload->cacheKey = result.cacheKey.value;
    payload->resultKind = result.resultKind;
    payload->qualityTier = result.qualityTier;
    payload->providerStatus = result.providerStatus;
    payload->drawMode = result.request.drawMode;
    payload->viewRevision = result.request.viewRevision;
    payload->policyRevision = result.request.policyRevision;
    payload->counts = result.counts;
    payload->bounds = result.bounds;
    payload->hasSnappedPoints = result.hasSnappedPoints;
    payload->hasNormals = result.hasNormals;
    payload->diagnostic = result.diagnostic;

    std::vector<std::string> keys;
    view_lod_source_keys(keys, source);
    view_lod_result_keys(keys, result);
    for (size_t i = 0; i < keys.size(); i++)
	this->cadBindings[keys[i]] = payload;

    return TRUE;
}

const BRLObolViewLodState::MeshPayload *
BRLObolViewLodState::findMesh(const SoBRLMeshShape *shape) const
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

const BRLObolViewLodState::MeshPayload *
BRLObolViewLodState::findMeshForResult(
    const BRLObolLodResult &result) const
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

const BRLObolViewLodState::ProxyPayload *
BRLObolViewLodState::findProxy(const SoBRLMeshShape *shape) const
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

const BRLObolViewLodState::ProxyPayload *
BRLObolViewLodState::findProxyForResult(
    const BRLObolLodResult &result) const
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

const BRLObolViewLodState::CadPayload *
BRLObolViewLodState::findCad(const SoBRLDatabaseSource *source) const
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

const BRLObolViewLodState::CadPayload *
BRLObolViewLodState::findCadForResult(
    const BRLObolLodResult &result) const
{
    if (this->cadBindings.empty())
	return NULL;

    std::vector<std::string> keys;
    view_lod_result_keys(keys, result);
    for (size_t i = 0; i < keys.size(); i++) {
	std::unordered_map<std::string, CadPayloadPtr>::const_iterator it =
	    this->cadBindings.find(keys[i]);
	if (it != this->cadBindings.end() && it->second &&
	    it->second->isValid())
	    return it->second.get();
    }

    return NULL;
}

size_t
BRLObolViewLodState::bindingCount(void) const
{
    return this->meshBindings.size() + this->proxyBindings.size() +
	   this->cadBindings.size();
}

static std::vector<BRLObolViewLodState::MeshPayloadPtr>
view_lod_unique_payloads(
    const std::unordered_map<std::string,
    BRLObolViewLodState::MeshPayloadPtr> &bindings)
{
    std::vector<BRLObolViewLodState::MeshPayloadPtr> payloads;
    for (std::unordered_map<std::string,
	 BRLObolViewLodState::MeshPayloadPtr>::const_iterator it =
	     bindings.begin(); it != bindings.end(); ++it) {
	if (!it->second)
	    continue;
	if (std::find(payloads.begin(), payloads.end(), it->second) ==
	    payloads.end())
	    payloads.push_back(it->second);
    }

    return payloads;
}

static std::vector<BRLObolViewLodState::ProxyPayloadPtr>
view_lod_unique_proxy_payloads(
    const std::unordered_map<std::string,
    BRLObolViewLodState::ProxyPayloadPtr> &bindings)
{
    std::vector<BRLObolViewLodState::ProxyPayloadPtr> payloads;
    for (std::unordered_map<std::string,
	 BRLObolViewLodState::ProxyPayloadPtr>::const_iterator it =
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
BRLObolViewLodState::payloadCount(void) const
{
    return view_lod_unique_payloads(this->meshBindings).size() +
	   view_lod_unique_proxy_payloads(this->proxyBindings).size() +
	   this->cadPayloadCount();
}

size_t
BRLObolViewLodState::meshPayloadCount(void) const
{
    return view_lod_unique_payloads(this->meshBindings).size();
}

size_t
BRLObolViewLodState::proxyPayloadCount(int proxyKind) const
{
    std::vector<ProxyPayloadPtr> payloads =
	view_lod_unique_proxy_payloads(this->proxyBindings);
    size_t count = 0;
    for (size_t i = 0; i < payloads.size(); i++) {
	if (!payloads[i] || !payloads[i]->isValid())
	    continue;
	if (proxyKind == BRLOBOL_LOD_PROXY_NONE ||
	    payloads[i]->proxy.kind == proxyKind)
	    count++;
    }

    return count;
}

size_t
BRLObolViewLodState::cadPayloadCount(void) const
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
BRLObolViewLodState::cadMeshPayloadCount(void) const
{
    std::vector<CadPayloadPtr> payloads;
    for (const auto &binding : this->cadBindings) {
	const CadPayloadPtr &payload = binding.second;
	if (!payload || !payload->isValid() ||
	    (payload->resultKind != BRLOBOL_LOD_RESULT_MESH &&
	     payload->resultKind != BRLOBOL_LOD_RESULT_FULL_DETAIL))
	    continue;
	if (std::find(payloads.begin(), payloads.end(), payload) ==
	    payloads.end())
	    payloads.push_back(payload);
    }
    return payloads.size();
}

size_t
BRLObolViewLodState::cadProxyPayloadCount(int proxyKind) const
{
    std::vector<CadPayloadPtr> payloads;
    for (const auto &binding : this->cadBindings) {
	const CadPayloadPtr &payload = binding.second;
	if (!payload || !payload->isValid() ||
	    (payload->resultKind != BRLOBOL_LOD_RESULT_AABB &&
	     payload->resultKind != BRLOBOL_LOD_RESULT_PROXY) ||
	    (proxyKind != BRLOBOL_LOD_PROXY_NONE &&
	     payload->proxy.kind != proxyKind))
	    continue;
	if (std::find(payloads.begin(), payloads.end(), payload) ==
	    payloads.end())
	    payloads.push_back(payload);
    }
    return payloads.size();
}

size_t
BRLObolViewLodState::estimateDisplayMeshBytes(void) const
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

size_t
BRLObolViewLodState::evictDisplayMeshes(unsigned int *evictedMeshCount)
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
			 const BRLObolViewLodState *newViewState)
{
    if (!state || !state->isElementEnabled(SoBRLViewLodElement::getClassStackIndex()))
	return;

    SoBRLViewLodElement *element =
	static_cast<SoBRLViewLodElement *>(
	    SoElement::getElement(state,
				  SoBRLViewLodElement::getClassStackIndex()));
    element->viewState = newViewState;
}

const BRLObolViewLodState *
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
    viewState(NULL)
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
SoBRLViewLodGroup::setViewLodState(BRLObolViewLodState *newViewState)
{
    this->viewState = newViewState;
}

BRLObolViewLodState *
SoBRLViewLodGroup::getViewLodState(void) const
{
    return this->viewState;
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

const BRLObolViewLodState::MeshPayload *
brlobol_view_lod_mesh_for_action(SoAction *action,
				 const SoBRLMeshShape *shape)
{
    if (!action || !shape)
	return NULL;

    const BRLObolViewLodState *viewState =
	SoBRLViewLodElement::get(action->getState());
    return viewState ? viewState->findMesh(shape) : NULL;
}

const BRLObolViewLodState::ProxyPayload *
brlobol_view_lod_proxy_for_action(SoAction *action,
				  const SoBRLMeshShape *shape)
{
    if (!action || !shape)
	return NULL;

    const BRLObolViewLodState *viewState =
	SoBRLViewLodElement::get(action->getState());
    return viewState ? viewState->findProxy(shape) : NULL;
}

const BRLObolViewLodState::CadPayload *
brlobol_view_lod_cad_for_action(SoAction *action,
				const SoBRLDatabaseSource *source)
{
    if (!action || !source)
	return NULL;

    const BRLObolViewLodState *viewState =
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
