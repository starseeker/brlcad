/*        M E S H _ L O D _ S U B M I T _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/database_source.h"
#include "brlobol/lod_service.h"
#include "brlobol/mesh_lod_submit_action.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/view_lod.h"

#include "raytrace.h"

#include <Inventor/SbString.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>

#include <errno.h>
#include <stdlib.h>
#include <string.h>

SO_ACTION_SOURCE(SoBRLMeshLodSubmitAction);

SoBRLMeshLodSubmitAction::SoBRLMeshLodSubmitAction(void) :
    service(NULL),
    dbip(NULL),
    databaseId(""),
    databaseRevision(0),
    generation(0),
    viewRevision(0),
    policyRevision(0),
    providerId("brlobol_mesh_lod"),
    providerVersion("brlobol-cache-v1"),
    qualityTier(BRLOBOL_LOD_QUALITY_FAST_DISPLAY),
    refreshMissing(TRUE),
    reset(0),
    useForcedLevel(FALSE),
    forcedLevel(0),
    requireLodBacked(TRUE),
    submitAabbProxyStage(FALSE),
    submitObbProxyStage(FALSE),
    viewState(NULL),
    visitedMeshCount(0),
    submittedTaskCount(0),
    skippedMeshCount(0),
    diagnostics("")
{
    SO_ACTION_CONSTRUCTOR(SoBRLMeshLodSubmitAction);
    bv_view_info_init(&this->view);
}

SoBRLMeshLodSubmitAction::~SoBRLMeshLodSubmitAction(void)
{
}

void
SoBRLMeshLodSubmitAction::initClass(void)
{
    SO_ACTION_INIT_CLASS(SoBRLMeshLodSubmitAction, SoAction);
    SO_ACTION_ADD_METHOD(SoNode, SoBRLMeshLodSubmitAction::nodeAction);
    SO_ACTION_ADD_METHOD(SoGroup, SoBRLMeshLodSubmitAction::nodeAction);
    SO_ACTION_ADD_METHOD(SoBRLDatabaseSource,
			 SoBRLMeshLodSubmitAction::databaseSourceAction);
    SO_ACTION_ADD_METHOD(SoBRLMeshShape,
			 SoBRLMeshLodSubmitAction::meshShapeAction);
}

void
SoBRLMeshLodSubmitAction::setService(BRLObolLodService *newService)
{
    this->service = newService;
}

BRLObolLodService *
SoBRLMeshLodSubmitAction::getService(void) const
{
    return this->service;
}

void
SoBRLMeshLodSubmitAction::setDatabase(struct db_i *newDbip,
				      const char *newDatabaseId, uint64_t newDatabaseRevision)
{
    this->dbip = newDbip;
    this->databaseId = newDatabaseId ? newDatabaseId : "";
    this->databaseRevision = newDatabaseRevision;
}

void
SoBRLMeshLodSubmitAction::setViewInfo(const struct bv_view_info *info)
{
    if (info)
	this->view = *info;
    else
	bv_view_info_init(&this->view);
    bv_view_info_sanitize(&this->view);
}

const struct bv_view_info &
SoBRLMeshLodSubmitAction::getViewInfo(void) const {
    return this->view;
}

void
SoBRLMeshLodSubmitAction::setGeneration(uint64_t newGeneration)
{
    this->generation = newGeneration;
}

uint64_t
SoBRLMeshLodSubmitAction::getGeneration(void) const
{
    return this->generation;
}

void
SoBRLMeshLodSubmitAction::setRevisions(uint64_t newViewRevision,
				       uint64_t newPolicyRevision)
{
    this->viewRevision = newViewRevision;
    this->policyRevision = newPolicyRevision;
}

void
SoBRLMeshLodSubmitAction::setProvider(const char *newProviderId,
				      const char *newProviderVersion)
{
    this->providerId = newProviderId ? newProviderId : "";
    this->providerVersion = newProviderVersion ? newProviderVersion : "";
}

void
SoBRLMeshLodSubmitAction::setQualityTier(int newQualityTier)
{
    this->qualityTier = newQualityTier;
}

void
SoBRLMeshLodSubmitAction::setRefreshMissing(SbBool newRefreshMissing)
{
    this->refreshMissing = newRefreshMissing ? TRUE : FALSE;
}

void
SoBRLMeshLodSubmitAction::setReset(int newReset)
{
    this->reset = newReset;
}

void
SoBRLMeshLodSubmitAction::setForcedLevel(int level)
{
    this->useForcedLevel = TRUE;
    this->forcedLevel = level < 0 ? 0 : level;
}

void
SoBRLMeshLodSubmitAction::clearForcedLevel(void)
{
    this->useForcedLevel = FALSE;
    this->forcedLevel = 0;
}

SbBool
SoBRLMeshLodSubmitAction::hasForcedLevel(void) const
{
    return this->useForcedLevel;
}

int
SoBRLMeshLodSubmitAction::getForcedLevel(void) const
{
    return this->forcedLevel;
}

void
SoBRLMeshLodSubmitAction::setRequireLodBacked(SbBool newRequireLodBacked)
{
    this->requireLodBacked = newRequireLodBacked ? TRUE : FALSE;
}

SbBool
SoBRLMeshLodSubmitAction::getRequireLodBacked(void) const
{
    return this->requireLodBacked;
}

void
SoBRLMeshLodSubmitAction::setProxyStages(SbBool submitAabb,
	SbBool submitObb)
{
    this->submitAabbProxyStage = submitAabb ? TRUE : FALSE;
    this->submitObbProxyStage = submitObb ? TRUE : FALSE;
}

SbBool
SoBRLMeshLodSubmitAction::getSubmitAabbProxyStage(void) const
{
    return this->submitAabbProxyStage;
}

SbBool
SoBRLMeshLodSubmitAction::getSubmitObbProxyStage(void) const
{
    return this->submitObbProxyStage;
}

void
SoBRLMeshLodSubmitAction::setViewLodState(
    const BRLObolViewLodState *newViewState)
{
    this->viewState = newViewState;
}

const BRLObolViewLodState *
SoBRLMeshLodSubmitAction::getViewLodState(void) const
{
    return this->viewState;
}

unsigned int
SoBRLMeshLodSubmitAction::getVisitedMeshCount(void) const
{
    return this->visitedMeshCount;
}

unsigned int
SoBRLMeshLodSubmitAction::getSubmittedTaskCount(void) const
{
    return this->submittedTaskCount;
}

unsigned int
SoBRLMeshLodSubmitAction::getSkippedMeshCount(void) const
{
    return this->skippedMeshCount;
}

const SbString &
SoBRLMeshLodSubmitAction::getDiagnostics(void) const
{
    return this->diagnostics;
}

void
SoBRLMeshLodSubmitAction::beginTraversal(SoNode *node)
{
    this->visitedMeshCount = 0;
    this->submittedTaskCount = 0;
    this->skippedMeshCount = 0;
    this->diagnostics = "";
    this->traverse(node);
}

void
SoBRLMeshLodSubmitAction::nodeAction(SoAction *action, SoNode *node)
{
    if (node->isOfType(SoGroup::getClassTypeId()))
	node->doAction(action);
}

static const char *
mesh_lod_source_leaf_name(const SoBRLDatabaseSource *source)
{
    if (!source)
	return "";
    const char *path = source->path.getValue().getString();
    if (!path || !path[0])
	return "";
    const char *slash = strrchr(path, '/');
    return (slash && slash[1]) ? slash + 1 : path;
}

static int
mesh_lod_source_is_threshold_bot(const SoBRLDatabaseSource *source)
{
    if (!source || source->lodBotThreshold.getValue() == 0)
	return 0;

    struct db_i *dbip = source->getDatabase();
    const char *name = mesh_lod_source_leaf_name(source);
    if (!dbip || !name || !name[0])
	return 0;

    struct directory *dp = db_lookup(dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL ||
	dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BOT)
	return 0;

    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, dbip, NULL) < 0)
	return 0;

    int ret = 0;
    if (intern.idb_type == ID_BOT && intern.idb_ptr) {
	const struct rt_bot_internal *bot =
	    static_cast<const struct rt_bot_internal *>(intern.idb_ptr);
	if (bot && bot->num_faces >=
	    static_cast<size_t>(source->lodBotThreshold.getValue()))
	    ret = 1;
    }
    rt_db_free_internal(&intern);
    return ret;
}

static void
mesh_lod_make_source_request(const SoBRLDatabaseSource *source,
			     BRLObolLodRequest &request,
			     const char *databaseId,
			     uint64_t databaseRevision,
			     uint64_t viewRevision,
			     uint64_t policyRevision,
			     int requestDrawMode,
			     const char *providerId,
			     const char *providerVersion,
			     int qualityTier)
{
    request.clear();
    request.databaseId = databaseId ? databaseId : "";
    request.databaseRevision = databaseRevision;
    request.sourceRevision = source ? source->sourceRevision.getValue() : 0;
    request.objectPath = source ? source->path.getValue() : SbString("");
    request.objectName = mesh_lod_source_leaf_name(source);
    request.viewRevision = viewRevision;
    request.policyRevision = policyRevision;
    request.drawMode = requestDrawMode;
    request.lodPolicy = 0;
    request.providerId = providerId ? providerId : "";
    request.providerVersion = providerVersion ? providerVersion : "";
    request.qualityTier = qualityTier;
    if (source) {
	SbBox3f bounds;
	if (source->getSourceBounds(bounds))
	    request.bounds = bounds;
    }
}

static uint32_t
mesh_lod_debug_delay_milliseconds(const char *env_name)
{
    const char *delay = getenv(env_name);
    if (!delay || !delay[0])
	return 0;

    errno = 0;
    char *end = NULL;
    unsigned long value = strtoul(delay, &end, 10);
    if (errno != 0 || end == delay || value > UINT32_MAX)
	return 0;

    return (uint32_t)value;
}

static int
mesh_lod_submit_proxy_task(BRLObolLodService *service,
			   struct db_i *dbip,
			   uint64_t generation,
			   const BRLObolLodRequest &request,
			   int proxyKind,
			   uint64_t dependencyTaskId,
			   const char *delayEnv,
			   unsigned int *submittedTaskCount,
			   uint64_t *taskIdOut)
{
    if (taskIdOut)
	*taskIdOut = 0;
    if (!service || !service->isRunning())
	return 0;
    if (service->hasActiveRequest(request))
	return 0;

    BRLObolRtProxyProvider *provider = new BRLObolRtProxyProvider;
    provider->dbip = dbip;
    provider->proxyKind = proxyKind;
    provider->useRequestBounds = TRUE;

    BRLObolLodTask task;
    task.generation = generation;
    task.request = request;
    if (dependencyTaskId != 0)
	task.addDependency(dependencyTaskId);
    task.realize = brlobol_rt_proxy_provider_task;
    task.realizeData = provider;
    task.realizeDataFree = brlobol_rt_proxy_provider_free;
    task.debugDelayMilliseconds =
	mesh_lod_debug_delay_milliseconds(delayEnv);

    uint64_t taskId = service->submitIfNotActive(task);
    if (taskId == 0) {
	brlobol_rt_proxy_provider_free(provider);
	return 0;
    }

    if (taskIdOut)
	*taskIdOut = taskId;
    if (submittedTaskCount)
	(*submittedTaskCount)++;
    return 1;
}

void
SoBRLMeshLodSubmitAction::databaseSourceAction(SoAction *action, SoNode *node)
{
    SoBRLMeshLodSubmitAction *submitAction =
	static_cast<SoBRLMeshLodSubmitAction *>(action);
    SoBRLDatabaseSource *source = static_cast<SoBRLDatabaseSource *>(node);

    if (!source->hasCompactInstanceIndex()) {
	source->doAction(action);
	return;
    }

    submitAction->visitedMeshCount++;
    SbString target = source->path.getValue();

    if (!submitAction->service || !submitAction->service->isRunning()) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target, "LoD service is not running");
	source->doAction(action);
	return;
    }

    if (!submitAction->dbip) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target, "LoD submit action has no database");
	source->doAction(action);
	return;
    }

    const int sourceDrawMode =
	source->drawMode.getValue() == SoBRLDatabaseSource::SHADED ?
	BRLOBOL_LOD_DRAW_SHADED : BRLOBOL_LOD_DRAW_WIRE;
    const int roleFlags = source->realizationRoleFlags.getValue();
    const int meshEligible =
	(roleFlags & SoBRLDatabaseSource::REALIZATION_ROLE_MESH) ||
	source->drawMode.getValue() == SoBRLDatabaseSource::SHADED ||
	source->representationMode.getValue() ==
	SoBRLDatabaseSource::REPRESENTATION_SHADED_BOTS ||
	source->representationMode.getValue() ==
	SoBRLDatabaseSource::REPRESENTATION_SHADED ||
	mesh_lod_source_is_threshold_bot(source);
    SbBox3f sourceBounds;
    const SbBool hasSourceBounds = source->getSourceBounds(sourceBounds);

    uint64_t dependencyTaskId = 0;
    if (hasSourceBounds &&
	(submitAction->submitAabbProxyStage ||
	 submitAction->submitObbProxyStage)) {
	BRLObolLodRequest aabbRequest;
	mesh_lod_make_source_request(source, aabbRequest,
	    submitAction->databaseId.getString(),
	    submitAction->databaseRevision,
	    submitAction->viewRevision,
	    submitAction->policyRevision,
	    sourceDrawMode, "rt_proxy_aabb", "rt-proxy-v1",
	    BRLOBOL_LOD_QUALITY_PROXY);

	if (submitAction->submitAabbProxyStage)
	    (void)mesh_lod_submit_proxy_task(submitAction->service,
		submitAction->dbip, submitAction->generation, aabbRequest,
		BRLOBOL_LOD_PROXY_AABB, 0,
		"BRLOBOL_LOD_AABB_TASK_DELAY_MS",
		&submitAction->submittedTaskCount, &dependencyTaskId);

	if (submitAction->submitObbProxyStage) {
	    BRLObolLodRequest obbRequest;
	    mesh_lod_make_source_request(source, obbRequest,
		submitAction->databaseId.getString(),
		submitAction->databaseRevision,
		submitAction->viewRevision,
		submitAction->policyRevision,
		sourceDrawMode, "rt_proxy_obb", "rt-proxy-v1",
		BRLOBOL_LOD_QUALITY_PROXY);
	    uint64_t obbTaskId = 0;
	    (void)mesh_lod_submit_proxy_task(submitAction->service,
		submitAction->dbip, submitAction->generation, obbRequest,
		BRLOBOL_LOD_PROXY_OBB, dependencyTaskId,
		"BRLOBOL_LOD_OBB_TASK_DELAY_MS",
		&submitAction->submittedTaskCount, &obbTaskId);
	    if (obbTaskId != 0)
		dependencyTaskId = obbTaskId;
	}
    } else if (submitAction->submitAabbProxyStage ||
	       submitAction->submitObbProxyStage) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target,
				       "compact source has no leaf-derived bounds for proxy LoD");
    }

    if (!meshEligible) {
	source->doAction(action);
	return;
    }
    if (submitAction->requireLodBacked &&
	source->lodBotThreshold.getValue() == 0) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target,
				       "compact source is not LoD-backed");
	source->doAction(action);
	return;
    }

    BRLObolLodRequest request;
    mesh_lod_make_source_request(source, request,
	submitAction->databaseId.getString(),
	submitAction->databaseRevision,
	submitAction->viewRevision,
	submitAction->policyRevision,
	sourceDrawMode,
	submitAction->providerId.getString(),
	submitAction->providerVersion.getString(),
	submitAction->qualityTier);

    const SbBool suppressActiveDuplicate =
	(!submitAction->useForcedLevel && submitAction->reset == 0) ?
	TRUE : FALSE;
    if (suppressActiveDuplicate &&
	submitAction->service->hasActiveRequest(request)) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target,
				       "current compact LoD request is already active");
	source->doAction(action);
	return;
    }

    BRLObolMeshLodProvider *provider = new BRLObolMeshLodProvider;
    provider->dbip = submitAction->dbip;
    provider->view = submitAction->view;
    provider->useView = TRUE;
    provider->refreshMissing = submitAction->refreshMissing;
    provider->useForcedLevel = submitAction->useForcedLevel;
    provider->forcedLevel = submitAction->forcedLevel;
    provider->shrinkAfterCopy = TRUE;
    provider->reset = submitAction->reset;

    BRLObolLodTask task;
    task.generation = submitAction->generation;
    task.request = request;
    task.realize = brlobol_mesh_lod_provider_task;
    task.realizeData = provider;
    task.realizeDataFree = brlobol_mesh_lod_provider_free;
    if (dependencyTaskId != 0)
	task.addDependency(dependencyTaskId);
    task.debugDelayMilliseconds =
	mesh_lod_debug_delay_milliseconds("BRLOBOL_LOD_TASK_DELAY_MS");

    uint64_t taskId = suppressActiveDuplicate ?
		      submitAction->service->submitIfNotActive(task) :
		      submitAction->service->submit(task);
    if (taskId == 0) {
	brlobol_mesh_lod_provider_free(provider);
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target,
				       "LoD service rejected compact mesh task");
    } else {
	submitAction->submittedTaskCount++;
    }

    source->doAction(action);
}

void
SoBRLMeshLodSubmitAction::meshShapeAction(SoAction *action, SoNode *node)
{
    SoBRLMeshLodSubmitAction *submitAction =
	static_cast<SoBRLMeshLodSubmitAction *>(action);
    SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);

    submitAction->visitedMeshCount++;

    SbString target = shape->sourcePath.getValue().getLength() > 0 ?
		      shape->sourcePath.getValue() : shape->sourceName.getValue();

    if (submitAction->requireLodBacked && !shape->isLodBackedMesh()) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target, "mesh is not LoD-backed");
	return;
    }

    if (!submitAction->service || !submitAction->service->isRunning()) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target, "LoD service is not running");
	return;
    }

    if (!submitAction->dbip) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target, "LoD submit action has no database");
	return;
    }

    BRLObolLodRequest request;
    shape->makeLodRequest(request,
			  submitAction->databaseId.getString(),
			  submitAction->databaseRevision,
			  submitAction->viewRevision,
			  submitAction->policyRevision,
			  BRLOBOL_LOD_DRAW_SHADED,
			  submitAction->providerId.getString(),
			  submitAction->providerVersion.getString(),
			  submitAction->qualityTier);

    BRLObolLodCacheKey requestKey = brlobol_lod_cache_key(request);
    const BRLObolViewLodState::MeshPayload *viewPayload =
	submitAction->viewState ?
	submitAction->viewState->findMesh(shape) : NULL;
    const SbBool viewPayloadResident =
	(viewPayload &&
	 viewPayload->resultKind == BRLOBOL_LOD_RESULT_MESH &&
	 viewPayload->providerStatus == BRLOBOL_LOD_PROVIDER_READY &&
	 requestKey.isValid() &&
	 viewPayload->cacheKey.getLength() > 0 &&
	 strcmp(viewPayload->cacheKey.getString(),
		requestKey.value.getString()) == 0) ? TRUE : FALSE;
    const SbBool shapePayloadResident =
	(shape->lodAvailable.getValue() &&
	 shape->lodResultKind.getValue() == BRLOBOL_LOD_RESULT_MESH &&
	 shape->lodProviderStatus.getValue() == BRLOBOL_LOD_PROVIDER_READY &&
	 requestKey.isValid() &&
	 shape->lodCacheKey.getValue().getLength() > 0 &&
	 strcmp(shape->lodCacheKey.getValue().getString(),
		requestKey.value.getString()) == 0) ? TRUE : FALSE;
    if (!submitAction->useForcedLevel &&
	submitAction->reset == 0 &&
	(viewPayloadResident || shapePayloadResident)) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target, "current LoD request is already resident");
	return;
    }

    const SbBool suppressActiveDuplicate =
	(!submitAction->useForcedLevel && submitAction->reset == 0) ?
	TRUE : FALSE;
    if (suppressActiveDuplicate &&
	submitAction->service->hasActiveRequest(request)) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target,
				       "current LoD request is already active");
	return;
    }

    uint64_t dependencyTaskId = 0;
    if (submitAction->submitAabbProxyStage ||
	submitAction->submitObbProxyStage) {
	BRLObolLodRequest aabbRequest;
	shape->makeLodRequest(aabbRequest,
			      submitAction->databaseId.getString(),
			      submitAction->databaseRevision,
			      submitAction->viewRevision,
			      submitAction->policyRevision,
			      BRLOBOL_LOD_DRAW_SHADED,
			      "rt_proxy_aabb",
			      "rt-proxy-v1",
			      BRLOBOL_LOD_QUALITY_PROXY);

	if (submitAction->submitAabbProxyStage &&
	    !submitAction->service->hasActiveRequest(aabbRequest)) {
	    BRLObolRtProxyProvider *provider = new BRLObolRtProxyProvider;
	    provider->dbip = submitAction->dbip;
	    provider->proxyKind = BRLOBOL_LOD_PROXY_AABB;
	    provider->useRequestBounds = TRUE;

	    BRLObolLodTask task;
	    task.generation = submitAction->generation;
	    task.request = aabbRequest;
	    task.realize = brlobol_rt_proxy_provider_task;
	    task.realizeData = provider;
	    task.realizeDataFree = brlobol_rt_proxy_provider_free;
	    task.debugDelayMilliseconds =
		mesh_lod_debug_delay_milliseconds(
		    "BRLOBOL_LOD_AABB_TASK_DELAY_MS");
	    uint64_t taskId = submitAction->service->submitIfNotActive(task);
	    if (taskId != 0) {
		dependencyTaskId = taskId;
		submitAction->submittedTaskCount++;
	    } else {
		brlobol_rt_proxy_provider_free(provider);
	    }
	}

	if (submitAction->submitObbProxyStage) {
	    BRLObolLodRequest obbRequest;
	    shape->makeLodRequest(obbRequest,
				  submitAction->databaseId.getString(),
				  submitAction->databaseRevision,
				  submitAction->viewRevision,
				  submitAction->policyRevision,
				  BRLOBOL_LOD_DRAW_SHADED,
				  "rt_proxy_obb",
				  "rt-proxy-v1",
				  BRLOBOL_LOD_QUALITY_PROXY);

	    if (!submitAction->service->hasActiveRequest(obbRequest)) {
		BRLObolRtProxyProvider *provider = new BRLObolRtProxyProvider;
		provider->dbip = submitAction->dbip;
		provider->proxyKind = BRLOBOL_LOD_PROXY_OBB;
		provider->useRequestBounds = TRUE;

		BRLObolLodTask task;
		task.generation = submitAction->generation;
		task.request = obbRequest;
		if (dependencyTaskId != 0)
		    task.addDependency(dependencyTaskId);
		task.realize = brlobol_rt_proxy_provider_task;
		task.realizeData = provider;
		task.realizeDataFree = brlobol_rt_proxy_provider_free;
		task.debugDelayMilliseconds =
		    mesh_lod_debug_delay_milliseconds(
			"BRLOBOL_LOD_OBB_TASK_DELAY_MS");
		uint64_t taskId = submitAction->service->submitIfNotActive(task);
		if (taskId != 0) {
		    dependencyTaskId = taskId;
		    submitAction->submittedTaskCount++;
		} else {
		    brlobol_rt_proxy_provider_free(provider);
		}
	    }
	}
    }

    BRLObolMeshLodProvider *provider = new BRLObolMeshLodProvider;
    provider->dbip = submitAction->dbip;
    provider->view = submitAction->view;
    provider->useView = TRUE;
    provider->refreshMissing = submitAction->refreshMissing;
    provider->useForcedLevel = submitAction->useForcedLevel;
    provider->forcedLevel = submitAction->forcedLevel;
    provider->shrinkAfterCopy = TRUE;
    provider->reset = submitAction->reset;

    BRLObolLodTask task;
    task.generation = submitAction->generation;
    task.request = request;
    task.realize = brlobol_mesh_lod_provider_task;
    task.realizeData = provider;
    task.realizeDataFree = brlobol_mesh_lod_provider_free;
    if (dependencyTaskId != 0)
	task.addDependency(dependencyTaskId);
    task.debugDelayMilliseconds =
	mesh_lod_debug_delay_milliseconds("BRLOBOL_LOD_TASK_DELAY_MS");

    uint64_t taskId = suppressActiveDuplicate ?
		      submitAction->service->submitIfNotActive(task) :
		      submitAction->service->submit(task);
    if (taskId == 0) {
	brlobol_mesh_lod_provider_free(provider);
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target,
				       suppressActiveDuplicate &&
				       submitAction->service->hasActiveRequest(request) ?
				       "current LoD request is already active" :
				       "LoD service rejected mesh task");
	return;
    }

    submitAction->submittedTaskCount++;
}

void
SoBRLMeshLodSubmitAction::appendDiagnostic(const SbString &target,
	const char *message)
{
    if (this->diagnostics.getLength() > 0)
	this->diagnostics += "\n";

    this->diagnostics += target.getLength() > 0 ? target : SbString("<unknown>");
    this->diagnostics += ": ";
    this->diagnostics += message ? message : "LoD submit diagnostic";
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
