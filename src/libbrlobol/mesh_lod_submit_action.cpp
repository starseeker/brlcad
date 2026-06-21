/*        M E S H _ L O D _ S U B M I T _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/lod_service.h"
#include "brlobol/mesh_lod_submit_action.h"
#include "brlobol/mesh_shape.h"

#include <Inventor/SbString.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>

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
    providerId("rt_mesh_lod"),
    providerVersion("rt-cache-v1"),
    qualityTier(BRLOBOL_LOD_QUALITY_FAST_DISPLAY),
    refreshMissing(TRUE),
    reset(0),
    useForcedLevel(FALSE),
    forcedLevel(0),
    requireLodBacked(TRUE),
    visitedMeshCount(0),
    submittedTaskCount(0),
    skippedMeshCount(0),
    diagnostics("")
{
    SO_ACTION_CONSTRUCTOR(SoBRLMeshLodSubmitAction);
    rt_view_info_init(&this->view);
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
SoBRLMeshLodSubmitAction::setViewInfo(const struct rt_view_info *info)
{
    if (info)
	this->view = *info;
    else
	rt_view_info_init(&this->view);
    rt_view_info_sanitize(&this->view);
}

const struct rt_view_info &
SoBRLMeshLodSubmitAction::getViewInfo(void) const
{
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
    if (!submitAction->useForcedLevel &&
	    submitAction->reset == 0 &&
	    shape->lodAvailable.getValue() &&
	    shape->lodResultKind.getValue() == BRLOBOL_LOD_RESULT_MESH &&
	    shape->lodProviderStatus.getValue() == BRLOBOL_LOD_PROVIDER_READY &&
	    requestKey.isValid() &&
	    shape->lodCacheKey.getValue().getLength() > 0 &&
	    strcmp(shape->lodCacheKey.getValue().getString(),
		requestKey.value.getString()) == 0) {
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

    BRLObolRtMeshLodProvider *provider = new BRLObolRtMeshLodProvider;
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
    task.realize = brlobol_rt_mesh_lod_provider_task;
    task.realizeData = provider;
    task.realizeDataFree = brlobol_rt_mesh_lod_provider_free;

    uint64_t taskId = suppressActiveDuplicate ?
	submitAction->service->submitIfNotActive(task) :
	submitAction->service->submit(task);
    if (taskId == 0) {
	brlobol_rt_mesh_lod_provider_free(provider);
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
