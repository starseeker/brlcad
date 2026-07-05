/*                 V I E W _ C O N T R O L L E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/edit_preview.h"
#include "brlobol/export_action.h"
#include "brlobol/hud_label_overlay.h"
#include "brlobol/line_layer_overlay.h"
#include "brlobol/lod_service.h"
#include "brlobol/lod_update_action.h"
#include "brlobol/measure_action.h"
#include "brlobol/mesh_lod_submit_action.h"
#include "brlobol/mesh_residency_action.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/pick_detail.h"
#include "brlobol/snap_action.h"
#include "brlobol/view_controller.h"
#include "brlobol/view_lod.h"
#include "brlobol/view_store.h"
#include "raytrace.h"
#include "rt/view.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string.h>
#include <vector>

#include <Inventor/SbName.h>
#include <Inventor/SoRenderManager.h>
#include <Inventor/SoViewport.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoEnvironment.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoLightModel.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/nodes/SoSeparator.h>

static const char *controller_database_id(const struct db_i *dbip);
static std::vector<SoBRLDatabaseSource *> controller_render_database_sources(
    const BRLObolViewController *controller);
static std::vector<SoBRLDatabaseSource *> controller_render_database_source_roots(
    const BRLObolViewController *controller);
static std::vector<SoBRLMeshShape *> controller_render_mesh_shapes(
    const BRLObolViewController *controller);

SbMatrix
brlobol_sbmatrix_from_brl_mat(const mat_t mat)
{
    if (!mat)
	return SbMatrix::identity();

    return SbMatrix(
	       static_cast<float>(mat[0]),  static_cast<float>(mat[4]),
	       static_cast<float>(mat[8]),  static_cast<float>(mat[12]),
	       static_cast<float>(mat[1]),  static_cast<float>(mat[5]),
	       static_cast<float>(mat[9]),  static_cast<float>(mat[13]),
	       static_cast<float>(mat[2]),  static_cast<float>(mat[6]),
	       static_cast<float>(mat[10]), static_cast<float>(mat[14]),
	       static_cast<float>(mat[3]),  static_cast<float>(mat[7]),
	       static_cast<float>(mat[11]), static_cast<float>(mat[15]));
}

SbRotation
brlobol_camera_orientation_from_brl_rotation(const mat_t rotation)
{
    if (!rotation)
	return SbRotation::identity();

    SbMatrix cameraAxes = SbMatrix::identity();
    cameraAxes[0][0] = static_cast<float>(rotation[0]);
    cameraAxes[0][1] = static_cast<float>(rotation[1]);
    cameraAxes[0][2] = static_cast<float>(rotation[2]);
    cameraAxes[1][0] = static_cast<float>(rotation[4]);
    cameraAxes[1][1] = static_cast<float>(rotation[5]);
    cameraAxes[1][2] = static_cast<float>(rotation[6]);
    cameraAxes[2][0] = static_cast<float>(rotation[8]);
    cameraAxes[2][1] = static_cast<float>(rotation[9]);
    cameraAxes[2][2] = static_cast<float>(rotation[10]);
    return SbRotation(cameraAxes);
}

static double
controller_aspect_from_region(const SbViewportRegion &region)
{
    SbVec2s window = region.getWindowSize();
    if (window[0] <= 0 || window[1] <= 0)
	return 0.0;

    return static_cast<double>(window[0]) / static_cast<double>(window[1]);
}

static const char *
controller_render_environment_name(void)
{
    return "BRLObolRenderEnvironment";
}

static SoGroup *
controller_find_render_environment(SoSeparator *root)
{
    if (!root)
	return NULL;

    const char *name = controller_render_environment_name();
    for (int i = 0; i < root->getNumChildren(); i++) {
	SoNode *child = root->getChild(i);
	if (child &&
	    child->isOfType(SoGroup::getClassTypeId()) &&
	    strcmp(child->getName().getString(), name) == 0)
	    return static_cast<SoGroup *>(child);
    }

    return NULL;
}

static void
controller_configure_render_environment(SoViewport *viewport)
{
    if (!viewport || !viewport->getRoot())
	return;

    SoSeparator *root = viewport->getRoot();
    SoGroup *renderEnvironment = controller_find_render_environment(root);
    if (renderEnvironment) {
	const int index = root->findChild(renderEnvironment);
	if (index > 0) {
	    renderEnvironment->ref();
	    root->removeChild(renderEnvironment);
	    root->insertChild(renderEnvironment, 0);
	    renderEnvironment->unref();
	}
	return;
    }

    renderEnvironment = new SoGroup;
    renderEnvironment->setName(SbName(controller_render_environment_name()));

    SoEnvironment *environment = new SoEnvironment;
    environment->ambientIntensity = 0.3f;
    environment->ambientColor = SbColor(1.0f, 1.0f, 1.0f);
    renderEnvironment->addChild(environment);

    SoLightModel *lightModel = new SoLightModel;
    lightModel->model = SoLightModel::PHONG;
    renderEnvironment->addChild(lightModel);

    SoDirectionalLight *headlight = new SoDirectionalLight;
    headlight->color = SbColor(1.0f, 1.0f, 1.0f);
    headlight->intensity = 1.0f;
    headlight->direction = SbVec3f(0.0f, 0.0f, -1.0f);
    renderEnvironment->addChild(headlight);

    root->insertChild(renderEnvironment, 0);
}

static void
controller_append_database_source_unique(
    std::vector<SoBRLDatabaseSource *> &sources,
    SoBRLDatabaseSource *source)
{
    if (!source)
	return;

    const SbString sourceInstance = source->instanceKey.getValue();
    const SbString sourcePath = source->path.getValue();
    const int sourceDrawMode = source->drawMode.getValue();
    struct db_i *sourceDbip = source->getDatabase();

    for (size_t i = 0; i < sources.size(); i++) {
	SoBRLDatabaseSource *existing = sources[i];
	if (!existing)
	    continue;
	if (existing == source)
	    return;

	const SbString existingInstance = existing->instanceKey.getValue();
	if (sourceInstance.getLength() > 0 &&
	    existingInstance.getLength() > 0 &&
	    strcmp(sourceInstance.getString(),
		   existingInstance.getString()) == 0)
	    return;

	if (sourcePath.getLength() > 0 &&
	    existing->path.getValue().getLength() > 0 &&
	    sourceDbip == existing->getDatabase() &&
	    sourceDrawMode == existing->drawMode.getValue() &&
	    strcmp(sourcePath.getString(),
		   existing->path.getValue().getString()) == 0)
	    return;
    }

    sources.push_back(source);
}

static void
controller_collect_database_sources(SoNode *node,
				    std::vector<SoBRLDatabaseSource *> &sources)
{
    if (!node)
	return;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	SoBRLDatabaseSource *source = static_cast<SoBRLDatabaseSource *>(node);
	if (std::find(sources.begin(), sources.end(), source) == sources.end())
	    sources.push_back(source);
    }

    if (!node->isOfType(SoGroup::getClassTypeId()))
	return;

    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++)
	controller_collect_database_sources(group->getChild(i), sources);
}

static void
controller_collect_database_source_roots(
    SoNode *node,
    std::vector<SoBRLDatabaseSource *> &sources)
{
    if (!node)
	return;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	SoBRLDatabaseSource *source = static_cast<SoBRLDatabaseSource *>(node);
	controller_append_database_source_unique(sources, source);
	return;
    }

    if (!node->isOfType(SoGroup::getClassTypeId()))
	return;

    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++)
	controller_collect_database_source_roots(group->getChild(i), sources);
}

static void
controller_collect_mesh_shapes(SoNode *node,
			       std::vector<SoBRLMeshShape *> &shapes)
{
    if (!node)
	return;

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	if (std::find(shapes.begin(), shapes.end(), shape) == shapes.end())
	    shapes.push_back(shape);
    }

    if (!node->isOfType(SoGroup::getClassTypeId()))
	return;

    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++)
	controller_collect_mesh_shapes(group->getChild(i), shapes);
}

static std::vector<SoBRLDatabaseSource *>
controller_render_database_source_roots(
    const BRLObolViewController *controller)
{
    std::vector<SoBRLDatabaseSource *> sources;
    if (!controller)
	return sources;

    SoNode *root = controller->getRenderSceneRoot();
    if (!root)
	root = controller->getSceneRoot();
    controller_collect_database_source_roots(root, sources);
    if (!sources.empty())
	return sources;

    const int sourceCount = controller->getDatabaseSourceCount();
    for (int i = 0; i < sourceCount; i++) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	controller_append_database_source_unique(sources, source);
    }

    return sources;
}

static std::vector<SoBRLDatabaseSource *>
controller_render_database_sources(const BRLObolViewController *controller)
{
    std::vector<SoBRLDatabaseSource *> sources;
    if (!controller)
	return sources;

    SoNode *root = controller->getRenderSceneRoot();
    if (!root)
	root = controller->getSceneRoot();
    controller_collect_database_sources(root, sources);
    if (!sources.empty())
	return sources;

    const int sourceCount = controller->getDatabaseSourceCount();
    for (int i = 0; i < sourceCount; i++) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	if (source && std::find(sources.begin(), sources.end(), source) ==
	    sources.end())
	    sources.push_back(source);
    }

    return sources;
}

static std::vector<SoBRLMeshShape *>
controller_render_mesh_shapes(const BRLObolViewController *controller)
{
    std::vector<SoBRLMeshShape *> shapes;
    if (!controller)
	return shapes;

    SoNode *root = controller->getRenderSceneRoot();
    if (!root)
	root = controller->getSceneRoot();
    controller_collect_mesh_shapes(root, shapes);
    return shapes;
}

static SbString
rt_pick_result_path(const BRLObolRtPickResult &pick)
{
    return pick.detail.getPath();
}

static SbBool
rt_pick_result_path_recorded(
    const std::vector<BRLObolRtPickResult> &results,
    const BRLObolRtPickResult &candidate)
{
    const SbString candidatePath = rt_pick_result_path(candidate);
    if (candidatePath.getLength() == 0)
	return FALSE;

    for (size_t i = 0; i < results.size(); i++) {
	if (strcmp(rt_pick_result_path(results[i]).getString(),
		   candidatePath.getString()) == 0)
	    return TRUE;
    }

    return FALSE;
}

static void
insert_rt_pick_result(std::vector<BRLObolRtPickResult> &results,
		      const BRLObolRtPickResult &pick,
		      SbBool pickAll)
{
    if (pickAll) {
	std::vector<BRLObolRtPickResult>::iterator it = results.begin();
	while (it != results.end() && it->distance <= pick.distance)
	    ++it;
	results.insert(it, pick);
	return;
    }

    if (results.empty() || pick.distance < results[0].distance)
	results.assign(1, pick);
}

static const char *
controller_path_skip_slashes(const char *path)
{
    if (!path)
	return "";
    while (*path == '/')
	path++;
    return path;
}

static const char *
controller_path_leaf(const char *path)
{
    const char *clean = controller_path_skip_slashes(path);
    const char *slash = strrchr(clean, '/');
    return slash && slash[1] ? slash + 1 : clean;
}

static SbBool
controller_path_contains_request(const char *sourcePath,
				 const char *requestPath)
{
    const char *source = controller_path_skip_slashes(sourcePath);
    const char *request = controller_path_skip_slashes(requestPath);
    if (!source[0] || !request[0])
	return FALSE;

    size_t sourceLen = strlen(source);
    if (strncmp(request, source, sourceLen) != 0)
	return FALSE;

    return request[sourceLen] == '\0' || request[sourceLen] == '/' ?
	   TRUE : FALSE;
}

static SbBool
controller_source_matches_request(SoBRLDatabaseSource *source,
				  const BRLObolSourceMeshRequest &request)
{
    if (!source)
	return FALSE;

    const char *sourcePath = source->path.getValue().getString();
    if (controller_path_contains_request(sourcePath,
					 request.path.getString()))
	return TRUE;

    if (request.path.getLength() > 0)
	return FALSE;

    const char *sourceName = request.sourceName.getString();
    if (!sourceName || !sourceName[0])
	return FALSE;

    const char *cleanSourcePath = controller_path_skip_slashes(sourcePath);
    return strcmp(sourceName, cleanSourcePath) == 0 ||
	   strcmp(sourceName, controller_path_leaf(sourcePath)) == 0 ? TRUE :
	   FALSE;
}

static SoBRLDatabaseSource *
controller_database_source_for_request(BRLObolViewController *controller,
				       const BRLObolSourceMeshRequest &request)
{
    SoBRLDatabaseSource *singleSource = NULL;
    SoBRLDatabaseSource *matchedSource = NULL;
    int sourceCountWithDatabase = 0;
    int matchedSourceCount = 0;

    if (!controller)
	return NULL;

    std::vector<SoBRLDatabaseSource *> sources =
	controller_render_database_sources(controller);
    for (size_t i = 0; i < sources.size(); i++) {
	SoBRLDatabaseSource *source = sources[i];
	if (!source || !source->getDatabase())
	    continue;

	sourceCountWithDatabase++;
	if (!singleSource)
	    singleSource = source;

	if (controller_source_matches_request(source, request)) {
	    matchedSource = source;
	    matchedSourceCount++;
	}
    }

    if (matchedSourceCount == 1)
	return matchedSource;
    if (matchedSourceCount > 1)
	return NULL;

    return sourceCountWithDatabase == 1 ? singleSource : NULL;
}

static int
controller_database_source_count(BRLObolViewController *controller)
{
    int count = 0;

    if (!controller)
	return 0;

    std::vector<SoBRLDatabaseSource *> sources =
	controller_render_database_sources(controller);
    for (size_t i = 0; i < sources.size(); i++) {
	SoBRLDatabaseSource *source = sources[i];
	if (source && source->getDatabase())
	    count++;
    }

    return count;
}

static void
controller_source_request_template(BRLObolLodRequest &requestTemplate,
				   SoBRLDatabaseSource *source)
{
    requestTemplate.clear();
    if (!source)
	return;

    requestTemplate.databaseId =
	controller_database_id(source->getDatabase());
    requestTemplate.sourceRevision = source->sourceRevision.getValue();
}

static SbBool
controller_consume_one_source_full_detail_result(
    SoBRLExportAction &action,
    const BRLObolSourceMeshRequest &sourceRequest,
    const BRLObolLodResult &result)
{
    return action.appendSourceBackedFullDetailResult(sourceRequest, result);
}

static SbBool
controller_consume_one_source_full_detail_result(
    SoBRLMeasureAction &action,
    const BRLObolSourceMeshRequest &sourceRequest,
    const BRLObolLodResult &result)
{
    return action.consumeSourceBackedFullDetailResult(sourceRequest, result);
}

static SbBool
controller_consume_one_source_full_detail_result(
    SoBRLSnapAction &action,
    const BRLObolSourceMeshRequest &sourceRequest,
    const BRLObolLodResult &result)
{
    return action.consumeSourceBackedFullDetailResult(sourceRequest, result);
}

template <typename Action>
static int
controller_consume_source_full_detail(BRLObolViewController *controller,
				      Action &action,
				      uint64_t generation,
				      int *submittedRequestCount)
{
    if (submittedRequestCount)
	*submittedRequestCount = 0;

    if (!controller)
	return 0;

    const int requestCount =
	action.getSourceBackedFullDetailRequestCount();
    if (requestCount <= 0)
	return 0;

    BRLObolLodService *service = controller->getLodService();
    if (!service || !service->isRunning())
	return 0;

    std::vector<BRLObolLodRequest> expectedRequests;
    std::vector<BRLObolSourceMeshRequest> expectedSourceRequests;
    std::vector<int> expectedRequestIndices;
    std::vector<BRLObolSourceMeshRequest> submitSourceRequests;
    std::vector<SoBRLDatabaseSource *> requestSources;
    std::vector<int> submitRequestIndices;
    const int databaseSourceCount = controller_database_source_count(controller);
    for (int i = 0; i < requestCount; i++) {
	const BRLObolSourceMeshRequest &sourceRequest =
	    action.getSourceBackedFullDetailRequest(i);
	SoBRLDatabaseSource *source =
	    controller_database_source_for_request(controller, sourceRequest);
	if (source) {
	    BRLObolLodRequest requestTemplate;
	    controller_source_request_template(requestTemplate, source);

	    BRLObolLodRequest request;
	    if (action.makeSourceBackedFullDetailLodRequest(i, request,
		    &requestTemplate)) {
		expectedRequests.push_back(request);
		expectedSourceRequests.push_back(sourceRequest);
		expectedRequestIndices.push_back(i);
		submitSourceRequests.push_back(sourceRequest);
		requestSources.push_back(source);
		submitRequestIndices.push_back(i);
	    }
	}

	if (databaseSourceCount <= 1) {
	    BRLObolLodRequest request;
	    if (action.makeSourceBackedFullDetailLodRequest(i, request)) {
		expectedRequests.push_back(request);
		expectedSourceRequests.push_back(sourceRequest);
		expectedRequestIndices.push_back(i);
	    }
	}
    }
    if (expectedRequests.empty() && submitSourceRequests.empty())
	return 0;

    std::vector<SbBool> requestMatched(
	static_cast<size_t>(requestCount), FALSE);
    int consumed = 0;
    if (!expectedRequests.empty()) {
	std::vector<BRLObolLodResult> sourceResults;
	service->drainMatchingResults(sourceResults, expectedRequests);
	if (!sourceResults.empty()) {
	    std::vector<SbBool> used(sourceResults.size(), FALSE);
	    std::vector<SbBool> requestConsumed(
		static_cast<size_t>(requestCount), FALSE);

	    for (size_t i = 0; i < expectedRequests.size(); i++) {
		const int requestIndex = expectedRequestIndices[i];
		if (requestIndex < 0 ||
		    static_cast<size_t>(requestIndex) >=
		    requestConsumed.size() ||
		    requestConsumed[static_cast<size_t>(requestIndex)])
		    continue;

		for (size_t j = 0; j < sourceResults.size(); j++) {
		    if (used[j] ||
			!brlobol_lod_result_matches_request(
			    sourceResults[j], expectedRequests[i]))
			continue;

		    used[j] = TRUE;
		    requestMatched[static_cast<size_t>(requestIndex)] =
			TRUE;
		    if (controller_consume_one_source_full_detail_result(action,
			    expectedSourceRequests[i], sourceResults[j])) {
			requestConsumed[static_cast<size_t>(requestIndex)] =
			    TRUE;
			consumed++;
		    }
		    break;
		}
	    }
	}
    }

    int submitted = 0;
    for (size_t i = 0; i < submitSourceRequests.size(); i++) {
	const int requestIndex = submitRequestIndices[i];
	if (requestIndex >= 0 &&
	    static_cast<size_t>(requestIndex) < requestMatched.size() &&
	    requestMatched[static_cast<size_t>(requestIndex)])
	    continue;

	BRLObolLodRequest requestTemplate;
	controller_source_request_template(requestTemplate, requestSources[i]);
	if (brlobol_lod_submit_rt_source_full_detail_request(service,
		generation, submitSourceRequests[i],
		requestSources[i]->getDatabase(), &requestTemplate,
		controller->getMaxExactFullDetailFaceCount(),
		controller->getMaxExactFullDetailPointCount()) != 0)
	    submitted++;
    }
    if (submittedRequestCount)
	*submittedRequestCount = submitted;

    return consumed;
}

BRLObolViewController::BRLObolViewController(void) :
    sceneController(),
    viewport(new SoViewport),
    renderLodRoot(NULL),
    viewLodState(new BRLObolViewLodState),
    renderManager(new SoRenderManager),
    activeCamera(NULL),
    viewportRegion(1, 1),
    renderRequested(FALSE),
    renderReason(""),
    lodService(NULL),
    lodResultSubscriberId(0),
    lodResultsPending(0),
    lodAutoSubmit(FALSE),
    lodActiveGeneration(0),
    lodLastSubmittedViewRevision(0),
    lodLastSubmittedPolicyRevision(0),
    lodLastSubmittedSourceSignature(""),
    lodViewSignature(""),
    lodViewRevision(1),
    lodPolicyRevision(1),
    lodUseForcedLevel(FALSE),
    lodForcedLevel(0),
    maxExactFullDetailFaceCount(0),
    maxExactFullDetailPointCount(0),
    meshResidencyBudgetEnabled(FALSE),
    maxResidentMeshBytes(0),
    meshResidencyEvictDisplayPayloads(TRUE),
    lastMeshBudgetInitialResidentBytes(0),
    lastMeshBudgetFinalResidentBytes(0),
    lastMeshBudgetFreedFullDetailBytes(0),
    lastMeshBudgetFreedDisplayBytes(0),
    lastMeshBudgetVisitedMeshCount(0),
    lastMeshBudgetEvictedFullDetailMeshCount(0),
    lastMeshBudgetEvictedDisplayMeshCount(0),
    lastLodVisitedMeshCount(0),
    lastLodSubmittedTaskCount(0),
    lastLodSkippedMeshCount(0),
    lastLodResultCount(0),
    lastLodMatchedResultCount(0),
    lastLodAppliedResultCount(0),
    lastLodRejectedResultCount(0),
    lastLodUnmatchedResultCount(0),
    lastLodDiagnostics(""),
    featureStore(new BRLObolFeatureStore(this)),
    polygonStore(new BRLObolPolygonStore(this)),
    selectionStore(new BRLObolSelectionStore)
{
    controller_configure_render_environment(this->viewport);
}

BRLObolViewController::BRLObolViewController(SoNode *root, SoCamera *camera) :
    sceneController(),
    viewport(new SoViewport),
    renderLodRoot(NULL),
    viewLodState(new BRLObolViewLodState),
    renderManager(new SoRenderManager),
    activeCamera(NULL),
    viewportRegion(1, 1),
    renderRequested(FALSE),
    renderReason(""),
    lodService(NULL),
    lodResultSubscriberId(0),
    lodResultsPending(0),
    lodAutoSubmit(FALSE),
    lodActiveGeneration(0),
    lodLastSubmittedViewRevision(0),
    lodLastSubmittedPolicyRevision(0),
    lodLastSubmittedSourceSignature(""),
    lodViewSignature(""),
    lodViewRevision(1),
    lodPolicyRevision(1),
    lodUseForcedLevel(FALSE),
    lodForcedLevel(0),
    maxExactFullDetailFaceCount(0),
    maxExactFullDetailPointCount(0),
    meshResidencyBudgetEnabled(FALSE),
    maxResidentMeshBytes(0),
    meshResidencyEvictDisplayPayloads(TRUE),
    lastMeshBudgetInitialResidentBytes(0),
    lastMeshBudgetFinalResidentBytes(0),
    lastMeshBudgetFreedFullDetailBytes(0),
    lastMeshBudgetFreedDisplayBytes(0),
    lastMeshBudgetVisitedMeshCount(0),
    lastMeshBudgetEvictedFullDetailMeshCount(0),
    lastMeshBudgetEvictedDisplayMeshCount(0),
    lastLodVisitedMeshCount(0),
    lastLodSubmittedTaskCount(0),
    lastLodSkippedMeshCount(0),
    lastLodResultCount(0),
    lastLodMatchedResultCount(0),
    lastLodAppliedResultCount(0),
    lastLodRejectedResultCount(0),
    lastLodUnmatchedResultCount(0),
    lastLodDiagnostics(""),
    featureStore(new BRLObolFeatureStore(this)),
    polygonStore(new BRLObolPolygonStore(this)),
    selectionStore(new BRLObolSelectionStore)
{
    controller_configure_render_environment(this->viewport);
    this->setSceneRoot(root);
    this->setCamera(camera);
}

BRLObolViewController::~BRLObolViewController(void)
{
    this->setLodService(NULL);
    this->clearRtPickCaches();
    delete this->featureStore;
    this->featureStore = NULL;
    delete this->polygonStore;
    this->polygonStore = NULL;
    delete this->selectionStore;
    this->selectionStore = NULL;
    this->setCamera(NULL);
    this->setSceneRoot(NULL);
    delete this->viewLodState;
    this->viewLodState = NULL;
    this->renderManager->setSceneGraph(NULL);
    delete this->renderManager;
    this->renderManager = NULL;
    delete this->viewport;
    this->viewport = NULL;
}

void
BRLObolViewController::setViewportSceneGraphWithLod(SoNode *root)
{
    if (this->renderLodRoot) {
	this->viewport->setSceneGraph(NULL);
	this->renderLodRoot->unref();
	this->renderLodRoot = NULL;
    }

    if (!root) {
	this->viewport->setSceneGraph(NULL);
	return;
    }

    SoBRLViewLodGroup *wrapper = new SoBRLViewLodGroup;
    wrapper->ref();
    wrapper->setViewLodState(this->viewLodState);
    wrapper->addChild(root);
    this->renderLodRoot = wrapper;
    this->viewport->setSceneGraph(wrapper);
}

void
BRLObolViewController::setSceneRoot(SoNode *root)
{
    this->clearRtPickCaches();
    if (this->viewLodState)
	this->viewLodState->clear();
    this->sceneController.setSceneRoot(root);
    this->setViewportSceneGraphWithLod(root);
    this->syncRenderManager();
    this->requestRender("scene-root");
}

SoNode *
BRLObolViewController::getSceneRoot(void) const
{
    return this->sceneController.getSceneRoot();
}

void
BRLObolViewController::setRenderSceneRoot(SoNode *root)
{
    this->clearRtPickCaches();
    if (this->viewLodState)
	this->viewLodState->clear();
    this->setViewportSceneGraphWithLod(root);
    this->syncRenderManager();
    this->requestRender("render-scene-root");
}

SoNode *
BRLObolViewController::getRenderSceneRoot(void) const
{
    return this->viewport->getSceneGraph();
}

SoNode *
BRLObolViewController::getRenderRoot(void) const
{
    return this->viewport->getRoot();
}

BRLObolViewLodState *
BRLObolViewController::getViewLodState(void) const
{
    return this->viewLodState;
}

void
BRLObolViewController::clearViewLodState(void)
{
    if (this->viewLodState)
	this->viewLodState->clear();
}

void
BRLObolViewController::setCamera(SoCamera *camera)
{
    if (camera)
	camera->ref();
    if (this->activeCamera)
	this->activeCamera->unref();
    this->activeCamera = camera;
    this->viewport->setCamera(camera);
    controller_configure_render_environment(this->viewport);
    this->syncRenderManager();
    this->syncLodViewSignature(TRUE);
    this->requestRender("camera");
}

SoCamera *
BRLObolViewController::getCamera(void) const
{
    return this->activeCamera;
}

void
BRLObolViewController::setViewportRegion(const SbViewportRegion &region)
{
    this->viewportRegion = region;
    this->viewport->setViewportRegion(region);
    this->syncRenderManager();
    this->syncLodViewSignature(TRUE);
    this->requestRender("viewport");
}

const SbViewportRegion &
BRLObolViewController::getViewportRegion(void) const
{
    return this->viewportRegion;
}

void
BRLObolViewController::setViewportSize(unsigned int width, unsigned int height)
{
    if (width > 32767)
	width = 32767;
    if (height > 32767)
	height = 32767;
    if (width == 0)
	width = 1;
    if (height == 0)
	height = 1;
    this->viewportRegion.setWindowSize((short)width, (short)height);
    this->viewport->setViewportRegion(this->viewportRegion);
    this->syncRenderManager();
    this->syncLodViewSignature(TRUE);
    this->requestRender("viewport-size");
}

SbBool
BRLObolViewController::syncCameraFromRtViewContext(const void *viewCtx,
	SbBool createCamera)
{
    if (!viewCtx)
	return FALSE;

    const double perspectiveDegrees =
	rt_view_context_perspective_get(viewCtx);
    const SbBool wantPerspective = perspectiveDegrees > SMALL_FASTF ?
				   TRUE : FALSE;

    SoCamera *camera = this->activeCamera;
    if (!camera ||
	(wantPerspective &&
	 !camera->isOfType(SoPerspectiveCamera::getClassTypeId())) ||
	(!wantPerspective &&
	 !camera->isOfType(SoOrthographicCamera::getClassTypeId()))) {
	if (!createCamera)
	    return FALSE;

	camera = wantPerspective ?
		 static_cast<SoCamera *>(new SoPerspectiveCamera) :
		 static_cast<SoCamera *>(new SoOrthographicCamera);
	this->setCamera(camera);
    }

    int viewWidth = rt_view_context_width_get(viewCtx);
    int viewHeight = rt_view_context_height_get(viewCtx);
    SbVec2s window = this->viewportRegion.getWindowSize();
    if (window[0] <= 1 && window[1] <= 1 &&
	viewWidth > 0 && viewHeight > 0) {
	this->setViewportSize(static_cast<unsigned int>(viewWidth),
			      static_cast<unsigned int>(viewHeight));
	window = this->viewportRegion.getWindowSize();
    }

    double aspect = controller_aspect_from_region(this->viewportRegion);
    if (aspect <= SMALL_FASTF && viewWidth > 0 && viewHeight > 0)
	aspect = static_cast<double>(viewWidth) /
		 static_cast<double>(viewHeight);
    if (aspect <= SMALL_FASTF)
	aspect = 1.0;

    mat_t viewRotation;
    MAT_IDN(viewRotation);
    (void)rt_view_context_rotation_get(viewRotation, viewCtx);

    mat_t viewCenter;
    MAT_IDN(viewCenter);
    (void)rt_view_context_center_get(viewCenter, viewCtx);

    vect_t center;
    MAT_DELTAS_GET_NEG(center, viewCenter);

    const double horizontalSizeRaw = rt_view_context_size_get(viewCtx);
    double horizontalSize = horizontalSizeRaw;
    if (horizontalSize <= SMALL_FASTF) {
	const double scale = rt_view_context_scale_get(viewCtx);
	horizontalSize = scale > SMALL_FASTF ? scale * 2.0 : 2.0;
    }
    const double verticalSize = horizontalSize / aspect;

    double heightAngle = perspectiveDegrees * DEG2RAD;
    if (heightAngle <= SMALL_FASTF)
	heightAngle = 2.0 * std::atan(0.1);
    if (heightAngle < 0.001)
	heightAngle = 0.001;
    if (heightAngle > 3.0)
	heightAngle = 3.0;

    const SbBool orthographic =
	camera->isOfType(SoOrthographicCamera::getClassTypeId());
    double distance = orthographic ?
		      horizontalSize * 5.0 :
		      (verticalSize * 0.5) / std::tan(heightAngle * 0.5);
    if (distance <= horizontalSize * 0.5)
	distance = horizontalSize * 0.5;
    if (distance <= SMALL_FASTF)
	distance = 1.0;

    const SbRotation orientation =
	brlobol_camera_orientation_from_brl_rotation(viewRotation);
    const double viewZ[3] = {
	viewRotation[8], viewRotation[9], viewRotation[10]
    };

    camera->viewportMapping = SoCamera::LEAVE_ALONE;
    camera->aspectRatio = static_cast<float>(aspect);
    camera->position = SbVec3f(
			   static_cast<float>(center[X] + viewZ[X] * distance),
			   static_cast<float>(center[Y] + viewZ[Y] * distance),
			   static_cast<float>(center[Z] + viewZ[Z] * distance));
    camera->orientation = orientation;
    camera->focalDistance = static_cast<float>(distance);
    camera->nearDistance = static_cast<float>(
			       std::max(distance * 0.001, 1.0e-6));
    camera->farDistance = static_cast<float>(
			      std::max(distance + horizontalSize * 100.0,
				       distance + 1.0));

    if (camera->isOfType(SoPerspectiveCamera::getClassTypeId())) {
	SoPerspectiveCamera *perspectiveCamera =
	    static_cast<SoPerspectiveCamera *>(camera);
	perspectiveCamera->heightAngle = static_cast<float>(heightAngle);
    } else if (orthographic) {
	SoOrthographicCamera *orthographicCamera =
	    static_cast<SoOrthographicCamera *>(camera);
	orthographicCamera->height = static_cast<float>(verticalSize);
    }

    this->syncLodViewSignature(TRUE);
    this->requestRender("rt-view-camera");
    return TRUE;
}

SbBool
BRLObolViewController::getRtViewInfo(struct rt_view_info *info) const
{
    if (!info)
	return FALSE;

    rt_view_info_init(info);

    SbVec2s window = this->viewportRegion.getWindowSize();
    info->width = window[0] > 0 ? window[0] : 1;
    info->height = window[1] > 0 ? window[1] : 1;

    if (!this->activeCamera) {
	rt_view_info_sanitize(info);
	return FALSE;
    }

    if (this->activeCamera->isOfType(SoOrthographicCamera::getClassTypeId())) {
	SoOrthographicCamera *camera =
	    static_cast<SoOrthographicCamera *>(this->activeCamera);
	double aspect = controller_aspect_from_region(this->viewportRegion);
	if (aspect <= SMALL_FASTF)
	    aspect = camera->aspectRatio.getValue();
	if (aspect <= SMALL_FASTF)
	    aspect = 1.0;
	info->size = camera->height.getValue() * aspect;
    } else if (this->activeCamera->isOfType(SoPerspectiveCamera::getClassTypeId())) {
	SoPerspectiveCamera *camera =
	    static_cast<SoPerspectiveCamera *>(this->activeCamera);
	double focal = this->activeCamera->focalDistance.getValue();
	double angle = camera->heightAngle.getValue();
	double aspect = controller_aspect_from_region(this->viewportRegion);
	if (aspect <= SMALL_FASTF)
	    aspect = this->activeCamera->aspectRatio.getValue();
	if (aspect <= SMALL_FASTF)
	    aspect = 1.0;
	if (focal <= 0.0)
	    focal = 1.0;
	if (angle <= 0.0)
	    angle = 2.0 * std::atan(0.1);
	info->size = 2.0 * focal * std::tan(angle * 0.5) * aspect;
    } else {
	info->size = this->activeCamera->focalDistance.getValue();
    }

    rt_view_info_sanitize(info);
    return TRUE;
}

SbBool
BRLObolViewController::realizePending(void)
{
    const SbBool ret = this->sceneController.realizePending();
    this->requestRender(ret ? "realize" : "realize-failed");
    return ret;
}

unsigned int
BRLObolViewController::getLastVisitedSourceCount(void) const
{
    return this->sceneController.getLastVisitedSourceCount();
}

unsigned int
BRLObolViewController::getLastRealizedSourceCount(void) const
{
    return this->sceneController.getLastRealizedSourceCount();
}

unsigned int
BRLObolViewController::getLastFailedSourceCount(void) const
{
    return this->sceneController.getLastFailedSourceCount();
}

const SbString &
BRLObolViewController::getLastDiagnostics(void) const
{
    return this->sceneController.getLastDiagnostics();
}

void
BRLObolViewController::requestRender(const char *reason)
{
    this->renderRequested = TRUE;
    this->renderReason = reason ? reason : "";
}

void
BRLObolViewController::clearRenderRequest(void)
{
    this->renderRequested = FALSE;
    this->renderReason = "";
}

SbBool
BRLObolViewController::consumeRenderRequest(SbString *reason)
{
    const SbBool ret = this->renderRequested;
    if (reason)
	*reason = this->renderReason;
    this->clearRenderRequest();
    return ret;
}

SbBool
BRLObolViewController::renderPending(SbBool clearWindow,
				     SbBool clearZBuffer,
				     SbString *reason)
{
    if (this->lodAutoSubmit)
	(void)this->submitLodRequestsIfNeeded();

    if (this->hasPendingLodResults())
	this->processPendingLodResults();

    if (!this->renderRequested || !this->renderManager ||
	!this->activeCamera || !this->getRenderRoot())
	return FALSE;

    SbString renderReasonCopy;
    if (!this->consumeRenderRequest(&renderReasonCopy))
	return FALSE;

    if (reason)
	*reason = renderReasonCopy;

    this->renderManager->render(clearWindow, clearZBuffer);
    return TRUE;
}

SbBool
BRLObolViewController::isRenderRequested(void) const
{
    return (this->renderRequested || this->hasPendingLodResults()) ?
	   TRUE : FALSE;
}

const SbString &
BRLObolViewController::getRenderReason(void) const
{
    return this->renderReason;
}

static int
size_to_int_saturated(size_t value)
{
    return value > static_cast<size_t>(std::numeric_limits<int>::max()) ?
	   std::numeric_limits<int>::max() : static_cast<int>(value);
}

static void
append_controller_lod_diagnostic(SbString &diagnostics,
				 const SbString &target,
				 const char *message)
{
    if (diagnostics.getLength() > 0)
	diagnostics += "\n";

    diagnostics += target.getLength() > 0 ? target : SbString("<unknown>");
    diagnostics += ": ";
    diagnostics += message ? message : "LoD controller diagnostic";
}

static const char *
controller_database_id(const struct db_i *dbip)
{
    return dbip && dbip->dbi_filename ? dbip->dbi_filename : "";
}

static void
append_signature_string(std::ostringstream &out, const char *value)
{
    size_t len = value ? strlen(value) : 0;
    out << len << ':' << (value ? value : "") << ';';
}

static SbString
controller_lod_view_signature(const struct rt_view_info &view,
			      SbBool haveCamera)
{
    std::ostringstream out;

    out << (haveCamera ? 1 : 0) << ';'
	<< view.width << ';'
	<< view.height << ';'
	<< std::setprecision(17)
	<< view.size << ';'
	<< view.lod.scale << ';'
	<< view.lod.curve_scale << ';'
	<< view.lod.point_scale << ';'
	<< view.lod.bot_threshold << ';';

    return SbString(out.str().c_str());
}

static SbString
controller_lod_source_signature(const BRLObolViewController *controller)
{
    std::ostringstream out;

    if (!controller)
	return SbString("");

    std::vector<SoBRLDatabaseSource *> sources =
	controller_render_database_sources(controller);
    for (size_t i = 0; i < sources.size(); i++) {
	SoBRLDatabaseSource *source = sources[i];
	if (!source || !source->getDatabase())
	    continue;
	if (source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
	    source->needsRealization() ||
	    source->getRealizedMeshCount() <= 0)
	    continue;

	out << "source;";
	append_signature_string(out, controller_database_id(source->getDatabase()));
	append_signature_string(out, source->path.getValue().getString());
	out << source->drawMode.getValue() << ';'
	    << source->lodBotThreshold.getValue() << ';'
	    << source->sourceRevision.getValue() << ';'
	    << source->inputsRevision.getValue() << ';'
	    << source->realizedSourceRevision.getValue() << ';'
	    << source->realizedInputsRevision.getValue() << ';'
	    << source->realizedViewRevision.getValue() << ';';
    }

    return SbString(out.str().c_str());
}

void
BRLObolViewController::lodResultReadyCB(
    BRLObolLodService *UNUSED(service), void *userData)
{
    BRLObolViewController *controller =
	static_cast<BRLObolViewController *>(userData);
    if (controller)
	controller->lodResultsPending.store(1);
}

void
BRLObolViewController::setLodService(BRLObolLodService *service)
{
    if (this->lodService == service)
	return;

    if (this->lodService && this->lodResultSubscriberId != 0)
	this->lodService->unsubscribeResultReady(this->lodResultSubscriberId);

    this->lodService = service;
    this->lodResultSubscriberId = 0;
    this->lodResultsPending.store(0);
    this->lodActiveGeneration = 0;
    this->lodLastSubmittedViewRevision = 0;
    this->lodLastSubmittedPolicyRevision = 0;
    this->lodLastSubmittedSourceSignature = "";

    if (this->lodService)
	this->lodResultSubscriberId =
	    this->lodService->subscribeResultReady(
		BRLObolViewController::lodResultReadyCB, this);
}

BRLObolLodService *
BRLObolViewController::getLodService(void) const
{
    return this->lodService;
}

void
BRLObolViewController::setLodAutoSubmit(SbBool enabled)
{
    this->lodAutoSubmit = enabled ? TRUE : FALSE;
    if (this->lodAutoSubmit)
	this->requestRender("lod-auto-submit");
}

SbBool
BRLObolViewController::isLodAutoSubmitEnabled(void) const
{
    return this->lodAutoSubmit;
}

void
BRLObolViewController::setLodForcedLevel(int level)
{
    if (level < 0)
	level = 0;

    if (this->lodUseForcedLevel && this->lodForcedLevel == level)
	return;

    this->lodUseForcedLevel = TRUE;
    this->lodForcedLevel = level;
    this->advanceLodPolicyRevision();
    this->requestRender("lod-policy");
}

void
BRLObolViewController::clearLodForcedLevel(void)
{
    if (!this->lodUseForcedLevel)
	return;

    this->lodUseForcedLevel = FALSE;
    this->lodForcedLevel = 0;
    this->advanceLodPolicyRevision();
    this->requestRender("lod-policy");
}

SbBool
BRLObolViewController::hasLodForcedLevel(void) const
{
    return this->lodUseForcedLevel;
}

int
BRLObolViewController::getLodForcedLevel(void) const
{
    return this->lodForcedLevel;
}

void
BRLObolViewController::setExactFullDetailBudget(uint64_t maxFaceCount,
						uint64_t maxPointCount)
{
    this->maxExactFullDetailFaceCount = maxFaceCount;
    this->maxExactFullDetailPointCount = maxPointCount;
}

uint64_t
BRLObolViewController::getMaxExactFullDetailFaceCount(void) const
{
    return this->maxExactFullDetailFaceCount;
}

uint64_t
BRLObolViewController::getMaxExactFullDetailPointCount(void) const
{
    return this->maxExactFullDetailPointCount;
}

int
BRLObolViewController::consumeExportSourceFullDetail(
    SoBRLExportAction &exportAction,
    uint64_t generation,
    int *submittedRequestCount)
{
    return controller_consume_source_full_detail(this, exportAction,
	    generation, submittedRequestCount);
}

int
BRLObolViewController::consumeMeasureSourceFullDetail(
    SoBRLMeasureAction &measureAction,
    uint64_t generation,
    int *submittedRequestCount)
{
    return controller_consume_source_full_detail(this, measureAction,
	    generation, submittedRequestCount);
}

int
BRLObolViewController::consumeSnapSourceFullDetail(
    SoBRLSnapAction &snapAction,
    uint64_t generation,
    int *submittedRequestCount)
{
    return controller_consume_source_full_detail(this, snapAction,
	    generation, submittedRequestCount);
}

void
BRLObolViewController::clearRtPickCaches(void)
{
    for (size_t i = 0; i < this->rtPickCaches.size(); i++)
	delete this->rtPickCaches[i];
    this->rtPickCaches.clear();
    this->rtPickCachePaths.clear();
    this->rtPickCacheDatabases.clear();
    this->rtPickCacheSourceRevisions.clear();
}

int
BRLObolViewController::prepareRtPickCaches(void)
{
    std::vector<SbString> sourcePaths;
    std::vector<struct db_i *> sourceDatabases;
    std::vector<uint32_t> sourceRevisions;

    std::vector<SoBRLDatabaseSource *> sources =
	controller_render_database_sources(this);
    for (size_t i = 0; i < sources.size(); i++) {
	SoBRLDatabaseSource *source = sources[i];
	if (!source || !source->getDatabase() ||
	    source->path.getValue().getLength() == 0)
	    continue;
	sourcePaths.push_back(source->path.getValue());
	sourceDatabases.push_back(source->getDatabase());
	sourceRevisions.push_back(
	    static_cast<uint32_t>(source->sourceRevision.getValue()));
    }

    SbBool sameSignature =
	sourcePaths.size() == this->rtPickCachePaths.size() &&
	sourceDatabases.size() == this->rtPickCacheDatabases.size() &&
	sourceRevisions.size() == this->rtPickCacheSourceRevisions.size() &&
	this->rtPickCaches.size() == this->rtPickCachePaths.size();
    if (sameSignature) {
	for (size_t i = 0; i < sourcePaths.size(); i++) {
	    if (sourceDatabases[i] != this->rtPickCacheDatabases[i] ||
		sourceRevisions[i] != this->rtPickCacheSourceRevisions[i] ||
		strcmp(sourcePaths[i].getString(),
		       this->rtPickCachePaths[i].getString()) != 0 ||
		!this->rtPickCaches[i] ||
		!this->rtPickCaches[i]->isReady()) {
		sameSignature = FALSE;
		break;
	    }
	}
    }

    if (sameSignature)
	return static_cast<int>(this->rtPickCaches.size());

    this->clearRtPickCaches();
    for (size_t i = 0; i < sourcePaths.size(); i++) {
	std::vector<SbString> objectPaths;
	objectPaths.push_back(sourcePaths[i]);

	BRLObolRtPickCache *cache = new BRLObolRtPickCache;
	if (!cache->prepare(sourceDatabases[i], objectPaths)) {
	    delete cache;
	    continue;
	}

	this->rtPickCaches.push_back(cache);
	this->rtPickCachePaths.push_back(sourcePaths[i]);
	this->rtPickCacheDatabases.push_back(sourceDatabases[i]);
	this->rtPickCacheSourceRevisions.push_back(sourceRevisions[i]);
    }

    return static_cast<int>(this->rtPickCaches.size());
}

int
BRLObolViewController::getRtPickCacheCount(void) const
{
    return static_cast<int>(this->rtPickCaches.size());
}

BRLObolRtPickCache *
BRLObolViewController::getRtPickCache(int index) const
{
    if (index < 0 || static_cast<size_t>(index) >= this->rtPickCaches.size())
	return NULL;
    return this->rtPickCaches[static_cast<size_t>(index)];
}

uint32_t
BRLObolViewController::getRtPickCacheSourceRevision(int index) const
{
    if (index < 0 ||
	static_cast<size_t>(index) >=
	this->rtPickCacheSourceRevisions.size())
	return 0;
    return this->rtPickCacheSourceRevisions[static_cast<size_t>(index)];
}

int
BRLObolViewController::pickSourceMeshExactRay(
    BRLObolSourceMeshPickResult &pick,
    const SbVec3f &rayOrigin,
    const SbVec3f &rayDirection,
    uint64_t generation,
    int *submittedRequestCount)
{
    pick.clear();
    if (submittedRequestCount)
	*submittedRequestCount = 0;

    if (!this->viewport || !this->viewport->getRoot())
	return 0;

    BRLObolLodService *service = this->getLodService();
    if (!service || !service->isRunning())
	return 0;

    SoBRLSourceMeshPickAction sourcePickAction;
    sourcePickAction.setRay(rayOrigin, rayDirection);
    sourcePickAction.apply(this->viewport->getRoot());

    const int requestCount =
	sourcePickAction.getSourceBackedFullDetailRequestCount();
    if (requestCount <= 0)
	return 0;

    std::vector<BRLObolLodRequest> expectedRequests;
    std::vector<BRLObolSourceMeshRequest> expectedSourceRequests;
    std::vector<int> expectedRequestIndices;
    std::vector<BRLObolSourceMeshRequest> submitSourceRequests;
    std::vector<SoBRLDatabaseSource *> requestSources;
    std::vector<int> submitRequestIndices;
    const int databaseSourceCount = controller_database_source_count(this);
    for (int i = 0; i < requestCount; i++) {
	const BRLObolSourceMeshRequest &sourceRequest =
	    sourcePickAction.getSourceBackedFullDetailRequest(i);
	SoBRLDatabaseSource *source =
	    controller_database_source_for_request(this, sourceRequest);
	if (source) {
	    BRLObolLodRequest requestTemplate;
	    controller_source_request_template(requestTemplate, source);

	    BRLObolLodRequest request;
	    if (sourcePickAction.makeSourceBackedFullDetailLodRequest(i,
		    request, &requestTemplate)) {
		expectedRequests.push_back(request);
		expectedSourceRequests.push_back(sourceRequest);
		expectedRequestIndices.push_back(i);
		submitSourceRequests.push_back(sourceRequest);
		requestSources.push_back(source);
		submitRequestIndices.push_back(i);
	    }
	}

	if (databaseSourceCount <= 1) {
	    BRLObolLodRequest request;
	    if (sourcePickAction.makeSourceBackedFullDetailLodRequest(i,
		    request)) {
		expectedRequests.push_back(request);
		expectedSourceRequests.push_back(sourceRequest);
		expectedRequestIndices.push_back(i);
	    }
	}
    }
    if (expectedRequests.empty() && submitSourceRequests.empty())
	return 0;

    std::vector<SbBool> requestMatched(
	static_cast<size_t>(requestCount), FALSE);
    if (!expectedRequests.empty()) {
	std::vector<BRLObolLodResult> sourceResults;
	service->drainMatchingResults(sourceResults, expectedRequests);
	if (!sourceResults.empty()) {
	    std::vector<SbBool> used(sourceResults.size(), FALSE);
	    std::vector<SbBool> requestConsumed(
		static_cast<size_t>(requestCount), FALSE);
	    for (size_t i = 0; i < expectedRequests.size(); i++) {
		const int requestIndex = expectedRequestIndices[i];
		if (requestIndex < 0 ||
		    static_cast<size_t>(requestIndex) >=
		    requestConsumed.size() ||
		    requestConsumed[static_cast<size_t>(requestIndex)])
		    continue;

		for (size_t j = 0; j < sourceResults.size(); j++) {
		    if (used[j] ||
			!brlobol_lod_result_matches_request(
			    sourceResults[j], expectedRequests[i]))
			continue;

		    used[j] = TRUE;
		    requestMatched[static_cast<size_t>(requestIndex)] =
			TRUE;
		    BRLObolSourceMeshPickResult candidate;
		    if (brlobol_pick_source_full_detail_result(candidate,
			    expectedSourceRequests[i], sourceResults[j],
			    sourcePickAction.getRayOrigin(),
			    sourcePickAction.getRayDirection())) {
			requestConsumed[static_cast<size_t>(requestIndex)] =
			    TRUE;
			if (!pick.hit || candidate.distance < pick.distance)
			    pick = candidate;
		    }
		    break;
		}
	    }
	}
    }

    int submitted = 0;
    for (size_t i = 0; i < submitSourceRequests.size(); i++) {
	const int requestIndex = submitRequestIndices[i];
	if (requestIndex >= 0 &&
	    static_cast<size_t>(requestIndex) < requestMatched.size() &&
	    requestMatched[static_cast<size_t>(requestIndex)])
	    continue;

	BRLObolLodRequest requestTemplate;
	controller_source_request_template(requestTemplate, requestSources[i]);
	if (brlobol_lod_submit_rt_source_full_detail_request(service,
		generation, submitSourceRequests[i],
		requestSources[i]->getDatabase(), &requestTemplate,
		this->getMaxExactFullDetailFaceCount(),
		this->getMaxExactFullDetailPointCount()) != 0)
	    submitted++;
    }
    if (submittedRequestCount)
	*submittedRequestCount = submitted;

    if (pick.hit)
	return 1;
    pick.clear();
    return 0;
}

int
BRLObolViewController::pickRtExactRay(
    std::vector<BRLObolRtPickResult> &results,
    const SbVec3f &rayOrigin,
    const SbVec3f &rayDirection,
    SbBool pickAll)
{
    results.clear();

    SbVec3f direction = rayDirection;
    if (direction.length() <= 0.0f)
	return 0;
    direction.normalize();

    const int cacheCount = this->prepareRtPickCaches();
    for (int i = 0; i < cacheCount; i++) {
	BRLObolRtPickCache *cache = this->getRtPickCache(i);
	if (!cache || !cache->isReady())
	    continue;

	BRLObolRtPickResult rtPick;
	if (!cache->pickRay(rtPick, rayOrigin, direction) || !rtPick.hit)
	    continue;
	if (rt_pick_result_path_recorded(results, rtPick))
	    continue;

	insert_rt_pick_result(results, rtPick, pickAll ? TRUE : FALSE);
    }

    return static_cast<int>(results.size());
}

void
BRLObolViewController::setMeshResidencyBudget(
    size_t maxBytes,
    SbBool evictDisplayPayloads)
{
    this->meshResidencyBudgetEnabled = TRUE;
    this->maxResidentMeshBytes = maxBytes;
    this->meshResidencyEvictDisplayPayloads =
	evictDisplayPayloads ? TRUE : FALSE;
}

void
BRLObolViewController::clearMeshResidencyBudget(void)
{
    this->meshResidencyBudgetEnabled = FALSE;
    this->maxResidentMeshBytes = 0;
    this->meshResidencyEvictDisplayPayloads = TRUE;
}

SbBool
BRLObolViewController::hasMeshResidencyBudget(void) const
{
    return this->meshResidencyBudgetEnabled;
}

size_t
BRLObolViewController::getMaxResidentMeshBytes(void) const
{
    return this->maxResidentMeshBytes;
}

SbBool
BRLObolViewController::isMeshResidencyDisplayEvictionEnabled(void) const
{
    return this->meshResidencyEvictDisplayPayloads;
}

size_t
BRLObolViewController::evictMeshPayloadsToBudget(
    size_t maxBytes,
    SbBool evictDisplayPayloads)
{
    this->lastMeshBudgetInitialResidentBytes = 0;
    this->lastMeshBudgetFinalResidentBytes = 0;
    this->lastMeshBudgetFreedFullDetailBytes = 0;
    this->lastMeshBudgetFreedDisplayBytes = 0;
    this->lastMeshBudgetVisitedMeshCount = 0;
    this->lastMeshBudgetEvictedFullDetailMeshCount = 0;
    this->lastMeshBudgetEvictedDisplayMeshCount = 0;

    SoNode *root = this->getRenderSceneRoot();
    if (!root)
	root = this->getSceneRoot();
    if (!root)
	return 0;

    SoBRLMeshResidencyAction action;
    action.setMaxResidentMeshBytes(maxBytes);
    action.setEvictDisplayPayloads(evictDisplayPayloads);
    action.apply(root);

    const size_t viewLodBytes = this->viewLodState ?
				this->viewLodState->estimateDisplayMeshBytes() : 0;
    this->lastMeshBudgetInitialResidentBytes =
	action.getInitialResidentMeshBytes() + viewLodBytes;
    this->lastMeshBudgetFinalResidentBytes =
	action.getFinalResidentMeshBytes() + viewLodBytes;
    this->lastMeshBudgetFreedFullDetailBytes =
	action.getFreedFullDetailBytes();
    this->lastMeshBudgetFreedDisplayBytes =
	action.getFreedDisplayBytes();
    this->lastMeshBudgetVisitedMeshCount = action.getVisitedMeshCount();
    this->lastMeshBudgetEvictedFullDetailMeshCount =
	action.getEvictedFullDetailMeshCount();
    this->lastMeshBudgetEvictedDisplayMeshCount =
	action.getEvictedDisplayMeshCount();

    if (this->lastMeshBudgetFinalResidentBytes > maxBytes &&
	this->viewLodState) {
	std::vector<SoBRLMeshShape *> shapes =
	    controller_render_mesh_shapes(this);
	for (size_t i = 0; i < shapes.size(); i++) {
	    if (this->lastMeshBudgetFinalResidentBytes <= maxBytes)
		break;
	    if (!this->viewLodState->findMesh(shapes[i]))
		continue;
	    size_t freed = shapes[i]->evictDisplayMeshPreservingSourceMetrics();
	    if (freed == 0)
		continue;
	    this->lastMeshBudgetFreedFullDetailBytes += freed;
	    this->lastMeshBudgetEvictedFullDetailMeshCount++;
	    this->lastMeshBudgetFinalResidentBytes =
		this->lastMeshBudgetFinalResidentBytes > freed ?
		this->lastMeshBudgetFinalResidentBytes - freed : 0;
	}
    }

    if (evictDisplayPayloads &&
	this->lastMeshBudgetFinalResidentBytes > maxBytes &&
	this->viewLodState) {
	unsigned int evicted = 0;
	size_t freed = this->viewLodState->evictDisplayMeshes(&evicted);
	if (freed > 0) {
	    this->lastMeshBudgetFreedDisplayBytes += freed;
	    this->lastMeshBudgetEvictedDisplayMeshCount += evicted;
	    this->lastMeshBudgetFinalResidentBytes =
		this->lastMeshBudgetFinalResidentBytes > freed ?
		this->lastMeshBudgetFinalResidentBytes - freed : 0;
	}
    }

    size_t freedBytes = this->getLastMeshBudgetFreedResidentBytes();
    if (freedBytes > 0)
	this->requestRender("lod-memory-budget");
    return freedBytes;
}

size_t
BRLObolViewController::enforceMeshResidencyBudget(void)
{
    if (!this->meshResidencyBudgetEnabled)
	return 0;

    return this->evictMeshPayloadsToBudget(
	       this->maxResidentMeshBytes,
	       this->meshResidencyEvictDisplayPayloads);
}

size_t
BRLObolViewController::getLastMeshBudgetInitialResidentBytes(void) const
{
    return this->lastMeshBudgetInitialResidentBytes;
}

size_t
BRLObolViewController::getLastMeshBudgetFinalResidentBytes(void) const
{
    return this->lastMeshBudgetFinalResidentBytes;
}

size_t
BRLObolViewController::getLastMeshBudgetFreedResidentBytes(void) const
{
    return this->lastMeshBudgetInitialResidentBytes >
	   this->lastMeshBudgetFinalResidentBytes ?
	   this->lastMeshBudgetInitialResidentBytes -
	   this->lastMeshBudgetFinalResidentBytes : 0;
}

size_t
BRLObolViewController::getLastMeshBudgetFreedFullDetailBytes(void) const
{
    return this->lastMeshBudgetFreedFullDetailBytes;
}

size_t
BRLObolViewController::getLastMeshBudgetFreedDisplayBytes(void) const
{
    return this->lastMeshBudgetFreedDisplayBytes;
}

unsigned int
BRLObolViewController::getLastMeshBudgetVisitedMeshCount(void) const
{
    return this->lastMeshBudgetVisitedMeshCount;
}

unsigned int
BRLObolViewController::getLastMeshBudgetEvictedFullDetailMeshCount(void) const
{
    return this->lastMeshBudgetEvictedFullDetailMeshCount;
}

unsigned int
BRLObolViewController::getLastMeshBudgetEvictedDisplayMeshCount(void) const
{
    return this->lastMeshBudgetEvictedDisplayMeshCount;
}

SbBool
BRLObolViewController::hasPendingLodResults(void) const
{
    return this->lodResultsPending.load() != 0 ? TRUE : FALSE;
}

size_t
BRLObolViewController::processPendingLodResults(size_t maxResults)
{
    if (!this->lodService)
	return 0;

    if (!this->hasPendingLodResults() &&
	this->lodService->queuedResultCountForDiagnostics() == 0)
	return 0;

    (void)this->applyLodResults(this->lodService, maxResults);
    return this->lastLodResultCount;
}

int
BRLObolViewController::submitLodRequestsIfNeeded(SbBool refreshMissing,
	int reset)
{
    if (!this->lodService || !this->lodService->isRunning())
	return 0;
    if (!this->activeCamera ||
	(!this->getSceneRoot() && !this->getRenderSceneRoot()))
	return 0;

    this->syncLodViewSignature(TRUE);

    SbString signature = controller_lod_source_signature(this);
    if (signature.getLength() == 0)
	return 0;

    if (this->lodLastSubmittedViewRevision == this->lodViewRevision &&
	this->lodLastSubmittedPolicyRevision == this->lodPolicyRevision &&
	strcmp(this->lodLastSubmittedSourceSignature.getString(),
	       signature.getString()) == 0)
	return 0;

    if (this->lodActiveGeneration != 0)
	this->lodService->cancelGeneration(this->lodActiveGeneration);

    uint64_t generation = this->lodService->beginGeneration();
    int submitted = this->submitLodRequests(this->lodService, generation,
					    refreshMissing, reset);
    if (submitted >= 0) {
	this->lodActiveGeneration = generation;
	this->lodLastSubmittedViewRevision = this->lodViewRevision;
	this->lodLastSubmittedPolicyRevision = this->lodPolicyRevision;
	this->lodLastSubmittedSourceSignature = signature;
    }

    return submitted;
}

int
BRLObolViewController::submitLodRequests(BRLObolLodService *service,
	uint64_t generation,
	SbBool refreshMissing,
	int reset)
{
    if (!service)
	service = this->lodService;

    this->lastLodVisitedMeshCount = 0;
    this->lastLodSubmittedTaskCount = 0;
    this->lastLodSkippedMeshCount = 0;
    this->lastLodDiagnostics = "";

    if (!service || !service->isRunning()) {
	this->lastLodDiagnostics = "LoD service is not running";
	return -1;
    }

    struct rt_view_info view = RT_VIEW_INFO_INIT;
    if (!this->getRtViewInfo(&view)) {
	this->lastLodDiagnostics = "LoD submission requires an active camera";
	return -1;
    }

    if (!this->getSceneRoot() && !this->getRenderSceneRoot()) {
	this->lastLodDiagnostics = "LoD submission requires a scene root";
	return -1;
    }

    std::vector<SoBRLDatabaseSource *> sources =
	controller_render_database_source_roots(this);
    for (size_t i = 0; i < sources.size(); i++) {
	SoBRLDatabaseSource *source = sources[i];
	if (!source)
	    continue;

	struct db_i *dbip = source->getDatabase();
	if (!dbip) {
	    append_controller_lod_diagnostic(this->lastLodDiagnostics,
					     source->path.getValue(),
					     "database source has no database for LoD submission");
	    continue;
	}

	SoBRLMeshLodSubmitAction action;
	action.setService(service);
	action.setDatabase(dbip, controller_database_id(dbip),
			   source->sourceRevision.getValue());
	action.setViewInfo(&view);
	action.setGeneration(generation);
	action.setRevisions(this->lodViewRevision, this->lodPolicyRevision);
	action.setRefreshMissing(refreshMissing);
	action.setReset(reset);
	action.setViewLodState(this->viewLodState);
	action.setProxyStages(TRUE, TRUE);
	if (this->lodUseForcedLevel)
	    action.setForcedLevel(this->lodForcedLevel);
	action.apply(source);

	this->lastLodVisitedMeshCount += action.getVisitedMeshCount();
	this->lastLodSubmittedTaskCount += action.getSubmittedTaskCount();
	this->lastLodSkippedMeshCount += action.getSkippedMeshCount();
	if (action.getDiagnostics().getLength() > 0)
	    append_controller_lod_diagnostic(this->lastLodDiagnostics,
					     source->path.getValue(),
					     action.getDiagnostics().getString());
    }

    if (this->lastLodSubmittedTaskCount > 0)
	this->requestRender("lod-submit");

    return size_to_int_saturated(
	       static_cast<size_t>(this->lastLodSubmittedTaskCount));
}

int
BRLObolViewController::applyLodResults(BRLObolLodService *service,
				       size_t maxResults)
{
    if (!service)
	service = this->lodService;

    this->lastLodResultCount = 0;
    this->lastLodMatchedResultCount = 0;
    this->lastLodAppliedResultCount = 0;
    this->lastLodRejectedResultCount = 0;
    this->lastLodUnmatchedResultCount = 0;
    this->lastLodDiagnostics = "";

    if (!service) {
	this->lastLodDiagnostics = "LoD service is not configured";
	return -1;
    }

    SoNode *root = this->getRenderSceneRoot();
    if (!root)
	root = this->getSceneRoot();
    if (!root) {
	this->lastLodDiagnostics = "LoD result application requires a scene root";
	return -1;
    }

    std::vector<BRLObolLodResult> drained;
    this->lastLodResultCount = service->drainResults(drained, maxResults);
    this->lodResultsPending.store(
	service->queuedResultCountForDiagnostics() > 0 ? 1 : 0);
    if (this->lastLodResultCount == 0)
	return 0;

    SoBRLLodUpdateAction update;
    update.setViewLodState(this->viewLodState);
    for (size_t i = 0; i < drained.size(); i++) {
	if (drained[i].request.viewRevision != this->lodViewRevision ||
	    drained[i].request.policyRevision != this->lodPolicyRevision) {
	    this->lastLodRejectedResultCount++;
	    append_controller_lod_diagnostic(this->lastLodDiagnostics,
					     drained[i].request.objectPath,
					     "stale LoD result revision rejected");
	    continue;
	}
	update.addResult(drained[i]);
    }

    update.apply(root);
    this->lastLodMatchedResultCount = update.getMatchedResultCount();
    this->lastLodAppliedResultCount = update.getAppliedResultCount();
    this->lastLodRejectedResultCount += update.getRejectedResultCount();
    this->lastLodUnmatchedResultCount = update.getUnmatchedResultCount();
    if (update.getDiagnostics().getLength() > 0) {
	if (this->lastLodDiagnostics.getLength() > 0)
	    this->lastLodDiagnostics += "\n";
	this->lastLodDiagnostics += update.getDiagnostics();
    }
    this->lodResultsPending.store(
	service->queuedResultCountForDiagnostics() > 0 ? 1 : 0);

    if (this->lastLodAppliedResultCount > 0) {
	this->requestRender("lod-result");
	(void)this->enforceMeshResidencyBudget();
    }

    return size_to_int_saturated(
	       static_cast<size_t>(this->lastLodAppliedResultCount));
}

uint64_t
BRLObolViewController::getLodViewRevision(void) const
{
    return this->lodViewRevision;
}

void
BRLObolViewController::setLodPolicyRevision(uint64_t revision)
{
    uint64_t newRevision = revision == 0 ? 1 : revision;
    if (this->lodPolicyRevision == newRevision)
	return;
    this->lodPolicyRevision = newRevision;
    this->requestRender("lod-policy");
}

uint64_t
BRLObolViewController::getLodPolicyRevision(void) const
{
    return this->lodPolicyRevision;
}

unsigned int
BRLObolViewController::getLastLodVisitedMeshCount(void) const
{
    return this->lastLodVisitedMeshCount;
}

unsigned int
BRLObolViewController::getLastLodSubmittedTaskCount(void) const
{
    return this->lastLodSubmittedTaskCount;
}

unsigned int
BRLObolViewController::getLastLodSkippedMeshCount(void) const
{
    return this->lastLodSkippedMeshCount;
}

size_t
BRLObolViewController::getLastLodResultCount(void) const
{
    return this->lastLodResultCount;
}

unsigned int
BRLObolViewController::getLastLodMatchedResultCount(void) const
{
    return this->lastLodMatchedResultCount;
}

unsigned int
BRLObolViewController::getLastLodAppliedResultCount(void) const
{
    return this->lastLodAppliedResultCount;
}

unsigned int
BRLObolViewController::getLastLodRejectedResultCount(void) const
{
    return this->lastLodRejectedResultCount;
}

unsigned int
BRLObolViewController::getLastLodUnmatchedResultCount(void) const
{
    return this->lastLodUnmatchedResultCount;
}

const SbString &
BRLObolViewController::getLastLodDiagnostics(void) const
{
    return this->lastLodDiagnostics;
}

size_t
BRLObolViewController::getActiveLodMeshPayloadCount(void) const
{
    return this->viewLodState ? this->viewLodState->meshPayloadCount() : 0;
}

size_t
BRLObolViewController::getActiveLodProxyPayloadCount(int proxyKind) const
{
    return this->viewLodState ?
	   this->viewLodState->proxyPayloadCount(proxyKind) : 0;
}

SoBRLSceneController *
BRLObolViewController::getSceneController(void)
{
    return &this->sceneController;
}

const SoBRLSceneController *
BRLObolViewController::getSceneController(void) const
{
    return &this->sceneController;
}

BRLObolFeatureStore &
BRLObolViewController::features(void)
{
    return *this->featureStore;
}

const BRLObolFeatureStore &
BRLObolViewController::features(void) const
{
    return *this->featureStore;
}

BRLObolPolygonStore &
BRLObolViewController::polygons(void)
{
    return *this->polygonStore;
}

const BRLObolPolygonStore &
BRLObolViewController::polygons(void) const
{
    return *this->polygonStore;
}

BRLObolSelectionStore &
BRLObolViewController::selection(void)
{
    return *this->selectionStore;
}

const BRLObolSelectionStore &
BRLObolViewController::selection(void) const
{
    return *this->selectionStore;
}

SoViewport *
BRLObolViewController::getViewport(void)
{
    return this->viewport;
}

const SoViewport *
BRLObolViewController::getViewport(void) const
{
    return this->viewport;
}

SoRenderManager *
BRLObolViewController::getRenderManager(void)
{
    return this->renderManager;
}

const SoRenderManager *
BRLObolViewController::getRenderManager(void) const
{
    return this->renderManager;
}

static int
find_edit_preview_child(SoGroup *group, const char *previewId)
{
    if (!group || !previewId)
	return -1;
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *node = group->getChild(i);
	if (!node || !node->isOfType(SoBRLEditPreview::getClassTypeId()))
	    continue;
	SoBRLEditPreview *preview = static_cast<SoBRLEditPreview *>(node);
	if (strcmp(preview->previewId.getValue().getString(), previewId) == 0)
	    return i;
    }
    return -1;
}

static int
find_hud_label_overlay_child(SoGroup *group, const char *labelId)
{
    if (!group || !labelId)
	return -1;
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *node = group->getChild(i);
	if (!node || !node->isOfType(SoBRLHUDLabelOverlay::getClassTypeId()))
	    continue;
	SoBRLHUDLabelOverlay *label = static_cast<SoBRLHUDLabelOverlay *>(node);
	if (strcmp(label->labelId.getValue().getString(), labelId) == 0)
	    return i;
    }
    return -1;
}

static int
find_line_layer_overlay_child(SoGroup *group, const char *overlayId)
{
    if (!group || !overlayId)
	return -1;
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *node = group->getChild(i);
	if (!node || !node->isOfType(SoBRLLineLayerOverlay::getClassTypeId()))
	    continue;
	SoBRLLineLayerOverlay *overlay = static_cast<SoBRLLineLayerOverlay *>(node);
	if (strcmp(overlay->overlayId.getValue().getString(), overlayId) == 0)
	    return i;
    }
    return -1;
}

int
BRLObolViewController::replaceEditPreview(const char *previewId,
	const char *identity,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	uint32_t sourceRevision,
	uint32_t inputsRevision)
{
    return this->replaceEditPreviewWithIntent(previewId, identity, NULL, NULL,
	    points, commands, count, sourceRevision, inputsRevision);
}

int
BRLObolViewController::replaceEditPreviewWithIntent(const char *previewId,
	const char *identity,
	const char *editIntentId,
	const char *editIntentRole,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	uint32_t sourceRevision,
	uint32_t inputsRevision)
{
    if (!previewId || !previewId[0] || !points || !commands || count <= 0)
	return -1;

    SoNode *root = this->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(root);
    const int childIndex = find_edit_preview_child(group, previewId);
    SoBRLEditPreview *preview = NULL;
    if (childIndex >= 0)
	preview = static_cast<SoBRLEditPreview *>(group->getChild(childIndex));
    else
	preview = new SoBRLEditPreview;

    if (sourceRevision == 0)
	sourceRevision = preview->sourceRevision.getValue() + 1;
    if (inputsRevision == 0)
	inputsRevision = preview->inputsRevision.getValue() + 1;

    preview->previewId = previewId;
    preview->setEditIntent(editIntentId ? editIntentId : "",
			   editIntentRole ? editIntentRole : "preview");
    preview->sourceRevision = sourceRevision;
    preview->inputsRevision = inputsRevision;
    SoBRLVListShape *shape = preview->setLineSet(
				 (identity && identity[0]) ? identity : previewId,
				 points, commands, count);
    if (!shape)
	return -1;

    if (childIndex < 0)
	group->addChild(preview);

    this->requestRender("edit-preview");
    return 1;
}

int
BRLObolViewController::removeEditPreview(const char *previewId)
{
    if (!previewId || !previewId[0])
	return 0;

    SoNode *root = this->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(root);
    const int childIndex = find_edit_preview_child(group, previewId);
    if (childIndex < 0)
	return 0;

    group->removeChild(childIndex);
    this->requestRender("edit-preview");
    return 1;
}

int
BRLObolViewController::replaceLineLayerOverlay(const char *overlayId,
	const struct bg_line_layer_builder *builder,
	uint32_t sourceId,
	SbBool selectable)
{
    if (!overlayId || !overlayId[0])
	return -1;

    SoNode *root = this->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(root);
    const int childIndex = find_line_layer_overlay_child(group, overlayId);
    if (!builder) {
	if (childIndex >= 0) {
	    group->removeChild(childIndex);
	    this->requestRender("line-layer-overlay");
	}
	return 0;
    }

    SoBRLLineLayerOverlay *overlay = new SoBRLLineLayerOverlay;
    overlay->overlayId = overlayId;
    overlay->sourceId = sourceId;
    overlay->selectable = selectable;
    const int realized = overlay->rebuildGeometry(builder);

    if (childIndex >= 0)
	group->replaceChild(childIndex, overlay);
    else
	group->addChild(overlay);

    this->requestRender("line-layer-overlay");
    return realized;
}

int
BRLObolViewController::replaceHUDLabelOverlay(const char *labelId,
	const char *text,
	const SbVec2f &position,
	const SbColor &color,
	float fontSize,
	uint32_t sourceId)
{
    if (!labelId || !labelId[0])
	return -1;
    if (!text || !text[0])
	return this->removeHUDLabelOverlay(labelId);

    SoNode *root = this->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(root);
    const int childIndex = find_hud_label_overlay_child(group, labelId);
    SoBRLHUDLabelOverlay *label = childIndex >= 0 ?
				  static_cast<SoBRLHUDLabelOverlay *>(group->getChild(childIndex)) :
				  new SoBRLHUDLabelOverlay;

    label->labelId = labelId;
    label->sourceId = sourceId;
    label->text = text;
    label->position = position;
    label->color = color;
    label->fontSize = fontSize;
    label->visible = TRUE;
    label->rebuildGeometry();

    if (childIndex < 0)
	group->addChild(label);

    this->requestRender("hud-label-overlay");
    return 1;
}

int
BRLObolViewController::removeHUDLabelOverlay(const char *labelId)
{
    if (!labelId || !labelId[0])
	return 0;

    SoNode *root = this->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoGroup *group = static_cast<SoGroup *>(root);
    const int childIndex = find_hud_label_overlay_child(group, labelId);
    if (childIndex < 0)
	return 0;

    group->removeChild(childIndex);
    this->requestRender("hud-label-overlay");
    return 1;
}

SoGroup *
BRLObolViewController::findGroup(const char *groupPath) const
{
    return this->sceneController.findGroup(groupPath);
}

SoGroup *
BRLObolViewController::ensureGroup(const char *groupPath)
{
    const uint64_t revision = this->sceneController.getStructuralRevision();
    SoGroup *group = this->sceneController.ensureGroup(groupPath);
    if (group && this->sceneController.getStructuralRevision() != revision)
	this->requestRender("scene-group");
    return group;
}

int
BRLObolViewController::setGroupDrawIntent(const char *groupPath,
	const char *intentPath,
	int drawMode,
	int fallbackDrawMode,
	SbBool overlayIntent,
	uint32_t revalidationRevision)
{
    const int changed = this->sceneController.setGroupDrawIntent(groupPath,
			intentPath, drawMode, fallbackDrawMode, overlayIntent,
			revalidationRevision);
    if (changed > 0)
	this->requestRender("scene-group");
    return changed;
}

int
BRLObolViewController::setGroupDisplayState(const char *groupPath,
	SbBool visible,
	SbBool selected,
	SbBool highlighted,
	int lineStyle,
	int lineWidth,
	float transparency,
	SbBool colorOverride,
	const SbColor &color,
	SbBool materialColorValid,
	const SbColor &materialColor,
	uint32_t materialRevision)
{
    const int changed = this->sceneController.setGroupDisplayState(
			    groupPath, visible, selected, highlighted, lineStyle,
			    lineWidth, transparency, colorOverride, color,
			    materialColorValid, materialColor, materialRevision);
    if (changed > 0)
	this->requestRender("scene-group");
    return changed;
}

int
BRLObolViewController::renameGroup(const char *groupPath,
				   const char *newLeafName)
{
    const int changed =
	this->sceneController.renameGroup(groupPath, newLeafName);
    if (changed > 0)
	this->requestRender("scene-group");
    return changed;
}

int
BRLObolViewController::appendChildToGroup(const char *groupPath,
	SoNode *child)
{
    const int changed =
	this->sceneController.appendChildToGroup(groupPath, child);
    if (changed > 0)
	this->requestRender("scene-group");
    return changed;
}

int
BRLObolViewController::removeChildFromGroup(const char *groupPath,
	SoNode *child)
{
    const int changed =
	this->sceneController.removeChildFromGroup(groupPath, child);
    if (changed > 0)
	this->requestRender("scene-group");
    return changed;
}

int
BRLObolViewController::eraseGroupSubpath(const char *parentGroupPath,
	const char *subpath)
{
    const int changed =
	this->sceneController.eraseGroupSubpath(parentGroupPath, subpath);
    if (changed > 0)
	this->requestRender("scene-group");
    return changed;
}

int
BRLObolViewController::removeGroup(const char *groupPath)
{
    const int removed = this->sceneController.removeGroup(groupPath);
    if (removed > 0)
	this->requestRender("scene-group");
    return removed;
}

int
BRLObolViewController::clearGroup(const char *groupPath)
{
    const int removed = this->sceneController.clearGroup(groupPath);
    if (removed > 0)
	this->requestRender("scene-group");
    return removed;
}

int
BRLObolViewController::getGroupChildCount(const char *groupPath) const
{
    return this->sceneController.getGroupChildCount(groupPath);
}

int
BRLObolViewController::getGroupDescendantGroupCount(
    const char *groupPath) const
{
    return this->sceneController.getGroupDescendantGroupCount(groupPath);
}

int
BRLObolViewController::getGroupDatabaseSourceCount(
    const char *groupPath) const
{
    return this->sceneController.getGroupDatabaseSourceCount(groupPath);
}

SoNode *
BRLObolViewController::findShape(const char *shapePath) const
{
    return this->sceneController.findShape(shapePath);
}

SoGroup *
BRLObolViewController::findShapeParent(const char *shapePath) const
{
    return this->sceneController.findShapeParent(shapePath);
}

int
BRLObolViewController::moveShapeToGroup(const char *shapePath,
					const char *groupPath)
{
    const int changed =
	this->sceneController.moveShapeToGroup(shapePath, groupPath);
    if (changed > 0)
	this->requestRender("scene-shape");
    return changed;
}

int
BRLObolViewController::removeShape(const char *shapePath)
{
    const int removed = this->sceneController.removeShape(shapePath);
    if (removed > 0)
	this->requestRender("scene-shape");
    return removed;
}

int
BRLObolViewController::setShapeDrawState(const char *shapePath,
	int drawMode,
	SbBool databaseIntent,
	SbBool overlayIntent,
	SbBool hudIntent)
{
    const int changed = this->sceneController.setShapeDrawState(shapePath,
			drawMode, databaseIntent, overlayIntent, hudIntent);
    if (changed > 0)
	this->requestRender("scene-shape");
    return changed;
}

int
BRLObolViewController::setShapeDisplayState(const char *shapePath,
	SbBool visible,
	SbBool selected,
	SbBool highlighted,
	int lineStyle,
	int lineWidth,
	float transparency,
	SbBool colorOverride,
	const SbColor &color,
	SbBool materialColorValid,
	const SbColor &materialColor,
	uint32_t materialRevision)
{
    const int changed = this->sceneController.setShapeDisplayState(
			    shapePath, visible, selected, highlighted, lineStyle, lineWidth,
			    transparency, colorOverride, color, materialColorValid,
			    materialColor, materialRevision);
    if (changed > 0)
	this->requestRender("scene-shape");
    return changed;
}

int
BRLObolViewController::setShapePlacementState(const char *shapePath,
	SbBool drawMatrixValid,
	const SbMatrix &drawMatrix,
	SbBool drawCenterValid,
	const SbVec3f &drawCenter,
	SbBool drawSizeValid,
	float drawSize)
{
    const int changed = this->sceneController.setShapePlacementState(
			    shapePath, drawMatrixValid, drawMatrix, drawCenterValid,
			    drawCenter, drawSizeValid, drawSize);
    if (changed > 0)
	this->requestRender("scene-shape");
    return changed;
}

int
BRLObolViewController::setShapeSourceState(const char *shapePath,
	const char *ownerSourcePath,
	uint32_t ownerSourceRevision,
	uint32_t ownerInputsRevision,
	uint32_t ownerViewRevision,
	uint32_t ownerRealizedRevision,
	uint32_t ownerRealizedSourceRevision,
	uint32_t ownerRealizedInputsRevision,
	uint32_t ownerRealizedViewRevision,
	int ownerRealizationStatus,
	const char *ownerRealizationDiagnostic,
	const char *ownerRealizationIdentity,
	SbBool ownerSourceStale,
	uint32_t ownerStaleReason)
{
    const int changed = this->sceneController.setShapeSourceState(
			    shapePath, ownerSourcePath, ownerSourceRevision,
			    ownerInputsRevision, ownerViewRevision, ownerRealizedRevision,
			    ownerRealizedSourceRevision, ownerRealizedInputsRevision,
			    ownerRealizedViewRevision, ownerRealizationStatus,
			    ownerRealizationDiagnostic, ownerRealizationIdentity,
			    ownerSourceStale, ownerStaleReason);
    if (changed > 0)
	this->requestRender("scene-shape");
    return changed;
}

int
BRLObolViewController::replaceDatabaseSource(const char *sourcePath,
	struct db_i *dbip,
	int drawMode,
	uint32_t sourceRevision)
{
    int changed = this->sceneController.replaceDatabaseSource(sourcePath,
		  dbip, drawMode, sourceRevision);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BRLObolViewController::replaceDatabaseSourceInstance(
    const char *sourceInstanceKey,
    const char *sourcePath,
    struct db_i *dbip,
    int drawMode,
    uint32_t sourceRevision)
{
    int changed = this->sceneController.replaceDatabaseSourceInstance(
		      sourceInstanceKey, sourcePath, dbip, drawMode, sourceRevision);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BRLObolViewController::setDatabaseSourceState(const char *sourcePath,
	SbBool sourceRevisionValid,
	uint32_t sourceRevision,
	uint32_t inputsRevision,
	SbBool visible,
	SbBool highlighted,
	int lineStyle,
	int lineWidth,
	float transparency,
	SbBool colorOverride,
	const SbColor &color,
	SbBool materialColorValid,
	const SbColor &materialColor,
	uint32_t materialRevision)
{
    const int changed = this->sceneController.setDatabaseSourceState(
			    sourcePath, sourceRevisionValid, sourceRevision, inputsRevision,
			    visible, highlighted, lineStyle, lineWidth, transparency,
			    colorOverride, color, materialColorValid, materialColor,
			    materialRevision);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BRLObolViewController::setDatabaseSourceInstanceState(
    const char *sourceInstanceKey,
    SbBool sourceRevisionValid,
    uint32_t sourceRevision,
    uint32_t inputsRevision,
    SbBool visible,
    SbBool highlighted,
    int lineStyle,
    int lineWidth,
    float transparency,
    SbBool colorOverride,
    const SbColor &color,
    SbBool materialColorValid,
    const SbColor &materialColor,
    uint32_t materialRevision)
{
    const int changed = this->sceneController.setDatabaseSourceInstanceState(
			    sourceInstanceKey, sourceRevisionValid, sourceRevision,
			    inputsRevision, visible, highlighted, lineStyle, lineWidth,
			    transparency, colorOverride, color, materialColorValid,
			    materialColor, materialRevision);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BRLObolViewController::setDatabaseSourceDisplayPatch(const char *sourcePath,
	const BRLObolDatabaseSourceDisplayPatch &patch)
{
    const int changed = this->sceneController.setDatabaseSourceDisplayPatch(
			    sourcePath, patch);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BRLObolViewController::setDatabaseSourceInstanceDisplayPatch(
    const char *sourceInstanceKey,
    const BRLObolDatabaseSourceDisplayPatch &patch)
{
    const int changed =
	this->sceneController.setDatabaseSourceInstanceDisplayPatch(
	    sourceInstanceKey, patch);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BRLObolViewController::setDatabaseSourceDisplayName(const char *sourcePath,
	const char *displayName)
{
    const int changed = this->sceneController.setDatabaseSourceDisplayName(
			    sourcePath, displayName);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BRLObolViewController::setDatabaseSourceInstanceDisplayName(
    const char *sourceInstanceKey,
    const char *displayName)
{
    const int changed =
	this->sceneController.setDatabaseSourceInstanceDisplayName(
	    sourceInstanceKey, displayName);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BRLObolViewController::setDatabaseSourceBoundsState(const char *sourcePath,
	SbBool boundsValid,
	const SbVec3f &boundsMin,
	const SbVec3f &boundsMax)
{
    const int changed = this->sceneController.setDatabaseSourceBoundsState(
			    sourcePath, boundsValid, boundsMin, boundsMax);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BRLObolViewController::setDatabaseSourceInstanceBoundsState(
    const char *sourceInstanceKey,
    SbBool boundsValid,
    const SbVec3f &boundsMin,
    const SbVec3f &boundsMax)
{
    const int changed =
	this->sceneController.setDatabaseSourceInstanceBoundsState(
	    sourceInstanceKey, boundsValid, boundsMin, boundsMax);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BRLObolViewController::setDatabaseSourceMaterialPolicy(
    const char *sourcePath,
    int materialPolicy)
{
    const int changed = this->sceneController.setDatabaseSourceMaterialPolicy(
			    sourcePath, materialPolicy);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BRLObolViewController::setDatabaseSourceInstanceMaterialPolicy(
    const char *sourceInstanceKey,
    int materialPolicy)
{
    const int changed =
	this->sceneController.setDatabaseSourceInstanceMaterialPolicy(
	    sourceInstanceKey, materialPolicy);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BRLObolViewController::markDatabaseSourceStale(const char *sourcePath,
	uint32_t staleReason)
{
    const int changed = this->sceneController.markDatabaseSourceStale(
			    sourcePath, staleReason);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BRLObolViewController::markDatabaseSourceInstanceStale(
    const char *sourceInstanceKey,
    uint32_t staleReason)
{
    const int changed = this->sceneController.markDatabaseSourceInstanceStale(
			    sourceInstanceKey, staleReason);
    if (changed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return changed;
}

int
BRLObolViewController::removeDatabaseSource(const char *sourcePath)
{
    int removed = this->sceneController.removeDatabaseSource(sourcePath);
    if (removed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return removed;
}

int
BRLObolViewController::removeDatabaseSourceInstance(
    const char *sourceInstanceKey)
{
    int removed = this->sceneController.removeDatabaseSourceInstance(
		      sourceInstanceKey);
    if (removed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return removed;
}

int
BRLObolViewController::moveDatabaseSourceToGroup(const char *sourcePath,
	const char *groupPath)
{
    int moved = this->sceneController.moveDatabaseSourceToGroup(sourcePath,
		groupPath);
    if (moved > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return moved;
}

int
BRLObolViewController::moveDatabaseSourceInstanceToGroup(
    const char *sourceInstanceKey,
    const char *groupPath)
{
    int moved = this->sceneController.moveDatabaseSourceInstanceToGroup(
		    sourceInstanceKey, groupPath);
    if (moved > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return moved;
}

int
BRLObolViewController::clearDatabaseSources(void)
{
    int removed = this->sceneController.clearDatabaseSources();
    if (removed > 0) {
	this->clearRtPickCaches();
	this->requestRender("database-source");
    }
    return removed;
}

SoBRLDatabaseSource *
BRLObolViewController::getDatabaseSource(int index) const
{
    return this->sceneController.getDatabaseSource(index);
}

int
BRLObolViewController::getDatabaseSourceCount(void) const
{
    return this->sceneController.getDatabaseSourceCount();
}

SoBRLDatabaseSource *
BRLObolViewController::findDatabaseSourceInstance(
    const char *sourceInstanceKey) const
{
    return this->sceneController.findDatabaseSourceInstance(
	       sourceInstanceKey);
}

SbBool
BRLObolViewController::getDatabaseSourceSummary(int index,
	BRLObolDatabaseSourceSummary &summary) const
{
    return this->sceneController.getDatabaseSourceSummary(index, summary);
}

void
BRLObolViewController::syncRenderManager(void)
{
    this->renderManager->setSceneGraph(this->viewport->getRoot());
    this->renderManager->setCamera(this->activeCamera);
    this->renderManager->setViewportRegion(this->viewportRegion);
}

void
BRLObolViewController::advanceLodViewRevision(void)
{
    this->lodViewRevision++;
    if (this->lodViewRevision == 0)
	this->lodViewRevision++;
}

void
BRLObolViewController::advanceLodPolicyRevision(void)
{
    this->lodPolicyRevision++;
    if (this->lodPolicyRevision == 0)
	this->lodPolicyRevision++;
}

void
BRLObolViewController::syncLodViewSignature(SbBool advanceOnChange)
{
    struct rt_view_info view = RT_VIEW_INFO_INIT;
    SbBool haveCamera = this->getRtViewInfo(&view);
    SbString signature = controller_lod_view_signature(view, haveCamera);

    if (strcmp(this->lodViewSignature.getString(), signature.getString()) == 0)
	return;

    this->lodViewSignature = signature;
    if (advanceOnChange)
	this->advanceLodViewRevision();
}
