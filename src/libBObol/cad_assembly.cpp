/*                    C A D _ A S S E M B L Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "./cad_assembly_private.h"
#include "cad_publication_private.h"
#include "BObol/BDatabaseSource.h"
#include "BObol/BViewLod.h"

#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/SbName.h>
#include <Inventor/misc/SoState.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoGroup.h>

#include <Obol/cad/SoCADViewState.h>

#include <map>
#include <optional>
#include <vector>

SO_NODE_SOURCE(SoBRLCadAssembly);
SO_NODE_SOURCE(SoBRLCadRenderBatch);

static thread_local const std::unordered_set<const SoBRLDatabaseSource *> *
    activeCadBatchSources = NULL;

namespace {

class CadViewPolicyScope {
public:
    CadViewPolicyScope(SoAction *action, Obol::CadDrawMode drawMode,
	Obol::CadPickMode pickMode) :
	state(action ? action->getState() : NULL)
    {
	if (!this->state)
	    return;
	this->state->push();
	if (!this->state->isElementEnabled(
		SoCADViewStateElement::getClassStackIndex()))
	    return;
	Obol::CadViewState view = SoCADViewStateElement::get(this->state);
	view.drawMode = drawMode;
	view.pickMode = pickMode;
	SoCADViewStateElement::set(this->state, view);
    }

    ~CadViewPolicyScope()
    {
	if (this->state)
	    this->state->pop();
    }

    CadViewPolicyScope(const CadViewPolicyScope &) = delete;
    CadViewPolicyScope &operator=(const CadViewPolicyScope &) = delete;

private:
    SoState *state;
};

} // namespace

int
bobol_cad_batch_source_suppressed(const SoBRLDatabaseSource *source)
{
    return (source && activeCadBatchSources &&
	activeCadBatchSources->find(source) != activeCadBatchSources->end()) ? 1 : 0;
}

static void
cad_batch_collect_sources(SoNode *node,
	std::vector<SoBRLDatabaseSource *> &sources)
{
    if (!node)
	return;
    SoBRLDatabaseSource *source = dynamic_cast<SoBRLDatabaseSource *>(node);
    if (source) {
	sources.push_back(source);
	return;
    }
    SoGroup *group = dynamic_cast<SoGroup *>(node);
    if (!group)
	return;
    for (int i = 0; i < group->getNumChildren(); i++)
	cad_batch_collect_sources(group->getChild(i), sources);
}

static uint64_t
cad_batch_instance_state_signature(const std::vector<Obol::InstanceId> &ids)
{
    uint64_t signature = 1469598103934665603ULL;
    for (const Obol::InstanceId &id : ids) {
	signature ^= id.w0;
	signature *= 1099511628211ULL;
	signature ^= id.w1;
	signature *= 1099511628211ULL;
    }
    return signature ? signature : 1;
}

SoBRLCadAssembly::SoBRLCadAssembly(void) :
    presentationDrawModeValue(Obol::CadDrawMode::Wireframe),
    presentationPickModeValue(Obol::CadPickMode::Automatic)
{
    SO_NODE_CONSTRUCTOR(SoBRLCadAssembly);
}

SoBRLCadAssembly::~SoBRLCadAssembly(void)
{
}

void
SoBRLCadAssembly::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLCadAssembly, SoCADAssembly, "SoCADAssembly");
}

void
SoBRLCadAssembly::clearSemanticMap(void)
{
    this->semantics.clear();
}

void
SoBRLCadAssembly::setInstanceSemantic(Obol::InstanceId id,
    const InstanceSemantic &semantic)
{
    this->semantics[id] = semantic;
}

void
SoBRLCadAssembly::lastUncollapsedStructuralProxyOccurrenceKeys(
    std::vector<SbString> &occurrenceKeys) const
{
    occurrenceKeys.clear();
    const std::vector<Obol::InstanceId> instances =
	SoCADAssembly::lastUncollapsedStructuralProxyInstances();
    occurrenceKeys.reserve(instances.size());
    for (const Obol::InstanceId instance : instances) {
	const auto semantic = this->semantics.find(instance);
	if (semantic == this->semantics.end() ||
	    semantic->second.sourceInstanceKey.getLength() == 0)
	    continue;
	occurrenceKeys.push_back(semantic->second.sourceInstanceKey);
    }
}

bool
SoBRLCadAssembly::lastStructuralProxyOccurrenceKeysAbovePixels(
    float pixels, std::vector<SbString> &occurrenceKeys) const
{
    occurrenceKeys.clear();
    const std::vector<Obol::InstanceId> instances =
	this->lastStructuralProxyInstancesAbovePixels(pixels);
    occurrenceKeys.reserve(instances.size());
    for (const Obol::InstanceId &instance : instances) {
	const auto semantic = this->semantics.find(instance);
	if (semantic == this->semantics.end() ||
		semantic->second.sourceInstanceKey.getLength() == 0)
	    return false;
	occurrenceKeys.push_back(semantic->second.sourceInstanceKey);
    }
    return true;
}

void
SoBRLCadAssembly::reserveCompactPresentationCapacity(
    size_t expectedOccurrences)
{
    if (!expectedOccurrences)
	return;
    this->reserveStreamingCapacity(expectedOccurrences);
    this->semantics.reserve(expectedOccurrences);
    this->compactInstancePresentations.reserve(expectedOccurrences);
    this->compactActivePartReferences.reserve(expectedOccurrences);
    this->compactPartChannels.reserve(expectedOccurrences);
    this->compactLodParts.reserve(expectedOccurrences);
}

void
SoBRLCadAssembly::render(SoGLRenderAction *action)
{
    this->GLRender(action);
}

void
SoBRLCadAssembly::setPresentationDrawMode(Obol::CadDrawMode mode)
{
    if (this->presentationDrawModeValue == mode)
	return;
    this->presentationDrawModeValue = mode;
    this->touch();
}

void
SoBRLCadAssembly::setPresentationPickMode(Obol::CadPickMode mode)
{
    if (this->presentationPickModeValue == mode)
	return;
    this->presentationPickModeValue = mode;
    this->touch();
}

void
SoBRLCadAssembly::GLRender(SoGLRenderAction *action)
{
    CadViewPolicyScope policy(action, this->presentationDrawModeValue,
	this->presentationPickModeValue);
    inherited::GLRender(action);
}

void
SoBRLCadAssembly::getBounds(SoGetBoundingBoxAction *action)
{
    this->getBoundingBox(action);
}

void
SoBRLCadAssembly::pickRay(SoRayPickAction *action)
{
    this->rayPick(action);
}

void
SoBRLCadAssembly::rayPick(SoRayPickAction *action)
{
    CadViewPolicyScope policy(action, this->presentationDrawModeValue,
	this->presentationPickModeValue);
    inherited::rayPick(action);
}

SoDetail *
SoBRLCadAssembly::createPickDetail(
    const Obol::CadPickDetailRecord &hit) const
{
    const auto found = this->semantics.find(hit.instance);
    if (found == this->semantics.end())
	return inherited::createPickDetail(hit);

    const InstanceSemantic &semantic = found->second;
    SoBRLPickDetail *detail = new SoBRLPickDetail;
    detail->setPath(semantic.path);
    detail->setSourceInstanceKey(semantic.sourceInstanceKey);
    detail->setSourceName(semantic.sourceName);
    detail->setSourceType(semantic.sourceType);
    detail->setSourceId(semantic.sourceId);
    detail->setRegionId(semantic.regionId);
    detail->setAirCode(semantic.airCode);
    detail->setMaterialId(semantic.materialId);
    detail->setLos(semantic.los);
    detail->setMaterialColor(semantic.materialColorValid,
			     semantic.materialColor);
    detail->setMaterialShader(semantic.materialShader);
    detail->setEditIntent(semantic.editIntentId,
			  semantic.editIntentRole);
    detail->setModelPoint(hit.point);
    detail->setPrimitive(semantic.primitiveKind,
	static_cast<int>(hit.primIndex0));
    return detail;
}

SoBRLCadRenderBatch::SoBRLCadRenderBatch(void) :
    assembly(new SoBRLCadAssembly),
    sourceRoot(NULL),
    cachedSourceSignature(0),
    cachedStructureSignature(0),
    cachedStyleSignature(0),
    cachedSemanticSignature(0),
    cachedHiddenSignature(0),
    cachedSelectedSignature(0),
    cachedUnpickableSignature(0),
    batchValid(FALSE),
    softwareWireMode(SoCADViewState::SOFTWARE_WIRE_AUTO)
{
    SO_NODE_CONSTRUCTOR(SoBRLCadRenderBatch);
    this->assembly->ref();
}

SoBRLCadRenderBatch::~SoBRLCadRenderBatch(void)
{
    this->setBatchSourceRoot(NULL);
    this->assembly->unref();
    this->assembly = NULL;
}

void
SoBRLCadRenderBatch::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLCadRenderBatch, SoSeparator, "Separator");
}

void
SoBRLCadRenderBatch::setBatchSourceRoot(SoNode *root)
{
    if (root == this->sourceRoot)
	return;
    if (root)
	root->ref();
    if (this->sourceRoot)
	this->sourceRoot->unref();
    this->sourceRoot = root;
    this->cachedSourceSignature = 0;
    this->cachedStructureSignature = 0;
    this->cachedStyleSignature = 0;
    this->cachedSemanticSignature = 0;
    this->cachedHiddenSignature = 0;
    this->cachedSelectedSignature = 0;
    this->cachedUnpickableSignature = 0;
    this->batchValid = FALSE;
    this->batchedSources.clear();
}

void
SoBRLCadRenderBatch::setSoftwareWireMode(int mode)
{
    if (mode < SoCADViewState::SOFTWARE_WIRE_AUTO ||
	mode > SoCADViewState::SOFTWARE_WIRE_FAST)
	mode = SoCADViewState::SOFTWARE_WIRE_AUTO;
    this->softwareWireMode = mode;
}

SbBool
SoBRLCadRenderBatch::syncBatch(const BObolViewLodState *viewState)
{
    if (this->getNumChildren() != 1 || !this->getChild(0))
	return FALSE;

    SoNode *root = this->sourceRoot ? this->sourceRoot : this->getChild(0);
    std::vector<SoBRLDatabaseSource *> sources;
    cad_batch_collect_sources(root, sources);
    /* A compact source already compiles its occurrences into one retained
     * assembly.  Cross-source batching only pays for itself when it combines
     * many independently compiled sources; decide that before collecting or
     * copying any part geometry. */
    if (sources.size() < 32) {
	this->cachedSourceSignature = 0;
	this->cachedStructureSignature = 0;
	this->cachedStyleSignature = 0;
	this->cachedSemanticSignature = 0;
	this->cachedHiddenSignature = 0;
	this->cachedSelectedSignature = 0;
	this->cachedUnpickableSignature = 0;
	this->batchValid = FALSE;
	this->batchedSources.clear();
	return FALSE;
    }
    uint64_t signature = 1469598103934665603ULL;
    uint64_t structureSignature = signature;
    uint64_t styleSignature = signature;
    uint64_t semanticSignature = signature;
    for (SoBRLDatabaseSource *source : sources) {
	signature ^= static_cast<uint64_t>(reinterpret_cast<uintptr_t>(source));
	signature *= 1099511628211ULL;
	signature ^= source->cadBatchRevisionGet();
	signature *= 1099511628211ULL;
	structureSignature ^=
	    source ? source->cadBatchStructureSignature() : 0;
	structureSignature *= 1099511628211ULL;
	styleSignature ^=
	    source ? source->cadBatchStyleSignature() : 0;
	styleSignature *= 1099511628211ULL;
	semanticSignature ^=
	    source ? source->cadBatchSemanticSignature() : 0;
	semanticSignature *= 1099511628211ULL;
	const BObolViewLodState::CadPayload *payload =
	    viewState ? viewState->findCad(source) : NULL;
	signature ^= payload ? 1ULL : 0ULL;
	signature *= 1099511628211ULL;
	structureSignature ^= payload ? 1ULL : 0ULL;
	structureSignature *= 1099511628211ULL;
	styleSignature ^= payload ? 1ULL : 0ULL;
	styleSignature *= 1099511628211ULL;
	semanticSignature ^= payload ? 1ULL : 0ULL;
	semanticSignature *= 1099511628211ULL;
    }
    signature ^= static_cast<uint64_t>(sources.size());
    structureSignature ^= static_cast<uint64_t>(sources.size());
    styleSignature ^= static_cast<uint64_t>(sources.size());
    semanticSignature ^= static_cast<uint64_t>(sources.size());
    if (signature == this->cachedSourceSignature)
	return this->batchValid;

    const SbBool structureChanged = !this->batchValid ||
	structureSignature != this->cachedStructureSignature;
    const SbBool styleChanged = !this->batchValid ||
	styleSignature != this->cachedStyleSignature;
    const SbBool semanticChanged = !this->batchValid ||
	semanticSignature != this->cachedSemanticSignature;
    const SbBool includeSemantics = structureChanged || semanticChanged;

    BObolCadBatchBuildState state;
    state.parts.reserve(sources.size());
    state.partIds.reserve(sources.size());
    state.instances.reserve(sources.size());
    this->batchedSources.clear();
    for (SoBRLDatabaseSource *source : sources) {
	if (viewState && viewState->findCad(source))
	    continue;
	if (source && source->appendCadRenderBatch(&state, structureChanged,
		includeSemantics))
	    this->batchedSources.insert(source);
    }
    if (!state.valid) {
	this->batchedSources.clear();
	this->batchValid = FALSE;
	return FALSE;
    }
    std::vector<Obol::InstanceStyleUpdate> styles;
    if (!structureChanged && styleChanged) {
	styles.reserve(state.instances.size());
	for (const Obol::InstanceUpdate &instance : state.instances) {
	    Obol::InstanceStyleUpdate style;
	    style.instance = instance.instance;
	    style.style = instance.record.style;
	    styles.push_back(style);
	}
    }
    if ((structureChanged &&
	    (!bobol_cad_validate_shared_parts(state.parts,
		"CAD render batch preflight") ||
	     !bobol_cad_validate_instances(state.instances,
		"CAD render batch preflight"))) ||
	(!structureChanged && styleChanged &&
	 !bobol_cad_validate_styles(styles, "CAD style batch preflight"))) {
	this->batchedSources.clear();
	this->batchValid = FALSE;
	return FALSE;
    }
    std::optional<SoCADAssembly::UpdateScope> updateScope;
    if (structureChanged) {
	updateScope.emplace(this->assembly->batchUpdate());
	if (!bobol_cad_replace_scene(this->assembly, state.parts,
		state.instances, "CAD render batch replacement")) {
	    this->batchedSources.clear();
	    this->batchValid = FALSE;
	    return FALSE;
	}
	this->assembly->clearSemanticMap();
	for (const auto &semantic : state.semantics)
	    this->assembly->setInstanceSemantic(
		semantic.first, semantic.second);
    } else if (styleChanged) {
	if (!bobol_cad_publish_styles(this->assembly, styles,
		"CAD style batch publication")) {
	    this->batchValid = FALSE;
	    return FALSE;
	}
    }
    const uint64_t hiddenSignature =
	cad_batch_instance_state_signature(state.hiddenInstances);
    if (structureChanged || hiddenSignature != this->cachedHiddenSignature)
	this->assembly->setHiddenInstances(state.hiddenInstances);
    const uint64_t selectedSignature =
	cad_batch_instance_state_signature(state.selectedInstances);
    if (structureChanged || selectedSignature != this->cachedSelectedSignature)
	this->assembly->setSelectedInstances(state.selectedInstances);
    const uint64_t unpickableSignature =
	cad_batch_instance_state_signature(state.unpickableInstances);
    if (structureChanged ||
	unpickableSignature != this->cachedUnpickableSignature)
	this->assembly->setUnpickableInstances(state.unpickableInstances);
    Obol::CadDrawMode drawMode = Obol::CadDrawMode::Wireframe;
    if (state.shadedCount > 0 && state.wireCount > 0)
	drawMode = Obol::CadDrawMode::ShadedWithEdges;
    else if (state.shadedCount > 0)
	drawMode = Obol::CadDrawMode::Shaded;
    this->assembly->setPresentationDrawMode(drawMode);
    if (updateScope)
	updateScope->finish();
    this->cachedSourceSignature = signature;
    this->cachedStructureSignature = structureSignature;
    this->cachedStyleSignature = styleSignature;
    this->cachedSemanticSignature = semanticSignature;
    this->cachedHiddenSignature = hiddenSignature;
    this->cachedSelectedSignature = selectedSignature;
    this->cachedUnpickableSignature = unpickableSignature;
    this->batchValid = this->batchedSources.size() >= 32 &&
	this->assembly->instanceCount() > 0 ? TRUE : FALSE;
    return this->batchValid;
}

void
SoBRLCadRenderBatch::renderBatch(SoGLRenderAction *action)
{
    if (!action)
	return;
    const BObolViewLodState *viewState =
	SoBRLViewLodElement::get(action->getState());
    if (viewState && (viewState->meshPayloadCount() > 0 ||
	    viewState->proxyPayloadCount() > 0)) {
	inherited::GLRenderBelowPath(action);
	return;
    }
    if (!this->syncBatch(viewState)) {
	inherited::GLRenderBelowPath(action);
	return;
    }

    SoState *state = action->getState();
    if (state) {
	state->push();
	Obol::CadViewState cadState = SoCADViewStateElement::get(state);
	cadState.softwareWireMode =
	    static_cast<Obol::CadSoftwareWireMode>(this->softwareWireMode);
	SoCADViewStateElement::set(state, cadState);
    }
    SoGroup *viewportRoot = dynamic_cast<SoGroup *>(this->getChild(0));
    if (viewportRoot) {
	const SbName environmentName("BObolRenderEnvironment");
	const SbName sceneLightsName("BObolSceneLights");
	for (int i = 0; i < viewportRoot->getNumChildren(); i++) {
	    SoNode *child = viewportRoot->getChild(i);
	    /* The compact renderer bypasses the source subtree, but still needs
	     * the controller's camera, depth, light, clip, and in-scene light
	     * state first.  Index order keeps BObolSceneLights after the camera
	     * so its world-space light positions render correctly. */
	    if (!dynamic_cast<SoCamera *>(child) &&
		child->getName() != environmentName &&
		child->getName() != sceneLightsName)
		continue;
	    action->pushCurPath(i, child);
	    action->traverse(child);
	    action->popCurPath();
	}
    }
    this->assembly->render(action);
    if (state)
	state->pop();
    const std::unordered_set<const SoBRLDatabaseSource *> *previous =
	activeCadBatchSources;
    activeCadBatchSources = &this->batchedSources;
    inherited::GLRenderBelowPath(action);
    activeCadBatchSources = previous;
}

void
SoBRLCadRenderBatch::GLRender(SoGLRenderAction *action)
{
    inherited::GLRender(action);
}

void
SoBRLCadRenderBatch::GLRenderBelowPath(SoGLRenderAction *action)
{
    this->renderBatch(action);
}
