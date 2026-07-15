/*                    C A D _ A S S E M B L Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "./cad_assembly_private.h"
#include "brlobol/database_source.h"
#include "brlobol/view_lod.h"

#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/SbName.h>
#include <Inventor/misc/SoState.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoGroup.h>

#include <obol/cad/SoCADViewState.h>

#include <map>
#include <vector>

SO_NODE_SOURCE(SoBRLCadAssembly);
SO_NODE_SOURCE(SoBRLCadRenderBatch);

static thread_local const std::unordered_set<const SoBRLDatabaseSource *> *
    activeCadBatchSources = NULL;

int
brlobol_cad_batch_source_suppressed(const SoBRLDatabaseSource *source)
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
cad_batch_instance_state_signature(const std::vector<obol::InstanceId> &ids)
{
    uint64_t signature = 1469598103934665603ULL;
    for (const obol::InstanceId &id : ids) {
	signature ^= id.w0;
	signature *= 1099511628211ULL;
	signature ^= id.w1;
	signature *= 1099511628211ULL;
    }
    return signature ? signature : 1;
}

SoBRLCadAssembly::SoBRLCadAssembly(void)
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
SoBRLCadAssembly::setInstanceSemantic(obol::InstanceId id,
	const InstanceSemantic &semantic)
{
    this->semantics[id] = semantic;
}

void
SoBRLCadAssembly::render(SoGLRenderAction *action)
{
    this->GLRender(action);
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

SoDetail *
SoBRLCadAssembly::createPickDetail(
    const obol::CadPickDetailRecord &hit) const
{
    std::map<obol::InstanceId, InstanceSemantic>::const_iterator found =
	this->semantics.find(hit.instance);
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
SoBRLCadRenderBatch::syncBatch(const BRLObolViewLodState *viewState)
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
    uint64_t semanticSignature = signature;
    for (SoBRLDatabaseSource *source : sources) {
	signature ^= static_cast<uint64_t>(reinterpret_cast<uintptr_t>(source));
	signature *= 1099511628211ULL;
	signature ^= source->cadBatchRevision;
	signature *= 1099511628211ULL;
	structureSignature ^=
	    source ? source->cadBatchStructureSignature() : 0;
	structureSignature *= 1099511628211ULL;
	semanticSignature ^=
	    source ? source->cadBatchSemanticSignature() : 0;
	semanticSignature *= 1099511628211ULL;
	const BRLObolViewLodState::CadPayload *payload =
	    viewState ? viewState->findCad(source) : NULL;
	signature ^= payload ? 1ULL : 0ULL;
	signature *= 1099511628211ULL;
	structureSignature ^= payload ? 1ULL : 0ULL;
	structureSignature *= 1099511628211ULL;
	semanticSignature ^= payload ? 1ULL : 0ULL;
	semanticSignature *= 1099511628211ULL;
    }
    signature ^= static_cast<uint64_t>(sources.size());
    structureSignature ^= static_cast<uint64_t>(sources.size());
    semanticSignature ^= static_cast<uint64_t>(sources.size());
    if (signature == this->cachedSourceSignature)
	return this->batchValid;

    const SbBool structureChanged = !this->batchValid ||
	structureSignature != this->cachedStructureSignature;
    const SbBool semanticChanged = !this->batchValid ||
	semanticSignature != this->cachedSemanticSignature;
    const SbBool includeSemantics = structureChanged || semanticChanged;

    BRLObolCadBatchBuildState state;
    state.assembly = this->assembly;
    state.parts.reserve(sources.size());
    state.partIds.reserve(sources.size());
    state.instances.reserve(sources.size());
    this->batchedSources.clear();
    if (structureChanged) {
	this->assembly->beginUpdate();
	this->assembly->clear();
	this->assembly->clearSemanticMap();
    }
    for (SoBRLDatabaseSource *source : sources) {
	if (viewState && viewState->findCad(source))
	    continue;
	if (source && source->appendCadRenderBatch(&state, structureChanged,
		includeSemantics))
	    this->batchedSources.insert(source);
    }
    if (structureChanged) {
	this->assembly->upsertSharedParts(state.parts);
	this->assembly->upsertInstances(state.instances);
    } else {
	std::vector<obol::InstanceStyleUpdate> styles;
	styles.reserve(state.instances.size());
	for (const obol::InstanceUpdate &instance : state.instances) {
	    obol::InstanceStyleUpdate style;
	    style.instance = instance.instance;
	    style.style = instance.record.style;
	    styles.push_back(style);
	}
	this->assembly->updateInstanceStyles(styles);
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
    if (state.shadedCount > 0 && state.wireCount > 0)
	this->assembly->drawMode = SoCADAssembly::SHADED_WITH_EDGES;
    else if (state.shadedCount > 0)
	this->assembly->drawMode = SoCADAssembly::SHADED;
    else
	this->assembly->drawMode = SoCADAssembly::WIREFRAME;
    if (structureChanged)
	this->assembly->endUpdate();

    this->cachedSourceSignature = signature;
    this->cachedStructureSignature = structureSignature;
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
    const BRLObolViewLodState *viewState =
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
	obol::CadViewState cadState = SoCADViewStateElement::get(state);
	cadState.softwareWireMode =
	    static_cast<obol::CadSoftwareWireMode>(this->softwareWireMode);
	SoCADViewStateElement::set(state, cadState);
    }
    SoGroup *viewportRoot = dynamic_cast<SoGroup *>(this->getChild(0));
    if (viewportRoot) {
	const SbName environmentName("BRLObolRenderEnvironment");
	for (int i = 0; i < viewportRoot->getNumChildren(); i++) {
	    SoNode *child = viewportRoot->getChild(i);
	    /* The compact renderer bypasses the source subtree, but still needs
	     * the controller's camera, depth, light, and clip state first. */
	    if (!dynamic_cast<SoCamera *>(child) &&
		child->getName() != environmentName)
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
