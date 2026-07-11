/*                    C A D _ A S S E M B L Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "./cad_assembly_private.h"
#include "brlobol/database_source.h"

#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/misc/SoState.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoGroup.h>

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
    batchValid(FALSE)
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
    this->batchValid = FALSE;
    this->batchedSources.clear();
}

SbBool
SoBRLCadRenderBatch::syncBatch(void)
{
    if (this->getNumChildren() != 1 || !this->getChild(0))
	return FALSE;

    SoNode *root = this->sourceRoot ? this->sourceRoot : this->getChild(0);
    std::vector<SoBRLDatabaseSource *> sources;
    cad_batch_collect_sources(root, sources);
    uint64_t signature = 1469598103934665603ULL;
    for (SoBRLDatabaseSource *source : sources) {
	signature ^= static_cast<uint64_t>(reinterpret_cast<uintptr_t>(source));
	signature *= 1099511628211ULL;
	signature ^= source->cadBatchRevision;
	signature *= 1099511628211ULL;
    }
    signature ^= static_cast<uint64_t>(sources.size());
    if (signature == this->cachedSourceSignature)
	return this->batchValid;

    BRLObolCadBatchBuildState state;
    state.assembly = this->assembly;
    this->batchedSources.clear();
    this->assembly->beginUpdate();
    this->assembly->clear();
    this->assembly->clearSemanticMap();
    for (SoBRLDatabaseSource *source : sources) {
	if (source && source->appendCadRenderBatch(&state))
	    this->batchedSources.insert(source);
    }
    this->assembly->upsertParts(state.parts);
    this->assembly->upsertInstances(state.instances);
    this->assembly->setHiddenInstances(state.hiddenInstances);
    this->assembly->setSelectedInstances(state.selectedInstances);
    this->assembly->setUnpickableInstances(state.unpickableInstances);
    if (state.shadedCount > 0 && state.wireCount > 0)
	this->assembly->drawMode = SoCADAssembly::SHADED_WITH_EDGES;
    else if (state.shadedCount > 0)
	this->assembly->drawMode = SoCADAssembly::SHADED;
    else
	this->assembly->drawMode = SoCADAssembly::WIREFRAME;
    this->assembly->endUpdate();

    this->cachedSourceSignature = signature;
    this->batchValid = this->batchedSources.size() >= 32 &&
	this->assembly->instanceCount() > 0 ? TRUE : FALSE;
    return this->batchValid;
}

void
SoBRLCadRenderBatch::renderBatch(SoGLRenderAction *action)
{
    if (!action)
	return;
    if (!this->syncBatch()) {
	inherited::GLRenderBelowPath(action);
	return;
    }

    SoState *state = action->getState();
    if (state)
	state->push();
    SoGroup *viewportRoot = dynamic_cast<SoGroup *>(this->getChild(0));
    if (viewportRoot) {
	for (int i = 0; i < viewportRoot->getNumChildren(); i++) {
	    SoNode *child = viewportRoot->getChild(i);
	    if (!dynamic_cast<SoCamera *>(child))
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
