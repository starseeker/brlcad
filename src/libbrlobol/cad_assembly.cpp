/*                    C A D _ A S S E M B L Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "./cad_assembly_private.h"

SO_NODE_SOURCE(SoBRLCadAssembly);

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
