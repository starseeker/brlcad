/*              C A D _ A S S E M B L Y _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file cad_assembly_private.h
 *
 * libbrlobol-private BRL-CAD semantic layer for Obol CAD assemblies.
 */

#ifndef BRLOBOL_CAD_ASSEMBLY_PRIVATE_H
#define BRLOBOL_CAD_ASSEMBLY_PRIVATE_H

#include "brlobol/pick_detail.h"

#include <obol/cad/SoCADAssembly.h>

#include <cstdint>
#include <map>

class SoGetBoundingBoxAction;
class SoGLRenderAction;
class SoRayPickAction;

class SoBRLCadAssembly : public SoCADAssembly {
    typedef SoCADAssembly inherited;

    SO_NODE_HEADER(SoBRLCadAssembly);

public:
    struct InstanceSemantic {
	InstanceSemantic(void) :
	    sourceId(0),
	    regionId(0),
	    airCode(0),
	    materialId(0),
	    los(0),
	    materialColorValid(FALSE),
	    primitiveKind(SoBRLPickDetail::UNKNOWN)
	{
	}

	SbString path;
	SbString sourceInstanceKey;
	SbString sourceName;
	SbString sourceType;
	SbString materialShader;
	SbString editIntentId;
	SbString editIntentRole;
	SbColor materialColor;
	uint32_t sourceId;
	int regionId;
	int airCode;
	int materialId;
	int los;
	SbBool materialColorValid;
	SoBRLPickDetail::PrimitiveKind primitiveKind;
    };

    SoBRLCadAssembly(void);
    static void initClass(void);

    void clearSemanticMap(void);
    void setInstanceSemantic(obol::InstanceId id,
	const InstanceSemantic &semantic);

    void render(SoGLRenderAction *action);
    void getBounds(SoGetBoundingBoxAction *action);
    void pickRay(SoRayPickAction *action);

protected:
    ~SoBRLCadAssembly(void) override;
    SoDetail *createPickDetail(
	const obol::CadPickDetailRecord &hit) const override;

private:
    std::map<obol::InstanceId, InstanceSemantic> semantics;
};

#endif /* BRLOBOL_CAD_ASSEMBLY_PRIVATE_H */
