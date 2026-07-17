/*                B E D I T P R E V I E W . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BEditPreview.h */

#ifndef BOBOL_BEDITPREVIEW_H
#define BOBOL_BEDITPREVIEW_H

#include "BObol/BDefines.h"

#include <Inventor/SbMatrix.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFEnum.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFUInt32.h>
#include <Inventor/nodes/SoSeparator.h>

class SoBRLVListShape;
class SoFieldSensor;
class SoSensor;

class BOBOL_EXPORT SoBRLEditPreview : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLEditPreview);

public:
    enum PreviewStatus {
	EMPTY = 0,
	CURRENT = 1,
	STALE = 2,
	FAILED = 3
    };

    SoSFString previewId;
    SoSFString editIntentId;
    SoSFString editIntentRole;
    SoSFUInt32 sourceRevision;
    SoSFUInt32 inputsRevision;
    SoSFUInt32 realizedSourceRevision;
    SoSFUInt32 realizedInputsRevision;
    SoSFEnum previewStatus;
    SoSFBool stale;

    SoBRLEditPreview(void);
    static void initClass(void);

    void setEditIntent(const SbString &id, const SbString &role);
    void markSourceRevision(uint32_t revision);
    void markInputsRevision(uint32_t revision);
    SbBool needsRealization(void) const;

    void clearPreview(void);
    SoBRLVListShape *setLineSet(const SbString &identity,
	    const SbVec3f *points,
	    const int32_t *commands,
	    int count);
    SoBRLVListShape *setTransformedLineSet(const SbString &identity,
	    const SbMatrix &matrix,
	    const SbVec3f *points,
	    const int32_t *commands,
	    int count);

protected:
    virtual ~SoBRLEditPreview(void);

private:
    static void fieldSensorCB(void *data, SoSensor *sensor);
    void attachFieldSensors(void);
    void detachFieldSensors(void);
    void markStale(void);
    SoBRLVListShape *appendLineSet(const SbString &identity,
	    const SbVec3f *points,
	    const int32_t *commands,
	    int count);
    void markCurrent(void);

    SoFieldSensor *previewIdSensor;
    SoFieldSensor *editIntentIdSensor;
    SoFieldSensor *editIntentRoleSensor;
    SoFieldSensor *sourceRevisionSensor;
    SoFieldSensor *inputsRevisionSensor;
};

#endif /* BOBOL_BEDITPREVIEW_H */
