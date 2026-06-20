/*                D A T A B A S E _ S O U R C E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/database_source.h */

#ifndef BRLOBOL_DATABASE_SOURCE_H
#define BRLOBOL_DATABASE_SOURCE_H

#include "brlobol/defines.h"

#include <stdint.h>

#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFEnum.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFUInt32.h>
#include <Inventor/nodes/SoSeparator.h>

class SoBRLVListShape;
class SoBRLMeshShape;
class SoFieldSensor;
class SoSensor;
struct db_i;

class BRLOBOL_EXPORT SoBRLDatabaseSource : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLDatabaseSource);

public:
    enum DrawMode {
	WIREFRAME = 0,
	SHADED = 1
    };

    enum RealizationStatus {
	UNREALIZED = 0,
	REALIZED = 1,
	FAILED = 2
    };

    enum StaleReason {
	STALE_NONE = 0,
	STALE_SOURCE = 1,
	STALE_INPUTS = 2,
	STALE_VIEW = 4,
	STALE_DATABASE = 8,
	STALE_DRAW = 16,
	STALE_TESSELLATION = 32
    };

    SoSFString path;
    SoSFEnum drawMode;
    SoSFFloat tessellationAbsTol;
    SoSFFloat tessellationRelTol;
    SoSFFloat tessellationNormTol;
    SoSFUInt32 lodBotThreshold;
    SoSFUInt32 sourceRevision;
    SoSFUInt32 inputsRevision;
    SoSFUInt32 viewRevision;
    SoSFUInt32 realizedRevision;
    SoSFUInt32 realizedSourceRevision;
    SoSFUInt32 realizedInputsRevision;
    SoSFUInt32 realizedViewRevision;
    SoSFEnum realizationStatus;
    SoSFString realizationDiagnostic;
    SoSFBool stale;
    SoSFUInt32 staleReason;

    SoBRLDatabaseSource(void);
    static void initClass(void);

    /**
     * Set the BRL-CAD database used for realization.
     *
     * The caller owns the db_i pointer.  SoBRLDatabaseSource stores a borrowed
     * pointer and does not close or otherwise manage the database lifetime.
     * The database must remain valid while this source can be realized.
     */
    void setDatabase(struct db_i *dbip);
    struct db_i *getDatabase(void) const;
    void configureDatabaseSource(const char *sourcePath,
	struct db_i *database,
	int mode,
	uint32_t revision);

    void markStale(void);
    void markStale(uint32_t reason);
    SbBool needsRealization(void) const;
    SbBool realizePrototypeWireframe(void);
    SbBool realizeDatabaseWireframe(void);
    SbBool realizeDatabaseMesh(void);
    SoBRLVListShape *getRealizedShape(void) const;
    SoBRLVListShape *getRealizedShape(int index) const;
    int getRealizedShapeCount(void) const;
    SoBRLMeshShape *getRealizedMesh(void) const;
    SoBRLMeshShape *getRealizedMesh(int index) const;
    int getRealizedMeshCount(void) const;

protected:
    virtual ~SoBRLDatabaseSource(void);

private:
    static void fieldSensorCB(void *data, SoSensor *sensor);
    void attachFieldSensors(void);
    void detachFieldSensors(void);

    struct db_i *dbip;
    SoFieldSensor *pathSensor;
    SoFieldSensor *drawModeSensor;
    SoFieldSensor *tessellationAbsTolSensor;
    SoFieldSensor *tessellationRelTolSensor;
    SoFieldSensor *tessellationNormTolSensor;
    SoFieldSensor *lodBotThresholdSensor;
    SoFieldSensor *sourceRevisionSensor;
    SoFieldSensor *inputsRevisionSensor;
    SoFieldSensor *viewRevisionSensor;
};

#endif /* BRLOBOL_DATABASE_SOURCE_H */
