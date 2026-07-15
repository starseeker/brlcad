/*          M E S H _ L O D _ S U B M I T _ A C T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/mesh_lod_submit_action.h */

#ifndef BRLOBOL_MESH_LOD_SUBMIT_ACTION_H
#define BRLOBOL_MESH_LOD_SUBMIT_ACTION_H

#include "brlobol/defines.h"
#include "brlobol/lod_realization.h"
#include "bv/view.h"

#include <Inventor/SbString.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/actions/SoSubAction.h>

#include <stdint.h>
#include <stddef.h>

class BRLObolLodService;
class BRLObolViewLodState;
struct db_i;

class BRLOBOL_EXPORT SoBRLMeshLodSubmitAction : public SoAction {
    typedef SoAction inherited;

    SO_ACTION_HEADER(SoBRLMeshLodSubmitAction);

public:
    SoBRLMeshLodSubmitAction(void);
    virtual ~SoBRLMeshLodSubmitAction(void);
    static void initClass(void);

    void setService(BRLObolLodService *service);
    BRLObolLodService *getService(void) const;
    void setDatabase(struct db_i *dbip, const char *databaseId = NULL,
	    uint64_t databaseRevision = 0);
    void setViewInfo(const struct bv_view_info *info);
    const struct bv_view_info &getViewInfo(void) const;
    void setGeneration(uint64_t generation);
    uint64_t getGeneration(void) const;
    void setRevisions(uint64_t viewRevision, uint64_t policyRevision);
    void setProvider(const char *providerId, const char *providerVersion);
    void setQualityTier(int qualityTier);
    void setRefreshMissing(SbBool refreshMissing);
    void setReset(int reset);
    void setForcedLevel(int level);
    void clearForcedLevel(void);
    SbBool hasForcedLevel(void) const;
    int getForcedLevel(void) const;
    void setRequireLodBacked(SbBool requireLodBacked);
    SbBool getRequireLodBacked(void) const;
    void setProxyStages(SbBool submitAabb, SbBool submitObb);
    SbBool getSubmitAabbProxyStage(void) const;
    SbBool getSubmitObbProxyStage(void) const;
    void setViewLodState(const BRLObolViewLodState *viewState);
    const BRLObolViewLodState *getViewLodState(void) const;
    void setCompactEntryRange(size_t first, size_t count);
    size_t getCompactEntryNext(void) const;
    size_t getCompactEntryTotal(void) const;
    SbBool hasDeferredCompactEntries(void) const;

    unsigned int getVisitedMeshCount(void) const;
    unsigned int getSubmittedTaskCount(void) const;
    unsigned int getSkippedMeshCount(void) const;
    const SbString &getDiagnostics(void) const;

protected:
    virtual void beginTraversal(SoNode *node);

private:
    static void nodeAction(SoAction *action, SoNode *node);
    static void databaseSourceAction(SoAction *action, SoNode *node);
    static void meshShapeAction(SoAction *action, SoNode *node);
    void appendDiagnostic(const SbString &target, const char *message);

    BRLObolLodService *service;
    struct db_i *dbip;
    SbString databaseId;
    uint64_t databaseRevision;
    struct bv_view_info view;
    uint64_t generation;
    uint64_t viewRevision;
    uint64_t policyRevision;
    SbString providerId;
    SbString providerVersion;
    int qualityTier;
    SbBool refreshMissing;
    int reset;
    SbBool useForcedLevel;
    int forcedLevel;
    SbBool requireLodBacked;
    SbBool submitAabbProxyStage;
    SbBool submitObbProxyStage;
    const BRLObolViewLodState *viewState;
    size_t compactEntryFirst;
    size_t compactEntryLimit;
    size_t compactEntryNext;
    size_t compactEntryTotal;
    SbBool deferredCompactEntries;
    unsigned int visitedMeshCount;
    unsigned int submittedTaskCount;
    unsigned int skippedMeshCount;
    SbString diagnostics;
};

#endif /* BRLOBOL_MESH_LOD_SUBMIT_ACTION_H */
