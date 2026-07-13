/*                  R E A L I Z E _ A C T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/realize_action.h */

#ifndef BRLOBOL_REALIZE_ACTION_H
#define BRLOBOL_REALIZE_ACTION_H

#include "brlobol/defines.h"

#include <Inventor/SbString.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/actions/SoSubAction.h>

class SoBRLDatabaseSource;
struct BRLObolDatabaseSourceRealizationCache;

class BRLOBOL_EXPORT BRLObolRealizationRepository {
public:
    BRLObolRealizationRepository(void);
    ~BRLObolRealizationRepository(void);
    BRLObolRealizationRepository(
	const BRLObolRealizationRepository &) = delete;
    BRLObolRealizationRepository &operator=(
	const BRLObolRealizationRepository &) = delete;

    void clear(void);
    void invalidateObject(const char *name);
    void renameObject(const char *oldName, const char *newName);
    void invalidateViewVariants(void);
    void seedSource(SoBRLDatabaseSource *source);
    void releaseSource(SoBRLDatabaseSource *source);

private:
    friend class SoBRLRealizeAction;
    struct Residency;
    BRLObolDatabaseSourceRealizationCache *cache;
    Residency *residency;
};

class BRLOBOL_EXPORT SoBRLRealizeAction : public SoAction {
    typedef SoAction inherited;

    SO_ACTION_HEADER(SoBRLRealizeAction);

public:
    SoBRLRealizeAction(void);
    virtual ~SoBRLRealizeAction(void);
    static void initClass(void);

    unsigned int getVisitedSourceCount(void) const;
    unsigned int getRealizedSourceCount(void) const;
    unsigned int getFailedSourceCount(void) const;
    const SbString &getDiagnostics(void) const;
    void setRetainRealizationCache(SbBool retain);
    void clearRealizationCache(void);
    void invalidateRealizationObject(const char *name);
    /* The repository is borrowed and must outlive this action. */
    void setRealizationRepository(BRLObolRealizationRepository *repository);

protected:
    virtual void beginTraversal(SoNode *node);

private:
    static void nodeAction(SoAction *action, SoNode *node);
    static void databaseSourceAction(SoAction *action, SoNode *node);
    void appendDiagnostic(const SoBRLDatabaseSource *source);

    unsigned int visitedSourceCount;
    unsigned int realizedSourceCount;
    unsigned int failedSourceCount;
    SbString diagnostics;
    BRLObolDatabaseSourceRealizationCache *realizationCache;
    BRLObolRealizationRepository *realizationRepository;
    SbBool ownsRealizationRepository;
    SbBool seedingCache;
    SbBool retainRealizationCache;
};

#endif /* BRLOBOL_REALIZE_ACTION_H */
