/*              B R E A L I Z E A C T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BRealizeAction.h */

#ifndef BOBOL_BREALIZEACTION_H
#define BOBOL_BREALIZEACTION_H

#include "BObol/BDefines.h"

#include <Inventor/SbString.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/actions/SoSubAction.h>

class SoBRLDatabaseSource;
struct BObolDatabaseSourceRealizationCache;

class BOBOL_EXPORT BObolRealizationRepository {
public:
    BObolRealizationRepository(void);
    ~BObolRealizationRepository(void);
    BObolRealizationRepository(
	const BObolRealizationRepository &) = delete;
    BObolRealizationRepository &operator=(
	const BObolRealizationRepository &) = delete;

    void clear(void);
    void invalidateObject(const char *name);
    void renameObject(const char *oldName, const char *newName);
    void invalidateViewVariants(void);
    void seedSource(SoBRLDatabaseSource *source);
    void releaseSource(SoBRLDatabaseSource *source);

private:
    friend class SoBRLRealizeAction;
    struct Residency;
    BObolDatabaseSourceRealizationCache *cache;
    Residency *residency;
};

class BOBOL_EXPORT SoBRLRealizeAction : public SoAction {
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
    void setRealizationRepository(BObolRealizationRepository *repository);

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
    BObolDatabaseSourceRealizationCache *realizationCache;
    BObolRealizationRepository *realizationRepository;
    SbBool ownsRealizationRepository;
    SbBool seedingCache;
    SbBool retainRealizationCache;
};

#endif /* BOBOL_BREALIZEACTION_H */
