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
};

#endif /* BRLOBOL_REALIZE_ACTION_H */
