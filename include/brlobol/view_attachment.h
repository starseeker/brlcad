/*                  V I E W _ A T T A C H M E N T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this file; see the file named COPYING for more information.
 */
/** @file brlobol/view_attachment.h */

#ifndef BRLOBOL_VIEW_ATTACHMENT_H
#define BRLOBOL_VIEW_ATTACHMENT_H

#include "brlobol/defines.h"

#include "bv/view.h"

#include <Inventor/SbBasic.h>

class BRLObolViewLodState;
class SoNode;

/**
 * Per-view Obol attachment state for a BRL-CAD view.
 *
 * This object is deliberately independent of GED.  Application-facing GED view
 * records may refer to it, but scene roots, view-local LoD realization state,
 * and passive LoD policy belong to the Obol/libbrlobol view attachment.
 */
class BRLOBOL_EXPORT BRLObolViewAttachment
{
public:
    BRLObolViewAttachment(void);
    ~BRLObolViewAttachment(void);

    void ref(void);
    void unref(void);
    int getRefCount(void) const;

    void copyHostStateFrom(const BRLObolViewAttachment *source);

    void setSceneRoot(SoNode *root);
    SoNode *getSceneRoot(void) const;
    SbBool hasSceneRoot(void) const;
    void setSceneRootToken(void *token);
    void *getSceneRootToken(void) const;
    SbBool hasSceneRootToken(void) const;

    void setIndependentScopeCreated(SbBool created);
    SbBool isIndependentScopeCreated(void) const;

    void setLodPolicy(const struct bv_lod_policy *policy);
    void getLodPolicy(struct bv_lod_policy *policy) const;
    void setLodBoundsCallbackSet(SbBool set);
    SbBool isLodBoundsCallbackSet(void) const;

    BRLObolViewLodState *getViewLodState(void) const;
    void clearViewLodState(void);

private:
    int refCount;
    SoNode *sceneRoot;
    void *sceneRootToken;
    BRLObolViewLodState *viewLodState;
    struct bv_lod_policy lodPolicy;
    SbBool lodBoundsCallbackSet;
    SbBool independentScopeCreated;
};

#endif /* BRLOBOL_VIEW_ATTACHMENT_H */
