/*                 V I E W _ A T T A C H M E N T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file libbrlobol/view_attachment.cpp */

#include "common.h"

#include "brlobol/view_attachment.h"
#include "brlobol/view_lod.h"

#include <Inventor/nodes/SoNode.h>

BRLObolViewAttachment::BRLObolViewAttachment(void) :
    refCount(0),
    sceneRoot(NULL),
    sceneRootToken(NULL),
    viewLodState(new BRLObolViewLodState),
    lodPolicy(BV_LOD_POLICY_INIT),
    lodBoundsCallbackSet(FALSE),
    independentScopeCreated(FALSE)
{
    bv_lod_policy_sanitize(&this->lodPolicy);
}

BRLObolViewAttachment::~BRLObolViewAttachment(void)
{
    this->setSceneRoot(NULL);
    delete this->viewLodState;
    this->viewLodState = NULL;
}

void
BRLObolViewAttachment::ref(void)
{
    this->refCount++;
}

void
BRLObolViewAttachment::unref(void)
{
    this->refCount--;
    if (this->refCount <= 0)
	delete this;
}

int
BRLObolViewAttachment::getRefCount(void) const
{
    return this->refCount;
}

void
BRLObolViewAttachment::copyHostStateFrom(
	const BRLObolViewAttachment *source)
{
    if (!source || source == this)
	return;

    this->sceneRootToken = source->sceneRootToken;
    this->lodPolicy = source->lodPolicy;
    bv_lod_policy_sanitize(&this->lodPolicy);
    this->lodBoundsCallbackSet = source->lodBoundsCallbackSet;
    this->independentScopeCreated = source->independentScopeCreated;
}

void
BRLObolViewAttachment::setSceneRoot(SoNode *root)
{
    if (this->sceneRoot == root)
	return;

    if (root)
	root->ref();
    if (this->sceneRoot)
	this->sceneRoot->unref();
    this->sceneRoot = root;

    this->clearViewLodState();
}

SoNode *
BRLObolViewAttachment::getSceneRoot(void) const
{
    return this->sceneRoot;
}

SbBool
BRLObolViewAttachment::hasSceneRoot(void) const
{
    return this->sceneRoot ? TRUE : FALSE;
}

void
BRLObolViewAttachment::setSceneRootToken(void *token)
{
    this->sceneRootToken = token;
}

void *
BRLObolViewAttachment::getSceneRootToken(void) const
{
    return this->sceneRootToken;
}

SbBool
BRLObolViewAttachment::hasSceneRootToken(void) const
{
    return this->sceneRootToken ? TRUE : FALSE;
}

void
BRLObolViewAttachment::setIndependentScopeCreated(SbBool created)
{
    this->independentScopeCreated = created ? TRUE : FALSE;
}

SbBool
BRLObolViewAttachment::isIndependentScopeCreated(void) const
{
    return this->independentScopeCreated;
}

void
BRLObolViewAttachment::setLodPolicy(const struct bv_lod_policy *policy)
{
    if (!policy)
	return;

    this->lodPolicy = *policy;
    bv_lod_policy_sanitize(&this->lodPolicy);
}

void
BRLObolViewAttachment::getLodPolicy(struct bv_lod_policy *policy) const
{
    if (!policy)
	return;

    *policy = this->lodPolicy;
    bv_lod_policy_sanitize(policy);
}

void
BRLObolViewAttachment::setLodBoundsCallbackSet(SbBool set)
{
    this->lodBoundsCallbackSet = set ? TRUE : FALSE;
}

SbBool
BRLObolViewAttachment::isLodBoundsCallbackSet(void) const
{
    return this->lodBoundsCallbackSet;
}

BRLObolViewLodState *
BRLObolViewAttachment::getViewLodState(void) const
{
    return this->viewLodState;
}

void
BRLObolViewAttachment::clearViewLodState(void)
{
    if (this->viewLodState)
	this->viewLodState->clear();
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
