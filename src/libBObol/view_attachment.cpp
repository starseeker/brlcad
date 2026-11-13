/*                 V I E W _ A T T A C H M E N T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file libBObol/view_attachment.cpp */

#include "common.h"

#include "BObol/BViewAttachment.h"
#include "BObol/BViewLod.h"

#include <Inventor/nodes/SoNode.h>

BObolViewAttachment::BObolViewAttachment(void) :
    refCount(0),
    sceneRoot(NULL),
    sceneRootToken(NULL),
    viewLodState(new BObolViewLodState),
    lodPolicy(BV_LOD_POLICY_INIT),
    lodBoundsCallbackSet(FALSE),
    independentScopeCreated(FALSE)
{
    bv_lod_policy_sanitize(&this->lodPolicy);
}

BObolViewAttachment::~BObolViewAttachment(void)
{
    this->setSceneRoot(NULL);
    delete this->viewLodState;
    this->viewLodState = NULL;
}

void
BObolViewAttachment::ref(void)
{
    this->refCount++;
}

void
BObolViewAttachment::unref(void)
{
    this->refCount--;
    if (this->refCount <= 0)
	delete this;
}

int
BObolViewAttachment::getRefCount(void) const
{
    return this->refCount;
}

void
BObolViewAttachment::copyHostStateFrom(
	const BObolViewAttachment *source)
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
BObolViewAttachment::setSceneRoot(SoNode *root)
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
BObolViewAttachment::getSceneRoot(void) const
{
    return this->sceneRoot;
}

SbBool
BObolViewAttachment::hasSceneRoot(void) const
{
    return this->sceneRoot ? TRUE : FALSE;
}

void
BObolViewAttachment::setSceneRootToken(void *token)
{
    this->sceneRootToken = token;
}

void *
BObolViewAttachment::getSceneRootToken(void) const
{
    return this->sceneRootToken;
}

SbBool
BObolViewAttachment::hasSceneRootToken(void) const
{
    return this->sceneRootToken ? TRUE : FALSE;
}

void
BObolViewAttachment::setIndependentScopeCreated(SbBool created)
{
    this->independentScopeCreated = created ? TRUE : FALSE;
}

SbBool
BObolViewAttachment::isIndependentScopeCreated(void) const
{
    return this->independentScopeCreated;
}

void
BObolViewAttachment::setLodPolicy(const struct bv_lod_policy *policy)
{
    if (!policy)
	return;

    this->lodPolicy = *policy;
    bv_lod_policy_sanitize(&this->lodPolicy);
}

void
BObolViewAttachment::getLodPolicy(struct bv_lod_policy *policy) const
{
    if (!policy)
	return;

    *policy = this->lodPolicy;
    bv_lod_policy_sanitize(policy);
}

void
BObolViewAttachment::setLodBoundsCallbackSet(SbBool set)
{
    this->lodBoundsCallbackSet = set ? TRUE : FALSE;
}

SbBool
BObolViewAttachment::isLodBoundsCallbackSet(void) const
{
    return this->lodBoundsCallbackSet;
}

BObolViewLodState *
BObolViewAttachment::getViewLodState(void) const
{
    return this->viewLodState;
}

void
BObolViewAttachment::clearViewLodState(void)
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
