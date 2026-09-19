/*        M E S H _ R E S I D E N C Y _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BMeshResidencyAction.h"
#include "BObol/BMeshShape.h"

#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>

#include <algorithm>
#include <limits>

SO_ACTION_SOURCE(SoBRLMeshResidencyAction);

SoBRLMeshResidencyAction::SoBRLMeshResidencyAction(void) :
    maxResidentMeshBytes(std::numeric_limits<size_t>::max()),
    evictDisplayPayloads(TRUE),
    entries(),
    visitedMeshCount(0),
    evictedFullDetailMeshCount(0),
    evictedDisplayMeshCount(0),
    skippedDisplayMeshCount(0),
    initialResidentMeshBytes(0),
    finalResidentMeshBytes(0),
    freedFullDetailBytes(0),
    freedDisplayBytes(0)
{
    SO_ACTION_CONSTRUCTOR(SoBRLMeshResidencyAction);
}

SoBRLMeshResidencyAction::~SoBRLMeshResidencyAction(void)
{
}

void
SoBRLMeshResidencyAction::initClass(void)
{
    SO_ACTION_INIT_CLASS(SoBRLMeshResidencyAction, SoAction);
    SO_ACTION_ADD_METHOD(SoNode, SoBRLMeshResidencyAction::nodeAction);
    SO_ACTION_ADD_METHOD(SoGroup, SoBRLMeshResidencyAction::nodeAction);
    SO_ACTION_ADD_METHOD(SoBRLMeshShape,
			 SoBRLMeshResidencyAction::meshShapeAction);
}

void
SoBRLMeshResidencyAction::setMaxResidentMeshBytes(size_t maxBytes)
{
    this->maxResidentMeshBytes = maxBytes;
}

size_t
SoBRLMeshResidencyAction::getMaxResidentMeshBytes(void) const
{
    return this->maxResidentMeshBytes;
}

void
SoBRLMeshResidencyAction::setEvictDisplayPayloads(SbBool enabled)
{
    this->evictDisplayPayloads = enabled ? TRUE : FALSE;
}

SbBool
SoBRLMeshResidencyAction::isEvictDisplayPayloadsEnabled(void) const
{
    return this->evictDisplayPayloads;
}

unsigned int
SoBRLMeshResidencyAction::getVisitedMeshCount(void) const
{
    return this->visitedMeshCount;
}

unsigned int
SoBRLMeshResidencyAction::getEvictedFullDetailMeshCount(void) const
{
    return this->evictedFullDetailMeshCount;
}

unsigned int
SoBRLMeshResidencyAction::getEvictedDisplayMeshCount(void) const
{
    return this->evictedDisplayMeshCount;
}

unsigned int
SoBRLMeshResidencyAction::getSkippedDisplayMeshCount(void) const
{
    return this->skippedDisplayMeshCount;
}

size_t
SoBRLMeshResidencyAction::getInitialResidentMeshBytes(void) const
{
    return this->initialResidentMeshBytes;
}

size_t
SoBRLMeshResidencyAction::getFinalResidentMeshBytes(void) const
{
    return this->finalResidentMeshBytes;
}

size_t
SoBRLMeshResidencyAction::getFreedFullDetailBytes(void) const
{
    return this->freedFullDetailBytes;
}

size_t
SoBRLMeshResidencyAction::getFreedDisplayBytes(void) const
{
    return this->freedDisplayBytes;
}

size_t
SoBRLMeshResidencyAction::getFreedResidentMeshBytes(void) const
{
    return this->initialResidentMeshBytes > this->finalResidentMeshBytes ?
	   this->initialResidentMeshBytes - this->finalResidentMeshBytes : 0;
}

void
SoBRLMeshResidencyAction::beginTraversal(SoNode *node)
{
    this->resetResults();
    this->traverse(node);
    this->finalResidentMeshBytes = this->initialResidentMeshBytes;
    this->evictToBudget();
    this->recomputeFinalResidentBytes();
}

void
SoBRLMeshResidencyAction::nodeAction(SoAction *action, SoNode *node)
{
    if (node->isOfType(SoGroup::getClassTypeId()))
	node->doAction(action);
}

void
SoBRLMeshResidencyAction::meshShapeAction(SoAction *action, SoNode *node)
{
    SoBRLMeshResidencyAction *budgetAction =
	static_cast<SoBRLMeshResidencyAction *>(action);
    SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);

    budgetAction->visitedMeshCount++;

    Entry entry;
    entry.shape = shape;
    entry.residentBytes = shape->estimateResidentMeshBytes();
    entry.fullDetailBytes = shape->estimateFullDetailMeshBytes();
    entry.displayBytes = shape->estimateDisplayMeshBytes();

    budgetAction->initialResidentMeshBytes += entry.residentBytes;
    if (entry.residentBytes > 0)
	budgetAction->entries.push_back(entry);
}

void
SoBRLMeshResidencyAction::resetResults(void)
{
    this->entries.clear();
    this->visitedMeshCount = 0;
    this->evictedFullDetailMeshCount = 0;
    this->evictedDisplayMeshCount = 0;
    this->skippedDisplayMeshCount = 0;
    this->initialResidentMeshBytes = 0;
    this->finalResidentMeshBytes = 0;
    this->freedFullDetailBytes = 0;
    this->freedDisplayBytes = 0;
}

void
SoBRLMeshResidencyAction::evictToBudget(void)
{
    if (this->finalResidentMeshBytes <= this->maxResidentMeshBytes)
	return;

    std::sort(this->entries.begin(), this->entries.end(),
    [](const Entry &a, const Entry &b) {
	return a.fullDetailBytes > b.fullDetailBytes;
    });
    for (size_t i = 0; i < this->entries.size(); i++) {
	if (this->finalResidentMeshBytes <= this->maxResidentMeshBytes)
	    return;
	SoBRLMeshShape *shape = this->entries[i].shape;
	if (!shape)
	    continue;
	size_t freedBytes = shape->evictFullDetailMesh();
	if (freedBytes == 0)
	    continue;
	this->freedFullDetailBytes += freedBytes;
	this->evictedFullDetailMeshCount++;
	this->finalResidentMeshBytes =
	    this->finalResidentMeshBytes > freedBytes ?
	    this->finalResidentMeshBytes - freedBytes : 0;
    }

    if (!this->evictDisplayPayloads ||
	this->finalResidentMeshBytes <= this->maxResidentMeshBytes)
	return;

    for (size_t i = 0; i < this->entries.size(); i++) {
	SoBRLMeshShape *shape = this->entries[i].shape;
	this->entries[i].displayBytes = shape ?
					shape->estimateDisplayMeshBytes() : 0;
    }
    std::sort(this->entries.begin(), this->entries.end(),
    [](const Entry &a, const Entry &b) {
	return a.displayBytes > b.displayBytes;
    });

    for (size_t i = 0; i < this->entries.size(); i++) {
	if (this->finalResidentMeshBytes <= this->maxResidentMeshBytes)
	    return;
	SoBRLMeshShape *shape = this->entries[i].shape;
	if (!shape || !shape->isLodDisplayActive())
	    continue;
	if (!shape->needsSourceBackedFullDetail()) {
	    this->skippedDisplayMeshCount++;
	    continue;
	}
	size_t freedBytes = shape->evictActiveDisplayMesh();
	if (freedBytes == 0)
	    continue;
	this->freedDisplayBytes += freedBytes;
	this->evictedDisplayMeshCount++;
	this->finalResidentMeshBytes =
	    this->finalResidentMeshBytes > freedBytes ?
	    this->finalResidentMeshBytes - freedBytes : 0;
    }
}

void
SoBRLMeshResidencyAction::recomputeFinalResidentBytes(void)
{
    this->finalResidentMeshBytes = 0;
    for (size_t i = 0; i < this->entries.size(); i++) {
	if (this->entries[i].shape)
	    this->finalResidentMeshBytes +=
		this->entries[i].shape->estimateResidentMeshBytes();
    }
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
