/*            C O M P A C T _ P R E S E N T A T I O N _ S T A G I N G . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "compact_presentation_staging_private.h"

#include <algorithm>
#include <utility>

BObolCompactPresentationStaging::BObolCompactPresentationStaging(
    SoBRLCadAssembly *assembly, bool reset) :
    assembly_(assembly), reset_(reset)
{
}

uint8_t
BObolCompactPresentationStaging::partChannel(Obol::PartId part) const
{
    const auto staged = this->partChannels_.find(part);
    if (staged != this->partChannels_.end())
	return staged->second;
    if (this->reset_)
	return 0;
    const auto retained = this->assembly_->compactPartChannels.find(part);
    return retained == this->assembly_->compactPartChannels.end() ? 0 :
	retained->second;
}

void
BObolCompactPresentationStaging::setPartChannel(
    Obol::PartId part, uint8_t channels)
{
    this->removedPartChannels_.erase(part);
    this->partChannels_[part] = channels;
}

void
BObolCompactPresentationStaging::addLodPart(Obol::PartId part)
{
    this->removedLodParts_.erase(part);
    this->lodParts_.insert(part);
}

bool
BObolCompactPresentationStaging::isSpatiallyCulled(
    Obol::InstanceId instance) const
{
    const auto staged = this->spatialCull_.find(instance);
    if (staged != this->spatialCull_.end())
	return staged->second;
    return !this->reset_ &&
	this->assembly_->compactSpatiallyCulledInstances.count(instance) != 0;
}

void
BObolCompactPresentationStaging::setSpatiallyCulled(
    Obol::InstanceId instance, bool culled)
{
    this->spatialCull_[instance] = culled;
}

const BObolCompactPresentationStaging::CompactPresentation *
BObolCompactPresentationStaging::findPresentation(
    Obol::InstanceId instance) const
{
    const auto staged = this->presentations_.find(instance);
    if (staged != this->presentations_.end())
	return &staged->second;
    if (this->reset_)
	return NULL;
    const auto retained =
	this->assembly_->compactInstancePresentations.find(instance);
    return retained == this->assembly_->compactInstancePresentations.end() ?
	NULL : &retained->second;
}

BObolCompactPresentationStaging::CompactPresentation &
BObolCompactPresentationStaging::editPresentation(
    Obol::InstanceId instance)
{
    auto staged = this->presentations_.find(instance);
    if (staged != this->presentations_.end())
	return staged->second;
    CompactPresentation initial;
    if (!this->reset_) {
	const auto retained =
	    this->assembly_->compactInstancePresentations.find(instance);
	if (retained != this->assembly_->compactInstancePresentations.end())
	    initial = retained->second;
    }
    return this->presentations_.emplace(
	instance, std::move(initial)).first->second;
}

const BObolCompactPresentationStaging::CompactLayers *
BObolCompactPresentationStaging::findLayers(Obol::InstanceId instance) const
{
    const auto staged = this->layerPresentations_.find(instance);
    if (staged != this->layerPresentations_.end())
	return &staged->second;
    if (this->removedLayerPresentations_.count(instance) || this->reset_)
	return NULL;
    const auto retained =
	this->assembly_->compactLayerPresentations.find(instance);
    return retained == this->assembly_->compactLayerPresentations.end() ?
	NULL : &retained->second;
}

void
BObolCompactPresentationStaging::setLayers(
    Obol::InstanceId instance, CompactLayers layers)
{
    if (layers.empty()) {
	this->layerPresentations_.erase(instance);
	this->removedLayerPresentations_.insert(instance);
	return;
    }
    this->removedLayerPresentations_.erase(instance);
    this->layerPresentations_[instance] = std::move(layers);
}

size_t
BObolCompactPresentationStaging::activePartReferences(Obol::PartId part) const
{
    ptrdiff_t count = 0;
    if (!this->reset_) {
	const auto retained =
	    this->assembly_->compactActivePartReferences.find(part);
	if (retained != this->assembly_->compactActivePartReferences.end())
	    count = static_cast<ptrdiff_t>(retained->second);
    }
    const auto staged = this->activePartReferenceDeltas_.find(part);
    if (staged != this->activePartReferenceDeltas_.end())
	count += staged->second;
    return count > 0 ? static_cast<size_t>(count) : 0;
}

void
BObolCompactPresentationStaging::addPartReference(Obol::PartId part)
{
    ++this->activePartReferenceDeltas_[part];
}

void
BObolCompactPresentationStaging::removePartReference(Obol::PartId part)
{
    if (this->activePartReferences(part) != 0)
	--this->activePartReferenceDeltas_[part];
}

std::vector<Obol::InstanceId>
BObolCompactPresentationStaging::expandInstanceSet(
    const std::vector<Obol::InstanceId> &baseSet) const
{
    std::vector<Obol::InstanceId> expanded = baseSet;
    const auto appendLayers = [&baseSet, &expanded](
	    Obol::InstanceId instance, const CompactLayers &layers) {
	if (!std::binary_search(baseSet.begin(), baseSet.end(), instance))
	    return;
	for (size_t layerIndex = 1; layerIndex < layers.size(); ++layerIndex)
	    expanded.push_back(layers[layerIndex].instance);
    };
    if (!this->reset_) {
	for (const auto &item : this->assembly_->compactLayerPresentations) {
	    if (this->layerPresentations_.count(item.first) ||
		this->removedLayerPresentations_.count(item.first))
		continue;
	    appendLayers(item.first, item.second);
	}
    }
    for (const auto &item : this->layerPresentations_)
	appendLayers(item.first, item.second);
    std::sort(expanded.begin(), expanded.end());
    expanded.erase(std::unique(expanded.begin(), expanded.end()),
	expanded.end());
    return expanded;
}

void
BObolCompactPresentationStaging::appendSpatiallyCulled(
    std::vector<Obol::InstanceId> &instances) const
{
    if (!this->reset_) {
	for (const Obol::InstanceId instance :
	     this->assembly_->compactSpatiallyCulledInstances) {
	    if (this->spatialCull_.find(instance) == this->spatialCull_.end())
		instances.push_back(instance);
	}
    }
    for (const auto &item : this->spatialCull_)
	if (item.second)
	    instances.push_back(item.first);
}

void
BObolCompactPresentationStaging::removePart(Obol::PartId part)
{
    this->partChannels_.erase(part);
    this->removedPartChannels_.insert(part);
    this->lodParts_.erase(part);
    this->removedLodParts_.insert(part);
}

void
BObolCompactPresentationStaging::commit()
{
    if (this->reset_) {
	this->assembly_->compactPartChannels.clear();
	this->assembly_->compactInstancePresentations.clear();
	this->assembly_->compactLayerPresentations.clear();
	this->assembly_->compactActivePartReferences.clear();
	this->assembly_->compactLodParts.clear();
	this->assembly_->compactSpatiallyCulledInstances.clear();
    }
    for (const Obol::PartId part : this->removedPartChannels_)
	this->assembly_->compactPartChannels.erase(part);
    for (const auto &item : this->partChannels_)
	this->assembly_->compactPartChannels[item.first] = item.second;
    for (const Obol::PartId part : this->removedLodParts_)
	this->assembly_->compactLodParts.erase(part);
    this->assembly_->compactLodParts.insert(
	this->lodParts_.begin(), this->lodParts_.end());
    for (const auto &item : this->spatialCull_) {
	if (item.second)
	    this->assembly_->compactSpatiallyCulledInstances.insert(item.first);
	else
	    this->assembly_->compactSpatiallyCulledInstances.erase(item.first);
    }
    for (const auto &item : this->presentations_)
	this->assembly_->compactInstancePresentations[item.first] = item.second;
    for (const Obol::InstanceId instance :
	 this->removedLayerPresentations_)
	this->assembly_->compactLayerPresentations.erase(instance);
    for (auto &item : this->layerPresentations_)
	this->assembly_->compactLayerPresentations[item.first] =
	    std::move(item.second);
    for (const auto &item : this->activePartReferenceDeltas_) {
	const size_t references = this->activePartReferences(item.first);
	if (references)
	    this->assembly_->compactActivePartReferences[item.first] = references;
	else
	    this->assembly_->compactActivePartReferences.erase(item.first);
    }
}
