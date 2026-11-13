/*              E D I T _ M A N I P U L A T O R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BEditManipulator.h"
#include "scalar_publication_private.h"

#include <Inventor/SbViewVolume.h>
#include <Inventor/misc/SoChildList.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoLightModel.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoPointSet.h>
#include <Inventor/tools/SbModernUtils.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

SO_NODE_SOURCE(SoBRLEditManipulator);
SO_NODE_SOURCE(SoBRLIndexedEditManipulator);

static constexpr size_t AXIS_MANIPULATOR_CHILD_COUNT = 4;
static constexpr size_t MAX_INDEXED_MANIPULATOR_CHILD_COUNT = 6;

struct ManipulatorGeometry {
    std::vector<SbModernUtils::SoNodeRef> owners;
    std::vector<SoNode *> children;
};

template <typename Notifications, typename Commit>
static void
publish_manipulator_geometry(SoSeparator &manipulator,
	const ManipulatorGeometry &geometry, Notifications &notifications,
	Commit commit)
{
    std::unique_ptr<SoChildList::Replacement> replacement;
    if (!geometry.children.empty() || manipulator.getNumChildren())
	replacement = manipulator.getChildren()->prepareReplacement(
	    geometry.children);

    if (replacement)
	replacement->commit();
    commit();
    notifications.restore();

    std::exception_ptr failure;
    if (replacement) {
	try { replacement->notify(); }
	catch (...) { failure = std::current_exception(); }
    }
    if (notifications.changed())
	notifications.notify(failure);
    if (failure)
	std::rethrow_exception(failure);
}

static ManipulatorGeometry
build_axis_manipulator_geometry(const SbVec3f &center,
	const SbVec3f (&axes)[3], SbBool visible, int hoverHandle,
	int activeHandle)
{
    ManipulatorGeometry geometry;
    if (!visible)
	return geometry;

    geometry.owners.reserve(AXIS_MANIPULATOR_CHILD_COUNT);
    geometry.children.reserve(AXIS_MANIPULATOR_CHILD_COUNT);

    SbModernUtils::SoNodeRef lightModelOwner(new SoLightModel);
    auto *lightModel = static_cast<SoLightModel *>(lightModelOwner.get());
    lightModel->model = SoLightModel::BASE_COLOR;
    geometry.children.push_back(lightModel);
    geometry.owners.push_back(std::move(lightModelOwner));

    const SbColor colors[3] = {
	SbColor(1.0f, 0.2f, 0.2f),
	SbColor(0.2f, 1.0f, 0.2f),
	SbColor(0.25f, 0.55f, 1.0f)
    };
    for (int i = 0; i < 3; i++) {
	SbModernUtils::SoNodeRef axisRootOwner(new SoSeparator);
	auto *axisRoot = static_cast<SoSeparator *>(axisRootOwner.get());
	SbModernUtils::SoNodeRef materialOwner(new SoMaterial);
	auto *material = static_cast<SoMaterial *>(materialOwner.get());
	SbColor color = colors[i];
	if (activeHandle == i)
	    color = SbColor(1.0f, 1.0f, 1.0f);
	else if (hoverHandle == i)
	    color = color + SbColor(0.25f, 0.25f, 0.25f);
	color[0] = std::min(color[0], 1.0f);
	color[1] = std::min(color[1], 1.0f);
	color[2] = std::min(color[2], 1.0f);
	material->diffuseColor = color;
	axisRoot->addChild(material);

	SbModernUtils::SoNodeRef styleOwner(new SoDrawStyle);
	auto *style = static_cast<SoDrawStyle *>(styleOwner.get());
	style->lineWidth = activeHandle == i ? 4.0f : 2.0f;
	style->pointSize = activeHandle == i ? 13.0f :
	    (hoverHandle == i ? 12.0f : 10.0f);
	axisRoot->addChild(style);

	SbModernUtils::SoNodeRef coordinatesOwner(new SoCoordinate3);
	auto *coordinates = static_cast<SoCoordinate3 *>(coordinatesOwner.get());
	const SbVec3f values[2] = {center, center + axes[i]};
	coordinates->point.setValues(0, 2, values);
	axisRoot->addChild(coordinates);

	SbModernUtils::SoNodeRef lineOwner(new SoLineSet);
	auto *line = static_cast<SoLineSet *>(lineOwner.get());
	line->numVertices.set1Value(0, 2);
	axisRoot->addChild(line);

	SbModernUtils::SoNodeRef pointOwner(new SoPointSet);
	auto *point = static_cast<SoPointSet *>(pointOwner.get());
	point->startIndex = 1;
	point->numPoints = 1;
	axisRoot->addChild(point);
	geometry.children.push_back(axisRoot);
	geometry.owners.push_back(std::move(axisRootOwner));
    }
    return geometry;
}

SoBRLEditManipulator::SoBRLEditManipulator(void) :
    editCenter(0.0f, 0.0f, 0.0f)
{
    SO_NODE_CONSTRUCTOR(SoBRLEditManipulator);
    SO_NODE_ADD_FIELD(manipulatorId, (""));
    SO_NODE_ADD_FIELD(sessionRevision, (0));
    SO_NODE_ADD_FIELD(visible, (TRUE));
    SO_NODE_ADD_FIELD(hoverHandle, (HANDLE_NONE));
    SO_NODE_ADD_FIELD(activeHandle, (HANDLE_NONE));

    editAxes[0] = SbVec3f(1.0f, 0.0f, 0.0f);
    editAxes[1] = SbVec3f(0.0f, 1.0f, 0.0f);
    editAxes[2] = SbVec3f(0.0f, 0.0f, 1.0f);
    this->rebuildGeometry();
}

SoBRLEditManipulator::~SoBRLEditManipulator(void)
{
}

void
SoBRLEditManipulator::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLEditManipulator, SoSeparator, "Separator");
}

void
SoBRLEditManipulator::setEllipsoidAxes(const SbVec3f &nextCenter,
	const SbVec3f &axisA, const SbVec3f &axisB, const SbVec3f &axisC)
{
    const SbVec3f nextAxes[3] = {axisA, axisB, axisC};
    const ManipulatorGeometry geometry = build_axis_manipulator_geometry(
	nextCenter, nextAxes, this->visible.getValue(),
	this->hoverHandle.getValue(), this->activeHandle.getValue());
    PreparedFieldNotifications<0> notifications(*this, {});
    publish_manipulator_geometry(*this, geometry, notifications, [&] {
	editCenter = nextCenter;
	for (int i = 0; i < 3; i++)
	    editAxes[i] = nextAxes[i];
    });
}

void
SoBRLEditManipulator::setVisible(SbBool value)
{
    if (this->visible.getValue() == value)
	return;
    const ManipulatorGeometry geometry = build_axis_manipulator_geometry(
	editCenter, editAxes, value, this->hoverHandle.getValue(),
	this->activeHandle.getValue());
    PreparedFieldNotifications<1> notifications(*this, {{{&this->visible,
	true}}});
    publish_manipulator_geometry(*this, geometry, notifications,
	[&] { this->visible = value; });
}

void
SoBRLEditManipulator::setHoverHandle(Handle handle)
{
    if (this->hoverHandle.getValue() == static_cast<int>(handle))
	return;
    const int nextHandle = static_cast<int>(handle);
    const ManipulatorGeometry geometry = build_axis_manipulator_geometry(
	editCenter, editAxes, this->visible.getValue(), nextHandle,
	this->activeHandle.getValue());
    PreparedFieldNotifications<1> notifications(*this, {{{&this->hoverHandle,
	true}}});
    publish_manipulator_geometry(*this, geometry, notifications,
	[&] { this->hoverHandle = nextHandle; });
}

void
SoBRLEditManipulator::setActiveHandle(Handle handle)
{
    if (this->activeHandle.getValue() == static_cast<int>(handle))
	return;
    const int nextHandle = static_cast<int>(handle);
    const ManipulatorGeometry geometry = build_axis_manipulator_geometry(
	editCenter, editAxes, this->visible.getValue(),
	this->hoverHandle.getValue(), nextHandle);
    PreparedFieldNotifications<1> notifications(*this, {{{&this->activeHandle,
	true}}});
    publish_manipulator_geometry(*this, geometry, notifications,
	[&] { this->activeHandle = nextHandle; });
}

SbVec3f
SoBRLEditManipulator::center(void) const
{
    return editCenter;
}

SbVec3f
SoBRLEditManipulator::axis(Handle handle) const
{
    const int index = static_cast<int>(handle);
    return index >= 0 && index < 3 ? editAxes[index] : SbVec3f(0, 0, 0);
}

void
SoBRLEditManipulator::rebuildGeometry(void)
{
    if (!this->visible.getValue() && !this->getNumChildren())
	return;
    const ManipulatorGeometry geometry = build_axis_manipulator_geometry(
	editCenter, editAxes, this->visible.getValue(),
	this->hoverHandle.getValue(), this->activeHandle.getValue());
    PreparedFieldNotifications<0> notifications(*this, {});
    publish_manipulator_geometry(*this, geometry, notifications, [] {});
}

SbBool
SoBRLEditManipulator::project(const SbVec3f &point, int width, int height,
	const SoCamera *camera, SbVec3f &pixel) const
{
    if (!camera || width <= 0 || height <= 0)
	return FALSE;
    const float aspect = static_cast<float>(width) /
	static_cast<float>(height);
    SbVec3f normalized;
    camera->getViewVolume(aspect).projectToScreen(point, normalized);
    if (!std::isfinite(normalized[0]) || !std::isfinite(normalized[1]) ||
	!std::isfinite(normalized[2]))
	return FALSE;
    pixel.setValue(normalized[0] * static_cast<float>(width),
	(1.0f - normalized[1]) * static_cast<float>(height), normalized[2]);
    return TRUE;
}

SoBRLEditManipulator::Handle
SoBRLEditManipulator::hitTest(int x, int y, int width, int height,
	const SoCamera *camera, float radiusPixels) const
{
    if (!this->visible.getValue() || !camera || radiusPixels <= 0.0f)
	return HANDLE_NONE;
    const float radiusSquared = radiusPixels * radiusPixels;
    float bestDistance = std::numeric_limits<float>::max();
    float bestDepth = std::numeric_limits<float>::max();
    Handle best = HANDLE_NONE;
    for (int i = 0; i < 3; i++) {
	SbVec3f endpoint;
	if (!this->project(editCenter + editAxes[i], width, height, camera,
		endpoint))
	    continue;
	const float dx = endpoint[0] - static_cast<float>(x);
	const float dy = endpoint[1] - static_cast<float>(y);
	const float distance = dx * dx + dy * dy;
	const float distanceDelta = std::fabs(distance - bestDistance);
	if (distance <= radiusSquared &&
	    (distance < bestDistance ||
	    (distanceDelta <= 1.0e-6f && endpoint[2] < bestDepth))) {
	    bestDistance = distance;
	    bestDepth = endpoint[2];
	    best = static_cast<Handle>(i);
	}
    }
    return best;
}

SbBool
SoBRLEditManipulator::projectedScale(Handle handle, int x, int y,
	int width, int height, const SoCamera *camera, float &factor) const
{
    factor = 1.0f;
    const int index = static_cast<int>(handle);
    if (index < 0 || index >= 3)
	return FALSE;
    SbVec3f start;
    SbVec3f end;
    if (!this->project(editCenter, width, height, camera, start) ||
	!this->project(editCenter + editAxes[index], width, height, camera, end))
	return FALSE;
    const float dx = end[0] - start[0];
    const float dy = end[1] - start[1];
    const float lengthSquared = dx * dx + dy * dy;
    if (lengthSquared < 1.0e-6f)
	return FALSE;
    factor = ((static_cast<float>(x) - start[0]) * dx +
	(static_cast<float>(y) - start[1]) * dy) / lengthSquared;
    return std::isfinite(factor) ? TRUE : FALSE;
}

SbBool
SoBRLEditManipulator::screenPosition(Handle handle, float factor, int width,
	int height, const SoCamera *camera, int &x, int &y) const
{
    x = 0;
    y = 0;
    const int index = static_cast<int>(handle);
    if (index < 0 || index >= 3 || !std::isfinite(factor))
	return FALSE;
    SbVec3f pixel;
    if (!this->project(editCenter + editAxes[index] * factor, width, height,
	    camera, pixel))
	return FALSE;
    x = static_cast<int>(std::lround(pixel[0]));
    y = static_cast<int>(std::lround(pixel[1]));
    return TRUE;
}


class SoBRLIndexedEditManipulator::Private {
public:
    struct Face {
	std::vector<int32_t> vertices;
    };

    std::vector<SbVec3f> points;
    std::vector<int32_t> edges;
    std::vector<int32_t> edgeFeatures;
    int edgeFeatureCount = 0;
    int pointFeatureCount = 0;
    std::vector<Face> faces;

    int pointCount(void) const { return pointFeatureCount; }
    int edgeCount(void) const { return edgeFeatureCount; }
    int faceCount(void) const { return static_cast<int>(faces.size()); }
    SbBool featurePosition(SoBRLIndexedEditManipulator::Domain domain,
	int index, SbVec3f &position) const;
    ManipulatorGeometry buildGeometry(SbBool visible,
	SoBRLIndexedEditManipulator::Domain domain, int selectedIndex,
	int hoverIndex, int activeIndex) const;
};


SoBRLIndexedEditManipulator::SoBRLIndexedEditManipulator(void) :
    d(new Private)
{
    SO_NODE_CONSTRUCTOR(SoBRLIndexedEditManipulator);
    SO_NODE_ADD_FIELD(manipulatorId, (""));
    SO_NODE_ADD_FIELD(sessionRevision, (0));
    SO_NODE_ADD_FIELD(visible, (TRUE));
    SO_NODE_ADD_FIELD(selectionDomain, (DOMAIN_VERTEX));
    SO_NODE_ADD_FIELD(selectedIndex, (-1));
    SO_NODE_ADD_FIELD(hoverIndex, (-1));
    SO_NODE_ADD_FIELD(activeIndex, (-1));
}


SoBRLIndexedEditManipulator::~SoBRLIndexedEditManipulator(void)
{
    delete d;
    d = nullptr;
}


void
SoBRLIndexedEditManipulator::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLIndexedEditManipulator, SoSeparator, "Separator");
}


void
SoBRLIndexedEditManipulator::setTopology(const SbVec3f *points,
	int nextPointCount, const int32_t *edgeIndices, int nextEdgeCount,
	const int32_t *faceIndices, const int32_t *faceVertexCounts,
	int nextFaceCount, const int32_t *edgeFeatureIndices,
	int nextVertexFeatureCount)
{
    auto next = std::make_unique<Private>();
    if (points && nextPointCount > 0)
	next->points.assign(points, points + nextPointCount);
    const int storedPointCount = static_cast<int>(next->points.size());
    next->pointFeatureCount = nextVertexFeatureCount < 0 ? storedPointCount :
	std::max(0, std::min(storedPointCount, nextVertexFeatureCount));
    if (edgeIndices && nextEdgeCount > 0) {
	for (int i = 0; i < nextEdgeCount; i++) {
	    const int32_t a = edgeIndices[i * 2];
	    const int32_t b = edgeIndices[i * 2 + 1];
	    if (a < 0 || b < 0 || a >= storedPointCount ||
		b >= storedPointCount ||
		a == b)
		continue;
	    next->edges.push_back(a);
	    next->edges.push_back(b);
	    const int32_t feature = edgeFeatureIndices ?
		edgeFeatureIndices[i] : static_cast<int32_t>(i);
	    if (feature < 0 || feature == std::numeric_limits<int32_t>::max()) {
		next->edges.resize(next->edges.size() - 2);
		continue;
	    }
	    next->edgeFeatures.push_back(feature);
	    next->edgeFeatureCount = std::max(next->edgeFeatureCount,
		static_cast<int>(feature) + 1);
	}
    }
    if (faceIndices && faceVertexCounts && nextFaceCount > 0) {
	size_t offset = 0;
	for (int fi = 0; fi < nextFaceCount; fi++) {
	    const int count = faceVertexCounts[fi];
	    Private::Face face;
	    if (count >= 3) {
		for (int vi = 0; vi < count; vi++) {
		    const int32_t vertex = faceIndices[offset +
			static_cast<size_t>(vi)];
		    if (vertex < 0 || vertex >= storedPointCount) {
			face.vertices.clear();
			break;
		    }
		    face.vertices.push_back(vertex);
		}
	    }
	    if (count > 0)
		offset += static_cast<size_t>(count);
	    if (face.vertices.size() >= 3)
		next->faces.push_back(std::move(face));
	}
    }

    const Domain domain = static_cast<Domain>(
	this->selectionDomain.getValue());
    int nextSelected = this->selectedIndex.getValue();
    if ((domain == DOMAIN_VERTEX && nextSelected >= next->pointCount()) ||
	(domain == DOMAIN_EDGE && nextSelected >= next->edgeCount()) ||
	(domain == DOMAIN_FACE && nextSelected >= next->faceCount()))
	nextSelected = -1;

    const ManipulatorGeometry geometry = next->buildGeometry(
	this->visible.getValue(), domain, nextSelected,
	this->hoverIndex.getValue(), this->activeIndex.getValue());
    const bool selectionChanged =
	this->selectedIndex.getValue() != nextSelected;
    PreparedFieldNotifications<1> notifications(*this,
	{{{&this->selectedIndex, selectionChanged}}});
    publish_manipulator_geometry(*this, geometry, notifications, [&] {
	Private *previous = d;
	d = next.release();
	next.reset(previous);
	if (selectionChanged)
	    this->selectedIndex = nextSelected;
    });
}


void
SoBRLIndexedEditManipulator::setVisible(SbBool value)
{
    if (this->visible.getValue() == value)
	return;
    const ManipulatorGeometry geometry = d->buildGeometry(value,
	static_cast<Domain>(this->selectionDomain.getValue()),
	this->selectedIndex.getValue(), this->hoverIndex.getValue(),
	this->activeIndex.getValue());
    PreparedFieldNotifications<1> notifications(*this,
	{{{&this->visible, true}}});
    publish_manipulator_geometry(*this, geometry, notifications,
	[&] { this->visible = value; });
}


void
SoBRLIndexedEditManipulator::setSelectionDomain(Domain domain)
{
    if (this->selectionDomain.getValue() == static_cast<int>(domain))
	return;
    const ManipulatorGeometry geometry = d->buildGeometry(
	this->visible.getValue(), domain, -1, -1, -1);
    PreparedFieldNotifications<4> notifications(*this, {{
	{&this->selectionDomain, true},
	{&this->selectedIndex, this->selectedIndex.getValue() != -1},
	{&this->hoverIndex, this->hoverIndex.getValue() != -1},
	{&this->activeIndex, this->activeIndex.getValue() != -1}
    }});
    publish_manipulator_geometry(*this, geometry, notifications, [&] {
	this->selectionDomain = static_cast<int>(domain);
	this->selectedIndex = -1;
	this->hoverIndex = -1;
	this->activeIndex = -1;
    });
}


void
SoBRLIndexedEditManipulator::setSelectedIndex(int index)
{
    if (this->selectedIndex.getValue() == index)
	return;
    const ManipulatorGeometry geometry = d->buildGeometry(
	this->visible.getValue(),
	static_cast<Domain>(this->selectionDomain.getValue()), index,
	this->hoverIndex.getValue(), this->activeIndex.getValue());
    PreparedFieldNotifications<1> notifications(*this,
	{{{&this->selectedIndex, true}}});
    publish_manipulator_geometry(*this, geometry, notifications,
	[&] { this->selectedIndex = index; });
}


void
SoBRLIndexedEditManipulator::setHoverIndex(int index)
{
    if (this->hoverIndex.getValue() == index)
	return;
    const ManipulatorGeometry geometry = d->buildGeometry(
	this->visible.getValue(),
	static_cast<Domain>(this->selectionDomain.getValue()),
	this->selectedIndex.getValue(), index, this->activeIndex.getValue());
    PreparedFieldNotifications<1> notifications(*this,
	{{{&this->hoverIndex, true}}});
    publish_manipulator_geometry(*this, geometry, notifications,
	[&] { this->hoverIndex = index; });
}


void
SoBRLIndexedEditManipulator::setActiveIndex(int index)
{
    if (this->activeIndex.getValue() == index)
	return;
    const ManipulatorGeometry geometry = d->buildGeometry(
	this->visible.getValue(),
	static_cast<Domain>(this->selectionDomain.getValue()),
	this->selectedIndex.getValue(), this->hoverIndex.getValue(), index);
    PreparedFieldNotifications<1> notifications(*this,
	{{{&this->activeIndex, true}}});
    publish_manipulator_geometry(*this, geometry, notifications,
	[&] { this->activeIndex = index; });
}


int
SoBRLIndexedEditManipulator::pointCount(void) const
{
    return d->pointCount();
}


int
SoBRLIndexedEditManipulator::edgeCount(void) const
{
    return d->edgeCount();
}


int
SoBRLIndexedEditManipulator::faceCount(void) const
{
    return d->faceCount();
}


ManipulatorGeometry
SoBRLIndexedEditManipulator::Private::buildGeometry(SbBool visible,
	SoBRLIndexedEditManipulator::Domain domain, int selectedIndex,
	int hoverIndex, int activeIndex) const
{
    ManipulatorGeometry geometry;
    if (!visible || points.empty())
	return geometry;

    geometry.owners.reserve(MAX_INDEXED_MANIPULATOR_CHILD_COUNT);
    geometry.children.reserve(MAX_INDEXED_MANIPULATOR_CHILD_COUNT);

    SbModernUtils::SoNodeRef lightModelOwner(new SoLightModel);
    auto *lightModel = static_cast<SoLightModel *>(lightModelOwner.get());
    lightModel->model = SoLightModel::BASE_COLOR;
    geometry.children.push_back(lightModel);
    geometry.owners.push_back(std::move(lightModelOwner));

    if (!edges.empty()) {
	SbModernUtils::SoNodeRef edgeRootOwner(new SoSeparator);
	auto *edgeRoot = static_cast<SoSeparator *>(edgeRootOwner.get());
	SbModernUtils::SoNodeRef materialOwner(new SoMaterial);
	auto *material = static_cast<SoMaterial *>(materialOwner.get());
	material->diffuseColor = SbColor(0.35f, 0.65f, 1.0f);
	edgeRoot->addChild(material);
	SbModernUtils::SoNodeRef styleOwner(new SoDrawStyle);
	auto *style = static_cast<SoDrawStyle *>(styleOwner.get());
	style->lineWidth = 2.0f;
	edgeRoot->addChild(style);
	std::vector<SbVec3f> edgePoints;
	edgePoints.reserve(edges.size());
	for (const int32_t index : edges)
	    edgePoints.push_back(points[static_cast<size_t>(index)]);
	SbModernUtils::SoNodeRef coordinatesOwner(new SoCoordinate3);
	auto *coordinates = static_cast<SoCoordinate3 *>(coordinatesOwner.get());
	coordinates->point.setValues(0, static_cast<int>(edgePoints.size()),
	    edgePoints.data());
	edgeRoot->addChild(coordinates);
	SbModernUtils::SoNodeRef linesOwner(new SoLineSet);
	auto *lines = static_cast<SoLineSet *>(linesOwner.get());
	std::vector<int32_t> counts(edges.size() / 2, 2);
	lines->numVertices.setValues(0, static_cast<int>(counts.size()),
	    counts.data());
	edgeRoot->addChild(lines);
	geometry.children.push_back(edgeRoot);
	geometry.owners.push_back(std::move(edgeRootOwner));
    }

    SbModernUtils::SoNodeRef pointRootOwner(new SoSeparator);
    auto *pointRoot = static_cast<SoSeparator *>(pointRootOwner.get());
    SbModernUtils::SoNodeRef pointMaterialOwner(new SoMaterial);
    auto *pointMaterial = static_cast<SoMaterial *>(pointMaterialOwner.get());
    pointMaterial->diffuseColor = SbColor(1.0f, 0.65f, 0.15f);
    pointRoot->addChild(pointMaterial);
    SbModernUtils::SoNodeRef pointStyleOwner(new SoDrawStyle);
    auto *pointStyle = static_cast<SoDrawStyle *>(pointStyleOwner.get());
    pointStyle->pointSize = 9.0f;
    pointRoot->addChild(pointStyle);
    SbModernUtils::SoNodeRef pointCoordinatesOwner(new SoCoordinate3);
    auto *pointCoordinates =
	static_cast<SoCoordinate3 *>(pointCoordinatesOwner.get());
    pointCoordinates->point.setValues(0, this->pointCount(), points.data());
    pointRoot->addChild(pointCoordinates);
    SbModernUtils::SoNodeRef pointsOwner(new SoPointSet);
    auto *pointSet = static_cast<SoPointSet *>(pointsOwner.get());
    pointSet->numPoints = this->pointCount();
    pointRoot->addChild(pointSet);
    geometry.children.push_back(pointRoot);
    geometry.owners.push_back(std::move(pointRootOwner));

    const int emphasis[3] = {
	selectedIndex, hoverIndex, activeIndex
    };
    const SbColor emphasisColors[3] = {
	SbColor(1.0f, 1.0f, 1.0f), SbColor(1.0f, 1.0f, 0.2f),
	SbColor(0.2f, 1.0f, 1.0f)
    };
    for (int pass = 0; pass < 3; pass++) {
	const int index = emphasis[pass];
	SbVec3f representative;
	if (index < 0 || !this->featurePosition(domain, index, representative))
	    continue;
	SbModernUtils::SoNodeRef rootOwner(new SoSeparator);
	auto *root = static_cast<SoSeparator *>(rootOwner.get());
	SbModernUtils::SoNodeRef materialOwner(new SoMaterial);
	auto *material = static_cast<SoMaterial *>(materialOwner.get());
	material->diffuseColor = emphasisColors[pass];
	root->addChild(material);
	SbModernUtils::SoNodeRef styleOwner(new SoDrawStyle);
	auto *style = static_cast<SoDrawStyle *>(styleOwner.get());
	style->lineWidth = pass == 2 ? 5.0f : 4.0f;
	style->pointSize = pass == 2 ? 15.0f : 13.0f;
	root->addChild(style);
	std::vector<SbVec3f> featurePoints;
	if (domain == DOMAIN_VERTEX) {
	    featurePoints.push_back(representative);
	} else if (domain == DOMAIN_EDGE) {
	    for (size_t edge = 0; edge < edgeFeatures.size(); edge++) {
		if (edgeFeatures[edge] != index)
		    continue;
		const size_t offset = edge * 2;
		featurePoints.push_back(points[edges[offset]]);
		featurePoints.push_back(points[edges[offset + 1]]);
	    }
	} else if (domain == DOMAIN_FACE) {
	    const Face &face = faces[static_cast<size_t>(index)];
	    for (const int32_t vertex : face.vertices)
		featurePoints.push_back(points[vertex]);
	    featurePoints.push_back(points[face.vertices.front()]);
	}
	SbModernUtils::SoNodeRef coordinatesOwner(new SoCoordinate3);
	auto *coordinates = static_cast<SoCoordinate3 *>(coordinatesOwner.get());
	coordinates->point.setValues(0, static_cast<int>(featurePoints.size()),
	    featurePoints.data());
	root->addChild(coordinates);
	if (domain == DOMAIN_VERTEX) {
	    SbModernUtils::SoNodeRef pointOwner(new SoPointSet);
	    auto *point = static_cast<SoPointSet *>(pointOwner.get());
	    point->numPoints = 1;
	    root->addChild(point);
	} else {
	    SbModernUtils::SoNodeRef lineOwner(new SoLineSet);
	    auto *line = static_cast<SoLineSet *>(lineOwner.get());
	    if (domain == DOMAIN_EDGE) {
		std::vector<int32_t> counts(featurePoints.size() / 2, 2);
		line->numVertices.setValues(0, static_cast<int>(counts.size()),
		    counts.data());
	    } else {
		line->numVertices.set1Value(0,
		    static_cast<int32_t>(featurePoints.size()));
	    }
	    root->addChild(line);
	}
	geometry.children.push_back(root);
	geometry.owners.push_back(std::move(rootOwner));
    }

    return geometry;
}


void
SoBRLIndexedEditManipulator::rebuildGeometry(void)
{
    if ((!this->visible.getValue() || d->points.empty()) &&
	!this->getNumChildren())
	return;
    const ManipulatorGeometry geometry = d->buildGeometry(
	this->visible.getValue(),
	static_cast<Domain>(this->selectionDomain.getValue()),
	this->selectedIndex.getValue(), this->hoverIndex.getValue(),
	this->activeIndex.getValue());
    PreparedFieldNotifications<0> notifications(*this, {});
    publish_manipulator_geometry(*this, geometry, notifications, [] {});
}


SbBool
SoBRLIndexedEditManipulator::project(const SbVec3f &point, int width,
	int height, const SoCamera *camera, SbVec3f &pixel) const
{
    if (!camera || width <= 0 || height <= 0)
	return FALSE;
    const float aspect = static_cast<float>(width) /
	static_cast<float>(height);
    SbVec3f normalized;
    camera->getViewVolume(aspect).projectToScreen(point, normalized);
    if (!std::isfinite(normalized[0]) || !std::isfinite(normalized[1]) ||
	!std::isfinite(normalized[2]))
	return FALSE;
    pixel.setValue(normalized[0] * static_cast<float>(width),
	(1.0f - normalized[1]) * static_cast<float>(height), normalized[2]);
    return TRUE;
}


SbBool
SoBRLIndexedEditManipulator::Private::featurePosition(
	SoBRLIndexedEditManipulator::Domain domain, int index,
	SbVec3f &position) const
{
    position.setValue(0.0f, 0.0f, 0.0f);
    if (domain == DOMAIN_VERTEX) {
	if (index < 0 || index >= this->pointCount())
	    return FALSE;
	position = points[static_cast<size_t>(index)];
	return TRUE;
    }
    if (domain == DOMAIN_EDGE) {
	if (index < 0 || index >= this->edgeCount())
	    return FALSE;
	int count = 0;
	for (size_t edge = 0; edge < edgeFeatures.size(); edge++) {
	    if (edgeFeatures[edge] != index)
		continue;
	    const size_t offset = edge * 2;
	    position += (points[edges[offset]] +
		points[edges[offset + 1]]) * 0.5f;
	    count++;
	}
	if (!count)
	    return FALSE;
	position /= static_cast<float>(count);
	return TRUE;
    }
    if (domain == DOMAIN_FACE) {
	if (index < 0 || index >= this->faceCount())
	    return FALSE;
	const Face &face = faces[static_cast<size_t>(index)];
	for (const int32_t vertex : face.vertices)
	    position += points[vertex];
	position /= static_cast<float>(face.vertices.size());
	return TRUE;
    }
    return FALSE;
}


SbBool
SoBRLIndexedEditManipulator::featurePosition(Domain domain, int index,
	SbVec3f &position) const
{
    return d->featurePosition(domain, index, position);
}


int
SoBRLIndexedEditManipulator::hitTest(Domain domain, int x, int y, int width,
	int height, const SoCamera *camera, float radiusPixels) const
{
    if (!this->visible.getValue() || !camera || radiusPixels <= 0.0f)
	return -1;
    std::vector<SbVec3f> projected(d->points.size());
    for (size_t i = 0; i < d->points.size(); i++) {
	if (!this->project(d->points[i], width, height, camera, projected[i]))
	    return -1;
    }
    const float radiusSq = radiusPixels * radiusPixels;
    float bestDistance = std::numeric_limits<float>::max();
    float bestDepth = std::numeric_limits<float>::max();
    int best = -1;
    if (domain == DOMAIN_VERTEX) {
	for (int i = 0; i < this->pointCount(); i++) {
	    const float dx = projected[i][0] - static_cast<float>(x);
	    const float dy = projected[i][1] - static_cast<float>(y);
	    const float distance = dx * dx + dy * dy;
	    if (distance <= radiusSq && (distance < bestDistance ||
		(std::fabs(distance - bestDistance) <= 1.0e-6f &&
		 projected[i][2] < bestDepth))) {
		best = i;
		bestDistance = distance;
		bestDepth = projected[i][2];
	    }
	}
	return best;
    }

    if (domain == DOMAIN_EDGE) {
	for (size_t i = 0; i < d->edgeFeatures.size(); i++) {
	    const SbVec3f &a = projected[d->edges[i * 2]];
	    const SbVec3f &b = projected[d->edges[i * 2 + 1]];
	    const float dx = b[0] - a[0];
	    const float dy = b[1] - a[1];
	    const float lengthSq = dx * dx + dy * dy;
	    float t = lengthSq > 1.0e-6f ?
		((static_cast<float>(x) - a[0]) * dx +
		 (static_cast<float>(y) - a[1]) * dy) / lengthSq : 0.0f;
	    t = std::max(0.0f, std::min(1.0f, t));
	    const float ex = a[0] + t * dx - static_cast<float>(x);
	    const float ey = a[1] + t * dy - static_cast<float>(y);
	    const float distance = ex * ex + ey * ey;
	    const float depth = a[2] + t * (b[2] - a[2]);
	    if (distance <= radiusSq && (distance < bestDistance ||
		(std::fabs(distance - bestDistance) <= 1.0e-6f &&
		 depth < bestDepth))) {
		best = d->edgeFeatures[i];
		bestDistance = distance;
		bestDepth = depth;
	    }
	}
	return best;
    }
    if (domain == DOMAIN_FACE) {
	for (int fi = 0; fi < this->faceCount(); fi++) {
	    const Private::Face &face = d->faces[static_cast<size_t>(fi)];
	    bool inside = false;
	    float depth = 0.0f;
	    for (size_t i = 0, j = face.vertices.size() - 1;
		i < face.vertices.size(); j = i++) {
		const SbVec3f &pi = projected[face.vertices[i]];
		const SbVec3f &pj = projected[face.vertices[j]];
		depth += pi[2];
		const bool crosses = ((pi[1] > static_cast<float>(y)) !=
		    (pj[1] > static_cast<float>(y))) &&
		    (static_cast<float>(x) < (pj[0] - pi[0]) *
		    (static_cast<float>(y) - pi[1]) /
		    (pj[1] - pi[1]) + pi[0]);
		if (crosses)
		    inside = !inside;
	    }
	    depth /= static_cast<float>(face.vertices.size());
	    if (inside && depth < bestDepth) {
		best = fi;
		bestDepth = depth;
	    }
	}
    }
    return best;
}


SbBool
SoBRLIndexedEditManipulator::screenPosition(Domain domain, int index,
	int width, int height, const SoCamera *camera, int &x, int &y) const
{
    x = 0;
    y = 0;
    SbVec3f position;
    SbVec3f pixel;
    if (!this->featurePosition(domain, index, position) ||
	!this->project(position, width, height, camera, pixel))
	return FALSE;
    x = static_cast<int>(std::lround(pixel[0]));
    y = static_cast<int>(std::lround(pixel[1]));
    return TRUE;
}
