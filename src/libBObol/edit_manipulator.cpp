/*              E D I T _ M A N I P U L A T O R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BEditManipulator.h"

#include <Inventor/SbViewVolume.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoLightModel.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoPointSet.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

SO_NODE_SOURCE(SoBRLEditManipulator);
SO_NODE_SOURCE(SoBRLIndexedEditManipulator);

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
    editCenter = nextCenter;
    editAxes[0] = axisA;
    editAxes[1] = axisB;
    editAxes[2] = axisC;
    this->rebuildGeometry();
}

void
SoBRLEditManipulator::setVisible(SbBool value)
{
    if (this->visible.getValue() == value)
	return;
    this->visible = value;
    this->rebuildGeometry();
}

void
SoBRLEditManipulator::setHoverHandle(Handle handle)
{
    if (this->hoverHandle.getValue() == static_cast<int>(handle))
	return;
    this->hoverHandle = static_cast<int>(handle);
    this->rebuildGeometry();
}

void
SoBRLEditManipulator::setActiveHandle(Handle handle)
{
    if (this->activeHandle.getValue() == static_cast<int>(handle))
	return;
    this->activeHandle = static_cast<int>(handle);
    this->rebuildGeometry();
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
    this->removeAllChildren();
    if (!this->visible.getValue())
	return;

    SoLightModel *lightModel = new SoLightModel;
    lightModel->model = SoLightModel::BASE_COLOR;
    this->addChild(lightModel);

    const SbColor colors[3] = {
	SbColor(1.0f, 0.2f, 0.2f),
	SbColor(0.2f, 1.0f, 0.2f),
	SbColor(0.25f, 0.55f, 1.0f)
    };
    for (int i = 0; i < 3; i++) {
	SoSeparator *axisRoot = new SoSeparator;
	SoMaterial *material = new SoMaterial;
	SbColor color = colors[i];
	if (this->activeHandle.getValue() == i)
	    color = SbColor(1.0f, 1.0f, 1.0f);
	else if (this->hoverHandle.getValue() == i)
	    color = color + SbColor(0.25f, 0.25f, 0.25f);
	color[0] = std::min(color[0], 1.0f);
	color[1] = std::min(color[1], 1.0f);
	color[2] = std::min(color[2], 1.0f);
	material->diffuseColor = color;
	axisRoot->addChild(material);

	SoDrawStyle *style = new SoDrawStyle;
	style->lineWidth = this->activeHandle.getValue() == i ? 4.0f : 2.0f;
	style->pointSize = this->activeHandle.getValue() == i ? 13.0f :
	    (this->hoverHandle.getValue() == i ? 12.0f : 10.0f);
	axisRoot->addChild(style);

	SoCoordinate3 *coordinates = new SoCoordinate3;
	const SbVec3f values[2] = {editCenter, editCenter + editAxes[i]};
	coordinates->point.setValues(0, 2, values);
	axisRoot->addChild(coordinates);

	SoLineSet *line = new SoLineSet;
	line->numVertices.set1Value(0, 2);
	axisRoot->addChild(line);

	SoPointSet *point = new SoPointSet;
	point->startIndex = 1;
	point->numPoints = 1;
	axisRoot->addChild(point);
	this->addChild(axisRoot);
    }
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
    d->points.clear();
    d->edges.clear();
    d->edgeFeatures.clear();
    d->edgeFeatureCount = 0;
    d->pointFeatureCount = 0;
    d->faces.clear();
    if (points && nextPointCount > 0)
	d->points.assign(points, points + nextPointCount);
    d->pointFeatureCount = nextVertexFeatureCount < 0 ? nextPointCount :
	std::max(0, std::min(nextPointCount, nextVertexFeatureCount));
    if (edgeIndices && nextEdgeCount > 0) {
	for (int i = 0; i < nextEdgeCount; i++) {
	    const int32_t a = edgeIndices[i * 2];
	    const int32_t b = edgeIndices[i * 2 + 1];
	    if (a < 0 || b < 0 || a >= nextPointCount || b >= nextPointCount ||
		a == b)
		continue;
	    d->edges.push_back(a);
	    d->edges.push_back(b);
	    const int32_t feature = edgeFeatureIndices ?
		edgeFeatureIndices[i] : static_cast<int32_t>(i);
	    if (feature < 0) {
		d->edges.resize(d->edges.size() - 2);
		continue;
	    }
	    d->edgeFeatures.push_back(feature);
	    d->edgeFeatureCount = std::max(d->edgeFeatureCount,
		static_cast<int>(feature) + 1);
	}
    }
    if (faceIndices && faceVertexCounts && nextFaceCount > 0) {
	int offset = 0;
	for (int fi = 0; fi < nextFaceCount; fi++) {
	    const int count = faceVertexCounts[fi];
	    Private::Face face;
	    if (count >= 3) {
		for (int vi = 0; vi < count; vi++) {
		    const int32_t vertex = faceIndices[offset + vi];
		    if (vertex < 0 || vertex >= nextPointCount) {
			face.vertices.clear();
			break;
		    }
		    face.vertices.push_back(vertex);
		}
	    }
	    offset += std::max(0, count);
	    if (face.vertices.size() >= 3)
		d->faces.push_back(face);
	}
    }
    if (this->selectedIndex.getValue() >= this->pointCount() &&
	this->selectionDomain.getValue() == DOMAIN_VERTEX)
	this->selectedIndex = -1;
    if (this->selectedIndex.getValue() >= this->edgeCount() &&
	this->selectionDomain.getValue() == DOMAIN_EDGE)
	this->selectedIndex = -1;
    if (this->selectedIndex.getValue() >= this->faceCount() &&
	this->selectionDomain.getValue() == DOMAIN_FACE)
	this->selectedIndex = -1;
    this->rebuildGeometry();
}


void
SoBRLIndexedEditManipulator::setVisible(SbBool value)
{
    if (this->visible.getValue() == value)
	return;
    this->visible = value;
    this->rebuildGeometry();
}


void
SoBRLIndexedEditManipulator::setSelectionDomain(Domain domain)
{
    if (this->selectionDomain.getValue() == static_cast<int>(domain))
	return;
    this->selectionDomain = static_cast<int>(domain);
    this->selectedIndex = -1;
    this->hoverIndex = -1;
    this->activeIndex = -1;
    this->rebuildGeometry();
}


void
SoBRLIndexedEditManipulator::setSelectedIndex(int index)
{
    if (this->selectedIndex.getValue() == index)
	return;
    this->selectedIndex = index;
    this->rebuildGeometry();
}


void
SoBRLIndexedEditManipulator::setHoverIndex(int index)
{
    if (this->hoverIndex.getValue() == index)
	return;
    this->hoverIndex = index;
    this->rebuildGeometry();
}


void
SoBRLIndexedEditManipulator::setActiveIndex(int index)
{
    if (this->activeIndex.getValue() == index)
	return;
    this->activeIndex = index;
    this->rebuildGeometry();
}


int
SoBRLIndexedEditManipulator::pointCount(void) const
{
    return d->pointFeatureCount;
}


int
SoBRLIndexedEditManipulator::edgeCount(void) const
{
    return d->edgeFeatureCount;
}


int
SoBRLIndexedEditManipulator::faceCount(void) const
{
    return static_cast<int>(d->faces.size());
}


void
SoBRLIndexedEditManipulator::rebuildGeometry(void)
{
    this->removeAllChildren();
    if (!this->visible.getValue() || d->points.empty())
	return;

    SoLightModel *lightModel = new SoLightModel;
    lightModel->model = SoLightModel::BASE_COLOR;
    this->addChild(lightModel);

    if (!d->edges.empty()) {
	SoSeparator *edgeRoot = new SoSeparator;
	SoMaterial *material = new SoMaterial;
	material->diffuseColor = SbColor(0.35f, 0.65f, 1.0f);
	edgeRoot->addChild(material);
	SoDrawStyle *style = new SoDrawStyle;
	style->lineWidth = 2.0f;
	edgeRoot->addChild(style);
	std::vector<SbVec3f> edgePoints;
	edgePoints.reserve(d->edges.size());
	for (const int32_t index : d->edges)
	    edgePoints.push_back(d->points[static_cast<size_t>(index)]);
	SoCoordinate3 *coordinates = new SoCoordinate3;
	coordinates->point.setValues(0, static_cast<int>(edgePoints.size()),
	    edgePoints.data());
	edgeRoot->addChild(coordinates);
	SoLineSet *lines = new SoLineSet;
	std::vector<int32_t> counts(d->edges.size() / 2, 2);
	lines->numVertices.setValues(0, static_cast<int>(counts.size()),
	    counts.data());
	edgeRoot->addChild(lines);
	this->addChild(edgeRoot);
    }

    SoSeparator *pointRoot = new SoSeparator;
    SoMaterial *pointMaterial = new SoMaterial;
    pointMaterial->diffuseColor = SbColor(1.0f, 0.65f, 0.15f);
    pointRoot->addChild(pointMaterial);
    SoDrawStyle *pointStyle = new SoDrawStyle;
    pointStyle->pointSize = 9.0f;
    pointRoot->addChild(pointStyle);
    SoCoordinate3 *pointCoordinates = new SoCoordinate3;
    pointCoordinates->point.setValues(0, this->pointCount(), d->points.data());
    pointRoot->addChild(pointCoordinates);
    SoPointSet *points = new SoPointSet;
    points->numPoints = this->pointCount();
    pointRoot->addChild(points);
    this->addChild(pointRoot);

    const Domain domain = static_cast<Domain>(this->selectionDomain.getValue());
    const int emphasis[3] = {
	this->selectedIndex.getValue(), this->hoverIndex.getValue(),
	this->activeIndex.getValue()
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
	SoSeparator *root = new SoSeparator;
	SoMaterial *material = new SoMaterial;
	material->diffuseColor = emphasisColors[pass];
	root->addChild(material);
	SoDrawStyle *style = new SoDrawStyle;
	style->lineWidth = pass == 2 ? 5.0f : 4.0f;
	style->pointSize = pass == 2 ? 15.0f : 13.0f;
	root->addChild(style);
	std::vector<SbVec3f> featurePoints;
	if (domain == DOMAIN_VERTEX) {
	    featurePoints.push_back(representative);
	} else if (domain == DOMAIN_EDGE) {
	    for (size_t edge = 0; edge < d->edgeFeatures.size(); edge++) {
		if (d->edgeFeatures[edge] != index)
		    continue;
		const size_t offset = edge * 2;
		featurePoints.push_back(d->points[d->edges[offset]]);
		featurePoints.push_back(d->points[d->edges[offset + 1]]);
	    }
	} else if (domain == DOMAIN_FACE) {
	    const Private::Face &face = d->faces[static_cast<size_t>(index)];
	    for (const int32_t vertex : face.vertices)
		featurePoints.push_back(d->points[vertex]);
	    featurePoints.push_back(d->points[face.vertices.front()]);
	}
	SoCoordinate3 *coordinates = new SoCoordinate3;
	coordinates->point.setValues(0, static_cast<int>(featurePoints.size()),
	    featurePoints.data());
	root->addChild(coordinates);
	if (domain == DOMAIN_VERTEX) {
	    SoPointSet *point = new SoPointSet;
	    point->numPoints = 1;
	    root->addChild(point);
	} else {
	    SoLineSet *line = new SoLineSet;
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
	this->addChild(root);
    }
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
SoBRLIndexedEditManipulator::featurePosition(Domain domain, int index,
	SbVec3f &position) const
{
    position.setValue(0.0f, 0.0f, 0.0f);
    if (domain == DOMAIN_VERTEX) {
	if (index < 0 || index >= this->pointCount())
	    return FALSE;
	position = d->points[static_cast<size_t>(index)];
	return TRUE;
    }
    if (domain == DOMAIN_EDGE) {
	if (index < 0 || index >= this->edgeCount())
	    return FALSE;
	int count = 0;
	for (size_t edge = 0; edge < d->edgeFeatures.size(); edge++) {
	    if (d->edgeFeatures[edge] != index)
		continue;
	    const size_t offset = edge * 2;
	    position += (d->points[d->edges[offset]] +
		d->points[d->edges[offset + 1]]) * 0.5f;
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
	const Private::Face &face = d->faces[static_cast<size_t>(index)];
	for (const int32_t vertex : face.vertices)
	    position += d->points[vertex];
	position /= static_cast<float>(face.vertices.size());
	return TRUE;
    }
    return FALSE;
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
