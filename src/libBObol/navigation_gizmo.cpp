/*              N A V I G A T I O N _ G I Z M O . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BNavigationGizmo.h"

#include <Inventor/SbViewportRegion.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/annex/HUD/nodekits/SoHUDKit.h>
#include <Inventor/elements/SoViewingMatrixElement.h>
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoFaceSet.h>
#include <Inventor/nodes/SoFont.h>
#include <Inventor/nodes/SoLightModel.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoShapeHints.h>
#include <Inventor/nodes/SoText2.h>
#include <Inventor/nodes/SoTransparencyType.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/sensors/SoFieldSensor.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <utility>
#include <vector>

SO_NODE_SOURCE(SoBRLNavigationGizmo);

namespace {

static const float gizmo_rad_to_deg = 57.2957795130823208768f;
static const float gizmo_pi = 3.14159265358979323846f;

struct GizmoFace {
    int normal[3];
    int vertices[4];
};

static const SbVec3f gizmo_vertices[8] = {
    SbVec3f(-1.0f, -1.0f, -1.0f),
    SbVec3f( 1.0f, -1.0f, -1.0f),
    SbVec3f( 1.0f,  1.0f, -1.0f),
    SbVec3f(-1.0f,  1.0f, -1.0f),
    SbVec3f(-1.0f, -1.0f,  1.0f),
    SbVec3f( 1.0f, -1.0f,  1.0f),
    SbVec3f( 1.0f,  1.0f,  1.0f),
    SbVec3f(-1.0f,  1.0f,  1.0f)
};

static const GizmoFace gizmo_faces[6] = {
    {{-1,  0,  0}, {0, 4, 7, 3}},
    {{ 1,  0,  0}, {1, 2, 6, 5}},
    {{ 0, -1,  0}, {0, 1, 5, 4}},
    {{ 0,  1,  0}, {3, 7, 6, 2}},
    {{ 0,  0, -1}, {0, 3, 2, 1}},
    {{ 0,  0,  1}, {4, 5, 6, 7}}
};

static SoBRLNavigationGizmo::Part
gizmo_part_from_direction(int x, int y, int z)
{
    if ((x < -1 || x > 1) || (y < -1 || y > 1) ||
	(z < -1 || z > 1) || (!x && !y && !z))
	return SoBRLNavigationGizmo::PART_NONE;
    const int code = 1 + (x + 1) * 9 + (y + 1) * 3 + (z + 1);
    return static_cast<SoBRLNavigationGizmo::Part>(code);
}

static SbColor
gizmo_scaled_color(const SbColor &color, float scale)
{
    return SbColor(std::min(1.0f, color[0] * scale),
	std::min(1.0f, color[1] * scale),
	std::min(1.0f, color[2] * scale));
}

static SbColor
gizmo_blend_color(const SbColor &a, const SbColor &b, float weight)
{
    const float aw = 1.0f - weight;
    return SbColor(std::min(1.0f, a[0] * aw + b[0] * weight),
	std::min(1.0f, a[1] * aw + b[1] * weight),
	std::min(1.0f, a[2] * aw + b[2] * weight));
}

static SbColor
gizmo_face_color(const SoBRLNavigationGizmo *gizmo,
	const GizmoFace &face, const SbVec3f &hover, const SbVec3f &active)
{
    SbColor color = face.normal[0] ? gizmo->xColor.getValue() :
	(face.normal[1] ? gizmo->yColor.getValue() : gizmo->zColor.getValue());
    const int sign = face.normal[0] ? face.normal[0] :
	(face.normal[1] ? face.normal[1] : face.normal[2]);
    color = gizmo_scaled_color(color, sign > 0 ? 0.90f : 0.58f);

    const int axis = face.normal[0] ? 0 : (face.normal[1] ? 1 : 2);
    const float faceSign = static_cast<float>(face.normal[axis]);
    if (active[axis] * faceSign > 0.5f)
	return gizmo_blend_color(color, gizmo->highlightColor.getValue(), 0.78f);
    if (hover[axis] * faceSign > 0.5f)
	return gizmo_blend_color(color, gizmo->highlightColor.getValue(), 0.48f);
    return color;
}

static SoSeparator *
gizmo_face_node(const SbVec3f points[4], const SbColor &color)
{
    SoSeparator *sep = new SoSeparator;
    SoBaseColor *base = new SoBaseColor;
    base->rgb = color;
    sep->addChild(base);

    SoShapeHints *hints = new SoShapeHints;
    hints->vertexOrdering = SoShapeHints::COUNTERCLOCKWISE;
    hints->shapeType = SoShapeHints::SOLID;
    sep->addChild(hints);

    SoCoordinate3 *coordinates = new SoCoordinate3;
    coordinates->point.setValues(0, 4, points);
    sep->addChild(coordinates);

    SoFaceSet *face = new SoFaceSet;
    face->numVertices.set1Value(0, 4);
    sep->addChild(face);
    return sep;
}

static SoSeparator *
gizmo_lines_node(const std::vector<SbVec3f> &points,
	const SbColor &color, float lineWidth)
{
    if (points.size() < 2 || points.size() % 2)
	return NULL;

    SoSeparator *sep = new SoSeparator;
    SoBaseColor *base = new SoBaseColor;
    base->rgb = color;
    sep->addChild(base);
    SoDrawStyle *style = new SoDrawStyle;
    style->lineWidth = lineWidth;
    sep->addChild(style);
    SoCoordinate3 *coordinates = new SoCoordinate3;
    coordinates->point.setValues(0, static_cast<int>(points.size()),
	points.data());
    sep->addChild(coordinates);
    SoLineSet *lines = new SoLineSet;
    for (size_t i = 0; i < points.size() / 2; i++)
	lines->numVertices.set1Value(static_cast<int>(i), 2);
    sep->addChild(lines);
    return sep;
}

static void
gizmo_add_ring(SoSeparator *root, const SbRotation &rotation,
	const SbVec3f &axisA, const SbVec3f &axisB, float radius,
	const SbColor &color, SbBool highlighted)
{
    if (!root)
	return;

    const int segments = 64;
    std::vector<SbVec3f> rear;
    std::vector<SbVec3f> front;
    for (int i = 0; i < segments; i++) {
	const float a0 = 2.0f * gizmo_pi * static_cast<float>(i) / segments;
	const float a1 = 2.0f * gizmo_pi * static_cast<float>(i + 1) /
	    segments;
	SbVec3f p0;
	SbVec3f p1;
	rotation.multVec((axisA * std::cos(a0) + axisB * std::sin(a0)) *
	    radius, p0);
	rotation.multVec((axisA * std::cos(a1) + axisB * std::sin(a1)) *
	    radius, p1);
	const float depth = (p0[2] + p1[2]) * 0.5f;
	p0[2] = p1[2] = depth >= 0.0f ? 0.08f : -0.08f;
	std::vector<SbVec3f> &target = depth >= 0.0f ? front : rear;
	target.push_back(p0);
	target.push_back(p1);
    }

    const float scale = highlighted ? 1.0f : 0.88f;
    if (!rear.empty())
	root->addChild(gizmo_lines_node(rear,
	    gizmo_scaled_color(color, scale * 0.48f),
	    highlighted ? 2.5f : 1.5f));
    if (!front.empty())
	root->addChild(gizmo_lines_node(front,
	    gizmo_scaled_color(color, scale),
	    highlighted ? 3.0f : 2.0f));
}

static SoSeparator *
gizmo_circle_panel(float radius, const SbColor &color)
{
    const int segments = 48;
    std::vector<SbVec3f> points;
    points.reserve(segments);
    for (int i = 0; i < segments; i++) {
	const float angle = 2.0f * gizmo_pi * static_cast<float>(i) /
	    segments;
	points.push_back(SbVec3f(radius * std::cos(angle),
	    radius * std::sin(angle), -0.6f));
    }

    SoSeparator *panel = new SoSeparator;
    SoTransparencyType *transparencyType = new SoTransparencyType;
    transparencyType->value = SoTransparencyType::BLEND;
    panel->addChild(transparencyType);
    SoMaterial *material = new SoMaterial;
    material->diffuseColor = color;
    material->transparency = 0.46f;
    panel->addChild(material);
    SoCoordinate3 *coordinates = new SoCoordinate3;
    coordinates->point.setValues(0, segments, points.data());
    panel->addChild(coordinates);
    SoFaceSet *face = new SoFaceSet;
    face->numVertices.set1Value(0, segments);
    panel->addChild(face);
    return panel;
}

static SoSeparator *
gizmo_endpoint_node(const SbVec3f &center, float radius,
	const SbColor &color)
{
    const int segments = 24;
    std::vector<SbVec3f> points;
    points.reserve(segments);
    for (int i = 0; i < segments; i++) {
	const float angle = 2.0f * gizmo_pi * static_cast<float>(i) /
	    segments;
	points.push_back(SbVec3f(center[0] + radius * std::cos(angle),
	    center[1] + radius * std::sin(angle), 0.22f));
    }

    SoSeparator *endpoint = new SoSeparator;
    SoBaseColor *base = new SoBaseColor;
    base->rgb = color;
    endpoint->addChild(base);
    SoCoordinate3 *coordinates = new SoCoordinate3;
    coordinates->point.setValues(0, segments, points.data());
    endpoint->addChild(coordinates);
    SoFaceSet *face = new SoFaceSet;
    face->numVertices.set1Value(0, segments);
    endpoint->addChild(face);
    return endpoint;
}

static void
gizmo_add_label(SoSeparator *root, const SbRotation &rotation,
	const SbVec3f &axis, float distance, float fontSize,
	const char *label, const SbColor &color)
{
    if (!root || !label)
	return;
    SbVec3f position;
    rotation.multVec(axis * distance, position);
    position[2] = 0.25f;

    SoSeparator *sep = new SoSeparator;
    SoTranslation *translation = new SoTranslation;
    translation->translation = position;
    sep->addChild(translation);
    SoBaseColor *base = new SoBaseColor;
    base->rgb = color;
    sep->addChild(base);
    SoFont *font = new SoFont;
    font->size = fontSize;
    sep->addChild(font);
    SoText2 *text = new SoText2;
    text->string.set1Value(0, label);
    text->justification = SoText2::CENTER;
    text->depthTest = FALSE;
    sep->addChild(text);
    root->addChild(sep);
}

} // namespace

SoBRLNavigationGizmo::SoBRLNavigationGizmo(void) :
    trackedCamera(NULL),
    cameraSensor(new SoFieldSensor(SoBRLNavigationGizmo::cameraChanged, this)),
    hudKit(NULL),
    anchorTranslation(NULL),
    rotation(SbRotation::identity())
{
    SO_NODE_CONSTRUCTOR(SoBRLNavigationGizmo);

    SO_NODE_DEFINE_ENUM_VALUE(Corner, LOWER_LEFT);
    SO_NODE_DEFINE_ENUM_VALUE(Corner, LOWER_RIGHT);
    SO_NODE_DEFINE_ENUM_VALUE(Corner, UPPER_LEFT);
    SO_NODE_DEFINE_ENUM_VALUE(Corner, UPPER_RIGHT);

    SO_NODE_DEFINE_ENUM_VALUE(Style, CUBE);
    SO_NODE_DEFINE_ENUM_VALUE(Style, CIRCLES);

    SO_NODE_ADD_FIELD(overlayId, ("faceplate::navigation_gizmo"));
    SO_NODE_ADD_FIELD(visible, (TRUE));
    SO_NODE_ADD_FIELD(style, (CUBE));
    SO_NODE_SET_SF_ENUM_TYPE(style, Style);
    SO_NODE_ADD_FIELD(corner, (UPPER_RIGHT));
    SO_NODE_SET_SF_ENUM_TYPE(corner, Corner);
    SO_NODE_ADD_FIELD(size, (96.0f));
    SO_NODE_ADD_FIELD(margin, (12.0f));
    SO_NODE_ADD_FIELD(xColor, (SbColor(0.88f, 0.22f, 0.20f)));
    SO_NODE_ADD_FIELD(yColor, (SbColor(0.28f, 0.78f, 0.30f)));
    SO_NODE_ADD_FIELD(zColor, (SbColor(0.25f, 0.48f, 0.94f)));
    SO_NODE_ADD_FIELD(panelColor, (SbColor(0.08f, 0.09f, 0.11f)));
    SO_NODE_ADD_FIELD(highlightColor, (SbColor(1.0f, 0.82f, 0.22f)));
    SO_NODE_ADD_FIELD(hoverPart, (PART_NONE));
    SO_NODE_ADD_FIELD(activePart, (PART_NONE));

    this->cameraSensor->setPriority(0);
}

SoBRLNavigationGizmo::~SoBRLNavigationGizmo(void)
{
    if (this->cameraSensor) {
	this->cameraSensor->detach();
	delete this->cameraSensor;
	this->cameraSensor = NULL;
    }
    if (this->trackedCamera) {
	this->trackedCamera->unref();
	this->trackedCamera = NULL;
    }
}

void
SoBRLNavigationGizmo::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLNavigationGizmo, SoSeparator, "Separator");
}

void
SoBRLNavigationGizmo::setCamera(SoCamera *camera)
{
    if (camera == this->trackedCamera)
	return;
    this->cameraSensor->detach();
    if (camera)
	camera->ref();
    if (this->trackedCamera)
	this->trackedCamera->unref();
    this->trackedCamera = camera;
    if (this->trackedCamera)
	this->cameraSensor->attach(&this->trackedCamera->orientation);
    this->syncCamera();
}

SoCamera *
SoBRLNavigationGizmo::getCamera(void) const
{
    return this->trackedCamera;
}

void
SoBRLNavigationGizmo::cameraChanged(void *userData, SoSensor *UNUSED(sensor))
{
    SoBRLNavigationGizmo *gizmo =
	static_cast<SoBRLNavigationGizmo *>(userData);
    if (gizmo)
	gizmo->syncCamera();
}

void
SoBRLNavigationGizmo::syncCamera(void)
{
    this->rotation = this->trackedCamera ?
	this->trackedCamera->orientation.getValue().inverse() :
	SbRotation::identity();
    if (this->hudKit)
	this->rebuildGeometry();
}

SbRotation
SoBRLNavigationGizmo::displayRotation(void) const
{
    return this->rotation;
}

SoHUDKit *
SoBRLNavigationGizmo::rebuildGeometry(void)
{
    this->hudKit = NULL;
    this->anchorTranslation = NULL;
    this->removeAllChildren();
    if (!this->visible.getValue())
	return NULL;

    const float controlSize = std::max(40.0f, this->size.getValue());
    const Style gizmoStyle = static_cast<Style>(this->style.getValue());
    const float cubeHalf = controlSize * 0.235f;
    const float axisLength = controlSize * 0.43f;
    const float panelHalf = controlSize * 0.50f;
    const SbRotation viewRotation = this->displayRotation();

    SoHUDKit *hud = new SoHUDKit;
    SoSeparator *widget = new SoSeparator;
    this->anchorTranslation = new SoTranslation;
    widget->addChild(this->anchorTranslation);

    SoLightModel *lightModel = new SoLightModel;
    lightModel->model = SoLightModel::BASE_COLOR;
    widget->addChild(lightModel);

    if (gizmoStyle == CIRCLES) {
	widget->addChild(gizmo_circle_panel(panelHalf,
	    this->panelColor.getValue()));
    } else {
	SbVec3f panelPoints[4] = {
	    SbVec3f(-panelHalf, -panelHalf, -0.6f),
	    SbVec3f( panelHalf, -panelHalf, -0.6f),
	    SbVec3f( panelHalf,  panelHalf, -0.6f),
	    SbVec3f(-panelHalf,  panelHalf, -0.6f)
	};
	SoSeparator *panel = new SoSeparator;
	SoTransparencyType *transparencyType = new SoTransparencyType;
	transparencyType->value = SoTransparencyType::BLEND;
	panel->addChild(transparencyType);
	SoMaterial *panelMaterial = new SoMaterial;
	panelMaterial->diffuseColor = this->panelColor.getValue();
	panelMaterial->transparency = 0.32f;
	panel->addChild(panelMaterial);
	SoCoordinate3 *panelCoordinates = new SoCoordinate3;
	panelCoordinates->point.setValues(0, 4, panelPoints);
	panel->addChild(panelCoordinates);
	SoFaceSet *panelFace = new SoFaceSet;
	panelFace->numVertices.set1Value(0, 4);
	panel->addChild(panelFace);
	widget->addChild(panel);
    }

    SbVec3f hover(0.0f, 0.0f, 0.0f);
    SbVec3f active(0.0f, 0.0f, 0.0f);
    (void)partDirection(static_cast<Part>(this->hoverPart.getValue()), hover);
    (void)partDirection(static_cast<Part>(this->activePart.getValue()), active);

    struct VisibleFace {
	int index;
	float depth;
    };
    if (gizmoStyle == CUBE) {
	std::vector<VisibleFace> visibleFaces;
	for (int i = 0; i < 6; i++) {
	    SbVec3f normal(static_cast<float>(gizmo_faces[i].normal[0]),
		static_cast<float>(gizmo_faces[i].normal[1]),
		static_cast<float>(gizmo_faces[i].normal[2]));
	    SbVec3f viewNormal;
	    viewRotation.multVec(normal, viewNormal);
	    if (viewNormal[2] <= 1.0e-5f)
		continue;
	    VisibleFace face = {i, viewNormal[2]};
	    visibleFaces.push_back(face);
	}
	std::sort(visibleFaces.begin(), visibleFaces.end(),
	    [](const VisibleFace &a, const VisibleFace &b) {
		return a.depth < b.depth;
	    });

	std::set<std::pair<int, int>> visibleEdges;
	for (const VisibleFace &visibleFace : visibleFaces) {
	    const GizmoFace &face = gizmo_faces[visibleFace.index];
	    SbVec3f points[4];
	    for (int v = 0; v < 4; v++) {
		viewRotation.multVec(gizmo_vertices[face.vertices[v]] *
		    cubeHalf, points[v]);
		points[v][2] = -0.25f;
		const int a = face.vertices[v];
		const int b = face.vertices[(v + 1) % 4];
		visibleEdges.insert(std::make_pair(std::min(a, b),
		    std::max(a, b)));
	    }
	    widget->addChild(gizmo_face_node(points,
		gizmo_face_color(this, face, hover, active)));
	}

	std::vector<SbVec3f> edgePoints;
	for (const std::pair<int, int> &edge : visibleEdges) {
	    SbVec3f a;
	    SbVec3f b;
	    viewRotation.multVec(gizmo_vertices[edge.first] * cubeHalf, a);
	    viewRotation.multVec(gizmo_vertices[edge.second] * cubeHalf, b);
	    a[2] = b[2] = 0.05f;
	    edgePoints.push_back(a);
	    edgePoints.push_back(b);
	}
	SoSeparator *edges = new SoSeparator;
	SoBaseColor *edgeColor = new SoBaseColor;
	edgeColor->rgb = (this->hoverPart.getValue() != PART_NONE ||
	    this->activePart.getValue() != PART_NONE) ?
	    gizmo_scaled_color(this->highlightColor.getValue(), 0.85f) :
	    SbColor(0.84f, 0.84f, 0.86f);
	edges->addChild(edgeColor);
	SoDrawStyle *edgeStyle = new SoDrawStyle;
	edgeStyle->lineWidth = 1.5f;
	edges->addChild(edgeStyle);
	SoCoordinate3 *edgeCoordinates = new SoCoordinate3;
	edgeCoordinates->point.setValues(0,
	    static_cast<int>(edgePoints.size()), edgePoints.data());
	edges->addChild(edgeCoordinates);
	SoLineSet *edgeLines = new SoLineSet;
	for (size_t i = 0; i < edgePoints.size() / 2; i++)
	    edgeLines->numVertices.set1Value(static_cast<int>(i), 2);
	edges->addChild(edgeLines);
	widget->addChild(edges);
    } else {
	const SbBool orbitHighlight =
	    this->hoverPart.getValue() == PART_ORBIT ||
	    this->activePart.getValue() == PART_ORBIT;
	gizmo_add_ring(widget, viewRotation,
	    SbVec3f(1.0f, 0.0f, 0.0f), SbVec3f(0.0f, 1.0f, 0.0f),
	    controlSize * 0.35f, this->zColor.getValue(), orbitHighlight);
	gizmo_add_ring(widget, viewRotation,
	    SbVec3f(1.0f, 0.0f, 0.0f), SbVec3f(0.0f, 0.0f, 1.0f),
	    controlSize * 0.35f, this->yColor.getValue(), orbitHighlight);
	gizmo_add_ring(widget, viewRotation,
	    SbVec3f(0.0f, 1.0f, 0.0f), SbVec3f(0.0f, 0.0f, 1.0f),
	    controlSize * 0.35f, this->xColor.getValue(), orbitHighlight);
    }

    SbVec3f axisPoints[6];
    const SbVec3f axes[3] = {
	SbVec3f(1.0f, 0.0f, 0.0f),
	SbVec3f(0.0f, 1.0f, 0.0f),
	SbVec3f(0.0f, 0.0f, 1.0f)
    };
    for (int i = 0; i < 3; i++) {
	viewRotation.multVec(axes[i] * (gizmoStyle == CUBE ?
	    cubeHalf * 1.08f : -axisLength), axisPoints[i * 2]);
	viewRotation.multVec(axes[i] * axisLength, axisPoints[i * 2 + 1]);
	axisPoints[i * 2][2] = axisPoints[i * 2 + 1][2] = 0.12f;
    }
    const SbColor axisColors[3] = {this->xColor.getValue(),
	this->yColor.getValue(), this->zColor.getValue()};
    for (int i = 0; i < 3; i++) {
	SoSeparator *axisSep = new SoSeparator;
	SoBaseColor *axisColor = new SoBaseColor;
	axisColor->rgb = axisColors[i];
	axisSep->addChild(axisColor);
	SoDrawStyle *axisStyle = new SoDrawStyle;
	axisStyle->lineWidth = 2.5f;
	axisSep->addChild(axisStyle);
	SoCoordinate3 *axisCoordinates = new SoCoordinate3;
	axisCoordinates->point.setValues(0, 2, &axisPoints[i * 2]);
	axisSep->addChild(axisCoordinates);
	SoLineSet *axisLine = new SoLineSet;
	axisLine->numVertices.set1Value(0, 2);
	axisSep->addChild(axisLine);
	widget->addChild(axisSep);
	if (gizmoStyle == CIRCLES) {
	    SbVec3f positive;
	    SbVec3f negative;
	    viewRotation.multVec(axes[i] * axisLength, positive);
	    viewRotation.multVec(axes[i] * -axisLength, negative);
	    const SbColor positiveColor =
		(active[i] > 0.5f || hover[i] > 0.5f) ?
		gizmo_blend_color(axisColors[i],
		    this->highlightColor.getValue(),
		    active[i] > 0.5f ? 0.80f : 0.52f) : axisColors[i];
	    const SbColor negativeBase = gizmo_scaled_color(axisColors[i], 0.48f);
	    const SbColor negativeColor =
		(active[i] < -0.5f || hover[i] < -0.5f) ?
		gizmo_blend_color(negativeBase,
		    this->highlightColor.getValue(),
		    active[i] < -0.5f ? 0.80f : 0.52f) : negativeBase;
	    widget->addChild(gizmo_endpoint_node(positive,
		controlSize * 0.045f, positiveColor));
	    widget->addChild(gizmo_endpoint_node(negative,
		controlSize * 0.035f, negativeColor));
	}
    }

    const float fontSize = std::max(10.0f, controlSize * 0.13f);
    gizmo_add_label(widget, viewRotation, axes[0], axisLength + fontSize * 0.35f,
	fontSize, "X", axisColors[0]);
    gizmo_add_label(widget, viewRotation, axes[1], axisLength + fontSize * 0.35f,
	fontSize, "Y", axisColors[1]);
    gizmo_add_label(widget, viewRotation, axes[2], axisLength + fontSize * 0.35f,
	fontSize, "Z", axisColors[2]);

    hud->addWidget(widget);
    this->addChild(hud);
    this->hudKit = hud;
    return hud;
}

SoHUDKit *
SoBRLNavigationGizmo::getHUDKit(void) const
{
    return this->hudKit;
}

void
SoBRLNavigationGizmo::setHoverPart(Part part)
{
    if (this->hoverPart.getValue() == static_cast<int>(part))
	return;
    this->hoverPart = static_cast<int>(part);
    this->rebuildGeometry();
}

void
SoBRLNavigationGizmo::setActivePart(Part part)
{
    if (this->activePart.getValue() == static_cast<int>(part))
	return;
    this->activePart = static_cast<int>(part);
    this->rebuildGeometry();
}

SbBool
SoBRLNavigationGizmo::partDirection(Part part, SbVec3f &direction)
{
    const int code = static_cast<int>(part);
    if (code < 1 || code > 27 || code == 14) {
	direction.setValue(0.0f, 0.0f, 0.0f);
	return FALSE;
    }
    const int value = code - 1;
    const int x = value / 9 - 1;
    const int y = (value % 9) / 3 - 1;
    const int z = value % 3 - 1;
    direction.setValue(static_cast<float>(x), static_cast<float>(y),
	static_cast<float>(z));
    return TRUE;
}

SbBool
SoBRLNavigationGizmo::partAet(Part part, float &azimuth, float &elevation)
{
    SbVec3f direction;
    if (!partDirection(part, direction))
	return FALSE;
    const float horizontal = std::hypot(direction[0], direction[1]);
    if (horizontal <= 1.0e-6f)
	azimuth = 270.0f;
    else {
	azimuth = std::atan2(direction[1], direction[0]) *
	    gizmo_rad_to_deg;
	if (azimuth < 0.0f)
	    azimuth += 360.0f;
    }
    elevation = std::atan2(direction[2], horizontal) *
	gizmo_rad_to_deg;
    return TRUE;
}

void
SoBRLNavigationGizmo::updateAnchor(int width, int height)
{
    if (!this->anchorTranslation || width <= 0 || height <= 0)
	return;
    const float controlSize = std::max(40.0f, this->size.getValue());
    const float offset = std::max(0.0f, this->margin.getValue()) +
	controlSize * 0.5f;
    float x = offset;
    float y = offset;
    switch (static_cast<Corner>(this->corner.getValue())) {
	case LOWER_RIGHT:
	    x = static_cast<float>(width) - offset;
	    break;
	case UPPER_LEFT:
	    y = static_cast<float>(height) - offset;
	    break;
	case UPPER_RIGHT:
	    x = static_cast<float>(width) - offset;
	    y = static_cast<float>(height) - offset;
	    break;
	case LOWER_LEFT:
	default:
	    break;
    }
    this->anchorTranslation->translation = SbVec3f(x, y, 0.0f);
}

void
SoBRLNavigationGizmo::GLRender(SoGLRenderAction *action)
{
    if (action) {
	/* A screen overlay is rendered after the CAD separator has restored its
	 * traversal state, so its viewing-matrix element may legitimately be the
	 * identity.  The assigned camera is authoritative when present.  The
	 * action matrix remains a useful fallback when the node is embedded
	 * directly in a camera-bearing scene. */
	const SbRotation renderRotation = this->trackedCamera ?
	    this->trackedCamera->orientation.getValue().inverse() :
	    SbRotation(SoViewingMatrixElement::get(action->getState()));
	if (!renderRotation.equals(this->rotation, 1.0e-6f)) {
	    this->rotation = renderRotation;
	    (void)this->rebuildGeometry();
	}
	const SbVec2s viewport =
	    action->getViewportRegion().getViewportSizePixels();
	this->updateAnchor(static_cast<int>(viewport[0]),
	    static_cast<int>(viewport[1]));
    }
    inherited::GLRender(action);
}

SoBRLNavigationGizmo::Part
SoBRLNavigationGizmo::hitTest(int x, int y, int width, int height) const
{
    if (!this->visible.getValue() || width <= 0 || height <= 0)
	return PART_NONE;

    const float controlSize = std::max(40.0f, this->size.getValue());
    const float offset = std::max(0.0f, this->margin.getValue()) +
	controlSize * 0.5f;
    float centerX = offset;
    float centerY = offset;
    switch (static_cast<Corner>(this->corner.getValue())) {
	case LOWER_RIGHT:
	    centerX = static_cast<float>(width) - offset;
	    break;
	case UPPER_LEFT:
	    centerY = static_cast<float>(height) - offset;
	    break;
	case UPPER_RIGHT:
	    centerX = static_cast<float>(width) - offset;
	    centerY = static_cast<float>(height) - offset;
	    break;
	case LOWER_LEFT:
	default:
	    break;
    }

    const float px = static_cast<float>(x);
    const float py = static_cast<float>(height - y);
    if (static_cast<Style>(this->style.getValue()) == CIRCLES) {
	const float localX = px - centerX;
	const float localY = py - centerY;
	const float axisLength = controlSize * 0.43f;
	const float pickRadius = controlSize * 0.105f;
	const SbVec3f axes[3] = {
	    SbVec3f(1.0f, 0.0f, 0.0f),
	    SbVec3f(0.0f, 1.0f, 0.0f),
	    SbVec3f(0.0f, 0.0f, 1.0f)
	};
	const Part positiveParts[3] = {PART_POS_X, PART_POS_Y, PART_POS_Z};
	const Part negativeParts[3] = {PART_NEG_X, PART_NEG_Y, PART_NEG_Z};
	float bestDistance = std::numeric_limits<float>::max();
	float bestDepth = -std::numeric_limits<float>::max();
	Part bestPart = PART_NONE;
	for (int sign = -1; sign <= 1; sign += 2) {
	    for (int axis = 0; axis < 3; axis++) {
		SbVec3f endpoint;
		this->displayRotation().multVec(axes[axis] *
		    axisLength * static_cast<float>(sign), endpoint);
		const float dx = localX - endpoint[0];
		const float dy = localY - endpoint[1];
		const float distance = std::hypot(dx, dy);
		if (distance <= pickRadius &&
		    (distance < bestDistance - 0.5f ||
		    (std::fabs(distance - bestDistance) <= 0.5f &&
		    endpoint[2] > bestDepth))) {
		    bestDistance = distance;
		    bestDepth = endpoint[2];
		    bestPart = sign > 0 ? positiveParts[axis] :
			negativeParts[axis];
		}
	    }
	}
	if (bestPart != PART_NONE)
	    return bestPart;
	return std::hypot(localX, localY) <= controlSize * 0.47f ?
	    PART_ORBIT : PART_NONE;
    }

    const float cubeHalf = controlSize * 0.235f;
    if (std::fabs(px - centerX) > cubeHalf * 1.75f ||
	std::fabs(py - centerY) > cubeHalf * 1.75f)
	return PART_NONE;

    SbVec3f rayOrigin((px - centerX) / cubeHalf,
	(py - centerY) / cubeHalf, 4.0f);
    SbVec3f rayDirection(0.0f, 0.0f, -1.0f);
    const SbRotation inverse = this->displayRotation().inverse();
    SbVec3f localOrigin;
    SbVec3f localDirection;
    inverse.multVec(rayOrigin, localOrigin);
    inverse.multVec(rayDirection, localDirection);

    float tMinimum = 0.0f;
    float tMaximum = std::numeric_limits<float>::max();
    for (int axis = 0; axis < 3; axis++) {
	if (std::fabs(localDirection[axis]) < 1.0e-7f) {
	    if (localOrigin[axis] < -1.0f || localOrigin[axis] > 1.0f)
		return PART_NONE;
	    continue;
	}
	float t1 = (-1.0f - localOrigin[axis]) / localDirection[axis];
	float t2 = ( 1.0f - localOrigin[axis]) / localDirection[axis];
	if (t1 > t2)
	    std::swap(t1, t2);
	tMinimum = std::max(tMinimum, t1);
	tMaximum = std::min(tMaximum, t2);
	if (tMaximum < tMinimum)
	    return PART_NONE;
    }

    const SbVec3f hit = localOrigin + localDirection * tMinimum;
    const float edgeThreshold = 0.68f;
    int direction[3] = {0, 0, 0};
    int largestAxis = 0;
    for (int axis = 0; axis < 3; axis++) {
	if (std::fabs(hit[axis]) > std::fabs(hit[largestAxis]))
	    largestAxis = axis;
	if (std::fabs(hit[axis]) >= edgeThreshold)
	    direction[axis] = hit[axis] < 0.0f ? -1 : 1;
    }
    if (!direction[0] && !direction[1] && !direction[2])
	direction[largestAxis] = hit[largestAxis] < 0.0f ? -1 : 1;
    return gizmo_part_from_direction(direction[0], direction[1],
	direction[2]);
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
