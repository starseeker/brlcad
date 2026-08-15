/*            T E S T _ E D I T _ M A N I P U L A T O R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BEditManipulator.h"
#include "BObol/BInit.h"

#include <Inventor/nodes/SoOrthographicCamera.h>

#include <cmath>
#include <cstdio>

static int failures = 0;

#define CHECK(_condition, _message) do { \
    if (!(_condition)) { \
	std::fprintf(stderr, "FAIL: %s\n", (_message)); \
	failures++; \
    } \
} while (0)

int
main(void)
{
    bobol_init(NULL);

    SoOrthographicCamera *camera = new SoOrthographicCamera;
    camera->ref();
    camera->position = SbVec3f(0.0f, 0.0f, 10.0f);
    camera->focalDistance = 10.0f;
    camera->height = 10.0f;

    SoBRLEditManipulator *manipulator = new SoBRLEditManipulator;
    manipulator->ref();
    manipulator->setEllipsoidAxes(SbVec3f(0.0f, 0.0f, 0.0f),
	SbVec3f(2.0f, 0.0f, 0.0f), SbVec3f(0.0f, 1.0f, 0.0f),
	SbVec3f(0.0f, 0.0f, 1.0f));

    CHECK(manipulator->getNumChildren() > 0,
	"visible manipulator realizes retained handle geometry");
    CHECK(manipulator->hitTest(520, 300, 800, 600, camera) ==
	SoBRLEditManipulator::HANDLE_AXIS_A,
	"screen hit testing identifies the ellipsoid A handle");
    CHECK(manipulator->hitTest(400, 240, 800, 600, camera) ==
	SoBRLEditManipulator::HANDLE_AXIS_B,
	"screen hit testing identifies the ellipsoid B handle");
    CHECK(manipulator->hitTest(40, 40, 800, 600, camera) ==
	SoBRLEditManipulator::HANDLE_NONE,
	"screen hit testing rejects points outside the handles");

    float factor = 0.0f;
    CHECK(manipulator->projectedScale(
	SoBRLEditManipulator::HANDLE_AXIS_A, 580, 300, 800, 600,
	camera, factor) && std::fabs(factor - 1.5f) < 1.0e-5f,
	"screen dragging projects to a stable semantic axis scale");
    int projectedX = 0;
    int projectedY = 0;
    CHECK(manipulator->screenPosition(
	SoBRLEditManipulator::HANDLE_AXIS_A, 1.5f, 800, 600, camera,
	projectedX, projectedY) && projectedX == 580 && projectedY == 300,
	"semantic handle positions project back to canonical input pixels");

    manipulator->setActiveHandle(SoBRLEditManipulator::HANDLE_AXIS_A);
    CHECK(manipulator->activeHandle.getValue() ==
	SoBRLEditManipulator::HANDLE_AXIS_A,
	"active handle state is retained on the Obol node");
    manipulator->setVisible(FALSE);
    CHECK(manipulator->getNumChildren() == 0 &&
	manipulator->hitTest(520, 300, 800, 600, camera) ==
	    SoBRLEditManipulator::HANDLE_NONE,
	"hidden manipulators neither render nor accept input");

    const SbVec3f topologyPoints[8] = {
	SbVec3f(-1.0f, -1.0f, 0.0f), SbVec3f(1.0f, -1.0f, 0.0f),
	SbVec3f(1.0f, 1.0f, 0.0f), SbVec3f(-1.0f, 1.0f, 0.0f),
	SbVec3f(-1.0f, -1.0f, 1.0f), SbVec3f(1.0f, -1.0f, 1.0f),
	SbVec3f(1.0f, 1.0f, 1.0f), SbVec3f(-1.0f, 1.0f, 1.0f)
    };
    const int32_t topologyEdges[24] = {
	0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6,
	6, 7, 7, 4, 0, 4, 1, 5, 2, 6, 3, 7
    };
    const int32_t topologyFaces[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    const int32_t topologyFaceCounts[2] = {4, 4};
    SoBRLIndexedEditManipulator *indexed =
	new SoBRLIndexedEditManipulator;
    indexed->ref();
    indexed->setTopology(topologyPoints, 8, topologyEdges, 12,
	topologyFaces, topologyFaceCounts, 2);
    CHECK(indexed->pointCount() == 8 && indexed->edgeCount() == 12 &&
	indexed->faceCount() == 2,
	"one retained indexed manipulator stores point, edge, and face sets");
    CHECK(indexed->hitTest(SoBRLIndexedEditManipulator::DOMAIN_VERTEX,
	460, 360, 800, 600, camera) == 5,
	"indexed vertex hit testing resolves coincident projections by depth");
    CHECK(indexed->hitTest(SoBRLIndexedEditManipulator::DOMAIN_EDGE,
	400, 360, 800, 600, camera) == 4,
	"indexed edge hit testing resolves the closest retained feature");
    CHECK(indexed->hitTest(SoBRLIndexedEditManipulator::DOMAIN_FACE,
	400, 300, 800, 600, camera) == 1,
	"indexed face hit testing resolves the frontmost projected polygon");
    indexed->setSelectionDomain(SoBRLIndexedEditManipulator::DOMAIN_FACE);
    indexed->setSelectedIndex(1);
    indexed->setHoverIndex(1);
    indexed->setActiveIndex(1);
    CHECK(indexed->selectionDomain.getValue() ==
	SoBRLIndexedEditManipulator::DOMAIN_FACE &&
	indexed->selectedIndex.getValue() == 1 &&
	indexed->getNumChildren() > 2,
	"indexed selection, hover, and active state are retained and presented");
    CHECK(indexed->screenPosition(SoBRLIndexedEditManipulator::DOMAIN_FACE,
	1, 800, 600, camera, projectedX, projectedY) &&
	projectedX == 400 && projectedY == 300,
	"indexed feature centroids project to canonical input pixels");

    const SbVec3f groupedPoints[5] = {
	SbVec3f(-2.0f, 0.0f, 0.0f), SbVec3f(-1.0f, 1.0f, 0.0f),
	SbVec3f(0.0f, 0.0f, 0.0f), SbVec3f(1.0f, -1.0f, 0.0f),
	SbVec3f(2.0f, 0.0f, 0.0f)
    };
    const int32_t groupedEdges[8] = {0, 1, 1, 2, 2, 3, 3, 4};
    const int32_t groupedFeatures[4] = {0, 0, 1, 1};
    indexed->setTopology(groupedPoints, 5, groupedEdges, 4, nullptr,
	nullptr, 0, groupedFeatures);
    CHECK(indexed->edgeCount() == 2,
	"piecewise edge pairs retain two semantic feature identities");
    CHECK(indexed->hitTest(SoBRLIndexedEditManipulator::DOMAIN_EDGE,
	310, 270, 800, 600, camera) == 0 &&
	indexed->hitTest(SoBRLIndexedEditManipulator::DOMAIN_EDGE,
	490, 330, 800, 600, camera) == 1,
	"piecewise curve hit testing returns semantic group indices");
    indexed->setSelectionDomain(SoBRLIndexedEditManipulator::DOMAIN_EDGE);
    indexed->setSelectedIndex(1);
    CHECK(indexed->screenPosition(SoBRLIndexedEditManipulator::DOMAIN_EDGE,
	1, 800, 600, camera, projectedX, projectedY) &&
	projectedX == 460 && projectedY == 330,
	"grouped feature positioning averages all retained edge pieces");
    indexed->setVisible(FALSE);
    CHECK(indexed->getNumChildren() == 0 &&
	indexed->hitTest(SoBRLIndexedEditManipulator::DOMAIN_VERTEX,
	460, 360, 800, 600, camera) == -1,
	"hidden indexed manipulators neither render nor accept input");

    indexed->unref();
    manipulator->unref();
    camera->unref();
    if (failures)
	std::fprintf(stderr, "%d edit manipulator test(s) failed\n", failures);
    return failures ? 1 : 0;
}
