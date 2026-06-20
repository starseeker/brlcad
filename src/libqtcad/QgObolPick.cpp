/*                  Q G O B O L P I C K . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolPick.cpp */

#include "common.h"

#include "qtcad/QgObolPick.h"

#include "brlobol/database_source.h"
#include "brlobol/lod_service.h"
#include "brlobol/pick_detail.h"
#include "brlobol/view_controller.h"
#include "qtcad/QgView.h"

#include <Inventor/SbLine.h>
#include <Inventor/SbViewVolume.h>
#include <Inventor/SbVec2f.h>
#include <Inventor/SoPath.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/SoViewport.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/lists/SoPickedPointList.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoNode.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <float.h>
#include <string>
#include <unordered_set>

QgObolPickRecord::QgObolPickRecord(void) :
    path(),
    sourceName(),
    sourceType(),
    materialShader(),
    point(0.0f, 0.0f, 0.0f),
    materialColor(1.0f, 1.0f, 1.0f),
    distance(FLT_MAX),
    sourceId(0),
    regionId(0),
    airCode(0),
    materialId(0),
    los(0),
    primitiveKind(UNKNOWN),
    primitiveIndex(-1),
    materialColorValid(false)
{
}

static int
qg_obol_display_extent(QgView *display, int &width, int &height)
{
    width = 0;
    height = 0;
    if (!display)
	return 0;

    const double dpr = display->devicePixelRatioF();
    width = std::max(1, static_cast<int>(std::ceil(display->width() * dpr)));
    height = std::max(1, static_cast<int>(std::ceil(display->height() * dpr)));
    return width > 0 && height > 0;
}

static const SoBRLPickDetail *
qg_obol_brl_pick_detail(const SoPickedPoint *pickedPoint)
{
    if (!pickedPoint)
	return NULL;

    const SoDetail *detail = pickedPoint->getDetail();
    if (detail && detail->isOfType(SoBRLPickDetail::getClassTypeId()))
	return static_cast<const SoBRLPickDetail *>(detail);

    SoPath *path = pickedPoint->getPath();
    if (!path)
	return NULL;

    for (int i = path->getLength() - 1; i >= 0; i--) {
	SoNode *node = path->getNode(i);
	if (!node)
	    continue;
	detail = pickedPoint->getDetail(node);
	if (detail && detail->isOfType(SoBRLPickDetail::getClassTypeId()))
	    return static_cast<const SoBRLPickDetail *>(detail);
    }

    return NULL;
}

static QgObolPickRecord
qg_obol_pick_record(const SoBRLPickDetail *detail,
	const SbVec3f &point,
	float distance)
{
    QgObolPickRecord record;
    record.point = point;
    record.distance = distance;
    if (!detail)
	return record;

    record.path = detail->getPath().getString();
    record.sourceName = detail->getSourceName().getString();
    record.sourceType = detail->getSourceType().getString();
    record.materialShader = detail->getMaterialShader().getString();
    record.sourceId = detail->getSourceId();
    record.regionId = detail->getRegionId();
    record.airCode = detail->getAirCode();
    record.materialId = detail->getMaterialId();
    record.los = detail->getLos();
    record.primitiveKind = static_cast<int>(detail->getPrimitiveKind());
    record.primitiveIndex = detail->getPrimitiveIndex();
    record.materialColorValid = detail->hasMaterialColor() ? true : false;
    record.materialColor = detail->getMaterialColor();
    return record;
}

static QgObolPickRecord
qg_obol_pick_record(const SoPickedPoint *pickedPoint,
	const SoBRLPickDetail *detail,
	const SbVec3f *rayOrigin)
{
    SbVec3f point(0.0f, 0.0f, 0.0f);
    float distance = FLT_MAX;

    if (pickedPoint) {
	point = pickedPoint->getPoint();
	if (rayOrigin)
	    distance = (point - *rayOrigin).length();
    }

    return qg_obol_pick_record(detail, point, distance);
}

static QgObolPickRecord
qg_obol_pick_record(const BRLObolSourceMeshPickResult &pick)
{
    return qg_obol_pick_record(&pick.detail, pick.point, pick.distance);
}

static QgObolPickRecord
qg_obol_pick_record(const BRLObolRtPickResult &pick)
{
    return qg_obol_pick_record(&pick.detail, pick.point, pick.distance);
}

static void
qg_obol_insert_pick_record(std::vector<QgObolPickRecord> &records,
	const QgObolPickRecord &record,
	bool pickAll)
{
    if (pickAll) {
	records.push_back(record);
	return;
    }

    if (records.empty() || record.distance < records[0].distance)
	records.assign(1, record);
}

static SoBRLDatabaseSource *
qg_obol_first_database_source(BRLObolViewController *controller)
{
    if (!controller)
	return NULL;

    for (int i = 0; i < controller->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	if (source && source->getDatabase())
	    return source;
    }

    return NULL;
}

static int
qg_obol_pick_source_full_detail(BRLObolViewController *controller,
	const SbLine &line,
	bool pickAll,
	std::vector<QgObolPickRecord> &records)
{
    if (!controller || !controller->getViewport() ||
	    !controller->getViewport()->getRoot())
	return 0;

    BRLObolLodService *service = controller->getLodService();
    if (!service || !service->isRunning())
	return 0;

    SoBRLSourceMeshPickAction sourcePickAction;
    sourcePickAction.setRay(line.getPosition(), line.getDirection());
    sourcePickAction.apply(controller->getViewport()->getRoot());

    const int requestCount =
	sourcePickAction.getSourceBackedFullDetailRequestCount();
    if (requestCount <= 0)
	return 0;

    std::vector<BRLObolLodRequest> expectedRequests;
    for (int i = 0; i < requestCount; i++) {
	BRLObolLodRequest request;
	if (sourcePickAction.makeSourceBackedFullDetailLodRequest(i, request))
	    expectedRequests.push_back(request);
    }
    if (expectedRequests.empty())
	return 0;

    std::vector<BRLObolLodResult> sourceResults;
    service->drainMatchingResults(sourceResults, expectedRequests);

    int added = 0;
    if (!sourceResults.empty()) {
	BRLObolSourceMeshPickResult sourcePick;
	if (sourcePickAction.consumeSourceBackedFullDetailResults(sourcePick,
		sourceResults) > 0 && sourcePick.hit) {
	    qg_obol_insert_pick_record(records,
		    qg_obol_pick_record(sourcePick), pickAll);
	    added = 1;
	}
    } else {
	SoBRLDatabaseSource *source = qg_obol_first_database_source(controller);
	if (source)
	    (void)sourcePickAction.submitSourceBackedFullDetailRequests(
		    service, 0, source->getDatabase(), NULL,
		    controller->getMaxExactFullDetailFaceCount(),
		    controller->getMaxExactFullDetailPointCount());
    }

    return added;
}

static SbBool
qg_obol_pick_camera_line(BRLObolViewController *controller,
	int vx,
	int vy,
	SbLine &line)
{
    if (!controller || !controller->getCamera())
	return FALSE;

    const SbViewportRegion &region = controller->getViewportRegion();
    SbVec2s size = region.getViewportSizePixels();
    if (size[0] <= 0 || size[1] <= 0)
	return FALSE;

    const float nx = size[0] > 1 ?
	static_cast<float>(vx) / static_cast<float>(size[0] - 1) : 0.5f;
    const float ny = size[1] > 1 ?
	static_cast<float>(vy) / static_cast<float>(size[1] - 1) : 0.5f;
    const float aspect = static_cast<float>(size[0]) /
	static_cast<float>(size[1]);

    SbViewVolume viewVolume = controller->getCamera()->getViewVolume(aspect);
    viewVolume.projectPointToLine(SbVec2f(nx, ny), line);
    return TRUE;
}

static int
qg_obol_pick_rt_exact(BRLObolViewController *controller,
	const SbLine &line,
	bool pickAll,
	std::vector<QgObolPickRecord> &records)
{
    if (!controller)
	return 0;

    int added = 0;
    for (int i = 0; i < controller->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	if (!source || !source->getDatabase() ||
		source->path.getValue().getLength() == 0)
	    continue;

	std::vector<SbString> paths;
	paths.push_back(source->path.getValue());
	BRLObolRtPickResult rtPick;
	if (!brlobol_pick_rt_ray(rtPick, source->getDatabase(), paths,
		line.getPosition(), line.getDirection()) || !rtPick.hit)
	    continue;

	qg_obol_insert_pick_record(records, qg_obol_pick_record(rtPick),
		pickAll);
	added++;
    }

    return added;
}

int
qg_obol_pick_point(QgView *display,
	int x,
	int y,
	float radiusPixels,
	bool pickAll,
	std::vector<QgObolPickRecord> &records)
{
    records.clear();
    if (!display)
	return 0;

    BRLObolViewController *controller = display->obolViewController();
    if (!controller || !controller->getViewport() ||
	    !controller->getViewport()->getRoot())
	return 0;

    int width = 0;
    int height = 0;
    if (qg_obol_display_extent(display, width, height))
	controller->setViewportSize(static_cast<unsigned int>(width),
		static_cast<unsigned int>(height));

    const SbViewportRegion &region = controller->getViewportRegion();
    SbVec2s size = region.getViewportSizePixels();
    if (size[0] <= 0 || size[1] <= 0)
	return 0;

    int vx = x;
    int vy = size[1] - 1 - y;
    if (vx < 0)
	vx = 0;
    if (vy < 0)
	vy = 0;
    if (vx >= size[0])
	vx = size[0] - 1;
    if (vy >= size[1])
	vy = size[1] - 1;

    SoRayPickAction pickAction(region);
    pickAction.setPoint(SbVec2s(static_cast<short>(vx),
	    static_cast<short>(vy)));
    pickAction.setRadius(radiusPixels > 0.0f ? radiusPixels : 1.0f);
    pickAction.setPickAll(pickAll ? TRUE : FALSE);
    pickAction.apply(controller->getViewport()->getRoot());

    const SbLine &line = pickAction.getLine();
    const SbVec3f &rayOrigin = line.getPosition();

    if (pickAll) {
	const SoPickedPointList &pickedPoints =
	    pickAction.getPickedPointList();
	for (int i = 0; i < pickedPoints.getLength(); i++) {
	    const SoPickedPoint *pickedPoint = pickedPoints[i];
	    const SoBRLPickDetail *detail =
		qg_obol_brl_pick_detail(pickedPoint);
	    if (detail)
		qg_obol_insert_pick_record(records,
			qg_obol_pick_record(pickedPoint, detail, &rayOrigin),
			pickAll);
	}
    } else {
	const SoPickedPoint *pickedPoint = pickAction.getPickedPoint();
	const SoBRLPickDetail *detail =
	    qg_obol_brl_pick_detail(pickedPoint);
	if (detail)
	    qg_obol_insert_pick_record(records,
		    qg_obol_pick_record(pickedPoint, detail, &rayOrigin),
		    pickAll);
    }

    qg_obol_pick_source_full_detail(controller, line, pickAll, records);

    if (records.empty()) {
	SbLine rtLine = line;
	(void)qg_obol_pick_camera_line(controller, vx, vy, rtLine);
	qg_obol_pick_rt_exact(controller, rtLine, pickAll, records);
    }

    return static_cast<int>(records.size());
}

static std::string
qg_obol_pick_unique_key(const QgObolPickRecord &record)
{
    char buffer[128] = {0};
    std::snprintf(buffer, sizeof(buffer), ":%u:%d:%d:%d",
	    record.sourceId,
	    record.primitiveKind,
	    record.primitiveIndex,
	    record.materialId);
    std::string key = record.path.empty() ? record.sourceName : record.path;
    key.append(buffer);
    return key;
}

int
qg_obol_pick_rect(QgView *display,
	int x0,
	int y0,
	int x1,
	int y1,
	float radiusPixels,
	bool firstOnly,
	std::vector<QgObolPickRecord> &records)
{
    records.clear();
    if (!display)
	return 0;

    int minX = std::min(x0, x1);
    int maxX = std::max(x0, x1);
    int minY = std::min(y0, y1);
    int maxY = std::max(y0, y1);

    int displayWidth = 0;
    int displayHeight = 0;
    if (qg_obol_display_extent(display, displayWidth, displayHeight)) {
	minX = std::max(0, std::min(minX, displayWidth - 1));
	maxX = std::max(0, std::min(maxX, displayWidth - 1));
	minY = std::max(0, std::min(minY, displayHeight - 1));
	maxY = std::max(0, std::min(maxY, displayHeight - 1));
    }

    int width = std::max(1, maxX - minX);
    int height = std::max(1, maxY - minY);
    int xSteps = std::max(1, std::min(6, width / 16));
    int ySteps = std::max(1, std::min(6, height / 16));

    std::unordered_set<std::string> seen;
    for (int yi = 0; yi <= ySteps; yi++) {
	int y = minY + (height * yi) / ySteps;
	for (int xi = 0; xi <= xSteps; xi++) {
	    int x = minX + (width * xi) / xSteps;
	    std::vector<QgObolPickRecord> sampledRecords;
	    qg_obol_pick_point(display, x, y, radiusPixels, !firstOnly,
		    sampledRecords);
	    for (const QgObolPickRecord &record : sampledRecords) {
		std::string key = qg_obol_pick_unique_key(record);
		if (!seen.insert(key).second)
		    continue;
		records.push_back(record);
		if (firstOnly)
		    return static_cast<int>(records.size());
	    }
	}
    }

    return static_cast<int>(records.size());
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
