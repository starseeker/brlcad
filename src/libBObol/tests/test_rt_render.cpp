/*                 T E S T _ R T _ R E N D E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BInit.h"
#include "BObol/BDisplayEndpoint.h"
#include "BObol/BFramebuffer.h"
#include "BObol/BRtRender.h"
#include "BObol/BViewController.h"
#include "BObol/BViewportImage.h"
#include "BObol/BWindowHost.h"

#include "bu/app.h"
#include "bu/file.h"
#include "bu/malloc.h"
#include "bu/str.h"
#include "raytrace.h"
#include "wdb.h"

#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoPointLight.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSpotLight.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/SbMatrix.h>

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <vector>

#define FAIL(_message) \
    do { \
	std::fprintf(stderr, "FAIL: %s\n", (_message)); \
	return 1; \
    } while (0)

static int
make_fixture(const char *path)
{
    struct rt_wdb *wdbp = wdb_fopen_v(path, 5);
    if (!wdbp)
	return 0;
    point_t center = VINIT_ZERO;
    point_t back_center = VINIT_ZERO;
    point_t target_center = VINIT_ZERO;
    back_center[Z] = -3.0;
    target_center[Z] = -5.0;
    int ret = mk_sph(wdbp, "ball.s", center, 2.0) == 0 &&
	mk_sph(wdbp, "back.s", back_center, 2.0) == 0 &&
	mk_sph(wdbp, "target.s", target_center, 1.25) == 0 ? 0 : -1;
    if (ret == 0) {
	struct wmember region;
	BU_LIST_INIT(&region.l);
	unsigned char red[3] = {255, 0, 0};
	if (!mk_addmember("ball.s", &region.l, NULL, WMOP_UNION) ||
	    mk_lrcomb(wdbp, "ball.r", &region, 1, "plastic",
		"di=0 sp=1 sh=4", red, 1, 0, 0, 0, 0) != 0)
	    ret = -1;
    }
    wdb_close(wdbp);
    return ret == 0;
}

static size_t
count_lit(const std::vector<unsigned char> &pixels)
{
    size_t count = 0;
    for (size_t i = 0; i + 2 < pixels.size(); i += 3)
	if (pixels[i] || pixels[i + 1] || pixels[i + 2])
	    ++count;
    return count;
}

static int
has_red_surface(const std::vector<unsigned char> &pixels)
{
    for (size_t i = 0; i + 2 < pixels.size(); i += 3) {
	if (pixels[i] > pixels[i + 1] + 16u &&
	    pixels[i] > pixels[i + 2] + 16u)
	    return 1;
    }
    return 0;
}

static int
has_green_surface(const std::vector<unsigned char> &pixels)
{
    for (size_t i = 0; i + 2 < pixels.size(); i += 3) {
	if (pixels[i + 1u] > pixels[i] + 16u &&
	    pixels[i + 1u] > pixels[i + 2u] + 16u)
	    return 1;
    }
    return 0;
}

static int
has_blue_surface(const std::vector<unsigned char> &pixels)
{
    for (size_t i = 0; i + 2 < pixels.size(); i += 3) {
	if (pixels[i + 2u] > pixels[i] + 16u &&
	    pixels[i + 2u] > pixels[i + 1u] + 16u)
	    return 1;
    }
    return 0;
}

static int
mostly_green(const std::vector<unsigned char> &pixels)
{
    size_t green = 0;
    size_t count = 0;
    for (size_t i = 0; i + 2 < pixels.size(); i += 3) {
	count++;
	if (pixels[i] < 32u && pixels[i + 1u] > 224u && pixels[i + 2u] < 32u)
	    green++;
    }
    return count && green * 10u >= count * 9u;
}

static int
capture_composite(bobol_display_endpoint_t *endpoint,
	std::vector<unsigned char> &pixels)
{
    unsigned char *capture = NULL;
    size_t size = 0;
    unsigned int width = 0;
    unsigned int height = 0;
    unsigned int components = 0;
    pixels.clear();
    if (!bobol_display_endpoint_capture(endpoint, &capture, &size, &width,
	    &height, &components) || !capture || width != 64u || height != 64u ||
	components != 3u || size != 64u * 64u * 3u) {
	if (capture)
	    bu_free(capture, "mixed RT/framebuffer capture");
	return 0;
    }
    pixels.assign(capture, capture + size);
    bu_free(capture, "mixed RT/framebuffer capture");
    return 1;
}

static int
has_unlit_red_surface(const std::vector<unsigned char> &pixels)
{
    for (size_t i = 0; i + 2 < pixels.size(); i += 3) {
	if (pixels[i] == 255u && pixels[i + 1u] == 0u &&
	    pixels[i + 2u] == 0u)
	    return 1;
    }
    return 0;
}

static int
has_neutral_highlight(const std::vector<unsigned char> &pixels)
{
    for (size_t i = 0; i + 2 < pixels.size(); i += 3) {
	if (pixels[i] > 192u && pixels[i + 1u] > 192u &&
	    pixels[i + 2u] > 192u)
	    return 1;
    }
    return 0;
}

static int
has_transparent_red_green_surface(const std::vector<unsigned char> &pixels)
{
    for (size_t i = 0; i + 2 < pixels.size(); i += 3) {
	if (pixels[i] > 48u && pixels[i + 1u] > 48u &&
	    pixels[i + 2u] < pixels[i] / 2u &&
	    pixels[i + 2u] < pixels[i + 1u] / 2u)
	    return 1;
    }
    return 0;
}

static int
center_pixel_green_dominant(const std::vector<unsigned char> &pixels,
	unsigned int width, unsigned int height)
{
    if (!width || !height ||
	pixels.size() != static_cast<size_t>(width) * height * 3u)
	return 0;
    const size_t offset =
	(static_cast<size_t>(height / 2u) * width + width / 2u) * 3u;
    return pixels[offset + 1u] > pixels[offset] + 32u &&
	pixels[offset + 1u] > pixels[offset + 2u] + 32u;
}

static int
center_pixel_red_dominant(const std::vector<unsigned char> &pixels,
	unsigned int width, unsigned int height)
{
    if (!width || !height ||
	pixels.size() != static_cast<size_t>(width) * height * 3u)
	return 0;
    const size_t offset =
	(static_cast<size_t>(height / 2u) * width + width / 2u) * 3u;
    return pixels[offset] > pixels[offset + 1u] + 32u &&
	pixels[offset] > pixels[offset + 2u] + 32u;
}

static int
center_pixel_blue_dominant(const std::vector<unsigned char> &pixels,
	unsigned int width, unsigned int height)
{
    if (!width || !height ||
	pixels.size() != static_cast<size_t>(width) * height * 3u)
	return 0;
    const size_t offset =
	(static_cast<size_t>(height / 2u) * width + width / 2u) * 3u;
    return pixels[offset + 2u] > pixels[offset] + 32u &&
	pixels[offset + 2u] > pixels[offset + 1u] + 32u;
}

static SoBRLViewportImage *
viewport_image(SoGroup *root)
{
    if (!root)
	return NULL;
    for (int i = 0; i < root->getNumChildren(); ++i) {
	SoNode *node = root->getChild(i);
	if (node && node->isOfType(SoBRLViewportImage::getClassTypeId()))
	    return static_cast<SoBRLViewportImage *>(node);
    }
    return NULL;
}

static int
test_reflective_refractive_materials(struct db_i *dbip)
{
    SoSeparator *root = new SoSeparator;
    root->ref();
    SoOrthographicCamera *camera = new SoOrthographicCamera;
    camera->position.setValue(0.0f, 0.0f, 10.0f);
    camera->orientation = SbRotation::identity();
    camera->focalDistance = 10.0f;
    camera->nearDistance = 0.1f;
    camera->farDistance = 100.0f;
    camera->height = 8.0f;
    root->addChild(camera);

    int ret = 0;
    {
	BObolViewController controller(root, camera);
	controller.setRenderContextManager(bobol_headless_context_manager());
	controller.setViewportSize(64, 64);
	controller.setLightingEnabled(FALSE);
	if (controller.replaceDatabaseSource("ball.s", dbip,
		SoBRLDatabaseSource::SHADED, 1) != 1) {
	    std::fprintf(stderr, "FAIL: retained rt secondary-ray fixture did not attach\n");
	    ret = 1;
	}

	SoBRLDatabaseSource *front = ret ? NULL : controller.getDatabaseSource(0);
	BObolRtRenderer renderer;
	BObolRtRenderSettings settings;
	settings.width = 64;
	settings.height = 64;
	settings.workers = 2;
	settings.samples = 1;
	settings.backgroundBottom = SbColor(0.0f, 1.0f, 0.0f);
	settings.backgroundTop = settings.backgroundBottom;
	if (!ret && front) {
	    front->colorOverride = TRUE;
	    front->color = SbColor(1.0f, 0.0f, 0.0f);
	    front->databaseMaterialShader =
		"plastic {diffuse 1 specular 0 reflect 0 transmit 0}";
	}
	if (!ret && (!front || !renderer.synchronize(&controller))) {
	    std::fprintf(stderr, "FAIL: retained rt secondary-ray fixture did not synchronize\n");
	    ret = 1;
	}

	const uint64_t geometryRevision = renderer.getGeometryRevision();
	const uint64_t plasticRevision = renderer.getPresentationRevision();
	std::vector<unsigned char> plastic;
	if (!ret && (!renderer.render(settings, plastic, NULL) ||
		!center_pixel_red_dominant(plastic, settings.width,
		    settings.height))) {
	    std::fprintf(stderr, "FAIL: retained rt secondary-ray plastic baseline is invalid\n");
	    ret = 1;
	}

	if (!ret) {
	    front->databaseMaterialShader =
		"mirror {diffuse 0.4 specular 0.6 reflect 0.75}";
	}
	std::vector<unsigned char> mirror;
	BObolRtRenderStatus mirrorStatus;
	if (!ret && (!renderer.synchronize(&controller) ||
		renderer.getGeometryRevision() != geometryRevision ||
		renderer.getPresentationRevision() <= plasticRevision ||
		!renderer.render(settings, mirror, &mirrorStatus) ||
		mirror == plastic ||
		!center_pixel_green_dominant(mirror, settings.width,
		    settings.height) ||
		mirrorStatus.raysShot <= 64u * 64u ||
		mirrorStatus.raysShot > 64u * 64u * 2u)) {
	    std::fprintf(stderr, "FAIL: retained rt mirror did not reflect the bounded environment\n");
	    ret = 1;
	}

	if (!ret && controller.replaceDatabaseSource("target.s", dbip,
		SoBRLDatabaseSource::SHADED, 1) != 1) {
	    std::fprintf(stderr, "FAIL: retained rt refraction target did not attach\n");
	    ret = 1;
	}
	SoBRLDatabaseSource *target = ret ? NULL : controller.getDatabaseSource(1);
	if (!ret && target) {
	    target->colorOverride = TRUE;
	    target->color = SbColor(0.0f, 1.0f, 0.0f);
	    target->databaseMaterialShader =
		"plastic {diffuse 1 specular 0 reflect 0 transmit 0}";
	    front->databaseMaterialShader =
		"glass {diffuse 0.3 specular 0.7 reflect 0.1 transmit 0.8 ri 1.0}";
	    settings.backgroundBottom = SbColor(0.0f, 0.0f, 1.0f);
	    settings.backgroundTop = settings.backgroundBottom;
	}
	if (!ret && (!target || !renderer.synchronize(&controller))) {
	    std::fprintf(stderr, "FAIL: retained rt refraction fixture did not synchronize\n");
	    ret = 1;
	}
	const uint64_t glassGeometryRevision = renderer.getGeometryRevision();
	const uint64_t unitIndexRevision = renderer.getPresentationRevision();
	std::vector<unsigned char> unitIndexGlass;
	if (!ret && (!renderer.render(settings, unitIndexGlass, NULL) ||
		!center_pixel_green_dominant(unitIndexGlass, settings.width,
		    settings.height))) {
	    std::fprintf(stderr, "FAIL: retained rt glass did not transmit geometry behind the solid\n");
	    ret = 1;
	}

	if (!ret)
	    front->databaseMaterialShader =
		"glass {diffuse 0.3 specular 0.7 reflect 0.1 transmit 0.8 ri 1.65}";
	std::vector<unsigned char> refracted;
	BObolRtRenderStatus refractedStatus;
	if (!ret && (!renderer.synchronize(&controller) ||
		renderer.getGeometryRevision() != glassGeometryRevision ||
		renderer.getPresentationRevision() <= unitIndexRevision ||
		!renderer.render(settings, refracted, &refractedStatus) ||
		refracted == unitIndexGlass ||
		!center_pixel_green_dominant(refracted, settings.width,
		    settings.height) ||
		refractedStatus.raysShot <= 64u * 64u * 2u ||
		refractedStatus.raysShot > 64u * 64u * 30u)) {
	    std::fprintf(stderr, "FAIL: retained rt glass did not apply bounded refraction\n");
	    ret = 1;
	}

	settings.workers = 1;
	std::vector<unsigned char> serialRefracted;
	if (!ret && (!renderer.render(settings, serialRefracted, NULL) ||
		serialRefracted != refracted)) {
	    std::fprintf(stderr, "FAIL: retained rt secondary rays depend on worker layout\n");
	    ret = 1;
	}

	if (!ret)
	    front->databaseMaterialShader =
		"glass {diffuse 1 specular 0 reflect 0 transmit 0 ri 1.65}";
	std::vector<unsigned char> opaqueGlass;
	if (!ret && (!renderer.synchronize(&controller) ||
		renderer.getGeometryRevision() != glassGeometryRevision ||
		!renderer.render(settings, opaqueGlass, NULL) ||
		!center_pixel_red_dominant(opaqueGlass, settings.width,
		    settings.height) || opaqueGlass == refracted)) {
	    std::fprintf(stderr, "FAIL: retained rt glass parameter update did not remain presentation-only\n");
	    ret = 1;
	}
    }

    root->unref();
    return ret;
}

static int
test_authored_retained_lights(struct db_i *dbip)
{
    SoSeparator *root = new SoSeparator;
    root->ref();
    SoOrthographicCamera *camera = new SoOrthographicCamera;
    camera->position.setValue(0.0f, 0.0f, 10.0f);
    camera->orientation = SbRotation::identity();
    camera->focalDistance = 10.0f;
    camera->nearDistance = 0.1f;
    camera->farDistance = 100.0f;
    camera->height = 8.0f;
    root->addChild(camera);

    SoDirectionalLight *directional = new SoDirectionalLight;
    directional->on = FALSE;
    directional->direction = SbVec3f(0.0f, 0.0f, -1.0f);
    directional->color = SbColor(0.0f, 0.0f, 1.0f);
    root->addChild(directional);

    SoSeparator *pointGroup = new SoSeparator;
    SoTranslation *pointTransform = new SoTranslation;
    pointTransform->translation = SbVec3f(0.0f, 0.0f, 10.0f);
    SoPointLight *point = new SoPointLight;
    point->on = FALSE;
    point->location = SbVec3f(0.0f, 0.0f, 0.0f);
    point->color = SbColor(0.0f, 1.0f, 0.0f);
    pointGroup->addChild(pointTransform);
    pointGroup->addChild(point);
    root->addChild(pointGroup);

    SoSeparator *spotGroup = new SoSeparator;
    SoTranslation *spotTransform = new SoTranslation;
    spotTransform->translation = SbVec3f(0.0f, 0.0f, 10.0f);
    SoSpotLight *spot = new SoSpotLight;
    spot->on = FALSE;
    spot->location = SbVec3f(0.0f, 0.0f, 0.0f);
    spot->direction = SbVec3f(0.0f, 0.0f, -1.0f);
    spot->cutOffAngle = 0.4f;
    spot->dropOffRate = 0.25f;
    spot->color = SbColor(1.0f, 0.0f, 0.0f);
    spotGroup->addChild(spotTransform);
    spotGroup->addChild(spot);
    root->addChild(spotGroup);

    int ret = 0;
    {
	BObolViewController controller(root, camera);
	controller.setViewportSize(64, 64);
	controller.setHeadlightIntensity(0.0f);
	if (controller.replaceDatabaseSource("ball.s", dbip,
		SoBRLDatabaseSource::SHADED, 1) != 1) {
	    std::fprintf(stderr, "FAIL: retained authored-light fixture did not attach\n");
	    ret = 1;
	}
	SoBRLDatabaseSource *source = ret ? NULL : controller.getDatabaseSource(0);
	if (!ret && source) {
	    source->colorOverride = TRUE;
	    source->color = SbColor(1.0f, 1.0f, 1.0f);
	    source->databaseMaterialShader =
		"plastic {diffuse 1 specular 0 reflect 0 transmit 0}";
	}

	BObolRtRenderer renderer;
	BObolRtRenderSettings settings;
	settings.width = 64;
	settings.height = 64;
	settings.workers = 2;
	settings.samples = 1;
	if (!ret && (!source || !renderer.synchronize(&controller))) {
	    std::fprintf(stderr, "FAIL: retained authored-light fixture did not synchronize\n");
	    ret = 1;
	}
	const uint64_t geometryRevision = renderer.getGeometryRevision();
	uint64_t presentationRevision = renderer.getPresentationRevision();
	std::vector<unsigned char> unlit;
	if (!ret && !renderer.render(settings, unlit, NULL)) {
	    std::fprintf(stderr, "FAIL: retained authored-light baseline did not render\n");
	    ret = 1;
	}

	directional->on = TRUE;
	std::vector<unsigned char> directionalImage;
	if (!ret && (!renderer.synchronize(&controller) ||
		renderer.getGeometryRevision() != geometryRevision ||
		renderer.getPresentationRevision() <= presentationRevision ||
		!renderer.render(settings, directionalImage, NULL) ||
		directionalImage == unlit ||
		!center_pixel_blue_dominant(directionalImage, settings.width,
		    settings.height))) {
	    std::fprintf(stderr, "FAIL: retained rt did not snapshot directional scene light\n");
	    ret = 1;
	}
	presentationRevision = renderer.getPresentationRevision();
	directional->direction = SbVec3f(0.0f, 0.0f, 1.0f);
	std::vector<unsigned char> reversedDirectional;
	if (!ret && (!renderer.synchronize(&controller) ||
		renderer.getGeometryRevision() != geometryRevision ||
		renderer.getPresentationRevision() <= presentationRevision ||
		!renderer.render(settings, reversedDirectional, NULL) ||
		reversedDirectional == directionalImage ||
		center_pixel_blue_dominant(reversedDirectional, settings.width,
		    settings.height))) {
	    std::fprintf(stderr, "FAIL: retained rt directional light edit was not presentation-only\n");
	    ret = 1;
	}

	presentationRevision = renderer.getPresentationRevision();
	directional->on = FALSE;
	point->on = TRUE;
	std::vector<unsigned char> pointImage;
	if (!ret && (!renderer.synchronize(&controller) ||
		renderer.getGeometryRevision() != geometryRevision ||
		renderer.getPresentationRevision() <= presentationRevision ||
		!renderer.render(settings, pointImage, NULL) ||
		!center_pixel_green_dominant(pointImage, settings.width,
		    settings.height))) {
	    std::fprintf(stderr, "FAIL: retained rt did not transform point scene light\n");
	    ret = 1;
	}
	settings.workers = 1;
	std::vector<unsigned char> serialPoint;
	if (!ret && (!renderer.render(settings, serialPoint, NULL) ||
		serialPoint != pointImage)) {
	    std::fprintf(stderr, "FAIL: retained authored light depends on worker layout\n");
	    ret = 1;
	}

	presentationRevision = renderer.getPresentationRevision();
	point->on = FALSE;
	spot->on = TRUE;
	std::vector<unsigned char> spotImage;
	if (!ret && (!renderer.synchronize(&controller) ||
		renderer.getGeometryRevision() != geometryRevision ||
		renderer.getPresentationRevision() <= presentationRevision ||
		!renderer.render(settings, spotImage, NULL) ||
		!center_pixel_red_dominant(spotImage, settings.width,
		    settings.height))) {
	    std::fprintf(stderr, "FAIL: retained rt did not transform spot scene light\n");
	    ret = 1;
	}
	presentationRevision = renderer.getPresentationRevision();
	spot->direction = SbVec3f(0.0f, 0.0f, 1.0f);
	std::vector<unsigned char> reversedSpot;
	if (!ret && (!renderer.synchronize(&controller) ||
		renderer.getGeometryRevision() != geometryRevision ||
		renderer.getPresentationRevision() <= presentationRevision ||
		!renderer.render(settings, reversedSpot, NULL) ||
		reversedSpot == spotImage ||
		center_pixel_red_dominant(reversedSpot, settings.width,
		    settings.height))) {
	    std::fprintf(stderr, "FAIL: retained rt spot cone edit was not applied\n");
	    ret = 1;
	}
    }

    root->unref();
    return ret;
}

/* Hosted applications such as MGED keep their authoritative database sources
 * in a shared render root rather than in each endpoint controller's local
 * scene.  Retained RT must discover that same generic scene contract. */
static int
test_shared_render_root(struct db_i *dbip)
{
    SoSeparator *localRoot = new SoSeparator;
    localRoot->ref();
    SoOrthographicCamera *camera = new SoOrthographicCamera;
    camera->position.setValue(0.0f, 0.0f, 10.0f);
    camera->orientation = SbRotation::identity();
    camera->focalDistance = 10.0f;
    camera->nearDistance = 0.1f;
    camera->farDistance = 100.0f;
    camera->height = 8.0f;
    localRoot->addChild(camera);

    SoSeparator *sharedRoot = new SoSeparator;
    sharedRoot->ref();
    SoBRLDatabaseSource *sharedSource = new SoBRLDatabaseSource;
    sharedSource->configureDatabaseSourceInstance("shared-retained-rt",
	"ball.r", dbip, SoBRLDatabaseSource::SHADED, 1);
    sharedRoot->addChild(sharedSource);

    int ret = 0;
    {
	BObolViewController controller(localRoot, camera);
	controller.setViewportSize(64, 64);
	controller.setBackgroundColors(SbColor(0.0f, 0.0f, 0.0f),
	    SbColor(0.0f, 0.0f, 0.0f));
	controller.setRenderSceneRoot(sharedRoot);

	const std::vector<SoBRLDatabaseSource *> sources =
	    controller.getRenderDatabaseSources();
	if (controller.getDatabaseSourceCount() != 0 || sources.size() != 1u ||
	    sources[0] != sharedSource) {
	    std::fprintf(stderr, "FAIL: shared render-root source discovery disagrees with rendered scene\n");
	    ret = 1;
	} else {
	    BObolRtRenderer renderer;
	    BObolRtRenderSettings settings;
	    settings.width = 64;
	    settings.height = 64;
	    settings.workers = 2;
	    settings.samples = 1;
	    std::vector<unsigned char> pixels;
	    if (!renderer.synchronize(&controller) ||
		renderer.getPreparedSourceCount() != 1u ||
		!renderer.getGeometryRevision() ||
		!renderer.render(settings, pixels, NULL) || count_lit(pixels) == 0) {
		std::fprintf(stderr, "FAIL: retained rt ignored authoritative shared render root\n");
		ret = 1;
	    }
	}
    }

    sharedRoot->unref();
    localRoot->unref();
    return ret;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    (void)argc;
    bobol_init(NULL);

    char dbpath[MAXPATHLEN] = {0};
    FILE *file = bu_temp_file(dbpath, MAXPATHLEN);
    if (!file)
	FAIL("failed to create retained rt fixture path");
    std::fclose(file);
    if (!make_fixture(dbpath)) {
	bu_file_delete(dbpath);
	FAIL("failed to create retained rt fixture database");
    }

    struct db_i *dbip = db_open(dbpath, DB_OPEN_READONLY);
    if (!dbip || db_dirbuild(dbip) < 0) {
	if (dbip)
	    db_close(dbip);
	bu_file_delete(dbpath);
	FAIL("failed to open retained rt fixture database");
    }

    int ret = 0;
    SoSeparator *root = new SoSeparator;
    root->ref();
    SoOrthographicCamera *camera = new SoOrthographicCamera;
    camera->position.setValue(0.0f, 0.0f, 10.0f);
    camera->orientation = SbRotation::identity();
    camera->focalDistance = 10.0f;
    camera->nearDistance = 0.1f;
    camera->farDistance = 100.0f;
    camera->height = 8.0f;
    root->addChild(camera);

    {
	BObolViewController controller(root, camera);
	controller.setRenderContextManager(bobol_headless_context_manager());
	controller.setViewportSize(64, 64);
	controller.setBackgroundColors(SbColor(0.0f, 0.0f, 0.0f),
	    SbColor(0.0f, 0.0f, 0.0f));
	std::vector<BObolSceneLightRealization> cameraLights;
	controller.getCameraLights(cameraLights);
	if (controller.getLightingProfile() !=
		BObolViewController::LIGHTING_STUDIO ||
	    cameraLights.size() != 3u ||
	    std::fabs(controller.getLightingAmbientIntensity() - 0.18f) >
		1.0e-6f) {
	    std::fprintf(stderr, "FAIL: retained rt did not inherit studio lighting defaults\n");
	    ret = 1;
	}
	/* Material-metadata assertions below use the historical straight-on
	 * highlight as a controlled fixture.  Studio rendering has a separate
	 * image-transition assertion after the material checks. */
	controller.setLightingProfile(BObolViewController::LIGHTING_MGED);
	controller.getCameraLights(cameraLights);
	if (!ret && (cameraLights.size() != 1u ||
	    std::fabs(controller.getLightingAmbientIntensity() - 0.30f) >
		1.0e-6f)) {
	    std::fprintf(stderr, "FAIL: retained rt MGED lighting profile is incomplete\n");
	    ret = 1;
	}
	if (controller.replaceDatabaseSource("ball.s", dbip,
		SoBRLDatabaseSource::SHADED, 1) != 1) {
	    std::fprintf(stderr, "FAIL: failed to attach retained rt source\n");
	    ret = 1;
	} else {
	    BObolRtRenderer renderer;
	    if (!renderer.synchronize(&controller) ||
		renderer.getPreparedSourceCount() != 1) {
		std::fprintf(stderr, "FAIL: retained rt renderer did not prepare source\n");
		ret = 1;
	    } else {
		uint64_t geometryRevision = renderer.getGeometryRevision();
		const uint64_t presentationRevision = renderer.getPresentationRevision();
		BObolRtRenderSettings settings;
		settings.width = 64;
		settings.height = 64;
		settings.workers = 2;
		settings.samples = 1;
		std::vector<unsigned char> first;
		BObolRtRenderStatus firstStatus;
		if (!renderer.render(settings, first, &firstStatus) ||
		    !firstStatus.complete || firstStatus.cancelled ||
		    firstStatus.preparedSources != 1 ||
		    firstStatus.raysShot != 64u * 64u ||
		    count_lit(first) == 0) {
		    std::fprintf(stderr, "FAIL: retained rt renderer did not produce image\n");
		    ret = 1;
		}

		/* The worker layout must not affect a fixed final frame. */
		settings.workers = 1;
		std::vector<unsigned char> serial;
		BObolRtRenderStatus serialStatus;
		if (!ret && (!renderer.render(settings, serial, &serialStatus) ||
		    serial != first || !serialStatus.complete)) {
		    std::fprintf(stderr, "FAIL: retained rt final frame is not deterministic\n");
		    ret = 1;
		}
		std::vector<unsigned char> tiled;
		for (unsigned int row = 0; !ret && row < settings.height; row += 7u) {
		    BObolRtRenderStatus tileStatus;
		    if (!renderer.renderRows(settings, tiled, row,
			std::min(7u, settings.height - row), &tileStatus) ||
			!tileStatus.complete) {
			std::fprintf(stderr, "FAIL: retained rt tiled pass did not complete\n");
			ret = 1;
		    }
		}
		if (!ret && tiled != serial) {
		    std::fprintf(stderr, "FAIL: retained rt tiled pass differs from final frame\n");
		    ret = 1;
		}
		BObolRtRenderPlanes planes;
		std::vector<unsigned char> planeRgb;
		if (!ret && (!renderer.renderWithPlanes(settings, planeRgb, planes,
		    NULL) || planeRgb != serial ||
		    planes.depth.size() != 64u * 64u ||
		    planes.sourceIdentity.size() != 64u * 64u)) {
		    std::fprintf(stderr, "FAIL: retained rt output planes are incomplete\n");
		    ret = 1;
		}
		uint32_t hitIdentity = 0u;
		float nearestDepth = 1.0f;
		for (size_t pixel = 0; !ret && pixel < planes.depth.size(); ++pixel) {
		    if (planes.sourceIdentity[pixel] != 0u) {
			hitIdentity = planes.sourceIdentity[pixel];
			nearestDepth = std::min(nearestDepth, planes.depth[pixel]);
		    }
		}
		BObolRtSourceIdentity sourceIdentity;
		if (!ret && (!hitIdentity || nearestDepth >= 1.0f ||
		    !renderer.getSourceIdentity(hitIdentity, sourceIdentity) ||
		    sourceIdentity.database != dbip ||
		    bu_strcmp(sourceIdentity.path.getString(), "ball.s") != 0)) {
		    std::fprintf(stderr, "FAIL: retained rt output identity did not resolve source\n");
		    ret = 1;
		}

		/* Camera state is presentation state.  Re-synchronizing it must retain
		 * the prepared source acceleration structure. */
		camera->position.setValue(2.0f, 0.0f, 10.0f);
		if (!ret && (!renderer.synchronize(&controller) ||
		    renderer.getGeometryRevision() != geometryRevision ||
		    renderer.getPresentationRevision() <= presentationRevision)) {
		    std::fprintf(stderr, "FAIL: camera update did not retain RT geometry\n");
		    ret = 1;
		}

		/* Appearance updates must reuse the retained acceleration structure. */
		SoBRLDatabaseSource *source = controller.getDatabaseSource(0);
		source->colorOverride = TRUE;
		source->color = SbColor(1.0f, 0.0f, 0.0f);
		source->databaseMaterialShader = "plastic {di 0 sp 1 sh 4}";
		if (!ret && (!renderer.synchronize(&controller) ||
		    renderer.getGeometryRevision() != geometryRevision ||
		    renderer.getPresentationRevision() == presentationRevision)) {
		    std::fprintf(stderr, "FAIL: appearance update rebuilt retained rt geometry\n");
		    ret = 1;
		}
		std::vector<unsigned char> specularRed;
		if (!ret && (!renderer.render(settings, specularRed, NULL) ||
		    !has_neutral_highlight(specularRed))) {
		    std::fprintf(stderr, "FAIL: retained rt did not apply plastic specular metadata\n");
		    ret = 1;
		}
		const uint64_t specularRevision = renderer.getPresentationRevision();
		source->databaseMaterialShader =
		    "plastic {diffuse 1 specular 0 shine 4}";
		std::vector<unsigned char> red;
		if (!ret && (!renderer.synchronize(&controller) ||
		    renderer.getGeometryRevision() != geometryRevision ||
		    renderer.getPresentationRevision() <= specularRevision ||
		    !renderer.render(settings, red, NULL) || red == specularRed ||
		    !has_red_surface(red) || has_neutral_highlight(red))) {
		    std::fprintf(stderr, "FAIL: retained rt material update did not reuse geometry and change lighting\n");
		    ret = 1;
		}

		/* Region shader metadata wins over the source fallback without
		 * introducing a renderer-specific scene representation. */
		if (!ret && (controller.removeDatabaseSource("ball.s") != 1 ||
		    controller.replaceDatabaseSource("ball.r", dbip,
			SoBRLDatabaseSource::SHADED, 1) != 1)) {
		    std::fprintf(stderr, "FAIL: retained rt primitive material fixture did not replace source\n");
		    ret = 1;
		}
		if (!ret) {
		    source = controller.getDatabaseSource(0);
		    source->colorOverride = TRUE;
		    source->color = SbColor(1.0f, 0.0f, 0.0f);
		    source->databaseMaterialShader =
			"plastic {diffuse 1 specular 0 shine 4}";
		    if (!renderer.synchronize(&controller)) {
			std::fprintf(stderr, "FAIL: retained rt primitive material fixture did not synchronize\n");
			ret = 1;
		    } else {
			geometryRevision = renderer.getGeometryRevision();
		    }
		}
		std::vector<unsigned char> primitiveSpecular;
		if (!ret && (!renderer.render(settings, primitiveSpecular, NULL) ||
		    !has_neutral_highlight(primitiveSpecular))) {
		    std::fprintf(stderr, "FAIL: retained rt did not prioritize primitive shader metadata\n");
		    ret = 1;
		}
		const uint64_t whiteLightRevision = renderer.getPresentationRevision();
		controller.setHeadlightColor(SbColor(0.0f, 0.0f, 1.0f));
		std::vector<unsigned char> blueLight;
		if (!ret && (!renderer.synchronize(&controller) ||
		    renderer.getGeometryRevision() != geometryRevision ||
		    renderer.getPresentationRevision() <= whiteLightRevision ||
		    !renderer.render(settings, blueLight, NULL) ||
		    !has_blue_surface(blueLight) || blueLight == primitiveSpecular)) {
		    std::fprintf(stderr, "FAIL: retained rt headlight update did not reuse geometry and change lighting\n");
		    ret = 1;
		}
		controller.setHeadlightColor(SbColor(1.0f, 1.0f, 1.0f));
		if (!ret && !renderer.synchronize(&controller)) {
		    std::fprintf(stderr, "FAIL: retained rt did not restore white headlight\n");
		    ret = 1;
		}
		std::vector<unsigned char> mgedLighting;
		if (!ret && !renderer.render(settings, mgedLighting, NULL)) {
		    std::fprintf(stderr, "FAIL: retained rt did not render MGED lighting profile\n");
		    ret = 1;
		}
		const uint64_t mgedLightingRevision =
		    renderer.getPresentationRevision();
		controller.setLightingProfile(BObolViewController::LIGHTING_STUDIO);
		std::vector<unsigned char> studioLighting;
		if (!ret && (!renderer.synchronize(&controller) ||
		    renderer.getGeometryRevision() != geometryRevision ||
		    renderer.getPresentationRevision() <= mgedLightingRevision ||
		    !renderer.render(settings, studioLighting, NULL) ||
		    studioLighting == mgedLighting)) {
		    std::fprintf(stderr, "FAIL: retained rt did not apply studio lighting as presentation policy\n");
		    ret = 1;
		}
		const uint64_t litRevision = renderer.getPresentationRevision();
		controller.setLightingEnabled(FALSE);
		if (!ret && (!renderer.synchronize(&controller) ||
		    renderer.getGeometryRevision() != geometryRevision ||
		    renderer.getPresentationRevision() <= litRevision)) {
		    std::fprintf(stderr, "FAIL: lighting update rebuilt retained rt geometry\n");
		    ret = 1;
		}
		std::vector<unsigned char> unlit;
		if (!ret && (!renderer.render(settings, unlit, NULL) ||
		    !has_unlit_red_surface(unlit))) {
		    std::fprintf(stderr, "FAIL: retained rt did not apply lighting policy\n");
		    ret = 1;
		}

		/* Semi-transparent sources must retain the opaque geometry behind
		 * them instead of blending only against the background. */
		if (!ret && controller.replaceDatabaseSource("back.s", dbip,
		    SoBRLDatabaseSource::SHADED, 1) != 1) {
		    std::fprintf(stderr, "FAIL: retained rt transparency fixture did not attach\n");
		    ret = 1;
		}
		SoBRLDatabaseSource *back_source = controller.getDatabaseSource(1);
		if (!ret && !back_source) {
		    std::fprintf(stderr, "FAIL: retained rt transparency fixture source is missing\n");
		    ret = 1;
		}
		if (!ret) {
		    source->transparency = 0.5f;
		    back_source->colorOverride = TRUE;
		    back_source->color = SbColor(0.0f, 1.0f, 0.0f);
		    if (!renderer.synchronize(&controller)) {
			std::fprintf(stderr, "FAIL: retained rt transparency fixture did not synchronize\n");
			ret = 1;
		    }
		}
		std::vector<unsigned char> transparent;
		if (!ret && (!renderer.render(settings, transparent, NULL) ||
		    !has_transparent_red_green_surface(transparent))) {
		    std::fprintf(stderr, "FAIL: retained rt did not composite geometry behind transparency\n");
		    ret = 1;
		}
		controller.setTransparencyEnabled(FALSE);
		std::vector<unsigned char> transparency_disabled;
		if (!ret && (!renderer.synchronize(&controller) ||
		    !renderer.render(settings, transparency_disabled, NULL) ||
		    has_transparent_red_green_surface(transparency_disabled) ||
		    !has_unlit_red_surface(transparency_disabled))) {
		    std::fprintf(stderr, "FAIL: retained rt transparency-disable policy did not keep front source opaque\n");
		    ret = 1;
		}
		if (back_source)
		    (void)controller.removeDatabaseSource("back.s");
		source->transparency = 0.0f;
		controller.setLightingEnabled(TRUE);
		controller.setTransparencyEnabled(TRUE);
		if (!ret && !renderer.synchronize(&controller)) {
		    std::fprintf(stderr, "FAIL: retained rt did not restore opaque fixture state\n");
		    ret = 1;
		}
		geometryRevision = renderer.getGeometryRevision();

		/* The endpoint owns presentation: a worker-produced librt frame must
		 * arrive through an Obol image node and be capturable deterministically. */
		bobol_display_endpoint_t *endpoint =
		    bobol_display_endpoint_create(&controller, 0);
		struct bv_display_property_value rtPolicy =
		    BV_DISPLAY_PROPERTY_VALUE_INIT;
		rtPolicy.type = BV_DISPLAY_PROPERTY_UINT;
		rtPolicy.uint_value = 1;
		if (!ret && (!endpoint || bobol_display_endpoint_property_set(
		    endpoint, "render.rt.frame_budget_ms", &rtPolicy) !=
		    BV_DISPLAY_PROPERTY_OK)) {
		    std::fprintf(stderr, "FAIL: retained rt frame-budget policy was rejected\n");
		    ret = 1;
		}
		rtPolicy = BV_DISPLAY_PROPERTY_VALUE_INIT;
		rtPolicy.type = BV_DISPLAY_PROPERTY_ENUM;
		rtPolicy.string_value = "interactive";
	if (!ret && bobol_display_endpoint_property_set(endpoint,
	    "render.rt.quality", &rtPolicy) != BV_DISPLAY_PROPERTY_OK) {
		    std::fprintf(stderr, "FAIL: retained rt quality policy was rejected\n");
	    ret = 1;
	}
	rtPolicy = BV_DISPLAY_PROPERTY_VALUE_INIT;
	rtPolicy.type = BV_DISPLAY_PROPERTY_ENUM;
	rtPolicy.string_value = "overlay";
	if (!ret && bobol_display_endpoint_property_set(endpoint,
	    "composition.rt.layer", &rtPolicy) != BV_DISPLAY_PROPERTY_OK) {
	    std::fprintf(stderr, "FAIL: retained RT image-layer policy was rejected\n");
	    ret = 1;
	}
		rtPolicy = BV_DISPLAY_PROPERTY_VALUE_INIT;
		if (!ret && (bobol_display_endpoint_property_get(endpoint,
		    "render.rt.frame_budget_ms", &rtPolicy) !=
		    BV_DISPLAY_PROPERTY_OK || rtPolicy.uint_value != 1u)) {
		    std::fprintf(stderr, "FAIL: retained rt frame-budget policy did not persist\n");
		    ret = 1;
		}
		rtPolicy = BV_DISPLAY_PROPERTY_VALUE_INIT;
	if (!ret && (bobol_display_endpoint_property_get(endpoint,
	    "render.rt.quality", &rtPolicy) != BV_DISPLAY_PROPERTY_OK ||
		    !rtPolicy.string_value ||
		    bu_strcmp(rtPolicy.string_value, "interactive") != 0)) {
		    std::fprintf(stderr, "FAIL: retained rt quality policy did not persist\n");
	    ret = 1;
	}
	rtPolicy = BV_DISPLAY_PROPERTY_VALUE_INIT;
	if (!ret && (bobol_display_endpoint_property_get(endpoint,
	    "composition.rt.layer", &rtPolicy) != BV_DISPLAY_PROPERTY_OK ||
	    !rtPolicy.string_value || bu_strcmp(rtPolicy.string_value, "overlay") != 0)) {
	    std::fprintf(stderr, "FAIL: retained RT image-layer policy did not persist\n");
	    ret = 1;
	}
		if (!ret && (!endpoint || !bobol_display_endpoint_render_engine_set(
		    endpoint, BOBOL_RENDER_ENGINE_RT))) {
		    std::fprintf(stderr, "FAIL: retained rt endpoint did not initialize\n");
		    ret = 1;
		}
		unsigned char *capture = NULL;
		size_t captureSize = 0;
		unsigned int captureWidth = 0;
		unsigned int captureHeight = 0;
		unsigned int captureComponents = 0;
		if (!ret && (!bobol_display_endpoint_capture(endpoint, &capture,
		    &captureSize, &captureWidth, &captureHeight, &captureComponents) ||
		    !capture || captureWidth != 64u || captureHeight != 64u ||
		    captureComponents != 3u || captureSize != 64u * 64u * 3u)) {
		    std::fprintf(stderr, "FAIL: retained rt endpoint did not capture frame\n");
		    ret = 1;
		}
		std::vector<unsigned char> endpointPixels;
		if (capture && captureSize)
		    endpointPixels.assign(capture, capture + captureSize);
		if (capture)
		    bu_free(capture, "retained rt endpoint capture");
	if (!ret && count_lit(endpointPixels) == 0) {
	    unsigned int maximumRed = 0;
	    unsigned int maximumRedDelta = 0;
	    for (size_t i = 0; i + 2 < endpointPixels.size(); i += 3) {
		maximumRed = std::max(maximumRed,
		    static_cast<unsigned int>(endpointPixels[i]));
		maximumRedDelta = std::max(maximumRedDelta,
		    static_cast<unsigned int>(endpointPixels[i]) -
		    std::min(static_cast<unsigned int>(endpointPixels[i]),
			std::max(static_cast<unsigned int>(endpointPixels[i + 1]),
			    static_cast<unsigned int>(endpointPixels[i + 2]))));
	    }
	    std::fprintf(stderr,
		"FAIL: retained rt endpoint did not present ray image (max red %u, delta %u, lit %zu)\n",
		maximumRed, maximumRedDelta, count_lit(endpointPixels));
	    ret = 1;
	}
	struct bv_display_property_value rtGeometryRevision =
	    BV_DISPLAY_PROPERTY_VALUE_INIT;
	struct bv_display_property_value rtPresentationRevision =
	    BV_DISPLAY_PROPERTY_VALUE_INIT;
	if (!ret && (bobol_display_endpoint_property_get(endpoint,
	    "render.rt.geometry_revision", &rtGeometryRevision) !=
	    BV_DISPLAY_PROPERTY_OK ||
	    bobol_display_endpoint_property_get(endpoint,
	    "render.rt.presentation_revision", &rtPresentationRevision) !=
	    BV_DISPLAY_PROPERTY_OK || !rtGeometryRevision.uint_value ||
	    !rtPresentationRevision.uint_value)) {
	    std::fprintf(stderr, "FAIL: retained rt endpoint did not expose renderer revisions\n");
	    ret = 1;
	}

	/* A view change arriving during an expensive frame must cancel that
	 * generation, clear its presentation, and publish only the replacement
	 * camera.  This is the controller-owned path used by QGED and MGED. */
	rtPolicy = BV_DISPLAY_PROPERTY_VALUE_INIT;
	rtPolicy.type = BV_DISPLAY_PROPERTY_UINT;
	rtPolicy.uint_value = 64;
	if (!ret && bobol_display_endpoint_property_set(endpoint,
	    "render.rt.samples", &rtPolicy) != BV_DISPLAY_PROPERTY_OK) {
	    std::fprintf(stderr, "FAIL: retained rt cancellation fixture did not start expensive frame\n");
	    ret = 1;
	}
	camera->position.setValue(20.0f, 0.0f, 10.0f);
	controller.requestLodCapacityRender("camera");
	unsigned char *cameraOwnerImage = NULL;
	if (!ret && controller.renderToImage(&cameraOwnerImage, 0, 0, NULL,
	    NULL, NULL) != BRLCAD_OK) {
	    std::fprintf(stderr, "FAIL: retained rt camera cancellation redraw failed\n");
	    ret = 1;
	}
	if (cameraOwnerImage)
	    bu_free(cameraOwnerImage, "retained rt camera cancellation redraw");
	std::vector<unsigned char> movedCameraPixels;
	if (!ret && (!capture_composite(endpoint, movedCameraPixels) ||
	    count_lit(movedCameraPixels) != 0 ||
	    has_red_surface(movedCameraPixels))) {
	    std::fprintf(stderr, "FAIL: retained rt published stale pixels after camera cancellation\n");
	    ret = 1;
	}
	struct bv_display_property_value movedGeometryRevision =
	    BV_DISPLAY_PROPERTY_VALUE_INIT;
	struct bv_display_property_value movedPresentationRevision =
	    BV_DISPLAY_PROPERTY_VALUE_INIT;
	if (!ret && (bobol_display_endpoint_property_get(endpoint,
	    "render.rt.geometry_revision", &movedGeometryRevision) !=
	    BV_DISPLAY_PROPERTY_OK ||
	    bobol_display_endpoint_property_get(endpoint,
	    "render.rt.presentation_revision", &movedPresentationRevision) !=
	    BV_DISPLAY_PROPERTY_OK ||
	    movedGeometryRevision.uint_value != rtGeometryRevision.uint_value ||
	    movedPresentationRevision.uint_value <=
		rtPresentationRevision.uint_value)) {
	    std::fprintf(stderr, "FAIL: retained rt camera cancellation rebuilt geometry or missed presentation revision\n");
	    ret = 1;
	}
	rtPolicy.uint_value = 1;
	if (!ret && bobol_display_endpoint_property_set(endpoint,
	    "render.rt.samples", &rtPolicy) != BV_DISPLAY_PROPERTY_OK) {
	    std::fprintf(stderr, "FAIL: retained rt cancellation fixture did not restore sample policy\n");
	    ret = 1;
	}
	camera->position.setValue(2.0f, 0.0f, 10.0f);
	controller.requestLodCapacityRender("camera");
	cameraOwnerImage = NULL;
	if (!ret && controller.renderToImage(&cameraOwnerImage, 0, 0, NULL,
	    NULL, NULL) != BRLCAD_OK) {
	    std::fprintf(stderr, "FAIL: retained rt camera restore redraw failed\n");
	    ret = 1;
	}
	if (cameraOwnerImage)
	    bu_free(cameraOwnerImage, "retained rt camera restore redraw");
	if (!ret && (!capture_composite(endpoint, endpointPixels) ||
	    count_lit(endpointPixels) == 0)) {
	    std::fprintf(stderr, "FAIL: retained rt did not publish replacement camera frame\n");
	    ret = 1;
	}
	SoBRLViewportImage *rtViewport = viewport_image(
	    controller.getFramebufferOverlayRoot());
	if (!ret && (!rtViewport ||
	    rtViewport->layer.getValue() != SoBRLViewportImage::OVERLAY ||
	    viewport_image(controller.getFramebufferInterlayRoot()) != NULL)) {
	    std::fprintf(stderr, "FAIL: retained RT overlay policy did not reparent image\n");
	    ret = 1;
	}
	rtPolicy = BV_DISPLAY_PROPERTY_VALUE_INIT;
	rtPolicy.type = BV_DISPLAY_PROPERTY_ENUM;
	rtPolicy.string_value = "underlay";
	if (!ret && bobol_display_endpoint_property_set(endpoint,
	    "composition.rt.layer", &rtPolicy) != BV_DISPLAY_PROPERTY_OK) {
	    std::fprintf(stderr, "FAIL: retained RT live layer update was rejected\n");
	    ret = 1;
	}
	rtViewport = viewport_image(controller.getFramebufferUnderlayRoot());
	if (!ret && (!rtViewport ||
	    rtViewport->layer.getValue() != SoBRLViewportImage::UNDERLAY ||
	    viewport_image(controller.getFramebufferOverlayRoot()) != NULL)) {
	    std::fprintf(stderr, "FAIL: retained RT live layer update did not reparent image\n");
	    ret = 1;
	}
	/* Renderer and host images may intentionally select the same named
	 * composition layer.  The renderer image must remain first in that layer,
	 * so host framebuffer pixels composite after it regardless of which image
	 * was attached or reparented most recently. */
	{
	    BObolWindowHost framebufferHost;
	    framebufferHost.attachController(&controller, FALSE);
	    BObolFramebufferStream framebuffer(&framebufferHost);
	    std::vector<unsigned char> greenPixels(64u * 64u * 3u, 0u);
	    for (size_t i = 0; i < greenPixels.size(); i += 3u)
		greenPixels[i + 1u] = 255u;
	    if (!ret && (framebuffer.configure(64, 64) != 0 ||
		framebuffer.ensure() != 0 ||
		framebuffer.writerect(0, 0, 64, 64, greenPixels.data()) !=
		    64 * 64 || framebuffer.present() != 0)) {
		std::fprintf(stderr, "FAIL: mixed-composition framebuffer did not initialize\n");
		ret = 1;
	    }

	    std::vector<unsigned char> composite;
	    if (!ret && (!capture_composite(endpoint, composite) ||
		!mostly_green(composite))) {
		std::fprintf(stderr, "FAIL: framebuffer overlay did not cover RT underlay\n");
		ret = 1;
	    }

	    rtPolicy.string_value = "overlay";
	    if (!ret && (bobol_display_endpoint_property_set(endpoint,
		    "composition.rt.layer", &rtPolicy) !=
		    BV_DISPLAY_PROPERTY_OK ||
		!capture_composite(endpoint, composite) ||
		!mostly_green(composite))) {
		std::fprintf(stderr, "FAIL: same-layer RT reparent changed framebuffer precedence\n");
		ret = 1;
	    }

	    if (!ret && (framebuffer.setComposition(
		    BOBOL_FRAMEBUFFER_COMPOSITION_UNDERLAY) != 0 ||
		!capture_composite(endpoint, composite) ||
		count_lit(composite) == 0 || mostly_green(composite))) {
		std::fprintf(stderr, "FAIL: RT overlay did not cover framebuffer underlay\n");
		ret = 1;
	    }

	    rtPolicy.string_value = "underlay";
	    if (!ret && (bobol_display_endpoint_property_set(endpoint,
		    "composition.rt.layer", &rtPolicy) !=
		    BV_DISPLAY_PROPERTY_OK ||
		!capture_composite(endpoint, composite) ||
		!has_green_surface(composite))) {
		std::fprintf(stderr, "FAIL: same-layer framebuffer did not follow renderer image\n");
		ret = 1;
	    }

	    if (!ret && (framebuffer.setComposition(
		    BOBOL_FRAMEBUFFER_COMPOSITION_OFF) != 0 ||
		!capture_composite(endpoint, composite) ||
		count_lit(composite) == 0)) {
		std::fprintf(stderr, "FAIL: disabled framebuffer remained in composite capture\n");
		ret = 1;
	    }
	}
		void *planeSamples = NULL;
		size_t planeSize = 0;
		unsigned int planeWidth = 0;
		unsigned int planeHeight = 0;
		if (!ret && (!bobol_display_endpoint_rt_plane_capture(endpoint,
		    BOBOL_RT_OUTPUT_DEPTH, &planeSamples, &planeSize, &planeWidth,
		    &planeHeight) || !planeSamples || planeWidth != 64u ||
		    planeHeight != 64u || planeSize != 64u * 64u * sizeof(float))) {
		    std::fprintf(stderr, "FAIL: retained rt endpoint did not capture depth plane\n");
		    ret = 1;
		}
		float endpointNearestDepth = 1.0f;
		if (planeSamples) {
		    const float *depth = static_cast<const float *>(planeSamples);
		    for (size_t pixel = 0; pixel < 64u * 64u; ++pixel)
			endpointNearestDepth = std::min(endpointNearestDepth, depth[pixel]);
		    bu_free(planeSamples, "retained rt endpoint depth plane");
		}
		planeSamples = NULL;
		planeSize = 0;
		planeWidth = 0;
		planeHeight = 0;
		if (!ret && (!bobol_display_endpoint_rt_plane_capture(endpoint,
		    BOBOL_RT_OUTPUT_SOURCE_ID, &planeSamples, &planeSize, &planeWidth,
		    &planeHeight) || !planeSamples || planeWidth != 64u ||
		    planeHeight != 64u || planeSize != 64u * 64u * sizeof(uint32_t))) {
		    std::fprintf(stderr, "FAIL: retained rt endpoint did not capture identity plane\n");
		    ret = 1;
		}
		uint32_t endpointIdentity = 0u;
		if (planeSamples) {
		    const uint32_t *identity = static_cast<const uint32_t *>(planeSamples);
		    for (size_t pixel = 0; pixel < 64u * 64u; ++pixel) {
			if (identity[pixel] != 0u) {
			    endpointIdentity = identity[pixel];
			    break;
			}
		    }
		    bu_free(planeSamples, "retained rt endpoint identity plane");
		}
		if (!ret && (!endpointIdentity || endpointNearestDepth >= 1.0f)) {
		    std::fprintf(stderr, "FAIL: retained rt endpoint planes missed visible source\n");
		    ret = 1;
		}
		struct bobol_rt_source_identity endpointSource =
		    BOBOL_RT_SOURCE_IDENTITY_INIT;
		if (!ret && (!bobol_display_endpoint_rt_source_identity_get(endpoint,
		    endpointIdentity, &endpointSource) || endpointSource.database != dbip ||
		    !endpointSource.path ||
		    bu_strcmp(endpointSource.path, "ball.r") != 0)) {
		    std::fprintf(stderr, "FAIL: retained rt endpoint identity did not resolve source\n");
		    ret = 1;
		}
		bobol_display_endpoint_rt_source_identity_clear(&endpointSource);

		/* Controller policy changes can arrive outside endpoint_property_set.
		 * The owning render pass must restart RT before capture observes them. */
		/* The fixture restored lighting before creating the endpoint.  Toggle
		 * it off here so this is an actual state transition rather than a
		 * no-op whose image comparison depends on incidental earlier state. */
		controller.setLightingEnabled(FALSE);
		unsigned char *lightingOwnerImage = NULL;
		if (!ret && controller.renderToImage(&lightingOwnerImage, 0, 0, NULL,
		    NULL, NULL) != BRLCAD_OK) {
		    std::fprintf(stderr, "FAIL: controller-owned retained RT lighting redraw failed\n");
		    ret = 1;
		}
		if (lightingOwnerImage)
		    bu_free(lightingOwnerImage, "retained rt lighting redraw");
		capture = NULL;
		captureSize = 0;
		captureWidth = 0;
		captureHeight = 0;
		captureComponents = 0;
		if (!ret && (!bobol_display_endpoint_capture(endpoint, &capture,
		    &captureSize, &captureWidth, &captureHeight, &captureComponents) ||
		    !capture || !captureSize)) {
		    std::fprintf(stderr, "FAIL: retained rt endpoint did not recapture lighting\n");
		    ret = 1;
		}
		std::vector<unsigned char> litEndpointPixels;
		if (capture && captureSize)
		    litEndpointPixels.assign(capture, capture + captureSize);
		if (capture)
		    bu_free(capture, "retained rt lighting capture");
		if (!ret && litEndpointPixels == endpointPixels) {
		    std::fprintf(stderr, "FAIL: controller lighting did not refresh retained rt\n");
		    ret = 1;
		}

		/* QGED and MGED can update controller state directly instead of
		 * entering through endpoint_view_sync.  The owning render thread must
		 * notice that request, restart the retained ray job, and present its
		 * result without any Coin access from the worker thread. */
		source->color = SbColor(0.0f, 1.0f, 0.0f);
		controller.requestLodCapacityRender("database-source");
		unsigned char *ownerImage = NULL;
		if (!ret && controller.renderToImage(&ownerImage, 0, 0, NULL,
		    NULL, NULL) != BRLCAD_OK) {
		    std::fprintf(stderr, "FAIL: controller-owned retained rt redraw failed\n");
		    ret = 1;
		}
		if (ownerImage)
		    bu_free(ownerImage, "retained rt owner redraw");
		capture = NULL;
		captureSize = 0;
		captureWidth = 0;
		captureHeight = 0;
		captureComponents = 0;
		if (!ret && (!bobol_display_endpoint_capture(endpoint, &capture,
		    &captureSize, &captureWidth, &captureHeight, &captureComponents) ||
		    !capture || !captureSize)) {
		    std::fprintf(stderr, "FAIL: retained rt endpoint did not recapture redraw\n");
		    ret = 1;
		}
		endpointPixels.clear();
		if (capture && captureSize)
		    endpointPixels.assign(capture, capture + captureSize);
		if (capture)
		    bu_free(capture, "retained rt endpoint redraw capture");
		if (!ret && !has_green_surface(endpointPixels)) {
		    std::fprintf(stderr, "FAIL: controller-owned redraw did not refresh retained rt\n");
		    ret = 1;
		}
		source->color = SbColor(1.0f, 0.0f, 0.0f);
		bobol_display_endpoint_destroy(endpoint);

		/* Selection state is presentation state and follows the same retained
		 * acceleration structure as ordinary source color updates. */
		source->highlighted = TRUE;
		source->highlightedColor = SbColor(0.0f, 1.0f, 0.0f);
		if (!ret && (!renderer.synchronize(&controller) ||
		    renderer.getGeometryRevision() != geometryRevision)) {
		    std::fprintf(stderr, "FAIL: highlight update rebuilt retained rt geometry\n");
		    ret = 1;
		}
		std::vector<unsigned char> highlighted;
		if (!ret && (!renderer.render(settings, highlighted, NULL) ||
		    !has_green_surface(highlighted))) {
		    std::fprintf(stderr, "FAIL: retained rt did not apply highlight color\n");
		    ret = 1;
		}
		source->highlighted = FALSE;

		/* Placement is source presentation state.  The retained tree remains
		 * valid, while librt must trace the translated source in scene space. */
		SbMatrix placement;
		placement.setTranslate(SbVec3f(20.0f, 0.0f, 0.0f));
		source->setPlacementState(TRUE, placement, FALSE,
		    SbVec3f(0.0f, 0.0f, 0.0f), FALSE, 0.0f);
		if (!ret && (!renderer.synchronize(&controller) ||
		    renderer.getGeometryRevision() != geometryRevision)) {
		    std::fprintf(stderr, "FAIL: placement update rebuilt retained rt geometry\n");
		    ret = 1;
		}
		std::vector<unsigned char> translated;
		if (!ret && (!renderer.render(settings, translated, NULL) ||
		    has_red_surface(translated))) {
		    std::fprintf(stderr, "FAIL: retained rt did not apply source placement\n");
		    ret = 1;
		}

		std::atomic_bool cancel(true);
		BObolRtRenderStatus cancelledStatus;
		if (!ret && (renderer.render(settings, serial, &cancelledStatus,
		    &cancel) || !cancelledStatus.cancelled ||
		    cancelledStatus.complete)) {
		    std::fprintf(stderr, "FAIL: cancelled retained rt frame reported completion\n");
		    ret = 1;
		}
	    }
	}
    }

    if (!ret)
	ret = test_reflective_refractive_materials(dbip);
    if (!ret)
	ret = test_authored_retained_lights(dbip);
    if (!ret)
	ret = test_shared_render_root(dbip);

    root->unref();
    db_close(dbip);
    bu_file_delete(dbpath);
    return ret;
}
