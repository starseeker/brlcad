/*                 T E S T _ R T _ R E N D E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/init.h"
#include "brlobol/display_endpoint.h"
#include "brlobol/rt_render.h"
#include "brlobol/view_controller.h"
#include "brlobol/viewport_image.h"

#include "bu/app.h"
#include "bu/file.h"
#include "bu/malloc.h"
#include "raytrace.h"
#include "wdb.h"

#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/SbMatrix.h>

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstring>
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
    back_center[Z] = -3.0;
    const int ret = mk_sph(wdbp, "ball.s", center, 2.0) == 0 &&
	mk_sph(wdbp, "back.s", back_center, 2.0) == 0 ? 0 : -1;
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

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    (void)argc;
    brlobol_init(NULL);

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
	BRLObolViewController controller(root, camera);
	controller.setViewportSize(64, 64);
	controller.setBackgroundColors(SbColor(0.0f, 0.0f, 0.0f),
	    SbColor(0.0f, 0.0f, 0.0f));
	if (controller.replaceDatabaseSource("ball.s", dbip,
		SoBRLDatabaseSource::SHADED, 1) != 1) {
	    std::fprintf(stderr, "FAIL: failed to attach retained rt source\n");
	    ret = 1;
	} else {
	    BRLObolRtRenderer renderer;
	    if (!renderer.synchronize(&controller) ||
		renderer.getPreparedSourceCount() != 1) {
		std::fprintf(stderr, "FAIL: retained rt renderer did not prepare source\n");
		ret = 1;
	    } else {
		uint64_t geometryRevision = renderer.getGeometryRevision();
		const uint64_t presentationRevision = renderer.getPresentationRevision();
		BRLObolRtRenderSettings settings;
		settings.width = 64;
		settings.height = 64;
		settings.workers = 2;
		settings.samples = 1;
		std::vector<unsigned char> first;
		BRLObolRtRenderStatus firstStatus;
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
		BRLObolRtRenderStatus serialStatus;
		if (!ret && (!renderer.render(settings, serial, &serialStatus) ||
		    serial != first || !serialStatus.complete)) {
		    std::fprintf(stderr, "FAIL: retained rt final frame is not deterministic\n");
		    ret = 1;
		}
		std::vector<unsigned char> tiled;
		for (unsigned int row = 0; !ret && row < settings.height; row += 7u) {
		    BRLObolRtRenderStatus tileStatus;
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
		BRLObolRtRenderPlanes planes;
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
		BRLObolRtSourceIdentity sourceIdentity;
		if (!ret && (!hitIdentity || nearestDepth >= 1.0f ||
		    !renderer.getSourceIdentity(hitIdentity, sourceIdentity) ||
		    sourceIdentity.database != dbip ||
		    std::strcmp(sourceIdentity.path.getString(), "ball.s") != 0)) {
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
		if (!ret && (!renderer.synchronize(&controller) ||
		    renderer.getGeometryRevision() != geometryRevision ||
		    renderer.getPresentationRevision() == presentationRevision)) {
		    std::fprintf(stderr, "FAIL: appearance update rebuilt retained rt geometry\n");
		    ret = 1;
		}
		std::vector<unsigned char> red;
		if (!ret && (!renderer.render(settings, red, NULL) ||
		    !has_red_surface(red))) {
		    std::fprintf(stderr, "FAIL: retained rt appearance update did not affect frame\n");
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
		brlobol_display_endpoint_t *endpoint =
		    brlobol_display_endpoint_create(&controller, 0);
		struct brlobol_endpoint_property_value rtPolicy =
		    BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
		rtPolicy.type = BRLOBOL_ENDPOINT_PROPERTY_UINT;
		rtPolicy.uint_value = 1;
		if (!ret && (!endpoint || brlobol_display_endpoint_property_set(
		    endpoint, "render.rt.frame_budget_ms", &rtPolicy) !=
		    BRLOBOL_ENDPOINT_PROPERTY_OK)) {
		    std::fprintf(stderr, "FAIL: retained rt frame-budget policy was rejected\n");
		    ret = 1;
		}
		rtPolicy = BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
		rtPolicy.type = BRLOBOL_ENDPOINT_PROPERTY_ENUM;
		rtPolicy.string_value = "interactive";
	if (!ret && brlobol_display_endpoint_property_set(endpoint,
	    "render.rt.quality", &rtPolicy) != BRLOBOL_ENDPOINT_PROPERTY_OK) {
		    std::fprintf(stderr, "FAIL: retained rt quality policy was rejected\n");
	    ret = 1;
	}
	rtPolicy = BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	rtPolicy.type = BRLOBOL_ENDPOINT_PROPERTY_ENUM;
	rtPolicy.string_value = "overlay";
	if (!ret && brlobol_display_endpoint_property_set(endpoint,
	    "composition.rt.layer", &rtPolicy) != BRLOBOL_ENDPOINT_PROPERTY_OK) {
	    std::fprintf(stderr, "FAIL: retained RT image-layer policy was rejected\n");
	    ret = 1;
	}
		rtPolicy = BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
		if (!ret && (brlobol_display_endpoint_property_get(endpoint,
		    "render.rt.frame_budget_ms", &rtPolicy) !=
		    BRLOBOL_ENDPOINT_PROPERTY_OK || rtPolicy.uint_value != 1u)) {
		    std::fprintf(stderr, "FAIL: retained rt frame-budget policy did not persist\n");
		    ret = 1;
		}
		rtPolicy = BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	if (!ret && (brlobol_display_endpoint_property_get(endpoint,
	    "render.rt.quality", &rtPolicy) != BRLOBOL_ENDPOINT_PROPERTY_OK ||
		    !rtPolicy.string_value ||
		    std::strcmp(rtPolicy.string_value, "interactive") != 0)) {
		    std::fprintf(stderr, "FAIL: retained rt quality policy did not persist\n");
	    ret = 1;
	}
	rtPolicy = BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	if (!ret && (brlobol_display_endpoint_property_get(endpoint,
	    "composition.rt.layer", &rtPolicy) != BRLOBOL_ENDPOINT_PROPERTY_OK ||
	    !rtPolicy.string_value || std::strcmp(rtPolicy.string_value, "overlay") != 0)) {
	    std::fprintf(stderr, "FAIL: retained RT image-layer policy did not persist\n");
	    ret = 1;
	}
		if (!ret && (!endpoint || !brlobol_display_endpoint_render_engine_set(
		    endpoint, BRLOBOL_RENDER_ENGINE_RT))) {
		    std::fprintf(stderr, "FAIL: retained rt endpoint did not initialize\n");
		    ret = 1;
		}
		unsigned char *capture = NULL;
		size_t captureSize = 0;
		unsigned int captureWidth = 0;
		unsigned int captureHeight = 0;
		unsigned int captureComponents = 0;
		if (!ret && (!brlobol_display_endpoint_capture(endpoint, &capture,
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
	if (!ret && !has_red_surface(endpointPixels)) {
		    std::fprintf(stderr, "FAIL: retained rt endpoint did not present ray image\n");
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
	rtPolicy = BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	rtPolicy.type = BRLOBOL_ENDPOINT_PROPERTY_ENUM;
	rtPolicy.string_value = "underlay";
	if (!ret && brlobol_display_endpoint_property_set(endpoint,
	    "composition.rt.layer", &rtPolicy) != BRLOBOL_ENDPOINT_PROPERTY_OK) {
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
		void *planeSamples = NULL;
		size_t planeSize = 0;
		unsigned int planeWidth = 0;
		unsigned int planeHeight = 0;
		if (!ret && (!brlobol_display_endpoint_rt_plane_capture(endpoint,
		    BRLOBOL_RT_OUTPUT_DEPTH, &planeSamples, &planeSize, &planeWidth,
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
		if (!ret && (!brlobol_display_endpoint_rt_plane_capture(endpoint,
		    BRLOBOL_RT_OUTPUT_SOURCE_ID, &planeSamples, &planeSize, &planeWidth,
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
		struct brlobol_rt_source_identity endpointSource =
		    BRLOBOL_RT_SOURCE_IDENTITY_INIT;
		if (!ret && (!brlobol_display_endpoint_rt_source_identity_get(endpoint,
		    endpointIdentity, &endpointSource) || endpointSource.database != dbip ||
		    !endpointSource.path ||
		    std::strcmp(endpointSource.path, "ball.s") != 0)) {
		    std::fprintf(stderr, "FAIL: retained rt endpoint identity did not resolve source\n");
		    ret = 1;
		}
		brlobol_display_endpoint_rt_source_identity_clear(&endpointSource);

		/* Controller policy changes can arrive outside endpoint_property_set.
		 * The owning render pass must restart RT before capture observes them. */
		controller.setLightingEnabled(TRUE);
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
		if (!ret && (!brlobol_display_endpoint_capture(endpoint, &capture,
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
		controller.requestRender("database-source");
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
		if (!ret && (!brlobol_display_endpoint_capture(endpoint, &capture,
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
		brlobol_display_endpoint_destroy(endpoint);

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
		BRLObolRtRenderStatus cancelledStatus;
		if (!ret && (renderer.render(settings, serial, &cancelledStatus,
		    &cancel) || !cancelledStatus.cancelled ||
		    cancelledStatus.complete)) {
		    std::fprintf(stderr, "FAIL: cancelled retained rt frame reported completion\n");
		    ret = 1;
		}
	    }
	}
    }

    root->unref();
    db_close(dbip);
    bu_file_delete(dbpath);
    return ret;
}
