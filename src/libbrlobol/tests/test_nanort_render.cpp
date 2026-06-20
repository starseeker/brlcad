/*                 T E S T _ N A N O R T _ R E N D E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol.h"

#include "bu/app.h"
#include "bu/file.h"
#include "raytrace.h"
#include "wdb.h"

#include <Inventor/SbColor.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/SoOffscreenRenderer.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/annex/HUD/nodekits/SoHUDKit.h>
#include <Inventor/annex/HUD/nodes/SoHUDLabel.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/nodes/SoSeparator.h>
#include <OSMesa/osmesa.h>

#include "nanort_context_manager.h"

#include <limits.h>
#include <memory>
#include <stdio.h>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

class BRLOBOLRenderContext {
public:
    BRLOBOLRenderContext(unsigned int w, unsigned int h) :
	context(NULL),
	buffer(NULL),
	width(w),
	height(h),
	previousContext(NULL),
	previousBuffer(NULL),
	previousWidth(0),
	previousHeight(0),
	previousFormat(0)
    {
	context = OSMesaCreateContextExt(OSMESA_RGBA, 24, 0, 0, NULL);
	if (context)
	    buffer.reset(new unsigned char[width * height * 4]);
    }

    ~BRLOBOLRenderContext(void)
    {
	if (context)
	    OSMesaDestroyContext(context);
    }

    bool isValid(void) const
    {
	return context != NULL && buffer.get() != NULL;
    }

    SbBool makeCurrent(void)
    {
	if (!this->isValid())
	    return FALSE;

	previousContext = OSMesaGetCurrentContext();
	previousBuffer = NULL;
	previousWidth = 0;
	previousHeight = 0;
	previousFormat = 0;
	if (previousContext) {
	    GLint fmt = 0;
	    OSMesaGetColorBuffer(previousContext, &previousWidth, &previousHeight,
		    &fmt, &previousBuffer);
	    previousFormat = (GLenum)fmt;
	}

	return OSMesaMakeCurrent(context, buffer.get(), GL_UNSIGNED_BYTE,
		static_cast<GLsizei>(width), static_cast<GLsizei>(height)) ? TRUE : FALSE;
    }

    void restorePrevious(void)
    {
	if (previousContext && previousBuffer) {
	    OSMesaMakeCurrent(previousContext, previousBuffer, GL_UNSIGNED_BYTE,
		    previousWidth, previousHeight);
	} else {
	    OSMesaMakeCurrent(NULL, NULL, 0, 0, 0);
	}
    }

    OSMesaContext context;
    std::unique_ptr<unsigned char[]> buffer;
    unsigned int width;
    unsigned int height;
    OSMesaContext previousContext;
    void *previousBuffer;
    GLsizei previousWidth;
    GLsizei previousHeight;
    GLenum previousFormat;
};

class BRLOBOLRenderContextManager : public SoDB::ContextManager {
public:
    virtual void *createOffscreenContext(unsigned int width, unsigned int height)
    {
	BRLOBOLRenderContext *ctx = new BRLOBOLRenderContext(width, height);
	if (ctx->isValid())
	    return ctx;
	delete ctx;
	return NULL;
    }

    virtual SbBool isOSMesaContext(void *UNUSED(context))
    {
	return TRUE;
    }

    virtual SbBool makeContextCurrent(void *context)
    {
	return context ? static_cast<BRLOBOLRenderContext *>(context)->makeCurrent() : FALSE;
    }

    virtual void restorePreviousContext(void *context)
    {
	if (context)
	    static_cast<BRLOBOLRenderContext *>(context)->restorePrevious();
    }

    virtual void destroyContext(void *context)
    {
	BRLOBOLRenderContext *ctx = static_cast<BRLOBOLRenderContext *>(context);
	if (ctx && ctx->context && OSMesaGetCurrentContext() == ctx->context)
	    OSMesaMakeCurrent(NULL, NULL, 0, 0, 0);
	delete ctx;
    }

    virtual void *getProcAddress(const char *funcName)
    {
	return reinterpret_cast<void *>(OSMesaGetProcAddress(funcName));
    }
};

static int
write_backend_test_db(char *dbpath, size_t dbpath_len)
{
    FILE *fp = bu_temp_file(dbpath, dbpath_len);
    if (!fp)
	return 0;
    fclose(fp);

    struct rt_wdb *wdbp = wdb_fopen(dbpath);
    if (!wdbp)
	return 0;

    point_t sph_center = {0.0, 0.0, 0.0};
    if (mk_sph(wdbp, "backend_ball.s", sph_center, 2.0) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    wdb_close(wdbp);
    return 1;
}

static int
open_database(const char *dbpath, struct db_i **dbipp)
{
    struct db_i *dbip = db_open(dbpath, DB_OPEN_READONLY);
    if (!dbip)
	return 0;
    if (db_dirbuild(dbip) < 0) {
	db_close(dbip);
	return 0;
    }
    *dbipp = dbip;
    return 1;
}

static int
count_lit_pixels(const unsigned char *buffer, int width, int height)
{
    int count = 0;
    for (int i = 0; i < width * height; i++) {
	const unsigned char *p = buffer + i * 3;
	if (p[0] > 32 || p[1] > 32 || p[2] > 32)
	    count++;
    }
    return count;
}

static int
count_red_pixels(const unsigned char *buffer, int width, int height)
{
    int count = 0;
    for (int i = 0; i < width * height; i++) {
	const unsigned char *p = buffer + i * 3;
	if (p[0] > 96 && p[1] < 80 && p[2] < 80)
	    count++;
    }
    return count;
}

static int
count_green_pixels(const unsigned char *buffer, int width, int height)
{
    int count = 0;
    for (int i = 0; i < width * height; i++) {
	const unsigned char *p = buffer + i * 3;
	if (p[1] > 96 && p[0] < 80 && p[2] < 80)
	    count++;
    }
    return count;
}

int
main(int UNUSED(argc), const char **UNUSED(argv))
{
    const int width = 180;
    const int height = 140;
    char dbpath[MAXPATHLEN] = {0};

    if (!write_backend_test_db(dbpath, MAXPATHLEN))
	FAIL("failed to create NanoRT backend test database");

    struct db_i *dbip = NULL;
    if (!open_database(dbpath, &dbip)) {
	bu_file_delete(dbpath);
	FAIL("failed to open NanoRT backend test database");
    }

    SoNanoRTContextManager nanortManager;
    brlobol_init(&nanortManager);
    if (SoDB::getContextManager() != &nanortManager)
	FAIL("brlobol_init should install the application-provided NanoRT context manager");

    SoSeparator *root = new SoSeparator;
    root->ref();

    SoPerspectiveCamera *camera = new SoPerspectiveCamera;
    camera->position = SbVec3f(0.0f, 0.0f, 8.0f);
    camera->nearDistance = 1.0f;
    camera->farDistance = 20.0f;
    root->addChild(camera);

    SoDirectionalLight *light = new SoDirectionalLight;
    light->direction = SbVec3f(-0.25f, -0.45f, -1.0f);
    light->intensity = 0.95f;
    root->addChild(light);

    SoMaterial *material = new SoMaterial;
    material->diffuseColor = SbColor(0.85f, 0.05f, 0.05f);
    material->emissiveColor = SbColor(0.25f, 0.0f, 0.0f);
    material->ambientColor = SbColor(0.20f, 0.0f, 0.0f);
    root->addChild(material);

    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->setDatabase(dbip);
    source->path = "backend_ball.s";
    source->drawMode = SoBRLDatabaseSource::SHADED;
    source->sourceRevision = 1;
    root->addChild(source);

    SoHUDKit *hud = new SoHUDKit;
    root->addChild(hud);

    SoHUDLabel *status = new SoHUDLabel;
    status->position = SbVec2f(8.0f, 116.0f);
    status->string.set1Value(0, "NanoRT backend");
    status->color = SbColor(0.0f, 1.0f, 0.0f);
    status->fontSize = 13.0f;
    hud->addWidget(status);

    SoBRLRealizeAction realize;
    realize.apply(root);
    if (realize.getRealizedSourceCount() != 1 ||
	    source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED)
	FAIL("NanoRT scene should realize database-backed shaded source");

    SoBRLMeshShape *mesh = source->getRealizedMesh();
    if (!mesh || mesh->getTriangleCount() <= 12)
	FAIL("NanoRT scene should expose tessellated BRL-CAD mesh geometry");
    mesh->colorOverride = TRUE;
    mesh->color = SbColor(0.85f, 0.05f, 0.05f);

    SbViewportRegion viewport(width, height);
    BRLOBOLRenderContextManager osmesaManager;
    SoOffscreenRenderer glRenderer(&osmesaManager, viewport);
    glRenderer.setComponents(SoOffscreenRenderer::RGB);
    glRenderer.setBackgroundColor(SbColor(0.0f, 0.0f, 0.0f));
    if (!glRenderer.render(root))
	FAIL("OSMesa should render the same database-backed BRL-CAD Obol scene");

    const unsigned char *buffer = glRenderer.getBuffer();
    if (!buffer)
	FAIL("OSMesa renderer should expose framebuffer readback");

    if (count_lit_pixels(buffer, width, height) < 150)
	FAIL("OSMesa backend should produce visible framebuffer pixels");
    if (count_red_pixels(buffer, width, height) < 60)
	FAIL("OSMesa backend should render material-colored database mesh pixels");
    if (count_green_pixels(buffer, width, height) < 20)
	FAIL("OSMesa backend should render direct Obol HUD overlays");

    SoOffscreenRenderer renderer(&nanortManager, viewport);
    renderer.setComponents(SoOffscreenRenderer::RGB);
    renderer.setBackgroundColor(SbColor(0.0f, 0.0f, 0.0f));
    if (!renderer.render(root))
	FAIL("NanoRT ContextManager should render database-backed BRL-CAD Obol scene");

    buffer = renderer.getBuffer();
    if (!buffer)
	FAIL("NanoRT renderer should expose framebuffer readback");

    if (count_lit_pixels(buffer, width, height) < 150)
	FAIL("NanoRT backend should produce visible framebuffer pixels");
    if (count_red_pixels(buffer, width, height) < 60)
	FAIL("NanoRT backend should render material-colored database mesh pixels");
    if (count_green_pixels(buffer, width, height) < 20)
	FAIL("NanoRT backend should composite direct Obol HUD overlays");

    SoRayPickAction backendPick(viewport);
    backendPick.setRay(SbVec3f(0.0f, 0.0f, 8.0f), SbVec3f(0.0f, 0.0f, -1.0f));
    backendPick.apply(root);
    const SoPickedPoint *pickedPoint = backendPick.getPickedPoint();
    if (!pickedPoint)
	FAIL("NanoRT-rendered scene should remain pickable through Obol traversal");
    const SoDetail *rawDetail = pickedPoint->getDetail(mesh);
    if (!rawDetail || !rawDetail->isOfType(SoBRLPickDetail::getClassTypeId()))
	FAIL("NanoRT-rendered scene pick should return a BRL-CAD Obol detail");
    const SoBRLPickDetail *pickDetail = static_cast<const SoBRLPickDetail *>(rawDetail);
    if (strcmp(pickDetail->getPath().getString(), "/backend_ball.s") != 0 ||
	    pickDetail->getPrimitiveKind() != SoBRLPickDetail::FACE)
	FAIL("NanoRT-rendered scene pick detail should preserve mesh face identity");

    char rgbpath[MAXPATHLEN] = {0};
    FILE *rgbfp = bu_temp_file(rgbpath, MAXPATHLEN);
    if (!rgbfp)
	FAIL("failed to create NanoRT diagnostic RGB output path");
    fclose(rgbfp);

    if (!renderer.writeToRGB(rgbpath)) {
	bu_file_delete(rgbpath);
	FAIL("NanoRT renderer should write diagnostic RGB output without the display-manager path");
    }
    int rgbsize = bu_file_size(rgbpath);
    bu_file_delete(rgbpath);
    if (rgbsize <= 512)
	FAIL("NanoRT diagnostic RGB output should contain image data");

    root->unref();
    db_close(dbip);
    bu_file_delete(dbpath);
    return 0;
}
