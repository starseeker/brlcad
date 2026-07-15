/*             Q G O B O L D I S P L A Y P R O V I D E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "qtcad/display_provider.h"

#include "brlobol/display_endpoint.h"
#include "brlobol/display_session.h"
#include "bu/app.h"
#include "qtcad/QgObolWindowHost.h"

#include <QApplication>
#include <QByteArray>
#include <QCoreApplication>
#include <QEventLoop>
#include <QThread>

#include <new>
#include <limits.h>

struct QgObolDisplayProvider {
    QApplication *application = NULL;
    bool owns_application = false;
    QByteArray program_name;
    int argc = 1;
    char *argv[2] = {NULL, NULL};
};

static QApplication *
qg_obol_display_application(QgObolDisplayProvider *provider)
{
    if (!provider)
	return NULL;
    QCoreApplication *existing = QCoreApplication::instance();
    if (existing)
	return qobject_cast<QApplication *>(existing);

    const char *program_name = bu_getprogname();
    provider->program_name = program_name && program_name[0] ?
	program_name : "brlcad";
    provider->argv[0] = provider->program_name.data();
    provider->application = new (std::nothrow) QApplication(provider->argc,
	provider->argv);
    provider->owns_application = provider->application != NULL;
    return provider->application;
}

static int
qg_obol_display_provider_open(brlobol_display_endpoint_t *endpoint,
	const imgstream_fb_spec_info_t *spec, size_t width, size_t height,
	const char *title, void **instance, void *UNUSED(user_data))
{
    if (!endpoint || !spec || !instance || !width || !height ||
	width > UINT_MAX || height > UINT_MAX)
	return 0;

    const char *factory_name = NULL;
    enum brlobol_render_engine engine = BRLOBOL_RENDER_ENGINE_AUTO;
    switch (spec->display) {
	case IMGSTREAM_FB_DISPLAY_X:
	case IMGSTREAM_FB_DISPLAY_SWRAST:
	    factory_name = "qt-sw";
	    engine = BRLOBOL_RENDER_ENGINE_SW;
	    break;
	case IMGSTREAM_FB_DISPLAY_QTGL:
	case IMGSTREAM_FB_DISPLAY_OGL:
	case IMGSTREAM_FB_DISPLAY_WGL:
	    factory_name = "qt-gl";
	    engine = BRLOBOL_RENDER_ENGINE_HW;
	    break;
	case IMGSTREAM_FB_DISPLAY_NONE:
	default:
	    return 0;
    }

    QgObolDisplayProvider *provider = new (std::nothrow) QgObolDisplayProvider;
    if (!provider)
	return 0;
    QApplication *application = qg_obol_display_application(provider);
    if (!application || QThread::currentThread() != application->thread() ||
	!qtcad_obol_host_factories_register() ||
	!brlobol_display_endpoint_render_engine_set(endpoint, engine)) {
	if (provider->owns_application)
	    delete provider->application;
	delete provider;
	return 0;
    }

    struct brlobol_host_desc desc = {};
    desc.struct_size = sizeof(desc);
    desc.mode = BRLOBOL_HOST_MODE_TOPLEVEL;
    desc.width = static_cast<unsigned int>(width);
    desc.height = static_cast<unsigned int>(height);
    desc.device_pixel_ratio = 1.0;
    desc.visible = 1;
    desc.title = title ? title : "BRL-CAD framebuffer";
    desc.vsync = BRLOBOL_HOST_VSYNC_AUTO;
    if (!brlobol_display_endpoint_host_open(endpoint, factory_name, &desc)) {
	if (provider->owns_application)
	    delete provider->application;
	delete provider;
	return 0;
    }

    *instance = provider;
    return 1;
}

static void
qg_obol_display_provider_close(void *instance, void *UNUSED(user_data))
{
    QgObolDisplayProvider *provider =
	static_cast<QgObolDisplayProvider *>(instance);
    if (!provider)
	return;
    if (provider->owns_application)
	delete provider->application;
    delete provider;
}

static int
qg_obol_display_provider_poll(void *instance, void *UNUSED(user_data))
{
    QgObolDisplayProvider *provider =
	static_cast<QgObolDisplayProvider *>(instance);
    QApplication *application = provider ?
	qobject_cast<QApplication *>(QCoreApplication::instance()) : NULL;
    if (!application || QThread::currentThread() != application->thread())
	return -1;
    QCoreApplication::processEvents(QEventLoop::AllEvents, 0);
    return 0;
}

static long
qg_obol_display_provider_poll_rate(const void *UNUSED(instance),
	void *UNUSED(user_data))
{
    return 16000;
}

extern "C" QTCAD_EXPORT int
qtcad_obol_display_provider_register(void)
{
    static const brlobol_display_provider_t provider = {
	BRLOBOL_DISPLAY_PROVIDER_ABI_VERSION,
	sizeof(brlobol_display_provider_t),
	"qt-obol",
	100,
	NULL,
	qg_obol_display_provider_open,
	qg_obol_display_provider_close,
	qg_obol_display_provider_poll,
	qg_obol_display_provider_poll_rate
    };
    return brlobol_display_provider_register(&provider);
}
