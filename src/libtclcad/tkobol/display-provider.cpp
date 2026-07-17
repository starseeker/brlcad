/*              D I S P L A Y - P R O V I D E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BDisplaySession.h"
#include "BObol/BHostFactory.h"
#include "bu/app.h"
#include "bu/log.h"
#include "bu/vls.h"
#include "tclcad/setup.h"

#include <limits.h>
#include <new>

struct TkObolDisplayProvider {
    Tcl_Interp *interp = NULL;
};

static int
tk_obol_display_provider_open(bobol_display_endpoint_t *endpoint,
	const imgstream_fb_spec_info_t *spec, size_t width, size_t height,
	const char *title, void **instance, void *UNUSED(user_data))
{
    if (!endpoint || !spec || !instance || !width || !height ||
	width > UINT_MAX || height > UINT_MAX)
	return 0;

    const char *factory_name = NULL;
    enum bobol_render_engine engine = BOBOL_RENDER_ENGINE_AUTO;
    switch (spec->display) {
	case IMGSTREAM_FB_DISPLAY_X:
	case IMGSTREAM_FB_DISPLAY_SWRAST:
	    factory_name = "tk-photo";
	    engine = BOBOL_RENDER_ENGINE_SW;
	    break;
	case IMGSTREAM_FB_DISPLAY_QTGL:
	case IMGSTREAM_FB_DISPLAY_OGL:
	case IMGSTREAM_FB_DISPLAY_WGL:
	    factory_name = "tk-gl";
	    engine = BOBOL_RENDER_ENGINE_HW;
	    break;
	case IMGSTREAM_FB_DISPLAY_NONE:
	default:
	    return 0;
    }

    TkObolDisplayProvider *provider = new (std::nothrow) TkObolDisplayProvider;
    if (!provider)
	return 0;

    struct bu_vls log = BU_VLS_INIT_ZERO;

    Tcl_FindExecutable(bu_getprogname());
    provider->interp = Tcl_CreateInterp();
    int opened = 0;
    const char *stage = "creating Tcl interpreter";
    if (provider->interp &&
	((stage = "initializing Tcl/Tk"),
	 tclcad_init(provider->interp, 1, &log) == TCL_OK) &&
	((stage = "registering Tk Obol host factories"),
	 tclcad_obol_host_factories_register()) &&
	((stage = "selecting the requested Obol renderer"),
	 bobol_display_endpoint_render_engine_set(endpoint, engine))) {
	struct bobol_host_desc desc = {};
	desc.struct_size = sizeof(desc);
	desc.mode = BOBOL_HOST_MODE_TOPLEVEL;
	desc.width = static_cast<unsigned int>(width);
	desc.height = static_cast<unsigned int>(height);
	desc.device_pixel_ratio = 1.0;
	desc.visible = 1;
	desc.title = title ? title : "BRL-CAD framebuffer";
	desc.application_context = provider->interp;
	desc.vsync = BOBOL_HOST_VSYNC_AUTO;
	stage = "opening the Tk Obol host";
	opened = bobol_display_endpoint_host_open(endpoint, factory_name, &desc);
    }
    if (opened) {
	bu_vls_free(&log);
	*instance = provider;
	return 1;
    }
    bu_log("Tk Obol display provider failed while %s%s%s\n", stage,
	bu_vls_strlen(&log) ? ": " : "",
	bu_vls_strlen(&log) ? bu_vls_cstr(&log) : "");
    bu_vls_free(&log);
    if (provider->interp)
	Tcl_DeleteInterp(provider->interp);
    delete provider;
    return 0;
}

static void
tk_obol_display_provider_close(void *instance, void *UNUSED(user_data))
{
    TkObolDisplayProvider *provider =
	static_cast<TkObolDisplayProvider *>(instance);
    if (!provider)
	return;
    if (provider->interp)
	Tcl_DeleteInterp(provider->interp);
    delete provider;
}

static int
tk_obol_display_provider_poll(void *instance, void *UNUSED(user_data))
{
    TkObolDisplayProvider *provider =
	static_cast<TkObolDisplayProvider *>(instance);
    if (!provider || !provider->interp)
	return -1;
    for (int count = 0; count < 64 && Tcl_DoOneEvent(TCL_DONT_WAIT); count++)
	;
    return 0;
}

static long
tk_obol_display_provider_poll_rate(const void *UNUSED(instance),
	void *UNUSED(user_data))
{
    return 16000;
}

extern "C" TCLCAD_EXPORT int
tclcad_obol_display_provider_register(void)
{
    static const bobol_display_provider_t provider = {
	BOBOL_DISPLAY_PROVIDER_ABI_VERSION,
	sizeof(bobol_display_provider_t),
	"tk-obol",
	100,
	NULL,
	tk_obol_display_provider_open,
	tk_obol_display_provider_close,
	tk_obol_display_provider_poll,
	tk_obol_display_provider_poll_rate
    };
    return bobol_display_provider_register(&provider);
}
