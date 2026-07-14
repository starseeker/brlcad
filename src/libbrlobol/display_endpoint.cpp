/*              D I S P L A Y _ E N D P O I N T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/str.h"

#include "brlobol/display_endpoint.h"
#include "brlobol/init.h"
#include "brlobol/scene_group.h"
#include "brlobol/view_controller.h"
#include "brlobol/window_host.h"

#include <cmath>
#include <cstring>
#include <new>
#include <string>

struct brlobol_display_endpoint {
    brlobol_display_endpoint(void) :
	controller(NULL),
	owns_controller(false),
	host(NULL),
	owns_host(false),
	factory(NULL),
	factory_instance(NULL),
	engine(BRLOBOL_RENDER_ENGINE_AUTO),
	host_mode(BRLOBOL_HOST_MODE_HEADLESS),
	width(1),
	height(1),
	device_pixel_ratio(1.0),
	visible(false),
	vsync(true),
	framebuffer_capture_callback(NULL),
	framebuffer_capture_user_data(NULL),
	factory_frame_callback_bound(false),
	property_get_callback(NULL),
	property_set_callback(NULL),
	property_user_data(NULL)
    {
	input.setProfile(brlobol_input_default_view_profile());
    }

    BRLObolViewController *controller;
    bool owns_controller;
    BRLObolWindowHost *host;
    bool owns_host;
    brlobol_host_factory_token_t *factory;
    void *factory_instance;
    enum brlobol_render_engine engine;
    enum brlobol_host_mode host_mode;
    unsigned int width;
    unsigned int height;
    double device_pixel_ratio;
    std::string title;
    bool visible;
    bool vsync;
    brlobol_endpoint_framebuffer_capture_callback framebuffer_capture_callback;
    void *framebuffer_capture_user_data;
    bool factory_frame_callback_bound;
    BRLObolInputContext input;
    brlobol_endpoint_property_get_callback property_get_callback;
    brlobol_endpoint_property_set_callback property_set_callback;
    void *property_user_data;
};

static bool valid_engine(enum brlobol_render_engine engine);

static void
endpoint_frame_requested(void *user_data, const char *reason)
{
    brlobol_display_endpoint_t *endpoint =
	static_cast<brlobol_display_endpoint_t *>(user_data);
    if (!endpoint || !endpoint->factory || !endpoint->factory_instance)
	return;
    (void)brlobol_host_factory_instance_request_frame(endpoint->factory,
	endpoint->factory_instance, reason);
}

static int
endpoint_input_dispatch(void *user_data, const BRLObolInputEvent *event)
{
    return brlobol_display_endpoint_input_dispatch(
	static_cast<brlobol_display_endpoint_t *>(user_data), event);
}

static void
endpoint_dimensions_refresh(brlobol_display_endpoint_t *endpoint)
{
    if (!endpoint || !endpoint->factory || !endpoint->factory_instance)
	return;
    unsigned int width = 0;
    unsigned int height = 0;
    double device_pixel_ratio = 0.0;
    if (!brlobol_host_factory_instance_dimensions(endpoint->factory,
	endpoint->factory_instance, &width, &height, &device_pixel_ratio))
	return;
    if (width)
	endpoint->width = width;
    if (height)
	endpoint->height = height;
    if (device_pixel_ratio > 0.0)
	endpoint->device_pixel_ratio = device_pixel_ratio;
}

#define FACEPLATE_AXES_STYLE_PROPERTIES(_axis) \
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate." _axis ".position.x", \
	BRLOBOL_ENDPOINT_PROPERTY_DOUBLE, \
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE, \
	0, -1.0e12, 1.0e12, NULL}, \
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate." _axis ".position.y", \
	BRLOBOL_ENDPOINT_PROPERTY_DOUBLE, \
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE, \
	0, -1.0e12, 1.0e12, NULL}, \
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate." _axis ".position.z", \
	BRLOBOL_ENDPOINT_PROPERTY_DOUBLE, \
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE, \
	0, -1.0e12, 1.0e12, NULL}, \
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate." _axis ".size", \
	BRLOBOL_ENDPOINT_PROPERTY_DOUBLE, \
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE, \
	0, 0.0, 1.0e12, NULL}, \
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate." _axis ".line_width", \
	BRLOBOL_ENDPOINT_PROPERTY_UINT, \
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE, \
	0, 1.0, 1048576.0, NULL}, \
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate." _axis ".color", \
	BRLOBOL_ENDPOINT_PROPERTY_COLOR3, \
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE, \
	0, 0.0, 1.0, NULL}, \
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate." _axis ".position_only", \
	BRLOBOL_ENDPOINT_PROPERTY_BOOL, \
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE, \
	0, 0.0, 1.0, NULL}, \
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate." _axis ".labels.visible", \
	BRLOBOL_ENDPOINT_PROPERTY_BOOL, \
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE, \
	0, 0.0, 1.0, NULL}, \
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate." _axis ".labels.color", \
	BRLOBOL_ENDPOINT_PROPERTY_COLOR3, \
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE, \
	0, 0.0, 1.0, NULL}, \
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate." _axis ".triple_color", \
	BRLOBOL_ENDPOINT_PROPERTY_BOOL, \
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE, \
	0, 0.0, 1.0, NULL}, \
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate." _axis ".ticks.visible", \
	BRLOBOL_ENDPOINT_PROPERTY_BOOL, \
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE, \
	0, 0.0, 1.0, NULL}, \
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate." _axis ".ticks.length", \
	BRLOBOL_ENDPOINT_PROPERTY_UINT, \
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE, \
	0, 1.0, 1048576.0, NULL}, \
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate." _axis ".ticks.major_length", \
	BRLOBOL_ENDPOINT_PROPERTY_UINT, \
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE, \
	0, 1.0, 1048576.0, NULL}, \
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate." _axis ".ticks.interval", \
	BRLOBOL_ENDPOINT_PROPERTY_DOUBLE, \
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE, \
	0, 0.0, 1.0e12, NULL}, \
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate." _axis ".ticks.per_major", \
	BRLOBOL_ENDPOINT_PROPERTY_UINT, \
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE, \
	0, 0.0, 1048576.0, NULL}, \
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate." _axis ".ticks.threshold", \
	BRLOBOL_ENDPOINT_PROPERTY_UINT, \
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE, \
	0, 0.0, 2147483647.0, NULL}, \
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate." _axis ".ticks.color", \
	BRLOBOL_ENDPOINT_PROPERTY_COLOR3, \
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE, \
	0, 0.0, 1.0, NULL}, \
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate." _axis ".ticks.major_color", \
	BRLOBOL_ENDPOINT_PROPERTY_COLOR3, \
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE, \
	0, 0.0, 1.0, NULL}

static const brlobol_endpoint_property_desc endpoint_properties[] = {
    {sizeof(brlobol_endpoint_property_desc), "endpoint.renderer",
	BRLOBOL_ENDPOINT_PROPERTY_ENUM,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 0.0, "auto,hw,sw,rt,none,diagnostic"},
    {sizeof(brlobol_endpoint_property_desc), "endpoint.width",
	BRLOBOL_ENDPOINT_PROPERTY_UINT,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 1.0, 32767.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "endpoint.height",
	BRLOBOL_ENDPOINT_PROPERTY_UINT,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 1.0, 32767.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "endpoint.device_pixel_ratio",
	BRLOBOL_ENDPOINT_PROPERTY_DOUBLE,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.01, 64.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "endpoint.host",
	BRLOBOL_ENDPOINT_PROPERTY_STRING, BRLOBOL_ENDPOINT_PROPERTY_READ,
	0, 0.0, 0.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "endpoint.title",
	BRLOBOL_ENDPOINT_PROPERTY_STRING,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	BRLOBOL_HOST_CAP_TOPLEVEL, 0.0, 0.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "endpoint.visible",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	BRLOBOL_HOST_CAP_TOPLEVEL, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "endpoint.vsync",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	BRLOBOL_HOST_CAP_PRESENT_VSYNC, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "controller.background.bottom",
	BRLOBOL_ENDPOINT_PROPERTY_COLOR3,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "controller.background.top",
	BRLOBOL_ENDPOINT_PROPERTY_COLOR3,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "controller.software_wire",
	BRLOBOL_ENDPOINT_PROPERTY_ENUM,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 0.0, "auto,quality,fast"},
    {sizeof(brlobol_endpoint_property_desc), "renderer.depth_test",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "renderer.lighting",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "renderer.transparency",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "renderer.antialiasing",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "renderer.clip.minimum",
	BRLOBOL_ENDPOINT_PROPERTY_DOUBLE,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, -1.0e12, 1.0e12, NULL},
    {sizeof(brlobol_endpoint_property_desc), "renderer.clip.maximum",
	BRLOBOL_ENDPOINT_PROPERTY_DOUBLE,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, -1.0e12, 1.0e12, NULL},
    {sizeof(brlobol_endpoint_property_desc), "renderer.depth_cue",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.perspective",
	BRLOBOL_ENDPOINT_PROPERTY_DOUBLE,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 179.999, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.zclip",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.navigation.min_delta",
	BRLOBOL_ENDPOINT_PROPERTY_DOUBLE,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, -10000.0, 0.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.navigation.max_delta",
	BRLOBOL_ENDPOINT_PROPERTY_DOUBLE,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 10000.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.navigation.rotate_scale",
	BRLOBOL_ENDPOINT_PROPERTY_DOUBLE,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.001, 10.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.navigation.scale_scale",
	BRLOBOL_ENDPOINT_PROPERTY_DOUBLE,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.001, 100.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.adc.visible",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.adc.line_color",
	BRLOBOL_ENDPOINT_PROPERTY_COLOR3,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.adc.tick_color",
	BRLOBOL_ENDPOINT_PROPERTY_COLOR3,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.adc.line_width",
	BRLOBOL_ENDPOINT_PROPERTY_UINT,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 1.0, 1048576.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.center_dot.visible",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.grid.visible",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.grid.adaptive",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.grid.snap",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.grid.anchor.x",
	BRLOBOL_ENDPOINT_PROPERTY_DOUBLE,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, -1.0e12, 1.0e12, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.grid.anchor.y",
	BRLOBOL_ENDPOINT_PROPERTY_DOUBLE,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, -1.0e12, 1.0e12, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.grid.anchor.z",
	BRLOBOL_ENDPOINT_PROPERTY_DOUBLE,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, -1.0e12, 1.0e12, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.grid.resolution.horizontal",
	BRLOBOL_ENDPOINT_PROPERTY_DOUBLE,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0e12, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.grid.resolution.vertical",
	BRLOBOL_ENDPOINT_PROPERTY_DOUBLE,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0e12, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.grid.major.horizontal",
	BRLOBOL_ENDPOINT_PROPERTY_UINT,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 2147483647.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.grid.major.vertical",
	BRLOBOL_ENDPOINT_PROPERTY_UINT,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 2147483647.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.grid.color",
	BRLOBOL_ENDPOINT_PROPERTY_COLOR3,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.model_axes.visible",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    FACEPLATE_AXES_STYLE_PROPERTIES("model_axes"),
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.scale.visible",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.view_axes.visible",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    FACEPLATE_AXES_STYLE_PROPERTIES("view_axes"),
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.params.visible",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.params.size",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.params.center",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.params.azimuth",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.params.elevation",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.params.twist",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.params.fps",
	BRLOBOL_ENDPOINT_PROPERTY_BOOL,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.params.font_size",
	BRLOBOL_ENDPOINT_PROPERTY_UINT,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 5.0, 96.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.center_dot.font_size",
	BRLOBOL_ENDPOINT_PROPERTY_UINT,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 5.0, 96.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.scale.font_size",
	BRLOBOL_ENDPOINT_PROPERTY_UINT,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 5.0, 96.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.params.color",
	BRLOBOL_ENDPOINT_PROPERTY_COLOR3,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.center_dot.color",
	BRLOBOL_ENDPOINT_PROPERTY_COLOR3,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "view.faceplate.scale.color",
	BRLOBOL_ENDPOINT_PROPERTY_COLOR3,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(brlobol_endpoint_property_desc), "composition.framebuffer.mode",
	BRLOBOL_ENDPOINT_PROPERTY_ENUM,
	BRLOBOL_ENDPOINT_PROPERTY_READ | BRLOBOL_ENDPOINT_PROPERTY_WRITE,
	0, 0.0, 0.0, "off,overlay,underlay,interlay"}
};

#undef FACEPLATE_AXES_STYLE_PROPERTIES

static const brlobol_endpoint_property_desc *
endpoint_property(const char *name)
{
    if (!name || !name[0])
	return NULL;
    for (const brlobol_endpoint_property_desc &property : endpoint_properties) {
	if (bu_strcmp(property.name, name) == 0)
	    return &property;
    }
    return NULL;
}

static int
endpoint_property_host_supported(const brlobol_display_endpoint_t *endpoint,
	const brlobol_endpoint_property_desc *property)
{
    if (!endpoint || !property)
	return 0;
    if (property->required_host_capabilities &&
	(brlobol_display_endpoint_host_capabilities(endpoint) &
	 property->required_host_capabilities) !=
	 property->required_host_capabilities)
	return 0;
    if ((bu_strcmp(property->name, "endpoint.title") == 0 ||
	 bu_strcmp(property->name, "endpoint.visible") == 0) &&
	endpoint->host_mode != BRLOBOL_HOST_MODE_TOPLEVEL)
	return 0;
    return 1;
}

static const char *
render_engine_name(enum brlobol_render_engine engine)
{
    static const char *names[] = {
	"auto", "hw", "sw", "rt", "none", "diagnostic"
    };
    return valid_engine(engine) ? names[static_cast<int>(engine)] : "auto";
}

static int
render_engine_from_name(const char *name, enum brlobol_render_engine *engine)
{
    if (!name || !engine)
	return 0;
    for (int i = BRLOBOL_RENDER_ENGINE_AUTO;
	    i <= BRLOBOL_RENDER_ENGINE_DIAGNOSTIC; i++) {
	if (bu_strcmp(name, render_engine_name(
		static_cast<enum brlobol_render_engine>(i))) == 0) {
	    *engine = static_cast<enum brlobol_render_engine>(i);
	    return 1;
	}
    }
    return 0;
}

static bool
valid_engine(enum brlobol_render_engine engine)
{
    return engine >= BRLOBOL_RENDER_ENGINE_AUTO &&
	engine <= BRLOBOL_RENDER_ENGINE_DIAGNOSTIC;
}

static bool
engine_host_compatible(enum brlobol_render_engine engine,
	uint64_t capabilities)
{
    switch (engine) {
	case BRLOBOL_RENDER_ENGINE_HW:
	    return (capabilities & BRLOBOL_HOST_CAP_SYSTEM_GL) != 0;
	case BRLOBOL_RENDER_ENGINE_SW:
	    return (capabilities & BRLOBOL_HOST_CAP_PIXEL_PRESENT) != 0;
	case BRLOBOL_RENDER_ENGINE_RT:
	    return false;
	default:
	    return true;
    }
}

extern "C" brlobol_display_endpoint_t *
brlobol_display_endpoint_create(void *controller, unsigned int flags)
{
    brlobol_init(NULL);

    brlobol_display_endpoint_t *endpoint =
	new (std::nothrow) brlobol_display_endpoint_t;
    if (!endpoint)
	return NULL;

    endpoint->controller = static_cast<BRLObolViewController *>(controller);
    endpoint->owns_controller =
	(flags & BRLOBOL_ENDPOINT_OWN_CONTROLLER) != 0;
    if (!endpoint->controller) {
	endpoint->controller = new (std::nothrow)
	    BRLObolViewController(new SoBRLSceneGroup);
	endpoint->owns_controller = true;
    }

    if (!endpoint->controller) {
	delete endpoint;
	return NULL;
    }

    const SbVec2s viewport =
	endpoint->controller->getViewportRegion().getWindowSize();
    endpoint->width = viewport[0] > 0 ?
	static_cast<unsigned int>(viewport[0]) : 1;
    endpoint->height = viewport[1] > 0 ?
	static_cast<unsigned int>(viewport[1]) : 1;

    return endpoint;
}

extern "C" void
brlobol_display_endpoint_host_detach(brlobol_display_endpoint_t *endpoint)
{
    if (!endpoint)
	return;

	/* The factory callback is endpoint-owned.  Once it is cleared, a second
	 * detach must not touch a controller that GED may already have released. */
    if (endpoint->factory_frame_callback_bound && endpoint->controller) {
	endpoint->controller->clearFrameRequestCallback(endpoint);
    }
    endpoint->factory_frame_callback_bound = false;

    if (endpoint->factory) {
	brlobol_host_factory_instance_destroy(endpoint->factory,
	    endpoint->factory_instance);
	brlobol_host_factory_release(endpoint->factory);
	endpoint->factory = NULL;
	endpoint->factory_instance = NULL;
    }

    endpoint->host_mode = BRLOBOL_HOST_MODE_HEADLESS;
    endpoint->title.clear();
    endpoint->visible = false;
    endpoint->vsync = true;

    if (!endpoint->host)
	return;

    BRLObolWindowHost *host = endpoint->host;
    const bool owns_host = endpoint->owns_host;
    endpoint->host = NULL;
    endpoint->owns_host = false;

    host->attachController(NULL, FALSE);
    if (owns_host)
	delete host;
}

extern "C" void
brlobol_display_endpoint_destroy(brlobol_display_endpoint_t *endpoint)
{
    if (!endpoint)
	return;

    brlobol_display_endpoint_host_detach(endpoint);
    if (endpoint->owns_controller)
	delete endpoint->controller;
    endpoint->controller = NULL;
    delete endpoint;
}

extern "C" void *
brlobol_display_endpoint_controller(const brlobol_display_endpoint_t *endpoint)
{
    return endpoint ? endpoint->controller : NULL;
}

extern "C" int
brlobol_display_endpoint_view_sync(brlobol_display_endpoint_t *endpoint,
	const void *view_ctx)
{
    return endpoint && endpoint->controller && view_ctx &&
	endpoint->controller->syncCameraFromViewContext(view_ctx) ? 1 : 0;
}

extern "C" int
brlobol_display_endpoint_host_bind(brlobol_display_endpoint_t *endpoint,
	void *host_ptr, unsigned int flags)
{
    BRLObolWindowHost *host = static_cast<BRLObolWindowHost *>(host_ptr);
    if (!endpoint || !endpoint->controller || !host)
	return 0;
    if (endpoint->engine == BRLOBOL_RENDER_ENGINE_HW ||
	endpoint->engine == BRLOBOL_RENDER_ENGINE_SW ||
	endpoint->engine == BRLOBOL_RENDER_ENGINE_RT)
	return 0;

    if (endpoint->host == host) {
	endpoint->owns_host = endpoint->owns_host ||
	    (flags & BRLOBOL_ENDPOINT_OWN_HOST) != 0;
	if (host->getController() != endpoint->controller)
	    host->attachController(endpoint->controller, FALSE);
	return 1;
    }

    brlobol_display_endpoint_host_detach(endpoint);
    host->attachController(endpoint->controller, FALSE);
    endpoint->host = host;
    endpoint->owns_host = (flags & BRLOBOL_ENDPOINT_OWN_HOST) != 0;
    return 1;
}

extern "C" void *
brlobol_display_endpoint_host(const brlobol_display_endpoint_t *endpoint)
{
    if (!endpoint)
	return NULL;
    return endpoint->factory_instance ? endpoint->factory_instance :
	endpoint->host;
}

extern "C" int
brlobol_display_endpoint_host_open(brlobol_display_endpoint_t *endpoint,
	const char *factory_name, const struct brlobol_host_desc *desc)
{
    if (!endpoint || !endpoint->controller)
	return 0;

    if (endpoint->engine == BRLOBOL_RENDER_ENGINE_RT)
	return 0;

    struct brlobol_host_desc required_desc = {};
    if (desc) {
	const size_t copy_size = desc->struct_size < sizeof(required_desc) ?
	    desc->struct_size : sizeof(required_desc);
	std::memcpy(&required_desc, desc, copy_size);
	required_desc.struct_size = sizeof(required_desc);
    } else {
	required_desc.struct_size = sizeof(required_desc);
	required_desc.mode = BRLOBOL_HOST_MODE_HEADLESS;
	required_desc.width = endpoint->width;
	required_desc.height = endpoint->height;
	required_desc.device_pixel_ratio = endpoint->device_pixel_ratio;
	required_desc.visible = 0;
    }
    if (endpoint->engine == BRLOBOL_RENDER_ENGINE_HW)
	required_desc.required_capabilities |= BRLOBOL_HOST_CAP_SYSTEM_GL;
    else if (endpoint->engine == BRLOBOL_RENDER_ENGINE_SW)
	required_desc.required_capabilities |= BRLOBOL_HOST_CAP_PIXEL_PRESENT;
    if (required_desc.vsync != BRLOBOL_HOST_VSYNC_AUTO)
	required_desc.required_capabilities |= BRLOBOL_HOST_CAP_PRESENT_VSYNC;
    required_desc.input_dispatch = endpoint_input_dispatch;
    required_desc.input_dispatch_data = endpoint;

    brlobol_host_factory_token_t *factory =
	brlobol_host_factory_acquire(factory_name, &required_desc);
    if (!factory)
	return 0;
    if (!engine_host_compatible(endpoint->engine,
	    brlobol_host_factory_capabilities(factory))) {
	brlobol_host_factory_release(factory);
	return 0;
    }

    void *instance = NULL;
	if (!brlobol_host_factory_instance_create(factory, &required_desc,
	endpoint->controller, &instance)) {
	brlobol_host_factory_release(factory);
	return 0;
    }

    brlobol_display_endpoint_host_detach(endpoint);
    endpoint->factory = factory;
    endpoint->factory_instance = instance;
    endpoint->host_mode = required_desc.mode;
    endpoint->width = required_desc.width ? required_desc.width :
	endpoint->width;
    endpoint->height = required_desc.height ? required_desc.height :
	endpoint->height;
    endpoint->device_pixel_ratio = required_desc.device_pixel_ratio > 0.0 ?
	required_desc.device_pixel_ratio : endpoint->device_pixel_ratio;
    endpoint->title = required_desc.title ? required_desc.title : "";
    endpoint->visible = required_desc.visible != 0;
    endpoint->vsync = required_desc.vsync != BRLOBOL_HOST_VSYNC_OFF;
    endpoint->controller->setFrameRequestCallback(endpoint_frame_requested,
	endpoint);
	endpoint->factory_frame_callback_bound = true;
    if (endpoint->controller->isRenderRequested() ||
	endpoint->controller->hasProgressiveWorkPending())
	endpoint_frame_requested(endpoint, "endpoint-host-open");
    return 1;
}

extern "C" const char *
brlobol_display_endpoint_host_factory_name(
	const brlobol_display_endpoint_t *endpoint)
{
    return endpoint && endpoint->factory ?
	brlobol_host_factory_name(endpoint->factory) : NULL;
}

extern "C" uint64_t
brlobol_display_endpoint_host_capabilities(
	const brlobol_display_endpoint_t *endpoint)
{
    return endpoint && endpoint->factory ?
	brlobol_host_factory_capabilities(endpoint->factory) : 0;
}

extern "C" int
brlobol_display_endpoint_request_frame(brlobol_display_endpoint_t *endpoint,
	const char *reason)
{
    if (!endpoint || !endpoint->controller)
	return 0;
    endpoint->controller->requestRender(reason);
    if (!endpoint->factory)
	return 1;
    return brlobol_host_factory_instance_request_frame(endpoint->factory,
	endpoint->factory_instance, reason);
}

extern "C" int
brlobol_display_endpoint_resize(brlobol_display_endpoint_t *endpoint,
	unsigned int width, unsigned int height, double device_pixel_ratio)
{
    if (!endpoint || !endpoint->controller || !width || !height ||
	device_pixel_ratio <= 0.0)
	return 0;
    if (endpoint->factory && !brlobol_host_factory_instance_resize(
	endpoint->factory, endpoint->factory_instance, width, height,
	device_pixel_ratio))
	return 0;
    endpoint->controller->setViewportSize(width, height);
    endpoint->width = width;
    endpoint->height = height;
    endpoint->device_pixel_ratio = device_pixel_ratio;
    return 1;
}

extern "C" int
brlobol_display_endpoint_input_profile_set(
	brlobol_display_endpoint_t *endpoint,
	const BRLObolInputProfile *profile)
{
    if (!endpoint)
	return 0;
    endpoint->input.setProfile(profile ? profile :
	brlobol_input_default_view_profile());
    return 1;
}

extern "C" int
brlobol_display_endpoint_input_action_handler_set(
	brlobol_display_endpoint_t *endpoint, BRLObolInputActionHandler handler,
	void *user_data)
{
    if (!endpoint)
	return 0;
    endpoint->input.setActionHandler(handler, user_data);
    return 1;
}

extern "C" int
brlobol_display_endpoint_input_action_handler_clear_if(
	brlobol_display_endpoint_t *endpoint, BRLObolInputActionHandler handler,
	void *user_data)
{
    return endpoint ? endpoint->input.clearActionHandlerIf(handler,
	user_data) : 0;
}

extern "C" int
brlobol_display_endpoint_input_dispatch(brlobol_display_endpoint_t *endpoint,
	const BRLObolInputEvent *event)
{
    return endpoint ? endpoint->input.dispatch(event) : -1;
}

extern "C" int
brlobol_display_endpoint_capture(brlobol_display_endpoint_t *endpoint,
	unsigned char **pixels, size_t *size, unsigned int *width,
	unsigned int *height, unsigned int *components)
{
    return brlobol_display_endpoint_capture_plane(endpoint,
	BRLOBOL_CAPTURE_COMPOSITE, pixels, size, width, height, components);
}

extern "C" int
brlobol_display_endpoint_capture_plane(brlobol_display_endpoint_t *endpoint,
	enum brlobol_capture_plane plane, unsigned char **pixels, size_t *size,
	unsigned int *width, unsigned int *height, unsigned int *components)
{
    if (pixels)
	*pixels = NULL;
    if (size)
	*size = 0;
    if (width)
	*width = 0;
    if (height)
	*height = 0;
    if (components)
	*components = 0;
    if (!endpoint || !endpoint->controller || !pixels || !size || !width ||
	!height || !components)
	return 0;

    if (plane == BRLOBOL_CAPTURE_FRAMEBUFFER) {
	return endpoint->framebuffer_capture_callback ?
	    endpoint->framebuffer_capture_callback(
		endpoint->framebuffer_capture_user_data, pixels, size, width,
		height, components) > 0 ? 1 : 0 : 0;
    }
    if (plane != BRLOBOL_CAPTURE_COMPOSITE)
	return 0;

    if (endpoint->factory)
	return brlobol_host_factory_instance_capture(endpoint->factory,
	    endpoint->factory_instance, pixels, size, width, height,
	    components) > 0 ? 1 : 0;

    if (endpoint->controller->renderToImage(pixels, 0, 0) != BRLCAD_OK ||
	!(*pixels))
	return 0;
    const SbVec2s viewport =
	endpoint->controller->getViewportRegion().getWindowSize();
    *width = viewport[0] > 0 ? (unsigned int)viewport[0] : 1;
    *height = viewport[1] > 0 ? (unsigned int)viewport[1] : 1;
    *components = 3;
    *size = (size_t)(*width) * (size_t)(*height) * (*components);
    return 1;
}

extern "C" int
brlobol_display_endpoint_framebuffer_capture_provider_set(
	brlobol_display_endpoint_t *endpoint,
	brlobol_endpoint_framebuffer_capture_callback callback,
	void *user_data)
{
    if (!endpoint || (callback && !user_data))
	return 0;
    if (!callback && endpoint->framebuffer_capture_user_data != user_data)
	return 0;
    endpoint->framebuffer_capture_callback = callback;
    endpoint->framebuffer_capture_user_data = callback ? user_data : NULL;
    return 1;
}

extern "C" int
brlobol_display_endpoint_render_engine_set(
	brlobol_display_endpoint_t *endpoint,
	enum brlobol_render_engine engine)
{
    if (!endpoint || !valid_engine(engine))
	return 0;
    if (engine == BRLOBOL_RENDER_ENGINE_RT)
	return 0;
    if (endpoint->factory && !engine_host_compatible(engine,
	brlobol_host_factory_capabilities(endpoint->factory)))
	return 0;
    if (endpoint->host && !endpoint->factory &&
	(engine == BRLOBOL_RENDER_ENGINE_HW ||
	 engine == BRLOBOL_RENDER_ENGINE_SW))
	return 0;
    endpoint->engine = engine;
    return 1;
}

extern "C" enum brlobol_render_engine
brlobol_display_endpoint_render_engine_get(
	const brlobol_display_endpoint_t *endpoint)
{
    return endpoint ? endpoint->engine : BRLOBOL_RENDER_ENGINE_AUTO;
}

extern "C" size_t
brlobol_display_endpoint_property_count(void)
{
    return sizeof(endpoint_properties) / sizeof(endpoint_properties[0]);
}

extern "C" int
brlobol_display_endpoint_property_descriptor(size_t index,
	struct brlobol_endpoint_property_desc *out)
{
    if (!out || out->struct_size < sizeof(*out) ||
	index >= brlobol_display_endpoint_property_count())
	return BRLOBOL_ENDPOINT_PROPERTY_INVALID;
    *out = endpoint_properties[index];
    return BRLOBOL_ENDPOINT_PROPERTY_OK;
}

static void
property_value_prepare(brlobol_endpoint_property_value *value,
	enum brlobol_endpoint_property_type type)
{
    const uint32_t struct_size = value->struct_size;
    std::memset(value, 0, sizeof(*value));
    value->struct_size = struct_size;
    value->type = type;
}

extern "C" int
brlobol_display_endpoint_property_get(
	const brlobol_display_endpoint_t *endpoint, const char *name,
	struct brlobol_endpoint_property_value *out)
{
    if (!endpoint || !endpoint->controller || !out ||
	out->struct_size < sizeof(*out))
	return BRLOBOL_ENDPOINT_PROPERTY_INVALID;
    const brlobol_endpoint_property_desc *property = endpoint_property(name);
    if (!property)
	return BRLOBOL_ENDPOINT_PROPERTY_UNKNOWN;
    if (!endpoint_property_host_supported(endpoint, property))
	return BRLOBOL_ENDPOINT_PROPERTY_UNSUPPORTED;
    property_value_prepare(out, property->type);

    if (bu_strcmp(name, "endpoint.width") == 0 ||
	bu_strcmp(name, "endpoint.height") == 0 ||
	bu_strcmp(name, "endpoint.device_pixel_ratio") == 0)
	endpoint_dimensions_refresh(
	    const_cast<brlobol_display_endpoint_t *>(endpoint));

    if (bu_strcmp(name, "endpoint.renderer") == 0) {
	out->string_value = render_engine_name(endpoint->engine);
    } else if (bu_strcmp(name, "endpoint.width") == 0) {
	out->uint_value = endpoint->width;
    } else if (bu_strcmp(name, "endpoint.height") == 0) {
	out->uint_value = endpoint->height;
    } else if (bu_strcmp(name, "endpoint.device_pixel_ratio") == 0) {
	out->double_value = endpoint->device_pixel_ratio;
    } else if (bu_strcmp(name, "endpoint.host") == 0) {
	const char *factory_name =
	    brlobol_display_endpoint_host_factory_name(endpoint);
	out->string_value = factory_name ? factory_name :
	    (endpoint->host ? "bound" : "");
    } else if (bu_strcmp(name, "endpoint.title") == 0) {
	out->string_value = endpoint->title.c_str();
    } else if (bu_strcmp(name, "endpoint.visible") == 0) {
	out->bool_value = endpoint->visible ? 1 : 0;
    } else if (bu_strcmp(name, "endpoint.vsync") == 0) {
	out->bool_value = endpoint->vsync ? 1 : 0;
    } else if (bu_strcmp(name, "controller.background.bottom") == 0 ||
	bu_strcmp(name, "controller.background.top") == 0) {
	const SbColor color = bu_strcmp(name,
	    "controller.background.bottom") == 0 ?
	    endpoint->controller->getBackgroundBottomColor() :
	    endpoint->controller->getBackgroundTopColor();
	out->color3[0] = color[0];
	out->color3[1] = color[1];
	out->color3[2] = color[2];
    } else if (bu_strcmp(name, "controller.software_wire") == 0) {
	switch (endpoint->controller->getSoftwareWireMode()) {
	    case BRLObolViewController::SOFTWARE_WIRE_QUALITY:
		out->string_value = "quality";
		break;
	    case BRLObolViewController::SOFTWARE_WIRE_FAST:
		out->string_value = "fast";
		break;
	    default:
		out->string_value = "auto";
		break;
	}
    } else if (bu_strcmp(name, "renderer.depth_test") == 0) {
	out->bool_value = endpoint->controller->isDepthTestEnabled() ? 1 : 0;
    } else if (bu_strcmp(name, "renderer.lighting") == 0) {
	out->bool_value = endpoint->controller->isLightingEnabled() ? 1 : 0;
    } else if (bu_strcmp(name, "renderer.transparency") == 0) {
	out->bool_value = endpoint->controller->isTransparencyEnabled() ? 1 : 0;
    } else if (bu_strcmp(name, "renderer.antialiasing") == 0) {
	out->bool_value = endpoint->controller->isAntialiasingEnabled() ? 1 : 0;
    } else if (bu_strcmp(name, "renderer.clip.minimum") == 0 ||
	bu_strcmp(name, "renderer.clip.maximum") == 0) {
	double minimum = 0.0;
	double maximum = 0.0;
	endpoint->controller->getClipBounds(minimum, maximum);
	out->double_value = bu_strcmp(name, "renderer.clip.minimum") == 0 ?
	    minimum : maximum;
    } else if (bu_strcmp(name, "renderer.depth_cue") == 0) {
	out->bool_value = endpoint->controller->isDepthCueEnabled() ? 1 : 0;
    } else if (endpoint->property_get_callback) {
	return endpoint->property_get_callback(endpoint->property_user_data,
	    name, out);
    } else {
	return BRLOBOL_ENDPOINT_PROPERTY_UNSUPPORTED;
    }
    return BRLOBOL_ENDPOINT_PROPERTY_OK;
}

static int
valid_color3(const double color[3])
{
    return color && std::isfinite(color[0]) && std::isfinite(color[1]) &&
	std::isfinite(color[2]) && color[0] >= 0.0 && color[0] <= 1.0 &&
	color[1] >= 0.0 && color[1] <= 1.0 && color[2] >= 0.0 &&
	color[2] <= 1.0;
}

extern "C" int
brlobol_display_endpoint_property_set(
	brlobol_display_endpoint_t *endpoint, const char *name,
	const struct brlobol_endpoint_property_value *value)
{
    if (!endpoint || !endpoint->controller || !value ||
	value->struct_size < sizeof(*value))
	return BRLOBOL_ENDPOINT_PROPERTY_INVALID;
    const brlobol_endpoint_property_desc *property = endpoint_property(name);
    if (!property)
	return BRLOBOL_ENDPOINT_PROPERTY_UNKNOWN;
    if (!(property->access & BRLOBOL_ENDPOINT_PROPERTY_WRITE))
	return BRLOBOL_ENDPOINT_PROPERTY_READ_ONLY;
    if (value->type != property->type)
	return BRLOBOL_ENDPOINT_PROPERTY_INVALID;
    if (property->type == BRLOBOL_ENDPOINT_PROPERTY_BOOL &&
	value->bool_value != 0 && value->bool_value != 1)
	return BRLOBOL_ENDPOINT_PROPERTY_INVALID;
    if (property->type == BRLOBOL_ENDPOINT_PROPERTY_DOUBLE &&
	(!std::isfinite(value->double_value) ||
	 (property->minimum < property->maximum &&
	  (value->double_value < property->minimum ||
	   value->double_value > property->maximum))))
	return BRLOBOL_ENDPOINT_PROPERTY_INVALID;
    if (property->type == BRLOBOL_ENDPOINT_PROPERTY_UINT &&
	property->minimum < property->maximum &&
	(value->uint_value < static_cast<uint64_t>(property->minimum) ||
	 value->uint_value > static_cast<uint64_t>(property->maximum)))
	return BRLOBOL_ENDPOINT_PROPERTY_INVALID;
    if (property->type == BRLOBOL_ENDPOINT_PROPERTY_COLOR3 &&
	!valid_color3(value->color3))
	return BRLOBOL_ENDPOINT_PROPERTY_INVALID;
    if (!endpoint_property_host_supported(endpoint, property))
	return BRLOBOL_ENDPOINT_PROPERTY_UNSUPPORTED;

    if (bu_strcmp(name, "endpoint.renderer") == 0) {
	enum brlobol_render_engine engine = BRLOBOL_RENDER_ENGINE_AUTO;
	if (!render_engine_from_name(value->string_value, &engine))
	    return BRLOBOL_ENDPOINT_PROPERTY_INVALID;
	if (engine == BRLOBOL_RENDER_ENGINE_RT)
	    return BRLOBOL_ENDPOINT_PROPERTY_UNSUPPORTED;
	if (!brlobol_display_endpoint_render_engine_set(endpoint, engine))
	    return BRLOBOL_ENDPOINT_PROPERTY_UNSUPPORTED;
    } else if (bu_strcmp(name, "endpoint.width") == 0) {
	if (value->uint_value < 1 || value->uint_value > 32767 ||
	    !brlobol_display_endpoint_resize(endpoint,
		static_cast<unsigned int>(value->uint_value), endpoint->height,
		endpoint->device_pixel_ratio))
	    return BRLOBOL_ENDPOINT_PROPERTY_INVALID;
    } else if (bu_strcmp(name, "endpoint.height") == 0) {
	if (value->uint_value < 1 || value->uint_value > 32767 ||
	    !brlobol_display_endpoint_resize(endpoint, endpoint->width,
		static_cast<unsigned int>(value->uint_value),
		endpoint->device_pixel_ratio))
	    return BRLOBOL_ENDPOINT_PROPERTY_INVALID;
    } else if (bu_strcmp(name, "endpoint.device_pixel_ratio") == 0) {
	if (!std::isfinite(value->double_value) ||
	    value->double_value < 0.01 || value->double_value > 64.0 ||
	    !brlobol_display_endpoint_resize(endpoint, endpoint->width,
		endpoint->height, value->double_value))
	    return BRLOBOL_ENDPOINT_PROPERTY_INVALID;
    } else if (bu_strcmp(name, "endpoint.title") == 0) {
	if (!value->string_value)
	    return BRLOBOL_ENDPOINT_PROPERTY_INVALID;
	if (!brlobol_host_factory_instance_set_title(endpoint->factory,
	    endpoint->factory_instance, value->string_value))
	    return BRLOBOL_ENDPOINT_PROPERTY_UNSUPPORTED;
	endpoint->title = value->string_value;
    } else if (bu_strcmp(name, "endpoint.visible") == 0) {
	if (!brlobol_host_factory_instance_set_visible(endpoint->factory,
	    endpoint->factory_instance, value->bool_value))
	    return BRLOBOL_ENDPOINT_PROPERTY_UNSUPPORTED;
	endpoint->visible = value->bool_value != 0;
    } else if (bu_strcmp(name, "endpoint.vsync") == 0) {
	if (!brlobol_host_factory_instance_set_vsync(endpoint->factory,
	    endpoint->factory_instance, value->bool_value))
	    return BRLOBOL_ENDPOINT_PROPERTY_UNSUPPORTED;
	endpoint->vsync = value->bool_value != 0;
    } else if (bu_strcmp(name, "controller.background.bottom") == 0 ||
	bu_strcmp(name, "controller.background.top") == 0) {
	if (!valid_color3(value->color3))
	    return BRLOBOL_ENDPOINT_PROPERTY_INVALID;
	const SbColor updated(static_cast<float>(value->color3[0]),
	    static_cast<float>(value->color3[1]),
	    static_cast<float>(value->color3[2]));
	const SbColor bottom = bu_strcmp(name,
	    "controller.background.bottom") == 0 ? updated :
	    endpoint->controller->getBackgroundBottomColor();
	const SbColor top = bu_strcmp(name,
	    "controller.background.top") == 0 ? updated :
	    endpoint->controller->getBackgroundTopColor();
	endpoint->controller->setBackgroundColors(bottom, top);
    } else if (bu_strcmp(name, "controller.software_wire") == 0) {
	BRLObolViewController::SoftwareWireMode mode =
	    BRLObolViewController::SOFTWARE_WIRE_AUTO;
	if (!value->string_value)
	    return BRLOBOL_ENDPOINT_PROPERTY_INVALID;
	if (bu_strcmp(value->string_value, "quality") == 0)
	    mode = BRLObolViewController::SOFTWARE_WIRE_QUALITY;
	else if (bu_strcmp(value->string_value, "fast") == 0)
	    mode = BRLObolViewController::SOFTWARE_WIRE_FAST;
	else if (bu_strcmp(value->string_value, "auto") != 0)
	    return BRLOBOL_ENDPOINT_PROPERTY_INVALID;
	endpoint->controller->setSoftwareWireMode(mode);
    } else if (bu_strcmp(name, "renderer.depth_test") == 0) {
	endpoint->controller->setDepthTestEnabled(
	    value->bool_value ? TRUE : FALSE);
    } else if (bu_strcmp(name, "renderer.lighting") == 0) {
	endpoint->controller->setLightingEnabled(
	    value->bool_value ? TRUE : FALSE);
    } else if (bu_strcmp(name, "renderer.transparency") == 0) {
	endpoint->controller->setTransparencyEnabled(
	    value->bool_value ? TRUE : FALSE);
    } else if (bu_strcmp(name, "renderer.antialiasing") == 0) {
	endpoint->controller->setAntialiasingEnabled(
	    value->bool_value ? TRUE : FALSE);
    } else if (bu_strcmp(name, "renderer.clip.minimum") == 0 ||
	bu_strcmp(name, "renderer.clip.maximum") == 0) {
	if (!std::isfinite(value->double_value) ||
	    value->double_value < -1.0e12 || value->double_value > 1.0e12)
	    return BRLOBOL_ENDPOINT_PROPERTY_INVALID;
	double minimum = 0.0;
	double maximum = 0.0;
	endpoint->controller->getClipBounds(minimum, maximum);
	if (bu_strcmp(name, "renderer.clip.minimum") == 0)
	    minimum = value->double_value;
	else
	    maximum = value->double_value;
	if (!endpoint->controller->setClipBounds(minimum, maximum))
	    return BRLOBOL_ENDPOINT_PROPERTY_INVALID;
    } else if (bu_strcmp(name, "renderer.depth_cue") == 0) {
	endpoint->controller->setDepthCueEnabled(
	    value->bool_value ? TRUE : FALSE);
    } else if (endpoint->property_set_callback) {
	const int ret = endpoint->property_set_callback(
	    endpoint->property_user_data, name, value);
	if (ret == BRLOBOL_ENDPOINT_PROPERTY_OK)
	    endpoint->controller->requestRender("external property");
	return ret;
    } else {
	return BRLOBOL_ENDPOINT_PROPERTY_UNSUPPORTED;
    }

    return BRLOBOL_ENDPOINT_PROPERTY_OK;
}

extern "C" int
brlobol_display_endpoint_property_provider_set(
	brlobol_display_endpoint_t *endpoint,
	brlobol_endpoint_property_get_callback get_callback,
	brlobol_endpoint_property_set_callback set_callback,
	void *user_data)
{
    if (!endpoint || ((get_callback || set_callback) && !user_data))
	return 0;
    endpoint->property_get_callback = get_callback;
    endpoint->property_set_callback = set_callback;
    endpoint->property_user_data = user_data;
    return 1;
}

BRLObolDisplayEndpoint::BRLObolDisplayEndpoint(
	BRLObolViewController *controller, bool takeOwnership) :
    endpoint(brlobol_display_endpoint_create(controller,
	takeOwnership ? BRLOBOL_ENDPOINT_OWN_CONTROLLER : 0))
{
}

BRLObolDisplayEndpoint::~BRLObolDisplayEndpoint(void)
{
    brlobol_display_endpoint_destroy(this->endpoint);
    this->endpoint = NULL;
}

bool
BRLObolDisplayEndpoint::isValid(void) const
{
    return this->endpoint != NULL;
}

BRLObolViewController *
BRLObolDisplayEndpoint::controller(void) const
{
    return static_cast<BRLObolViewController *>(
	brlobol_display_endpoint_controller(this->endpoint));
}

int
BRLObolDisplayEndpoint::bindHost(BRLObolWindowHost *host, bool takeOwnership)
{
    return brlobol_display_endpoint_host_bind(this->endpoint, host,
	takeOwnership ? BRLOBOL_ENDPOINT_OWN_HOST : 0);
}

void
BRLObolDisplayEndpoint::detachHost(void)
{
    brlobol_display_endpoint_host_detach(this->endpoint);
}

BRLObolWindowHost *
BRLObolDisplayEndpoint::host(void) const
{
    return static_cast<BRLObolWindowHost *>(
	brlobol_display_endpoint_host(this->endpoint));
}

int
BRLObolDisplayEndpoint::setRenderEngine(enum brlobol_render_engine engine)
{
    return brlobol_display_endpoint_render_engine_set(this->endpoint, engine);
}

enum brlobol_render_engine
BRLObolDisplayEndpoint::renderEngine(void) const
{
    return brlobol_display_endpoint_render_engine_get(this->endpoint);
}

int
BRLObolDisplayEndpoint::propertyGet(const char *name,
	struct brlobol_endpoint_property_value *value) const
{
    return brlobol_display_endpoint_property_get(this->endpoint, name, value);
}

int
BRLObolDisplayEndpoint::propertySet(const char *name,
	const struct brlobol_endpoint_property_value *value)
{
    return brlobol_display_endpoint_property_set(this->endpoint, name, value);
}

brlobol_display_endpoint_t *
BRLObolDisplayEndpoint::release(void)
{
    brlobol_display_endpoint_t *released = this->endpoint;
    this->endpoint = NULL;
    return released;
}

brlobol_display_endpoint_t *
BRLObolDisplayEndpoint::get(void) const
{
    return this->endpoint;
}
