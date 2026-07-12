/*              T E S T _ W I N D O W _ H O S T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol.h"

#include "bu/log.h"
#include "bu/malloc.h"

#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wshadow"
#endif
#include "imgstream/fb_compat.h"
#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

#include <Inventor/SbViewportRegion.h>
#include <Inventor/nodes/SoClipPlane.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoTexture2.h>

#include <math.h>
#include <string.h>

#define CHECK(_expr, _msg) do { \
    if (!(_expr)) { \
	bu_log("FAIL: %s\n", _msg); \
	return 1; \
    } \
} while (0)

static int
float_equal(float a, float b)
{
    return fabs((double)a - (double)b) <= 1.0e-6;
}

static SoClipPlane *
find_clip_plane(SoNode *node, const char *name)
{
    if (!node || !name)
	return NULL;
    if (node->isOfType(SoClipPlane::getClassTypeId()) &&
	strcmp(node->getName().getString(), name) == 0)
	return static_cast<SoClipPlane *>(node);
    if (!node->isOfType(SoGroup::getClassTypeId()))
	return NULL;
    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoClipPlane *plane = find_clip_plane(group->getChild(i), name);
	if (plane)
	    return plane;
    }
    return NULL;
}

static void
property_value_init(struct brlobol_endpoint_property_value *value)
{
    memset(value, 0, sizeof(*value));
    value->struct_size = sizeof(*value);
    value->type = BRLOBOL_ENDPOINT_PROPERTY_BOOL;
}

static int
test_property_provider_get(void *data, const char *name,
	struct brlobol_endpoint_property_value *value)
{
    if (!data || !name || !value || strcmp(name, "view.zclip") != 0)
	return BRLOBOL_ENDPOINT_PROPERTY_UNSUPPORTED;
    value->bool_value = *(int *)data;
    return BRLOBOL_ENDPOINT_PROPERTY_OK;
}

static int
test_property_provider_set(void *data, const char *name,
	const struct brlobol_endpoint_property_value *value)
{
    if (!data || !name || !value || strcmp(name, "view.zclip") != 0)
	return BRLOBOL_ENDPOINT_PROPERTY_UNSUPPORTED;
    *(int *)data = value->bool_value ? 1 : 0;
    return BRLOBOL_ENDPOINT_PROPERTY_OK;
}

static int
test_window_host_contract(void)
{
    BRLObolWindowHost host;
    BRLObolWindowDesc desc;

    desc.mode = BRLOBOL_WINDOW_HEADLESS;
    desc.backend = BRLOBOL_WINDOW_BACKEND_OFFSCREEN;
    desc.width = 6;
    desc.height = 4;
    desc.title = "window-host-test";
    desc.display = "";
    desc.nativeIdHint = "";
    desc.visible = FALSE;

    CHECK(host.open(&desc) == 0, "window host opens neutral controller");
    CHECK(host.isOpen(), "window host records open state");
    CHECK(host.getController() != NULL, "window host owns a view controller");
    CHECK(host.getController()->getSceneRoot() != NULL, "window host creates scene root");
    CHECK(host.getController()->getSceneRoot()->isOfType(SoGroup::getClassTypeId()),
	  "window host scene root is a group");
    CHECK(host.getController()->getViewportRegion().getWindowSize()[0] == 6 &&
	  host.getController()->getViewportRegion().getWindowSize()[1] == 4,
	  "window host applies requested viewport size");

    BRLObolInputBinding binding;
    binding.eventType = BRLOBOL_INPUT_KEY;
    binding.key = 'f';
    binding.button = 0;
    binding.modifiers = 0;
    binding.action = BRLOBOL_ACTION_VIEW_FRONT;

    BRLObolInputProfile profile;
    profile.name = "test";
    profile.bindings = &binding;
    profile.bindingCount = 1;

    BRLObolInputEvent event;
    memset(&event, 0, sizeof(event));
    event.type = BRLOBOL_INPUT_KEY;
    event.key = 'f';
    CHECK(host.handleInputEvent(&event, &profile) == 1,
	  "window host applies semantic input profile action");
    CHECK(host.getController()->isRenderRequested(),
	  "semantic input action requests render");
    host.getController()->clearRenderRequest();

    host.close();
    CHECK(!host.isOpen(), "window host closes");
    return 0;
}

class CountingWindowHost : public BRLObolWindowHost {
public:
    ~CountingWindowHost(void) override
    {
	destroyed++;
    }

    static int destroyed;
};

int CountingWindowHost::destroyed = 0;

struct FactoryTestState {
    FactoryTestState(void) :
	probe_result(1), open_result(1), creates(0), destroys(0), binds(0),
	detaches(0), opens(0), closes(0), frames(0), resizes(0)
    {
    }

    int probe_result;
    int open_result;
    int creates;
    int destroys;
    int binds;
    int detaches;
    int opens;
    int closes;
    int frames;
    int resizes;
    int captures = 0;
    int dimension_queries = 0;
    int framebuffer_captures = 0;
};

struct FactoryTestInstance {
    FactoryTestState *state;
    void *controller;
};

static int
factory_test_probe(const struct brlobol_host_desc *UNUSED(desc), void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    return state->probe_result;
}

static void *
factory_test_create(const struct brlobol_host_desc *UNUSED(desc), void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    FactoryTestInstance *instance = new FactoryTestInstance;
    instance->state = state;
    instance->controller = NULL;
    state->creates++;
    return instance;
}

static void
factory_test_destroy(void *instance_ptr, void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    state->destroys++;
    delete static_cast<FactoryTestInstance *>(instance_ptr);
}

static int
factory_test_bind(void *instance_ptr, void *controller, void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    FactoryTestInstance *instance =
	static_cast<FactoryTestInstance *>(instance_ptr);
    instance->controller = controller;
    if (controller)
	state->binds++;
    else
	state->detaches++;
    return 1;
}

static int
factory_test_open(void *UNUSED(instance),
	const struct brlobol_host_desc *UNUSED(desc), void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    state->opens++;
    return state->open_result;
}

static void
factory_test_close(void *UNUSED(instance), void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    state->closes++;
}

static int
factory_test_frame(void *UNUSED(instance), const char *UNUSED(reason),
	void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    state->frames++;
    return 1;
}

static int
factory_test_resize(void *UNUSED(instance), unsigned int UNUSED(width),
	unsigned int UNUSED(height), double UNUSED(device_pixel_ratio), void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    state->resizes++;
    return 1;
}

static int
factory_test_capture(void *UNUSED(instance), unsigned char **pixels,
	size_t *size, unsigned int *width, unsigned int *height,
	unsigned int *components, void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    *width = 2;
    *height = 2;
    *components = 3;
    *size = 12;
    *pixels = static_cast<unsigned char *>(bu_calloc(*size, 1,
	"factory test capture"));
    state->captures++;
    return 1;
}

static int
factory_test_dimensions(void *UNUSED(instance), unsigned int *width,
	unsigned int *height, double *device_pixel_ratio, void *data)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    *width = 13;
    *height = 9;
    *device_pixel_ratio = 2.0;
    state->dimension_queries++;
    return 1;
}

static int
factory_test_framebuffer_capture(void *data, unsigned char **pixels,
	size_t *size, unsigned int *width, unsigned int *height,
	unsigned int *components)
{
    FactoryTestState *state = static_cast<FactoryTestState *>(data);
    if (!state || !pixels || !size || !width || !height || !components)
	return 0;
    *width = 1;
    *height = 1;
    *components = 3;
    *size = 3;
    *pixels = static_cast<unsigned char *>(bu_malloc(*size,
	"factory framebuffer capture"));
    (*pixels)[0] = 17;
    (*pixels)[1] = 34;
    (*pixels)[2] = 51;
    state->framebuffer_captures++;
    return 1;
}

static struct brlobol_host_factory
factory_test_desc(const char *name, int priority, uint64_t capabilities,
	FactoryTestState *state)
{
    struct brlobol_host_factory factory;
    memset(&factory, 0, sizeof(factory));
    factory.abi_version = BRLOBOL_HOST_FACTORY_ABI_VERSION;
    factory.struct_size = sizeof(factory);
    factory.name = name;
    factory.priority = priority;
    factory.capabilities = capabilities;
    factory.user_data = state;
    factory.probe = factory_test_probe;
    factory.create = factory_test_create;
    factory.destroy = factory_test_destroy;
    factory.bind_controller = factory_test_bind;
    factory.open = factory_test_open;
    factory.close = factory_test_close;
    factory.request_frame = factory_test_frame;
    factory.resize = factory_test_resize;
    factory.capture = factory_test_capture;
    factory.dimensions = factory_test_dimensions;
    return factory;
}

static int
test_host_factory_contract(void)
{
    const size_t initial_count = brlobol_host_factory_registry_count();
    FactoryTestState low_state;
    FactoryTestState high_state;
    FactoryTestState failed_state;
    failed_state.open_result = 0;

    struct brlobol_host_factory invalid =
	factory_test_desc("endpoint-test-invalid", 0, 0, &low_state);
    invalid.abi_version++;
    CHECK(brlobol_host_factory_register(&invalid) == NULL,
	  "host factory rejects an unsupported ABI");

    struct brlobol_host_factory low = factory_test_desc(
	"endpoint-test-low", 10, BRLOBOL_HOST_CAP_PIXEL_PRESENT, &low_state);
    struct brlobol_host_factory high = factory_test_desc(
	"endpoint-test-high", 20,
	BRLOBOL_HOST_CAP_PIXEL_PRESENT | BRLOBOL_HOST_CAP_READBACK,
	&high_state);
    struct brlobol_host_factory failed = factory_test_desc(
	"endpoint-test-failed", 30, BRLOBOL_HOST_CAP_PIXEL_PRESENT,
	&failed_state);

    brlobol_host_factory_token_t *low_token =
	brlobol_host_factory_register(&low);
    brlobol_host_factory_token_t *high_token =
	brlobol_host_factory_register(&high);
    brlobol_host_factory_token_t *failed_token =
	brlobol_host_factory_register(&failed);
    CHECK(low_token && high_token && failed_token,
	  "host factories register copied descriptors");
    CHECK(brlobol_host_factory_register(&high) == NULL,
	  "host factory rejects duplicate stable names");
    CHECK(brlobol_host_factory_registry_count() == initial_count + 3,
	  "host factory registry reports registrations");
    int found_high_capabilities = 0;
    for (size_t i = 0; i < brlobol_host_factory_registry_count(); i++) {
	char name[64] = {0};
	if (brlobol_host_factory_registry_name(i, name, sizeof(name)) &&
	    strcmp(name, "endpoint-test-high") == 0) {
	    found_high_capabilities =
		brlobol_host_factory_registry_capabilities(i) ==
		(BRLOBOL_HOST_CAP_PIXEL_PRESENT |
		 BRLOBOL_HOST_CAP_READBACK);
	}
    }
    CHECK(found_high_capabilities,
	  "host factory registry exposes immutable capabilities");

    struct brlobol_host_desc desc;
    memset(&desc, 0, sizeof(desc));
    desc.struct_size = sizeof(desc);
    desc.mode = BRLOBOL_HOST_MODE_HEADLESS;
    desc.width = 8;
    desc.height = 6;
    desc.device_pixel_ratio = 1.0;
    desc.required_capabilities = BRLOBOL_HOST_CAP_READBACK;

    brlobol_display_endpoint_t *endpoint =
	brlobol_display_endpoint_create(NULL, 0);
    CHECK(endpoint != NULL, "factory test creates display endpoint");
    CHECK(brlobol_display_endpoint_host_open(endpoint, NULL, &desc),
	  "endpoint selects a compatible registered factory");
    CHECK(strcmp(brlobol_display_endpoint_host_factory_name(endpoint),
	  "endpoint-test-high") == 0,
	  "endpoint selection honors capabilities before priority");
    CHECK(brlobol_display_endpoint_host_capabilities(endpoint) ==
	  (BRLOBOL_HOST_CAP_PIXEL_PRESENT | BRLOBOL_HOST_CAP_READBACK),
	  "endpoint exposes its active host capabilities");
    FactoryTestInstance *instance = static_cast<FactoryTestInstance *>(
	brlobol_display_endpoint_host(endpoint));
    CHECK(instance && instance->controller ==
	  brlobol_display_endpoint_controller(endpoint),
	  "factory host instance binds the endpoint controller");
    CHECK(!brlobol_host_factory_unregister(high_token),
	  "live endpoint prevents host factory unregister");

    struct brlobol_endpoint_property_value host_dimension =
	BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    CHECK(brlobol_display_endpoint_property_get(endpoint, "endpoint.width",
	  &host_dimension) == BRLOBOL_ENDPOINT_PROPERTY_OK &&
	  host_dimension.uint_value == 13 && high_state.dimension_queries == 1,
	  "endpoint dimensions refresh from the active toolkit host");

    CHECK(brlobol_display_endpoint_request_frame(endpoint, "test-frame") &&
	  high_state.frames == 1,
	  "factory dispatches frame requests");
    CHECK(brlobol_display_endpoint_resize(endpoint, 10, 7, 1.5) &&
	  high_state.resizes == 1 &&
	  instance->controller == brlobol_display_endpoint_controller(endpoint),
	  "factory dispatches resize requests");
    unsigned char *capture = NULL;
    size_t capture_size = 0;
    unsigned int capture_width = 0;
    unsigned int capture_height = 0;
    unsigned int capture_components = 0;
    CHECK(brlobol_display_endpoint_capture(endpoint, &capture, &capture_size,
	  &capture_width, &capture_height, &capture_components) &&
	  capture && capture_size == 12 && capture_width == 2 &&
	  capture_height == 2 && capture_components == 3 &&
	  high_state.captures == 1,
	  "endpoint dispatches capture through a readback factory");
    bu_free(capture, "factory test capture");

    CHECK(!brlobol_display_endpoint_capture_plane(endpoint,
	  BRLOBOL_CAPTURE_FRAMEBUFFER, &capture, &capture_size,
	  &capture_width, &capture_height, &capture_components),
	  "framebuffer-plane capture fails without an explicit provider");
    CHECK(brlobol_display_endpoint_framebuffer_capture_provider_set(endpoint,
	  factory_test_framebuffer_capture, &high_state),
	  "endpoint accepts a per-instance framebuffer capture provider");
    CHECK(brlobol_display_endpoint_capture_plane(endpoint,
	  BRLOBOL_CAPTURE_FRAMEBUFFER, &capture, &capture_size,
	  &capture_width, &capture_height, &capture_components) && capture &&
	  capture_size == 3 && capture_width == 1 && capture_height == 1 &&
	  capture_components == 3 && capture[0] == 17 && capture[1] == 34 &&
	  capture[2] == 51 && high_state.framebuffer_captures == 1,
	  "endpoint dispatches framebuffer-plane capture to its provider");
    bu_free(capture, "factory framebuffer capture");
    CHECK(!brlobol_display_endpoint_framebuffer_capture_provider_set(endpoint,
	  NULL, &low_state),
	  "unrelated owner cannot clear an endpoint framebuffer provider");
    CHECK(brlobol_display_endpoint_framebuffer_capture_provider_set(endpoint,
	  NULL, &high_state),
	  "framebuffer provider owner can clear its endpoint binding");

    brlobol_display_endpoint_destroy(endpoint);
    CHECK(high_state.binds == 1 && high_state.opens == 1 &&
	  high_state.closes == 1 && high_state.detaches == 1 &&
	  high_state.destroys == 1,
	  "endpoint closes, detaches, and destroys its factory host");

    desc.required_capabilities = BRLOBOL_HOST_CAP_PIXEL_PRESENT;
    endpoint = brlobol_display_endpoint_create(NULL, 0);
    CHECK(!brlobol_display_endpoint_host_open(endpoint,
	  "endpoint-test-failed", &desc),
	  "endpoint reports explicit factory open failure");
    CHECK(failed_state.creates == 1 && failed_state.opens == 1 &&
	  failed_state.detaches == 1 && failed_state.destroys == 1 &&
	  failed_state.closes == 0,
	  "failed factory open detaches and destroys without close");
    brlobol_display_endpoint_destroy(endpoint);

    CHECK(brlobol_host_factory_unregister(failed_token) &&
	  brlobol_host_factory_unregister(high_token) &&
	  brlobol_host_factory_unregister(low_token),
	  "unused host factories unregister");
    CHECK(brlobol_host_factory_registry_count() == initial_count,
	  "host factory test restores registry state");
    return 0;
}

static int
test_display_endpoint_contract(void)
{
    brlobol_display_endpoint_t *endpoint =
	brlobol_display_endpoint_create(NULL, 0);
    CHECK(endpoint != NULL, "display endpoint creates an owned controller");
    CHECK(brlobol_display_endpoint_controller(endpoint) != NULL,
	  "display endpoint exposes its controller");
    CHECK(brlobol_display_endpoint_render_engine_get(endpoint) ==
	  BRLOBOL_RENDER_ENGINE_AUTO,
	  "display endpoint defaults to automatic renderer selection");
    CHECK(brlobol_display_endpoint_render_engine_set(endpoint,
	  BRLOBOL_RENDER_ENGINE_SW),
	  "display endpoint accepts a typed renderer policy");
    CHECK(brlobol_display_endpoint_render_engine_get(endpoint) ==
	  BRLOBOL_RENDER_ENGINE_SW,
	  "display endpoint retains renderer policy");
    CHECK(!brlobol_display_endpoint_render_engine_set(endpoint,
	  (enum brlobol_render_engine)99),
	  "display endpoint rejects an invalid renderer policy");

    CHECK(brlobol_display_endpoint_property_count() == 15,
	  "display endpoint exposes the typed property registry");
    int found_renderer_property = 0;
    for (size_t i = 0; i < brlobol_display_endpoint_property_count(); i++) {
	struct brlobol_endpoint_property_desc property = {};
	property.struct_size = sizeof(property);
	CHECK(brlobol_display_endpoint_property_descriptor(i, &property) ==
	      BRLOBOL_ENDPOINT_PROPERTY_OK && property.name,
	      "display endpoint property descriptors are readable");
	if (strcmp(property.name, "endpoint.renderer") == 0) {
	    found_renderer_property = property.type ==
		BRLOBOL_ENDPOINT_PROPERTY_ENUM &&
		(property.access & BRLOBOL_ENDPOINT_PROPERTY_WRITE) &&
		property.allowed_values &&
		strcmp(property.allowed_values,
		    "auto,hw,sw,rt,none,diagnostic") == 0;
	}
    }
    CHECK(found_renderer_property,
	  "renderer property declares its type, access, and allowed values");

    struct brlobol_endpoint_property_value property_value =
	BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    CHECK(brlobol_display_endpoint_property_get(endpoint,
	  "endpoint.renderer", &property_value) ==
	  BRLOBOL_ENDPOINT_PROPERTY_OK &&
	  property_value.type == BRLOBOL_ENDPOINT_PROPERTY_ENUM &&
	  strcmp(property_value.string_value, "sw") == 0,
	  "typed renderer property reads endpoint policy");
    property_value_init(&property_value);
    property_value.type = BRLOBOL_ENDPOINT_PROPERTY_ENUM;
    property_value.string_value = "quality";
    CHECK(brlobol_display_endpoint_property_set(endpoint,
	  "controller.software_wire", &property_value) ==
	  BRLOBOL_ENDPOINT_PROPERTY_OK,
	  "typed software wire property accepts QUALITY");
    BRLObolViewController *endpoint_controller =
	static_cast<BRLObolViewController *>(
	    brlobol_display_endpoint_controller(endpoint));
    CHECK(endpoint_controller->getSceneRoot() &&
	  endpoint_controller->getSceneRoot()->isOfType(
	      SoBRLSceneGroup::getClassTypeId()),
	  "endpoint-created controller starts with an authoritative scene root");
    CHECK(endpoint_controller->getSoftwareWireMode() ==
	  BRLObolViewController::SOFTWARE_WIRE_QUALITY,
	  "typed software wire property updates the controller");

    property_value_init(&property_value);
    property_value.type = BRLOBOL_ENDPOINT_PROPERTY_BOOL;
    property_value.bool_value = 0;
    CHECK(brlobol_display_endpoint_property_set(endpoint,
	  "renderer.depth_test", &property_value) ==
	  BRLOBOL_ENDPOINT_PROPERTY_OK &&
	  !endpoint_controller->isDepthTestEnabled(),
	  "typed depth property updates the Obol render environment");
    CHECK(brlobol_display_endpoint_property_set(endpoint,
	  "renderer.lighting", &property_value) ==
	  BRLOBOL_ENDPOINT_PROPERTY_OK &&
	  !endpoint_controller->isLightingEnabled(),
	  "typed lighting property updates the Obol render environment");
    CHECK(brlobol_display_endpoint_property_set(endpoint,
	  "renderer.transparency", &property_value) ==
	  BRLOBOL_ENDPOINT_PROPERTY_OK &&
	  !endpoint_controller->isTransparencyEnabled(),
	  "typed transparency property updates the Obol render action");
    property_value_init(&property_value);
    CHECK(brlobol_display_endpoint_property_get(endpoint,
	  "renderer.transparency", &property_value) ==
	  BRLOBOL_ENDPOINT_PROPERTY_OK && !property_value.bool_value,
	  "typed transparency property round trips render policy");
    property_value_init(&property_value);
    property_value.type = BRLOBOL_ENDPOINT_PROPERTY_DOUBLE;
    property_value.double_value = -4.0;
    CHECK(brlobol_display_endpoint_property_set(endpoint,
	  "renderer.clip.minimum", &property_value) ==
	  BRLOBOL_ENDPOINT_PROPERTY_OK,
	  "typed clip minimum updates controller camera policy");
    property_value.double_value = 5.0;
    CHECK(brlobol_display_endpoint_property_set(endpoint,
	  "renderer.clip.maximum", &property_value) ==
	  BRLOBOL_ENDPOINT_PROPERTY_OK,
	  "typed clip maximum updates controller camera policy");
    double clip_minimum = 0.0;
    double clip_maximum = 0.0;
    endpoint_controller->getClipBounds(clip_minimum, clip_maximum);
    CHECK(fabs(clip_minimum + 4.0) < 0.0001 &&
	  fabs(clip_maximum - 5.0) < 0.0001,
	  "typed clip bounds preserve the ordered controller range");
    struct bv_context *clip_view = bv_context_create();
    CHECK(clip_view != NULL &&
	  bv_context_dimensions_set(clip_view, 320, 240) &&
	  bv_zclip_set(bv_context_view(clip_view), 1) &&
	  endpoint_controller->syncCameraFromViewContext(clip_view),
	  "controller synchronizes enabled retained clipping");
    SoClipPlane *clip_minimum_node = find_clip_plane(
	endpoint_controller->getRenderRoot(), "BRLObolClipMinimum");
    SoClipPlane *clip_maximum_node = find_clip_plane(
	endpoint_controller->getRenderRoot(), "BRLObolClipMaximum");
    CHECK(clip_minimum_node && clip_maximum_node &&
	  clip_minimum_node->on.getValue() &&
	  clip_maximum_node->on.getValue(),
	  "zclip enables the retained Obol clipping-plane pair");
    CHECK(bv_zclip_set(bv_context_view(clip_view), 0) &&
	  endpoint_controller->syncCameraFromViewContext(clip_view) &&
	  !clip_minimum_node->on.getValue() &&
	  !clip_maximum_node->on.getValue(),
	  "disabling zclip preserves but deactivates retained clip planes");
    bv_context_destroy(clip_view);
    property_value_init(&property_value);
    property_value.type = BRLOBOL_ENDPOINT_PROPERTY_BOOL;
    property_value.bool_value = 1;
    CHECK(brlobol_display_endpoint_property_set(endpoint,
	  "renderer.depth_cue", &property_value) ==
	  BRLOBOL_ENDPOINT_PROPERTY_OK &&
	  endpoint_controller->isDepthCueEnabled(),
	  "typed depth-cue property updates the Obol render environment");

    CHECK(brlobol_display_endpoint_property_get(endpoint, "view.zclip",
	  &property_value) == BRLOBOL_ENDPOINT_PROPERTY_UNSUPPORTED,
	  "external view property reports unsupported without an owner");
    int zclip = 0;
    CHECK(brlobol_display_endpoint_property_provider_set(endpoint,
	  test_property_provider_get, test_property_provider_set, &zclip),
	  "endpoint accepts an external property owner");
    property_value.bool_value = 1;
    CHECK(brlobol_display_endpoint_property_set(endpoint, "view.zclip",
	  &property_value) == BRLOBOL_ENDPOINT_PROPERTY_OK && zclip == 1,
	  "external property setter updates owner state");
    property_value.bool_value = 0;
    CHECK(brlobol_display_endpoint_property_get(endpoint, "view.zclip",
	  &property_value) == BRLOBOL_ENDPOINT_PROPERTY_OK &&
	  property_value.bool_value == 1,
	  "external property getter reads owner state");
    CHECK(brlobol_display_endpoint_property_provider_set(endpoint,
	  NULL, NULL, NULL),
	  "endpoint clears its external property owner");

    property_value_init(&property_value);
    property_value.type = BRLOBOL_ENDPOINT_PROPERTY_COLOR3;
    property_value.color3[0] = 0.125;
    property_value.color3[1] = 0.25;
    property_value.color3[2] = 0.5;
    CHECK(brlobol_display_endpoint_property_set(endpoint,
	  "controller.background.bottom", &property_value) ==
	  BRLOBOL_ENDPOINT_PROPERTY_OK,
	  "typed background property accepts normalized RGB");
    property_value_init(&property_value);
    CHECK(brlobol_display_endpoint_property_get(endpoint,
	  "controller.background.bottom", &property_value) ==
	  BRLOBOL_ENDPOINT_PROPERTY_OK &&
	  fabs(property_value.color3[0] - 0.125) < 0.0001 &&
	  fabs(property_value.color3[1] - 0.25) < 0.0001 &&
	  fabs(property_value.color3[2] - 0.5) < 0.0001,
	  "typed background property round trips controller color");

    property_value_init(&property_value);
    property_value.type = BRLOBOL_ENDPOINT_PROPERTY_UINT;
    property_value.uint_value = 320;
    CHECK(brlobol_display_endpoint_property_set(endpoint, "endpoint.width",
	  &property_value) == BRLOBOL_ENDPOINT_PROPERTY_OK,
	  "typed width property resizes the endpoint");
    property_value.uint_value = 240;
    CHECK(brlobol_display_endpoint_property_set(endpoint, "endpoint.height",
	  &property_value) == BRLOBOL_ENDPOINT_PROPERTY_OK,
	  "typed height property resizes the endpoint");
    const SbVec2s property_viewport =
	endpoint_controller->getViewportRegion().getWindowSize();
    CHECK(property_viewport[0] == 320 && property_viewport[1] == 240,
	  "typed size properties update the controller viewport");

    property_value_init(&property_value);
    property_value.type = BRLOBOL_ENDPOINT_PROPERTY_STRING;
    property_value.string_value = "invalid";
    CHECK(brlobol_display_endpoint_property_set(endpoint, "endpoint.host",
	  &property_value) == BRLOBOL_ENDPOINT_PROPERTY_READ_ONLY,
	  "typed property API rejects writes to read-only host identity");
    CHECK(brlobol_display_endpoint_property_get(endpoint,
	  "does.not.exist", &property_value) ==
	  BRLOBOL_ENDPOINT_PROPERTY_UNKNOWN,
	  "typed property API distinguishes unknown properties");
    property_value.type = BRLOBOL_ENDPOINT_PROPERTY_DOUBLE;
    CHECK(brlobol_display_endpoint_property_set(endpoint, "endpoint.width",
	  &property_value) == BRLOBOL_ENDPOINT_PROPERTY_INVALID,
	  "typed property API rejects mismatched value types");

    property_value_init(&property_value);
    property_value.type = BRLOBOL_ENDPOINT_PROPERTY_ENUM;
    property_value.string_value = "rt";
    CHECK(brlobol_display_endpoint_property_set(endpoint,
	  "endpoint.renderer", &property_value) ==
	  BRLOBOL_ENDPOINT_PROPERTY_UNSUPPORTED,
	  "typed renderer property rejects the unavailable rt engine");

    BRLObolWindowHost borrowed_host;
    CHECK(!brlobol_display_endpoint_host_bind(endpoint, &borrowed_host, 0),
	  "explicit software renderer rejects an untyped direct host");
    CHECK(brlobol_display_endpoint_render_engine_set(endpoint,
	  BRLOBOL_RENDER_ENGINE_AUTO),
	  "automatic renderer policy permits transitional direct host binding");
    CHECK(brlobol_display_endpoint_host_bind(endpoint, &borrowed_host, 0),
	  "display endpoint binds a borrowed host");
    CHECK(brlobol_display_endpoint_host(endpoint) == &borrowed_host,
	  "display endpoint reports its host");
    CHECK(borrowed_host.getController() ==
	  brlobol_display_endpoint_controller(endpoint),
	  "bound host borrows the endpoint controller");
    brlobol_display_endpoint_host_detach(endpoint);
    CHECK(brlobol_display_endpoint_host(endpoint) == NULL &&
	  borrowed_host.getController() == NULL,
	  "display endpoint detaches a borrowed host");

    CountingWindowHost::destroyed = 0;
    CountingWindowHost *owned_host = new CountingWindowHost;
    CHECK(brlobol_display_endpoint_host_bind(endpoint, owned_host,
	  BRLOBOL_ENDPOINT_OWN_HOST),
	  "display endpoint binds an owned host");
    CHECK(brlobol_display_endpoint_host_bind(endpoint, owned_host, 0),
	  "display endpoint permits idempotent host binding");
    brlobol_display_endpoint_destroy(endpoint);
    CHECK(CountingWindowHost::destroyed == 1,
	  "idempotent binding does not relinquish owned host lifetime");

    BRLObolViewController borrowed_controller;
    BRLObolWindowHost wrapper_host;
    {
	BRLObolDisplayEndpoint wrapper(&borrowed_controller, false);
	CHECK(wrapper.isValid() && wrapper.controller() == &borrowed_controller,
	      "C++ endpoint wrapper borrows an explicit controller");
	CHECK(wrapper.bindHost(&wrapper_host) && wrapper.host() == &wrapper_host,
	      "C++ endpoint wrapper binds a typed host");
	CHECK(!wrapper.setRenderEngine(BRLOBOL_RENDER_ENGINE_HW),
	      "explicit hardware policy rejects an untyped direct host");
	CHECK(wrapper.setRenderEngine(BRLOBOL_RENDER_ENGINE_NONE) &&
	      wrapper.renderEngine() == BRLOBOL_RENDER_ENGINE_NONE,
	      "C++ endpoint wrapper exposes renderer policy");
    }
    CHECK(wrapper_host.getController() == NULL,
	  "C++ endpoint wrapper detaches its borrowed host on destruction");
    return 0;
}

static int
texture_matches(SoTexture2 *texture, int width, int height, int channels)
{
    if (!texture)
	return 0;
    int texWidth = 0;
    int texHeight = 0;
    int texChannels = 0;
    const unsigned char *pixels = texture->getImageData(texWidth, texHeight, texChannels);
    return pixels && texWidth == width && texHeight == height && texChannels == channels;
}

static int
test_imgstream_display_host_bridge(void)
{
    BRLObolWindowHost host;

    CHECK(brlobol_window_host_open_display_framebuffer(NULL, "/dev/qtgl", 5, 4) == NULL,
	  "null Obol display host is rejected");

    imgstream_fb_t *fb = brlobol_window_host_open_display_framebuffer(&host, "/dev/qtgl", 5, 4);
    CHECK(fb != NULL, "display framebuffer opens through Obol host");
    CHECK(host.isOpen(), "display framebuffer opens host");
    CHECK(host.getFramebufferCount() == 1, "host records one framebuffer");

    SoBRLImageSource *source = host.getFramebufferImageSource(fb);
    SoBRLViewportImage *viewport = host.getFramebufferViewportImage(fb);
    CHECK(source != NULL, "display host creates image source");
    CHECK(viewport != NULL, "display host creates viewport image");
    CHECK(source->getStream() == imgstream_fb_stream(fb),
	  "display host image source borrows framebuffer stream");
    CHECK(viewport->getImageSource() == source,
	  "display host viewport image references source");
    CHECK(viewport->getTextureNode() != NULL &&
	  texture_matches(viewport->getTextureNode(), 5, 4, 3),
	  "display host realizes framebuffer texture");
    CHECK(host.getController()->getSceneRoot()->isOfType(SoGroup::getClassTypeId()),
	  "display host controller root is a group");
    CHECK(static_cast<SoGroup *>(host.getController()->getSceneRoot())->getNumChildren() >= 2,
	  "display host attaches image nodes to controller root");

    unsigned char red[3] = {255, 0, 0};
    CHECK(imgstream_fb_writerect(fb, 2, 1, 1, 1, red) == 1,
	  "framebuffer write updates stream");
    CHECK(imgstream_fb_flush(fb) == 0, "display host flush syncs stream");
    CHECK(viewport->realizedDirtyRevision.getValue() == source->dirtyRevision.getValue(),
	  "display host flush refreshes viewport dirty revision");

    CHECK(imgstream_fb_view(fb, 3, 2, 4, 4) == 0,
	  "display host accepts framebuffer view transform");
    CHECK(float_equal(viewport->sourceCenter.getValue()[0], 3.0f) &&
	  float_equal(viewport->sourceCenter.getValue()[1], 2.0f) &&
	  float_equal(viewport->sourceZoom.getValue(), 4.0f),
	  "display host maps framebuffer view to viewport image");

    CHECK(imgstream_fb_cursor(fb, 1, 4, 3) == 0,
	  "display host accepts cursor state");
    CHECK(viewport->cursorVisible.getValue() == TRUE &&
	  float_equal(viewport->cursorImagePosition.getValue()[0], 4.0f) &&
	  float_equal(viewport->cursorImagePosition.getValue()[1], 3.0f),
	  "display host maps cursor state to viewport image");

    CHECK(imgstream_fb_viewport(fb, 1, 2, 11, 8) == 0,
	  "display host accepts viewport state");
    CHECK(host.getController()->getViewportRegion().getWindowSize()[0] == 10 &&
	  host.getController()->getViewportRegion().getWindowSize()[1] == 6,
	  "display host maps viewport state to controller size");

    unsigned char cursorBits[1] = {0xff};
    CHECK(imgstream_fb_setcursor(fb, cursorBits, 8, 1, 0, 0) == 0,
	  "display host accepts custom cursor shape");
    CHECK(viewport->cursorShape.getValue() == SoBRLViewportImage::CURSOR_CUSTOM,
	  "display host records custom cursor policy");

    CHECK(imgstream_fb_poll(fb) == 0, "display host poll succeeds");
    CHECK(imgstream_fb_poll_rate(fb) == 0, "display host poll rate defaults to zero");

    imgstream_fb_close(fb);
    CHECK(host.getFramebufferCount() == 0, "display host closes framebuffer attachment");

    CHECK(imgstream_fb_open("/dev/qtgl", 2, 2) == NULL,
	  "display framebuffer requires an explicit host");
    return 0;
}

static int
test_framebuffer_stream_helper(void)
{
    BRLObolWindowHost host;
    BRLObolFramebufferStream stream(&host);

    CHECK(stream.configure(4, 3) == 0, "framebuffer stream accepts size");

    BRLObolFramebufferInfo info;
    CHECK(stream.info(&info) == 0, "framebuffer stream reports dimensions");
    CHECK(info.width == 4 && info.height == 3 &&
	  info.max_width == 4 && info.max_height == 3,
	  "framebuffer stream dimensions match configuration");
    CHECK(stream.framebuffer() != NULL, "framebuffer stream opens imgstream framebuffer");
    CHECK(host.getFramebufferCount() == 1,
	  "framebuffer stream attaches one image to Obol host");

    imgstream_fb_t *fb = stream.framebuffer();
    SoBRLImageSource *source = host.getFramebufferImageSource(fb);
    SoBRLViewportImage *viewport = host.getFramebufferViewportImage(fb);
    CHECK(source != NULL && viewport != NULL,
	  "framebuffer stream exposes host image nodes");

    unsigned char blue[3] = {0, 0, 255};
    CHECK(stream.writerect(1, 1, 1, 1, blue) == 1,
	  "framebuffer stream writes pixels");
    CHECK(source->hasPendingStreamUpdate() == TRUE,
	  "framebuffer stream write marks source update pending");
    CHECK(stream.present() == 0,
	  "framebuffer stream presents pending pixels to Obol image");
    CHECK(viewport->realizedDirtyRevision.getValue() ==
	  source->dirtyRevision.getValue(),
	  "framebuffer stream present refreshes Obol image");
    CHECK(source->hasPendingStreamUpdate() == FALSE,
	  "framebuffer stream present clears pending source update");

    CHECK(stream.view(2, 1, 3, 3) == 0,
	  "framebuffer stream records view state");
    CHECK(stream.present() == 0,
	  "framebuffer stream presents view state");
    CHECK(float_equal(viewport->sourceCenter.getValue()[0], 2.0f) &&
	  float_equal(viewport->sourceCenter.getValue()[1], 1.0f) &&
	  float_equal(viewport->sourceZoom.getValue(), 3.0f),
	  "framebuffer stream maps view state to viewport image");

    CHECK(stream.cursor(1, 3, 2) == 0,
	  "framebuffer stream records cursor state");
    CHECK(stream.present() == 0,
	  "framebuffer stream presents cursor state");
    CHECK(viewport->cursorVisible.getValue() == TRUE &&
	  float_equal(viewport->cursorImagePosition.getValue()[0], 3.0f) &&
	  float_equal(viewport->cursorImagePosition.getValue()[1], 2.0f),
	  "framebuffer stream maps cursor state to viewport image");

    stream.close();
    CHECK(host.getFramebufferCount() == 0,
	  "framebuffer stream close detaches host image nodes");

    CHECK(stream.configure(2, 2) == 0 && stream.ensure() == 0,
	  "framebuffer stream reopens after close");
    CHECK(host.getFramebufferCount() == 1,
	  "framebuffer stream reattaches after close");
    stream.setHost(NULL);
    CHECK(host.getFramebufferCount() == 0,
	  "framebuffer stream host change closes attachment");
    return 0;
}

int
main(int ac, char **av)
{
    (void)ac;
    (void)av;

    brlobol_init(NULL);

    if (test_window_host_contract())
	return 1;
    if (test_display_endpoint_contract())
	return 1;
    if (test_host_factory_contract())
	return 1;
    if (test_imgstream_display_host_bridge())
	return 1;
    if (test_framebuffer_stream_helper())
	return 1;

    return 0;
}
