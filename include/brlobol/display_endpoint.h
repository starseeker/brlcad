/*              D I S P L A Y _ E N D P O I N T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/display_endpoint.h */

#ifndef BRLOBOL_DISPLAY_ENDPOINT_H
#define BRLOBOL_DISPLAY_ENDPOINT_H

#include "brlobol/defines.h"
#include "brlobol/host_factory.h"
#include "brlobol/input.h"

struct brlobol_display_endpoint;
typedef struct brlobol_display_endpoint brlobol_display_endpoint_t;

enum brlobol_render_engine {
    BRLOBOL_RENDER_ENGINE_AUTO = 0,
    BRLOBOL_RENDER_ENGINE_HW = 1,
    BRLOBOL_RENDER_ENGINE_SW = 2,
    BRLOBOL_RENDER_ENGINE_RT = 3,
    BRLOBOL_RENDER_ENGINE_NONE = 4,
    BRLOBOL_RENDER_ENGINE_DIAGNOSTIC = 5
};

enum brlobol_capture_plane {
    BRLOBOL_CAPTURE_COMPOSITE = 0,
    BRLOBOL_CAPTURE_FRAMEBUFFER = 1
};

/** Typed retained-RT planes.  DEPTH captures float samples in [0, 1];
 * SOURCE_ID captures uint32_t samples where zero denotes background. */
enum brlobol_rt_output_plane {
    BRLOBOL_RT_OUTPUT_DEPTH = 0,
    BRLOBOL_RT_OUTPUT_SOURCE_ID = 1
};

/** Source metadata addressed by BRLOBOL_RT_OUTPUT_SOURCE_ID.  The caller owns
 * instance_key and path after a successful query and releases them with
 * brlobol_display_endpoint_rt_source_identity_clear(). */
struct brlobol_rt_source_identity {
    uint32_t struct_size;
    void *database;
    char *instance_key;
    char *path;
    uint32_t source_revision;
};

#define BRLOBOL_RT_SOURCE_IDENTITY_INIT { \
    sizeof(struct brlobol_rt_source_identity), NULL, NULL, NULL, 0u }

typedef int (*brlobol_endpoint_framebuffer_capture_callback)(void *user_data,
	unsigned char **pixels, size_t *size, unsigned int *width,
	unsigned int *height, unsigned int *components);

enum brlobol_endpoint_property_type {
    BRLOBOL_ENDPOINT_PROPERTY_BOOL = 1,
    BRLOBOL_ENDPOINT_PROPERTY_INT = 2,
    BRLOBOL_ENDPOINT_PROPERTY_UINT = 3,
    BRLOBOL_ENDPOINT_PROPERTY_DOUBLE = 4,
    BRLOBOL_ENDPOINT_PROPERTY_STRING = 5,
    BRLOBOL_ENDPOINT_PROPERTY_COLOR3 = 6,
    BRLOBOL_ENDPOINT_PROPERTY_ENUM = 7
};

#define BRLOBOL_ENDPOINT_PROPERTY_READ  0x01u
#define BRLOBOL_ENDPOINT_PROPERTY_WRITE 0x02u

enum brlobol_endpoint_property_result {
    BRLOBOL_ENDPOINT_PROPERTY_OK = 1,
    BRLOBOL_ENDPOINT_PROPERTY_UNKNOWN = 0,
    BRLOBOL_ENDPOINT_PROPERTY_INVALID = -1,
    BRLOBOL_ENDPOINT_PROPERTY_READ_ONLY = -2,
    BRLOBOL_ENDPOINT_PROPERTY_UNSUPPORTED = -3
};

struct brlobol_endpoint_property_desc {
    uint32_t struct_size;
    const char *name;
    enum brlobol_endpoint_property_type type;
    unsigned int access;
    uint64_t required_host_capabilities;
    double minimum;
    double maximum;
    const char *allowed_values;
};

struct brlobol_endpoint_property_value {
    uint32_t struct_size;
    enum brlobol_endpoint_property_type type;
    int bool_value;
    int64_t int_value;
    uint64_t uint_value;
    double double_value;
    double color3[3];
    const char *string_value;
};

typedef int (*brlobol_endpoint_property_get_callback)(void *user_data,
	const char *name, struct brlobol_endpoint_property_value *out);
typedef int (*brlobol_endpoint_property_set_callback)(void *user_data,
	const char *name,
	const struct brlobol_endpoint_property_value *value);

#define BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT { \
    sizeof(struct brlobol_endpoint_property_value), \
    BRLOBOL_ENDPOINT_PROPERTY_BOOL, 0, 0, 0, 0.0, {0.0, 0.0, 0.0}, NULL }

#define BRLOBOL_ENDPOINT_OWN_CONTROLLER 0x01u
#define BRLOBOL_ENDPOINT_OWN_HOST       0x01u

__BEGIN_DECLS

/**
 * Create an endpoint around an optional BRLObolViewController.
 *
 * A NULL controller creates an endpoint-owned controller.  A non-NULL
 * controller is borrowed unless BRLOBOL_ENDPOINT_OWN_CONTROLLER is set.
 */
BRLOBOL_EXPORT brlobol_display_endpoint_t *
brlobol_display_endpoint_create(void *controller, unsigned int flags);

BRLOBOL_EXPORT void
brlobol_display_endpoint_destroy(brlobol_display_endpoint_t *endpoint);

BRLOBOL_EXPORT void *
brlobol_display_endpoint_controller(const brlobol_display_endpoint_t *endpoint);

/** Synchronize the endpoint camera from a libbv view context. */
BRLOBOL_EXPORT int
brlobol_display_endpoint_view_sync(brlobol_display_endpoint_t *endpoint,
	const void *view_ctx);

/**
 * Bind an existing BRLObolWindowHost to the endpoint controller.
 *
 * Binding closes the host before changing its controller.  The host is
 * borrowed unless BRLOBOL_ENDPOINT_OWN_HOST is set.
 */
BRLOBOL_EXPORT int
brlobol_display_endpoint_host_bind(brlobol_display_endpoint_t *endpoint,
	void *host, unsigned int flags);

BRLOBOL_EXPORT void
brlobol_display_endpoint_host_detach(brlobol_display_endpoint_t *endpoint);

BRLOBOL_EXPORT void *
brlobol_display_endpoint_host(const brlobol_display_endpoint_t *endpoint);

/**
 * Return a borrowed BRLObolWindowHost offered by the active factory for
 * retained framebuffer attachment.  Toolkit hosts without a compatible
 * retained attachment return NULL; callers must then use controller-owned
 * presentation instead.  The result remains opaque to C callers.
 */
BRLOBOL_EXPORT void *
brlobol_display_endpoint_framebuffer_window_host(
	const brlobol_display_endpoint_t *endpoint);

/** Select a compatible registered factory and create its host instance. */
BRLOBOL_EXPORT int
brlobol_display_endpoint_host_open(brlobol_display_endpoint_t *endpoint,
	const char *factory_name, const struct brlobol_host_desc *desc);

BRLOBOL_EXPORT const char *
brlobol_display_endpoint_host_factory_name(
	const brlobol_display_endpoint_t *endpoint);

/** Return capabilities of the active factory host, or zero if unbound. */
BRLOBOL_EXPORT uint64_t brlobol_display_endpoint_host_capabilities(
	const brlobol_display_endpoint_t *endpoint);

BRLOBOL_EXPORT int
brlobol_display_endpoint_request_frame(brlobol_display_endpoint_t *endpoint,
	const char *reason);

BRLOBOL_EXPORT int
brlobol_display_endpoint_resize(brlobol_display_endpoint_t *endpoint,
	unsigned int width, unsigned int height, double device_pixel_ratio);

/** Set an endpoint-local input profile.  A NULL profile restores BRL-CAD's
 * standard view bindings. */
BRLOBOL_EXPORT int
brlobol_display_endpoint_input_profile_set(
	brlobol_display_endpoint_t *endpoint,
	const BRLObolInputProfile *profile);

/** Set an endpoint-local action handler.  A NULL handler leaves events
 * unhandled so the embedding application may retain its native bindings. */
BRLOBOL_EXPORT int
brlobol_display_endpoint_input_action_handler_set(
	brlobol_display_endpoint_t *endpoint,
	BRLObolInputActionHandler handler, void *user_data);

/** Clear an action handler only when it is still owned by the specified
 * adapter.  This protects application-installed handlers during teardown. */
BRLOBOL_EXPORT int
brlobol_display_endpoint_input_action_handler_clear_if(
	brlobol_display_endpoint_t *endpoint,
	BRLObolInputActionHandler handler, void *user_data);

/**
 * Install a scoped application semantic-gesture overlay.  The overlay takes
 * precedence over equal-score default bindings but does not replace default
 * view navigation.  Only its owner may clear or replace it.
 */
BRLOBOL_EXPORT int
brlobol_display_endpoint_input_semantic_profile_set(
	brlobol_display_endpoint_t *endpoint,
	const BRLObolInputProfile *profile, void *owner);

BRLOBOL_EXPORT int
brlobol_display_endpoint_input_semantic_profile_clear_if(
	brlobol_display_endpoint_t *endpoint, void *owner);

/**
 * Set the endpoint-local application semantic-action handler.  Semantic
 * actions return BRLOBOL_INPUT_RESULT_*; view navigation remains owned by
 * brlobol_display_endpoint_input_action_handler_set.
 */
BRLOBOL_EXPORT int
brlobol_display_endpoint_input_semantic_action_handler_set(
	brlobol_display_endpoint_t *endpoint,
	BRLObolInputActionHandler handler, void *user_data);

BRLOBOL_EXPORT int
brlobol_display_endpoint_input_semantic_action_handler_clear_if(
	brlobol_display_endpoint_t *endpoint,
	BRLObolInputActionHandler handler, void *user_data);

/** Match and deliver a normalized event using this endpoint's profile. */
BRLOBOL_EXPORT int
brlobol_display_endpoint_input_dispatch(brlobol_display_endpoint_t *endpoint,
	const BRLObolInputEvent *event);

/**
 * Apply one common faceplate-toggle action to an endpoint-backed view.  The
 * action may be BRLOBOL_ACTION_TOGGLE_ADC,
 * BRLOBOL_ACTION_TOGGLE_MODEL_AXES, or BRLOBOL_ACTION_TOGGLE_VIEW_AXES.
 *
 * An endpoint property owner is preferred.  If that property is explicitly
 * unsupported, the supplied bv_context receives the equivalent passive state
 * update for a standalone view.  Other actions and property failures are
 * left unhandled.  @p visible receives the resulting visibility when non-NULL.
 */
BRLOBOL_EXPORT int
brlobol_display_endpoint_input_faceplate_toggle_apply(
	brlobol_display_endpoint_t *endpoint, void *view_ctx,
	BRLObolInputAction action, int *visible);

/** Capture RGB/RGBA pixels in bottom-left row order.  A successful call
 * returns storage freed with bu_free. */
BRLOBOL_EXPORT int
brlobol_display_endpoint_capture(brlobol_display_endpoint_t *endpoint,
	unsigned char **pixels, size_t *size, unsigned int *width,
	unsigned int *height, unsigned int *components);

/** Capture a selected retained composition plane. */
BRLOBOL_EXPORT int
brlobol_display_endpoint_capture_plane(brlobol_display_endpoint_t *endpoint,
	enum brlobol_capture_plane plane, unsigned char **pixels, size_t *size,
	unsigned int *width, unsigned int *height, unsigned int *components);

/** Capture the final retained-RT output plane.  The returned storage is
 * released with bu_free.  It is unavailable unless the endpoint is using the
 * retained RT engine and has produced a final frame. */
BRLOBOL_EXPORT int
brlobol_display_endpoint_rt_plane_capture(brlobol_display_endpoint_t *endpoint,
	enum brlobol_rt_output_plane plane, void **samples, size_t *size,
	unsigned int *width, unsigned int *height);

/** Resolve a nonzero RT source-ID sample from the currently displayed RT
 * snapshot.  Initialize @p out with BRLOBOL_RT_SOURCE_IDENTITY_INIT or clear
 * a prior result before reuse. */
BRLOBOL_EXPORT int
brlobol_display_endpoint_rt_source_identity_get(
	const brlobol_display_endpoint_t *endpoint, uint32_t identifier,
	struct brlobol_rt_source_identity *out);

BRLOBOL_EXPORT void
brlobol_display_endpoint_rt_source_identity_clear(
	struct brlobol_rt_source_identity *identity);

/** Bind or conditionally clear the endpoint's framebuffer-plane provider.
 * A NULL callback clears only a provider whose data pointer matches user_data. */
BRLOBOL_EXPORT int
brlobol_display_endpoint_framebuffer_capture_provider_set(
	brlobol_display_endpoint_t *endpoint,
	brlobol_endpoint_framebuffer_capture_callback callback,
	void *user_data);

BRLOBOL_EXPORT int
brlobol_display_endpoint_render_engine_set(
	brlobol_display_endpoint_t *endpoint,
	enum brlobol_render_engine engine);

BRLOBOL_EXPORT enum brlobol_render_engine
brlobol_display_endpoint_render_engine_get(
	const brlobol_display_endpoint_t *endpoint);

/** Enumerate the stable typed endpoint property registry. */
BRLOBOL_EXPORT size_t brlobol_display_endpoint_property_count(void);
BRLOBOL_EXPORT int brlobol_display_endpoint_property_descriptor(
	size_t index, struct brlobol_endpoint_property_desc *out);

/** Get or set one typed property.  Returns brlobol_endpoint_property_result. */
BRLOBOL_EXPORT int brlobol_display_endpoint_property_get(
	const brlobol_display_endpoint_t *endpoint, const char *name,
	struct brlobol_endpoint_property_value *out);
BRLOBOL_EXPORT int brlobol_display_endpoint_property_set(
	brlobol_display_endpoint_t *endpoint, const char *name,
	const struct brlobol_endpoint_property_value *value);

/** Bind optional storage callbacks for properties owned outside libbrlobol. */
BRLOBOL_EXPORT int brlobol_display_endpoint_property_provider_set(
	brlobol_display_endpoint_t *endpoint,
	brlobol_endpoint_property_get_callback get_callback,
	brlobol_endpoint_property_set_callback set_callback,
	void *user_data);

__END_DECLS

#ifdef __cplusplus

class BRLObolViewController;
class BRLObolWindowHost;

/** RAII wrapper for the C display-endpoint lifetime. */
class BRLOBOL_EXPORT BRLObolDisplayEndpoint {
public:
    explicit BRLObolDisplayEndpoint(BRLObolViewController *controller = NULL,
	bool takeOwnership = false);
    ~BRLObolDisplayEndpoint(void);

    bool isValid(void) const;
    BRLObolViewController *controller(void) const;
    int bindHost(BRLObolWindowHost *host, bool takeOwnership = false);
    void detachHost(void);
    BRLObolWindowHost *host(void) const;
    int setRenderEngine(enum brlobol_render_engine engine);
    enum brlobol_render_engine renderEngine(void) const;
    int propertyGet(const char *name,
	struct brlobol_endpoint_property_value *value) const;
    int propertySet(const char *name,
	const struct brlobol_endpoint_property_value *value);

    brlobol_display_endpoint_t *release(void);
    brlobol_display_endpoint_t *get(void) const;

private:
    BRLObolDisplayEndpoint(const BRLObolDisplayEndpoint &) = delete;
    BRLObolDisplayEndpoint &operator=(const BRLObolDisplayEndpoint &) = delete;

    brlobol_display_endpoint_t *endpoint;
};

#endif

#endif /* BRLOBOL_DISPLAY_ENDPOINT_H */
