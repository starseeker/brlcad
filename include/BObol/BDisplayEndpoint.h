/*            B D I S P L A Y E N D P O I N T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BDisplayEndpoint.h */

#ifndef BOBOL_BDISPLAYENDPOINT_H
#define BOBOL_BDISPLAYENDPOINT_H

#include "BObol/BDefines.h"
#include "BObol/BHostFactory.h"
#include "BObol/BInput.h"

struct bobol_display_endpoint;
typedef struct bobol_display_endpoint bobol_display_endpoint_t;

enum bobol_render_engine {
    BOBOL_RENDER_ENGINE_AUTO = 0,
    BOBOL_RENDER_ENGINE_HW = 1,
    BOBOL_RENDER_ENGINE_SW = 2,
    BOBOL_RENDER_ENGINE_RT = 3,
    BOBOL_RENDER_ENGINE_NONE = 4,
    BOBOL_RENDER_ENGINE_DIAGNOSTIC = 5
};

#define BOBOL_RENDER_ENGINE_CAP_AUTO       UINT64_C(0x00000001)
#define BOBOL_RENDER_ENGINE_CAP_HW         UINT64_C(0x00000002)
#define BOBOL_RENDER_ENGINE_CAP_SW         UINT64_C(0x00000004)
#define BOBOL_RENDER_ENGINE_CAP_RT         UINT64_C(0x00000008)
#define BOBOL_RENDER_ENGINE_CAP_NONE       UINT64_C(0x00000010)
#define BOBOL_RENDER_ENGINE_CAP_DIAGNOSTIC UINT64_C(0x00000020)

enum bobol_capture_plane {
    BOBOL_CAPTURE_COMPOSITE = 0,
    BOBOL_CAPTURE_FRAMEBUFFER = 1
};

/** Typed retained-RT planes.  DEPTH captures float samples in [0, 1];
 * SOURCE_ID captures uint32_t samples where zero denotes background. */
enum bobol_rt_output_plane {
    BOBOL_RT_OUTPUT_DEPTH = 0,
    BOBOL_RT_OUTPUT_SOURCE_ID = 1
};

/** Source metadata addressed by BOBOL_RT_OUTPUT_SOURCE_ID.  The caller owns
 * instance_key and path after a successful query and releases them with
 * bobol_display_endpoint_rt_source_identity_clear(). */
struct bobol_rt_source_identity {
    uint32_t struct_size;
    void *database;
    char *instance_key;
    char *path;
    uint32_t source_revision;
};

#define BOBOL_RT_SOURCE_IDENTITY_INIT { \
    sizeof(struct bobol_rt_source_identity), NULL, NULL, NULL, 0u }

typedef int (*bobol_endpoint_framebuffer_capture_callback)(void *user_data,
	unsigned char **pixels, size_t *size, unsigned int *width,
	unsigned int *height, unsigned int *components);

enum bobol_endpoint_property_type {
    BOBOL_ENDPOINT_PROPERTY_BOOL = 1,
    BOBOL_ENDPOINT_PROPERTY_INT = 2,
    BOBOL_ENDPOINT_PROPERTY_UINT = 3,
    BOBOL_ENDPOINT_PROPERTY_DOUBLE = 4,
    BOBOL_ENDPOINT_PROPERTY_STRING = 5,
    BOBOL_ENDPOINT_PROPERTY_COLOR3 = 6,
    BOBOL_ENDPOINT_PROPERTY_ENUM = 7
};

#define BOBOL_ENDPOINT_PROPERTY_READ  0x01u
#define BOBOL_ENDPOINT_PROPERTY_WRITE 0x02u

enum bobol_endpoint_property_result {
    BOBOL_ENDPOINT_PROPERTY_OK = 1,
    BOBOL_ENDPOINT_PROPERTY_UNKNOWN = 0,
    BOBOL_ENDPOINT_PROPERTY_INVALID = -1,
    BOBOL_ENDPOINT_PROPERTY_READ_ONLY = -2,
    BOBOL_ENDPOINT_PROPERTY_UNSUPPORTED = -3
};

struct bobol_endpoint_property_desc {
    uint32_t struct_size;
    const char *name;
    enum bobol_endpoint_property_type type;
    unsigned int access;
    uint64_t required_host_capabilities;
    double minimum;
    double maximum;
    const char *allowed_values;
};

struct bobol_endpoint_property_value {
    uint32_t struct_size;
    enum bobol_endpoint_property_type type;
    int bool_value;
    int64_t int_value;
    uint64_t uint_value;
    double double_value;
    double color3[3];
    const char *string_value;
};

typedef int (*bobol_endpoint_property_get_callback)(void *user_data,
	const char *name, struct bobol_endpoint_property_value *out);
typedef int (*bobol_endpoint_property_set_callback)(void *user_data,
	const char *name,
	const struct bobol_endpoint_property_value *value);

#define BOBOL_ENDPOINT_PROPERTY_VALUE_INIT { \
    sizeof(struct bobol_endpoint_property_value), \
    BOBOL_ENDPOINT_PROPERTY_BOOL, 0, 0, 0, 0.0, {0.0, 0.0, 0.0}, NULL }

#define BOBOL_ENDPOINT_OWN_CONTROLLER 0x01u
#define BOBOL_ENDPOINT_OWN_HOST       0x01u

__BEGIN_DECLS

/**
 * Create an endpoint around an optional BObolViewController.
 *
 * A NULL controller creates an endpoint-owned controller.  A non-NULL
 * controller is borrowed unless BOBOL_ENDPOINT_OWN_CONTROLLER is set.
 */
BOBOL_EXPORT bobol_display_endpoint_t *
bobol_display_endpoint_create(void *controller, unsigned int flags);

BOBOL_EXPORT void
bobol_display_endpoint_destroy(bobol_display_endpoint_t *endpoint);

BOBOL_EXPORT void *
bobol_display_endpoint_controller(const bobol_display_endpoint_t *endpoint);

/** Synchronize the endpoint camera from a libbv view context. */
BOBOL_EXPORT int
bobol_display_endpoint_view_sync(bobol_display_endpoint_t *endpoint,
	const void *view_ctx);

/**
 * Bind an existing BObolWindowHost to the endpoint controller.
 *
 * Binding closes the host before changing its controller.  The host is
 * borrowed unless BOBOL_ENDPOINT_OWN_HOST is set.
 */
BOBOL_EXPORT int
bobol_display_endpoint_host_bind(bobol_display_endpoint_t *endpoint,
	void *host, unsigned int flags);

BOBOL_EXPORT void
bobol_display_endpoint_host_detach(bobol_display_endpoint_t *endpoint);

BOBOL_EXPORT void *
bobol_display_endpoint_host(const bobol_display_endpoint_t *endpoint);

/**
 * Return a borrowed BObolWindowHost offered by the active factory for
 * retained framebuffer attachment.  Toolkit hosts without a compatible
 * retained attachment return NULL; callers must then use controller-owned
 * presentation instead.  The result remains opaque to C callers.
 */
BOBOL_EXPORT void *
bobol_display_endpoint_framebuffer_window_host(
	const bobol_display_endpoint_t *endpoint);

/** Select a compatible registered factory and create its host instance. */
BOBOL_EXPORT int
bobol_display_endpoint_host_open(bobol_display_endpoint_t *endpoint,
	const char *factory_name, const struct bobol_host_desc *desc);

BOBOL_EXPORT const char *
bobol_display_endpoint_host_factory_name(
	const bobol_display_endpoint_t *endpoint);

/** Return capabilities of the active factory host, or zero if unbound. */
BOBOL_EXPORT uint64_t bobol_display_endpoint_host_capabilities(
	const bobol_display_endpoint_t *endpoint);

BOBOL_EXPORT int
bobol_display_endpoint_request_frame(bobol_display_endpoint_t *endpoint,
	const char *reason);

BOBOL_EXPORT int
bobol_display_endpoint_resize(bobol_display_endpoint_t *endpoint,
	unsigned int width, unsigned int height, double device_pixel_ratio);

/** Set an endpoint-local input profile.  A NULL profile restores BRL-CAD's
 * standard view bindings. */
BOBOL_EXPORT int
bobol_display_endpoint_input_profile_set(
	bobol_display_endpoint_t *endpoint,
	const BObolInputProfile *profile);

/** Set an endpoint-local action handler.  A NULL handler leaves events
 * unhandled so the embedding application may retain its native bindings. */
BOBOL_EXPORT int
bobol_display_endpoint_input_action_handler_set(
	bobol_display_endpoint_t *endpoint,
	BObolInputActionHandler handler, void *user_data);

/** Clear an action handler only when it is still owned by the specified
 * adapter.  This protects application-installed handlers during teardown. */
BOBOL_EXPORT int
bobol_display_endpoint_input_action_handler_clear_if(
	bobol_display_endpoint_t *endpoint,
	BObolInputActionHandler handler, void *user_data);

/**
 * Atomically install a scoped, application-defined action layer.  Its
 * bindings take precedence over equal-score default bindings, and winning
 * actions dispatch only to the handler stored in the layer.  Action IDs are
 * opaque to libBObol and local to that handler.  An UNHANDLED result falls
 * through to the matching default action.  Only the owner token may replace
 * or clear an installed layer.
 */
BOBOL_EXPORT int
bobol_display_endpoint_input_action_layer_set(
	bobol_display_endpoint_t *endpoint,
	const BObolInputActionLayer *layer, void *owner, void *user_data);

BOBOL_EXPORT int
bobol_display_endpoint_input_action_layer_clear_if(
	bobol_display_endpoint_t *endpoint, void *owner);

/** Match and deliver a normalized event using this endpoint's profile. */
BOBOL_EXPORT int
bobol_display_endpoint_input_dispatch(bobol_display_endpoint_t *endpoint,
	const BObolInputEvent *event);

/**
 * Apply one common faceplate-toggle action to an endpoint-backed view.  The
 * action may be BOBOL_ACTION_TOGGLE_ADC,
 * BOBOL_ACTION_TOGGLE_MODEL_AXES, or BOBOL_ACTION_TOGGLE_VIEW_AXES.
 *
 * An endpoint property owner is preferred.  If that property is explicitly
 * unsupported, the supplied bv_context receives the equivalent passive state
 * update for a standalone view.  Other actions and property failures are
 * left unhandled.  @p visible receives the resulting visibility when non-NULL.
 */
BOBOL_EXPORT int
bobol_display_endpoint_input_faceplate_toggle_apply(
	bobol_display_endpoint_t *endpoint, void *view_ctx,
	BObolInputAction action, int *visible);

/** Capture RGB/RGBA pixels in bottom-left row order.  A successful call
 * returns storage freed with bu_free. */
BOBOL_EXPORT int
bobol_display_endpoint_capture(bobol_display_endpoint_t *endpoint,
	unsigned char **pixels, size_t *size, unsigned int *width,
	unsigned int *height, unsigned int *components);

/** Capture a selected retained composition plane. */
BOBOL_EXPORT int
bobol_display_endpoint_capture_plane(bobol_display_endpoint_t *endpoint,
	enum bobol_capture_plane plane, unsigned char **pixels, size_t *size,
	unsigned int *width, unsigned int *height, unsigned int *components);

/** Capture the final retained-RT output plane.  The returned storage is
 * released with bu_free.  It is unavailable unless the endpoint is using the
 * retained RT engine and has produced a final frame. */
BOBOL_EXPORT int
bobol_display_endpoint_rt_plane_capture(bobol_display_endpoint_t *endpoint,
	enum bobol_rt_output_plane plane, void **samples, size_t *size,
	unsigned int *width, unsigned int *height);

/** Resolve a nonzero RT source-ID sample from the currently displayed RT
 * snapshot.  Initialize @p out with BOBOL_RT_SOURCE_IDENTITY_INIT or clear
 * a prior result before reuse. */
BOBOL_EXPORT int
bobol_display_endpoint_rt_source_identity_get(
	const bobol_display_endpoint_t *endpoint, uint32_t identifier,
	struct bobol_rt_source_identity *out);

BOBOL_EXPORT void
bobol_display_endpoint_rt_source_identity_clear(
	struct bobol_rt_source_identity *identity);

/** Bind or conditionally clear the endpoint's framebuffer-plane provider.
 * A NULL callback clears only a provider whose data pointer matches user_data. */
BOBOL_EXPORT int
bobol_display_endpoint_framebuffer_capture_provider_set(
	bobol_display_endpoint_t *endpoint,
	bobol_endpoint_framebuffer_capture_callback callback,
	void *user_data);

BOBOL_EXPORT int
bobol_display_endpoint_render_engine_set(
	bobol_display_endpoint_t *endpoint,
	enum bobol_render_engine engine);

/** Report whether the current endpoint/host binding can honor an engine.
 * An unbound endpoint reports policies that may be selected before opening a
 * compatible host. */
BOBOL_EXPORT int
bobol_display_endpoint_render_engine_supported(
	const bobol_display_endpoint_t *endpoint,
	enum bobol_render_engine engine);

/** Return BOBOL_RENDER_ENGINE_CAP_* bits for the current binding. */
BOBOL_EXPORT uint64_t
bobol_display_endpoint_render_engine_capabilities(
	const bobol_display_endpoint_t *endpoint);

/** Perform one non-graphical diagnostic traversal.  This is available only
 * while the diagnostic engine is selected; typed diagnostic properties retain
 * the resulting counts and summary. */
BOBOL_EXPORT int
bobol_display_endpoint_diagnostic_refresh(
	bobol_display_endpoint_t *endpoint);

BOBOL_EXPORT enum bobol_render_engine
bobol_display_endpoint_render_engine_get(
	const bobol_display_endpoint_t *endpoint);

/** Return the renderer selected by the endpoint policy and active factory.
 * AUTO resolves to HW before SW for a capable factory host.  It remains AUTO
 * while unbound or when a transitional direct host has no capability record.
 * Retained RT is never selected implicitly. */
BOBOL_EXPORT enum bobol_render_engine
bobol_display_endpoint_render_engine_resolved_get(
	const bobol_display_endpoint_t *endpoint);

/** Enumerate the stable typed endpoint property registry. */
BOBOL_EXPORT size_t bobol_display_endpoint_property_count(void);
BOBOL_EXPORT int bobol_display_endpoint_property_descriptor(
	size_t index, struct bobol_endpoint_property_desc *out);

/** Get or set one typed property.  Returns bobol_endpoint_property_result. */
BOBOL_EXPORT int bobol_display_endpoint_property_get(
	const bobol_display_endpoint_t *endpoint, const char *name,
	struct bobol_endpoint_property_value *out);
BOBOL_EXPORT int bobol_display_endpoint_property_set(
	bobol_display_endpoint_t *endpoint, const char *name,
	const struct bobol_endpoint_property_value *value);

/** Bind optional storage callbacks for properties owned outside libBObol. */
BOBOL_EXPORT int bobol_display_endpoint_property_provider_set(
	bobol_display_endpoint_t *endpoint,
	bobol_endpoint_property_get_callback get_callback,
	bobol_endpoint_property_set_callback set_callback,
	void *user_data);

__END_DECLS

#ifdef __cplusplus

class BObolViewController;
class BObolWindowHost;

/** RAII wrapper for the C display-endpoint lifetime. */
class BOBOL_EXPORT BObolDisplayEndpoint {
public:
    explicit BObolDisplayEndpoint(BObolViewController *controller = NULL,
	bool takeOwnership = false);
    ~BObolDisplayEndpoint(void);

    bool isValid(void) const;
    BObolViewController *controller(void) const;
    int bindHost(BObolWindowHost *host, bool takeOwnership = false);
    void detachHost(void);
    BObolWindowHost *host(void) const;
    int setRenderEngine(enum bobol_render_engine engine);
    int supportsRenderEngine(enum bobol_render_engine engine) const;
    uint64_t renderEngineCapabilities(void) const;
    int refreshDiagnostics(void);
    enum bobol_render_engine renderEngine(void) const;
    enum bobol_render_engine resolvedRenderEngine(void) const;
    int propertyGet(const char *name,
	struct bobol_endpoint_property_value *value) const;
    int propertySet(const char *name,
	const struct bobol_endpoint_property_value *value);

    bobol_display_endpoint_t *release(void);
    bobol_display_endpoint_t *get(void) const;

private:
    BObolDisplayEndpoint(const BObolDisplayEndpoint &) = delete;
    BObolDisplayEndpoint &operator=(const BObolDisplayEndpoint &) = delete;

    bobol_display_endpoint_t *endpoint;
};

#endif

#endif /* BOBOL_BDISPLAYENDPOINT_H */
