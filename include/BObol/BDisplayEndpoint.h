/*            B D I S P L A Y E N D P O I N T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BDisplayEndpoint.h
 * @ingroup bobol_display
 *
 * Stable C ABI for the sole retained-render attachment boundary.  Unless a
 * declaration explicitly says otherwise, endpoint and host operations run on
 * the endpoint owner thread, pointer results are borrowed, and integer
 * operations return nonzero on success and zero on failure.  Callbacks are
 * synchronous and must not destroy or recursively mutate their endpoint.
 */

#ifndef BOBOL_BDISPLAYENDPOINT_H
#define BOBOL_BDISPLAYENDPOINT_H

#include "BObol/BDefines.h"
#include "BObol/BHostFactory.h"
#include "BObol/BInput.h"
#include "bv/display.h"

/** @addtogroup bobol_display
 * @{ */

struct bobol_display_endpoint;
/** Opaque endpoint handle.  It owns or borrows its controller and host as
 * selected by creation/bind flags and is released only with
 * bobol_display_endpoint_destroy(). */
typedef struct bobol_display_endpoint bobol_display_endpoint_t;

/** Retained rendering policy selected for one endpoint. */
enum bobol_render_engine {
    BOBOL_RENDER_ENGINE_AUTO = 0,
    BOBOL_RENDER_ENGINE_HW = 1,
    BOBOL_RENDER_ENGINE_SW = 2,
    BOBOL_RENDER_ENGINE_RT = 3,
    BOBOL_RENDER_ENGINE_NONE = 4,
    BOBOL_RENDER_ENGINE_DIAGNOSTIC = 5
};

/** Capability bit for automatic renderer selection. */
#define BOBOL_RENDER_ENGINE_CAP_AUTO       UINT64_C(0x00000001)
/** Capability bit for system/hardware OpenGL rendering. */
#define BOBOL_RENDER_ENGINE_CAP_HW         UINT64_C(0x00000002)
/** Capability bit for software OpenGL rendering. */
#define BOBOL_RENDER_ENGINE_CAP_SW         UINT64_C(0x00000004)
/** Capability bit for retained ray-traced rendering. */
#define BOBOL_RENDER_ENGINE_CAP_RT         UINT64_C(0x00000008)
/** Capability bit for an endpoint with presentation disabled. */
#define BOBOL_RENDER_ENGINE_CAP_NONE       UINT64_C(0x00000010)
/** Capability bit for non-graphical diagnostic traversal. */
#define BOBOL_RENDER_ENGINE_CAP_DIAGNOSTIC UINT64_C(0x00000020)

/** Pixel composition plane requested from endpoint capture. */
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
    uint32_t struct_size; /**< Set to sizeof(struct bobol_rt_source_identity). */
    void *database; /**< Borrowed database identity; never dereference from C. */
    char *instance_key; /**< Owned instance key freed by the clear function. */
    char *path; /**< Owned database path freed by the clear function. */
    uint32_t source_revision; /**< Revision represented by the captured sample. */
};

/** Initializer for bobol_rt_source_identity before its first query. */
#define BOBOL_RT_SOURCE_IDENTITY_INIT { \
    sizeof(struct bobol_rt_source_identity), NULL, NULL, NULL, 0u }

/** Synchronous framebuffer-plane provider.  On success it transfers a
 * bu_free-compatible pixel allocation to the caller and fills byte size,
 * dimensions, and component count. */
typedef int (*bobol_endpoint_framebuffer_capture_callback)(void *user_data,
	unsigned char **pixels, size_t *size, unsigned int *width,
	unsigned int *height, unsigned int *components);

/** Creation flag transferring controller ownership to an endpoint. */
#define BOBOL_ENDPOINT_OWN_CONTROLLER 0x01u
/** Bind flag transferring host ownership to an endpoint. */
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

/** Destroy an endpoint, its owned controller/host, and all callbacks.  NULL is
 * accepted.  Borrowed objects are detached but not destroyed. */
BOBOL_EXPORT void
bobol_display_endpoint_destroy(bobol_display_endpoint_t *endpoint);

/** Return the borrowed BObolViewController as an opaque pointer, or NULL. */
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

/** Close and detach the current host, destroying it only when endpoint-owned. */
BOBOL_EXPORT void
bobol_display_endpoint_host_detach(bobol_display_endpoint_t *endpoint);

/** Return the borrowed opaque BObolWindowHost, or NULL when unbound. */
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

/** Return the borrowed active factory name, or NULL when no factory host is
 * bound.  The pointer remains valid until host detach/replacement. */
BOBOL_EXPORT const char *
bobol_display_endpoint_host_factory_name(
	const bobol_display_endpoint_t *endpoint);

/** Return capabilities of the active factory host, or zero if unbound. */
BOBOL_EXPORT uint64_t bobol_display_endpoint_host_capabilities(
	const bobol_display_endpoint_t *endpoint);

/** Schedule a frame for an open host.  @p reason is borrowed for the call. */
BOBOL_EXPORT int
bobol_display_endpoint_request_frame(bobol_display_endpoint_t *endpoint,
	const char *reason);

/** Resize endpoint viewport and host.  Dimensions are physical pixels and
 * device_pixel_ratio must be positive. */
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
 * Atomically add or replace a scoped, application-defined action layer.
 * Multiple owners may coexist.  Matching layers are tried by priority until
 * one consumes the event; equal-score layer bindings precede the default
 * profile and newer layers precede older layers.  Action IDs are opaque to
 * libBObol and local to their handler.  Only the owner token may replace or
 * clear its layer.
 */
BOBOL_EXPORT int
bobol_display_endpoint_input_action_layer_set(
	bobol_display_endpoint_t *endpoint,
	const BObolInputActionLayer *layer, void *owner, void *user_data);

/** Remove an installed action layer only if @p owner is the current identity.
 * Returns nonzero when cleared and zero for mismatch or absence. */
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

/** Clear and release allocations in an RT source identity.  The structure is
 * reset and may be reused; NULL is accepted. */
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

/** Select an explicit engine policy.  The operation fails without changing
 * policy when the current or prospective host cannot support it. */
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

/** Return the selected engine policy, or BOBOL_RENDER_ENGINE_AUTO for a NULL
 * endpoint. */
BOBOL_EXPORT enum bobol_render_engine
bobol_display_endpoint_render_engine_get(
	const bobol_display_endpoint_t *endpoint);

/** Return the renderer selected by the endpoint policy and active factory.
 * AUTO resolves to HW before SW for a capable factory host.  It remains AUTO
 * while unbound or when a directly bound custom host has no factory capability
 * record.
 * Retained RT is never selected implicitly. */
BOBOL_EXPORT enum bobol_render_engine
bobol_display_endpoint_render_engine_resolved_get(
	const bobol_display_endpoint_t *endpoint);

/** Enumerate the stable typed endpoint property registry. */
BOBOL_EXPORT size_t bobol_display_endpoint_property_count(void);
/** Copy descriptor @p index into a correctly sized output structure.  Returns
 * bv_display_property_result. */
BOBOL_EXPORT int bobol_display_endpoint_property_descriptor(
	size_t index, struct bv_display_property_desc *out);

/** Get or set one typed property.  Returns bv_display_property_result. */
BOBOL_EXPORT int bobol_display_endpoint_property_get(
	const bobol_display_endpoint_t *endpoint, const char *name,
	struct bv_display_property_value *out);
/** Set one typed property; the input and any string payload are borrowed only
 * for the call.  Returns bv_display_property_result. */
BOBOL_EXPORT int bobol_display_endpoint_property_set(
	bobol_display_endpoint_t *endpoint, const char *name,
	const struct bv_display_property_value *value);

/** Bind optional storage callbacks for properties owned outside libBObol. */
BOBOL_EXPORT int bobol_display_endpoint_property_provider_set(
	bobol_display_endpoint_t *endpoint,
	bv_display_property_get_callback get_callback,
	bv_display_property_set_callback set_callback,
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
	struct bv_display_property_value *value) const;
    int propertySet(const char *name,
	const struct bv_display_property_value *value);

    bobol_display_endpoint_t *release(void);
    bobol_display_endpoint_t *get(void) const;

private:
    BObolDisplayEndpoint(const BObolDisplayEndpoint &) = delete;
    BObolDisplayEndpoint &operator=(const BObolDisplayEndpoint &) = delete;

    bobol_display_endpoint_t *endpoint;
};

#endif

/** @} */

#endif /* BOBOL_BDISPLAYENDPOINT_H */
