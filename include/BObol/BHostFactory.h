/*                B H O S T F A C T O R Y . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BHostFactory.h
 * @ingroup bobol_display
 *
 * Stable C plug-in contract for endpoint host factories.  Registry functions
 * are synchronized; host-instance calls and every factory callback run on the
 * owning endpoint thread.  Descriptors are copied at registration, while host
 * instances and toolkit objects remain owned by their factory.
 */

#ifndef BOBOL_BHOSTFACTORY_H
#define BOBOL_BHOSTFACTORY_H

#include "BObol/BDefines.h"
#include "BObol/BInput.h"

#include <stddef.h>
#include <stdint.h>

/** @addtogroup bobol_display
 * @{ */

/** bobol_host_factory descriptor ABI implemented by this release. */
#define BOBOL_HOST_FACTORY_ABI_VERSION 1u

/** Factory can create independent top-level windows. */
#define BOBOL_HOST_CAP_TOPLEVEL             UINT64_C(0x00000001)
/** Factory can bind a caller-supplied native canvas. */
#define BOBOL_HOST_CAP_EMBEDDED             UINT64_C(0x00000002)
/** Factory offers system/hardware OpenGL presentation. */
#define BOBOL_HOST_CAP_SYSTEM_GL            UINT64_C(0x00000004)
/** Factory accepts software-rendered pixel frames. */
#define BOBOL_HOST_CAP_PIXEL_PRESENT        UINT64_C(0x00000008)
/** Factory can present progressive updates. */
#define BOBOL_HOST_CAP_PROGRESSIVE_PRESENT  UINT64_C(0x00000010)
/** Factory forwards normalized input events. */
#define BOBOL_HOST_CAP_INPUT                UINT64_C(0x00000020)
/** Factory can capture its displayed pixels. */
#define BOBOL_HOST_CAP_READBACK             UINT64_C(0x00000040)
/** Factory accepts retained framebuffer composition. */
#define BOBOL_HOST_CAP_FRAMEBUFFER_PRESENT  UINT64_C(0x00000080)
/** Factory may host multiple independently bound views. */
#define BOBOL_HOST_CAP_MULTI_VIEW           UINT64_C(0x00000100)
/** Factory and instances must be called on their creation thread. */
#define BOBOL_HOST_CAP_THREAD_AFFINE        UINT64_C(0x00000200)
/** Factory supports explicit presentation-vsync policy. */
#define BOBOL_HOST_CAP_PRESENT_VSYNC        UINT64_C(0x00000400)

/** Requested presentation-vsync policy. */
enum bobol_host_vsync {
    BOBOL_HOST_VSYNC_AUTO = 0,
    BOBOL_HOST_VSYNC_OFF = 1,
    BOBOL_HOST_VSYNC_ON = 2
};

/** Native hosting mode requested from a factory. */
enum bobol_host_mode {
    BOBOL_HOST_MODE_TOPLEVEL = 0,
    BOBOL_HOST_MODE_EMBEDDED = 1,
    BOBOL_HOST_MODE_HEADLESS = 2,
    BOBOL_HOST_MODE_DIAGNOSTIC = 3
};

/** Immutable host request borrowed by probe/create/open callbacks.  Strings,
 * application context, and input callback storage need remain valid only for
 * the synchronous call unless a factory explicitly copies them. */
struct bobol_host_desc {
    uint32_t struct_size; /**< Set to sizeof(struct bobol_host_desc). */
    enum bobol_host_mode mode; /**< Required hosting mode. */
    unsigned int width; /**< Initial physical width in pixels. */
    unsigned int height; /**< Initial physical height in pixels. */
    double device_pixel_ratio; /**< Logical-to-physical pixel ratio; positive. */
    int visible; /**< Nonzero requests initial visibility. */
    uint64_t required_capabilities; /**< Required BOBOL_HOST_CAP_* bit mask. */
    const char *title; /**< Borrowed UTF-8 title, or NULL. */
    const char *display; /**< Borrowed platform display name, or NULL. */
    const char *native_id_hint; /**< Borrowed platform-native identifier hint. */
    /** Factory-defined toolkit/application context.  Embedded factories use
     * this to receive the native canvas or interpreter they must borrow. */
    void *application_context;
    enum bobol_host_vsync vsync; /**< Initial vsync policy. */
    /** Endpoint-owned normalized input dispatcher for this host instance. */
    BObolInputEventHandler input_dispatch;
    void *input_dispatch_data; /**< Opaque first argument for input_dispatch. */
};

/** Opaque retained registry reference.  Acquire/register returns one owned
 * reference and bobol_host_factory_release/unregister consumes it. */
typedef struct bobol_host_factory_token bobol_host_factory_token_t;

/** Versioned factory callback table copied at registration.  All callbacks
 * are synchronous owner-thread calls.  Callback return values use nonzero for
 * success and zero for unsupported/failure; create returns NULL on failure. */
struct bobol_host_factory {
    uint32_t abi_version; /**< BOBOL_HOST_FACTORY_ABI_VERSION. */
    uint32_t struct_size; /**< sizeof(struct bobol_host_factory). */
    const char *name; /**< Unique factory name copied by registration. */
    int priority; /**< Selection priority; larger values are preferred. */
    uint64_t capabilities; /**< BOBOL_HOST_CAP_* mask. */
    void *user_data; /**< Factory-owned callback context. */

    /** Test a request without allocating; nonzero means compatible. */
    int (*probe)(const struct bobol_host_desc *desc, void *user_data);
    /** Create a closed factory-owned instance, or NULL on failure. */
    void *(*create)(const struct bobol_host_desc *desc, void *user_data);
    /** Destroy a closed instance. */
    void (*destroy)(void *instance, void *user_data);
    /** Bind the controller and its concrete rendering provider.  Factories
     * advertising GL, pixel, or progressive presentation must install a
     * non-NULL provider before open succeeds. */
    int (*bind_controller)(void *instance, void *controller, void *user_data);
    /** Open native resources for an existing instance. */
    int (*open)(void *instance, const struct bobol_host_desc *desc,
	void *user_data);
    /** Close native resources without destroying the instance. */
    void (*close)(void *instance, void *user_data);
    /** Queue or perform one presentation; reason is borrowed. */
    int (*request_frame)(void *instance, const char *reason, void *user_data);
    /** Resize to physical pixels and a positive pixel ratio. */
    int (*resize)(void *instance, unsigned int width, unsigned int height,
	double device_pixel_ratio, void *user_data);
    /** Return bottom-left-order RGB/RGBA pixels in bu_free-compatible storage. */
    int (*capture)(void *instance, unsigned char **pixels, size_t *size,
	unsigned int *width, unsigned int *height, unsigned int *components,
	void *user_data);
    /** Query physical dimensions and logical-to-physical pixel ratio. */
    int (*dimensions)(void *instance, unsigned int *width,
	unsigned int *height, double *device_pixel_ratio, void *user_data);
    /** Set a borrowed UTF-8 title. */
    int (*set_title)(void *instance, const char *title, void *user_data);
    /** Change native visibility. */
    int (*set_visible)(void *instance, int visible, void *user_data);
    /** Enable or disable presentation synchronization. */
    int (*set_vsync)(void *instance, int enabled, void *user_data);
    /** Optional.  Return a borrowed BObolWindowHost suitable for retained
     * framebuffer attachment, or NULL when the factory host is opaque. */
    void *(*framebuffer_window_host)(void *instance, void *user_data);
};

__BEGIN_DECLS

/** Register a copied factory descriptor.  Duplicate names are rejected. */
BOBOL_EXPORT bobol_host_factory_token_t *
bobol_host_factory_register(const struct bobol_host_factory *factory);

/**
 * Unregister and destroy a registration token.
 *
 * Returns zero while an endpoint or acquired reference still uses the token.
 */
BOBOL_EXPORT int
bobol_host_factory_unregister(bobol_host_factory_token_t *token);

/** Return a snapshot count of registered factories. */
BOBOL_EXPORT size_t bobol_host_factory_registry_count(void);
/** Copy registry entry @p index name including a terminating NUL when space
 * permits; return required bytes, or zero for an invalid index. */
BOBOL_EXPORT size_t bobol_host_factory_registry_name(
	size_t index, char *name, size_t name_size);
/** Return BOBOL_HOST_CAP_* bits for registry entry @p index, or zero. */
BOBOL_EXPORT uint64_t bobol_host_factory_registry_capabilities(
	size_t index);
/** Return a borrowed process-lifetime name for a retained token, or NULL. */
BOBOL_EXPORT const char *bobol_host_factory_name(
	const bobol_host_factory_token_t *token);
/** Return BOBOL_HOST_CAP_* bits for a retained token, or zero. */
BOBOL_EXPORT uint64_t bobol_host_factory_capabilities(
	const bobol_host_factory_token_t *token);

/** Acquire a named or highest-priority compatible factory. */
BOBOL_EXPORT bobol_host_factory_token_t *
bobol_host_factory_acquire(const char *name,
	const struct bobol_host_desc *desc);
/** Release one acquired token reference.  NULL is accepted. */
BOBOL_EXPORT void bobol_host_factory_release(
	bobol_host_factory_token_t *token);

/** Create/open and close/destroy one host instance through a retained token. */
BOBOL_EXPORT int bobol_host_factory_instance_create(
	bobol_host_factory_token_t *token,
	const struct bobol_host_desc *desc,
	void *controller, void **instance);
/** Explicitly replace or reassert an instance's controller binding. */
BOBOL_EXPORT int bobol_host_factory_instance_bind_controller(
	bobol_host_factory_token_t *token, void *instance, void *controller);
/** Close and destroy a factory instance; token remains retained by caller. */
BOBOL_EXPORT void bobol_host_factory_instance_destroy(
	bobol_host_factory_token_t *token, void *instance);

/** Request one frame from an open instance. */
BOBOL_EXPORT int bobol_host_factory_instance_request_frame(
	bobol_host_factory_token_t *token, void *instance,
	const char *reason);
/** Resize an open instance in physical pixels. */
BOBOL_EXPORT int bobol_host_factory_instance_resize(
	bobol_host_factory_token_t *token, void *instance,
	unsigned int width, unsigned int height, double device_pixel_ratio);
/** Capture bottom-left-order RGB/RGBA pixels into bu_free-compatible storage. */
BOBOL_EXPORT int bobol_host_factory_instance_capture(
	bobol_host_factory_token_t *token, void *instance,
	unsigned char **pixels, size_t *size, unsigned int *width,
	unsigned int *height, unsigned int *components);
/** Query physical dimensions and device pixel ratio. */
BOBOL_EXPORT int bobol_host_factory_instance_dimensions(
	bobol_host_factory_token_t *token, void *instance,
	unsigned int *width, unsigned int *height, double *device_pixel_ratio);
/** Set an open instance's UTF-8 title. */
BOBOL_EXPORT int bobol_host_factory_instance_set_title(
	bobol_host_factory_token_t *token, void *instance, const char *title);
/** Set an open instance's native visibility. */
BOBOL_EXPORT int bobol_host_factory_instance_set_visible(
	bobol_host_factory_token_t *token, void *instance, int visible);
/** Set presentation-vsync enabled state. */
BOBOL_EXPORT int bobol_host_factory_instance_set_vsync(
	bobol_host_factory_token_t *token, void *instance, int enabled);
/** Return a factory-provided borrowed BObolWindowHost, if supported. */
BOBOL_EXPORT void *bobol_host_factory_instance_framebuffer_window_host(
	bobol_host_factory_token_t *token, void *instance);

__END_DECLS

/** @} */

#endif /* BOBOL_BHOSTFACTORY_H */
