/*                B H O S T F A C T O R Y . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BHostFactory.h */

#ifndef BOBOL_BHOSTFACTORY_H
#define BOBOL_BHOSTFACTORY_H

#include "BObol/BDefines.h"
#include "BObol/BInput.h"

#include <stddef.h>
#include <stdint.h>

#define BOBOL_HOST_FACTORY_ABI_VERSION 1u

#define BOBOL_HOST_CAP_TOPLEVEL             UINT64_C(0x00000001)
#define BOBOL_HOST_CAP_EMBEDDED             UINT64_C(0x00000002)
#define BOBOL_HOST_CAP_SYSTEM_GL            UINT64_C(0x00000004)
#define BOBOL_HOST_CAP_PIXEL_PRESENT        UINT64_C(0x00000008)
#define BOBOL_HOST_CAP_PROGRESSIVE_PRESENT  UINT64_C(0x00000010)
#define BOBOL_HOST_CAP_INPUT                UINT64_C(0x00000020)
#define BOBOL_HOST_CAP_READBACK             UINT64_C(0x00000040)
#define BOBOL_HOST_CAP_FRAMEBUFFER_PRESENT  UINT64_C(0x00000080)
#define BOBOL_HOST_CAP_MULTI_VIEW           UINT64_C(0x00000100)
#define BOBOL_HOST_CAP_THREAD_AFFINE        UINT64_C(0x00000200)
#define BOBOL_HOST_CAP_PRESENT_VSYNC        UINT64_C(0x00000400)

enum bobol_host_vsync {
    BOBOL_HOST_VSYNC_AUTO = 0,
    BOBOL_HOST_VSYNC_OFF = 1,
    BOBOL_HOST_VSYNC_ON = 2
};

enum bobol_host_mode {
    BOBOL_HOST_MODE_TOPLEVEL = 0,
    BOBOL_HOST_MODE_EMBEDDED = 1,
    BOBOL_HOST_MODE_HEADLESS = 2,
    BOBOL_HOST_MODE_DIAGNOSTIC = 3
};

struct bobol_host_desc {
    uint32_t struct_size;
    enum bobol_host_mode mode;
    unsigned int width;
    unsigned int height;
    double device_pixel_ratio;
    int visible;
    uint64_t required_capabilities;
    const char *title;
    const char *display;
    const char *native_id_hint;
    /** Factory-defined toolkit/application context.  Embedded factories use
     * this to receive the native canvas or interpreter they must borrow. */
    void *application_context;
    enum bobol_host_vsync vsync;
    /** Endpoint-owned normalized input dispatcher for this host instance. */
    BObolInputEventHandler input_dispatch;
    void *input_dispatch_data;
};

typedef struct bobol_host_factory_token bobol_host_factory_token_t;

struct bobol_host_factory {
    uint32_t abi_version;
    uint32_t struct_size;
    const char *name;
    int priority;
    uint64_t capabilities;
    void *user_data;

    int (*probe)(const struct bobol_host_desc *desc, void *user_data);
    void *(*create)(const struct bobol_host_desc *desc, void *user_data);
    void (*destroy)(void *instance, void *user_data);
    /** Bind the controller and its concrete rendering provider.  Factories
     * advertising GL, pixel, or progressive presentation must install a
     * non-NULL provider before open succeeds. */
    int (*bind_controller)(void *instance, void *controller, void *user_data);
    int (*open)(void *instance, const struct bobol_host_desc *desc,
	void *user_data);
    void (*close)(void *instance, void *user_data);
    int (*request_frame)(void *instance, const char *reason, void *user_data);
    int (*resize)(void *instance, unsigned int width, unsigned int height,
	double device_pixel_ratio, void *user_data);
    int (*capture)(void *instance, unsigned char **pixels, size_t *size,
	unsigned int *width, unsigned int *height, unsigned int *components,
	void *user_data);
    int (*dimensions)(void *instance, unsigned int *width,
	unsigned int *height, double *device_pixel_ratio, void *user_data);
    int (*set_title)(void *instance, const char *title, void *user_data);
    int (*set_visible)(void *instance, int visible, void *user_data);
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

BOBOL_EXPORT size_t bobol_host_factory_registry_count(void);
BOBOL_EXPORT size_t bobol_host_factory_registry_name(
	size_t index, char *name, size_t name_size);
BOBOL_EXPORT uint64_t bobol_host_factory_registry_capabilities(
	size_t index);
BOBOL_EXPORT const char *bobol_host_factory_name(
	const bobol_host_factory_token_t *token);
BOBOL_EXPORT uint64_t bobol_host_factory_capabilities(
	const bobol_host_factory_token_t *token);

/** Acquire a named or highest-priority compatible factory. */
BOBOL_EXPORT bobol_host_factory_token_t *
bobol_host_factory_acquire(const char *name,
	const struct bobol_host_desc *desc);
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
BOBOL_EXPORT void bobol_host_factory_instance_destroy(
	bobol_host_factory_token_t *token, void *instance);

BOBOL_EXPORT int bobol_host_factory_instance_request_frame(
	bobol_host_factory_token_t *token, void *instance,
	const char *reason);
BOBOL_EXPORT int bobol_host_factory_instance_resize(
	bobol_host_factory_token_t *token, void *instance,
	unsigned int width, unsigned int height, double device_pixel_ratio);
BOBOL_EXPORT int bobol_host_factory_instance_capture(
	bobol_host_factory_token_t *token, void *instance,
	unsigned char **pixels, size_t *size, unsigned int *width,
	unsigned int *height, unsigned int *components);
BOBOL_EXPORT int bobol_host_factory_instance_dimensions(
	bobol_host_factory_token_t *token, void *instance,
	unsigned int *width, unsigned int *height, double *device_pixel_ratio);
BOBOL_EXPORT int bobol_host_factory_instance_set_title(
	bobol_host_factory_token_t *token, void *instance, const char *title);
BOBOL_EXPORT int bobol_host_factory_instance_set_visible(
	bobol_host_factory_token_t *token, void *instance, int visible);
BOBOL_EXPORT int bobol_host_factory_instance_set_vsync(
	bobol_host_factory_token_t *token, void *instance, int enabled);
/** Return a factory-provided borrowed BObolWindowHost, if supported. */
BOBOL_EXPORT void *bobol_host_factory_instance_framebuffer_window_host(
	bobol_host_factory_token_t *token, void *instance);

__END_DECLS

#endif /* BOBOL_BHOSTFACTORY_H */
