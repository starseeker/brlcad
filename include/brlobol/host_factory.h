/*                  H O S T _ F A C T O R Y . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/host_factory.h */

#ifndef BRLOBOL_HOST_FACTORY_H
#define BRLOBOL_HOST_FACTORY_H

#include "brlobol/defines.h"
#include "brlobol/input.h"

#include <stddef.h>
#include <stdint.h>

#define BRLOBOL_HOST_FACTORY_ABI_VERSION 1u

#define BRLOBOL_HOST_CAP_TOPLEVEL             UINT64_C(0x00000001)
#define BRLOBOL_HOST_CAP_EMBEDDED             UINT64_C(0x00000002)
#define BRLOBOL_HOST_CAP_SYSTEM_GL            UINT64_C(0x00000004)
#define BRLOBOL_HOST_CAP_PIXEL_PRESENT        UINT64_C(0x00000008)
#define BRLOBOL_HOST_CAP_PROGRESSIVE_PRESENT  UINT64_C(0x00000010)
#define BRLOBOL_HOST_CAP_INPUT                UINT64_C(0x00000020)
#define BRLOBOL_HOST_CAP_READBACK             UINT64_C(0x00000040)
#define BRLOBOL_HOST_CAP_FRAMEBUFFER_PRESENT  UINT64_C(0x00000080)
#define BRLOBOL_HOST_CAP_MULTI_VIEW           UINT64_C(0x00000100)
#define BRLOBOL_HOST_CAP_THREAD_AFFINE        UINT64_C(0x00000200)
#define BRLOBOL_HOST_CAP_PRESENT_VSYNC        UINT64_C(0x00000400)

enum brlobol_host_vsync {
    BRLOBOL_HOST_VSYNC_AUTO = 0,
    BRLOBOL_HOST_VSYNC_OFF = 1,
    BRLOBOL_HOST_VSYNC_ON = 2
};

enum brlobol_host_mode {
    BRLOBOL_HOST_MODE_TOPLEVEL = 0,
    BRLOBOL_HOST_MODE_EMBEDDED = 1,
    BRLOBOL_HOST_MODE_HEADLESS = 2,
    BRLOBOL_HOST_MODE_DIAGNOSTIC = 3
};

struct brlobol_host_desc {
    uint32_t struct_size;
    enum brlobol_host_mode mode;
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
    enum brlobol_host_vsync vsync;
    /** Endpoint-owned normalized input dispatcher for this host instance. */
    BRLObolInputEventHandler input_dispatch;
    void *input_dispatch_data;
};

typedef struct brlobol_host_factory_token brlobol_host_factory_token_t;

struct brlobol_host_factory {
    uint32_t abi_version;
    uint32_t struct_size;
    const char *name;
    int priority;
    uint64_t capabilities;
    void *user_data;

    int (*probe)(const struct brlobol_host_desc *desc, void *user_data);
    void *(*create)(const struct brlobol_host_desc *desc, void *user_data);
    void (*destroy)(void *instance, void *user_data);
    int (*bind_controller)(void *instance, void *controller, void *user_data);
    int (*open)(void *instance, const struct brlobol_host_desc *desc,
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
};

__BEGIN_DECLS

/** Register a copied factory descriptor.  Duplicate names are rejected. */
BRLOBOL_EXPORT brlobol_host_factory_token_t *
brlobol_host_factory_register(const struct brlobol_host_factory *factory);

/**
 * Unregister and destroy a registration token.
 *
 * Returns zero while an endpoint or acquired reference still uses the token.
 */
BRLOBOL_EXPORT int
brlobol_host_factory_unregister(brlobol_host_factory_token_t *token);

BRLOBOL_EXPORT size_t brlobol_host_factory_registry_count(void);
BRLOBOL_EXPORT size_t brlobol_host_factory_registry_name(
	size_t index, char *name, size_t name_size);
BRLOBOL_EXPORT uint64_t brlobol_host_factory_registry_capabilities(
	size_t index);
BRLOBOL_EXPORT const char *brlobol_host_factory_name(
	const brlobol_host_factory_token_t *token);
BRLOBOL_EXPORT uint64_t brlobol_host_factory_capabilities(
	const brlobol_host_factory_token_t *token);

/** Acquire a named or highest-priority compatible factory. */
BRLOBOL_EXPORT brlobol_host_factory_token_t *
brlobol_host_factory_acquire(const char *name,
	const struct brlobol_host_desc *desc);
BRLOBOL_EXPORT void brlobol_host_factory_release(
	brlobol_host_factory_token_t *token);

/** Create/open and close/destroy one host instance through a retained token. */
BRLOBOL_EXPORT int brlobol_host_factory_instance_create(
	brlobol_host_factory_token_t *token,
	const struct brlobol_host_desc *desc,
	void *controller, void **instance);
BRLOBOL_EXPORT void brlobol_host_factory_instance_destroy(
	brlobol_host_factory_token_t *token, void *instance);

BRLOBOL_EXPORT int brlobol_host_factory_instance_request_frame(
	brlobol_host_factory_token_t *token, void *instance,
	const char *reason);
BRLOBOL_EXPORT int brlobol_host_factory_instance_resize(
	brlobol_host_factory_token_t *token, void *instance,
	unsigned int width, unsigned int height, double device_pixel_ratio);
BRLOBOL_EXPORT int brlobol_host_factory_instance_capture(
	brlobol_host_factory_token_t *token, void *instance,
	unsigned char **pixels, size_t *size, unsigned int *width,
	unsigned int *height, unsigned int *components);
BRLOBOL_EXPORT int brlobol_host_factory_instance_dimensions(
	brlobol_host_factory_token_t *token, void *instance,
	unsigned int *width, unsigned int *height, double *device_pixel_ratio);
BRLOBOL_EXPORT int brlobol_host_factory_instance_set_title(
	brlobol_host_factory_token_t *token, void *instance, const char *title);
BRLOBOL_EXPORT int brlobol_host_factory_instance_set_visible(
	brlobol_host_factory_token_t *token, void *instance, int visible);
BRLOBOL_EXPORT int brlobol_host_factory_instance_set_vsync(
	brlobol_host_factory_token_t *token, void *instance, int enabled);

__END_DECLS

#endif /* BRLOBOL_HOST_FACTORY_H */
