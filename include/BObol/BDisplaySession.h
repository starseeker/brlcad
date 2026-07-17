/*             B D I S P L A Y S E S S I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BDisplaySession.h */

#ifndef BOBOL_BDISPLAYSESSION_H
#define BOBOL_BDISPLAYSESSION_H

#include "BObol/BDefines.h"
#include "BObol/BDisplayEndpoint.h"
#include "imgstream/fb_compat.h"

#include <stddef.h>
#include <stdint.h>

#define BOBOL_DISPLAY_PROVIDER_ABI_VERSION 1u

struct bobol_display_provider;
typedef struct bobol_display_provider bobol_display_provider_t;
struct bobol_display_session;
typedef struct bobol_display_session bobol_display_session_t;
typedef int (*bobol_display_session_task_t)(void *data);

/** A toolkit provider opens an endpoint host and owns only its toolkit state.
 * The core session owns the endpoint, imgstream framebuffer, and retained
 * framebuffer presentation. */
struct bobol_display_provider {
    uint32_t abi_version;
    uint32_t struct_size;
    const char *name;
    int priority;
    void *user_data;
    int (*open)(bobol_display_endpoint_t *endpoint,
	const imgstream_fb_spec_info_t *spec, size_t width, size_t height,
	const char *title, void **instance, void *user_data);
    void (*close)(void *instance, void *user_data);
    int (*poll)(void *instance, void *user_data);
    long (*poll_rate)(const void *instance, void *user_data);
};

__BEGIN_DECLS

/** Register a process-lifetime toolkit provider.  Registration is idempotent
 * for the same provider name. */
BOBOL_EXPORT int bobol_display_provider_register(
	const bobol_display_provider_t *provider);

/** Open a legacy display specification through a registered Obol provider.
 * Only IMGSTREAM_FB_SPEC_DISPLAY specifications are accepted. */
BOBOL_EXPORT bobol_display_session_t *bobol_display_session_open(
	const char *spec, size_t width, size_t height, const char *title);
BOBOL_EXPORT void bobol_display_session_close(
	bobol_display_session_t *session);
BOBOL_EXPORT imgstream_fb_t *bobol_display_session_framebuffer(
	const bobol_display_session_t *session);
BOBOL_EXPORT int bobol_display_session_present(
	bobol_display_session_t *session);
BOBOL_EXPORT int bobol_display_session_poll(
	bobol_display_session_t *session);
BOBOL_EXPORT long bobol_display_session_poll_rate(
	const bobol_display_session_t *session);

/** Run work on a background thread while the caller continues to own and poll
 * the native display session.  The task must not access toolkit or retained
 * presentation state. */
BOBOL_EXPORT int bobol_display_session_run(
	bobol_display_session_t *session,
	bobol_display_session_task_t task, void *task_data);

__END_DECLS

#endif /* BOBOL_BDISPLAYSESSION_H */
