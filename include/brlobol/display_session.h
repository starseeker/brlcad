/*              D I S P L A Y _ S E S S I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/display_session.h */

#ifndef BRLOBOL_DISPLAY_SESSION_H
#define BRLOBOL_DISPLAY_SESSION_H

#include "brlobol/defines.h"
#include "brlobol/display_endpoint.h"
#include "imgstream/fb_compat.h"

#include <stddef.h>
#include <stdint.h>

#define BRLOBOL_DISPLAY_PROVIDER_ABI_VERSION 1u

struct brlobol_display_provider;
typedef struct brlobol_display_provider brlobol_display_provider_t;
struct brlobol_display_session;
typedef struct brlobol_display_session brlobol_display_session_t;
typedef int (*brlobol_display_session_task_t)(void *data);

/** A toolkit provider opens an endpoint host and owns only its toolkit state.
 * The core session owns the endpoint, imgstream framebuffer, and retained
 * framebuffer presentation. */
struct brlobol_display_provider {
    uint32_t abi_version;
    uint32_t struct_size;
    const char *name;
    int priority;
    void *user_data;
    int (*open)(brlobol_display_endpoint_t *endpoint,
	const imgstream_fb_spec_info_t *spec, size_t width, size_t height,
	const char *title, void **instance, void *user_data);
    void (*close)(void *instance, void *user_data);
    int (*poll)(void *instance, void *user_data);
    long (*poll_rate)(const void *instance, void *user_data);
};

__BEGIN_DECLS

/** Register a process-lifetime toolkit provider.  Registration is idempotent
 * for the same provider name. */
BRLOBOL_EXPORT int brlobol_display_provider_register(
	const brlobol_display_provider_t *provider);

/** Open a legacy display specification through a registered Obol provider.
 * Only IMGSTREAM_FB_SPEC_DISPLAY specifications are accepted. */
BRLOBOL_EXPORT brlobol_display_session_t *brlobol_display_session_open(
	const char *spec, size_t width, size_t height, const char *title);
BRLOBOL_EXPORT void brlobol_display_session_close(
	brlobol_display_session_t *session);
BRLOBOL_EXPORT imgstream_fb_t *brlobol_display_session_framebuffer(
	const brlobol_display_session_t *session);
BRLOBOL_EXPORT int brlobol_display_session_present(
	brlobol_display_session_t *session);
BRLOBOL_EXPORT int brlobol_display_session_poll(
	brlobol_display_session_t *session);
BRLOBOL_EXPORT long brlobol_display_session_poll_rate(
	const brlobol_display_session_t *session);

/** Run work on a background thread while the caller continues to own and poll
 * the native display session.  The task must not access toolkit or retained
 * presentation state. */
BRLOBOL_EXPORT int brlobol_display_session_run(
	brlobol_display_session_t *session,
	brlobol_display_session_task_t task, void *task_data);

__END_DECLS

#endif /* BRLOBOL_DISPLAY_SESSION_H */
