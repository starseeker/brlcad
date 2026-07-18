/*             B D I S P L A Y S E S S I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BDisplaySession.h
 * @ingroup bobol_display
 *
 * Stable C convenience API for polling a native display while non-toolkit
 * work runs separately.  Session and provider callbacks are owner-thread
 * operations.  Provider registration is process-lifetime and synchronized.
 */

#ifndef BOBOL_BDISPLAYSESSION_H
#define BOBOL_BDISPLAYSESSION_H

#include "BObol/BDefines.h"
#include "BObol/BDisplayEndpoint.h"
#include "imgstream/fb_compat.h"

#include <stddef.h>
#include <stdint.h>

/** @addtogroup bobol_display
 * @{ */

/** bobol_display_provider descriptor ABI implemented by this release. */
#define BOBOL_DISPLAY_PROVIDER_ABI_VERSION 1u

struct bobol_display_provider;
/** Immutable provider descriptor copied into the process registry. */
typedef struct bobol_display_provider bobol_display_provider_t;
struct bobol_display_session;
/** Opaque owned display session closed by bobol_display_session_close(). */
typedef struct bobol_display_session bobol_display_session_t;
/** Background task called once with borrowed caller data.  It must not access
 * toolkit, endpoint, framebuffer, or retained scene state; return zero for
 * failure and nonzero for success. */
typedef int (*bobol_display_session_task_t)(void *data);

/** A toolkit provider opens an endpoint host and owns only its toolkit state.
 * The core session owns the endpoint, imgstream framebuffer, and retained
 * framebuffer presentation. */
struct bobol_display_provider {
    uint32_t abi_version; /**< BOBOL_DISPLAY_PROVIDER_ABI_VERSION. */
    uint32_t struct_size; /**< sizeof(struct bobol_display_provider). */
    const char *name; /**< Unique provider name copied at registration. */
    int priority; /**< Selection priority; larger values are preferred. */
    void *user_data; /**< Provider-owned callback context. */
    /** Open toolkit state and return the provider-owned instance through
     * @p instance.  All inputs are borrowed for the call. */
    int (*open)(bobol_display_endpoint_t *endpoint,
	const imgstream_fb_spec_info_t *spec, size_t width, size_t height,
	const char *title, void **instance, void *user_data);
    /** Close and destroy a provider-owned instance. */
    void (*close)(void *instance, void *user_data);
    /** Process pending toolkit events; return nonzero while usable. */
    int (*poll)(void *instance, void *user_data);
    /** Return recommended poll interval in microseconds, or a negative value
     * when the provider supplies no recommendation. */
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
/** Close a session and its endpoint, framebuffer, and provider state.  NULL is
 * accepted. */
BOBOL_EXPORT void bobol_display_session_close(
	bobol_display_session_t *session);
/** Return the borrowed imgstream framebuffer owned by @p session, or NULL. */
BOBOL_EXPORT imgstream_fb_t *bobol_display_session_framebuffer(
	const bobol_display_session_t *session);
/** Publish pending framebuffer pixels and request presentation. */
BOBOL_EXPORT int bobol_display_session_present(
	bobol_display_session_t *session);
/** Process one owner-thread toolkit event batch. */
BOBOL_EXPORT int bobol_display_session_poll(
	bobol_display_session_t *session);
/** Return the provider poll interval in microseconds, or a negative value when
 * unavailable. */
BOBOL_EXPORT long bobol_display_session_poll_rate(
	const bobol_display_session_t *session);

/** Run work on a background thread while the caller continues to own and poll
 * the native display session.  The task must not access toolkit or retained
 * presentation state. */
BOBOL_EXPORT int bobol_display_session_run(
	bobol_display_session_t *session,
	bobol_display_session_task_t task, void *task_data);

__END_DECLS

/** @} */

#endif /* BOBOL_BDISPLAYSESSION_H */
