/*                         I N I T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/init.h */

#ifndef BRLOBOL_INIT_H
#define BRLOBOL_INIT_H

#include "brlobol/defines.h"

#include <Inventor/SoDB.h>

/**
 * Initialize Obol and register the BRL-CAD Obol node/action/detail classes.
 *
 * Rendering applications should pass the application-owned Obol context
 * manager that will service GL or renderScene() backend requests.  The caller
 * owns the manager and must keep it alive for the life of any Obol renderer
 * that can use the global SoDB context manager.
 *
 * Passing NULL installs libbrlobol's no-op fallback manager.  That is useful
 * for non-rendering scene/action tests, but any renderer that relies on the
 * global context manager will fail unless a real manager is supplied.
 *
 * If libbrlobol is already initialized, passing a non-NULL manager updates
 * Obol's active context manager without re-registering classes.
 */
BRLOBOL_EXPORT void brlobol_init(SoDB::ContextManager *contextManager);

/** Return libbrlobol's process-lifetime OSMesa context manager. */
BRLOBOL_EXPORT SoDB::ContextManager *brlobol_headless_context_manager(void);

#endif /* BRLOBOL_INIT_H */
