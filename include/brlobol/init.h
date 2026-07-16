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
 * The optional manager initializes Coin/Obol's compatibility process-global
 * context manager.  BRL-CAD view controllers and hosts do not select their
 * rendering provider from that global state: hosts must bind a concrete
 * provider with BRLObolViewController::setRenderContextManager().
 *
 * Passing NULL initializes Coin with libbrlobol's process-lifetime OSMesa
 * manager for compatibility.  It does not implicitly enable controller
 * rendering.
 *
 * If libbrlobol is already initialized, passing a non-NULL manager updates
 * Obol's active context manager without re-registering classes.
 */
BRLOBOL_EXPORT void brlobol_init(SoDB::ContextManager *contextManager);

/** Return libbrlobol's process-lifetime OSMesa context manager. */
BRLOBOL_EXPORT SoDB::ContextManager *brlobol_headless_context_manager(void);

#endif /* BRLOBOL_INIT_H */
