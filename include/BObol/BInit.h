/*                       B I N I T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BInit.h */

#ifndef BOBOL_BINIT_H
#define BOBOL_BINIT_H

#include "BObol/BDefines.h"

#include <Inventor/SoDB.h>

/**
 * Initialize Obol and register the BRL-CAD Obol node/action/detail classes.
 *
 * The optional manager initializes Coin/Obol's compatibility process-global
 * context manager.  BRL-CAD view controllers and hosts do not select their
 * rendering provider from that global state: hosts must bind a concrete
 * provider with BObolViewController::setRenderContextManager().
 *
 * Passing NULL initializes Coin with libBObol's process-lifetime OSMesa
 * manager for compatibility.  It does not implicitly enable controller
 * rendering.
 *
 * If libBObol is already initialized, passing a non-NULL manager updates
 * Obol's active context manager without re-registering classes.
 */
BOBOL_EXPORT void bobol_init(SoDB::ContextManager *contextManager);

/** Return libBObol's process-lifetime OSMesa context manager. */
BOBOL_EXPORT SoDB::ContextManager *bobol_headless_context_manager(void);

#endif /* BOBOL_BINIT_H */
