/*             Q T C A D / D I S P L A Y _ P R O V I D E R . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file qtcad/display_provider.h */

#ifndef QTCAD_DISPLAY_PROVIDER_H
#define QTCAD_DISPLAY_PROVIDER_H

#include "common.h"
#include "qtcad/defines.h"

__BEGIN_DECLS

/** Register libqtcad's top-level Qt Obol display-session provider.  Clients
 * use libBObol's toolkit-neutral display-session API after this selection. */
QTCAD_EXPORT int qtcad_obol_display_provider_register(void);

__END_DECLS

#endif /* QTCAD_DISPLAY_PROVIDER_H */
