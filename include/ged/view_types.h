/*                 V I E W _ T Y P E S . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file ged/view_types.h
 *
 * Opaque GED view ownership types.  The concrete implementation is private to
 * libged; installed callers must not cast these pointers to renderer or libbv
 * implementation objects.
 */

#ifndef GED_VIEW_TYPES_H
#define GED_VIEW_TYPES_H

#include "bu/defines.h"

#ifndef GED_EXPORT
#  if defined(GED_DLL_EXPORTS) && defined(GED_DLL_IMPORTS)
#    error "Only GED_DLL_EXPORTS or GED_DLL_IMPORTS can be defined, not both."
#  elif defined(GED_DLL_EXPORTS)
#    define GED_EXPORT COMPILER_DLLEXPORT
#  elif defined(GED_DLL_IMPORTS)
#    define GED_EXPORT COMPILER_DLLIMPORT
#  else
#    define GED_EXPORT
#  endif
#endif

__BEGIN_DECLS

struct ged_view_context;
struct ged_view_set;
struct bv_context;
struct bv_lod_policy;

typedef struct ged_view_context ged_view_context;
typedef struct ged_view_set ged_view_set;
typedef struct bv_lod_policy ged_view_lod_policy;

/**
 * Explicit compatibility boundary between GED's owned view context and its
 * borrowed libbv record.  Returned pointers alias the input, are never owned
 * by the caller, and remain valid only for the GED context's lifetime.  These
 * functions do not accept arbitrary renderer or application data.
 */
GED_EXPORT extern struct bv_context *ged_view_context_bv(struct ged_view_context *view);
GED_EXPORT extern const struct bv_context *ged_view_context_bv_const(const struct ged_view_context *view);
GED_EXPORT extern struct ged_view_context *ged_view_context_from_bv(struct bv_context *view);
GED_EXPORT extern const struct ged_view_context *ged_view_context_from_bv_const(const struct bv_context *view);

__END_DECLS

#endif /* GED_VIEW_TYPES_H */
