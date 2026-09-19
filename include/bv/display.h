/*                    D I S P L A Y . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file bv/display.h
 *
 * Renderer-neutral display property vocabulary.
 *
 * This value ABI is shared by view owners and retained-render backends.  It
 * deliberately contains no renderer, scene-graph, window-system, or toolkit
 * handles.  Property names define the semantic contract; endpoint
 * implementations only transport the typed values.
 */

#ifndef BV_DISPLAY_H
#define BV_DISPLAY_H

#include "common.h"

#include <stdint.h>

__BEGIN_DECLS

/** Discriminator for one value in a typed display property registry. */
enum bv_display_property_type {
    BV_DISPLAY_PROPERTY_BOOL = 1,
    BV_DISPLAY_PROPERTY_INT = 2,
    BV_DISPLAY_PROPERTY_UINT = 3,
    BV_DISPLAY_PROPERTY_DOUBLE = 4,
    BV_DISPLAY_PROPERTY_STRING = 5,
    BV_DISPLAY_PROPERTY_COLOR3 = 6,
    BV_DISPLAY_PROPERTY_ENUM = 7
};

/** Property descriptor access bit permitting reads. */
#define BV_DISPLAY_PROPERTY_READ  0x01u
/** Property descriptor access bit permitting writes. */
#define BV_DISPLAY_PROPERTY_WRITE 0x02u

/** Complete result vocabulary for typed display-property operations. */
enum bv_display_property_result {
    BV_DISPLAY_PROPERTY_OK = 1,
    BV_DISPLAY_PROPERTY_UNKNOWN = 0,
    BV_DISPLAY_PROPERTY_INVALID = -1,
    BV_DISPLAY_PROPERTY_READ_ONLY = -2,
    BV_DISPLAY_PROPERTY_UNSUPPORTED = -3
};

/** Static metadata for one registered display property.  Returned strings are
 * borrowed process-lifetime registry storage. */
struct bv_display_property_desc {
    uint32_t struct_size; /**< Caller sets this to sizeof the structure. */
    const char *name; /**< Borrowed canonical property name. */
    enum bv_display_property_type type; /**< Required value discriminator. */
    unsigned int access; /**< BV_DISPLAY_PROPERTY_READ/WRITE bits. */
    uint64_t required_host_capabilities; /**< Backend-defined capability bits. */
    double minimum; /**< Inclusive numeric minimum when applicable. */
    double maximum; /**< Inclusive numeric maximum when applicable. */
    const char *allowed_values; /**< Borrowed enum spelling list, or NULL. */
};

/** Tagged display property value.  Only the field selected by @ref type is
 * meaningful.  Returned strings are borrowed until the next endpoint
 * mutation; setter strings need remain valid only during the call. */
struct bv_display_property_value {
    uint32_t struct_size; /**< Caller sets this to sizeof the structure. */
    enum bv_display_property_type type; /**< Active payload discriminator. */
    int bool_value; /**< Boolean payload, normalized to zero or one. */
    int64_t int_value; /**< Signed integer or enum payload. */
    uint64_t uint_value; /**< Unsigned integer payload. */
    double double_value; /**< Floating-point payload in property-defined units. */
    double color3[3]; /**< Normalized RGB channels in the range [0,1]. */
    const char *string_value; /**< Borrowed/setter-input UTF-8 string payload. */
};

/** External-property getter.  It runs synchronously on the display owner
 * thread and returns bv_display_property_result. */
typedef int (*bv_display_property_get_callback)(void *user_data,
	const char *name, struct bv_display_property_value *out);

/** External-property setter.  The value is borrowed for the synchronous
 * owner-thread callback and returns bv_display_property_result. */
typedef int (*bv_display_property_set_callback)(void *user_data,
	const char *name, const struct bv_display_property_value *value);

/** Initializer for a bool-valued bv_display_property_value. */
#define BV_DISPLAY_PROPERTY_VALUE_INIT { \
    sizeof(struct bv_display_property_value), \
    BV_DISPLAY_PROPERTY_BOOL, 0, 0, 0, 0.0, {0.0, 0.0, 0.0}, NULL }

__END_DECLS

#endif /* BV_DISPLAY_H */
