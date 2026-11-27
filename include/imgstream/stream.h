/*                        S T R E A M . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file imgstream/stream.h */
#ifndef IMGSTREAM_STREAM_H
#define IMGSTREAM_STREAM_H

#include "common.h"
#include <stddef.h>
#include <stdint.h>

#include "icv.h"

__BEGIN_DECLS

#ifndef IMGSTREAM_EXPORT
#  if defined(IMGSTREAM_DLL_EXPORTS) && defined(IMGSTREAM_DLL_IMPORTS)
#    error "Only IMGSTREAM_DLL_EXPORTS or IMGSTREAM_DLL_IMPORTS can be defined, not both."
#  elif defined(IMGSTREAM_DLL_EXPORTS)
#    define IMGSTREAM_EXPORT COMPILER_DLLEXPORT
#  elif defined(IMGSTREAM_DLL_IMPORTS)
#    define IMGSTREAM_EXPORT COMPILER_DLLIMPORT
#  else
#    define IMGSTREAM_EXPORT
#  endif
#endif

struct imgstream;
typedef struct imgstream imgstream_t;

enum imgstream_pixel_format {
    IMGSTREAM_PIXEL_RGB8 = 3,
    IMGSTREAM_PIXEL_RGBA8 = 4
};

struct imgstream_rect {
    size_t x;
    size_t y;
    size_t width;
    size_t height;
};

struct imgstream_info {
    size_t width;
    size_t height;
    enum imgstream_pixel_format format;
    size_t channels;
    size_t stride;
    uint64_t generation;
    int dirty;
    struct imgstream_rect dirty_rect;
    int producer_active;
};

typedef int imgstream_subscriber_id;
typedef void (*imgstream_dirty_cb)(void *ctx, const struct imgstream_rect *rect, uint64_t generation);

IMGSTREAM_EXPORT imgstream_t *imgstream_create(size_t width, size_t height, enum imgstream_pixel_format format);
IMGSTREAM_EXPORT imgstream_t *imgstream_create_from_icv(const icv_image_t *image);
IMGSTREAM_EXPORT void imgstream_destroy(imgstream_t *stream);

IMGSTREAM_EXPORT size_t imgstream_width(const imgstream_t *stream);
IMGSTREAM_EXPORT size_t imgstream_height(const imgstream_t *stream);
IMGSTREAM_EXPORT enum imgstream_pixel_format imgstream_format(const imgstream_t *stream);
IMGSTREAM_EXPORT size_t imgstream_channels(const imgstream_t *stream);
IMGSTREAM_EXPORT uint64_t imgstream_generation(const imgstream_t *stream);
IMGSTREAM_EXPORT int imgstream_producer_active(const imgstream_t *stream);
IMGSTREAM_EXPORT int imgstream_get_info(const imgstream_t *stream, struct imgstream_info *info);

IMGSTREAM_EXPORT int imgstream_resize(imgstream_t *stream, size_t width, size_t height);
IMGSTREAM_EXPORT int imgstream_clear(imgstream_t *stream, const unsigned char *pixel);
IMGSTREAM_EXPORT int imgstream_write(imgstream_t *stream, const unsigned char *pixels, size_t stride);
IMGSTREAM_EXPORT int imgstream_read(const imgstream_t *stream, unsigned char *pixels, size_t stride);
IMGSTREAM_EXPORT int imgstream_write_rect(imgstream_t *stream, size_t x, size_t y, size_t width, size_t height, const unsigned char *pixels, size_t stride);
IMGSTREAM_EXPORT int imgstream_read_rect(const imgstream_t *stream, size_t x, size_t y, size_t width, size_t height, unsigned char *pixels, size_t stride);

IMGSTREAM_EXPORT int imgstream_dirty_rect(const imgstream_t *stream, struct imgstream_rect *rect);
IMGSTREAM_EXPORT void imgstream_dirty_clear(imgstream_t *stream);

IMGSTREAM_EXPORT imgstream_subscriber_id imgstream_subscribe(imgstream_t *stream, imgstream_dirty_cb callback, void *ctx);
IMGSTREAM_EXPORT int imgstream_unsubscribe(imgstream_t *stream, imgstream_subscriber_id id);

IMGSTREAM_EXPORT int imgstream_producer_begin(imgstream_t *stream);
IMGSTREAM_EXPORT int imgstream_producer_end(imgstream_t *stream);

IMGSTREAM_EXPORT icv_image_t *imgstream_snapshot_icv(const imgstream_t *stream);
IMGSTREAM_EXPORT icv_image_t *imgstream_to_icv(const imgstream_t *stream);

__END_DECLS

#endif /* IMGSTREAM_STREAM_H */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
