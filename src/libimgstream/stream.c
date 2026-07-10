/*                         S T R E A M . C
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
/** @file libimgstream/stream.c */

#include "common.h"

#include <string.h>

#include "bu/malloc.h"
#include "imgstream.h"


struct imgstream_subscriber {
    imgstream_subscriber_id id;
    imgstream_dirty_cb callback;
    void *ctx;
};

struct imgstream {
    size_t width;
    size_t height;
    enum imgstream_pixel_format format;
    size_t channels;
    unsigned char *pixels;
    uint64_t generation;
    int dirty;
    struct imgstream_rect dirty_rect;
    struct imgstream_subscriber *subscribers;
    size_t subscriber_count;
    size_t subscriber_capacity;
    imgstream_subscriber_id next_subscriber_id;
    size_t producer_depth;
};


static int
imgstream_format_channels(enum imgstream_pixel_format format, size_t *channels)
{
    switch (format) {
	case IMGSTREAM_PIXEL_RGB8:
	    *channels = 3;
	    return 0;
	case IMGSTREAM_PIXEL_RGBA8:
	    *channels = 4;
	    return 0;
	default:
	    return -1;
    }
}


static int
imgstream_byte_count(size_t width, size_t height, size_t channels, size_t *bytes)
{
    if (!bytes || !width || !height || !channels)
	return -1;

    if (width > ((size_t)-1) / height)
	return -1;

    size_t pixels = width * height;
    if (pixels > ((size_t)-1) / channels)
	return -1;

    *bytes = pixels * channels;
    return 0;
}


static int
imgstream_rect_valid(const imgstream_t *stream, size_t x, size_t y, size_t width, size_t height)
{
    if (!stream || !width || !height)
	return 0;

    if (x >= stream->width || y >= stream->height)
	return 0;

    if (width > stream->width - x || height > stream->height - y)
	return 0;

    return 1;
}


static void
imgstream_notify(imgstream_t *stream, const struct imgstream_rect *rect)
{
    if (!stream || !rect)
	return;

    for (size_t i = 0; i < stream->subscriber_count; i++) {
	if (stream->subscribers[i].callback)
	    stream->subscribers[i].callback(stream->subscribers[i].ctx, rect, stream->generation);
    }
}


static void
imgstream_mark_dirty(imgstream_t *stream, const struct imgstream_rect *rect)
{
    if (!stream || !rect || !rect->width || !rect->height)
	return;

    if (!stream->dirty) {
	stream->dirty_rect = *rect;
	stream->dirty = 1;
    } else {
	size_t x0 = (rect->x < stream->dirty_rect.x) ? rect->x : stream->dirty_rect.x;
	size_t y0 = (rect->y < stream->dirty_rect.y) ? rect->y : stream->dirty_rect.y;
	size_t x1 = ((rect->x + rect->width) > (stream->dirty_rect.x + stream->dirty_rect.width)) ?
	    (rect->x + rect->width) : (stream->dirty_rect.x + stream->dirty_rect.width);
	size_t y1 = ((rect->y + rect->height) > (stream->dirty_rect.y + stream->dirty_rect.height)) ?
	    (rect->y + rect->height) : (stream->dirty_rect.y + stream->dirty_rect.height);

	stream->dirty_rect.x = x0;
	stream->dirty_rect.y = y0;
	stream->dirty_rect.width = x1 - x0;
	stream->dirty_rect.height = y1 - y0;
    }

    stream->generation++;
    imgstream_notify(stream, rect);
}


imgstream_t *
imgstream_create(size_t width, size_t height, enum imgstream_pixel_format format)
{
    size_t channels = 0;
    size_t bytes = 0;

    if (imgstream_format_channels(format, &channels) != 0)
	return NULL;

    if (imgstream_byte_count(width, height, channels, &bytes) != 0)
	return NULL;

    imgstream_t *stream = NULL;
    BU_ALLOC(stream, struct imgstream);
    stream->width = width;
    stream->height = height;
    stream->format = format;
    stream->channels = channels;
    stream->pixels = (unsigned char *)bu_calloc(bytes, sizeof(unsigned char), "imgstream pixels");
    stream->next_subscriber_id = 1;

    return stream;
}


imgstream_t *
imgstream_create_from_icv(const icv_image_t *image)
{
    if (!ICV_IMAGE_IS_INITIALIZED(image) || !image->data || !image->width || !image->height)
	return NULL;

    enum imgstream_pixel_format format = (image->channels >= 4) ? IMGSTREAM_PIXEL_RGBA8 : IMGSTREAM_PIXEL_RGB8;
    imgstream_t *stream = imgstream_create(image->width, image->height, format);
    if (!stream)
	return NULL;

    size_t src_channels = image->channels;
    for (size_t y = 0; y < stream->height; y++) {
	for (size_t x = 0; x < stream->width; x++) {
	    size_t dst_idx = (y * stream->width + x) * stream->channels;
	    size_t src_idx = (y * image->width + x) * src_channels;
	    for (size_t c = 0; c < stream->channels; c++) {
		double value = 1.0;
		if (c < src_channels)
		    value = image->data[src_idx + c];
		if (value < 0.0)
		    value = 0.0;
		if (value > 1.0)
		    value = 1.0;
		stream->pixels[dst_idx + c] = (unsigned char)(value * 255.0 + 0.5);
	    }
	}
    }

    struct imgstream_rect rect = {0, 0, stream->width, stream->height};
    imgstream_mark_dirty(stream, &rect);
    return stream;
}


void
imgstream_destroy(imgstream_t *stream)
{
    if (!stream)
	return;

    if (stream->pixels)
	bu_free(stream->pixels, "imgstream pixels");
    if (stream->subscribers)
	bu_free(stream->subscribers, "imgstream subscribers");
    BU_FREE(stream, struct imgstream);
}


size_t
imgstream_width(const imgstream_t *stream)
{
    return stream ? stream->width : 0;
}


size_t
imgstream_height(const imgstream_t *stream)
{
    return stream ? stream->height : 0;
}


enum imgstream_pixel_format
imgstream_format(const imgstream_t *stream)
{
    return stream ? stream->format : (enum imgstream_pixel_format)0;
}


size_t
imgstream_channels(const imgstream_t *stream)
{
    return stream ? stream->channels : 0;
}


uint64_t
imgstream_generation(const imgstream_t *stream)
{
    return stream ? stream->generation : 0;
}


int
imgstream_producer_active(const imgstream_t *stream)
{
    return (stream && stream->producer_depth) ? 1 : 0;
}


int
imgstream_get_info(const imgstream_t *stream, struct imgstream_info *info)
{
    if (!stream || !info)
	return -1;

    memset(info, 0, sizeof(*info));
    info->width = stream->width;
    info->height = stream->height;
    info->format = stream->format;
    info->channels = stream->channels;
    info->stride = stream->width * stream->channels;
    info->generation = stream->generation;
    info->dirty = stream->dirty ? 1 : 0;
    if (stream->dirty)
	info->dirty_rect = stream->dirty_rect;
    info->producer_active = stream->producer_depth ? 1 : 0;

    return 0;
}


int
imgstream_resize(imgstream_t *stream, size_t width, size_t height)
{
    size_t bytes = 0;

    if (!stream || stream->producer_depth)
	return -1;

    if (imgstream_byte_count(width, height, stream->channels, &bytes) != 0)
	return -1;

    unsigned char *pixels = (unsigned char *)bu_calloc(bytes, sizeof(unsigned char), "imgstream resized pixels");
    if (stream->pixels)
	bu_free(stream->pixels, "imgstream pixels");

    stream->pixels = pixels;
    stream->width = width;
    stream->height = height;
    stream->dirty = 0;

    struct imgstream_rect rect = {0, 0, width, height};
    imgstream_mark_dirty(stream, &rect);
    return 0;
}


int
imgstream_clear(imgstream_t *stream, const unsigned char *pixel)
{
    if (!stream || !stream->pixels)
	return -1;

    if (!pixel) {
	size_t bytes = 0;
	if (imgstream_byte_count(stream->width, stream->height, stream->channels, &bytes) != 0)
	    return -1;
	memset(stream->pixels, 0, bytes);
    } else {
	for (size_t y = 0; y < stream->height; y++) {
	    unsigned char *row = stream->pixels + y * stream->width * stream->channels;
	    for (size_t x = 0; x < stream->width; x++)
		memcpy(row + x * stream->channels, pixel, stream->channels);
	}
    }

    struct imgstream_rect rect = {0, 0, stream->width, stream->height};
    imgstream_mark_dirty(stream, &rect);
    return 0;
}


int
imgstream_write(imgstream_t *stream, const unsigned char *pixels, size_t stride)
{
    if (!stream)
	return -1;

    return imgstream_write_rect(stream, 0, 0, stream->width, stream->height, pixels, stride);
}


int
imgstream_read(const imgstream_t *stream, unsigned char *pixels, size_t stride)
{
    if (!stream)
	return -1;

    return imgstream_read_rect(stream, 0, 0, stream->width, stream->height, pixels, stride);
}


int
imgstream_write_rect(imgstream_t *stream, size_t x, size_t y, size_t width, size_t height, const unsigned char *pixels, size_t stride)
{
    if (!imgstream_rect_valid(stream, x, y, width, height) || !pixels)
	return -1;

    size_t row_bytes = width * stream->channels;
    if (stride < row_bytes)
	return -1;

    for (size_t row = 0; row < height; row++) {
	unsigned char *dst = stream->pixels + ((y + row) * stream->width + x) * stream->channels;
	const unsigned char *src = pixels + row * stride;
	memcpy(dst, src, row_bytes);
    }

    struct imgstream_rect rect = {x, y, width, height};
    imgstream_mark_dirty(stream, &rect);
    return 0;
}


int
imgstream_read_rect(const imgstream_t *stream, size_t x, size_t y, size_t width, size_t height, unsigned char *pixels, size_t stride)
{
    if (!imgstream_rect_valid(stream, x, y, width, height) || !pixels)
	return -1;

    size_t row_bytes = width * stream->channels;
    if (stride < row_bytes)
	return -1;

    for (size_t row = 0; row < height; row++) {
	const unsigned char *src = stream->pixels + ((y + row) * stream->width + x) * stream->channels;
	unsigned char *dst = pixels + row * stride;
	memcpy(dst, src, row_bytes);
    }

    return 0;
}


int
imgstream_dirty_rect(const imgstream_t *stream, struct imgstream_rect *rect)
{
    if (!stream || !rect || !stream->dirty)
	return 0;

    *rect = stream->dirty_rect;
    return 1;
}


void
imgstream_dirty_clear(imgstream_t *stream)
{
    if (!stream)
	return;

    stream->dirty = 0;
    memset(&stream->dirty_rect, 0, sizeof(stream->dirty_rect));
}


imgstream_subscriber_id
imgstream_subscribe(imgstream_t *stream, imgstream_dirty_cb callback, void *ctx)
{
    if (!stream || !callback)
	return -1;

    if (stream->subscriber_count == stream->subscriber_capacity) {
	size_t new_capacity = stream->subscriber_capacity ? stream->subscriber_capacity * 2 : 4;
	struct imgstream_subscriber *new_subscribers = (struct imgstream_subscriber *)bu_calloc(
	    new_capacity, sizeof(struct imgstream_subscriber), "imgstream subscribers");
	if (stream->subscribers) {
	    memcpy(new_subscribers, stream->subscribers, stream->subscriber_count * sizeof(struct imgstream_subscriber));
	    bu_free(stream->subscribers, "imgstream subscribers");
	}
	stream->subscribers = new_subscribers;
	stream->subscriber_capacity = new_capacity;
    }

    imgstream_subscriber_id id = stream->next_subscriber_id++;
    if (id < 0) {
	stream->next_subscriber_id = 1;
	id = stream->next_subscriber_id++;
    }

    stream->subscribers[stream->subscriber_count].id = id;
    stream->subscribers[stream->subscriber_count].callback = callback;
    stream->subscribers[stream->subscriber_count].ctx = ctx;
    stream->subscriber_count++;

    return id;
}


int
imgstream_unsubscribe(imgstream_t *stream, imgstream_subscriber_id id)
{
    if (!stream || id < 0)
	return -1;

    for (size_t i = 0; i < stream->subscriber_count; i++) {
	if (stream->subscribers[i].id == id) {
	    if (i + 1 < stream->subscriber_count) {
		memmove(&stream->subscribers[i], &stream->subscribers[i + 1],
		    (stream->subscriber_count - i - 1) * sizeof(struct imgstream_subscriber));
	    }
	    stream->subscriber_count--;
	    return 0;
	}
    }

    return -1;
}


int
imgstream_producer_begin(imgstream_t *stream)
{
    if (!stream)
	return -1;

    stream->producer_depth++;
    return 0;
}


int
imgstream_producer_end(imgstream_t *stream)
{
    if (!stream || !stream->producer_depth)
	return -1;

    stream->producer_depth--;
    return 0;
}


icv_image_t *
imgstream_snapshot_icv(const imgstream_t *stream)
{
    if (!stream || !stream->pixels)
	return NULL;

    icv_image_t *image = icv_image_create_with_channels(stream->width, stream->height,
	    ICV_COLOR_SPACE_RGB, stream->channels);
    if (!image)
	return NULL;

    for (size_t y = 0; y < stream->height; y++) {
	for (size_t x = 0; x < stream->width; x++) {
	    size_t idx = (y * stream->width + x) * stream->channels;
	    for (size_t c = 0; c < stream->channels; c++)
		image->data[idx + c] = ((double)stream->pixels[idx + c]) / 255.0;
	}
    }

    return image;
}


icv_image_t *
imgstream_to_icv(const imgstream_t *stream)
{
    return imgstream_snapshot_icv(stream);
}
