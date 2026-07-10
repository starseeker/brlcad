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
/** @file libimgstream/tests/stream.c */

#include "common.h"

#include <string.h>

#include "bu/log.h"
#include "bu/malloc.h"
#include "icv.h"
#include "imgstream.h"


struct callback_state {
    int count;
    uint64_t generation;
    struct imgstream_rect rect;
};


static void
dirty_cb(void *ctx, const struct imgstream_rect *rect, uint64_t generation)
{
    struct callback_state *state = (struct callback_state *)ctx;
    state->count++;
    state->generation = generation;
    state->rect = *rect;
}


#define CHECK(_expr, _msg) do { \
    if (!(_expr)) { \
	bu_log("FAIL: %s\n", _msg); \
	return 1; \
    } \
} while (0)


static int
test_lifecycle_resize_and_dirty(void)
{
    imgstream_t *stream = imgstream_create(4, 3, IMGSTREAM_PIXEL_RGB8);
    CHECK(stream != NULL, "created RGB stream");
    CHECK(imgstream_width(stream) == 4, "width recorded");
    CHECK(imgstream_height(stream) == 3, "height recorded");
    CHECK(imgstream_channels(stream) == 3, "RGB channel count recorded");
    CHECK(imgstream_generation(stream) == 0, "new stream is generation zero");

    struct imgstream_info info;
    CHECK(imgstream_get_info(stream, &info) == 0, "stream info query accepted");
    CHECK(info.width == 4 && info.height == 3, "stream info dimensions recorded");
    CHECK(info.format == IMGSTREAM_PIXEL_RGB8 && info.channels == 3, "stream info format recorded");
    CHECK(info.stride == 12, "stream info stride recorded");
    CHECK(info.generation == 0 && info.dirty == 0, "stream info initial state recorded");

    struct callback_state cb = {0, 0, {0, 0, 0, 0}};
    imgstream_subscriber_id id = imgstream_subscribe(stream, dirty_cb, &cb);
    CHECK(id > 0, "subscriber registered");

    unsigned char pixel[3] = {10, 20, 30};
    CHECK(imgstream_clear(stream, pixel) == 0, "clear accepted");
    CHECK(imgstream_generation(stream) == 1, "clear advanced generation");
    CHECK(cb.count == 1, "subscriber received clear notification");
    CHECK(cb.rect.x == 0 && cb.rect.y == 0 && cb.rect.width == 4 && cb.rect.height == 3, "clear notification covers full image");
    CHECK(imgstream_get_info(stream, &info) == 0, "stream info query after clear accepted");
    CHECK(info.generation == 1 && info.dirty == 1, "stream info reports dirty generation");
    CHECK(info.dirty_rect.width == 4 && info.dirty_rect.height == 3, "stream info reports dirty rect");

    struct imgstream_rect dirty;
    CHECK(imgstream_dirty_rect(stream, &dirty) == 1, "dirty rect present");
    CHECK(dirty.width == 4 && dirty.height == 3, "dirty rect covers full image");
    imgstream_dirty_clear(stream);
    CHECK(imgstream_dirty_rect(stream, &dirty) == 0, "dirty clear removed dirty rect");

    CHECK(imgstream_producer_begin(stream) == 0, "producer begin accepted");
    CHECK(imgstream_producer_active(stream) == 1, "producer marked active");
    CHECK(imgstream_get_info(stream, &info) == 0 && info.producer_active == 1, "stream info reports active producer");
    CHECK(imgstream_resize(stream, 2, 2) == -1, "resize rejected during active producer");
    CHECK(imgstream_producer_end(stream) == 0, "producer end accepted");
    CHECK(imgstream_resize(stream, 2, 2) == 0, "resize accepted after producer end");
    CHECK(imgstream_width(stream) == 2 && imgstream_height(stream) == 2, "resized dimensions recorded");

    CHECK(imgstream_unsubscribe(stream, id) == 0, "subscriber removed");
    cb.count = 0;
    CHECK(imgstream_clear(stream, NULL) == 0, "zero clear accepted");
    CHECK(cb.count == 0, "removed subscriber not called");

    imgstream_destroy(stream);
    return 0;
}


static int
test_read_write_rect(void)
{
    imgstream_t *stream = imgstream_create(4, 4, IMGSTREAM_PIXEL_RGBA8);
    CHECK(stream != NULL, "created RGBA stream");
    CHECK(imgstream_channels(stream) == 4, "RGBA channel count recorded");

    unsigned char full[4 * 4 * 4];
    for (size_t i = 0; i < sizeof(full); i++)
	full[i] = (unsigned char)i;

    CHECK(imgstream_write(stream, full, 4 * 4) == 0, "full write accepted");

    unsigned char rect[2 * 2 * 4];
    memset(rect, 200, sizeof(rect));
    CHECK(imgstream_write_rect(stream, 1, 1, 2, 2, rect, 2 * 4) == 0, "rect write accepted");

    unsigned char out[2 * 2 * 4];
    memset(out, 0, sizeof(out));
    CHECK(imgstream_read_rect(stream, 1, 1, 2, 2, out, 2 * 4) == 0, "rect read accepted");
    CHECK(memcmp(rect, out, sizeof(rect)) == 0, "rect read returns written pixels");

    struct imgstream_rect dirty;
    CHECK(imgstream_dirty_rect(stream, &dirty) == 1, "dirty rect present after writes");
    CHECK(dirty.x == 0 && dirty.y == 0 && dirty.width == 4 && dirty.height == 4, "dirty rect is union of full and rect writes");

    CHECK(imgstream_write_rect(stream, 3, 3, 2, 2, rect, 2 * 4) == -1, "out-of-bounds write rejected");
    CHECK(imgstream_read_rect(stream, 0, 0, 2, 2, out, 1) == -1, "short stride rejected");

    imgstream_destroy(stream);
    return 0;
}


static int
test_icv_conversion(void)
{
    imgstream_t *stream = imgstream_create(2, 1, IMGSTREAM_PIXEL_RGB8);
    CHECK(stream != NULL, "created RGB stream for ICV conversion");

    unsigned char pixels[6] = {0, 127, 255, 255, 0, 64};
    CHECK(imgstream_write(stream, pixels, 6) == 0, "wrote RGB pixels");

    icv_image_t *image = imgstream_snapshot_icv(stream);
    CHECK(image != NULL, "snapshotted stream to ICV image");
    CHECK(image->width == 2 && image->height == 1 && image->channels == 3, "ICV image dimensions/channels match");
    CHECK(image->data[0] > -0.000001 && image->data[0] < 0.000001, "first ICV red channel near zero");
    CHECK(image->data[2] > 0.999999 && image->data[2] < 1.000001, "first ICV blue channel near one");

    imgstream_t *roundtrip = imgstream_create_from_icv(image);
    CHECK(roundtrip != NULL, "created stream from ICV image");
    unsigned char out[6] = {0, 0, 0, 0, 0, 0};
    CHECK(imgstream_read(roundtrip, out, 6) == 0, "read roundtrip stream");
    CHECK(memcmp(pixels, out, sizeof(pixels)) == 0, "roundtrip preserved RGB bytes");

    icv_destroy(image);
    imgstream_destroy(roundtrip);
    imgstream_destroy(stream);

    icv_image_t *rgba = icv_image_create_with_channels(1, 1, ICV_COLOR_SPACE_RGB, 4);
    CHECK(rgba != NULL, "created RGBA ICV image");
    rgba->data[0] = 1.0;
    rgba->data[1] = 0.5;
    rgba->data[2] = 0.0;
    rgba->data[3] = 0.25;

    imgstream_t *rgba_stream = imgstream_create_from_icv(rgba);
    CHECK(rgba_stream != NULL, "created RGBA stream from ICV image");
    CHECK(imgstream_format(rgba_stream) == IMGSTREAM_PIXEL_RGBA8, "RGBA stream format selected from four-channel ICV image");

    icv_image_t *rgba_snapshot = imgstream_snapshot_icv(rgba_stream);
    CHECK(rgba_snapshot != NULL, "converted RGBA stream to ICV image");
    CHECK(rgba_snapshot->channels == 4 && rgba_snapshot->alpha_channel == 1, "RGBA ICV snapshot preserves alpha channel");

    icv_destroy(rgba_snapshot);
    imgstream_destroy(rgba_stream);
    icv_destroy(rgba);
    return 0;
}


int
main(int ac, char **av)
{
    (void)ac;
    (void)av;

    if (test_lifecycle_resize_and_dirty())
	return 1;
    if (test_read_write_rect())
	return 1;
    if (test_icv_conversion())
	return 1;

    return 0;
}
