/*                    F B _ C O M P A T . C
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
/** @file libimgstream/fb_compat.c */

#include "common.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bu/log.h"
#include "bu/malloc.h"
#include "bu/str.h"
#include "imgstream.h"


#define IMGSTREAM_FB_DEFAULT_SIZE 512
#define IMGSTREAM_FB_RGB_CHANNELS 3
#define IMGSTREAM_FB_FANOUT_MAX 32
#define IMGSTREAM_FB_DIAGNOSTIC_LOG_CMAP 2
#define IMGSTREAM_FB_DIAGNOSTIC_LOG_PIXELS 4

struct imgstream_fb {
    char *name;
    char *file_path;
    imgstream_t *stream;
    int owns_stream;
    int file_backed;
    int display_backed;
    enum imgstream_fb_diagnostic_kind diagnostic_kind;
    unsigned int diagnostic_flags;
    struct imgstream_fb_display_host display_host;
    void *display_host_data;
    struct imgstream_fb_colormap colormap;
    struct imgstream_fb **children;
    size_t child_count;
    struct imgstream_fb_view view;
    struct imgstream_fb_cursor cursor;
};


static struct imgstream_fb_display_host fb_display_host;
static void *fb_display_host_data = NULL;
static int fb_display_host_active = 0;


static int
fb_display_host_registered(void)
{
    return fb_display_host_active && fb_display_host.open;
}


static int
fb_spec_numeric(const char *spec, size_t len)
{
    if (!spec || !len)
	return 0;

    for (size_t i = 0; i < len; i++) {
	if (spec[i] < '0' || spec[i] > '9')
	    return 0;
    }

    return 1;
}


static int
fb_spec_int(const char *spec, size_t len, int *value)
{
    if (!value || !fb_spec_numeric(spec, len))
	return 0;

    long parsed = strtol(spec, NULL, 10);
    if (parsed < 0 || parsed > INT_MAX)
	return 0;

    *value = (int)parsed;
    return 1;
}


static int
fb_remote_port_map(int port)
{
    return port < 1024 ? port + 5559 : port;
}


static const char *
fb_localhost(size_t *len)
{
    if (len)
	*len = 9;
    return "localhost";
}


static int
fb_spec_empty_or_memory(const char *spec)
{
    if (!spec || !spec[0])
	return 1;

    if (BU_STR_EQUAL(spec, "mem") || BU_STR_EQUAL(spec, "memory"))
	return 1;

    if (BU_STR_EQUAL(spec, "/dev/mem"))
	return 1;

    if (bu_strncmp(spec, "/dev/mem", 8) == 0)
	return 1;

    return 0;
}


static int
fb_spec_is_numeric_remote(const char *spec)
{
    return spec && fb_spec_numeric(spec, strlen(spec));
}


static int
fb_spec_has_remote_shape(const char *spec)
{
    if (!spec || !spec[0])
	return 0;

    if (bu_strncmp(spec, "ipc:", 4) == 0)
	return 1;

    if (bu_strncmp(spec, "tcp:", 4) == 0)
	return 1;

    if (strchr(spec, ':'))
	return 1;

    return 0;
}


static enum imgstream_fb_display_kind
fb_display_kind_from_spec(const char *spec)
{
    if (!spec || !spec[0])
	return IMGSTREAM_FB_DISPLAY_NONE;

    if (BU_STR_EQUAL(spec, "/dev/qtgl"))
	return IMGSTREAM_FB_DISPLAY_QTGL;
    if (BU_STR_EQUAL(spec, "/dev/ogl"))
	return IMGSTREAM_FB_DISPLAY_OGL;
    if (BU_STR_EQUAL(spec, "/dev/wgl"))
	return IMGSTREAM_FB_DISPLAY_WGL;
    if (BU_STR_EQUAL(spec, "/dev/X"))
	return IMGSTREAM_FB_DISPLAY_X;
    if (BU_STR_EQUAL(spec, "/dev/swrast"))
	return IMGSTREAM_FB_DISPLAY_SWRAST;

    return IMGSTREAM_FB_DISPLAY_NONE;
}


static const char *
fb_file_path_from_spec(const char *spec)
{
    if (!spec || !spec[0])
	return NULL;

    if (BU_STR_EQUAL(spec, "/dev/disk"))
	return NULL;

    if (bu_strncmp(spec, "/dev/disk:", 10) == 0 ||
	    bu_strncmp(spec, "/dev/disk=", 10) == 0)
	return spec + 10;

    if (bu_strncmp(spec, "/dev/disk/", 10) == 0)
	return spec + 9;

    if (bu_strncmp(spec, "/dev/", 5) != 0 && !fb_spec_has_remote_shape(spec))
	return spec;

    return NULL;
}


static int
fb_remote_info_from_legacy_spec(const char *spec, struct imgstream_fb_spec_info *info)
{
    const char *file = spec;
    const char *colon = NULL;
    const char *rest = NULL;
    int port = 0;

    if (!file || !info)
	return -1;

    if (fb_spec_int(file, strlen(file), &port)) {
	info->remote = IMGSTREAM_FB_REMOTE_TCP;
	info->host = fb_localhost(&info->host_len);
	info->port = fb_remote_port_map(port);
	return 0;
    }

    colon = strchr(file, ':');
    if (!colon)
	return -1;

    rest = colon + 1;
    size_t prefix_len = (size_t)(colon - file);
    if (fb_spec_int(file, prefix_len, &port)) {
	info->remote = IMGSTREAM_FB_REMOTE_TCP;
	info->host = fb_localhost(&info->host_len);
	info->port = fb_remote_port_map(port);
	info->device = rest;
	info->device_len = strlen(rest);
	return 0;
    }

    info->remote = IMGSTREAM_FB_REMOTE_TCP;
    info->host = file;
    info->host_len = prefix_len;
    if (info->host_len == 0 ||
	    (info->host_len == 4 && bu_strncmp(info->host, "unix", 4) == 0))
	info->host = fb_localhost(&info->host_len);

    if (fb_spec_int(rest, strlen(rest), &port)) {
	info->port = fb_remote_port_map(port);
	return 0;
    }

    const char *port_colon = strchr(rest, ':');
    if (port_colon && fb_spec_int(rest, (size_t)(port_colon - rest), &port)) {
	info->port = fb_remote_port_map(port);
	info->device = port_colon + 1;
	info->device_len = strlen(info->device);
	return 0;
    }

    info->port = 5558;
    info->device = rest;
    info->device_len = strlen(rest);
    return 0;
}


static int
fb_remote_info_from_spec(const char *spec, struct imgstream_fb_spec_info *info)
{
    if (!spec || !info)
	return -1;

    info->target = spec;
    info->target_len = strlen(spec);

    if (bu_strncmp(spec, "ipc:", 4) == 0 && spec[4] != '\0') {
	info->remote = IMGSTREAM_FB_REMOTE_IPC;
	info->target = spec + 4;
	info->target_len = strlen(spec + 4);
	return 0;
    }

    if (bu_strncmp(spec, "tcp:", 4) == 0 && spec[4] != '\0') {
	info->target = spec + 4;
	info->target_len = strlen(spec + 4);
	return fb_remote_info_from_legacy_spec(spec + 4, info);
    }

    if (BU_STR_EQUAL(spec, "/dev/remote")) {
	info->remote = IMGSTREAM_FB_REMOTE_TCP;
	info->host = fb_localhost(&info->host_len);
	info->port = 5558;
	return 0;
    }

    if (bu_strncmp(spec, "/dev/remote:", 12) == 0 ||
	    bu_strncmp(spec, "/dev/remote=", 12) == 0) {
	info->target = spec + 12;
	info->target_len = strlen(spec + 12);
	return fb_remote_info_from_legacy_spec(spec + 12, info);
    }

    return fb_remote_info_from_legacy_spec(spec, info);
}


static enum imgstream_fb_diagnostic_kind
fb_diagnostic_kind_from_spec(const char *spec)
{
    if (!spec || !spec[0])
	return IMGSTREAM_FB_DIAGNOSTIC_NONE;

    if (BU_STR_EQUAL(spec, "/dev/null"))
	return IMGSTREAM_FB_DIAGNOSTIC_NULL;

    if (BU_STR_EQUAL(spec, "/dev/txt"))
	return IMGSTREAM_FB_DIAGNOSTIC_TXT;

    if (BU_STR_EQUAL(spec, "/dev/debug") ||
	    (bu_strncmp(spec, "/dev/debug", 10) == 0 && spec[10] >= '0' && spec[10] <= '9'))
	return IMGSTREAM_FB_DIAGNOSTIC_DEBUG;

    return IMGSTREAM_FB_DIAGNOSTIC_NONE;
}


static unsigned int
fb_diagnostic_flags_from_spec(const char *spec)
{
    if (!spec)
	return 0;

    const char *cp = spec;
    while (*cp && (*cp < '0' || *cp > '9'))
	cp++;
    return *cp ? (unsigned int)strtoul(cp, NULL, 10) : 0;
}


static int
fb_is_diagnostic(const imgstream_fb_t *fb)
{
    return fb && fb->diagnostic_kind != IMGSTREAM_FB_DIAGNOSTIC_NONE;
}


void
imgstream_fb_colormap_linear(struct imgstream_fb_colormap *cmap)
{
    if (!cmap)
	return;

    for (size_t i = 0; i < 256; i++) {
	uint16_t value = (uint16_t)(i << 8);
	cmap->red[i] = value;
	cmap->green[i] = value;
	cmap->blue[i] = value;
    }
}


int
imgstream_fb_colormap_is_linear(const struct imgstream_fb_colormap *cmap)
{
    if (!cmap)
	return 0;

    for (size_t i = 0; i < 256; i++) {
	uint16_t value = (uint16_t)(i << 8);
	if (cmap->red[i] != value || cmap->green[i] != value ||
		cmap->blue[i] != value)
	    return 0;
    }

    return 1;
}


static char *
fb_trimmed_copy(const char *start, size_t len)
{
    if (!start)
	return NULL;

    while (len && (*start == ' ' || *start == '\t')) {
	start++;
	len--;
    }
    while (len && (start[len - 1] == ' ' || start[len - 1] == '\t'))
	len--;

    if (!len)
	return NULL;

    char *copy = (char *)bu_malloc(len + 1, "imgstream fb spec token");
    memcpy(copy, start, len);
    copy[len] = '\0';
    return copy;
}


static void
fb_child_specs_free(char **children, size_t child_count)
{
    if (!children)
	return;

    for (size_t i = 0; i < child_count; i++) {
	if (children[i])
	    bu_free(children[i], "imgstream fb child spec");
    }
    bu_free(children, "imgstream fb child spec list");
}


static int
fb_fanout_parse(const char *spec, char ***children_out, size_t *child_count_out)
{
    if (!children_out || !child_count_out)
	return -1;

    *children_out = NULL;
    *child_count_out = 0;

    if (!spec || bu_strncmp(spec, "/dev/stack", 10) != 0)
	return -1;

    const char *cp = spec + 10;
    while (*cp && *cp != ' ' && *cp != '\t')
	cp++;
    if (!*cp)
	return -1;

    char **children = (char **)bu_calloc(IMGSTREAM_FB_FANOUT_MAX, sizeof(char *),
	    "imgstream fb child spec list");
    size_t child_count = 0;

    while (*cp && child_count < IMGSTREAM_FB_FANOUT_MAX) {
	while (*cp == ' ' || *cp == '\t' || *cp == ';')
	    cp++;
	if (!*cp)
	    break;

	const char *start = cp;
	while (*cp && *cp != ';')
	    cp++;

	char *child = fb_trimmed_copy(start, (size_t)(cp - start));
	if (child)
	    children[child_count++] = child;
    }

    if (!child_count) {
	fb_child_specs_free(children, child_count);
	return -1;
    }

    *children_out = children;
    *child_count_out = child_count;
    return 0;
}


enum imgstream_fb_spec_kind
imgstream_fb_spec_kind(const char *spec)
{
    if (fb_spec_empty_or_memory(spec))
	return IMGSTREAM_FB_SPEC_MEMORY;

    if (!spec || !spec[0])
	return IMGSTREAM_FB_SPEC_MEMORY;

    if (fb_spec_is_numeric_remote(spec))
	return IMGSTREAM_FB_SPEC_REMOTE;

    if (BU_STR_EQUAL(spec, "/dev/disk") || bu_strncmp(spec, "/dev/disk", 9) == 0)
	return IMGSTREAM_FB_SPEC_FILE;

    if (BU_STR_EQUAL(spec, "/dev/stack") || bu_strncmp(spec, "/dev/stack", 10) == 0)
	return IMGSTREAM_FB_SPEC_FANOUT;

    if (BU_STR_EQUAL(spec, "/dev/remote") || bu_strncmp(spec, "/dev/remote", 11) == 0 || fb_spec_has_remote_shape(spec))
	return IMGSTREAM_FB_SPEC_REMOTE;

    if (fb_diagnostic_kind_from_spec(spec) != IMGSTREAM_FB_DIAGNOSTIC_NONE)
	return IMGSTREAM_FB_SPEC_DIAGNOSTIC;

    if (fb_display_kind_from_spec(spec) != IMGSTREAM_FB_DISPLAY_NONE)
	return IMGSTREAM_FB_SPEC_DISPLAY;

    if (bu_strncmp(spec, "/dev/", 5) != 0 && !fb_spec_has_remote_shape(spec))
	return IMGSTREAM_FB_SPEC_FILE;

    return IMGSTREAM_FB_SPEC_UNSUPPORTED;
}


int
imgstream_fb_spec_info(const char *spec, struct imgstream_fb_spec_info *info)
{
    if (!info)
	return -1;

    memset(info, 0, sizeof(*info));
    info->port = -1;
    info->kind = imgstream_fb_spec_kind(spec);

    if (!spec || !spec[0])
	return 0;

    switch (info->kind) {
	case IMGSTREAM_FB_SPEC_FILE:
	    info->target = fb_file_path_from_spec(spec);
	    if (info->target)
		info->target_len = strlen(info->target);
	    break;
	case IMGSTREAM_FB_SPEC_FANOUT:
	    if (bu_strncmp(spec, "/dev/stack", 10) == 0) {
		const char *cp = spec + 10;
		while (*cp == ' ' || *cp == '\t')
		    cp++;
		info->target = cp;
		info->target_len = strlen(cp);
	    }
	    break;
	case IMGSTREAM_FB_SPEC_REMOTE:
	    return fb_remote_info_from_spec(spec, info);
	case IMGSTREAM_FB_SPEC_DISPLAY:
	    info->display = fb_display_kind_from_spec(spec);
	    if (bu_strncmp(spec, "/dev/", 5) == 0) {
		info->target = spec + 5;
		info->target_len = strlen(spec + 5);
	    }
	    break;
	case IMGSTREAM_FB_SPEC_DIAGNOSTIC:
	    info->diagnostic = fb_diagnostic_kind_from_spec(spec);
	    if (bu_strncmp(spec, "/dev/", 5) == 0) {
		info->target = spec + 5;
		info->target_len = strlen(spec + 5);
	    }
	    break;
	case IMGSTREAM_FB_SPEC_MEMORY:
	case IMGSTREAM_FB_SPEC_UNSUPPORTED:
	    break;
    }

    return 0;
}


int
imgstream_fb_display_host_set(const struct imgstream_fb_display_host *host, void *data)
{
    if (!host || !host->open)
	return -1;

    fb_display_host = *host;
    fb_display_host_data = data;
    fb_display_host_active = 1;
    return 0;
}


void
imgstream_fb_display_host_clear(void)
{
    memset(&fb_display_host, 0, sizeof(fb_display_host));
    fb_display_host_data = NULL;
    fb_display_host_active = 0;
}


const char *
imgstream_fb_spec_kind_name(enum imgstream_fb_spec_kind kind)
{
    switch (kind) {
	case IMGSTREAM_FB_SPEC_MEMORY:
	    return "memory";
	case IMGSTREAM_FB_SPEC_FILE:
	    return "file";
	case IMGSTREAM_FB_SPEC_FANOUT:
	    return "fanout";
	case IMGSTREAM_FB_SPEC_REMOTE:
	    return "remote";
	case IMGSTREAM_FB_SPEC_DISPLAY:
	    return "display";
	case IMGSTREAM_FB_SPEC_DIAGNOSTIC:
	    return "diagnostic";
	case IMGSTREAM_FB_SPEC_UNSUPPORTED:
	    return "unsupported";
    }

    return "unsupported";
}


static int
fb_spec_supported_no_fanout(const char *spec)
{
    enum imgstream_fb_spec_kind kind = imgstream_fb_spec_kind(spec);
    if (kind == IMGSTREAM_FB_SPEC_MEMORY)
	return 1;

    if (kind == IMGSTREAM_FB_SPEC_FILE) {
	const char *path = fb_file_path_from_spec(spec);
	return (path && path[0] && !BU_STR_EQUAL(path, "-")) ? 1 : 0;
    }

    if (kind == IMGSTREAM_FB_SPEC_DIAGNOSTIC)
	return 1;

    if (kind == IMGSTREAM_FB_SPEC_DISPLAY)
	return fb_display_host_registered() ? 1 : 0;

    return 0;
}


int
imgstream_fb_spec_supported(const char *spec)
{
    enum imgstream_fb_spec_kind kind = imgstream_fb_spec_kind(spec);
    if (kind == IMGSTREAM_FB_SPEC_FANOUT) {
	char **children = NULL;
	size_t child_count = 0;
	int supported = 1;

	if (fb_fanout_parse(spec, &children, &child_count) != 0)
	    return 0;

	for (size_t i = 0; i < child_count; i++) {
	    if (!fb_spec_supported_no_fanout(children[i])) {
		supported = 0;
		break;
	    }
	}
	fb_child_specs_free(children, child_count);
	return supported;
    }

    return fb_spec_supported_no_fanout(spec);
}


static int
fb_rgb_byte_count(size_t width, size_t height, size_t *bytes)
{
    if (!bytes || !width || !height)
	return -1;
    if (width > ((size_t)-1) / height)
	return -1;
    size_t pixels = width * height;
    if (pixels > ((size_t)-1) / IMGSTREAM_FB_RGB_CHANNELS)
	return -1;

    *bytes = pixels * IMGSTREAM_FB_RGB_CHANNELS;
    return 0;
}


static int
fb_rgb_pixel_byte_count(size_t pixels, size_t *bytes)
{
    if (!bytes)
	return -1;
    if (pixels > ((size_t)-1) / IMGSTREAM_FB_RGB_CHANNELS)
	return -1;

    *bytes = pixels * IMGSTREAM_FB_RGB_CHANNELS;
    return 0;
}


static void
fb_diagnostic_log_pixels(const imgstream_fb_t *fb, const unsigned char *rgb, size_t count)
{
    if (!fb || fb->diagnostic_kind != IMGSTREAM_FB_DIAGNOSTIC_DEBUG ||
	    !(fb->diagnostic_flags & IMGSTREAM_FB_DIAGNOSTIC_LOG_PIXELS) || !rgb)
	return;

    for (size_t i = 0; i < count; i++) {
	if (i % 4 == 0)
	    bu_log("%4zu:", i);
	bu_log("  [%3d, %3d, %3d]",
		rgb[(i * IMGSTREAM_FB_RGB_CHANNELS) + 0],
		rgb[(i * IMGSTREAM_FB_RGB_CHANNELS) + 1],
		rgb[(i * IMGSTREAM_FB_RGB_CHANNELS) + 2]);
	if (i % 4 == 3)
	    bu_log("\n");
    }
    if (count % 4 != 0)
	bu_log("\n");
}


static void
fb_diagnostic_log(const imgstream_fb_t *fb, const char *op)
{
    if (fb && fb->diagnostic_kind == IMGSTREAM_FB_DIAGNOSTIC_DEBUG)
	bu_log("imgstream_fb_%s(%s)\n", op, fb->name ? fb->name : "/dev/debug");
}


static int
fb_stream_load_file(imgstream_t *stream, const char *path)
{
    if (!stream || !path || !path[0])
	return -1;

    FILE *fp = fopen(path, "rb");
    if (!fp)
	return 0;

    size_t width = imgstream_width(stream);
    size_t row_bytes = width * IMGSTREAM_FB_RGB_CHANNELS;
    unsigned char *row = (unsigned char *)bu_calloc(row_bytes, sizeof(unsigned char), "imgstream fb file row");
    int ret = 0;

    for (size_t y = 0; y < imgstream_height(stream); y++) {
	size_t got = fread(row, 1, row_bytes, fp);
	if (got < row_bytes)
	    memset(row + got, 0, row_bytes - got);
	if (got || y == 0) {
	    if (imgstream_write_rect(stream, 0, y, width, 1, row, row_bytes) != 0) {
		ret = -1;
		break;
	    }
	}
	if (got < row_bytes)
	    break;
    }

    if (ferror(fp))
	ret = -1;

    fclose(fp);
    bu_free(row, "imgstream fb file row");
    imgstream_dirty_clear(stream);
    return ret;
}


static int
fb_stream_flush_file(const imgstream_fb_t *fb)
{
    if (!fb || !fb->stream || !fb->file_backed || !fb->file_path)
	return 0;

    size_t byte_count = 0;
    if (fb_rgb_byte_count(imgstream_fb_width(fb), imgstream_fb_height(fb), &byte_count) != 0)
	return -1;
    if (!byte_count)
	return -1;

    FILE *fp = fopen(fb->file_path, "wb");
    if (!fp)
	return -1;

    size_t width = imgstream_fb_width(fb);
    size_t row_bytes = width * IMGSTREAM_FB_RGB_CHANNELS;
    unsigned char *row = (unsigned char *)bu_malloc(row_bytes, "imgstream fb file row");
    int ret = 0;

    for (size_t y = 0; y < imgstream_fb_height(fb); y++) {
	if (imgstream_read_rect(fb->stream, 0, y, width, 1, row, row_bytes) != 0) {
	    ret = -1;
	    break;
	}
	if (fwrite(row, 1, row_bytes, fp) != row_bytes) {
	    ret = -1;
	    break;
	}
    }

    bu_free(row, "imgstream fb file row");
    if (fclose(fp) != 0)
	ret = -1;

    return ret;
}


static imgstream_fb_t *
fb_wrap_stream(const char *name, imgstream_t *stream, const char *file_path)
{
    if (!stream)
	return NULL;

    size_t width = imgstream_width(stream);
    size_t height = imgstream_height(stream);
    imgstream_fb_t *fb = NULL;
    BU_ALLOC(fb, struct imgstream_fb);
    fb->name = bu_strdup((name && name[0]) ? name : "/dev/mem");
    if (file_path && file_path[0]) {
	fb->file_path = bu_strdup(file_path);
	fb->file_backed = 1;
    }
    fb->stream = stream;
    fb->owns_stream = 1;
    fb->view.xcenter = (int)(width / 2);
    fb->view.ycenter = (int)(height / 2);
    fb->view.xzoom = 1;
    fb->view.yzoom = 1;
    fb->cursor.mode = 0;
    fb->cursor.x = 0;
    fb->cursor.y = 0;
    imgstream_fb_colormap_linear(&fb->colormap);

    return fb;
}


static imgstream_fb_t *
fb_open_fanout(const char *spec, size_t width, size_t height)
{
    char **child_specs = NULL;
    size_t child_count = 0;

    if (fb_fanout_parse(spec, &child_specs, &child_count) != 0)
	return NULL;

    struct imgstream_fb **children = (struct imgstream_fb **)bu_calloc(child_count,
	    sizeof(struct imgstream_fb *), "imgstream fb fanout children");
    size_t opened = 0;
    for (size_t i = 0; i < child_count; i++) {
	children[i] = imgstream_fb_open(child_specs[i], width, height);
	if (!children[i])
	    break;
	opened++;
    }

    fb_child_specs_free(child_specs, child_count);

    if (opened != child_count) {
	for (size_t i = 0; i < opened; i++)
	    imgstream_fb_close(children[i]);
	bu_free(children, "imgstream fb fanout children");
	return NULL;
    }

    imgstream_fb_t *fb = NULL;
    BU_ALLOC(fb, struct imgstream_fb);
    fb->name = bu_strdup((spec && spec[0]) ? spec : "/dev/stack");
    fb->children = children;
    fb->child_count = child_count;
    fb->stream = children[0]->stream;
    fb->owns_stream = 0;
    fb->view.xcenter = (int)(imgstream_fb_width(children[0]) / 2);
    fb->view.ycenter = (int)(imgstream_fb_height(children[0]) / 2);
    fb->view.xzoom = 1;
    fb->view.yzoom = 1;
    fb->cursor.mode = 0;
    fb->cursor.x = 0;
    fb->cursor.y = 0;
    imgstream_fb_colormap_linear(&fb->colormap);

    return fb;
}


imgstream_fb_t *
imgstream_fb_open(const char *spec, size_t width, size_t height)
{
    if (!imgstream_fb_spec_supported(spec))
	return NULL;

    enum imgstream_fb_spec_kind kind = imgstream_fb_spec_kind(spec);

    if (kind == IMGSTREAM_FB_SPEC_FANOUT)
	return fb_open_fanout(spec, width, height);

    if (!width)
	width = IMGSTREAM_FB_DEFAULT_SIZE;
    if (!height)
	height = IMGSTREAM_FB_DEFAULT_SIZE;

    const char *file_path = NULL;
    if (kind == IMGSTREAM_FB_SPEC_FILE)
	file_path = fb_file_path_from_spec(spec);
    enum imgstream_fb_diagnostic_kind diagnostic_kind = fb_diagnostic_kind_from_spec(spec);

    imgstream_t *stream = imgstream_create(width, height, IMGSTREAM_PIXEL_RGB8);
    if (!stream)
	return NULL;

    if (file_path) {
	FILE *fp = fopen(file_path, "r+b");
	if (!fp)
	    fp = fopen(file_path, "w+b");
	if (!fp) {
	    imgstream_destroy(stream);
	    return NULL;
	}
	fclose(fp);
	if (fb_stream_load_file(stream, file_path) != 0) {
	    imgstream_destroy(stream);
	    return NULL;
	}
    }

    imgstream_fb_t *fb = fb_wrap_stream(spec, stream, file_path);
    if (fb) {
	fb->diagnostic_kind = diagnostic_kind;
	fb->diagnostic_flags = fb_diagnostic_flags_from_spec(spec);
	if (kind == IMGSTREAM_FB_SPEC_DISPLAY) {
	    struct imgstream_fb_spec_info info;
	    struct imgstream_fb_display_host host = fb_display_host;
	    void *host_data = fb_display_host_data;
	    if (!fb_display_host_registered() ||
		    imgstream_fb_spec_info(fb->name, &info) != 0 ||
		    host.open(fb, &info, host_data) != 0) {
		imgstream_fb_close(fb);
		return NULL;
	    }
	    fb->display_host = host;
	    fb->display_host_data = host_data;
	    fb->display_backed = 1;
	}
	fb_diagnostic_log(fb, "open");
    }
    return fb;
}


void
imgstream_fb_close(imgstream_fb_t *fb)
{
    if (!fb)
	return;

    if (fb->child_count) {
	for (size_t i = 0; i < fb->child_count; i++)
	    imgstream_fb_close(fb->children[i]);
	bu_free(fb->children, "imgstream fb fanout children");
    } else if (fb->file_backed) {
	(void)fb_stream_flush_file(fb);
    } else if (fb->display_backed && fb->display_host.close) {
	fb->display_host.close(fb, fb->display_host_data);
    }
    fb_diagnostic_log(fb, "close");
    if (fb->owns_stream && fb->stream)
	imgstream_destroy(fb->stream);
    if (fb->file_path)
	bu_free(fb->file_path, "imgstream fb file path");
    if (fb->name)
	bu_free(fb->name, "imgstream fb name");
    BU_FREE(fb, struct imgstream_fb);
}


imgstream_t *
imgstream_fb_stream(imgstream_fb_t *fb)
{
    return fb ? fb->stream : NULL;
}


const imgstream_t *
imgstream_fb_cstream(const imgstream_fb_t *fb)
{
    return fb ? fb->stream : NULL;
}


const char *
imgstream_fb_name(const imgstream_fb_t *fb)
{
    return fb ? fb->name : NULL;
}


size_t
imgstream_fb_width(const imgstream_fb_t *fb)
{
    return fb ? imgstream_width(fb->stream) : 0;
}


size_t
imgstream_fb_height(const imgstream_fb_t *fb)
{
    return fb ? imgstream_height(fb->stream) : 0;
}


int
imgstream_fb_clear(imgstream_fb_t *fb, const unsigned char *rgb)
{
    if (!fb || !fb->stream)
	return -1;

    if (fb->child_count) {
	for (size_t i = 0; i < fb->child_count; i++) {
	    if (imgstream_fb_clear(fb->children[i], rgb) != 0)
		return -1;
	}
	return 0;
    }

    if (fb_is_diagnostic(fb)) {
	fb_diagnostic_log(fb, "clear");
	return 0;
    }

    return imgstream_clear(fb->stream, rgb);
}


static int
fb_coords_valid(const imgstream_fb_t *fb, int x, int y)
{
    if (!fb || !fb->stream || x < 0 || y < 0)
	return 0;

    if ((size_t)x >= imgstream_fb_width(fb) || (size_t)y >= imgstream_fb_height(fb))
	return 0;

    return 1;
}


ssize_t
imgstream_fb_read(const imgstream_fb_t *fb, int x, int y, unsigned char *rgb, size_t count)
{
    if (fb && fb->child_count)
	return imgstream_fb_read(fb->children[0], x, y, rgb, count);

    if (fb_is_diagnostic(fb)) {
	fb_diagnostic_log(fb, "read");
	if (rgb) {
	    size_t bytes = 0;
	    if (fb_rgb_pixel_byte_count(count, &bytes) != 0)
		return -1;
	    memset(rgb, 0, bytes);
	}
	return (ssize_t)count;
    }

    if (!fb_coords_valid(fb, x, y) || !rgb)
	return -1;

    size_t width = imgstream_fb_width(fb);
    size_t height = imgstream_fb_height(fb);
    size_t offset = (size_t)y * width + (size_t)x;
    size_t pixels_to_end = width * height - offset;
    if (count > pixels_to_end)
	count = pixels_to_end;

    size_t done = 0;
    size_t cx = (size_t)x;
    size_t cy = (size_t)y;
    while (done < count) {
	size_t run = width - cx;
	if (run > count - done)
	    run = count - done;
	if (imgstream_read_rect(fb->stream, cx, cy, run, 1, rgb + done * 3, run * 3) != 0)
	    return -1;
	done += run;
	cx = 0;
	cy++;
    }

    return (ssize_t)done;
}


ssize_t
imgstream_fb_write(imgstream_fb_t *fb, int x, int y, const unsigned char *rgb, size_t count)
{
    if (fb && fb->child_count) {
	ssize_t ret = imgstream_fb_write(fb->children[0], x, y, rgb, count);
	if (ret < 0)
	    return -1;
	for (size_t i = 1; i < fb->child_count; i++) {
	    if (imgstream_fb_write(fb->children[i], x, y, rgb, count) != ret)
		return -1;
	}
	return ret;
    }

    if (fb_is_diagnostic(fb)) {
	(void)x;
	(void)y;
	fb_diagnostic_log(fb, "write");
	fb_diagnostic_log_pixels(fb, rgb, count);
	return (ssize_t)count;
    }

    if (!fb_coords_valid(fb, x, y) || !rgb)
	return -1;

    size_t width = imgstream_fb_width(fb);
    size_t height = imgstream_fb_height(fb);
    size_t offset = (size_t)y * width + (size_t)x;
    size_t pixels_to_end = width * height - offset;
    if (count > pixels_to_end)
	count = pixels_to_end;

    size_t done = 0;
    size_t cx = (size_t)x;
    size_t cy = (size_t)y;
    while (done < count) {
	size_t run = width - cx;
	if (run > count - done)
	    run = count - done;
	if (imgstream_write_rect(fb->stream, cx, cy, run, 1, rgb + done * 3, run * 3) != 0)
	    return -1;
	done += run;
	cx = 0;
	cy++;
    }

    return (ssize_t)done;
}


static int
fb_rect_valid(const imgstream_fb_t *fb, int x, int y, int width, int height)
{
    if (!fb || !fb->stream || x < 0 || y < 0 || width <= 0 || height <= 0)
	return 0;

    if ((size_t)x >= imgstream_fb_width(fb) || (size_t)y >= imgstream_fb_height(fb))
	return 0;

    if ((size_t)width > imgstream_fb_width(fb) - (size_t)x)
	return 0;

    if ((size_t)height > imgstream_fb_height(fb) - (size_t)y)
	return 0;

    return 1;
}


int
imgstream_fb_readrect(const imgstream_fb_t *fb, int xmin, int ymin, int width, int height, unsigned char *rgb)
{
    if (fb && fb->child_count)
	return imgstream_fb_readrect(fb->children[0], xmin, ymin, width, height, rgb);

    if (fb_is_diagnostic(fb)) {
	if (width < 0 || height < 0)
	    return -1;
	if (width == 0 || height == 0)
	    return 0;
	if ((size_t)width > (size_t)INT_MAX / (size_t)height)
	    return -1;
	int pixels = width * height;
	fb_diagnostic_log(fb, "readrect");
	if (rgb) {
	    size_t bytes = 0;
	    if (fb_rgb_pixel_byte_count((size_t)pixels, &bytes) != 0)
		return -1;
	    memset(rgb, 0, bytes);
	}
	(void)xmin;
	(void)ymin;
	return pixels;
    }

    if (!fb_rect_valid(fb, xmin, ymin, width, height) || !rgb)
	return -1;

    if (imgstream_read_rect(fb->stream, (size_t)xmin, (size_t)ymin, (size_t)width, (size_t)height, rgb, (size_t)width * 3) != 0)
	return -1;

    return width * height;
}


int
imgstream_fb_writerect(imgstream_fb_t *fb, int xmin, int ymin, int width, int height, const unsigned char *rgb)
{
    if (fb && fb->child_count) {
	int ret = imgstream_fb_writerect(fb->children[0], xmin, ymin, width, height, rgb);
	if (ret < 0)
	    return -1;
	for (size_t i = 1; i < fb->child_count; i++) {
	    if (imgstream_fb_writerect(fb->children[i], xmin, ymin, width, height, rgb) != ret)
		return -1;
	}
	return ret;
    }

    if (fb_is_diagnostic(fb)) {
	if (width < 0 || height < 0)
	    return -1;
	if (width == 0 || height == 0)
	    return 0;
	if ((size_t)width > (size_t)INT_MAX / (size_t)height)
	    return -1;
	int pixels = width * height;
	fb_diagnostic_log(fb, "writerect");
	fb_diagnostic_log_pixels(fb, rgb, (size_t)pixels);
	(void)xmin;
	(void)ymin;
	return pixels;
    }

    if (!fb_rect_valid(fb, xmin, ymin, width, height) || !rgb)
	return -1;

    if (imgstream_write_rect(fb->stream, (size_t)xmin, (size_t)ymin, (size_t)width, (size_t)height, rgb, (size_t)width * 3) != 0)
	return -1;

    return width * height;
}


static int
fb_rect_pixel_count(int width, int height, int *pixels)
{
    if (!pixels || width < 0 || height < 0)
	return -1;
    if (width == 0 || height == 0) {
	*pixels = 0;
	return 0;
    }
    if ((size_t)width > (size_t)INT_MAX / (size_t)height)
	return -1;

    *pixels = width * height;
    return 0;
}


int
imgstream_fb_bwreadrect(const imgstream_fb_t *fb, int xmin, int ymin, int width, int height, unsigned char *bw)
{
    if (fb && fb->child_count)
	return imgstream_fb_bwreadrect(fb->children[0], xmin, ymin, width, height, bw);

    if (fb_is_diagnostic(fb)) {
	int pixels = 0;
	if (fb_rect_pixel_count(width, height, &pixels) != 0)
	    return -1;
	fb_diagnostic_log(fb, "bwreadrect");
	if (bw && pixels)
	    memset(bw, 0, (size_t)pixels);
	return pixels;
    }

    if (!fb_rect_valid(fb, xmin, ymin, width, height) || !bw)
	return -1;

    size_t row_bytes = 0;
    if (fb_rgb_pixel_byte_count((size_t)width, &row_bytes) != 0)
	return -1;

    unsigned char *row = (unsigned char *)bu_malloc(row_bytes, "imgstream fb bw read row");
    for (int y = 0; y < height; y++) {
	if (imgstream_read_rect(fb->stream, (size_t)xmin, (size_t)(ymin + y),
		    (size_t)width, 1, row, row_bytes) != 0) {
	    bu_free(row, "imgstream fb bw read row");
	    return -1;
	}
	for (int x = 0; x < width; x++)
	    bw[(size_t)y * (size_t)width + (size_t)x] = row[(size_t)x * 3 + 1];
    }
    bu_free(row, "imgstream fb bw read row");

    return width * height;
}


int
imgstream_fb_bwwriterect(imgstream_fb_t *fb, int xmin, int ymin, int width, int height, const unsigned char *bw)
{
    if (fb && fb->child_count) {
	int ret = imgstream_fb_bwwriterect(fb->children[0], xmin, ymin, width, height, bw);
	if (ret < 0)
	    return -1;
	for (size_t i = 1; i < fb->child_count; i++) {
	    if (imgstream_fb_bwwriterect(fb->children[i], xmin, ymin, width, height, bw) != ret)
		return -1;
	}
	return ret;
    }

    if (fb_is_diagnostic(fb)) {
	int pixels = 0;
	if (fb_rect_pixel_count(width, height, &pixels) != 0)
	    return -1;
	fb_diagnostic_log(fb, "bwwriterect");
	(void)xmin;
	(void)ymin;
	(void)bw;
	return pixels;
    }

    if (!fb_rect_valid(fb, xmin, ymin, width, height) || !bw)
	return -1;

    size_t row_bytes = 0;
    if (fb_rgb_pixel_byte_count((size_t)width, &row_bytes) != 0)
	return -1;

    unsigned char *row = (unsigned char *)bu_malloc(row_bytes, "imgstream fb bw write row");
    for (int y = 0; y < height; y++) {
	for (int x = 0; x < width; x++) {
	    unsigned char value = bw[(size_t)y * (size_t)width + (size_t)x];
	    row[(size_t)x * 3 + 0] = value;
	    row[(size_t)x * 3 + 1] = value;
	    row[(size_t)x * 3 + 2] = value;
	}
	if (imgstream_write_rect(fb->stream, (size_t)xmin, (size_t)(ymin + y),
		    (size_t)width, 1, row, row_bytes) != 0) {
	    bu_free(row, "imgstream fb bw write row");
	    return -1;
	}
    }
    bu_free(row, "imgstream fb bw write row");

    return width * height;
}


int
imgstream_fb_flush(imgstream_fb_t *fb)
{
    if (!fb)
	return -1;

    if (fb->child_count) {
	for (size_t i = 0; i < fb->child_count; i++) {
	    if (imgstream_fb_flush(fb->children[i]) != 0)
		return -1;
	}
	return 0;
    }

    fb_diagnostic_log(fb, "flush");
    if (fb->display_backed && fb->display_host.flush)
	return fb->display_host.flush(fb, fb->display_host_data);
    return fb_stream_flush_file(fb);
}


int
imgstream_fb_reset(imgstream_fb_t *fb)
{
    if (!fb)
	return 0;

    if (fb->child_count) {
	for (size_t i = 0; i < fb->child_count; i++) {
	    if (imgstream_fb_reset(fb->children[i]) != 0)
		return -1;
	}
    }

    fb_diagnostic_log(fb, "reset");
    if (fb->display_backed && fb->display_host.reset)
	return fb->display_host.reset(fb, fb->display_host_data);
    return 0;
}


int
imgstream_fb_viewport(imgstream_fb_t *fb, int left, int top, int right, int bottom)
{
    if (!fb)
	return 0;

    if (fb->child_count) {
	for (size_t i = 0; i < fb->child_count; i++) {
	    if (imgstream_fb_viewport(fb->children[i], left, top, right, bottom) != 0)
		return -1;
	}
    }

    fb_diagnostic_log(fb, "viewport");
    if (fb->display_backed && fb->display_host.viewport)
	return fb->display_host.viewport(fb, left, top, right, bottom,
		fb->display_host_data);
    return 0;
}


int
imgstream_fb_view(imgstream_fb_t *fb, int xcenter, int ycenter, int xzoom, int yzoom)
{
    if (!fb)
	return -1;

    if (fb->child_count) {
	for (size_t i = 0; i < fb->child_count; i++) {
	    if (imgstream_fb_view(fb->children[i], xcenter, ycenter, xzoom, yzoom) != 0)
		return -1;
	}
    }
    fb->view.xcenter = xcenter;
    fb->view.ycenter = ycenter;
    fb->view.xzoom = xzoom;
    fb->view.yzoom = yzoom;
    if (fb->display_backed && fb->display_host.view)
	return fb->display_host.view(fb, &fb->view, fb->display_host_data);
    return 0;
}


int
imgstream_fb_getview(const imgstream_fb_t *fb, int *xcenter, int *ycenter, int *xzoom, int *yzoom)
{
    if (!fb || !xcenter || !ycenter || !xzoom || !yzoom)
	return -1;

    *xcenter = fb->view.xcenter;
    *ycenter = fb->view.ycenter;
    *xzoom = fb->view.xzoom;
    *yzoom = fb->view.yzoom;
    return 0;
}


int
imgstream_fb_window(imgstream_fb_t *fb, int xcenter, int ycenter)
{
    if (!fb)
	return -1;

    return imgstream_fb_view(fb, xcenter, ycenter, fb->view.xzoom, fb->view.yzoom);
}


int
imgstream_fb_zoom(imgstream_fb_t *fb, int xzoom, int yzoom)
{
    if (!fb)
	return -1;

    return imgstream_fb_view(fb, fb->view.xcenter, fb->view.ycenter, xzoom, yzoom);
}


int
imgstream_fb_cursor(imgstream_fb_t *fb, int mode, int x, int y)
{
    if (!fb)
	return -1;

    if (fb->child_count) {
	for (size_t i = 0; i < fb->child_count; i++) {
	    if (imgstream_fb_cursor(fb->children[i], mode, x, y) != 0)
		return -1;
	}
    }
    fb->cursor.mode = mode;
    fb->cursor.x = x;
    fb->cursor.y = y;
    if (fb->display_backed && fb->display_host.cursor)
	return fb->display_host.cursor(fb, &fb->cursor, fb->display_host_data);
    return 0;
}


int
imgstream_fb_getcursor(const imgstream_fb_t *fb, int *mode, int *x, int *y)
{
    if (!fb || !mode || !x || !y)
	return -1;

    *mode = fb->cursor.mode;
    *x = fb->cursor.x;
    *y = fb->cursor.y;
    return 0;
}


int
imgstream_fb_scursor(imgstream_fb_t *fb, int mode, int x, int y)
{
    if (!fb)
	return -1;

    if (fb->child_count) {
	for (size_t i = 0; i < fb->child_count; i++) {
	    if (imgstream_fb_scursor(fb->children[i], mode, x, y) != 0)
		return -1;
	}
    }

    if (fb->display_backed && fb->display_host.scursor)
	return fb->display_host.scursor(fb, mode, x, y, fb->display_host_data);
    return 0;
}


int
imgstream_fb_setcursor(imgstream_fb_t *fb, const unsigned char *bits, int xbits, int ybits, int xorig, int yorig)
{
    if (!fb)
	return -1;

    if (fb->child_count) {
	for (size_t i = 0; i < fb->child_count; i++) {
	    if (imgstream_fb_setcursor(fb->children[i], bits, xbits, ybits, xorig, yorig) != 0)
		return -1;
	}
    }

    fb_diagnostic_log(fb, "setcursor");
    if (fb->display_backed && fb->display_host.setcursor)
	return fb->display_host.setcursor(fb, bits, xbits, ybits, xorig,
		yorig, fb->display_host_data);
    return 0;
}


int
imgstream_fb_rmap(const imgstream_fb_t *fb, struct imgstream_fb_colormap *cmap)
{
    if (!fb || !cmap)
	return -1;

    if (fb->child_count)
	return imgstream_fb_rmap(fb->children[0], cmap);

    fb_diagnostic_log(fb, "rmap");
    *cmap = fb->colormap;
    return 0;
}


static void
fb_diagnostic_log_colormap(const imgstream_fb_t *fb, const struct imgstream_fb_colormap *cmap)
{
    if (!fb || fb->diagnostic_kind != IMGSTREAM_FB_DIAGNOSTIC_DEBUG ||
	    !(fb->diagnostic_flags & IMGSTREAM_FB_DIAGNOSTIC_LOG_CMAP) || !cmap)
	return;

    for (size_t i = 0; i < 256; i++) {
	bu_log("%3zu: [ 0x%4x, 0x%4x, 0x%4x ]\n", i,
		(unsigned int)cmap->red[i],
		(unsigned int)cmap->green[i],
		(unsigned int)cmap->blue[i]);
    }
}


int
imgstream_fb_wmap(imgstream_fb_t *fb, const struct imgstream_fb_colormap *cmap)
{
    if (!fb)
	return -1;

    if (fb->child_count) {
	for (size_t i = 0; i < fb->child_count; i++) {
	    if (imgstream_fb_wmap(fb->children[i], cmap) != 0)
		return -1;
	}
	if (cmap)
	    fb->colormap = *cmap;
	else
	    imgstream_fb_colormap_linear(&fb->colormap);
	return 0;
    }

    fb_diagnostic_log(fb, "wmap");
    fb_diagnostic_log_colormap(fb, cmap);
    if (cmap)
	fb->colormap = *cmap;
    else
	imgstream_fb_colormap_linear(&fb->colormap);
    return 0;
}


int
imgstream_fb_poll(imgstream_fb_t *fb)
{
    if (!fb)
	return -1;

    if (fb->child_count) {
	for (size_t i = 0; i < fb->child_count; i++) {
	    if (imgstream_fb_poll(fb->children[i]) != 0)
		return -1;
	}
    }

    if (fb->display_backed && fb->display_host.poll)
	return fb->display_host.poll(fb, fb->display_host_data);

    return 0;
}


long
imgstream_fb_poll_rate(const imgstream_fb_t *fb)
{
    if (fb && fb->display_backed && fb->display_host.poll_rate)
	return fb->display_host.poll_rate(fb, fb->display_host_data);

    return 0;
}
