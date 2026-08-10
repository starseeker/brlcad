/*                    F B _ C O M P A T . H
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
/** @file imgstream/fb_compat.h */
#ifndef IMGSTREAM_FB_COMPAT_H
#define IMGSTREAM_FB_COMPAT_H

#include "common.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "imgstream/stream.h"

__BEGIN_DECLS

struct imgstream_fb;
typedef struct imgstream_fb imgstream_fb_t;
struct icv_image;

struct imgstream_fb_import_options {
    int file_xoff;
    int file_yoff;
    int screen_xoff;
    int screen_yoff;
    int clear;
    int zoom;
    int inverse;
};

#define IMGSTREAM_FB_IMPORT_OPTIONS_INIT {0, 0, 0, 0, 0, 0, 0}

enum imgstream_fb_spec_kind {
    IMGSTREAM_FB_SPEC_MEMORY,
    IMGSTREAM_FB_SPEC_FILE,
    IMGSTREAM_FB_SPEC_FANOUT,
    IMGSTREAM_FB_SPEC_REMOTE,
    IMGSTREAM_FB_SPEC_DISPLAY,
    IMGSTREAM_FB_SPEC_DIAGNOSTIC,
    IMGSTREAM_FB_SPEC_UNSUPPORTED
};

enum imgstream_fb_remote_kind {
    IMGSTREAM_FB_REMOTE_NONE,
    IMGSTREAM_FB_REMOTE_TCP,
    IMGSTREAM_FB_REMOTE_IPC
};

enum imgstream_fb_display_kind {
    IMGSTREAM_FB_DISPLAY_NONE,
    IMGSTREAM_FB_DISPLAY_QTGL,
    IMGSTREAM_FB_DISPLAY_OGL,
    IMGSTREAM_FB_DISPLAY_WGL,
    IMGSTREAM_FB_DISPLAY_X,
    IMGSTREAM_FB_DISPLAY_SWRAST
};

enum imgstream_fb_diagnostic_kind {
    IMGSTREAM_FB_DIAGNOSTIC_NONE = 0,
    IMGSTREAM_FB_DIAGNOSTIC_NULL,
    IMGSTREAM_FB_DIAGNOSTIC_DEBUG,
    IMGSTREAM_FB_DIAGNOSTIC_TXT
};

struct imgstream_fb_spec_info_s {
    enum imgstream_fb_spec_kind kind;
    enum imgstream_fb_remote_kind remote;
    enum imgstream_fb_display_kind display;
    enum imgstream_fb_diagnostic_kind diagnostic;
    /* Borrowed string spans; valid while the input spec string is valid. */
    const char *target;
    size_t target_len;
    const char *host;
    size_t host_len;
    int port;
    const char *device;
    size_t device_len;
};
typedef struct imgstream_fb_spec_info_s imgstream_fb_spec_info_t;

struct imgstream_fb_view_s {
    int xcenter;
    int ycenter;
    int xzoom;
    int yzoom;
};
typedef struct imgstream_fb_view_s imgstream_fb_view_t;

struct imgstream_fb_cursor_s {
    int mode;
    int x;
    int y;
};
typedef struct imgstream_fb_cursor_s imgstream_fb_cursor_t;

struct imgstream_fb_display_host {
    int (*open)(imgstream_fb_t *fb, const imgstream_fb_spec_info_t *info, void *data);
    void (*close)(imgstream_fb_t *fb, void *data);
    int (*flush)(imgstream_fb_t *fb, void *data);
    int (*reset)(imgstream_fb_t *fb, void *data);
    int (*viewport)(imgstream_fb_t *fb, int left, int top, int right, int bottom, void *data);
    int (*view)(imgstream_fb_t *fb, const imgstream_fb_view_t *view, void *data);
    int (*cursor)(imgstream_fb_t *fb, const imgstream_fb_cursor_t *cursor, void *data);
    int (*scursor)(imgstream_fb_t *fb, int mode, int x, int y, void *data);
    int (*setcursor)(imgstream_fb_t *fb, const unsigned char *bits, int xbits, int ybits, int xorig, int yorig, void *data);
    int (*poll)(imgstream_fb_t *fb, void *data);
    long (*poll_rate)(const imgstream_fb_t *fb, void *data);
};

struct imgstream_fb_remote_options {
    uint32_t struct_size;
    const char *auth_token;
    int use_tls;
    /** Lowercase/uppercase SHA-256 fingerprint of the expected TLS certificate. */
    const char *tls_server_sha256;
    /** Explicitly permit unpinned TLS or non-loopback plaintext TCP. */
    int allow_insecure;
};

#define IMGSTREAM_FB_REMOTE_OPTIONS_INIT { \
	sizeof(struct imgstream_fb_remote_options), NULL, 0, NULL, 0 }

struct imgstream_fb_colormap {
    uint16_t red[256];
    uint16_t green[256];
    uint16_t blue[256];
};

IMGSTREAM_EXPORT enum imgstream_fb_spec_kind imgstream_fb_spec_kind(const char *spec);
IMGSTREAM_EXPORT int imgstream_fb_spec_info(const char *spec, imgstream_fb_spec_info_t *info);
IMGSTREAM_EXPORT const char *imgstream_fb_spec_kind_name(enum imgstream_fb_spec_kind kind);
IMGSTREAM_EXPORT int imgstream_fb_spec_supported(const char *spec);
IMGSTREAM_EXPORT void imgstream_fb_colormap_linear(struct imgstream_fb_colormap *cmap);
IMGSTREAM_EXPORT int imgstream_fb_colormap_is_linear(const struct imgstream_fb_colormap *cmap);
IMGSTREAM_EXPORT int imgstream_image_name_size(size_t *width, size_t *height,
	const char *name);
IMGSTREAM_EXPORT int imgstream_image_file_size(size_t *width, size_t *height,
	const char *filename, int bytes_per_pixel);
IMGSTREAM_EXPORT int imgstream_image_size(size_t *width, size_t *height,
	size_t pixel_count);

IMGSTREAM_EXPORT imgstream_fb_t *imgstream_fb_open(const char *spec, size_t width, size_t height);
IMGSTREAM_EXPORT imgstream_fb_t *imgstream_fb_open_remote(const char *spec,
	size_t width, size_t height,
	const struct imgstream_fb_remote_options *options);
IMGSTREAM_EXPORT imgstream_fb_t *imgstream_fb_open_display(const char *spec, size_t width, size_t height,
	const struct imgstream_fb_display_host *host, void *data);
/*
 * Disconnect the display host, invoking its close callback once now while the
 * pixel stream remains usable.  Later framebuffer close does not call it.
 */
IMGSTREAM_EXPORT int imgstream_fb_detach_display_host(imgstream_fb_t *fb);
IMGSTREAM_EXPORT void imgstream_fb_close(imgstream_fb_t *fb);
/* /dev/null is a zero-storage diagnostic sink and returns NULL here. */
IMGSTREAM_EXPORT imgstream_t *imgstream_fb_stream(imgstream_fb_t *fb);
IMGSTREAM_EXPORT const imgstream_t *imgstream_fb_cstream(const imgstream_fb_t *fb);
IMGSTREAM_EXPORT const char *imgstream_fb_name(const imgstream_fb_t *fb);

IMGSTREAM_EXPORT size_t imgstream_fb_width(const imgstream_fb_t *fb);
IMGSTREAM_EXPORT size_t imgstream_fb_height(const imgstream_fb_t *fb);

IMGSTREAM_EXPORT int imgstream_fb_clear(imgstream_fb_t *fb, const unsigned char *rgb);
IMGSTREAM_EXPORT ssize_t imgstream_fb_read(const imgstream_fb_t *fb, int x, int y, unsigned char *rgb, size_t count);
IMGSTREAM_EXPORT ssize_t imgstream_fb_write(imgstream_fb_t *fb, int x, int y, const unsigned char *rgb, size_t count);
IMGSTREAM_EXPORT int imgstream_fb_readrect(const imgstream_fb_t *fb, int xmin, int ymin, int width, int height, unsigned char *rgb);
IMGSTREAM_EXPORT int imgstream_fb_writerect(imgstream_fb_t *fb, int xmin, int ymin, int width, int height, const unsigned char *rgb);
IMGSTREAM_EXPORT int imgstream_fb_bwreadrect(const imgstream_fb_t *fb, int xmin, int ymin, int width, int height, unsigned char *bw);
IMGSTREAM_EXPORT int imgstream_fb_bwwriterect(imgstream_fb_t *fb, int xmin, int ymin, int width, int height, const unsigned char *bw);
IMGSTREAM_EXPORT int imgstream_fb_flush(imgstream_fb_t *fb);
IMGSTREAM_EXPORT int imgstream_fb_export_pix_fp(const imgstream_fb_t *fb,
	FILE *fp, size_t width, size_t height, int crunch, int inverse);
IMGSTREAM_EXPORT int imgstream_fb_import_pix_fd(imgstream_fb_t *fb, int fd,
	const char *filename, size_t width, size_t height, int autosize,
	const struct imgstream_fb_import_options *options);
IMGSTREAM_EXPORT int imgstream_fb_import_icv(imgstream_fb_t *fb,
	const struct icv_image *image,
	const struct imgstream_fb_import_options *options);

IMGSTREAM_EXPORT int imgstream_fb_reset(imgstream_fb_t *fb);
IMGSTREAM_EXPORT int imgstream_fb_viewport(imgstream_fb_t *fb, int left, int top, int right, int bottom);
IMGSTREAM_EXPORT int imgstream_fb_view(imgstream_fb_t *fb, int xcenter, int ycenter, int xzoom, int yzoom);
IMGSTREAM_EXPORT int imgstream_fb_getview(const imgstream_fb_t *fb, int *xcenter, int *ycenter, int *xzoom, int *yzoom);
IMGSTREAM_EXPORT int imgstream_fb_window(imgstream_fb_t *fb, int xcenter, int ycenter);
IMGSTREAM_EXPORT int imgstream_fb_zoom(imgstream_fb_t *fb, int xzoom, int yzoom);
IMGSTREAM_EXPORT int imgstream_fb_cursor(imgstream_fb_t *fb, int mode, int x, int y);
IMGSTREAM_EXPORT int imgstream_fb_getcursor(const imgstream_fb_t *fb, int *mode, int *x, int *y);
IMGSTREAM_EXPORT int imgstream_fb_scursor(imgstream_fb_t *fb, int mode, int x, int y);
IMGSTREAM_EXPORT int imgstream_fb_setcursor(imgstream_fb_t *fb, const unsigned char *bits, int xbits, int ybits, int xorig, int yorig);
IMGSTREAM_EXPORT int imgstream_fb_rmap(const imgstream_fb_t *fb, struct imgstream_fb_colormap *cmap);
IMGSTREAM_EXPORT int imgstream_fb_wmap(imgstream_fb_t *fb, const struct imgstream_fb_colormap *cmap);
IMGSTREAM_EXPORT int imgstream_fb_poll(imgstream_fb_t *fb);
IMGSTREAM_EXPORT long imgstream_fb_poll_rate(const imgstream_fb_t *fb);

__END_DECLS

#endif /* IMGSTREAM_FB_COMPAT_H */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
