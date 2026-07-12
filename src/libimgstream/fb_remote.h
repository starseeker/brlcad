/* Private libimgstream remote framebuffer client contract. */

#ifndef LIBIMGSTREAM_FB_REMOTE_H
#define LIBIMGSTREAM_FB_REMOTE_H

#include "imgstream/fb_compat.h"

struct imgstream_fb_remote;

int imgstream_fb_remote_open(const imgstream_fb_spec_info_t *info,
    size_t requested_width, size_t requested_height,
    const struct imgstream_fb_remote_options *options,
    struct imgstream_fb_remote **remote, size_t *width, size_t *height);
void imgstream_fb_remote_close(struct imgstream_fb_remote *remote);
int imgstream_fb_remote_clear(struct imgstream_fb_remote *remote,
    const unsigned char *rgb);
ssize_t imgstream_fb_remote_read(struct imgstream_fb_remote *remote,
    int x, int y, unsigned char *rgb, size_t count);
ssize_t imgstream_fb_remote_write(struct imgstream_fb_remote *remote,
    int x, int y, const unsigned char *rgb, size_t count);
int imgstream_fb_remote_readrect(struct imgstream_fb_remote *remote,
    int x, int y, int width, int height, unsigned char *rgb);
int imgstream_fb_remote_writerect(struct imgstream_fb_remote *remote,
    int x, int y, int width, int height, const unsigned char *rgb);
int imgstream_fb_remote_bwreadrect(struct imgstream_fb_remote *remote,
    int x, int y, int width, int height, unsigned char *bw);
int imgstream_fb_remote_bwwriterect(struct imgstream_fb_remote *remote,
    int x, int y, int width, int height, const unsigned char *bw);
int imgstream_fb_remote_view(struct imgstream_fb_remote *remote,
    int xcenter, int ycenter, int xzoom, int yzoom);
int imgstream_fb_remote_getview(struct imgstream_fb_remote *remote,
    int *xcenter, int *ycenter, int *xzoom, int *yzoom);
int imgstream_fb_remote_cursor(struct imgstream_fb_remote *remote,
    int mode, int x, int y);
int imgstream_fb_remote_getcursor(struct imgstream_fb_remote *remote,
    int *mode, int *x, int *y);
int imgstream_fb_remote_scursor(struct imgstream_fb_remote *remote,
    int mode, int x, int y);
int imgstream_fb_remote_setcursor(struct imgstream_fb_remote *remote,
    const unsigned char *bits, int xbits, int ybits, int xorig, int yorig);
int imgstream_fb_remote_rmap(struct imgstream_fb_remote *remote,
    struct imgstream_fb_colormap *cmap);
int imgstream_fb_remote_wmap(struct imgstream_fb_remote *remote,
    const struct imgstream_fb_colormap *cmap);
int imgstream_fb_remote_flush(struct imgstream_fb_remote *remote);
int imgstream_fb_remote_poll(struct imgstream_fb_remote *remote);

#endif /* LIBIMGSTREAM_FB_REMOTE_H */
