/*                    F R A M E B U F F E R . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/framebuffer.h */

#ifndef BRLOBOL_FRAMEBUFFER_H
#define BRLOBOL_FRAMEBUFFER_H

#include "brlobol/defines.h"

#include <stddef.h>
#include <sys/types.h>

class BRLObolWindowHost;

struct imgstream_fb;
typedef struct imgstream_fb imgstream_fb_t;
struct imgstream_fb_colormap;

struct BRLOBOL_EXPORT BRLObolFramebufferInfo {
    int max_width;
    int max_height;
    int width;
    int height;
};

struct BRLObolFramebufferStreamPrivate;

class BRLOBOL_EXPORT BRLObolFramebufferStream {
public:
    explicit BRLObolFramebufferStream(BRLObolWindowHost *host = NULL);
    ~BRLObolFramebufferStream(void);

    void setHost(BRLObolWindowHost *host);
    BRLObolWindowHost *host(void) const;

    int configure(int width, int height, const char *spec = NULL);
    int ensure(void);
    void close(void);

    imgstream_fb_t *framebuffer(void) const;
    int info(BRLObolFramebufferInfo *info);

    int clear(const unsigned char rgb[3]);
    ssize_t read(int x, int y, unsigned char *rgb, size_t count);
    ssize_t write(int x, int y, const unsigned char *rgb, size_t count);
    int readrect(int xmin, int ymin, int width, int height, unsigned char *rgb);
    int writerect(int xmin, int ymin, int width, int height,
	const unsigned char *rgb);
    int bwreadrect(int xmin, int ymin, int width, int height,
	unsigned char *bw);
    int bwwriterect(int xmin, int ymin, int width, int height,
	const unsigned char *bw);

    int cursor(int mode, int x, int y);
    int getcursor(int *mode, int *x, int *y);
    int setcursor(const unsigned char *bits, int xbits, int ybits,
	int xorig, int yorig);
    int scursor(int mode, int x, int y);
    int window(int xcenter, int ycenter);
    int zoom(int xzoom, int yzoom);
    int view(int xcenter, int ycenter, int xzoom, int yzoom);
    int getview(int *xcenter, int *ycenter, int *xzoom, int *yzoom);

    int rmap(struct imgstream_fb_colormap *cmap);
    int wmap(const struct imgstream_fb_colormap *cmap);
    int present(void);
    int flush(void);
    int poll(void);

private:
    friend class BRLObolWindowHost;

    void detachHost(BRLObolWindowHost *host);
    BRLObolFramebufferStreamPrivate *p;
};

#endif /* BRLOBOL_FRAMEBUFFER_H */
