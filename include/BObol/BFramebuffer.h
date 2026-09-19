/*                B F R A M E B U F F E R . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BFramebuffer.h */

#ifndef BOBOL_BFRAMEBUFFER_H
#define BOBOL_BFRAMEBUFFER_H

#include "BObol/BDefines.h"

#include <stddef.h>
#include <sys/types.h>

class BObolWindowHost;

struct imgstream_fb;
typedef struct imgstream_fb imgstream_fb_t;
struct imgstream_fb_colormap;

enum BObolFramebufferComposition {
    BOBOL_FRAMEBUFFER_COMPOSITION_OFF = 0,
    BOBOL_FRAMEBUFFER_COMPOSITION_OVERLAY = 1,
    BOBOL_FRAMEBUFFER_COMPOSITION_UNDERLAY = 2,
    /* Between the CAD scene and view-local screen features. */
    BOBOL_FRAMEBUFFER_COMPOSITION_INTERLAY = 3
};

struct BOBOL_EXPORT BObolFramebufferInfo {
    int max_width;
    int max_height;
    int width;
    int height;
};

struct BObolFramebufferStreamPrivate;

class BOBOL_EXPORT BObolFramebufferStream {
public:
    explicit BObolFramebufferStream(BObolWindowHost *host = NULL);
    ~BObolFramebufferStream(void);

    /**
     * Move this stream's live framebuffer attachment to another host.
     *
     * The imgstream framebuffer and its pixel/view/cursor state are preserved.
     * When the target rejects the attachment, the old host remains active.
     * A NULL target detaches presentation without closing the stream.
     *
     * @return 0 on success, -1 when the target host rejects the attachment.
     */
    int setHost(BObolWindowHost *host);
    int setHost(BObolWindowHost *host,
	BObolFramebufferComposition composition);
    BObolWindowHost *host(void) const;

    /**
     * Set the retained image layer.  A live update is transactional: failure
     * leaves both the previous host attachment and policy unchanged.
     */
    int setComposition(BObolFramebufferComposition composition);
    BObolFramebufferComposition composition(void) const;

    int configure(int width, int height, const char *spec = NULL);
    int ensure(void);
    void close(void);

    imgstream_fb_t *framebuffer(void) const;
    int info(BObolFramebufferInfo *info);

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
    friend class BObolWindowHost;

    void detachHost(BObolWindowHost *host);
    BObolFramebufferStreamPrivate *p;
};

#endif /* BOBOL_BFRAMEBUFFER_H */
