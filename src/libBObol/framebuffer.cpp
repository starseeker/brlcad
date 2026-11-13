/*                 F R A M E B U F F E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BFramebuffer.h"
#include "BObol/BWindowHost.h"

#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wshadow"
#endif
#include "imgstream/fb_compat.h"
#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

#include <string>

struct BObolFramebufferStreamPrivate {
    BObolFramebufferStreamPrivate(BObolWindowHost *h) :
	host(h),
	fb(NULL),
	width(0),
	height(0),
	spec("/dev/mem"),
	composition(BOBOL_FRAMEBUFFER_COMPOSITION_OVERLAY)
    {
    }

    BObolWindowHost *host;
    imgstream_fb_t *fb;
    int width;
    int height;
    std::string spec;
    BObolFramebufferComposition composition;
};

BObolFramebufferStream::BObolFramebufferStream(BObolWindowHost *host) :
    p(new BObolFramebufferStreamPrivate(host))
{
    if (host)
	host->registerFramebufferStream(this);
}

BObolFramebufferStream::~BObolFramebufferStream(void)
{
    this->close();
    if (this->p->host)
	this->p->host->unregisterFramebufferStream(this);
    delete this->p;
    this->p = NULL;
}

int
BObolFramebufferStream::setHost(BObolWindowHost *host)
{
    return this->setHost(host, this->p->composition);
}

int
BObolFramebufferStream::setHost(BObolWindowHost *host,
	BObolFramebufferComposition composition)
{
    if (composition < BOBOL_FRAMEBUFFER_COMPOSITION_OFF ||
	composition > BOBOL_FRAMEBUFFER_COMPOSITION_INTERLAY)
	return -1;
    if (this->p->host == host)
	return this->setComposition(composition);

    BObolWindowHost *old_host = this->p->host;
    if (this->p->fb && host) {
	imgstream_fb_spec_info_t info;
	imgstream_fb_spec_info_t *ip = NULL;
	if (imgstream_fb_spec_info(this->p->spec.c_str(), &info) == 0)
	    ip = &info;
	/* Attach first: a failed replacement must leave the live host and its
	 * retained image nodes untouched. */
	if (host->openFramebuffer(this->p->fb, ip) != 0)
	    return -1;
	if (host->setFramebufferComposition(this->p->fb, composition) != 0) {
	    host->closeFramebuffer(this->p->fb);
	    return -1;
	}
    }

    if (old_host && this->p->fb)
	old_host->closeFramebuffer(this->p->fb);
    if (old_host)
	old_host->unregisterFramebufferStream(this);
    this->p->host = host;
    this->p->composition = composition;
    if (host)
	host->registerFramebufferStream(this);
    return 0;
}

BObolWindowHost *
BObolFramebufferStream::host(void) const
{
    return this->p->host;
}

int
BObolFramebufferStream::setComposition(
	BObolFramebufferComposition composition)
{
    if (composition < BOBOL_FRAMEBUFFER_COMPOSITION_OFF ||
	composition > BOBOL_FRAMEBUFFER_COMPOSITION_INTERLAY)
	return -1;
    if (this->p->fb && this->p->host &&
	this->p->host->setFramebufferComposition(this->p->fb,
	    composition) != 0)
	return -1;
    this->p->composition = composition;
    return 0;
}

BObolFramebufferComposition
BObolFramebufferStream::composition(void) const
{
    return this->p->composition;
}

void
BObolFramebufferStream::detachHost(BObolWindowHost *host)
{
    if (this->p->host == host)
	this->p->host = NULL;
}

int
BObolFramebufferStream::configure(int width, int height, const char *spec)
{
    if (width <= 0 || height <= 0)
	return -1;

    std::string nextSpec = (spec && spec[0]) ? spec : "/dev/mem";
    if (this->p->fb && this->p->width == width &&
	this->p->height == height && this->p->spec == nextSpec)
	return 0;

    this->close();
    this->p->width = width;
    this->p->height = height;
    this->p->spec = nextSpec;
    return 0;
}

int
BObolFramebufferStream::ensure(void)
{
    if (!this->p->host || this->p->width <= 0 || this->p->height <= 0)
	return -1;
    if (this->p->fb)
	return 0;

    imgstream_fb_t *fb = imgstream_fb_open(this->p->spec.c_str(),
	(size_t)this->p->width, (size_t)this->p->height);
    if (!fb)
	return -1;

    imgstream_fb_spec_info_t info;
    imgstream_fb_spec_info_t *ip = NULL;
    if (imgstream_fb_spec_info(this->p->spec.c_str(), &info) == 0)
	ip = &info;

    if (this->p->host->openFramebuffer(fb, ip) != 0) {
	imgstream_fb_close(fb);
	return -1;
    }

    if (this->p->host->setFramebufferComposition(fb,
	    this->p->composition) != 0) {
	this->p->host->closeFramebuffer(fb);
	imgstream_fb_close(fb);
	return -1;
    }

    this->p->fb = fb;
    return 0;
}

void
BObolFramebufferStream::close(void)
{
    if (!this->p->fb)
	return;
    if (this->p->host)
	this->p->host->closeFramebuffer(this->p->fb);
    imgstream_fb_close(this->p->fb);
    this->p->fb = NULL;
}

imgstream_fb_t *
BObolFramebufferStream::framebuffer(void) const
{
    return this->p->fb;
}

int
BObolFramebufferStream::info(BObolFramebufferInfo *info)
{
    if (!info || this->ensure() != 0)
	return -1;

    info->max_width = this->p->width;
    info->max_height = this->p->height;
    info->width = this->p->width;
    info->height = this->p->height;
    return 0;
}

int
BObolFramebufferStream::clear(const unsigned char rgb[3])
{
    if (this->ensure() != 0)
	return -1;
    return imgstream_fb_clear(this->p->fb, rgb);
}

ssize_t
BObolFramebufferStream::read(int x, int y, unsigned char *rgb, size_t count)
{
    return this->ensure() == 0 ?
	imgstream_fb_read(this->p->fb, x, y, rgb, count) : -1;
}

ssize_t
BObolFramebufferStream::write(int x, int y, const unsigned char *rgb,
				size_t count)
{
    if (this->ensure() != 0)
	return -1;
    return imgstream_fb_write(this->p->fb, x, y, rgb, count);
}

int
BObolFramebufferStream::readrect(int xmin, int ymin, int width, int height,
				   unsigned char *rgb)
{
    return this->ensure() == 0 ?
	imgstream_fb_readrect(this->p->fb, xmin, ymin, width, height, rgb) : -1;
}

int
BObolFramebufferStream::writerect(int xmin, int ymin, int width, int height,
				    const unsigned char *rgb)
{
    if (this->ensure() != 0)
	return -1;
    return imgstream_fb_writerect(this->p->fb, xmin, ymin, width, height,
	rgb);
}

int
BObolFramebufferStream::bwreadrect(int xmin, int ymin, int width, int height,
				     unsigned char *bw)
{
    return this->ensure() == 0 ?
	imgstream_fb_bwreadrect(this->p->fb, xmin, ymin, width, height, bw) : -1;
}

int
BObolFramebufferStream::bwwriterect(int xmin, int ymin, int width, int height,
				      const unsigned char *bw)
{
    if (this->ensure() != 0)
	return -1;
    return imgstream_fb_bwwriterect(this->p->fb, xmin, ymin, width, height,
	bw);
}

int
BObolFramebufferStream::cursor(int mode, int x, int y)
{
    if (this->ensure() != 0)
	return -1;
    return imgstream_fb_cursor(this->p->fb, mode, x, y);
}

int
BObolFramebufferStream::getcursor(int *mode, int *x, int *y)
{
    return this->ensure() == 0 ?
	imgstream_fb_getcursor(this->p->fb, mode, x, y) : -1;
}

int
BObolFramebufferStream::setcursor(const unsigned char *bits, int xbits,
				    int ybits, int xorig, int yorig)
{
    if (this->ensure() != 0)
	return -1;
    return imgstream_fb_setcursor(this->p->fb, bits, xbits, ybits,
	xorig, yorig);
}

int
BObolFramebufferStream::scursor(int mode, int x, int y)
{
    if (this->ensure() != 0)
	return -1;
    return imgstream_fb_scursor(this->p->fb, mode, x, y);
}

int
BObolFramebufferStream::window(int xcenter, int ycenter)
{
    if (this->ensure() != 0)
	return -1;
    return imgstream_fb_window(this->p->fb, xcenter, ycenter);
}

int
BObolFramebufferStream::zoom(int xzoom, int yzoom)
{
    if (this->ensure() != 0)
	return -1;
    return imgstream_fb_zoom(this->p->fb, xzoom, yzoom);
}

int
BObolFramebufferStream::view(int xcenter, int ycenter, int xzoom, int yzoom)
{
    if (this->ensure() != 0)
	return -1;
    return imgstream_fb_view(this->p->fb, xcenter, ycenter, xzoom, yzoom);
}

int
BObolFramebufferStream::getview(int *xcenter, int *ycenter, int *xzoom,
				  int *yzoom)
{
    return this->ensure() == 0 ?
	imgstream_fb_getview(this->p->fb, xcenter, ycenter, xzoom, yzoom) : -1;
}

int
BObolFramebufferStream::rmap(struct imgstream_fb_colormap *cmap)
{
    return this->ensure() == 0 ? imgstream_fb_rmap(this->p->fb, cmap) : -1;
}

int
BObolFramebufferStream::wmap(const struct imgstream_fb_colormap *cmap)
{
    return this->ensure() == 0 ? imgstream_fb_wmap(this->p->fb, cmap) : -1;
}

int
BObolFramebufferStream::present(void)
{
    if (this->ensure() != 0)
	return -1;

    int ret = 0;
    if (this->p->host->flushFramebuffer(this->p->fb) != 0)
	ret = -1;

    imgstream_fb_view_t view_state;
    if (imgstream_fb_getview(this->p->fb, &view_state.xcenter,
	    &view_state.ycenter, &view_state.xzoom,
	    &view_state.yzoom) == 0) {
	if (this->p->host->setFramebufferView(this->p->fb, &view_state) != 0)
	    ret = -1;
    }

    imgstream_fb_cursor_t cursor_state;
    if (imgstream_fb_getcursor(this->p->fb, &cursor_state.mode,
	    &cursor_state.x, &cursor_state.y) == 0) {
	if (this->p->host->setFramebufferCursor(this->p->fb,
		&cursor_state) != 0)
	    ret = -1;
    }

    return ret;
}

int
BObolFramebufferStream::flush(void)
{
    if (this->present() != 0)
	return -1;
    return imgstream_fb_flush(this->p->fb);
}

int
BObolFramebufferStream::poll(void)
{
    return this->p->host ? this->p->host->poll(NULL) : -1;
}
