/*               Q G O B O L C O N T E X T M A N A G E R . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolContextManager.h
 *
 * Private qtcad Obol ContextManager backed by Qt-provided OpenGL contexts.
 * QgGL renders directly through the current QOpenGLWidget context; this
 * manager supplies matching Qt offscreen contexts for viewport image readback.
 * QgSW pins its separate manager to OSMesa; hardware readback may fall back to
 * OSMesa when Qt cannot create an offscreen context.
 */

#ifndef QGOBOLCONTEXTMANAGER_H
#define QGOBOLCONTEXTMANAGER_H

#include "common.h"

#include "bu/str.h"

#include <Inventor/SoDB.h>

#include <cstdlib>

#include <QOffscreenSurface>
#include <QOpenGLContext>
#include <QSurface>
#include <QSurfaceFormat>

static inline bool
qg_obol_context_force_osmesa(void)
{
    const char *value = std::getenv("QTCAD_OBOL_FORCE_OSMESA");
    return value ? bu_str_true(value) != 0 : false;
}

class QgObolContextManager : public SoDB::ContextManager {
private:
    struct QgObolContext {
	enum ContextKind {
	    QT_CONTEXT,
	    FALLBACK_CONTEXT
	};

	ContextKind kind = QT_CONTEXT;
	QOffscreenSurface *surface = NULL;
	QOpenGLContext *context = NULL;
	QOpenGLContext *previousContext = NULL;
	QSurface *previousSurface = NULL;
	void *fallbackContext = NULL;
    };

public:
    explicit QgObolContextManager(bool softwareOnly = false) :
	forceSoftware(softwareOnly),
	fallbackManager(SoDB::createOSMesaContextManager())
    {
    }

    virtual ~QgObolContextManager(void)
    {
	delete this->fallbackManager;
	this->fallbackManager = NULL;
    }

    virtual void *createOffscreenContext(unsigned int width, unsigned int height)
    {
	QgObolContext *ctx = new QgObolContext;

	if (this->forceSoftware || qg_obol_context_force_osmesa())
	    return this->createFallbackContext(ctx, width, height);

	QSurfaceFormat format = QSurfaceFormat::defaultFormat();
	if (format.renderableType() == QSurfaceFormat::DefaultRenderableType)
	    format.setRenderableType(QSurfaceFormat::OpenGL);
	if (format.depthBufferSize() < 24)
	    format.setDepthBufferSize(24);
	if (format.stencilBufferSize() < 8)
	    format.setStencilBufferSize(8);

	ctx->context = new QOpenGLContext;
	ctx->context->setFormat(format);
	QOpenGLContext *shareContext = QOpenGLContext::globalShareContext();
	if (!shareContext)
	    shareContext = QOpenGLContext::currentContext();
	if (shareContext)
	    ctx->context->setShareContext(shareContext);
	if (!ctx->context->create()) {
	    delete ctx->context;
	    ctx->context = NULL;
	    return this->createFallbackContext(ctx, width, height);
	}

	ctx->surface = new QOffscreenSurface;
	ctx->surface->setFormat(ctx->context->format());
	ctx->surface->create();
	if (!ctx->surface->isValid()) {
	    delete ctx->context;
	    ctx->context = NULL;
	    delete ctx->surface;
	    ctx->surface = NULL;
	    return this->createFallbackContext(ctx, width, height);
	}

	ctx->kind = QgObolContext::QT_CONTEXT;
	return ctx;
    }

    virtual SbBool makeContextCurrent(void *context)
    {
	QgObolContext *ctx = static_cast<QgObolContext *>(context);
	if (ctx && ctx->kind == QgObolContext::FALLBACK_CONTEXT)
	    return this->fallbackManager ?
		this->fallbackManager->makeContextCurrent(ctx->fallbackContext) : FALSE;
	if (!ctx || !ctx->context || !ctx->surface)
	    return FALSE;

	ctx->previousContext = QOpenGLContext::currentContext();
	ctx->previousSurface = ctx->previousContext ? ctx->previousContext->surface() : NULL;
	return ctx->context->makeCurrent(ctx->surface) ? TRUE : FALSE;
    }

    virtual void restorePreviousContext(void *context)
    {
	QgObolContext *ctx = static_cast<QgObolContext *>(context);
	if (ctx && ctx->kind == QgObolContext::FALLBACK_CONTEXT) {
	    if (this->fallbackManager)
		this->fallbackManager->restorePreviousContext(ctx->fallbackContext);
	    return;
	}
	if (!ctx || !ctx->context)
	    return;

	if (ctx->previousContext && ctx->previousSurface)
	    ctx->previousContext->makeCurrent(ctx->previousSurface);
	else
	    ctx->context->doneCurrent();
    }

    virtual void destroyContext(void *context)
    {
	QgObolContext *ctx = static_cast<QgObolContext *>(context);
	if (!ctx)
	    return;

	if (ctx->kind == QgObolContext::FALLBACK_CONTEXT) {
	    if (this->fallbackManager)
		this->fallbackManager->destroyContext(ctx->fallbackContext);
	    delete ctx;
	    return;
	}

	if (ctx->context)
	    ctx->context->doneCurrent();
	delete ctx->context;
	delete ctx->surface;
	delete ctx;
    }

    virtual SbBool isOSMesaContext(void *context)
    {
	QgObolContext *ctx = static_cast<QgObolContext *>(context);
	if (ctx && ctx->kind == QgObolContext::FALLBACK_CONTEXT &&
		this->fallbackManager)
	    return this->fallbackManager->isOSMesaContext(ctx->fallbackContext);
	return FALSE;
    }

    virtual void maxOffscreenDimensions(unsigned int &width,
	    unsigned int &height) const
    {
	if (this->fallbackManager) {
	    this->fallbackManager->maxOffscreenDimensions(width, height);
	    if (width > 0 && height > 0)
		return;
	}
	width = 16384;
	height = 16384;
    }

    virtual void getActualSurfaceSize(void *context,
	    unsigned int &width,
	    unsigned int &height) const
    {
	QgObolContext *ctx = static_cast<QgObolContext *>(context);
	if (ctx && ctx->kind == QgObolContext::FALLBACK_CONTEXT &&
		this->fallbackManager) {
	    this->fallbackManager->getActualSurfaceSize(ctx->fallbackContext,
		    width, height);
	    return;
	}
	width = 0;
	height = 0;
    }

    virtual void *getProcAddress(const char *funcName)
    {
	QOpenGLContext *ctx = QOpenGLContext::currentContext();
	if (ctx && funcName)
	    return reinterpret_cast<void *>(ctx->getProcAddress(funcName));
	return this->fallbackManager ? this->fallbackManager->getProcAddress(funcName) : NULL;
    }

    virtual SbBool getCurrentSoftwareFramebuffer(unsigned char *&pixels,
	    unsigned int &width,
	    unsigned int &height,
	    unsigned int &components)
    {
	if (this->fallbackManager)
	    return this->fallbackManager->getCurrentSoftwareFramebuffer(
		pixels, width, height, components);
	pixels = NULL;
	width = height = components = 0;
	return FALSE;
    }

private:
    void *createFallbackContext(QgObolContext *ctx,
	    unsigned int width,
	    unsigned int height)
    {
	if (!this->fallbackManager) {
	    delete ctx;
	    return NULL;
	}

	ctx->fallbackContext = this->fallbackManager->createOffscreenContext(width, height);
	if (!ctx->fallbackContext) {
	    delete ctx;
	    return NULL;
	}

	ctx->kind = QgObolContext::FALLBACK_CONTEXT;
	return ctx;
    }

    bool forceSoftware = false;
    SoDB::ContextManager *fallbackManager = NULL;
};

#endif /* QGOBOLCONTEXTMANAGER_H */
