/*              Q G O B O L W I N D O W H O S T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file qtcad/QgObolWindowHost.h */

#ifndef QGOBOLWINDOWHOST_H
#define QGOBOLWINDOWHOST_H

#include "common.h"

#include "brlobol/window_host.h"
#include "brlobol/host_factory.h"
#include "qtcad/defines.h"

class QImage;
class QgCanvasBase;
struct QgObolWindowHostPrivate;

class QTCAD_EXPORT QgObolWindowHost : public BRLObolWindowHost {
public:
    explicit QgObolWindowHost(QgCanvasBase *canvas = NULL,
	bool takeCanvasOwnership = false);
    virtual ~QgObolWindowHost(void);

    void setCanvas(QgCanvasBase *canvas);
    QgCanvasBase *canvas(void) const;
    void bindController(BRLObolViewController *controller);

    virtual int open(const BRLObolWindowDesc *desc = NULL);
    virtual void close(void);
    virtual int poll(const BRLObolInputProfile *profile = NULL);
    virtual long pollRate(void) const;

    void setPollRate(long rate);

    const QImage &lastFrameImage(void) const;
    int lastFrameWidth(void) const;
    int lastFrameHeight(void) const;
    int renderCount(void) const;
    const SbString &lastRenderReason(void) const;

private:
    QgObolWindowHostPrivate *qp;
};

/**
 * Register Qt software and, when available, system-GL host factories.
 *
 * For BRLOBOL_HOST_MODE_EMBEDDED, brlobol_host_desc::application_context
 * must point to the QgCanvasBase the factory will borrow.  The canvas type
 * must match the selected qt-sw or qt-gl factory.
 */
QTCAD_EXPORT int qtcad_obol_host_factories_register(void);

#endif /* QGOBOLWINDOWHOST_H */
