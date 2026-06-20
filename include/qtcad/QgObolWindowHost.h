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
#include "qtcad/defines.h"

class QImage;
class QgCanvasBase;
struct QgObolWindowHostPrivate;

class QTCAD_EXPORT QgObolWindowHost : public BRLObolWindowHost {
public:
    explicit QgObolWindowHost(QgCanvasBase *canvas = NULL);
    virtual ~QgObolWindowHost(void);

    void setCanvas(QgCanvasBase *canvas);
    QgCanvasBase *canvas(void) const;

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

#endif /* QGOBOLWINDOWHOST_H */
