/*          H E A D L E S S _ W I N D O W _ H O S T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/headless_window_host.h */

#ifndef BRLOBOL_HEADLESS_WINDOW_HOST_H
#define BRLOBOL_HEADLESS_WINDOW_HOST_H

#include "brlobol/window_host.h"

#include <Inventor/SbColor.h>
#include <Inventor/SbString.h>

#include <stddef.h>

struct BRLObolHeadlessWindowHostPrivate;

class BRLOBOL_EXPORT BRLObolHeadlessWindowHost : public BRLObolWindowHost {
public:
    BRLObolHeadlessWindowHost(void);
    virtual ~BRLObolHeadlessWindowHost(void);

    virtual int open(const BRLObolWindowDesc *desc = NULL);
    virtual int poll(const BRLObolInputProfile *profile = NULL);
    virtual long pollRate(void) const;

    void setPollRate(long rate);
    int renderPending(void);

    void setBackgroundColor(const SbColor &color);
    const SbColor &getBackgroundColor(void) const;

    void setOutputComponents(int components);
    int getOutputComponents(void) const;

    const unsigned char *getLastFrameBuffer(void) const;
    size_t getLastFrameBufferSize(void) const;
    unsigned int getLastFrameWidth(void) const;
    unsigned int getLastFrameHeight(void) const;
    int getLastFrameComponents(void) const;
    int getRenderCount(void) const;
    const SbString &getLastRenderReason(void) const;

private:
    BRLObolHeadlessWindowHostPrivate *hp;
};

#endif /* BRLOBOL_HEADLESS_WINDOW_HOST_H */
