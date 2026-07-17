/*         B H E A D L E S S W I N D O W H O S T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BHeadlessWindowHost.h */

#ifndef BOBOL_BHEADLESSWINDOWHOST_H
#define BOBOL_BHEADLESSWINDOWHOST_H

#include "BObol/BWindowHost.h"
#include "BObol/BHostFactory.h"

#include <Inventor/SbColor.h>
#include <Inventor/SbString.h>
#include <Inventor/SoDB.h>

#include <stddef.h>

struct BObolHeadlessWindowHostPrivate;

class BOBOL_EXPORT BObolHeadlessWindowHost : public BObolWindowHost {
public:
    BObolHeadlessWindowHost(void);
    virtual ~BObolHeadlessWindowHost(void);

    virtual int open(const BObolWindowDesc *desc = NULL);
    virtual void close(void);
    virtual int poll(const BObolInputProfile *profile = NULL);
    virtual long pollRate(void) const;

    void setPollRate(long rate);
    int renderPending(void);
    void setContextManager(SoDB::ContextManager *manager);
    SoDB::ContextManager *getContextManager(void) const;

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
    BObolHeadlessWindowHostPrivate *hp;
};

/** Return the process-lifetime built-in headless host factory token. */
BOBOL_EXPORT bobol_host_factory_token_t *
bobol_headless_host_factory_register(void);

#endif /* BOBOL_BHEADLESSWINDOWHOST_H */
