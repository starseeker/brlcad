/*                  I M A G E _ S O U R C E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/image_source.h */

#ifndef BRLOBOL_IMAGE_SOURCE_H
#define BRLOBOL_IMAGE_SOURCE_H

#include "brlobol/defines.h"

#include "imgstream/stream.h"

#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFEnum.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFUInt32.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoSubNode.h>

class BRLOBOL_EXPORT SoBRLImageSource : public SoNode {
    typedef SoNode inherited;

    SO_NODE_HEADER(SoBRLImageSource);

public:
    enum SourceKind {
	SOURCE_EMPTY = 0,
	SOURCE_STATIC_IMAGE = 1,
	SOURCE_IMAGE_STREAM = 2
    };

    enum Status {
	STATUS_EMPTY = 0,
	STATUS_READY = 1,
	STATUS_STREAMING = 2,
	STATUS_FAILED = 3
    };

    enum PixelFormat {
	PIXEL_UNKNOWN = 0,
	PIXEL_RGB8 = IMGSTREAM_PIXEL_RGB8,
	PIXEL_RGBA8 = IMGSTREAM_PIXEL_RGBA8
    };

    SoSFString imageId;
    SoSFString sourceUri;
    SoSFEnum sourceKind;
    SoSFEnum status;
    SoSFString diagnostic;
    SoSFUInt32 pixelWidth;
    SoSFUInt32 pixelHeight;
    SoSFEnum pixelFormat;
    SoSFString colorSpace;
    SoSFString alphaMode;
    SoSFUInt32 sourceRevision;
    SoSFUInt32 dataRevision;
    SoSFUInt32 dirtyRevision;
    SoSFUInt32 dirtyX;
    SoSFUInt32 dirtyY;
    SoSFUInt32 dirtyWidth;
    SoSFUInt32 dirtyHeight;
    SoSFBool streamConnected;
    SoSFBool producerActive;

    SoBRLImageSource(void);
    static void initClass(void);

    int setImage(const icv_image_t *image);
    int setStream(imgstream_t *stream);
    void clearSource(void);
    imgstream_t *getStream(void) const;
    int refreshFromStream(void);
    SbBool ownsStream(void) const;

protected:
    virtual ~SoBRLImageSource(void);

private:
    void releaseStream(void);
    int attachStream(imgstream_t *stream, SbBool owned, SourceKind kind);
    void resetStreamFields(Status statusValue, const char *diagnosticValue);
    static void dirtyCB(void *ctx, const struct imgstream_rect *rect, uint64_t generation);

    imgstream_t *stream;
    imgstream_subscriber_id subscriberId;
    SbBool streamOwned;
};

#endif /* BRLOBOL_IMAGE_SOURCE_H */
