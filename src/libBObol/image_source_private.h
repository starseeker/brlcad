/* BRL-CAD
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
#ifndef LIBBOBOL_IMAGE_SOURCE_PRIVATE_H
#define LIBBOBOL_IMAGE_SOURCE_PRIVATE_H

#include "BObol/BImageSource.h"

#include <exception>
#include <memory>

struct BObolImageSourceSubscription;
struct bobol_image_payload;

enum class BObolImageSourceRefresh {
    Failed,
    Current,
    Required
};

/* A source publication may be composed with a display publication.  All
 * potentially failing successor work happens before commit; notifications
 * are explicitly drained after every participant is live. */
class BObolImageSourcePublication {
public:
    explicit BObolImageSourcePublication(SoBRLImageSource &target);
    ~BObolImageSourcePublication();

    BObolImageSourcePublication(const BObolImageSourcePublication &) = delete;
    BObolImageSourcePublication &operator=(const BObolImageSourcePublication &) = delete;

    /* Adopt a stream for a retained graph whose producer owner may disappear
     * before the graph itself. On failure the stream is released here. */
    static int adoptStream(SoBRLImageSource &target, imgstream_t *stream);
    static BObolImageSourceRefresh refreshRequired(SoBRLImageSource &target);

    SoBRLImageSource &next();
    const SoBRLImageSource &successor() const;
    void replaceSubscription(
	std::unique_ptr<BObolImageSourceSubscription> replacement,
	uint64_t realizedGeneration);
    void setRealizedGeneration(uint64_t generation);
    BObolImageSourceRefresh prepareRefresh();
    int loadPreparedPayload(struct bobol_image_payload &payload) const;
    void prepare();
    void commit();
    void restore();
    void notify(std::exception_ptr &failure);
    void notify();
    void publish();

private:
    static BObolImageSourceRefresh queryRefresh(SoBRLImageSource &target,
	struct imgstream_info &info);

    struct Impl;
    std::unique_ptr<Impl> impl;
};

#endif /* LIBBOBOL_IMAGE_SOURCE_PRIVATE_H */
