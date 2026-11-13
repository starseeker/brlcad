/*    S E R I A L I Z E D _ B O T _ S O U R C E _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_SERIALIZED_BOT_SOURCE_PRIVATE_H
#define LIBBOBOL_SERIALIZED_BOT_SOURCE_PRIVATE_H

#include <stddef.h>

struct bu_mapped_file;
struct db_i;
struct directory;

/* A borrowed, validated view of the fixed V5 BoT records needed by display
 * discovery and checked cache preparation.  It never owns database bytes. */
struct BObolSerializedBotView {
    const unsigned char *vertices = NULL;
    const unsigned char *faces = NULL;
    size_t vertexCount = 0;
    size_t faceCount = 0;
    size_t vertexStride = 0;
    size_t faceStride = 0;
    unsigned char mode = 0;
    unsigned char orientation = 0;
    unsigned char flags = 0;
};

bool bobol_serialized_bot_view(struct db_i *dbip, struct directory *dp,
	BObolSerializedBotView &view,
	const struct bu_mapped_file *mappedFile = NULL);

#endif /* LIBBOBOL_SERIALIZED_BOT_SOURCE_PRIVATE_H */
