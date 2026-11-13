/*    S E R I A L I Z E D _ B O T _ S O U R C E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "serialized_bot_source_private.h"

#include "bu/mapped_file.h"
#include "raytrace.h"
#include "rt/db5.h"

bool
bobol_serialized_bot_view(struct db_i *dbip, struct directory *dp,
	BObolSerializedBotView &view, const struct bu_mapped_file *mappedFile)
{
    view = BObolSerializedBotView();
    if (!dbip || !dp || db_version(dbip) != 5 ||
	dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BOT)
	return false;

    size_t serializedBytes = 0;
    const unsigned char *serialized = db_external_view(dbip, dp,
	&serializedBytes);
    if (!serialized && mappedFile && mappedFile->buf && dp->d_addr >= 0 &&
	static_cast<size_t>(dp->d_addr) <= mappedFile->buflen &&
	dp->d_len <= mappedFile->buflen - static_cast<size_t>(dp->d_addr)) {
	serialized = static_cast<const unsigned char *>(mappedFile->buf) +
	    static_cast<size_t>(dp->d_addr);
	serializedBytes = dp->d_len;
    }
    const size_t headerBytes = sizeof(struct db5_ondisk_header);
    if (!serialized || serializedBytes < headerBytes)
	return false;
    const int width = (serialized[1] & DB5HDR_HFLAGS_OBJECT_WIDTH_MASK) >>
	DB5HDR_HFLAGS_OBJECT_WIDTH_SHIFT;
    if (width < DB5HDR_WIDTHCODE_8BIT || width > DB5HDR_WIDTHCODE_64BIT)
	return false;
    const size_t widthBytes = static_cast<size_t>(1) << width;
    if (widthBytes > serializedBytes - headerBytes)
	return false;
    size_t encodedObjectUnits = 0;
    (void)db5_decode_length(&encodedObjectUnits, serialized + headerBytes,
	width);
    if (encodedObjectUnits > SIZE_MAX / 8)
	return false;
    const size_t objectLengthBytes = encodedObjectUnits * 8;
    if (objectLengthBytes < headerBytes || objectLengthBytes > serializedBytes)
	return false;

    struct db5_raw_internal raw;
    if (!db5_get_raw_internal_ptr(&raw, serialized) ||
	raw.object_length != objectLengthBytes ||
	raw.major_type != DB5_MAJORTYPE_BRLCAD ||
	raw.minor_type != DB5_MINORTYPE_BRLCAD_BOT || !raw.body.ext_buf ||
	raw.body.ext_nbytes < 2 * SIZEOF_NETWORK_LONG + 3 ||
	raw.body.ext_buf < serialized)
	return false;
    const size_t bodyOffset =
	static_cast<size_t>(raw.body.ext_buf - serialized);
    if (bodyOffset >= raw.object_length ||
	raw.body.ext_nbytes > raw.object_length - bodyOffset - 1)
	return false;

    const unsigned char *body = raw.body.ext_buf;
    const size_t vertices = static_cast<size_t>(BU_GLONG(body));
    const size_t faces = static_cast<size_t>(
	BU_GLONG(body + SIZEOF_NETWORK_LONG));
    const size_t fixedBytes = 2 * SIZEOF_NETWORK_LONG + 3;
    view.vertexStride = SIZEOF_NETWORK_DOUBLE * ELEMENTS_PER_POINT;
    view.faceStride = 3 * SIZEOF_NETWORK_LONG;
    if (!vertices || !faces ||
	vertices > (SIZE_MAX - fixedBytes) / view.vertexStride)
	return false;
    const size_t vertexEnd = fixedBytes + vertices * view.vertexStride;
    if (vertexEnd > raw.body.ext_nbytes ||
	faces > (raw.body.ext_nbytes - vertexEnd) / view.faceStride)
	return false;

    view.vertices = body + fixedBytes;
    view.faces = body + vertexEnd;
    view.vertexCount = vertices;
    view.faceCount = faces;
    view.orientation = body[2 * SIZEOF_NETWORK_LONG];
    view.mode = body[2 * SIZEOF_NETWORK_LONG + 1];
    view.flags = body[2 * SIZEOF_NETWORK_LONG + 2];
    return true;
}
