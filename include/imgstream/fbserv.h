/*                        F B S E R V . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file imgstream/fbserv.h */
#ifndef IMGSTREAM_FBSERV_H
#define IMGSTREAM_FBSERV_H

#include "common.h"

__BEGIN_DECLS

#ifndef IMGSTREAM_EXPORT
#  if defined(IMGSTREAM_DLL_EXPORTS) && defined(IMGSTREAM_DLL_IMPORTS)
#    error "Only IMGSTREAM_DLL_EXPORTS or IMGSTREAM_DLL_IMPORTS can be defined, not both."
#  elif defined(IMGSTREAM_DLL_EXPORTS)
#    define IMGSTREAM_EXPORT COMPILER_DLLEXPORT
#  elif defined(IMGSTREAM_DLL_IMPORTS)
#    define IMGSTREAM_EXPORT COMPILER_DLLIMPORT
#  else
#    define IMGSTREAM_EXPORT
#  endif
#endif

#define FBSERV_MSG_FBOPEN        1
#define FBSERV_MSG_FBCLOSE       2
#define FBSERV_MSG_FBCLEAR       3
#define FBSERV_MSG_FBREAD        4
#define FBSERV_MSG_FBWRITE       5
#define FBSERV_MSG_FBCURSOR      6
#define FBSERV_MSG_FBWINDOW      7
#define FBSERV_MSG_FBZOOM        8
#define FBSERV_MSG_FBSCURSOR     9
#define FBSERV_MSG_FBVIEW        10
#define FBSERV_MSG_FBGETVIEW     11
#define FBSERV_MSG_FBRMAP        12
#define FBSERV_MSG_FBWMAP        13
#define FBSERV_MSG_FBHELP        14
#define FBSERV_MSG_FBREADRECT    15
#define FBSERV_MSG_FBWRITERECT   16
#define FBSERV_MSG_FBFLUSH       17
#define FBSERV_MSG_FBFREE        18
#define FBSERV_MSG_FBGETCURSOR   19
#define FBSERV_MSG_DATA          20
#define FBSERV_MSG_RETURN        21
#define FBSERV_MSG_CLOSE         22
#define FBSERV_MSG_ERROR         23
#define FBSERV_MSG_FBPOLL        30
#define FBSERV_MSG_FBSETCURSOR   31
#define FBSERV_MSG_FBBWREADRECT  32
#define FBSERV_MSG_FBBWWRITERECT 33
#define FBSERV_MSG_FBAUTH        34
#define FBSERV_MSG_FBTLS         35
#define FBSERV_MSG_NORETURN      100

#define FBSERV_NET_LONG_LEN 4
#define FBSERV_AUTH_TOKEN_LEN 64
#define FBSERV_AUTH_TOKEN_ENVVAR "FBSERV_AUTH_TOKEN"
#define FBSERV_TLS_SERVER_SHA256_ENVVAR "FBSERV_TLS_SHA256"
#define FBSERV_INSECURE_ENVVAR "FBSERV_INSECURE"

IMGSTREAM_EXPORT void fbserv_generate_token(char *buffer);
IMGSTREAM_EXPORT int fbserv_verify_token(const char *provided,
	const char *expected);

__END_DECLS

#endif /* IMGSTREAM_FBSERV_H */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
