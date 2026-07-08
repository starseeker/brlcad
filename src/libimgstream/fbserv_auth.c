/*                     F B S E R V _ A U T H . C
 * BRL-CAD
 *
 * Copyright (c) 2025-2026 United States Government as represented by
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
/** @file libimgstream/fbserv_auth.c */

#include "common.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "bu/process.h"
#include "imgstream/fbserv.h"

#ifdef HAVE_SYS_TIME_H
#  include <sys/time.h>
#endif

#ifdef HAVE_OPENSSL_RAND_H
#  include <openssl/rand.h>
#endif

#ifdef HAVE_WINDOWS_H
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#  include <wincrypt.h>
#endif


void
fbserv_generate_token(char *buf)
{
    unsigned char raw[32];
    int i;
    int ok = 0;

    if (!buf)
	return;

    memset(raw, 0, sizeof(raw));

#ifdef HAVE_OPENSSL_RAND_H
    if (!ok && RAND_bytes(raw, (int)sizeof(raw)) == 1)
	ok = 1;
#endif

#ifdef HAVE_WINDOWS_H
    if (!ok) {
	HCRYPTPROV hProv = 0;
	if (CryptAcquireContextW(&hProv, NULL, NULL,
		PROV_RSA_FULL, CRYPT_VERIFYCONTEXT)) {
	    if (CryptGenRandom(hProv, (DWORD)sizeof(raw), raw))
		ok = 1;
	    CryptReleaseContext(hProv, 0);
	}
    }
#else
    if (!ok) {
	FILE *fp = fopen("/dev/urandom", "rb");
	if (fp) {
	    if (fread(raw, 1, sizeof(raw), fp) == sizeof(raw))
		ok = 1;
	    fclose(fp);
	}
    }
#endif

    if (!ok) {
#ifdef HAVE_SYS_TIME_H
	struct timeval tv;
	unsigned int seed;
	gettimeofday(&tv, NULL);
	seed = (unsigned int)tv.tv_sec
	    ^ (unsigned int)tv.tv_usec
	    ^ (unsigned int)bu_pid();
#else
	unsigned int seed = (unsigned int)time(NULL)
	    ^ (unsigned int)bu_pid();
#endif
	srand(seed);
	for (i = 0; i < 32; i++)
	    raw[i] = (unsigned char)(rand() & 0xff);
    }

    for (i = 0; i < 32; i++)
	snprintf(buf + i * 2, 3, "%02x", (unsigned int)raw[i]);
    buf[FBSERV_AUTH_TOKEN_LEN] = '\0';
}


int
fbserv_verify_token(const char *provided, const char *expected)
{
    int diff = 0;
    size_t i;

    if (!provided || !expected)
	return 0;
    if (strlen(provided) != FBSERV_AUTH_TOKEN_LEN
	|| strlen(expected) != FBSERV_AUTH_TOKEN_LEN)
	return 0;

    for (i = 0; i < FBSERV_AUTH_TOKEN_LEN; i++)
	diff |= (unsigned char)provided[i] ^ (unsigned char)expected[i];

    return (diff == 0) ? 1 : 0;
}


/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
