/*                   F B S E R V _ A U T H . C
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
/** @file libimgstream/tests/fbserv_auth.c */

#include "common.h"

#include <stdio.h>
#include <string.h>

#include "bu/app.h"
#include "imgstream/fbserv.h"


static int
token_is_hex(const char *token)
{
    if (!token || strlen(token) != FBSERV_AUTH_TOKEN_LEN)
	return 0;

    for (size_t i = 0; i < FBSERV_AUTH_TOKEN_LEN; i++) {
	if (!((token[i] >= '0' && token[i] <= '9') ||
	    (token[i] >= 'a' && token[i] <= 'f')))
	    return 0;
    }

    return 1;
}


int
main(int argc, char **argv)
{
    char token[FBSERV_AUTH_TOKEN_LEN + 1] = {0};
    char mismatch[FBSERV_AUTH_TOKEN_LEN + 1] = {0};

    bu_setprogname(argv[0]);
    (void)argc;

    fbserv_generate_token(token);
    if (!token_is_hex(token)) {
	fprintf(stderr, "generated token is not lowercase hexadecimal\n");
	return 1;
    }
    if (!fbserv_verify_token(token, token)) {
	fprintf(stderr, "generated token did not verify\n");
	return 1;
    }

    memcpy(mismatch, token, sizeof(mismatch));
    mismatch[0] = mismatch[0] == '0' ? '1' : '0';
    if (fbserv_verify_token(mismatch, token) ||
	fbserv_verify_token("short", token) ||
	fbserv_verify_token(NULL, token) ||
	fbserv_verify_token(token, NULL)) {
	fprintf(stderr, "invalid token verified\n");
	return 1;
    }

    return 0;
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
