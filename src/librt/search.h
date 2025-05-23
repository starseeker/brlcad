/*                        S E A R C H . H
 * BRL-CAD
 *
 * Copyright (c) 2008-2025 United States Government as represented by
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

/* Includes code from OpenBSD's find command:
 *
 * Copyright (c) 1990, 1993
 *	The Regents of the University of California.  All rights reserved.
 *
 * This code is derived from software contributed to Berkeley by
 * Cimarron D. Taylor of the University of California, Berkeley.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */

#ifndef LIBRT_SEARCH_H
#define LIBRT_SEARCH_H

#include "common.h"

#include <sys/types.h> /* for gid_t */

#include "bu/ptbl.h"
#include "raytrace.h"

// Execution context for the -exec filter.
struct _db_search_ctx {
    bu_clbk_t clbk; /**< @brief A function that evaluates an array of strings and returns a boolean. */
    void *u1; /**< @brief A pointer that will be passed to the callback. */
    void *u2; /**< @brief A pointer that will be passed to the callback. */
};

/* node struct - holds data specific to each node under consideration */
struct db_node_t {
    struct db_full_path *path;
    struct bu_ptbl *full_paths;
    int flags;
    int matched_filters;
};

/* search node type */
enum db_search_ntype {
    N_ABOVE = 1,                        /* must start > 0 */
    N_AND, N_ATTR, N_BELOW,
    N_CLOSEPAREN, N_DEPTH, N_EXEC, N_EXECDIR, N_EXPR,
    N_FLAGS, N_INAME, N_IREGEX, N_LS, N_MATRIX, N_MAXDEPTH,
    N_MINDEPTH, N_NAME, N_NNODES, N_NOT, N_OK, N_OPENPAREN, N_OR, N_PATH,
    N_PRINT, N_PRUNE, N_REGEX, N_SIZE, N_STDATTR, N_TYPE, N_BOOL, N_PARAM
};


struct db_plan_t {
    struct db_plan_t *next;			/* next node */
    int (*eval)(struct db_plan_t *, struct db_node_t *, struct db_i *dbip, struct bu_ptbl *results);
    /* node evaluation function */
#define F_EQUAL 1 /* [acm]time inum links size */
#define F_LESSTHAN 2
#define F_GREATER 3
#define F_NEEDOK 1 /* exec ok */
#define F_MTFLAG 1 /* fstype */
#define F_MTTYPE 2
#define F_ATLEAST 1 /* perm */
    int min_depth;
    int max_depth;
    mat_t m;
    int flags;				/* private flags */
    enum db_search_ntype type;		/* plan node type */
    struct bu_ptbl *plans;              /* set of all allocated plans */
    union {
	gid_t _g_data;			/* gid */
	struct {
	    unsigned int _f_flags;
	    unsigned int _f_mask;
	} fl;
	struct db_plan_t *_p_data[2];	/* PLAN trees */
	struct _ex {
	    bu_clbk_t _e_callback;
	    void *_e_userdata1;
	    void *_e_userdata2;
	    char **_e_argv;		/* argv array */
	    int _e_argc;		/* allocated length */
	    int *_e_holes;
	    int _e_nholes;
	} ex;
	struct db_plan_t *_ab_data[2];	/* PLAN trees */
	struct db_plan_t *_bl_data[2];  /* PLAN trees */
	char *_a_data[2];		/* array of char pointers */
	char *_c_data;			/* char pointer */
	char *_ci_data;			/* char pointer */
	char *_path_data;		/* char pointer */
	char *_attr_data;		/* char pointer */
	char *_param_data;		/* char pointer */
	char *_depth_data;		/* char pointer */
	char *_node_data;		/* char pointer */
	char *_type_data;
	char *_regex_pattern;		/* regexp pattern */
	int _bool_data;
	int _max_data;			/* tree depth */
	int _min_data;			/* tree depth */
    } p_un;
};


typedef struct _option {
    const char *name;				/* option name */
    enum db_search_ntype token;			/* token type */
    int (*create)(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);	/* create function */
#define O_NONE		0x01			/* no call required */
#define O_ZERO		0x02			/* pass: nothing */
#define O_ARGV		0x04			/* pass: argv, increment argv */
#define O_ARGVP		0x08		/* pass: *argv, N_OK || N_EXEC || N_EXECDIR */
    int flags;
} OPTION;

__BEGIN_DECLS

static int c_attr(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_objparam(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_iname(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_maxdepth(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_mindepth(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_depth(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_name(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_nnodes(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_regex(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_iregex(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_path(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_print(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_stdattr(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_matrix(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_type(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_bool(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_openparen(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_closeparen(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_not(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_or(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_above(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_below(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_exec(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);
static int c_size(char *, char ***, int, struct db_plan_t **, int *, struct bu_ptbl *, struct _db_search_ctx *);

__END_DECLS

#endif /* LIBRT_SEARCH_H */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
