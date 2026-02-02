/*                        P L U G I N . H
 * BRL-CAD
 *
 * Copyright (c) 2025 United States Government as represented by
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
/** @file plugin.h
 *
 */

#ifndef LIBGED_PLUGIN_H
#define LIBGED_PLUGIN_H

#include "../ged_private.h"
#include "brlcad_version.h"

#define BU_PLUGIN_NAME ged
#define BU_PLUGIN_CMD_RET int
#define BU_PLUGIN_CMD_ARGS struct ged *, int, const char **
#include "../../libbu/bu_plugin.h"

#define GED_API (2*1000000 + (BRLCAD_VERSION_MAJOR*10000) + (BRLCAD_VERSION_MINOR*100) + BRLCAD_VERSION_PATCH)

#define GED_CMD_DEFAULT       0
#define GED_CMD_INTERACTIVE   (1ULL << 0)
#define GED_CMD_UPDATE_SCENE  (1ULL << 1)
#define GED_CMD_UPDATE_VIEW   (1ULL << 2)
#define GED_CMD_AUTOVIEW      (1ULL << 3)
#define GED_CMD_ALL_VIEWS     (1ULL << 4)
#define GED_CMD_VIEW_CALLBACK (1ULL << 5)

struct ged_cmd_impl {
    const char *cname;
    ged_func_ptr cmd;
    unsigned long long opts;
};

struct ged_cmd_process_impl {
    ged_process_ptr func;
};

#ifdef __cplusplus
extern "C" {
#endif

    /* Registry API: ensure all callers have prototypes in scope */
    int ged_register_command(const struct ged_cmd *cmd);
    int ged_command_exists(const char *name);
    size_t ged_registered_count(void);
    void ged_list_command_names(struct bu_vls *out_csv);
    void ged_list_command_array(const char * const **cl, size_t *cnt);
    void ged_ensure_initialized();
    void libged_shutdown(void);
    const char *ged_init_msgs(void);

#ifdef __cplusplus
}
#endif

#if defined(LIBGED_STATIC_CORE) && !defined(GED_PLUGIN_ONLY)

/* We only need two things:
 *  - the canonical struct ged_cmd object (so there is a stable data symbol)
 *  - a stable anchor pointer __ged_cmd_ptr_<cmd> the scanner can reference
 * No constructors, no linker sections, no CRT $XCU.
 */
#if defined(__GNUC__) || defined(__clang__)
#define GED_CMD_USED __attribute__((used))
#else
#define GED_CMD_USED
#endif

#ifdef __cplusplus
#define REGISTER_GED_COMMAND(cmd_symbol)                                        \
    GED_CMD_USED const struct ged_cmd cmd_symbol##_rcmd = { &cmd_symbol##_impl };      \
    extern "C" GED_CMD_USED const struct ged_cmd * const                        \
        __ged_cmd_ptr_##cmd_symbol = &cmd_symbol##_rcmd
#else
#define REGISTER_GED_COMMAND(cmd_symbol)                                        \
    GED_CMD_USED const struct ged_cmd cmd_symbol##_rcmd = { &cmd_symbol##_impl };      \
    GED_CMD_USED const struct ged_cmd * const                                   \
        __ged_cmd_ptr_##cmd_symbol = &cmd_symbol##_rcmd
#endif

/* If we have odd characters in command names (for example the ? command) we
 * can't use REGISTER_GED_COMMAND.  In those cases we need to construct the
 * REGISTER_GED_COMMAND structures manually using legal names.  However, we
 * also still need the scanner to pick up the command - to achieve this, we add
 * a no-op marker to identify commands for the scanner in such cases.
 *
 * See the ? command in the help command's code for a worked example.
 *
 * Argument MUST be the sanitized C identifier (e.g., questionmark, threeptarb)
 * rather than the cmd string.
 */
#define LABEL_GED_COMMAND(cmd_symbol)

#else
#define REGISTER_GED_COMMAND(cmd_symbol) /* static registration disabled */
#endif /* LIBGED_STATIC_CORE && !GED_PLUGIN_ONLY */


/* =======================================================================================
 * Canonical command list helpers (token-based entries)
 *
 * These macros allow a single per-TU list to drive:
 *   - static registration (when LIBGED_STATIC_CORE is enabled)
 *   - plugin manifest export (when GED_PLUGIN is defined)
 *   - scanner extraction of command names and static symbols
 *
 * Expected list signature:
 *   #define GED_SOMETHING_COMMANDS(X, XID) \
 *     X(token, fn, opts) \
 *     XID(symbol, "cmdname", fn, opts)
 *
 * Where:
 *   - token is a C identifier token representing the command name, e.g. bot_condense
 *   - fn is a ged_func_ptr-compatible function pointer
 *   - opts is a GED_CMD_* flags expression
 *   - symbol is an explicit C identifier for odd/aliased commands
 *
 * From X(token,...), we derive:
 *   - command string: "token"
 *   - symbol id: token##_cmd
 *
 * From XID(symbol,"cmd",...), we use:
 *   - command string: "cmd"
 *   - symbol id: symbol
 * ======================================================================================= */
#define GED_STR_IMPL(x) #x
#define GED_STR(x) GED_STR_IMPL(x)

#define GED_CAT2_IMPL(a,b) a##b
#define GED_CAT2(a,b) GED_CAT2_IMPL(a,b)

/* token -> token_cmd (must expand before use in REGISTER_GED_COMMAND) */
#define GED_CMD_SYM(token) GED_CAT2(token,_cmd)

/* Wrapper to force expansion before REGISTER_GED_COMMAND applies ## internally */
#define GED_REGISTER_SYM(sym) REGISTER_GED_COMMAND(sym)

/* Emit ged_cmd_impl + static registration anchor (token case) */
#define GED__DECL_STATIC_X(token, fn, opts)                                           \
    struct ged_cmd_impl GED_CAT2(GED_CMD_SYM(token), _impl) = { GED_STR(token), fn, opts }; \
    GED_REGISTER_SYM(GED_CMD_SYM(token));

/* Emit ged_cmd_impl + static registration anchor (explicit-id case) */
#define GED__DECL_STATIC_XID(sym, cmdstr, fn, opts)                                   \
    struct ged_cmd_impl GED_CAT2(sym, _impl) = { cmdstr, fn, opts };                  \
    GED_REGISTER_SYM(sym);

/* Emit plugin command array entry (token case) */
#define GED__DECL_PCMD_X(token, fn, opts) { GED_STR(token), fn },

/* Emit plugin command array entry (explicit-id case) */
#define GED__DECL_PCMD_XID(sym, cmdstr, fn, opts) { cmdstr, fn },

/* Public: declare a command set in a TU */
#if defined(LIBGED_STATIC_CORE) && !defined(GED_PLUGIN_ONLY)
#  define GED_DECLARE_COMMAND_SET(LIST_MACRO) \
      LIST_MACRO(GED__DECL_STATIC_X, GED__DECL_STATIC_XID)
#else
  /* File-scope safe no-op (no variables, no statements) */
#  define GED_DECLARE_COMMAND_SET(LIST_MACRO) \
      /* static registration disabled */
#endif

/* Public: declare plugin manifest in a TU */
#ifdef GED_PLUGIN
#  define GED_DECLARE_PLUGIN_MANIFEST(plugin_name_str, plugin_version_u32, LIST_MACRO) \
    static bu_plugin_cmd pcommands[] = {                                           \
        LIST_MACRO(GED__DECL_PCMD_X, GED__DECL_PCMD_XID)                           \
    };                                                                             \
    static bu_plugin_manifest pinfo = {                                            \
        plugin_name_str,                                                           \
        (unsigned int)(plugin_version_u32),                                        \
        (unsigned int)(sizeof(pcommands)/sizeof(pcommands[0])),                    \
        pcommands,                                                                 \
        BU_PLUGIN_ABI_VERSION,                                                     \
        sizeof(bu_plugin_manifest)                                                 \
    };                                                                             \
    BU_PLUGIN_DECLARE_MANIFEST(pinfo)
#else
/* File-scope safe no-op (no variables, no statements) */
#  define GED_DECLARE_PLUGIN_MANIFEST(plugin_name_str, plugin_version_u32, LIST_MACRO) \
    /* not building plugin */
#endif

#endif /* LIBGED_PLUGIN_H */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8 cino=N-s
 */
