/*                          M A N . C
 * BRL-CAD
 *
 * Copyright (c) 1988-2026 United States Government as represented by
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
/** @addtogroup libged */
/** @{ */
/** @file libged/man.c
 *
 * Define a command for launching the man page viewer.
 */
/** @} */

#include "common.h"

#include <ctype.h> /* for isdigit */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif
#include "bresource.h"

#include "bu/app.h"
#include "bu/cmd.h"
#include "bu/env.h"
#include "bu/file.h"
#include "bu/process.h"

#include "../ged_private.h"


/**
 * Invoke brlman with the current environmental options, or
 * (if supplied) user defined options
 */
int
ged_man_core(struct ged *gedp, int argc, const char *argv[])
{
    int ret = BRLCAD_OK;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    const char *brlman = NULL;
    char brlmancmd[MAXPATHLEN] = {0};
    brlman = bu_dir(brlmancmd, MAXPATHLEN, BU_DIR_BIN, "brlman", BU_DIR_EXT, NULL);
    if (!bu_file_exists(brlman, NULL)) {
	bu_vls_printf(gedp->ged_result_str, "ERROR: Unable to find 'brlman' viewer: %s.\n", brlman);
	return BRLCAD_ERROR;
    }

    int enable_gui = 0;
    int disable_gui = 0;
    int print_help = 0;
    int dump = 0;
    enum { default_text_width = 80 };
    int width = default_text_width;
    char man_section = 'n';
    struct bu_vls lang = BU_VLS_INIT_ZERO;
    struct bu_vls optparse_msg = BU_VLS_INIT_ZERO;
    struct bu_opt_desc d[8];
    BU_OPT(d[0], "h", "help",        "",                NULL, &print_help,  "Print help and exit");
    BU_OPT(d[1], "g", "gui",         "",                NULL, &enable_gui,  "Enable GUI");
    BU_OPT(d[2], "",  "no-gui",      "",                NULL, &disable_gui, "Disable GUI");
    BU_OPT(d[3], "L", "language",  "lg",        &bu_opt_lang, &lang,        "Set language");
    BU_OPT(d[4], "S", "section",    "#", &bu_opt_man_section, &man_section, "Set section");
    BU_OPT(d[5], "", "dump", "", NULL, &dump, "Return formatted text in the command result");
    BU_OPT(d[6], "", "width", "columns", &bu_opt_int, &width, "Set formatted text width (default: 80)");
    BU_OPT_NULL(d[7]);

    // Parse options (if any)
    argv++; argc--;
    int uac = bu_opt_parse(&optparse_msg, argc, argv, d);
    if (uac == -1) {
	bu_vls_printf(gedp->ged_result_str, "%s", bu_vls_cstr(&optparse_msg));
	bu_vls_free(&optparse_msg);
	bu_vls_free(&lang);
	return BRLCAD_ERROR;
    }
    bu_vls_free(&optparse_msg);

    /* If we want help, print help */
    if (print_help) {
	char *option_help = bu_opt_describe(d, NULL);
	bu_vls_printf(gedp->ged_result_str, "Usage: brlman [options] [man_page]\n");
	if (option_help) {
	    bu_vls_printf(gedp->ged_result_str, "Options:\n%s\n", option_help);
	}
	bu_free(option_help, "help str");
	bu_vls_free(&lang);
	return BRLCAD_OK;
    }

    /* If we only have one non-option arg, assume it's the man name */
    const char *man_name = NULL;
    if (uac == 1) {
	man_name = argv[0];
    } else {
	/* If we've got more, check for a section number. */
	for (int i = 0; i < uac; i++) {
	    if (strlen(argv[i]) == 1 && (isdigit(argv[i][0]) || argv[i][0] == 'n')) {

		/* Record the section */
		man_section = argv[i][0];

		/* If we have a section identifier and it's not already the last
		 * element in argv, adjust the array */
		if (i < uac - 1) {
		    const char *tmp = argv[uac-1];
		    argv[uac-1] = argv[i];
		    argv[i] = tmp;
		}

		/* Found a number and processed - one less argv entry. */
		uac--;

		/* First number wins. */
		break;
	    }
	}

	/* For now, only support specifying one man page at a time */
	if (uac > 1) {
	    bu_vls_printf(gedp->ged_result_str, "Error - need a single man page name");
	    bu_vls_free(&lang);
	    return BRLCAD_ERROR;
	}

	/* If we have one, use it.  Zero is allowed if we're in graphical
	 * mode, so uac == at this point isn't necessarily a failure */
	if (uac == 1) {
	    man_name = argv[0];
	}
    }

    // Next, we check environment variables for any options.  If command line values are
    // set they override them, but if not let the environment tell us what we need
    if (dump)
	disable_gui = 1;

    if (!enable_gui && !disable_gui) {
	char *cmd_flags;
	cmd_flags = getenv("GED_MAN_CMDLINE_MODE");
	if (cmd_flags && bu_str_true(cmd_flags)) {
	    disable_gui = 1;
	}
    }
    if (!enable_gui && !disable_gui)
	enable_gui = 1;

    if (!bu_vls_strlen(&lang)) {
	char *lang_flags;
	lang_flags = getenv("GED_MAN_LANG_MODE");
	if (lang_flags) {
	    struct bu_vls langparse_msg = BU_VLS_INIT_ZERO;
	    const char *av[2] = {NULL};
	    av[0] = (const char *)lang_flags;
	    int lac = bu_opt_lang(&langparse_msg, 1, av, &lang);
	    if (lac == -1) {
		bu_vls_printf(gedp->ged_result_str, "getenv(GED_MAN_LANG_MODE) failure: %s", bu_vls_cstr(&langparse_msg));
		bu_vls_free(&lang);
		bu_vls_free(&langparse_msg);
		return BRLCAD_ERROR;
	    }
	    bu_vls_free(&langparse_msg);
	}
    }

    if (width <= 0) {
	bu_vls_printf(gedp->ged_result_str, "Error - width must be positive.\n");
	bu_vls_free(&lang);
	return BRLCAD_ERROR;
    }
    if (disable_gui && !man_name) {
	bu_vls_printf(gedp->ged_result_str, "Error - please specify a man page.\n");
	bu_vls_free(&lang);
	return BRLCAD_ERROR;
    }
    if (enable_gui && disable_gui) {
	bu_vls_printf(gedp->ged_result_str, "Error - both GUI and command line man viewer modes specified.");
	bu_vls_free(&lang);
	return BRLCAD_ERROR;
    }

    const char *brlman_exec_cmd[11] = {NULL};
    struct bu_vls man_sect = BU_VLS_INIT_ZERO;
    struct bu_vls man_width = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&man_sect, "%c", man_section);
    bu_vls_sprintf(&man_width, "%d", width);
    brlman_exec_cmd[0] = brlmancmd;
    brlman_exec_cmd[1] = "-L";
    brlman_exec_cmd[2] = bu_vls_strlen(&lang) ? bu_vls_cstr(&lang) : "en";
    brlman_exec_cmd[3] = "-S";
    brlman_exec_cmd[4] = bu_vls_cstr(&man_sect);
    brlman_exec_cmd[5] = disable_gui ? "--dump" : "--gui";
    int next_arg = 6;
    if (disable_gui) {
	brlman_exec_cmd[next_arg++] = "--width";
	brlman_exec_cmd[next_arg++] = bu_vls_cstr(&man_width);
    }
    if (man_name) {
	brlman_exec_cmd[next_arg++] = "--";
	brlman_exec_cmd[next_arg++] = man_name;
    }

    struct bu_process *process = NULL;
    int process_options = BU_PROCESS_HIDE_WINDOW;
    if (disable_gui)
	process_options |= BU_PROCESS_OUT_EQ_ERR;
    bu_process_create(&process, brlman_exec_cmd, process_options);
    if (!process || bu_process_pid(process) == -1) {
	bu_vls_printf(gedp->ged_result_str, "Unable to launch brlman viewer.\n");
	if (process)
	    bu_process_wait_n(&process, 0);
	ret = BRLCAD_ERROR;
    } else if (disable_gui) {
	char buffer[BUFSIZ];
	int count;
	while ((count = bu_process_read_n(process, BU_PROCESS_STDOUT,
		(int)sizeof(buffer), buffer)) > 0)
	    bu_vls_strncat(gedp->ged_result_str, buffer, (size_t)count);
	if (count < 0) {
	    bu_vls_printf(gedp->ged_result_str, "Unable to read brlman output.\n");
	    bu_process_terminate(process);
	}
	int exit_status = bu_process_wait_n(&process, 0);
	if (count < 0 || exit_status != 0)
	    ret = BRLCAD_ERROR;
    }

    bu_vls_free(&lang);
    bu_vls_free(&man_sect);
    bu_vls_free(&man_width);

    return ret;
}


#include "../include/plugin.h"

#define GED_MAN_COMMANDS(X, XID) \
    X(man, ged_man_core, GED_CMD_DEFAULT) \

GED_DECLARE_COMMAND_SET(GED_MAN_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_man", 1, GED_MAN_COMMANDS)

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
