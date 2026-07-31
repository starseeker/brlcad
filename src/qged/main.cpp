/*                        M A I N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2014-2026 United States Government as represented by
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
/** @file main.cpp
 *
 * Command-line parsing and main application launching for qged.
 */

#include "common.h"

#include <cstdlib>

#include <QSurfaceFormat>

#include "brlcad_version.h"
#include "bu/app.h"
#include "bu/log.h"
#include "bu/opt.h"
#include "bu/vls.h"
#include "QgEdApp.h"
#include "QgGuiTestDriver.h"

static void
qged_usage(const char *cmd, struct bu_opt_desc *descriptions)
{
    struct bu_vls str = BU_VLS_INIT_ZERO;
    char *optionHelp = bu_opt_describe(descriptions, NULL);
    bu_vls_sprintf(&str, "Usage: %s [options] [file.g]\n", cmd);
    if (optionHelp)
	bu_vls_printf(&str, "Options:\n%s\n", optionHelp);
    bu_free(optionHelp, "help str");
    bu_log("%s", bu_vls_cstr(&str));
    bu_vls_free(&str);
}

#ifdef HAVE_WINDOWS_H
int APIENTRY
WinMain(HINSTANCE hInstance,
	HINSTANCE hPrevInstance,
	LPSTR lpszCmdLine,
	int nCmdShow)
{
    int argc = __argc;
    char **argv = __argv;
#else
int
main(int argc, char **argv)
{
#endif
    int consoleMode = 0;
    int softwareRasterMode = 0;
    int quadMode = 0;
    int printHelp = 0;
    struct bu_vls message = BU_VLS_INIT_ZERO;
    const char *startupCommands = NULL;
    const char *testScript = NULL;
    const char *testReport = NULL;
    const char *executableName = argv[0];

    bu_setprogname(argv[0]);

    if (argc > 0) {
	argc--;
	argv++;
    }

    struct bu_opt_desc descriptions[9];
    BU_OPT(descriptions[0], "h", "help", "", NULL, &printHelp,
	"Print help and exit");
    BU_OPT(descriptions[1], "?", "", "", NULL, &printHelp, "");
    BU_OPT(descriptions[2], "c", "no-gui", "", NULL, &consoleMode,
	"Run without GUI");
    BU_OPT(descriptions[3], "s", "swrast", "", NULL, &softwareRasterMode,
	"Use offscreen rendering for 3D view");
    BU_OPT(descriptions[4], "4", "quad", "", NULL, &quadMode,
	"Launch using quad view");
    BU_OPT(descriptions[5], "e", "exec", "commands", &bu_opt_str,
	&startupCommands,
	"Run semicolon-separated GED commands after GUI initialization");
    BU_OPT(descriptions[6], "", "test-script", "events.json", &bu_opt_str,
	&testScript, "Replay a qtcad GUI event stream and exit with test status");
    BU_OPT(descriptions[7], "", "test-report", "report.json", &bu_opt_str,
	&testReport, "Write timing and LoD samples while replaying a GUI test");
    BU_OPT_NULL(descriptions[8]);

    int optionArgumentCount = 0;
    for (int i = 0; i < argc; i++) {
	if (argv[i][0] != '-')
	    break;
	optionArgumentCount++;
	if ((BU_STR_EQUAL(argv[i], "-e") ||
		BU_STR_EQUAL(argv[i], "--exec") ||
		BU_STR_EQUAL(argv[i], "--test-script") ||
		BU_STR_EQUAL(argv[i], "--test-report")) &&
	    i + 1 < argc) {
	    optionArgumentCount++;
	    i++;
	}
    }

    int parsedCount = bu_opt_parse(
	&message, optionArgumentCount, (const char **)argv, descriptions);
    if (parsedCount < 0) {
	bu_log("%s\n", bu_vls_cstr(&message));
	bu_vls_free(&message);
	return BRLCAD_ERROR;
    }

    int remainingCount = parsedCount;
    for (int i = optionArgumentCount; i < argc; i++) {
	argv[parsedCount + (i - optionArgumentCount)] = argv[i];
	remainingCount++;
    }
    argc = remainingCount;

    if (printHelp) {
	qged_usage(executableName, descriptions);
	bu_vls_free(&message);
	return BRLCAD_OK;
    }

    if (argc > 1 && !consoleMode) {
	bu_log("For qged GUI mode need either zero or one .g files specified\n");
	bu_vls_free(&message);
	return BRLCAD_ERROR;
    }

    if (consoleMode) {
	bu_log("Unimplemented\n");
	bu_vls_free(&message);
	return BRLCAD_ERROR;
    }

    if (startupCommands)
	bu_log("qged: queued startup commands: %s\n", startupCommands);

#ifdef BRLCAD_OPENGL
    if (!softwareRasterMode) {
	QSurfaceFormat format;
	format.setRenderableType(QSurfaceFormat::OpenGL);
	format.setDepthBufferSize(24);
	format.setStencilBufferSize(8);

	/* QOpenGLWidget renders into a Qt-owned FBO.  The compositor owns final
	 * presentation pacing; applying another v-sync gate to the intermediate
	 * context severely reduces interactive retained-view throughput. */
	int swapInterval = 0;
	const char *swapOverride = std::getenv("QGED_SWAP_INTERVAL");
	if (swapOverride && swapOverride[0]) {
	    char *end = NULL;
	    const long value = std::strtol(swapOverride, &end, 10);
	    if (end && end[0] == '\0' && value >= -1 && value <= 1)
		swapInterval = static_cast<int>(value);
	    else
		bu_log("qged: ignoring invalid QGED_SWAP_INTERVAL=%s "
		       "(expected -1, 0, or 1)\n", swapOverride);
	}
	format.setSwapInterval(swapInterval);
	QSurfaceFormat::setDefaultFormat(format);
    }
#endif

    const QString startup = startupCommands ?
	QString::fromLocal8Bit(startupCommands) : QString();

    /* QApplication requires a valid executable argv[0].  The qged option
     * parser has removed it from the optional database arguments. */
    const char *databaseFile = argc ? argv[0] : NULL;
    int qtArgc = 1;
    char *qtArgv[] = {const_cast<char *>(executableName), NULL};
    QgEdApp app(
	qtArgc, qtArgv, databaseFile, softwareRasterMode, quadMode);

    if (!startup.isEmpty()) {
	const QStringList commands = startup.split(';', Qt::SkipEmptyParts);
	for (const QString &command : commands) {
	    const QString trimmed = command.trimmed();
	    bu_log("qged: executing startup command: %s\n",
		trimmed.toLocal8Bit().constData());
	    app.run_qcmd(trimmed);
	}
    }

    if (testScript) {
	const QString script = QString::fromLocal8Bit(testScript);
	const QString reportFile = testReport ?
	    QString::fromLocal8Bit(testReport) : QString();
	qged_schedule_gui_test(
	    app, script, reportFile, softwareRasterMode != 0);
    }
    bu_vls_free(&message);
    return app.exec();
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
