/*                        M A I N . C X X
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
/** @file main.cxx
 *
 * Command line parsing and main application launching for qbrlcad
 *
 */

#include "common.h"

#include <cstdlib>
#include <iostream>

#include <QApplication>
#include <QImage>
#include <QOpenGLWidget>
#include <QPixmap>
#include <QSurfaceFormat>
#include <QTimer>

#include "bu/app.h"
#include "bu/log.h"
#include "bu/opt.h"
#include "brlcad_version.h"

#include "QgEdApp.h"
#include "fbserv.h"
#include "qtcad/QgTestEvent.h"

static void
qged_usage(const char *cmd, struct bu_opt_desc *d) {
    struct bu_vls str = BU_VLS_INIT_ZERO;
    char *option_help = bu_opt_describe(d, NULL);
    bu_vls_sprintf(&str, "Usage: %s [options] [file.g]\n", cmd);
    if (option_help) {
	bu_vls_printf(&str, "Options:\n%s\n", option_help);
    }
    bu_free(option_help, "help str");
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

    int console_mode = 0;
    int swrast_mode = 0;
    int quad_mode = 0;
    int print_help = 0;
    struct bu_vls msg = BU_VLS_INIT_ZERO;
    const char *startup_commands = NULL;
    const char *test_script = NULL;
    const char *exec_name = argv[0];

    // All BRL-CAD programs need to set this in order for relative path lookups
    // to work reliably.
    bu_setprogname(argv[0]);

    /* Done with command name argv[0] */
    argc-=(argc>0); argv+=(argc>0);

    /* Handle top level application options */
    struct bu_opt_desc d[8];
    BU_OPT(d[0],  "h", "help",   "", NULL, &print_help,    "Print help and exit");
    BU_OPT(d[1],  "?", "",       "", NULL, &print_help,    "");
    BU_OPT(d[2],  "c", "no-gui", "", NULL, &console_mode,  "Run without GUI");
    BU_OPT(d[3],  "s", "swrast", "", NULL, &swrast_mode,   "Use offscreen rendering for 3D view");
    BU_OPT(d[4],  "4", "quad",   "", NULL, &quad_mode,     "Launch using quad view");
    BU_OPT(d[5],  "e", "exec",   "commands", &bu_opt_str, &startup_commands,
	    "Run semicolon-separated GED commands after GUI initialization");
    BU_OPT(d[6],  "", "test-script", "events.json", &bu_opt_str, &test_script,
	    "Replay a qtcad GUI event stream and exit with test status");
    BU_OPT_NULL(d[7]);

    // High level options are only defined prior to the file argument (if there
    // is one).  See if we need to limit our processing
    int acmax = 0;
    for (int i = 0; i < argc; i++) {
	if (argv[i][0] == '-') {
	    acmax++;
	    if ((BU_STR_EQUAL(argv[i], "-e") ||
		    BU_STR_EQUAL(argv[i], "--exec") ||
		    BU_STR_EQUAL(argv[i], "--test-script")) && i + 1 < argc) {
		acmax++;
		i++;
	    }
	} else {
	    break;
	}
    }

    // Process high level args
    int opt_ac = bu_opt_parse(&msg, acmax, (const char **)argv, d);
    if (opt_ac < 0) {
	bu_log("%s\n", bu_vls_cstr(&msg));
	bu_vls_free(&msg);
	return BRLCAD_ERROR;
    }

    // Shift everything not processed by bu_opt_parse down in argv to the last
    // option left behind by bu_opt_parse (or the beginning of the array if
    // opt_ac == 0
    int opt_rem = opt_ac;
    for (int i = acmax; i < argc; i++) {
	argv[opt_ac + (i - acmax)] = argv[i];
	opt_rem++;
    }

    // Set argc to the full count of whatever is left
    argc = opt_rem;

    if (print_help) {
	qged_usage(exec_name, d);
	bu_vls_free(&msg);
	return BRLCAD_OK;
    }

    if (argc > 1 && !console_mode) {
	bu_log("For qged GUI mode need either zero or one .g files specified\n");
	return BRLCAD_ERROR;
    }

    if (console_mode) {
	bu_log("Unimplemented\n");
	return BRLCAD_ERROR;
    }

    if (startup_commands)
	bu_log("qged: queued startup commands: %s\n", startup_commands);

    // Qt6 requires QSurfaceFormat::setDefaultFormat() to be called BEFORE
    // QApplication is constructed.  Without specifying QSurfaceFormat::OpenGL
    // here, Qt's RHI layer may fall back to OpenGL ES (QRhiGles2) even on a
    // desktop system, causing context-creation failures at runtime.
#ifdef BRLCAD_OPENGL
    if (!swrast_mode) {
	QSurfaceFormat fmt;
	fmt.setRenderableType(QSurfaceFormat::OpenGL);
	fmt.setDepthBufferSize(24);
	fmt.setStencilBufferSize(8);
	/* QOpenGLWidget is composited, but the platform presentation contracts are
	 * not equivalent.  DWM already paces the composed desktop and historically
	 * incurs an extra multi-vblank wait here, while unpaced GLX presentation
	 * can expose partially updated/compositor-transition frames during zoom.
	 * Keep Windows latency-oriented and Linux/Unix presentation-safe by
	 * default.  The narrow override makes 0/1/adaptive experiments repeatable
	 * without baking a driver-specific result into global policy. */
#ifdef HAVE_WINDOWS_H
	int swapInterval = 0;
#else
	int swapInterval = 1;
#endif
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
	fmt.setSwapInterval(swapInterval);
	QSurfaceFormat::setDefaultFormat(fmt);
    }
#endif

    const QString startup = startup_commands ?
	QString::fromLocal8Bit(startup_commands) : QString();

    // QApplication requires argv[0] to be the executable name.  The qged
    // option parser above has removed that entry and left only the optional
    // database argument, so give Qt a valid argument vector and pass the
    // database path separately.
    const char *db_file = argc ? argv[0] : NULL;
    int qt_argc = 1;
    char *qt_argv[] = {const_cast<char *>(exec_name), NULL};
    QgEdApp app(qt_argc, qt_argv, db_file, swrast_mode, quad_mode);
    if (!startup.isEmpty()) {
	const QStringList commands = startup.split(';', Qt::SkipEmptyParts);
	for (const QString &command : commands) {
	    const QString trimmed = command.trimmed();
	    bu_log("qged: executing startup command: %s\n",
		    trimmed.toLocal8Bit().constData());
	    app.run_qcmd(trimmed);
	}
    }
    if (test_script) {
	const QString script = QString::fromLocal8Bit(test_script);
	QTimer::singleShot(0, &app, [&app, script]() {
	    QVector<QgTestEvent> events;
	    QString error;
	    QgEventPlayer player(app.w);
	    player.setCheckpointHandler([](QWidget *widget,
		    const QString &name, QString *checkpointError) {
		if (!widget || name.isEmpty()) {
		    if (checkpointError)
			*checkpointError =
			    QStringLiteral("checkpoint needs a widget and output path");
		    return false;
		}
		bool saved = false;
		if (QOpenGLWidget *glWidget =
			qobject_cast<QOpenGLWidget *>(widget)) {
		    const QImage frame = glWidget->grabFramebuffer();
		    saved = !frame.isNull() && frame.save(name);
		} else {
		    const QPixmap frame = widget->grab();
		    saved = !frame.isNull() && frame.save(name);
		}
		if (!saved) {
		    if (checkpointError)
			*checkpointError =
			    QStringLiteral("unable to save checkpoint image: %1")
			    .arg(name);
		    return false;
		}
		return true;
	    });
	    if (!QgEventRecorder::load(script, events, &error) ||
		!player.play(events, &error)) {
		bu_log("qged: GUI test replay failed: %s\n",
		    error.toLocal8Bit().constData());
		app.exit(BRLCAD_ERROR);
		return;
	    }
	    app.exit(BRLCAD_OK);
	});
    }
    bu_vls_free(&msg);

    // Setup complete - time to enter the interactive event loop
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
