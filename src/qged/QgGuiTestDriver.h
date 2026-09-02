/*                    Q G G U I T E S T D R I V E R . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */

#ifndef QGGUITESTDRIVER_H
#define QGGUITESTDRIVER_H

#include <QJsonObject>
#include <QString>

#include "qtcad/QgTestEvent.h"

class QgEdApp;
class BObolViewController;
struct BObolLodConvergenceStatus;

BObolViewController *qged_progressive_event_controller(
    QgEdApp &app, const QgTestEvent &event);

/**
 * Append the compact controller fields used by the progressive-control
 * refinement checker.  Keeping this projection in one function prevents the
 * per-event report and the denser wait-loop trace from acquiring different
 * meanings for the same controller state.
 */
void qged_append_progressive_control_diagnostics(
    QJsonObject &sample, const BObolViewController &controller,
    const BObolLodConvergenceStatus &status);

/**
 * Collect one progressive-display diagnostic sample.  Kept separate from the
 * event driver so normal GUI replay can avoid deep scene inspection except at
 * explicit structural checkpoints.
 */
QJsonObject qged_collect_progressive_sample(
    QgEdApp &app, int eventIndex, const QgTestEvent &event,
    qint64 elapsedMilliseconds, qint64 eventMicroseconds,
    bool collectDeepLodDiagnostics);

/**
 * Schedule a GUI event stream after application initialization.  Completion
 * exits @p app with the replay status.
 */
void qged_schedule_gui_test(
    QgEdApp &app, const QString &script, const QString &reportFile,
    bool softwareRenderer);

#endif
