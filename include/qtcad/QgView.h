/*                     Q G V I E W . H
 * BRL-CAD
 *
 * Copyright (c) 2021-2026 United States Government as represented by
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
/** @file QgView.h
 *
 */

#ifndef QGVIEW_H
#define QGVIEW_H

#include "common.h"

#include <vector>
#include <QBoxLayout>
#include <QWidget>
#include "qtcad/defines.h"
#include "qtcad/QgTypes.h"

class QgCanvasBase;
class QImage;
class QObject;
class QgViewFilter;
class BRLObolViewController;
struct bv_context;
struct brlobol_display_endpoint;

class QTCAD_EXPORT QgView : public QWidget {
Q_OBJECT
Q_DISABLE_COPY_MOVE(QgView)


public:
explicit QgView(QWidget *parent = nullptr,
	QgViewType type = QgViewType::Auto);
~QgView();

QgViewType view_type() const;
void set_current(int);
int current();

void stash_hashes(); // Store current view hash values
bool diff_hashes();  // Request refresh if current hashes != stashed hashes.  (Does not update stored hash values - use stash_hashes for that operation.)

void save_image(int quad = 0);
void render_to_file(const QString &filename);
/* Render the current Obol view and return the pixel data. */
void get_viewport_image(QImage &img);
/* Render the Obol-canonical scene and return the raw pixel data. */
void get_obol_viewport_image(QImage &img);

bool isValid();

/* Widget-owned libbv view context used for GED and endpoint coordination. */
struct bv_context *viewContext();
const struct bv_context *viewContext() const;
BRLObolViewController *obolViewController();
struct brlobol_display_endpoint *displayEndpoint();
QgCanvasBase *canvasBase();

void aet(double a, double e, double t);

QObject *active_event_filter() const
{
return curr_event_filter;
}

// Wrappers around Qt's facility for adding eventFilter objects to
// widgets.  This is how custom key binding modes are enabled and
// disabled in QgView windows.
void add_event_filter(QObject *);
void installFilter(QgViewFilter *);

// If a filter object is supplied, remove just that filter.  If nullptr is
// passed in, remove all filters added using add_event_filter.  Does
// not clear all event filters of any sort (i.e. internal filters used
// by Qt), just those managed using these methods.
void clear_event_filter(QObject *);
void clearFilter(QgViewFilter *);

void enableDefaultKeyBindings();
void disableDefaultKeyBindings();

void enableDefaultMouseBindings();
void disableDefaultMouseBindings();

signals:
void changed(QgView *);
void init_done();

public slots:
void need_update(QgViewUpdateFlags);
void do_view_changed();
void do_init_done();
void set_lmouse_move_default(int);

private:
QBoxLayout  *l = nullptr;
QgCanvasBase *canvas = nullptr;
struct brlobol_display_endpoint *endpoint = nullptr;
QObject     *curr_event_filter = nullptr;
std::vector<QObject *> filters;
};

#endif /* QGVIEW_H */

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
