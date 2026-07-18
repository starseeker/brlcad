/*                      Q G V I E W . C P P
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
/** @file QgView.cpp
 *
 * Wrapper widget that handles the various widget types which may
 * constitute a Qt based geometry view.
 *
 */

#include "common.h"

#include "BObol/BDisplayEndpoint.h"
#include "qtcad/QgCanvasBase.h"
#include "qtcad/QgGL.h"
#include "qtcad/QgSW.h"
#include "qtcad/QgView.h"
#include "qtcad/QgViewFilter.h"
#include "qtcad/QgSignalFlags.h"
#include "qtcad/QgObolWindowHost.h"
#include "bv.h"
#include "ged/view.h"

extern "C" {
#include "bu/malloc.h"
}

#include <algorithm>
#include <cmath>
#include <cstring>

static uint32_t
qg_refresh_flags(QgViewUpdateFlags flags)
{
    uint32_t refresh_flags = 0;

    if (!flags)
	return BV_REFRESH_ALL;
    if (flags & QG_VIEW_REFRESH)
	refresh_flags |= BV_REFRESH_VIEW;
    if (flags & QG_VIEW_DRAWN)
	refresh_flags |= BV_REFRESH_DRAW;
    if (flags & QG_VIEW_SELECT)
	refresh_flags |= BV_REFRESH_OVERLAY;
    if (flags & QG_VIEW_MODE)
	refresh_flags |= BV_REFRESH_EDIT;
    if (flags & QG_VIEW_DB)
	refresh_flags |= BV_REFRESH_DRAW;

    return refresh_flags ? refresh_flags : BV_REFRESH_ALL;
}

/* ---------------------------------------------------------------------- */
/* Factory: create the appropriate canvas widget for the requested type.  */
/* The single BRLCAD_OPENGL guard lives here rather than being duplicated  */
/* in every QgView method body.                                            */
/* ---------------------------------------------------------------------- */
static QgCanvasBase *
make_canvas(QWidget *parent, QgViewType type)
{
    switch (type) {
#ifdef BRLCAD_OPENGL
    case QgViewType::GL:
return new QgGL(parent, nullptr, false);
#endif
    case QgViewType::SW:
return new QgSW(parent, nullptr, false);
    default:
/* Auto prefers hardware GL and otherwise uses software presentation. */
#ifdef BRLCAD_OPENGL
return new QgGL(parent, nullptr, false);
#else
return new QgSW(parent, nullptr, false);
#endif
    }
}

/* ---------------------------------------------------------------------- */

/* The view, rather than a native host, initializes a freshly created
 * endpoint from its passive GED state.  A host may be replaced while the
 * endpoint survives, so it must never own this policy transfer. */
static bool
qg_endpoint_background_seed(bobol_display_endpoint_t *endpoint,
	struct bv_context *view_ctx)
{
    if (!endpoint || !view_ctx)
	return false;

    struct bv_background_state background = BV_BACKGROUND_STATE_INIT;
    if (!bv_background_state_get(&background,
	    bv_context_view_const(view_ctx)))
	return false;

    struct bobol_endpoint_property_value value =
	BOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    value.type = BOBOL_ENDPOINT_PROPERTY_COLOR3;
    for (int i = 0; i < 3; i++)
	value.color3[i] = background.bottom[i] / 255.0;
    if (bobol_display_endpoint_property_set(endpoint,
	    "controller.background.bottom", &value) !=
	BOBOL_ENDPOINT_PROPERTY_OK)
	return false;
    for (int i = 0; i < 3; i++)
	value.color3[i] = background.top[i] / 255.0;
    return bobol_display_endpoint_property_set(endpoint,
	"controller.background.top", &value) ==
	BOBOL_ENDPOINT_PROPERTY_OK;
}

/* ---------------------------------------------------------------------- */

QgView::QgView(QWidget *parent, QgViewType type)
    : QWidget(parent)
{
    this->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    l = new QBoxLayout(QBoxLayout::LeftToRight, this);
    l->setSpacing(0);
    l->setContentsMargins(0, 0, 0, 0);

    canvas = make_canvas(this, type);
    if (!canvas)
return;

    QWidget *w = canvas->canvasWidget();
    w->setMinimumSize(50, 50);
    w->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    l->addWidget(w);

    if (qtcad_obol_host_factories_register()) {
	struct bobol_host_desc desc;
	std::memset(&desc, 0, sizeof(desc));
	desc.struct_size = sizeof(desc);
	desc.mode = BOBOL_HOST_MODE_EMBEDDED;
	const qreal dpr = w->devicePixelRatioF();
	desc.width = static_cast<unsigned int>(std::max(1,
	    static_cast<int>(std::ceil(w->width() * dpr))));
	desc.height = static_cast<unsigned int>(std::max(1,
	    static_cast<int>(std::ceil(w->height() * dpr))));
	desc.device_pixel_ratio = dpr > 0.0 ? dpr : 1.0;
	desc.visible = w->isVisible() ? 1 : 0;
	desc.application_context = canvas;
#ifdef BRLCAD_OPENGL
	const bool use_gl = dynamic_cast<QgGL *>(canvas) != nullptr;
#else
	const bool use_gl = false;
#endif
	desc.required_capabilities = BOBOL_HOST_CAP_EMBEDDED |
	    BOBOL_HOST_CAP_READBACK | (use_gl ?
	    BOBOL_HOST_CAP_SYSTEM_GL : BOBOL_HOST_CAP_PIXEL_PRESENT);
	endpoint = bobol_display_endpoint_create(NULL, 0);
	if (!endpoint || !bobol_display_endpoint_render_engine_set(endpoint,
		use_gl ? BOBOL_RENDER_ENGINE_HW : BOBOL_RENDER_ENGINE_SW) ||
	    !bobol_display_endpoint_host_open(endpoint,
		use_gl ? "qt-gl" : "qt-sw", &desc) ||
	    !qg_endpoint_background_seed(endpoint, viewContext())) {
	    canvas->setObolViewController(nullptr);
	    bobol_display_endpoint_destroy(endpoint);
	    endpoint = nullptr;
	}
    } else {
	canvas->setObolViewController(nullptr);
    }

    /* Connect canvas signals via old-style macros (QgCanvasBase is not a
     * QObject, so we obtain the QObject* via asQObject()). */
    QObject::connect(canvas->asQObject(), SIGNAL(changed()),
     this, SLOT(do_view_changed()));
    QObject::connect(canvas->asQObject(), SIGNAL(init_done()),
     this, SLOT(do_init_done()));
}

QgView::~QgView()
{
    struct bv_context *view_ctx = viewContext();
    if (endpoint && view_ctx &&
	ged_view_context_display_endpoint_get(
	    ged_view_context_from_bv(view_ctx)) == endpoint)
	(void)ged_view_context_display_endpoint_set(
	    ged_view_context_from_bv(view_ctx), nullptr, 0);
    if (canvas) {
	canvas->setObolInputEndpoint(nullptr);
	canvas->setObolViewController(nullptr);
    }
    bobol_display_endpoint_destroy(endpoint);
    endpoint = nullptr;
    delete canvas;
}

bool
QgView::isValid()
{
    if (!canvas || !endpoint)
return false;
    return canvas->isValid();
}

QgViewType
QgView::view_type() const
{
    if (!canvas)
	return QgViewType::Auto;
#ifdef BRLCAD_OPENGL
    if (dynamic_cast<QgGL *>(canvas))
	return QgViewType::GL;
#endif
    return QgViewType::SW;
}


void
QgView::save_image(int UNUSED(quad))
{
}

void
QgView::render_to_file(const QString &filename)
{
    if (canvas)
canvas->render_to_file(filename);
}

void
QgView::get_viewport_image(QImage &img)
{
    if (canvas)
canvas->get_viewport_image(img);
    else
img = QImage();
}

void
QgView::get_obol_viewport_image(QImage &img)
{
    if (canvas)
canvas->get_obol_viewport_image(img);
    else
img = QImage();
}

void
QgView::do_view_changed()
{
    QTCAD_SLOT("QgView::do_view_changed", 1);
    emit changed(this);
}

void
QgView::need_update(QgViewUpdateFlags flags)
{
    QTCAD_SLOT("QgView::need_update", 1);
    uint32_t refresh_flags = qg_refresh_flags(flags);
    if (struct bv_context *view_ctx = viewContext())
	bv_refresh_request(bv_context_view(view_ctx), refresh_flags);
    if (canvas)
canvas->request_update(refresh_flags);
}

struct bv_context *
QgView::viewContext()
{
    return canvas ? canvas->viewContext() : nullptr;
}

const struct bv_context *
QgView::viewContext() const
{
	return canvas ? canvas->viewContext() : nullptr;
}

BObolViewController *
QgView::obolViewController()
{
    return canvas ? canvas->obolViewController() : nullptr;
}

struct bobol_display_endpoint *
QgView::displayEndpoint()
{
    return endpoint;
}

QgCanvasBase *
QgView::canvasBase()
{
    return canvas;
}

void
QgView::stash_hashes()
{
    if (canvas)
canvas->stash_hashes();
}

bool
QgView::diff_hashes()
{
    return canvas ? canvas->diff_hashes() : false;
}

void
QgView::aet(double a, double e, double t)
{
    if (canvas)
canvas->aet(a, e, t);
}

void
QgView::set_current(int i)
{
    if (canvas)
canvas->set_current(i);
}

int
QgView::current()
{
    return canvas ? canvas->currentView() : 0;
}

void
QgView::add_event_filter(QObject *o)
{
    curr_event_filter = o;
    filters.push_back(o);
    if (canvas)
canvas->canvasWidget()->installEventFilter(o);
}

void
QgView::installFilter(QgViewFilter *f)
{
    if (!f)
	return;
    f->set_view_widget(this);
    add_event_filter(f);
}

void
QgView::clear_event_filter(QObject *o)
{
    if (!canvas)
return;

    QWidget *w = canvas->canvasWidget();
    if (o) {
w->removeEventFilter(o);
auto fit = std::find(filters.begin(), filters.end(), o);
if (fit != filters.end())
    filters.erase(fit);
    } else {
for (size_t i = 0; i < filters.size(); i++)
    w->removeEventFilter(filters[i]);
filters.clear();
    }

    /* Passing nullptr is the documented "clear all managed filters" mode. */
    if (!o || curr_event_filter == o)
curr_event_filter = nullptr;
}

void
QgView::clearFilter(QgViewFilter *f)
{
    if (!f)
	return;
    clear_event_filter(f);
    f->set_view_widget(nullptr);
}

void
QgView::enableDefaultKeyBindings()
{
    if (canvas)
canvas->enableDefaultKeyBindings();
}

void
QgView::disableDefaultKeyBindings()
{
    if (canvas)
canvas->disableDefaultKeyBindings();
}

void
QgView::enableDefaultMouseBindings()
{
    if (canvas)
canvas->enableDefaultMouseBindings();
}

void
QgView::disableDefaultMouseBindings()
{
    if (canvas)
canvas->disableDefaultMouseBindings();
}


void
QgView::set_lmouse_move_default(int mm)
{
    QTCAD_SLOT("QgView::set_lmouse_move_default", 1);
    if (canvas)
canvas->set_lmouse_move_default(mm);
}


void
QgView::do_init_done()
{
    QTCAD_SLOT("QgView::do_init_done", 1);
    emit init_done();
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
