/*                      Q T C A D Q U A D . C P P
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
/** @file QgQuadView.cpp
 *
 * TODO - initialize non-current views to the standard views (check MGED for
 * what the defaults should be.)  Also, need to implement an event filter for
 * this widget (I think there is an example in qged...)  Events should pass
 * through to the current selected widget when they are either key based or
 * take place with the xy coordinates matching the current widget.  However, a
 * mouse click over the quad widget but NOT located with xy coordinates over
 * the currently selected view should change the selected view (updating
 * the session active view, perhaps changing the background or some other
 * visual signature of which widget is currently active.
 *
 * One open question is whether the faceplate settings of the previous
 * selection should be copied/transferred to the new current selection (in
 * effect, making the faceplate settings independent of the specific view
 * selected.)  Maybe this should be a widget setting, since there are arguments
 * that could be made either way...  that we we wouldn't be locked into one
 * approach at the app level.
 *
 */

#include "common.h"

#include "ged/display_obol_private.h"

#include "bv.h"

#include <QGridLayout>
#include <QMouseEvent>
#include <QtGlobal>

#include "bu/str.h"
#include "ged.h"
#include "ged/defines.h"
#include "ged/commands.h"
#include "ged/draw.h"
#include "ged/view.h"
#include "qtcad/QgQuadView.h"
#include "qtcad/QgSession.h"
#include "qtcad/QgView.h"

static const char *VIEW_NAMES[] = {"Q1", "Q2", "Q3", "Q4"};

static constexpr int
qg_quadrant_index(QgQuadrantId quadrant)
{
	return static_cast<int>(quadrant);
}

static constexpr int kUpperRightQuadrant =
	qg_quadrant_index(QgQuadrantId::UpperRight);
static constexpr int kUpperLeftQuadrant =
	qg_quadrant_index(QgQuadrantId::UpperLeft);
static constexpr int kLowerLeftQuadrant =
	qg_quadrant_index(QgQuadrantId::LowerLeft);
static constexpr int kLowerRightQuadrant =
	qg_quadrant_index(QgQuadrantId::LowerRight);
static constexpr int kQuadrantCount = kLowerRightQuadrant + 1;

static struct bv_context *
qg_quad_view_context(QgView *view)
{
    return view ? view->viewContext() : NULL;
}

static const struct bv_context *
qg_quad_view_context_const(const QgView *view)
{
    return view ? view->viewContext() : NULL;
}

static int
qg_quad_view_set_add(struct ged *gedp, QgView *view)
{
    struct ged_view_context *view_ctx = view ?
	ged_view_context_from_bv(view->viewContext()) : NULL;
    if (!gedp || !view_ctx ||
	!ged_view_set_context_add(ged_view_set_ctx(gedp), view_ctx))
	return 0;
    return ged_view_context_host_attach(gedp, view_ctx);
}

static int
qg_quad_view_set_attach(struct ged *gedp, QgView *view)
{
    struct ged_view_context *view_ctx = view ?
	ged_view_context_from_bv(view->viewContext()) : NULL;
    if (!gedp || !view_ctx ||
	!ged_view_context_view_set_attach(view_ctx, ged_view_set_ctx(gedp)))
	return 0;
    return ged_view_context_host_attach(gedp, view_ctx);
}

static int
qg_quad_view_set_remove(struct ged *gedp, QgView *view)
{
    return gedp && view ? ged_view_set_context_remove(ged_view_set_ctx(gedp),
	ged_view_context_from_bv(view->viewContext())) : 0;
}

static void
qg_quad_view_destroy(struct ged *gedp, QgView *&view)
{
    if (!view)
	return;
    (void)qg_quad_view_set_remove(gedp, view);
    delete view;
    view = nullptr;
}

static int
qg_quad_attach_endpoint(struct ged *gedp, QgView *view)
{
    if (!gedp || !view || !view->displayEndpoint())
	return 0;
    struct ged_view_context *view_ctx =
	ged_view_context_from_bv(view->viewContext());
    return ged_view_context_obol_endpoint_set(view_ctx,
	view->displayEndpoint(), 0);
}

/**
 * @brief Construct a new Qt C A D Quad:: Qt C A D Quad object
 *
 * @param parent   Parent Qt widget
 * @param session  Associated qtcad/GED session
 * @param type     Requesting either a GL or SWRAST display mechanism
 */
QgQuadView::QgQuadView(QWidget *parent, QgSession *session,
	QgViewType type) : QWidget(parent)
{
	m_session = session;
	graphicsType = type;

	views[kUpperRightQuadrant] = createView(kUpperRightQuadrant);
	struct ged *gedp = m_session ? m_session->ged() : nullptr;
	if (gedp)
		qg_quad_view_set_add(gedp, views[kUpperRightQuadrant]);
	if (m_session)
		m_session->setActiveViewContext(
		    views[kUpperRightQuadrant]->viewContext());

	views[kUpperRightQuadrant]->set_current(1);
	currentView = views[kUpperRightQuadrant];

	default_views(1);
}

QgQuadView::~QgQuadView()
{
	struct ged *gedp = m_session ? m_session->ged() : nullptr;
	if (m_session) {
		struct bv_context *active_view = m_session->activeViewContext();
		for (int i = kUpperRightQuadrant;
			i <= kLowerRightQuadrant; i++) {
			if (views[i] && views[i]->viewContext() == active_view) {
				m_session->setActiveViewContext(nullptr);
				break;
			}
		}
	}
	for (int i = 0; i < 4; i++) {
		qg_quad_view_destroy(gedp, views[i]);
	}

	if (spacerTop)
		delete spacerTop;
	if (spacerBottom)
		delete spacerBottom;
	if (spacerLeft)
		delete spacerLeft;
	if (spacerRight)
		delete spacerRight;
	if (spacerCenter)
		delete spacerCenter;
}

QgView *
QgQuadView::curr_view()
{
	return get(get_selected());
}

/**
 * @brief Creates a view for the viewport. Convenience method of common things that need to be done to the view.
 *
 * @param index of the view names to use from the constant list of names
 * @return QgView*
 */
QgView *
QgQuadView::createView(unsigned int index)
{
	QgView *view = new QgView(this, graphicsType);
	static const char *test_ids[] = {
	    "view-upper-right",
	    "view-upper-left",
	    "view-lower-left",
	    "view-lower-right"
	};
	if (index < sizeof(test_ids) / sizeof(test_ids[0]))
	    view->setProperty("qgTestId",
		QString::fromLatin1(test_ids[index]));
	bv_context_name_set(qg_quad_view_context(view), VIEW_NAMES[index]);
	view->set_current(0);
	view->installEventFilter(this);

	struct ged *gedp = m_session ? m_session->ged() : nullptr;
	qg_quad_view_set_attach(gedp, view);
	(void)qg_quad_attach_endpoint(gedp, view);
	/* Each viewport owns independent libbv view state. */

	QObject::connect(view, &QgView::changed, this, &QgQuadView::do_view_changed);
	QObject::connect(view, &QgView::init_done, this, &QgQuadView::do_init_done);
	return view;
}

/**
 * @brief Convenience method to create the layout and set common parameters.
 *
 * @return QGridLayout*
 */
QGridLayout *
QgQuadView::createLayout()
{
	QGridLayout *layout = new QGridLayout(this);
	layout->setSpacing(0);
	layout->setContentsMargins(0, 0, 0, 0);
	layout->setAlignment(Qt::AlignTop | Qt::AlignLeft);

	this->setLayout(layout);

	if (currentLayout != nullptr) {
		delete currentLayout;
	}
	currentLayout = layout;

	return layout;
}

/**
 * @brief Changes the viewport layout to only have the single view.  We destroy the other views if needed and the flag is set
 *
 */
void
QgQuadView::changeToSingleFrame()
{
	struct ged *gedp = m_session ? m_session->ged() : nullptr;
	QGridLayout *layout = (QGridLayout *)this->layout();
	if (layout == nullptr) {
		layout = createLayout();
	}
	QLayoutItem *layout_item = nullptr;
	while ((layout_item = layout->takeAt(0)) != nullptr)
		delete layout_item;
	layout->addWidget(views[kUpperRightQuadrant], 0, 2);

	/* A secondary quadrant may have been active.  Select the surviving view
	 * before releasing its peers so GED never retains a dangling active context. */
	if (m_session)
		m_session->setActiveViewContext(
			views[kUpperRightQuadrant]->viewContext());

	for (int i = 1; i < 4; i++) {
		// Don't want use cpu for views that are not visible
		if (views[i] != nullptr) {
			views[i]->disconnect();
			qg_quad_view_destroy(gedp, views[i]);
		}
	}

	views[kUpperRightQuadrant]->set_current(1);
	currentView = views[kUpperRightQuadrant];

	// No need to indicate active quad
	delete spacerTop;
	delete spacerBottom;
	delete spacerLeft;
	delete spacerRight;
	delete spacerCenter;
	spacerTop = nullptr;
	spacerBottom = nullptr;
	spacerLeft = nullptr;
	spacerRight = nullptr;
	spacerCenter = nullptr;

	default_views(0);
}

/**
 * @brief Changes the viewport layout to have 4 views the views will be equal size and not resizeable.  This will create the extra view if needed.
 *
 */
void
QgQuadView::changeToQuadFrame()
{
	struct ged *gedp = m_session ? m_session->ged() : nullptr;
	for (int i = kUpperRightQuadrant + 1; i < kLowerRightQuadrant + 1; i++) {
		if (views[i] == nullptr) {
			// Make the new view
			views[i] = createView(i);

			// Out of the gate, have the new view units match the first view's
			// units (which should usually be based on the database units)
			const struct bv_context *source_view_ctx =
				qg_quad_view_context_const(views[0]);
			bv_unit_conversion_set(bv_context_view(qg_quad_view_context(views[i])),
				bv_local2base_get(bv_context_view_const(source_view_ctx)),
				bv_base2local_get(bv_context_view_const(source_view_ctx)));

			// For initial layout calculations, we need to set a screen width
			// and height.  This won't be right in the end, but it gives
			// the LoD bounds update something to work with
			const struct bv_context *layout_view_ctx =
				qg_quad_view_context_const(views[kUpperRightQuadrant]);
			bv_context_dimensions_set(qg_quad_view_context(views[i]),
				bv_width_get(bv_context_view_const(layout_view_ctx)),
				bv_height_get(bv_context_view_const(layout_view_ctx)));
		}
		// Copy the LoD source policy so all quadrants use the same
		// source-selection behavior.
		ged_view_lod_policy lod_policy = BV_LOD_POLICY_INIT;
		if (ged_view_lod_policy_get(&lod_policy,
				ged_view_context_from_bv(
				    views[kUpperRightQuadrant]->viewContext()))) {
			ged_view_lod_policy_apply(
				ged_view_context_from_bv(views[i]->viewContext()),
				&lod_policy);
		}
		if (gedp)
			qg_quad_view_set_add(gedp, views[i]);
	}

	// Define the spacers
	spacerTop = new QWidget;
	spacerTop->setMinimumWidth(1);
	spacerTop->setMaximumWidth(1);
	spacerTop->setStyleSheet("");
	spacerBottom = new QWidget;
	spacerBottom->setMinimumWidth(1);
	spacerBottom->setMaximumWidth(1);
	spacerBottom->setStyleSheet("");
	spacerLeft = new QWidget;
	spacerLeft->setMinimumHeight(1);
	spacerLeft->setMaximumHeight(1);
	spacerLeft->setStyleSheet("");
	spacerRight = new QWidget;
	spacerRight->setMinimumHeight(1);
	spacerRight->setMaximumHeight(1);
	spacerRight->setStyleSheet("");
	spacerCenter = new QWidget;
	spacerCenter->setMinimumSize(1,1);
	spacerCenter->setMaximumSize(1,1);
	// Something is always selected, so the center widget is always colored
	// accordingly.
	spacerCenter->setStyleSheet("background-color:yellow;");

	QGridLayout *layout = (QGridLayout *)this->layout();
	if (layout == nullptr) {
		layout = createLayout();
	}
	QLayoutItem *layout_item = nullptr;
	while ((layout_item = layout->takeAt(0)) != nullptr)
		delete layout_item;

	layout->addWidget(views[kUpperLeftQuadrant], 0, 0);
	layout->addWidget(spacerTop, 0, 1);
	layout->addWidget(views[kUpperRightQuadrant], 0, 2);
	layout->addWidget(spacerLeft, 1, 0);
	layout->addWidget(spacerCenter, 1, 1);
	layout->addWidget(spacerRight, 1, 2);
	layout->addWidget(views[kLowerLeftQuadrant], 2, 0);
	layout->addWidget(spacerBottom, 2, 1);
	layout->addWidget(views[kLowerRightQuadrant], 2, 2);

	default_views(0);

	/* Every non-independent view is attached to GED's shared semantic scene
	 * root when its endpoint is registered.  Updating the view bounds is
	 * sufficient to seed each pane's independent LoD policy; reconstructing
	 * semantic draw intent from rendered geometry would be both incorrect and
	 * O(number of realized primitives). */
	for (int i = kUpperRightQuadrant + 1; i < kLowerRightQuadrant + 1; i++) {
		struct ged_view_context *view_ctx =
		    ged_view_context_from_bv(views[i]->viewContext());
		ged_view_lod_bounds_update(view_ctx);
	}

	const struct bv_context *layout_view_ctx =
		qg_quad_view_context_const(views[kUpperRightQuadrant]);
	for (int i = kUpperRightQuadrant + 1; i < kLowerRightQuadrant + 1; i++) {
		bv_context_dimensions_set(qg_quad_view_context(views[i]),
			bv_width_get(bv_context_view_const(layout_view_ctx)),
			bv_height_get(bv_context_view_const(layout_view_ctx)));
	}

	// Current view selection pieces
	select(QgQuadrantId::UpperRight);
	if (m_session)
		m_session->setActiveViewContext(
		    views[kUpperRightQuadrant]->viewContext());
	views[kUpperRightQuadrant]->set_current(1);
	currentView = views[kUpperRightQuadrant];
}

void
QgQuadView::do_view_changed()
{
	QTCAD_SLOT("QgQuadView::do_view_changed", 1);
	emit changed(currentView);
}

bool
QgQuadView::isValid()
{
	for (int i = kUpperRightQuadrant; i < kLowerRightQuadrant + 1; i++) {
		if (views[i] != nullptr && !views[i]->isValid())
			return false;
	}
	return true;
}

bool
QgQuadView::eventFilter(QObject *t, QEvent *e)
{
	if (e->type() == QEvent::KeyPress || e->type() == QEvent::MouseButtonPress) {
		for (int i = kUpperRightQuadrant; i < kLowerRightQuadrant + 1; i++) {
			if (views[i] != nullptr && t == views[i]) {
				select(static_cast<QgQuadrantId>(i));
				break;
			}
		}
	}
	return false;
}

void
QgQuadView::default_views(int all_views)
{
	if (all_views && views[kUpperRightQuadrant] != nullptr) {
		if (views[kUpperLeftQuadrant] == nullptr) {
			views[kUpperRightQuadrant]->aet(270, 90, 0);
		}
		else {
			views[kUpperRightQuadrant]->aet(35, 25, 0);
		}
	}
	if (views[kUpperLeftQuadrant] != nullptr) {
		views[kUpperLeftQuadrant]->aet(0, 90, 0);
	}
	if (views[kLowerLeftQuadrant] != nullptr) {
		views[kLowerLeftQuadrant]->aet(0, 0, 0);
	}
	if (views[kLowerRightQuadrant] != nullptr) {
		views[kLowerRightQuadrant]->aet(90, 0, 0);
	}
}

QgView *
QgQuadView::get(QgQuadrantId quadrant)
{
	const int quadrant_index = qg_quadrant_index(quadrant);
	if (quadrant_index < kUpperRightQuadrant ||
		quadrant_index >= kQuadrantCount)
		return currentView;
	if (views[quadrant_index] != nullptr) {
		return views[quadrant_index];
	}

	return currentView;
}

QgView *
QgQuadView::get(const QPoint &gpos)
{
	QgView *retv = nullptr;
	for (int i = kUpperRightQuadrant; i < kLowerRightQuadrant + 1; i++) {
		QgView *cv = views[i];
		if (cv == nullptr)
			continue;
		QWidget *cw = (QWidget *)cv;
		QRect br = cw->geometry();
		QWidget *pcw = (QWidget *)cw->parent();
		QPoint lp = pcw->mapFromGlobal(gpos);
		if (br.contains(lp)) {
			retv = cv;
			break;
		}
	}

	return retv;
}

QgView *
QgQuadView::get(QEvent *e)
{
	if (e->type() != QEvent::MouseButtonPress)
		return nullptr;
	QMouseEvent *m_e = (QMouseEvent *)e;
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
	QPoint gpos = m_e->globalPos();
#else
	QPoint gpos = m_e->globalPosition().toPoint();
#endif
	return get(gpos);
}

void
QgQuadView::select(QgQuadrantId quadrant)
{
	const int quadrant_index = qg_quadrant_index(quadrant);
	if (quadrant_index < kUpperRightQuadrant ||
		quadrant_index >= kQuadrantCount)
		return;

	QgView *oc = currentView;

	// Set new selection
	if (views[quadrant_index] != nullptr) {
		currentView = views[quadrant_index];
	}

	// Clear any old selections
	if (spacerTop)
		spacerTop->setStyleSheet("");
	if (spacerBottom)
		spacerBottom->setStyleSheet("");
	if (spacerLeft)
		spacerLeft->setStyleSheet("");
	if (spacerRight)
		spacerRight->setStyleSheet("");

	// If we're not in Quad mode, done
	if (views[kUpperLeftQuadrant] == nullptr)
		return;

	// If we're in quad mode, more work to do
	currentView->set_current(1);

	for (int i = kUpperRightQuadrant; i < kLowerRightQuadrant + 1; i++) {
		if (i == quadrant_index)
			continue;
		if (views[i] != nullptr) {
			views[i]->set_current(0);
		}
	}

	if (quadrant == QgQuadrantId::UpperRight) {
		spacerTop->setStyleSheet("background-color:yellow;");
		spacerRight->setStyleSheet("background-color:yellow;");
	}

	if (quadrant == QgQuadrantId::UpperLeft) {
		spacerTop->setStyleSheet("background-color:yellow;");
		spacerLeft->setStyleSheet("background-color:yellow;");
	}

	if (quadrant == QgQuadrantId::LowerLeft) {
		spacerBottom->setStyleSheet("background-color:yellow;");
		spacerLeft->setStyleSheet("background-color:yellow;");
	}

	if (quadrant == QgQuadrantId::LowerRight) {
		spacerBottom->setStyleSheet("background-color:yellow;");
		spacerRight->setStyleSheet("background-color:yellow;");
	}

	if (oc != currentView)
		emit selected(currentView);
}

void
QgQuadView::select(const char *quadrant_id)
{
	if (BU_STR_EQUIV(quadrant_id, "ur")) {
		select(QgQuadrantId::UpperRight);
		return;
	}
	if (BU_STR_EQUIV(quadrant_id, "ul")) {
		select(QgQuadrantId::UpperLeft);
		return;
	}
	if (BU_STR_EQUIV(quadrant_id, "ll")) {
		select(QgQuadrantId::LowerLeft);
		return;
	}
	if (BU_STR_EQUIV(quadrant_id, "lr")) {
		select(QgQuadrantId::LowerRight);
		return;
	}
}


QgQuadrantId
QgQuadView::get_selected()
{
	if (currentView == views[kUpperRightQuadrant]) {
		return QgQuadrantId::UpperRight;
	}
	if (currentView == views[kUpperLeftQuadrant]) {
		return QgQuadrantId::UpperLeft;
	}
	if (currentView == views[kLowerLeftQuadrant]) {
		return QgQuadrantId::LowerLeft;
	}
	if (currentView == views[kLowerRightQuadrant]) {
		return QgQuadrantId::LowerRight;
	}

	return QgQuadrantId::UpperRight;
}

void
QgQuadView::do_view_update(QgViewUpdateFlags flags)
{
	QTCAD_SLOT("QgQuadView::do_view_update", 1);
	for (int i = kUpperRightQuadrant; i < kLowerRightQuadrant + 1; i++) {
		if (views[i] != nullptr) {
			views[i]->need_update(flags);
		}
	}
}

void
QgQuadView::stash_hashes()
{
	for (int i = kUpperRightQuadrant; i < kLowerRightQuadrant + 1; i++) {
		if (views[i] != nullptr) {
			views[i]->stash_hashes();
		}
	}
}

bool
QgQuadView::diff_hashes()
{
	bool ret = false;
	for (int i = kUpperRightQuadrant; i < kLowerRightQuadrant + 1; i++) {
		if (views[i] != nullptr) {
			if (views[i]->diff_hashes()) {
				ret = true;
			}
		}
	}

	return ret;
}

void
QgQuadView::enableDefaultKeyBindings()
{
	for (int i = kUpperRightQuadrant; i < kLowerRightQuadrant + 1; i++) {
		if (views[i] != nullptr) {
			views[i]->enableDefaultKeyBindings();
		}
	}
}

void
QgQuadView::disableDefaultKeyBindings()
{
	for (int i = kUpperRightQuadrant; i < kLowerRightQuadrant + 1; i++) {
		if (views[i] != nullptr) {
			views[i]->disableDefaultKeyBindings();
		}
	}
}

void
QgQuadView::enableDefaultMouseBindings()
{
	for (int i = kUpperRightQuadrant; i < kLowerRightQuadrant + 1; i++) {
		if (views[i] != nullptr) {
			views[i]->enableDefaultMouseBindings();
		}
	}
}

void
QgQuadView::disableDefaultMouseBindings()
{
	for (int i = kUpperRightQuadrant; i < kLowerRightQuadrant + 1; i++) {
		if (views[i] != nullptr) {
			views[i]->disableDefaultMouseBindings();
		}
	}
}

void
QgQuadView::set_lmouse_move_default(int mm)
{
	QTCAD_SLOT("QgQuadView::set_lmouse_move_default", 1);
	for (int i = kUpperRightQuadrant; i < kLowerRightQuadrant + 1; i++) {
		if (views[i] != nullptr) {
			views[i]->set_lmouse_move_default(mm);
		}
	}
}

void
QgQuadView::do_init_done()
{
	QTCAD_SLOT("QgQuadView::do_init_done", 1);
	if (!init_done_flag) {
		init_done_flag = true;
		emit init_done();
	}
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
