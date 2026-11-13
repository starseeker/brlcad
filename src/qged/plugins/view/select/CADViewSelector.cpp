/*             C A D V I E W S E L E C T O R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2023-2026 United States Government as represented by
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
/** @file CADViewSelector.cpp
 *
 * Brief description
 *
 */

#include "common.h"

#include "BObol/BDisplayEndpoint.h"
#include <QVBoxLayout>
#include <string>
#include <vector>
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgSignalFlags.h"

#include "bu/vls.h"
#include "bg/aabb_ray.h"
#include "bg/plane.h"
#include "ged.h"
#include "ged/scene.h"
#include "ged/selection.h"

#include "qtcad/QgSelectFilter.h"
#include "qtcad/QgView.h"
#include "./CADViewSelector.h"

static struct ged_view_context *
qged_selector_view(const QgPluginContext *ctx)
{
    return ctx ? ged_view_context_from_bv(ctx->activeViewContext()) : nullptr;
}

static enum ged_scene_status
qged_scene_draw_paths(struct ged *gedp, struct ged_view_context *view,
	const std::vector<std::string> &paths, int *changed)
{
    if (changed)
	*changed = 0;
    if (!gedp || !view || paths.empty())
	return GED_SCENE_INVALID;

    std::vector<const char *> path_values;
    path_values.reserve(paths.size());
    for (const std::string &path : paths)
	path_values.push_back(path.c_str());

    struct ged_scene_draw_request request;
    ged_scene_draw_request_init(&request);
    request.view = view;
    request.paths = path_values.data();
    request.path_count = path_values.size();

    struct ged_scene_result *result = ged_scene_result_create();
    if (!result)
	return GED_SCENE_ERROR;
    const enum ged_scene_status status = ged_scene_draw(gedp, &request,
	result);
    if (changed)
	*changed = ged_scene_result_changed(result);
    ged_scene_result_destroy(result);
    return status;
}

static enum ged_scene_status
qged_scene_erase_paths(struct ged *gedp, struct ged_view_context *view,
	const std::vector<std::string> &paths, int *changed)
{
    if (changed)
	*changed = 0;
    if (!gedp || !view || paths.empty())
	return GED_SCENE_INVALID;

    struct ged_scene_result *result = ged_scene_result_create();
    if (!result)
	return GED_SCENE_ERROR;

    enum ged_scene_status status = GED_SCENE_OK;
    int scene_changed = 0;
    for (const std::string &path : paths) {
	struct ged_scene_erase_request request;
	ged_scene_erase_request_init(&request);
	request.view = view;
	request.path = path.c_str();
	status = ged_scene_erase(gedp, &request, result);
	if (status != GED_SCENE_OK)
	    break;
	scene_changed |= ged_scene_result_changed(result);
    }
    if (changed)
	*changed = scene_changed;
    ged_scene_result_destroy(result);
    return status;
}

const BObolInputActionLayer *
CADViewSelector::inputActionLayer()
{
    static const unsigned int modifier_mask = BOBOL_INPUT_MOD_SHIFT |
	BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_ALT |
	BOBOL_INPUT_MOD_META;
    static const BObolInputBinding bindings[] = {
	{BOBOL_INPUT_POINTER_PRESS, BOBOL_INPUT_ANY, 0, 0, modifier_mask,
	 10, QG_SELECT_INPUT_BEGIN},
	{BOBOL_INPUT_POINTER_MOTION, BOBOL_INPUT_ANY, 0, 0, modifier_mask,
	 10, QG_SELECT_INPUT_UPDATE},
	{BOBOL_INPUT_POINTER_RELEASE, BOBOL_INPUT_ANY, 0, 0, modifier_mask,
	 10, QG_SELECT_INPUT_COMMIT}
    };
    static const BObolInputActionLayer layer = {
	"qged-selection", bindings, sizeof(bindings) / sizeof(bindings[0]),
	CADViewSelector::inputActionDispatch
    };
    return &layer;
}

static std::vector<std::string>
qged_selection_paths(struct ged *gedp)
{
    std::vector<std::string> ret;
    if (!ged_selection_state_available(gedp))
	return ret;

    struct bu_vls paths = BU_VLS_INIT_ZERO;
    if (!ged_selection_list_paths(gedp, nullptr, &paths)) {
	bu_vls_free(&paths);
	return ret;
    }

    const char *pstr = bu_vls_cstr(&paths);
    const char *start = pstr;
    for (const char *c = pstr; c && *c; c++) {
	if (*c != '\n')
	    continue;
	if (c > start)
	    ret.push_back(std::string(start, (size_t)(c - start)));
	start = c + 1;
    }
    if (start && *start)
	ret.push_back(std::string(start));
    bu_vls_free(&paths);

    return ret;
}

CADViewSelector::CADViewSelector(QWidget *)
{
    const QString testPrefix = QStringLiteral("org.brlcad.qged.view.select");
    this->setProperty("qgTestId", testPrefix + QStringLiteral(".controls"));
    QVBoxLayout *wl = new QVBoxLayout;
    wl->setAlignment(Qt::AlignTop);

    QGroupBox *sstyle_box = new QGroupBox("Selection Style");
    QButtonGroup *sstyle_grp = new QButtonGroup();
    QVBoxLayout *sstyle_gl = new QVBoxLayout;
    sstyle_gl->setAlignment(Qt::AlignTop);

    use_pnt_select_button = new QRadioButton("Select Under Point");
    use_pnt_select_button->setProperty("qgTestId",
	testPrefix + QStringLiteral(".point"));
    use_pnt_select_button->setChecked(true);
    sstyle_grp->addButton(use_pnt_select_button);
    sstyle_gl->addWidget(use_pnt_select_button);

    use_rect_select_button = new QRadioButton("Select Under Rectangle");
    use_rect_select_button->setProperty("qgTestId",
	testPrefix + QStringLiteral(".rectangle"));
    sstyle_grp->addButton(use_rect_select_button);
    sstyle_gl->addWidget(use_rect_select_button);

    sstyle_box->setLayout(sstyle_gl);

    select_all_depth_ckbx = new QCheckBox("Use All Intersections");
    select_all_depth_ckbx->setProperty("qgTestId",
	testPrefix + QStringLiteral(".all-intersections"));
    use_ray_test_ckbx = new QCheckBox("Test with Raytracing");
    use_ray_test_ckbx->setProperty("qgTestId",
	testPrefix + QStringLiteral(".raytrace"));

    QObject::connect(use_pnt_select_button, &QRadioButton::clicked, this, &CADViewSelector::enable_raytrace_opt);
    QObject::connect(use_pnt_select_button, &QRadioButton::clicked, this, &CADViewSelector::enable_useall_opt);
    QObject::connect(use_rect_select_button, &QRadioButton::clicked, this, &CADViewSelector::disable_raytrace_opt);
    QObject::connect(use_rect_select_button, &QRadioButton::clicked, this, &CADViewSelector::disable_useall_opt);

    QGroupBox *smode_box = new QGroupBox("Selection Mode");
    QVBoxLayout *smode_gl = new QVBoxLayout;
    smode_gl->setAlignment(Qt::AlignTop);

    QButtonGroup *smode_grp = new QButtonGroup();
    erase_from_scene_button = new QRadioButton("Erase from Scene");
    erase_from_scene_button->setProperty("qgTestId",
	testPrefix + QStringLiteral(".erase"));
    erase_from_scene_button->setChecked(true);
    smode_grp->addButton(erase_from_scene_button);
    smode_gl->addWidget(erase_from_scene_button);

    add_to_group_button = new QRadioButton("Add to Current Set");
    add_to_group_button->setProperty("qgTestId",
	testPrefix + QStringLiteral(".add"));
    smode_grp->addButton(add_to_group_button);
    smode_gl->addWidget(add_to_group_button);

    rm_from_group_button = new QRadioButton("Remove from Current Set");
    rm_from_group_button->setProperty("qgTestId",
	testPrefix + QStringLiteral(".remove"));
    smode_grp->addButton(rm_from_group_button);
    smode_gl->addWidget(rm_from_group_button);

    QGroupBox *wgroups = new QGroupBox("Current Selection Set");
    QVBoxLayout *sgrp_gl = new QVBoxLayout;
    sgrp_gl->setAlignment(Qt::AlignTop);


    draw_selections = new QPushButton("Draw selected");
    draw_selections->setProperty("qgTestId",
	testPrefix + QStringLiteral(".draw-selected"));
    sgrp_gl->addWidget(draw_selections);
    QObject::connect(draw_selections, &QPushButton::clicked, this, &CADViewSelector::do_draw_selections);

    erase_selections = new QPushButton("Erase selected");
    erase_selections->setProperty("qgTestId",
	testPrefix + QStringLiteral(".erase-selected"));
    sgrp_gl->addWidget(erase_selections);
    QObject::connect(erase_selections, &QPushButton::clicked, this, &CADViewSelector::do_erase_selections);


    QWidget *sgrp = new QWidget();
    QHBoxLayout *groups_gl = new QHBoxLayout;
    groups_gl->setAlignment(Qt::AlignLeft);
    //current_group = new QComboBox();
    //groups_gl->addWidget(current_group);
    //add_new_group = new QPushButton("+");
    //groups_gl->addWidget(add_new_group);
    //rm_group = new QPushButton("-");
    //groups_gl->addWidget(rm_group);
    sgrp->setLayout(groups_gl);
    sgrp_gl->addWidget(sgrp);

    group_contents = new QListWidget();
    group_contents->setProperty("qgTestId",
	testPrefix + QStringLiteral(".selection-list"));
    sgrp_gl->addWidget(group_contents);

    disable_groups(false);

    QObject::connect(erase_from_scene_button , &QRadioButton::clicked, this, &CADViewSelector::disable_groups);
    QObject::connect(add_to_group_button , &QRadioButton::clicked, this, &CADViewSelector::enable_groups);
    QObject::connect(rm_from_group_button , &QRadioButton::clicked, this, &CADViewSelector::enable_groups);

    wgroups->setLayout(sgrp_gl);
    smode_gl->addWidget(wgroups);

    smode_box->setLayout(smode_gl);

    wl->addWidget(sstyle_box);
    wl->addWidget(select_all_depth_ckbx);
    wl->addWidget(use_ray_test_ckbx);
    wl->addWidget(smode_box );

    this->setLayout(wl);

    pf = new QgSelectPntFilter();
    bf = new QgSelectBoxFilter();
    rf = new QgSelectRayFilter();
    cf = pf;

    QObject::connect(pf, &QgViewFilter::view_updated, this,
	&CADViewSelector::view_changed);
    QObject::connect(bf, &QgViewFilter::view_updated, this,
	&CADViewSelector::view_changed);
    QObject::connect(rf, &QgViewFilter::view_updated, this,
	&CADViewSelector::view_changed);
}

CADViewSelector::~CADViewSelector()
{
	detachFromView(nullptr);
	delete pf;
	delete bf;
	delete rf;
}

void
CADViewSelector::attachToView(QgView *view)
{
    if (!view)
	return;
    if (m_input_endpoint && m_input_view != view)
	detachFromView(m_input_view);

    struct bobol_display_endpoint *endpoint = view->displayEndpoint();
    if (!endpoint || !bobol_display_endpoint_input_action_layer_set(
	endpoint, CADViewSelector::inputActionLayer(), this, this))
	return;
    m_input_view = view;
    m_input_endpoint = endpoint;
    QObject::connect(view, &QObject::destroyed, this, [this, view]() {
	if (m_input_view != view)
	    return;
	m_input_view = nullptr;
	m_input_endpoint = nullptr;
	if (pf)
	    pf->set_view_widget(nullptr);
	if (bf)
	    bf->set_view_widget(nullptr);
	if (rf)
	    rf->set_view_widget(nullptr);
    });
}

void
CADViewSelector::detachFromView(QgView *view)
{
    if (!m_input_endpoint || (view && view != m_input_view))
	return;
    (void)bobol_display_endpoint_input_action_layer_clear_if(
	m_input_endpoint, this);
    if (pf)
	pf->set_view_widget(nullptr);
    if (bf)
	bf->set_view_widget(nullptr);
    if (rf)
	rf->set_view_widget(nullptr);
    m_input_view = nullptr;
    m_input_endpoint = nullptr;
}

void
CADViewSelector::enable_groups(bool)
{
    //current_group->setEnabled(true);
    group_contents->setEnabled(true);
    //add_new_group->setEnabled(true);
    //rm_group->setEnabled(true);
    draw_selections->setEnabled(true);
    erase_selections->setEnabled(true);
}

void
CADViewSelector::disable_groups(bool)
{
    //current_group->setEnabled(false);
    group_contents->setEnabled(false);
    //add_new_group->setEnabled(false);
    //rm_group->setEnabled(false);
    draw_selections->setEnabled(false);
    erase_selections->setEnabled(false);
}

void
CADViewSelector::enable_raytrace_opt(bool)
{
    use_ray_test_ckbx->setEnabled(true);
}

void
CADViewSelector::disable_raytrace_opt(bool)
{
    use_ray_test_ckbx->setEnabled(false);
}

void
CADViewSelector::enable_useall_opt(bool)
{
    select_all_depth_ckbx->setEnabled(true);
}

void
CADViewSelector::disable_useall_opt(bool)
{
    select_all_depth_ckbx->setEnabled(false);
}

void
CADViewSelector::do_view_update(QgViewUpdateFlags flags)
{
    struct ged *gedp = m_ctx ? m_ctx->getGed() : nullptr;
    if (!ged_selection_state_available(gedp))
	return;

    unsigned long long chash = ged_selection_state_hash(gedp, nullptr);
    if ((flags & QG_VIEW_SELECT) || chash != ohash) {
	group_contents->setUpdatesEnabled(false);
	group_contents->clear();
	ohash = chash;

	std::vector<std::string> paths = qged_selection_paths(gedp);
	QStringList items;
	items.reserve(static_cast<int>(paths.size()));
	for (const std::string &path : paths)
	    items.push_back(QString::fromStdString(path));
	group_contents->addItems(items);
	group_contents->setUpdatesEnabled(true);
    }
}

void
CADViewSelector::select_objs()
{
    struct ged *gedp = m_ctx ? m_ctx->getGed() : nullptr;
    if (!ged_selection_state_available(gedp))
	return;

    const std::vector<std::string> &paths = cf->selected_paths();
    if (paths.empty())
	return;

    std::vector<const char *> path_values;
    path_values.reserve(paths.size());
    for (const std::string &path : paths)
	path_values.push_back(path.c_str());
    (void)ged_selection_select_paths(gedp, nullptr, path_values.data(),
	path_values.size(), 0);
}

void
CADViewSelector::deselect_objs()
{
    struct ged *gedp = m_ctx ? m_ctx->getGed() : nullptr;
    if (!ged_selection_state_available(gedp))
	return;

    const std::vector<std::string> &paths = cf->selected_paths();
    if (paths.empty())
	return;

    std::vector<const char *> path_values;
    path_values.reserve(paths.size());
    for (const std::string &path : paths)
	path_values.push_back(path.c_str());
    (void)ged_selection_deselect_paths(gedp, nullptr, path_values.data(),
	path_values.size(), 0);
}


int
CADViewSelector::erase_objs()
{
    struct ged *gedp = m_ctx ? m_ctx->getGed() : nullptr;
    struct ged_view_context *view = qged_selector_view(m_ctx);
    if (!gedp || !view)
	return 0;

    const std::vector<std::string> &paths = cf->selected_paths();
    if (paths.empty())
	return 0;

    int changed = 0;
    return qged_scene_erase_paths(gedp, view, paths, &changed) ==
	GED_SCENE_OK && changed;
}

void
CADViewSelector::do_draw_selections()
{
    struct ged *gedp = m_ctx ? m_ctx->getGed() : nullptr;
    struct ged_view_context *view = qged_selector_view(m_ctx);
    if (!gedp || !view)
	return;

    std::vector<std::string> paths = qged_selection_paths(gedp);
    if (paths.empty())
	return;

    int changed = 0;
    if (qged_scene_draw_paths(gedp, view, paths, &changed) == GED_SCENE_OK &&
	changed)
	emit view_changed(QG_VIEW_DRAWN|QG_VIEW_SELECT);
}

void
CADViewSelector::do_erase_selections()
{
    struct ged *gedp = m_ctx ? m_ctx->getGed() : nullptr;
    struct ged_view_context *view = qged_selector_view(m_ctx);
    if (!gedp || !view)
	return;

    std::vector<std::string> paths = qged_selection_paths(gedp);
    if (paths.empty())
	return;

    int changed = 0;
    if (qged_scene_erase_paths(gedp, view, paths, &changed) ==
	GED_SCENE_OK && changed)
	emit view_changed(QG_VIEW_DRAWN|QG_VIEW_SELECT);
}

int
CADViewSelector::inputActionDispatch(void *user_data,
	BObolInputAction action, const BObolInputEvent *event)
{
    CADViewSelector *selector = static_cast<CADViewSelector *>(user_data);
    return selector ? selector->applyInputAction(action, event) :
	BOBOL_INPUT_RESULT_ERROR;
}

int
CADViewSelector::applyInputAction(BObolInputAction action,
	const BObolInputEvent *event)
{
    if (!event || !m_input_view ||
	(action != QG_SELECT_INPUT_BEGIN &&
	 action != QG_SELECT_INPUT_UPDATE &&
	 action != QG_SELECT_INPUT_COMMIT &&
	 action != QG_SELECT_INPUT_CANCEL))
	return BOBOL_INPUT_RESULT_UNHANDLED;

    struct ged *gedp = m_ctx ? m_ctx->getGed() : nullptr;
    if (!gedp || !m_input_view->viewContext())
	return BOBOL_INPUT_RESULT_UNHANDLED;

    cf = pf;
    if (use_rect_select_button->isChecked())
	cf = bf;
    if (use_ray_test_ckbx->isChecked()) {
	rf->dbip = gedp->dbip;
	cf = rf;
    }
    /* A rectangle is a set-selection gesture: with its depth option disabled,
     * it must collect all paths in the box rather than silently degrade to
     * the single nearest hit.  Point/ray selection retains the explicit
     * all-intersections choice. */
    cf->first_only = use_rect_select_button->isChecked() ? false :
	(select_all_depth_ckbx->isChecked() ? false : true);
    cf->set_view_widget(m_input_view);
    if (!cf->semanticInput(action, event))
	return BOBOL_INPUT_RESULT_UNHANDLED;

    /* Gesture updates only draw the selection affordance.  Applying an
     * operation before commit can act on a stale result from an earlier
     * gesture. */
    if (action != QG_SELECT_INPUT_COMMIT)
	return action == QG_SELECT_INPUT_CANCEL ?
	    BOBOL_INPUT_RESULT_CANCELLED : BOBOL_INPUT_RESULT_HANDLED;

    if (erase_from_scene_button->isChecked()) {
	if (erase_objs())
	    emit view_changed(QG_VIEW_DRAWN);
	return BOBOL_INPUT_RESULT_HANDLED;
    }

    if (add_to_group_button->isChecked()) {
	select_objs();
	emit view_changed(QG_VIEW_DRAWN|QG_VIEW_SELECT);
	return BOBOL_INPUT_RESULT_HANDLED;
    }

    if (rm_from_group_button->isChecked()) {
	deselect_objs();
	emit view_changed(QG_VIEW_DRAWN|QG_VIEW_SELECT);
	return BOBOL_INPUT_RESULT_HANDLED;
    }

    return BOBOL_INPUT_RESULT_HANDLED;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
