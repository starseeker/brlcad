/*                       Q B O T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
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
/** @file QBot.cpp
 *
 * Edit-preview BOT editor.
 *
 * The transient preview and pickable surface route BoT edits through the
 * neutral GED edit transaction and view-feature APIs.
 */

#include "common.h"

#include "bv.h"
#include <QEvent>
#include <QMouseEvent>
#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <QGroupBox>
#include <QVBoxLayout>
#include "ged.h"
#include "ged/draw.h"
#include "ged/view_feature_batch.h"
#include "rt/db_io.h"
#include "rt/directory.h"
#include "rt/primitives/bot.h"
#include "qtcad/QgGedEventBatch.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgSignalFlags.h"
#include "../qged_edit_preview_util.h"
#include "QBot.h"

#include <algorithm>

/* ---- QBot constructor --------------------------------------------------- */

QBot::QBot()
    : QWidget()
{
    QVBoxLayout *l = new QVBoxLayout;

    QLabel *name_label = new QLabel("Object name:");
    l->addWidget(name_label);
    bot_name = new QLineEdit();
    l->addWidget(bot_name);

    QGroupBox *mode_box = new QGroupBox("Edit mode");
    QVBoxLayout *mbl = new QVBoxLayout;
    edit_mode = new QComboBox();
    edit_mode->addItem("Vertex");
    edit_mode->addItem("Face");
    edit_mode->addItem("Edge");
    mbl->addWidget(edit_mode);
    mode_box->setLayout(mbl);
    l->addWidget(mode_box);

    QGroupBox *ac_box = new QGroupBox("Actions");
    QVBoxLayout *acl = new QVBoxLayout;
    write_edit = new QPushButton("Apply");
    acl->addWidget(write_edit);
    reset_values = new QPushButton("Reset");
    acl->addWidget(reset_values);
    ac_box->setLayout(acl);
    l->addWidget(ac_box);

    l->setAlignment(Qt::AlignTop);
    this->setLayout(l);

    QObject::connect(bot_name, &QLineEdit::textChanged,
		     this, &QBot::update_viewobj_name);
    QObject::connect(write_edit, &QPushButton::clicked,
		     this, &QBot::write_to_db);
    QObject::connect(reset_values, &QPushButton::clicked,
		     this, &QBot::read_from_db);
}

QBot::~QBot()
{
    if (!qged_edit_feature_ref_is_null(p))
	qged_edit_preview_publish_event(m_ctx, p, "_bot_edit",
	    QGED_EDIT_PREVIEW_CANCEL, bu_vls_cstr(&oname));
    qged_edit_feature_clear_geometry(p);
    if (!qged_edit_feature_ref_is_null(p) && m_ctx) {
	qged_edit_feature_remove(m_ctx, "_bot_edit");
	(void)ged_view_feature_remove(qged_edit_ged_view_context(m_ctx),
	    "_bot_edit_surface");
	(void)ged_view_feature_remove(qged_edit_ged_view_context(m_ctx),
	    "_bot_edit_handle");
	p = QGED_EDIT_FEATURE_REF_NULL;
    }
    clear_edit_state();
    bu_vls_free(&oname);
}

struct ged *
QBot::getGed() const
{
    if (!m_ctx)
	return nullptr;
    return m_ctx->getGed();
}

void
QBot::clear_edit_state()
{
    if (bot) {
	struct rt_db_internal intern = RT_DB_INTERNAL_INIT_ZERO;
	intern.idb_major_type = DB5_MAJORTYPE_BRLCAD;
	intern.idb_type = ID_BOT;
	intern.idb_minor_type = DB5_MINORTYPE_BRLCAD_BOT;
	intern.idb_ptr = bot;
	intern.idb_meth = &OBJ[ID_BOT];
	rt_db_free_internal(&intern);
	bot = NULL;
    }
    selected_vertices.clear();
    drag_vertex_positions.clear();
    dragging = false;
    dirty = false;
}

bool
QBot::load_bot()
{
    struct ged *gedp = getGed();
    if (!gedp || !gedp->dbip || !bu_vls_strlen(&oname))
	return false;

    struct directory *next_dp = db_lookup(gedp->dbip,
	bu_vls_cstr(&oname), LOOKUP_QUIET);
    if (!next_dp || next_dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BOT)
	return false;

    struct rt_db_internal intern = RT_DB_INTERNAL_INIT_ZERO;
    if (rt_db_get_internal(&intern, next_dp, gedp->dbip, NULL) < 0)
	return false;
    if (intern.idb_minor_type != DB5_MINORTYPE_BRLCAD_BOT ||
	!intern.idb_ptr) {
	rt_db_free_internal(&intern);
	return false;
    }

    struct rt_bot_internal *next_bot = rt_bot_dup(
	(const struct rt_bot_internal *)intern.idb_ptr);
    rt_db_free_internal(&intern);
    if (!next_bot)
	return false;

    clear_edit_state();
    bot = next_bot;
    dp = next_dp;
    return true;
}

void
QBot::read_from_db()
{
    struct ged *gedp = getGed();
    if (!gedp || !gedp->dbip || !bu_vls_strlen(&oname))
	return;


    qged_edit_preview_publish_event(m_ctx, p, "_bot_edit",
	QGED_EDIT_PREVIEW_CANCEL, bu_vls_cstr(&oname));
    p = QGED_EDIT_FEATURE_REF_NULL;
    if (!load_bot())
	return;
    update_obj_wireframe();
    emit view_updated(QG_VIEW_REFRESH);
}

void
QBot::write_to_db()
{
    struct ged *gedp = getGed();
    if (!gedp || !gedp->dbip || !dp || !bot)
	return;

    struct rt_bot_internal *write_bot = rt_bot_dup(bot);
    if (!write_bot)
	return;
    struct rt_db_internal intern = RT_DB_INTERNAL_INIT_ZERO;
    intern.idb_major_type = DB5_MAJORTYPE_BRLCAD;
    intern.idb_type = ID_BOT;
    intern.idb_minor_type = DB5_MINORTYPE_BRLCAD_BOT;
    intern.idb_ptr = write_bot;
    intern.idb_meth = &OBJ[ID_BOT];
    {
	QgGedEventBatch event_batch(gedp);
	if (rt_db_put_internal(dp, gedp->dbip, &intern) < 0) {
	    rt_db_free_internal(&intern);
	    return;
	}
	(void)ged_event_notify_object_modified(gedp, bu_vls_cstr(&oname),
	    1, NULL);
    }
    dirty = false;
    qged_edit_preview_publish_event(m_ctx, p, "_bot_edit",
	    QGED_EDIT_PREVIEW_COMMIT,
	    bu_vls_cstr(&oname));
    p = QGED_EDIT_FEATURE_REF_NULL;
    emit view_updated(QG_VIEW_DB);
}

void
QBot::update_obj_wireframe()
{
    struct ged *gedp = getGed();
    if (!gedp)
	return;

    p = qged_edit_feature_overlay_ensure(m_ctx, "_bot_edit",
	    bu_vls_cstr(&oname));
    if (qged_edit_feature_ref_is_null(p))
	return;

    if (!gedp->dbip || !bu_vls_strlen(&oname)) {
	qged_edit_feature_clear_geometry(p);
	qged_edit_feature_set_visible(p, 0);
	(void)ged_view_feature_remove(qged_edit_ged_view_context(m_ctx),
	    "_bot_edit_surface");
	(void)ged_view_feature_remove(qged_edit_ged_view_context(m_ctx),
	    "_bot_edit_handle");
	return;
    }

    dp = db_lookup(gedp->dbip, bu_vls_cstr(&oname), LOOKUP_QUIET);
    if (!dp || dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BOT) {
	qged_edit_feature_clear_geometry(p);
	qged_edit_feature_set_visible(p, 0);
	(void)ged_view_feature_remove(qged_edit_ged_view_context(m_ctx),
	    "_bot_edit_surface");
	(void)ged_view_feature_remove(qged_edit_ged_view_context(m_ctx),
	    "_bot_edit_handle");
	return;
    }

    if (!bot && !load_bot())
	return;

    qged_edit_feature_clear_geometry(p);
    qged_edit_feature_set_visible(p, 1);

    if (qged_edit_feature_replace_bot_face_lines(p,
	    QGED_EDIT_FEATURE_TRANSIENT_PREVIEW, bot) &&
	    !qged_edit_feature_ref_is_null(p)) {
	qged_edit_preview_publish_event(m_ctx, p, "_bot_edit",
		QGED_EDIT_PREVIEW_UPDATE,
		bu_vls_cstr(&oname));
    }

    struct ged_view_context *view_ctx = qged_edit_ged_view_context(m_ctx);
    if (view_ctx) {
	point_t *surface_points = (point_t *)bu_calloc(bot->num_vertices,
		sizeof(point_t), "BOT edit surface points");
	for (size_t i = 0; i < bot->num_vertices; i++)
	    VMOVE(surface_points[i], &bot->vertices[i * 3]);
	struct ged_view_feature_style surface_style =
	    ged_view_feature_style_default();
	surface_style.visible = 1;
	surface_style.selectable = 1;
	surface_style.color_valid = 1;
	surface_style.color[0] = 72;
	surface_style.color[1] = 126;
	surface_style.color[2] = 168;
	struct ged_view_feature_batch_desc desc =
	    ged_view_feature_batch_desc_default();
	desc.owner_id = "qged-bot-edit";
	desc.owner_role = "edit-handle";
	desc.overlay_class = GED_VIEW_FEATURE_OVERLAY_CLASS_EDIT_HANDLE;
	desc.lifecycle = GED_VIEW_FEATURE_LIFECYCLE_PER_TOOL;
	desc.local = 1;
	struct ged_view_feature_batch *batch =
	    ged_view_feature_batch_begin(view_ctx, &desc);
	if (batch && ged_view_feature_batch_indexed_face_set_replace(batch,
		"_bot_edit_surface", (const point_t *)surface_points,
		bot->num_vertices, NULL, 0, bot->faces,
		bot->num_faces * 3, &surface_style))
	    (void)ged_view_feature_batch_commit(batch);
	else if (batch)
	    ged_view_feature_batch_abort(batch);
	bu_free(surface_points, "BOT edit surface points");
    }

    const char *wcolor = "255/255/255";
    const char *av[2] = {wcolor, NULL};
    struct bu_color cval;
    bu_opt_color(NULL, 1, (const char **)&av[0], (void *)&cval);
    unsigned char rgb[3] = {0, 0, 0};
    bu_color_to_rgb_chars(&cval, rgb);
    qged_edit_feature_set_color(p, rgb[0], rgb[1], rgb[2]);

    qged_edit_feature_touch(p);
    publish_selection_handle();
}

void
QBot::update_viewobj_name(const QString &ostr)
{
    bu_vls_sprintf(&oname, "%s", ostr.toLocal8Bit().data());
    if (!load_bot()) {
	clear_edit_state();
	update_obj_wireframe();
	return;
    }
    update_obj_wireframe();
    emit view_updated(QG_VIEW_REFRESH);
}

void
QBot::publish_selection_handle()
{
    struct ged_view_context *view_ctx = qged_edit_ged_view_context(m_ctx);
    if (!view_ctx || !bot || selected_vertices.empty()) {
	if (view_ctx)
	    (void)ged_view_feature_remove(view_ctx,
		"_bot_edit_handle");
	return;
    }

    const size_t count = selected_vertices.size();
    point_t *points = (point_t *)bu_calloc(count, sizeof(point_t),
	"BOT edit selection handles");
    int *commands = (int *)bu_calloc(count, sizeof(int),
	"BOT edit selection handle commands");
    for (size_t i = 0; i < count; i++) {
	const int vertex = selected_vertices[i];
	if (vertex >= 0 && (size_t)vertex < bot->num_vertices)
	    VMOVE(points[i], &bot->vertices[(size_t)vertex * 3]);
	commands[i] = GED_DRAW_VIEW_LINE_POINT_DRAW;
    }

    struct ged_view_feature_style style =
	ged_view_feature_style_default();
    style.visible = 1;
    style.selectable = 1;
    style.color_valid = 1;
    style.color[0] = 255;
    style.color[1] = 196;
    style.color[2] = 32;
    struct ged_view_feature_batch_desc desc =
	ged_view_feature_batch_desc_default();
    desc.owner_id = "qged-bot-edit";
    desc.owner_role = "edit-handle";
    desc.overlay_class = GED_VIEW_FEATURE_OVERLAY_CLASS_EDIT_HANDLE;
    desc.lifecycle = GED_VIEW_FEATURE_LIFECYCLE_PER_TOOL;
    desc.local = 1;
    struct ged_view_feature_batch *batch =
	ged_view_feature_batch_begin(view_ctx, &desc);
    if (batch && ged_view_feature_batch_line_set_replace(batch,
	    "_bot_edit_handle", (const point_t *)points, commands, count,
	    &style))
	(void)ged_view_feature_batch_commit(batch);
    else if (batch)
	ged_view_feature_batch_abort(batch);
    bu_free(commands, "BOT edit selection handle commands");
    bu_free(points, "BOT edit selection handles");
}

static fastf_t
bot_point_segment_distance_sq(const point_t point, const fastf_t *a,
	const fastf_t *b)
{
    vect_t ab;
    vect_t ap;
    VSUB2(ab, b, a);
    VSUB2(ap, point, a);
    const fastf_t length_sq = MAGSQ(ab);
    fastf_t t = length_sq > SMALL_FASTF ? VDOT(ap, ab) / length_sq : 0.0;
    t = std::max((fastf_t)0.0, std::min((fastf_t)1.0, t));
    point_t nearest;
    VJOIN1(nearest, a, t, ab);
    vect_t delta;
    VSUB2(delta, point, nearest);
    return MAGSQ(delta);
}

bool
QBot::eventFilter(QObject *, QEvent *event)
{
    if (!event || !m_ctx || !bot)
	return false;
    QMouseEvent *mouse = static_cast<QMouseEvent *>(event);

    if (event->type() == QEvent::MouseMove && dragging) {
	point_t current;
	const struct bv *view = bv_context_view_const(
	    static_cast<const struct bv_context *>(qged_edit_view_context(m_ctx)));
	if (!view || !bv_screen_to_model(current, view,
		mouse->pos().x(), mouse->pos().y()))
	    return true;
	vect_t delta;
	VSUB2(delta, current, drag_start);
	for (size_t i = 0; i < selected_vertices.size(); i++) {
	    const int vertex = selected_vertices[i];
	    if (vertex < 0 || (size_t)vertex >= bot->num_vertices)
		continue;
	    fastf_t *dst = &bot->vertices[(size_t)vertex * 3];
	    const fastf_t *src = &drag_vertex_positions[i * 3];
	    VADD2(dst, src, delta);
	}
	dirty = true;
	update_obj_wireframe();
	emit view_updated(QG_VIEW_REFRESH);
	return true;
    }

    if (event->type() == QEvent::MouseButtonRelease && dragging &&
	mouse->button() == Qt::LeftButton) {
	dragging = false;
	qged_edit_preview_publish_event(m_ctx, p, "_bot_edit",
	    QGED_EDIT_PREVIEW_UPDATE, bu_vls_cstr(&oname));
	emit view_updated(QG_VIEW_REFRESH);
	return true;
    }

    if (event->type() != QEvent::MouseButtonPress ||
	mouse->button() != Qt::LeftButton)
	return false;

    struct ged_view_context *view_ctx = qged_edit_ged_view_context(m_ctx);
    if (!view_ctx)
	return false;
    struct ged_pick_result *result =
	ged_pick_nearest(view_ctx,
	    mouse->pos().x(), mouse->pos().y());
    struct ged_pick_detail detail = ged_pick_detail_default();
    struct bu_vls path = BU_VLS_INIT_ZERO;
    const int picked = result && ged_pick_result_count(result) > 0 &&
	ged_pick_result_path(result, 0, &path) &&
	ged_pick_result_detail(result, 0, &detail) &&
	BU_STR_EQUAL(bu_vls_cstr(&path), "_bot_edit_surface") &&
	detail.primitive_kind == 3 && detail.primitive_index >= 0;
    if (!picked) {
	bu_vls_free(&path);
	ged_pick_result_free(result);
	return false;
    }

    const int face = detail.primitive_index;
    if (face < 0 || (size_t)face >= bot->num_faces) {
	bu_vls_free(&path);
	ged_pick_result_free(result);
	return false;
    }
    (void)ged_view_feature_set_selection(
	view_ctx, "_bot_edit_surface", &face, 1);

    selected_vertices.clear();
    if (edit_mode->currentIndex() == 0) {
	int vertex = detail.nearest_face_vertex_index;
	if (vertex < 0)
	    vertex = bot->faces[(size_t)face * 3];
	selected_vertices.push_back(vertex);
    } else if (edit_mode->currentIndex() == 1) {
	for (int i = 0; i < 3; i++)
	    selected_vertices.push_back(bot->faces[(size_t)face * 3 + i]);
    } else {
	int edge = 0;
	if (detail.model_point_valid) {
	    fastf_t best = INFINITY;
	    for (int i = 0; i < 3; i++) {
		const int a = bot->faces[(size_t)face * 3 + i];
		const int b = bot->faces[(size_t)face * 3 + ((i + 1) % 3)];
		const fastf_t distance = bot_point_segment_distance_sq(
		    detail.model_point, &bot->vertices[(size_t)a * 3],
		    &bot->vertices[(size_t)b * 3]);
		if (distance < best) {
		    best = distance;
		    edge = i;
		}
	    }
	}
	selected_vertices.push_back(bot->faces[(size_t)face * 3 + edge]);
	selected_vertices.push_back(
	    bot->faces[(size_t)face * 3 + ((edge + 1) % 3)]);
    }

    std::sort(selected_vertices.begin(), selected_vertices.end());
    selected_vertices.erase(std::unique(selected_vertices.begin(),
	selected_vertices.end()), selected_vertices.end());
    drag_vertex_positions.clear();
    drag_vertex_positions.reserve(selected_vertices.size() * 3);
    for (const int vertex : selected_vertices) {
	const fastf_t *point = &bot->vertices[(size_t)vertex * 3];
	drag_vertex_positions.insert(drag_vertex_positions.end(), point,
	    point + 3);
    }

    const struct bv *view = bv_context_view_const(
	qged_edit_view_context(m_ctx));
    if (!view || !bv_screen_to_model(drag_start, view,
	    mouse->pos().x(), mouse->pos().y())) {
	const int vertex = selected_vertices.front();
	VMOVE(drag_start, &bot->vertices[(size_t)vertex * 3]);
    }
    dragging = true;
    publish_selection_handle();
    qged_edit_preview_publish_event(m_ctx, p, "_bot_edit",
	QGED_EDIT_PREVIEW_BEGIN, bu_vls_cstr(&oname));

    bu_vls_free(&path);
    ged_pick_result_free(result);
    emit view_updated(QG_VIEW_REFRESH);
    return true;
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
