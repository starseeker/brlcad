/*                          Q B O T . H
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
/** @file QBot.h
 *
 * Bag-of-triangles (BOT) editor using the neutral GED edit transaction and
 * Obol primitive-pick details.
 *
 * Workflow
 * --------
 * 1. `_bot_edit` carries the wire preview and `_bot_edit_surface` carries
 *    pickable face geometry.
 * 2. GED edit transactions own preview lifecycle and revision updates.
 * 3. Obol pick details identify faces and nearest vertices for edit handles.
 * 4. Vertex, edge, and face drags update an owned BoT copy; Apply commits a
 *    duplicate through the database event batch and Reset restores the DB.
 * 5. Preview teardown is owned by neutral view features.
 */

#ifndef QBOT_H
#define QBOT_H

#include <QCheckBox>
#include <QComboBox>
#include <QLineEdit>
#include <QPushButton>
#include <QGroupBox>
#include <vector>
#include "raytrace.h"
#include "qtcad/QgTypes.h"
#include "../qged_edit_preview_util.h"

class QgPluginContext;

class QBot : public QWidget
{
    Q_OBJECT

    public:
	QBot();
	~QBot();

	void setContext(QgPluginContext *ctx) { m_ctx = ctx; }

	/* Primitive name */
	QLineEdit *bot_name;
	/* Mode selector (vertex / face / edge) */
	QComboBox *edit_mode;
	QPushButton *write_edit;
	QPushButton *reset_values;

    signals:
	void view_updated(QgViewUpdateFlags);

    private slots:
	void read_from_db();
	void write_to_db();
	void update_obj_wireframe();
	void update_viewobj_name(const QString &);

    protected:
	bool eventFilter(QObject *, QEvent *);

    private:
	struct directory *dp = NULL;
	struct rt_bot_internal *bot = NULL;
	std::vector<int> selected_vertices;
	std::vector<fastf_t> drag_vertex_positions;
	point_t drag_start = VINIT_ZERO;
	bool dragging = false;
	bool dirty = false;
	/* Edit-preview overlay feature. */
	qged_edit_feature_ref p = QGED_EDIT_FEATURE_REF_NULL;
	struct bu_vls oname = BU_VLS_INIT_ZERO;
	QgPluginContext *m_ctx = nullptr;

	struct ged *getGed() const;
	void clear_edit_state();
	bool load_bot();
	void publish_selection_handle();
};

#endif /* QBOT_H */

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
