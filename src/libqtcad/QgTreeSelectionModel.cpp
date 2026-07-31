/*            Q G T R E E S E L E C T I O N M O D E L . C P P
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
/** @file QgTreeSelectionModel.cpp
 *
 */

#include "common.h"
#include <algorithm>
#include <queue>
#include <unordered_set>
#include <stack>
#include <vector>
#include <QGuiApplication>
#include "qtcad/QgUtil.h"
#include "qtcad/QgModel.h"
#include "qtcad/QgTreeSelectionModel.h"
#include "qtcad/QgSignalFlags.h"

#include "ged/selection_state.h"

static void
qg_commit_selection_change(QgTreeView *treeview, QgModel *model)
{
	if (!treeview || !model || !model->ged())
		return;

	/* All path mutations in one Qt selection request are complete before
	 * deriving ancestor/descendant metadata or touching retained
	 * presentation.  Large range selections therefore pay for exactly one
	 * hierarchy recompute, one sparse draw sync, and one indexed row
	 * notification pass. */
	ged_selection_recompute(model->ged(), nullptr);

	QgViewUpdateFlags flags = QG_VIEW_SELECT;
	if (ged_selection_draw_sync(model->ged(), nullptr))
		flags |= QG_VIEW_REFRESH;

	emit treeview->view_changed(flags);
	model->notifySelectionItemsChanged();
}

void
QgTreeSelectionModel::clear_all()
{
	QgModel *m = treeview->cadModel();

	ged_selection_clear(m->ged(), nullptr);
}

void
QgTreeSelectionModel::select(const QItemSelection &selection, QItemSelectionModel::SelectionFlags flags)
{
	QTCAD_SLOT("QgTreeSelectionModel::select QItemSelection", 1);
	QgModel *m = treeview->cadModel();
	/* QItemSelectionModel adjusts its private ranges while the source model is
	 * being reset or rows are changing.  Those maintenance calls are not CAD
	 * selection requests and must not clear or redraw the GED scene. */
	if (m->structureChangeInProgress()) {
		QItemSelectionModel::select(selection, flags);
		return;
	}
	struct ged *gedp = m->ged();

	QModelIndexList dl = selection.indexes();
	std::unordered_set<QgItem *> visited;
	for (long int i = 0; i < dl.size(); i++) {
		QgItem *snode = static_cast<QgItem *>(dl.at(i).internalPointer());
		/* QItemSelection::indexes() can contain one index per column for
		 * the same tree item.  A CAD path is a row identity, so never
		 * repeat its hash lookup and mutation. */
		if (!snode || !visited.insert(snode).second)
			continue;

		// If we are selecting an already selected node, clear it
		if (flags & QItemSelectionModel::Select &&
			ged_selection_is_path_selected(gedp, nullptr, snode->path_hash())) {
			if (!(QGuiApplication::keyboardModifiers().testFlag(Qt::ShiftModifier))) {
				if (flags & QItemSelectionModel::Clear &&
					ged_selection_count(gedp, nullptr) > 1) {
					ged_selection_clear(gedp, nullptr);
					std::vector<unsigned long long> path_hashes = snode->path_items();
					ged_selection_select_path_ids(gedp, nullptr,
						path_hashes.data(), path_hashes.size(), 0);
				}
				else {
					std::vector<unsigned long long> path_hashes = snode->path_items();
					ged_selection_deselect_path_ids(gedp, nullptr,
						path_hashes.data(), path_hashes.size(), 0);
				}
			}
		}
		else {
			if (flags & QItemSelectionModel::Clear)
				ged_selection_clear(gedp, nullptr);

			std::vector<unsigned long long> path_hashes = snode->path_items();
			ged_selection_select_path_ids(gedp, nullptr,
				path_hashes.data(), path_hashes.size(), 0);
		}
	}

	qg_commit_selection_change(treeview, m);
}

void
QgTreeSelectionModel::select(const QModelIndex &index, QItemSelectionModel::SelectionFlags flags)
{
	QTCAD_SLOT("QgTreeSelectionModel::select QModelIndex", 1);
	QgModel *m = treeview->cadModel();
	/* See the QItemSelection overload above. */
	if (m->structureChangeInProgress()) {
		QItemSelectionModel::select(index, flags);
		return;
	}
	struct ged *gedp = m->ged();

	if (flags & QItemSelectionModel::Clear)
		ged_selection_clear(gedp, nullptr);

	QgItem *snode = static_cast<QgItem *>(index.internalPointer());
	if (!snode) {
		ged_selection_clear(gedp, nullptr);
		qg_commit_selection_change(treeview, m);
		return;
	}

	if (!(flags & QItemSelectionModel::Deselect)) {

		// If we are selecting an already selected node, clear it
		if (flags & QItemSelectionModel::Select &&
			ged_selection_is_path_selected(gedp, nullptr, snode->path_hash())) {
			std::vector<unsigned long long> path_hashes = snode->path_items();
			ged_selection_deselect_path_ids(gedp, nullptr,
				path_hashes.data(), path_hashes.size(), 0);
			qg_commit_selection_change(treeview, m);
			return;
		}

		std::vector<unsigned long long> path_hashes = snode->path_items();
		ged_selection_select_path_ids(gedp, nullptr,
			path_hashes.data(), path_hashes.size(), 0);

	}
	else {

		std::vector<unsigned long long> path_hashes = snode->path_items();
		ged_selection_deselect_path_ids(gedp, nullptr,
			path_hashes.data(), path_hashes.size(), 0);

	}

	qg_commit_selection_change(treeview, m);
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
