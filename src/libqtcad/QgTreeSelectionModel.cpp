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
#include <unordered_set>
#include <QGuiApplication>
#include "qtcad/QgModel.h"
#include "qtcad/QgTreeSelectionModel.h"
#include "qtcad/QgTreeView.h"

static void
qg_commit_selection_change(QgTreeView *tree_view, QgModel *model)
{
	if (!tree_view || !model || !model->ged())
		return;

	/* All path mutations in one Qt selection request are complete before
	 * deriving ancestor/descendant metadata or touching retained
	 * presentation.  Large range selections therefore pay for exactly one
	 * hierarchy recompute, one sparse selection commit, and one indexed row
	 * notification pass. */
	QgViewUpdateFlags flags = model->selectionCommit();

	emit tree_view->view_changed(flags);
}

void
QgTreeSelectionModel::clear_all()
{
	QgModel *model = m_tree_view->cadModel();

	model->selectionClear();
	qg_commit_selection_change(m_tree_view, model);
}

void
QgTreeSelectionModel::select(const QItemSelection &selection, QItemSelectionModel::SelectionFlags flags)
{
	QTCAD_SLOT("QgTreeSelectionModel::select QItemSelection", 1);
	QgModel *model = m_tree_view->cadModel();
	/* QItemSelectionModel adjusts its private ranges while the source model is
	 * being reset or rows are changing.  Those maintenance calls are not CAD
	 * selection requests and must not clear or redraw the GED scene. */
	if (model->structureChangeInProgress()) {
		QItemSelectionModel::select(selection, flags);
		return;
	}
	QModelIndexList indexes = selection.indexes();
	std::unordered_set<QgItem *> visited;
	for (const QModelIndex &index : indexes) {
		QgItem *snode = static_cast<QgItem *>(index.internalPointer());
		/* QItemSelection::indexes() can contain one index per column for
		 * the same tree item.  A CAD path is a row identity, so never
		 * repeat its hash lookup and mutation. */
		if (!snode || !visited.insert(snode).second)
			continue;

		// If we are selecting an already selected node, clear it
		if (flags & QItemSelectionModel::Select &&
			model->itemSelected(snode)) {
			if (!(QGuiApplication::keyboardModifiers().testFlag(Qt::ShiftModifier))) {
				if (flags & QItemSelectionModel::Clear &&
					model->selectionCount() > 1) {
					model->selectionClear();
					model->itemSelect(snode);
				}
				else
					model->itemDeselect(snode);
			}
		}
		else {
			if (flags & QItemSelectionModel::Clear)
				model->selectionClear();

			model->itemSelect(snode);
		}
	}

	qg_commit_selection_change(m_tree_view, model);
}

void
QgTreeSelectionModel::select(const QModelIndex &index, QItemSelectionModel::SelectionFlags flags)
{
	QTCAD_SLOT("QgTreeSelectionModel::select QModelIndex", 1);
	QgModel *model = m_tree_view->cadModel();
	/* See the QItemSelection overload above. */
	if (model->structureChangeInProgress()) {
		QItemSelectionModel::select(index, flags);
		return;
	}
	if (flags & QItemSelectionModel::Clear)
		model->selectionClear();

	QgItem *snode = static_cast<QgItem *>(index.internalPointer());
	if (!snode) {
		model->selectionClear();
		qg_commit_selection_change(m_tree_view, model);
		return;
	}

	if (!(flags & QItemSelectionModel::Deselect)) {

		// If we are selecting an already selected node, clear it
		if (flags & QItemSelectionModel::Select &&
			model->itemSelected(snode)) {
			model->itemDeselect(snode);
			qg_commit_selection_change(m_tree_view, model);
			return;
		}

		model->itemSelect(snode);

	}
	else {

		model->itemDeselect(snode);

	}

	qg_commit_selection_change(m_tree_view, model);
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
