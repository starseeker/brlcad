/*                    Q G B R E P F I L T E R . H
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this file; see the file named COPYING for more information.
 */
/** @file QgBrepFilter.h
 *
 * Reusable Qt view filters for topology-aware B-rep NURBS CV editing.
 */

#ifndef QGBREPFILTER_H
#define QGBREPFILTER_H

#include "common.h"

extern "C" {
#include "bu/ptbl.h"
#include "rt/edit.h"
}

#include <QEvent>
#include <QString>

#include "qtcad/defines.h"
#include "qtcad/QgViewFilter.h"

/**
 * Common B-rep editing filter support.
 *
 * The caller owns es and must keep it valid while a filter is installed.
 * Picking can select both editable and locked CVs so the UI can report the
 * applicable B-rep constraint policy.
 */
class QTCAD_EXPORT QgBrepFilter : public QgViewFilter {
	Q_OBJECT
	Q_DISABLE_COPY_MOVE(QgBrepFilter)

public:
	explicit QgBrepFilter(QObject *parent = nullptr);

	/** Select the closest projected surface CV within max_px pixels. */
	bool pick_cv(int sx, int sy, fastf_t max_px = 12.0);

	/** Read the current rt_edit B-rep CV selection. */
	bool selected_cv(int *face, int *cv_u, int *cv_v) const;

	/** Return whether the current selection can move without changing a trim. */
	bool selected_cv_is_topology_safe() const;

	/** Return whether the current selection has a supported translation path. */
	bool selected_cv_can_translate() const;

	/** Convert a screen-space motion into primitive-local model units. */
	bool screen_delta_to_local(int from_x, int from_y, int to_x, int to_y,
		vect_t delta) const;

	struct rt_edit *es = nullptr;

signals:
	void brep_selection_changed(int face, int cv_u, int cv_v,
		bool topology_safe, bool can_translate);
	void brep_changed();
	void edit_rejected(const QString &reason);
};

/** Left-click CV selection by projected screen-space proximity. */
class QTCAD_EXPORT QgBrepPickCVFilter : public QgBrepFilter {
	Q_OBJECT
	Q_DISABLE_COPY_MOVE(QgBrepPickCVFilter)

public:
	using QgBrepFilter::QgBrepFilter;
	bool eventFilter(QObject *, QEvent *e) override;

	fastf_t pick_radius_px = 12.0;
};

/**
 * Left-drag the selected CV in the current view plane.
 *
 * If no CV is selected, the press first picks the nearest one.  Interior CVs
 * and eligible C0-coupled isoparametric boundary CVs may be dragged; other
 * trim-influencing CVs remain selectable but locked.
 */
class QTCAD_EXPORT QgBrepMoveCVFilter : public QgBrepFilter {
	Q_OBJECT
	Q_DISABLE_COPY_MOVE(QgBrepMoveCVFilter)

public:
	using QgBrepFilter::QgBrepFilter;
	bool eventFilter(QObject *, QEvent *e) override;

	fastf_t pick_radius_px = 12.0;
	bool pick_on_press = true;

private:
	bool m_dragging = false;
	int m_prev_x = 0;
	int m_prev_y = 0;
};

#endif /* QGBREPFILTER_H */

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
