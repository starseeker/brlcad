/*                    Q G P O L Y F I L T E R . H
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
/** @file QgPolyFilter.h
 *
 *  Qt mouse filters for polygon editing modes
 */

#ifndef QGPOLYFILTER_H
#define QGPOLYFILTER_H

#include "common.h"

extern "C" {
#include "bg/polygon.h"
#include "bu/color.h"
#include "ged/draw.h"
}

#include <vector>
#include <QBoxLayout>
#include <QEvent>
#include <QMouseEvent>
#include <QWidget>
#include "qtcad/defines.h"
#include "qtcad/QgViewFilter.h"

enum qg_polygon_type {
	QG_POLYGON_GENERAL = 0,
	QG_POLYGON_CIRCLE = 1,
	QG_POLYGON_ELLIPSE = 2,
	QG_POLYGON_RECTANGLE = 3,
	QG_POLYGON_SQUARE = 4
};

enum qg_polygon_update_mode {
	QG_POLYGON_UPDATE_DEFAULT = 0,
	QG_POLYGON_UPDATE_PROPS_ONLY = 1,
	QG_POLYGON_UPDATE_PT_SELECT = 2,
	QG_POLYGON_UPDATE_PT_SELECT_CLEAR = 3,
	QG_POLYGON_UPDATE_PT_MOVE = 4,
	QG_POLYGON_UPDATE_PT_APPEND = 5
};

typedef ged_draw_view_polygon_ref qg_polygon_ref;
#define QG_POLYGON_REF_NULL_INIT GED_DRAW_VIEW_POLYGON_REF_NULL_INIT
typedef struct ged_draw_view_polygon_record qg_polygon_record;

// Filters designed for specific editing modes
class QTCAD_EXPORT QgPolyFilter : public QgViewFilter {
	Q_OBJECT
	Q_DISABLE_COPY_MOVE(QgPolyFilter)


    public:
	QgPolyFilter() = default;

	// We want to be able to swap derived QgPolyFilter classes in
	// parent calling code - make eventFilter virtual to help
	// simplify doing so.
	bool eventFilter(QObject *, QEvent *) override
	{
		return false;
	};

signals:
	void finalized(bool);

public:
	bool close_polygon();

	bg_clip_t op = bg_None;
	qg_polygon_ref polygon = QG_POLYGON_REF_NULL_INIT;
	int ptype = QG_POLYGON_CIRCLE;
	bool close_general_poly = true; // set to false if application wants to allow non-closed polygons
	struct bu_color fill_color = BU_COLOR_BLUE;
	struct bu_color edge_color = BU_COLOR_YELLOW;
	bool fill_poly = false;
	double fill_slope_x = 1.0;
	double fill_slope_y = 1.0;
	double fill_density = 10.0;
	double vZ = 0.0;
	std::string vname;
};

class QTCAD_EXPORT QgPolyCreateFilter : public QgPolyFilter {
	Q_OBJECT
	Q_DISABLE_COPY_MOVE(QgPolyCreateFilter)


public:
	QgPolyCreateFilter() = default;
	bool eventFilter(QObject *, QEvent *e) override;
	void finalize(bool);

	std::vector<qg_polygon_ref> bool_objs;
};

class QTCAD_EXPORT QgPolyUpdateFilter : public QgPolyFilter {
	Q_OBJECT
	Q_DISABLE_COPY_MOVE(QgPolyUpdateFilter)


public:
	QgPolyUpdateFilter() = default;
	bool eventFilter(QObject *, QEvent *e) override;

	std::vector<qg_polygon_ref> bool_objs;
};

class QTCAD_EXPORT QgPolySelectFilter : public QgPolyFilter {
	Q_OBJECT
	Q_DISABLE_COPY_MOVE(QgPolySelectFilter)


public:
	QgPolySelectFilter() = default;
	bool eventFilter(QObject *, QEvent *e) override;
};

class QTCAD_EXPORT QgPolyPointFilter : public QgPolyFilter {
	Q_OBJECT
	Q_DISABLE_COPY_MOVE(QgPolyPointFilter)


public:
	QgPolyPointFilter() = default;
	bool eventFilter(QObject *, QEvent *e) override;
};

class QTCAD_EXPORT QgPolyMoveFilter : public QgPolyFilter {
	Q_OBJECT
	Q_DISABLE_COPY_MOVE(QgPolyMoveFilter)


public:
	QgPolyMoveFilter() = default;
	bool eventFilter(QObject *, QEvent *e) override;
	std::vector<qg_polygon_ref> move_objs;

private:
	point_t m_prev_point = VINIT_ZERO;
	bool m_prev_point_valid = false;
};


using QPolyCreateFilter = QgPolyCreateFilter;
using QPolyUpdateFilter = QgPolyUpdateFilter;
using QPolySelectFilter = QgPolySelectFilter;
using QPolyPointFilter = QgPolyPointFilter;
using QPolyMoveFilter = QgPolyMoveFilter;

#endif /* QGPOLYFILTER_H */

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
