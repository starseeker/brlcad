/*                Q G P L U G I N C O N T E X T . C P P
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

#include "common.h"

#include "qtcad/QgModel.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgView.h"

QgPluginNotifier::QgPluginNotifier(QObject *parent)
    : QObject(parent)
{
}

QgPluginNotifier::~QgPluginNotifier() = default;

qg_legacy_view *
QgPluginContext::activeView() const
{
    QgView *view = this->getViewWidget();
    if (view)
	return view->view();
    return (this->model && this->model->session()) ?
	this->model->session()->activeView() : nullptr;
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
