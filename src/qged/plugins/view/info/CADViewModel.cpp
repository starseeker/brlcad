/*                  C A D V I E W M O D E L . C P P
 * BRL-CAD
 *
 * Copyright (c) 2020-2026 United States Government as represented by
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
/** @file CADViewModel.cpp
 *
 * Brief description
 *
 */

#include <QPainter>
#include <QString>

#include "bu/vls.h"
#include "qtcad/QgLegacyView.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgSignalFlags.h"
#include "rt/view.h"
#include "CADViewModel.h"

static qg_legacy_view *
qged_view_model_view(const QgPluginContext *ctx)
{
    return ctx ? ctx->activeView() : nullptr;
}

CADViewModel::CADViewModel(QObject *parentobj)
    : QgKeyValModel(parentobj)
{
    m_root = NULL;
    refresh(QG_VIEW_REFRESH);
}

CADViewModel::~CADViewModel()
{
    if (m_root)
	delete m_root;
}

void
CADViewModel::update()
{
    refresh(QG_VIEW_REFRESH);
}

void
CADViewModel::refresh(unsigned long long)
{
    qg_legacy_view *v = qged_view_model_view(m_ctx);
    if (!v)
	return;

    struct bu_vls val = BU_VLS_INIT_ZERO;
    QMap<QString, QgKeyValNode*> standard_nodes;
    int i = 0;
    struct rt_view_info view_info;
    vect_t aet;
    mat_t view_center;
    if (m_root)
	delete m_root;
    m_root = new QgKeyValNode();
    beginResetModel();

    qg_legacy_view_info_get(v, &view_info);
    qg_legacy_view_aet_get(v, aet);
    qg_legacy_view_center_get(v, view_center);

    const char *view_name = qg_legacy_view_name_get(v);
    standard_nodes.insert("Name", add_pair("Name", view_name ? view_name : "", m_root, i));
    bu_vls_sprintf(&val, "%g", view_info.size);
    standard_nodes.insert("Size", add_pair("Size", bu_vls_cstr(&val), m_root, i));
    bu_vls_sprintf(&val, "%d", view_info.width);
    standard_nodes.insert("Width", add_pair("Width", bu_vls_cstr(&val), m_root, i));
    bu_vls_sprintf(&val, "%d", view_info.height);
    standard_nodes.insert("Height", add_pair("Height", bu_vls_cstr(&val), m_root, i));
    bu_vls_sprintf(&val, "%g", aet[0]);
    standard_nodes.insert("Az", add_pair("Az", bu_vls_cstr(&val), m_root, i));
    bu_vls_sprintf(&val, "%g", aet[1]);
    standard_nodes.insert("El", add_pair("El", bu_vls_cstr(&val), m_root, i));
    bu_vls_sprintf(&val, "%g", aet[2]);
    standard_nodes.insert("Tw", add_pair("Tw", bu_vls_cstr(&val), m_root, i));

    vect_t center;
    MAT_DELTAS_GET_NEG(center, view_center);
    bu_vls_sprintf(&val, "%g", center[0]);
    standard_nodes.insert("Center[X]", add_pair("Center[X]", bu_vls_cstr(&val), m_root, i));
    bu_vls_sprintf(&val, "%g", center[1]);
    standard_nodes.insert("Center[Y]", add_pair("Center[Y]", bu_vls_cstr(&val), m_root, i));
    bu_vls_sprintf(&val, "%g", center[2]);
    standard_nodes.insert("Center[Z]", add_pair("Center[Z]", bu_vls_cstr(&val), m_root, i));

    endResetModel();

    bu_vls_free(&val);
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
