/*                       B R E P . C P P
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

#include <QIcon>
#include <QMetaObject>
#include <QSizePolicy>

#include "qtcad/QgPluginContext.h"
#include "qtcad/QgToolPalette.h"
#include "qtcad/QgView.h"
#include "EditBrepFactory.h"

EditBrepTool::EditBrepTool(QgPluginContext *ctx, QObject *parent)
    : QgToolBase(ctx, parent)
{
}

QgToolPaletteElement *
EditBrepTool::paletteElement()
{
    QgToolPaletteElement *element = QgToolBase::paletteElement();
    if (m_brep && element && !m_extra_wired) {
	QObject::connect(m_brep, &QBrep::view_updated,
		element, &QgToolPaletteElement::element_view_changed);
	m_extra_wired = true;
    }
    return element;
}

QWidget *
EditBrepTool::createWidget()
{
    m_brep = new QBrep;
    m_brep->setContext(m_ctx);
    m_brep->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    return m_brep;
}

QIcon *
EditBrepTool::createIcon()
{
    return new QIcon(":/images/primitives/brep.png");
}

void
EditBrepTool::refresh()
{
    if (m_brep)
	QMetaObject::invokeMethod(m_brep, "refresh", Qt::DirectConnection);
}

void
EditBrepTool::onDbChanged()
{
    refresh();
}

void
EditBrepTool::onViewChanged()
{
    refresh();
}

void
EditBrepTool::attachToView(QgView *view)
{
    if (m_brep)
	m_brep->attachToView(view);
}

void
EditBrepTool::detachFromView(QgView *view)
{
    if (m_brep)
	m_brep->detachFromView(view);
}

QgPluginDescriptor
EditBrepFactory::descriptor() const
{
    QgPluginDescriptor descriptor;
    descriptor.id = "org.brlcad.qged.edit.brep";
    descriptor.displayName = "B-rep NURBS Edit";
    descriptor.category = "qged.object";
    descriptor.iconName = ":/images/primitives/brep.png";
    descriptor.sortKey = 1400;
    descriptor.requiresOpenDb = true;
    descriptor.requiresView = true;
    descriptor.description =
	"Topology-aware control-point editing for B-rep NURBS surfaces.";
    descriptor.vendor = "BRL-CAD";
    descriptor.version = "0.1";
    return descriptor;
}

QgToolBase *
EditBrepFactory::create(QgPluginContext *ctx, QObject *parent)
{
    return new EditBrepTool(ctx, parent);
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
