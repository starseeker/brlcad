/*             Q G P A L E T T E C O N T R O L L E R . C P P
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

#include <QSet>

#include "qtcad/QgPaletteController.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgPluginInterfaces.h"
#include "qtcad/QgPluginManager.h"
#include "qtcad/QgToolBase.h"
#include "qtcad/QgToolPalette.h"
#include "qtcad/QgView.h"

QgPaletteController::QgPaletteController(QgToolPalette *palette,
					 QgPluginManager *manager,
					 const QString &category,
					 QgPluginContext *ctx,
					 QObject *parent)
    : QObject(parent)
    , m_palette(palette)
    , m_manager(manager)
    , m_category(category)
    , m_ctx(ctx)
{
    if (m_palette) {
	connect(m_palette, &QgToolPalette::palette_element_selected,
		this, &QgPaletteController::onElementSelected);
	connect(m_palette, &QObject::destroyed, this,
		&QgPaletteController::onPaletteDestroyed,
		Qt::DirectConnection);
    }
    if (m_manager) {
	connect(m_manager, &QgPluginManager::pluginsAboutToUnload,
		this, &QgPaletteController::clear, Qt::DirectConnection);
	connect(m_manager, &QgPluginManager::catalogChanged,
		this, &QgPaletteController::populate);
    }
}

QgPaletteController::~QgPaletteController()
{
    clear();
}

QgToolPalette *
QgPaletteController::palette() const
{
    return m_palette.data();
}

QgView *
QgPaletteController::activeView() const
{
    return m_view.data();
}

void
QgPaletteController::removeTool(const QString &id)
{
    QgToolBase *tool = m_tools.take(id);
    if (!tool)
	return;

    if (tool == m_currentTool) {
	tool->setActive(false);
	detachFilter(tool);
	m_currentTool = nullptr;
	emit currentToolChanged(nullptr);
    }

    QgToolPaletteElement *element = m_elementToId.key(id, nullptr);
    if (element) {
	m_elementToId.remove(element);
	if (m_palette)
	    m_palette->deleteElement(element, false);
    }
    delete tool;
}

void
QgPaletteController::clear()
{
    const QStringList ids = m_tools.keys();
    for (const QString &id : ids)
	removeTool(id);
    m_elementToId.clear();
}

void
QgPaletteController::populate()
{
    if (!m_palette || !m_manager)
	return;
    const QList<IQgToolFactory *> factories =
	m_manager->factories<IQgToolFactory>(m_category);
    QSet<QString> availableIds;
    for (IQgToolFactory *factory : factories) {
	if (!factory)
	    continue;
	const QgPluginDescriptor plugin = factory->descriptor();
	if (!plugin.id.isEmpty())
	    availableIds.insert(plugin.id);
    }

    const QStringList existingIds = m_tools.keys();
    for (const QString &id : existingIds) {
	if (!availableIds.contains(id))
	    removeTool(id);
    }

    for (IQgToolFactory *factory : factories) {
	if (!factory)
	    continue;
	const QgPluginDescriptor plugin = factory->descriptor();
	if (plugin.id.isEmpty() || m_tools.contains(plugin.id))
	    continue;

	QgToolBase *tool = factory->create(m_ctx, this);
	if (!tool)
	    continue;
	QgToolPaletteElement *element = tool->paletteElement();
	if (!element) {
	    delete tool;
	    continue;
	}
	/* Tool palettes are intentionally free to change their layout.  Give the
	 * semantic tool, activation button, and controls stable ids so recorded
	 * GUI workflows do not depend on QScrollArea's private QObject hierarchy
	 * or plugin discovery order. */
	tool->setProperty("qgTestId",
		plugin.id + QStringLiteral(".tool"));
	element->setProperty("qgTestId",
		plugin.id + QStringLiteral(".element"));
	if (element->buttonWidget())
	    element->buttonWidget()->setProperty("qgTestId",
		    plugin.id + QStringLiteral(".activate"));
	if (element->controlsWidget())
	    element->controlsWidget()->setProperty("qgTestId",
		    plugin.id + QStringLiteral(".controls"));

	m_tools.insert(plugin.id, tool);
	m_elementToId.insert(element, plugin.id);
	m_palette->addElement(element);
    }
}

void
QgPaletteController::setActiveView(QgView *view)
{
    if (m_view == view)
	return;
    if (m_currentTool) {
	detachFilter(m_currentTool);
    }
    m_view = view;
    if (m_currentTool) {
	attachFilter(m_currentTool);
    }
}

void
QgPaletteController::attachFilter(QgToolBase *tool)
{
    if (!tool || !m_view)
	return;
    QObject *obj = tool->viewEventFilter();
    if (!obj)
	return;
    IQgViewEventFilter *f = qobject_cast<IQgViewEventFilter *>(obj);
    if (f)
	f->attachToView(m_view);
}

void
QgPaletteController::detachFilter(QgToolBase *tool)
{
    if (!tool || !m_view)
	return;
    QObject *obj = tool->viewEventFilter();
    if (!obj)
	return;
    IQgViewEventFilter *f = qobject_cast<IQgViewEventFilter *>(obj);
    if (f)
	f->detachFromView(m_view);
}

void
QgPaletteController::onElementSelected(QgToolPaletteElement *el)
{
    QString id = m_elementToId.value(el);
    QgToolBase *tool = id.isEmpty() ? nullptr : m_tools.value(id, nullptr);
    if (tool == m_currentTool)
	return;
    if (m_currentTool) {
	m_currentTool->setActive(false);
	detachFilter(m_currentTool);
    }
    m_currentTool = tool;
    if (m_currentTool) {
	m_currentTool->setActive(true);
	attachFilter(m_currentTool);
    }
    emit currentToolChanged(m_currentTool);
}

void
QgPaletteController::onPaletteDestroyed()
{
    /* QObject emits destroyed() before deleting its children.  The QPointer
     * is already clear, but tool controls and event filters are still alive;
     * tear down the tools now so no adapter keeps a raw widget pointer after
     * the palette's child destruction. */
    if (m_currentTool) {
	m_currentTool->setActive(false);
	detachFilter(m_currentTool);
    }
    m_currentTool = nullptr;
    m_elementToId.clear();
    const QList<QgToolBase *> tools = m_tools.values();
    m_tools.clear();
    for (QgToolBase *tool : tools)
	delete tool;
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
