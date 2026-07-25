/*               E D I T B R E P F A C T O R Y . H
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

#ifndef EDITBREPFACTORY_H
#define EDITBREPFACTORY_H

#include <QIcon>

#include "qtcad/QgPluginDescriptor.h"
#include "qtcad/QgPluginInterfaces.h"
#include "qtcad/QgToolBase.h"
#include "qtcad/QgToolPalette.h"
#include "QBrep.h"

class QgPluginContext;
class QgView;

class EditBrepTool : public QgToolBase, public IQgViewEventFilter
{
    Q_OBJECT
    Q_INTERFACES(IQgViewEventFilter)

public:
    explicit EditBrepTool(QgPluginContext *ctx, QObject *parent = nullptr);

    QgToolPaletteElement *paletteElement() override;
    QWidget *createWidget() override;
    QIcon *createIcon() override;

    QObject *viewEventFilter() override { return this; }
    void attachToView(QgView *view) override;
    void detachFromView(QgView *view) override;

public slots:
    void refresh() override;
    void onDbChanged() override;
    void onViewChanged() override;

private:
    QBrep *m_brep = nullptr;
    bool m_extra_wired = false;
};

class EditBrepFactory : public QObject, public IQgToolFactory
{
    Q_OBJECT
    Q_INTERFACES(IQgToolFactory)
    Q_PLUGIN_METADATA(IID "org.brlcad.qtcad.IQgToolFactory/1" FILE "brep.json")

public:
    QgPluginDescriptor descriptor() const override;
    QgToolBase *create(QgPluginContext *ctx, QObject *parent = nullptr) override;
};

#endif /* EDITBREPFACTORY_H */

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
