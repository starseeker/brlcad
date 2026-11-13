/*                 E D I T A R B F A C T O R Y . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef EDITARBFACTORY_H
#define EDITARBFACTORY_H

#include <QIcon>

#include "qtcad/QgPluginDescriptor.h"
#include "qtcad/QgPluginInterfaces.h"
#include "qtcad/QgToolBase.h"
#include "qtcad/QgToolPalette.h"
#include "QArb.h"

class QgPluginContext;
class QgView;

class EditArbTool : public QgToolBase, public IQgViewEventFilter
{
    Q_OBJECT
    Q_INTERFACES(IQgViewEventFilter)

    public:
	explicit EditArbTool(QgPluginContext *ctx, QObject *parent = nullptr);

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
	QArb *m_arb = nullptr;
	bool m_extra_wired = false;
};

class EditArbFactory : public QObject, public IQgToolFactory
{
    Q_OBJECT
    Q_INTERFACES(IQgToolFactory)
    Q_PLUGIN_METADATA(IID "org.brlcad.qtcad.IQgToolFactory/1" FILE "arb.json")

    public:
	QgPluginDescriptor descriptor() const override;
	QgToolBase *create(QgPluginContext *ctx, QObject *parent = nullptr) override;
};

#endif /* EDITARBFACTORY_H */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
