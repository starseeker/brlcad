/*                        A R B . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include <QIcon>
#include <QMetaObject>
#include <QPixmap>
#include <QSizePolicy>

#include "qtcad/QgPluginContext.h"
#include "qtcad/QgToolPalette.h"
#include "qtcad/QgView.h"
#include "EditArbFactory.h"

EditArbTool::EditArbTool(QgPluginContext *ctx, QObject *parent)
    : QgToolBase(ctx, parent)
{
}

QgToolPaletteElement *
EditArbTool::paletteElement()
{
    QgToolPaletteElement *element = QgToolBase::paletteElement();
    if (m_arb && element && !m_extra_wired) {
	QObject::connect(m_arb, &QArb::view_updated,
	    element, &QgToolPaletteElement::element_view_changed);
	if (m_ctx && m_ctx->notifier) {
	    QObject::connect(m_ctx->notifier, &QgPluginNotifier::viewUpdated,
		this, [this](QgViewUpdateFlags flags) {
		    if (m_arb && (flags & QG_VIEW_SELECT))
			QMetaObject::invokeMethod(m_arb, "sync_selection",
			    Qt::DirectConnection);
		});
	}
	m_extra_wired = true;
    }
    return element;
}

QWidget *
EditArbTool::createWidget()
{
    m_arb = new QArb();
    m_arb->setContext(m_ctx);
    m_arb->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    return m_arb;
}

QIcon *
EditArbTool::createIcon()
{
    return new QIcon(QPixmap(":arb.svg"));
}

void
EditArbTool::refresh()
{
    if (!m_arb)
	return;
    m_arb->blockSignals(true);
    QMetaObject::invokeMethod(m_arb, "refresh_preview", Qt::DirectConnection);
    m_arb->blockSignals(false);
}

void
EditArbTool::onDbChanged()
{
    if (!m_arb)
	return;
    m_arb->blockSignals(true);
    QMetaObject::invokeMethod(m_arb, "reset_for_database",
	Qt::DirectConnection);
    m_arb->blockSignals(false);
}

void
EditArbTool::onViewChanged()
{
    refresh();
}

void
EditArbTool::attachToView(QgView *view)
{
    if (m_arb && view)
	view->add_event_filter(m_arb);
}

void
EditArbTool::detachFromView(QgView *view)
{
    if (m_arb && view)
	view->clear_event_filter(m_arb);
}

QgPluginDescriptor
EditArbFactory::descriptor() const
{
    QgPluginDescriptor descriptor;
    descriptor.id = "org.brlcad.qged.edit.arb";
    descriptor.displayName = "ARB Edit";
    descriptor.category = "qged.object";
    descriptor.iconName = ":arb.svg";
    descriptor.sortKey = 1050;
    descriptor.requiresOpenDb = true;
    descriptor.requiresView = true;
    descriptor.description =
	"Edit ARB4 through ARB8 primitives with topology-aware viewport handles.";
    descriptor.vendor = "BRL-CAD";
    descriptor.version = "1.0";
    return descriptor;
}

QgToolBase *
EditArbFactory::create(QgPluginContext *ctx, QObject *parent)
{
    return new EditArbTool(ctx, parent);
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
