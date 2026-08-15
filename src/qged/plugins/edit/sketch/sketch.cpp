/*                     S K E T C H . C P P
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
#include "EditSketchFactory.h"

EditSketchTool::EditSketchTool(QgPluginContext *ctx, QObject *parent) :
    QgToolBase(ctx, parent)
{
}

QgToolPaletteElement *
EditSketchTool::paletteElement()
{
    QgToolPaletteElement *element = QgToolBase::paletteElement();
    if (m_sketch && element && !m_extra_wired) {
	QObject::connect(m_sketch, &QSketch::view_updated,
	    element, &QgToolPaletteElement::element_view_changed);
	if (m_ctx && m_ctx->notifier) {
	    QObject::connect(m_ctx->notifier, &QgPluginNotifier::viewUpdated,
		this, [this](QgViewUpdateFlags flags) {
		    if (m_sketch && (flags & QG_VIEW_SELECT))
			QMetaObject::invokeMethod(m_sketch, "sync_selection",
			    Qt::DirectConnection);
		});
	}
	m_extra_wired = true;
    }
    return element;
}

QWidget *
EditSketchTool::createWidget()
{
    m_sketch = new QSketch();
    m_sketch->setContext(m_ctx);
    m_sketch->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    return m_sketch;
}

QIcon *
EditSketchTool::createIcon()
{
    return new QIcon(QPixmap(":sketch_edit.svg"));
}

void
EditSketchTool::refresh()
{
    if (!m_sketch)
	return;
    m_sketch->blockSignals(true);
    QMetaObject::invokeMethod(m_sketch, "refresh_preview",
	Qt::DirectConnection);
    m_sketch->blockSignals(false);
}

void
EditSketchTool::onDbChanged()
{
    if (!m_sketch)
	return;
    m_sketch->blockSignals(true);
    QMetaObject::invokeMethod(m_sketch, "reset_for_database",
	Qt::DirectConnection);
    m_sketch->blockSignals(false);
}

void
EditSketchTool::onViewChanged()
{
    refresh();
}

void
EditSketchTool::attachToView(QgView *view)
{
    if (m_sketch && view)
	view->add_event_filter(m_sketch);
}

void
EditSketchTool::detachFromView(QgView *view)
{
    if (m_sketch && view)
	view->clear_event_filter(m_sketch);
}

QgPluginDescriptor
EditSketchFactory::descriptor() const
{
    QgPluginDescriptor descriptor;
    descriptor.id = "org.brlcad.qged.edit.sketch";
    descriptor.displayName = "Sketch Edit";
    descriptor.category = "qged.object";
    descriptor.iconName = ":sketch_edit.svg";
    descriptor.sortKey = 1060;
    descriptor.requiresOpenDb = true;
    descriptor.requiresView = true;
    descriptor.description =
	"Edit sketch vertices and curve segments through one shared session.";
    descriptor.vendor = "BRL-CAD";
    descriptor.version = "1.0";
    return descriptor;
}

QgToolBase *
EditSketchFactory::create(QgPluginContext *ctx, QObject *parent)
{
    return new EditSketchTool(ctx, parent);
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
