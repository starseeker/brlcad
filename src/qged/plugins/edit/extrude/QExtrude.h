/*                      Q E X T R U D E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QExtrude.h */

#ifndef QEXTRUDE_H
#define QEXTRUDE_H

#include <QString>
#include <QWidget>

#include "qtcad/QgTypes.h"
class QgPluginContext;
class QgPrimitiveEdit;

/**
 * Extrusion adapter for the shared descriptor-generated primitive editor.
 *
 * The widget owns no primitive copy.  Its preview is rebuilt from a copied
 * snapshot of the path-scoped GED edit session, so CLI and GUI changes share
 * exactly one transaction and revision stream.
 */
class QExtrude : public QWidget
{
    Q_OBJECT

    public:
	QExtrude();
	~QExtrude() override = default;

	void setContext(QgPluginContext *ctx);

    signals:
	void view_updated(QgViewUpdateFlags);

    private slots:
	void sync_selection();
	void reset_for_database();
	void refresh_preview();
	void update_preview(int kind, qulonglong revision);

    protected:
	bool eventFilter(QObject *, QEvent *) override;

    private:
	struct ged *getGed() const;

	QgPrimitiveEdit *editor = nullptr;
	QgPluginContext *m_ctx = nullptr;
	QString selection_path;
};

#endif /* QEXTRUDE_H */

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
