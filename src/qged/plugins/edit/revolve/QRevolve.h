/*                       Q R E V O L V E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QRevolve.h */

#ifndef QREVOLVE_H
#define QREVOLVE_H

#include <QString>
#include <QWidget>

#include "qtcad/QgTypes.h"
class QgPluginContext;
class QgPrimitiveEdit;

/** Shared-session adapter for descriptor-generated revolve editing. */
class QRevolve : public QWidget
{
    Q_OBJECT

    public:
	QRevolve();
	~QRevolve() override = default;

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

#endif /* QREVOLVE_H */

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
