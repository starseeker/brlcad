/*                        Q E L L . H
 * BRL-CAD
 *
 * Copyright (c) 2022-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file QEll.h */

#ifndef QELL_H
#define QELL_H

#include <QString>
#include <QWidget>
#include <vector>

#include "BObol/BInput.h"
#include "qtcad/QgTypes.h"
#include "vmath.h"

class QgPluginContext;
class QgPrimitiveEdit;
struct QEllManipulatorState;

class QEll : public QWidget
{
    Q_OBJECT

    public:
	QEll();
	~QEll() override;

	void setContext(QgPluginContext *ctx);

    signals:
	void view_updated(QgViewUpdateFlags);

    private slots:
	void sync_selection();
	void reset_for_database();
	void refresh_preview();
	void update_preview(int kind, qulonglong revision);
	void target_changed(const QString &path);

    protected:
	bool eventFilter(QObject *, QEvent *) override;

    private:
	void clear_preview();
	void clear_labels();
	void clear_manipulator();
	void update_manipulator(const point_t center, const vect_t axis_a,
	    const vect_t axis_b, const vect_t axis_c,
	    const fastf_t local_lengths[3], uint64_t revision);
	int handle_manipulator_input(QEllManipulatorState *,
	    BObolInputAction action,
	    const BObolInputEvent *event);
	static int manipulator_input(void *user_data, BObolInputAction action,
	    const BObolInputEvent *event);
	struct ged *getGed() const;

	QgPrimitiveEdit *editor = nullptr;
	QgPluginContext *m_ctx = nullptr;
	QString selection_path;
	QString preview_path;
	std::vector<QEllManipulatorState *> manipulator_states;
};

#endif /* QELL_H */

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
