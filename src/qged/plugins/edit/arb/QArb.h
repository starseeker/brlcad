/*                         Q A R B . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QArb.h
 *
 * Topology-aware ARB adapter for the shared GED primitive-edit session.
 */

#ifndef QARB_H
#define QARB_H

#include <QString>
#include <QWidget>
#include <vector>

#include "BObol/BInput.h"
#include "ged/edit.h"
#include "qtcad/QgTypes.h"

class QComboBox;
class QgPluginContext;
class QgPrimitiveEdit;
struct QArbManipulatorState;
struct QArbPresentationState;

class QArb : public QWidget
{
    Q_OBJECT

    public:
	QArb();
	~QArb() override;

	void setContext(QgPluginContext *ctx);

    signals:
	void view_updated(QgViewUpdateFlags);

    private slots:
	void sync_selection();
	void reset_for_database();
	void refresh_preview();
	void update_preview(int kind, qulonglong revision);
	void target_changed(const QString &path);
	void mode_changed(int mode);

    protected:
	bool eventFilter(QObject *, QEvent *) override;

    private:
	struct ged *getGed() const;
	int command_id(const char *name) const;
	void clear_preview();
	void clear_labels();
	void clear_manipulators();
	bool rebuild_preview(uint64_t revision);
	void update_manipulators(uint64_t revision);
	void sync_session_selection(int command_id);
	int handle_manipulator_input(QArbManipulatorState *,
	    BObolInputAction action, const BObolInputEvent *event);
	static int manipulator_input(void *user_data, BObolInputAction action,
	    const BObolInputEvent *event);

	QgPrimitiveEdit *editor = nullptr;
	QComboBox *edit_mode = nullptr;
	QgPluginContext *m_ctx = nullptr;
	QString selection_path;
	QString preview_path;
	QArbPresentationState *state = nullptr;
};

#endif /* QARB_H */

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
