/*                       Q S K E T C H . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QSketch.h
 *
 * Session-backed 2-D sketch editing with retained vertex/curve handles.
 */

#ifndef QSKETCH_H
#define QSKETCH_H

#include <QString>
#include <QWidget>

#include "BObol/BInput.h"
#include "ged/edit.h"
#include "qtcad/QgTypes.h"

class QComboBox;
class QgPluginContext;
class QgPrimitiveEdit;
class QgSketchEdit;
struct QSketchManipulatorState;
struct QSketchPresentationState;

class QSketch : public QWidget
{
    Q_OBJECT

    public:
	QSketch();
	~QSketch() override;

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
	void clear_preview();
	void clear_manipulators();
	bool rebuild_preview(uint64_t revision);
	void update_manipulators(uint64_t revision);
	void sync_session_selection(int command_id);
	int handle_manipulator_input(QSketchManipulatorState *,
	    BObolInputAction action, const BObolInputEvent *event);
	static int manipulator_input(void *user_data, BObolInputAction action,
	    const BObolInputEvent *event);

	QgPrimitiveEdit *editor = nullptr;
	QgSketchEdit *interaction = nullptr;
	QComboBox *edit_mode = nullptr;
	QgPluginContext *m_ctx = nullptr;
	QString selection_path;
	QString preview_path;
	QSketchPresentationState *state = nullptr;
};

#endif /* QSKETCH_H */

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
