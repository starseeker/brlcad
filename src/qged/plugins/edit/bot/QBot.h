/*                          Q B O T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QBot.h
 *
 * Scalable BoT adapter for the shared GED primitive-edit session.
 */

#ifndef QBOT_H
#define QBOT_H

#include <QString>
#include <QWidget>

#include "BObol/BViewStore.h"
#include "ged/edit.h"
#include "qtcad/QgTypes.h"

class QComboBox;
class BObolViewController;
class QgPluginContext;
class QgPrimitiveEdit;
struct QBotPresentationState;

class QBot : public QWidget
{
    Q_OBJECT

    public:
	QBot();
	~QBot() override;

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
	struct ged *getGed() const;
	void clear_preview();
	bool rebuild_preview(uint64_t revision);
	bool sync_view_presentations(bool replace_all);
	void clear_view_presentations();
	bool sync_session_selection(int command_id,
	    bool publish_presentation = true);
	bool patch_session_move(int command_id);
	bool publish_presentations(uint32_t source_id,
	    BObolFeaturePublicationCommit committed = nullptr,
	    void *context = nullptr);
	bool publish_surfaces(uint32_t source_id,
	    BObolFeaturePublicationCommit committed, void *context);
	bool publish_surface(BObolViewController *controller);
	bool publish_handle(struct ged_view_context *view_context);
	BObolFeaturePublication surface_publication(SoNode *node);
	int command_id(const char *name) const;

	QgPrimitiveEdit *editor = nullptr;
	QComboBox *edit_mode = nullptr;
	QgPluginContext *m_ctx = nullptr;
	QString selection_path;
	QString preview_path;
	QBotPresentationState *state = nullptr;
};

#endif /* QBOT_H */

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
