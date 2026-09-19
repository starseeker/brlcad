/*                    Q E X T R U D E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file QExtrude.cpp */

#include "common.h"

#include <QEvent>
#include <QVBoxLayout>

#include "ged.h"
#include "ged/edit.h"
#include "ged/selection.h"
#include "rt/db_fullpath.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgPrimitiveEdit.h"
#include "qtcad/QgSignalFlags.h"
#include "QExtrude.h"


QExtrude::QExtrude()
    : QWidget()
{
    QVBoxLayout *layout = new QVBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);
    editor = new QgPrimitiveEdit(this);
    editor->setObjectName(QStringLiteral("extrude.sharedPrimitiveEditor"));
    layout->addWidget(editor);

    connect(editor, &QgPrimitiveEdit::sessionEvent,
	this, &QExtrude::update_preview);
}


void
QExtrude::setContext(QgPluginContext *ctx)
{
    m_ctx = ctx;
    editor->setGed(getGed());
}


struct ged *
QExtrude::getGed() const
{
    return m_ctx ? m_ctx->getGed() : nullptr;
}


void
QExtrude::refresh_preview()
{
    editor->refreshFromSession();
    emit view_updated(QG_VIEW_REFRESH);
}


void
QExtrude::update_preview(int kindValue, qulonglong UNUSED(revision))
{
    const enum ged_edit_session_event_kind kind =
	static_cast<enum ged_edit_session_event_kind>(kindValue);
    emit view_updated(kind == GED_EDIT_SESSION_COMMIT ? QG_VIEW_DB :
	QG_VIEW_REFRESH);
}


void
QExtrude::sync_selection()
{
    struct ged *gedp = getGed();
    if (!gedp || !gedp->dbip)
	return;

    QString selectedPath;
    if (ged_selection_count(gedp, nullptr) == 1) {
	struct bu_vls paths = BU_VLS_INIT_ZERO;
	(void)ged_selection_list_paths(gedp, nullptr, &paths);
	selectedPath = QString::fromUtf8(bu_vls_cstr(&paths)).trimmed();
	bu_vls_free(&paths);

	struct db_full_path fullPath;
	db_full_path_init(&fullPath);
	const QByteArray pathBytes = selectedPath.toUtf8();
	const bool isExtrude = !selectedPath.isEmpty() &&
	    db_string_to_path(&fullPath, gedp->dbip,
		pathBytes.constData()) == 0 &&
	    DB_FULL_PATH_CUR_DIR(&fullPath) &&
	    DB_FULL_PATH_CUR_DIR(&fullPath)->d_minor_type ==
		DB5_MINORTYPE_BRLCAD_EXTRUDE;
	db_free_full_path(&fullPath);
	if (!isExtrude)
	    selectedPath.clear();
    }

    if (selectedPath.isEmpty()) {
	if (!selection_path.isEmpty() &&
	    editor->targetPath() == selection_path)
	    editor->setTargetPath(QString());
	selection_path.clear();
	return;
    }
    selection_path = selectedPath;
    editor->setTargetPath(selectedPath);
}


void
QExtrude::reset_for_database()
{
    selection_path.clear();
    editor->setGed(nullptr);
    editor->setTargetPath(QString());
    editor->setGed(getGed());
    emit view_updated(QG_VIEW_REFRESH);
}


bool
QExtrude::eventFilter(QObject *, QEvent *)
{
    return false;
}

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
