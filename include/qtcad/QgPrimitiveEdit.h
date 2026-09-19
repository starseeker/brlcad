/*                 Q G P R I M I T I V E E D I T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file QgPrimitiveEdit.h
 *
 * Descriptor-generated Qt controls for a GED primitive edit session.
 *
 * This widget is deliberately geometry-neutral.  Librt descriptors define
 * the available operations and values, while GED owns the path-scoped edit
 * transaction.  Commands, other widgets, and retained scene manipulators may
 * therefore update the same state and this widget will follow their revisions.
 */

#ifndef QGPRIMITIVEEDIT_H
#define QGPRIMITIVEEDIT_H

#include <QWidget>

#include "ged/edit.h"
#include "qtcad/defines.h"

class QString;

struct ged;

class QTCAD_EXPORT QgPrimitiveEdit : public QWidget
{
    Q_OBJECT

    public:
	explicit QgPrimitiveEdit(QWidget *parent = nullptr);
	~QgPrimitiveEdit() override;

	/** Change the non-owning GED endpoint observed by this widget. */
	void setGed(struct ged *gedp);
	struct ged *ged() const;

	/** Select and begin/join an edit session for a database instance path. */
	void setTargetPath(const QString &path);
	QString targetPath() const;

	/** Current opaque session, or GED_EDIT_SESSION_REF_NULL. */
	ged_edit_session_ref session() const;

    public slots:
	void refreshFromSession();
	void checkpoint();
	void revert();
	void commit();
	void cancel();

    signals:
	/** Target text changed after it was resolved by the editor. */
	void targetChanged(const QString &path);

	/**
	 * A shared session transition was observed.  This is suitable for
	 * refreshing a retained preview; consumers must query GED for state.
	 */
	void sessionEvent(int kind, qulonglong revision);

	/** A user-facing validation or session error occurred. */
	void errorMessage(const QString &message);

    private:
	class Private;
	Private *d;

	static void sessionObserver(struct ged *,
	    const struct ged_edit_session_event *, void *);
	void attachTarget();
	void rebuildOperations();
	void rebuildParameters();
	void selectSessionOperation();
	void applyCurrentOperation();
	void detachSession(bool cancelOwned);
	void handleSessionEvent(enum ged_edit_session_event_kind,
	    ged_edit_session_ref, const QString &, const QString &,
	    enum ged_edit_session_invalidation_reason, uint64_t);
};

#endif /* QGPRIMITIVEEDIT_H */

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
