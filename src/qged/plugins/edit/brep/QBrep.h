/*                        Q B R E P . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this file; see the file named COPYING for more information.
 */
/** @file QBrep.h
 *
 * Transactional, topology-aware NURBS control-vertex editor for B-reps.
 */

#ifndef QBREP_H
#define QBREP_H

#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QDoubleSpinBox>

#include "bn/tol.h"
#include "qtcad/QgBrepFilter.h"
#include "qtcad/QgTypes.h"
#include "../qged_edit_preview_util.h"

class QgPluginContext;
class QgView;

class QBrep : public QWidget
{
    Q_OBJECT

public:
    QBrep();
    ~QBrep() override;

    void setContext(QgPluginContext *ctx) { m_ctx = ctx; }
    void attachToView(QgView *view);
    void detachFromView(QgView *view);
    QgBrepMoveCVFilter *viewFilter() const { return m_filter; }

signals:
    void view_updated(QgViewUpdateFlags);

public slots:
    void refresh();

private slots:
    void load_from_db();
    void apply_to_db();
    void reset_from_db();
    void set_selected_weight();
    void selection_changed(int face, int cv_u, int cv_v, bool topology_safe);
    void edit_changed();
    void edit_rejected(const QString &reason);

private:
    struct ged *getGed() const;
    void clear_edit_state();
    void clear_preview(bool cancel);
    void update_preview();
    void update_cv_overlay();
    void set_status(const QString &message, bool error = false);

    QgPluginContext *m_ctx = nullptr;
    QgView *m_view = nullptr;
    struct directory *m_dp = nullptr;
    struct rt_edit *m_es = nullptr;
    struct bn_tol m_tol;
    bool m_baseline_valid = false;
    bool m_dirty = false;

    QgBrepMoveCVFilter *m_filter = nullptr;
    struct qged_edit_feature_ref m_preview = QGED_EDIT_FEATURE_REF_NULL;
    struct qged_edit_feature_ref m_cv_overlay = QGED_EDIT_FEATURE_REF_NULL;

    QLineEdit *m_name = nullptr;
    QLabel *m_selection = nullptr;
    QLabel *m_status = nullptr;
    QDoubleSpinBox *m_weight = nullptr;
    QPushButton *m_set_weight = nullptr;
    QPushButton *m_apply = nullptr;
    QPushButton *m_reset = nullptr;
};

#endif /* QBREP_H */

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
