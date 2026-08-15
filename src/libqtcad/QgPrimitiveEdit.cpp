/*              Q G P R I M I T I V E E D I T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file QgPrimitiveEdit.cpp */

#include "common.h"

#include <algorithm>
#include <cfloat>
#include <vector>

#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMetaObject>
#include <QPushButton>
#include <QRegularExpression>
#include <QSignalBlocker>
#include <QSpinBox>
#include <QStringList>
#include <QThread>
#include <QTimer>
#include <QVBoxLayout>

#include "bu/str.h"
#include "rt/edit.h"
#include "qtcad/QgPrimitiveEdit.h"


namespace {

static bool
same_session(ged_edit_session_ref a, ged_edit_session_ref b)
{
    return a.owner == b.owner && a.id == b.id &&
	a.generation == b.generation;
}


static QString
canonical_event_path(QString path)
{
    path = path.trimmed();
    while (path.size() > 1 && path.endsWith(QLatin1Char('/')))
	path.chop(1);
    if (!path.isEmpty() && !path.startsWith(QLatin1Char('/')))
	path.prepend(QLatin1Char('/'));
    return path;
}


static bool
same_target_path(const QString &targetPath, const QString &eventPath)
{
    return !targetPath.trimmed().isEmpty() &&
	canonical_event_path(targetPath) == canonical_event_path(eventPath);
}


static bool
no_limit(fastf_t value)
{
    return value <= -DBL_MAX;
}


static int
param_width(int type)
{
    switch (type) {
	case RT_EDIT_PARAM_POINT2:
	case RT_EDIT_PARAM_VECTOR2:
	    return 2;
	case RT_EDIT_PARAM_POINT:
	case RT_EDIT_PARAM_VECTOR:
	case RT_EDIT_PARAM_COLOR:
	    return 3;
	case RT_EDIT_PARAM_MATRIX:
	    return 16;
	case RT_EDIT_PARAM_STRING:
	case RT_EDIT_PARAM_SCALAR_LIST:
	case RT_EDIT_PARAM_INTEGER_LIST:
	    return 0;
	default:
	    return 1;
    }
}

}


class QgPrimitiveEdit::Private
{
    public:
	struct ParameterControl {
	    const struct rt_edit_param_desc *descriptor = nullptr;
	    std::vector<QWidget *> values;
	};

	struct ged *gedp = nullptr;
	ged_edit_observer_token observer = 0;
	ged_edit_session_ref session = GED_EDIT_SESSION_REF_NULL;
	const struct rt_edit_prim_desc *descriptor = nullptr;
	const struct rt_edit_cmd_desc *command = nullptr;
	uint64_t revision = 0;
	bool ownsSession = false;
	bool refreshing = false;
	bool restarting = false;
	bool formValueError = false;

	QLineEdit *target = nullptr;
	QComboBox *operation = nullptr;
	QGroupBox *parametersBox = nullptr;
	QFormLayout *parameters = nullptr;
	QLabel *status = nullptr;
	QPushButton *apply = nullptr;
	QPushButton *checkpoint = nullptr;
	QPushButton *revert = nullptr;
	QPushButton *commit = nullptr;
	QPushButton *cancel = nullptr;
	std::vector<ParameterControl> controls;
};


QgPrimitiveEdit::QgPrimitiveEdit(QWidget *parent)
    : QWidget(parent), d(new Private)
{
    QVBoxLayout *layout = new QVBoxLayout(this);

    QFormLayout *identity = new QFormLayout;
    d->target = new QLineEdit(this);
    d->target->setObjectName(QStringLiteral("primitiveEdit.target"));
    d->target->setProperty("qgTestId", QStringLiteral("primitiveEdit.target"));
    d->target->setPlaceholderText(QStringLiteral("database object or instance path"));
    identity->addRow(tr("Target"), d->target);

    d->operation = new QComboBox(this);
    d->operation->setObjectName(QStringLiteral("primitiveEdit.operation"));
    d->operation->setProperty("qgTestId", QStringLiteral("primitiveEdit.operation"));
    identity->addRow(tr("Operation"), d->operation);
    layout->addLayout(identity);

    d->parametersBox = new QGroupBox(tr("Parameters"), this);
    d->parameters = new QFormLayout(d->parametersBox);
    layout->addWidget(d->parametersBox);

    d->status = new QLabel(tr("No edit target"), this);
    d->status->setObjectName(QStringLiteral("primitiveEdit.status"));
    d->status->setProperty("qgTestId", QStringLiteral("primitiveEdit.status"));
    d->status->setWordWrap(true);
    layout->addWidget(d->status);

    QGridLayout *actions = new QGridLayout;
    d->apply = new QPushButton(tr("Apply"), this);
    d->apply->setObjectName(QStringLiteral("primitiveEdit.apply"));
    d->apply->setProperty("qgTestId", QStringLiteral("primitiveEdit.apply"));
    d->checkpoint = new QPushButton(tr("Checkpoint"), this);
    d->checkpoint->setObjectName(QStringLiteral("primitiveEdit.checkpoint"));
    d->checkpoint->setProperty("qgTestId", QStringLiteral("primitiveEdit.checkpoint"));
    d->revert = new QPushButton(tr("Revert"), this);
    d->revert->setObjectName(QStringLiteral("primitiveEdit.revert"));
    d->revert->setProperty("qgTestId", QStringLiteral("primitiveEdit.revert"));
    d->commit = new QPushButton(tr("Commit"), this);
    d->commit->setObjectName(QStringLiteral("primitiveEdit.commit"));
    d->commit->setProperty("qgTestId", QStringLiteral("primitiveEdit.commit"));
    d->cancel = new QPushButton(tr("Cancel"), this);
    d->cancel->setObjectName(QStringLiteral("primitiveEdit.cancel"));
    d->cancel->setProperty("qgTestId", QStringLiteral("primitiveEdit.cancel"));
    actions->addWidget(d->apply, 0, 0);
    actions->addWidget(d->checkpoint, 0, 1);
    actions->addWidget(d->revert, 0, 2);
    actions->addWidget(d->commit, 1, 0);
    actions->addWidget(d->cancel, 1, 1);
    layout->addLayout(actions);
    layout->addStretch(1);

    connect(d->target, &QLineEdit::editingFinished,
	this, &QgPrimitiveEdit::attachTarget);
    connect(d->operation,
	QOverload<int>::of(&QComboBox::currentIndexChanged),
	this, [this](int) { rebuildParameters(); });
    connect(d->apply, &QPushButton::clicked,
	this, &QgPrimitiveEdit::applyCurrentOperation);
    connect(d->checkpoint, &QPushButton::clicked,
	this, &QgPrimitiveEdit::checkpoint);
    connect(d->revert, &QPushButton::clicked,
	this, &QgPrimitiveEdit::revert);
    connect(d->commit, &QPushButton::clicked,
	this, &QgPrimitiveEdit::commit);
    connect(d->cancel, &QPushButton::clicked,
	this, &QgPrimitiveEdit::cancel);

    rebuildOperations();
}


QgPrimitiveEdit::~QgPrimitiveEdit()
{
    if (d->gedp && d->observer)
	(void)ged_edit_observer_remove(d->gedp, d->observer);
    detachSession(true);
    delete d;
}


void
QgPrimitiveEdit::setGed(struct ged *gedp)
{
    if (gedp == d->gedp)
	return;
    detachSession(true);
    if (d->gedp && d->observer)
	(void)ged_edit_observer_remove(d->gedp, d->observer);
    d->observer = 0;
    d->gedp = gedp;
    d->session = GED_EDIT_SESSION_REF_NULL;
    d->ownsSession = false;
    d->descriptor = nullptr;
    d->command = nullptr;
    d->revision = 0;
    rebuildOperations();
    if (d->gedp) {
	d->observer = ged_edit_observer_add(d->gedp,
	    QgPrimitiveEdit::sessionObserver, this);
	attachTarget();
    } else {
	d->status->setText(tr("No open database"));
    }
}


struct ged *
QgPrimitiveEdit::ged() const
{
    return d->gedp;
}


void
QgPrimitiveEdit::setTargetPath(const QString &path)
{
    if (d->target->text() != path) {
	QSignalBlocker blocker(d->target);
	d->target->setText(path);
    }
    attachTarget();
}


QString
QgPrimitiveEdit::targetPath() const
{
    return d->target->text().trimmed();
}


ged_edit_session_ref
QgPrimitiveEdit::session() const
{
    return d->session;
}


void
QgPrimitiveEdit::attachTarget()
{
    const QString path = targetPath();
    emit targetChanged(path);

    if (d->gedp && !path.isEmpty() &&
	!ged_edit_session_ref_is_null(d->session)) {
	ged_edit_session_ref targetSession = GED_EDIT_SESSION_REF_NULL;
	const QByteArray pathBytes = path.toUtf8();
	if (ged_edit_session_find(d->gedp, pathBytes.constData(),
		&targetSession) == GED_EDIT_OK &&
	    same_session(targetSession, d->session)) {
	    refreshFromSession();
	    emit sessionEvent(static_cast<int>(GED_EDIT_SESSION_UPDATE),
		static_cast<qulonglong>(d->revision));
	    return;
	}
    }

    detachSession(true);
    const struct rt_edit_prim_desc *previousDescriptor = d->descriptor;
    d->revision = 0;

    if (!d->gedp || path.isEmpty()) {
	d->descriptor = nullptr;
	d->command = nullptr;
	rebuildOperations();
	d->status->setText(path.isEmpty() ? tr("No edit target") :
	    tr("No open database"));
	return;
    }

    const QByteArray pathBytes = path.toUtf8();
    ged_edit_session_ref existing = GED_EDIT_SESSION_REF_NULL;
    const bool joined = ged_edit_session_find(d->gedp,
	pathBytes.constData(), &existing) == GED_EDIT_OK;
    ged_edit_session_ref ref = GED_EDIT_SESSION_REF_NULL;
    enum ged_edit_status result = ged_edit_session_begin(d->gedp,
	pathBytes.constData(), nullptr, &ref);
    if (result != GED_EDIT_OK) {
	d->descriptor = nullptr;
	d->command = nullptr;
	rebuildOperations();
	d->status->setText(tr("Unable to begin edit session for %1").arg(path));
	emit errorMessage(d->status->text());
	return;
    }

    d->session = ref;
    d->ownsSession = !joined;
    struct ged_edit_session_info info;
    if (ged_edit_session_info_get(d->gedp, ref, &info) == GED_EDIT_OK) {
	d->revision = info.revision;
    }
    const struct rt_edit_prim_desc *sessionDescriptor = nullptr;
    if (ged_edit_session_descriptor_get(d->gedp, ref,
	&sessionDescriptor) != GED_EDIT_OK)
	sessionDescriptor = nullptr;
    d->descriptor = sessionDescriptor;
    if (sessionDescriptor != previousDescriptor)
	rebuildOperations();
    selectSessionOperation();
    refreshFromSession();

    /* A newly created session reports BEGIN synchronously from inside
     * ged_edit_session_begin, before this widget can install the returned
     * opaque reference.  A joined session has no lifecycle event at all.
     * Publish one attached-snapshot signal in either case so retained scene
     * controls are initialized from the same authoritative state as widgets.
     * This is especially important after commit/cancel automatically replaces
     * a persistent editor's closed session. */
    emit sessionEvent(static_cast<int>(joined ? GED_EDIT_SESSION_UPDATE :
	GED_EDIT_SESSION_BEGIN), static_cast<qulonglong>(d->revision));
}


void
QgPrimitiveEdit::detachSession(bool cancelOwned)
{
    struct ged *sessionGed = d->gedp;
    const ged_edit_session_ref oldSession = d->session;
    const bool owned = d->ownsSession;
    d->session = GED_EDIT_SESSION_REF_NULL;
    d->ownsSession = false;
    d->revision = 0;
    if (cancelOwned && owned && sessionGed &&
	!ged_edit_session_ref_is_null(oldSession))
	(void)ged_edit_session_cancel(sessionGed, oldSession);
}


void
QgPrimitiveEdit::rebuildOperations()
{
    const int oldCommand = d->operation->currentData().toInt();
    QSignalBlocker blocker(d->operation);
    d->operation->clear();
    if (d->descriptor) {
	for (int i = 0; i < d->descriptor->ncmd; ++i) {
	    const struct rt_edit_cmd_desc *cmd = &d->descriptor->cmds[i];
	    if (rt_edit_cmd_control_class(d->descriptor, cmd) ==
		    RT_EDIT_CONTROL_UNSUPPORTED)
		continue;
	    d->operation->addItem(QString::fromUtf8(cmd->label), cmd->cmd_id);
	}
    }
    int restore = d->operation->findData(oldCommand);
    d->operation->setCurrentIndex(restore >= 0 ? restore :
	(d->operation->count() ? 0 : -1));
    d->operation->setEnabled(d->operation->count() > 0);
    rebuildParameters();
    if (d->descriptor && !d->operation->count())
	d->status->setText(tr("%1 has no operations available in this editor")
	    .arg(QString::fromUtf8(d->descriptor->prim_label)));
}


void
QgPrimitiveEdit::rebuildParameters()
{
    while (d->parameters->rowCount())
	d->parameters->removeRow(0);
    d->controls.clear();
    d->command = nullptr;
    d->formValueError = false;
    d->parametersBox->setTitle(tr("Parameters"));

    if (!d->descriptor || d->operation->currentIndex() < 0) {
	d->parametersBox->setEnabled(false);
	d->apply->setEnabled(false);
	return;
    }

    const int commandId = d->operation->currentData().toInt();
    for (int i = 0; i < d->descriptor->ncmd; ++i) {
	if (d->descriptor->cmds[i].cmd_id == commandId) {
	    d->command = &d->descriptor->cmds[i];
	    break;
	}
    }
    if (!d->command)
	return;

    const enum rt_edit_control_class controlClass =
	rt_edit_cmd_control_class(d->descriptor, d->command);
    const bool formEditable = controlClass == RT_EDIT_CONTROL_GENERATED ||
	controlClass == RT_EDIT_CONTROL_ACTION;

    for (int pi = 0; pi < d->command->nparam; ++pi) {
	const struct rt_edit_param_desc *p = &d->command->params[pi];
	fastf_t rangeMinimum = p->range_min;
	fastf_t rangeMaximum = p->range_max;
	struct rt_edit_param_bounds bounds = {};
	if (d->gedp && !ged_edit_session_ref_is_null(d->session) &&
	    ged_edit_session_parameter_bounds_get(d->gedp, d->session,
		d->command->cmd_id, pi, &bounds) == GED_EDIT_OK) {
	    rangeMinimum = bounds.has_minimum ? bounds.minimum :
		RT_EDIT_PARAM_NO_LIMIT;
	    rangeMaximum = bounds.has_maximum ? bounds.maximum :
		RT_EDIT_PARAM_NO_LIMIT;
	}
	Private::ParameterControl control;
	control.descriptor = p;
	QWidget *field = nullptr;

	if (p->type == RT_EDIT_PARAM_BOOLEAN) {
	    QCheckBox *check = new QCheckBox(d->parametersBox);
	    check->setObjectName(QStringLiteral("primitiveEdit.param.%1")
		.arg(QString::fromUtf8(p->name)));
	    check->setProperty("qgTestId", check->objectName());
	    control.values.push_back(check);
	    field = check;
	} else if (p->type == RT_EDIT_PARAM_ENUM) {
	    QComboBox *combo = new QComboBox(d->parametersBox);
	    combo->setObjectName(QStringLiteral("primitiveEdit.param.%1")
		.arg(QString::fromUtf8(p->name)));
	    combo->setProperty("qgTestId", combo->objectName());
	    for (int ei = 0; ei < p->nenum; ++ei) {
		const int id = p->enum_ids ? p->enum_ids[ei] : ei;
		combo->addItem(QString::fromUtf8(p->enum_labels[ei]), id);
	    }
	    control.values.push_back(combo);
	    field = combo;
	} else if (p->type == RT_EDIT_PARAM_STRING ||
	    p->type == RT_EDIT_PARAM_SCALAR_LIST ||
	    p->type == RT_EDIT_PARAM_INTEGER_LIST) {
	    QLineEdit *edit = new QLineEdit(d->parametersBox);
	    edit->setObjectName(QStringLiteral("primitiveEdit.param.%1")
		.arg(QString::fromUtf8(p->name)));
	    edit->setProperty("qgTestId", edit->objectName());
	    control.values.push_back(edit);
	    field = edit;
	} else {
	    const int width = param_width(p->type);
	    QWidget *container = new QWidget(d->parametersBox);
	    QHBoxLayout *row = new QHBoxLayout(container);
	    row->setContentsMargins(0, 0, 0, 0);
	    for (int vi = 0; vi < width; ++vi) {
		if (p->type == RT_EDIT_PARAM_INTEGER ||
		    p->type == RT_EDIT_PARAM_COLOR) {
		    QSpinBox *spin = new QSpinBox(container);
		    int lower = p->type == RT_EDIT_PARAM_COLOR ? 0 : -1000000000;
		    int upper = p->type == RT_EDIT_PARAM_COLOR ? 255 : 1000000000;
		    if (!no_limit(rangeMinimum))
			lower = static_cast<int>(rangeMinimum);
		    if (!no_limit(rangeMaximum))
			upper = static_cast<int>(rangeMaximum);
		    spin->setRange(lower, upper);
		    spin->setObjectName(QStringLiteral("primitiveEdit.param.%1.%2")
			.arg(QString::fromUtf8(p->name)).arg(vi));
		    spin->setProperty("qgTestId", spin->objectName());
		    control.values.push_back(spin);
		    row->addWidget(spin);
		} else {
		    QDoubleSpinBox *spin = new QDoubleSpinBox(container);
		    spin->setDecimals(9);
		    spin->setRange(no_limit(rangeMinimum) ? -1.0e100 : rangeMinimum,
			no_limit(rangeMaximum) ? 1.0e100 : rangeMaximum);
		    spin->setObjectName(QStringLiteral("primitiveEdit.param.%1.%2")
			.arg(QString::fromUtf8(p->name)).arg(vi));
		    spin->setProperty("qgTestId", spin->objectName());
		    control.values.push_back(spin);
		    row->addWidget(spin);
		}
	    }
	    field = container;
	}

	QString label = QString::fromUtf8(p->label);
	if (p->units && p->units[0])
	    label += QStringLiteral(" (%1)").arg(QString::fromUtf8(p->units));
	d->parameters->addRow(label, field);
	d->controls.push_back(control);
    }

    d->parametersBox->setEnabled(formEditable);
    d->apply->setEnabled(formEditable);
    if (controlClass == RT_EDIT_CONTROL_CUSTOM)
	d->parametersBox->setTitle(tr("Parameters (managed by custom editor)"));
    else
	d->parametersBox->setTitle(tr("Parameters"));
    refreshFromSession();
}


void
QgPrimitiveEdit::selectSessionOperation()
{
    if (!d->gedp || ged_edit_session_ref_is_null(d->session) ||
	!d->descriptor)
	return;

    struct ged_edit_session_info info = {};
    if (ged_edit_session_info_get(d->gedp, d->session, &info) != GED_EDIT_OK ||
	!info.command_id || (d->command &&
	    d->command->cmd_id == info.command_id))
	return;

    const int index = d->operation->findData(info.command_id);
    if (index < 0)
	return;
    {
	QSignalBlocker blocker(d->operation);
	d->operation->setCurrentIndex(index);
    }
    rebuildParameters();
}


void
QgPrimitiveEdit::refreshFromSession()
{
    if (!d->gedp || ged_edit_session_ref_is_null(d->session) ||
	!d->command || d->refreshing)
	return;

    d->refreshing = true;
    struct rt_edit_cmd_values values;
    const enum rt_edit_control_class controlClass =
	rt_edit_cmd_control_class(d->descriptor, d->command);
    const enum ged_edit_status result = controlClass == RT_EDIT_CONTROL_ACTION ?
	GED_EDIT_UNSUPPORTED : ged_edit_session_command_values_get(d->gedp,
	    d->session, d->command->cmd_id, &values);
    d->formValueError = controlClass == RT_EDIT_CONTROL_GENERATED &&
	result != GED_EDIT_OK;
    if (result == GED_EDIT_OK) {
	for (const Private::ParameterControl &control : d->controls) {
	    const struct rt_edit_param_desc *p = control.descriptor;
	    if (p->type == RT_EDIT_PARAM_STRING) {
		if (p->index < 0 || p->index >= RT_EDIT_MAXSTR ||
		    !values.string_valid[p->index] || control.values.empty())
		    continue;
		if (QLineEdit *edit =
			qobject_cast<QLineEdit *>(control.values[0])) {
		    QSignalBlocker blocker(edit);
		    edit->setText(QString::fromUtf8(values.strings[p->index]));
		}
		continue;
	    }
	    if (p->type == RT_EDIT_PARAM_SCALAR_LIST ||
		p->type == RT_EDIT_PARAM_INTEGER_LIST) {
		if (p->index < 0 ||
		    (size_t)p->index >= values.value_count ||
		    control.values.empty())
		    continue;
		QStringList fields;
		for (size_t vi = (size_t)p->index;
		    vi < values.value_count && vi < RT_EDIT_MAXPARA; vi++) {
		    if (values.value_valid[vi])
			fields.append(QString::number(values.values[vi], 'g', 17));
		}
		if (QLineEdit *edit =
			qobject_cast<QLineEdit *>(control.values[0])) {
		    QSignalBlocker blocker(edit);
		    edit->setText(fields.join(QLatin1Char(' ')));
		}
		continue;
	    }
	    for (size_t vi = 0; vi < control.values.size(); ++vi) {
		const size_t index = static_cast<size_t>(p->index) + vi;
		if (index >= values.value_count ||
		    !values.value_valid[index])
		    continue;
		QSignalBlocker blocker(control.values[vi]);
		if (QDoubleSpinBox *doubleSpin =
			qobject_cast<QDoubleSpinBox *>(control.values[vi]))
		    doubleSpin->setValue(values.values[index]);
		else if (QSpinBox *integerSpin =
			qobject_cast<QSpinBox *>(control.values[vi]))
		    integerSpin->setValue(
			static_cast<int>(values.values[index]));
		else if (QCheckBox *check =
			qobject_cast<QCheckBox *>(control.values[vi]))
		    check->setChecked(!NEAR_ZERO(values.values[index],
			SMALL_FASTF));
		else if (QComboBox *combo =
			qobject_cast<QComboBox *>(control.values[vi])) {
		    const int ci = combo->findData(
			static_cast<int>(values.values[index]));
		    if (ci >= 0)
			combo->setCurrentIndex(ci);
		}
	    }
	}
    }

    struct ged_edit_session_info info;
    if (ged_edit_session_info_get(d->gedp, d->session, &info) == GED_EDIT_OK) {
	d->revision = info.revision;
	QString presentation;
	if (controlClass == RT_EDIT_CONTROL_CUSTOM)
	    presentation = tr("; use the topology-aware controls or viewport handles");
	else if (d->formValueError)
	    presentation = tr("; current values unavailable");
	d->status->setText(tr("Active %1 edit, revision %2%3%4")
	    .arg(d->descriptor ? QString::fromUtf8(d->descriptor->prim_label) :
		tr("primitive"))
	    .arg(static_cast<qulonglong>(d->revision))
	    .arg(info.dirty ? tr(" (modified)") : QString())
	    .arg(presentation));
    }
    d->apply->setEnabled((controlClass == RT_EDIT_CONTROL_GENERATED &&
	result == GED_EDIT_OK) || controlClass == RT_EDIT_CONTROL_ACTION);
    d->refreshing = false;
}


void
QgPrimitiveEdit::applyCurrentOperation()
{
    if (!d->gedp || ged_edit_session_ref_is_null(d->session) || !d->command)
	return;
    const enum rt_edit_control_class controlClass =
	rt_edit_cmd_control_class(d->descriptor, d->command);
    if (controlClass != RT_EDIT_CONTROL_GENERATED &&
	controlClass != RT_EDIT_CONTROL_ACTION)
	return;

    fastf_t values[RT_EDIT_MAXPARA] = {0.0};
    std::vector<QByteArray> stringStorage(RT_EDIT_MAXSTR);
    const char *strings[RT_EDIT_MAXSTR] = {nullptr};
    size_t valueCount = 0;
    size_t stringCount = 0;

    for (const Private::ParameterControl &control : d->controls) {
	const struct rt_edit_param_desc *p = control.descriptor;
	if (p->type == RT_EDIT_PARAM_STRING) {
	    if (p->index < 0 || p->index >= RT_EDIT_MAXSTR ||
		control.values.empty())
		continue;
	    QLineEdit *edit = qobject_cast<QLineEdit *>(control.values[0]);
	    stringStorage[p->index] = edit ? edit->text().toUtf8() : QByteArray();
	    strings[p->index] = stringStorage[p->index].constData();
	    stringCount = std::max(stringCount, static_cast<size_t>(p->index + 1));
	    continue;
	}
	if (p->type == RT_EDIT_PARAM_SCALAR_LIST ||
	    p->type == RT_EDIT_PARAM_INTEGER_LIST) {
	    QLineEdit *edit = control.values.empty() ? nullptr :
		qobject_cast<QLineEdit *>(control.values[0]);
	    const QStringList fields = edit ? edit->text().split(
		QRegularExpression(QStringLiteral("[,\\s]+")),
		Qt::SkipEmptyParts) : QStringList();
	    const int maxCount = RT_EDIT_MAXPARA - p->index;
	    if (fields.size() < p->nenum || fields.size() > maxCount ||
		p->index < 0 || p->index + fields.size() > RT_EDIT_MAXPARA) {
		d->status->setText(tr("Parameter %1 requires %2 to %3 values")
		    .arg(QString::fromUtf8(p->label))
		    .arg(p->nenum).arg(maxCount));
		emit errorMessage(d->status->text());
		return;
	    }
	    for (int li = 0; li < fields.size(); li++) {
		bool ok = false;
		const double value = fields[li].toDouble(&ok);
		if (!ok) {
		    d->status->setText(tr(
			"Parameter %1 contains a non-numeric value")
			.arg(QString::fromUtf8(p->label)));
		    emit errorMessage(d->status->text());
		    return;
		}
		values[p->index + li] = value;
	    }
	    valueCount = std::max(valueCount,
		(size_t)(p->index + fields.size()));
	    continue;
	}

	for (size_t vi = 0; vi < control.values.size(); ++vi) {
	    const int index = p->index + static_cast<int>(vi);
	    if (index < 0 || index >= RT_EDIT_MAXPARA)
		continue;
	    QWidget *widget = control.values[vi];
	    if (QDoubleSpinBox *doubleSpin = qobject_cast<QDoubleSpinBox *>(widget))
		values[index] = doubleSpin->value();
	    else if (QSpinBox *integerSpin = qobject_cast<QSpinBox *>(widget))
		values[index] = integerSpin->value();
	    else if (QCheckBox *check = qobject_cast<QCheckBox *>(widget))
		values[index] = check->isChecked() ? 1.0 : 0.0;
	    else if (QComboBox *combo = qobject_cast<QComboBox *>(widget))
		values[index] = combo->currentData().toInt();
	    valueCount = std::max(valueCount, static_cast<size_t>(index + 1));
	}
    }

    struct ged_edit_command_input input;
    input.command_id = d->command->cmd_id;
    input.values = values;
    input.value_count = valueCount;
    input.strings = strings;
    input.string_count = stringCount;
    input.view = nullptr;
    const enum ged_edit_status result = ged_edit_session_apply(d->gedp,
	d->session, &input);
    if (result != GED_EDIT_OK) {
	d->status->setText(tr("Edit operation was rejected (%1)").arg(result));
	emit errorMessage(d->status->text());
    }
}


void
QgPrimitiveEdit::checkpoint()
{
    if (d->gedp && !ged_edit_session_ref_is_null(d->session) &&
	ged_edit_session_checkpoint(d->gedp, d->session) != GED_EDIT_OK)
	emit errorMessage(tr("Unable to checkpoint edit session"));
}


void
QgPrimitiveEdit::revert()
{
    if (d->gedp && !ged_edit_session_ref_is_null(d->session) &&
	ged_edit_session_revert(d->gedp, d->session) != GED_EDIT_OK)
	emit errorMessage(tr("Unable to revert edit session"));
}


void
QgPrimitiveEdit::commit()
{
    if (d->gedp && !ged_edit_session_ref_is_null(d->session) &&
	ged_edit_session_commit(d->gedp, d->session) != GED_EDIT_OK)
	emit errorMessage(tr("Unable to commit edit session"));
}


void
QgPrimitiveEdit::cancel()
{
    if (d->gedp && !ged_edit_session_ref_is_null(d->session) &&
	ged_edit_session_cancel(d->gedp, d->session) != GED_EDIT_OK)
	emit errorMessage(tr("Unable to cancel edit session"));
}


void
QgPrimitiveEdit::sessionObserver(struct ged *gedp,
    const struct ged_edit_session_event *event, void *clientData)
{
    QgPrimitiveEdit *self = static_cast<QgPrimitiveEdit *>(clientData);
    if (!self || !event || self->d->gedp != gedp)
	return;
    const enum ged_edit_session_event_kind kind = event->kind;
    const ged_edit_session_ref session = event->session;
    const QString path = QString::fromUtf8(event->path ? event->path : "");
    const QString replacementPath = QString::fromUtf8(
	event->replacement_path ? event->replacement_path : "");
    const enum ged_edit_session_invalidation_reason invalidationReason =
	event->invalidation_reason;
    const uint64_t revision = event->revision;
    if (self->thread() == QThread::currentThread()) {
	self->handleSessionEvent(kind, session, path, replacementPath,
	    invalidationReason, revision);
	return;
    }
    QMetaObject::invokeMethod(self,
	[self, kind, session, path, replacementPath, invalidationReason,
	    revision]() {
	self->handleSessionEvent(kind, session, path, replacementPath,
	    invalidationReason, revision);
    }, Qt::QueuedConnection);
}


void
QgPrimitiveEdit::handleSessionEvent(enum ged_edit_session_event_kind kind,
    ged_edit_session_ref session, const QString &eventPath,
    const QString &replacementPath,
    enum ged_edit_session_invalidation_reason invalidationReason,
    uint64_t revision)
{
    bool relevant = same_session(session, d->session);
    if (!relevant && d->gedp && same_target_path(targetPath(), eventPath)) {
	/* The event already carries GED's canonical session identity.  Adopt a
	 * matching externally-created session directly instead of resolving the
	 * target through the database for every unrelated notification.  The
	 * latter is both expensive and wrong after erase/remove has made a
	 * persistent widget's displayed target path stale. */
	d->session = session;
	d->ownsSession = false;
	const struct rt_edit_prim_desc *sessionDescriptor = nullptr;
	if (ged_edit_session_descriptor_get(d->gedp, session,
		&sessionDescriptor) != GED_EDIT_OK)
	    sessionDescriptor = nullptr;
	if (d->descriptor != sessionDescriptor) {
	    d->descriptor = sessionDescriptor;
	    rebuildOperations();
	}
	relevant = true;
    }
    if (!relevant)
	return;

    d->revision = revision;
    /* An external command is also a change of the active edit operation.
     * Follow it before notifying scene presenters so every frontend observes
     * one coherent revision: operation selector, parameter widgets, CLI
     * readback, and retained scene controls all describe the same command. */
    if (kind == GED_EDIT_SESSION_BEGIN || kind == GED_EDIT_SESSION_UPDATE)
	selectSessionOperation();
    emit sessionEvent(static_cast<int>(kind),
	static_cast<qulonglong>(revision));
    if (kind == GED_EDIT_SESSION_COMMIT || kind == GED_EDIT_SESSION_CANCEL ||
	kind == GED_EDIT_SESSION_INVALIDATE) {
	d->session = GED_EDIT_SESSION_REF_NULL;
	d->ownsSession = false;
	if (kind == GED_EDIT_SESSION_COMMIT)
	    d->status->setText(tr("Edit committed"));
	else if (kind == GED_EDIT_SESSION_CANCEL)
	    d->status->setText(tr("Edit canceled"));
	else
	    d->status->setText(tr("Edit invalidated by a database change"));
	if (kind == GED_EDIT_SESSION_INVALIDATE && !replacementPath.isEmpty() &&
	    d->target->text() != replacementPath) {
	    QSignalBlocker blocker(d->target);
	    d->target->setText(replacementPath);
	}
	const bool invalidationCanRestart =
	    invalidationReason == GED_EDIT_INVALIDATION_SOURCE_CHANGED ||
	    invalidationReason == GED_EDIT_INVALIDATION_SOURCE_RENAMED ||
	    invalidationReason == GED_EDIT_INVALIDATION_DATABASE_REBUILT;
	if (kind == GED_EDIT_SESSION_INVALIDATE && !invalidationCanRestart) {
	    d->descriptor = nullptr;
	    d->command = nullptr;
	    rebuildOperations();
	}
	/* A persistent editor follows the now-current database state.  Delay the
	 * replacement session until the closing event dispatch is complete. */
	if (!d->restarting && d->gedp && !targetPath().isEmpty() &&
	    (kind != GED_EDIT_SESSION_INVALIDATE || invalidationCanRestart)) {
	    d->restarting = true;
	    QTimer::singleShot(0, this, [this]() {
		d->restarting = false;
		attachTarget();
	    });
	}
	return;
    }
    refreshFromSession();
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
