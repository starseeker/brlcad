/*               T E S T _ Q T C A D _ E V E N T _ P L A Y E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/app.h"
#include "qtcad/QgTestEvent.h"

#include <QApplication>
#include <QComboBox>
#include <QJsonArray>
#include <QLineEdit>
#include <QMouseEvent>
#include <QPushButton>
#include <QStandardItemModel>
#include <QTemporaryDir>
#include <QTreeView>
#include <QVBoxLayout>

#include <cstdio>

#define CHECK(_condition, _message) do { \
    if (!(_condition)) { \
	std::fprintf(stderr, "FAIL: %s\n", (_message)); \
	return 1; \
    } \
} while (0)

class CaptureWidget : public QWidget {
public:
    QPointF lastPress;
    int pressCount = 0;

protected:
    void mousePressEvent(QMouseEvent *event) override
    {
#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
	lastPress = event->position();
#else
	lastPress = event->localPos();
#endif
	pressCount++;
	QWidget::mousePressEvent(event);
    }
};

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    QApplication app(argc, argv);

    QWidget root;
    root.setObjectName(QStringLiteral("testRoot"));
    QVBoxLayout layout(&root);
    QPushButton button(QStringLiteral("Apply"), &root);
    button.setObjectName(QStringLiteral("applyButton"));
    button.setProperty("qgTestId", QStringLiteral("apply-action"));
    QLineEdit edit(&root);
    edit.setObjectName(QStringLiteral("nameEdit"));
    QComboBox combo(&root);
    combo.setObjectName(QStringLiteral("modeCombo"));
    combo.addItems({QStringLiteral("Wire"), QStringLiteral("Shaded"),
	QStringLiteral("Hidden line")});
    QTreeView tree(&root);
    tree.setObjectName(QStringLiteral("sceneTree"));
    QStandardItemModel model;
    QStandardItem *assembly = new QStandardItem(QStringLiteral("assembly"));
    assembly->appendRow(new QStandardItem(QStringLiteral("wheel")));
    model.appendRow(assembly);
    tree.setModel(&model);
    CaptureWidget canvas;
    canvas.setObjectName(QStringLiteral("canvas"));
    canvas.setFixedSize(200, 100);
    canvas.setParent(&root);
    layout.addWidget(&button);
    layout.addWidget(&edit);
    layout.addWidget(&combo);
    layout.addWidget(&tree);
    layout.addWidget(&canvas);

    int clicks = 0;
    QObject::connect(&button, &QPushButton::clicked, [&clicks]() { clicks++; });

    QgEventRecorder recorder(&root);
    recorder.start();
    button.click();
    edit.setText(QStringLiteral("Generic_Twin"));
    QMetaObject::invokeMethod(&edit, "editingFinished", Qt::DirectConnection);
    combo.setCurrentIndex(2);
    const QModelIndex wheel = model.index(0, 0, model.index(0, 0));
    tree.setCurrentIndex(wheel);
    QCoreApplication::processEvents();
    CHECK(recorder.events().size() >= 4,
	"semantic recorder did not capture button/edit/combo/tree changes");

    const QString buttonPath = QgEventRecorder::objectPath(&root, &button);
    CHECK(buttonPath.contains(QLatin1String("i:apply-action")) &&
	QgEventRecorder::resolveObject(&root, buttonPath) == &button,
	"stable QObject path did not round-trip");
    CHECK(QgEventRecorder::resolveObject(&root,
	QStringLiteral("i:apply-action")) == &button,
	"unique global test-id lookup did not resolve the control");

    QTemporaryDir temporary;
    CHECK(temporary.isValid(), "could not create temporary event directory");
    const QString stream = temporary.filePath(QStringLiteral("events.json"));
    QString error;
    CHECK(recorder.save(stream, &error), error.toLocal8Bit().constData());
    QVector<QgTestEvent> loaded;
    CHECK(QgEventRecorder::load(stream, loaded, &error),
	error.toLocal8Bit().constData());
    CHECK(loaded.size() == recorder.events().size(),
	"saved event stream did not round-trip");

    recorder.stop();
    clicks = 0;
    edit.clear();
    combo.setCurrentIndex(0);
    tree.setCurrentIndex(QModelIndex());
    QgEventPlayer player(&root);
    CHECK(player.play(loaded, &error), error.toLocal8Bit().constData());
    CHECK(clicks == 1, "player did not activate the recorded button exactly once");
    CHECK(edit.text() == QLatin1String("Generic_Twin"),
	"player did not restore line-edit text");
    CHECK(combo.currentIndex() == 2, "player did not restore combo index");
    CHECK(tree.currentIndex() == wheel, "player did not restore tree current item");

    QgTestEvent assertCombo;
    assertCombo.target = QgEventRecorder::objectPath(&root, &combo);
    assertCombo.action = QStringLiteral("assert_state");
    assertCombo.arguments.insert(QStringLiteral("text"),
	QStringLiteral("Hidden line"));
    assertCombo.arguments.insert(QStringLiteral("index"), 2);
    assertCombo.arguments.insert(QStringLiteral("count"), 3);
    assertCombo.arguments.insert(QStringLiteral("enabled"), true);
    assertCombo.arguments.insert(QStringLiteral("items"), QJsonArray{
	QStringLiteral("Wire"), QStringLiteral("Shaded"),
	QStringLiteral("Hidden line")});
    CHECK(player.play(assertCombo, &error), error.toLocal8Bit().constData());

    QgTestEvent assertEdit;
    assertEdit.target = QgEventRecorder::objectPath(&root, &edit);
    assertEdit.action = QStringLiteral("assert_state");
    assertEdit.arguments.insert(QStringLiteral("text"),
	QStringLiteral("Generic_Twin"));
    CHECK(player.play(assertEdit, &error), error.toLocal8Bit().constData());

    const QString treePath = QgEventRecorder::objectPath(&root, &tree);
    QgTestEvent collapse;
    collapse.target = treePath;
    collapse.action = QStringLiteral("set_expanded");
    collapse.arguments.insert(QStringLiteral("rows"), QJsonArray{0});
    collapse.arguments.insert(QStringLiteral("expanded"), false);
    CHECK(player.play(collapse, &error) && !tree.isExpanded(model.index(0, 0)),
	"player did not collapse a tree item");

    QgTestEvent expand = collapse;
    expand.arguments.remove(QStringLiteral("rows"));
    expand.arguments.insert(QStringLiteral("labels"),
	QJsonArray{QStringLiteral("assembly")});
    expand.arguments.insert(QStringLiteral("expanded"), true);
    CHECK(player.play(expand, &error) && tree.isExpanded(model.index(0, 0)),
	"player did not expand a tree item resolved by label");

    QgTestEvent clearSelection;
    clearSelection.target = treePath;
    clearSelection.action = QStringLiteral("clear_selection");
    CHECK(player.play(clearSelection, &error),
	error.toLocal8Bit().constData());
    CHECK(tree.selectionModel()->selectedIndexes().isEmpty() &&
	!tree.currentIndex().isValid(),
	"player did not clear tree selection and current item");

    QgTestEvent resize;
    resize.target = QStringLiteral(".");
    resize.action = QStringLiteral("resize");
    resize.arguments.insert(QStringLiteral("width"), 640);
    resize.arguments.insert(QStringLiteral("height"), 480);
    CHECK(player.play(resize, &error) && root.size() == QSize(640, 480),
	"player did not deterministically resize its root widget");

    QgTestEvent windowState;
    windowState.target = QStringLiteral(".");
    windowState.action = QStringLiteral("window_state");
    windowState.arguments.insert(QStringLiteral("state"),
	QStringLiteral("minimized"));
    CHECK(player.play(windowState, &error) && root.isMinimized(),
	"player did not minimize its root widget");
    windowState.arguments.insert(QStringLiteral("state"),
	QStringLiteral("normal"));
    CHECK(player.play(windowState, &error) && !root.isMinimized() &&
	!root.isMaximized() && !root.isFullScreen(),
	"player did not restore its root widget");

    const QPointF sourcePos(50.0, 25.0);
    QMouseEvent sourcePress(QEvent::MouseButtonPress, sourcePos,
	canvas.mapToGlobal(sourcePos.toPoint()),
	Qt::LeftButton, Qt::LeftButton, Qt::NoModifier);
    QgTestEvent translated;
    CHECK(QgEventTranslator::translate(&canvas, &sourcePress, translated),
	"mouse translator rejected a widget mouse event");
    translated.target = QgEventRecorder::objectPath(&root, &canvas);
    CHECK(player.play(translated, &error), error.toLocal8Bit().constData());
    CHECK(canvas.pressCount == 1 && qAbs(canvas.lastPress.x() - 50.0) < 1.0 &&
	qAbs(canvas.lastPress.y() - 25.0) < 1.0,
	"normalized mouse coordinates did not replay at the expected location");

    bool checkpointSeen = false;
    player.setCheckpointHandler([&checkpointSeen](QWidget *, const QString &name,
	QString *) {
	checkpointSeen = name == QLatin1String("after-draw");
	return checkpointSeen;
    });
    QgTestEvent checkpoint;
    checkpoint.target = translated.target;
    checkpoint.action = QStringLiteral("checkpoint");
    checkpoint.arguments.insert(QStringLiteral("name"),
	QStringLiteral("after-draw"));
    CHECK(player.play(checkpoint, &error) && checkpointSeen,
	"checkpoint callback was not invoked");

    return 0;
}
