/*          T E S T _ Q T C A D _ P R I M I T I V E _ E D I T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include <cmath>
#include <cstdio>

#include <QApplication>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QGroupBox>
#include <QLineEdit>
#include <QPushButton>

#include "bu/app.h"
#include "bu/file.h"
#include "bu/malloc.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "ged.h"
#include "qtcad/QgPrimitiveEdit.h"
#include "qtcad/QgSketchEdit.h"
#include "rt/edit.h"
#include "rt/geom.h"
#include "rt/primitives/sketch.h"
#include "wdb.h"

#define CHECK(_condition, _message) do { \
    if (!(_condition)) { \
	std::fprintf(stderr, "FAIL: %s\n", (_message)); \
	return 1; \
    } \
} while (0)


static int
make_edit_db(const char *path)
{
    struct rt_wdb *wdbp = wdb_fopen(path);
    if (!wdbp)
	return 0;
    point_t vertex = {0.0, 0.0, 0.0};
    vect_t axisA = {2.0, 0.0, 0.0};
    vect_t axisB = {0.0, 3.0, 0.0};
    vect_t axisC = {0.0, 0.0, 4.0};
    int ok = mk_ell(wdbp, "ell.s", vertex, axisA, axisB, axisC) == 0;

    vect_t hypH = {0.0, 0.0, 10.0};
    vect_t hypA = {4.0, 0.0, 0.0};
    ok = ok && mk_hyp(wdbp, "hyp.s", vertex, hypH, hypA, 3.0, 0.5) == 0;

    point_t arb[8] = {
	{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
	{1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
	{0.0, 0.0, 1.0}, {1.0, 0.0, 1.0},
	{1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}
    };
    ok = ok && mk_arb8(wdbp, "arb8.s", &arb[0][0]) == 0;

    struct rt_sketch_internal *sketch = nullptr;
    BU_ALLOC(sketch, struct rt_sketch_internal);
    sketch->magic = RT_SKETCH_INTERNAL_MAGIC;
    VSET(sketch->V, 1.0, 2.0, 3.0);
    VSET(sketch->u_vec, 2.0, 0.0, 0.0);
    VSET(sketch->v_vec, 1.0, 4.0, 0.0);
    sketch->vert_count = 4;
    sketch->verts = static_cast<point2d_t *>(bu_malloc(
	4 * sizeof(point2d_t), "qtcad session sketch vertices"));
    V2SET(sketch->verts[0], 0.0, 0.0);
    V2SET(sketch->verts[1], 10.0, 0.0);
    V2SET(sketch->verts[2], 10.0, 10.0);
    V2SET(sketch->verts[3], 0.0, 10.0);
    struct line_seg *line = nullptr;
    BU_ALLOC(line, struct line_seg);
    line->magic = CURVE_LSEG_MAGIC;
    line->start = 0;
    line->end = 1;
    sketch->curve.count = 1;
    sketch->curve.segment = static_cast<void **>(bu_malloc(sizeof(void *),
	"qtcad session sketch segments"));
    sketch->curve.reverse = static_cast<int *>(bu_calloc(1, sizeof(int),
	"qtcad session sketch reverse"));
    sketch->curve.segment[0] = line;
    ok = ok && wdb_export(wdbp, "sketch.s", sketch, ID_SKETCH, 1.0) == 0;

    struct rt_annot_internal annot = {};
    annot.magic = RT_ANNOT_INTERNAL_MAGIC;
    VSET(annot.V, 1.0, 2.0, 0.0);
    annot.vert_count = 1;
    annot.verts = static_cast<point2d_t *>(
	bu_calloc(1, sizeof(point2d_t), "annotation vertices"));
    V2SET(annot.verts[0], 0.1, 0.2);
    annot.ant.count = 1;
    annot.ant.reverse = static_cast<int *>(
	bu_calloc(1, sizeof(int), "annotation reverse"));
    annot.ant.segments = static_cast<void **>(
	bu_calloc(1, sizeof(void *), "annotation segments"));
    struct txt_seg *segment = nullptr;
    BU_ALLOC(segment, struct txt_seg);
    segment->magic = ANN_TSEG_MAGIC;
    segment->ref_pt = 0;
    segment->rel_pos = RT_TXT_POS_BL;
    bu_vls_init(&segment->label);
    bu_vls_strcpy(&segment->label, "Initial label");
    segment->txt_size = 10.0;
    segment->txt_rot_angle = 0.0;
    annot.ant.segments[0] = segment;
    ok = ok && mk_annot(wdbp, "annot.s", &annot) == 0;
    bu_vls_free(&segment->label);
    BU_PUT(segment, struct txt_seg);
    bu_free(annot.ant.segments, "annotation segments");
    bu_free(annot.ant.reverse, "annotation reverse");
    bu_free(annot.verts, "annotation vertices");
    wdb_close(wdbp);
    return ok;
}


static fastf_t
database_axis_a(struct ged *gedp)
{
    struct directory *dp = db_lookup(gedp->dbip, "ell.s", LOOKUP_QUIET);
    struct rt_db_internal intern = RT_DB_INTERNAL_INIT_ZERO;
    if (!dp || rt_db_get_internal(&intern, dp, gedp->dbip, nullptr) < 0)
	return -1.0;
    struct rt_ell_internal *ell =
	static_cast<struct rt_ell_internal *>(intern.idb_ptr);
    const fastf_t magnitude = MAGNITUDE(ell->a);
    rt_db_free_internal(&intern);
    return magnitude;
}


int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    QApplication app(argc, argv);

    const char *dbPath = "qtcad_primitive_edit_tmp.g";
    (void)bu_file_delete(dbPath);
    CHECK(make_edit_db(dbPath), "could not create primitive-edit database");
    struct ged *gedp = ged_open("db", dbPath, 1);
    CHECK(gedp != nullptr, "could not open primitive-edit database");

    {
	ged_edit_session_ref session = GED_EDIT_SESSION_REF_NULL;
	CHECK(ged_edit_session_begin(gedp, "sketch.s", nullptr, &session) ==
	    GED_EDIT_OK, "could not begin the sketch interaction session");
	QgSketchEdit sketch;
	CHECK(sketch.refresh(gedp, session, QStringLiteral("sketch.s")) &&
	    sketch.vertexCount() == 4 && sketch.segmentCount() == 1 &&
	    sketch.linePointCount() >= 2,
	    "session-backed sketch topology did not initialize");
	point_t point;
	point_t expected = {31.0, 42.0, 3.0};
	CHECK(sketch.vertexDisplayPoint(2, point) &&
	    VNEAR_EQUAL(point, expected, SMALL_FASTF),
	    "sketch occurrence presentation lost its non-unit UV basis");
	fastf_t u = 0.0;
	fastf_t v = 0.0;
	CHECK(sketch.displayPointToUV(point, &u, &v) &&
	    NEAR_EQUAL(u, 10.0, SMALL_FASTF) &&
	    NEAR_EQUAL(v, 10.0, SMALL_FASTF),
	    "display-to-sketch conversion did not invert the UV basis");

	CHECK(sketch.selectFeature(QgSketchEdit::Vertex, 2) == GED_EDIT_OK,
	    "retained interaction could not select a sketch vertex");
	point_t moved = {63.0, 90.0, 3.0};
	CHECK(sketch.moveVertexToDisplayPoint(moved) == GED_EDIT_OK &&
	    sketch.refresh(gedp, session, QStringLiteral("sketch.s")) &&
	    sketch.vertexDisplayPoint(2, point) &&
	    VNEAR_EQUAL(point, moved, SMALL_FASTF),
	    "retained vertex drag did not update the authoritative session");

	point_t added = {101.0, 162.0, 3.0};
	CHECK(sketch.addVertexAtDisplayPoint(added) == GED_EDIT_OK &&
	    sketch.refresh(gedp, session, QStringLiteral("sketch.s")) &&
	    sketch.vertexCount() == 5 && sketch.vertexDisplayPoint(4, point) &&
	    VNEAR_EQUAL(point, added, SMALL_FASTF),
	    "retained add-vertex interaction did not update topology");

	CHECK(sketch.selectFeature(QgSketchEdit::Segment, 0) == GED_EDIT_OK,
	    "retained interaction could not select a sketch segment");
	vect_t delta = {13.0, 12.0, 0.0};
	CHECK(sketch.moveSegmentByDisplayDelta(delta) == GED_EDIT_OK &&
	    sketch.refresh(gedp, session, QStringLiteral("sketch.s")),
	    "retained segment drag did not update the authoritative session");
	point_t firstExpected = {14.0, 14.0, 3.0};
	point_t secondExpected = {34.0, 14.0, 3.0};
	CHECK(sketch.vertexDisplayPoint(0, point) &&
	    VNEAR_EQUAL(point, firstExpected, SMALL_FASTF) &&
	    sketch.vertexDisplayPoint(1, point) &&
	    VNEAR_EQUAL(point, secondExpected, SMALL_FASTF),
	    "segment drag used display coordinates as raw sketch UV values");
	CHECK(ged_edit_session_cancel(gedp, session) == GED_EDIT_OK,
	    "could not cancel the sketch interaction session");
    }

    {
	QgPrimitiveEdit editor;
	editor.setGed(gedp);
	editor.setTargetPath(QStringLiteral("hyp.s"));
	QComboBox *operation = editor.findChild<QComboBox *>(
	    QStringLiteral("primitiveEdit.operation"));
	QDoubleSpinBox *height = editor.findChild<QDoubleSpinBox *>(
	    QStringLiteral("primitiveEdit.param.h.0"));
	QPushButton *apply = editor.findChild<QPushButton *>(
	    QStringLiteral("primitiveEdit.apply"));
	CHECK(operation && height && apply &&
	    operation->currentText() == QStringLiteral("Set H") &&
	    std::fabs(height->value() - 10.0) < 1.0e-9,
	    "generated HYP editor initializes from absolute descriptor values");
	height->setValue(12.0);
	apply->click();
	const char *getHeight[] = {"edit", "hyp.s", "get", "h", nullptr};
	CHECK(ged_exec(gedp, 4, getHeight) == BRLCAD_OK &&
	    QString::fromUtf8(bu_vls_cstr(gedp->ged_result_str)).contains(
		QStringLiteral("12")),
	    "CLI observes an absolute HYP update made by the generated widget");
	const char *setB[] = {"edit", "-i", "hyp.s", "b", "5", nullptr};
	CHECK(ged_exec(gedp, 5, setB) == BRLCAD_OK &&
	    operation->currentText() == QStringLiteral("Set B"),
	    "HYP command operation change updates the generated widget");
	QDoubleSpinBox *axisB = editor.findChild<QDoubleSpinBox *>(
	    QStringLiteral("primitiveEdit.param.b.0"));
	CHECK(axisB && std::fabs(axisB->value() - 5.0) < 1.0e-9,
	    "HYP command value refreshes the generated widget absolutely");
	editor.cancel();
	QApplication::processEvents();
	editor.setGed(nullptr);
    }

    {
	const char *copy[] = {"cp", "ell.s", "resurrect.s", nullptr};
	CHECK(ged_exec(gedp, 3, copy) == BRLCAD_OK,
	    "could not create the recreate-after-remove fixture");
	QgPrimitiveEdit editor;
	editor.setGed(gedp);
	editor.setTargetPath(QStringLiteral("resurrect.s"));
	const char *remove[] = {"kill", "resurrect.s", nullptr};
	CHECK(ged_exec(gedp, 2, remove) == BRLCAD_OK,
	    "could not remove the recreate-after-remove fixture");
	QApplication::processEvents();
	CHECK(ged_edit_session_ref_is_null(editor.session()),
	    "removed target left a live generated-editor session");
	CHECK(ged_exec(gedp, 3, copy) == BRLCAD_OK,
	    "could not recreate a removed editor target");
	const char *cliSet[] = {
	    "edit", "-i", "resurrect.s", "a", "13", nullptr
	};
	CHECK(ged_exec(gedp, 5, cliSet) == BRLCAD_OK,
	    "CLI could not begin an edit for a recreated widget target");
	QApplication::processEvents();
	QDoubleSpinBox *axisA = editor.findChild<QDoubleSpinBox *>(
	    QStringLiteral("primitiveEdit.param.a.0"));
	CHECK(!ged_edit_session_ref_is_null(editor.session()) && axisA &&
	    std::fabs(axisA->value() - 13.0) < 1.0e-9,
	    "widget did not adopt and display an externally begun replacement session");
	const char *cliGet[] = {
	    "edit", "resurrect.s", "get", "a", nullptr
	};
	CHECK(ged_exec(gedp, 4, cliGet) == BRLCAD_OK &&
	    QString::fromUtf8(bu_vls_cstr(gedp->ged_result_str)).contains(
		QStringLiteral("13")),
	    "command readback disagreed with the adopted replacement session");
	editor.cancel();
	editor.setGed(nullptr);
	const char *removeCopy[] = {"kill", "resurrect.s", nullptr};
	CHECK(ged_exec(gedp, 2, removeCopy) == BRLCAD_OK,
	    "could not remove the recreate-after-remove fixture copy");
    }

    {
	QgPrimitiveEdit editor;
	editor.setGed(gedp);
	editor.setTargetPath(QStringLiteral("ell.s"));
	CHECK(!ged_edit_session_ref_is_null(editor.session()),
	    "generated editor did not attach to an edit session");

	QComboBox *operation = editor.findChild<QComboBox *>(
	    QStringLiteral("primitiveEdit.operation"));
	QDoubleSpinBox *axisA = editor.findChild<QDoubleSpinBox *>(
	    QStringLiteral("primitiveEdit.param.a.0"));
	QPushButton *apply = editor.findChild<QPushButton *>(
	    QStringLiteral("primitiveEdit.apply"));
	QPushButton *commit = editor.findChild<QPushButton *>(
	    QStringLiteral("primitiveEdit.commit"));
	CHECK(operation && axisA && apply && commit,
	    "generated editor did not expose its stable controls");
	CHECK(operation->currentText() == QStringLiteral("Set A") &&
	    std::fabs(axisA->value() - 2.0) < 1.0e-9,
	    "generated editor did not initialize from the shared session");

	QgPrimitiveEdit peerEditor;
	peerEditor.setGed(gedp);
	peerEditor.setTargetPath(QStringLiteral("ell.s"));
	QComboBox *peerOperation = peerEditor.findChild<QComboBox *>(
	    QStringLiteral("primitiveEdit.operation"));
	QDoubleSpinBox *peerAxisA = peerEditor.findChild<QDoubleSpinBox *>(
	    QStringLiteral("primitiveEdit.param.a.0"));
	QPushButton *peerApply = peerEditor.findChild<QPushButton *>(
	    QStringLiteral("primitiveEdit.apply"));
	CHECK(peerOperation && peerAxisA && peerApply &&
	    peerEditor.session().id == editor.session().id,
	    "a second generated editor joins the same authoritative session");

	const char *cliSet[] = {"edit", "-i", "ell.s", "a", "7", nullptr};
	CHECK(ged_exec(gedp, 5, cliSet) == BRLCAD_OK,
	    "CLI could not update the editor's shared session");
	CHECK(std::fabs(axisA->value() - 7.0) < 1.0e-9 &&
	    std::fabs(peerAxisA->value() - 7.0) < 1.0e-9,
	    "command-originated edit did not refresh every generated widget");

	const char *cliSetB[] = {"edit", "-i", "ell.s", "b", "8", nullptr};
	CHECK(ged_exec(gedp, 5, cliSetB) == BRLCAD_OK &&
	    operation->currentText() == QStringLiteral("Set B") &&
	    peerOperation->currentText() == QStringLiteral("Set B"),
	    "command-originated operation did not update every widget selector");
	QDoubleSpinBox *axisB = editor.findChild<QDoubleSpinBox *>(
	    QStringLiteral("primitiveEdit.param.b.0"));
	QDoubleSpinBox *peerAxisB = peerEditor.findChild<QDoubleSpinBox *>(
	    QStringLiteral("primitiveEdit.param.b.0"));
	CHECK(axisB && peerAxisB && std::fabs(axisB->value() - 8.0) < 1.0e-9 &&
	    std::fabs(peerAxisB->value() - 8.0) < 1.0e-9,
	    "command-originated operation did not rebuild all parameter widgets");
	axisB->setValue(10.0);
	apply->click();
	const char *cliGetB[] = {"edit", "ell.s", "get", "b", nullptr};
	CHECK(ged_exec(gedp, 4, cliGetB) == BRLCAD_OK &&
	    QString::fromUtf8(bu_vls_cstr(gedp->ged_result_str)).contains(
		QStringLiteral("10")) &&
	    std::fabs(peerAxisB->value() - 10.0) < 1.0e-9,
	    "CLI and peer widget did not observe a widget-originated update");
	peerAxisB->setValue(10.5);
	peerApply->click();
	CHECK(std::fabs(axisB->value() - 10.5) < 1.0e-9 &&
	    ged_exec(gedp, 4, cliGetB) == BRLCAD_OK &&
	    QString::fromUtf8(bu_vls_cstr(gedp->ged_result_str)).contains(
		QStringLiteral("10.5")),
	    "primary widget and CLI did not observe a peer-widget update");
	peerEditor.setGed(nullptr);
	operation->setCurrentIndex(0);
	axisA = editor.findChild<QDoubleSpinBox *>(
	    QStringLiteral("primitiveEdit.param.a.0"));
	CHECK(axisA && std::fabs(axisA->value() - 7.0) < 1.0e-9,
	    "widget operation switch did not restore authoritative session state");

	axisA->setValue(9.0);
	apply->click();
	const char *cliGet[] = {"edit", "ell.s", "get", "a", nullptr};
	CHECK(ged_exec(gedp, 4, cliGet) == BRLCAD_OK &&
	    QString::fromUtf8(bu_vls_cstr(gedp->ged_result_str)).contains(
		QStringLiteral("9")),
	    "CLI did not observe a generated-widget update");

	const char *checkpoint[] = {"edit", "ell.s", "checkpoint", nullptr};
	const char *setEleven[] = {"edit", "-i", "ell.s", "a", "11", nullptr};
	const char *revert[] = {"edit", "-i", "ell.s", "revert", nullptr};
	CHECK(ged_exec(gedp, 3, checkpoint) == BRLCAD_OK &&
	    ged_exec(gedp, 5, setEleven) == BRLCAD_OK &&
	    std::fabs(axisA->value() - 11.0) < 1.0e-9 &&
	    ged_exec(gedp, 4, revert) == BRLCAD_OK &&
	    std::fabs(axisA->value() - 9.0) < 1.0e-9,
	    "checkpoint/revert did not stay synchronized across CLI and widget");

	commit->click();
	CHECK(std::fabs(database_axis_a(gedp) - 9.0) < 1.0e-9,
	    "widget commit did not write the shared session state");
	QApplication::processEvents();
	CHECK(!ged_edit_session_ref_is_null(editor.session()) &&
	    std::fabs(axisA->value() - 9.0) < 1.0e-9,
	    "persistent editor did not reattach to committed database state");

	const char *adjust[] = {"adjust", "ell.s", "A", "12 0 0", nullptr};
	CHECK(ged_exec(gedp, 4, adjust) == BRLCAD_OK,
	    "external CLI adjust could not replace the editor source geometry");
	QApplication::processEvents();
	CHECK(!ged_edit_session_ref_is_null(editor.session()) &&
	    editor.targetPath() == QStringLiteral("ell.s") &&
	    std::fabs(axisA->value() - 12.0) < 1.0e-9,
	    "persistent editor did not reload after external geometry replacement");

	const char *rename[] = {"move", "ell.s", "ell_renamed.s", nullptr};
	CHECK(ged_exec(gedp, 3, rename) == BRLCAD_OK,
	    "external CLI rename could not rename the editor source");
	QApplication::processEvents();
	CHECK(editor.targetPath() == QStringLiteral("/ell_renamed.s") &&
	    !ged_edit_session_ref_is_null(editor.session()) &&
	    std::fabs(axisA->value() - 12.0) < 1.0e-9,
	    "persistent editor did not follow the rename replacement path");
	const char *renamedGet[] = {
	    "edit", "ell_renamed.s", "get", "a", nullptr
	};
	CHECK(ged_exec(gedp, 4, renamedGet) == BRLCAD_OK &&
	    QString::fromUtf8(bu_vls_cstr(gedp->ged_result_str)).contains(
		QStringLiteral("12")),
	    "CLI did not observe the reattached renamed editor state");

	const char *renameBack[] = {"move", "ell_renamed.s", "ell.s", nullptr};
	CHECK(ged_exec(gedp, 3, renameBack) == BRLCAD_OK,
	    "external CLI rename could not restore the editor source name");
	QApplication::processEvents();
	CHECK(editor.targetPath() == QStringLiteral("/ell.s") &&
	    !ged_edit_session_ref_is_null(editor.session()),
	    "persistent editor did not follow a second rename");
	editor.cancel();
	editor.setGed(nullptr);
    }

    {
	QgPrimitiveEdit editor;
	editor.setGed(gedp);
	editor.setTargetPath(QStringLiteral("ell.s"));
	const ged_edit_session_ref owned = editor.session();
	CHECK(!ged_edit_session_ref_is_null(owned),
	    "target-switch test begins a widget-owned session");
	QDoubleSpinBox *axisA = editor.findChild<QDoubleSpinBox *>(
	    QStringLiteral("primitiveEdit.param.a.0"));
	QPushButton *apply = editor.findChild<QPushButton *>(
	    QStringLiteral("primitiveEdit.apply"));
	CHECK(axisA && apply, "target-switch test finds generated controls");
	axisA->setValue(14.0);
	apply->click();
	editor.setTargetPath(QStringLiteral("annot.s"));
	struct ged_edit_session_info oldInfo = {};
	CHECK(ged_edit_session_info_get(gedp, owned, &oldInfo) == GED_EDIT_STALE &&
	    editor.targetPath() == QStringLiteral("annot.s") &&
	    !ged_edit_session_ref_is_null(editor.session()),
	    "switching targets cancels only the widget-owned intermediate");
	CHECK(std::fabs(database_axis_a(gedp) - 12.0) < 1.0e-9,
	    "target switching does not commit abandoned widget geometry");
	editor.cancel();
    }

    {
	ged_edit_session_ref cliSession = GED_EDIT_SESSION_REF_NULL;
	CHECK(ged_edit_session_begin(gedp, "ell.s", nullptr, &cliSession) ==
	    GED_EDIT_OK, "CLI-owned target-switch fixture begins a session");
	QgPrimitiveEdit editor;
	editor.setGed(gedp);
	editor.setTargetPath(QStringLiteral("ell.s"));
	CHECK(editor.session().id == cliSession.id,
	    "generated editor joins an independently owned session");
	editor.setTargetPath(QStringLiteral("annot.s"));
	struct ged_edit_session_info cliInfo = {};
	CHECK(ged_edit_session_info_get(gedp, cliSession, &cliInfo) == GED_EDIT_OK,
	    "switching widgets leaves an independently owned CLI session active");
	editor.cancel();
	CHECK(ged_edit_session_cancel(gedp, cliSession) == GED_EDIT_OK,
	    "CLI-owned session remains explicitly controllable by its owner");
    }

    {
	QgPrimitiveEdit editor;
	editor.setGed(gedp);
	editor.setTargetPath(QStringLiteral("annot.s"));
	QComboBox *operation = editor.findChild<QComboBox *>(
	    QStringLiteral("primitiveEdit.operation"));
	QLineEdit *text = editor.findChild<QLineEdit *>(
	    QStringLiteral("primitiveEdit.param.text"));
	QPushButton *apply = editor.findChild<QPushButton *>(
	    QStringLiteral("primitiveEdit.apply"));
	CHECK(operation && text && apply &&
	    operation->currentText() == QStringLiteral("Set Text") &&
	    text->text() == QStringLiteral("Initial label"),
	    "generated string editor did not initialize from session state");

	const char *cliSet[] = {
	    "edit", "-i", "annot.s", "set_text", "CLI label", nullptr
	};
	CHECK(ged_exec(gedp, 5, cliSet) == BRLCAD_OK &&
	    text->text() == QStringLiteral("CLI label"),
	    "CLI string update did not refresh the generated widget");

	text->setText(QStringLiteral("Widget label"));
	apply->click();
	struct rt_edit_cmd_values current;
	CHECK(ged_edit_session_command_values_get(gedp, editor.session(),
	    42010, &current) == GED_EDIT_OK && current.string_valid[0] &&
	    bu_strcmp(current.strings[0], "Widget label") == 0,
	    "shared session did not retain the generated-widget string update");
	const char *cliGet[] = {
	    "edit", "annot.s", "get", "set_text", nullptr
	};
	CHECK(ged_exec(gedp, 4, cliGet) == BRLCAD_OK &&
	    QString::fromUtf8(bu_vls_cstr(gedp->ged_result_str)).contains(
		QStringLiteral("Widget")) &&
	    QString::fromUtf8(bu_vls_cstr(gedp->ged_result_str)).contains(
		QStringLiteral("label")),
	    "CLI did not observe a generated-widget string update");
	const char *remove[] = {"kill", "annot.s", nullptr};
	CHECK(ged_exec(gedp, 2, remove) == BRLCAD_OK,
	    "external CLI remove could not delete the editor source");
	QApplication::processEvents();
	CHECK(ged_edit_session_ref_is_null(editor.session()) &&
	    operation->count() == 0,
	    "persistent editor retained controls or a session for a removed source");
	editor.setGed(nullptr);
    }

    {
	QgPrimitiveEdit editor;
	editor.setGed(gedp);
	editor.setTargetPath(QStringLiteral("arb8.s"));
	QComboBox *operation = editor.findChild<QComboBox *>(
	    QStringLiteral("primitiveEdit.operation"));
	QGroupBox *parameters = editor.findChild<QGroupBox *>();
	QPushButton *apply = editor.findChild<QPushButton *>(
	    QStringLiteral("primitiveEdit.apply"));
	CHECK(operation && parameters && apply && operation->count() == 7 &&
	    operation->currentText() == QStringLiteral("Select Vertex") &&
	    !parameters->isEnabled() && !apply->isEnabled(),
	    "custom ARB operations are visible state but not unsafe generated controls");
	const char *selectEdge[] = {
	    "edit", "-i", "arb8.s", "select_edge", "0", nullptr
	};
	CHECK(ged_exec(gedp, 5, selectEdge) == BRLCAD_OK &&
	    operation->currentText() == QStringLiteral("Select Edge") &&
	    !parameters->isEnabled() && !apply->isEnabled(),
	    "command-originated ARB selection updates the custom widget state");
	const char *moveVertex[] = {
	    "edit", "-i", "arb8.s", "move_vertex", "0", "2", "2", "2",
	    nullptr
	};
	CHECK(ged_exec(gedp, 8, moveVertex) == BRLCAD_OK &&
	    operation->currentText() == QStringLiteral("Move Vertex") &&
	    !parameters->isEnabled() && !apply->isEnabled(),
	    "command-originated custom edit updates the widget without enabling a false form");
	editor.cancel();
	editor.setGed(nullptr);
    }

    ged_close(gedp);
    (void)bu_file_delete(dbPath);
    std::fprintf(stdout,
	"PASS: generated/custom primitive editor CLI-widget synchronization\n");
    return 0;
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
