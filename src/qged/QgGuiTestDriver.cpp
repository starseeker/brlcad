/*                  Q G G U I T E S T D R I V E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */

#include "common.h"

#include <algorithm>
#include <cstdlib>
#include <memory>
#include <vector>

#include <QDir>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QFile>
#include <QHash>
#include <QImage>
#include <QJsonArray>
#include <QJsonDocument>
#include <QOpenGLWidget>
#include <QPixmap>
#include <QSaveFile>
#include <QTimer>

#include "BObol/BLodService.h"
#include "BObol/BDatabaseSource.h"
#include "BObol/BViewController.h"
#include "BObol/BViewLod.h"
#include "brlcad_version.h"
#include "bu/log.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "qtcad/QgCanvasBase.h"
#include "qtcad/QgGL.h"
#include "qtcad/QgSW.h"
#include "qtcad/QgView.h"
#include "QgEdApp.h"
#include "QgGuiTestDriver.h"

static std::vector<BObolViewController *>
qged_test_all_controllers(QgEdApp &app)
{
    std::vector<BObolViewController *> controllers;
    if (!app.w)
	return controllers;

    const QList<QgView *> views = app.w->findChildren<QgView *>();
    for (QgView *view : views) {
	BObolViewController *controller =
	    view ? view->obolViewController() : nullptr;
	if (controller &&
	    std::find(controllers.begin(), controllers.end(), controller) ==
		controllers.end())
	    controllers.push_back(controller);
    }
    QgView *current = app.w->CurrentDisplay();
    BObolViewController *controller =
	current ? current->obolViewController() : nullptr;
    if (controller &&
	std::find(controllers.begin(), controllers.end(), controller) ==
	    controllers.end())
	controllers.push_back(controller);
    return controllers;
}

/* A controller render request records renderer intent, but it is deliberately
 * not a Qt paint primitive.  In particular, a request made while another one
 * is already latched has no false-to-true callback edge with which to wake a
 * quiescent canvas.  Test presentation barriers need an actual endpoint frame,
 * so route the wake through the canvas abstraction without claiming a semantic
 * camera or scene change. */
static void
qged_test_present_controller_frame(QgEdApp &app,
	BObolViewController *controller)
{
    if (!app.w || app.w->isMinimized() || !app.w->isVisible() ||
	!controller)
	return;

    /* A quad/switchable display may retain hidden QgView instances which
     * share the endpoint controller.  update() on such a widget is legally
     * ignored by Qt, so prefer the active display before considering the
     * child list. */
    QgView *current = app.w->CurrentDisplay();
    if (current && current->isVisible() &&
	current->obolViewController() == controller) {
	QgCanvasBase *canvas = current->canvasBase();
	if (canvas) {
	    canvas->present_frame();
	    QWidget *widget = canvas->canvasWidget();
	    if (widget && widget->isVisible())
		widget->repaint();
	    return;
	}
    }

    const QList<QgView *> views = app.w->findChildren<QgView *>();
    for (QgView *view : views) {
	if (!view || !view->isVisible() ||
	    view->obolViewController() != controller)
	    continue;
	QgCanvasBase *canvas = view->canvasBase();
	if (canvas) {
	    canvas->present_frame();
	    QWidget *widget = canvas->canvasWidget();
	    if (widget && widget->isVisible())
		widget->repaint();
	}
	return;
    }
}

static bool
qged_test_controller_has_visible_canvas(QgEdApp &app,
	BObolViewController *controller)
{
    if (!app.w || app.w->isMinimized() || !app.w->isVisible() ||
	!controller)
	return false;

    QgView *current = app.w->CurrentDisplay();
    if (current && current->isVisible() &&
	current->obolViewController() == controller)
	return true;

    const QList<QgView *> views = app.w->findChildren<QgView *>();
    for (QgView *view : views) {
	if (view && view->isVisible() &&
	    view->obolViewController() == controller)
	    return true;
    }
    return false;
}

static bool
qged_write_test_report(const QString &fileName, const QJsonObject &report,
    QString *error)
{
    if (fileName.isEmpty())
	return true;
    QSaveFile file(fileName);
    if (!file.open(QIODevice::WriteOnly)) {
	if (error)
	    *error = file.errorString();
	return false;
    }
    const QByteArray data = QJsonDocument(report).toJson(
	QJsonDocument::Indented);
    if (file.write(data) != data.size() || !file.commit()) {
	if (error)
	    *error = file.errorString();
	return false;
    }
    return true;
}

static bool
qged_test_command_expect(QgEdApp &app, const QJsonObject &arguments,
    QHash<QString, QString> *captures, QString *actualOutput, QString *error)
{
    const QString command =
	arguments.value(QStringLiteral("command")).toString();
    if (command.isEmpty()) {
	if (error)
	    *error = QStringLiteral(
		"qged_command_expect requires a command argument");
	return false;
    }

    const QByteArray encoded = command.toLocal8Bit();
    std::vector<char> input(static_cast<size_t>(encoded.size()) + 1, '\0');
    if (!encoded.isEmpty())
	std::copy(encoded.constBegin(), encoded.constEnd(), input.begin());
    std::vector<char *> argv(input.size() + 1, nullptr);
    const int argc = static_cast<int>(bu_argv_from_string(
	argv.data(), input.size(), input.data()));
    if (argc <= 0) {
	if (error)
	    *error = QStringLiteral("qged_command_expect could not parse: %1")
		.arg(command);
	return false;
    }
    std::vector<const char *> commandArgv(static_cast<size_t>(argc), nullptr);
    for (int i = 0; i < argc; ++i)
	commandArgv[static_cast<size_t>(i)] = argv[static_cast<size_t>(i)];

    struct bu_vls result = BU_VLS_INIT_ZERO;
    const int ret = app.run_cmd(&result, argc,
	commandArgv.data());
    const QString output = QString::fromLocal8Bit(bu_vls_cstr(&result));
    bu_vls_free(&result);
    if (actualOutput)
	*actualOutput = output;
    if (ret & BRLCAD_ERROR) {
	if (error)
	    *error = QStringLiteral("qged command failed: %1 (%2)")
		.arg(command, output);
	return false;
    }

    const auto checkTerms = [&](const QString &key, bool wanted) {
	const QJsonValue value = arguments.value(key);
	QJsonArray terms;
	if (value.isString())
	    terms.append(value);
	else
	    terms = value.toArray();
	for (const QJsonValue &termValue : terms) {
	    const QString term = termValue.toString();
	    if (term.isEmpty())
		continue;
	    if (output.contains(term) != wanted) {
		if (error)
		    *error = wanted ?
			QStringLiteral(
			    "qged command output did not contain '%1': %2")
			    .arg(term, output) :
			QStringLiteral(
			    "qged command output unexpectedly contained '%1': %2")
			    .arg(term, output);
		return false;
	    }
	}
	return true;
    };
    if (!checkTerms(QStringLiteral("contains"), true) ||
	!checkTerms(QStringLiteral("not_contains"), false))
	return false;

    const QString numericKeys[] = {
	QStringLiteral("numeric_gt"), QStringLiteral("numeric_ge"),
	QStringLiteral("numeric_lt"), QStringLiteral("numeric_le")
    };
    bool numericOk = false;
    const double numericOutput = output.trimmed().toDouble(&numericOk);
    for (const QString &key : numericKeys) {
	if (!arguments.contains(key))
	    continue;
	if (!numericOk) {
	    if (error)
		*error = QStringLiteral(
		    "qged command output is not numeric for %1: %2")
		    .arg(key, output);
	    return false;
	}
	const double limit = arguments.value(key).toDouble();
	const bool pass = key == QLatin1String("numeric_gt") ?
	    numericOutput > limit : key == QLatin1String("numeric_ge") ?
	    numericOutput >= limit : key == QLatin1String("numeric_lt") ?
	    numericOutput < limit : numericOutput <= limit;
	if (!pass) {
	    if (error)
		*error = QStringLiteral(
		    "qged command numeric assertion %1 %2 failed: %3")
		    .arg(key).arg(limit, 0, 'g', 17)
		    .arg(numericOutput, 0, 'g', 17);
	    return false;
	}
    }

    const QString captureKeys[] = {
	QStringLiteral("numeric_gt_capture"),
	QStringLiteral("numeric_lt_capture"),
	QStringLiteral("numeric_near_capture")
    };
    for (const QString &key : captureKeys) {
	if (!arguments.contains(key))
	    continue;
	const QString captureName = arguments.value(key).toString();
	if (!captures || !captures->contains(captureName)) {
	    if (error)
		*error = QStringLiteral("unknown qged command capture: %1")
		    .arg(captureName);
	    return false;
	}
	bool capturedOk = false;
	const double captured = captures->value(captureName).trimmed().toDouble(
	    &capturedOk);
	if (!numericOk || !capturedOk) {
	    if (error)
		*error = QStringLiteral(
		    "non-numeric qged command capture comparison: %1")
		    .arg(captureName);
	    return false;
	}
	bool pass = key == QLatin1String("numeric_gt_capture") ?
	    numericOutput > captured : key == QLatin1String("numeric_lt_capture") ?
	    numericOutput < captured :
	    qAbs(numericOutput - captured) <=
		arguments.value(QStringLiteral("numeric_tolerance")).toDouble(1.0e-6);
	if (!pass) {
	    if (error)
		*error = QStringLiteral(
		    "qged command capture assertion %1 %2 failed: %3 versus %4")
		    .arg(key, captureName)
		    .arg(numericOutput, 0, 'g', 17)
		    .arg(captured, 0, 'g', 17);
	    return false;
	}
    }

    const QString captureName =
	arguments.value(QStringLiteral("capture")).toString();
    if (!captureName.isEmpty() && captures)
	captures->insert(captureName, output);
    return true;
}

class QgedPresentedFrameCapture {
public:
    explicit QgedPresentedFrameCapture(const QString &directory)
	: directory_(QDir(directory).absolutePath())
    {
	bool limitOk = false;
	const QByteArray configuredLimit = qgetenv("QGED_TEST_FRAME_LIMIT");
	const int configured = configuredLimit.toInt(&limitOk);
	if (limitOk && configured > 0)
	    limit_ = configured;
	if (!QDir().mkpath(directory_)) {
	    error_ = QStringLiteral("unable to create frame directory: %1")
		.arg(directory_);
	    return;
	}
	manifest_.setFileName(
	    QDir(directory_).filePath(QStringLiteral("frames.tsv")));
	stateManifest_.setFileName(
	    QDir(directory_).filePath(QStringLiteral("states.tsv")));
	if (!manifest_.open(QIODevice::WriteOnly | QIODevice::Truncate) ||
	    !stateManifest_.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
	    error_ = manifest_.errorString();
	    if (error_.isEmpty())
		error_ = stateManifest_.errorString();
	    return;
	}
	clock_.start();
	valid_ = true;
    }

    void capture(const QImage &presentedFrame,
	BObolViewController *controller)
    {
	if (!valid_ || presentedFrame.isNull() || capturing_)
	    return;
	if (captured_ >= limit_) {
	    ++dropped_;
	    return;
	}
	capturing_ = true;
	const QImage rgba = presentedFrame.convertToFormat(
	    QImage::Format_RGBA8888);
	if (rgba.isNull() || rgba.bytesPerLine() != rgba.width() * 4) {
	    ++failed_;
	    capturing_ = false;
	    return;
	}
	const QString baseName = QStringLiteral("frame-%1.rgba")
	    .arg(captured_, 6, 10, QLatin1Char('0'));
	const QString framePath = QDir(directory_).filePath(baseName);
	QFile output(framePath);
	const qint64 bytes = static_cast<qint64>(rgba.width()) *
	    static_cast<qint64>(rgba.height()) * 4;
	if (!output.open(QIODevice::WriteOnly | QIODevice::Truncate) ||
	    output.write(reinterpret_cast<const char *>(rgba.constBits()),
		bytes) != bytes) {
	    ++failed_;
	    output.remove();
	    capturing_ = false;
	    return;
	}
	output.close();
	const qint64 elapsedUsec = clock_.nsecsElapsed() / 1000;
	const QByteArray row = QStringLiteral("%1\t%2\t%3\t%4\n")
	    .arg(framePath)
	    .arg(rgba.width())
	    .arg(rgba.height())
	    .arg(elapsedUsec)
	    .toUtf8();
	if (manifest_.write(row) != row.size()) {
	    ++failed_;
	    QFile::remove(framePath);
	} else {
	    QByteArray reason = controller ?
		controller->getRenderReason().getString() : "";
	    reason.replace('\t', ' ');
	    reason.replace('\n', ' ');
	    BObolViewLodState *lodState = controller ?
		controller->getViewLodState() : nullptr;
	    size_t presentedPrimitives = 0;
	    size_t presentedRenderCost = 0;
	    const bool exactCadFrame = lodState &&
		lodState->lastCadPresentationFrameExact();
	    const bool havePresentedPrimitives = lodState &&
		lodState->lastCadPresentedPrimitiveCount(
		    presentedPrimitives);
	    const bool havePresentedRenderCost = lodState &&
		lodState->lastCadPresentedRenderCost(presentedRenderCost);
	    const QByteArray stateRow = QStringLiteral(
		"%1\t%2\t%3\t%4\t%5\t%6\t%7\t%8\t%9\t%10\t%11\t%12\t%13\t%14\n")
		.arg(captured_)
		.arg(elapsedUsec)
		.arg(controller ? static_cast<qulonglong>(
		    controller->getActiveLodMeshPayloadCount()) : 0)
		.arg(controller ? static_cast<qulonglong>(
		    controller->getActiveLodCadPayloadCount()) : 0)
		.arg(controller ? static_cast<qulonglong>(
		    controller->getLodViewRevision()) : 0)
		.arg(controller ? static_cast<qulonglong>(
		    controller->getLodPolicyRevision()) : 0)
		.arg(controller ?
		    controller->getLodInteractiveProgressiveCeiling() : -1)
		.arg(controller && controller->isLodInteractionActive() ? 1 : 0)
		.arg(exactCadFrame ? 1 : 0)
		.arg(havePresentedPrimitives ?
		    static_cast<qulonglong>(presentedPrimitives) : 0)
		.arg(havePresentedRenderCost ?
		    static_cast<qulonglong>(presentedRenderCost) : 0)
		.arg(controller ? static_cast<qulonglong>(
		    controller->getInterruptedPresentationFrameCount()) : 0)
		.arg(controller ? static_cast<qulonglong>(
		    controller->getLastInterruptedPresentationTimeNanoseconds()) : 0)
		.arg(QString::fromLocal8Bit(reason))
		.toUtf8();
	    if (stateManifest_.write(stateRow) != stateRow.size())
		++failed_;
	    ++captured_;
	}
	capturing_ = false;
    }

    void finish()
    {
	if (manifest_.isOpen())
	    manifest_.flush();
	if (stateManifest_.isOpen())
	    stateManifest_.flush();
    }

    QJsonObject summary() const
    {
	QJsonObject result;
	result.insert(QStringLiteral("directory"), directory_);
	result.insert(QStringLiteral("manifest"), manifest_.fileName());
	result.insert(QStringLiteral("state_manifest"),
	    stateManifest_.fileName());
	result.insert(QStringLiteral("captured"), captured_);
	result.insert(QStringLiteral("dropped"), dropped_);
	result.insert(QStringLiteral("failed"), failed_);
	result.insert(QStringLiteral("limit"), limit_);
	result.insert(QStringLiteral("valid"), valid_);
	if (!error_.isEmpty())
	    result.insert(QStringLiteral("error"), error_);
	return result;
    }

private:
    QString directory_;
    QFile manifest_;
    QFile stateManifest_;
    QElapsedTimer clock_;
    int captured_ = 0;
    int dropped_ = 0;
    int failed_ = 0;
    /*
     * libicv deliberately owns lossless image copies while assembling an
     * APNG.  Keep an accidental long-running diagnostic from consuming
     * unbounded memory; callers may raise this explicit test-only limit.
     */
    int limit_ = 300;
    bool valid_ = false;
    bool capturing_ = false;
    QString error_;
};

static bool
qged_test_wait_progressive_idle(QgEdApp &app, int timeoutMilliseconds,
    int quietMilliseconds, QString *error)
{
    const std::vector<BObolViewController *> controllers =
	qged_test_all_controllers(app);
    if (controllers.empty()) {
	if (error)
	    *error = QStringLiteral(
		"wait_progressive_idle requires an Obol view controller");
	return false;
    }

    timeoutMilliseconds = std::max(0, timeoutMilliseconds);
    quietMilliseconds = std::max(0, quietMilliseconds);
    QElapsedTimer elapsed;
    QElapsedTimer quiet;
    elapsed.start();

    while (elapsed.elapsed() <= timeoutMilliseconds) {
	bool idle = true;
	for (BObolViewController *controller : controllers) {
	    BObolLodService *service = controller->getLodService();
	    const bool serviceIdle =
		!service ||
		(service->pendingTaskCountForDiagnostics() == 0 &&
		 service->delayedTaskCountForDiagnostics() == 0 &&
		 service->inFlightCount() == 0 &&
		 service->activeRequestCountForDiagnostics() == 0 &&
		 service->queuedResultCountForDiagnostics() == 0 &&
		 service->queuedCacheWriteCountForDiagnostics() == 0);
	    if (controller->hasProgressiveWorkPending() ||
		controller->hasPendingLodResults() ||
		controller->hasPendingLodSubmissions() ||
		controller->hasPendingLodRefinementFrame() ||
		controller->isRenderRequested() ||
		!serviceIdle) {
		idle = false;
		break;
	    }
	}

	if (idle) {
	    if (!quiet.isValid())
		quiet.start();
	    if (quiet.elapsed() >= quietMilliseconds)
		return true;
	} else {
	    quiet.invalidate();
	}

	/* Let queued paints, worker wakeups, and the progressive frame pump
	 * advance while retaining a hard deadline for a defective pipeline. */
	QEventLoop loop;
	const int remaining = timeoutMilliseconds -
	    static_cast<int>(elapsed.elapsed());
	QTimer::singleShot(std::max(1, std::min(16, remaining)),
	    &loop, &QEventLoop::quit);
	loop.exec(QEventLoop::AllEvents);
    }

    if (error) {
	BObolViewController *controller = controllers.front();
	for (BObolViewController *candidate : controllers) {
	    BObolLodService *candidateService = candidate->getLodService();
	    if (candidate->hasProgressiveWorkPending() ||
		candidate->hasPendingLodResults() ||
		candidate->hasPendingLodSubmissions() ||
		candidate->hasPendingLodRefinementFrame() ||
		candidate->isRenderRequested() ||
		(candidateService &&
		 (candidateService->pendingTaskCountForDiagnostics() != 0 ||
		  candidateService->delayedTaskCountForDiagnostics() != 0 ||
		  candidateService->inFlightCount() != 0 ||
		  candidateService->activeRequestCountForDiagnostics() != 0 ||
		  candidateService->queuedResultCountForDiagnostics() != 0 ||
		  candidateService->queuedCacheWriteCountForDiagnostics() != 0))) {
		controller = candidate;
		break;
	    }
	}
	BObolLodService *service = controller->getLodService();
	*error = QStringLiteral(
	    "progressive pipeline did not become idle within %1 ms "
	    "across %2 controller(s) "
	    "(progressive=%3 results=%4 submissions=%5 refinement_frame=%6 "
	    "render=%7 pending=%8 delayed=%9 in_flight=%10 active=%11 "
	    "queued=%12 cache_writes=%13 request_serial=%14 "
	    "completion_serial=%15 presented_serial=%16 settle_after=%17 "
	    "refinement_after=%18)")
	    .arg(timeoutMilliseconds)
	    .arg(controllers.size())
	    .arg(controller->hasProgressiveWorkPending() ? 1 : 0)
	    .arg(controller->hasPendingLodResults() ? 1 : 0)
	    .arg(controller->hasPendingLodSubmissions() ? 1 : 0)
	    .arg(controller->hasPendingLodRefinementFrame() ? 1 : 0)
	    .arg(controller->isRenderRequested() ? 1 : 0)
	    .arg(service ?
		static_cast<qulonglong>(
		    service->pendingTaskCountForDiagnostics()) : 0)
	    .arg(service ?
		static_cast<qulonglong>(
		    service->delayedTaskCountForDiagnostics()) : 0)
	    .arg(service ?
		static_cast<qulonglong>(service->inFlightCount()) : 0)
	    .arg(service ?
		static_cast<qulonglong>(
		    service->activeRequestCountForDiagnostics()) : 0)
	    .arg(service ?
		static_cast<qulonglong>(
		    service->queuedResultCountForDiagnostics()) : 0)
	    .arg(service ?
		static_cast<qulonglong>(
		    service->queuedCacheWriteCountForDiagnostics()) : 0)
	    .arg(static_cast<qulonglong>(
		controller->getRenderRequestSerial()))
	    .arg(static_cast<qulonglong>(
		controller->getRenderCompletionSerial()))
	    .arg(static_cast<qulonglong>(
		controller->getPresentedFrameSerial()))
	    .arg(static_cast<qulonglong>(
		controller->getLodSettleAfterRenderSerial()))
	    .arg(static_cast<qulonglong>(
		controller->getLodRefinementResumeAfterRenderSerial()));
    }
    return false;
}

static bool
qged_test_wait_progressive_scope_ready(QgEdApp &app,
    int timeoutMilliseconds, int quietMilliseconds, QString *error)
{
    QgView *currentView = app.w ? app.w->CurrentDisplay() : nullptr;
    BObolViewController *currentController = currentView ?
	currentView->obolViewController() : nullptr;
    if (!currentController) {
	if (error)
	    *error = QStringLiteral(
		"wait_progressive_scope_ready requires an Obol view controller");
	return false;
    }
    /* Event targets named "." address the active canvas.  A QgQuadView also
     * retains inactive controllers whose source-visible flags remain true;
     * those widgets are intentionally not painted and cannot satisfy a
     * presentation barrier.  Whole-target readiness here is therefore the
     * target view's scope, not every dormant view owned by the window. */
    const std::vector<BObolViewController *> controllers(1,
	currentController);

    timeoutMilliseconds = std::max(0, timeoutMilliseconds);
    quietMilliseconds = std::max(0, quietMilliseconds);
    QElapsedTimer elapsed;
    QElapsedTimer quiet;
    std::vector<std::pair<BObolViewController *, uint64_t>>
	scopePresentationBarriers;
    size_t lastVisibleSourceCount = 0;
    size_t lastInexactSourceCount = 0;
    QString lastInexactSource;
    int lastInexactSourceStatus = -1;
    QString lastInexactSourceDiagnostic;
    elapsed.start();

    while (elapsed.elapsed() <= timeoutMilliseconds) {
	bool foundVisibleSource = false;
	bool ready = true;
	size_t visibleSourceCount = 0;
	size_t inexactSourceCount = 0;
	QString inexactSource;
	int inexactSourceStatus = -1;
	QString inexactSourceDiagnostic;
	std::vector<BObolViewController *> visibleControllers;
	for (BObolViewController *controller : controllers) {
	    const bool endpointVisible =
		qged_test_controller_has_visible_canvas(app, controller);
	    bool controllerVisible = false;
	    const std::vector<SoBRLDatabaseSource *> sources =
		controller->getRenderDatabaseSources();
	    for (SoBRLDatabaseSource *source : sources) {
		if (!endpointVisible || !source ||
		    !source->visible.getValue())
		    continue;
		foundVisibleSource = true;
		controllerVisible = true;
		visibleSourceCount++;
		if (source->realizationStatus.getValue() ==
			SoBRLDatabaseSource::FAILED) {
		    if (error)
			*error = QStringLiteral(
			    "progressive whole-target source failed "
			    "(source=%1 diagnostic=%2)")
			    .arg(QString::fromLocal8Bit(
				source->path.getValue().getString()))
			    .arg(QString::fromLocal8Bit(
				source->realizationDiagnostic.getValue().
				getString()));
		    return false;
		}
		SbBox3f bounds;
		if (!source->hasExactSourceBounds() ||
		    !source->getSourceBounds(bounds) || bounds.isEmpty()) {
		    inexactSourceCount++;
		    if (inexactSource.isEmpty())
			inexactSource = QString::fromLocal8Bit(
			    source->path.getValue().getString());
		    if (inexactSourceStatus < 0) {
			inexactSourceStatus =
			    source->realizationStatus.getValue();
			inexactSourceDiagnostic = QString::fromLocal8Bit(
			    source->realizationDiagnostic.getValue().getString());
		    }
		    ready = false;
		    break;
		}
	    }
	    if (controllerVisible)
		visibleControllers.push_back(controller);
	    if (!ready)
		break;
	}
	ready = ready && foundVisibleSource;
	lastVisibleSourceCount = visibleSourceCount;
	lastInexactSourceCount = inexactSourceCount;
	lastInexactSource = inexactSource;
	lastInexactSourceStatus = inexactSourceStatus;
	lastInexactSourceDiagnostic = inexactSourceDiagnostic;

	if (ready) {
	    if (scopePresentationBarriers.empty()) {
		/* Exact source bounds and the one-shot deferred autoview are
		 * published on the owner thread before the resulting frame is
		 * necessarily presented.  Capturing immediately produced a stale
		 * pre-autoview framebuffer on sufficiently large cold scenes.  Ask
		 * every visible endpoint for one deterministic scope frame and wait
		 * for its presentation serial before starting the quiet interval. */
		for (BObolViewController *controller : visibleControllers) {
		    scopePresentationBarriers.emplace_back(controller,
			controller->getPresentedFrameSerial());
		    controller->requestRender("qged-test-scope-ready");
		    qged_test_present_controller_frame(app, controller);
		}
		quiet.invalidate();
	    }
	    bool scopePresented = !scopePresentationBarriers.empty();
	    for (const auto &barrier : scopePresentationBarriers) {
		if (!barrier.first ||
		    barrier.first->getPresentedFrameSerial() <= barrier.second) {
		    scopePresented = false;
		    break;
		}
	    }
	    if (scopePresented && !quiet.isValid())
		quiet.start();
	    if (scopePresented && quiet.elapsed() >= quietMilliseconds)
		return true;
	} else {
	    scopePresentationBarriers.clear();
	    quiet.invalidate();
	}

	QEventLoop loop;
	const int remaining = timeoutMilliseconds -
	    static_cast<int>(elapsed.elapsed());
	QTimer::singleShot(std::max(1, std::min(16, remaining)),
	    &loop, &QEventLoop::quit);
	loop.exec(QEventLoop::AllEvents);
    }

    if (error)
	*error = QStringLiteral(
	    "progressive whole-target scope did not become ready within %1 ms "
	    "(visible_sources=%2 inexact_sources=%3 inexact_source=%4 "
	    "inexact_status=%5 inexact_diagnostic=%6 request_serial=%7 "
	    "completion_serial=%8 presented_serial=%9 settle_after=%10 "
	    "refinement_after=%11)")
	    .arg(timeoutMilliseconds)
	    .arg(static_cast<qulonglong>(lastVisibleSourceCount))
	    .arg(static_cast<qulonglong>(lastInexactSourceCount))
	    .arg(lastInexactSource)
	    .arg(lastInexactSourceStatus)
	    .arg(lastInexactSourceDiagnostic)
	    .arg(static_cast<qulonglong>(
		controllers.front()->getRenderRequestSerial()))
	    .arg(static_cast<qulonglong>(
		controllers.front()->getRenderCompletionSerial()))
	    .arg(static_cast<qulonglong>(
		controllers.front()->getPresentedFrameSerial()))
	    .arg(static_cast<qulonglong>(
		controllers.front()->getLodSettleAfterRenderSerial()))
	    .arg(static_cast<qulonglong>(
		controllers.front()->getLodRefinementResumeAfterRenderSerial()));
    return false;
}

static bool
qged_test_wait_progressive_discovery_complete(QgEdApp &app,
    int timeoutMilliseconds, int quietMilliseconds, QString *error)
{
    QgView *currentView = app.w ? app.w->CurrentDisplay() : nullptr;
    BObolViewController *controller = currentView ?
	currentView->obolViewController() : nullptr;
    if (!controller) {
	if (error)
	    *error = QStringLiteral(
		"wait_progressive_discovery_complete requires an Obol view controller");
	return false;
    }

    timeoutMilliseconds = std::max(0, timeoutMilliseconds);
    quietMilliseconds = std::max(0, quietMilliseconds);
    QElapsedTimer elapsed;
    QElapsedTimer quiet;
    size_t lastAvailable = 0;
    size_t lastExpected = 0;
    int lastPhase = -1;
    elapsed.start();
    while (elapsed.elapsed() <= timeoutMilliseconds) {
	BObolLodConvergenceStatus status;
	controller->getLodConvergenceStatus(status);
	lastAvailable = status.availableLeafCount;
	lastExpected = status.expectedLeafCount;
	lastPhase = status.phase;
	const bool complete = lastExpected > 0 &&
	    lastAvailable >= lastExpected;
	if (complete) {
	    if (!quiet.isValid())
		quiet.start();
	    if (quiet.elapsed() >= quietMilliseconds)
		return true;
	} else {
	    quiet.invalidate();
	}

	QEventLoop loop;
	const int remaining = timeoutMilliseconds -
	    static_cast<int>(elapsed.elapsed());
	QTimer::singleShot(std::max(1, std::min(16, remaining)),
	    &loop, &QEventLoop::quit);
	loop.exec(QEventLoop::AllEvents);
    }

    if (error)
	*error = QStringLiteral(
	    "progressive discovery did not publish all leaves within %1 ms "
	    "(available=%2 expected=%3 phase=%4 request_serial=%5 "
	    "completion_serial=%6 presented_serial=%7)")
	    .arg(timeoutMilliseconds)
	    .arg(static_cast<qulonglong>(lastAvailable))
	    .arg(static_cast<qulonglong>(lastExpected))
	    .arg(lastPhase)
	    .arg(static_cast<qulonglong>(
		controller->getRenderRequestSerial()))
	    .arg(static_cast<qulonglong>(
		controller->getRenderCompletionSerial()))
	    .arg(static_cast<qulonglong>(
		controller->getPresentedFrameSerial()));
    return false;
}

void
qged_schedule_gui_test(QgEdApp &app, const QString &script,
    const QString &reportFile, bool softwareRenderer)
{
    QTimer::singleShot(0, &app,
	[&app, script, reportFile, softwareRenderer]() {
	    QVector<QgTestEvent> events;
	    QString error;
	    QgEventPlayer player(app.w);
	    QElapsedTimer elapsed;
	    QJsonArray samples;
	    QHash<QString, QString> commandCaptures;
	    QJsonObject report;
	    std::shared_ptr<QgedPresentedFrameCapture> frameCapture;
	    report.insert(QStringLiteral("schema"),
		QStringLiteral("brlcad.qged.gui-report"));
	    /* Version 2 names producer-authored progressive cuts explicitly and
	     * removes the former *_level telemetry fields. */
	    report.insert(QStringLiteral("version"), 2);
	    report.insert(QStringLiteral("event_script"), script);
	    report.insert(QStringLiteral("backend"),
		softwareRenderer ? QStringLiteral("osmesa") :
		QStringLiteral("system_gl"));
	    report.insert(QStringLiteral("cache_directory"),
		QString::fromLocal8Bit(std::getenv("BU_DIR_CACHE") ?
		    std::getenv("BU_DIR_CACHE") : ""));
	    report.insert(QStringLiteral("swap_interval"),
		QString::fromLocal8Bit(std::getenv("QGED_SWAP_INTERVAL") ?
		    std::getenv("QGED_SWAP_INTERVAL") : "default"));
	    const QString frameDirectory = QString::fromLocal8Bit(
		std::getenv("QGED_TEST_FRAME_DIR") ?
		std::getenv("QGED_TEST_FRAME_DIR") : "");
	    if (!frameDirectory.isEmpty() && app.w) {
		frameCapture = std::make_shared<QgedPresentedFrameCapture>(
		    frameDirectory);
		QgView *view = app.w->CurrentDisplay();
		QgCanvasBase *canvas = view ? view->canvasBase() : nullptr;
		QWidget *canvasWidget = canvas ? canvas->canvasWidget() : nullptr;
		BObolViewController *captureController =
		    view ? view->obolViewController() : nullptr;
		if (QgGL *glWidget = qobject_cast<QgGL *>(canvasWidget)) {
		    QObject::connect(glWidget, &QgGL::frame_presented, &app,
			[frameCapture, captureController](const QImage &image) {
			    frameCapture->capture(image, captureController);
			});
		} else if (QgSW *softwareWidget =
			qobject_cast<QgSW *>(canvasWidget)) {
		    QObject::connect(softwareWidget, &QgSW::frame_presented,
			&app, [frameCapture, captureController](const QImage &image) {
			    frameCapture->capture(image, captureController);
			});
		}
	    }
	    player.setCheckpointHandler([](QWidget *widget,
		    const QString &name, QString *checkpointError) {
		if (!widget || name.isEmpty()) {
		    if (checkpointError)
			*checkpointError =
			    QStringLiteral(
				"checkpoint needs a widget and output path");
		    return false;
		}
		bool saved = false;
	if (QOpenGLWidget *glWidget =
			qobject_cast<QOpenGLWidget *>(widget)) {
		    const QImage frame = glWidget->grabFramebuffer();
		    saved = !frame.isNull() && frame.save(name);
		} else if (QgSW *softwareWidget =
			qobject_cast<QgSW *>(widget)) {
		    /* QWidget::grab() invokes QgSW::paintEvent(), turning a
		     * diagnostic readback into a second deadline-governed
		     * presentation.  On software GL the readback itself can cross
		     * that deadline and make checkpoint frequency alter the retained
		     * PoP cut.  The canvas image API performs the same traversal as an
		     * observational export: no progressive advance, capacity sample,
		     * or deadline feedback. */
		    QImage frame;
		    softwareWidget->get_viewport_image(frame);
		    saved = !frame.isNull() && frame.save(name);
		} else {
		    const QPixmap frame = widget->grab();
		    saved = !frame.isNull() && frame.save(name);
		}
		if (!saved) {
		    if (checkpointError)
			*checkpointError = QStringLiteral(
			    "unable to save checkpoint image: %1").arg(name);
		    return false;
		}
		return true;
	    });
	    elapsed.start();
	    bool success = QgEventRecorder::load(script, events, &error);
	    for (int i = 0; success && i < events.size(); ++i) {
		QElapsedTimer eventElapsed;
		eventElapsed.start();
		QString commandOutput;
		if (events[i].action ==
		    QLatin1String("qged_command_expect")) {
		    success = qged_test_command_expect(
			app, events[i].arguments, &commandCaptures,
			&commandOutput, &error);
		} else if (events[i].action == QLatin1String("qged_command")) {
		    const QString command = events[i].arguments.value(
			QStringLiteral("command")).toString();
		    if (command.isEmpty()) {
			error = QStringLiteral(
			    "qged_command requires a command argument");
			success = false;
		    } else {
			app.run_qcmd(command);
		    }
		} else if (events[i].action ==
		    QLatin1String("qged_command_batch")) {
		    const QJsonArray commands = events[i].arguments.value(
			QStringLiteral("commands")).toArray();
		    if (commands.isEmpty()) {
			error = QStringLiteral(
			    "qged_command_batch requires a commands array");
			success = false;
		    } else {
			/* Framing commonly requires an orientation followed by an
			 * autoview.  Execute that transaction without exposing the
			 * deliberately invalid intermediate camera as a one-frame GUI
			 * flash; controller state and final rendering remain ordinary. */
			const bool updatesEnabled = app.w &&
			    app.w->updatesEnabled();
			if (app.w && updatesEnabled)
			    app.w->setUpdatesEnabled(false);
			for (const QJsonValue &value : commands) {
			    const QString command = value.toString();
			    if (command.isEmpty()) {
				error = QStringLiteral(
				    "qged_command_batch contains an empty command");
				success = false;
				break;
			    }
			    app.run_qcmd(command);
			}
			if (app.w && updatesEnabled) {
			    app.w->setUpdatesEnabled(true);
			    app.w->update();
			}
		    }
		} else if (events[i].action ==
		    QLatin1String("wait_progressive_idle")) {
		    success = qged_test_wait_progressive_idle(app,
			events[i].arguments.value(
			    QStringLiteral("timeout_ms")).toInt(30000),
			events[i].arguments.value(
			    QStringLiteral("quiet_ms")).toInt(100),
			&error);
		} else if (events[i].action ==
		    QLatin1String("wait_progressive_scope_ready")) {
		    success = qged_test_wait_progressive_scope_ready(app,
			events[i].arguments.value(
			    QStringLiteral("timeout_ms")).toInt(30000),
			events[i].arguments.value(
			    QStringLiteral("quiet_ms")).toInt(50),
			&error);
		} else if (events[i].action ==
		    QLatin1String("wait_progressive_discovery_complete")) {
		    success = qged_test_wait_progressive_discovery_complete(app,
			events[i].arguments.value(
			    QStringLiteral("timeout_ms")).toInt(30000),
			events[i].arguments.value(
			    QStringLiteral("quiet_ms")).toInt(50),
			&error);
		} else {
		    success = player.play(events[i], &error);
		}

		/* Commands and synthetic input dispatch synchronously.  Explicit
		 * wait events, rather than a processEvents drain after every
		 * input, provide real event-loop intervals without consuming a
		 * progressive-frame backlog. */
		const qint64 eventMicroseconds =
		    eventElapsed.nsecsElapsed() / 1000;
		const char *deepReportMode =
		    std::getenv("QGED_TEST_DEEP_LOD_REPORT");
		const bool deepReportExplicitlyDisabled =
		    deepReportMode &&
		    (BU_STR_EQUAL(deepReportMode, "0") ||
		     BU_STR_EQUAL(deepReportMode, "false") ||
		     BU_STR_EQUAL(deepReportMode, "off"));
		const bool deepReportEverySample =
		    deepReportMode && !deepReportExplicitlyDisabled;
		const bool deepLodDiagnostics =
		    deepReportEverySample ||
		    (!deepReportExplicitlyDisabled &&
		     i + 1 == events.size());
		QElapsedTimer sampleElapsed;
		sampleElapsed.start();
		QJsonObject sample = qged_collect_progressive_sample(
		    app, i, events[i], elapsed.elapsed(), eventMicroseconds,
		    deepLodDiagnostics);
		if (!commandOutput.isNull())
		    sample.insert(QStringLiteral("command_output"),
			commandOutput);
		sample.insert(QStringLiteral("sample_duration_us"),
		    sampleElapsed.nsecsElapsed() / 1000);
		samples.append(sample);
	    }
	    report.insert(QStringLiteral("success"), success);
	    report.insert(QStringLiteral("elapsed_ms"), elapsed.elapsed());
	    report.insert(QStringLiteral("samples"), samples);
	    if (frameCapture) {
		frameCapture->finish();
		report.insert(QStringLiteral("presented_frame_capture"),
		    frameCapture->summary());
	    }
	    if (!success)
		report.insert(QStringLiteral("error"), error);
	    QString reportError;
	    if (!qged_write_test_report(reportFile, report, &reportError)) {
		bu_log("qged: unable to write GUI test report: %s\n",
		    reportError.toLocal8Bit().constData());
		success = false;
	    }
	    if (!success) {
		bu_log("qged: GUI test replay failed: %s\n",
		    error.toLocal8Bit().constData());
		app.exit(BRLCAD_ERROR);
		return;
	    }
	    app.exit(BRLCAD_OK);
	});
}
