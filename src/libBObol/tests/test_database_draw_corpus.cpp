/*          T E S T _ D A T A B A S E _ D R A W _ C O R P U S . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol.h"
#include "BObol/BDatabaseSource.h"

#include "bu/app.h"
#include "raytrace.h"
#include "rt/display_bounds.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <string>
#include <system_error>
#include <vector>

namespace {

constexpr uint32_t kProgressiveBotThreshold = 1;
constexpr float kViewScale = 100.0f;
constexpr float kLodScale = 1.0f;
constexpr int kViewWidth = 1024;
constexpr int kViewHeight = 1024;
constexpr size_t kDrainBatchSize = 4096;

enum class DrawMode {
    Wire,
    Shaded
};

struct DatabaseObject {
    std::string name;
};

struct AuditTotals {
    size_t databaseCount = 0;
    size_t objectCount = 0;
    size_t modeCount = 0;
    size_t fallbackCount = 0;
    size_t emptyCount = 0;
    size_t failureCount = 0;
};

const char *
draw_mode_name(DrawMode mode)
{
    return mode == DrawMode::Wire ? "wire" : "shaded";
}

void
report_failure(const std::filesystem::path &databasePath,
	const DatabaseObject &object, DrawMode mode, const char *failure,
	const char *diagnostic, AuditTotals &totals)
{
    std::fprintf(stderr, "FAIL\t%s\t%s\t%s\t%s",
	databasePath.string().c_str(), object.name.c_str(),
	draw_mode_name(mode), failure ? failure : "unspecified failure");
    if (diagnostic && diagnostic[0])
	std::fprintf(stderr, "\t%s", diagnostic);
    std::fprintf(stderr, "\n");
    totals.failureCount++;
}

bool
database_object_has_display_bounds(struct db_i *dbip,
	const DatabaseObject &object)
{
    point_t minimum;
    point_t maximum;
    struct bu_vls diagnostic = BU_VLS_INIT_ZERO;
    const char *objects[1] = {object.name.c_str()};
    const int result = rt_display_bounds(&diagnostic, dbip, 1, objects,
	minimum, maximum);
    bu_vls_free(&diagnostic);
    if (result != BRLCAD_OK)
	return false;
    for (int axis = 0; axis < 3; axis++) {
	if (!std::isfinite(minimum[axis]) || !std::isfinite(maximum[axis]) ||
	    minimum[axis] > maximum[axis])
	    return false;
    }
    return true;
}

bool
source_mesh_identity_is_valid(struct db_i *dbip,
	const BObolCompactManifestOccurrence &occurrence)
{
    if (!occurrence.sourceMeshRequestValid)
	return true;

    const char *assetPath = occurrence.meshAssetPath.getString();
    const char *assetName = occurrence.meshAssetName.getString();
    if (!assetPath || !assetPath[0] || !assetName || !assetName[0])
	return false;
    if (db_lookup(dbip, assetName, LOOKUP_QUIET) == RT_DIR_NULL)
	return false;

    const size_t pathLength = std::char_traits<char>::length(assetPath);
    const size_t nameLength = std::char_traits<char>::length(assetName);
    if (pathLength < nameLength)
	return false;
    return std::string(assetPath + pathLength - nameLength) == assetName;
}

bool
source_mesh_identity_is_valid(struct db_i *dbip,
	const BObolCompactOccurrence &occurrence)
{
    if (!occurrence.sourceMeshRequestValid)
	return true;

    BObolCompactManifestOccurrence manifest;
    manifest.sourceMeshRequestValid = TRUE;
    manifest.meshAssetPath = occurrence.sourceMeshRequest.meshAssetPath;
    manifest.meshAssetName = occurrence.sourceMeshRequest.meshAssetName;
    return source_mesh_identity_is_valid(dbip, manifest);
}

bool
audit_object_mode(struct db_i *dbip,
	const std::filesystem::path &databasePath, const DatabaseObject &object,
	DrawMode mode, uint32_t revision, AuditTotals &totals)
{
    totals.modeCount++;

    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->ref();
    source->setDatabase(dbip);
    source->path = object.name.c_str();
    source->drawMode = mode == DrawMode::Wire ?
	SoBRLDatabaseSource::WIREFRAME : SoBRLDatabaseSource::SHADED;
    source->lodBotThreshold = kProgressiveBotThreshold;
    source->sourceRevision = revision;
    source->setRealizationViewPolicy(TRUE, TRUE, TRUE, kViewScale,
	kLodScale, kViewWidth, kViewHeight, kProgressiveBotThreshold,
	0.0f, 0.0f);

    BObolCompactOccurrenceStream stream;
    const bool meshRealization = source->usesMeshRealization() ? true : false;
    SbBool realized = meshRealization ?
	source->realizeDatabaseMesh(&stream) :
	source->realizeDatabaseWireframe(&stream);
    if (!realized && meshRealization) {
	/* This is the same deliberate compatibility fallback used by the GED
	 * drawing owner for primitives without a shaded representation. */
	realized = source->realizeDatabaseWireframe(&stream);
	if (realized)
	    totals.fallbackCount++;
    }

    const char *diagnostic = source->realizationDiagnostic.getValue().getString();
    if (!realized || source->realizationStatus.getValue() !=
	    SoBRLDatabaseSource::REALIZED) {
	if (!database_object_has_display_bounds(dbip, object)) {
	    totals.emptyCount++;
	    source->unref();
	    return true;
	}
	report_failure(databasePath, object, mode,
	    "realization did not reach REALIZED", diagnostic, totals);
	source->unref();
	return false;
    }

    const size_t expectedCount = stream.getExpectedCount();
    if (expectedCount == 0) {
	report_failure(databasePath, object, mode,
	    "successful realization published no expected occurrence count",
	    diagnostic, totals);
	source->unref();
	return false;
    }

    size_t preparationCompleted = 0;
    size_t preparationTotal = 0;
    stream.getPreparationWorkCount(preparationCompleted, preparationTotal);
    if (preparationCompleted != preparationTotal) {
	report_failure(databasePath, object, mode,
	    "preparation work did not close", diagnostic, totals);
	source->unref();
	return false;
    }

    BObolCompactSourceProfile profile;
    if (!stream.getSourceProfile(profile) || !profile.valid ||
	profile.occurrenceCount != expectedCount ||
	profile.uniqueAssetCount == 0 ||
	profile.uniqueAssetCount > profile.occurrenceCount ||
	profile.encodedSourceBytes == 0 || profile.largestAssetBytes == 0) {
	report_failure(databasePath, object, mode,
	    "source profile is absent or inconsistent", diagnostic, totals);
	source->unref();
	return false;
    }

    std::vector<BObolCompactManifestOccurrence> manifest;
    if (!stream.takeManifest(manifest) || manifest.size() != expectedCount) {
	report_failure(databasePath, object, mode,
	    "authoritative manifest did not seal at the expected count",
	    diagnostic, totals);
	source->unref();
	return false;
    }
    for (const BObolCompactManifestOccurrence &occurrence : manifest) {
	if (occurrence.path.getLength() == 0 ||
	    !source_mesh_identity_is_valid(dbip, occurrence)) {
	    report_failure(databasePath, object, mode,
		"manifest contains an invalid path or mesh asset identity",
		diagnostic, totals);
	    source->unref();
	    return false;
	}
    }

    size_t occurrenceCount = 0;
    std::vector<BObolCompactOccurrence> occurrences;
    for (;;) {
	occurrences.clear();
	const size_t drained = stream.drain(occurrences, kDrainBatchSize);
	if (drained == 0)
	    break;
	occurrenceCount += drained;
	for (const BObolCompactOccurrence &occurrence : occurrences) {
	    if (!occurrence.geometry ||
		!source_mesh_identity_is_valid(dbip, occurrence)) {
		report_failure(databasePath, object, mode,
		    "stream contains null geometry or an invalid mesh asset identity",
		    diagnostic, totals);
		source->unref();
		return false;
	    }
	}
    }
    if (occurrenceCount < expectedCount) {
	report_failure(databasePath, object, mode,
	    "stream contains fewer occurrences than its expected population",
	    diagnostic, totals);
	source->unref();
	return false;
    }

    SbBox3f sourceBounds;
    if (!source->hasExactSourceBounds() ||
	!source->getSourceBounds(sourceBounds) || sourceBounds.isEmpty()) {
	report_failure(databasePath, object, mode,
	    "successful realization has no exact nonempty source bound",
	    diagnostic, totals);
	source->unref();
	return false;
    }

    source->unref();
    return true;
}

bool
audit_database(const std::filesystem::path &databasePath,
	const std::vector<DrawMode> &modes, const std::string &objectFilter,
	bool reportProgress, AuditTotals &totals)
{
    struct db_i *dbip = db_open(databasePath.string().c_str(),
	DB_OPEN_READONLY);
    if (!dbip) {
	std::fprintf(stderr, "FAIL\t%s\t-\t-\tcould not open database\n",
	    databasePath.string().c_str());
	totals.failureCount++;
	return false;
    }
    if (db_dirbuild(dbip) < 0) {
	std::fprintf(stderr,
	    "FAIL\t%s\t-\t-\tcould not build database directory\n",
	    databasePath.string().c_str());
	totals.failureCount++;
	db_close(dbip);
	return false;
    }

    std::vector<DatabaseObject> objects;
    struct directory *dp = RT_DIR_NULL;
    FOR_ALL_DIRECTORY_START(dp, dbip) {
	if (!dp->d_namep || !(dp->d_flags & (RT_DIR_SOLID | RT_DIR_COMB)))
	    continue;
	DatabaseObject object;
	object.name = dp->d_namep;
	if (!objectFilter.empty() && object.name != objectFilter)
	    continue;
	objects.push_back(std::move(object));
    }
    FOR_ALL_DIRECTORY_END;
    std::sort(objects.begin(), objects.end(),
	[](const DatabaseObject &left, const DatabaseObject &right) {
	    return left.name < right.name;
	});

    if (!objectFilter.empty() && objects.empty()) {
	std::fprintf(stderr, "FAIL\t%s\t%s\t-\tobject was not found\n",
	    databasePath.string().c_str(), objectFilter.c_str());
	totals.databaseCount++;
	totals.failureCount++;
	db_close(dbip);
	return false;
    }

    totals.databaseCount++;
    totals.objectCount += objects.size();
    size_t databaseFailuresBefore = totals.failureCount;
    uint32_t revision = 1;
    for (const DatabaseObject &object : objects) {
	for (DrawMode mode : modes) {
	    const size_t failuresBefore = totals.failureCount;
	    const size_t emptyBefore = totals.emptyCount;
	    const auto start = std::chrono::steady_clock::now();
	    const bool passed = audit_object_mode(dbip, databasePath, object,
		mode, revision++, totals);
	    if (reportProgress) {
		const auto elapsed = std::chrono::duration_cast<
		    std::chrono::milliseconds>(
		    std::chrono::steady_clock::now() - start);
		const char *result = totals.emptyCount > emptyBefore ? "empty" :
		    (passed && totals.failureCount == failuresBefore ?
		    "pass" : "fail");
		std::printf("OBJECT\t%s\t%s\t%s\t%s\t%lldms\n",
		    databasePath.string().c_str(), object.name.c_str(),
		    draw_mode_name(mode), result,
		    static_cast<long long>(elapsed.count()));
		std::fflush(stdout);
	    }
	}
    }

    std::printf("DATABASE\t%s\tobjects=%zu\tmodes=%zu\tfailures=%zu\n",
	databasePath.string().c_str(), objects.size(), modes.size(),
	totals.failureCount - databaseFailuresBefore);
    std::fflush(stdout);
    db_close(dbip);
    return totals.failureCount == databaseFailuresBefore;
}

void
collect_database_paths(const std::filesystem::path &input,
	std::vector<std::filesystem::path> &databases)
{
    std::error_code error;
    if (std::filesystem::is_regular_file(input, error)) {
	if (input.extension() == ".g")
	    databases.push_back(input);
	return;
    }
    if (!std::filesystem::is_directory(input, error))
	return;

    const auto options =
	std::filesystem::directory_options::skip_permission_denied;
    for (std::filesystem::recursive_directory_iterator iterator(input,
	     options, error), end;
	 iterator != end; iterator.increment(error)) {
	if (error) {
	    error.clear();
	    continue;
	}
	if (iterator->is_regular_file(error) &&
	    iterator->path().extension() == ".g")
	    databases.push_back(iterator->path());
    }
}

void
usage(const char *program)
{
    std::fprintf(stderr,
	"Usage: %s [--wire|--shaded] [--object name] "
	"[--progress] database-or-directory [...]\n",
	program);
}

} // namespace

int
main(int argc, const char **argv)
{
    bu_setprogname(argv[0]);
    bobol_init(NULL);

    std::vector<DrawMode> modes = {DrawMode::Wire, DrawMode::Shaded};
    std::vector<std::filesystem::path> databases;
    std::string objectFilter;
    bool reportProgress = false;
    for (int argumentIndex = 1; argumentIndex < argc; argumentIndex++) {
	const std::string argument = argv[argumentIndex];
	if (argument == "--wire") {
	    modes = {DrawMode::Wire};
	    continue;
	}
	if (argument == "--shaded") {
	    modes = {DrawMode::Shaded};
	    continue;
	}
	if (argument == "--object" && argumentIndex + 1 < argc) {
	    objectFilter = argv[++argumentIndex];
	    continue;
	}
	if (argument == "--progress") {
	    reportProgress = true;
	    continue;
	}
	if (!argument.empty() && argument[0] == '-') {
	    usage(argv[0]);
	    return 2;
	}
	collect_database_paths(argument, databases);
    }
    if (databases.empty()) {
	usage(argv[0]);
	return 2;
    }

    std::sort(databases.begin(), databases.end());
    databases.erase(std::unique(databases.begin(), databases.end()),
	databases.end());

    AuditTotals totals;
    for (const std::filesystem::path &database : databases)
	(void)audit_database(database, modes, objectFilter, reportProgress,
	    totals);

    std::printf("TOTAL\tdatabases=%zu\tobjects=%zu\tmodes=%zu\t"
	"fallbacks=%zu\tempty=%zu\tfailures=%zu\n", totals.databaseCount,
	totals.objectCount, totals.modeCount, totals.fallbackCount,
	totals.emptyCount,
	totals.failureCount);
    return totals.failureCount == 0 ? 0 : 1;
}
