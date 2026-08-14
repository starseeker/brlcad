/*            M E S H _ L O D _ C A C H E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * Based off of code from https://github.com/bhaettasch/pop-buffer-demo
 * Copyright (c) 2016 Benjamin Hattasch and X3DOM
 * The MIT License (MIT)
 */

#include "common.h"

#include "bv.h"
#include "BObol/BDrawCache.h"
#include "BObol/BMeshLodCache.h"

#include "bg/trimesh.h"
#include "bu/app.h"
#include "bu/cache.h"
#include "bu/file.h"
#include "bu/hash.h"
#include "bu/malloc.h"
#include "bu/parallel.h"
#include "bu/path.h"
#include "bu/str.h"
#include "bu/datetime.h"
#include "bu/vls.h"
#include "raytrace.h"

#include <algorithm>
#include <array>
#include <cfloat>
#include <charconv>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <mutex>
#include <shared_mutex>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#define POP_CUT_COUNT_MAX BOBOL_MESH_LOD_CUT_COUNT_MAX
#define POP_CACHEDIR BOBOL_DRAW_CACHE_DIR
#define CACHE_CURRENT_FORMAT 18

#define CACHE_POP_MAX_CUT "max"
#define CACHE_POP_MIN_CUT "min"
#define CACHE_VERTEX_COUNT "vc"
#define CACHE_TRI_COUNT "tc"
#define CACHE_OBJ_BOUNDS "bb"
#define CACHE_SHADED_CULL_BACKFACES "cb"
#define CACHE_HAS_NORMALS "hn"
#define CACHE_CUT_QUANTIZATION "qb"
#define CACHE_VERT_CUT "v"
#define CACHE_VERTNORM_CUT "vn"
#define CACHE_TRI_CUT "t"
#define CACHE_CLUSTER_GRID "cg"
#define CACHE_CLUSTER_BOUNDS "gb"
#define CACHE_CLUSTER_RANGE_OFFSETS "go"
#define CACHE_CLUSTER_RANGES "gr"
#define CACHE_CHUNK_COUNT "kc"
#define CACHE_CHUNK_BOUNDS "kb"
#define CACHE_CHUNK_MINMAX "km"
#define CACHE_CHUNK_FACE_COUNTS "kf"
#define CACHE_CHUNK_POINT_COUNTS "kp"
#define CACHE_CHUNK_RESIDENT_BYTES "kx"
#define CACHE_CHUNK_DATA_PREFIX "k"

static const uint8_t meshLodChunkMagic[8] = {
    'B', 'O', 'B', 'C', 'H', 'N', 'K', '1'
};

/* Serialized cache records must not inherit compiler ABI padding or size_t
 * width.  Keep this format explicitly fixed so a hierarchy produced on one
 * supported host cannot turn a range offset into a wild mesh index on
 * another. */
struct BObolMeshLodClusterRangeDisk {
    uint32_t firstIndex;
    uint32_t indexCount;
    uint8_t activationCut;
    uint8_t reserved[3];
};
static_assert(sizeof(BObolMeshLodClusterRangeDisk) == 12,
    "cluster range cache record must have a fixed width");

struct BObolMeshLodContextInternal {
    struct bu_cache *lodCache;
    struct bu_cache *nameCache;
    struct bu_vls *fname;
    char *registryKey;
    std::shared_mutex *accessMutex;
    std::shared_mutex *nameMutex;
    std::unordered_map<std::string, unsigned long long> *nameKeys;
};

struct BObolMeshLodContext {
    struct BObolMeshLodContextInternal *i;
    size_t refs;
};

class BObolPopState;

struct BObolMeshLod {
    struct BObolMeshLodContext *context;
    BObolPopState *state;

    int fcnt;
    const uint32_t *faces;
    int pcnt;
    const point_t *points;
    int porig_cnt;
    const point_t *points_orig;
    const vect_t *normals;
    point_t bmin;
    point_t bmax;
};

class BObolPopRec
{
public:
    unsigned short x = 0;
    unsigned short y = 0;
    unsigned short z = 0;
};

struct BObolPopCut {
    uint8_t bits[3] = {0, 0, 0};
    double objectError = 0.0;
};

static int
mesh_lod_cache_write_semaphore(void)
{
    static int sem = 0;
    if (!sem)
	sem = bu_semaphore_register("BOBOL_MESH_LOD_CACHE_WRITE");
    return sem;
}

static std::mutex &
mesh_lod_context_registry_mutex(void)
{
    static std::mutex m;
    return m;
}

/*
 * A PoP cache for one real vehicle database can legitimately exceed libbu's
 * generic 4 GiB cache default.  LMDB's map size is an address-space ceiling,
 * not an eager disk allocation, so give mesh payloads a production-scale
 * ceiling while retaining an explicit deployment override.  Name->asset-key
 * records remain in the small generic cache below.
 */
static size_t
mesh_lod_cache_max_bytes(void)
{
    static const size_t bytes = []() {
	const size_t gibibyte = 1024ULL * 1024ULL * 1024ULL;
	const char *configured = getenv("BOBOL_MESH_LOD_CACHE_GB");
	if (configured && configured[0]) {
	    char *end = NULL;
	    const unsigned long long value = strtoull(configured, &end, 10);
	    if (end && end != configured && *end == '\0' && value > 0) {
		if (value > SIZE_MAX / gibibyte)
		    return SIZE_MAX - (SIZE_MAX % 4096u);
		return static_cast<size_t>(value) * gibibyte;
	    }
	}
	const unsigned long long defaultGibibytes = 64;
	if (defaultGibibytes > SIZE_MAX / gibibyte)
	    return SIZE_MAX - (SIZE_MAX % 4096u);
	return static_cast<size_t>(defaultGibibytes) * gibibyte;
    }();
    return bytes;
}

static std::map<std::string, struct BObolMeshLodContext *> &
mesh_lod_context_registry(void)
{
    static std::map<std::string, struct BObolMeshLodContext *> registry;
    return registry;
}

struct BObolMeshLodDatabaseContextHint {
    std::string databaseIdentity;
    std::string registryKey;
};

/*
 * These hints do not own contexts; the registry remains authoritative.  They
 * only avoid normalizing, hashing, and basename-formatting the same database
 * path for every asset opened by a many-leaf worker stream.  Bound the table
 * so applications which cycle through many transient databases cannot grow a
 * process-lifetime optimization cache without limit.
 */
static std::unordered_map<const struct db_i *,
    BObolMeshLodDatabaseContextHint> &
mesh_lod_database_context_hints(void)
{
    static std::unordered_map<const struct db_i *,
	BObolMeshLodDatabaseContextHint> hints;
    return hints;
}

static int
mesh_lod_size_to_int(size_t count, int *out)
{
    if (!out || count > static_cast<size_t>(std::numeric_limits<int>::max()))
	return 0;

    *out = static_cast<int>(count);
    return 1;
}

static int
mesh_lod_arrays_validate(const int *faces,
			 size_t faceCountIn,
			 const point_t *points,
			 size_t pointCountIn,
			 const point_t *pointsOrig,
			 size_t pointOrigCountIn,
			 int *faceCount,
			 int *pointCount,
			 int *pointOrigCount)
{
    int fcnt = 0;
    int pcnt = 0;
    int porigCnt = 0;

    if (!mesh_lod_size_to_int(faceCountIn, &fcnt) ||
	!mesh_lod_size_to_int(pointCountIn, &pcnt) ||
	!mesh_lod_size_to_int(pointOrigCountIn, &porigCnt))
	return 0;
    if (!faces || fcnt <= 0 || !points || pcnt <= 0 ||
	!pointsOrig || porigCnt <= 0)
	return 0;
    if (faceCountIn > ((size_t)-1) / 3)
	return 0;

    for (size_t pointIndex = 0; pointIndex < pointCountIn; ++pointIndex) {
	for (size_t axis = 0; axis < 3; ++axis) {
	    if (!std::isfinite(points[pointIndex][axis]))
		return 0;
	}
    }
    for (size_t pointIndex = 0; pointIndex < pointOrigCountIn;
	 pointIndex++) {
	for (size_t axis = 0; axis < 3; ++axis) {
	    if (!std::isfinite(pointsOrig[pointIndex][axis]))
		return 0;
	}
    }

    const size_t indexCount = faceCountIn * 3;
    for (size_t index = 0; index < indexCount; index++) {
	if (faces[index] < 0 || faces[index] >= pcnt ||
	    faces[index] >= porigCnt)
	    return 0;
    }

    if (faceCount)
	*faceCount = fcnt;
    if (pointCount)
	*pointCount = pcnt;
    if (pointOrigCount)
	*pointOrigCount = porigCnt;
    return 1;
}

static bool
mesh_lod_bot_authored_normals(std::vector<fastf_t> &normalStorage,
			      const struct rt_bot_internal *bot,
			      const int *displayFaces,
			      bool flippedWinding)
{
    normalStorage.clear();
    if (!bot || !(bot->bot_flags & RT_BOT_HAS_SURFACE_NORMALS) ||
	!(bot->bot_flags & RT_BOT_USE_NORMALS) || !bot->normals ||
	!bot->face_normals || bot->num_face_normals < bot->num_faces ||
	!displayFaces)
	return false;

    normalStorage.reserve(bot->num_faces * 9);
    for (size_t faceIndex = 0; faceIndex < bot->num_faces; ++faceIndex) {
	vect_t edgeA;
	vect_t edgeB;
	vect_t faceNormal;
	const int *face = &displayFaces[faceIndex * 3];
	VSUB2(edgeA, &bot->vertices[face[1] * 3],
	    &bot->vertices[face[0] * 3]);
	VSUB2(edgeB, &bot->vertices[face[2] * 3],
	    &bot->vertices[face[0] * 3]);
	VCROSS(faceNormal, edgeA, edgeB);
	if (MAGNITUDE(faceNormal) <= SMALL_FASTF) {
	    normalStorage.clear();
	    return false;
	}
	VUNITIZE(faceNormal);
	for (size_t corner = 0; corner < 3; ++corner) {
	    const size_t sourceCorner =
		(flippedWinding && corner > 0) ?
		3 - corner : corner;
	    const int normalIndex =
		bot->face_normals[faceIndex * 3 + sourceCorner];
	    if (normalIndex < 0 ||
		static_cast<size_t>(normalIndex) >= bot->num_normals) {
		normalStorage.clear();
		return false;
	    }
	    vect_t normal;
	    VMOVE(normal, &bot->normals[static_cast<size_t>(normalIndex) * 3]);
	    if (MAGNITUDE(normal) <= SMALL_FASTF) {
		normalStorage.clear();
		return false;
	    }
	    VUNITIZE(normal);
	    /* Normalize authored normals to the exterior-CCW display winding.
	     * A declared-CW BoT commonly stores normals in its source winding,
	     * while imported unoriented solids are less predictable. */
	    if (VDOT(normal, faceNormal) < 0.0)
		VREVERSE(normal, normal);
	    normalStorage.push_back(normal[X]);
	    normalStorage.push_back(normal[Y]);
	    normalStorage.push_back(normal[Z]);
	}
    }
    return true;
}

static void
mesh_lod_active_data_clear(struct BObolMeshLod *lod)
{
    if (!lod)
	return;

    lod->fcnt = 0;
    lod->faces = NULL;
    lod->pcnt = 0;
    lod->points = NULL;
    lod->porig_cnt = 0;
    lod->points_orig = NULL;
    lod->normals = NULL;
}

static void
mesh_lod_context_close(struct BObolMeshLodContext *context)
{
    if (!context)
	return;

    if (context->i) {
	if (context->i->nameCache)
	    bu_cache_close(context->i->nameCache);
	if (context->i->lodCache)
	    bu_cache_close(context->i->lodCache);
	if (context->i->fname) {
	    bu_vls_free(context->i->fname);
	    BU_PUT(context->i->fname, struct bu_vls);
	}
	if (context->i->registryKey)
	    bu_free(context->i->registryKey, "mesh lod context registry key");
	delete context->i->nameKeys;
	delete context->i->nameMutex;
	delete context->i->accessMutex;
	BU_PUT(context->i, struct BObolMeshLodContextInternal);
    }
    BU_PUT(context, struct BObolMeshLodContext);
}

static void
mesh_lod_context_destroy(struct BObolMeshLodContext *context)
{
    if (!context)
	return;

    std::lock_guard<std::mutex> guard(mesh_lod_context_registry_mutex());
    if (context->refs > 1) {
	context->refs--;
	return;
    }
    if (context->i && context->i->registryKey)
	mesh_lod_context_registry().erase(context->i->registryKey);
    context->refs = 0;
    /* Keep creation serialized until the final LMDB handles are fully
     * closed.  Releasing the registry lock first lets another worker open
     * the same cache while the previous context is still closing, which is
     * an intermittent failure on Windows. */
    mesh_lod_context_close(context);
}

static void
mesh_lod_cache_clear_context(struct BObolMeshLodContext *context,
			     unsigned long long key);

static struct BObolMeshLodContext *
mesh_lod_context_create(const char *name)
{
    if (!name)
	return NULL;

    /* bu_path_normalize returns one process-global scratch buffer.  Context
     * creation is already registry-serialized; normalize and copy inside that
     * critical section so concurrent LoD workers cannot corrupt one another's
     * cache identity before the registry lookup. */
    std::lock_guard<std::mutex> guard(mesh_lod_context_registry_mutex());
    const std::string normalizedPath(bu_path_normalize(name));
    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s", normalizedPath.c_str());
    if (bu_vls_strlen(&fname) < 10)
	bu_vls_printf(&fname, "GGGGGGGGGGGGG");
    unsigned long long hash = bu_data_hash(bu_vls_cstr(&fname),
					   bu_vls_strlen(&fname) * sizeof(char));
    bu_path_component(&fname, normalizedPath.c_str(),
	BU_PATH_BASENAME_EXTLESS);
    bu_vls_printf(&fname, "_%llu", hash);

    std::string registryKey(bu_vls_cstr(&fname));
    {
	auto it = mesh_lod_context_registry().find(registryKey);
	if (it != mesh_lod_context_registry().end()) {
	    it->second->refs++;
	    bu_vls_free(&fname);
	    return it->second;
	}
    }

    struct BObolMeshLodContext *context;
    BU_GET(context, struct BObolMeshLodContext);
    BU_GET(context->i, struct BObolMeshLodContextInternal);
    context->refs = 1;
    struct BObolMeshLodContextInternal *internal = context->i;
    BU_GET(internal->fname, struct bu_vls);
    bu_vls_init(internal->fname);
    bu_vls_sprintf(internal->fname, "%s", bu_vls_cstr(&fname));
    internal->registryKey = bu_strdup(registryKey.c_str());
    internal->lodCache = NULL;
    internal->nameCache = NULL;
    internal->nameKeys = NULL;
    internal->nameMutex = NULL;
    internal->accessMutex = new (std::nothrow) std::shared_mutex;
    if (!internal->accessMutex) {
	mesh_lod_context_close(context);
	bu_vls_free(&fname);
	return NULL;
    }
    internal->nameMutex = new (std::nothrow) std::shared_mutex;
    if (!internal->nameMutex) {
	mesh_lod_context_close(context);
	bu_vls_free(&fname);
	return NULL;
    }
    internal->nameKeys =
	new (std::nothrow) std::unordered_map<std::string,
	    unsigned long long>;
    if (!internal->nameKeys) {
	mesh_lod_context_close(context);
	bu_vls_free(&fname);
	return NULL;
    }

    char dir[MAXPATHLEN];
    bu_dir(dir, MAXPATHLEN, BU_DIR_CACHE, NULL);
    if (dir[0] && !bu_file_exists(dir, NULL))
	bu_mkdir(dir);
    bu_dir(dir, MAXPATHLEN, BU_DIR_CACHE, POP_CACHEDIR, NULL);
    if (!bu_file_exists(dir, NULL))
	bu_mkdir(dir);

    {
	char formatPath[MAXPATHLEN];
	bu_dir(formatPath, MAXPATHLEN, BU_DIR_CACHE, POP_CACHEDIR,
	       "mesh_lod.format", NULL);
	long diskFormatVersion = -1;
	{
	    std::ifstream formatFile(formatPath);
	    if (formatFile.is_open())
		formatFile >> diskFormatVersion;
	}
	if (diskFormatVersion > 0 && diskFormatVersion != CACHE_CURRENT_FORMAT) {
	    bu_log("Old mesh lod cache version (%ld) found in format file at %s - clearing\n",
		   diskFormatVersion, formatPath);
	    mesh_lod_cache_clear_context(NULL, 0);
	    bu_dir(dir, MAXPATHLEN, BU_DIR_CACHE, POP_CACHEDIR, NULL);
	    if (!bu_file_exists(dir, NULL))
		bu_mkdir(dir);
	}
	FILE *fp = fopen(formatPath, "w");
	if (fp) {
	    fprintf(fp, "%d\n", CACHE_CURRENT_FORMAT);
	    fclose(fp);
	}
    }

    struct bu_vls lodCachePath = BU_VLS_INIT_ZERO;
    struct bu_vls nameCachePath = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&lodCachePath, "%s/%s", POP_CACHEDIR, bu_vls_cstr(&fname));
    bu_vls_sprintf(&nameCachePath, "%s/%s_namekey", POP_CACHEDIR,
		   bu_vls_cstr(&fname));

    /*
     * View scheduling opens unrelated asset prefixes in screen-value order,
     * not LMDB page order.  Kernel readahead on a multi-gigabyte PoP cache
     * therefore loads speculative pages and can evict anonymous scene data
     * while a large database is also being streamed.  Both databases are
     * reconstructible caches, so make commits visible immediately but defer
     * physical durability to their clean close.  A cold many-asset build
     * otherwise forces two filesystem syncs per asset (payload plus name
     * mapping), serializing an otherwise parallel preparation pipeline.
     */
    internal->lodCache = bu_cache_open_with_options(
	bu_vls_cstr(&lodCachePath), 1, mesh_lod_cache_max_bytes(),
	BU_CACHE_OPEN_NORDAHEAD | BU_CACHE_OPEN_DEFER_SYNC);
    internal->nameCache = bu_cache_open_with_options(
	bu_vls_cstr(&nameCachePath), 1, 0, BU_CACHE_OPEN_DEFER_SYNC);

    bu_vls_free(&lodCachePath);
    bu_vls_free(&nameCachePath);
    bu_vls_free(&fname);

    if (!internal->lodCache || !internal->nameCache) {
	mesh_lod_context_close(context);
	return NULL;
    }

    mesh_lod_context_registry()[registryKey] = context;
    return context;
}

static struct BObolMeshLodContext *
mesh_lod_context_create_for_db(struct db_i *dbip)
{
    if (!dbip)
	return NULL;

    char ctxName[MAXPATHLEN] = {0};
    if (dbip->dbi_filename) {
	bu_strlcpy(ctxName, dbip->dbi_filename, sizeof(ctxName));
    } else {
	snprintf(ctxName, sizeof(ctxName), "bobol_inmem_mesh_lod_%p",
		 (void *)dbip);
    }

    {
	std::lock_guard<std::mutex> guard(mesh_lod_context_registry_mutex());
	const auto hint = mesh_lod_database_context_hints().find(dbip);
	if (hint != mesh_lod_database_context_hints().end() &&
	    hint->second.databaseIdentity == ctxName) {
	    const auto context =
		mesh_lod_context_registry().find(hint->second.registryKey);
	    if (context != mesh_lod_context_registry().end()) {
		context->second->refs++;
		return context->second;
	    }
	}
    }

    struct BObolMeshLodContext *context =
	mesh_lod_context_create(ctxName);
    if (!context || !context->i || !context->i->registryKey)
	return context;
    {
	std::lock_guard<std::mutex> guard(mesh_lod_context_registry_mutex());
	auto &hints = mesh_lod_database_context_hints();
	if (hints.size() >= 256 && hints.find(dbip) == hints.end())
	    hints.clear();
	BObolMeshLodDatabaseContextHint &hint = hints[dbip];
	hint.databaseIdentity = ctxName;
	hint.registryKey = context->i->registryKey;
    }
    return context;
}

static bool
mesh_lod_name_cache_key(char *buffer, size_t bufferSize, const char *name)
{
    if (!buffer || !bufferSize || !name)
	return false;

    const size_t nameLength = strlen(name);
    unsigned long long hash = 0;
    if (nameLength >= 10) {
	hash = bu_data_hash(name, nameLength);
    } else {
	char padded[32] = {0};
	static const char padding[] = "GGGGGGGGGGGGG";
	memcpy(padded, name, nameLength);
	memcpy(padded + nameLength, padding, sizeof(padding) - 1);
	hash = bu_data_hash(
	    padded, nameLength + sizeof(padding) - 1);
    }

    const std::to_chars_result converted =
	std::to_chars(buffer, buffer + bufferSize, hash);
    static const char suffix[] = ":namekey";
    if (converted.ec != std::errc() ||
	static_cast<size_t>(buffer + bufferSize - converted.ptr) <
	    sizeof(suffix))
	return false;
    memcpy(converted.ptr, suffix, sizeof(suffix));
    return true;
}

static unsigned long long
mesh_lod_key_get(struct BObolMeshLodContext *context, const char *name)
{
    if (!context || !name)
	return 0;

    {
	std::shared_lock<std::shared_mutex> lock(
	    *context->i->nameMutex);
	const auto found = context->i->nameKeys->find(name);
	if (found != context->i->nameKeys->end())
	    return found->second;
    }

    char keystr[32] = {0};
    if (!mesh_lod_name_cache_key(keystr, sizeof(keystr), name))
	return 0;

    void *data = NULL;
    std::shared_lock<std::shared_mutex> lock(*context->i->accessMutex);
    size_t dsize = bu_cache_get(&data, keystr,
				context->i->nameCache, NULL);
    lock.unlock();
    if (dsize != sizeof(uint64_t) || !data) {
	if (data)
	    bu_free(data, "invalid lod key data");
	return 0;
    }

    uint64_t diskKey = 0;
    memcpy(&diskKey, data, sizeof(diskKey));
    bu_free(data, "lod key data");

    const unsigned long long meshKey =
	static_cast<unsigned long long>(diskKey);
    /*
     * The on-disk name map is authoritative.  Caching a first lookup in the
     * process map is optional; never stall seven other cache readers merely
     * to install one hint while a 100k-asset warm scene is fanning out.
     */
    std::unique_lock<std::shared_mutex> memoryLock(
	*context->i->nameMutex, std::try_to_lock);
    if (memoryLock.owns_lock())
	context->i->nameKeys->emplace(name, meshKey);
    return meshKey;
}

static int
mesh_lod_key_put(struct BObolMeshLodContext *context,
		 const char *name,
		 unsigned long long key)
{
    if (!context || !name || !key)
	return -1;

    char keystr[32] = {0};
    if (!mesh_lod_name_cache_key(keystr, sizeof(keystr), name))
	return -1;

    const uint64_t diskKey = static_cast<uint64_t>(key);
    std::unique_lock<std::shared_mutex> lock(*context->i->accessMutex);
    size_t wsize = bu_cache_write((void *)&diskKey, sizeof(diskKey),
				  keystr,
				  context->i->nameCache, NULL);
    if (wsize > 0) {
	std::unique_lock<std::shared_mutex> memoryLock(
	    *context->i->nameMutex);
	(*context->i->nameKeys)[name] = key;
    }
    lock.unlock();
    return (wsize > 0) ? 0 : -1;
}

class BObolPopState
{
public:
    BObolPopState(struct BObolMeshLodContext *ctx,
		    const point_t *vertices,
		    size_t vertexCount,
		    const vect_t *normals,
		    const int *faces,
		    size_t faceCount,
		    unsigned long long userKey,
		    bool cullBackfaces);
    BObolPopState(struct BObolMeshLodContext *ctx,
		    unsigned long long key);
    BObolPopState(struct BObolMeshLodContext *ctx,
		    unsigned long long key,
		    bool retainHeaderSnapshot);
    ~BObolPopState();

    bool setCut(int cut, bool materializeSnapped = true);
    bool materializeInitialPrefix(int cut);
    void releaseGenerationScratch(void);
    bool readSuffix(int residentCut, int targetCut,
	BObolMeshLodSuffixCallback callback, void *callbackData);
    bool readChunks(const uint32_t *chunkIds, size_t chunkCount, int cut,
	BObolMeshLodChunkCallback callback, void *callbackData);
    void shrinkMemory(void);
    size_t residentBytes(void) const;
    size_t residentPrefixBytes(void) const;
    void cutPoint(point_t *out, const point_t *point, int cut);
    void hierarchyInfo(struct BObolMeshLodHierarchyInfo *info) const;

    std::vector<uint32_t> lodTris;
    std::vector<fastf_t> lodTriPoints;
    std::vector<fastf_t> lodTriPointsSnapped;
    std::vector<fastf_t> lodTriNormals;

    int currCut = -1;
    bool forceUpdate = false;
    int minPopCut = 0;
    int maxPopCut = 0;
    uint32_t cutCount = 0;
    bool isValid = false;
    bool shadedCullBackfaces = false;
    bool hasNormals = false;
    unsigned long long hash = 0;
    point_t bbmin = VINIT_ZERO;
    point_t bbmax = VINIT_ZERO;
    struct BObolMeshLod *lod = NULL;

private:
    void buildCutSchedule(void);
    void triProcess(void);
    bool buildClusters(void);
    bool buildChunks(void);
    bool cacheChunks(void);
    fastf_t snap(fastf_t value, fastf_t min, fastf_t max, uint8_t bits);
    void cache(void);
    bool cacheTri(void);
    bool cacheWrite(const char *component, std::stringstream &stream);
    bool cacheWriteData(const char *component, const void *data, size_t size);
    size_t cacheGet(void **data, const char *component);
    void cacheDone(void);
    bool triPopLoad(int startCut, int cut, bool materializeSnapped);
    void triPopTrim(int cut, bool materializeSnapped);
    void updateSnappedPoints(int cut);
    bool loadCachedHeader(bool retainHeaderSnapshot);
    void initializeCacheKeyPrefix(void);
    bool cacheComponentKey(char *buffer, size_t bufferSize,
	const char *component) const;

    fastf_t minx = std::numeric_limits<fastf_t>::max();
    fastf_t miny = std::numeric_limits<fastf_t>::max();
    fastf_t minz = std::numeric_limits<fastf_t>::max();
    fastf_t maxx = -std::numeric_limits<fastf_t>::max();
    fastf_t maxy = -std::numeric_limits<fastf_t>::max();
    fastf_t maxz = -std::numeric_limits<fastf_t>::max();
    BObolPopCut cuts[POP_CUT_COUNT_MAX];
    size_t cutVertexCount[POP_CUT_COUNT_MAX + 1] = {0};
    size_t cutTriangleCount[POP_CUT_COUNT_MAX + 1] = {0};

    /* Every source vertex and face belongs to exactly one activation cut.
     * Dense cut vectors are substantially cheaper than the old
     * map<cut, unordered_set<vertex>> representation on scan meshes: Lucy's
     * 14 million vertices otherwise spent more memory on hash nodes than on
     * the source coordinates themselves.  BoT indices are int-sized by
     * contract, so 32-bit remap entries are sufficient. */
    std::vector<uint32_t> triIndexMap;
    std::vector<uint8_t> vertexTriMinCut;
    std::vector<uint8_t> faceActivationCut;
    std::vector<std::vector<uint32_t>> cutTriVerts;
    std::vector<std::vector<uint32_t>> cutTris;
    std::vector<struct BObolMeshLodClusterInfo> clusterInfos;
    std::vector<struct BObolMeshLodClusterRange> clusterRanges;
    std::vector<std::vector<uint32_t>> chunkFaces;
    std::vector<struct BObolMeshLodChunkInfo> chunkInfos;

    size_t vertexCount = 0;
    const point_t *vertexArray = NULL;
    const vect_t *normalArray = NULL;
    size_t faceCount = 0;
    const int *faceArray = NULL;
    /* Authored BoTs may contain zero-area faces with a repeated vertex
     * index.  They contribute no pixels, but retaining them in a progressive
     * prefix can make the immutable mesh payload reject its own index stream.
     * Allocate normalized storage only for that exceptional input; ordinary
     * production meshes continue to borrow the staged arrays directly. */
    std::vector<int> normalizedFaceArray;
    std::vector<fastf_t> normalizedNormalArray;

    struct BObolMeshLodContext *context = NULL;
    struct bu_cache_txn *readTxn = NULL;
    struct bu_cache_txn *writeTxn = NULL;
    std::shared_lock<std::shared_mutex> readLock;
    char cacheKeyPrefix[32] = {0};
    size_t cacheKeyPrefixLength = 0;
};

void
BObolPopState::buildCutSchedule(void)
{
    const double extent[3] = {
	static_cast<double>(maxx - minx),
	static_cast<double>(maxy - miny),
	static_cast<double>(maxz - minz)
    };
    uint8_t bits[3] = {1, 1, 1};
    for (int axis = 0; axis < 3; ++axis)
	if (!(extent[axis] > 0.0))
	    bits[axis] = BOBOL_MESH_LOD_QUANTIZATION_BITS;

    cutCount = 0;
    while (cutCount < POP_CUT_COUNT_MAX) {
	BObolPopCut &cut = cuts[cutCount++];
	for (int axis = 0; axis < 3; ++axis)
	    cut.bits[axis] = bits[axis];
	if (bits[0] >= BOBOL_MESH_LOD_QUANTIZATION_BITS &&
	    bits[1] >= BOBOL_MESH_LOD_QUANTIZATION_BITS &&
	    bits[2] >= BOBOL_MESH_LOD_QUANTIZATION_BITS) {
	    cut.objectError = 0.0;
	    break;
	}

	double squaredError = 0.0;
	for (int axis = 0; axis < 3; ++axis) {
	    const double mask = std::ldexp(1.0,
		BOBOL_MESH_LOD_QUANTIZATION_BITS - bits[axis]);
	    const double axisError =
		0.5 * mask * extent[axis] / 65535.0;
	    squaredError += axisError * axisError;
	}
	cut.objectError = std::sqrt(squaredError);

	int refineAxis = -1;
	double largestAxisError = -1.0;
	for (int axis = 0; axis < 3; ++axis) {
	    if (bits[axis] >= BOBOL_MESH_LOD_QUANTIZATION_BITS)
		continue;
	    const double axisError =
		extent[axis] / std::ldexp(1.0, bits[axis]);
	    if (axisError > largestAxisError) {
		largestAxisError = axisError;
		refineAxis = axis;
	    }
	}
	if (refineAxis < 0)
	    break;
	++bits[refineAxis];
    }
    maxPopCut = cutCount ? static_cast<int>(cutCount) - 1 : 0;
}

void
BObolPopState::triProcess(void)
{
    vertexTriMinCut.assign(vertexCount,
	UINT8_MAX);
    cutTriVerts.resize(cutCount);
    cutTris.resize(cutCount);

    /* Face activation is independent for every triangle.  The old loop did
     * nine floating-point divisions and all activation work serially before
     * any cold PoP cut could become usable.  Classify large meshes in bounded
     * disjoint ranges, recording only one byte per source face.  The merge
     * below remains serial and source ordered, preserving byte-for-byte
     * deterministic cache topology regardless of worker count. */
    faceActivationCut.assign(faceCount, UINT8_MAX);
    const fastf_t extent[3] = {
	maxx - minx, maxy - miny, maxz - minz
    };
    const fastf_t minimum[3] = {minx, miny, minz};
    fastf_t quantizeScale[3] = {0.0, 0.0, 0.0};
    for (int axis = 0; axis < 3; ++axis) {
	if (extent[axis] > 0.0)
	    quantizeScale[axis] =
		static_cast<fastf_t>(USHRT_MAX) / extent[axis];
    }
    const auto quantize = [&minimum, &extent, &quantizeScale](
	    fastf_t value, int axis) -> unsigned short {
	if (!(extent[axis] > 0.0))
	    return 0;
	if (value <= minimum[axis])
	    return 0;
	if (value >= minimum[axis] + extent[axis])
	    return USHRT_MAX;
	return static_cast<unsigned short>(
	    floor((value - minimum[axis]) * quantizeScale[axis]));
    };
    uint8_t firstCutForBits[3][BOBOL_MESH_LOD_QUANTIZATION_BITS + 1];
    memset(firstCutForBits, UINT8_MAX, sizeof(firstCutForBits));
    for (int axis = 0; axis < 3; ++axis) {
	firstCutForBits[axis][0] = 0;
	for (uint32_t cut = 0; cut < cutCount; ++cut)
	    for (int bits = 1;
		    bits <= cuts[cut].bits[axis]; ++bits)
		if (firstCutForBits[axis][bits] == UINT8_MAX)
		    firstCutForBits[axis][bits] = static_cast<uint8_t>(cut);
    }
    const auto differingBits = [](unsigned short a,
	    unsigned short b) -> int {
	const unsigned int differing = static_cast<unsigned int>(a ^ b);
	if (!differing)
	    return 0;
#if defined(__GNUC__) || defined(__clang__)
	const int mostSignificantBit = 31 - __builtin_clz(differing);
#else
	int mostSignificantBit = 0;
	for (unsigned int bits = differing; bits >>= 1;)
	    ++mostSignificantBit;
#endif
	return BOBOL_MESH_LOD_QUANTIZATION_BITS - mostSignificantBit;
    };
    const auto pairCut = [this, &firstCutForBits, &differingBits](
	    const BObolPopRec &a, const BObolPopRec &b) -> int {
	const unsigned short av[3] = {a.x, a.y, a.z};
	const unsigned short bv[3] = {b.x, b.y, b.z};
	int cut = maxPopCut;
	bool distinguishable = false;
	for (int axis = 0; axis < 3; ++axis) {
	    const int bits = differingBits(av[axis], bv[axis]);
	    if (!bits)
		continue;
	    distinguishable = true;
	    cut = std::min(cut,
		static_cast<int>(firstCutForBits[axis][bits]));
	}
	return distinguishable ? cut : maxPopCut;
    };
    const auto classifyRange = [this, &quantize, &pairCut](
	    size_t begin, size_t end) {
	for (size_t faceIndex = begin; faceIndex < end; ++faceIndex) {
	BObolPopRec triangle[3];
	bool badFace = false;
	for (size_t cornerIndex = 0; cornerIndex < 3; cornerIndex++) {
	    int faceVertex = faceArray[3 * faceIndex + cornerIndex];
	    if (faceVertex < 0 || static_cast<size_t>(faceVertex) >= vertexCount) {
		badFace = true;
		break;
	    }
	    triangle[cornerIndex].x =
		quantize(vertexArray[faceVertex][X], X);
	    triangle[cornerIndex].y =
		quantize(vertexArray[faceVertex][Y], Y);
	    triangle[cornerIndex].z =
		quantize(vertexArray[faceVertex][Z], Z);
	}
	if (badFace)
	    continue;

	const int cut = std::max(pairCut(triangle[0], triangle[1]),
	    std::max(pairCut(triangle[1], triangle[2]),
		     pairCut(triangle[0], triangle[2])));
	faceActivationCut[faceIndex] = static_cast<uint8_t>(cut);
	}
    };

    size_t workerCount = 1;
    if (faceCount >= 1000000) {
	const size_t usefulWorkers =
	    std::max<size_t>(1, (faceCount + 999999) / 1000000);
	workerCount = std::min<size_t>(
	    std::max<size_t>(1, bu_avail_cpus()), usefulWorkers);
    }
    if (workerCount <= 1) {
	classifyRange(0, faceCount);
    } else {
	std::vector<std::thread> workers;
	workers.reserve(workerCount);
	const size_t chunk = (faceCount + workerCount - 1) / workerCount;
	for (size_t workerIndex = 0; workerIndex < workerCount; ++workerIndex) {
	    const size_t begin = workerIndex * chunk;
	    const size_t end = std::min(faceCount, begin + chunk);
	    if (begin >= end)
		break;
	    workers.emplace_back(classifyRange, begin, end);
	}
	for (std::thread &worker : workers)
	    worker.join();
    }

    size_t cutFaceCounts[POP_CUT_COUNT_MAX] = {0};
    size_t badFaceCount = 0;
    for (uint8_t cut : faceActivationCut) {
	if (cut >= cutCount)
	    badFaceCount++;
	else
	    cutFaceCounts[cut]++;
    }
    for (uint32_t cut = 0; cut < cutCount; ++cut)
	cutTris[cut].reserve(cutFaceCounts[cut]);

    for (size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex) {
	const uint8_t cut = faceActivationCut[faceIndex];
	if (cut >= cutCount)
	    continue;
	cutTris[cut].push_back(static_cast<uint32_t>(faceIndex));
	for (size_t cornerIndex = 0; cornerIndex < 3; cornerIndex++) {
	    int faceVertex = faceArray[3 * faceIndex + cornerIndex];
	    if (vertexTriMinCut[faceVertex] > cut)
		vertexTriMinCut[faceVertex] = cut;
	}
    }
    if (badFaceCount)
	bu_log("%zu invalid mesh faces skipped during PoP classification\n",
	       badFaceCount);

    for (size_t vertexIndex = 0; vertexIndex < vertexTriMinCut.size();
	 vertexIndex++) {
	const uint8_t activation = vertexTriMinCut[vertexIndex];
	/* Unreferenced vertices do not participate in any rendered triangle.
	 * Assigning them to the terminal global prefix while spatial chunks omit
	 * them violates the cache population conservation proof and wastes memory
	 * on authored construction debris. */
	if (activation < cutCount)
	    cutTriVerts[activation].push_back(
		static_cast<uint32_t>(vertexIndex));
    }

    for (uint32_t cut = 0; cut < cutCount; ++cut) {
	cutVertexCount[cut] = cutTriVerts[cut].size();
	cutTriangleCount[cut] = cutTris[cut].size();
    }

    triIndexMap.resize(vertexCount);

    uint32_t outputVertexIndex = 0;
    for (const std::vector<uint32_t> &cutVertices : cutTriVerts) {
	for (uint32_t sourceVertexIndex : cutVertices) {
	    triIndexMap[sourceVertexIndex] = outputVertexIndex;
	    ++outputVertexIndex;
	}
    }

    int firstPopulatedCut = -1;
    for (size_t cut = 0; cut < cutTris.size(); ++cut) {
	if (!cutTris[cut].empty()) {
	    if (firstPopulatedCut < 0)
		firstPopulatedCut = static_cast<int>(cut);
	}
    }

    /* Face activation and coordinate precision are independent progress axes.
     * A mesh may have every face active at a low cut while its vertices are
     * still snapped to a visibly coarse grid.  The last topology cut is
     * therefore not a valid display terminal.  Persist the cumulative
     * topology through the full 16-bit precision range; cuts with no new
     * faces or vertices are cheap and still re-snap the resident coordinates
     * to a finer grid. */
    minPopCut = firstPopulatedCut >= 0 ? firstPopulatedCut : 0;
}

bool
BObolPopState::buildClusters(void)
{
    constexpr size_t side = BOBOL_MESH_LOD_CLUSTER_GRID_RESOLUTION;
    constexpr size_t clusterCount = BOBOL_MESH_LOD_CLUSTER_COUNT;
    static_assert(side > 1 && side * side * side == clusterCount,
	"invalid fixed PoP cluster grid");

    clusterInfos.clear();
    clusterRanges.clear();
    if (!vertexArray || !faceArray || !faceCount ||
	faceCount > UINT32_MAX / 3 || cutTris.size() != cutCount)
	return false;

    struct ClusterBounds {
	double minimum[3] = {
	    std::numeric_limits<double>::infinity(),
	    std::numeric_limits<double>::infinity(),
	    std::numeric_limits<double>::infinity()
	};
	double maximum[3] = {
	    -std::numeric_limits<double>::infinity(),
	    -std::numeric_limits<double>::infinity(),
	    -std::numeric_limits<double>::infinity()
	};
    };
    std::array<ClusterBounds, clusterCount> bounds;
    std::array<std::vector<BObolMeshLodClusterRange>, clusterCount>
	rangesByCluster;
    const double minimum[3] = {minx, miny, minz};
    const double extent[3] = {maxx - minx, maxy - miny, maxz - minz};

    const auto faceCluster = [this, &minimum, &extent](
	    uint32_t sourceFace) -> size_t {
	size_t cell[3] = {0, 0, 0};
	for (int axis = 0; axis < 3; ++axis) {
	    double centroid = 0.0;
	    for (size_t corner = 0; corner < 3; ++corner) {
		const int vertex = faceArray[
		    3 * static_cast<size_t>(sourceFace) + corner];
		centroid += static_cast<double>(vertexArray[vertex][axis]);
	    }
	    centroid /= 3.0;
	    if (extent[axis] > 0.0) {
		const double normalized =
		    (centroid - minimum[axis]) / extent[axis];
		const double scaled = floor(normalized * side);
		cell[axis] = scaled <= 0.0 ? 0 :
		    static_cast<size_t>(std::min<double>(side - 1, scaled));
	    }
	}
	return cell[X] + side * (cell[Y] + side * cell[Z]);
    };

    uint64_t cumulativeFaces = 0;
    for (uint32_t cut = 0; cut < cutCount; ++cut) {
	std::vector<uint32_t> &sourceFaces = cutTris[cut];
	std::array<size_t, clusterCount> counts = {};
	for (uint32_t sourceFace : sourceFaces) {
	    if (sourceFace >= faceCount)
		return false;
	    const size_t cluster = faceCluster(sourceFace);
	    ++counts[cluster];
	    for (size_t corner = 0; corner < 3; ++corner) {
		const int vertex = faceArray[
		    3 * static_cast<size_t>(sourceFace) + corner];
		if (vertex < 0 || static_cast<size_t>(vertex) >= vertexCount)
		    return false;
		for (int axis = 0; axis < 3; ++axis) {
		    const double value = vertexArray[vertex][axis];
		    bounds[cluster].minimum[axis] =
			std::min(bounds[cluster].minimum[axis], value);
		    bounds[cluster].maximum[axis] =
			std::max(bounds[cluster].maximum[axis], value);
		}
	    }
	}

	std::array<size_t, clusterCount> offsets = {};
	size_t next = 0;
	for (size_t cluster = 0; cluster < clusterCount; ++cluster) {
	    offsets[cluster] = next;
	    if (counts[cluster]) {
		const uint64_t first = (cumulativeFaces + next) * 3;
		const uint64_t count = counts[cluster] * 3;
		if (first > UINT32_MAX || count > UINT32_MAX ||
		    first + count > UINT32_MAX)
		    return false;
		rangesByCluster[cluster].push_back({
		    static_cast<uint32_t>(first),
		    static_cast<uint32_t>(count),
		    static_cast<uint8_t>(cut)
		});
	    }
	    next += counts[cluster];
	}
	if (next != sourceFaces.size())
	    return false;

	/* Counting-sort by cell within this activation cut.  Cut order, point
	 * activation order, and every cumulative prefix remain unchanged. */
	std::vector<uint32_t> ordered(sourceFaces.size());
	std::array<size_t, clusterCount> cursor = offsets;
	for (uint32_t sourceFace : sourceFaces) {
	    const size_t cluster = faceCluster(sourceFace);
	    ordered[cursor[cluster]++] = sourceFace;
	}
	sourceFaces.swap(ordered);
	cumulativeFaces += sourceFaces.size();
    }

    size_t totalRangeCount = 0;
    for (const auto &ranges : rangesByCluster)
	totalRangeCount += ranges.size();
    if (!totalRangeCount || totalRangeCount > clusterCount * cutCount)
	return false;
    clusterRanges.reserve(totalRangeCount);
    clusterInfos.resize(clusterCount);
    for (size_t cluster = 0; cluster < clusterCount; ++cluster) {
	BObolMeshLodClusterInfo &info = clusterInfos[cluster];
	const size_t firstRange = clusterRanges.size();
	clusterRanges.insert(clusterRanges.end(),
	    rangesByCluster[cluster].begin(), rangesByCluster[cluster].end());
	info.range_count = static_cast<uint32_t>(
	    clusterRanges.size() - firstRange);
	if (!info.range_count) {
	    VSETALL(info.bmin, 0.0);
	    VSETALL(info.bmax, 0.0);
	    info.ranges = NULL;
	    continue;
	}
	VSET(info.bmin, bounds[cluster].minimum[X],
	    bounds[cluster].minimum[Y], bounds[cluster].minimum[Z]);
	VSET(info.bmax, bounds[cluster].maximum[X],
	    bounds[cluster].maximum[Y], bounds[cluster].maximum[Z]);
    }
    /* Assign borrowed pointers only after the flat vector has reached its
	 * final capacity. */
    size_t rangeOffset = 0;
    for (BObolMeshLodClusterInfo &info : clusterInfos) {
	if (info.range_count) {
	    info.ranges = clusterRanges.data() + rangeOffset;
	    rangeOffset += info.range_count;
	}
    }
    return rangeOffset == clusterRanges.size();
}

bool
BObolPopState::buildChunks(void)
{
    constexpr size_t side = BOBOL_MESH_LOD_CLUSTER_GRID_RESOLUTION;
    constexpr size_t cellCount = BOBOL_MESH_LOD_CLUSTER_COUNT;
    constexpr size_t faceTarget = BOBOL_MESH_LOD_CHUNK_FACE_TARGET;

    chunkFaces.clear();
    chunkInfos.clear();
    if (!vertexArray || !faceArray || !faceCount || !faceTarget ||
	faceActivationCut.size() != faceCount || cutCount == 0)
	return false;

    const double minimum[3] = {minx, miny, minz};
    const double extent[3] = {maxx - minx, maxy - miny, maxz - minz};
    const auto faceCell = [this, &minimum, &extent](
	    uint32_t sourceFace) -> size_t {
	size_t cell[3] = {0, 0, 0};
	for (int axis = 0; axis < 3; ++axis) {
	    double centroid = 0.0;
	    for (size_t corner = 0; corner < 3; ++corner) {
		const int vertex = faceArray[
		    3 * static_cast<size_t>(sourceFace) + corner];
		centroid += static_cast<double>(vertexArray[vertex][axis]);
	    }
	    centroid /= 3.0;
	    if (extent[axis] > 0.0) {
		const double normalized =
		    (centroid - minimum[axis]) / extent[axis];
		const double scaled = floor(normalized * side);
		cell[axis] = scaled <= 0.0 ? 0 :
		    static_cast<size_t>(std::min<double>(side - 1, scaled));
	    }
	}
	return cell[X] + side * (cell[Y] + side * cell[Z]);
    };

    std::array<std::vector<uint32_t>, cellCount> cells;
    size_t validFaces = 0;
    for (size_t face = 0; face < faceCount; ++face) {
	if (faceActivationCut[face] >= cutCount)
	    continue;
	const uint32_t sourceFace = static_cast<uint32_t>(face);
	cells[faceCell(sourceFace)].push_back(sourceFace);
	++validFaces;
    }
    if (!validFaces)
	return false;

    /* Preserve spatial cell order while packing bounded pages.  A fixed
     * grid alone creates one enormous page for a scan concentrated in one
     * cell; one page per nonempty cell creates hundreds of tiny pages for a
     * modest mesh.  Packing and splitting by face count provides both a
     * deterministic spatial order and a hard preparation/IO bound. */
    for (std::vector<uint32_t> &cell : cells) {
	size_t offset = 0;
	while (offset < cell.size()) {
	    if (chunkFaces.empty() ||
		chunkFaces.back().size() >= faceTarget)
		chunkFaces.emplace_back();
	    std::vector<uint32_t> &chunk = chunkFaces.back();
	    const size_t available = faceTarget - chunk.size();
	    const size_t take = std::min(available, cell.size() - offset);
	    chunk.insert(chunk.end(), cell.begin() + offset,
		cell.begin() + offset + take);
	    offset += take;
	}
	std::vector<uint32_t>().swap(cell);
    }
    if (chunkFaces.empty() || chunkFaces.size() > UINT32_MAX)
	return false;

    chunkInfos.resize(chunkFaces.size());
    for (size_t chunkIndex = 0; chunkIndex < chunkFaces.size();
	 chunkIndex++) {
	BObolMeshLodChunkInfo &info = chunkInfos[chunkIndex];
	memset(&info, 0, sizeof(info));
	info.chunk_id = static_cast<uint32_t>(chunkIndex);
	info.min_cut = maxPopCut;
	info.max_cut = maxPopCut;
	double chunkMin[3] = {
	    std::numeric_limits<double>::infinity(),
	    std::numeric_limits<double>::infinity(),
	    std::numeric_limits<double>::infinity()
	};
	double chunkMax[3] = {
	    -std::numeric_limits<double>::infinity(),
	    -std::numeric_limits<double>::infinity(),
	    -std::numeric_limits<double>::infinity()
	};
	for (uint32_t sourceFace : chunkFaces[chunkIndex]) {
	    if (sourceFace >= faceCount ||
		faceActivationCut[sourceFace] >= cutCount)
		return false;
	    info.min_cut = std::min(info.min_cut,
		static_cast<int>(faceActivationCut[sourceFace]));
	    for (size_t corner = 0; corner < 3; ++corner) {
		const int vertex = faceArray[
		    3 * static_cast<size_t>(sourceFace) + corner];
		if (vertex < 0 || static_cast<size_t>(vertex) >= vertexCount)
		    return false;
		for (int axis = 0; axis < 3; ++axis) {
		    const double value = vertexArray[vertex][axis];
		    chunkMin[axis] = std::min(chunkMin[axis], value);
		    chunkMax[axis] = std::max(chunkMax[axis], value);
		}
	    }
	}
	if (!std::isfinite(chunkMin[X]) || !std::isfinite(chunkMax[X]))
	    return false;
	VSET(info.bmin, chunkMin[X], chunkMin[Y], chunkMin[Z]);
	VSET(info.bmax, chunkMax[X], chunkMax[Y], chunkMax[Z]);
    }
    return true;
}

BObolPopState::BObolPopState(
    struct BObolMeshLodContext *ctx,
    const point_t *vertices,
    size_t inputVertexCount,
    const vect_t *normals,
    const int *faces,
    size_t inputFaceCount,
    unsigned long long userKey,
    bool cullBackfaces)
{
    context = ctx;
    shadedCullBackfaces = cullBackfaces;
    hasNormals = normals != NULL;

    vertexCount = inputVertexCount;
    vertexArray = vertices;
    normalArray = normals;
    faceCount = inputFaceCount;
    faceArray = faces;

    size_t drawableFaceCount = 0;
    for (size_t faceIndex = 0; faces && faceIndex < inputFaceCount;
	 faceIndex++) {
	const int a = faces[3 * faceIndex];
	const int b = faces[3 * faceIndex + 1];
	const int c = faces[3 * faceIndex + 2];
	if (a != b && b != c && a != c)
	    drawableFaceCount++;
    }
    if (faces && drawableFaceCount != inputFaceCount) {
	normalizedFaceArray.reserve(drawableFaceCount * 3);
	if (normals)
	    normalizedNormalArray.reserve(drawableFaceCount * 9);
	for (size_t faceIndex = 0; faceIndex < inputFaceCount; faceIndex++) {
	    const int a = faces[3 * faceIndex];
	    const int b = faces[3 * faceIndex + 1];
	    const int c = faces[3 * faceIndex + 2];
	    if (a == b || b == c || a == c)
		continue;
	    normalizedFaceArray.push_back(a);
	    normalizedFaceArray.push_back(b);
	    normalizedFaceArray.push_back(c);
	    if (normals) {
		for (size_t corner = 0; corner < 3; corner++)
		    for (size_t axis = 0; axis < 3; axis++)
			normalizedNormalArray.push_back(
			    normals[3 * faceIndex + corner][axis]);
	    }
	}
	faceCount = drawableFaceCount;
	faceArray = normalizedFaceArray.empty() ? NULL :
	    normalizedFaceArray.data();
	if (normals)
	    normalArray = normalizedNormalArray.empty() ? NULL :
		reinterpret_cast<const vect_t *>(
		    normalizedNormalArray.data());
	if (getenv("BOBOL_DRAW_TIMING"))
	    bu_log("[obol-timing] pop cache: discarded %zu "
		   "topologically degenerate faces\n",
		   inputFaceCount - drawableFaceCount);
    }
    if (!vertexArray || !vertexCount || !faceArray || !faceCount)
	return;

    if (!userKey) {
	struct bu_data_hash_state *state = bu_data_hash_create();
	static const char semantics[] = "BObol-chunked-PoP-format-19";
	bu_data_hash_update(state, semantics, sizeof(semantics));
	bu_data_hash_update(state, vertexArray, vertexCount * sizeof(point_t));
	bu_data_hash_update(state, faceArray, 3 * faceCount * sizeof(int));
	const unsigned char normalFlag = normals ? 1u : 0u;
	bu_data_hash_update(state, &normalFlag, sizeof(normalFlag));
	if (normalArray)
	    bu_data_hash_update(state, normalArray,
		3 * faceCount * sizeof(vect_t));
	const unsigned char cullFlag = cullBackfaces ? 1u : 0u;
	bu_data_hash_update(state, &cullFlag, sizeof(cullFlag));
	hash = bu_data_hash_val(state);
	bu_data_hash_destroy(state);
    } else {
	hash = userKey;
    }
    initializeCacheKeyPrefix();

    /* A content-identical hierarchy may already have been generated for a
     * different database name.  A marker proves only that some record is
     * present; it is not a usable hierarchy state.  Validate the complete
     * immutable header before reusing it so generate-and-open can materialize
     * its first prefix without reclassifying the mesh or publishing an
     * uninitialized shell. */
    void *cacheData = NULL;
    const size_t cacheSize = cacheGet(&cacheData, CACHE_POP_MAX_CUT);
    const bool cachedMarkerPresent = cacheSize && cacheData;
    cacheDone();
    if (cachedMarkerPresent && loadCachedHeader(false))
	return;

    /* A partial same-key cache must not leak metadata into the new producer
     * run.  Complete generation below rewrites the immutable records. */
    isValid = false;
    currCut = -1;
    minPopCut = 0;
    maxPopCut = 0;
    cutCount = 0;
    memset(cutVertexCount, 0, sizeof(cutVertexCount));
    memset(cutTriangleCount, 0, sizeof(cutTriangleCount));
    const int64_t buildStarted = bu_gettime();
    bg_trimesh_aabb(&bbmin, &bbmax, const_cast<int *>(faceArray), faceCount,
		    vertexArray, vertexCount);

    for (size_t vertexIndex = 0; vertexIndex < inputVertexCount;
	 vertexIndex++) {
	minx = (vertices[vertexIndex][X] < minx) ? vertices[vertexIndex][X] : minx;
	miny = (vertices[vertexIndex][Y] < miny) ? vertices[vertexIndex][Y] : miny;
	minz = (vertices[vertexIndex][Z] < minz) ? vertices[vertexIndex][Z] : minz;
	maxx = (vertices[vertexIndex][X] > maxx) ? vertices[vertexIndex][X] : maxx;
	maxy = (vertices[vertexIndex][Y] > maxy) ? vertices[vertexIndex][Y] : maxy;
	maxz = (vertices[vertexIndex][Z] > maxz) ? vertices[vertexIndex][Z] : maxz;
    }

    if (!std::isfinite(maxx - minx) || !std::isfinite(maxy - miny) ||
	!std::isfinite(maxz - minz))
	return;

    const int64_t classifyStarted = bu_gettime();
    buildCutSchedule();
    currCut = maxPopCut;
    triProcess();
    if (!buildClusters())
	return;
    if (!buildChunks())
	return;
    const int64_t cacheStarted = bu_gettime();

    isValid = true;
    cache();
    if (getenv("BOBOL_DRAW_TIMING")) {
	const int64_t finished = bu_gettime();
	bu_log("[obol-timing] pop cache: bounds scan %8.1f ms "
	       "classify %8.1f ms write %8.1f ms total %8.1f ms "
	       "(faces=%zu points=%zu)\n",
	       (classifyStarted - buildStarted) / 1000.0,
	       (cacheStarted - classifyStarted) / 1000.0,
	       (finished - cacheStarted) / 1000.0,
	       (finished - buildStarted) / 1000.0,
	       faceCount, vertexCount);
    }
}

BObolPopState::BObolPopState(struct BObolMeshLodContext *ctx,
				 unsigned long long key) :
    BObolPopState(ctx, key, false)
{
}

BObolPopState::BObolPopState(struct BObolMeshLodContext *ctx,
				 unsigned long long key,
				 bool retainHeaderSnapshot)
{
    context = ctx;
    if (!key)
	return;

    currCut = -1;
    hash = key;
    initializeCacheKeyPrefix();

    loadCachedHeader(retainHeaderSnapshot);
}

bool
BObolPopState::loadCachedHeader(bool retainHeaderSnapshot)
{
    isValid = false;

    /* These records form one immutable hierarchy header.  Read them through
     * one LMDB snapshot: besides cutting transaction churn, this prevents a
     * handle from accepting a mixture of records if an invalidation races
     * with construction. */
    const auto readMetadata = [this](const char *component, void *output,
	    size_t expected) -> bool {
	void *buffer = NULL;
	const size_t size = cacheGet(&buffer, component);
	if (size != expected || !buffer)
	    return false;
	memcpy(output, buffer, expected);
	return true;
    };
    const auto readCutCounts = [this](const char *component,
	size_t *output) -> bool {
	void *buffer = NULL;
	const size_t size = cacheGet(&buffer, component);
	if (!buffer || !output)
	    return false;
	const size_t count = POP_CUT_COUNT_MAX + 1;
	if (size == count * sizeof(uint64_t)) {
	    const unsigned char *bytes =
		static_cast<const unsigned char *>(buffer);
	    for (size_t cut = 0; cut < count; ++cut) {
		uint64_t diskCount = 0;
		memcpy(&diskCount, bytes + cut * sizeof(diskCount),
		    sizeof(diskCount));
		if (diskCount > static_cast<uint64_t>(SIZE_MAX))
		    return false;
		output[cut] = static_cast<size_t>(diskCount);
	    }
	    return true;
	}
	return false;
    };
    int32_t diskMaxCut = -1;
    int32_t diskMinCut = -1;
    uint8_t diskCutBits[POP_CUT_COUNT_MAX][3] = {{0}};
    uint8_t diskCullBackfaces = 0;
    uint8_t diskHasNormals = 0;
    if (!readMetadata(CACHE_POP_MAX_CUT, &diskMaxCut,
	    sizeof(diskMaxCut)) ||
	!readMetadata(CACHE_POP_MIN_CUT, &diskMinCut,
	    sizeof(diskMinCut)) ||
	!readMetadata(CACHE_CUT_QUANTIZATION, diskCutBits,
	    sizeof(diskCutBits)) ||
	!readMetadata(CACHE_SHADED_CULL_BACKFACES,
	    &diskCullBackfaces, sizeof(diskCullBackfaces)) ||
	!readMetadata(CACHE_HAS_NORMALS,
	    &diskHasNormals, sizeof(diskHasNormals)) ||
	!readCutCounts(CACHE_VERTEX_COUNT, cutVertexCount) ||
	!readCutCounts(CACHE_TRI_COUNT, cutTriangleCount)) {
	cacheDone();
	return false;
    }
    if (diskCullBackfaces > 1 || diskHasNormals > 1 ||
	diskMaxCut < 0 || diskMaxCut >= POP_CUT_COUNT_MAX ||
	diskMinCut < 0 || diskMinCut > diskMaxCut) {
	cacheDone();
	return false;
    }
    shadedCullBackfaces = diskCullBackfaces != 0;
    hasNormals = diskHasNormals != 0;

    {
	double diskBounds[12] = {};
	void *buffer = NULL;
	size_t bufferSize = cacheGet(&buffer, CACHE_OBJ_BOUNDS);
	if (bufferSize != sizeof(diskBounds) || !buffer) {
	    cacheDone();
	    return false;
	}
	memcpy(diskBounds, buffer, sizeof(diskBounds));
	VSET(bbmin, diskBounds[0], diskBounds[1], diskBounds[2]);
	VSET(bbmax, diskBounds[3], diskBounds[4], diskBounds[5]);
	minx = static_cast<fastf_t>(diskBounds[6]);
	miny = static_cast<fastf_t>(diskBounds[7]);
	minz = static_cast<fastf_t>(diskBounds[8]);
	maxx = static_cast<fastf_t>(diskBounds[9]);
	maxy = static_cast<fastf_t>(diskBounds[10]);
	maxz = static_cast<fastf_t>(diskBounds[11]);
    }
    {
	const bool validBounds =
	    std::isfinite(minx) && std::isfinite(miny) &&
	    std::isfinite(minz) && std::isfinite(maxx) &&
	    std::isfinite(maxy) && std::isfinite(maxz) &&
	    minx <= maxx && miny <= maxy && minz <= maxz &&
	    std::isfinite(bbmin[X]) && std::isfinite(bbmin[Y]) &&
	    std::isfinite(bbmin[Z]) && std::isfinite(bbmax[X]) &&
	    std::isfinite(bbmax[Y]) && std::isfinite(bbmax[Z]) &&
	    bbmin[X] <= bbmax[X] && bbmin[Y] <= bbmax[Y] &&
	    bbmin[Z] <= bbmax[Z] &&
	    std::isfinite(maxx - minx) &&
	    std::isfinite(maxy - miny) &&
	    std::isfinite(maxz - minz);
	if (!validBounds) {
	    cacheDone();
	    return false;
	}

	/* This format has one canonical anisotropic schedule.  Rebuild it from
	 * the validated quantization domain and require the stored cut records
	 * to match exactly.  This rejects truncated, stale, or corrupt metadata
	 * before any ordinal can index fixed-size storage. */
	buildCutSchedule();
	bool validCuts = maxPopCut == static_cast<int>(diskMaxCut) &&
	    cutCount == static_cast<uint32_t>(diskMaxCut + 1);
	for (uint32_t cut = 0; validCuts && cut < cutCount; ++cut)
	    for (int axis = 0; axis < 3; ++axis)
		if (cuts[cut].bits[axis] != diskCutBits[cut][axis])
		    validCuts = false;
	minPopCut = static_cast<int>(diskMinCut);

	size_t totalPoints = 0;
	size_t totalFaces = 0;
	bool validCounts = cutVertexCount[POP_CUT_COUNT_MAX] == 0 &&
	    cutTriangleCount[POP_CUT_COUNT_MAX] == 0;
	for (int cut = 0; validCounts && cut < POP_CUT_COUNT_MAX; ++cut) {
	    if (cutVertexCount[cut] > SIZE_MAX - totalPoints ||
		cutTriangleCount[cut] > SIZE_MAX - totalFaces) {
		validCounts = false;
		break;
	    }
	    totalPoints += cutVertexCount[cut];
	    totalFaces += cutTriangleCount[cut];
	    if (cut < minPopCut &&
		(cutVertexCount[cut] || cutTriangleCount[cut]))
		validCounts = false;
	    if (cut > maxPopCut &&
		(cutVertexCount[cut] || cutTriangleCount[cut]))
		validCounts = false;
	}
	if (!validCuts || !validCounts ||
	    !totalPoints || !totalFaces ||
	    !cutVertexCount[minPopCut] || !cutTriangleCount[minPopCut] ||
	    totalPoints > static_cast<size_t>(
		std::numeric_limits<int>::max()) ||
	    totalFaces > static_cast<size_t>(
		std::numeric_limits<int>::max())) {
	    cacheDone();
	    return false;
	}
    }
    {
	uint16_t diskGrid = 0;
	if (!readMetadata(CACHE_CLUSTER_GRID, &diskGrid,
		sizeof(diskGrid)) ||
	    diskGrid != BOBOL_MESH_LOD_CLUSTER_GRID_RESOLUTION) {
	    cacheDone();
	    return false;
	}

	void *boundsBytes = NULL;
	void *offsetBytes = NULL;
	void *rangeBytes = NULL;
	const size_t clusterCount = BOBOL_MESH_LOD_CLUSTER_COUNT;
	const size_t boundsSize = cacheGet(&boundsBytes,
	    CACHE_CLUSTER_BOUNDS);
	const size_t offsetsSize = cacheGet(&offsetBytes,
	    CACHE_CLUSTER_RANGE_OFFSETS);
	const size_t rangesSize = cacheGet(&rangeBytes,
	    CACHE_CLUSTER_RANGES);
	if (!boundsBytes || !offsetBytes || !rangeBytes ||
	    boundsSize != clusterCount * 6 * sizeof(double) ||
	    offsetsSize != (clusterCount + 1) * sizeof(uint32_t) ||
	    rangesSize == 0 ||
	    rangesSize % sizeof(BObolMeshLodClusterRangeDisk) != 0) {
	    cacheDone();
	    return false;
	}
	const size_t rangeCount =
	    rangesSize / sizeof(BObolMeshLodClusterRangeDisk);
	if (rangeCount > clusterCount * cutCount) {
	    cacheDone();
	    return false;
	}

	std::vector<uint32_t> offsets(clusterCount + 1);
	memcpy(offsets.data(), offsetBytes, offsetsSize);
	if (offsets[0] != 0 || offsets.back() != rangeCount) {
	    cacheDone();
	    return false;
	}
	for (size_t cluster = 0; cluster < clusterCount; ++cluster) {
	    if (offsets[cluster] > offsets[cluster + 1]) {
		cacheDone();
		return false;
	    }
	}

	size_t totalFaces = 0;
	std::array<uint64_t, POP_CUT_COUNT_MAX + 1> cutEnds = {};
	for (uint32_t cut = 0; cut < cutCount; ++cut) {
	    totalFaces += cutTriangleCount[cut];
	    cutEnds[cut + 1] = static_cast<uint64_t>(totalFaces) * 3;
	}
	const uint64_t totalIndices = static_cast<uint64_t>(totalFaces) * 3;
	clusterRanges.resize(rangeCount);
	const unsigned char *diskRanges =
	    static_cast<const unsigned char *>(rangeBytes);
	for (size_t range = 0; range < rangeCount; ++range) {
	    BObolMeshLodClusterRangeDisk disk = {};
	    memcpy(&disk,
		diskRanges + range * sizeof(disk), sizeof(disk));
	    const uint64_t end = static_cast<uint64_t>(disk.firstIndex) +
		disk.indexCount;
	    if (disk.reserved[0] || disk.reserved[1] || disk.reserved[2] ||
		disk.activationCut >= cutCount ||
		disk.firstIndex % 3 || disk.indexCount < 3 ||
		disk.indexCount % 3 || end > totalIndices ||
		disk.firstIndex < cutEnds[disk.activationCut] ||
		end > cutEnds[disk.activationCut + 1]) {
		cacheDone();
		return false;
	    }
	    clusterRanges[range] = {
		disk.firstIndex, disk.indexCount, disk.activationCut
	    };
	}

	clusterInfos.resize(clusterCount);
	const unsigned char *diskBounds =
	    static_cast<const unsigned char *>(boundsBytes);
	for (size_t cluster = 0; cluster < clusterCount; ++cluster) {
	    double values[6] = {};
	    memcpy(values, diskBounds + cluster * sizeof(values),
		sizeof(values));
	    BObolMeshLodClusterInfo &info = clusterInfos[cluster];
	    info.range_count = offsets[cluster + 1] - offsets[cluster];
	    info.ranges = info.range_count ?
		clusterRanges.data() + offsets[cluster] : NULL;
	    if (!info.range_count) {
		VSETALL(info.bmin, 0.0);
		VSETALL(info.bmax, 0.0);
		continue;
	    }
	    bool valid = true;
	    for (double value : values)
		valid = valid && std::isfinite(value);
	    valid = valid && values[0] <= values[3] &&
		values[1] <= values[4] && values[2] <= values[5] &&
		values[0] >= minx && values[1] >= miny &&
		values[2] >= minz && values[3] <= maxx &&
		values[4] <= maxy && values[5] <= maxz;
	    if (!valid) {
		cacheDone();
		return false;
	    }
	    VSET(info.bmin, values[0], values[1], values[2]);
	    VSET(info.bmax, values[3], values[4], values[5]);
	}
    }
    {
	const auto chunkHeaderFailure = [this](const char *reason) {
	    if (getenv("BOBOL_CHUNK_DEBUG"))
		bu_log("BObol chunk header validation failed: %s\n", reason);
	    cacheDone();
	    return false;
	};
	uint32_t diskChunkCount = 0;
	if (!readMetadata(CACHE_CHUNK_COUNT, &diskChunkCount,
		sizeof(diskChunkCount)) || !diskChunkCount ||
	    diskChunkCount > static_cast<uint32_t>(
		std::numeric_limits<int>::max())) {
	    return chunkHeaderFailure("count");
	}
	void *boundsBytes = NULL;
	void *minmaxBytes = NULL;
	void *faceBytes = NULL;
	void *pointBytes = NULL;
	void *residentBytes = NULL;
	const size_t chunkCount = diskChunkCount;
	if (chunkCount > SIZE_MAX / (POP_CUT_COUNT_MAX * sizeof(uint64_t)) ||
	    cacheGet(&boundsBytes, CACHE_CHUNK_BOUNDS) !=
		chunkCount * 6 * sizeof(double) ||
	    cacheGet(&minmaxBytes, CACHE_CHUNK_MINMAX) !=
		chunkCount * 2 * sizeof(uint8_t) ||
	    cacheGet(&faceBytes, CACHE_CHUNK_FACE_COUNTS) !=
		chunkCount * POP_CUT_COUNT_MAX * sizeof(uint32_t) ||
	    cacheGet(&pointBytes, CACHE_CHUNK_POINT_COUNTS) !=
		chunkCount * POP_CUT_COUNT_MAX * sizeof(uint32_t) ||
	    cacheGet(&residentBytes, CACHE_CHUNK_RESIDENT_BYTES) !=
		chunkCount * POP_CUT_COUNT_MAX * sizeof(uint64_t) ||
	    !boundsBytes || !minmaxBytes || !faceBytes || !pointBytes ||
	    !residentBytes) {
	    return chunkHeaderFailure("arrays");
	}
	chunkInfos.resize(chunkCount);
	std::array<uint64_t, POP_CUT_COUNT_MAX> summedFaces = {};
	std::array<uint64_t, POP_CUT_COUNT_MAX> summedPoints = {};
	const unsigned char *diskBounds =
	    static_cast<const unsigned char *>(boundsBytes);
	const uint8_t *diskMinmax = static_cast<const uint8_t *>(minmaxBytes);
	const unsigned char *diskFaces =
	    static_cast<const unsigned char *>(faceBytes);
	const unsigned char *diskPoints =
	    static_cast<const unsigned char *>(pointBytes);
	const unsigned char *diskResident =
	    static_cast<const unsigned char *>(residentBytes);
	for (size_t chunk = 0; chunk < chunkCount; ++chunk) {
	    BObolMeshLodChunkInfo &info = chunkInfos[chunk];
	    memset(&info, 0, sizeof(info));
	    info.chunk_id = static_cast<uint32_t>(chunk);
	    info.min_cut = diskMinmax[chunk * 2 + 0];
	    info.max_cut = diskMinmax[chunk * 2 + 1];
	    double values[6] = {};
	    memcpy(values, diskBounds + chunk * sizeof(values),
		sizeof(values));
	    bool valid = info.min_cut >= minPopCut &&
		info.min_cut <= info.max_cut && info.max_cut == maxPopCut;
	    for (double value : values)
		valid = valid && std::isfinite(value);
	    valid = valid && values[0] <= values[3] &&
		values[1] <= values[4] && values[2] <= values[5] &&
		values[0] >= minx && values[1] >= miny &&
		values[2] >= minz && values[3] <= maxx &&
		values[4] <= maxy && values[5] <= maxz;
	    if (!valid) {
		return chunkHeaderFailure("bounds");
	    }
	    VSET(info.bmin, values[0], values[1], values[2]);
	    VSET(info.bmax, values[3], values[4], values[5]);
	    uint32_t priorFaces = 0;
	    uint32_t priorPoints = 0;
	    uint64_t priorResident = 0;
	    for (uint32_t cut = 0; cut < POP_CUT_COUNT_MAX; ++cut) {
		const size_t offset = chunk * POP_CUT_COUNT_MAX + cut;
		uint32_t faces = 0;
		uint32_t points = 0;
		uint64_t resident = 0;
		memcpy(&faces, diskFaces + offset * sizeof(faces),
		    sizeof(faces));
		memcpy(&points, diskPoints + offset * sizeof(points),
		    sizeof(points));
		memcpy(&resident, diskResident + offset * sizeof(resident),
		    sizeof(resident));
		if ((cut <= static_cast<uint32_t>(info.max_cut) &&
		     (faces < priorFaces || points < priorPoints ||
		      resident < priorResident)) ||
		    (cut < static_cast<uint32_t>(info.min_cut) &&
		     (faces || points || resident)) ||
		    (cut > static_cast<uint32_t>(info.max_cut) &&
		     (faces || points || resident))) {
		    if (getenv("BOBOL_CHUNK_DEBUG"))
			bu_log("chunk=%zu cut=%u min=%d max=%d "
			       "faces=%u/%u points=%u/%u bytes=%llu/%llu\n",
			       chunk, cut, info.min_cut, info.max_cut,
			       faces, priorFaces, points, priorPoints,
			       static_cast<unsigned long long>(resident),
			       static_cast<unsigned long long>(priorResident));
		    return chunkHeaderFailure("monotonic counts");
		}
		info.cuts[cut].face_count = faces;
		info.cuts[cut].point_count = points;
		info.cuts[cut].resident_bytes = resident;
		if (cut <= static_cast<uint32_t>(info.max_cut)) {
		    priorFaces = faces;
		    priorPoints = points;
		    priorResident = resident;
		}
		if (cut < cutCount) {
		    summedFaces[cut] += faces;
		    summedPoints[cut] += points;
		}
	    }
	    if (!priorFaces || !priorPoints) {
		return chunkHeaderFailure("empty terminal");
	    }
	}
	uint64_t globalFaces = 0;
	uint64_t globalPoints = 0;
	for (uint32_t cut = 0; cut < cutCount; ++cut) {
	    globalFaces += cutTriangleCount[cut];
	    globalPoints += cutVertexCount[cut];
	    /* Vertices shared by chunks are duplicated by design, so only the
	     * face population is conserved exactly. */
	    if (summedFaces[cut] != globalFaces ||
		summedPoints[cut] < globalPoints) {
		return chunkHeaderFailure("global population sums");
	    }
	}
    }
    /*
     * Retained drawing loads its first cumulative prefix immediately.  Let
     * that explicit path reuse this immutable header snapshot; ordinary
     * handles close it here so an idle API consumer never pins an LMDB reader
     * indefinitely.
     */
    if (!retainHeaderSnapshot)
	cacheDone();

    /* Metadata validation must not also read geometry.  The caller supplies
     * the view-selected target immediately after opening the handle, so
     * eagerly loading the minimum cut here doubled cache I/O for every
     * request above that cut. */
    isValid = true;
    return true;
}

BObolPopState::~BObolPopState()
{
    if (readTxn || readLock.owns_lock())
	cacheDone();
}

bool
BObolPopState::materializeInitialPrefix(int cut)
{
    if (!isValid)
	return false;
    /* Generate-and-open may discover an immutable hierarchy already produced
     * for identical content under another source name.  Reuse its cache;
     * only a genuinely new hierarchy has producer scratch for direct
     * materialization. */
    if (!vertexArray || !faceArray || cutTriVerts.empty() ||
	cutTris.empty() || triIndexMap.size() != vertexCount)
	return setCut(cut);
    cut = std::max(minPopCut, std::min(maxPopCut, cut));

    size_t pointCount = 0;
    size_t triangleCount = 0;
    for (int cutIndex = 0; cutIndex <= cut; ++cutIndex) {
	pointCount += cutVertexCount[cutIndex];
	triangleCount += cutTriangleCount[cutIndex];
    }
    if (!pointCount || !triangleCount ||
	pointCount > SIZE_MAX / 3 || triangleCount > SIZE_MAX / 3)
	return false;

    lodTriPoints.clear();
    lodTriPoints.reserve(pointCount * 3);
    for (int cutIndex = 0; cutIndex <= cut; ++cutIndex) {
	for (uint32_t sourceVertex : cutTriVerts[cutIndex]) {
	    if (sourceVertex >= vertexCount)
		return false;
	    lodTriPoints.push_back(vertexArray[sourceVertex][X]);
	    lodTriPoints.push_back(vertexArray[sourceVertex][Y]);
	    lodTriPoints.push_back(vertexArray[sourceVertex][Z]);
	}
    }

    lodTris.clear();
    lodTris.reserve(triangleCount * 3);
    if (normalArray) {
	if (triangleCount > SIZE_MAX / 9)
	    return false;
	lodTriNormals.clear();
	lodTriNormals.reserve(triangleCount * 9);
    }
    for (int cutIndex = 0; cutIndex <= cut; ++cutIndex) {
	for (uint32_t sourceFace : cutTris[cutIndex]) {
	    if (sourceFace >= faceCount)
		return false;
	    for (size_t corner = 0; corner < 3; ++corner) {
		const int sourceVertex =
		    faceArray[3 * static_cast<size_t>(sourceFace) + corner];
		if (sourceVertex < 0 ||
		    static_cast<size_t>(sourceVertex) >=
			triIndexMap.size())
		    return false;
		lodTris.push_back(triIndexMap[sourceVertex]);
		if (normalArray) {
		    const size_t normalIndex =
			3 * static_cast<size_t>(sourceFace) + corner;
		    lodTriNormals.push_back(normalArray[normalIndex][X]);
		    lodTriNormals.push_back(normalArray[normalIndex][Y]);
		    lodTriNormals.push_back(normalArray[normalIndex][Z]);
		}
	    }
	}
    }
    lodTriPointsSnapped.clear();
    currCut = cut;
    return true;
}

void
BObolPopState::releaseGenerationScratch(void)
{
    triIndexMap.clear();
    triIndexMap.shrink_to_fit();
    vertexTriMinCut.clear();
    vertexTriMinCut.shrink_to_fit();
    faceActivationCut.clear();
    faceActivationCut.shrink_to_fit();
    cutTriVerts.clear();
    cutTriVerts.shrink_to_fit();
    cutTris.clear();
    cutTris.shrink_to_fit();
    chunkFaces.clear();
    chunkFaces.shrink_to_fit();
    vertexCount = 0;
    vertexArray = NULL;
    normalArray = NULL;
    faceCount = 0;
    faceArray = NULL;
}

void
BObolPopState::updateSnappedPoints(int cut)
{
    lodTriPointsSnapped.clear();
    lodTriPointsSnapped.reserve(lodTriPoints.size());
    for (size_t pointIndex = 0; pointIndex < lodTriPoints.size() / 3;
	 pointIndex++) {
	point_t point;
	point_t snapped;
	VSET(point, lodTriPoints[3 * pointIndex + 0],
	     lodTriPoints[3 * pointIndex + 1],
	     lodTriPoints[3 * pointIndex + 2]);
	cutPoint(&snapped, &point, cut);
	for (int coordIndex = 0; coordIndex < 3; coordIndex++)
	    lodTriPointsSnapped.push_back(snapped[coordIndex]);
    }
}

bool
BObolPopState::triPopLoad(int startCut, int cut,
			  bool materializeSnapped)
{
    char component[8] = {0};
    const auto cutComponent = [&component](const char *prefix,
	    int cutIndex) -> const char * {
	char *out = component;
	const size_t prefixLength = strlen(prefix);
	memcpy(out, prefix, prefixLength);
	out += prefixLength;
	const std::to_chars_result converted =
	    std::to_chars(out, component + sizeof(component) - 1, cutIndex);
	if (converted.ec != std::errc()) {
	    component[0] = '\0';
	    return component;
	}
	*converted.ptr = '\0';
	return component;
    };

    /* A requested PoP prefix is one immutable cache snapshot.  Keep the LMDB
     * read transaction open while copying every vertex, triangle, and normal
     * chunk in the prefix.  Finishing a transaction after every component
     * made a 50k-asset warm scene perform hundreds of thousands of tiny read
     * transactions and reduced an otherwise memory-resident cache to roughly
     * a hundred meshes per second. */
    for (int cutIndex = startCut + 1; cutIndex <= cut; ++cutIndex) {
	if (!cutVertexCount[cutIndex])
	    continue;
	void *buffer = NULL;
	const size_t bufferSize = cacheGet(&buffer,
	    cutComponent(CACHE_VERT_CUT, cutIndex));
	const size_t scalarCount = cutVertexCount[cutIndex] * 3;
	if (bufferSize != scalarCount * sizeof(fastf_t) || !buffer ||
	    scalarCount > lodTriPoints.max_size() - lodTriPoints.size()) {
	    cacheDone();
	    return false;
	}
	const size_t offset = lodTriPoints.size();
	lodTriPoints.resize(offset + scalarCount);
	memcpy(lodTriPoints.data() + offset, buffer, bufferSize);
	const fastf_t minimum[3] = {minx, miny, minz};
	const fastf_t maximum[3] = {maxx, maxy, maxz};
	for (size_t scalar = offset; scalar < offset + scalarCount; ++scalar) {
	    const int axis = static_cast<int>(scalar % 3);
	    if (!std::isfinite(lodTriPoints[scalar]) ||
		lodTriPoints[scalar] < minimum[axis] ||
		lodTriPoints[scalar] > maximum[axis]) {
		cacheDone();
		return false;
	    }
	}
    }

    if (materializeSnapped)
	updateSnappedPoints(cut);
    else
	lodTriPointsSnapped.clear();

    for (int cutIndex = startCut + 1; cutIndex <= cut; ++cutIndex) {
	if (!cutTriangleCount[cutIndex])
	    continue;
	void *buffer = NULL;
	const size_t bufferSize = cacheGet(&buffer,
	    cutComponent(CACHE_TRI_CUT, cutIndex));
	const size_t indexCount = cutTriangleCount[cutIndex] * 3;
	if (bufferSize != indexCount * sizeof(uint32_t) || !buffer ||
	    indexCount > lodTris.max_size() - lodTris.size()) {
	    cacheDone();
	    return false;
	}
	const size_t offset = lodTris.size();
	lodTris.resize(offset + indexCount);
	memcpy(lodTris.data() + offset, buffer, bufferSize);
	const size_t cumulativePointCount = lodTriPoints.size() / 3;
	for (size_t index = offset; index < offset + indexCount; ++index) {
	    if (static_cast<size_t>(lodTris[index]) >= cumulativePointCount) {
		cacheDone();
		return false;
	    }
	}
    }

    for (int cutIndex = startCut + 1; cutIndex <= cut; ++cutIndex) {
	if (!cutTriangleCount[cutIndex])
	    continue;
	void *buffer = NULL;
	const size_t bufferSize = cacheGet(&buffer,
	    cutComponent(CACHE_VERTNORM_CUT, cutIndex));
	const size_t scalarCount = cutTriangleCount[cutIndex] * 9;
	const size_t expectedSize = hasNormals ?
	    scalarCount * sizeof(fastf_t) : 0;
	if (bufferSize != expectedSize ||
	    (expectedSize && (!buffer || scalarCount >
		lodTriNormals.max_size() - lodTriNormals.size()))) {
	    cacheDone();
	    return false;
	}
	if (bufferSize) {
	    const size_t offset = lodTriNormals.size();
	    lodTriNormals.resize(offset + scalarCount);
	    memcpy(lodTriNormals.data() + offset, buffer, bufferSize);
	    for (size_t scalar = offset; scalar < offset + scalarCount; ++scalar) {
		if (!std::isfinite(lodTriNormals[scalar])) {
		    cacheDone();
		    return false;
		}
	    }
	}
    }

    cacheDone();
    return true;
}

bool
BObolPopState::readSuffix(int residentCut, int targetCut,
	BObolMeshLodSuffixCallback callback, void *callbackData)
{
    if (!isValid || !callback || residentCut < minPopCut ||
	targetCut <= residentCut || targetCut > maxPopCut)
	return false;

    char component[16] = {0};
    const auto cutComponent = [&component](const char *prefix,
	int cutIndex) -> const char * {
	const size_t prefixLength = strlen(prefix);
	if (prefixLength >= sizeof(component)) {
	    component[0] = '\0';
	    return component;
	}
	memcpy(component, prefix, prefixLength);
	char *out = component + prefixLength;
	const std::to_chars_result converted =
	    std::to_chars(out, component + sizeof(component) - 1, cutIndex);
	if (converted.ec != std::errc()) {
	    component[0] = '\0';
	    return component;
	}
	*converted.ptr = '\0';
	return component;
    };

    /* LMDB values are byte records and need not satisfy point_t, vect_t, or
     * uint32_t alignment.  Copy one cut at a time into reusable aligned
     * scratch before exposing typed callback pointers.  This keeps transient
     * memory bounded by the largest activation chunk rather than the whole
     * cumulative prefix. */
    std::vector<fastf_t> pointStorage;
    std::vector<uint32_t> faceStorage;
    std::vector<fastf_t> normalStorage;
    const fastf_t minimum[3] = {minx, miny, minz};
    const fastf_t maximum[3] = {maxx, maxy, maxz};
    for (int cutIndex = residentCut + 1;
	 cutIndex <= targetCut; ++cutIndex) {
	const size_t pointCount = cutVertexCount[cutIndex];
	const size_t chunkFaceCount = cutTriangleCount[cutIndex];
	if (pointCount > SIZE_MAX / sizeof(point_t) ||
	    chunkFaceCount > SIZE_MAX / 3 ||
	    chunkFaceCount * 3 > SIZE_MAX / sizeof(uint32_t) ||
	    (hasNormals &&
	     chunkFaceCount * 3 > SIZE_MAX / sizeof(vect_t))) {
	    cacheDone();
	    return false;
	}
	pointStorage.resize(pointCount * 3);
	faceStorage.resize(chunkFaceCount * 3);
	normalStorage.resize(hasNormals ? chunkFaceCount * 9 : 0);
	if (pointCount) {
	    void *bytesData = NULL;
	    const size_t bytes = cacheGet(&bytesData,
		cutComponent(CACHE_VERT_CUT, cutIndex));
	    if (bytes != pointStorage.size() * sizeof(fastf_t) ||
		!bytesData) {
		cacheDone();
		return false;
	    }
	    memcpy(pointStorage.data(), bytesData, bytes);
	    for (size_t point = 0; point < pointCount; ++point) {
		for (int axis = 0; axis < 3; ++axis) {
		    const fastf_t value = pointStorage[point * 3 + axis];
		    if (!std::isfinite(value) || value < minimum[axis] ||
			value > maximum[axis]) {
			cacheDone();
			return false;
		    }
		}
	    }
	}
	if (chunkFaceCount) {
	    void *bytesData = NULL;
	    size_t bytes = cacheGet(&bytesData,
		cutComponent(CACHE_TRI_CUT, cutIndex));
	    if (bytes != faceStorage.size() * sizeof(uint32_t) ||
		!bytesData) {
		cacheDone();
		return false;
	    }
	    memcpy(faceStorage.data(), bytesData, bytes);
	    const size_t finalPointCount =
		static_cast<size_t>(cutVertexCount[0]);
	    size_t cumulativePointCount = finalPointCount;
	    for (int priorCut = 1; priorCut <= cutIndex; ++priorCut)
		cumulativePointCount += cutVertexCount[priorCut];
	    for (uint32_t index : faceStorage) {
		if (static_cast<size_t>(index) >= cumulativePointCount) {
		    cacheDone();
		    return false;
		}
	    }
	    if (hasNormals) {
		bytesData = NULL;
		bytes = cacheGet(&bytesData,
		    cutComponent(CACHE_VERTNORM_CUT, cutIndex));
		if (bytes != normalStorage.size() * sizeof(fastf_t) ||
		    !bytesData) {
		    cacheDone();
		    return false;
		}
		memcpy(normalStorage.data(), bytesData, bytes);
		for (fastf_t value : normalStorage) {
		    if (!std::isfinite(value)) {
			cacheDone();
			return false;
		    }
		}
	    }
	}
	const point_t *points = pointStorage.empty() ? NULL :
	    reinterpret_cast<const point_t *>(pointStorage.data());
	const uint32_t *faces = faceStorage.empty() ? NULL :
	    faceStorage.data();
	const vect_t *normals = normalStorage.empty() ? NULL :
	    reinterpret_cast<const vect_t *>(normalStorage.data());
	if (!callback(cutIndex, points, pointCount, faces, chunkFaceCount,
		normals, normals ? chunkFaceCount * 3 : 0, callbackData)) {
	    cacheDone();
	    return false;
	}
    }
    cacheDone();
    return true;
}

bool
BObolPopState::readChunks(const uint32_t *chunkIds, size_t requestedCount,
	int cut, BObolMeshLodChunkCallback callback, void *callbackData)
{
    if (!isValid || !chunkIds || !requestedCount || !callback ||
	cut < minPopCut || cut > maxPopCut || chunkInfos.empty())
	return false;
    for (size_t i = 0; i < requestedCount; ++i) {
	if (chunkIds[i] >= chunkInfos.size() ||
	    (i && chunkIds[i] <= chunkIds[i - 1]))
	    return false;
    }

    std::vector<fastf_t> pointStorage;
    std::vector<uint32_t> faceStorage;
    std::vector<fastf_t> normalStorage;
    for (size_t requested = 0; requested < requestedCount; ++requested) {
	const uint32_t chunkId = chunkIds[requested];
	const BObolMeshLodChunkInfo &info = chunkInfos[chunkId];
	if (cut < info.min_cut) {
	    if (!callback(chunkId, cut, NULL, 0, NULL, 0, NULL, 0,
		    callbackData)) {
		cacheDone();
		return false;
	    }
	    continue;
	}
	const BObolMeshLodChunkCutInfo &selected = info.cuts[cut];
	const BObolMeshLodChunkCutInfo &terminal = info.cuts[info.max_cut];
	if (!selected.point_count || !selected.face_count ||
	    selected.point_count > terminal.point_count ||
	    selected.face_count > terminal.face_count)
	    return false;

	char component[32] = {0};
	snprintf(component, sizeof(component), "%s%08x",
	    CACHE_CHUNK_DATA_PREFIX, chunkId);
	void *recordData = NULL;
	const size_t recordSize = cacheGet(&recordData, component);
	if (!recordData || recordSize < 88u +
		POP_CUT_COUNT_MAX * 2u * sizeof(uint32_t)) {
	    cacheDone();
	    return false;
	}
	const unsigned char *record =
	    static_cast<const unsigned char *>(recordData);
	size_t offset = 0;
	const auto read = [&record, recordSize, &offset](void *value,
		size_t size) -> bool {
	    if (!value || size > recordSize - std::min(recordSize, offset))
		return false;
	    memcpy(value, record + offset, size);
	    offset += size;
	    return true;
	};
	uint8_t magic[8] = {};
	uint32_t version = 0;
	uint32_t diskChunk = UINT32_MAX;
	uint32_t diskCutCount = 0;
	uint32_t flags = 0;
	uint32_t fullPoints = 0;
	uint32_t fullFaces = 0;
	uint32_t reservedA = 1;
	uint32_t reservedB = 1;
	double bounds[6] = {};
	if (!read(magic, sizeof(magic)) ||
	    !read(&version, sizeof(version)) ||
	    !read(&diskChunk, sizeof(diskChunk)) ||
	    !read(&diskCutCount, sizeof(diskCutCount)) ||
	    !read(&flags, sizeof(flags)) ||
	    !read(&fullPoints, sizeof(fullPoints)) ||
	    !read(&fullFaces, sizeof(fullFaces)) ||
	    !read(&reservedA, sizeof(reservedA)) ||
	    !read(&reservedB, sizeof(reservedB)) ||
	    !read(bounds, sizeof(bounds)) ||
	    memcmp(magic, meshLodChunkMagic, sizeof(magic)) != 0 ||
	    version != 1 || diskChunk != chunkId ||
	    diskCutCount != cutCount || flags > 1 ||
	    (flags != 0) != hasNormals || reservedA || reservedB ||
	    fullPoints != terminal.point_count ||
	    fullFaces != terminal.face_count)
	    goto chunk_read_failed;
	{
	    uint32_t diskPoints[POP_CUT_COUNT_MAX] = {};
	    uint32_t diskFaces[POP_CUT_COUNT_MAX] = {};
	    if (!read(diskPoints, sizeof(diskPoints)) ||
		!read(diskFaces, sizeof(diskFaces)))
		goto chunk_read_failed;
	    for (uint32_t c = 0; c < POP_CUT_COUNT_MAX; ++c) {
		if (diskPoints[c] != info.cuts[c].point_count ||
		    diskFaces[c] != info.cuts[c].face_count)
		    goto chunk_read_failed;
	    }
	}
	{
	    const uint64_t pointBytes =
		static_cast<uint64_t>(fullPoints) * 3u * sizeof(double);
	    const uint64_t indexBytes =
		static_cast<uint64_t>(fullFaces) * 3u * sizeof(uint32_t);
	    const uint64_t normalBytes = hasNormals ?
		static_cast<uint64_t>(fullFaces) * 9u * sizeof(double) : 0;
	    const uint64_t expected = static_cast<uint64_t>(offset) +
		pointBytes + indexBytes + normalBytes;
	    if (expected != recordSize || expected > SIZE_MAX)
		goto chunk_read_failed;
	    const size_t pointOffset = offset;
	    const size_t indexOffset = pointOffset +
		static_cast<size_t>(pointBytes);
	    const size_t normalOffset = indexOffset +
		static_cast<size_t>(indexBytes);
	    pointStorage.resize(static_cast<size_t>(selected.point_count) * 3);
	    faceStorage.resize(static_cast<size_t>(selected.face_count) * 3);
	    normalStorage.resize(hasNormals ?
		static_cast<size_t>(selected.face_count) * 9 : 0);
	    for (size_t scalar = 0; scalar < pointStorage.size(); ++scalar) {
		double value = 0.0;
		memcpy(&value, record + pointOffset + scalar * sizeof(value),
		    sizeof(value));
		const int axis = static_cast<int>(scalar % 3);
		if (!std::isfinite(value) || value < info.bmin[axis] ||
		    value > info.bmax[axis])
		    goto chunk_read_failed;
		pointStorage[scalar] = static_cast<fastf_t>(value);
	    }
	    memcpy(faceStorage.data(), record + indexOffset,
		faceStorage.size() * sizeof(faceStorage[0]));
	    for (uint32_t index : faceStorage)
		if (index >= selected.point_count)
		    goto chunk_read_failed;
	    for (size_t scalar = 0; scalar < normalStorage.size(); ++scalar) {
		double value = 0.0;
		memcpy(&value, record + normalOffset +
		    scalar * sizeof(value), sizeof(value));
		if (!std::isfinite(value))
		    goto chunk_read_failed;
		normalStorage[scalar] = static_cast<fastf_t>(value);
	    }
	}
	if (!callback(chunkId, cut,
		reinterpret_cast<const point_t *>(pointStorage.data()),
		selected.point_count, faceStorage.data(), selected.face_count,
		normalStorage.empty() ? NULL :
		    reinterpret_cast<const vect_t *>(normalStorage.data()),
		normalStorage.empty() ? 0 :
		    static_cast<size_t>(selected.face_count) * 3,
		callbackData)) {
	    cacheDone();
	    return false;
	}
	continue;

chunk_read_failed:
	cacheDone();
	return false;
    }
    cacheDone();
    return true;
}

void
BObolPopState::shrinkMemory(void)
{
    lodTriPoints.clear();
    lodTriPoints.shrink_to_fit();
    lodTriNormals.clear();
    lodTriNormals.shrink_to_fit();
    lodTris.clear();
    lodTris.shrink_to_fit();
    lodTriPointsSnapped.clear();
    lodTriPointsSnapped.shrink_to_fit();
}

static size_t
mesh_lod_saturating_bytes_add(size_t current, size_t count, size_t stride)
{
    if (count && stride > SIZE_MAX / count)
	return SIZE_MAX;
    const size_t amount = count * stride;
    return amount > SIZE_MAX - current ? SIZE_MAX : current + amount;
}

size_t
BObolPopState::residentPrefixBytes(void) const
{
    size_t bytes = 0;
    bytes = mesh_lod_saturating_bytes_add(
	bytes, lodTris.capacity(), sizeof(lodTris[0]));
    bytes = mesh_lod_saturating_bytes_add(
	bytes, lodTriPoints.capacity(), sizeof(lodTriPoints[0]));
    bytes = mesh_lod_saturating_bytes_add(
	bytes, lodTriPointsSnapped.capacity(),
	sizeof(lodTriPointsSnapped[0]));
    bytes = mesh_lod_saturating_bytes_add(
	bytes, lodTriNormals.capacity(), sizeof(lodTriNormals[0]));
    return bytes;
}

size_t
BObolPopState::residentBytes(void) const
{
    size_t bytes = sizeof(*this);
    bytes = mesh_lod_saturating_bytes_add(
	bytes, triIndexMap.capacity(), sizeof(triIndexMap[0]));
    bytes = mesh_lod_saturating_bytes_add(
	bytes, vertexTriMinCut.capacity(), sizeof(vertexTriMinCut[0]));
    bytes = mesh_lod_saturating_bytes_add(
	bytes, faceActivationCut.capacity(), sizeof(faceActivationCut[0]));
    bytes = mesh_lod_saturating_bytes_add(
	bytes, cutTriVerts.capacity(), sizeof(cutTriVerts[0]));
    bytes = mesh_lod_saturating_bytes_add(
	bytes, cutTris.capacity(), sizeof(cutTris[0]));
    bytes = mesh_lod_saturating_bytes_add(
	bytes, clusterInfos.capacity(), sizeof(clusterInfos[0]));
    bytes = mesh_lod_saturating_bytes_add(
	bytes, clusterRanges.capacity(), sizeof(clusterRanges[0]));
    bytes = mesh_lod_saturating_bytes_add(
	bytes, chunkInfos.capacity(), sizeof(chunkInfos[0]));
    bytes = mesh_lod_saturating_bytes_add(
	bytes, chunkFaces.capacity(), sizeof(chunkFaces[0]));
    for (const std::vector<uint32_t> &cut : cutTriVerts)
	bytes = mesh_lod_saturating_bytes_add(
	    bytes, cut.capacity(), sizeof(cut[0]));
    for (const std::vector<uint32_t> &cut : cutTris)
	bytes = mesh_lod_saturating_bytes_add(
	    bytes, cut.capacity(), sizeof(cut[0]));
    for (const std::vector<uint32_t> &chunk : chunkFaces)
	bytes = mesh_lod_saturating_bytes_add(
	    bytes, chunk.capacity(), sizeof(chunk[0]));
    const size_t prefix = residentPrefixBytes();
    return prefix > SIZE_MAX - bytes ? SIZE_MAX : bytes + prefix;
}

void
BObolPopState::triPopTrim(int cut, bool materializeSnapped)
{
    size_t keepVertexCount = 0;
    size_t keepFaceCount = 0;
    for (size_t cutIndex = 0; cutIndex <= static_cast<size_t>(cut);
	 ++cutIndex) {
	keepVertexCount += cutVertexCount[cutIndex];
	keepFaceCount += cutTriangleCount[cutIndex];
    }

    lodTriPoints.resize(keepVertexCount * 3);
    lodTriPoints.shrink_to_fit();
    lodTriNormals.resize(keepFaceCount * 3 * 3);
    lodTriNormals.shrink_to_fit();
    lodTris.resize(keepFaceCount * 3);
    lodTris.shrink_to_fit();

    if (materializeSnapped)
	updateSnappedPoints(cut);
    else
	lodTriPointsSnapped.clear();
}

void
BObolPopState::hierarchyInfo(struct BObolMeshLodHierarchyInfo *info) const
{
    if (!info)
	return;
    info->min_cut = minPopCut;
    info->max_cut = maxPopCut;
    info->resident_cut = currCut;
    info->cut_count = cutCount;
    info->has_normals = hasNormals ? 1 : 0;
    info->shaded_cull_backfaces = shadedCullBackfaces ? 1 : 0;
    VSET(info->quantization_min, minx, miny, minz);
    VSET(info->quantization_max, maxx, maxy, maxz);
    info->cluster_grid_resolution = clusterInfos.empty() ? 0 :
	BOBOL_MESH_LOD_CLUSTER_GRID_RESOLUTION;
    info->cluster_count = static_cast<uint32_t>(clusterInfos.size());
    info->clusters = clusterInfos.empty() ? NULL : clusterInfos.data();
    info->chunk_count = static_cast<uint32_t>(chunkInfos.size());
    info->chunks = chunkInfos.empty() ? NULL : chunkInfos.data();
    uint64_t points = 0;
    uint64_t faces = 0;
    for (uint32_t cut = 0; cut < cutCount; ++cut) {
	points += static_cast<uint64_t>(cutVertexCount[cut]);
	faces += static_cast<uint64_t>(cutTriangleCount[cut]);
	struct BObolMeshLodCutInfo &cutInfo = info->cuts[cut];
	cutInfo.point_count = points;
	cutInfo.face_count = faces;
	uint64_t bytes = points > UINT64_MAX / sizeof(point_t) ? UINT64_MAX :
	    points * sizeof(point_t);
	const uint64_t indexBytes = faces >
	    UINT64_MAX / (3 * sizeof(uint32_t)) ? UINT64_MAX :
	    faces * 3 * sizeof(uint32_t);
	bytes = indexBytes > UINT64_MAX - bytes ? UINT64_MAX :
	    bytes + indexBytes;
	if (hasNormals) {
	    const uint64_t normalBytes = faces >
		UINT64_MAX / (3 * sizeof(vect_t)) ? UINT64_MAX :
		faces * 3 * sizeof(vect_t);
	    bytes = normalBytes > UINT64_MAX - bytes ? UINT64_MAX :
		bytes + normalBytes;
	}
	cutInfo.resident_bytes = bytes;
	cutInfo.object_error = cuts[cut].objectError;
	for (int axis = 0; axis < 3; ++axis)
	    cutInfo.quantization_bits[axis] = cuts[cut].bits[axis];
	cutInfo.exact = cut == static_cast<uint32_t>(maxPopCut) ? 1 : 0;
    }
}

bool
BObolPopState::setCut(int cut, bool materializeSnapped)
{
    cut = std::max(minPopCut, std::min(maxPopCut, cut));

    if (cut == currCut && !forceUpdate) {
	if (materializeSnapped && lodTriPointsSnapped.empty() &&
	    !lodTriPoints.empty())
	    updateSnappedPoints(cut);
	return true;
    }

    if (forceUpdate) {
	forceUpdate = false;
	mesh_lod_active_data_clear(lod);
	shrinkMemory();
	currCut = -1;
	return setCut(cut, materializeSnapped);
    }

    if (cut > currCut && cut <= maxPopCut) {
	if (!lodTriPoints.size()) {
	    if (!triPopLoad(-1, cut, materializeSnapped))
		return false;
	} else if (!triPopLoad(currCut, cut, materializeSnapped)) {
	    return false;
	}
	/* A retained exact-coordinate load may be followed by a snapped-display
	 * request at the same or a higher cut.  Populate that derived array on demand
	 * without reopening any cache chunks. */
	if (materializeSnapped && lodTriPointsSnapped.empty())
	    updateSnappedPoints(cut);
    }

    if (cut < currCut && cut <= maxPopCut &&
	currCut <= maxPopCut) {
	if (!lodTriPoints.size()) {
	    if (!triPopLoad(-1, cut, materializeSnapped))
		return false;
	} else {
	    triPopTrim(cut, materializeSnapped);
	}
    }

    currCut = cut;
    return true;
}

bool
BObolPopState::cacheWrite(const char *component, std::stringstream &stream)
{
    std::string buffer = stream.str();
    return cacheWriteData(component, buffer.data(), buffer.size());
}

bool
BObolPopState::cacheWriteData(const char *component, const void *data,
			      size_t size)
{
    if (!component || !data || !size)
	return false;
    char keystr[64] = {0};
    if (!cacheComponentKey(keystr, sizeof(keystr), component))
	return false;
    size_t wsize = bu_cache_write(const_cast<void *>(data), size,
				  keystr, context->i->lodCache,
				  &writeTxn);
    return (wsize > 0);
}

void
BObolPopState::initializeCacheKeyPrefix(void)
{
    cacheKeyPrefixLength = 0;
    const std::to_chars_result converted =
	std::to_chars(cacheKeyPrefix,
	    cacheKeyPrefix + sizeof(cacheKeyPrefix) - 2, hash);
    if (converted.ec != std::errc())
	return;
    *converted.ptr = ':';
    cacheKeyPrefixLength =
	static_cast<size_t>(converted.ptr - cacheKeyPrefix) + 1;
    cacheKeyPrefix[cacheKeyPrefixLength] = '\0';
}

bool
BObolPopState::cacheComponentKey(char *buffer, size_t bufferSize,
				 const char *component) const
{
    if (!buffer || !component || !cacheKeyPrefixLength)
	return false;
    const size_t componentLength = strlen(component);
    if (cacheKeyPrefixLength >= bufferSize ||
	componentLength >= bufferSize - cacheKeyPrefixLength)
	return false;
    memcpy(buffer, cacheKeyPrefix, cacheKeyPrefixLength);
    memcpy(buffer + cacheKeyPrefixLength, component, componentLength + 1);
    return true;
}

size_t
BObolPopState::cacheGet(void **data, const char *component)
{
    if (context && context->i && context->i->accessMutex &&
	!readLock.owns_lock())
	readLock =
	    std::shared_lock<std::shared_mutex>(*context->i->accessMutex);
    char keystr[64] = {0};
    if (!cacheComponentKey(keystr, sizeof(keystr), component))
	return 0;
    return bu_cache_get(data, keystr, context->i->lodCache,
	&readTxn);
}

void
BObolPopState::cacheDone(void)
{
    bu_cache_get_done(&readTxn);
    if (readLock.owns_lock())
	readLock.unlock();
}

static bool
mesh_lod_chunk_append(std::vector<unsigned char> &bytes,
	const void *data, size_t size)
{
    if (!data || !size || size > SIZE_MAX - bytes.size())
	return false;
    const size_t offset = bytes.size();
    bytes.resize(offset + size);
    memcpy(bytes.data() + offset, data, size);
    return true;
}

bool
BObolPopState::cacheChunks(void)
{
    if (chunkFaces.empty() || chunkFaces.size() != chunkInfos.size() ||
	faceActivationCut.size() != faceCount)
	return false;

    std::vector<uint8_t> localVertexCut(vertexCount, UINT8_MAX);
    std::vector<uint32_t> localVertexIndex(vertexCount, UINT32_MAX);
    std::vector<uint32_t> touchedVertices;
    std::array<std::vector<uint32_t>, POP_CUT_COUNT_MAX> cutVertices;
    std::array<std::vector<uint32_t>, POP_CUT_COUNT_MAX> chunkCutFaces;

    for (size_t chunkIndex = 0; chunkIndex < chunkFaces.size();
	 chunkIndex++) {
	BObolMeshLodChunkInfo &info = chunkInfos[chunkIndex];
	for (std::vector<uint32_t> &values : cutVertices)
	    values.clear();
	for (std::vector<uint32_t> &values : chunkCutFaces)
	    values.clear();
	touchedVertices.clear();

	for (uint32_t sourceFace : chunkFaces[chunkIndex]) {
	    if (sourceFace >= faceCount)
		return false;
	    const uint8_t activation = faceActivationCut[sourceFace];
	    if (activation >= cutCount)
		return false;
	    chunkCutFaces[activation].push_back(sourceFace);
	    for (size_t corner = 0; corner < 3; ++corner) {
		const int vertex = faceArray[
		    3 * static_cast<size_t>(sourceFace) + corner];
		if (vertex < 0 || static_cast<size_t>(vertex) >= vertexCount)
		    return false;
		const size_t sourceVertex = static_cast<size_t>(vertex);
		if (localVertexCut[sourceVertex] == UINT8_MAX) {
		    localVertexCut[sourceVertex] = activation;
		    touchedVertices.push_back(static_cast<uint32_t>(sourceVertex));
		} else {
		    localVertexCut[sourceVertex] = std::min(
			localVertexCut[sourceVertex], activation);
		}
	    }
	}
	for (uint32_t sourceVertex : touchedVertices)
	    cutVertices[localVertexCut[sourceVertex]].push_back(sourceVertex);

	uint64_t cumulativePoints = 0;
	uint64_t cumulativeFaces = 0;
	uint32_t nextLocalVertex = 0;
	for (uint32_t cut = 0; cut < cutCount; ++cut) {
	    for (uint32_t sourceVertex : cutVertices[cut])
		localVertexIndex[sourceVertex] = nextLocalVertex++;
	    cumulativePoints += cutVertices[cut].size();
	    cumulativeFaces += chunkCutFaces[cut].size();
	    if (cumulativePoints > UINT32_MAX || cumulativeFaces > UINT32_MAX)
		return false;
	    BObolMeshLodChunkCutInfo &cutInfo = info.cuts[cut];
	    cutInfo.point_count = static_cast<uint32_t>(cumulativePoints);
	    cutInfo.face_count = static_cast<uint32_t>(cumulativeFaces);
	    uint64_t bytes = cumulativePoints * 3u * sizeof(float) +
		cumulativeFaces * 3u * sizeof(uint32_t);
	    if (hasNormals)
		bytes += cumulativeFaces * 9u * sizeof(float);
	    cutInfo.resident_bytes = bytes;
	}
	if (!cumulativePoints || !cumulativeFaces ||
	    cumulativePoints != touchedVertices.size() ||
	    cumulativeFaces != chunkFaces[chunkIndex].size())
	    return false;

	/* Fixed-width record.  Header size is deliberately a multiple of eight:
	 * cache readers may validate/copy the IEEE doubles without relying on
	 * compiler packing or size_t width. */
	std::vector<unsigned char> bytes;
	const uint64_t pointScalarCount = cumulativePoints * 3u;
	const uint64_t indexCount = cumulativeFaces * 3u;
	const uint64_t normalScalarCount = hasNormals ? cumulativeFaces * 9u : 0;
	const uint64_t expected = 88u +
	    static_cast<uint64_t>(POP_CUT_COUNT_MAX) * 2u * sizeof(uint32_t) +
	    pointScalarCount * sizeof(double) + indexCount * sizeof(uint32_t) +
	    normalScalarCount * sizeof(double);
	if (expected > SIZE_MAX)
	    return false;
	bytes.reserve(static_cast<size_t>(expected));
	const uint32_t diskVersion = 1;
	const uint32_t diskChunk = static_cast<uint32_t>(chunkIndex);
	const uint32_t diskCutCount = cutCount;
	const uint32_t diskFlags = hasNormals ? 1u : 0u;
	const uint32_t diskPointCount = static_cast<uint32_t>(cumulativePoints);
	const uint32_t diskFaceCount = static_cast<uint32_t>(cumulativeFaces);
	const uint32_t reserved = 0;
	const double bounds[6] = {
	    static_cast<double>(info.bmin[X]),
	    static_cast<double>(info.bmin[Y]),
	    static_cast<double>(info.bmin[Z]),
	    static_cast<double>(info.bmax[X]),
	    static_cast<double>(info.bmax[Y]),
	    static_cast<double>(info.bmax[Z])
	};
	if (!mesh_lod_chunk_append(bytes, meshLodChunkMagic,
		sizeof(meshLodChunkMagic)) ||
	    !mesh_lod_chunk_append(bytes, &diskVersion, sizeof(diskVersion)) ||
	    !mesh_lod_chunk_append(bytes, &diskChunk, sizeof(diskChunk)) ||
	    !mesh_lod_chunk_append(bytes, &diskCutCount, sizeof(diskCutCount)) ||
	    !mesh_lod_chunk_append(bytes, &diskFlags, sizeof(diskFlags)) ||
	    !mesh_lod_chunk_append(bytes, &diskPointCount,
		sizeof(diskPointCount)) ||
	    !mesh_lod_chunk_append(bytes, &diskFaceCount,
		sizeof(diskFaceCount)) ||
	    !mesh_lod_chunk_append(bytes, &reserved, sizeof(reserved)) ||
	    !mesh_lod_chunk_append(bytes, &reserved, sizeof(reserved)) ||
	    !mesh_lod_chunk_append(bytes, bounds, sizeof(bounds)))
	    return false;
	for (uint32_t cut = 0; cut < POP_CUT_COUNT_MAX; ++cut) {
	    const uint32_t value = info.cuts[cut].point_count;
	    if (!mesh_lod_chunk_append(bytes, &value, sizeof(value)))
		return false;
	}
	for (uint32_t cut = 0; cut < POP_CUT_COUNT_MAX; ++cut) {
	    const uint32_t value = info.cuts[cut].face_count;
	    if (!mesh_lod_chunk_append(bytes, &value, sizeof(value)))
		return false;
	}
	for (uint32_t cut = 0; cut < cutCount; ++cut) {
	    for (uint32_t sourceVertex : cutVertices[cut]) {
		const double point[3] = {
		    static_cast<double>(vertexArray[sourceVertex][X]),
		    static_cast<double>(vertexArray[sourceVertex][Y]),
		    static_cast<double>(vertexArray[sourceVertex][Z])
		};
		if (!mesh_lod_chunk_append(bytes, point, sizeof(point)))
		    return false;
	    }
	}
	for (uint32_t cut = 0; cut < cutCount; ++cut) {
	    for (uint32_t sourceFace : chunkCutFaces[cut]) {
		uint32_t indices[3] = {};
		for (size_t corner = 0; corner < 3; ++corner) {
		    const int sourceVertex = faceArray[
			3 * static_cast<size_t>(sourceFace) + corner];
		    indices[corner] = localVertexIndex[sourceVertex];
		    if (indices[corner] >= cumulativePoints)
			return false;
		}
		if (!mesh_lod_chunk_append(bytes, indices, sizeof(indices)))
		    return false;
	    }
	}
	if (hasNormals) {
	    if (!normalArray)
		return false;
	    for (uint32_t cut = 0; cut < cutCount; ++cut) {
		for (uint32_t sourceFace : chunkCutFaces[cut]) {
		    for (size_t corner = 0; corner < 3; ++corner) {
			const vect_t &source = normalArray[
			    3 * static_cast<size_t>(sourceFace) + corner];
			const double normal[3] = {
			    static_cast<double>(source[X]),
			    static_cast<double>(source[Y]),
			    static_cast<double>(source[Z])
			};
			if (!mesh_lod_chunk_append(bytes, normal,
				sizeof(normal)))
			    return false;
		    }
		}
	    }
	}
	if (bytes.size() != expected)
	    return false;

	char component[32] = {0};
	snprintf(component, sizeof(component), "%s%08x",
	    CACHE_CHUNK_DATA_PREFIX, diskChunk);
	if (!cacheWriteData(component, bytes.data(), bytes.size()))
	    return false;

	for (uint32_t sourceVertex : touchedVertices) {
	    localVertexCut[sourceVertex] = UINT8_MAX;
	    localVertexIndex[sourceVertex] = UINT32_MAX;
	}
    }

    const uint32_t diskChunkCount = static_cast<uint32_t>(chunkInfos.size());
    if (!cacheWriteData(CACHE_CHUNK_COUNT, &diskChunkCount,
	    sizeof(diskChunkCount)))
	return false;
    std::vector<double> bounds(chunkInfos.size() * 6);
    std::vector<uint8_t> minmax(chunkInfos.size() * 2);
    std::vector<uint32_t> faceCounts(
	chunkInfos.size() * POP_CUT_COUNT_MAX);
    std::vector<uint32_t> pointCounts(
	chunkInfos.size() * POP_CUT_COUNT_MAX);
    std::vector<uint64_t> residentBytes(
	chunkInfos.size() * POP_CUT_COUNT_MAX);
    for (size_t chunk = 0; chunk < chunkInfos.size(); ++chunk) {
	const BObolMeshLodChunkInfo &info = chunkInfos[chunk];
	bounds[chunk * 6 + 0] = info.bmin[X];
	bounds[chunk * 6 + 1] = info.bmin[Y];
	bounds[chunk * 6 + 2] = info.bmin[Z];
	bounds[chunk * 6 + 3] = info.bmax[X];
	bounds[chunk * 6 + 4] = info.bmax[Y];
	bounds[chunk * 6 + 5] = info.bmax[Z];
	minmax[chunk * 2 + 0] = static_cast<uint8_t>(info.min_cut);
	minmax[chunk * 2 + 1] = static_cast<uint8_t>(info.max_cut);
	for (uint32_t cut = 0; cut < POP_CUT_COUNT_MAX; ++cut) {
	    const size_t offset = chunk * POP_CUT_COUNT_MAX + cut;
	    faceCounts[offset] = info.cuts[cut].face_count;
	    pointCounts[offset] = info.cuts[cut].point_count;
	    residentBytes[offset] = info.cuts[cut].resident_bytes;
	}
    }
    return cacheWriteData(CACHE_CHUNK_BOUNDS, bounds.data(),
	       bounds.size() * sizeof(bounds[0])) &&
	cacheWriteData(CACHE_CHUNK_MINMAX, minmax.data(),
	       minmax.size() * sizeof(minmax[0])) &&
	cacheWriteData(CACHE_CHUNK_FACE_COUNTS, faceCounts.data(),
	       faceCounts.size() * sizeof(faceCounts[0])) &&
	cacheWriteData(CACHE_CHUNK_POINT_COUNTS, pointCounts.data(),
	       pointCounts.size() * sizeof(pointCounts[0])) &&
	cacheWriteData(CACHE_CHUNK_RESIDENT_BYTES, residentBytes.data(),
	       residentBytes.size() * sizeof(residentBytes[0]));
}

bool
BObolPopState::cacheTri(void)
{
    {
	const int32_t diskMaxCut = static_cast<int32_t>(maxPopCut);
	std::stringstream stream;
	stream.write(reinterpret_cast<const char *>(&diskMaxCut),
		     sizeof(diskMaxCut));
	if (!cacheWrite(CACHE_POP_MAX_CUT, stream))
	    return false;
    }

    {
	const int32_t diskMinCut = static_cast<int32_t>(minPopCut);
	std::stringstream stream;
	stream.write(reinterpret_cast<const char *>(&diskMinCut),
		     sizeof(diskMinCut));
	if (!cacheWrite(CACHE_POP_MIN_CUT, stream))
	    return false;
    }

    {
	uint8_t diskCutBits[POP_CUT_COUNT_MAX][3] = {{0}};
	for (uint32_t cut = 0; cut < cutCount; ++cut)
	    for (int axis = 0; axis < 3; ++axis)
		diskCutBits[cut][axis] = cuts[cut].bits[axis];
	if (!cacheWriteData(CACHE_CUT_QUANTIZATION, diskCutBits,
		sizeof(diskCutBits)))
	    return false;
    }

    {
	const uint8_t diskCullBackfaces = shadedCullBackfaces ? 1u : 0u;
	std::stringstream stream;
	stream.write(reinterpret_cast<const char *>(&diskCullBackfaces),
		     sizeof(diskCullBackfaces));
	if (!cacheWrite(CACHE_SHADED_CULL_BACKFACES, stream))
	    return false;
    }

    {
	const uint8_t diskHasNormals = hasNormals ? 1u : 0u;
	std::stringstream stream;
	stream.write(reinterpret_cast<const char *>(&diskHasNormals),
		     sizeof(diskHasNormals));
	if (!cacheWrite(CACHE_HAS_NORMALS, stream))
	    return false;
    }

    if (clusterInfos.size() != BOBOL_MESH_LOD_CLUSTER_COUNT ||
	clusterRanges.empty())
	return false;
    {
	const uint16_t diskGrid = BOBOL_MESH_LOD_CLUSTER_GRID_RESOLUTION;
	if (!cacheWriteData(CACHE_CLUSTER_GRID, &diskGrid,
		sizeof(diskGrid)))
	    return false;
    }
    {
	std::vector<double> diskBounds(clusterInfos.size() * 6, 0.0);
	for (size_t cluster = 0; cluster < clusterInfos.size(); ++cluster) {
	    const BObolMeshLodClusterInfo &info = clusterInfos[cluster];
	    if (!info.range_count)
		continue;
	    diskBounds[cluster * 6 + 0] = info.bmin[X];
	    diskBounds[cluster * 6 + 1] = info.bmin[Y];
	    diskBounds[cluster * 6 + 2] = info.bmin[Z];
	    diskBounds[cluster * 6 + 3] = info.bmax[X];
	    diskBounds[cluster * 6 + 4] = info.bmax[Y];
	    diskBounds[cluster * 6 + 5] = info.bmax[Z];
	}
	if (!cacheWriteData(CACHE_CLUSTER_BOUNDS, diskBounds.data(),
		diskBounds.size() * sizeof(diskBounds[0])))
	    return false;
    }
    {
	std::vector<uint32_t> offsets(clusterInfos.size() + 1, 0);
	for (size_t cluster = 0; cluster < clusterInfos.size(); ++cluster) {
	    const uint64_t next = static_cast<uint64_t>(offsets[cluster]) +
		clusterInfos[cluster].range_count;
	    if (next > UINT32_MAX)
		return false;
	    offsets[cluster + 1] = static_cast<uint32_t>(next);
	}
	if (offsets.back() != clusterRanges.size() ||
	    !cacheWriteData(CACHE_CLUSTER_RANGE_OFFSETS, offsets.data(),
		offsets.size() * sizeof(offsets[0])))
	    return false;
    }
    {
	std::vector<BObolMeshLodClusterRangeDisk> diskRanges(
	    clusterRanges.size());
	for (size_t range = 0; range < clusterRanges.size(); ++range) {
	    diskRanges[range].firstIndex = clusterRanges[range].first_index;
	    diskRanges[range].indexCount = clusterRanges[range].index_count;
	    diskRanges[range].activationCut =
		clusterRanges[range].activation_cut;
	    memset(diskRanges[range].reserved, 0,
		sizeof(diskRanges[range].reserved));
	}
	if (!cacheWriteData(CACHE_CLUSTER_RANGES, diskRanges.data(),
		diskRanges.size() * sizeof(diskRanges[0])))
	    return false;
    }

    {
	std::stringstream stream;
	for (size_t cutIndex = 0; cutIndex <= POP_CUT_COUNT_MAX; ++cutIndex) {
	    uint64_t count = 0;
	    if (cutIndex >= cutTriVerts.size()) {
		stream.write(reinterpret_cast<const char *>(&count),
			     sizeof(count));
		continue;
	    }
	    if (static_cast<int>(cutIndex) > maxPopCut ||
		cutTriVerts[cutIndex].empty()) {
		stream.write(reinterpret_cast<const char *>(&count),
			     sizeof(count));
		continue;
	    }
	    count = static_cast<uint64_t>(cutTriVerts[cutIndex].size());
	    stream.write(reinterpret_cast<const char *>(&count), sizeof(count));
	}
	if (!cacheWrite(CACHE_VERTEX_COUNT, stream))
	    return false;
    }

    {
	std::stringstream stream;
	for (size_t cutIndex = 0; cutIndex <= POP_CUT_COUNT_MAX; ++cutIndex) {
	    uint64_t count = 0;
	    if (cutIndex >= cutTris.size() ||
		static_cast<int>(cutIndex) > maxPopCut ||
		!cutTris[cutIndex].size()) {
		stream.write(reinterpret_cast<const char *>(&count),
			     sizeof(count));
		continue;
	    }
	    count = static_cast<uint64_t>(cutTris[cutIndex].size());
	    stream.write(reinterpret_cast<const char *>(&count), sizeof(count));
	}
	if (!cacheWrite(CACHE_TRI_COUNT, stream))
	    return false;
    }

    struct bu_vls keyBuffer = BU_VLS_INIT_ZERO;

    for (int cutIndex = 0; cutIndex <= maxPopCut; ++cutIndex) {
	if (cutIndex >= static_cast<int>(cutTriVerts.size()) ||
	    cutTriVerts[cutIndex].empty())
	    continue;
	std::vector<fastf_t> points;
	points.resize(cutTriVerts[cutIndex].size() * 3);
	size_t outputIndex = 0;
	for (uint32_t sourceVertexIndex : cutTriVerts[cutIndex]) {
	    points[outputIndex++] = vertexArray[sourceVertexIndex][X];
	    points[outputIndex++] = vertexArray[sourceVertexIndex][Y];
	    points[outputIndex++] = vertexArray[sourceVertexIndex][Z];
	}
	bu_vls_sprintf(&keyBuffer, "%s%d", CACHE_VERT_CUT, cutIndex);
	if (!cacheWriteData(bu_vls_cstr(&keyBuffer), points.data(),
		points.size() * sizeof(fastf_t))) {
	    bu_vls_free(&keyBuffer);
	    return false;
	}
    }

    for (int cutIndex = 0; cutIndex <= maxPopCut; ++cutIndex) {
	if (cutTris[cutIndex].empty())
	    continue;
	std::vector<uint32_t> triangles;
	triangles.resize(cutTris[cutIndex].size() * 3);
	size_t outputIndex = 0;
	for (uint32_t sourceFaceIndex : cutTris[cutIndex]) {
	    triangles[outputIndex++] =
		triIndexMap[faceArray[3 * sourceFaceIndex + 0]];
	    triangles[outputIndex++] =
		triIndexMap[faceArray[3 * sourceFaceIndex + 1]];
	    triangles[outputIndex++] =
		triIndexMap[faceArray[3 * sourceFaceIndex + 2]];
	}
	bu_vls_sprintf(&keyBuffer, "%s%d", CACHE_TRI_CUT, cutIndex);
	if (!cacheWriteData(bu_vls_cstr(&keyBuffer), triangles.data(),
		triangles.size() * sizeof(uint32_t))) {
	    bu_vls_free(&keyBuffer);
	    return false;
	}
    }

    if (normalArray) {
	for (int cutIndex = 0; cutIndex <= maxPopCut; ++cutIndex) {
	    if (cutTris[cutIndex].empty())
		continue;
	    std::vector<fastf_t> normals;
	    normals.resize(cutTris[cutIndex].size() * 9);
	    size_t outputIndex = 0;
	    for (uint32_t sourceFaceIndex : cutTris[cutIndex]) {
		for (int cornerIndex = 0; cornerIndex < 3; cornerIndex++) {
		    const size_t normalIndex =
			3 * static_cast<size_t>(sourceFaceIndex) + cornerIndex;
		    normals[outputIndex++] = normalArray[normalIndex][X];
		    normals[outputIndex++] = normalArray[normalIndex][Y];
		    normals[outputIndex++] = normalArray[normalIndex][Z];
		}
	    }
	    bu_vls_sprintf(&keyBuffer, "%s%d", CACHE_VERTNORM_CUT,
			   cutIndex);
	    if (!cacheWriteData(bu_vls_cstr(&keyBuffer), normals.data(),
		    normals.size() * sizeof(fastf_t))) {
		bu_vls_free(&keyBuffer);
		return false;
	    }
	}
    }

    bu_vls_free(&keyBuffer);
    return cacheChunks();
}

void
BObolPopState::cache(void)
{
    if (!hash) {
	isValid = false;
	return;
    }

    int writeSem = mesh_lod_cache_write_semaphore();
    bu_semaphore_acquire(writeSem);
    std::unique_lock<std::shared_mutex> cacheLock(
	*context->i->accessMutex);

    {
	const double diskBounds[12] = {
	    static_cast<double>(bbmin[X]),
	    static_cast<double>(bbmin[Y]),
	    static_cast<double>(bbmin[Z]),
	    static_cast<double>(bbmax[X]),
	    static_cast<double>(bbmax[Y]),
	    static_cast<double>(bbmax[Z]),
	    static_cast<double>(minx),
	    static_cast<double>(miny),
	    static_cast<double>(minz),
	    static_cast<double>(maxx),
	    static_cast<double>(maxy),
	    static_cast<double>(maxz)
	};
	std::stringstream stream;
	stream.write(reinterpret_cast<const char *>(diskBounds),
	    sizeof(diskBounds));
	isValid = cacheWrite(CACHE_OBJ_BOUNDS, stream);
    }

    if (!isValid) {
	bu_cache_write_abort(&writeTxn);
	cacheLock.unlock();
	bu_semaphore_release(writeSem);
	return;
    }

    isValid = cacheTri();

    if (writeTxn) {
	if (isValid) {
	    if (bu_cache_write_commit(context->i->lodCache, &writeTxn) !=
		BRLCAD_OK)
		isValid = false;
	} else {
	    bu_cache_write_abort(&writeTxn);
	}
    }

    cacheLock.unlock();
    bu_semaphore_release(writeSem);
}

fastf_t
BObolPopState::snap(fastf_t value, fastf_t min, fastf_t max, uint8_t bits)
{
    if (!(max > min) || !bits)
	return value;
    const double mask = std::ldexp(
	1.0, BOBOL_MESH_LOD_QUANTIZATION_BITS -
	    std::min<int>(BOBOL_MESH_LOD_QUANTIZATION_BITS, bits));
    const double scaled = std::max(0.0, std::min(
	static_cast<double>(USHRT_MAX),
	(static_cast<double>(value) - min) /
	    (static_cast<double>(max) - min) * USHRT_MAX));
    /* triProcess assigns activation from the high-bit prefix of the floored
     * 16-bit code.  Render exactly that cell's representative.  The former
     * floor/ceil average selected a cell boundary when a source coordinate
     * lay exactly on one; the activation classifier selected the upper cell,
     * so an omitted face could nevertheless have distinct displayed
     * vertices.  That disagreement produced transient cracks at coarse
     * cuts. */
    const double code = floor(scaled);
    const double cell = floor(code / mask);
    const double snapped = std::min(
	static_cast<double>(USHRT_MAX), (cell + 0.5) * mask);
    return ((snapped / USHRT_MAX) * (max - min)) + min;
}

void
BObolPopState::cutPoint(point_t *out, const point_t *point, int cut)
{
    /* The terminal cut is the stable exact state.  All source vertices and
     * faces are already present in the cumulative cache at this point, so
     * retaining the original coordinate costs no additional I/O and avoids
     * imposing a permanent 16-bit quantization error on very large extents. */
    if (cut >= maxPopCut) {
	VMOVE(*out, *point);
	return;
    }
    cut = std::max(0, std::min(maxPopCut, cut));
    fastf_t nx = snap((*point)[X], minx, maxx, cuts[cut].bits[X]);
    fastf_t ny = snap((*point)[Y], miny, maxy, cuts[cut].bits[Y]);
    fastf_t nz = snap((*point)[Z], minz, maxz, cuts[cut].bits[Z]);
    VSET(*out, nx, ny, nz);
}

static void
mesh_lod_active_pop_data_publish(struct BObolMeshLod *lod,
				 BObolPopState *state)
{
    if (!lod || !state)
	return;

    lod->fcnt = static_cast<int>(state->lodTris.size() / 3);
    lod->faces = state->lodTris.data();
    lod->points_orig = reinterpret_cast<const point_t *>(
			   state->lodTriPoints.data());
    lod->porig_cnt = static_cast<int>(state->lodTriPoints.size() / 3);
    lod->normals = (state->lodTriNormals.size() >=
		    state->lodTris.size() * 3) ?
		   reinterpret_cast<const vect_t *>(state->lodTriNormals.data()) : NULL;
    const std::vector<fastf_t> &displayPoints =
	state->lodTriPointsSnapped.empty() ?
	state->lodTriPoints : state->lodTriPointsSnapped;
    lod->points = reinterpret_cast<const point_t *>(displayPoints.data());
    lod->pcnt = static_cast<int>(displayPoints.size() / 3);
}

static unsigned long long
mesh_lod_cache_generate(struct BObolMeshLodContext *context,
			const point_t *vertices,
			size_t vertexCount,
			const vect_t *normals,
			const int *faces,
			size_t faceCount,
			unsigned long long userKey,
			bool shadedCullBackfaces)
{
    if (!context || !vertices || !vertexCount || !faces || !faceCount)
	return 0;

    BObolPopState state(context, vertices, vertexCount, normals, faces,
			  faceCount, userKey, shadedCullBackfaces);
    return state.isValid ? state.hash : 0;
}

static struct BObolMeshLod *
mesh_lod_create(struct BObolMeshLodContext *context,
		unsigned long long key,
		bool retainHeaderSnapshot = false)
{
    if (!context || !key)
	return NULL;

    BObolPopState *state = new (std::nothrow) BObolPopState(
	context, key, retainHeaderSnapshot);
    if (!state)
	return NULL;

    if (!state->isValid) {
	delete state;
	return NULL;
    }

    struct BObolMeshLod *lod;
    BU_GET(lod, struct BObolMeshLod);
    lod->context = context;
    lod->state = state;
    mesh_lod_active_data_clear(lod);
    VMOVE(lod->bmin, state->bbmin);
    VMOVE(lod->bmax, state->bbmax);
    state->lod = lod;
    return lod;
}

static struct BObolMeshLod *
mesh_lod_adopt_generated_state(struct BObolMeshLodContext *context,
			       BObolPopState *state)
{
    if (!context || !state || !state->isValid)
	return NULL;
    struct BObolMeshLod *lod;
    BU_GET(lod, struct BObolMeshLod);
    lod->context = context;
    lod->state = state;
    mesh_lod_active_data_clear(lod);
    VMOVE(lod->bmin, state->bbmin);
    VMOVE(lod->bmax, state->bbmax);
    state->lod = lod;
    return lod;
}

static int
mesh_lod_cut(struct BObolMeshLod *lod, int cut, int reset,
	       bool materializeSnapped = true)
{
    if (!lod || !lod->state)
	return -1;
    BObolPopState *state = lod->state;
    if (cut < 0)
	return state->currCut;

    if (reset)
	state->forceUpdate = true;
    if (!state->setCut(cut, materializeSnapped)) {
	mesh_lod_active_data_clear(lod);
	return -1;
    }

    if (state->currCut <= state->maxPopCut)
	mesh_lod_active_pop_data_publish(lod, state);

    return state->currCut;
}

static int
mesh_lod_payload_available(struct BObolMeshLodContext *context,
			   unsigned long long key)
{
    if (!context || !key)
	return 0;

    struct BObolMeshLod *lod = mesh_lod_create(context, key);
    if (!lod)
	return 0;

    lod->context = NULL;
    bobol_mesh_lod_destroy(lod);
    return 1;
}

static void
mesh_lod_status_current(struct db_i *dbip,
			struct BObolMeshLodContext *context,
			const char *name,
			struct BObolMeshLodCacheStatus *status)
{
    if (!status)
	return;

    bobol_mesh_lod_cache_status_init(status);
    if (!dbip || !name)
	return;

    struct directory *dp = db_lookup(dbip, name, LOOKUP_QUIET);
    status->directory_found = (dp != RT_DIR_NULL) ? 1 : 0;
    status->is_bot = (dp != RT_DIR_NULL &&
		      dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) ? 1 : 0;

    if (!context)
	return;

    status->cache_key = mesh_lod_key_get(context, name);
    status->has_cache_key = status->cache_key ? 1 : 0;
    status->has_cached_payload = mesh_lod_payload_available(
				     context, status->cache_key);
    status->stale_cache_entry = (status->has_cache_key &&
				 !status->has_cached_payload) ? 1 : 0;
}

static void
mesh_lod_cache_del(struct BObolMeshLodContext *context,
		   unsigned long long hash,
		   const char *component)
{
    std::string keystr = std::to_string(hash) + std::string(":") +
			 std::string(component);
    bu_cache_clear(keystr.c_str(), context->i->lodCache, NULL);
}

static void
mesh_lod_cache_clear_context(struct BObolMeshLodContext *context,
			     unsigned long long key)
{
    char dir[MAXPATHLEN];
    std::unique_lock<std::shared_mutex> lock;
    std::unique_lock<std::shared_mutex> nameLock;

    if (context && context->i && context->i->accessMutex) {
	lock =
	    std::unique_lock<std::shared_mutex>(*context->i->accessMutex);
	if (context->i->nameMutex)
	    nameLock =
		std::unique_lock<std::shared_mutex>(*context->i->nameMutex);
    }

    if (context && key) {
	mesh_lod_cache_del(context, key, CACHE_POP_MAX_CUT);
	mesh_lod_cache_del(context, key, CACHE_VERTEX_COUNT);
	mesh_lod_cache_del(context, key, CACHE_TRI_COUNT);
	mesh_lod_cache_del(context, key, CACHE_OBJ_BOUNDS);

	char **keysv = NULL;
	int nkeys = bu_cache_keys(&keysv, context->i->lodCache);
	std::string prefix = std::to_string(key) + std::string(":");
	for (int keyIndex = 0; keyIndex < nkeys; keyIndex++) {
	    if (bu_strncmp(keysv[keyIndex], prefix.c_str(),
			   prefix.length()) == 0)
		bu_cache_clear(keysv[keyIndex], context->i->lodCache, NULL);
	}
	if (nkeys)
	    bu_argv_free(static_cast<size_t>(nkeys), keysv);

	keysv = NULL;
	nkeys = bu_cache_keys(&keysv, context->i->nameCache);
	for (int keyIndex = 0; keyIndex < nkeys; keyIndex++) {
	    void *data = NULL;
	    size_t dsize = bu_cache_get(&data, keysv[keyIndex],
					context->i->nameCache, NULL);
	    if (dsize == sizeof(unsigned long long) && data) {
		unsigned long long storedKey = *(unsigned long long *)data;
		bu_free(data, "name cache data");
		if (storedKey == key)
		    bu_cache_clear(keysv[keyIndex], context->i->nameCache, NULL);
	    } else if (data) {
		bu_free(data, "name cache data");
	    }
	}
	if (nkeys)
	    bu_argv_free(static_cast<size_t>(nkeys), keysv);
	for (auto it = context->i->nameKeys->begin();
	     it != context->i->nameKeys->end();) {
	    if (it->second == key)
		it = context->i->nameKeys->erase(it);
	    else
		++it;
	}
	return;
    }

    if (context && !key) {
	context->i->nameKeys->clear();
	char **keysv = NULL;
	int nkeys = bu_cache_keys(&keysv, context->i->lodCache);
	for (int keyIndex = 0; keyIndex < nkeys; keyIndex++)
	    bu_cache_clear(keysv[keyIndex], context->i->lodCache, NULL);
	if (nkeys)
	    bu_argv_free(static_cast<size_t>(nkeys), keysv);

	keysv = NULL;
	nkeys = bu_cache_keys(&keysv, context->i->nameCache);
	for (int keyIndex = 0; keyIndex < nkeys; keyIndex++)
	    bu_cache_clear(keysv[keyIndex], context->i->nameCache, NULL);
	if (nkeys)
	    bu_argv_free(static_cast<size_t>(nkeys), keysv);
	return;
    }

    bu_dir(dir, MAXPATHLEN, BU_DIR_CACHE, POP_CACHEDIR, NULL);
    bu_dirclear((const char *)dir);
}

void
bobol_mesh_lod_cache_init(struct db_i *dbip, int verbose)
{
    if (!dbip)
	return;

    int completed = 0;
    int target = 0;
    struct directory *dp = RT_DIR_NULL;
    FOR_ALL_DIRECTORY_START(dp, dbip)
    if (dp->d_addr == RT_DIR_PHONY_ADDR)
	continue;
    if (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT)
	target++;
    FOR_ALL_DIRECTORY_END;

    int64_t start = bu_gettime();
    int64_t overallStart = start;
    FOR_ALL_DIRECTORY_START(dp, dbip)
    if (dp->d_addr == RT_DIR_PHONY_ADDR)
	continue;
    if (dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BOT)
	continue;

    struct BObolMeshLodCacheStatus status =
	    BOBOL_MESH_LOD_CACHE_STATUS_INIT;
    if (bobol_mesh_lod_cache_status(dbip, dp->d_namep, &status) == BRLCAD_OK &&
	status.has_cache_key && status.has_cached_payload)
	continue;

    if (verbose > 1)
	bu_log("Processing(%d):  %s\n", completed + 1, dp->d_namep);

    (void)bobol_mesh_lod_cache_update(dbip, dp->d_namep);
    completed++;

    int64_t elapsed = bu_gettime() - start;
    fastf_t seconds = elapsed / 1000000.0;
    if (verbose > 1)
	bu_log("Completed. (%g seconds)", seconds);
    if (seconds > 5.0) {
	if (verbose) {
	    elapsed = bu_gettime() - overallStart;
	    seconds = elapsed / 1000000.0;
	    bu_log("LoD cache processing (%g seconds): completed %d of %d BoTs\n",
		   seconds, completed, target);
	}
	start = bu_gettime();
    }
    FOR_ALL_DIRECTORY_END;

    int64_t elapsed = bu_gettime() - overallStart;
    int totalSeconds = static_cast<int>(elapsed / 1000000);
    int totalMinutes = totalSeconds / 60;
    int totalHours = totalMinutes / 60;
    totalMinutes = totalMinutes % 60;
    totalSeconds = totalSeconds % 60;
    bu_log("Mesh LoD caching complete (Elapsed time: %02d:%02d:%02d)\n",
	   totalHours, totalMinutes, totalSeconds);
}

void
bobol_mesh_lod_cache_clear_database(struct db_i *dbip)
{
    struct BObolMeshLodContext *context = mesh_lod_context_create_for_db(dbip);
    if (!context)
	return;
    mesh_lod_cache_clear_context(context, 0);
    mesh_lod_context_destroy(context);
}

void
bobol_mesh_lod_cache_clear_all(void)
{
    mesh_lod_cache_clear_context(NULL, 0);
}

int
bobol_mesh_lod_cache_update(struct db_i *dbip, const char *name)
{
    return bobol_mesh_lod_cache_refresh(dbip, name, NULL);
}

void
bobol_mesh_lod_cache_status_init(struct BObolMeshLodCacheStatus *status)
{
    struct BObolMeshLodCacheStatus defaults =
	    BOBOL_MESH_LOD_CACHE_STATUS_INIT;
    if (status)
	*status = defaults;
}

int
bobol_mesh_lod_cache_status(struct db_i *dbip,
			      const char *name,
			      struct BObolMeshLodCacheStatus *status)
{
    if (!dbip || !name || !status)
	return BRLCAD_ERROR;

    struct BObolMeshLodContext *context = mesh_lod_context_create_for_db(dbip);
    if (!context)
	return BRLCAD_ERROR;

    mesh_lod_status_current(dbip, context, name, status);
    mesh_lod_context_destroy(context);
    return BRLCAD_OK;
}

int
bobol_mesh_lod_cache_invalidate(struct db_i *dbip,
				  const char *name,
				  struct BObolMeshLodCacheStatus *status)
{
    struct BObolMeshLodCacheStatus current =
	    BOBOL_MESH_LOD_CACHE_STATUS_INIT;

    if (status)
	bobol_mesh_lod_cache_status_init(status);
    if (!dbip || !name)
	return BRLCAD_ERROR;

    struct BObolMeshLodContext *context = mesh_lod_context_create_for_db(dbip);
    if (!context)
	return BRLCAD_ERROR;

    mesh_lod_status_current(dbip, context, name, &current);
    if (current.has_cache_key) {
	current.cleared_cache_entry = 1;
	current.cleared_cache_key = current.cache_key;
	mesh_lod_cache_clear_context(context, current.cache_key);
	current.cache_key = 0;
	current.has_cache_key = 0;
	current.has_cached_payload = 0;
	current.stale_cache_entry = 0;
    }

    if (status)
	*status = current;
    mesh_lod_context_destroy(context);
    return BRLCAD_OK;
}

static int
mesh_lod_cache_refresh_impl(struct db_i *dbip, const char *name,
			    const struct rt_bot_internal *stagedBot,
			    bool forceRefresh,
			    struct BObolMeshLodCacheStatus *status,
			    struct BObolMeshLod **generatedLod)
{
    struct BObolMeshLodCacheStatus current =
	    BOBOL_MESH_LOD_CACHE_STATUS_INIT;

    if (status)
	bobol_mesh_lod_cache_status_init(status);
    if (generatedLod)
	*generatedLod = NULL;
    if (!dbip || !name)
	return BRLCAD_ERROR;

    struct BObolMeshLodContext *context = mesh_lod_context_create_for_db(dbip);
    if (!context)
	return BRLCAD_ERROR;

    mesh_lod_status_current(dbip, context, name, &current);
    if (!forceRefresh && current.has_cache_key &&
	current.has_cached_payload && !current.stale_cache_entry) {
	if (status)
	    *status = current;
	mesh_lod_context_destroy(context);
	return BRLCAD_OK;
    }
    if (current.has_cache_key) {
	current.cleared_cache_entry = 1;
	current.cleared_cache_key = current.cache_key;
	mesh_lod_cache_clear_context(context, current.cache_key);
	current.cache_key = 0;
	current.has_cache_key = 0;
	current.has_cached_payload = 0;
	current.stale_cache_entry = 0;
    }

    struct directory *dp = db_lookup(dbip, name, LOOKUP_QUIET);
    current.directory_found = (dp != RT_DIR_NULL) ? 1 : 0;
    current.is_bot = (dp != RT_DIR_NULL &&
		      dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) ? 1 : 0;
    if (dp == RT_DIR_NULL || dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BOT) {
	if (status)
	    *status = current;
	mesh_lod_context_destroy(context);
	return BRLCAD_OK;
    }

    struct rt_db_internal dbintern;
    RT_DB_INTERNAL_INIT(&dbintern);
    bool ownsInternal = false;
    const struct rt_bot_internal *bot = stagedBot;
    if (!bot) {
	if (rt_db_get_internal(&dbintern, dp, dbip, NULL) < 0) {
	    mesh_lod_context_destroy(context);
	    return BRLCAD_ERROR;
	}
	ownsInternal = true;
	if (dbintern.idb_minor_type != DB5_MINORTYPE_BRLCAD_BOT) {
	    rt_db_free_internal(&dbintern);
	    mesh_lod_context_destroy(context);
	    return BRLCAD_ERROR;
	}
	bot = static_cast<const struct rt_bot_internal *>(dbintern.idb_ptr);
    }
    RT_BOT_CK_MAGIC(bot);

    if (!mesh_lod_arrays_validate(bot->faces, bot->num_faces,
				  reinterpret_cast<const point_t *>(bot->vertices),
				  bot->num_vertices,
				  reinterpret_cast<const point_t *>(bot->vertices),
				  bot->num_vertices, NULL, NULL, NULL)) {
	if (ownsInternal)
	    rt_db_free_internal(&dbintern);
	if (status)
	    *status = current;
	mesh_lod_context_destroy(context);
	return BRLCAD_ERROR;
    }

    /* Obol's cull-safe contract is exterior CCW.  The BoT's authored
     * mode/orientation is the source of that semantic fact; proving it again
     * with bg_trimesh_solid2 sorts every half edge and used roughly 40% of
     * Lucy's true-cold CPU before any usable PoP prefix could be published.
     *
     * An explicitly oriented solid may therefore opt into culling directly.
     * An unoriented or non-solid BoT remains two-sided.  Topology validation
     * belongs in authoring/validation workflows (and may populate these
     * properties ahead of drawing), not on the interactive display critical
     * path. */
    std::vector<int> normalizedFaces;
    const int *cacheFaces = bot->faces;
    bool flippedWinding = bot->orientation == RT_BOT_CW;
    const auto applyWinding = [&]() {
	if (!flippedWinding) {
	    normalizedFaces.clear();
	    cacheFaces = bot->faces;
	    return;
	}
	normalizedFaces.resize(bot->num_faces * 3);
	for (size_t faceIndex = 0; faceIndex < bot->num_faces; faceIndex++) {
	    normalizedFaces[3 * faceIndex] = bot->faces[3 * faceIndex];
	    normalizedFaces[3 * faceIndex + 1] =
		bot->faces[3 * faceIndex + 2];
	    normalizedFaces[3 * faceIndex + 2] =
		bot->faces[3 * faceIndex + 1];
	}
	cacheFaces = normalizedFaces.data();
    };
    applyWinding();

    const bool cullBackfaces =
	bot->mode == RT_BOT_SOLID &&
	bot->orientation != RT_BOT_UNORIENTED;

    std::vector<fastf_t> authoredNormals;
    const point_t *botVertices =
	reinterpret_cast<const point_t *>(bot->vertices);
    const vect_t *botNormals =
	mesh_lod_bot_authored_normals(
	    authoredNormals, bot, cacheFaces, flippedWinding) ?
	reinterpret_cast<const vect_t *>(authoredNormals.data()) : NULL;

    BObolPopState *generatedState = NULL;
    unsigned long long key = 0;
    if (generatedLod) {
	generatedState = new (std::nothrow) BObolPopState(
	    context, botVertices, bot->num_vertices,
	    botNormals, cacheFaces, bot->num_faces, 0,
	    cullBackfaces);
	key = generatedState && generatedState->isValid ?
	    generatedState->hash : 0;
    } else {
	key = mesh_lod_cache_generate(
	    context, botVertices, bot->num_vertices,
	    botNormals, cacheFaces, bot->num_faces, 0,
	    cullBackfaces);
    }
    if (!key || mesh_lod_key_put(context, dp->d_namep, key) != 0) {
	delete generatedState;
	if (ownsInternal)
	    rt_db_free_internal(&dbintern);
	mesh_lod_context_destroy(context);
	return BRLCAD_ERROR;
    }

    if (generatedState) {
	struct BObolMeshLod *lod =
	    mesh_lod_adopt_generated_state(context, generatedState);
	if (!lod ||
	    !generatedState->materializeInitialPrefix(
		generatedState->minPopCut)) {
	    if (lod) {
		lod->context = NULL;
		bobol_mesh_lod_destroy(lod);
	    } else {
		delete generatedState;
	    }
	    if (ownsInternal)
		rt_db_free_internal(&dbintern);
	    mesh_lod_context_destroy(context);
	    return BRLCAD_ERROR;
	}
	mesh_lod_active_pop_data_publish(lod, generatedState);
	generatedState->releaseGenerationScratch();
	*generatedLod = lod;
    }

    if (ownsInternal)
	rt_db_free_internal(&dbintern);

    current.cache_key = key;
    current.has_cache_key = 1;
    current.has_cached_payload = generatedLod && *generatedLod ?
	1 : mesh_lod_payload_available(context, key);
    current.stale_cache_entry = current.has_cached_payload ? 0 : 1;
    current.generated_cache_entry = current.has_cached_payload ? 1 : 0;
    if (status)
	*status = current;

    if (!generatedLod || !*generatedLod)
	mesh_lod_context_destroy(context);
    return current.has_cached_payload ? BRLCAD_OK : BRLCAD_ERROR;
}

int
bobol_mesh_lod_cache_refresh(struct db_i *dbip,
			       const char *name,
			       struct BObolMeshLodCacheStatus *status)
{
    return mesh_lod_cache_refresh_impl(
	dbip, name, NULL, true, status, NULL);
}

int
bobol_mesh_lod_cache_refresh_from_bot(
    struct db_i *dbip,
    const char *name,
    const struct rt_bot_internal *bot,
    struct BObolMeshLodCacheStatus *status)
{
    if (!bot)
	return BRLCAD_ERROR;
    return mesh_lod_cache_refresh_impl(
	dbip, name, bot, false, status, NULL);
}

struct BObolMeshLod *
bobol_mesh_lod_cache_refresh_from_bot_open(
    struct db_i *dbip,
    const char *name,
    const struct rt_bot_internal *bot,
    struct BObolMeshLodCacheStatus *status)
{
    if (!bot)
	return NULL;
    struct BObolMeshLod *lod = NULL;
    if (mesh_lod_cache_refresh_impl(
	    dbip, name, bot, false, status, &lod) != BRLCAD_OK) {
	if (lod)
	    bobol_mesh_lod_destroy(lod);
	return NULL;
    }
    if (lod)
	return lod;
    if (status && status->has_cache_key)
	return bobol_mesh_lod_get_cached_prefix(
	    dbip, status->cache_key);
    return bobol_mesh_lod_get_named_cached_prefix(dbip, name);
}

static int
mesh_lod_cache_store_mesh_impl(
    struct db_i *dbip,
    const char *name,
    const point_t *vertices,
    size_t vertexCount,
    const vect_t *normals,
    const int *faces,
    size_t faceCount,
    unsigned long long userKey,
    int shadedCullBackfaces,
    struct BObolMeshLodCacheStatus *status,
    bool preserveVariants)
{
    struct BObolMeshLodCacheStatus current =
	    BOBOL_MESH_LOD_CACHE_STATUS_INIT;

    if (status)
	bobol_mesh_lod_cache_status_init(status);
    if (!dbip || !name)
	return BRLCAD_ERROR;

    struct BObolMeshLodContext *context = mesh_lod_context_create_for_db(dbip);
    if (!context)
	return BRLCAD_ERROR;

    mesh_lod_status_current(dbip, context, name, &current);
    if (!current.directory_found ||
	!vertices || !vertexCount || !faces || !faceCount ||
	!mesh_lod_arrays_validate(faces, faceCount, vertices, vertexCount,
				  vertices, vertexCount, NULL, NULL, NULL)) {
	if (status)
	    *status = current;
	mesh_lod_context_destroy(context);
	return BRLCAD_ERROR;
    }

    if (!preserveVariants && current.has_cache_key &&
	(current.cache_key != userKey || !current.has_cached_payload)) {
	current.cleared_cache_entry = 1;
	current.cleared_cache_key = current.cache_key;
	mesh_lod_cache_clear_context(context, current.cache_key);
	current.cache_key = 0;
	current.has_cache_key = 0;
	current.has_cached_payload = 0;
	current.stale_cache_entry = 0;
    }

    if (current.has_cache_key && current.has_cached_payload &&
	(!preserveVariants || current.cache_key == userKey)) {
	if (status)
	    *status = current;
	mesh_lod_context_destroy(context);
	return BRLCAD_OK;
    }

    /* The caller owns the cull-safety decision.  Re-running a complete
     * half-edge manifold proof here makes cache population scale as a second
     * full topology build and defeats progressive time-to-first-content. */
    const bool cullBackfaces = shadedCullBackfaces != 0;

    unsigned long long key = mesh_lod_cache_generate(
				 context, vertices, vertexCount, normals, faces,
				 faceCount, userKey, cullBackfaces);
    if (!key || mesh_lod_key_put(context, name, key) != 0) {
	if (status)
	    *status = current;
	mesh_lod_context_destroy(context);
	return BRLCAD_ERROR;
    }

    current.cache_key = key;
    current.has_cache_key = 1;
    current.has_cached_payload = mesh_lod_payload_available(context, key);
    current.stale_cache_entry = current.has_cached_payload ? 0 : 1;
    current.generated_cache_entry = current.has_cached_payload ? 1 : 0;
    if (status)
	*status = current;

    mesh_lod_context_destroy(context);
    return current.has_cached_payload ? BRLCAD_OK : BRLCAD_ERROR;
}

int
bobol_mesh_lod_cache_store_mesh(
    struct db_i *dbip,
    const char *name,
    const point_t *vertices,
    size_t vertexCount,
    const vect_t *normals,
    const int *faces,
    size_t faceCount,
    unsigned long long userKey,
    int shadedCullBackfaces,
    struct BObolMeshLodCacheStatus *status)
{
    return mesh_lod_cache_store_mesh_impl(dbip, name, vertices, vertexCount,
	normals, faces, faceCount, userKey, shadedCullBackfaces, status,
	false);
}

int
bobol_mesh_lod_cache_store_mesh_variant(
    struct db_i *dbip,
    const char *name,
    const point_t *vertices,
    size_t vertexCount,
    const vect_t *normals,
    const int *faces,
    size_t faceCount,
    unsigned long long userKey,
    int shadedCullBackfaces,
    struct BObolMeshLodCacheStatus *status)
{
    return mesh_lod_cache_store_mesh_impl(dbip, name, vertices, vertexCount,
	normals, faces, faceCount, userKey, shadedCullBackfaces, status,
	true);
}

struct BObolMeshLod *
bobol_mesh_lod_get(struct db_i *dbip, const char *name)
{
    if (!dbip || !name)
	return NULL;

    struct BObolMeshLodContext *context = mesh_lod_context_create_for_db(dbip);
    if (!context)
	return NULL;

    struct directory *dp = db_lookup(dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL) {
	mesh_lod_context_destroy(context);
	return NULL;
    }

    unsigned long long key = mesh_lod_key_get(context, name);
    if (!key) {
	mesh_lod_context_destroy(context);
	return NULL;
    }

    struct BObolMeshLod *lod = mesh_lod_create(context, key);
    if (!lod) {
	mesh_lod_context_destroy(context);
	return NULL;
    }

    return lod;
}

struct BObolMeshLod *
bobol_mesh_lod_get_named_cached_prefix(struct db_i *dbip,
				       const char *name)
{
    if (!dbip || !name || !name[0])
	return NULL;

    struct BObolMeshLodContext *context =
	mesh_lod_context_create_for_db(dbip);
    if (!context)
	return NULL;

    const unsigned long long key = mesh_lod_key_get(context, name);
    if (!key) {
	mesh_lod_context_destroy(context);
	return NULL;
    }

    struct BObolMeshLod *lod = mesh_lod_create(context, key, true);
    if (!lod) {
	mesh_lod_context_destroy(context);
	return NULL;
    }
    return lod;
}

struct BObolMeshLod *
bobol_mesh_lod_get_cached_prefix(struct db_i *dbip,
				 unsigned long long cacheKey)
{
    if (!dbip || !cacheKey)
	return NULL;

    struct BObolMeshLodContext *context =
	mesh_lod_context_create_for_db(dbip);
    if (!context)
	return NULL;

    struct BObolMeshLod *lod =
	mesh_lod_create(context, cacheKey, true);
    if (!lod) {
	mesh_lod_context_destroy(context);
	return NULL;
    }
    return lod;
}

unsigned long long
bobol_mesh_lod_cache_key_get(const struct BObolMeshLod *lod)
{
    return lod && lod->state ? lod->state->hash : 0;
}

int
bobol_mesh_lod_load_cut(struct BObolMeshLod *lod, int cut, int reset)
{
    if (!lod || !lod->state || cut < lod->state->minPopCut ||
	cut > lod->state->maxPopCut)
	return -1;
    return mesh_lod_cut(lod, cut, reset);
}

int
bobol_mesh_lod_load_resident_cut(struct BObolMeshLod *lod, int cut,
	int reset)
{
    if (!lod || !lod->state || cut < lod->state->minPopCut ||
	cut > lod->state->maxPopCut)
	return -1;
    return mesh_lod_cut(lod, cut, reset, false);
}

int
bobol_mesh_lod_read_resident_suffix(struct BObolMeshLod *lod,
	int residentCut, int targetCut,
	BObolMeshLodSuffixCallback callback, void *callbackData)
{
    if (!lod || !lod->state || !callback)
	return 0;
    residentCut = std::max(residentCut, lod->state->minPopCut);
    targetCut = std::min(targetCut, lod->state->maxPopCut);
    return lod->state->readSuffix(
	residentCut, targetCut, callback, callbackData) ? 1 : 0;
}

int
bobol_mesh_lod_read_chunk_prefixes(struct BObolMeshLod *lod,
	const uint32_t *chunkIds, size_t chunkCount, int cut,
	BObolMeshLodChunkCallback callback, void *callbackData)
{
    if (!lod || !lod->state || !chunkIds || !chunkCount || !callback ||
	cut < lod->state->minPopCut || cut > lod->state->maxPopCut)
	return 0;
    return lod->state->readChunks(
	chunkIds, chunkCount, cut, callback, callbackData) ? 1 : 0;
}

int
bobol_mesh_lod_current_cut(const struct BObolMeshLod *lod)
{
    return (lod && lod->state) ? lod->state->currCut : -1;
}

int
bobol_mesh_lod_has_active_data(const struct BObolMeshLod *lod)
{
    return (lod && lod->faces && lod->fcnt > 0 &&
	    lod->points && lod->pcnt > 0 &&
	    lod->points_orig && lod->porig_cnt > 0) ? 1 : 0;
}

int
bobol_mesh_lod_data_get(const struct BObolMeshLod *lod,
			  struct BObolMeshLodData *data)
{
    if (!data)
	return 0;

    memset(data, 0, sizeof(*data));
    if (!lod)
	return 0;

    data->faces = lod->faces;
    data->face_count = (lod->fcnt > 0) ? static_cast<size_t>(lod->fcnt) : 0;
    data->points = lod->points;
    data->point_count = (lod->pcnt > 0) ? static_cast<size_t>(lod->pcnt) : 0;
    data->points_orig = lod->points_orig;
    data->point_orig_count = (lod->porig_cnt > 0) ?
			     static_cast<size_t>(lod->porig_cnt) : 0;
    data->normals = lod->normals;
    data->normal_count = (lod->normals && data->face_count) ?
			 data->face_count * 3 : 0;
    VMOVE(data->bmin, lod->bmin);
    VMOVE(data->bmax, lod->bmax);

    return bobol_mesh_lod_has_active_data(lod);
}

void
bobol_mesh_lod_info_init(struct BObolMeshLodInfo *info)
{
    struct BObolMeshLodInfo defaults = BOBOL_MESH_LOD_INFO_INIT;
    if (info)
	*info = defaults;
}

int
bobol_mesh_lod_info_get(const struct BObolMeshLod *lod,
			  struct BObolMeshLodInfo *info)
{
    if (!info)
	return 0;

    bobol_mesh_lod_info_init(info);
    if (!lod || !lod->state)
	return 0;

    info->active_cut = lod->state->currCut;
    info->face_count = (lod->fcnt > 0) ? static_cast<size_t>(lod->fcnt) : 0;
    info->point_count = (lod->pcnt > 0) ? static_cast<size_t>(lod->pcnt) : 0;
    info->point_orig_count = (lod->porig_cnt > 0) ?
			     static_cast<size_t>(lod->porig_cnt) : 0;
    info->normal_count = (lod->normals && info->face_count) ?
			 info->face_count * 3 : 0;
    info->has_faces = (lod->faces && info->face_count) ? 1 : 0;
    info->has_points = (lod->points && info->point_count) ? 1 : 0;
    info->has_original_points =
	(lod->points_orig && info->point_orig_count) ? 1 : 0;
    info->has_snapped_points =
	(info->has_points && info->has_original_points &&
	 lod->points != lod->points_orig) ? 1 : 0;
    info->has_normals = (lod->normals && info->normal_count) ? 1 : 0;
    info->shaded_cull_backfaces =
	(lod->state->shadedCullBackfaces &&
	 lod->state->currCut <= lod->state->maxPopCut) ? 1 : 0;
    VMOVE(info->bmin, lod->bmin);
    VMOVE(info->bmax, lod->bmax);

    return info->has_faces && info->has_points && info->has_original_points;
}

void
bobol_mesh_lod_hierarchy_info_init(
    struct BObolMeshLodHierarchyInfo *info)
{
    struct BObolMeshLodHierarchyInfo defaults =
	BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
    if (info)
	*info = defaults;
}

int
bobol_mesh_lod_hierarchy_info_get(
    const struct BObolMeshLod *lod,
    struct BObolMeshLodHierarchyInfo *info)
{
    if (!info)
	return 0;
    bobol_mesh_lod_hierarchy_info_init(info);
    if (!lod || !lod->state)
	return 0;

    lod->state->hierarchyInfo(info);
    return info->min_cut >= 0 && info->max_cut >= info->min_cut;
}

int
bobol_mesh_lod_select_cut(
    const struct BObolMeshLodHierarchyInfo *hierarchy,
    double projectedPixelDiameter,
    double targetPixelError)
{
    if (!hierarchy || hierarchy->min_cut < 0 ||
	hierarchy->max_cut < hierarchy->min_cut ||
	hierarchy->cut_count !=
	    static_cast<uint32_t>(hierarchy->max_cut + 1) ||
	!std::isfinite(projectedPixelDiameter) ||
	!std::isfinite(targetPixelError) ||
	projectedPixelDiameter <= 0.0 || targetPixelError <= 0.0)
	return -1;
    double extentSquared = 0.0;
    for (int axis = 0; axis < 3; ++axis) {
	const double extent = hierarchy->quantization_max[axis] -
	    hierarchy->quantization_min[axis];
	if (!std::isfinite(extent) || extent < 0.0)
	    return -1;
	extentSquared += extent * extent;
    }
    const double extent = std::sqrt(extentSquared);
    if (!(extent > 0.0))
	return hierarchy->max_cut;
    const double maximumObjectError =
	targetPixelError * extent / projectedPixelDiameter;
    int low = hierarchy->min_cut;
    int high = hierarchy->max_cut;
    while (low < high) {
	const int middle = low + (high - low) / 2;
	if (hierarchy->cuts[middle].object_error <= maximumObjectError)
	    high = middle;
	else
	    low = middle + 1;
    }
    return low;
}

void
bobol_mesh_lod_memshrink(struct BObolMeshLod *lod)
{
    if (!lod || !lod->state)
	return;

    if (lod->state->currCut > lod->state->maxPopCut)
	return;

    lod->state->shrinkMemory();
    lod->state->forceUpdate = true;
    mesh_lod_active_data_clear(lod);
}

size_t
bobol_mesh_lod_resident_bytes(const struct BObolMeshLod *lod)
{
    if (!lod || !lod->state)
	return 0;
    const size_t stateBytes = lod->state->residentBytes();
    return stateBytes > SIZE_MAX - sizeof(*lod) ?
	SIZE_MAX : sizeof(*lod) + stateBytes;
}

size_t
bobol_mesh_lod_resident_prefix_bytes(const struct BObolMeshLod *lod)
{
    return lod && lod->state ?
	lod->state->residentPrefixBytes() : 0;
}

void
bobol_mesh_lod_destroy(struct BObolMeshLod *lod)
{
    if (!lod)
	return;

    BObolMeshLodContext *context = lod->context;
    delete lod->state;
    lod->state = NULL;
    lod->context = NULL;
    BU_PUT(lod, struct BObolMeshLod);
    mesh_lod_context_destroy(context);
}
